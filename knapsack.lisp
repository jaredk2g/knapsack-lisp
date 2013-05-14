;;; Initialize globals
(defvar *dataset-base-dir* nil)
(defvar *num-generations* nil)
(defvar *population-size* nil)
(defvar *elitism?* nil)
(defvar *mutation-rate* nil)
(defvar *crossover-rate* nil)
(defvar *selection-operator* nil)
(defvar *crossover-operator* nil)
(defvar *mutation-operator* nil)
(defvar *tournament-size* nil)
(defvar *sa-initial-temperature* nil)
(defvar *sa-initial-iterations* nil)
(defvar *sa-perturb-operator* nil)
(defvar *sa-time-limit* nil)
(defvar *sa-alpha* nil)
(defvar *sa-beta* nil)
(defvar *hc-perturb-operator* nil)
(defvar *hc-time-limit* nil)

;;; Set globals
(setf *dataset-base-dir* "/Volumes/HDD/Users/jared/Documents/Code/LISP/knapsack/datasets/")
; GA
(setf *num-generations* 200)
(setf *population-size* 100)
(setf *elitism?* t)
(setf *mutation-rate* 0.03)
(setf *crossover-rate* 0.7)
(setf *selection-operator* "roulette")
(setf *crossover-operator* "singlepoint")
(setf *mutation-operator* "random")
(setf *tournament-size* 5)
; SA
(setf *sa-time-limit* 100)
(setf *sa-initial-temperature* 250)
(setf *sa-initial-iterations* 5)
(setf *sa-perturb-operator* "swap")
(setf *sa-alpha* 0.85)
(setf *sa-beta* 1)
; Hill Climbing
(setf *hc-time-limit* 1000)
(setf *hc-perturb-operator* "swap")

;;; Checks if a file exists
(defun file-exists-p (filename)
  (probe-file filename))

;;; Tokenizes a line into a list of numbers
(defun tokenize (str)
  (read-from-string (concatenate 'string "(" str ")")))

;;; Tokenizes a dataset file for 0/1 knapsack
;;; The dataset file has the following format:
;;; - lines beginning with # = comment
;;; - first non-comment line = knapsack capacity
;;; - remaining non-comment lines = value weight (separated by space)
;;; Returns a list of knapsack capacity and a list of objects
;;; with value and weight.
;;; The result will look like:
;;;   (capacity (value weight) (value2 weight2) ... (valueN weightN))
(defun tokenize-file (file &aux result (first-line t)
  (infile (open file :direction :input :if-does-not-exist :error)))
  (loop
    ; read each line
    (let ((temp (read-line infile nil 'end-of-file)))
        ; stop at the end of file
        (cond
          ; reached the end of file
          ((eq temp 'end-of-file)
            (close infile)
            (return))
          ; skip comments
          ((equal (subseq temp 0 1) "#") nil)
          ; first line, the number of puzzles we are supposed to discover
          (first-line
            (push (parse-integer temp) result)
            (setf first-line nil))
          ; tokenize the line
          (t
            (let ((object nil))
              (dolist (token (tokenize temp))
              (push token object))
              (push (reverse object) result))))))
  (reverse result))

(defun knapsack (filename)
  ; append base directory
  (setf filename (concatenate 'string *dataset-base-dir* filename))
  ; check if file exists
  (unless (file-exists-p filename)
      (format t "Please supply a valid filename.~%")
      (return-from knapsack))
  (format t "--- Solving 0/1 Knapsack problem ---~%")
  (format t "--- Created by Jared King (jaredtking.com) ---~%")
  ; let's get this party started
  (let* ((knapsack (tokenize-file filename))
    (capacity (car knapsack)) (weights (cdr knapsack)))
      (genetic-algorithm capacity weights)
      (simulated-annealing capacity weights)
      (hill-climbing capacity weights)
      ))

;;;
;;;
;;; Genetic Algorithm
;;;
;;; Each chromosome is represented by a list
;;; with a '0' or '1' corresponding with each item
;;; from `weights`. '0' = not present, '1' = present
;;;
;;; The fitness is calculated by the following format:
;;;

(defun genetic-algorithm (capacity weights)
  (format t "~%--- Running genetic algorithm~%")
  (format t "--- Knapsack has capacity of ~a~%" capacity)
  (format t "--- Using ~a items~&" (length weights))
  (format t "--- Creating initial population pool~%")
  (let ((population (generate-initial-population capacity weights *population-size*))
    (past-runs nil) (i 0) (converged nil))
    (format t "--- Building generations")
    ; stop after N generations or if the population converges
    (loop while (and (< i *num-generations*) t) do
      (incf i)
      ; Reproduce
      (let ((new-population nil))
        ; carry over the best chromosomes when elitism enabled
        (when *elitism?*
          (setf new-population (select-best 2 population t)))
        ; produce the remaining children
        (dotimes (i (/ (- *population-size* (length new-population)) 2))
          ; Selection
          (let* ((selection (select population *selection-operator*))
            (a (car selection)) (b (cadr selection)))
            ; Crossover
            (dolist (chromosome (crossover a b *crossover-operator*))
              ; Mutate
              (setf chromosome (mutate chromosome *mutation-operator*))
              ; add to next generation
              (push (list chromosome (fitness capacity weights chromosome)) new-population))))
        ; replace old population with new generation
        (setf population new-population)
        ; calculate the variance between fitness values
        (push (mean-variance (mapcar #'(lambda (chromosome) (caadr chromosome)) population)) past-runs)
        (when (>= i 5)
          ; checks for convergence of the population
          (let* ((recent-n-runs (subseq past-runs 0 5)) (base-run (car recent-n-runs)))
            (mapcar #'(lambda(run)
              (when (and (< (abs (- (car base-run) (car run))) 1)
                  (< (abs (- (cadr base-run) (cadr run))) 1))
                (setf converged t))
              ) (cdr recent-n-runs))))))
    ; find the most fit feasible chromosome
    (let ((most-fit (car (select-best 1 population nil))))
      (format t "~%--- Result from genetic algorithm:~%")
      (format t "    Generations: ~a~%" i)
      (format t "    Value: ~a~%" (nth 1 (cadr most-fit)))
      (format t "    Weight: ~a~%" (nth 2 (cadr most-fit)))
      (format t "    Chromosome: ~a~%" (car most-fit))
      most-fit)))

;;; Generates a random population given the weights.
;;; Returns a population of N chromosomes.
(defun generate-initial-population (capacity weights population-size)
  (let ((population nil))
    (dotimes (i population-size)
      (push (generate-random-chromosome capacity weights) population))
    population))
  
;;; Selects the most fit, feasible N chromosomes
;;; The result looks like:
;;; ((chromosome-1 fitness) (chromosome-2 fitness) ... (chromosome-N fitness))
(defun select-best (n population allow-infeasibles?)
  (subseq (sort (copy-list population) #'(lambda(a b)
    (and (or allow-infeasibles? (nth 3 (cadr a))) (> (caadr a) (caadr b))))) 0 n))

;;; Selects 2 chromosomes from the population, using the
;;; specified selection technique.
;;; The result looks like: ((chromosome-1 fitness) (chromosome-2 fitness))
(defun select (population method)
  ; seed the random number generator
  (setf *random-state* (make-random-state t))  
  (cond
    ; Tournament Selection
    ((equal method "tournament")
      ; select N chromosomes from population at random
      (let* ((temp-population (copy-list population)) (tournament-population nil))
        (dotimes (i *tournament-size*)
          (let ((index (random (length temp-population))))
            (push (nth index temp-population) tournament-population)
            (setf temp-population (remove (nth index temp-population) temp-population))))
        ; return best 2
        (return-from select (select-best 2 tournament-population t))))
    ; Roulette Selection
    ((equal method "roulette")
      (let ((roulette nil) (max 0))
        ; generate the bounds for each chromosome
        (mapcar #'(lambda(chromosome)
          (push (incf max (caadr chromosome)) roulette)) population)
        (setf roulette (reverse roulette))
        (let ((selected nil))
          (dotimes (i 2)
            ; generate a random number between 0 and max
            (let ((rand (random max)) (index 0))
              ; loop through each fitness until rand is < 0
              (loop while (> rand 0) do
                (decf rand (nth index roulette))
                (incf index))
              (push (nth index population) selected)))
          (return-from select selected))))))

;;; Performs the specified crossover method on 2 chromosomes.
;;; Returns a list of 2 chromosomes.
(defun crossover (a b method)
  ; seed the random number generator
  (setf *random-state* (make-random-state t))  
  (cond
    ; Uniform Crossover
    ((equal method "uniform")
      (let ((child-a nil) (child-b nil))
        (dotimes (i (length a))
          (let ((rand (random 1)))
            (push (if (eq rand 1) (nth i child-a) (nth i child-b)) child-b)
            (push (if (eq rand 0) (nth i child-a) (nth i child-b)) child-b)))
        (return-from crossover (list child-a child-b))))
    ; Single-Point Crossover
    ((equal method "singlepoint")
      (return-from crossover (if (> (random 1.0) *crossover-rate*) (list a b)
        (let ((n (random (length a))))
          (list
            (append (subseq a 0 n) (nthcdr n b))
            (append (subseq b 0 n) (nthcdr n a)))))))))

;;; Mutates a chromosome according to the specified method.
;;; Returns a chromosome.
(defun mutate (chromosome method)
  ; seed the random number generator
  (setf *random-state* (make-random-state t))
  (cond
    ; Random Bit
    ((equal method "random")
      (mapcar #'(lambda (x)
           (if (> (random 1.0) *mutation-rate*) x (random 1))) chromosome)
       (return-from mutate chromosome))
    ; Flip Bit
    ((equal method "flip")
      (mapcar #'(lambda (x)
           (if (> (random 1.0) *mutation-rate*) x (if (= x 1) 0 1))) chromosome)
       (return-from mutate chromosome))))

;;; Calculates the fitness of a given capacity, weights, and chromosome
;;; The result looks like:
;;; (fitness value weight feasible?)
(defun fitness (capacity weights chromosome)
  (let ((value 0) (weight 0))
    ; calculate the value and weight of chromosome
    (dotimes (index (length weights))
      (when (eq (nth index chromosome) 1)
        (let ((item (nth index weights)))
          (setf weight (+ weight (cadr item)))
          (setf value (+ value (car item))))))
;    (format t "weight: ~a value: ~a~%" weight value)
    (let ((fitness value))
      ; total weight > capacity = infeasible chromosome
      (when (> weight capacity)
        (decf fitness (- weight capacity)))
      (setf fitness (max 1 fitness))
      (list fitness value weight (<= weight capacity)))))

;;;
;;;
;;; Simulated Annealing
;;;
;;;

(defun simulated-annealing (capacity weights)
  ; seed the random number generator
  (setf *random-state* (make-random-state t))
  (format t "~%--- Running simulated annealing~%")
  (format t "--- Knapsack has capacity of ~a~%" capacity)
  (format t "--- Using ~a items~&" (length weights))
  ; initialize solution, temperature, and iterations
  (let* ((init-chromosome (generate-random-chromosome capacity weights))
      (solution (car init-chromosome)) (solution-fitness (cadr init-chromosome))
      (temperature *sa-initial-temperature*) (iterations *sa-initial-iterations*)
      (time 0) (i 0) (num-pertubations 0))
    ; stop when we reach the time limit
    (loop while (< time *sa-time-limit*) do
      ; loop 'iterations' times
      (loop while (< i iterations) do 
        ; newS = perturb(S)
        (let* ((new-solution (perturb capacity weights solution *sa-perturb-operator*))
          (new-solution-fitness (fitness capacity weights new-solution)))
          (incf num-pertubations)
          ; when h(newS) < h(S) or random < e^((h(S) - h(newS))/T)
          (when (or (< (car new-solution-fitness) (car solution-fitness))
              (< (random 1.0) (exp (/ (- (car solution-fitness) (car new-solution-fitness)) temperature))))
            ; accept new solution
            (setf solution new-solution)
            (setf solution-fitness new-solution-fitness)))
        (incf i))
        ; T = alpha * T
        (setf temperature (* *sa-alpha* temperature))
        ; iterations = beta * iterations
        (setf iterations (* *sa-beta* iterations))
        ; reset inner loop
        (setf i 0)
        ; increment time
        (incf time))
    (format t "~%--- Result from simulated annealing:~%")
    (format t "    Pertubations: ~a~%" num-pertubations)
    (format t "    Value: ~a~%" (nth 1 solution-fitness))
    (format t "    Weight: ~a~%" (nth 2 solution-fitness))
    (format t "    Chromosome: ~a~%" solution)))

;;;
;;;
;;; Hill Climbing
;;;
;;;

(defun hill-climbing (capacity weights)
  ; seed the random number generator
  (setf *random-state* (make-random-state t))
  (format t "~%--- Running hill climbing~%")
  (format t "--- Knapsack has capacity of ~a~%" capacity)  
  (format t "--- Using ~a items~&" (length weights))
; initialize solution, temperature, and iterations
  (let* ((init-chromosome (generate-random-chromosome capacity weights))
      (solution (car init-chromosome)) (solution-fitness (cadr init-chromosome))
      (time 0) (num-pertubations 0))
    ; stop when we reach the time limit
    (loop while (< time *hc-time-limit*) do
      ; newS = perturb(S)
      (let* ((new-solution (perturb capacity weights solution *hc-perturb-operator*))
        (new-solution-fitness (fitness capacity weights new-solution)))
        (incf num-pertubations)
        ; when h(newS) < h(S)
        (when (< (car new-solution-fitness) (car solution-fitness))
          ; accept new solution
          (setf solution new-solution)
          (setf solution-fitness new-solution-fitness)))
      ; increment time
      (incf time))
    (format t "~%--- Result from hill climbing:~%")
    (format t "    Pertubations: ~a~%" num-pertubations)
    (format t "    Value: ~a~%" (nth 1 solution-fitness))
    (format t "    Weight: ~a~%" (nth 2 solution-fitness))
    (format t "    Chromosome: ~a~%" solution)))

;;; Perturbation operator for a given solution.
;;; Returns a feasible solution.
(defun perturb (capacity weights solution method)
  ; seed the random number generator
  (setf *random-state* (make-random-state t))
  (let ((new-solution (copy-list solution)))
    (cond
      ; Best Fit
      ((equal method "bestfit")
        (let ((weight 0))
          ; calculate the weight of current solution
          (dotimes (index (length weights))
            (when (eq (nth index new-solution) 1)
              (let ((item (nth index weights)))
                (setf weight (+ weight (cadr item)))
                )))
          ; sort the items according to weight and value
          (let ((items (sort (iota (length weights)) #'(lambda(a b) 
              (let ((item-a (nth a weights)) (item-b (nth b weights)))
                (or (and (eq (cadr item-a) (cadr item-b)) (< (car item-a) (car item-b)))
                  (> (cadr item-a) (cadr item-b))))))))
            (dolist (index items)
              ; add the item if it fits
              (when (<= (+ weight (cadr (nth index weights))) capacity)
                (setf (nth index new-solution) 1))))))
      ; Swap
      ((equal method "swap")
        ; swap out an object for one of similar size
        (let ((index (random (length solution))) (item nil)
            (i 0) (equal-min nil) (equal-index 0))
          ; pick a random item in knapsack
          (loop while (eq 0 (nth index solution)) do
            (setf index (random (length solution))))
          (setf item (nth index weights))
          (dolist (test weights)
            (unless (and (eq (car test) (car item)) (eq (cadr test) (cadr item)))
              (let ((delta (abs (- (car item) (car test)))))
                ; find an item of closest to equal value
                (when (or (not equal-min) (< delta equal-min))
                  (setf equal-min delta)
                  (setf equal-index i))))
            (incf i))
          ; perform swap
          (setf (nth index new-solution) 0)
          (setf (nth equal-index new-solution) 1))))
      ; do not allow infeasibles
      (make-feasible capacity weights new-solution)))

;;;
;;; Utility Functions
;;;

;;; Generates a random chromosome given the weights.
;;; Uses a greedy algorithm to generate more optimal chromosomes
;;; than pure randomness.
;;; The return looks like: (chromosome fitness)
(defun generate-random-chromosome (capacity weights)
  ; seed the random number generator
  (setf *random-state* (make-random-state t))
  (let* ((length (length weights))
    (chromosome (make-list length :initial-element 0))
    (weight 0) (tries nil))
    ; loop until tries = length:
    (loop while (and (< (length tries) length) (not (eq weight capacity))) do
      ; pick a random number between 0 and length
      (let* ((index (random length)) (item (nth index weights)))
        (unless (find index tries)
          ; if adding this item does not go over our capacity, add it
          (unless (> (+ weight (cadr item)) capacity)
            ; set the bit at that index to 1
            (setf (nth index chromosome) 1)
            ; add the result to the weight
            (setf weight (+ weight (cadr item))))
          ; document the index that was attempted
          (push index tries))))
    (list chromosome (fitness capacity weights chromosome))))

;;; Makes a chromosome feasible
;;; Returns a chromosom
(defun make-feasible (capacity weights chromosome)
  (let ((value 0) (weight 0) (new-solution (copy-list chromosome)))
    ; calculate the value and weight of chromosome
    (dotimes (index (length weights))
      (when (eq (nth index chromosome) 1)
        (let ((item (nth index weights)))
          (setf weight (+ weight (cadr item)))
          (setf value (+ value (car item))))))
    ; sort the items according to weight and value
    (let ((items (sort (iota (length weights)) #'(lambda(a b) 
        (let ((item-a (nth a weights)) (item-b (nth b weights)))
          (or (and (eq (cadr item-a) (cadr item-b)) (< (car item-a) (car item-b)))
            (< (cadr item-a) (cadr item-b))))))))
      ; while infeasible, take out the smallest objects until feasible
      (loop while (> weight capacity) do
        (let ((index (pop items)))
          (setf (nth index new-solution) 0)
          (setf weight (- weight (cadr (nth index weights)))))))
    new-solution))

; props to http://nklein.com/2011/02/calculating-the-mean-and-variance-with-one-pass/
; returns (mean variance)
(defun mean-variance (data)
  (flet ((square (x) (* x x)))
    (destructuring-bind (n xs x2s)
        (reduce #'(lambda (accum xi)
                    (list (1+ (first accum))
                          (+ (second accum) xi)
                          (+ (third accum) (square xi))))
                data :initial-value '(0 0 0))
      (let ((mu (/ xs n)))
        (list mu (- (/ x2s n) (square mu)))))))

;;; Generates a list of natural numbers from 0...N-1
(defun iota (stop)
  (let ((numbers nil))
    (dotimes (i stop) (push i numbers))
    (reverse numbers)))

;;; Deletes the nth element of a list
(defun delete-nth (n list)
  (if (zerop n)
    (cdr list)
    (let ((cons (nthcdr (1- n) list)))
      (if cons
        (setf (cdr cons) (cddr cons))
        cons))))