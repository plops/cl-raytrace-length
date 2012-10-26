#.(require :raytrace)

(in-package :raytrace)

(defun focal-length-thin-lens (c1 c2 n)
  (/
   (* (- n 1) (- c1 c2))))

(defun focal-length-thick-lens-in-air (c1 c2 n d)
 (/
  (+
   (* (- n 1) (- c1 c2))
   (* (expt (- n 1) 2) d c1 c2 (/ n)))))

#+nil
(focal-length-thick-lens-in-air 100 -100 1.515 12)

(defun accum-thickness (sys)
  (let ((start 0)
	(thick-old 0))
    (loop for (curv ind thick) in sys and i from 0 collect
	 (prog1
	     (list curv ind (incf start thick-old))
	   (setf thick-old thick)))))

#+nil
(unsplit-optical-system `((0 50 -30 0) (1 1.5 1) (100 10 100)))
#+nil
(accum-thickness (unsplit-optical-system `((0 50 -30 0) (1 1.5 1) (100 10 100))))

(defun calc-centers (sys)
  (loop for (curv ind vert) in sys collect
       (list curv ind (if (= curv 0)
			  vert
			  (+ vert (/ curv))))))


#+nil
(intersect-sphere (v 0 10 0) (v 0 0 1) (v 0 0 200) 100d0)

#+nil
(let* ((center (v 0 0 200))
       (radius 100d0)
       (i (intersect-sphere (v 0 10 0) (v 0 0 1) (v 0 0 200) 100d0)))
  (refract-plane (v 0 0 1)
		 (normalize (.- i center))
		 (/ 1d0 1.5)))


(defun refract-sphere (intersect dir center n1/n2 radius)
  (refract-plane dir
		 (.s (signum radius) (normalize (.- intersect center)))
		 n1/n2))

#+nil
(let* ((start (v 10 0 0))
       (dir (v 0 0 1))
       (center (v 0 0 200))
       (radius 100d0)
       (n 1.5d0)
       (i (intersect-sphere start dir
			    center radius)))
  (list i
   (refract-sphere i dir center (/ 1 n))))

(defun trace-ray (start dir center radius n1/n2)
  (let* ((i (intersect-sphere start dir
			      center radius)))
    (values i
	    (refract-sphere i dir center n1/n2 radius))))

(defun principal-plane (sys &optional (h 0.0001))
  (let* ((start (v h 0 0))
	 (dir (v 0 0 1))
	 (old-index 1d0))
    (loop for (curv index center-z) in
	 (calc-centers (accum-thickness sys))
       collect
	 (prog1
	     (multiple-value-bind (start-new dir-new)
		 (trace-ray start dir (v 0 0 center-z)
			    (/ curv)
			    (/ old-index index))
	       (setf dir dir-new
		     start start-new))
	   (setf old-index index)))
    (let* ((xd (aref dir 0))
	   (zd (aref dir 2))
	   (xt (aref start 0))
	   (zt (aref start 2))
	   (alpha (/ (- h xt) xd))
	   (beta (/ (- xt) xd))
	   (z0 (+ zt (* beta zd)))
	   (zp (+ zt (* alpha zd))))
      (values (- z0 zp) zp))))


(defun split-optical-system (sys)
  (let ((curvs (loop for (curv ind thick) in sys collect curv))
	(inds (loop for (curv ind thick) in sys collect ind))
	(thicks (loop for (curv ind thick) in sys collect thick)))
    (list curvs
	  (subseq inds 0 (1- (length curvs)))
	  (subseq thicks 0 (1- (length curvs))))))

(defun unsplit-optical-system (sys3 &key (outside-index 1d0))
  (destructuring-bind (curvs inds thicks) sys3
    (unless (= (1- (length curvs))
	       (length inds)
	       (length thicks))
      (error "system doesn't contain a correct number of parameters"))
    (loop for i below (length curvs) collect
	 (if (< i (length inds))
	  (list (elt curvs i) (elt inds i) (elt thicks i))
	  (list (elt curvs i) (elt inds (1- i)) 0d0)))))

#+nil
(unsplit-optical-system `((0 50 -30 0) (1 1.5 1) (100 10 100)))

#+nil
(unsplit-optical-system ;; zeiss 40x 0..2mm coverslip
 `(,(mapcar #'/ '(-3.566 -3.839 -659.4 -10.96 48.448 7.515 -8.485 -24.36
			17.013 -25.925 8.593 -29.03 5.975))
    (1.5205 1.5891 1.6134 1.4866 1.8807 1.5285 1.4585 1.6967)
    (3.24 1.43 2.50 0.62 1.0 4.7 1.0 0.2 3.0 15.75 4.0 1.5)))


(defun reverse-optical-system (sys)
  (let ((s3 (split-optical-system sys)))
    (destructuring-bind (c i d) s3
     (unsplit-optical-system (list (mapcar #'- (reverse c))
				   (reverse i)
				   (reverse d))))))
(defun bend-lens (f n d x)
  ;;Solve[(n-1)*(c1-c2)+(n-1)^2/n*d*c1*c2=1/f,c2=(x-1)/(x+1)*c1,c2]
  (if (= x 1)
      (values (/ (* f (- n 1)))
	      0d0)
      (if (= d 0)
	  (values (/ (+ 1 x)
		     (* 2 f (- n 1)))
		  (/ (- x 1)
		     (* 2 f (- n 1))))
	  (if (= d (* f n))
	      (values (/ (* f (- n 1)))
		      (/ (- x 1)
			 (* f (- n 1) (+ x 1))))
	      (let* ((sn (expt (- n 1) 2))
		     (discriminant (* f sn n (+ (* d x x) (* f n) (- d))))
		     (sd (sqrt discriminant))
		     (div1 (+ sd (* -1 f n n x) (* f n x)))
		     (div2 (+ (- sd) (* -1 f n n x) (* f n x))))
		(if (/= 0 div1)
		    (values
		     (/ (+ sd (* f n) (* -1 n n f))
			(* d f sn (- x 1)))
		     (- (/ (* n (+ (- sd) (* d n x) (* -1 d n) (* -1 d x) d (* f n n) (* -1 f n)))
			   (* d (- n 1) (- div1)))))
		    (if (/= 0 div2)
			(values
			 (/ (+ (- sd) (* f n) (* -1 n n f))
			    (* d f sn (- x 1)))
			 (- (/ (* n (+ sd (* d n x) (* -1 d n) (* -1 d x) d (* f n n) (* -1 f n)))
			       (* d (- n 1) (- div2))))))))))))
#+nil
(bend-lens 100d0 1.5d0 10d0 1d0)

(defvar *bla* nil)


(defvar *img* nil)

(defun .minmax (img)
  (let* ((i1 (sb-ext:array-storage-vector img))
	 (mi (reduce #'min i1))
	 (ma (reduce #'max i1))
	 (n (length i1))
	 (avg (loop for i below n sum
		(* (aref i1 i) (/ 1d0 n)))))
    (list :avg-mm avg :min-nm (* 1e9 (- avg mi)) :ma-nm (* 1e9 (- ma avg)))))

(defun .norm (img)
  (let* ((a (make-array (array-dimensions img)
			:element-type 'double-float))
	 (a1 (sb-ext:array-storage-vector a))
	 (i1 (sb-ext:array-storage-vector img))
	 (mi (reduce #'min i1))
	 (ma (reduce #'max i1)))
    (dotimes (i (length a1))
      (setf (aref a1 i) (/ (- (aref i1 i) mi)
			   (- ma mi))))
    a))

(defun .linear (a2)
  (sb-ext:array-storage-vector a2))

(defun .floor (a)
  (let* ((a1 (.linear a))
         (b (make-array (array-dimensions a) :element-type 'fixnum))
         (b1 (.linear b)))
    (dotimes (i (length a1))
      (setf (aref b1 i) (floor (aref a1 i))))
    b))
(defun .scale (val a)
  (let* ((a1 (.linear a))
         (b (make-array (array-dimensions a) :element-type 'double-float))
         (b1 (.linear b)))
    (dotimes (i (length a1))
      (setf (aref b1 i) (* val (aref a1 i))))
    b))

(defun .add (val a)
  (let* ((a1 (.linear a))
         (b (make-array (array-dimensions a) :element-type 'double-float))
         (b1 (.linear b)))
    (dotimes (i (length a1))
      (setf (aref b1 i) (+ val (aref a1 i))))
    b))

(defun .ub8 (a)
  (let* ((a1 (.linear a))
         (b (make-array (array-dimensions a) :element-type '(unsigned-byte 8)))
         (b1 (.linear b)))
    (dotimes (i (length a1))
      (setf (aref b1 i) (min 255 (max 0 (round (aref a1 i))))))
    b))

(defun write-pgm (filename img)
  (declare (simple-string filename)
           ((or (array (unsigned-byte 16) 2)
                (array (unsigned-byte 8) 2)) img)
           (values null &optional))
  (let ((ty (array-element-type img)))
   (destructuring-bind (h w)
       (array-dimensions img)
     (declare ((integer 0 65535) w h))
     (with-open-file (s filename
                        :direction :output
                        :if-exists :supersede
                        :if-does-not-exist :create)
       (declare (stream s))
       (if (equal '(unsigned-byte 16) ty)
           (format s "P5~%~D ~D~%65535~%" w h)
           (if (equal '(unsigned-byte 8) ty)
               (format s "P5~%~D ~D~%255~%" w h)
               (error "unsupported type"))))
     (with-open-file (s filename 
                        :element-type ty
                        :direction :output
                        :if-exists :append)
       (let ((data-1d (make-array 
                       (* h w)
                       :element-type ty
                       :displaced-to img)))
         (write-sequence data-1d s)))
     nil)))


#+nil
(let ((c1 (/ 100d0))
      (n 1.5d0)
      (d 10d0))
  (multiple-value-bind (f s)
      (principal-plane `((,c1 ,n ,d) (,(- c1) 1 0)))
    (- (* 2 f) s)))

(defvar *sys* nil)


(let* ((d 10d0)
       (f 100d0)
       (n 1.5d0))
  (multiple-value-bind (c1 c2) (bend-lens f n d .0d0)
   ;; curvature index-behind thickness-behind
   (let ((lens `((,c1 ,n ,d) (,c2 1 0))))
     (multiple-value-bind (f sa)
	 (principal-plane lens)
       (multiple-value-bind (f2 sb)
	   (principal-plane (reverse-optical-system lens))
	 (format t "f=~f mm, df=~f nm~%" f (* 1e9 (- f f2)))
	 (let ((sys `((0 1 ,(- (* 2 f) (- d sa)))
		      ,(list c1 n d)
		      ,(list c2 1 (- (* 2 f) (- d sb)))
		      (0 1 0))))
	   (defparameter *sys* sys)
	   (defparameter *bla*
	     (let* ((beam-radius .6d0)
		    (ap (/ beam-radius (* 2 f)))
		    (num 120))
	       (loop for dir-y from (- ap) upto ap by (/ ap num) collect
		    (loop for dir-x from (- ap) upto ap by (/ ap num) collect
			 (let ((start (v 0 0 0))
			       (dir (v dir-x dir-y (sqrt (- 1 (+ (expt dir-x 2)
								 (expt dir-y 2))))))
			       (old-index 1d0))
			   (loop for (curv index center-z) in
				(calc-centers (accum-thickness sys))
			      collect
			   
				(prog1
				    (if (= curv 0)
					(intersect-plane start dir (v 0 0 center-z)
							 (v 0 0 -1))
					(multiple-value-bind (start-new dir-new)
					    (trace-ray start dir (v 0 0 center-z)
						       (/ curv)
						       (/ old-index index))
					  (setf dir dir-new
						start start-new)))
				  (setf old-index index))))))))
	   (with-open-file (s "/dev/shm/o.dat" :direction :output
			      :if-does-not-exist :create
			      :if-exists :supersede)
	     (dolist (e (nth
			 (floor (length *bla*) 2)
			 *bla*))
	       (dolist (p e)
		 (format s "~f ~f~%" (aref p 2) (aref p 0)))
	       (terpri s)))
	   (with-open-file (s "/dev/shm/spot.dat" :direction :output
			      :if-does-not-exist :create
			      :if-exists :supersede)
	     (dolist (e *bla*)
	       (dolist (f e)
		 (let ((p (first (last f))))
		
		   (format s "~f ~f~%" (aref p 0) (aref p 1))))
	       (terpri s)))
	   (let* ((w (length *bla*))
       (h (length (nth 0 *bla*)))
       (img (make-array (list h w) :element-type 'double-float)))
  (let ((inds (second (split-optical-system *sys*))))
    (dotimes (i w)
      (dotimes (j h)
	(let* ((points (nth i (nth j *bla*)))
	       (start (first points)))
	  (setf (aref img i j)
		(loop for n in inds and p in (rest points) sum
		     (prog1
			 (* n (norm (.- start p)))
		       (setf start p)))
		))))

    (format t "minmax: ~a~%" (.minmax img))
    (defparameter *img* img)
    (write-pgm "/dev/shm/phase.pgm" (.ub8 (.scale 255 (.norm *img*))))))))))))


#+nil 
(* 180 (/ pi) (cos (
asin (/ 1.25 1.525)))) ;; aperture half-angle

#+nil
(let ((s *zeiss100*))
  (calc-centers
   (accum-thickness s)))


#.(require :cl-who)

(accum-thickness *zeiss100*)

(calc-centers (accum-thickness *zeiss100*))

(calc-centers
 (accum-thickness (reverse-optical-system *zeiss100*)))

(defun trace-sys (sys start dir &optional (old-index 1d0))
  (loop for (curv index center-z) in (calc-centers (accum-thickness sys)) collect	   
       (prog1
	   (if (= curv 0)
	       (let* ((n (v 0 0 -1))
		      (start-new (intersect-plane start dir (v 0 0 center-z)
						  n))
		      (dir-new (refract-plane dir n (/ old-index index))))
		 (setf dir dir-new
		       start start-new))
	       (multiple-value-bind (start-new dir-new)
		   (trace-ray start dir (v 0 0 center-z)
			      (/ curv)
			      (/ old-index index))
		 (setf dir dir-new
		       start start-new)))
	 (setf old-index index))))

#+nil
(progn

 #+nil (defparameter *zeiss100* 
    (unsplit-optical-system ;; zeiss US2009/0284841 cheap planachromat
     ;; tubelength 200, 100x NA1.25
     ;; visual field factor 20, .28 mm working distance
     ;; curv index thick
     `(,(append '(0d0 0d0 0d0)
		(mapcar #'/ '(-1.9845 -32.1451 -5.8498 -87.1808 10.0650
			      -16.4805 16.4805 -10.0650 87.1808 8.8700
			      4.4577 2.4528 -2.4528 -4.4577))
		'(0d0))
	,(mapcar #'(lambda (x) (+ 1d0 (/ x 1000d0))) 
		 '(525 518 517 0 758 0 762 667 0 667 762 0 489 813 0 813 0 0))
	(.17 .281 2.770 .2 2.2 7.127 3. 4. .2 4. 3. 4.823 6.5 4. 1.8 4. .5 100))))

  (defparameter *zeiss100* 
    (unsplit-optical-system ;; curv index thick
     `(,(append '(0d0)
		(mapcar #'/ '(20d0))
		'(0d0 0d0))
	,(mapcar #'(lambda (x) (+ 1d0 (/ x 1000d0))) 
		 '(0 500 0))
	(100 3 100))))


  (defparameter *bla* nil)
  (defparameter *bla2* nil)
  (defparameter *bla3* nil)

  (let ((sys *zeiss100*))
    (defparameter *sys* sys)
    
    (defparameter *bla* ;; margin ray
      (let* ((start (v 0 5d0 0.0d0))
	     (dir-x 0d0)
	     (dir-y 0d0 ;(/ 1.25d0 1.525)
	       )
	     (dir (v dir-x dir-y (sqrt (- 1 (+ (expt dir-x 2)
					       (expt dir-y 2))))))
	     (old-index 1d0))
	(loop for (curv index center-z) in
	     (calc-centers (accum-thickness sys))
	   collect	   
	     (prog1
		 (if (= curv 0)
		     (intersect-plane start dir (v 0 0 center-z)
				      (v 0 0 -1))
		     (multiple-value-bind (start-new dir-new)
			 (trace-ray start dir (v 0 0 center-z)
				    (/ curv)
				    (/ old-index index))
		       (setf dir dir-new
			     start start-new)))
	       (setf old-index index)))))

    #+nil
    (defparameter *bla3* ;; coma ray
      (let* ((start (v 0 .1596d0 0.0d0))
	     (dir-x 0d0)
	     (dir-y (/ -.5d0 1.525)
	       )
	     (dir (v dir-x dir-y (sqrt (- 1 (+ (expt dir-x 2)
					       (expt dir-y 2))))))
	     (old-index 1.525d0))
	(loop for (curv index center-z) in
	     (calc-centers (accum-thickness sys))
	   collect	   
	     (prog1
		 (if (= curv 0)
		     (intersect-plane start dir (v 0 0 center-z)
				      (v 0 0 -1))
		     (multiple-value-bind (start-new dir-new)
			 (trace-ray start dir (v 0 0 center-z)
				    (/ curv)
				    (/ old-index index))
		       (setf dir dir-new
			     start start-new)))
	       (setf old-index index)))))

    
    (defparameter *bla2* ;; chief ray
      (let* ((start (v 0 -2d0 0))
	     (dir-x 0d0)
	     (dir-y 0d0 ;(/ 10.d0 200d0) ;; y=12.5 mm in intermediate is 202um in sample
	       )
	     (dir (v dir-x dir-y (sqrt (- 1 (+ (expt dir-x 2)
					       (expt dir-y 2))))))
	     (syss (reverse-optical-system *zeiss100*))
	     (old-index (second (first syss))))
	(defparameter *sysss* syss)
	(trace-sys syss start dir old-index)))
    #+nil
    (with-open-file (s "/dev/shm/o.dat" :direction :output
		       :if-does-not-exist :create
		       :if-exists :supersede)
      (dolist (e (nth
		  (floor (length *bla*) 2)
		  *bla*))
	(dolist (p e)
	  (format s "~f ~f~%" (aref p 2) (aref p 0)))
	(terpri s))))





  (with-open-file (st "/dev/shm/o.svg" :direction :output
		      :if-does-not-exist :create
		      :if-exists :supersede)
    (format st "~a" "<?xml version='1.0' encoding='iso-8859-1'?>")
    (cl-who:with-html-output (s st :indent t)
      (cl-who:htm 
       (:svg :width "5000px" :height "5000px"
	     :version "1.1" :xmlns "http://www.w3.org/2000/svg" :|xmlns:xlink| "http://www.w3.org/1999/xlink" 
	     :|viewBox| "0 0 1200 600" ;:style "enable-background new 0 0 1200 600;" :|xml:space| "preserve"
	     (cl-who:str (loop for  (curv index center-z) in
			      (calc-centers (accum-thickness  *zeiss100*))
			    do
			      (let* ((sc 12)
				     (h (* sc 7)))
				(let ((xe (+ 20 (* sc 48.571)))
				      (rpupil (let* ((m 100) ;; pupil radius is 2.5
						     (ftl 200)
						     (f (/ ftl m))
						     (na 1.25))
						(* f na))))
				  (cl-who:htm (:line :x1 20 :y1 200 ;; optical axis
						     :x2 5000 :y2 200
						     :style "stroke:black; stroke-width:.1px;")
					      (:line :x1 100 :y1 (+ 200 (* sc 2.5)) ;; axis at 2.5mm
						     :x2 300 :y2 (+ 200 (* sc 2.5))
						     :style "stroke:black; stroke-width:.1px;")
					      (:line :x1 xe :y1 (+ 200 (* sc h)) ;; top stop
						     :x2 xe :y2 (+ 200 (* sc rpupil))
						     :style "stroke:black; stroke-width:2px;")
					      (:line :x1 xe :y1 (- 200 (* sc h)) ;; bottom stop
						     :x2 xe :y2 (- 200 (* sc rpupil))
						     :style "stroke:black; stroke-width:2px;")))
				
				(loop for (color f) in `(("green" ,*bla*) ("blue" ,*bla3*)) do
				     (when f
				       (loop for i from 1 below (length f) and e in f do
					    (let ((x1 (format nil "~a" (+ 20 (* sc (aref e 2)))))
						  (y1 (format nil "~a" (+ 200 (* sc (aref e 1)))))
						  (x2 (format nil "~a" (+ 20 (* sc (aref (elt f i) 2)))))
						 (y2 (format nil "~a" (+ 200 (* sc (aref (elt f i) 1))))))
					     (cl-who:htm (:line :x1 x1 :y1 y1
								:x2 x2 :y2 y2
								:style (format nil "stroke:~a;stroke-width:.1px"
									       color)))))))
				
				(when *bla2*
				 (loop for i from 1 below (length *bla2*) and e in *bla2* do
				      (let* ((len (aref (first (last *bla2*)) 2))
					     (x1 (format nil "~a" (+ 20 (* sc (- len (aref e 2))))))
					     (y1 (format nil "~a" (+ 200 (* sc (aref e 1)))))
					     (x2 (format nil "~a" (+ 20 (* sc (- len (aref (elt *bla2* i) 2))))))
					    (y2 (format nil "~a" (+ 200 (* sc (aref (elt *bla2* i) 1))))))
					(cl-who:htm (:line :x1 x1 :y1 y1
							   :x2 x2 :y2 y2
							   :style "stroke:red;stroke-width:.1px")))))
				(if (= 0 curv)
				    (let ((x (format nil "~a" (+ 20 (* sc center-z)))))
				      (cl-who:htm (:line :x1 x :y1 (format nil "~a" (- 200 h))
							 :x2 x :y2 (format nil "~a" (+ 200 h))
							 :style "stroke:black;")))
				    
				    
				    (cl-who:htm 
				     
				     #+nil (:circle :cx (+ 200 center-z)
						    :cy 200
						    :r (/ (abs curv))
						    :style "stroke: red; stroke-width: 3; fill:none;")
				     (:g :stroke "black" :fill "none" 
					 (:path :d (let* ((r (* sc (-  (/ curv))))
							  (cz (* sc center-z))
							  ;; (alpha (* (- 180 20) pi (/ 180)))
							  (alpha (if (< h (abs r))
								 (- pi (sin (/ h (abs r))))
								 (* (- 180 90) pi (/ 180))
								 ))
							  (x (- (+ r (* r (cos alpha)))))
							  (y (* -1 r (sin alpha))))
						     (concatenate 'string
								  ;; (rx ry x-axis-rotation large-arc-flag sweep-flag x y)+
								  (format nil "M ~a 200 a ~a ~a 0 0 0 ~a ~a"
									  (+ 20 r cz)
									  r r
									  x y)
								  (format nil "M ~a 200 a ~a ~a 0 0 1 ~a ~a"
									  (+ 20 r cz)
									  r r
									  x (- y))))))))))))))))



    (with-open-file (s "/dev/shm/index.html" :direction :output
		       :if-does-not-exist :create
		       :if-exists :supersede)
      (cl-who:with-html-output (s)
	(cl-who:htm (:html (:body (:img :src "file:///dev/shm/o.svg")
				  (:p "hallo"))))))



    ;; <svg width="200px" height="200px" viewBox="0 0 200 200">
    ;;    <circle cx="30" cy="30" r="20" style="stroke: black; fill: none;"/>
    ;;    <circle cx="80" cy="30" r="20"
    ;;       style="stroke-width: 5; stroke: black; fill: none;"/>

    ;;    <ellipse cx="30" cy="80" rx="10" ry="20" 
    ;;       style="stroke: black; fill: none;"/>
    ;;    <ellipse cx="80" cy="80" rx="20" ry="10" 
    ;;       style="stroke: black; fill: none;"/>
    ;; </svg>

;; (require :clack)

;; (defpackage simple-app
;;   (:use :cl
;; 	:clack))
;; (in-package :simple-app)

;; (defvar *handler*
;;   (clackup
;;    #'(lambda (env)
;;        '(200 (:content-type "text/plain") ("Hello, Clack!")))))

