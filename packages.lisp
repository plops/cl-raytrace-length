(defpackage :base
  (:use :cl)
  (:export 
   #:num #:vec #:with-arrays #:v #:copy-vec
   #:.+ #:.- #:.* #:./ #:dot #:norm #:normalize
   #:cross #:.s #:vx #:vy #:vz #:v-spherical
   #:check-unit-vector #:check-range
   #:req #:+pi+ #:+2*pi+ #:+pi/2+
   #:mat
   #:m
   #:rotation-matrix
   #:determinant
   #:transpose
   #:m*
   #:chop
   #:*read-default-float-format*))

(defpackage :raytrace
  (:use :cl :base)
  (:export
   #:refract-plane
   #:intersect-plane
   #:intersect-sphere
   #:refract-objective-detection
   #:refract-objective-illumination
   #:refract-thin-lens
   #:aberrate-index-plane))
