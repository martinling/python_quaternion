Quaternion module for Python
----------------------------

Version 0.8, 19th August 2012

Copyright (C) 2012 Martin Ling <martin@earth.li>
Released under a BSD-style license (see LICENSE.txt)

This Python module provides a Quaternion type that provides easy access to
fast quaternion math operations implemented in C.

It is designed as a standalone alternative to two other Python quaternion math
implementations that I have contributed to previously: firstly that in the
IMUSim package (http://www.imusim.org/), and secondly one designed to
integrate into NumPy (http://github.com/martinling/numpy_quaternion/). While
those versions offer some additional capabilities, this one aims to provide a
clean interface with no external dependencies.

Full support is provided for all basic arithmetic operations, as well as less
common ones including the quaternion logarithm and exponent. All arithmetic
expressions are polymorphic, accepting both quaternion and scalar arguments.

Brief usage example:

 >>> from quaternion import Quaternion
 >>> Quaternion(1,0,0,0)
 Quaternion(1.0, 0.0, 0.0, 0.0)
 >>> a = Quaternion(1,2,3,4)
 >>> b = Quaternion(5,6,7,8)
 >>> a * b
 Quaternion(-60.0, 12.0, 30.0, 24.0)
 >>> (a * b**2 + ~a).log()
 Quaternion(5.15495943019, 1.16003608989, 1.35337543821, 1.54671478653)

Implementation notes:

- Quaternion components are stored as doubles.

- Comparison operations follow the same lexicographic ordering as tuples.

- The inversion operator (~) can be used to obtain the quaternion conjugate.

- The math implementation in quaternion.{c, h} is designed to be independently
  reusable for other projects.

Building:

 The Python wrapper is written in Cython (http://www.cython.org) and requires
 the Cython compiler to to generate the C source for the extension:

 $ cython quaternion/__init__.pyx

 This step is only required for a fresh source checkout or if making changes
 to the .pyx source. Once done, the result is a standard extension package
 which can be distributed normally and requires only a working C compiler &
 Python development environment to build:

 $ python setup.py build

 To install (may require administrator rights):

 # python setup.py install

Bugs & Problems:

 Please report at: https://github.com/martinling/python_quaternion/issues
