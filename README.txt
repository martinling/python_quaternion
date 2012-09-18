Quaternion module for Python
----------------------------

Copyright (C) 2012 Martin Ling <martin@earth.li>
Released under a BSD-style license (see LICENSE.txt)

This Python module provides a Quaternion type that provides easy access to
fast quaternion math operations implemented in C. The aim is to provide a
comprehensive set of capabilities in a single library with only NumPy as a
prerequisite.

The types provided are a Quaternion class for representing single quaternion
values, and a QuaternionArray class for representing sequences of quaternions
with efficient array operations.

Full support is provided for all basic arithmetic operations, as well as less
common ones including the quaternion logarithm and exponent. All arithmetic
expressions are polymorphic: arguments may be scalars or quaternions, and
either single values or arrays.

In addition, some methods are included to aid the use of quaternions as a
representation of 3D rotations. Currently these include rotation operations
on vectors and conversions to and from Euler angle sequences.

Some brief example usage::

 >>> from quaternion import Quaternion, QuaternionArray
 >>> Quaternion(1,0,0,0)
 Quaternion(1.0, 0.0, 0.0, 0.0)
 >>> a = Quaternion(1,2,3,4)
 >>> b = Quaternion(5,6,7,8)
 >>> a * b
 Quaternion(-60.0, 12.0, 30.0, 24.0)
 >>> (a * b**2 + ~a).log()
 Quaternion(5.15495943019, 1.16003608989, 1.35337543821, 1.54671478653)
 >>> qa = QuaternionArray([a, b])
 >>> qa
 QuaternionArray(
 array([[ 1.,  2.,  3.,  4.],
        [ 5.,  6.,  7.,  8.]]))
 >>> qa * a
 QuaternionArray(
 array([[-28.,   4.,   6.,   8.],
        [-60.,  20.,  14.,  32.]]))

Implementation notes:

- Quaternions are mutable.

- Quaternion components are stored as doubles.

- Comparison operations follow the same lexicographic ordering as tuples.

- The inversion operator (~) can be used to obtain the quaternion conjugate.

- The math implementation in quaternion.{c, h} is designed to be independent
  of the Python wrapping code and reusable for other projects.

Background:

This library was designed as an alternative to two other Python quaternion
math implementations that I have contributed to previously: firstly that in
the IMUSim package (http://www.imusim.org/), and secondly an experimental one
designed to integrate quaternions more deeply into NumPy as a custom dtype
(http://github.com/martinling/numpy_quaternion/).

The IMUSim quaternion library provides similar capabilities to this one, but
currently requires installation of the whole IMUSim package and also has some
idiosyncracies. Eventually this library may replace it in a future version.

The NumPy integration approach was promising, but to finish the job properly
would also require some changes to numpy itself.

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
