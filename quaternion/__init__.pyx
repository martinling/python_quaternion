# Quaternion math for Python
#
# Copyright (c) 2012 Martin Ling <martin@earth.li>
#
# All rights reserved.
# 
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#
#     * Redistributions of source code must retain the above copyright
#       notice, this list of conditions and the following disclaimer.
#
#     * Redistributions in binary form must reproduce the above copyright
#       notice, this list of conditions and the following disclaimer in the
#       documentation and/or other materials provided with the distribution.
#
#     * The name of the author(s) may not be used to endorse or promote
#       products derived from this software without specific prior written
#       permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
# ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTERS BE
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
# POSSIBILITY OF SUCH DAMAGE.

from __future__ import division

cdef extern from "quaternion.h":
    ctypedef struct quaternion:
        double w
        double x
        double y
        double z
    int quaternion_equal(quaternion q1, quaternion q2)
    int quaternion_not_equal(quaternion q1, quaternion q2)
    int quaternion_less(quaternion q1, quaternion q2)
    int quaternion_less_equal(quaternion q1, quaternion q2)
    double quaternion_magnitude(quaternion q)
    double quaternion_dot(quaternion q1, quaternion q2)
    quaternion quaternion_negative(quaternion q)
    quaternion quaternion_conjugate(quaternion q)
    quaternion quaternion_add(quaternion q1, quaternion q2)
    quaternion quaternion_add_scalar(quaternion q, double s)
    quaternion quaternion_subtract(quaternion q1, quaternion q2)
    quaternion quaternion_subtract_scalar(quaternion q, double s)
    quaternion quaternion_multiply(quaternion q1, quaternion q2)
    quaternion quaternion_multiply_scalar(quaternion q, double s)
    quaternion quaternion_divide(quaternion q1, quaternion q2)
    quaternion quaternion_divide_scalar(quaternion q, double s)
    quaternion quaternion_log(quaternion q)
    quaternion quaternion_exp(quaternion q)
    quaternion quaternion_power(quaternion q1, quaternion q2)
    quaternion quaternion_power_scalar(quaternion q, double s)

cdef class Quaternion:

    cdef quaternion _value

    def __cinit__(self, double w=0, double x=0, double y=0, double z=0):
        self._value = quaternion(w, x, y, z)

    # Properties

    property w:
        """ Real component of this quaternion. """
        def __get__(self):
            return self._value.w
        def __set__(self, double value):
            self._value.w = value

    property x:
        """ First imaginary component of this quaternion. """
        def __get__(self):
            return self._value.x
        def __set__(self, double value):
            self._value.x = value

    property y:
        """ Second imaginary component of this quaternion. """
        def __get__(self):
            return self._value.y
        def __set__(self, double value):
            self._value.y = value

    property z:
        """ Third imaginary component of this quaternion. """
        def __get__(self):
            return self._value.z
        def __set__(self, double value):
            self._value.z = value

    property components:
        """ The (w, x, y, z) components of this quaternion as a tuple. """
        def __get__(Quaternion self):
            return (self.w, self.x, self.y, self.z)

    property magnitude:
        """ Magnitude of this quaternion. """
        def __get__(Quaternion self):
            return quaternion_magnitude(self._value)

    property conjugate:
        """ Conjugate of this quaternion. """
        def __get__(Quaternion self):
            result = Quaternion()
            result._value = quaternion_conjugate(self._value)
            return result

    property vector:
        """ The (x, y, z) components of this quaternion as a tuple. """
        def __get__(Quaternion self):
            return (self.x, self.y, self.z)

    # Identity, comparison and serialization

    def __repr__(Quaternion self):
        return "Quaternion(%r, %r, %r, %r)" % (self.w, self.x, self.y, self.z)

    __hash__ = None  # Not hashable because instances are mutable.

    def copy(Quaternion self):
        """ Obtain a copy of this quaternion. """
        return Quaternion(self.w, self.x, self.y, self.z)

    __copy__ = copy

    def __richcmp__(Quaternion self, Quaternion other, int op):
        if op == 0:
            return quaternion_less(self._value, other._value)
        elif op == 1:
            return quaternion_less_equal(self._value, other._value)
        elif op == 2:
            return quaternion_equal(self._value, other._value)
        elif op == 3:
            return quaternion_not_equal(self._value, other._value)
        else:
            return NotImplemented

    def __getstate__(self):
        return (self.w, self.x, self.y, self.z)

    def __setstate__(self, state):
        self.__init__(*state)

    def __reduce__(Quaternion self):
        return (Quaternion, (self.w, self.x, self.y, self.z))

    # General arithmetic (creates new instances)

    def __neg__(Quaternion self):
        cdef Quaternion result = Quaternion()
        result._value = quaternion_negative(self._value)
        return result

    def __invert__(Quaternion self):
        cdef Quaternion result = Quaternion()
        result._value = quaternion_conjugate(self._value)
        return result

    def __add__(a, b):
        cdef Quaternion result = Quaternion()
        cdef Quaternion qa, qb
        cdef double sa, sb
        if isinstance(a, Quaternion) and isinstance(b, Quaternion):
            qa, qb = a, b
            result._value = quaternion_add(qa._value, qb._value)
        elif isinstance(a, Quaternion):
            qa, sb = a, b
            result._value = quaternion_add_scalar(qa._value, sb)
        elif isinstance(b, Quaternion):
            sa, qb = a, b
            result._value = quaternion_add_scalar(qb._value, sa)
        else:
            return NotImplemented
        return result

    def __sub__(a, b):
        cdef Quaternion result = Quaternion()
        cdef Quaternion qa, qb
        cdef double sa, sb
        if isinstance(a, Quaternion) and isinstance(b, Quaternion):
            qa, qb = a, b
            result._value = quaternion_subtract(qa._value, qb._value)
        elif isinstance(a, Quaternion):
            qa, sb = a, b
            result._value = quaternion_subtract_scalar(qa._value, sb)
        elif isinstance(b, Quaternion):
            sa, qb = a, b
            result._value = quaternion_subtract_scalar(qb._value, sa)
        else:
            return NotImplemented
        return result

    def __mul__(a, b):
        cdef Quaternion result = Quaternion()
        cdef Quaternion qa, qb
        cdef double sa, sb
        if isinstance(a, Quaternion) and isinstance(b, Quaternion):
            qa, qb = a, b
            result._value = quaternion_multiply(qa._value, qb._value)
        elif isinstance(a, Quaternion):
            qa, sb = a, b
            result._value = quaternion_multiply_scalar(qa._value, sb)
        elif isinstance(b, Quaternion):
            sa, qb = a, b
            result._value = quaternion_multiply_scalar(qb._value, sa)
        else:
            return NotImplemented
        return result

    def __div__(a, b):
        cdef Quaternion result = Quaternion()
        cdef Quaternion qa, qb
        cdef double sa, sb
        if isinstance(a, Quaternion) and isinstance(b, Quaternion):
            qa, qb = a, b
            result._value = quaternion_divide(qa._value, qb._value)
        elif isinstance(a, Quaternion):
            qa, sb = a, b
            result._value = quaternion_divide_scalar(qa._value, sb)
        elif isinstance(b, Quaternion):
            sa, qb = a, b
            result._value = quaternion_divide_scalar(qb._value, sa)
        else:
            return NotImplemented
        return result

    __truediv__ = __div__

    def __pow__(a, b, c):
        cdef Quaternion result = Quaternion()
        cdef Quaternion qa, qb
        cdef double sa, sb
        if isinstance(a, Quaternion) and isinstance(b, Quaternion):
            qa, qb = a, b
            result._value = quaternion_power(qa._value, qb._value)
        elif isinstance(a, Quaternion):
            qa, sb = a, b
            result._value = quaternion_power_scalar(qa._value, sb)
        elif isinstance(b, Quaternion):
            sa, qb = a, b
            result._value = quaternion_power_scalar(qb._value, sa)
        else:
            return NotImplemented
        return result

    def log(Quaternion self):
        """ Natural logarithm of this quaternion. """
        cdef Quaternion result = Quaternion()
        result._value = quaternion_log(self._value)
        return result

    def exp(Quaternion self):
        """ Exponential of this quaternion. """
        cdef Quaternion result = Quaternion()
        result._value = quaternion_exp(self._value)
        return result

    def dot(Quaternion self, Quaternion other):
        """ Dot prodct of this quaternion with another. """
        return quaternion_dot(self._value, other._value)

    # In-place arithmetic (modifies existing instance)

    def __iadd__(Quaternion self, other):
        cdef Quaternion qo
        if isinstance(other, Quaternion):
            qo = other
            self._value = quaternion_add(self._value, qo._value)
        else:
            self._value = quaternion_add_scalar(self._value, other)
        return self

    def __isub__(Quaternion self, other):
        cdef Quaternion qo
        if isinstance(other, Quaternion):
            qo = other
            self._value = quaternion_subtract(self._value, qo._value)
        else:
            self._value = quaternion_subtract_scalar(self._value, other)
        return self

    def __imul__(Quaternion self, other):
        cdef Quaternion qo
        if isinstance(other, Quaternion):
            qo = other
            self._value = quaternion_multiply(self._value, qo._value)
        else:
            self._value = quaternion_multiply_scalar(self._value, other)
        return self

    def __idiv__(Quaternion self, other):
        cdef Quaternion qo
        if isinstance(other, Quaternion):
            qo = other
            self._value = quaternion_divide(self._value, qo._value)
        else:
            self._value = quaternion_divide_scalar(self._value, other)
        return self

    __itruediv__ = __div__

    def __ipow__(Quaternion self, other):
        cdef Quaternion qo
        if isinstance(other, Quaternion):
            qo = other
            self._value = quaternion_power(self._value, qo._value)
        else:
            self._value = quaternion_power_scalar(self._value, other)
        return self

    def normalise(Quaternion self):
        """ Normalise this quaternion to unit length. """
        self._value = quaternion_multiply_scalar(self._value,
            1 / quaternion_magnitude(self._value))
        return self

    def negate(Quaternion self):
        """ Negate this quaternion. """
        self._value = quaternion_negative(self._value)
        return self
