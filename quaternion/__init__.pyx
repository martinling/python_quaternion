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

import numpy as np
cimport numpy as np

cdef extern from "quaternion.h":
    ctypedef struct quaternion:
        double w
        double x
        double y
        double z
    int quaternion_equal(quaternion *q1, quaternion *q2)
    int quaternion_not_equal(quaternion *q1, quaternion *q2)
    int quaternion_less(quaternion *q1, quaternion *q2)
    int quaternion_less_equal(quaternion *q1, quaternion *q2)
    double quaternion_magnitude(quaternion *q)
    double quaternion_dot(quaternion *q1, quaternion *q2)
    void quaternion_negative(quaternion *q, quaternion *r)
    void quaternion_conjugate(quaternion *q, quaternion *r)
    void quaternion_add(quaternion *q1, quaternion *q2, quaternion *r)
    void quaternion_add_scalar(quaternion *q, double s, quaternion *r)
    void quaternion_subtract(quaternion *q1, quaternion *q2, quaternion *r)
    void quaternion_subtract_scalar(quaternion *q, double s, quaternion *r)
    void quaternion_multiply(quaternion *q1, quaternion *q2, quaternion *r)
    void quaternion_multiply_scalar(quaternion *q, double s, quaternion *r)
    void quaternion_divide(quaternion *q1, quaternion *q2, quaternion *r)
    void quaternion_divide_scalar(quaternion *q, double s, quaternion *r)
    void quaternion_log(quaternion *q, quaternion *r)
    void quaternion_exp(quaternion *q, quaternion *r)
    void quaternion_power(quaternion *q1, quaternion *q2, quaternion *r)
    void quaternion_power_scalar(quaternion *q, double s, quaternion *r)
    void quaternion_rotate_vector(quaternion *q, double v[3], double r[3])
    void quaternion_rotate_frame(quaternion *q, double v[3], double r[3])
    void quaternion_from_euler(char *order, double angles[3], quaternion *r)
    void quaternion_to_euler(quaternion *q, char *order, double r[3])

cdef extern from "stdlib.h":
    void *malloc(size_t size)
    void free(void *ptr)

cdef extern from "string.h":
    size_t strlen(char *s)

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
            return quaternion_magnitude(&self._value)

    property conjugate:
        """ Conjugate of this quaternion. """
        def __get__(Quaternion self):
            result = Quaternion()
            quaternion_conjugate(&self._value, &result._value)
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
            return quaternion_less(&self._value, &other._value)
        elif op == 1:
            return quaternion_less_equal(&self._value, &other._value)
        elif op == 2:
            return quaternion_equal(&self._value, &other._value)
        elif op == 3:
            return quaternion_not_equal(&self._value, &other._value)
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
        quaternion_negative(&self._value, &result._value)
        return result

    def __invert__(Quaternion self):
        cdef Quaternion result = Quaternion()
        quaternion_conjugate(&self._value, &result._value)
        return result

    def __add__(a, b):
        cdef Quaternion result = Quaternion()
        cdef Quaternion qa, qb
        cdef double sa, sb
        if isinstance(a, Quaternion) and isinstance(b, Quaternion):
            qa, qb = a, b
            quaternion_add(&qa._value, &qb._value, &result._value)
        elif isinstance(a, Quaternion):
            qa, sb = a, b
            quaternion_add_scalar(&qa._value, sb, &result._value)
        elif isinstance(b, Quaternion):
            sa, qb = a, b
            quaternion_add_scalar(&qb._value, sa, &result._value)
        else:
            return NotImplemented
        return result

    def __sub__(a, b):
        cdef Quaternion result = Quaternion()
        cdef Quaternion qa, qb
        cdef double sa, sb
        if isinstance(a, Quaternion) and isinstance(b, Quaternion):
            qa, qb = a, b
            quaternion_subtract(&qa._value, &qb._value, &result._value)
        elif isinstance(a, Quaternion):
            qa, sb = a, b
            quaternion_subtract_scalar(&qa._value, sb, &result._value)
        elif isinstance(b, Quaternion):
            sa, qb = a, b
            quaternion_negative(&qb._value, &result._value)
            quaternion_add_scalar(&result._value, sa, &result._value)
        else:
            return NotImplemented
        return result

    def __mul__(a, b):
        cdef Quaternion result = Quaternion()
        cdef Quaternion qa, qb
        cdef double sa, sb
        if isinstance(a, Quaternion) and isinstance(b, Quaternion):
            qa, qb = a, b
            quaternion_multiply(&qa._value, &qb._value, &result._value)
        elif isinstance(a, Quaternion):
            qa, sb = a, b
            quaternion_multiply_scalar(&qa._value, sb, &result._value)
        elif isinstance(b, Quaternion):
            sa, qb = a, b
            quaternion_multiply_scalar(&qb._value, sa, &result._value)
        else:
            return NotImplemented
        return result

    def __div__(a, b):
        cdef Quaternion result = Quaternion()
        cdef Quaternion qa, qb
        cdef double sa, sb
        if isinstance(a, Quaternion) and isinstance(b, Quaternion):
            qa, qb = a, b
            quaternion_divide(&qa._value, &qb._value, &result._value)
        elif isinstance(a, Quaternion):
            qa, sb = a, b
            quaternion_divide_scalar(&qa._value, sb, &result._value)
        elif isinstance(b, Quaternion):
            sa, qb = a, b
            quaternion_power_scalar(&qb._value, -1, &result._value)
            quaternion_multiply_scalar(&result._value, sa, &result._value)
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
            quaternion_power(&qa._value, &qb._value, &result._value)
        elif isinstance(a, Quaternion):
            qa, sb = a, b
            quaternion_power_scalar(&qa._value, sb, &result._value)
        elif isinstance(b, Quaternion):
            qa, qb = Quaternion(a, 0, 0, 0), b
            quaternion_power(&qa._value, &qb._value, &result._value)
        else:
            return NotImplemented
        return result

    def log(Quaternion self):
        """ Natural logarithm of this quaternion. """
        cdef Quaternion result = Quaternion()
        quaternion_log(&self._value, &result._value)
        return result

    def exp(Quaternion self):
        """ Exponential of this quaternion. """
        cdef Quaternion result = Quaternion()
        quaternion_exp(&self._value, &result._value)
        return result

    def dot(Quaternion self, Quaternion other):
        """ Dot prodct of this quaternion with another. """
        return quaternion_dot(&self._value, &other._value)

    # In-place arithmetic (modifies existing instance)

    def __iadd__(Quaternion self, other):
        cdef Quaternion qo
        if isinstance(other, Quaternion):
            qo = other
            quaternion_add(&self._value, &qo._value, &self._value)
        else:
            quaternion_add_scalar(&self._value, other, &self._value)
        return self

    def __isub__(Quaternion self, other):
        cdef Quaternion qo
        if isinstance(other, Quaternion):
            qo = other
            quaternion_subtract(&self._value, &qo._value, &self._value)
        else:
            quaternion_subtract_scalar(&self._value, other, &self._value)
        return self

    def __imul__(Quaternion self, other):
        cdef Quaternion qo
        if isinstance(other, Quaternion):
            qo = other
            quaternion_multiply(&self._value, &qo._value, &self._value)
        else:
            quaternion_multiply_scalar(&self._value, other, &self._value)
        return self

    def __idiv__(Quaternion self, other):
        cdef Quaternion qo
        if isinstance(other, Quaternion):
            qo = other
            quaternion_divide(&self._value, &qo._value, &self._value)
        else:
            quaternion_divide_scalar(&self._value, other, &self._value)
        return self

    __itruediv__ = __idiv__

    def __ipow__(Quaternion self, other):
        cdef Quaternion qo
        if isinstance(other, Quaternion):
            qo = other
            quaternion_power(&self._value, &qo._value, &self._value)
        else:
            quaternion_power_scalar(&self._value, other, &self._value)
        return self

    def normalise(Quaternion self):
        """ Normalise this quaternion to unit length. """
        quaternion_multiply_scalar(&self._value,
            1 / quaternion_magnitude(&self._value), &self._value)
        return self

    def negate(Quaternion self):
        """ Negate this quaternion. """
        quaternion_negative(&self._value, &self._value)
        return self

    # Rotations on vectors.

    def rotate_vector(Quaternion self, v):
        """
        Use this quaternion to rotate a vector.
        """
        cdef np.ndarray[np.float_t, ndim=1] vec1d
        cdef np.ndarray[np.float_t, ndim=1] res1d
        cdef np.ndarray[np.float_t, ndim=2] vec2d
        cdef np.ndarray[np.float_t, ndim=2] res2d
        cdef double[3] vec, res
        if isinstance(v, np.ndarray):
            if v.ndim == 2:
                vec2d = v
                res2d = np.empty_like(v)
                for i in range(len(v)):
                    vec1d = vec2d[i]
                    res1d = res2d[i]
                    quaternion_rotate_vector(&self._value,
                        <double *> vec1d.data, <double *> res1d.data)
                result = res2d
            else:
                vec1d = v
                res1d = np.empty_like(v)
                quaternion_rotate_vector(&self._value,
                    <double *> vec1d.data, <double *> res1d.data)
                result = res1d
        else:
            for i in range(3):
                vec[i] = v[i]
            quaternion_rotate_vector(&self._value, vec, res)
            result = type(v)(res[0], res[1], res[2])
        return result

    def rotate_frame(Quaternion self, v):
        """
        Use this quaternion to rotate the co-ordinate frame of a vector.
        """
        cdef np.ndarray[np.float_t, ndim=1] vec1d
        cdef np.ndarray[np.float_t, ndim=1] res1d
        cdef np.ndarray[np.float_t, ndim=2] vec2d
        cdef np.ndarray[np.float_t, ndim=2] res2d
        cdef double[3] vec, res
        if isinstance(v, np.ndarray):
            if v.ndim == 2:
                vec2d = v
                res2d = np.empty_like(v)
                for i in range(len(v)):
                    vec1d = vec2d[i]
                    res1d = res2d[i]
                    quaternion_rotate_frame(&self._value,
                        <double *> vec1d.data, <double *> res1d.data)
                result = res2d
            else:
                vec1d = v
                res1d = np.empty_like(v)
                quaternion_rotate_frame(&self._value,
                    <double *> vec1d.data, <double *> res1d.data)
                result = res1d
        else:
            for i in range(3):
                vec[i] = v[i]
            quaternion_rotate_frame(&self._value, vec, res)
            result = type(v)(res[0], res[1], res[2])
        return result

    # Conversions to and from other rotation formats.

    @classmethod
    def from_euler(cls, char *order, angles):
        cdef Quaternion result = Quaternion()
        cdef double *ang = <double *> malloc(strlen(order) * sizeof(double))
        for i in range(strlen(order)):
            ang[i] = angles[i]
        quaternion_from_euler(order, ang, &result._value)
        free(ang)
        return result

    def to_euler(Quaternion self, char *order):
        cdef double[3] res
        quaternion_to_euler(&self._value, order, res)
        return res[0], res[1], res[2]

cdef class QuaternionArray:

    cdef quaternion *_values
    cdef np.ndarray _components
    cdef unsigned int length

    def __cinit__(QuaternionArray self, values):
        cdef np.ndarray[np.float_t, ndim=2] components
        assert np.ndim(values) == 2 and np.shape(values)[1] == 4, \
            "Values must have shape (N, 4)"
        components = np.asarray(values, dtype=np.float)
        self._components = components
        self._values = <quaternion *> &components[0, 0]
        self.length = len(components)

    # Sequence operations

    def __len__(QuaternionArray self):
        return self.length

    def __getitem__(QuaternionArray self, key):
        if isinstance(key, slice):
            return QuaternionArray(self.components[key])
        else:
            return Quaternion(*self.components[key])
    
    def __setitem__(QuaternionArray self, key, value):
        if isinstance(value, (Quaternion, QuaternionArray)):
            self.components[key] = value.components
        else:
            self.components[key] = value

    # Properties

    property components:
        """ N x 4 array of quaternion [w, x, y, z] components. """
        def __get__(QuaternionArray self):
            cdef np.ndarray[np.float_t, ndim=2] components = self._components
            return components
        def __set__(QuaternionArray self, value):
            self._components[:] = value

    property w:
        """ Real components of these quaternions. """
        def __get__(QuaternionArray self):
            return self.components[:,0]
        def __set__(QuaternionArray self, value):
            self.components[:,0] = value

    property x:
        """ First imaginary components of these quaternions. """
        def __get__(QuaternionArray self):
            return self.components[:,1]
        def __set__(QuaternionArray self, value):
            self.components[:,1] = value

    property y:
        """ Second imaginary components of these quaternions. """
        def __get__(QuaternionArray self):
            return self.components[:,2]
        def __set__(QuaternionArray self, value):
            self.components[:,2] = value

    property z:
        """ Third imaginary components of these quaternions. """
        def __get__(QuaternionArray self):
            return self.components[:,3]
        def __set__(QuaternionArray self, value):
            self.components[:,3] = value

    property magnitude:
        """ Magnitudes of these quaternions. """
        def __get__(QuaternionArray self):
            cdef np.ndarray[np.float_t, ndim=1] result = np.empty(self.length)
            for i in range(self.length):
                result[i] = quaternion_magnitude(&self._values[i])
            return result

    property conjugate:
        """ Conjugates of these quaternions. """
        def __get__(QuaternionArray self):
            cdef QuaternionArray result = QuaternionArray.like(self)
            for i in range(self.length):
                quaternion_conjugate(&self._values[i], &result._values[i])
            return result

    property vector:
        """ N x 3 array of the [x, y, z] compnents of these quaternions. """
        def __get__(QuaternionArray self):
            return self.components[:,1:]
        def __set__(QuaternionArray self, value):
            self.components[:,1:] = value

    # Identity, comparison and serialization

    def __repr__(QuaternionArray self):
        return "QuaternionArray(\n%s)" % repr(self.components)

    __hash__ = None  # Not hashable because instances are mutable.

    def copy(QuaternionArray self):
        """ Obtain a copy of these quaternions. """
        return QuaternionArray(self.components.copy())

    @classmethod
    def like(cls, ref):
        """ Create an empty QuaternionArray with the length of an existing sequence. """
        return QuaternionArray(np.empty((len(ref), 4)))

    __copy__ = copy

    def __getstate__(QuaternionArray self):
        return (self.components,)

    def __setstate__(QuaternionArray self, state):
        self.__init__(*state)

    def __reduce__(QuaternionArray self):
        return (QuaternionArray, (self.components))

    # General arithmetic (creates new instances)

    def __neg__(QuaternionArray self):
        cdef QuaternionArray result = QuaternionArray.like(self)
        for i in range(self.length):
            quaternion_negative(&self._values[i], &result._values[i])
        return result

    def __invert__(QuaternionArray self):
        cdef QuaternionArray result = QuaternionArray.like(self)
        for i in range(self.length):
            quaternion_conjugate(&self._values[i], &result._values[i])
        return result

    def __add__(a, b):
        cdef QuaternionArray result
        cdef QuaternionArray qqa, qqb
        cdef Quaternion qa, qb
        cdef np.ndarray[np.float_t, ndim=1] ssa, ssb
        cdef double sa, sb
        if isinstance(a, QuaternionArray) and isinstance(b, QuaternionArray):
            qqa, qqb = a, b
            assert qqa.length == qqb.length, "Arrays must have equal length"
            result = QuaternionArray.like(qqa)
            for i in range(qqa.length):
                quaternion_add(&qqa._values[i], &qqb._values[i], &result._values[i])
        elif isinstance(a, QuaternionArray):
            qqa = a
            result = QuaternionArray.like(qqa)
            if isinstance(b, Quaternion):
                qb = b
                for i in range(qqa.length):
                    quaternion_add(&qqa._values[i], &qb._value, &result._values[i])
            elif np.shape(b) == (qqa.length,):
                ssb = np.asarray(b)
                for i in range(qqa.length):
                    quaternion_add_scalar(&qqa._values[i], ssb[i], &result._values[i])
            else:
                sb = b
                for i in range(qqa.length):
                    quaternion_add_scalar(&qqa._values[i], sb, &result._values[i])
        elif isinstance(b, QuaternionArray):
            qqb = b
            result = QuaternionArray.like(qqb)
            if isinstance(a, Quaternion):
                qa = a
                for i in range(qqb.length):
                    quaternion_add(&qa._value, &qqb._values[i], &result._values[i])
            elif np.shape(a) == (qqb.length,):
                ssa = np.asarray(a)
                for i in range(qqb.length):
                    quaternion_add_scalar(&qqb._values[i], ssa[i], &result._values[i])
            else:
                sa = a
                for i in range(qqb.length):
                    quaternion_add_scalar(&qqb._values[i], sa, &result._values[i])
        else:
            return NotImplemented
        return result

    def __sub__(a, b):
        cdef QuaternionArray result
        cdef QuaternionArray qqa, qqb
        cdef Quaternion qa, qb
        cdef np.ndarray[np.float_t, ndim=1] ssa, ssb
        cdef double sa, sb
        if isinstance(a, QuaternionArray) and isinstance(b, QuaternionArray):
            qqa, qqb = a, b
            assert qqa.length == qqb.length, "Arrays must have equal length"
            result = QuaternionArray.like(qqa)
            for i in range(qqa.length):
                quaternion_subtract(&qqa._values[i], &qqb._values[i], &result._values[i])
        elif isinstance(a, QuaternionArray):
            qqa = a
            result = QuaternionArray.like(qqa)
            if isinstance(b, Quaternion):
                qb = b
                for i in range(qqa.length):
                    quaternion_subtract(&qqa._values[i], &qb._value, &result._values[i])
            elif np.shape(b) == (qqa.length,):
                ssb = np.asarray(b)
                for i in range(qqa.length):
                    quaternion_subtract_scalar(&qqa._values[i], ssb[i], &result._values[i])
            else:
                sb = b
                for i in range(qqa.length):
                    quaternion_subtract_scalar(&qqa._values[i], sb, &result._values[i])
        elif isinstance(b, QuaternionArray):
            qqb = b
            result = QuaternionArray.like(qqb)
            if isinstance(a, Quaternion):
                qa = a
                for i in range(qqb.length):
                    quaternion_subtract(&qa._value, &qqb._values[i], &result._values[i])
            elif np.shape(a) == (qqb.length,):
                ssa = np.asarray(a)
                for i in range(qqb.length):
                    quaternion_negative(&qqb._values[i], &result._values[i])
                    quaternion_add_scalar(&result._values[i], ssa[i], &result._values[i])
            else:
                sa = a
                for i in range(qqb.length):
                    quaternion_negative(&qqb._values[i], &result._values[i])
                    quaternion_add_scalar(&result._values[i], sa, &result._values[i])
        else:
            return NotImplemented
        return result

    def __mul__(a, b):
        cdef QuaternionArray result
        cdef QuaternionArray qqa, qqb
        cdef Quaternion qa, qb
        cdef np.ndarray[np.float_t, ndim=1] ssa, ssb
        cdef double sa, sb
        if isinstance(a, QuaternionArray) and isinstance(b, QuaternionArray):
            qqa, qqb = a, b
            assert qqa.length == qqb.length, "Arrays must have equal length"
            result = QuaternionArray.like(qqa)
            for i in range(qqa.length):
                quaternion_multiply(&qqa._values[i], &qqb._values[i], &result._values[i])
        elif isinstance(a, QuaternionArray):
            qqa = a
            result = QuaternionArray.like(qqa)
            if isinstance(b, Quaternion):
                qb = b
                for i in range(qqa.length):
                    quaternion_multiply(&qqa._values[i], &qb._value, &result._values[i])
            elif np.shape(b) == (qqa.length,):
                ssb = np.asarray(b)
                for i in range(qqa.length):
                    quaternion_multiply_scalar(&qqa._values[i], ssb[i], &result._values[i])
            else:
                sb = b
                for i in range(qqa.length):
                    quaternion_multiply_scalar(&qqa._values[i], sb, &result._values[i])
        elif isinstance(b, QuaternionArray):
            qqb = b
            result = QuaternionArray.like(qqb)
            if isinstance(a, Quaternion):
                qa = a
                for i in range(qqb.length):
                    quaternion_multiply(&qa._value, &qqb._values[i], &result._values[i])
            elif np.shape(a) == (qqb.length,):
                ssa = np.asarray(a)
                for i in range(qqb.length):
                    quaternion_multiply_scalar(&qqb._values[i], ssa[i], &result._values[i])
            else:
                sa = a
                for i in range(qqb.length):
                    quaternion_multiply_scalar(&qqb._values[i], sa, &result._values[i])
        else:
            return NotImplemented
        return result

    def __div__(a, b):
        cdef QuaternionArray result
        cdef QuaternionArray qqa, qqb
        cdef Quaternion qa, qb
        cdef np.ndarray[np.float_t, ndim=1] ssa, ssb
        cdef double sa, sb
        if isinstance(a, QuaternionArray) and isinstance(b, QuaternionArray):
            qqa, qqb = a, b
            assert qqa.length == qqb.length, "Arrays must have equal length"
            result = QuaternionArray.like(qqa)
            for i in range(qqa.length):
                quaternion_divide(&qqa._values[i], &qqb._values[i], &result._values[i])
        elif isinstance(a, QuaternionArray):
            qqa = a
            result = QuaternionArray.like(qqa)
            if isinstance(b, Quaternion):
                qb = b
                for i in range(qqa.length):
                    quaternion_divide(&qqa._values[i], &qb._value, &result._values[i])
            elif np.shape(b) == (qqa.length,):
                ssb = np.asarray(b)
                for i in range(qqa.length):
                    quaternion_divide_scalar(&qqa._values[i], ssb[i], &result._values[i])
            else:
                sb = b
                for i in range(qqa.length):
                    quaternion_divide_scalar(&qqa._values[i], sb, &result._values[i])
        elif isinstance(b, QuaternionArray):
            qqb = b
            result = QuaternionArray.like(qqb)
            if isinstance(a, Quaternion):
                qa = a
                for i in range(qqb.length):
                    quaternion_divide(&qa._value, &qqb._values[i], &result._values[i])
            elif np.shape(a) == (qqb.length,):
                ssa = np.asarray(a)
                for i in range(qqb.length):
                    quaternion_power_scalar(&qqb._values[i], -1, &result._values[i])
                    quaternion_multiply_scalar(&result._values[i], ssa[i], &result._values[i])
            else:
                sa = a
                for i in range(qqb.length):
                    quaternion_power_scalar(&qqb._values[i], -1, &result._values[i])
                    quaternion_multiply_scalar(&result._values[i], sa, &result._values[i])
        else:
            return NotImplemented
        return result

    __truediv__ = __div__

    def __pow__(a, b, c):
        cdef QuaternionArray result
        cdef QuaternionArray qqa, qqb
        cdef Quaternion qa, qb
        cdef np.ndarray[np.float_t, ndim=1] ssb
        cdef double sa, sb
        if isinstance(a, QuaternionArray) and isinstance(b, QuaternionArray):
            qqa, qqb = a, b
            assert qqa.length == qqb.length, "Arrays must have equal length"
            result = QuaternionArray.like(qqa)
            for i in range(qqa.length):
                quaternion_power(&qqa._values[i], &qqb._values[i], &result._values[i])
        elif isinstance(a, QuaternionArray):
            qqa = a
            result = QuaternionArray.like(qqa)
            if isinstance(b, Quaternion):
                qb = b
                for i in range(qqa.length):
                    quaternion_power(&qqa._values[i], &qb._value, &result._values[i])
            elif np.shape(b) == (qqa.length,):
                ssb = np.asarray(b)
                for i in range(qqa.length):
                    quaternion_power_scalar(&qqa._values[i], ssb[i], &result._values[i])
            else:
                sb = b
                for i in range(qqa.length):
                    quaternion_power_scalar(&qqa._values[i], sb, &result._values[i])
        elif isinstance(b, QuaternionArray):
            qqb = b
            result = QuaternionArray.like(qqb)
            if isinstance(a, Quaternion):
                qa = a
                for i in range(qqb.length):
                    quaternion_power(&qa._value, &qqb._values[i], &result._values[i])
            elif np.shape(a) == (qqb.length,):
                qqa = QuaternionArray.like(qqb)
                qqa.w = a
                qqa.vector = 0
                for i in range(qqb.length):
                    quaternion_power(&qqa._values[i], &qqb._values[i], &result._values[i])
            else:
                qa = Quaternion(a, 0, 0, 0)
                for i in range(qqb.length):
                    quaternion_power(&qa._value, &qqb._values[i], &result._values[i])
        else:
            return NotImplemented
        return result

    def log(QuaternionArray self):
        """ Natural logarithms of these quaternions. """
        cdef QuaternionArray result = QuaternionArray.like(self)
        for i in range(self.length):
            quaternion_log(&self._values[i], &result._values[i])
        return result

    def exp(QuaternionArray self):
        """ Exponentials of these quaternions. """
        cdef QuaternionArray result = QuaternionArray.like(self)
        for i in range(self.length):
            quaternion_exp(&self._values[i], &result._values[i])
        return result

    def dot(QuaternionArray self, other):
        """ Dot prodct of these quaternions with another quaternion or quaternion array. """
        cdef np.ndarray[np.float_t, ndim=1] result = np.empty(self.length)
        cdef QuaternionArray qqo
        cdef Quaternion qo
        if isinstance(other, QuaternionArray):
            qqo = other
            for i in range(self.length):
                result[i] = quaternion_dot(&self._values[i], &qqo._values[i])
        elif isinstance(other, Quaternion):
            qo = other
            for i in range(self.length):
                result[i] = quaternion_dot(&self._values[i], &qo._value)
        else:
            return NotImplemented
        return result

    # In-place arithmetic (modifies existing instance)

    def __iadd__(QuaternionArray self, other):
        cdef QuaternionArray qqo
        cdef Quaternion qo
        cdef np.ndarray[np.float_t, ndim=1] sso
        if isinstance(other, QuaternionArray):
            qqo = other
            for i in range(self.length):
                quaternion_add(&self._values[i], &qqo._values[i], &self._values[i])
        elif isinstance(other, Quaternion):
            qo = other
            for i in range(self.length):
                quaternion_add(&self._values[i], &qo._value, &self._values[i])
        elif np.shape(other) == (self.length,):
            sso = np.asarray(other)
            for i in range(self.length):
                quaternion_add_scalar(&self._values[i], sso[i], &self._values[i])
        else:
            for i in range(self.length):
                quaternion_add_scalar(&self._values[i], other, &self._values[i])
        return self

    def __isub__(QuaternionArray self, other):
        cdef QuaternionArray qqo
        cdef Quaternion qo
        cdef np.ndarray[np.float_t, ndim=1] sso
        if isinstance(other, QuaternionArray):
            qqo = other
            for i in range(self.length):
                quaternion_subtract(&self._values[i], &qqo._values[i], &self._values[i])
        elif isinstance(other, Quaternion):
            qo = other
            for i in range(self.length):
                quaternion_subtract(&self._values[i], &qo._value, &self._values[i])
        elif np.shape(other) == (self.length,):
            sso = np.asarray(other)
            for i in range(self.length):
                quaternion_subtract_scalar(&self._values[i], sso[i], &self._values[i])
        else:
            for i in range(self.length):
                quaternion_subtract_scalar(&self._values[i], other, &self._values[i])
        return self

    def __imul__(QuaternionArray self, other):
        cdef QuaternionArray qqo
        cdef Quaternion qo
        cdef np.ndarray[np.float_t, ndim=1] sso
        if isinstance(other, QuaternionArray):
            qqo = other
            for i in range(self.length):
                quaternion_multiply(&self._values[i], &qqo._values[i], &self._values[i])
        elif isinstance(other, Quaternion):
            qo = other
            for i in range(self.length):
                quaternion_multiply(&self._values[i], &qo._value, &self._values[i])
        elif np.shape(other) == (self.length,):
            sso = np.asarray(other)
            for i in range(self.length):
                quaternion_multiply_scalar(&self._values[i], sso[i], &self._values[i])
        else:
            for i in range(self.length):
                quaternion_multiply_scalar(&self._values[i], other, &self._values[i])
        return self

    def __idiv__(QuaternionArray self, other):
        cdef QuaternionArray qqo
        cdef Quaternion qo
        cdef np.ndarray[np.float_t, ndim=1] sso
        if isinstance(other, QuaternionArray):
            qqo = other
            for i in range(self.length):
                quaternion_divide(&self._values[i], &qqo._values[i], &self._values[i])
        elif isinstance(other, Quaternion):
            qo = other
            for i in range(self.length):
                quaternion_divide(&self._values[i], &qo._value, &self._values[i])
        elif np.shape(other) == (self.length,):
            sso = np.asarray(other)
            for i in range(self.length):
                quaternion_divide_scalar(&self._values[i], sso[i], &self._values[i])
        else:
            for i in range(self.length):
                quaternion_divide_scalar(&self._values[i], other, &self._values[i])
        return self

    __itruediv__ = __idiv__

    def __ipow__(QuaternionArray self, other):
        cdef QuaternionArray qqo
        cdef Quaternion qo
        cdef np.ndarray[np.float_t, ndim=1] sso
        if isinstance(other, QuaternionArray):
            qqo = other
            for i in range(self.length):
                quaternion_power(&self._values[i], &qqo._values[i], &self._values[i])
        elif isinstance(other, Quaternion):
            qo = other
            for i in range(self.length):
                quaternion_power(&self._values[i], &qo._value, &self._values[i])
        elif np.shape(other) == (self.length,):
            sso = np.asarray(other)
            for i in range(self.length):
                quaternion_power_scalar(&self._values[i], sso[i], &self._values[i])
        else:
            for i in range(self.length):
                quaternion_power_scalar(&self._values[i], other, &self._values[i])
        return self

    def normalise(QuaternionArray self):
        """ Normalise these quaternions to unit length. """
        for i in range(self.length):
            quaternion_multiply_scalar(&self._values[i],
                1 / quaternion_magnitude(&self._values[i]), &self._values[i])
        return self

    def negate(Quaternion self):
        """ Negate these quaternions. """
        for i in range(self.length):
            quaternion_negative(&self._values[i], &self._values[i])
        return self
