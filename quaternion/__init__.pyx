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

cdef class Quaternion:

    cdef quaternion *value
    cdef np.ndarray _components

    def __cinit__(self, double w=0, double x=0, double y=0, double z=0, value=None):
        cdef np.ndarray[np.float_t, ndim=1, mode='c'] components
        if value is None:
            components = np.array((w, x, y, z))
        else:
            assert np.shape(value) == (4,), "Component array must have length 4"
            components = np.asarray(value, dtype=np.float, order='C')
        self._components = components
        self.value = <quaternion *> &components[0]

    # Properties

    property w:
        """ Real component of this quaternion. """
        def __get__(self):
            return self.value.w
        def __set__(self, double value):
            self.value.w = value

    property x:
        """ First imaginary component of this quaternion. """
        def __get__(self):
            return self.value.x
        def __set__(self, double value):
            self.value.x = value

    property y:
        """ Second imaginary component of this quaternion. """
        def __get__(self):
            return self.value.y
        def __set__(self, double value):
            self.value.y = value

    property z:
        """ Third imaginary component of this quaternion. """
        def __get__(self):
            return self.value.z
        def __set__(self, double value):
            self.value.z = value

    property components:
        """ The (w, x, y, z) components of this quaternion as a tuple. """
        def __get__(Quaternion self):
            return self._components
        def __set__(Quaternion self, value):
            self._components[:] = value

    property magnitude:
        """ Magnitude of this quaternion. """
        def __get__(Quaternion self):
            return quaternion_magnitude(self.value)

    property conjugate:
        """ Conjugate of this quaternion. """
        def __get__(Quaternion self):
            result = Quaternion()
            quaternion_conjugate(self.value, result.value)
            return result

    property vector:
        """ The (x, y, z) components of this quaternion as a tuple. """
        def __get__(Quaternion self):
            return self.components[1:]
        def __set__(Quaternion self, value):
            self.components[1:] = value

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
            return quaternion_less(self.value, other.value)
        elif op == 1:
            return quaternion_less_equal(self.value, other.value)
        elif op == 2:
            return quaternion_equal(self.value, other.value)
        elif op == 3:
            return quaternion_not_equal(self.value, other.value)
        else:
            return NotImplemented

    def __getstate__(self):
        return self.components

    def __setstate__(self, state):
        self.__init__(value=state)

    def __reduce__(Quaternion self):
        return (Quaternion, (self.w, self.x, self.y, self.z))

    # General arithmetic (creates new instances)

    def __neg__(Quaternion self):
        cdef Quaternion result = Quaternion()
        quaternion_negative(self.value, result.value)
        return result

    def __invert__(Quaternion self):
        cdef Quaternion result = Quaternion()
        quaternion_conjugate(self.value, result.value)
        return result

    def __add__(a, b):
        cdef Quaternion result = Quaternion()
        cdef Quaternion qa, qb
        cdef double sa, sb
        if isinstance(a, Quaternion) and isinstance(b, Quaternion):
            qa, qb = a, b
            quaternion_add(qa.value, qb.value, result.value)
        elif isinstance(a, Quaternion) and np.isscalar(b):
            qa, sb = a, b
            quaternion_add_scalar(qa.value, sb, result.value)
        elif np.isscalar(a) and isinstance(b, Quaternion):
            sa, qb = a, b
            quaternion_add_scalar(qb.value, sa, result.value)
        else:
            return NotImplemented
        return result

    def __sub__(a, b):
        cdef Quaternion result = Quaternion()
        cdef Quaternion qa, qb
        cdef double sa, sb
        if isinstance(a, Quaternion) and isinstance(b, Quaternion):
            qa, qb = a, b
            quaternion_subtract(qa.value, qb.value, result.value)
        elif isinstance(a, Quaternion) and np.isscalar(b):
            qa, sb = a, b
            quaternion_subtract_scalar(qa.value, sb, result.value)
        elif np.isscalar(a) and isinstance(b, Quaternion):
            sa, qb = a, b
            quaternion_negative(qb.value, result.value)
            quaternion_add_scalar(result.value, sa, result.value)
        else:
            return NotImplemented
        return result

    def __mul__(a, b):
        cdef Quaternion result = Quaternion()
        cdef Quaternion qa, qb
        cdef double sa, sb
        if isinstance(a, Quaternion) and isinstance(b, Quaternion):
            qa, qb = a, b
            quaternion_multiply(qa.value, qb.value, result.value)
        elif isinstance(a, Quaternion) and np.isscalar(b):
            qa, sb = a, b
            quaternion_multiply_scalar(qa.value, sb, result.value)
        elif np.isscalar(a) and isinstance(b, Quaternion):
            sa, qb = a, b
            quaternion_multiply_scalar(qb.value, sa, result.value)
        else:
            return NotImplemented
        return result

    def __div__(a, b):
        cdef Quaternion result = Quaternion()
        cdef Quaternion qa, qb
        cdef double sa, sb
        if isinstance(a, Quaternion) and isinstance(b, Quaternion):
            qa, qb = a, b
            quaternion_divide(qa.value, qb.value, result.value)
        elif isinstance(a, Quaternion) and np.isscalar(b):
            qa, sb = a, b
            quaternion_divide_scalar(qa.value, sb, result.value)
        elif np.isscalar(a) and isinstance(b, Quaternion):
            sa, qb = a, b
            quaternion_power_scalar(qb.value, -1, result.value)
            quaternion_multiply_scalar(result.value, sa, result.value)
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
            quaternion_power(qa.value, qb.value, result.value)
        elif isinstance(a, Quaternion) and np.isscalar(b):
            qa, sb = a, b
            quaternion_power_scalar(qa.value, sb, result.value)
        elif np.isscalar(a) and isinstance(b, Quaternion):
            qa, qb = Quaternion(a, 0, 0, 0), b
            quaternion_power(qa.value, qb.value, result.value)
        else:
            return NotImplemented
        return result

    def log(Quaternion self):
        """ Natural logarithm of this quaternion. """
        cdef Quaternion result = Quaternion()
        quaternion_log(self.value, result.value)
        return result

    def exp(Quaternion self):
        """ Exponential of this quaternion. """
        cdef Quaternion result = Quaternion()
        quaternion_exp(self.value, result.value)
        return result

    def dot(Quaternion self, Quaternion other):
        """ Dot prodct of this quaternion with another. """
        return quaternion_dot(self.value, other.value)

    # In-place arithmetic (modifies existing instance)

    def __iadd__(Quaternion self, other):
        cdef Quaternion qo
        if isinstance(other, Quaternion):
            qo = other
            quaternion_add(self.value, qo.value, self.value)
        else:
            quaternion_add_scalar(self.value, other, self.value)
        return self

    def __isub__(Quaternion self, other):
        cdef Quaternion qo
        if isinstance(other, Quaternion):
            qo = other
            quaternion_subtract(self.value, qo.value, self.value)
        else:
            quaternion_subtract_scalar(self.value, other, self.value)
        return self

    def __imul__(Quaternion self, other):
        cdef Quaternion qo
        if isinstance(other, Quaternion):
            qo = other
            quaternion_multiply(self.value, qo.value, self.value)
        else:
            quaternion_multiply_scalar(self.value, other, self.value)
        return self

    def __idiv__(Quaternion self, other):
        cdef Quaternion qo
        if isinstance(other, Quaternion):
            qo = other
            quaternion_divide(self.value, qo.value, self.value)
        else:
            quaternion_divide_scalar(self.value, other, self.value)
        return self

    __itruediv__ = __idiv__

    def __ipow__(Quaternion self, other):
        cdef Quaternion qo
        if isinstance(other, Quaternion):
            qo = other
            quaternion_power(self.value, qo.value, self.value)
        else:
            quaternion_power_scalar(self.value, other, self.value)
        return self

    def normalise(Quaternion self):
        """ Normalise this quaternion to unit length. """
        quaternion_multiply_scalar(self.value,
            1 / quaternion_magnitude(self.value), self.value)
        return self

    def negate(Quaternion self):
        """ Negate this quaternion. """
        quaternion_negative(self.value, self.value)
        return self

    # Rotations on vectors.

    def rotate_vector(Quaternion self, v):
        """
        Use this quaternion to rotate a vector.
        """
        cdef np.ndarray[np.float_t, ndim=1, mode='c'] vec1d
        cdef np.ndarray[np.float_t, ndim=1, mode='c'] res1d
        cdef np.ndarray[np.float_t, ndim=2, mode='c'] vec2d
        cdef np.ndarray[np.float_t, ndim=2, mode='c'] res2d
        if np.ndim(v) == 2:
            vec2d = np.asarray(v, dtype=np.float)
            assert vec2d.shape[1] == 3, "Vectors must have shape (N, 3)"
            res2d = np.empty_like(vec2d)
            for i in range(len(v)):
                quaternion_rotate_vector(self.value,
                    <double *> &vec2d[i, 0], <double *> &res2d[i, 0])
            return res2d
        else:
            assert np.shape(v) == (3,), "Vector must have length 3"
            vec1d = np.asarray(v, dtype=float)
            res1d = np.empty_like(vec1d)
            quaternion_rotate_vector(self.value,
                <double *> &vec1d[0], <double *> &res1d[0])
            return res1d

    def rotate_frame(Quaternion self, v):
        """
        Use this quaternion to rotate the co-ordinate frame of a vector.
        """
        cdef np.ndarray[np.float_t, ndim=1, mode='c'] vec1d
        cdef np.ndarray[np.float_t, ndim=1, mode='c'] res1d
        cdef np.ndarray[np.float_t, ndim=2, mode='c'] vec2d
        cdef np.ndarray[np.float_t, ndim=2, mode='c'] res2d
        if np.ndim(v) == 2:
            vec2d = np.asarray(v, dtype=np.float)
            assert vec2d.shape[1] == 3, "Vectors must have shape (N, 3)"
            res2d = np.empty_like(vec2d)
            for i in range(len(v)):
                quaternion_rotate_frame(self.value,
                    <double *> &vec2d[i, 0], <double *> &res2d[i, 0])
            return res2d
        else:
            assert np.shape(v) == (3,), "Vector must have length 3"
            vec1d = np.asarray(v, dtype=np.float)
            res1d = np.empty_like(vec1d)
            quaternion_rotate_frame(self.value,
                <double *> &vec1d[0], <double *> &res1d[0])
            return res1d

    # Conversions to and from other rotation formats.

    @classmethod
    def from_euler(cls, char *order, angles):
        cdef Quaternion result = Quaternion()
        assert len(order) == len(angles), "Order and angles must have same length"
        cdef np.ndarray[np.float_t, ndim=1, mode='c'] ang = np.asarray(angles, dtype=np.float)
        quaternion_from_euler(order, <double *> &ang[0], result.value)
        return result

    def to_euler(Quaternion self, char *order):
        cdef np.ndarray[np.float_t, ndim=1, mode='c'] result = np.zeros(len(order))
        quaternion_to_euler(self.value, order, <double *> &result[0])
        return result

cdef class QuaternionArray:

    cdef quaternion *values
    cdef np.ndarray _components
    cdef unsigned int length

    def __cinit__(QuaternionArray self, values):
        cdef np.ndarray[np.float_t, ndim=2, mode='c'] components
        if np.ndim(values) == 2:
            assert np.shape(values)[1] == 4, \
                "Expected components of shape (N, 4) or a sequence of quaternions"
            components = np.asarray(values, dtype=np.float, order='C')
        else:
            assert np.ndim(values) == 1, \
                "Expected components of shape (N, 4) or a sequence of quaternions"
            length = len(values)
            components = np.empty((length, 4))
            for i in range(length):
                components[i] = values[i].components
        self._components = components
        self.values = <quaternion *> &components[0, 0]
        self.length = len(components)

    # Sequence operations

    def __len__(QuaternionArray self):
        return self.length

    def __getitem__(QuaternionArray self, key):
        if isinstance(key, slice):
            return QuaternionArray(self.components[key])
        else:
            return Quaternion(value=self.components[key])
    
    def __setitem__(QuaternionArray self, key, value):
        if isinstance(value, (Quaternion, QuaternionArray)):
            self.components[key] = value.components
        else:
            self.components[key] = value

    # Properties

    property components:
        """ N x 4 array of quaternion [w, x, y, z] components. """
        def __get__(QuaternionArray self):
            cdef np.ndarray[np.float_t, ndim=2, mode='c'] components = self._components
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
            cdef np.ndarray[np.float_t, ndim=1, mode='c'] result = np.empty(self.length)
            for i in range(self.length):
                result[i] = quaternion_magnitude(&self.values[i])
            return result

    property conjugate:
        """ Conjugates of these quaternions. """
        def __get__(QuaternionArray self):
            cdef QuaternionArray result = QuaternionArray.like(self)
            for i in range(self.length):
                quaternion_conjugate(&self.values[i], &result.values[i])
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
    def empty(cls, length):
        """ Create an empty QuaternionArray with a given length. """
        return QuaternionArray(np.empty((length, 4)))

    @classmethod
    def like(cls, ref):
        """ Create an empty QuaternionArray with the length of an existing sequence. """
        return QuaternionArray.empty(len(ref))

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
            quaternion_negative(&self.values[i], &result.values[i])
        return result

    def __invert__(QuaternionArray self):
        cdef QuaternionArray result = QuaternionArray.like(self)
        for i in range(self.length):
            quaternion_conjugate(&self.values[i], &result.values[i])
        return result

    def __add__(a, b):
        cdef QuaternionArray result
        cdef QuaternionArray qqa, qqb
        cdef Quaternion qa, qb
        cdef np.ndarray[np.float_t, ndim=1, mode='c'] ssa, ssb
        cdef double sa, sb
        if isinstance(a, QuaternionArray) and isinstance(b, QuaternionArray):
            qqa, qqb = a, b
            assert qqa.length == qqb.length, "Arrays must have equal length"
            result = QuaternionArray.like(qqa)
            for i in range(qqa.length):
                quaternion_add(&qqa.values[i], &qqb.values[i], &result.values[i])
        elif isinstance(a, QuaternionArray):
            qqa = a
            result = QuaternionArray.like(qqa)
            if isinstance(b, Quaternion):
                qb = b
                for i in range(qqa.length):
                    quaternion_add(&qqa.values[i], qb.value, &result.values[i])
            elif np.shape(b) == (qqa.length,):
                ssb = np.asarray(b)
                for i in range(qqa.length):
                    quaternion_add_scalar(&qqa.values[i], ssb[i], &result.values[i])
            else:
                sb = b
                for i in range(qqa.length):
                    quaternion_add_scalar(&qqa.values[i], sb, &result.values[i])
        elif isinstance(b, QuaternionArray):
            qqb = b
            result = QuaternionArray.like(qqb)
            if isinstance(a, Quaternion):
                qa = a
                for i in range(qqb.length):
                    quaternion_add(qa.value, &qqb.values[i], &result.values[i])
            elif np.shape(a) == (qqb.length,):
                ssa = np.asarray(a)
                for i in range(qqb.length):
                    quaternion_add_scalar(&qqb.values[i], ssa[i], &result.values[i])
            else:
                sa = a
                for i in range(qqb.length):
                    quaternion_add_scalar(&qqb.values[i], sa, &result.values[i])
        else:
            return NotImplemented
        return result

    def __sub__(a, b):
        cdef QuaternionArray result
        cdef QuaternionArray qqa, qqb
        cdef Quaternion qa, qb
        cdef np.ndarray[np.float_t, ndim=1, mode='c'] ssa, ssb
        cdef double sa, sb
        if isinstance(a, QuaternionArray) and isinstance(b, QuaternionArray):
            qqa, qqb = a, b
            assert qqa.length == qqb.length, "Arrays must have equal length"
            result = QuaternionArray.like(qqa)
            for i in range(qqa.length):
                quaternion_subtract(&qqa.values[i], &qqb.values[i], &result.values[i])
        elif isinstance(a, QuaternionArray):
            qqa = a
            result = QuaternionArray.like(qqa)
            if isinstance(b, Quaternion):
                qb = b
                for i in range(qqa.length):
                    quaternion_subtract(&qqa.values[i], qb.value, &result.values[i])
            elif np.shape(b) == (qqa.length,):
                ssb = np.asarray(b)
                for i in range(qqa.length):
                    quaternion_subtract_scalar(&qqa.values[i], ssb[i], &result.values[i])
            else:
                sb = b
                for i in range(qqa.length):
                    quaternion_subtract_scalar(&qqa.values[i], sb, &result.values[i])
        elif isinstance(b, QuaternionArray):
            qqb = b
            result = QuaternionArray.like(qqb)
            if isinstance(a, Quaternion):
                qa = a
                for i in range(qqb.length):
                    quaternion_subtract(qa.value, &qqb.values[i], &result.values[i])
            elif np.shape(a) == (qqb.length,):
                ssa = np.asarray(a)
                for i in range(qqb.length):
                    quaternion_negative(&qqb.values[i], &result.values[i])
                    quaternion_add_scalar(&result.values[i], ssa[i], &result.values[i])
            else:
                sa = a
                for i in range(qqb.length):
                    quaternion_negative(&qqb.values[i], &result.values[i])
                    quaternion_add_scalar(&result.values[i], sa, &result.values[i])
        else:
            return NotImplemented
        return result

    def __mul__(a, b):
        cdef QuaternionArray result
        cdef QuaternionArray qqa, qqb
        cdef Quaternion qa, qb
        cdef np.ndarray[np.float_t, ndim=1, mode='c'] ssa, ssb
        cdef double sa, sb
        if isinstance(a, QuaternionArray) and isinstance(b, QuaternionArray):
            qqa, qqb = a, b
            assert qqa.length == qqb.length, "Arrays must have equal length"
            result = QuaternionArray.like(qqa)
            for i in range(qqa.length):
                quaternion_multiply(&qqa.values[i], &qqb.values[i], &result.values[i])
        elif isinstance(a, QuaternionArray):
            qqa = a
            result = QuaternionArray.like(qqa)
            if isinstance(b, Quaternion):
                qb = b
                for i in range(qqa.length):
                    quaternion_multiply(&qqa.values[i], qb.value, &result.values[i])
            elif np.shape(b) == (qqa.length,):
                ssb = np.asarray(b)
                for i in range(qqa.length):
                    quaternion_multiply_scalar(&qqa.values[i], ssb[i], &result.values[i])
            else:
                sb = b
                for i in range(qqa.length):
                    quaternion_multiply_scalar(&qqa.values[i], sb, &result.values[i])
        elif isinstance(b, QuaternionArray):
            qqb = b
            result = QuaternionArray.like(qqb)
            if isinstance(a, Quaternion):
                qa = a
                for i in range(qqb.length):
                    quaternion_multiply(qa.value, &qqb.values[i], &result.values[i])
            elif np.shape(a) == (qqb.length,):
                ssa = np.asarray(a)
                for i in range(qqb.length):
                    quaternion_multiply_scalar(&qqb.values[i], ssa[i], &result.values[i])
            else:
                sa = a
                for i in range(qqb.length):
                    quaternion_multiply_scalar(&qqb.values[i], sa, &result.values[i])
        else:
            return NotImplemented
        return result

    def __div__(a, b):
        cdef QuaternionArray result
        cdef QuaternionArray qqa, qqb
        cdef Quaternion qa, qb
        cdef np.ndarray[np.float_t, ndim=1, mode='c'] ssa, ssb
        cdef double sa, sb
        if isinstance(a, QuaternionArray) and isinstance(b, QuaternionArray):
            qqa, qqb = a, b
            assert qqa.length == qqb.length, "Arrays must have equal length"
            result = QuaternionArray.like(qqa)
            for i in range(qqa.length):
                quaternion_divide(&qqa.values[i], &qqb.values[i], &result.values[i])
        elif isinstance(a, QuaternionArray):
            qqa = a
            result = QuaternionArray.like(qqa)
            if isinstance(b, Quaternion):
                qb = b
                for i in range(qqa.length):
                    quaternion_divide(&qqa.values[i], qb.value, &result.values[i])
            elif np.shape(b) == (qqa.length,):
                ssb = np.asarray(b)
                for i in range(qqa.length):
                    quaternion_divide_scalar(&qqa.values[i], ssb[i], &result.values[i])
            else:
                sb = b
                for i in range(qqa.length):
                    quaternion_divide_scalar(&qqa.values[i], sb, &result.values[i])
        elif isinstance(b, QuaternionArray):
            qqb = b
            result = QuaternionArray.like(qqb)
            if isinstance(a, Quaternion):
                qa = a
                for i in range(qqb.length):
                    quaternion_divide(qa.value, &qqb.values[i], &result.values[i])
            elif np.shape(a) == (qqb.length,):
                ssa = np.asarray(a)
                for i in range(qqb.length):
                    quaternion_power_scalar(&qqb.values[i], -1, &result.values[i])
                    quaternion_multiply_scalar(&result.values[i], ssa[i], &result.values[i])
            else:
                sa = a
                for i in range(qqb.length):
                    quaternion_power_scalar(&qqb.values[i], -1, &result.values[i])
                    quaternion_multiply_scalar(&result.values[i], sa, &result.values[i])
        else:
            return NotImplemented
        return result

    __truediv__ = __div__

    def __pow__(a, b, c):
        cdef QuaternionArray result
        cdef QuaternionArray qqa, qqb
        cdef Quaternion qa, qb
        cdef np.ndarray[np.float_t, ndim=1, mode='c'] ssb
        cdef double sa, sb
        if isinstance(a, QuaternionArray) and isinstance(b, QuaternionArray):
            qqa, qqb = a, b
            assert qqa.length == qqb.length, "Arrays must have equal length"
            result = QuaternionArray.like(qqa)
            for i in range(qqa.length):
                quaternion_power(&qqa.values[i], &qqb.values[i], &result.values[i])
        elif isinstance(a, QuaternionArray):
            qqa = a
            result = QuaternionArray.like(qqa)
            if isinstance(b, Quaternion):
                qb = b
                for i in range(qqa.length):
                    quaternion_power(&qqa.values[i], qb.value, &result.values[i])
            elif np.shape(b) == (qqa.length,):
                ssb = np.asarray(b)
                for i in range(qqa.length):
                    quaternion_power_scalar(&qqa.values[i], ssb[i], &result.values[i])
            else:
                sb = b
                for i in range(qqa.length):
                    quaternion_power_scalar(&qqa.values[i], sb, &result.values[i])
        elif isinstance(b, QuaternionArray):
            qqb = b
            result = QuaternionArray.like(qqb)
            if isinstance(a, Quaternion):
                qa = a
                for i in range(qqb.length):
                    quaternion_power(qa.value, &qqb.values[i], &result.values[i])
            elif np.shape(a) == (qqb.length,):
                qqa = QuaternionArray.like(qqb)
                qqa.w = a
                qqa.vector = 0
                for i in range(qqb.length):
                    quaternion_power(&qqa.values[i], &qqb.values[i], &result.values[i])
            else:
                qa = Quaternion(a, 0, 0, 0)
                for i in range(qqb.length):
                    quaternion_power(qa.value, &qqb.values[i], &result.values[i])
        else:
            return NotImplemented
        return result

    def log(QuaternionArray self):
        """ Natural logarithms of these quaternions. """
        cdef QuaternionArray result = QuaternionArray.like(self)
        for i in range(self.length):
            quaternion_log(&self.values[i], &result.values[i])
        return result

    def exp(QuaternionArray self):
        """ Exponentials of these quaternions. """
        cdef QuaternionArray result = QuaternionArray.like(self)
        for i in range(self.length):
            quaternion_exp(&self.values[i], &result.values[i])
        return result

    def dot(QuaternionArray self, other):
        """ Dot prodct of these quaternions with another quaternion or quaternion array. """
        cdef np.ndarray[np.float_t, ndim=1, mode='c'] result = np.empty(self.length)
        cdef QuaternionArray qqo
        cdef Quaternion qo
        if isinstance(other, QuaternionArray):
            qqo = other
            for i in range(self.length):
                result[i] = quaternion_dot(&self.values[i], &qqo.values[i])
        elif isinstance(other, Quaternion):
            qo = other
            for i in range(self.length):
                result[i] = quaternion_dot(&self.values[i], qo.value)
        else:
            return NotImplemented
        return result

    # In-place arithmetic (modifies existing instance)

    def __iadd__(QuaternionArray self, other):
        cdef QuaternionArray qqo
        cdef Quaternion qo
        cdef np.ndarray[np.float_t, ndim=1, mode='c'] sso
        if isinstance(other, QuaternionArray):
            qqo = other
            for i in range(self.length):
                quaternion_add(&self.values[i], &qqo.values[i], &self.values[i])
        elif isinstance(other, Quaternion):
            qo = other
            for i in range(self.length):
                quaternion_add(&self.values[i], qo.value, &self.values[i])
        elif np.shape(other) == (self.length,):
            sso = np.asarray(other)
            for i in range(self.length):
                quaternion_add_scalar(&self.values[i], sso[i], &self.values[i])
        else:
            for i in range(self.length):
                quaternion_add_scalar(&self.values[i], other, &self.values[i])
        return self

    def __isub__(QuaternionArray self, other):
        cdef QuaternionArray qqo
        cdef Quaternion qo
        cdef np.ndarray[np.float_t, ndim=1, mode='c'] sso
        if isinstance(other, QuaternionArray):
            qqo = other
            for i in range(self.length):
                quaternion_subtract(&self.values[i], &qqo.values[i], &self.values[i])
        elif isinstance(other, Quaternion):
            qo = other
            for i in range(self.length):
                quaternion_subtract(&self.values[i], qo.value, &self.values[i])
        elif np.shape(other) == (self.length,):
            sso = np.asarray(other)
            for i in range(self.length):
                quaternion_subtract_scalar(&self.values[i], sso[i], &self.values[i])
        else:
            for i in range(self.length):
                quaternion_subtract_scalar(&self.values[i], other, &self.values[i])
        return self

    def __imul__(QuaternionArray self, other):
        cdef QuaternionArray qqo
        cdef Quaternion qo
        cdef np.ndarray[np.float_t, ndim=1, mode='c'] sso
        if isinstance(other, QuaternionArray):
            qqo = other
            for i in range(self.length):
                quaternion_multiply(&self.values[i], &qqo.values[i], &self.values[i])
        elif isinstance(other, Quaternion):
            qo = other
            for i in range(self.length):
                quaternion_multiply(&self.values[i], qo.value, &self.values[i])
        elif np.shape(other) == (self.length,):
            sso = np.asarray(other)
            for i in range(self.length):
                quaternion_multiply_scalar(&self.values[i], sso[i], &self.values[i])
        else:
            for i in range(self.length):
                quaternion_multiply_scalar(&self.values[i], other, &self.values[i])
        return self

    def __idiv__(QuaternionArray self, other):
        cdef QuaternionArray qqo
        cdef Quaternion qo
        cdef np.ndarray[np.float_t, ndim=1, mode='c'] sso
        if isinstance(other, QuaternionArray):
            qqo = other
            for i in range(self.length):
                quaternion_divide(&self.values[i], &qqo.values[i], &self.values[i])
        elif isinstance(other, Quaternion):
            qo = other
            for i in range(self.length):
                quaternion_divide(&self.values[i], qo.value, &self.values[i])
        elif np.shape(other) == (self.length,):
            sso = np.asarray(other)
            for i in range(self.length):
                quaternion_divide_scalar(&self.values[i], sso[i], &self.values[i])
        else:
            for i in range(self.length):
                quaternion_divide_scalar(&self.values[i], other, &self.values[i])
        return self

    __itruediv__ = __idiv__

    def __ipow__(QuaternionArray self, other):
        cdef QuaternionArray qqo
        cdef Quaternion qo
        cdef np.ndarray[np.float_t, ndim=1, mode='c'] sso
        if isinstance(other, QuaternionArray):
            qqo = other
            for i in range(self.length):
                quaternion_power(&self.values[i], &qqo.values[i], &self.values[i])
        elif isinstance(other, Quaternion):
            qo = other
            for i in range(self.length):
                quaternion_power(&self.values[i], qo.value, &self.values[i])
        elif np.shape(other) == (self.length,):
            sso = np.asarray(other)
            for i in range(self.length):
                quaternion_power_scalar(&self.values[i], sso[i], &self.values[i])
        else:
            for i in range(self.length):
                quaternion_power_scalar(&self.values[i], other, &self.values[i])
        return self

    def normalise(QuaternionArray self):
        """ Normalise these quaternions to unit length. """
        for i in range(self.length):
            quaternion_multiply_scalar(&self.values[i],
                1 / quaternion_magnitude(&self.values[i]), &self.values[i])
        return self

    def negate(Quaternion self):
        """ Negate these quaternions. """
        for i in range(self.length):
            quaternion_negative(&self.values[i], &self.values[i])
        return self

    def rotate_vector(QuaternionArray self, v):
        """
        Use these quaternions to rotate a vector or array of vectors.
        """
        cdef np.ndarray[np.float_t, ndim=1, mode='c'] vec1d
        cdef np.ndarray[np.float_t, ndim=2, mode='c'] vec2d
        cdef np.ndarray[np.float_t, ndim=2, mode='c'] result = np.empty((self.length, 3))
        if np.ndim(v) == 2:
            assert np.shape(v) == (self.length, 3), \
                "Vectors must have shape (N, 3) matching array length"
            vec2d = np.asarray(v, dtype=np.float)
            for i in range(self.length):
                quaternion_rotate_vector(&self.values[i],
                    <double *> &vec2d[i, 0], <double *> &result[i, 0])
        else:
            assert np.shape(v) == (3,), "Vector must have length 3"
            vec1d = np.asarray(v, dtype=float)
            for i in range(self.length):
                quaternion_rotate_vector(&self.values[i],
                    <double *> &vec1d[0], <double *> &result[i, 0])
        return result

    def rotate_frame(QuaternionArray self, v):
        """
        Use these quaternions to rotate the co-ordinate frame of a vector or array of vectors.
        """
        cdef np.ndarray[np.float_t, ndim=1, mode='c'] vec1d
        cdef np.ndarray[np.float_t, ndim=2, mode='c'] vec2d
        cdef np.ndarray[np.float_t, ndim=2, mode='c'] result = np.empty((self.length, 3))
        if np.ndim(v) == 2:
            assert np.shape(v) == (self.length, 3), \
                "Vectors must have shape (N, 3) matching array length"
            vec2d = np.asarray(v, dtype=np.float)
            for i in range(self.length):
                quaternion_rotate_frame(&self.values[i],
                    <double *> &vec2d[i, 0], <double *> &result[i, 0])
        else:
            assert np.shape(v) == (3,), "Vector must have length 3"
            vec1d = np.asarray(v, dtype=float)
            for i in range(self.length):
                quaternion_rotate_frame(&self.values[i],
                    <double *> &vec1d[0], <double *> &result[i, 0])
        return result

    # Conversions to and from other rotation formats.

    @classmethod
    def from_euler(cls, char *order, angles):
        shape = np.shape(angles)
        assert len(shape) == 2 and shape[1] == len(order), \
            "Angles must have shape (N, len(order))"
        cdef QuaternionArray result = QuaternionArray.empty(shape[0])
        cdef np.ndarray[np.float_t, ndim=2, mode='c'] ang = np.asarray(angles, dtype=np.float)
        for i in range(len(angles)):
            quaternion_from_euler(order, <double *> &ang[i, 0], &result.values[i])
        return result

    def to_euler(QuaternionArray self, char *order):
        cdef np.ndarray[np.float_t, ndim=2, mode='c'] result = np.zeros((self.length, len(order)))
        for i in range(self.length):
            quaternion_to_euler(&self.values[i], order, <double *> &result[i, 0])
        return result
