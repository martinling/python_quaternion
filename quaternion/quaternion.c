/*
 * Quaternion math routines
 *
 * Copyright (c) 2011-2012 Martin Ling <martin@earth.li>
 *
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 *     * Redistributions of source code must retain the above copyright
 *       notice, this list of conditions and the following disclaimer.
 *
 *     * Redistributions in binary form must reproduce the above copyright
 *       notice, this list of conditions and the following disclaimer in the
 *       documentation and/or other materials provided with the distribution.
 *
 *     * The name of the author(s) may not be used to endorse or promote
 *       products derived from this software without specific prior written
 *       permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTERS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 */

#include "quaternion.h"
#include "math.h"
#include "ctype.h"

#define _QUAT_EPS 1e-6

int
quaternion_isnonzero(const quaternion *q)
{
    return q->w != 0 && q->x != 0 && q->y != 0 && q->z != 0;
}

int
quaternion_isnan(const quaternion *q)
{
    return isnan(q->w) || isnan(q->x) || isnan(q->y) || isnan(q->z);
}

int
quaternion_isinf(const quaternion *q)
{
    return isinf(q->w) || isinf(q->x) || isinf(q->y) || isnan(q->z);
}

int
quaternion_isfinite(const quaternion *q)
{
    return isfinite(q->w) && isfinite(q->x) && isfinite(q->y) && isfinite(q->z);
}

int
quaternion_equal(const quaternion *q1, const quaternion *q2)
{
    return
        !quaternion_isnan(q1) &&
        !quaternion_isnan(q2) &&
        q1->w == q2->w &&
        q1->x == q2->x &&
        q1->y == q2->y &&
        q1->z == q2->z;
}

int
quaternion_not_equal(const quaternion *q1, const quaternion *q2)
{
    return !quaternion_equal(q1, q2);
}

int
quaternion_less(const quaternion *q1, const quaternion *q2)
{
    return
        (!quaternion_isnan(q1) &&
        !quaternion_isnan(q2)) && (
            q1->w != q2->w ? q1->w < q2->w :
            q1->x != q2->x ? q1->x < q2->x :
            q1->y != q2->y ? q1->y < q2->y :
            q1->z != q2->z ? q1->z < q2->z : 0);
}

int
quaternion_less_equal(const quaternion *q1, const quaternion *q2)
{
   return
        (!quaternion_isnan(q1) &&
        !quaternion_isnan(q2)) && (
            q1->w != q2->w ? q1->w < q2->w :
            q1->x != q2->x ? q1->x < q2->x :
            q1->y != q2->y ? q1->y < q2->y :
            q1->z != q2->z ? q1->z < q2->z : 1);
}

double
quaternion_magnitude(const quaternion *q)
{
   return sqrt(q->w*q->w + q->x*q->x + q->y*q->y + q->z*q->z);
}

double
quaternion_dot(const quaternion *q1, const quaternion *q2)
{
   return q1->w*q2->w + q1->x*q2->x + q1->y*q2->y + q1->z*q2->z;
}

void
quaternion_negative(const quaternion *q, quaternion *r)
{
   r->w = -q->w,
   r->x = -q->x,
   r->y = -q->y,
   r->z = -q->z;
}

void
quaternion_conjugate(const quaternion *q, quaternion *r)
{
   r->w = q->w,
   r->x = -q->x,
   r->y = -q->y,
   r->z = -q->z;
}

void
quaternion_copysign(const quaternion *q1, const quaternion *q2, quaternion *r)
{
   r->w = copysign(q1->w, q2->w),
   r->x = copysign(q1->x, q2->x),
   r->y = copysign(q1->y, q2->y),
   r->z = copysign(q1->z, q2->z);
}

void
quaternion_add(const quaternion *q1, const quaternion *q2, quaternion *r)
{
   r->w = q1->w+q2->w,
   r->x = q1->x+q2->x,
   r->y = q1->y+q2->y,
   r->z = q1->z+q2->z;
}

void
quaternion_add_scalar(const quaternion *q, double s, quaternion *r)
{
   r->w = q->w+s,
   r->x = q->x,
   r->y = q->y,
   r->z = q->z;
}

void
quaternion_subtract(const quaternion *q1, const quaternion *q2, quaternion *r)
{
   r->w = q1->w-q2->w,
   r->x = q1->x-q2->x,
   r->y = q1->y-q2->y,
   r->z = q1->z-q2->z;
}

void
quaternion_subtract_scalar(const quaternion *q, double s, quaternion *r)
{
   r->w = q->w-s,
   r->x = q->x,
   r->y = q->y,
   r->z = q->z;
}

void
quaternion_multiply(const quaternion *q1, const quaternion *q2, quaternion *r)
{
   r->w = q1->w*q2->w - q1->x*q2->x - q1->y*q2->y - q1->z*q2->z,
   r->x = q1->w*q2->x + q1->x*q2->w + q1->y*q2->z - q1->z*q2->y,
   r->y = q1->w*q2->y - q1->x*q2->z + q1->y*q2->w + q1->z*q2->x,
   r->z = q1->w*q2->z + q1->x*q2->y - q1->y*q2->x + q1->z*q2->w;
}

void
quaternion_multiply_scalar(const quaternion *q, double s, quaternion *r)
{
   r->w = s*q->w,
   r->x = s*q->x,
   r->y = s*q->y,
   r->z = s*q->z;
}

void
quaternion_divide(const quaternion *q1, const quaternion *q2, quaternion *r)
{
   double s = q2->w*q2->w + q2->x*q2->x + q2->y*q2->y + q2->z*q2->z;
   r->w = (  q1->w*q2->w + q1->x*q2->x + q1->y*q2->y + q1->z*q2->z) / s,
   r->x = (- q1->w*q2->x + q1->x*q2->w + q1->y*q2->z - q1->z*q2->y) / s,
   r->y = (- q1->w*q2->y - q1->x*q2->z + q1->y*q2->w + q1->z*q2->x) / s,
   r->z = (- q1->w*q2->z + q1->x*q2->y - q1->y*q2->x + q1->z*q2->w) / s;
}

void
quaternion_divide_scalar(const quaternion *q, double s, quaternion *r)
{
   r->w = q->w/s,
   r->x = q->x/s,
   r->y = q->y/s,
   r->z = q->z/s;
}

void
quaternion_log(const quaternion *q, quaternion *r)
{
   double sumvsq = q->x*q->x + q->y*q->y + q->z*q->z;
   double vnorm = sqrt(sumvsq);
   if (vnorm > _QUAT_EPS) {
      double m = sqrt(q->w*q->w + sumvsq);
      double s = acos(q->w/m) / vnorm;
      r->w = log(m),
      r->x = s*q->x,
      r->y = s*q->y,
      r->z = s*q->z;
   } else {
      r->w = 0,
      r->x = 0,
      r->y = 0,
      r->z = 0;
   }
}

void
quaternion_exp(const quaternion *q, quaternion *r)
{
   double vnorm = sqrt(q->x*q->x + q->y*q->y + q->z*q->z);
   if (vnorm > _QUAT_EPS) {
      double s = sin(vnorm) / vnorm;
      double e = exp(q->w);
      r->w = e*cos(vnorm),
      r->x = e*s*q->x,
      r->y = e*s*q->y,
      r->z = e*s*q->z;
   } else {
      r->w = exp(q->w),
      r->x = 0,
      r->y = 0,
      r->z = 0;
   }
}

void
quaternion_power(const quaternion *q1, const quaternion *q2, quaternion *r)
{
   quaternion_log(q1, r);
   quaternion_multiply(r, q2, r);
   quaternion_exp(r, r);
}

void
quaternion_power_scalar(const quaternion *q, double s, quaternion *r)
{
   quaternion_log(q, r);
   quaternion_multiply_scalar(r, s, r);
   quaternion_exp(r, r);
}

void
quaternion_rotate_vector(const quaternion *q, const double v[3], double r[3])
{
   double W = q->x * v[0] + q->y * v[1] + q->z * v[2];
   double X = q->w * v[0] + q->y * v[2] - q->z * v[1];
   double Y = q->w * v[1] - q->x * v[2] + q->z * v[0];
   double Z = q->w * v[2] + q->x * v[1] - q->y * v[0];

   r[0] = W * q->x + X * q->w - Y * q->z + Z * q->y,
   r[1] = W * q->y + X * q->z + Y * q->w - Z * q->x,
   r[2] = W * q->z - X * q->y + Y * q->x + Z * q->w;
}

void
quaternion_rotate_frame(const quaternion *q, const double v[3], double r[3])
{
   double W = q->x * v[0] + q->y * v[1] - q->z * v[2];
   double X = q->w * v[0] - q->y * v[2] + q->z * v[1];
   double Y = q->w * v[1] + q->x * v[2] - q->z * v[0];
   double Z = q->w * v[2] - q->x * v[1] + q->y * v[0];

   r[0] = W * q->x + X * q->w + Y * q->z - Z * q->y,
   r[1] = W * q->y - X * q->z + Y * q->w + Z * q->x,
   r[2] = W * q->z + X * q->y - Y * q->x + Z * q->w;
}


static inline void
quaternion_from_angle(int axis, double angle, quaternion *r)
{
   int i;
   r->w = cos(angle) / 2;
   for (i = 0; i < 3; i++)
      ((double *) &r->x)[i] = i == axis ? sin(angle / 2) : 0;
}

void
quaternion_from_euler(const char *order, double angles[], quaternion *r)
{
   int i;
   quaternion q;
   for (i = 0; order[i] != '\0'; i++)
   {
      char axis = tolower(*order);
      if (axis < 'x' || axis > 'z') {
         r->w = NAN, 
         r->x = NAN,
         r->y = NAN,
         r->z = NAN;
         return;
      } else {
         if (i == 0) {
            quaternion_from_angle(axis - 'x', angles[i], r);
         } else {
            quaternion_from_angle(axis - 'x', angles[i], &q);
            quaternion_multiply(r, &q, r);
         }
      }
   }
}

void
quaternion_to_euler(const quaternion *q, const char *order, double r[3])
{
   const double axes[3][3] = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
   int indices[3], i;
   for (i = 0; i < 3; i++) {
      char axis = tolower(order[i]);
      if (axis < 'x' || axis > 'z') {
         r[0] = NAN;
         r[1] = NAN;
         r[2] = NAN;
         return;
      }
      indices[i] = axis - 'x';
   }
   double v[3];
   quaternion_rotate_vector(q, axes[indices[2]], v);
   int a = indices[0];
   int b = (indices[0] + 1) % 3;
   int c = (indices[0] + 2) % 3;
   int non_circular = indices[1] != b;
   int repeated_axis = indices[0] == indices[2];
   if (non_circular) {
      if (repeated_axis) {
         r[0] = atan2(v[c], v[b]);
         r[1] = acos(v[a]);
      } else {
         r[0] = atan2(v[c], v[b]);
         r[1] = -asin(v[a]);
      }
   } else {
      if (repeated_axis) {
         r[0] = atan2(v[b], -v[c]);
         r[1] = acos(a);
      } else {
         r[0] = atan2(-v[b], v[c]);
         r[1] = asin(v[a]);
      }
   }
   quaternion qr[3];
   for (i = 0; i < 2; i++)
      quaternion_from_angle(indices[i], r[i], &qr[i]);
   quaternion_multiply(&qr[0], &qr[1], &qr[2]);
   quaternion_conjugate(&qr[2], &qr[2]);
   quaternion_multiply(&qr[2], q, &qr[2]);
   r[2] = 2 * acos(qr[2].w);
   if (r[2] > M_PI)
      r[2] -= 2 * M_PI;
}
