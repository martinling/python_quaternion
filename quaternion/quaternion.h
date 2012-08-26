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

#ifndef __QUATERNION_H__
#define __QUATERNION_H__

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
	double w;
	double x;
	double y;
	double z;
} quaternion;

int quaternion_isnonzero(const quaternion *q);
int quaternion_isnan(const quaternion *q);
int quaternion_isinf(const quaternion *q);
int quaternion_isfinite(const quaternion *q);
int quaternion_equal(const quaternion *q1, const quaternion *q2);
int quaternion_not_equal(const quaternion *q1, const quaternion *q2);
int quaternion_less(const quaternion *q1, const quaternion *q2);
int quaternion_less_equal(const quaternion *q1, const quaternion *q2);
double quaternion_magnitude(const quaternion *q);
double quaternion_dot(const quaternion *q1, const quaternion *q2);
void quaternion_negative(const quaternion *q, quaternion *r);
void quaternion_conjugate(const quaternion *q, quaternion *r);
void quaternion_copysign(const quaternion *q1, const quaternion *q2, quaternion *r);
void quaternion_add(const quaternion *q1, const quaternion *q2, quaternion *r);
void quaternion_add_scalar(const quaternion *q, double s, quaternion *r);
void quaternion_subtract(const quaternion *q1, const quaternion *q2, quaternion *r);
void quaternion_subtract_scalar(const quaternion *q, double s, quaternion *r);
void quaternion_multiply(const quaternion *q1, const quaternion *q2, quaternion *r);
void quaternion_multiply_scalar(const quaternion *q, double s, quaternion *r);
void quaternion_divide(const quaternion *q1, const quaternion *q2, quaternion *r);
void quaternion_divide_scalar(const quaternion *q, double s, quaternion *r);
void quaternion_log(const quaternion *q, quaternion *r);
void quaternion_exp(const quaternion *q, quaternion *r);
void quaternion_power(const quaternion *q1, const quaternion *q2, quaternion *r);
void quaternion_power_scalar(const quaternion *q, double s, quaternion *r);
void quaternion_rotate_vector(const quaternion *q, const double v[3], double r[3]);
void quaternion_rotate_frame(const quaternion *q, const double v[3], double r[3]);
void quaternion_from_euler(const char *order, double angles[3], quaternion *r);
void quaternion_to_euler(const quaternion *q, const char *order, double r[3]);

#ifdef __cplusplus
}
#endif

#endif
