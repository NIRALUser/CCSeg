#ifndef _FOURIER_FUNCS_H_
#define _FOURIER_FUNCS_H_

#include "CCsegtool_parameters.h"
#include <assert.h>
#include <memory.h>

typedef double Point4[4];
typedef double Point2[2];

/* Functions defined in, exported from fourier_funcs.c */
void sincos_table(int n_steps);
void parse_coefs(CCsegtool_parameters *parameters, int num_coefs, Point4* coefs, Point2 centroid);
void format_coefs(Point4* coefs, Point2 centroid, int num_coefs, CCsegtool_parameters *parameters);
void fourier_reconst(Point4* coefs, int num_coefs, Point2 centroid, int num_samples, Point2* sample_pts);
void fourier_interp(Point2* sample_pts, int num_samples, int num_coefs, Point4* coefs, Point2 centroid);
void compute_normals(Point4* coefs, int num_coefs, int num_pts, Point2* normals);
void normalize_coefs(Point4* coefs, int num_coefs, double *rot, double *scale);
void unnormalize_coefs(CCsegtool_parameters *parameters, Point4* _coefs,unsigned int num_coefs, double rotation, double scale);
void rot_coefs(Point4* coefs, int num_coefs, double rot);
void scale_coefs(Point4* coefs, int num_coefs, double scale);

#endif /* _FOURIER_FUNCS_H_ */
