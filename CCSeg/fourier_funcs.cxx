#include "fourier_funcs.h"

static double sin_t[MAX_STEPS];
static double cos_t[MAX_STEPS];

/* **************************************************************************/
/* sincos_table fills in the global sin_t and cos_t tables                  */
/* **************************************************************************/
void sincos_table(int n_steps)
{
	int i;
	for (i= 0; i < n_steps; i++) {
		sin_t[i]= sin(2.0*M_PI*i/n_steps);
		cos_t[i]= cos(2.0*M_PI*i/n_steps);
	}
}

/* ************************************************************************** 
 * parse_coefs splits off the centroid info from the rest of the coefs and
 * translates from a0c0 (real, sin/cos) format to w00z (complex) format.
 * Assumes input field is in row-major order (e.g. AVS field 1D 4-vector)
 * **************************************************************************/
void parse_coefs( CCsegtool_parameters *parameters, int num_coefs, Point4* coefs, Point2 centroid)
{
	int i,j;
	/* First line holds the centroid information */
	centroid[0] = parameters->GetSMMeanval(0,0);
	centroid[1] = parameters->GetSMMeanval(1,0);
	for (i = 1; i <= num_coefs; i++)
	{
		for (j = 0; j < 4; j++)
		{
			coefs[i-1][j] = (double) parameters->GetSMMeanval(j,i);
		}
	}
	for (i = num_coefs+1; i <= MAX_COEFS; i++)
	{
		for (j = 0; j < 4; j++)
		{
			coefs[i-1][j] = 0.0;
		}
	}
}

/* **************************************************************************/
/* format_coefs packages a Fourier contour from our internal data structure */
/* to a standard AVS field 1D 4-vector uniform float (row-major).           */
/* **************************************************************************/
void format_coefs(Point4* coefs,	/* input Fourier contour */
		  Point2 centroid, 	/* input centroid of Fourier contour */
		  int num_coefs,	/* input size of coefs[] */
		  CCsegtool_parameters *parameters)	/* output data field of AVS field */
{
	int i, j;
	std::vector<float> out_fld;
	out_fld.push_back(centroid[0]);
	out_fld.push_back(centroid[1]);
	out_fld.push_back(0);
	out_fld.push_back(0);
	for (i=1; i<=num_coefs; i++)
		for (j=0; j<4; j++)
			out_fld.push_back(coefs[i-1][j]);
	parameters->SetSSCoefs(out_fld);
}

/* constants used to select op in fourier_backend */
#define RECONST 0
#define TANGENT 1

/* **************************************************************************/
/* Double sum over points around the contour, and over the coefficients.    */
/* Both contour reconstruction (coeffs -> pts) and calculation of tangent   */
/* vectors (coeffs -> vectors) are of this form.  Internal function only.   */
/* **************************************************************************/
static void fourier_backend(Point4* coefs,		/* input contour */
			    int num_coefs,		/* input size of coefs[] */
			    Point2 centroid, 		/* input centroid of contour */
			    int num_samples, 		/* input number of points to extract */
			    int op,			/* input type of Fourier operation */
			    Point2* pts)		/* output points in 2D image space */
{
	int h, xy, phi;
	int idx;
	double c, s;
	for (phi = 0; phi < num_samples; phi++) {
		for (xy = 0; xy < 2; xy++) {
			pts[phi][xy] = centroid[xy];
			for (h = 1; h <= num_coefs; h++) {
				idx = (h*phi)%num_samples;
				c = coefs[h-1][2*xy  ];
				s = coefs[h-1][2*xy+1];
				pts[phi][xy] += (op == RECONST) ?   cos_t[idx]*c +   sin_t[idx]*s
					: h*cos_t[idx]*s - h*sin_t[idx]*c;
			}
		}
	}
}

/* **************************************************************************/
/* fourier_reconst calculates the (x,y) positions of the specified number   */
/* of sample points evenly spaced around the given Fourier contour.         */
/* This is just a proxy to the internal function fourier_backend.           */
/* **************************************************************************/
void fourier_reconst(Point4* coefs,	/* input contour */
		     int num_coefs,	/* input size of coefs[] */
		     Point2 centroid, 	/* input centroid of contour */
		     int num_samples, 	/* input number of points to extract */
		     Point2* pts)	/* output polygonal approximation */
{
	fourier_backend(coefs, num_coefs, centroid, num_samples, RECONST, pts);
}

/* **************************************************************************/
/* compute_normals calculates the unit normal vectors to the specified      */
/* Fourier contour at a number of points evenly spaced around the contour.  */
/* This is just a proxy to the internal function fourier_backend.           */
/* **************************************************************************/
void compute_normals(Point4* coefs,	/* input contour */
		     int num_coefs,	/* input size of coefs[] */
		     int num_pts,	/* input number of sample points */
		     Point2* normals)	/* output unit normal vectors */
{
	double len, nx, ny;
	Point2 centroid;
	centroid[0] = centroid[1] = 0.0;	/* don't use centroid here */
	
	fourier_backend(coefs, num_coefs, centroid, num_pts, TANGENT, normals);
	
	for (int i=0; i<num_pts; i++) {		/* normalize vectors */
		len = sqrt(normals[i][0]*normals[i][0] + normals[i][1]*normals[i][1]);
		nx = -normals[i][1]/len;		/* and take perpendicular */
		ny =  normals[i][0]/len;
		normals[i][0] = nx;
		normals[i][1] = ny;
	}
}

/* ************************************************************************** 
 * fourier_interp fits a Fourier contour through the list of 2D points.
 * The coefficients are not explicitly normalized for parameterization 
 * starting point.
 * **************************************************************************/
void fourier_interp(Point2* sample_pts,	/* input PDM of contour points */
		    int num_samples,	/* input size of sample_pts[] */
		    int num_coefs,		/* input size of coefs[] */
		    Point4* coefs, 		/* output Fourier coefficients */
		    Point2 centroid)	/* output centroid of input points */
{
	int i, idx, order;
	Point4 t_coef;
	double x, y;
	
	/* calculate centroid */
	centroid[0] = centroid[1] = 0.0;
	for (i=0; i<num_samples; i++) {
		centroid[0] += sample_pts[i][0];
		centroid[1] += sample_pts[i][1];
	}
	centroid[0] /= num_samples;
	centroid[1] /= num_samples;
	
	for (order=1; order<=num_coefs; order++) {	/* find one coef */
		t_coef[0] = t_coef[1] = t_coef[2] = t_coef[3] = 0.0;
		for (i=0; i<num_samples; i++) {		/* loop around contour */
			idx = (order*i)%num_samples;
			x = sample_pts[i][0] - centroid[0];
			y = sample_pts[i][1] - centroid[1];
			t_coef[0] += cos_t[idx]*x;
			t_coef[1] += sin_t[idx]*x;
			t_coef[2] += cos_t[idx]*y;
			t_coef[3] += sin_t[idx]*y;
		}
		for (i=0; i<4; i++)
			coefs[order-1][i] = t_coef[i]*2/num_samples;
	}
}

/*****************************************************************************
 * normalize rotation and scale by first ellipse (modify in-place)
 *****************************************************************************/
void normalize_coefs(Point4* coefs,	/* input/output contour description */
		     int num_coefs,	/* input size of coefs[] */
		     double *rot,	/* output original rotation */
		     double *scale)	/* output original scaling */
{
	/* need to norm starting point of parameterization first ? */
	*rot   = -atan2( coefs[0][2], coefs[0][0] );
	
	#define hypot(x,y) sqrt((x)*(x)+(y)*(y))
	*scale = hypot(coefs[0][2], coefs[0][0]);
	
	rot_coefs  (coefs, num_coefs, -(*rot));
	scale_coefs(coefs, num_coefs, 1.0/(*scale));
}

/***************************************************************************** 
 * unnormalize coefs
 *****************************************************************************/
void unnormalize_coefs(CCsegtool_parameters *parameters, Point4* _coefs,
		       unsigned int num_coefs,double rotation, double scale)
{	
	rot_coefs(_coefs, num_coefs, rotation);
	scale_coefs(_coefs, num_coefs, scale);
}

/***************************************************************************** 
 * Modify coefs in-place to achieve a rotation of the contour
 *****************************************************************************/
void rot_coefs(Point4* coefs,	/* input/output contour description */
		int num_coefs,	/* input size of coefs[] */
		double rot)	/* input rotation to apply */
{
	int i;
	Point4 tmp;
	for (i=0; i<num_coefs; i++) {
		tmp[0] = cos(rot)*coefs[i][0] + sin(rot)*coefs[i][2];
		tmp[1] = cos(rot)*coefs[i][1] + sin(rot)*coefs[i][3];
		tmp[2] = cos(rot)*coefs[i][2] - sin(rot)*coefs[i][0];
		tmp[3] = cos(rot)*coefs[i][3] - sin(rot)*coefs[i][1];
		coefs[i][0] = tmp[0];
		coefs[i][1] = tmp[1];
		coefs[i][2] = tmp[2];
		coefs[i][3] = tmp[3];
	}
}

/***************************************************************************** 
 * Modify coefs in-place to achieve a scaling  of the contour
 *****************************************************************************/
void scale_coefs(Point4* coefs,	/* input/output contour description */
		int num_coefs,	/* input size of coefs[] */
		double scale)	/* input scaling to apply */
{
	int i, j;
	for (i=0; i<num_coefs; i++) 
		for (j=0; j<4; j++)
			coefs[i][j] *= scale;
}

