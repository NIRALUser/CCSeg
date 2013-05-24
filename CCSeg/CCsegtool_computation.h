#ifndef _CCSEGTOOL_COMPUTATION_H_
#define _CCSEGTOOL_COMPUTATION_H_

#include "fourier_funcs.h"
#include <sstream>
#include <itkIndex.h>

class CCsegtool_computation
{
	public :
		/* constructor */
		CCsegtool_computation(CCsegtool_parameters * parameters, std::string nameofproject, std::string path_output,
				      bool Debug=false,double A=0.0, double r=0.0);
		/* set the repulsive points */
		void setRepulPoints(std::vector<double> *RepulPoints);
		
		/* methods */
		/* Main method : execution */
		float execution(CCsegtool_parameters * parameters, int step, int lambdamax,int iteration, bool update=false);
		
		/* get the profiles */
		void get_profiles();
		float get_shade(Point2 *pt);
		ExtractImageType::IndexType getindex(int x, int y);
		
		/* calcul of the shift */
		int profile_shift( float* imgprof, double* GOFvalue, int proflen, float grey_scaling,int pointID);
		double goodness_of_fit( float* profile, int proflen,int pointID);
		std::vector<float>::iterator getiterator(std::vector<float> vect,int i);
		double computeRepulsionExponentialPenalty(int pointID, double shift);
		double exponential(double X, double Y, int pointRepulsivID);
		void calcshift();
		
		/* functions on the coefficients */
		void subtractmean( Point4* _coefs);
		void addmean( Point4* _coefs);
		void proj_eigvec( Point4* _coefs);
		void restrict_eigloads(int lambdamax);
		void dotpord_eigvec( Point4* _coefs);
		
		/* functions to save the ouput and to reconstruct the curve */
		void reconstructCurve(Point2* pts, int number_pts);
		void writefileasc( int step, int iteration);
		void writeoutput(double Xsize, double Ysize);
		void writedp(float coefoptim1, float coefoptim2, bool unconstrained);
		
		/* calcul the optimization of the loop */
		float calc_dp(Point2* sample_pts, Point2* old_sample_pts);
		
		/* get or set the center, scale and rotation */
		void getparam(double *lastX, double *lastY, double *lastScale, double *lastRot);
		void setparam(double X, double Y, double Scale, double Rot);
		
		//accesd to private attributs
		void Setnum_coefs(unsigned int _num_coefs);
		unsigned int Getnum_coefs();
		void Setnum_samples(unsigned int _num_samples);
		unsigned int Getnum_samples();
		void Setnum_modes(unsigned int _num_modes);
		unsigned int Getnum_modes();
		unsigned int Getproflen();
		unsigned int Getext_proflen();
		double GetcenterX();
		double GetcenterY();
		void Setcoefs(Point4* _coefs);
		Point4* Getcoefs();
		void Setsample_pts(Point2* _sample_pts);
		Point2* Getsample_pts();
		double Getsample_ptsval(int i, int j);
		void Setnormals(Point2* _normals);
		Point2* Getnormals();
		
		//parameter
		Point2 m_centroid;
		
	private :
		/* Debug */
		bool debug;
		std::string m_path_output;
		std::string m_nameofproject;
		
		/*Pointer on the parameters*/
		CCsegtool_parameters *m_parameters;

		/* misc array dimensions */
		unsigned int m_num_coefs, m_num_samples, m_num_modes, m_proflen, m_ext_proflen;
		double       m_rot, m_scale;
		double       m_A, m_r, m_bestfitmean;
		bool         m_firstTime;
		
		/* local data structures */
		std::vector<float>  m_dp;
		std::vector<double> m_RepulPoints;
		Point4* m_coefs;
		Point2* m_sample_pts;
		Point2* m_old_sample_pts;
		Point2* m_normals;
};

#endif

