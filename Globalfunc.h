#ifndef _GLOBALFUNC_H_
#define _GLOBALFUNC_H_

#include "CCcurveProba.h"
#include "Imageview.h"
#include "CCsegtool_computation.h"

/* Functions defined in, exported from Globalfunc.cxx */
int initialize(std::string Image_filename, std::string Seg_filename, std::string CCAtlasDirectory,
		           std::string Path_output, std::string nameofproject, bool interpolationlinear,
			   bool vesselRemoveOn, int segLabel, int averageNum, bool permute_x_y, bool reflectXOn, 
			   bool reflectYOn, bool openOn, bool doubleOn, int sliceDir, std::string MidPlaneSliceNumber,
			   int number_pts, bool FSXForm, double PSDistance, bool Unconstrained, 
			   std::string WMintensity, std::string MPSDisplacement, std::string Number_iteration, 
			   int Lambdamax, std::string Coefofoptim, bool debug, bool withgui, Imageview * &imageview,
			   CCsegtool_parameters* &last_parameters, double *lastX, double *lastY, double *lastScale,
			   double *lastRot, float scalefactor, bool othercompo, int angle, 
			   bool rot90=false, bool rot180=false, bool rot270=false,
			   QWidget * parent=NULL, QLabel* parametersinitview=NULL, QLabel* Arealabel=NULL, 
			   QLabel* Midsaggitalplaneslicevalue=NULL);
void compute( std::string Path_output, std::string nameofproject, int number_pts, bool Unconstrained, 
	      std::string WMintensity, std::string MPSDisplacement, std::string Number_iteration, int Lambdamax,
	      std::string Coefofoptim, bool debug, bool withgui, Imageview * &imageview, 
	      CCsegtool_parameters* &last_parameters, double *lastX, double *lastY, double *lastScale, double *lastRot,
	      std::string A="default", double r=0.0, bool rot90=false, bool rot180=false, bool rot270=false,
	      QWidget * parent=NULL, QLabel* instructionsview=NULL, QLabel* parametersinitview=NULL,QLabel* Arealabel=NULL);
int readstringWMI(std::string str, int position, int initWMvalue);
float readstringinputf(std::string str, int position);
int readstringinputi(std::string str, int position, int type);
QString readparam(std::string filename);
void showinstruction(QLabel* view, int step, int iter, bool Unconstrained, int Number_iteration);
int run_step(CCsegtool_parameters* parameters, CCsegtool_computation* computation, Imageview* imageview, QWidget* parent,
	     int step, int WMvalue, int MPSDisplacement, int Number_iteration, float coef, Point2 *pts, int number_pts,
	     int Lambdamax, float disp, bool debug, bool withgui, bool preview=false);
void updatedrawing(double X,double Y,double Scale,double Rotation, QWidget* parent,CCsegtool_parameters *parameters, 			   Imageview* &imageview, int number_pts, std::string Path_output,std::string nameofproject, int Lambdamax,
		   bool debug);
double CalculBestfitmean(std::string Path_output,  std::string nameofproject);
double contourArea(Polygon2D poly);
void saveArea(double area, std::string Path_output, std::string nameofproject);

#endif
