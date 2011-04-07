#ifndef _CCCURVEPROBA_H_
#define _CCCURVEPROBA_H_

#include "CCsegtool_parameters.h"
#include <vtkPolyLine.h>
#include <string.h>
#include <sys/types.h>
#include <assert.h>

#include <itkImageRegionIterator.h>
#include <itkImageRegionIteratorWithIndex.h>
#include <itkImageSliceIteratorWithIndex.h>
#include <itkPolygonSpatialObject.h>
#include <itkSpatialObjectToImageFilter.h>

typedef struct{ /* A structure for a 2D point*/
        double x,y; /* x and y coordinates of the point*/
}       Point2D;

typedef struct{ /* A structure for a 2D polygon*/
	int numPoints; /* Number of points in the polygon*/
	Point2D *pointArray; /* An array to store the 2D points*/
}       Polygon2D;


typedef itk::PolygonSpatialObject<ImageIDimension>                                 PolygonSpatialObjectType;
typedef PolygonSpatialObjectType::PointType                                        PointType;
typedef itk::SpatialObjectToImageFilter<PolygonSpatialObjectType, InputImageType>  DigitizePolygonFilterType;
typedef itk::ImageDuplicator<InputImageType>                                       DuplicatorType;
typedef itk::ImageRegionIteratorWithIndex<InputImageType>                          IteratorType;

class CCcurveProba
{
	public:
		CCcurveProba(InputImageType::Pointer Image3D, vtkSmartPointer<vtkPoints> points, 
			     std::string CCAtlasDirectory, std::string outputfile, std::string nameofproject, 
			     bool vesselRemoveOn, int segLabel, int averageNum, bool permute_x_y, bool reflectXOn, 
			     bool reflectYOn, bool openOn,bool doubleOn, int sliceDir, int sliceNum, int curvenumpts,
			     bool Debug=false);
		void compute_proba();
		
	protected:
		int loadProbaMaps(std::string probFile);
		double squaredDistance (const double * point1, const double * point2);
		
	private:
		bool debug;
		/* Index */
		InputImageType::Pointer    m_Image3D;
		vtkSmartPointer<vtkPoints> m_points;
		
		/* Changement on the curve */
		bool m_vesselRemoveOn;
		int  m_segLabel;
		int  m_averageNum;
		bool m_permute_x_y;
		bool m_reflectXOn;
		bool m_reflectYOn;
		bool m_openOn;
		bool m_doubleOn;
		int  m_sliceDir;
		int  m_sliceNum;
		
		/* string for outputfile */
		std::string m_outputfile;
		std::string m_nameofproject;
		/* string for Atlas directory */
		std::string m_CCAtlasDirectory;
		
		/* Variables */
		double *m_probmaps;
		int     m_curvenumpts;
		int     m_numPts;
		int     m_numParc;
};

#endif
