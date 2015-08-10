#ifndef _CCSEGTOOL_INITIALIZATION_H_
#define _CCSEGTOOL_INITIALIZATION_H_


#include "VesselRemover.h"
#include <iostream>
#include <time.h>
#include <math.h>
#include <fstream>
#include <ostream>

//ITK libraries for Image
#include <itkObject.h>
#include <itkImage.h>
#include <itkObjectFactory.h>
#include <itkImageFileReader.h> 
#include <itkImageFileWriter.h>
#include <itkFlipImageFilter.h>
#include <itkImageMomentsCalculator.h>
#include <itkCenteredAffineTransform.h>
#include <itkAffineTransform.h>
#include <itkResampleImageFilter.h>
#include <itkBSplineInterpolateImageFunction.h>
#include <itkNearestNeighborInterpolateImageFunction.h>
#include <itkAddImageFilter.h>
#include <itkShiftScaleImageFilter.h>
#include <itkRelabelComponentImageFilter.h>
#include <itkConnectedComponentImageFilter.h>
#include <itkBinaryBallStructuringElement.h> 
#include <itkBinaryDilateImageFilter.h> 
#include <itkBinaryErodeImageFilter.h> 
#include <itkExtractImageFilter.h> 
#include <itkPermuteAxesImageFilter.h> 
#include <itkCastImageFilter.h> 
#include <itkMinimumMaximumImageFilter.h>
#include "itkIntensityWindowingImageFilter.h"
#include "itkScalarImageToHistogramGenerator.h"

typedef double Point2[2];
typedef double Point3[3];
typedef double Point4[4];

// itk definitions
/** Pixel type used for reading the input image */
typedef   short                         InputPixelType;

/** Pixel type to be used internally */
typedef   unsigned char                 BinaryPixelType;
typedef   unsigned char                 OutputPixelType; // the segmentation as an image

/** Dimension */
/** I for input, in 3D and E for extract, in 2D **/
enum { ImageIDimension = 3 };
enum { ImageEDimension = 2 };

/** Dimension 3 **/
/** Type of the input image */
typedef itk::Image<InputPixelType,ImageIDimension>     InputImageType;
typedef itk::Image<BinaryPixelType,ImageIDimension>    BinaryImageType;
typedef InputImageType::Pointer                                ImagePointer;
typedef InputImageType::RegionType                             ImageRegionType;
typedef InputImageType::SpacingType                            SpacingType;
typedef InputImageType::SizeType                               SizeType;
typedef itk::ContinuousIndex<float,ImageIDimension>            ContinuousIndexType;
typedef itk::Point<float,ImageIDimension>                      ImagePointType;
typedef itk::Index<ImageIDimension>                            InputImageIndex;

/** Dimension 2 **/
/** Type of the internal image */
typedef itk::Image<InputPixelType,ImageEDimension>          ExtractImageType;
/** Const Pointer type for the image. */
typedef ExtractImageType::ConstPointer                              ExtractImageConstPointer;
/** Type of the internal image */
typedef itk::Image<BinaryPixelType,ImageEDimension>         ExtractBinaryImageType;
typedef ExtractBinaryImageType::ConstPointer                        ExtractBinaryImageConstPointer;
/** Type of the output image */
typedef itk::Image<OutputPixelType,ImageEDimension>         OutputImageType;
/** Use in  compute_parameters **/
typedef itk::Statistics::Histogram<double>                  HistogramType;
typedef ExtractImageType::RegionType                        ExtractImageRegionType;
typedef itk::ImageRegionIterator< ExtractImageType >        Iterator;

typedef itk::ImageFileReader< InputImageType >    VolumeReaderType;
typedef itk::ImageFileReader< BinaryImageType >   BinaryVolumeReaderType;
typedef itk::ImageFileWriter< InputImageType >    ImageWriterType;
typedef itk::ImageFileWriter< ExtractImageType >  ImageWriterType2;
typedef itk::ImageFileWriter< OutputImageType >   ImageWriterType3;
typedef itk::ImageFileWriter< BinaryImageType >   BinaryWriterType;

typedef itk::CastImageFilter<ExtractImageType,InputImageType>              Cast2d3dFilterType;
typedef itk::CastImageFilter<ExtractBinaryImageType,BinaryImageType>       BinaryCast2d3dFilterType;
typedef itk::CastImageFilter<ExtractImageType, OutputImageType>            CastPNGFilterType;

typedef itk::CastImageFilter<BinaryImageType,InputImageType>               CastImageFilterType;
typedef itk::ExtractImageFilter< InputImageType, ExtractImageType >        ExtractFilterType ;
typedef itk::ExtractImageFilter< BinaryImageType, ExtractBinaryImageType > BinaryExtractFilterType ;
typedef itk::FlipImageFilter< ExtractImageType>                            ReflectFilterType ;
typedef itk::FlipImageFilter< ExtractBinaryImageType>                      BinaryReflectFilterType ;
typedef itk::PermuteAxesImageFilter< ExtractImageType>                     PermuteAxesImageFilterType;
typedef itk::PermuteAxesImageFilter< ExtractBinaryImageType>               BinaryPermuteAxesImageFilterType ;

typedef itk::BinaryThresholdImageFilter<ExtractBinaryImageType,ExtractBinaryImageType> BinaryThresholdFilterType;
typedef itk::BinaryBallStructuringElement<BinaryPixelType,2>                           StructuringElementType;
typedef itk::BinaryDilateImageFilter<ExtractBinaryImageType,ExtractBinaryImageType,StructuringElementType> dilateFilterType;
typedef itk::BinaryErodeImageFilter<ExtractBinaryImageType,ExtractBinaryImageType,StructuringElementType>  erodeFilterType;

typedef itk::AffineTransform<double, 2>                                             TransformType;
typedef itk::BSplineInterpolateImageFunction<ExtractImageType,double,double>        SplineType;
typedef itk::NearestNeighborInterpolateImageFunction<ExtractBinaryImageType,double> NNType;
typedef itk::ResampleImageFilter<ExtractImageType,ExtractImageType>                 ResampleFilterType;
typedef itk::ResampleImageFilter<ExtractBinaryImageType,ExtractBinaryImageType>     BinaryResampleFilterType;
typedef itk::AddImageFilter<ExtractImageType,ExtractImageType,ExtractImageType>     AddFilterType;
typedef itk::ShiftScaleImageFilter<ExtractImageType,ExtractImageType>               ScaleFilterType;
typedef itk::ImageMomentsCalculator< BinaryImageType >                              MomentsCalcType;

typedef itk::MinimumMaximumImageFilter< InputImageType >			        MinimumMaximumType;
typedef itk::Statistics::ScalarImageToHistogramGenerator< InputImageType > 		HistogramFilterType;
typedef itk::IntensityWindowingImageFilter< InputImageType, InputImageType > 		WindowingFilterType;

class CCsegtool_initialization
{
	public:
		/* constructor */
		CCsegtool_initialization(bool interpolationlinear, float scalefactor=1.35);
		
		static const int BG_VALUE = 0;
		static const int LABEL_VALUE = 1;
		static const int NUM_COMPONENTS = 5;
		
		int compute_initialization(std::string inputFileName, std::string segFile, bool vesselRemoveOn, 
					    bool segLabel, int averageNum, bool permute_x_y, bool reflectXOn, 
					    bool reflectYOn, bool openOn, bool doubleOn, int sliceDir, 
					    std::string outfilebase, std::string nameofproject,
					    std::string MidPlaneSliceNumber, bool othercompo, int angle, bool debug=false);
		void SetLambdaMax(int lambdaMax){m_LambdaMax=lambdaMax;}
		int compute_parameters(void);
		
		//Get and Set the parameters
		ExtractImageType::Pointer  GetImage();
		int GetCenterX();
		int GetCenterY();
		double GetRotation();
		double GetScale();
		InputPixelType GetWMvalue();
		ExtractImageConstPointer GetOutputImage();
		float GetvoxelsizeX();
		float GetvoxelsizeY();
		int   GetImageSSize();
		SizeType Get3DImageSize();
		SpacingType Get3DImageSpacing();
		ImagePointer Get3DImage();
		
	protected:
		void loadinginputimage(std::string inputFileName, std::string segFile);
		void vesselremoval(bool segLabel);
		void extract_Midsagtital_planes(int sliceDir, std::string MidPlaneSliceNumber);
		void averaging(int averageNum, int sliceDir, std::string MidPlaneSliceNumber);
		void dopermute_x_y();
		void reflectX();
		void reflectY();
		void extractLabel(bool segLabel);
		void closingop();
		void dodoubleOn();
		void writeInput(std::string outfilebase, std::string nameofproject);
		void writeoutput(std::string outfilebase, std::string nameofproject);
		

	private:
		/** Use in compute_parameters **/
		typedef   itk::Image<short,ImageEDimension>     CCLOutputImageType;
		typedef   itk::ConnectedComponentImageFilter< ExtractBinaryImageType, CCLOutputImageType> CCLFilterType;
		typedef   itk::RelabelComponentImageFilter<CCLOutputImageType, CCLOutputImageType> RelabelFilterType;
		typedef   itk::BinaryThresholdImageFilter< CCLOutputImageType,  ExtractBinaryImageType> threshCCLFilterType;
		typedef   itk::ImageMomentsCalculator<  ExtractBinaryImageType > MomentsCalcType;
		typedef   itk::MaskImageFilter< ExtractImageType, ExtractBinaryImageType, ExtractImageType> MaskImageType;
		typedef   itk::MinimumMaximumImageCalculator< ExtractImageType > minMaxCalcType;
		
		// data structures
		ImagePointer                       m_loadImage;
		ImageRegionType                    m_extractRegion;
		BinaryImageType::Pointer           m_loadMask;
		ImagePointer                       m_preProcImage ;
		ExtractImageType::Pointer          m_inputImage ;
		ExtractBinaryImageType::Pointer    m_mask ;
		
		ExtractImageConstPointer           m_extractImage;
		ExtractBinaryImageConstPointer     m_extractMask;
		ExtractImageRegionType             m_Region;
		BinaryThresholdFilterType::Pointer m_threshFilter;
		dilateFilterType::Pointer          m_dilateFilter;
		erodeFilterType::Pointer           m_erodeFilter;
		CastImageFilterType::Pointer       m_castfilter;
		ExtractFilterType::Pointer         m_extractFilter;
		Cast2d3dFilterType::Pointer        m_cast3DFilter;
		BinaryCast2d3dFilterType::Pointer  m_bincast3DFilter;
		threshCCLFilterType::Pointer       m_thresh1Filter;
		threshCCLFilterType::Pointer       m_thresh2Filter;
		threshCCLFilterType::Pointer       m_thresh3Filter;
		MaskImageType::Pointer             m_maskFilter;
		CastPNGFilterType::Pointer         m_uncharoutput;
		
		bool                 m_OnlyConstrainedOn;
		bool                 m_interpolationlinear;
		float                m_scalefactor;
		
		double               m_LambdaMax;
		InputPixelType       m_WMvalue;
		double               m_Rotation;
		double               m_Scale;
		int                  m_CenterX;
		int                  m_CenterY;
		int                  m_numSlices;
		
		float                m_voxelsizeX;
		float                m_voxelsizeY;
		
		SizeType             m_imagesize;
		SpacingType          m_imagespacing;
		
		int                  m_Xsize;
		int                  m_Ysize;
		
		bool                 m_permute_x_y;
		bool                 m_reflectXOn;
		bool                 m_reflectYOn;
		
		bool                 m_othercompo;
		int                  m_angle;
};

#endif
