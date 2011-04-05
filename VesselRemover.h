/*=========================================================================

  Header file for Vessel removal M. Styner 

=========================================================================*/
#ifndef _VESSEL_REMOVAL_H_
#define _VESSEL_REMOVAL_H_

#include <itkImage.h>
#include <itkOrientedImage.h>
#include <itkMinimumMaximumImageCalculator.h>
#include <itkHistogram.h> 
#include <itkMaskImageFilter.h>
#include <itkImageDuplicator.h>
#include <itkBinaryThresholdImageFilter.h> 

/** /class VesselRemover
 *
 * 
 *  Removes high intensity Vessel based on EMS segmentation (for IR-prepped SPGR only)
 *  Everything with higher intensity than mean WM and no segmentation label is set to 
 *  the mean of CSF,GM,WM masked image
 */

class ITK_EXPORT  VesselRemover : public itk::Object
{
	public:
	
		/** Standard class typedefs. */
		typedef VesselRemover Self;
		typedef itk::Object  Superclass;
		typedef itk::SmartPointer<Self>   Pointer;
		typedef itk::SmartPointer<const Self>  ConstPointer;
		
		/** Method for creation through the object factory. */
		itkNewMacro(Self);
		
		/** Run-time type information (and related methods). */
		itkTypeMacro(Self, Object);
		
		/** Dimension of the images to be segmented*/ 
		enum { ImageDimension = 3 };
		
		/** Pixel type used for reading the input image */
		typedef   short                        PixelType;
		
		/** Type of the input image */
		typedef   itk::OrientedImage<PixelType,ImageDimension>    ImageType;
		typedef   ImageType::Pointer  ImagePointer;
		
		/** Const Pointer type for the image. */
		typedef   ImageType::ConstPointer ImageConstPointer;
		typedef   itk::Statistics::Histogram<> HistogramType;
		typedef   ImageType::RegionType ImageRegionType;
		typedef   itk::ImageRegionIterator< ImageType > Iterator;
		typedef   itk::ImageRegionConstIterator< ImageType >   ConstIterator;
		
		itkSetConstObjectMacro(Image,ImageType);
		itkSetConstObjectMacro(EMSseg,ImageType);
		itkGetObjectMacro(ResultImage,ImageType);
		
		itkSetMacro(WMlabel, int);
		itkGetMacro(WMlabel, int);
		
		static const int BG_VALUE = 0;
		static const int LABEL_VALUE = 1;
		
		void RemoveVessels(void);
	
	protected:
		VesselRemover();
		virtual ~VesselRemover() {};
		void PrintSelf(std::ostream& os, itk::Indent indent) const;
	
	private:
		
		ImageType::ConstPointer m_Image;
		ImageType::ConstPointer m_EMSseg;
		ImageType::Pointer m_ResultImage;
		
		int       m_WMlabel;
		
		typedef   itk::MinimumMaximumImageCalculator< ImageType > minMaxCalcType;
		typedef   itk::MaskImageFilter< ImageType, ImageType, ImageType> MaskImageType;
		typedef   itk::ImageDuplicator< ImageType >       ImageDuplicatorType;
		typedef   itk::BinaryThresholdImageFilter< ImageType, ImageType > threshFilterType;

};



#endif /* _VESSEL_REMOVAL_H_ */
