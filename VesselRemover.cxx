 // *************************************************************************
// *************************************************************************
//  Removes high intensity Vessel based on EMS segmentation (for IR-prepped SPGR only)
//  Everything with higher intensity than mean WM and no segmentation label is set to
//  the mean of CSF,GM,WM masked image
// *************************************************************************
// *************************************************************************
#include "VesselRemover.h"

#include <ctime>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <fstream>
#include <iostream>

/*
 * Constructor
 */

VesselRemover::VesselRemover()
{
	m_Image = ImageType::New();
	m_EMSseg = ImageType::New();
}

void VesselRemover::PrintSelf( std::ostream& os, itk::Indent indent ) const
{
	Superclass::PrintSelf(os,indent);
}

void VesselRemover::RemoveVessels(void)
{
	// extract the label from the segmentation 
	threshFilterType::Pointer threshFilter = threshFilterType::New();
	threshFilter->SetInput(m_EMSseg);
	threshFilter->SetUpperThreshold(m_WMlabel);
	threshFilter->SetLowerThreshold(m_WMlabel);
	threshFilter->SetOutsideValue(BG_VALUE);
	threshFilter->SetInsideValue(LABEL_VALUE);
	threshFilter->Update();
	
	// Mask the input image
	MaskImageType::Pointer maskFilter = MaskImageType::New();
	maskFilter->SetInput1(m_Image);
	maskFilter->SetInput2(threshFilter->GetOutput());
	maskFilter->Update();
	
	double WMavg = 0;
	double WMnum = 0;
	{
		ImageType::Pointer maskedImage = maskFilter->GetOutput();
		// average value of masked image
		ConstIterator iterImage (maskedImage, maskedImage->GetBufferedRegion());
		// compute average
		while ( !iterImage.IsAtEnd() )  {
			double value =  iterImage.Get();
			if (value > 0.0) {
				WMavg += value;
				WMnum++;
			}
			++iterImage;
		}
		WMavg /= WMnum;
	}
	
	// Mask the input image
	MaskImageType::Pointer maskFilter2 = MaskImageType::New();
	maskFilter2->SetInput1(m_Image);
	maskFilter2->SetInput2(m_EMSseg);
	maskFilter2->Update();
	
	double avg = 0;
	double num = 0;
	{
		ImageType::Pointer maskedImage = maskFilter2->GetOutput();
		// average value of full image
		ConstIterator iterImage (maskedImage, maskedImage->GetBufferedRegion());
		// compute average
		while ( !iterImage.IsAtEnd() )  {
			double value =  iterImage.Get();
			if (value > 0.0) {
				avg += value;
				num++;
			}
			++iterImage;
		}
		avg /= num;
	}
	
	
	
	ImageDuplicatorType::Pointer duplicator = ImageDuplicatorType::New();
	duplicator->SetInputImage(m_Image);
	duplicator->Update();
	m_ResultImage = duplicator->GetOutput();
	
	// all pixels > WMavg and not in EMSsegmentation -> = avg
	{
		ConstIterator iterImage (m_Image, m_Image->GetBufferedRegion());
		ConstIterator iterMask  (m_EMSseg, m_EMSseg->GetBufferedRegion());
		Iterator iterRes (m_ResultImage, m_ResultImage->GetBufferedRegion());
		while ( !iterImage.IsAtEnd() )  {
			double value =  iterImage.Get();
			if (value > WMavg && iterMask.Get() == 0)
				iterRes.Set((short) avg);
			++iterRes;
			++iterImage;
			++iterMask;
		}
	}
}
