// ************************************************************************
// *************************************************************************
// Segment the Corpus callosum from a 2D image
// *************************************************************************
// *************************************************************************

#include "CCsegtool_initialization.h"
#include <iostream>
#include <fstream>

/***************************************************************************
 * Constructor
 ***************************************************************************/

CCsegtool_initialization::CCsegtool_initialization(bool interpolationlinear, float scalefactor)
{
	m_loadImage = InputImageType::New();
	m_loadMask = BinaryImageType::New();
	m_inputImage = ExtractImageType::New();
	m_mask = ExtractBinaryImageType::New();
	m_extractImage = ExtractImageType::New();
	m_extractMask = ExtractBinaryImageType::New();
	m_WMvalue = 100;
	m_Rotation = 1.9;
	m_Scale = 27;
	m_CenterX = 116;
	m_CenterY = 105;
	m_LambdaMax = 3;
	m_OnlyConstrainedOn = false;
	m_interpolationlinear = interpolationlinear;
	m_scalefactor = scalefactor;
}


/***************************************************************************
 * Compute the initialization, transform the input image (3D to 2D) and
 * several other operation ( opening, rotation ....)
 ***************************************************************************/
int CCsegtool_initialization::compute_initialization(std::string inputFileName,std::string segFile, bool vesselRemoveOn,
		 bool segLabel, int averageNum, bool permute_x_y, bool reflectXOn, bool reflectYOn, bool openOn, 
		 bool doubleOn, int sliceDir, std::string outfileBase, std::string nameofproject, 
		 std::string MidPlaneSliceNumber, bool othercompo, int angle, bool debug)
{
	//* flip/permutation//othercomponant and supplement angle */
	m_permute_x_y = permute_x_y;
	m_reflectXOn = reflectXOn;
	m_reflectYOn = reflectYOn;
	m_othercompo = othercompo;
	m_angle = angle;
	
	// load
	loadinginputimage(inputFileName, segFile);
	
	//Vessel removal
	if(vesselRemoveOn)
		vesselremoval(segLabel);
	else
		m_preProcImage = m_loadImage;

	// Extract Midsagtital planes
	extract_Midsagtital_planes(sliceDir, MidPlaneSliceNumber);
	
	//Averaging?
	if (averageNum > 0)
		averaging(averageNum, sliceDir, MidPlaneSliceNumber);
	else
	{
		try 
		{
			ExtractFilterType::Pointer extractFilter = ExtractFilterType::New();
			extractFilter->SetInput(m_preProcImage);
			extractFilter->SetExtractionRegion(m_extractRegion);
#if  ITK_VERSION_MAJOR >=4
			extractFilter->SetDirectionCollapseToIdentity();
#endif
			extractFilter->Update();
			m_inputImage = extractFilter->GetOutput();
		}
		catch( itk::ExceptionObject & e )
		{
			std::cerr << "Error: extract filter fail" << std::endl;
			exit(1);
		}
	}
	
	//double the size of the image
	if (doubleOn)
		dodoubleOn();
	
	//Permutation of x and y axes
	if (m_permute_x_y)
		dopermute_x_y();

	//Reflection on X
	if (m_reflectXOn)
		reflectX();

	//Reflection on Y
	if (m_reflectYOn)
		reflectY();
	
	// extract the label from the segmentation
	extractLabel(segLabel);
	
	// closing operator
	closingop();
	if (openOn) { //opening
		try 
		{
			m_erodeFilter->SetInput(m_threshFilter->GetOutput());
			m_dilateFilter->SetInput(m_erodeFilter->GetOutput());
			m_dilateFilter->Update();
			m_mask = m_dilateFilter->GetOutput();
		}
		catch( itk::ExceptionObject & e )
		{
			std::cerr << "Error: opening fail" << std::endl;
			exit(1);
		}
	}
	else { //closing
		try 
		{
			m_dilateFilter->SetInput(m_threshFilter->GetOutput());
			m_erodeFilter->SetInput(m_dilateFilter->GetOutput());
			m_erodeFilter->Update();
			m_mask = m_erodeFilter->GetOutput();
		}
		catch( itk::ExceptionObject & e )
		{
			std::cerr << "Error: clothing fail" << std::endl;
			exit(1);
		}
	}
	
	
	
	// write input to segmentation
	writeInput(outfileBase,nameofproject);
	
	//Preparing for the segmentation
	m_extractMask=m_mask;
	m_extractImage=m_inputImage;
	
	//compute the parameters for X/Y/Scale/Rotation and WMvalue 
	if(compute_parameters()==1)
		return 1;
	
	// Voxelsize value and size of image
	itk::Vector<double,2> buf;
	buf = m_inputImage->GetSpacing();
	m_voxelsizeX = buf[0];  //x
	m_voxelsizeY = buf[1];  //y
	m_Xsize = m_Region.GetSize(0);
	m_Ysize = m_Region.GetSize(1);
	
	//Read the resulting and write output files
	writeoutput(outfileBase,nameofproject);
	
	if(debug)
		std::cout<<"  --compute_initialization done--"<<std::endl;
	
	return 0;
}

/***************************************************************************
 * Load the input and Mask image
 ***************************************************************************/
void CCsegtool_initialization::loadinginputimage(std::string inputFileName, std::string segFile)
{
	try 
	{
		// load image
		VolumeReaderType::Pointer imageReader = VolumeReaderType::New();
		imageReader->SetFileName(inputFileName.c_str()) ;
		imageReader->Update() ;
		m_loadImage=imageReader->GetOutput();
		// Take the Size and spacing of the 3D input image for Proba
		m_imagesize    = m_loadImage->GetLargestPossibleRegion().GetSize();
		m_imagespacing = m_loadImage->GetSpacing();
		MinimumMaximumType::Pointer minmax = MinimumMaximumType::New();
		minmax->SetInput(imageReader->GetOutput());
		minmax->Update();

		if (minmax->GetMinimum() < 0 || minmax->GetMaximum() > 255) {
			HistogramFilterType::Pointer histogramFilter = HistogramFilterType::New();

			histogramFilter->SetNumberOfBins( 255 );
			histogramFilter->SetMarginalScale( 1 );
			histogramFilter->SetHistogramMin( 0 );
			histogramFilter->SetHistogramMax( 32000 );
			histogramFilter->SetInput( imageReader->GetOutput() );
			histogramFilter->Compute();
			const HistogramFilterType::HistogramType   *histo=histogramFilter->GetOutput();

			double lower;
			double upper;
			lower = histo->Quantile(0,0.02);
			upper = histo->Quantile(0,0.995);

			WindowingFilterType::Pointer windowingFilter = WindowingFilterType::New();
			windowingFilter->SetInput( imageReader->GetOutput() );
			windowingFilter->SetOutputMinimum( 0 );
			windowingFilter->SetOutputMaximum( 255 );
			windowingFilter->SetWindowMinimum( lower );
			windowingFilter->SetWindowMaximum( upper );
			windowingFilter->Update();
			m_loadImage = windowingFilter->GetOutput();
		}
	}	
	catch( itk::ExceptionObject & e )
	{
		std::cerr << "Error: loading image fail : " << e << std::endl;
		exit(1);
	}

	// threshold stuff below 0
	// scale 0..max to 0..255 to work

	try 
	{
		// load segmentation
		BinaryVolumeReaderType::Pointer binaryImageReader = BinaryVolumeReaderType::New();
		binaryImageReader->SetFileName(segFile.c_str()) ;
		binaryImageReader->Update();
		m_loadMask=binaryImageReader->GetOutput();
	}
	catch( itk::ExceptionObject & e )
	{
		std::cerr << "Error: loading mask fail : " << e << std::endl;
		exit(1);
	}
}

/***************************************************************************
 * Do a vesselremoval
 ***************************************************************************/
void CCsegtool_initialization::vesselremoval(bool segLabel)
{
	try 
	{
		//castfilter
		m_castfilter = CastImageFilterType::New();
		m_castfilter->SetInput(m_loadMask);
		m_castfilter->Update();
	}
	catch( itk::ExceptionObject & e )
	{
		std::cerr << "Error: castfilter fail : " << e << std::endl;
		exit(1);
	}
	
	try 
	{
		//Vessel removal
		int seglabel=1;
		if(segLabel)
			seglabel=0;
		VesselRemover::Pointer VesselFilter = VesselRemover::New();
		VesselFilter->SetImage(m_loadImage);
		VesselFilter->SetEMSseg(m_castfilter->GetOutput());
		VesselFilter->SetWMlabel(seglabel);
		VesselFilter->RemoveVessels();
		m_preProcImage = VesselFilter->GetResultImage();
	}
	catch( itk::ExceptionObject & e )
	{
		std::cerr << "Error: vessel removal fail : " << e << std::endl;
		exit(1);
	}
}


/***************************************************************************
 * Extract the Midsagtital planes and set it for the mask
 ***************************************************************************/
void CCsegtool_initialization::extract_Midsagtital_planes(int sliceDir, std::string MidPlaneSliceNumber)
{
	m_extractRegion = m_preProcImage->GetLargestPossibleRegion();
	m_extractRegion.SetIndex(0,0);
	m_extractRegion.SetIndex(1,0);
	m_extractRegion.SetIndex(2,0);
	m_numSlices = m_extractRegion.GetSize(sliceDir);
	if(MidPlaneSliceNumber.compare("default")==0)
		m_extractRegion.SetIndex(sliceDir,(m_numSlices-1) / 2);
	else
		m_extractRegion.SetIndex(sliceDir, atoi(MidPlaneSliceNumber.c_str()));
	m_extractRegion.SetSize(sliceDir, 0);
 
	try 
	{
		BinaryExtractFilterType::Pointer binaryExtractFilter = BinaryExtractFilterType::New();
		binaryExtractFilter->SetInput(m_loadMask);
		binaryExtractFilter->SetExtractionRegion(m_extractRegion);
#if  ITK_VERSION_MAJOR >=4
		binaryExtractFilter->SetDirectionCollapseToIdentity();
#endif
		binaryExtractFilter->Update();
		m_mask = binaryExtractFilter->GetOutput();
	}
	catch( itk::ExceptionObject & e )
	{
		std::cerr << "Error: extract_Midsagtital_planes fail : " << e << std::endl;
		exit(1);
	}
}

/***************************************************************************
 * Make an averaging
 ***************************************************************************/
void CCsegtool_initialization::averaging(int averageNum, int sliceDir,  std::string MidPlaneSliceNumber)
{
	try 
	{
		//get first slice
		int slice;
		if(MidPlaneSliceNumber.compare("default")==0)
			slice = (m_numSlices-1) / 2 - averageNum;
		else
			slice = atoi(MidPlaneSliceNumber.c_str()) - averageNum;
		m_extractRegion.SetIndex(sliceDir,slice);
		ExtractFilterType::Pointer extractFilter = ExtractFilterType::New();
		extractFilter->SetInput(m_preProcImage);
		extractFilter->SetExtractionRegion(m_extractRegion);
#if  ITK_VERSION_MAJOR >=4
		extractFilter->SetDirectionCollapseToIdentity();
#endif
		extractFilter->Update();
		m_inputImage = extractFilter->GetOutput();
		m_Region = m_inputImage->GetLargestPossibleRegion();
	}
	catch( itk::ExceptionObject & e )
	{
		std::cerr << "Error: get first slice fail : " << e << std::endl;
		exit(1);
	}
	
	//add other slices
	for (int i = - averageNum + 1; i <= averageNum; i++)
	{
		try 
		{
			if(MidPlaneSliceNumber.compare("default")==0)
				m_extractRegion.SetIndex(sliceDir,(m_numSlices-1) / 2 + i);
			else
				m_extractRegion.SetIndex(sliceDir, atoi(MidPlaneSliceNumber.c_str()) + i);
			m_extractFilter = ExtractFilterType::New();
			m_extractFilter->SetInput(m_preProcImage);
			m_extractFilter->SetExtractionRegion(m_extractRegion);
#if  ITK_VERSION_MAJOR >=4
			m_extractFilter->SetDirectionCollapseToIdentity();
#endif
			m_extractFilter->Update();
		}
		catch( itk::ExceptionObject & e )
		{
			std::cerr << "Error: second slices fail : " << e << std::endl;
			exit(1);
		}
		try 
		{
			AddFilterType::Pointer addFilter = AddFilterType::New();
			addFilter->SetInput1( m_inputImage );
			addFilter->SetInput2( m_extractFilter->GetOutput());
			addFilter->Update();
			m_inputImage = addFilter->GetOutput();
		}
		catch( itk::ExceptionObject & e )
		{
			std::cerr << "Error: add other slices fail : " << e << std::endl;
			exit(1);
		}
	}
	
	//divide summation slice with 2* averageNum + 1
	try 
	{
		ScaleFilterType::Pointer scaleFilter = ScaleFilterType::New();
		scaleFilter->SetInput(m_inputImage);
		scaleFilter->SetShift(0); // no shift
		scaleFilter->SetScale(1.0/(1.0 + 2*averageNum)); // division by averageNum +1
		scaleFilter->Update();
		m_inputImage = scaleFilter->GetOutput();
	}
	catch( itk::ExceptionObject & e )
	{
		std::cerr << "Error: divide summation slice with 2* averageNum + 1 fail : " << e << std::endl;
		exit(1);
	}
}

/***************************************************************************
 * Rotation 90 of the input image
 ***************************************************************************/
void CCsegtool_initialization::dopermute_x_y()
{
	PermuteAxesImageFilterType::PermuteOrderArrayType order;
	order[0] = 1;
	order[1] = 0;

	try 
	{
		PermuteAxesImageFilterType::Pointer permuteFilter = PermuteAxesImageFilterType::New();
		permuteFilter->SetInput(m_inputImage);
		permuteFilter->SetOrder(order);
		permuteFilter->Update();
		m_inputImage = permuteFilter->GetOutput();
	}
	catch( itk::ExceptionObject & e )
	{
		std::cerr << "Error: Permut axes image fail : " << e  << std::endl;
		exit(1);
	}
	
	try 
	{
		BinaryPermuteAxesImageFilterType::Pointer binaryPermuteFilter = BinaryPermuteAxesImageFilterType::New();
		binaryPermuteFilter->SetInput(m_mask);
		binaryPermuteFilter->SetOrder(order);
		binaryPermuteFilter->Update();
		m_mask = binaryPermuteFilter->GetOutput();
	}
	catch( itk::ExceptionObject & e )
	{
		std::cerr << "Error: Permut axes mask fail : "<< e << std::endl;
		exit(1);
	}
}

/***************************************************************************
 * Reflection on X 
 ***************************************************************************/
void CCsegtool_initialization::reflectX()
{
	ReflectFilterType::FlipAxesArrayType axes;
	axes[0] = true;
	axes[1] = false;
	
	try 
	{
		ReflectFilterType::Pointer reflectFilter = ReflectFilterType::New();
		reflectFilter->SetInput(m_inputImage);
		reflectFilter->SetFlipAxes(axes);
		reflectFilter->Update();
		m_inputImage = reflectFilter->GetOutput();
	}
	catch( itk::ExceptionObject & e )
	{
		std::cerr << "Error: Reflect filter image fail : " << e << std::endl;
		exit(1);
	}
	
	try 
	{
		BinaryReflectFilterType::Pointer binaryReflectFilter = BinaryReflectFilterType::New();
		binaryReflectFilter->SetInput(m_mask);
		binaryReflectFilter->SetFlipAxes(axes);
		binaryReflectFilter->Update();
		m_mask = binaryReflectFilter->GetOutput();
	}
	catch( itk::ExceptionObject & e )
	{
		std::cerr << "Error: Reflect filter mask fail : " << e << std::endl;
		exit(1);
	}
}

/***************************************************************************
 * Reflection on Y
 ***************************************************************************/
void CCsegtool_initialization::reflectY()
{
	ReflectFilterType::FlipAxesArrayType axes;
	axes[0] = false;
	axes[1] = true;
	
	try 
	{
		ReflectFilterType::Pointer reflectFilter = ReflectFilterType::New();
		reflectFilter->SetInput(m_inputImage);
		reflectFilter->SetFlipAxes(axes);
		reflectFilter->Update();
		m_inputImage = reflectFilter->GetOutput();
	}
	catch( itk::ExceptionObject & e )
	{
		std::cerr << "Error: Reflect filter image fail : " << e << std::endl;
		exit(1);
	}
	
	try 
	{
		BinaryReflectFilterType::Pointer binaryReflectFilter = BinaryReflectFilterType::New();
		binaryReflectFilter->SetInput(m_mask);
		binaryReflectFilter->SetFlipAxes(axes);
		binaryReflectFilter->Update();
		m_mask = binaryReflectFilter->GetOutput();
	}
	catch( itk::ExceptionObject & e )
	{
		std::cerr << "Error: Reflect filter mask fail : " << e << std::endl;
		exit(1);
	}
}

/***************************************************************************
 * Extract the Label
 ***************************************************************************/
void CCsegtool_initialization::extractLabel(bool segLabel)
{
	try 
	{
		int seglabel=1;
		if(segLabel)
			seglabel=0;
		m_threshFilter = BinaryThresholdFilterType::New();
		m_threshFilter->SetInput(m_mask);
		m_threshFilter->SetUpperThreshold(seglabel);
		m_threshFilter->SetLowerThreshold(seglabel);
		m_threshFilter->SetOutsideValue(BG_VALUE);
		m_threshFilter->SetInsideValue(LABEL_VALUE);
		m_threshFilter->Update();
	}
	catch( itk::ExceptionObject & e )
	{
		std::cerr << "Error: threshfilter fail : " << e << std::endl;
		exit(1);
	}
}

/***************************************************************************
 * Do an closing
 ***************************************************************************/
void CCsegtool_initialization::closingop()
{
	int elemSize = 1;
	StructuringElementType structuringElement;
	structuringElement.SetRadius( elemSize );  // 3x3x3 structuring element
	structuringElement.CreateStructuringElement( );
	m_dilateFilter = dilateFilterType::New();
	m_erodeFilter = erodeFilterType::New();
	m_dilateFilter->SetKernel( structuringElement );
	m_erodeFilter->SetKernel( structuringElement ); 
	m_dilateFilter->SetDilateValue (LABEL_VALUE);
	m_erodeFilter->SetErodeValue (LABEL_VALUE); 
}

/***************************************************************************
 * Double the size of the image
 ***************************************************************************/
void CCsegtool_initialization::dodoubleOn()
{
	TransformType::Pointer transform = TransformType::New();
	transform->SetIdentity();
	ExtractImageRegionType extractregion;
	extractregion = m_inputImage->GetLargestPossibleRegion();
	extractregion.SetSize(0, extractregion.GetSize(0) * 2 );
	extractregion.SetSize(1, extractregion.GetSize(1) * 2 );
	double spacing[2] ;
	spacing[0] = (m_inputImage->GetSpacing())[0] / 2;
	spacing[1] = (m_inputImage->GetSpacing())[1] / 2;
	double maskSpacing[2] ;
	maskSpacing[0] = (m_mask->GetSpacing())[0] / 2;
	maskSpacing[1] = (m_mask->GetSpacing())[1] / 2;

	SplineType::Pointer splineInterpolator = SplineType::New();
	splineInterpolator->SetSplineOrder(5);
	NNType::Pointer NNInterpolator = NNType::New();

	try 
	{
		ResampleFilterType::Pointer resampleFilter = ResampleFilterType::New();
		resampleFilter->SetInput(m_inputImage);
		resampleFilter->SetTransform(transform);
		if(!m_interpolationlinear)
			resampleFilter->SetInterpolator(splineInterpolator);
		resampleFilter->SetDefaultPixelValue( BG_VALUE);  
		resampleFilter->SetOutputOrigin(m_inputImage->GetOrigin());
		resampleFilter->SetSize(extractregion.GetSize());
		resampleFilter->SetOutputSpacing(spacing);
		resampleFilter->Update();
		m_inputImage = resampleFilter->GetOutput();
	}
	catch( itk::ExceptionObject & e )
	{
		std::cerr << "Error: Resample filter image fail : " << e << std::endl;
		exit(1);
	}

	try 
	{
		BinaryResampleFilterType::Pointer binaryResampleFilter = BinaryResampleFilterType::New();
		binaryResampleFilter->SetInput(m_mask);
		binaryResampleFilter->SetTransform(transform);
		binaryResampleFilter->SetInterpolator(NNInterpolator); 
		binaryResampleFilter->SetDefaultPixelValue(BG_VALUE);  
		binaryResampleFilter->SetOutputOrigin(m_mask->GetOrigin());
		binaryResampleFilter->SetSize(extractregion.GetSize());
		binaryResampleFilter->SetOutputSpacing(maskSpacing);
		binaryResampleFilter->Update();
		m_mask = binaryResampleFilter->GetOutput();
	}
	catch( itk::ExceptionObject & e )
	{
		std::cerr << "Error: Resample filter mask fail : " << e << std::endl;
		exit(1);
	}
}

/***************************************************************************
 * Write input: inSlice.gipl is the slice (2D-image) of the input image
 *              inMask.gipl is the 2D-image of the mask image
 ***************************************************************************/
void CCsegtool_initialization::writeInput(std::string outfileBase, std::string nameofproject )
{
	std::string outfileName;
	try 
	{
		m_cast3DFilter = Cast2d3dFilterType::New();
		m_cast3DFilter->SetInput(m_inputImage);
		m_cast3DFilter->Update();
	}
	catch( itk::ExceptionObject & e )
	{
		std::cerr << "Error: Cast2d3dFilter mask fail : " << e << std::endl;
		exit(1);
	}
	ImageWriterType::Pointer writer = ImageWriterType::New();
	outfileName = outfileBase + "/" + nameofproject +"_inSlice.gipl";
	writer->SetFileName(outfileName.c_str()); 
	writer->SetInput(m_cast3DFilter->GetOutput());
	writer->Write();
	
	try 
	{
		m_bincast3DFilter = BinaryCast2d3dFilterType::New();
		m_bincast3DFilter->SetInput(m_mask);
		m_bincast3DFilter->Update();
	}
	catch( itk::ExceptionObject & e )
	{
		std::cerr << "Error: Cast2d3dFilter mask fail : " << e << std::endl;
		exit(1);
	}
	
	BinaryWriterType::Pointer binaryWriter = BinaryWriterType::New();
	outfileName = outfileBase + "/" + nameofproject + "_inMask.gipl";
	binaryWriter->SetFileName(outfileName.c_str()); 
	binaryWriter->SetInput(m_bincast3DFilter->GetOutput());
	binaryWriter->Write();
}

/***************************************************************************
 * Write output: param.txt is a file with the value of the parameters
 *               CCslice.png is the input image which could be read by Qt
 ***************************************************************************/
void CCsegtool_initialization::writeoutput(std::string outfileBase, std::string nameofproject)
{
	std::string outfileParamName;
	outfileParamName = outfileBase + "/" + nameofproject + "_param.txt";
	std::ofstream cfile(outfileParamName.c_str(), std::ios::out);
	cfile << "CenterX:" << GetCenterX() << std::endl;
	cfile << "CenterY:" << GetCenterY() << std::endl;
	cfile << "Scale:" << GetScale() << std::endl;
	cfile << "Rotation:" << GetRotation() << std::endl;
	cfile << "WMvalue:" << GetWMvalue() << std::endl;
	cfile.close();
	
	try 
	{
		m_uncharoutput = CastPNGFilterType::New();
		m_uncharoutput->SetInput(m_extractImage);
		m_uncharoutput->Update();
	}
	catch( itk::ExceptionObject & e )
	{
		std::cerr << "Error: Cast2d3dFilter mask fail : " << e << std::endl;
		exit(1);
	}
	
	std::string imageFileName;
	ImageWriterType3::Pointer writer = ImageWriterType3::New();
	imageFileName = outfileBase + "/" + nameofproject +  "_CCslice.png";
	writer->SetFileName(imageFileName.c_str()); 
	writer->SetInput(m_uncharoutput->GetOutput());
	writer->UseCompressionOff();
	writer->Write();
}

/***************************************************************************
 * Comput parameters: calcul and save the value of the parameters
 ***************************************************************************/
int CCsegtool_initialization::compute_parameters(void)
{
	// get largest component in mask only? this can cause troubles if
	// the brain-stem or the cerebellum is larger than the CC, the alternative is
	// to take the first 3 largest components and select the one,
	// which is closer to the image center (heuristic...)
	
	CCLFilterType::Pointer CCLFilter = CCLFilterType::New();
	CCLFilter->SetInput(m_extractMask);
	
	//do erosion with star element
	RelabelFilterType::Pointer RelabelFilter = RelabelFilterType::New();
	RelabelFilter->SetInput(CCLFilter->GetOutput());
	ExtractBinaryImageType::Pointer comp1,comp2,comp3;
	MomentsCalcType::VectorType firstMomentsC1,firstMomentsC2,firstMomentsC3;

	
	try 
	{
		m_thresh1Filter= threshCCLFilterType::New();
		m_thresh1Filter->SetInput(RelabelFilter->GetOutput());
		m_thresh1Filter->SetUpperThreshold(1); // only the largest one
		m_thresh1Filter->SetLowerThreshold(1);
		m_thresh1Filter->SetOutsideValue( BG_VALUE );
		m_thresh1Filter->SetInsideValue( LABEL_VALUE );
		m_thresh1Filter->Update();

		comp1 = m_thresh1Filter->GetOutput();
		MomentsCalcType::Pointer momentCalcC1 = MomentsCalcType::New();
		momentCalcC1->SetImage(comp1);
		momentCalcC1->Compute();
		firstMomentsC1 = momentCalcC1->GetFirstMoments();
		
		m_thresh2Filter = threshCCLFilterType::New();
		m_thresh2Filter->SetInput(RelabelFilter->GetOutput());
		m_thresh2Filter->SetUpperThreshold(2); // only the largest one
		m_thresh2Filter->SetLowerThreshold(2);
		m_thresh2Filter->SetOutsideValue( BG_VALUE );
		m_thresh2Filter->SetInsideValue( LABEL_VALUE );
		m_thresh2Filter->Update();

		comp2 = m_thresh2Filter->GetOutput();
		MomentsCalcType::Pointer momentCalcC2 = MomentsCalcType::New();
		momentCalcC2->SetImage(comp2);
		momentCalcC2->Compute();
		firstMomentsC2 = momentCalcC2->GetFirstMoments();

		m_thresh3Filter = threshCCLFilterType::New();
		m_thresh3Filter->SetInput(RelabelFilter->GetOutput());
		m_thresh3Filter->SetUpperThreshold(3); // only the largest one
		m_thresh3Filter->SetLowerThreshold(3);
		m_thresh3Filter->SetOutsideValue( BG_VALUE );
		m_thresh3Filter->SetInsideValue( LABEL_VALUE );
		m_thresh3Filter->Update();

		comp3 = m_thresh3Filter->GetOutput();
		MomentsCalcType::Pointer momentCalcC3 = MomentsCalcType::New();
		momentCalcC3->SetImage(comp3);
		momentCalcC3->Compute();
		firstMomentsC3 = momentCalcC3->GetFirstMoments();
	}
	catch( itk::ExceptionObject & e )
	{
		std::cerr << "Error: Moment calculation error, empty binary segmentation?? : " << e << std::endl;
		return(1);
	}
	
	// compare centers
	double imageCenter[2];
	ExtractImageType::RegionType extractregion = m_extractImage->GetLargestPossibleRegion();
	//ExtractImageType::PointType origin =  m_extractImage->GetOrigin();
	// Assuming ACPC alignment the CC is slightly superiorly to the image center
	imageCenter[0] = (extractregion.GetSize())[0] / 2;
	imageCenter[1] = (extractregion.GetSize())[1] / 2;
	double distC1 = sqrt((imageCenter[0] - firstMomentsC1[0]) * 
				(imageCenter[0] - firstMomentsC1[0]) + 
				(imageCenter[1] - firstMomentsC1[1]) * 
				(imageCenter[1] - firstMomentsC1[1]));
	double distC2 = sqrt((imageCenter[0] - firstMomentsC2[0]) * 
				(imageCenter[0] - firstMomentsC2[0]) + 
				(imageCenter[1] - firstMomentsC2[1]) * 
				(imageCenter[1] - firstMomentsC2[1]));
	double distC3 = sqrt((imageCenter[0] - firstMomentsC3[0]) * 
				(imageCenter[0] - firstMomentsC3[0]) + 
				(imageCenter[1] - firstMomentsC3[1]) * 
				(imageCenter[1] - firstMomentsC3[1]));
	//std::cout << "dists: " << distC1 << ", " << distC2 << ", "<< distC3 << std::endl;
	//std::cout << (extractregion.GetSize())[0] << ", " << (extractregion.GetSize())[1] << std::endl;
	//std::cout << imageCenter[0] << ", " << imageCenter[1] << std::endl;
	if (distC1 < distC2) {
	  if (distC1 < distC3) {
	    if (!m_othercompo) {
		m_extractMask = comp1;
		std::cout << "choosing first" << std::endl;
	    } else {
	      if (distC2 < distC3) {
	    	m_extractMask = comp2;
		std::cout << "choosing second" << std::endl;
	      } else {
		m_extractMask = comp3;
		std::cout << "choosing third" << std::endl;
	      }
	    }
	  } else {
	    if (!m_othercompo) {
		m_extractMask = comp3;
		std::cout << "choosing third" << std::endl;
	    } else {
	    	m_extractMask = comp1;
		std::cout << "choosing first" << std::endl;
	    }
	  }
	} else {
	  if (distC2 < distC3) {
	    if (!m_othercompo) {
	      m_extractMask = comp2;
		std::cout << "choosing second" << std::endl;
	    } else {
	      if (distC1 < distC3) {
		m_extractMask = comp1;
		std::cout << "choosing first" << std::endl;
	      } else {
		m_extractMask = comp3;
		std::cout << "choosing third" << std::endl;
	      }
	    }
	  }else {
	    if (!m_othercompo) {
		m_extractMask = comp3;
		std::cout << "choosing third" << std::endl;
	    } else {
	    	m_extractMask = comp2;
		std::cout << "choosing second" << std::endl;
	    }
	  }
	}

	
	// Moments are needed for initialization of pose
	MomentsCalcType::Pointer momentCalc = MomentsCalcType::New();
	momentCalc->SetImage(m_extractMask);
	momentCalc->Compute();
	MomentsCalcType::VectorType center = momentCalc->GetCenterOfGravity();
	MomentsCalcType::ScalarType mass = sqrt(momentCalc->GetTotalMass());
	MomentsCalcType::VectorType firstMoments = momentCalc->GetFirstMoments();

	// COG is in physical coordinates, first moment is COG in index coordinates
	MomentsCalcType::MatrixType principalAxes = momentCalc->GetPrincipalAxes();
	// angle is angle between x-axis and 2nd principal axes (or y-axis and 1st)
	// principal axes are normalized
	// With a flip, change the rotation value by 2*Pi - init value
	double anglePrincAxes =  0;
	if(m_reflectXOn || m_reflectYOn || m_permute_x_y)
	{
		if(m_reflectXOn)
		{
			anglePrincAxes = acos (principalAxes[0][1]);
			anglePrincAxes = 2*3.14159265 - anglePrincAxes;
		}
		if(m_reflectYOn)
		{
			anglePrincAxes = acos (principalAxes[0][0]);
			anglePrincAxes = 2*3.14159265 - anglePrincAxes;
		}
		if(m_permute_x_y)
			anglePrincAxes = -acos (principalAxes[0][0]);
	}
	else
		anglePrincAxes = -acos (principalAxes[0][1]);
	
	m_CenterX = (int) firstMoments[0];
	m_CenterY = (int)  firstMoments[1];
	m_Scale =  mass * m_scalefactor; // sum of pixel values (all pixels==1) -> pixel area
	
	// Rotation with an additional value
	switch(m_angle)
	{
		case 0:
			m_Rotation = anglePrincAxes;
			break;
		case 1:
			m_Rotation = anglePrincAxes - 3.14159265/2.0;
			break;
		case 2:
			m_Rotation = anglePrincAxes + 3.14159265;
			break;
		case 3 :
			m_Rotation = anglePrincAxes + 3.14159265/2.0;
			break;
	}
	
	// Compute histogram for initialization of White matter mean value
	// Mask the input image
	try 
	{
		m_maskFilter = MaskImageType::New();
		m_maskFilter->SetInput1(m_extractImage);
		m_maskFilter->SetInput2(m_extractMask);
		m_maskFilter->Update();
	}
	catch( itk::ExceptionObject & e )
	{
		std::cerr << "Error: Mask the input image fail : " << e << std::endl;
		exit(1);
	}
	ExtractImageType::Pointer maskedImage = m_maskFilter->GetOutput();
	
	// MinMax computation
	minMaxCalcType::Pointer minMaxCalc = minMaxCalcType::New();
	minMaxCalc->SetImage(maskedImage);
	minMaxCalc->Compute();
	InputPixelType maxval = minMaxCalc->GetMaximum();
	InputPixelType minval = minMaxCalc->GetMinimum();
	int numBins =  maxval - minval + 1;
	
	// Histogram computation
	HistogramType::SizeType size;
	size.SetSize(1);
	size[0] = numBins;

	HistogramType::MeasurementVectorType minValVector, maxValVector;
	minValVector.SetSize(1);
	maxValVector.SetSize(1);
	minValVector[0] = minval;
	maxValVector[0] = maxval + 1;
	
	HistogramType::Pointer histogram = HistogramType::New();
	histogram->SetMeasurementVectorSize(1);
	histogram->Initialize( size, minValVector, maxValVector );
	
	// put each image pixel into the histogram
	HistogramType::MeasurementVectorType measurement;
	measurement.SetSize(1);
	Iterator iter (maskedImage, maskedImage->GetBufferedRegion());
	while ( !iter.IsAtEnd() )
	{
		InputPixelType value = iter.Get();
		if (value > BG_VALUE) { // is pixel in mask
			measurement[0] = value;
#if  ITK_VERSION_MAJOR >=4
			histogram->IncreaseFrequencyOfMeasurement( measurement , 1 );
#else
			histogram->IncreaseFrequency( measurement , 1 );
#endif
		}
		++iter;
	}
	
	//get mean value from histogram
	ExtractImageType::RegionType imageRegion = maskedImage->GetBufferedRegion();
	double numVoxels = 0;
	double sumIntensities = 0;
	HistogramType::Iterator histoIter(histogram);
	HistogramType::IndexType index;
	
	for (histoIter = histogram->Begin() ; histoIter != histogram->End() ; ++histoIter) 
	{
		numVoxels = numVoxels + histoIter.GetFrequency();
		index = histogram->GetIndex(histoIter.GetInstanceIdentifier());
		sumIntensities = sumIntensities + (histogram->GetHistogramMaxFromIndex(index)[0] +
				histogram->GetHistogramMinFromIndex(index)[0]) / 2 *histoIter.GetFrequency(); 
	}
	double meanVal = sumIntensities / numVoxels;
	// get stdev from histogram
	
	
	m_WMvalue = (InputPixelType) meanVal;
	
	return 0;
}


//intput
ExtractImageType::Pointer  CCsegtool_initialization::GetImage()
{
	return m_inputImage;
}

//output
int CCsegtool_initialization::GetCenterX()
{
	return m_CenterX;
}

int CCsegtool_initialization::GetCenterY()
{
	return m_CenterY;
}

double CCsegtool_initialization::GetRotation()
{
	return m_Rotation;
}

double CCsegtool_initialization::GetScale()
{
	return m_Scale;
}

InputPixelType CCsegtool_initialization::GetWMvalue()
{
	return m_WMvalue;
}

ExtractImageConstPointer CCsegtool_initialization::GetOutputImage()
{
	return m_extractImage;
}

float CCsegtool_initialization::GetvoxelsizeX()
{
	return m_voxelsizeX;
}

float CCsegtool_initialization::GetvoxelsizeY()
{
	return m_voxelsizeY;
}

int CCsegtool_initialization::GetImageSSize()
{
	return m_numSlices;
}

SizeType CCsegtool_initialization::Get3DImageSize()
{
	return m_imagesize;
}

SpacingType CCsegtool_initialization::Get3DImageSpacing()
{
	return m_imagespacing;
}

InputImageType::Pointer  CCsegtool_initialization::Get3DImage()
{
	return m_loadImage;
}

