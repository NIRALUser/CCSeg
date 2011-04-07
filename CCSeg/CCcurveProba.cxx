// ************************************************************************
// *************************************************************************
// Probability of the Corpus Collosum
// *************************************************************************
// *************************************************************************

#include "CCcurveProba.h"

CCcurveProba::CCcurveProba(InputImageType::Pointer Image3D, vtkSmartPointer<vtkPoints> points ,
			   std::string CCAtlasDirectory, std::string outputfile, std::string nameofproject, 
			   bool vesselRemoveOn, int segLabel, int averageNum, bool permute_x_y, bool reflectXOn, 
			   bool reflectYOn, bool openOn, bool doubleOn, int sliceDir, int sliceNum, int curvenumpts, 
			   bool Debug)
{
	debug = Debug;
	m_Image3D = Image3D;
	m_curvenumpts = curvenumpts;
	m_points = points;
	m_outputfile = outputfile;
	m_CCAtlasDirectory = CCAtlasDirectory;
	m_nameofproject = nameofproject;
	m_vesselRemoveOn = vesselRemoveOn;
	m_segLabel = segLabel;
	m_averageNum = averageNum;
	m_permute_x_y = permute_x_y;
	m_reflectXOn = reflectXOn;
	m_reflectYOn = reflectYOn;
	m_openOn = openOn;
	m_doubleOn = doubleOn;
	m_sliceDir = sliceDir;
	m_sliceNum = sliceNum;
}

void CCcurveProba::compute_proba()
{
	if(debug)
		std::cout<<" ---Compute Probability Model---"<<std::endl;
	
	//load probability maps
	std::string probFile = m_CCAtlasDirectory + "/meanCC_probModel.txt";
	
	if ( loadProbaMaps(probFile)){
		std::cout << "Error loading probabiblity maps" << std::endl;
		exit(-1);
	}
	
	if (m_numPts != m_curvenumpts) {
		std::cout<<"Mismatch number of contour points between contour and probability map "<<m_curvenumpts
				<< "," << m_numPts << std::endl;
		exit(-1);
	}
	ImagePointer  * probImage  = new  ImagePointer[m_numParc];
	
	try{
		//create the binary image
		PolygonSpatialObjectType::Pointer poly = PolygonSpatialObjectType::New(); 
		InputImageType::Pointer digiImage;
		 // create a PolygonSpatialObject from the curve
		poly->GetProperty()->SetName("sliceBasedPolygonIn3D");
		poly->SetId(1);
		//Thickness
		double thickness = 0;
		for(int i=0;i<ImageIDimension;i++)
		{
			if(m_Image3D->GetSpacing()[ImageIDimension]>thickness)
				thickness = m_Image3D->GetSpacing()[ImageIDimension];
		}
		
		ImagePointType point3Dphysical;
		ContinuousIndexType contindex;
		for (int index = 0; index < m_curvenumpts; index ++) {
			PointType point(m_points->GetPoint(index)); // Is not a spatial object point !
			poly->AddPoint(point);
		}
		PointType point(m_points->GetPoint(0));
		// close polygon by adding again first point
		poly->AddPoint(point);
		
		/* Print point */
		if (debug) {
			std::cout<<"Points :"<< std::endl;
			for (int index = 0; index < m_curvenumpts; index ++)
				std::cout<< poly->GetPoint(index)->GetPosition ()<<std::endl;
		}
		
		DigitizePolygonFilterType::Pointer digiFilter = DigitizePolygonFilterType::New();
		digiFilter->SetInput(poly);
		digiFilter->SetSpacing(m_Image3D->GetSpacing());
		digiFilter->SetSize(m_Image3D->GetLargestPossibleRegion().GetSize());
		digiFilter->SetOrigin(m_Image3D->GetOrigin());
		digiFilter->SetDirection(m_Image3D->GetDirection());
		digiFilter->Update();
		digiImage = digiFilter->GetOutput();
		
		if (debug) {
			std::cout<<"DigiFilter finished"<< std::endl;
		}
		
		//Initialize each slide to 0 except the numslice
		IteratorType iter(digiImage, digiImage->GetLargestPossibleRegion());
		iter.GoToBegin() ;
		while (!iter.IsAtEnd()) {
			if ( iter.GetIndex()[m_sliceDir] != m_sliceNum) {
				digiImage->SetPixel(iter.GetIndex(),0);
			}
			++iter;
		}
		
		ImageWriterType::Pointer writer = ImageWriterType::New();
		std::string BinaryImagefilename = m_outputfile + "/" + m_nameofproject + "_Digicurve.nrrd";
		writer->SetFileName(BinaryImagefilename.c_str());
		writer->SetUseCompression(true);
		writer->SetInput(digiImage);
		writer->Write();
		
		//duplicate the first image in 4 new images
		// copy the filled contour numParc images
		for (int parc = 0; parc < m_numParc; parc++) {
			DuplicatorType::Pointer duplicator = DuplicatorType::New();
			
			duplicator->SetInputImage(digiImage);
			duplicator->Update();
			probImage[parc] = duplicator->GetOutput(); 
		}
		
		 // For all points inside the filled contour, find closest point on curve 
		iter.GoToBegin() ;
		while (!iter.IsAtEnd()) {
			if ( iter.Get() > 0) {
				// find closest point to current position
				PointType point;
				const InputImageIndex imgIndex = iter.GetIndex();
				digiImage->TransformIndexToPhysicalPoint(imgIndex,point);
				double pointD[ImageIDimension];
				for (int dim = 0; dim < ImageIDimension; dim++ ) {
					pointD[dim] = point[dim];
				}
				int closestIndex = 0;
				double closestDistance = squaredDistance(m_points->GetPoint(closestIndex), pointD);
				for (int index = 0; index < m_numPts; index ++) {
					double newDist = squaredDistance(m_points->GetPoint(index), pointD);
					if (closestDistance > newDist) {
						closestDistance = newDist;
						closestIndex = index;
					}
				}
				// for each subdivision, set the image value to the corresponding probability
				const InputImageIndex imgIndex2 = iter.GetIndex();
				for (int parc = 0; parc < m_numParc; parc++) {
					double probVal = m_probmaps [ parc * m_numPts + closestIndex ];
					(probImage[parc])->SetPixel(imgIndex2, static_cast<short int>(probVal*1000) );
					 // imprint curve outline, multiply by 1000 because the pixeltype is short int
				}
			}
			++iter;
		}
		
		ImageWriterType::Pointer Probawriter = ImageWriterType::New();
		/* write the probability maps */
		for (int parc = 0; parc < m_numParc; parc++) {
			
			std::stringstream sparc;
			sparc << parc + 1;
			std::string Probafilename = m_outputfile+"/"+m_nameofproject+"_probCurve_"+sparc.str()+".nrrd";
			Probawriter->SetFileName(Probafilename.c_str());
			Probawriter->SetUseCompression(true);
			Probawriter->SetInput(probImage[parc]);
			Probawriter->Write();
		}
		
	}
	catch (itk::ExceptionObject e)  {
		std::cerr << "Error: Probability model fail : " << e << std::endl;
		exit(-1) ;
	}
	
	// initiliaze area values
	double * area = new double [m_numParc];
	for (int parc = 0; parc < m_numParc; parc++) {
		area [parc] = 0.0;
	}
	
	// sum the probability maps
	for (int parc = 0; parc < m_numParc; parc++) {
		IteratorType iter(probImage[parc], (probImage[parc])->GetLargestPossibleRegion()); 
		iter.GoToBegin() ;
		while (!iter.IsAtEnd()) {
			area[parc] = area[parc] + iter.Get();
			++iter;
		}
	}
	
	// scale the areas to correct size
	double real_scale;
	real_scale = 1.0;
	for (int dim = 0; dim < ImageIDimension; dim++ ) {
		real_scale = real_scale * m_Image3D->GetSpacing()[dim];
	}
	for (int parc = 0; parc < m_numParc; parc++) {
		area [parc] = area[parc] * real_scale / 1000; // divided by 1000 to go from short int to float
	}
	
	if (debug) cout << "areas scaled by " << real_scale << endl;
	
	// write probability areas
	std::string Areaprobfilename = m_outputfile+"/"+m_nameofproject+"_areaprob.txt";
	ofstream efile(Areaprobfilename.c_str(), ios::out);
	if (!efile) {
		cerr << "Error: open of file \"" << Areaprobfilename << "\" failed." << endl;
		exit(-1);
	}
	efile << "area = ";
	for (int parc = 0; parc < m_numParc; parc++) {
		efile << area[parc] << " " ;
	}
	efile << endl;
	efile.close();
	
	if(debug)
		std::cout << " ---Fin Probability Model---" << std::endl;
}

int CCcurveProba::loadProbaMaps(std::string probFile)
// loads the probabilities as a single linear double array
// addressing: probmaps [ mapindex * numPts + ptsIndex ]
{
	std::string buffer;
	size_t found;
	
	std::ifstream file(probFile.c_str() , ios::in);  // open in reading
	
	if (!file) {
		std::cout << "probfile " << probFile <<" does not exist" << std::endl;
		return -1;
	}
	
	// Read the first line
	getline(file, buffer);
	found=buffer.find_first_of("=");
	m_numPts = atoi((buffer.substr(found+2,buffer.size())).c_str());
	if (m_numPts == 0)  {
		std::cout << "No points in file ??" << std::endl;
		return -1;
	}
	
	//Second line
	getline(file, buffer);
	found=buffer.find_first_of("=");
	m_numParc = atoi((buffer.substr(found+2,buffer.size())).c_str());
	if (m_numParc == 0)  {
		std::cout << "No subdivisions in file ??" << std::endl;
		return -1;
	}
	
	//fil the probmaps with the values
	m_probmaps = new double [ m_numPts * m_numParc ];
	for (int i = 0; i < m_numParc; i ++) {
		buffer.clear();
		getline(file, buffer);
		// Get the probmaps values after the "="
		// +2 because there is a space after the "="
		found=buffer.find_first_of("=");
		buffer = buffer.substr(found + 2, buffer.size());
		// Find the first space and copy the value before this space,
		// Then erase the last value to buffer and do it again
		for (int p = 0; p < m_numPts; p ++) {
			found=buffer.find_first_of(" ");
			m_probmaps[i*m_numPts + p]  = atof((buffer.substr(0, found)).c_str());
			buffer = buffer.substr(found + 1, buffer.size());
		}
	}
	file.close();
	
	if(debug)
	{
		//Num of pts and parc
		std::cout << "Numpts " << m_numPts << std::endl;
		std::cout << "NumParc " << m_numParc << std::endl;
		//Print ProbMap values
		for (int i = 0; i < m_numParc; i ++) {
			for (int p = 0; p < m_numPts; p ++) {
				std::cout << i << "/" << p << " : " << m_probmaps[i*m_numPts + p] << std::endl;
			}
		}
	}
	
	return 0;
}

/****************************************************************************************
* Distance between two points in 3D
*****************************************************************************************/
double CCcurveProba::squaredDistance (const double * point1, const double * point2)
{
	double distance = 0.0;
	for (int dim = 0; dim < 3; dim++) {
		distance = distance +  (point1[dim] - point2[dim]) * (point1[dim] - point2[dim]);
	}
	return distance;
}

