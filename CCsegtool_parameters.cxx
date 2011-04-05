#include "CCsegtool_parameters.h"

/***************************************************************************** 
 * constructor
 *****************************************************************************/
CCsegtool_parameters::CCsegtool_parameters (ExtractImageConstPointer Image,	//filename of the Image
					    std::string  Way,
					    ImagePointer Image3D,
					    bool         permute_x_y,
					    bool         reflectXOn,
					    bool         reflectYOn,
					    bool         doubleOn,
					    int          sliceDir,
					    int          numSlice,
					    int          X,
					    int          Y, 
					    float        Rotation,
					    float        Scale,
					    bool         FSXForm,
					    double       PSDistance,
					    int          WMIntensity,
					    int          MPSDisplacement,
					    int          Oiter,
					    SizeType     ImageSize,
					    bool         Debug)
{
	/* Init param */
	debug=Debug;
	m_permute_x_y = permute_x_y;
	m_reflectXOn = reflectXOn;
	m_reflectYOn = reflectYOn;
	m_doubleOn = doubleOn;
	m_sliceDir = sliceDir;
	m_numSlice = numSlice;
	m_3DImage = Image3D;
	m_ImageSize = ImageSize;
	
	if(debug)
		std::cout<<"  --Initialization of Param--"<<endl;
	//load the Image
	SetImage(Image);
	
	//Each files have a number
	//SSMean
	readf4(Way,1);
	//Eigenvectors
	readf4(Way,2);
	//Bounds
	readf(Way,3);
	//PMMean
	readf(Way,4);
	//PMSigmaInv
	readf(Way,5);
	
	if(debug)
	{
		//Print input vector size and the first and last elements
		std::cout<<"  SMMean vector size:"<<m_SMMean.size()<<" ,first element "<<m_SMMean[0]<<" ,last element "
				<<m_SMMean[m_SMMean.size()-1]<<endl;
		std::cout<<"  Eigenvectors vector size:"<<m_Eigenvectors.size()<<" ,first element "<<m_Eigenvectors[0]
				<<" ,last element "<<m_Eigenvectors[m_Eigenvectors.size()-1]<<endl;
		std::cout<<"  Bound vector size:"<<m_Bound.size()<<" ,first element "<<m_Bound[0]<<" ,last element "
				<<m_Bound[m_Bound.size()-1]<<endl;
		std::cout<<"  PMMean vector size:"<<m_PMMean.size()<<" ,first element "<<m_PMMean[0]<<" ,last element "
				<<m_PMMean[m_PMMean.size()-1]<<endl;
		std::cout<<"  SigmaInv vector size:"<<m_SigmaInv.size()<<" ,first element "<<m_SigmaInv[0]
				<<" ,last element "<<m_SigmaInv[m_SigmaInv.size()-1]<<endl;
	}
	
	
	//Initialization
	m_X = X;
	m_Y = Y;
	m_Rotation = Rotation;
	m_Scale = Scale;
	m_FSXForm = FSXForm;
	m_Unconstrained = false;
	m_PSDistance = PSDistance;
	m_WMIntensity = WMIntensity;
	m_MPSDisplacement = MPSDisplacement;
	m_Oiter = Oiter;
	
	/* Initialization of coef at 0 */
	m_coefs = new Point4[MAX_COEFS];
	for (unsigned int i=0; i<MAX_COEFS; i++) 
		for (unsigned int j=0; j<4; j++)
			m_coefs[i][j] = 0;
	
	coefsgiven = false;
	
	if(debug)
		std::cout<<"  --Constructor param done--"<<endl;
}



/********************************************************************************* 
 * Fill a vector with the CC atlas directory
 ********************************************************************************/
void CCsegtool_parameters::readf(std::string way, int nbfile)
{
	int i[3];
	float j;
	char c;
	std::string buffer,filename;
	
	switch(nbfile)
	{
		case 3:
			filename=way + "/Shape_eigvals.asc";
			break;
		case 4:
			filename=way + "/Profile_mean.asc";
			break;
		case 5:
			filename=way + "/Profile_sigmainv.asc";
			break;
	}
	
	std::ifstream file(filename.c_str() , ios::in);  // open in reading
	if(debug)
		std::cout<< "  " << filename <<endl;
	
	if(file)  // if open
	{
		//Ignoring the first line
		getline(file, buffer);
		
		switch (nbfile)
		{
			case 3:
				//set the size
				file>>i[0]>>c;
				getsizeof(i, nbfile);
				//set the value of the vector
				for(int k=0; k<m_Bound_xsize;k++)
				{
					file>>j;
					m_Bound.push_back(j);
				}
				break;
			case 4:
				//set the size
				file>>i[0]>>i[1]>>c;
				getsizeof(i, nbfile);
				//set the value of the vector
				for(int k=0; k<m_PMMean_xsize*m_PMMean_ysize;k++)
				{
					file>>j;
					m_PMMean.push_back(j);
				}
				break;
			case 5:
				//set the size
				file>>i[0]>>i[1]>>i[2]>>c;
				getsizeof(i, nbfile);
				//set the value of the vector
				for(int k=0; k<m_SigmaInv_xsize*m_SigmaInv_ysize*m_SigmaInv_zsize;k++)
				{
					file>>j;
					m_SigmaInv.push_back(j);
				}
				break;	
		}
		
		file.close();
	}
	else cerr << "Error: The opening of the file failed" << endl;
}


/********************************************************************************* 
 * Fill a vector with the data of a file where there are 4 value on a line
 * like in CC atlas directory
 ********************************************************************************/
void CCsegtool_parameters::readf4(std::string way, int nbfile)
{
	char c;
	int i[2];
	std::string buffer,filename;
	float j[4];
	
	switch (nbfile)
	{
		case 1:
			filename = way + "/Shape_mean.asc";
			break;
		case 2:
			filename = way + "/Shape_eigvecs.asc";
			break;
	}
	
	std::ifstream file(filename.c_str() , ios::in);  // open in reading
	if(debug)
		std::cout<< "  " << filename <<endl;
	
	if(file)  // if open
	{
		//Ignoring the first line
		getline(file, buffer);
		
		switch (nbfile)
		{
			case 1:
				//set the size
				file>>i[0]>>c;
				getsizeof(i, nbfile);
				//set the value of the vector
				for(int k=0; k<m_SMMean_ysize;k++)
				{
					file>>j[0]>>j[1]>>j[2]>>j[3];
					for(int l=0;l<4;l++)
						m_SMMean.push_back(j[l]);
				}
				break;
			case 2:
				//set the size
				file>>i[0]>>i[1]>>c;
				getsizeof(i, nbfile);
				//set the value of the vector
				for(int k=0; k<m_Eigenvectors_ysize*m_Eigenvectors_zsize;k++)
				{
					file>>j[0]>>j[1]>>j[2]>>j[3];
					for(int l=0;l<4;l++)
						m_Eigenvectors.push_back(j[l]);
				}
				break;
		}
		
		file.close();
	}
	
	else cerr << "Error: The opening of the file failed" << endl;
}



/********************************************************************************* 
 * Get the size of Input vector 
 ********************************************************************************/
void CCsegtool_parameters::getsizeof(int *i, int nbfile)
{
	switch ( nbfile )
	{
		case 1:
			m_SMMean_xsize=4;
			m_SMMean_ysize=i[0];
			break;
		case 2:
			m_Eigenvectors_xsize=4;
			m_Eigenvectors_ysize=i[0];
			m_Eigenvectors_zsize=i[1];
			break;
		case 3:
			m_Bound_xsize=i[0];
			break;
		case 4:
			m_PMMean_xsize=i[0];
			m_PMMean_ysize=i[1];
			break;
		case 5:
			m_SigmaInv_xsize=i[0];
			m_SigmaInv_ysize=i[1];
			m_SigmaInv_zsize=i[2];
			break;
	}
}


/***************************************************************************** 
 * Set the parameters
 *****************************************************************************/
void CCsegtool_parameters::setparam(double X, double Y, double Scale, double Rot)
{
	m_X = (int)X;
	m_Y = (int)Y;
	m_Scale = Scale;
	m_Rotation = Rot;
}

/********************************************************************************* 
 * Transform 2D index in 3D index
 ********************************************************************************/
void CCsegtool_parameters::transform2DIndex_3Dindex(Point2* pts, int numpts, std::string path_output,
		 						    std::string nameofproject,double Xsize, double Ysize)
{
	/* Init */
	m_numpts = numpts;
	ContinuousIndexType * index3D = new ContinuousIndexType[numpts];
	ImagePointType physicalPoint3D;
	double buf=0.0;
	
	/* Reflect on Y */
	if(m_reflectYOn)
	{
		for(int i=0;i<numpts;i++)
			pts[i][1] = static_cast<double>(Ysize) - pts[i][1];
	}
	
	/* Reflect on X */
	if(m_reflectXOn)
	{
		for(int i=0;i<numpts;i++)
			pts[i][0] = static_cast<double>(Xsize) - pts[i][0];
	}
	
	/* Permute X and Y */
	if(m_permute_x_y)
	{
		for(int i=0;i<numpts;i++)
		{
			buf = pts[i][0];
			pts[i][0] = pts[i][1];
			pts[i][1] = buf;
		}
	}
	
	/* Double */
	if(m_doubleOn)
	{
		for(int i=0;i<numpts;i++)
		{
			pts[i][0] =  pts[i][0]/2;
			pts[i][1] =  pts[i][1]/2;
		}
	}
	
	/* 2D to 3D */
	if(m_sliceDir==0)
	{
		for(int i=0;i<numpts;i++)
		{
			index3D[i][0] = m_numSlice;
			index3D[i][1] = pts[i][0];
			index3D[i][2] = pts[i][1];
		}
	}
	else if(m_sliceDir==1)
	{
		for(int i=0;i<numpts;i++)
		{
			index3D[i][0] = pts[i][0];
			index3D[i][1] = m_numSlice;
			index3D[i][2] = pts[i][1];
		}
	}
	else if(m_sliceDir==2)
	{
		for(int i=0;i<numpts;i++)
		{
			index3D[i][0] = pts[i][0];
			index3D[i][1] = pts[i][1];
			index3D[i][2] = m_numSlice;
		}
	}
	
	/* 3DWorldSpace */
	//create a vtkPoints object and store the points in it
	m_points = vtkSmartPointer<vtkPoints>::New();
	vtkSmartPointer<vtkPoints> Points = vtkSmartPointer<vtkPoints>::New();
	for(int i=0;i<numpts;i++)
	{
		ContinuousIndexType contindex;
		for (int j = 0; j < 3; j++)
		{
			contindex[j] = index3D[i][j];
		}
		//index3D is in index space -> conversion needed
		m_3DImage->TransformContinuousIndexToPhysicalPoint(contindex,physicalPoint3D);
		m_points->InsertPoint(i,physicalPoint3D[0],physicalPoint3D[1],physicalPoint3D[2]);
		/* Reflect on Slice X */
		physicalPoint3D[0] = -physicalPoint3D[0];
		/* Reflect on Slice Z */
		physicalPoint3D[1] = -physicalPoint3D[1];
		Points->InsertPoint(i,physicalPoint3D[0],physicalPoint3D[1],physicalPoint3D[2]);
	}
	
	//Polyline
	vtkSmartPointer<vtkPolyLine> polyLine = vtkSmartPointer<vtkPolyLine>::New();
	polyLine->GetPointIds()->SetNumberOfIds(numpts+1);
	for(int i=0;i<numpts;i++)
	{
		polyLine->GetPointIds()->SetId(i,i);
	}
	//add the first point at the end to close the curve
	polyLine->GetPointIds()->SetId(numpts,0);
	
	//Create a cell array to store the lines in and add the lines to it
	vtkSmartPointer<vtkCellArray> cells = vtkSmartPointer<vtkCellArray>::New();
	cells->InsertNextCell(polyLine);
 
	//Create a polydata to store everything in
	vtkSmartPointer<vtkPolyData> polyData = vtkSmartPointer<vtkPolyData>::New();
	//add the points to the dataset
	polyData->SetPoints(Points);
	//add the lines to the dataset
	polyData->SetLines(cells);
	
	//Write the file
	std::string outfilePDMName;
	outfilePDMName = path_output + "/" + nameofproject + "_3D.vtk";
	vtkSmartPointer<vtkPolyDataWriter> writer = vtkSmartPointer<vtkPolyDataWriter>::New();
	writer->SetInput(polyData);
	writer->SetFileName(outfilePDMName.c_str());
	writer->SetFileTypeToBinary();
	writer->Write();
}

/**********************************************************************************************************/

/********************************************************************************* 
 * Functions to access to the parameters
 ********************************************************************************/
//attribut X
void CCsegtool_parameters::SetX(int X)
{
	m_X=X;
}
int CCsegtool_parameters::GetX()
{
	return m_X;
}

//attribut Y
void CCsegtool_parameters::SetY(int Y)
{
	m_Y=Y;
}
int CCsegtool_parameters::GetY()
{
	return m_Y;
}

//attribut Rotation
void CCsegtool_parameters::SetRotation(float Rotation)
{
	m_Rotation=Rotation;
}
float CCsegtool_parameters::GetRotation()
{
	return m_Rotation;
}

//attribut Scale
void CCsegtool_parameters::SetScale(float Scale)
{
	m_Scale=Scale;
}
float CCsegtool_parameters::GetScale()
{
	return m_Scale;
}

//attribut FSXForm
void CCsegtool_parameters::SetFSXForm(bool FSXForm)
{
	m_FSXForm=FSXForm;
}
bool CCsegtool_parameters::GetFSXForm()
{
	return m_FSXForm;
}

//attribut Unconstrained
void CCsegtool_parameters::SetUnconstrained(bool Unconstrained)
{
	m_Unconstrained = Unconstrained;
}
bool CCsegtool_parameters::GetUnconstrained()
{
	return m_Unconstrained;
}

//attribut PSDistance
void CCsegtool_parameters::SetPSDistance(double PSDistance)
{
	m_PSDistance=PSDistance;
}
double CCsegtool_parameters::GetPSDistance()
{
	return m_PSDistance;
}

//attribut WMIntensity
void CCsegtool_parameters::SetWMIntensity(int WMIntensity)
{
	m_WMIntensity=WMIntensity;
}
int CCsegtool_parameters::GetWMIntensity()
{
	return m_WMIntensity;
}

//attribut MPSDisplacement
void CCsegtool_parameters::SetMPSDisplacement(int MPSDisplacement)
{
	m_MPSDisplacement=MPSDisplacement;
}
int CCsegtool_parameters::GetMPSDisplacement()
{
	return m_MPSDisplacement;
}

//attribut Osteps
void CCsegtool_parameters::SetOiter(int Oiter)
{
	m_Oiter=Oiter;
}
int CCsegtool_parameters::GetOiter()
{
	return m_Oiter;
}

//attribut Resetman
void CCsegtool_parameters::SetResetMean(bool _ResetMean)
{
	m_ResetMean=_ResetMean;
}
bool CCsegtool_parameters::GetResetMean()
{
	return m_ResetMean;
}



//Input image
void CCsegtool_parameters::SetImage(ExtractImageConstPointer Image )
{
	m_Image = Image;
}
ExtractImageConstPointer CCsegtool_parameters::GetImage()
{
	return m_Image;
}
/*Input vector*/
//Shape Model Mean
vector<float> CCsegtool_parameters::GetSMMean()
{
	return m_SMMean;
}
float CCsegtool_parameters::GetSMMeanval(int x,int y)
{
	if(x<m_SMMean_xsize && y<m_SMMean_ysize)
		return m_SMMean[x + y*m_SMMean_xsize];
	else
	{
		std::cout<<"Error: error between x and y in GetSMMeanval(x,y)"<<endl;
		return 0.0;
	}
}

//Shape Model Eigenvectors
vector<float> CCsegtool_parameters::GetEigenvectors()
{
	return m_Eigenvectors;
}
float CCsegtool_parameters::GetEigenvectorsval(int x,int y,int z)
{
	if(x<m_Eigenvectors_xsize && y<m_Eigenvectors_ysize && z<m_Eigenvectors_zsize)
		return m_Eigenvectors[x + y*m_Eigenvectors_xsize + z*m_Eigenvectors_xsize*m_Eigenvectors_ysize];
	else
	{
		std::cout<<"Error: error between x and y and z in GetEigenvectorsval(x,y)"<<endl;
		return 0.0;
	}
}

//Shape Model Deformation Bounds
vector<float> CCsegtool_parameters::GetBound()
{
	return m_Bound;
}


//Profile Model Mean
vector<float> CCsegtool_parameters::GetPMMean()
{
	return m_PMMean;
}
float CCsegtool_parameters::GetPMMeanval(int x,int y)
{
	if(x<m_PMMean_xsize && y<m_PMMean_ysize)
		return m_PMMean[x + y*m_PMMean_xsize];
	else
	{
		std::cout<<"Error: error between x and y in GetPMMeanval(x,y)"<<endl;
		return 0.0;
	}
}


//Profile Model SigmaInv
vector<float> CCsegtool_parameters::GetSigmaInv()
{
	return m_SigmaInv;
}
float CCsegtool_parameters::GetSigmaInvval(int x,int y, int z)
{
	if(x < m_SigmaInv_xsize && y < m_SigmaInv_ysize && z < m_SigmaInv_zsize)
		return m_SigmaInv[x + y*m_SigmaInv_xsize + z*m_SigmaInv_xsize*m_SigmaInv_ysize];
	else
	{
		if(x > m_SigmaInv_xsize)
			std::cout<<"Error: error x in GetSigmaInvval(x,y,z)"<<endl;
		if(y > m_SigmaInv_ysize)
			std::cout<<"Error: error y in GetSigmaInvval(x,y,z)"<<endl;
		if(z > m_SigmaInv_zsize)
			std::cout<<"Error: error z in GetSigmaInvval(x,y,z)"<<endl;
		return 0.0;
	}
}


/*Output vector*/
//Segmented Shape Coefs
void CCsegtool_parameters::SetSSCoefs(vector<float> _SSCoefs)
{
	m_SSCoefs.clear();
	for(unsigned int i=0;i<_SSCoefs.size();i++)
		m_SSCoefs.push_back(_SSCoefs[i]);
}
vector<float> CCsegtool_parameters::GetSSCoefs()
{
	return m_SSCoefs;
}

//Image Profiles
void CCsegtool_parameters::SetIProfiles(vector<float> _IProfiles)
{
	m_IProfiles.clear();
	for(unsigned int i=0;i<_IProfiles.size();i++)
		m_IProfiles.push_back(_IProfiles[i]);
}
vector<float> CCsegtool_parameters::GetIProfiles()
{
	return m_IProfiles;
}

//Goodness of Fit
void CCsegtool_parameters::SetGoodnessoffit(vector<double> _Goodnessoffit)
{
	m_Goodnessoffit.clear();
	for(unsigned int i=0;i<_Goodnessoffit.size();i++)
		m_Goodnessoffit.push_back(_Goodnessoffit[i]);
}
vector<double> CCsegtool_parameters::GetGoodnessoffit()
{
	return m_Goodnessoffit;
}

//Per-Profile Shifts
void CCsegtool_parameters::SetPPShifts(vector<float> _PPShifts)
{
	m_PPShifts.clear();
	for(unsigned int i=0;i<_PPShifts.size();i++)
		m_PPShifts.push_back(_PPShifts[i]);
}
vector<float> CCsegtool_parameters::GetPPShifts()
{
	return m_PPShifts;
}

//Eigenloadings
void CCsegtool_parameters::Seteigloads(vector<float> _eigloads)
{
	m_eigloads.clear();
	for(unsigned int i=0;i<_eigloads.size();i++)
		m_eigloads.push_back(_eigloads[i]);
}
void CCsegtool_parameters::Seteigloadsval(int i, float value)
{
	m_eigloads[i]=value;
}
vector<float> CCsegtool_parameters::Geteigloads()
{
	return m_eigloads;
}

/*Size of the input vector*/
int CCsegtool_parameters::GetSMMeansize(int number)
{
	switch(number)
	{
		case 1:
			return m_SMMean_xsize;
		case 2:
			return m_SMMean_ysize;
	}
	return -1;
}
int CCsegtool_parameters::GetEigenvectorssize(int number)
{
	switch(number)
	{
		case 1:
			return m_Eigenvectors_xsize;
		case 2:
			return m_Eigenvectors_ysize;
		case 3:
			return m_Eigenvectors_zsize;
	}
	return -1;
}
int CCsegtool_parameters::GetBoundsize()
{
	return m_Bound_xsize;
}
int CCsegtool_parameters::GetPMMeansize(int number)
{
	switch(number)
	{
		case 1:
			return m_PMMean_xsize;
		case 2:
			return m_PMMean_ysize;
	}
	return -1;
}
int CCsegtool_parameters::GetSigmaInvsize(int number)
{
	switch(number)
	{
		case 1:
			return m_SigmaInv_xsize;
		case 2:
			return m_SigmaInv_ysize;
		case 3:
			return m_SigmaInv_zsize;
	}
	return -1;
}
int CCsegtool_parameters::GetNumPointsProba()
{
	return m_numpts;
}
int CCsegtool_parameters::GetNumSlice()
{
	return m_numSlice;
}

/********************************************************************************* 
 * Write output: filecoef.coef is a file with the coefficient at the end of the
 * segmentation, godnessoffit.gof is a file with the value of the Goodness of fit
 ********************************************************************************/
void CCsegtool_parameters::writeoutput(std::string path_output, std::string _nameofproject)
{
	// Write output
	std::string outfileCoefName;
	outfileCoefName = path_output + "/" + _nameofproject + "_filecoef.coef";
	ofstream afile(outfileCoefName.c_str(), ios::out);
	int num_degrees = GetSSCoefs().size()/4;
	afile << "field 1D 4-vector uniform float"<<endl;
	afile << num_degrees << " ;" << endl;
	for (int i = 0; i < num_degrees; i++) {
		afile << GetSSCoefs()[i*4 + 0] << " " 
				<< GetSSCoefs()[i*4 + 1] << " " 
				<< GetSSCoefs()[i*4 + 2] << " "
				<< GetSSCoefs()[i*4 + 3] << " " << endl;
	}
	afile.close();
	
	std::string outfileGOFName;
	outfileGOFName = path_output + "/" + _nameofproject + "_goodnessoffit.gof";
	ofstream file(outfileGOFName.c_str(), ios::out);
	int num = GetGoodnessoffit().size();
	file << num << endl;
	for (int i = 0; i < num; i++) {
		file << GetGoodnessoffit()[i] <<endl;
	}
	file.close();
}

/********************************************************************************* 
 * Set the coefficient for the actuel curve
 ********************************************************************************/
void CCsegtool_parameters::SetCoefs(Point4* coefs)
{
	coefsgiven = true;
	for(unsigned int i=0;i<m_num_coefs;i++)
	{
		for(unsigned int j=0;j<4;j++)
			m_coefs[i][j]=coefs[i][j];
	}
}

/********************************************************************************* 
 * Get the coefficient for the actuel curve
 ********************************************************************************/
Point4* CCsegtool_parameters::GetCoefs()
{
	return m_coefs;
}

/********************************************************************************* 
 * boolean : true if the coef changed
 ********************************************************************************/
bool CCsegtool_parameters::Getcoefsgiven()
{
	return coefsgiven;
}

/********************************************************************************* 
 * Number of coefficients
 ********************************************************************************/
void CCsegtool_parameters::SetNumcoefs(unsigned int num_coefs)
{
	m_num_coefs = num_coefs;
}

unsigned int CCsegtool_parameters::GetNumcoefs()
{
	return m_num_coefs;
}

/********************************************************************************* 
 * Vtk points
 ********************************************************************************/
vtkSmartPointer<vtkPoints> CCsegtool_parameters::GetPointsvtk()
{
	return m_points;
}

/********************************************************************************* 
 * itk 3D Input Image
 ********************************************************************************/
ImagePointer CCsegtool_parameters::GetImage3D()
{
	return m_3DImage;
}

/********************************************************************************* 
 * Global value from the initialization
 ********************************************************************************/
void CCsegtool_parameters::SetGlobalvalue(short WMvalue, float voxelsizeX, float voxelsizeY)
{
	m_WMvalue = WMvalue;
	m_voxelsizeX = voxelsizeX;
	m_voxelsizeY = voxelsizeY;
}

/********************************************************************************* 
 * White matter intensity value from initialization
 ********************************************************************************/
short CCsegtool_parameters::GetWM()
{
	return m_WMvalue;
}

/********************************************************************************* 
 * X voxel size value from initialization
 ********************************************************************************/
float CCsegtool_parameters::GetVoxelSizeX()
{
	return m_voxelsizeX;
}

/********************************************************************************* 
 * Y voxel size
 ********************************************************************************/
float CCsegtool_parameters::GetVoxelSizeY()
{
	return m_voxelsizeY;
}


