#include "Globalfunc.h"

/**************************************************************************************
 * Initialization of the parameters : calculate the CCsegtool_initialization, parameters
 * and computation, then run the first step for one iteration. 
 * Draw the image with the model.
 *************************************************************************************/
int initialize( std::string Image_filename, std::string Seg_filename, std::string CCAtlasDirectory,
		std::string Path_output, std::string nameofproject, bool interpolationlinear,
		bool vesselRemoveOn, int segLabel, int averageNum, bool permute_x_y, bool reflectXOn, 
		bool reflectYOn, bool openOn, bool doubleOn, int sliceDir, std::string MidPlaneSliceNumber,
		int number_pts, bool FSXForm, double PSDistance, bool Unconstrained, 
		std::string WMintensity, std::string MPSDisplacement, std::string Number_iteration, 
		int Lambdamax, std::string Coefofoptim, bool debug, bool withgui, Imageview * &imageview,
		CCsegtool_parameters* &last_parameters, double *lastX, double *lastY, double *lastScale, double *lastRot,
		float scalefactor, bool othercompo, int angle, bool rot90, bool rot180, bool rot270, 
		QWidget * parent, QLabel* parametersinitview, QLabel* Arealabel, QLabel* Midsaggitalplaneslicevalue)
{
	/* create local variables */
	float disp = 100000000;
	int iter;
	int numberSlice=0;
	Point2 *pts=new Point2[number_pts];
	SizeType ImageSize;
	SpacingType ImageSpacing;
	for(int i=0;i<number_pts;i++)
		pts[i][0]=pts[i][1]=0;
	/* initiate WMvalue, MPSdisplacement, number of iteration and Coefficient of optimization*/
	int WMvalue1,WMvalue2,WMvalue3;
	int MPSDisplacement1,MPSDisplacement2,MPSDisplacement3;
	int Number_iteration1,Number_iteration2,Number_iteration3;
	float coef1,coef2,coef3;
	
	/* Initialization */
	CCsegtool_initialization * initialization = new CCsegtool_initialization(interpolationlinear, scalefactor);
	if(initialization->compute_initialization(Image_filename, Seg_filename, vesselRemoveOn, segLabel, averageNum,
	   permute_x_y, reflectXOn, reflectYOn, openOn, doubleOn, sliceDir, Path_output,
	   nameofproject, MidPlaneSliceNumber, othercompo, angle, debug) == 1)
		return 1;
	
	//Get the number of slices for the display, and set the Probacurve
	if(withgui)
	{
		/* Number of Slice */
		numberSlice = initialization->GetImageSSize();
		std::stringstream SnumberSlice;
		SnumberSlice << numberSlice - 1;
		QString text = "(0 to ";
		text.append((SnumberSlice.str()).c_str());
		text.append(")");
		Midsaggitalplaneslicevalue->setText(text);
		
		/* Proba parameters setting */
		ImageSize = initialization->Get3DImageSize();
		ImageSpacing = initialization->Get3DImageSpacing();
		
		/* If with GUI, with preview and the imageview exists, change the PixelValueLabel */
		if(imageview!=NULL )
		{
			QString imageFileName( (Path_output + "/" + nameofproject +  "_CCslice.png").c_str());
			imageview->setGrayPixmap(imageFileName);
			imageview->show();
		}
		
		/* Show the value of the Init parameters */
		QString param = readparam(Path_output + "/" + nameofproject + "_param.txt");
		parametersinitview->setText( param );
	
		if(imageview==NULL)
		{
			/*creat the pixmap with the output image */
			QString imageFileName( (Path_output + "/" + nameofproject +  "_CCslice.png").c_str());
			Imageview *newimageview = new Imageview(imageFileName, rot90, rot180, rot270, debug,
								parent);
			if((newimageview->Getsize()).width()>(newimageview->Getsize()).height())
				newimageview->move(40,145);
			else	
				newimageview->move(110,85);
			imageview = newimageview;
			imageview->show();
			/* Install the filter */
			imageview->installEventFilter(parent);
		}
	}
	
	/*WMintensity*/
	WMvalue1= (int)(initialization->GetWMvalue()*1.3);
	WMvalue2= (int)(initialization->GetWMvalue());
	WMvalue3= (int)(initialization->GetWMvalue());
	if(WMintensity.compare("default,default,default")!=0)
	{
		WMvalue1=readstringWMI(WMintensity, 1, initialization->GetWMvalue());
		WMvalue2=readstringWMI(WMintensity, 2, initialization->GetWMvalue());
		WMvalue3=readstringWMI(WMintensity, 3, initialization->GetWMvalue());
	}
	
	/*MPSDisplacement*/
	MPSDisplacement1= 10;
	MPSDisplacement2= 3;
	MPSDisplacement3= 2;
	if(MPSDisplacement.compare("10,3,2")!=0)
	{
		MPSDisplacement1=readstringinputi(MPSDisplacement, 1, 2);
		MPSDisplacement2=readstringinputi(MPSDisplacement, 2, 2);
		MPSDisplacement3=readstringinputi(MPSDisplacement, 3, 2);
	}
	
	/*Number of iteration*/
	Number_iteration1= 50;
	Number_iteration2= 5;
	Number_iteration3= 3;
	if(Number_iteration.compare("50,15,3")!=0)
	{
		Number_iteration1=readstringinputi(Number_iteration, 1, 3);
		Number_iteration2=readstringinputi(Number_iteration, 2, 3);
		Number_iteration3=readstringinputi(Number_iteration, 3, 3);
	}
	
	/*Coefficient of optimization*/
	coef1= 0.25* initialization->GetvoxelsizeX();
	coef2= 0.1* initialization->GetvoxelsizeX();
	coef3= 0.05* initialization->GetvoxelsizeX();
	if(Coefofoptim.compare("0.25,0.1,0.05")!=0)
	{
		coef1=readstringinputf(Coefofoptim, 1)* initialization->GetvoxelsizeX();
		coef2=readstringinputf(Coefofoptim, 2)* initialization->GetvoxelsizeX();
		coef3=readstringinputf(Coefofoptim, 3)* initialization->GetvoxelsizeX();
	}
	
	/* MidPlaneSlice Number */
	int numslice = initialization->GetImageSSize();
	if(MidPlaneSliceNumber.compare("default")==0)
		numslice = (numslice-1)/2;
	
	CCsegtool_parameters* parameters = new CCsegtool_parameters(initialization->GetOutputImage(), CCAtlasDirectory,
					      initialization->Get3DImage(), permute_x_y, reflectXOn, reflectYOn, doubleOn,
					      sliceDir, numslice ,initialization->GetCenterX(),initialization->GetCenterY(),
					      initialization->GetRotation(),initialization->GetScale(), FSXForm, PSDistance,
					      WMvalue1,10,50,initialization->Get3DImageSize(), debug);
	
	// Save the value from the initialization in parameters
	parameters->SetGlobalvalue(initialization->GetWMvalue(), initialization->GetvoxelsizeX(),
				    initialization->GetvoxelsizeY());
	
	// computation
	CCsegtool_computation * computation = new CCsegtool_computation(parameters, nameofproject, Path_output, debug,
									0.0, 0.0);
	
	iter=run_step(parameters, computation, imageview, parent, 1, WMvalue1, MPSDisplacement1,
		      Number_iteration1,coef1, pts, number_pts, Lambdamax, disp, debug, withgui, true);
			
	/*write ouput file*/
	parameters->writeoutput(Path_output, nameofproject);
	
	/* Save parameters and image */
	if(withgui)
		computation->getparam(lastX, lastY, lastScale, lastRot);
	
	/* save the parameters and computation */
	/* Create the polygon */
	double ContourArea;
	//calcul the pts
	computation->reconstructCurve(pts,number_pts);
	// the Array
	Point2D *ContourArray = new Point2D[number_pts];;
	for(int i=0;i<number_pts;i++)
	{
		ContourArray[i].x = pts[i][0];
		ContourArray[i].y = pts[i][1];
	}
	// the polygon and add the number of points with the array
	Polygon2D Contour;
	Contour.numPoints = number_pts;
	Contour.pointArray = ContourArray;
	ContourArea = contourArea(Contour);
	
	// Multiply by the size of a voxel
	ContourArea *= initialization->GetvoxelsizeX()*initialization->GetvoxelsizeY()*1.0;
	saveArea(ContourArea,Path_output, nameofproject);
	
	// Set parameters and coef
	last_parameters = parameters;
	last_parameters->SetNumcoefs(computation->Getnum_coefs());
	last_parameters->SetCoefs(computation->Getcoefs());
	
	if(withgui)
	{
		// Set the label
		QString str;
		Arealabel->setText(str.setNum(ContourArea));
	}
	
	if(debug)
		std::cout<<" --Preview/Init--"<<std::endl;
	return 0;
}


/**************************************************************************************
 * Main function : calcule the segmentation of the corpus collosum and draw it for each
 * iteration. Creat just an object, the CCsegtool_computation
 *************************************************************************************/
void compute( std::string Path_output, std::string nameofproject, int number_pts, bool Unconstrained, 
	      std::string WMintensity, std::string MPSDisplacement, std::string Number_iteration, int Lambdamax,
	      std::string Coefofoptim, bool debug, bool withgui, Imageview * &imageview, 
	      CCsegtool_parameters* &last_parameters, double *lastX, double *lastY, double *lastScale, double *lastRot,
	      std::string A, double r, bool rot90, bool rot180, bool rot270,
	      QWidget * parent, QLabel* instructionsview, QLabel* parametersinitview,QLabel* Arealabel)
{
	/* create local variables */
	int iter;
	Point2 *pts=new Point2[number_pts];
	for(int i=0;i<number_pts;i++)
		pts[i][0]=pts[i][1]=0;
	
	/* Parameter */
	CCsegtool_parameters *parameters;
	
	/* Use the last parameters */
	parameters = last_parameters;
	
	/*coef A*/
	double coefA=0.0;
	if(A.compare("default")==0)
		coefA = 100.0*CalculBestfitmean(Path_output, nameofproject);
	else if(A.compare("")==0)
		coefA = 0.0;
	else
		coefA = atof(A.c_str());
	
	// computation 
	CCsegtool_computation * computation = new CCsegtool_computation(parameters, nameofproject, Path_output, debug,
									coefA, r);
	
	//Set the repulsive points
	if(withgui)
		computation->setRepulPoints(imageview->GetRepulPoints());
	
	/* global variable */
	float disp = 100000000;
	/* initiate WMvalue, MPSdisplacement, number of iteration and Coefficient of optimization*/
	int WMvalue1,WMvalue2,WMvalue3;
	int MPSDisplacement1,MPSDisplacement2,MPSDisplacement3;
	int Number_iteration1,Number_iteration2,Number_iteration3;
	float coef1,coef2,coef3;
	
	/*WMintensity*/
	WMvalue1= (int)(parameters->GetWM()*1.3);
	WMvalue2= (int)(parameters->GetWM());
	WMvalue3= (int)(parameters->GetWM());
	if(WMintensity.compare("default,default,default")!=0)
	{
		WMvalue1=readstringWMI(WMintensity, 1, parameters->GetWM());
		WMvalue2=readstringWMI(WMintensity, 2, parameters->GetWM());
		WMvalue3=readstringWMI(WMintensity, 3, parameters->GetWM());
	}
	
	/*MPSDisplacement*/
	MPSDisplacement1= 10;
	MPSDisplacement2= 3;
	MPSDisplacement3= 2;
	if(MPSDisplacement.compare("10,3,2")!=0)
	{
		MPSDisplacement1=readstringinputi(MPSDisplacement, 1, 2);
		MPSDisplacement2=readstringinputi(MPSDisplacement, 2, 2);
		MPSDisplacement3=readstringinputi(MPSDisplacement, 3, 2);
	}
	
	/*Number of iteration*/
	Number_iteration1= 50;
	Number_iteration2= 5;
	Number_iteration3= 3;
	if(Number_iteration.compare("50,15,3")!=0)
	{
		Number_iteration1=readstringinputi(Number_iteration, 1, 3);
		Number_iteration2=readstringinputi(Number_iteration, 2, 3);
		Number_iteration3=readstringinputi(Number_iteration, 3, 3);
	}
	
	/*Coefficient of optimization*/
	coef1= 0.25* parameters->GetVoxelSizeX();
	coef2= 0.1* parameters->GetVoxelSizeX();
	coef3= 0.05* parameters->GetVoxelSizeX();
	if(Coefofoptim.compare("0.25,0.1,0.05")!=0)
	{
		coef1=readstringinputf(Coefofoptim, 1)* parameters->GetVoxelSizeX();
		coef2=readstringinputf(Coefofoptim, 2)* parameters->GetVoxelSizeX();
		coef3=readstringinputf(Coefofoptim, 3)* parameters->GetVoxelSizeX();
	}
	
	/* Check the value of WMvalue/MPSDisplacement parameters */
	if( WMvalue1>0  && MPSDisplacement1>0 && MPSDisplacement1<22)
	{
		if(Number_iteration1!=0)
		{
			/*First Step */
			if(withgui)
			{
				QString instruct="Running Step 1";
				instructionsview->setText( instruct );
			}
			
			iter=run_step(parameters, computation, imageview, parent, 1, WMvalue1, MPSDisplacement1,
				Number_iteration1, coef1, pts, number_pts, Lambdamax, disp, debug,withgui);
			
			/* Show the instructions "Converge iteration x" or "Not converge" */
			if(withgui) 
				showinstruction(instructionsview,1,iter,Unconstrained,Number_iteration1);
		}
		
		/* More Step */
		if(WMvalue2>0 && MPSDisplacement2>0 && MPSDisplacement2<22)
		{
			if(Number_iteration2!=0)
			{
				/* second step */
				if(withgui)
				{
					QString instruct = instructionsview->text();
					instruct.append("\n");
					instruct.append("Running Step 2");
					instructionsview->setText( instruct );
				}
				
				iter=run_step(parameters, computation, imageview, parent,2, WMvalue2, MPSDisplacement2,
					      Number_iteration2, coef2, pts, number_pts, Lambdamax, disp, debug, withgui);
				
				/* Show the instructions "Converge iteration x" or "Not converge" */
				if(withgui)
					showinstruction(instructionsview,2,iter,Unconstrained,Number_iteration2);
			}
			
			/*Unconstrained step*/
			if(Unconstrained && Number_iteration3!=0)
			{
				if(WMvalue3>0 && MPSDisplacement3>0 && MPSDisplacement3<22)
				{
					/* third step */
					if(withgui)
					{
						QString instruct = instructionsview->text();
						instruct.append("\n");
						instruct.append("Running Step 3");
						instructionsview->setText( instruct );
					}
					iter=run_step(parameters, computation, imageview, parent, 3, WMvalue3,
							MPSDisplacement3, Number_iteration3,coef3, pts, number_pts,
							Lambdamax, disp, debug, withgui);
					/* Show the instructions "Converge iteration x" or "Not converge" */
					if(withgui)
						showinstruction(instructionsview,3,iter,Unconstrained,
								Number_iteration3);
				}
				else
					std::cout<<"Error: wrong value for White Matter intensity or MPSDisplacement for second step."<<std::endl;
			}
		}
		else
			std::cout<<"Error: wrong value for White Matter intensity or MPSDisplacement for second step."<<std::endl;
		
		/* Save the output image */
		if(debug)
			std::cout<<" Save output image named finalImage.png\n"<<std::endl;
		/*creat the pixmap with the output image */
		QString imageFileName( (Path_output + "/" + nameofproject +  "_CCslice.png").c_str());
		QString file((Path_output + "/" + nameofproject +  "_finalImage.png").c_str());
		
		/* reconstruct the contour */
		computation->reconstructCurve(pts,number_pts);
		
		/* Save parameters and image */
		computation->getparam(lastX, lastY, lastScale, lastRot);
		
		/* Save image with or without gui */
		if(withgui)
		{
			(imageview->getPixmap())->save(file);
			/*write ouput file*/
			computation->writeoutput((imageview->Getsize()).width(), (imageview->Getsize()).height());
		}
		else
		{
			Imageview imageviewer(imageFileName,false,false,false,debug,parent);
			/* Draw the contour on the image */
			imageviewer.SaveImage(file, pts,number_pts, computation->GetcenterX(), computation->GetcenterY());
			/* Save the image */
			(imageviewer.getsavePixmap())->save(file);
			/*write ouput file*/
			computation->writeoutput((imageviewer.Getsize()).width(), (imageviewer.Getsize()).height());
		}
		
		/*write ouput file*/
		parameters->writeoutput(Path_output, nameofproject);
		computation->writedp(coef1, coef2,Unconstrained);
		
		/* Save parameters */
		computation->getparam(lastX, lastY, lastScale, lastRot);
	}
	else
		std::cout << "Error: wrong value for White Matter intensity or MPSDisplacement for first step."
				<< std::endl;
	
	/* save the parameters and computation */
	/* Create the polygon */
	double ContourArea;
	// the Array
	Point2D *ContourArray = new Point2D[number_pts];
	for(int i=0;i<number_pts;i++)
	{
		ContourArray[i].x = pts[i][0];
		ContourArray[i].y = pts[i][1];
	}
	// the polygon and add the number of points with the array
	Polygon2D Contour;
	Contour.numPoints = number_pts;
	Contour.pointArray = ContourArray;
	ContourArea = contourArea(Contour);
	
	// Multiply by the size of a voxel
	ContourArea *= parameters->GetVoxelSizeX()*parameters->GetVoxelSizeY()*1.0;
	saveArea(ContourArea,Path_output, nameofproject);
	
	// Set parameters and coef
	last_parameters = parameters;
	last_parameters->SetNumcoefs(computation->Getnum_coefs());
	last_parameters->SetCoefs(computation->Getcoefs());
	
	if(withgui)
	{
		// Set the label
		QString str;
		Arealabel->setText(str.setNum(ContourArea));
	}
	
	if(debug)
		std::cout<<" --Main done--"<<std::endl;
}


/********************************************************************************* 
 * Change WMI string in WMI int.
 ********************************************************************************/
int readstringWMI(std::string str, int position, int initWMvalue)
{
	std::string str2,buf1,buf2;
	int WMvalue = 0;
	size_t found;
	
	/* initialization */
	found=str.find_first_of(",");
	if(found!=std::string::npos)
	{
		switch(position)
		{
			case 1:
				str2= str.substr (0,found);
				if(str2.compare("default")==0)
					WMvalue = static_cast<int>(1.3*initWMvalue);
				else
					WMvalue = atoi(str2.c_str());
				break;
			case 2:
				str2 = str.substr(str.find_first_of(",")+1,
						str.find_last_of(",")-str.find_first_of(",")-1);
				if(str2.compare("default")==0)
					WMvalue = initWMvalue;
				else
					WMvalue = atoi(str2.c_str());
				break;
			case 3:
				if(str.find_last_of(",")==str.find_first_of(","))
					WMvalue = initWMvalue;
				else
				{
					str2 = str.substr (str.find_last_of(",")+1,str.size());
					if(str2.compare("default")==0)
						WMvalue = initWMvalue;
					else
						WMvalue = atoi(str2.c_str());
				}
				break;
		}
	}
	else
	{
		switch(position)
		{
			case 1:
				if(str.compare("default")==0)
					WMvalue = static_cast<int>(1.3*initWMvalue);
				else
					WMvalue = atoi(str.c_str());
				break;
			case 2:
				WMvalue = initWMvalue;
				break;
			case 3:
				WMvalue = initWMvalue;
				break;
		}
	}
	return WMvalue;
}


/********************************************************************************* 
 * Change the string in int. String like "100,20,3"
 ********************************************************************************/
int readstringinputi(std::string str, int position, int type)
{
	std::string str2,buf1,buf2;
	size_t found;
	
	/* initialization */
	found=str.find_first_of(",");
	if(found!=std::string::npos)
	{
		switch(position)
		{
			case 1:
				str2= str.substr (0,found);
				break;
			case 2:
				str2 = str.substr(str.find_first_of(",")+1,
						str.find_last_of(",")-str.find_first_of(",")-1);
				break;
			case 3:
				if(str.find_last_of(",")==str.find_first_of(","))
				{
					if(type==1)
						str2 = "90";
					if(type==2)
						str2 = "2";
					else
						str2 = "3";
				}
				else
					str2 = str.substr (str.find_last_of(",")+1,str.size());
				break;
		}
	}
	else
	{
		switch(position)
		{
			case 1:
				str2= str;
				break;
			case 2:
				if(type==1)
					str2 = "90";
				if(type==2)
					str2 = "3";
				else
					str2 = "5";
				break;
			case 3:
				if(type==1)
					str2 = "90";
				if(type==2)
					str2 = "2";
				else
					str2 = "3";
				break;
		}
	}
	return atoi(str2.c_str());
}


/********************************************************************************* 
 * Change the string in float. String like "0.5,0.01,0.006"
 ********************************************************************************/
float readstringinputf(std::string str, int position)
{
	std::string str2,buf1,buf2;
	size_t found;
	
	/* initialization */
	found=str.find_first_of(",");
	if(found!=std::string::npos)
	{
		switch(position)
		{
			case 1:
				str2= str.substr (0,found);
				break;
			case 2:
				str2 = str.substr(str.find_first_of(",")+1,
						str.find_last_of(",")-str.find_first_of(",")-1);
				break;
			case 3:
				if(str.find_last_of(",")==str.find_first_of(","))
				{
					str2 = "0.5";
				}
				else
					str2 = str.substr (str.find_last_of(",")+1,str.size());
				break;
		}
	}
	else
	{
		switch(position)
		{
			case 1:
				str2= str;
				break;
			case 2:
				str2 = "0.8";
				break;
			case 3:
				str2 = "0.5";
				break;
		}
	}
	return atof(str2.c_str());
}


/********************************************************************************* 
 * Show the init Parameters value on the QLabel
 ********************************************************************************/
QString readparam(std::string filename)
{
	std::string buf, param;
	param = "";
	std::ifstream file(filename.c_str() , ios::in);  // open in reading
	
	if(file)  // if open
	{
		for(int i=0; i<5;i++)
		{
			getline(file, buf);
			param+=buf + "\n" ;
		}
		file.close();
	}
	else cerr << "Error: The opening of the file failed" << endl;
	
	return param.c_str();
}


/********************************************************************************* 
 * Show the instructions on the QLabel
 ********************************************************************************/
void showinstruction(QLabel* view, int step, int iter, bool Unconstrained, int Number_iteration)
{
	QString instruct = view->text();
	
	if( iter < 1)
	{
		instruct.append("\n");
		instruct.append("No convergence");
	}
	else
	{
		std::stringstream siter;
		siter << Number_iteration-iter;
		instruct.append("\n");
		instruct.append("Converge iteration ");
		instruct.append((siter.str()).c_str());
	}
	
	view->setText( instruct );
}


/********************************************************************************* 
 * Run_step execute the computation for one step.
 ********************************************************************************/
int run_step(CCsegtool_parameters* parameters, CCsegtool_computation* computation, Imageview* imageview, QWidget* parent,
	     int step, int WMvalue, int MPSDisplacement, int Number_iteration, float coef, Point2 *pts, int number_pts,
	     int Lambdamax, float disp, bool debug, bool withgui, bool preview)
{
	/* local variables */
	int iter;
	
	/* execution */
	parameters->SetWMIntensity(WMvalue);
	parameters->SetMPSDisplacement(MPSDisplacement);
	parameters->SetOiter(Number_iteration);
	iter = Number_iteration;
	if(step==3)
		parameters->SetUnconstrained(true);
	else
		parameters->SetUnconstrained(false);
	while ( (iter > 0) && (disp > coef) ) {
		if(debug)
			std::cout<<" Step : "<<step<<  ", number of iteration : "<<iter<<std::endl;
		disp = computation->execution(parameters, step, Lambdamax, iter);
		if(withgui)
		{
			/* Paint the contour in red for each iteration*/
			if(debug)
				std::cout<<" Paint the contour on input image"<<std::endl;
			computation->reconstructCurve(pts,number_pts);
			imageview->GetPts(pts,number_pts, computation->GetcenterX(), computation->GetcenterY());
			parent->repaint();
			parent->update();
		}
		if(preview)
			return 0;
		iter--;
	}
	return iter;
}

/********************************************************************************* 
 * Updatedrawing draws the segmentation with the new parameters
 ********************************************************************************/
void updatedrawing(double X,double Y,double Scale,double Rotation, QWidget* parent, CCsegtool_parameters *parameters,
		   Imageview* &imageview, int number_pts, std::string Path_output,  std::string nameofproject, 
		   int Lambdamax, bool debug)
{
	/* local variable */
	Point2 *pts=new Point2[number_pts];
	for(int i=0;i<number_pts;i++)
		pts[i][0]=pts[i][1]=0;
	/* Set the param in the CCsegtoolparameters */
	parameters->setparam(X,Y,Scale,Rotation);
	/* Create the computation */
	CCsegtool_computation computation(parameters, nameofproject, Path_output);
	float disp = computation.execution(parameters, 0, Lambdamax, 1,true);
	if(disp!=0)
		std::cout<<"Error in update drawing"<<std::endl;
	computation.reconstructCurve(pts,number_pts);
	imageview->GetPts(pts,number_pts, computation.GetcenterX(), computation.GetcenterY());
	parent->repaint();
	parent->update();
}

/********************************************************************************* 
 * Calcul the bestfitmean if A is default
 ********************************************************************************/
double CalculBestfitmean(std::string Path_output,  std::string nameofproject)
{
	std::string filename = Path_output + "/" + nameofproject + "_goodnessoffit.gof";
	std::ifstream file(filename.c_str() , ios::in);  // open in reading the file from Run
	std::string str;
	double sum = 0.0;
	double i=0.0;
	
	if(file)  // if open
	{
		//the first line
		getline(file, str);
		
		while(!file.eof())
		{
			getline(file, str);
			sum += atof(str.c_str());
			i=i+1.0;
		}
		file.close();
	}
	else cerr << "ERROR: No godnessoffit.gof file found" << endl;
	
	return sum/i;
}

/********************************************************************************* 
 * Calcul the area of the contour
 ********************************************************************************/
double contourArea(Polygon2D poly)
{
	double area = 0.0;
	double xi,xiplus1,yi,yiplus1;
	int i;
 
  	/* Use the sum of trapezoids to obtain the area of a polygon. Area = sigma((y(i)+y(i+1))*(x(i)-x(i+1))/2)  */
	for(i=0;i<poly.numPoints;i++){
		xi = poly.pointArray[i].x;
		yi = poly.pointArray[i].y;
		xiplus1 = poly.pointArray[(i+1)%poly.numPoints].x;
		yiplus1 = poly.pointArray[(i+1)%poly.numPoints].y;
		area = area + (yi+yiplus1)*(xi-xiplus1)/2.0;
	}
	return area;
}

/********************************************************************************* 
 * Save the area of the contour
 ********************************************************************************/
void saveArea(double area, std::string Path_output, std::string nameofproject)
{
	std::string outfile;
	outfile = Path_output + "/" + nameofproject + "_area.txt";
	ofstream bfile(outfile.c_str(), ios::out);
	bfile << "Area of the contour :" <<endl;
	bfile << area << endl;
	bfile.close();
}
