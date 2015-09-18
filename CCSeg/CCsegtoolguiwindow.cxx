#include "CCsegtoolguiwindow.h"

/********************************************************************************* 
 * Constructor
 ********************************************************************************/
CCsegtoolguiwindow::CCsegtoolguiwindow(float scalefactor, bool Debug, QWidget * parent , Qt::WFlags f  ): QMainWindow(parent, f)
{
	setupUi(this);
	
	this->setMouseTracking(true);
	
	/* set the main layout */
	/* ScrollArea */
	QScrollArea *const scroll(new QScrollArea);
	centralwidget->setLayout(GLayoutPrincipal);
	scroll->setWidget(centralwidget);
	setCentralWidget(scroll);
	
	/* set the Debug */
	m_Debug = Debug;
	if(Debug)
		debug->setCheckState(Qt::Checked);
	
	/* MessageBox */
	m_information = new QMessageBox(this);
	m_information->setWindowTitle("Help");
	/* Pixmap for the output image */
	m_outputimage= new QPixmap();
	/* Saveimage disable */
	actionSaveimage->setEnabled(false);
	
	/* Change display : size of texts */
	title->setFont(QFont("Courrier", 20,4));
	Paramtitle->setFont(QFont("Courrier", 14,4));
	pathforinputitle->setFont(QFont("Courrier", 14,4));
	CCinittile->setFont(QFont("Courrier", 14,4));
	instructionstitle->setFont(QFont("Courrier", 12,4));
	Initparamframetitle->setFont(QFont("Courrier", 12,4));
	pixeltitle->setFont(QFont("Courrier", 12,4));
	GOFmeantitle->setFont(QFont("Courrier", 12,4));
	Areatitle->setFont(QFont("Courrier", 12,4));
	IsolineLabel->setFont(QFont("Courrier", 12,4));
	visuflip->setFont(QFont("Courrier", 11,4));
	
	/* Initialization of the parameters */
	m_OldInputImage = "";
	m_old_permute_x_y=false;
	m_oldreflectX = false;
	m_oldreflectY = false;
	m_Scrollexist = false;
	m_scaleFactor = 1.0;
	m_oldfilename = "";
	// Parameters 
	m_last_parameters = NULL;
	// Imageview
	m_viewimage = NULL;
	//preview by default
	m_scalefactor = scalefactor;
	m_runcalled = false;
	
	/*connection of slots*/
	//Run
	connect(this->Run, SIGNAL(clicked()), this, SLOT(running()));
	
	//Continue
	connect(this->initialization, SIGNAL(clicked()), this, SLOT(InitPreview()));
	
	//Compute probability model
	connect(this->ComputeProbaModel, SIGNAL(clicked()), this, SLOT(ComputeProbabilityModel()));
	
	//help
	connect(this->helpbutton, SIGNAL(clicked()), this, SLOT(help()));
	
	//default Tab Run
	connect(this->defaultparam, SIGNAL(clicked()), this, SLOT(defaultparaminput()));
	
	//default Tab Parameters for computation
	connect(this->defaultcomputation, SIGNAL(clicked()), this, SLOT(defaultparamcomputation()));
	
	//browsers for input data
	connect(this->browImage, SIGNAL(clicked()), this, SLOT(browserImage()));
	connect(this->browmask, SIGNAL(clicked()), this, SLOT(browserMask()));
	connect(this->browCCAtlasDirectory, SIGNAL(clicked()), this, SLOT(browserCCAtlasDirectory()));
	connect(this->browpathoutput, SIGNAL(clicked()), this, SLOT(browserOutput()));
	CCAtlasDirectorybox->setText(getenv("CCSEG_ATLASDIR")?getenv("CCSEG_ATLASDIR"):"");
	
	//"File" action
	connect(this->actionNew, SIGNAL(triggered()), SLOT(newparam()));
	connect(this->actionOpen, SIGNAL(triggered()), SLOT(openparam()));
	connect(this->actionSave, SIGNAL(triggered()), SLOT(saveparam()));
	connect(this->actionSaveimage, SIGNAL(triggered()), SLOT(saveimage()));
	
	//Debug
	connect(this->debug, SIGNAL(stateChanged(int)), this, SLOT(setDebug()));
	
	//Preview image drawing connection bouton/signal/slot
	connect(this->permute_x_y, SIGNAL(stateChanged(int)), this, SLOT(InitPreview()));
	connect(this->ReflectXOn, SIGNAL(stateChanged(int)), this, SLOT(InitPreview()));
	connect(this->ReflectYOn, SIGNAL(stateChanged(int)), this, SLOT(InitPreview()));
	connect(this->DoubleOn, SIGNAL(stateChanged(int)), this, SLOT(ConPreview()));
	connect(this->interpolationlinear, SIGNAL(stateChanged(int)), this, SLOT(ConPreview()));
	connect(this->imagemask, SIGNAL(editingFinished()), this, SLOT(InitPreview()));
	connect(this->Midsaggitalplaneslice, SIGNAL(editingFinished()), this, SLOT(ConPreview()));
	connect(this->averagenum, SIGNAL(editingFinished()), this, SLOT(ConPreview()));
	connect(this->SliceDir, SIGNAL(editingFinished()), this, SLOT(InitPreview()));
	connect(this->othercomponant, SIGNAL(stateChanged(int)), this, SLOT(ConPreview()));
	connect(this->ContourRot0, SIGNAL(clicked()), this, SLOT(ConPreview()));
	connect(this->ContourRot90, SIGNAL(clicked()), this, SLOT(ConPreview()));
	connect(this->ContourRot180, SIGNAL(clicked()), this, SLOT(ConPreview()));
	connect(this->ContourRot270, SIGNAL(clicked()), this, SLOT(ConPreview()));

	//Update the drawing
	connect(this->updatebutton, SIGNAL(clicked()), this, SLOT(update()));
	connect(this->VisuRot0, SIGNAL(clicked()), this, SLOT(updateVisuImage()));
	connect(this->VisuRot90, SIGNAL(clicked()), this, SLOT(updateVisuImage()));
	connect(this->VisuRot180, SIGNAL(clicked()), this, SLOT(updateVisuImage()));
	connect(this->VisuRot270, SIGNAL(clicked()), this, SLOT(updateVisuImage()));
	
	//clear pts button
	connect(this->clearpts, SIGNAL(clicked()), this, SLOT(resetpts()));
	
	//isoline value
	connect(this->isoline, SIGNAL(valueChanged(int)), this, SLOT(changeisoline()));
	
	//coef r and A
	connect(this->A, SIGNAL(textChanged(const QString)), this, SLOT(changecoef()));
	connect(this->r, SIGNAL(textChanged(const QString)), this, SLOT(changecoef()));
}


/*********************************************************************************
 * Two different preview
 ********************************************************************************/
void CCsegtoolguiwindow::InitPreview()
{
	/* reset the Visualization rotation */
	Preview(0);
}

void CCsegtoolguiwindow::ConPreview()
{
	/* Keep the Visualization rotation */
	Preview(1);
}

/*********************************************************************************
 * Show Preview Image slot
 ********************************************************************************/
void CCsegtoolguiwindow::Preview(int resetVisu)
{
	//Create the ScrollArea
	if(!m_Scrollexist)
	{
		/* ScrollArea */
		m_scrollArea = new QScrollArea(this);
		m_scrollArea->setBackgroundRole(QPalette::Base);
	}
	
	/* Mouse Event when running */
	QApplication::setOverrideCursor(QCursor(Qt::WaitCursor));
	
	if(m_Debug)
		std::cout<<"-----PREVIEW-----"<<std::endl;
	/* Call compute */
	if(inputImage->text().compare("")==0)
		QMessageBox::information(this, "Warning", "Path for the image no set");
	else if  (imagemask->text().compare("")==0)
		QMessageBox::information(this, "Warning", "Path for the tissue segmentation mask no set");
	else if  (CCAtlasDirectorybox->text().compare("")==0)
		QMessageBox::information(this, "Warning", "Path for the vector file no set");
	else
	{
		/* Set parameters */
		setvalue();
		if(m_viewimage!=NULL)
		{
			if(resetVisu==0)
			{
				/* reset the rotation of the visualization image */
				//VisuRot0
				VisuRot0->setChecked ( true );
				/* reset the rotation boolean value at false */
				m_viewimage->setVisuFlip(false, false, false);
			}
			/* reset point */
			resetpts();
			normalSize();
		}
		else
		{
			/* reset the rotation of the visualization image */
			//VisuRot0
			VisuRot0->setChecked ( true );
		}
		
		/* CCsegtool execution */
		if(m_InputImage.compare("")==0 || m_MaskImage.compare("")==0 || m_CCAtlasDirectory.compare("")==0)
			std::cout<<" Error: Input Image or Binary Mask or CCAtlas directory path not specified."<<std::endl;
		/* Check the value of the parameters SliceDir, PSDistance */
		else if(m_PSDistance<0.0001 || m_PSDistance>1.0 || m_slicedir<0 || m_slicedir>2 || m_averagenum>100)
			std::cout<<" Error: PSDistance or Slice Direction or average number has a wrong value."<<std::endl;
		/* Initiate CCsegtool*/
		else
		{
			instructionsview->clear ();
			if(initialize( m_InputImage, m_MaskImage, m_CCAtlasDirectory, m_outputfolder, m_nameofproject, 
			   m_interpolationlinear, m_vesselremove, m_seglabel, m_averagenum, m_permute_x_y, m_reflectX, 
			   m_reflectY, m_opening, m_double, m_slicedir, m_MidPlaneSliceNumber, m_number_pts, m_FSXForm,
			   m_PSDistance, m_unconstrained, m_WMvalue, m_MPSDisplacement, m_loop, m_Lambdamax, m_Coefofoptim,
			   m_Debug, true, m_viewimage, m_last_parameters, &m_lastX, &m_lastY, &m_lastScale, &m_lastRot,
			   m_scalefactor, m_othercomponant, m_angle, m_rot90, m_rot180, m_rot270, 
			   this, this->parametersinitview, this->Arealabel, this->Midsaggitalplaneslicevalue)!=0)
				QMessageBox::information(this, "Warning", "Problem with the images or with the Slice direction. Total Mass of the image was zero. Aborting here to prevent division by zero later on.");
			else
			{
				/* Init the scroll */
				if(!m_Scrollexist)
				{
					/* Show the scroolArea */
					m_scrollArea->setFixedSize((m_viewimage->Getsize()).width(),
							(m_viewimage->Getsize()).height());
					m_scrollArea->setWidgetResizable(true);
					/* Set the widget for the scrollArea */
					if((m_viewimage->Getsize()).width()>(m_viewimage->Getsize()).height())
					{
						m_scrollArea->setFixedSize((m_viewimage->Getsize()).width(),
								(m_viewimage->Getsize()).height());
						/* Set on the Imageview */
						m_scrollArea->setWidget(m_viewimage);
						m_scrollArea->move(40,145);
					}
					else
					{
						m_scrollArea->setFixedSize((m_viewimage->Getsize()).width(),
								(m_viewimage->Getsize()).height());
						/* Set on the Imageview */
						m_scrollArea->setWidget(m_viewimage);
						m_scrollArea->move(110,85);
					}
					GLayoutPrincipal->addWidget(m_scrollArea,2,0);
					m_Scrollexist = true;
				}
				else if(resetVisu==0)
				{
					if((m_viewimage->Getsize()).width()>(m_viewimage->Getsize()).height())
					{
						m_scrollArea->setFixedSize((m_viewimage->Getsize()).width(),
								(m_viewimage->Getsize()).height());
						/* Set on the Imageview */
						m_scrollArea->setWidget(m_viewimage);
						m_scrollArea->move(40,145);
					}
					else
					{
						m_scrollArea->setFixedSize((m_viewimage->Getsize()).width(),
								(m_viewimage->Getsize()).height());
						/* Set on the Imageview */
						m_scrollArea->setWidget(m_viewimage);
						m_scrollArea->move(110,85);
					}
					
					if(m_Debug)
					{
						std::cout<<"QscrollArea size : "<<m_scrollArea->width()<<
								" "<<m_scrollArea->height()<<std::endl;
						std::cout<<"QLabel size : "<<m_viewimage->width()<<
								"  "<<m_viewimage->height()<<std::endl;
					}
					/* Install the filter */
					m_viewimage->installEventFilter(this);
				}
				
				/* Show parameters */
				if(m_Debug)
					std::cout<<"Show Parameters"<<std::endl;
				/* Show the value of the four variables */
				showparameters(m_lastX, m_lastY, m_lastScale, m_lastRot);
			}
		}
	}
	
	/* Enable saveimage */
	actionSaveimage->setEnabled(true);
	/* Restore the mouse */
	QApplication::restoreOverrideCursor();
}


/*********************************************************************************
 * When the button Run is clicked
 ********************************************************************************/
void CCsegtoolguiwindow::running()
{
	/* Mouse Event when running */
	QApplication::setOverrideCursor(QCursor(Qt::WaitCursor));
	if(m_viewimage!=NULL)
	{
		if(m_Debug)
			std::cout<<"-----RUN-----"<<std::endl;
		/* Call compute */
		if(inputImage->text().compare("")==0)
			QMessageBox::information(this, "Warning", "Path for the image no set");
		else if  (imagemask->text().compare("")==0)
			QMessageBox::information(this, "Warning", "Path for the tissue segmentation mask no set");
		else if  (CCAtlasDirectorybox->text().compare("")==0)
			QMessageBox::information(this, "Warning", "Path for the vector file no set");
		else
		{
			/* Set parameters */
			setvalue();
			if(m_viewimage!=NULL)
				normalSize();
	
			/* Use X/Y/Scale/Rotation from the last segmentation (initialization or run)*/
			m_last_parameters->setparam(m_lastX, m_lastY, m_lastScale, m_lastRot);
	
			if(m_PSDistance<0.0001 || m_PSDistance>1.0 || m_slicedir<0 || m_slicedir>2 || m_averagenum>100)
				std::cout<<" Error: PSDistance or Slice Direction has a wrong value."<<std::endl;
			/* Compute CCsegtool*/
			else
			{
				instructionsview->clear ();
				compute(m_outputfolder, m_nameofproject, m_number_pts, m_unconstrained, m_WMvalue,
					m_MPSDisplacement, m_loop, m_Lambdamax, m_Coefofoptim, m_Debug,true,
					m_viewimage, m_last_parameters, &m_lastX, &m_lastY, &m_lastScale,
					&m_lastRot, m_A, m_r, m_rot90, m_rot180, m_rot270,
					this, this->instructionsview, this->parametersinitview,this->Arealabel);
				m_runcalled = true;
			}
			
			/* Show parameters */
			if(m_Debug)
				std::cout<<"Show Parameters"<<std::endl;
			/* Show the value of the four variables */
			showparameters(m_lastX, m_lastY, m_lastScale, m_lastRot);
			/* Enable saveimage */
			actionSaveimage->setEnabled(true);
			/* set the isoline */
			setisoline();
		}
	}
	else
	{
		QMessageBox::information(this, "Warning", "Call Preview first.");
	}
	/* Restore the mouse */
	QApplication::restoreOverrideCursor();
}


/*********************************************************************************
 * When the button Compute Probability Model is clicked
 ********************************************************************************/
void CCsegtoolguiwindow::ComputeProbabilityModel()
{
	/* Mouse Event when running */
	QApplication::setOverrideCursor(QCursor(Qt::WaitCursor));
	if(m_Debug)
		std::cout<<"-----COMPUTE PROBABILITY MODEL-----"<<std::endl;
	/* Call compute */
	if(m_runcalled)
	{
		// Create the Probacurve and set the parameters
		CCcurveProba *Probacurve;
		Probacurve = new CCcurveProba(m_last_parameters->GetImage3D(), m_last_parameters->GetPointsvtk(),
					      m_CCAtlasDirectory, m_outputfolder, m_nameofproject, m_vesselremove,
					      m_seglabel, m_averagenum, m_permute_x_y, m_reflectX, m_reflectY,
					      m_opening, m_double, m_slicedir, m_last_parameters->GetNumSlice(),
					      m_last_parameters->GetNumPointsProba(),m_Debug);
		
		/* Compute */
		Probacurve->compute_proba();
	}
	else
		QMessageBox::information(this, "Warning", "Call Run before Probability model");
	/* Restore the mouse */
	QApplication::restoreOverrideCursor();
}


/*********************************************************************************
 * Set the value of the parameters
 ********************************************************************************/
void CCsegtoolguiwindow::setvalue()
{
	/* First TAB*/
	//MPSDisplacement
	m_MPSDisplacement = MPSDisplacement->text().toStdString();
	//Numberloop
	m_loop = Numberloop->text().toStdString();
	//WMintensity
	m_WMvalue = WMintensity->text().toStdString();
	//Lambdamax
	m_Lambdamax = Lambdamax->text().toInt();
	//m_Coefofoptim
	m_Coefofoptim = Coefofoptim->text().toStdString();
	//A
	m_A = A->text().toStdString();
	//r
	m_r = r->text().toDouble();
	//X
	m_lastX = X->text().toDouble();
	//Y
	m_lastY = Y->text().toDouble();
	//Rot
	m_lastRot = Rotation->text().toDouble();
	//Scale
	m_lastScale = Scale->text().toDouble();
	//visualisation rotation of 0 degres
	if(VisuRot0->isChecked())
		m_noRotation=true;
	else
		m_noRotation=false;
	//visualisation rotation of 90 degres
	if(VisuRot90->isChecked())
		m_rot90=true;
	else
		m_rot90=false;
	//visualisation rotation of 180 degres
	if(VisuRot180->isChecked())
		m_rot180=true;
	else
		m_rot180=false;
	//visualisation rotation of 180 degres
	if(VisuRot270->isChecked())
		m_rot270=true;
	else
		m_rot270=false;
	
	/* Second TAB*/
	/* Protection */
	//inputImage
	m_InputImage = inputImage->text().toStdString();
	//CCAtlasDirectorybox
	m_CCAtlasDirectory = CCAtlasDirectorybox->text().toStdString();
	//imagemask
	m_MaskImage = imagemask->text().toStdString();
	//pathoutput
	if(pathoutput->text().isEmpty())
		m_outputfolder= m_InputImage.substr(0, m_InputImage.find_last_of("/"));
	else
		m_outputfolder = pathoutput->text().toStdString();
	//Name of project
	if(nameofproject->text().isEmpty())
		m_nameofproject = m_InputImage.substr(m_InputImage.find_last_of("/")+1,
			m_InputImage.find_last_of(".")-m_InputImage.find_last_of("/")-1);
	else
		m_nameofproject = nameofproject->text().toStdString();
	
	/* Third TAB*/
	//Vesselremove
	if(VesselRemoveOn->isChecked())
		m_vesselremove=true;
	else
		m_vesselremove=false;
	//Rotation90
	if(permute_x_y->isChecked())
		m_permute_x_y=true;
	else
		m_permute_x_y=false;
	//ReflectX
	if(ReflectXOn->isChecked())
		m_reflectX=true;
	else
		m_reflectX=false;
	//ReflectYOn
	if(ReflectYOn->isChecked())
		m_reflectY=true;
	else
		m_reflectY=false;
	//DoubleOn
	if(DoubleOn->isChecked())
		m_double=true;
	else
		m_double=false;
	//OpeningOn
	if(OpeningOn->isChecked())
		m_opening=true;
	else
		m_opening=false;
	//UnConstrained
	if(UnConstrained->isChecked())
		m_unconstrained=true;
	else
		m_unconstrained=false;
	//SegLabel
	if(SegLabel->isChecked())
		m_seglabel=true;
	else
		m_seglabel=false;
	//FSXform
	if(FSXform->isChecked())
		m_FSXForm=true;
	else
		m_FSXForm=false;
	//interpolationlinear
	if(interpolationlinear->isChecked())
		m_interpolationlinear=true;
	else
		m_interpolationlinear=false;
	//other componant
	if(othercomponant->isChecked())
		m_othercomponant=true;
	else
		m_othercomponant=false;
	//Slice direction
	m_slicedir = SliceDir->text().toInt() ;
	//MidPlaneSliceNumber
	m_MidPlaneSliceNumber = Midsaggitalplaneslice->text().toStdString();
	//averagenum
	m_averagenum = averagenum->text().toInt() ;
	//PSDistance
	m_PSDistance = PSDistance->text().toFloat();
	//Number of points for drawing
	m_number_pts = number_pts->text().toInt();
	//Sup Angle
	if(ContourRot0->isChecked())
		m_angle = 0;
	if(ContourRot90->isChecked())
		m_angle = 1;
	if(ContourRot180->isChecked())
		m_angle = 2;
	if(ContourRot270->isChecked())
		m_angle = 3;
}


/*********************************************************************************
 * Set the default value of the parameters in the first Tab
 ********************************************************************************/
void CCsegtoolguiwindow::setDebug()
{
	if(debug->isChecked())
	{
		m_Debug=true;
		std::cout<<"Debug activated"<<std::endl;
	}
	else
	{
		m_Debug=false;
		std::cout<<"Debug desactivated"<<std::endl;
	}
}


/********************************************************************************* 
 * Set the default value of the parameters in the first Tab
 ********************************************************************************/
void CCsegtoolguiwindow::defaultparaminput()
{
	//MPSDisplacement
	m_MPSDisplacement = "10,3,2";
	MPSDisplacement->setText(m_MPSDisplacement.c_str());
	//Numberloop
	m_loop = "50,5,3";
	Numberloop->setText(m_loop.c_str());
	//WMintensity
	m_WMvalue = "default,default,default";
	WMintensity->setText(m_WMvalue.c_str());
	//Lambdamax
	m_Lambdamax = 100;
	QString qstr;
	Lambdamax->setText(qstr.setNum(m_Lambdamax));
	//Coefofoptim
	m_Coefofoptim = "0.25,0.1,0.05";
	Coefofoptim->setText(m_Coefofoptim.c_str());
	//A
	m_A="default";
	A->setText(m_A.c_str());
	//r
	m_r=0.1;
	r->setText("0.1");
}


/********************************************************************************* 
 * Set the default value of the parameters in the third Tab
 ********************************************************************************/
void CCsegtoolguiwindow::defaultparamcomputation()
{
	//Vesselremove
	VesselRemoveOn->setCheckState(Qt::Unchecked);
	//Rotation90
	permute_x_y->setCheckState(Qt::Unchecked);
	//ReflectX
	ReflectXOn->setCheckState(Qt::Unchecked);
	//ReflectYOn
	ReflectYOn->setCheckState(Qt::Unchecked);
	//DoubleOn
	DoubleOn->setCheckState(Qt::Checked);
	//linear interpolation
	interpolationlinear->setCheckState(Qt::Unchecked);
	//OpeningOn
	OpeningOn->setCheckState(Qt::Unchecked);
	//UnConstrained
	UnConstrained->setCheckState(Qt::Unchecked);
	//SegLabel
	SegLabel->setCheckState(Qt::Unchecked);
	//FSXform
	FSXform->setCheckState(Qt::Unchecked);
	//other componant
	othercomponant->setCheckState(Qt::Unchecked);
	//Slice direction
	QString qstr;
	m_slicedir = 0;
	SliceDir->setText(qstr.setNum(m_slicedir));
	//MidPlaneSliceNumber
	m_MidPlaneSliceNumber = "default";
	Midsaggitalplaneslice->setText(m_MidPlaneSliceNumber.c_str());
	//averagenum
	m_averagenum = 2;
	averagenum->setText(qstr.setNum(m_averagenum));
	//PSDistance
	m_PSDistance = 0.5;
	PSDistance->setText(qstr.setNum(m_PSDistance));
	//Number of points for drawing
	m_number_pts = 256;
	number_pts->setText(qstr.setNum(m_number_pts));
	//Sup Angle
	ContourRot0->setChecked ( true );
}


/********************************************************************************* 
 * Show the instructions on the QLabel
 ********************************************************************************/
void CCsegtoolguiwindow::showparameters(double lastX, double lastY, double lastScale, double lastRot)
{
	QString qstr1;
	X->setText(qstr1.setNum(lastX));
	Y->setText(qstr1.setNum(lastY));
	Scale->setText(qstr1.setNum(lastScale));
	Rotation->setText(qstr1.setNum(lastRot));
}


/********************************************************************************* 
 * Open a dialog to search the name of the image file and write it in the GUI
 * Widget
 ********************************************************************************/
void CCsegtoolguiwindow::browserImage()
{
	QString filename,type;
	std::string image;
	filename = QFileDialog::getOpenFileName(this, "Open File input image", "/", 
			"Images (*.gipl *.gipl.gz *.mhd *.mha *.img *.hdr *.nhdr *.nrrd *.nii *.nii.gz)",&type);
	
	if(m_Debug)
		std::cout<<"Filename : "<< (filename.toStdString()).c_str() <<std::endl;
	
	/* Keep the value if cancel */
	if(filename !=NULL)
	{
		inputImage->setText(filename);
		if(m_Debug)
			std::cout<<"Path saved"<<std::endl;
	}
	else if(m_Debug)
		std::cout<<"Error : filename is NULL."<<std::endl;
}

/********************************************************************************* 
 * Open a dialog to search the name of the mask file and write it in the GUI
 * Widget
 ********************************************************************************/
void CCsegtoolguiwindow::browserMask()
{
	QString filename,type;
	std::string image;
	image = inputImage->text().toStdString();
	if(image.compare("")!=0)
	{
		QString directoryPath = (image.substr(0, image.find_last_of("/"))).c_str();
		if(m_Debug)
			std::cout<<"directory path : "<<(directoryPath.toStdString()).c_str()<<std::endl;
		filename = QFileDialog::getOpenFileName(this, "Open File input image", directoryPath, 
				"Images (*.gipl *.gipl.gz *.mhd *.mha *.img *.hdr *.nhdr *.nrrd *.nii *.nii.gz)",&type);
	}
	else
	{
		filename = QFileDialog::getOpenFileName(this, "Open File binary mask","/",
			"Images (*.gipl *.gipl.gz *.mhd *.mha *.img *.hdr *.nhdr *.nrrd *.nii *.nii.gz)",&type);
	}
	
	if(m_Debug)
		std::cout<<"Filename : "<< (filename.toStdString()).c_str() <<std::endl;
	
	/* Keep the value if cancel */
	if(filename !=NULL)
	{
		m_oldfilename = filename;
		imagemask->setText(filename);
		// Clear the display
		Arealabel->clear();
		GOFmean->clear();
		pixellabel->clear();
		X->setText("");
		Y->setText("");
		Rotation->setText("");
		Scale->setText("");
		ConPreview();
		if(m_Debug)
			std::cout<<"Path saved"<<std::endl;
	}
	else if(m_Debug)
		std::cout<<"Error : filename is NULL."<<std::endl;
}

/********************************************************************************* 
 * Open a dialog to search the name of the vector's path and write it in the GUI
 * Widget
 ********************************************************************************/
void CCsegtoolguiwindow::browserCCAtlasDirectory()
{
	QString path = QFileDialog::getExistingDirectory(this, tr("Select Atlas Directory"), "/tmp");
	if(path!=NULL)
		CCAtlasDirectorybox->setText(path);
}

/********************************************************************************* 
 * Open a dialog to search the name of the file for the output and write it in the GUI
 * Widget. It's optionnal, Input filename if there is no precision.
 ********************************************************************************/
void CCsegtoolguiwindow::browserOutput()
{
	QString path = QFileDialog::getExistingDirectory(this);
	if(path!=NULL)
		pathoutput->setText(path);
}

/********************************************************************************* 
 * Dialog Help
 ********************************************************************************/
void CCsegtoolguiwindow::help()
{
	QString help;
	help =  "CCsegtool is a software which does the segmentation of the Corpus Collosum. \n\n";
	help += "The sofware can't start without four parameters, which are : \n";
	help += "  - 3D inputfile in ITK readable format (Analyze, gipl, mha, mhd)\n";
	help += "  - a binary mask for automatic computation of the initial parameters\n";
	help += "  - the path where there are the files to initialize input vectors\n";
	help += "  - the output path where the output files will be save (optional, if you don't tell it, it will be in the input image file).\n\n";
	help += "Furthermore, you can configurate the initialization/computation parameters.\n";
	help += "Here are each options :\n";
	help += "  - USEFUL PARAMETERS :\n";
	help += "  - MPSDisplacement : Max Per-Step Displacement, 3 values for the each step. DEFAULT=10,3,2\n";
	help += "  - Number of loop : Number of iteration, 3 values for each step. DEFAULT=50,5,3\n";
	help += "  - White Mater intensity : 3 values for each step. DEFAULT=set by initialization\n";
	help += "  - Lambda max : factor for allowing deformations within PCA shape space, DEFAULT=100\n";
	help += "  - Phase stopping coefficient : average distance across curve per voxel, 3 values for each step. DEFAULT=0.25,0.1,0\n.";
	help += "  - A and r are used for the repulsive points when you call continue.\n";
	help += "  - X,Y,Rotation and Scale are parameters used when you call 'continue'. You can set this value.\n\n";
	help += "  - PARAMETERS FOR COMPUTATION :\n";
	help += "  - Vessel Remove : Remove vessels as initial step, DEFAULT=off\n";
	help += "  - Permutation X/Y : Slice Image is rotated with 90 degrees, DEFAULT=off\n";
	help += "  - Flip X : Slice Image is reflected at y-axis (flip horizontal), DEFAULT=off\n";
	help += "  - Flip Y : Slice Image is reflected at x-axis (flip vertical),  DEFAULT=off\n";
	help += "  - Opening : Apply open operation instead of closing in preprocessing of segfile, DEFAULT=off\n";
	help += "  - use the other component during the preprocessing, DEFAULT=off\n";
	help += "  - Increase the size by two, DEFAULT=on\n";
	help += "  - Linear interpolation when you double the size of the image, DEFAULT=off\n";
	help += "  - Unconstrained : add a unconstrained segmentation step at the end of the process, DEFAULT=off\n";
	help += "  - Segmentation Label white matter label in the segmentation file is assumed to be '1', otherwise use this option, DEFAULT=off\n";
	help += "  - Fix Similarity Form, DEFAULT=off\n";
	help += "  - debug, Default=off\n";
	help += "  - PSDistance : Profile Sampling Distance, DEFAULT=0.5\n";
	help += "  - Slice Direction [0|1|2] slicing direction (0=x,1=y,2=z), DEFAULT = 0\n";
	help += "  - Midsagittal Plane Slice is the number of slice in the direction of SliceDir(0 to max), DEFAULT = middle\n";
	help += "  - Average number : number of slices for averaging the gray values, if 0, no averaging is done, DEFAULT=2\n";
	help += "  - Number of points for the drawing, Default=256\n";
	help += "  - rotate the contour by adding an angle, Default=0 degrees\n";
	help += "\nInfo : Left mouse's button to have the pixel value, right mouse's button for repulsive point and middle button for zoom\n";
	help += "Warning : When there is three values separated with a virgule, you can choose to give just one/two or three values.\n";
	
	/* Set the Text */
	m_information->setIcon(QMessageBox::Information);
	m_information->setInformativeText(help);
	/* Open the message */
	m_information->update();
	m_information->exec();
}

/********************************************************************************* 
 * Open a new window with default parameters
 ********************************************************************************/
void CCsegtoolguiwindow::newparam()
{
	if(actionSaveimage->isEnabled())
	{
		m_runcalled = false;
		actionSaveimage->setEnabled(false);
		parametersinitview->clear();
		instructionsview->clear();
		Arealabel->clear();
		Midsaggitalplaneslicevalue->clear();
		m_viewimage->hide();
		m_viewimage = NULL;
		defaultparaminput();
		defaultparamcomputation();
		pixellabel->clear();
		m_scrollArea = NULL;
		m_Scrollexist = false;
		this->repaint();
	}
	else
		std::cout<<"WARNING: Already new window"<<std::endl;
}

/********************************************************************************* 
 * Save the value of the parameters
 ********************************************************************************/
void CCsegtoolguiwindow::saveparam()
{
	QString file = QFileDialog::getSaveFileName(this, "Enregistrer un fichier", QString(), "Text (*.txt)");
	if(file!=NULL)
	{
		std::ofstream savefile((file.toStdString()).c_str(), ios::out);
		
		if(savefile)
		{
			savefile << "Parameters for CCsegtool : " <<endl;
			
			savefile << "MPSDisplacement : " << MPSDisplacement->text().toStdString() <<endl;
			savefile << "Number of iteration : " << Numberloop->text().toStdString() <<endl;
			savefile << "White Mater Value : " << WMintensity->text().toStdString() <<endl;
			savefile << "Lambda max : " << Lambdamax->text().toInt() <<endl;
			savefile << "Phase stopping coefficient average distance across curve per voxel : "
					<< Coefofoptim->text().toStdString() <<endl;
			
			savefile << "Vessel remove : " << VesselRemoveOn->isChecked() <<endl;
			savefile << "Rotation 90 : " << permute_x_y->isChecked() <<endl;
			savefile << "Reflect X : " << ReflectXOn->isChecked() <<endl;
			savefile << "Reflect Y : " << ReflectYOn->isChecked() <<endl;
			savefile << "Opening : " << OpeningOn->isChecked() <<endl;
			savefile << "Double the size : " << DoubleOn->isChecked() <<endl;
			savefile << "Use other componant : " << othercomponant->isChecked() <<endl;
			savefile << "Linear interpolation : " << interpolationlinear->isChecked() <<endl;
			savefile << "Unconstrained : " << UnConstrained->isChecked() <<endl;
			savefile << "Segmentation Label : " << SegLabel->isChecked() <<endl;
			savefile << "Slice Direction : " << SliceDir->text().toInt() <<endl;
			savefile << "Midsagittal Plane Slice : " << Midsaggitalplaneslice->text().toStdString() <<endl;
			savefile << "Number of average : " << averagenum->text().toInt() <<endl;
			savefile << "Fix Similarity XForm : " << FSXform->isChecked() <<endl;
			savefile << "PSDistance : " << PSDistance->text().toFloat() <<endl;
			savefile << "Number for drawing : " << number_pts->text().toInt() <<endl;
			savefile << "Add an angle : " << m_angle <<endl;
			//A
			if((A->text().toStdString()).compare("default")==0)
				savefile << "Coef A : " << "default" <<endl;
			else
				savefile << "Coef A : " << A->text().toDouble() <<endl;
			//r
			savefile << "Coef r : " << r->text().toDouble() <<endl;
			
			if(m_Debug)
				std::cout<<"-----Save parameters done-----"<<std::endl;
			savefile.close();
		}
		else
			std::cout<<"ERROR: Problem to open the file for saving parameters"<<std::endl;
	}
}


/********************************************************************************* 
 * load the parameters 
 ********************************************************************************/
void CCsegtoolguiwindow::openparam()
{
	QString filename = QFileDialog::getOpenFileName(this, "Open File", QString(), "Text (*.txt)");
	
	if(filename!=NULL)
	{
		std::ifstream file((filename.toStdString()).c_str() , ios::in);  // open in reading
		std::string str,buf1,buf2;
	
		if(file)  // if open
		{
			//the first line
			getline(file, buf1);
			if(buf1.compare(0,27,"Parameters for CCsegtool : ")==0)
			{
				/* Loop for reading the file and setting the parameters values */
				while(!file.eof())
				{
					getline(file, buf1);
					if(buf1.compare(0,18,"MPSDisplacement : ")==0)
						MPSDisplacement->setText((buf1.substr(18,buf1.size()-18)).c_str());
					
					else if(buf1.compare(0,22,"Number of iteration : ")==0)
						Numberloop->setText((buf1.substr(22,buf1.size()-22)).c_str());
					
					else if(buf1.compare(0,20,"White Mater Value : ")==0)
						WMintensity->setText((buf1.substr(20,buf1.size()-20)).c_str());
					
					else if(buf1.compare(0,13,"Lambda max : ")==0)
						Lambdamax->setText((buf1.substr(13,buf1.size()-13)).c_str());
					
					else if(buf1.compare(0,69,"Phase stopping coefficient average distance across curve per voxel : ")==0)
						Coefofoptim->setText((buf1.substr(69,buf1.size()-69)).c_str());
					
					else if(buf1.compare(0,16,"Vessel remove : ")==0)
					{
						buf2=buf1.substr(16,buf1.size()-16);
						if(buf2.compare("0")==0)
							VesselRemoveOn->setCheckState(Qt::Unchecked);
						else
							VesselRemoveOn->setCheckState(Qt::Checked);
					}
					
					else if(buf1.compare(0,25,"Permutation of axe X/Y : ")==0)
					{
						buf2=buf1.substr(25,buf1.size()-25);
						if(buf2.compare("0")==0)
							permute_x_y->setCheckState(Qt::Unchecked);
						else
							permute_x_y->setCheckState(Qt::Checked);
					}
					
					else if(buf1.compare(0,12,"Reflect X : ")==0)
					{
						buf2=buf1.substr(12,buf1.size()-12);
						if(buf2.compare("0")==0)
							ReflectXOn->setCheckState(Qt::Unchecked);
						else
							ReflectXOn->setCheckState(Qt::Checked);
					}
					
					else if(buf1.compare(0,12,"Reflect Y : ")==0)
					{
						buf2=buf1.substr(12,buf1.size()-12);
						if(buf2.compare("0")==0)
							ReflectYOn->setCheckState(Qt::Unchecked);
						else
							ReflectYOn->setCheckState(Qt::Checked);
					}
					
					else if(buf1.compare(0,10,"Opening : ")==0)
					{
						buf2=buf1.substr(10,buf1.size()-10);
						if(buf2.compare("0")==0)
							OpeningOn->setCheckState(Qt::Unchecked);
						else
							OpeningOn->setCheckState(Qt::Checked);
					}
					
					else if(buf1.compare(0,18,"Double the size : ")==0)
					{
						buf2=buf1.substr(18,buf1.size()-18);
						if(buf2.compare("0")==0)
							DoubleOn->setCheckState(Qt::Unchecked);
						else
							DoubleOn->setCheckState(Qt::Checked);
					}
					
					else if(buf1.compare(0,23,"Linear interpolation : ")==0)
					{
						buf2=buf1.substr(23,buf1.size()-23);
						if(buf2.compare("0")==0)
							interpolationlinear->setCheckState(Qt::Unchecked);
						else
							interpolationlinear->setCheckState(Qt::Checked);
					}
					
					else if(buf1.compare(0,16,"Unconstrained : ")==0)
					{
						buf2=buf1.substr(16,buf1.size()-16);
						if(buf2.compare("0")==0)
							UnConstrained->setCheckState(Qt::Unchecked);
						else
							UnConstrained->setCheckState(Qt::Checked);
					}
					
					else if(buf1.compare(0,21,"Segmentation Label : ")==0)
					{
						buf2=buf1.substr(21,buf1.size()-21);
						if(buf2.compare("0")==0)
							SegLabel->setCheckState(Qt::Unchecked);
						else
							SegLabel->setCheckState(Qt::Checked);
					}
					
					else if(buf1.compare(0,18,"Slice Direction : ")==0)
						SliceDir->setText((buf1.substr(18,buf1.size()-18)).c_str());
					
					else if(buf1.compare(0,25,"Midsagittal Plane Slice : ")==0)
						Midsaggitalplaneslice->setText((buf1.substr(25,buf1.size()-25)).c_str());
					
					else if(buf1.compare(0,20,"Number of average : ")==0)
						averagenum->setText((buf1.substr(20,buf1.size()-20)).c_str());
					
					else if(buf1.compare(0,23,"Fix Similarity XForm : ")==0)
					{
						buf2=buf1.substr(23,buf1.size()-23);
						if(buf2.compare("0")==0)
							FSXform->setCheckState(Qt::Unchecked);
						else
							FSXform->setCheckState(Qt::Checked);
					}
					
					else if(buf1.compare(0,13,"PSDistance : ")==0)
						PSDistance->setText((buf1.substr(13,buf1.size()-13)).c_str());
					
					else if(buf1.compare(0,21,"Number for drawing : ")==0)
						number_pts->setText((buf1.substr(21,buf1.size()-21)).c_str());
					
					else if(buf1.compare(0,9,"Coef A : ")==0)
						A->setText((buf1.substr(9,buf1.size()-9)).c_str());
					
					else if(buf1.compare(0,9,"Coef r : ")==0)
						r->setText((buf1.substr(9,buf1.size()-9)).c_str());
					
					else if(buf1.compare(0,22,"Use other componant : ")==0)
						r->setText((buf1.substr(22,buf1.size()-22)).c_str());
					
					else if(buf1.compare(0,15,"Add an angle : ")==0)
						r->setText((buf1.substr(15,buf1.size()-15)).c_str());
					
				}
			}
			else
				std::cout<<"ERROR: Wrong file for parameters"<<std::endl;
			
			file.close();
		}
		else cerr << "ERROR: No parameters file found" << endl;
	}
}

/********************************************************************************* 
 * Save an output image
 ********************************************************************************/
void CCsegtoolguiwindow::saveimage()
{
	QString type;
	QString file = QFileDialog::getSaveFileName(this, "Save an Image", QString(),
						"Images (*.png *.xpm *.jpg)", &type, QFileDialog:: DontUseNativeDialog);
	if((m_viewimage->getPixmap())->isNull())
		std::cout<<"ERROR: No Image"<<std::endl;
	else
		if((m_viewimage->getPixmap())->save(file) && m_Debug)
			std::cout<<"-----Save image done-----"<<std::endl;
}


/********************************************************************************* 
 * Update the drawing with the Parameters
 ********************************************************************************/
void CCsegtoolguiwindow::update()
{
	setvalue();
	normalSize();
	
	if( (X->text()).compare("")!=0 || (Y->text()).compare("")!=0 || 
		    (Rotation->text()).compare("")!=0 || (Scale->text()).compare("")!=0 )
		updatedrawing(m_lastX, m_lastY, m_lastScale, m_lastRot, this, m_last_parameters, m_viewimage,
			      m_number_pts, m_outputfolder, m_nameofproject, m_Lambdamax, m_Debug);
}


/********************************************************************************* 
 * Update Visualization 
 ********************************************************************************/
void CCsegtoolguiwindow::updateVisuImage()
{
	if(m_viewimage!=NULL)
	{
		setvalue();
		normalSize();
		m_viewimage->setVisuFlip(m_rot90 ,m_rot180,m_rot270);
		
		if((m_viewimage->Getsize()).width()>(m_viewimage->Getsize()).height())
		{
			if(m_rot90 || m_rot270)
			{
				m_scrollArea->setFixedSize((m_viewimage->Getsize()).height(),
						(m_viewimage->Getsize()).width());
				/* Set on the Imageview */
				m_scrollArea->setWidget(m_viewimage);
				m_scrollArea->move(110,85);
			}
			else
			{
				m_scrollArea->setFixedSize((m_viewimage->Getsize()).width(),
						(m_viewimage->Getsize()).height());
				/* Set on the Imageview */
				m_scrollArea->setWidget(m_viewimage);
				m_scrollArea->move(40,145);
			}
		}
		else
		{
			if(m_rot90 || m_rot270)
			{
				m_scrollArea->setFixedSize((m_viewimage->Getsize()).height(),
						(m_viewimage->Getsize()).width());
				/* Set on the Imageview */
				m_scrollArea->setWidget(m_viewimage);
				m_scrollArea->move(40,145);
			}
			else
			{
				m_scrollArea->setFixedSize((m_viewimage->Getsize()).width(),
						(m_viewimage->Getsize()).height());
				/* Set on the Imageview */
				m_scrollArea->setWidget(m_viewimage);
				m_scrollArea->move(110,85);
			}
		}
		this->repaint();
		this->update();
	}
	else
		QMessageBox::information(this, "Warning", "No image");
}

/********************************************************************************* 
 * Reset points for repulsive and attractive effect
 ********************************************************************************/
void CCsegtoolguiwindow::resetpts()
{
	normalSize();
	m_viewimage->resetR_Apts();
}


/********************************************************************************* 
 * Reset points for repulsive and attractive effect
 ********************************************************************************/
void CCsegtoolguiwindow::setisoline()
{
	m_isoline = isoline->sliderPosition();
	QString qstr;
	qstr.setNum(m_isoline);
	qstr.append("%GOFmean");
	isolinevalue->setText(qstr);
	GOFmean->setText(qstr.setNum(CalculBestfitmean(m_outputfolder, m_nameofproject)));
}


/********************************************************************************* 
 * calcul contour points and draw them
 ********************************************************************************/
void CCsegtoolguiwindow::changeisoline()
{
	setisoline();
	if(m_runcalled)
		m_viewimage->RepulsionZone(m_isoline, m_A, m_r, CalculBestfitmean(m_outputfolder, m_nameofproject));
}


/********************************************************************************* 
 * change contour points if r or A change and draw them
 ********************************************************************************/
void CCsegtoolguiwindow::changecoef()
{
	//A
	m_A = A->text().toStdString();
	//r
	m_r = r->text().toDouble();
	if(m_runcalled)
		m_viewimage->RepulsionZone(m_isoline, m_A, m_r, CalculBestfitmean(m_outputfolder, m_nameofproject));
}

/********************************************************************************* 
 * Functions for Event : filter of MouseEvent + MouseEvent
 ********************************************************************************/
bool CCsegtoolguiwindow::eventFilter(QObject *object, QEvent *event)
{
	if (object == m_viewimage && event->type()==QEvent::MouseButtonRelease)
	{
		QMouseEvent *mouseEvent = static_cast<QMouseEvent *>(event);
		mouseMoveEvent(mouseEvent);
	}
	else if(object == m_viewimage && event->type()==QEvent::Paint)
	{
		QPaintEvent *paintEvent = static_cast<QPaintEvent *>(event);
		m_viewimage->paintEvent(paintEvent);
	}
	else if(object == m_viewimage && event->type()==QEvent::Wheel)
	{
		QWheelEvent *wheelEvent = static_cast<QWheelEvent *>(event);
		mouseWheelEvent(wheelEvent);
	}
	return false;
}

void CCsegtoolguiwindow::mouseMoveEvent(QMouseEvent* event)
{
	/* Show pixel value */
	if(event->button() == Qt::LeftButton)
	{
		if(m_Debug)
			std::cout<<"MouseLeftButton clicked"<<endl;
		QString qstr,qstr2,qstr3;
		QString instruct = "Pixel (";
		int value;
		value = m_viewimage->CalculPixelvalue(event->x(), event->y());
		instruct.append(qstr.setNum(static_cast<int>(event->x()/m_scaleFactor))+","+
				qstr2.setNum(static_cast<int>(event->y()/m_scaleFactor))+") = ");
		instruct.append(qstr3.setNum(value));
		pixellabel->setText(instruct);
	}
	if(event->button() == Qt::RightButton && m_runcalled)
	{
		m_viewimage->drawrepul_attrac_point(event->x(),event->y(), isoline->sliderPosition(), m_A, m_r,
				CalculBestfitmean(m_outputfolder, m_nameofproject));
	}
}

void CCsegtoolguiwindow::mouseWheelEvent(QWheelEvent *event)
{
	/* Zoom */
	if(event->delta() > 0)
	{
		if(m_Debug)
			std::cout<<"MouseRightButton clicked"<<endl;
		zoomIn();
	}
	if(event->delta() < 0)
	{
		if(m_Debug)
			std::cout<<"MouseMiddleButton clicked"<<endl;
		zoomOut();
	}
}


/********************************************************************************* 
 * Functions for the zoom, action with the ScrollArea
 ********************************************************************************/
void CCsegtoolguiwindow::zoomIn()
{
	scaleImage(1.25);
}

void CCsegtoolguiwindow::zoomOut()
{
	scaleImage(0.8);
}

void CCsegtoolguiwindow::scaleImage(double factor)
{
	m_scaleFactor *= factor;
	if( (static_cast<int>(m_scaleFactor * ( m_viewimage->Getsize() ).width()) < 5000) && 
		   (static_cast<int> ( m_scaleFactor * ( m_viewimage->Getsize() ).height() ) < 5000)
		&& (static_cast<int>(m_scaleFactor * ( m_viewimage->Getsize() ).width()) >=
		   (m_viewimage->Getsize()).width()) && 
		   (static_cast<int> ( m_scaleFactor * ( m_viewimage->Getsize() ).height() ) >= 
		   (m_viewimage->Getsize()).height()))
	{
		m_viewimage->Zoom(m_scaleFactor);
		
		adjustScrollBar(m_scrollArea->horizontalScrollBar(), factor);
		adjustScrollBar(m_scrollArea->verticalScrollBar(), factor);
	}
	
	if(m_Debug)
		std::cout<<"ZoomFactor : "<<m_scaleFactor<<std::endl;
}

void CCsegtoolguiwindow::adjustScrollBar(QScrollBar *scrollBar, double factor)
{
	scrollBar->setValue(static_cast<int>(factor * scrollBar->value() + ((factor - 1) * scrollBar->pageStep()/2)));
}

void CCsegtoolguiwindow::normalSize()
{
	m_viewimage->normalsize();
	m_scaleFactor = 1.0;
}
