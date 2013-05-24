#include "CCsegtool_computation.h"

#define MIN(x,y) ((x) < (y)) ? (x) : (y)

/***************************************************************************** 
 * Constructor
 *****************************************************************************/
CCsegtool_computation::CCsegtool_computation(CCsegtool_parameters * parameters, std::string nameofproject,
                                             std::string path_output, bool Debug, double A, double r)
{
	/* Init debug */
	debug = Debug;
	m_nameofproject=nameofproject;
	m_path_output=path_output;
	m_A=A;
	m_r=r;
	m_rot=-1;
	m_scale=-1;
	m_firstTime=true;
	/* Memory allocation */
	m_coefs = new Point4[MAX_COEFS];
	m_sample_pts = new Point2[MAX_PROFS];
	m_old_sample_pts = new Point2[MAX_PROFS];
	m_normals = new Point2[MAX_PROFS];
	
	/*Start of initialization*/
	m_parameters = parameters;
	
	m_num_coefs   = MIN(MAX_COEFS, m_parameters->GetSMMeansize(2)-1 );
	m_num_modes   = MIN(MAX_MODES, m_parameters->GetEigenvectorssize(2));
	m_num_samples =  m_parameters->GetPMMeansize(1);
	m_proflen     = MIN(MAX_PROFLEN, m_parameters->GetPMMeansize(2));
	
	m_num_coefs = MIN( (int)m_num_coefs, m_parameters->GetEigenvectorssize(3));
	m_parameters->SetMPSDisplacement( MIN(m_parameters->GetMPSDisplacement(), 
					  (m_parameters->GetPMMeansize(2)-2)));
	if (m_proflen%2 == 0) m_proflen -= 1;		/* ensure proflen is odd */
	
	/* Precalculate sine/cosine tables (fourier_funcs) */
	sincos_table(m_num_samples);
	
	/* we will fetch profiles of this length out of the image */
	m_ext_proflen = 3*m_proflen;
	m_num_modes = MIN((int)m_num_modes, m_parameters->GetBoundsize());
	
	if( !parameters->Getcoefsgiven() )
	{
		parse_coefs(m_parameters, m_num_coefs, m_coefs, m_centroid);
	}
	else
	{
		m_centroid[0] = 0.0;
		m_centroid[1] = 0.0;
		Setcoefs(parameters->GetCoefs());
		Setnum_coefs(parameters->GetNumcoefs());
	}
	/* remove global transform */
	normalize_coefs(m_coefs, m_num_coefs, &m_rot, &m_scale);
	
	/* Setup global similarity transform */
	m_centroid[0] += m_parameters->GetX();
	m_centroid[1] += m_parameters->GetY();
	m_rot = m_parameters->GetRotation();
	m_scale = m_parameters->GetScale();
	
	/* apply global transform to coefficients */
	unnormalize_coefs(m_parameters, m_coefs, m_num_coefs, m_rot, m_scale);

	if(debug)
	{
		std::cout<<"  ------Variables------"<<std::endl;
		std::cout<<"Ext_proflen : "<<m_ext_proflen<<std::endl;
		std::cout<<"Number of coefs : "<<m_num_coefs<<std::endl;
		std::cout<<"Number of modes : "<<m_num_modes<<std::endl;
		std::cout<<"Number of samples :"<<m_num_samples<<std::endl;
		std::cout<<"Prolen : "<<m_proflen<<std::endl;
		std::cout<<"Centroid : "<<m_centroid[0]<<" "<<m_centroid[1]<<std::endl;
		
		std::cout<<"  ------Parameters------"<<std::endl;
		std::cout<<"  Center X :"<<m_parameters->GetX()<<std::endl;
		std::cout<<"  Center Y :"<<m_parameters->GetY()<<std::endl;
		std::cout<<"  Rotation :"<<m_parameters->GetRotation()<<std::endl;
		std::cout<<"  Scale    :"<<m_parameters->GetScale()<<std::endl;
		std::cout<<"  --Constructor computation done--"<<endl;
	}
	
}

void CCsegtool_computation::setRepulPoints(std::vector<double> *RepulPoints)
{
	if(RepulPoints!=NULL)
		m_RepulPoints.clear();
	
	 for(unsigned int i = 0 ; i < RepulPoints->size() ; i++)
		 	 m_RepulPoints.push_back((*RepulPoints)[i]);
}


/***************************************************************************** 
 * main fonction
 *****************************************************************************/
float CCsegtool_computation::execution(CCsegtool_parameters * parameters, int step,int lambdamax, int iteration, bool update)
{
	/*if this is a second or more step, reinitialization of the pointeur parameters to keep the changement*/
	if(step>1)
	{
		m_parameters=parameters;
	}
	
	/* Reset the table with just 100 pts because we change it in reconstructCurve */
	sincos_table(m_num_samples);
	
	/* Step One: deform contour in image space */
	if(debug)
		std::cout<<"  -fourier_reconst"<<std::endl;
	fourier_reconst(m_coefs, m_num_coefs, m_centroid, m_num_samples, m_sample_pts);
	
	if(debug)
		std::cout<<"  -compute_normals"<<std::endl;
	compute_normals(m_coefs, m_num_coefs, m_num_samples, m_normals);
	
	if(debug)
		std::cout<<"  -get_profiles"<<std::endl;
	get_profiles();
	
	if(debug)
		std::cout<<"  -calc shift"<<std::endl;
	calcshift();
	
	if(debug)
		std::cout<<"  -fourier_interp"<<std::endl;
	fourier_interp(m_sample_pts, m_num_samples, m_num_coefs, m_coefs, m_centroid);
	
	/* Step Two: restrict deformation in eigenspace */
	if (m_parameters->GetUnconstrained()==false) {
		if(debug)
			std::cout<<"  -normalize_coefs"<<std::endl;
		if(m_rot==-1 && m_scale==-1)
		{
			m_rot=0;
			m_scale=0;
		}
		normalize_coefs(m_coefs, m_num_coefs, &m_rot, &m_scale);
		
		/*Reset if fix_simil*/
		if (m_parameters->GetFSXForm())
		{
			m_centroid[0] = m_parameters->GetX();
			m_centroid[1] = m_parameters->GetY();
			m_rot=m_parameters->GetRotation();
			m_scale=m_parameters->GetScale();
		}
		
		if(debug)
			std::cout<<"  -subtractmean"<<std::endl;
		subtractmean(m_coefs);
		
		if(debug)
			std::cout<<"  -proj_eigvec"<<std::endl;
		proj_eigvec(m_coefs);
		
		if(debug)
			std::cout<<"  -restrict_eigloads"<<std::endl;
		restrict_eigloads(lambdamax);
		
		if(debug)
			std::cout<<"  -dotprod_eigvec"<<std::endl;
		dotpord_eigvec( m_coefs);
		
		if(debug)
			std::cout<<"  -addmean"<<std::endl;
		addmean(m_coefs);
		
		if(debug)
			std::cout<<"  -unnormalize coefs"<<std::endl;
		unnormalize_coefs(m_parameters, m_coefs, m_num_coefs, m_rot, m_scale);
	}
	format_coefs(m_coefs, m_centroid, m_num_coefs, m_parameters);
	
	if(debug)
		std::cout<<"  -fourier_reconst"<<std::endl;
	fourier_reconst(m_coefs, m_num_coefs, m_centroid, m_num_samples, m_sample_pts);
	
	if(!update)
	{
		/* Calcul for the optimization */
		if(debug)
			std::cout<<"  -calcul displacement"<<std::endl;
		float disp;
		if(step==1 && iteration==50)
			disp = 10000000;
		else
		{
			disp = calc_dp(m_sample_pts, m_old_sample_pts);
			m_dp.push_back(disp);
		}
		
		if(debug)
		{
			/* Write the output */
			std::cout<<"  -write file .asc"<<std::endl;
			writefileasc(step,iteration);
		}
		return disp;
	}
	
	if(debug)
		std::cout<<"  --Computation done--"<<std::endl;
	return 0;
}


/*****************************************************************************
* Get 1D grayscale profiles from an image.  Returns a 2D array (row-major
* order) of floats.  Image can be in colmaj or rowmaj; set the flag.
 *****************************************************************************/
void CCsegtool_computation::get_profiles()
{
  std::vector<float> profiles;
	Point2 shade_pt;
	for (unsigned int i=0; i<m_num_samples; i++) {
		for (unsigned int j = 0; j<m_ext_proflen; j++) {
			for (int xy=0; xy<2; xy++) {
				shade_pt[xy] = m_sample_pts[i][xy] + (j-m_ext_proflen/2.0) * (m_normals[i][xy]) *
						m_parameters->GetPSDistance();
			}
			profiles.push_back(get_shade(&shade_pt));
		}
	}
	m_parameters->SetIProfiles(profiles);
}

/* Internal functions used by get_profiles */
float CCsegtool_computation::get_shade(Point2 *pt)
{
	double shade;
	int x,y;
	ExtractImageType::PixelType f00,f10,f01,f11;
	double dx,dy;
	x = (int) (*pt)[0];
	y = (int) (*pt)[1];
	dx = (*pt)[0]-x;
	dy = (*pt)[1]-y;
	f00 = (m_parameters->GetImage())->GetPixel(getindex(x, y));
	f10 = (m_parameters->GetImage())->GetPixel(getindex(x+1, y));
	f01 = (m_parameters->GetImage())->GetPixel(getindex(x ,y+1));
	f11 = (m_parameters->GetImage())->GetPixel(getindex(x+1, y+1));
	shade = (float) f00 + (f10 - f00)*dx + (f01 - f00)*dy + (f11 + f00 - f10 - f01)*dx*dy;
	return shade;
}

/* Internal functions used by get_shade */
ExtractImageType::IndexType CCsegtool_computation::getindex(int x, int y)
{
	ExtractImageType::IndexType pixelIndex;
	pixelIndex[0] = x; // x position
	pixelIndex[1] = y; // y position
	return pixelIndex;
}

/***************************************************************************** 
 * Internal routine to evaluate how well a given profile matches its model.
 * Uses the square of the Mahalanobis distance.
 * Returns a small value for a good fit.
 ****************************************************************************/
double CCsegtool_computation::goodness_of_fit(float* profile,
					      int proflen, 
					      int pointID)
{
	int i,j;
	double maha=0.0;
	for (i=0; i<proflen; i++)
	#if 0
		maha += profile[i]*profile[i];
	#else
		for (j=0; j<proflen; j++)
				maha += abs((int)(m_parameters->GetSigmaInvval(pointID,j,i) * profile[i] * profile[j]));
	#endif
	return maha;
}

/***************************************************************************** 
 * Find optimal shift amount for a single given profile.
 * Uses the most simpleminded approach: step through integer shift amounts.
 * The input test profile should be longer than the model mean, to allow room
 * for shifting the model back and forth.  
 * More precisely, it should be (proflen + 2*max_shift) long.
 *****************************************************************************/
int CCsegtool_computation::profile_shift(float *imgprof,	/* input test profile */
					 double *GOFvalue,  	/* Value of Output goodness-of-fit*/
					 int proflen,		/* input size of profile[] */
					 float grey_scaling, 	/* input scale factor on profiles */
					 int pointID)		/* ID of the point for the profile */
{
	int i, firsttime=1;
	int shift, bestshift;
	double testfit, bestfit;
	float testprof[MAX_PROFLEN];
	
	/* Initialization */
	testfit=0;
	bestfit=0;
	bestshift=0;
	bool recordMatch = false;
	
	for (shift = -m_parameters->GetMPSDisplacement(); shift <= m_parameters->GetMPSDisplacement(); shift++) {
		/* Read in the test profile at this shift */
		for (i=0; i<proflen; i++) testprof[i] = imgprof[i +  (m_ext_proflen - m_proflen)/2 + shift];
	
		for (i=0; i<proflen; i++) testprof[i] *= grey_scaling;
		
		/* Subtract the model mean */
		for (i=0; i<proflen; i++) testprof[i] -= m_parameters->GetPMMeanval(pointID,i);
		/* Record how well this test profile matches */
		double gof = goodness_of_fit(testprof, proflen, pointID);
		double penalty = computeRepulsionExponentialPenalty(pointID,shift);
		testfit = gof + penalty;
		if (penalty > 1 && firsttime ) {
			recordMatch = true;
		}
		if(debug)
		{
			if (recordMatch) {
				std::cout  << "PointGOF" << pointID << "[" << shift << "] = " << gof   << "; ";
				std::cout  << "PointPenalty" << pointID << "[" << shift << "] = "  << penalty  << "; ";
				std::cout  << "PointTotal" << pointID << "[" << shift << "] = " << testfit   << "; " <<
				std::endl;
			}
		}
		if (firsttime || (testfit < bestfit)) {
			firsttime = 0;
			bestfit = testfit;
			bestshift = shift;
		}
	}
	*GOFvalue=bestfit;
	return bestshift;
}

std::vector<float>::iterator CCsegtool_computation::getiterator(std::vector<float> vect, int i)
{
	std::vector<float>::iterator it;
	int counter=0;
	for(it=vect.begin();it<vect.end();it++)
	{
		if(counter==i)
			return it;
		counter++;
	}
	return vect.end();
}

/***************************************************************************** 
 * functions to calcul the penalty of repulsive points
 *****************************************************************************/
double CCsegtool_computation::computeRepulsionExponentialPenalty(int pointID, double shift)
{
	double X,Y;
	double penalty = 0.0;
	if(!m_RepulPoints.empty())
	{
		/* Calcul of the X,Y for the pointID */
		X = m_sample_pts[pointID][0];
		Y = m_sample_pts[pointID][1];
		X += m_normals[pointID][0]*shift*(m_parameters->GetPSDistance());
		Y += m_normals[pointID][1]*shift*(m_parameters->GetPSDistance());
		/* Sum of the exponential for each repulsive point */
		for(unsigned int i=0;i<m_RepulPoints.size();i=i+2)
			penalty+=exponential(X,Y,i);
	}
	else if(debug)
		std::cout<<" No repulsion points"<<std::endl;
	
	return penalty;
}

double CCsegtool_computation::exponential(double X, double Y, int pointRepulsivID)
{
	double distance;
	/* Calcul of the distance */
	distance = sqrt((m_RepulPoints[pointRepulsivID]-X)*(m_RepulPoints[pointRepulsivID]-X) 
			+(m_RepulPoints[pointRepulsivID+1]-Y)*(m_RepulPoints[pointRepulsivID+1]-Y));
	//std::cout<<"Distance : "<<distance<<std::endl;
	return m_A*exp(-distance*m_r);
}

/*****************************************************************************
 * calcul shift
 *****************************************************************************/
void CCsegtool_computation::calcshift()
{
	int shift;
	double value;
	std::vector<float>  prof_shifts;
	std::vector<double> match;
	/* Keep the old value */
	for (unsigned int i=0; i<m_num_samples; i++) {
		m_old_sample_pts[i][0]=m_sample_pts[i][0];
		m_old_sample_pts[i][1]=m_sample_pts[i][1];
	}
	/* change the nw one with the new value with the calcul shift */
	for (unsigned int i=0; i<m_num_samples; i++) {
		shift = profile_shift(&m_parameters->GetIProfiles()[i*m_ext_proflen], &value, m_proflen,
				       255.0 / m_parameters->GetWMIntensity(), i);
		match.push_back(value);
		m_sample_pts[i][0] += m_normals[i][0]*shift*(m_parameters->GetPSDistance());
		m_sample_pts[i][1] += m_normals[i][1]*shift*(m_parameters->GetPSDistance());
		prof_shifts.push_back(shift*(m_parameters->GetPSDistance()));
	}
	m_parameters->SetGoodnessoffit(match);
	m_parameters->SetPPShifts(prof_shifts);
}


/*****************************************************************************
 * subtract mean.  Note i+1: first row of mean_coefs is 
 *****************************************************************************/
void CCsegtool_computation::subtractmean(Point4* coefs)
{
	for (unsigned int i=0; i<m_num_coefs; i++)
		for (unsigned int j=0; j<4; j++)
			coefs[i][j] -= m_parameters->GetSMMeanval(j,i+1);
}

/***************************************************************************** 
 * add mean.  Note i+1: first row of mean_coefs is 
 *****************************************************************************/
void CCsegtool_computation::addmean(Point4* coefs)
{
	for (unsigned int i=0; i<m_num_coefs; i++) 
		for (unsigned int j=0; j<4; j++)
			coefs[i][j] += m_parameters->GetSMMeanval(j,i+1);
}

/***************************************************************************** 
 * proj_eigvec
 *****************************************************************************/
void CCsegtool_computation::proj_eigvec(Point4* _coefs)
{
	std::vector<float> p;
	float buf[m_num_modes];
	for (unsigned int i=0; i<m_num_modes; i++) {
		buf[i]=0.;
		for (unsigned int j=0; j<m_num_coefs; j++) 
			for (unsigned int k=0; k<4; k++)
				buf[i] += m_parameters->GetEigenvectorsval(k,i,j) * _coefs[j][k];
	}
	for (unsigned int j=0; j<m_num_modes; j++) {
		p.push_back(buf[j]);
	}
	
	m_parameters->Seteigloads(p);
}


/***************************************************************************** 
 * Clamp the given eigload vector to the given magnitude bounds
 *****************************************************************************/
void CCsegtool_computation::restrict_eigloads(int lambdamax)
{
	for (unsigned int i=0; i<m_num_modes; i++) {
		if (m_parameters->Geteigloads()[i] > m_parameters->GetBound()[i] * lambdamax) 
		{
			if(debug)
				fprintf(stderr, "  mode %d: %f > %f\n", i,m_parameters->Geteigloads()[i],
					m_parameters->GetBound()[i] * lambdamax);
			m_parameters->Seteigloadsval(i, m_parameters->GetBound()[i] * lambdamax);
		}
		else if (m_parameters->Geteigloads()[i] < -(m_parameters->GetBound()[i] * lambdamax) ) 
		{
			if(debug)
				fprintf(stderr, "  mode %d: %f > %f\n", i, -(m_parameters->Geteigloads()[i]),
					m_parameters->GetBound()[i] * lambdamax);
			m_parameters->Seteigloadsval(i,-(m_parameters->GetBound()[i] * lambdamax));
		}
	}
}

/***************************************************************************** 
 * dotprod_eigvec
 *****************************************************************************/
void CCsegtool_computation::dotpord_eigvec(Point4* _coefs)
{
	for (unsigned int j=0; j<m_num_coefs; j++)
		for (unsigned int k=0; k<4; k++) {
			_coefs[j][k] = 0.;
			for (unsigned int i=0; i<m_num_modes; i++) 
				_coefs[j][k] += m_parameters->GetEigenvectorsval(k,i,j) * m_parameters->Geteigloads()[i];
		}
}

/***************************************************************************** 
 * calcul the distance for a point between two iteration
 *****************************************************************************/
float CCsegtool_computation::calc_dp(Point2* sample_pts, Point2* old_sample_pts)
{
	float sum=0;
	/*initialization of sum*/
	for(unsigned int j=0; j<m_num_samples; j++)
		sum += sqrt((old_sample_pts[j][0]-sample_pts[j][0])*(old_sample_pts[j][0]-sample_pts[j][0]) +
				(old_sample_pts[j][1]-sample_pts[j][1])*(old_sample_pts[j][1]-sample_pts[j][1]) );
	return sum/m_num_samples;
}

/***************************************************************************** 
 * Reconstruction of the curve
 *****************************************************************************/
void CCsegtool_computation::reconstructCurve(Point2* pts, int number_pts)
{
	if(debug)
		std::cout<<"  --ReconstructCurve--"<<endl;
	if(number_pts<=256)
	{
		/* Do the reconstruction */
		sincos_table(number_pts);
		fourier_reconst(m_coefs, m_num_coefs, m_centroid, number_pts, pts);
	}
	else
		std::cout<<"Error: To much point for drawing the contour, 256 is the maximun."<<std::endl;
}

/***************************************************************************** 
 * acced to private parameters
 *****************************************************************************/
//num_coefs
void CCsegtool_computation::Setnum_coefs(unsigned int _num_coefs)
{
	m_num_coefs=_num_coefs;
}
unsigned int CCsegtool_computation::Getnum_coefs()
{
	return m_num_coefs;
}
//num_samples
void CCsegtool_computation::Setnum_samples(unsigned int _num_samples)
{
	m_num_samples=_num_samples;
}
unsigned int CCsegtool_computation::Getnum_samples()
{
	return m_num_samples;
}
//num_modes
void CCsegtool_computation::Setnum_modes(unsigned int _num_modes)
{
	m_num_modes=_num_modes;
}
unsigned int CCsegtool_computation::Getnum_modes()
{
	return m_num_modes;
}//proflen
unsigned int CCsegtool_computation::Getproflen()
{
	return m_proflen;
}
//ext_proflen
unsigned int CCsegtool_computation::Getext_proflen()
{
	return m_ext_proflen;
}
//Center
double CCsegtool_computation::GetcenterX()
{
	return m_centroid[0];
}
double CCsegtool_computation::GetcenterY()
{
	return m_centroid[1];
}
//coefs
void CCsegtool_computation::Setcoefs(Point4* _coefs)
{
	for(unsigned int i=0;i<m_num_coefs;i++)
		for(unsigned int j=0;j<4;j++)
			m_coefs[i][j]=_coefs[i][j];
	for(unsigned int i=m_num_coefs;i<MAX_COEFS;i++)
		for(unsigned int j=0;j<4;j++)
			m_coefs[i][j]=0.0;
}
Point4* CCsegtool_computation::Getcoefs()
{
	return m_coefs;
}
//sample_pts
void CCsegtool_computation::Setsample_pts(Point2* _sample_pts)
{
	for(unsigned int i=0;i<sizeof(_sample_pts);i++)
		for(unsigned int j=0;j<4;j++)
			m_sample_pts[i][j]=_sample_pts[i][j];
}
Point2* CCsegtool_computation::Getsample_pts()
{
	return m_sample_pts;
}
double CCsegtool_computation::Getsample_ptsval(int i, int j)
{
	return m_sample_pts[i][j];
}
//normals
void CCsegtool_computation::Setnormals(Point2* _normals)
{
	for(unsigned int i=0;i<sizeof(_normals);i++)
		for(unsigned int j=0;j<4;j++)
			m_normals[i][j]=_normals[i][j];
}
Point2* CCsegtool_computation::Getnormals()
{
	return m_normals;
}


/***************************************************************************** 
 * Write asc if debug call : coefficient for each iteration
 *****************************************************************************/
void CCsegtool_computation::writefileasc(int step, int iteration)
{
	std::string outputfilename;
	std::stringstream iter,sstep;
	/* put some 0 in front of the iteration number */	
	if((m_parameters->GetOiter()- iteration) <10)
	{
		iter.put('0');
		iter << m_parameters->GetOiter()- iteration;
	}
	else
		iter<< m_parameters->GetOiter()- iteration;
	
	sstep<<step;
	outputfilename=m_path_output +  "/" + m_nameofproject + "_step" + sstep.str() + "_iter" + iter.str() + ".asc";
	std::ofstream ofile(outputfilename.c_str() , ios::out | ios::trunc);
	
	if(ofile)
	{
		ofile<<"field 1D 4-vector uniform float"<<endl;
		ofile<<"11 ;"<<endl;
		for(unsigned int i=0;i<11;i++)
		{
			ofile<<m_parameters->GetSSCoefs()[i*4]
					<<" "<<m_parameters->GetSSCoefs()[i*4+1]
					<<" "<<m_parameters->GetSSCoefs()[i*4+2]
					<<" "<<m_parameters->GetSSCoefs()[i*4+3]<<endl;
		}
		ofile.close();
	}
	else
		cerr << "Error: opening output file !" << endl;

}

/***************************************************************************** 
 * Write Output: 2D.pts is a file with all of the points (100 pts)
 *****************************************************************************/
void CCsegtool_computation::writeoutput(double Xsize, double Ysize)
{
	/* 2D index */
	std::string outfilePDMName;
	outfilePDMName = m_path_output + "/" + m_nameofproject + "_2D.pts";
	ofstream bfile(outfilePDMName.c_str(), ios::out);
	bfile << m_num_samples << endl;
	for (unsigned int i = 0; i < m_num_samples; i++)
		bfile << m_sample_pts[i][0] << " " << m_sample_pts[i][1] << " " << endl;
	bfile.close();
	
	/* 3D index */
	m_parameters->transform2DIndex_3Dindex(m_sample_pts, m_num_samples, m_path_output, m_nameofproject, Xsize,Ysize);
}

/*****************************************************************************
 * Write Output: DP.pts is a file with the distance between the same point for
 * two differents iterations.
 *****************************************************************************/
void CCsegtool_computation::writedp(float coefoptim1, float coefoptim2, bool unconstrained)
{
	std::string outfilePDMName;
	int k = 0;
	outfilePDMName = m_path_output + "/" + m_nameofproject + "_DP.pts";
	ofstream bfile(outfilePDMName.c_str(), ios::out);
	bfile << "Variations between two iterations" <<endl;
	for (unsigned int i = 0; i < m_dp.size(); i++)
	{
		bfile << m_dp[i] << endl;
		if (m_dp[i] < coefoptim1)
		{
			if (k==0)
			{
				bfile << "Step 2 :" << endl;
				k++;
			}
			if(m_dp[i] < coefoptim2 && unconstrained==true)
				bfile << "Step 3 :" << endl;
		}
	}
	bfile.close();
}

/***************************************************************************** 
 * Get the parameters for the last iteration last step
 *****************************************************************************/
void CCsegtool_computation::getparam(double *lastX, double *lastY, double *lastScale, double *lastRot)
{
	/* Rotation and Scale */
	*lastScale = m_scale;
	*lastRot = m_rot;
	/* Centroid */
	*lastX = m_centroid[0];
	*lastY = m_centroid[1];
}

/***************************************************************************** 
 * Get the parameters for the last iteration last step
 *****************************************************************************/
void CCsegtool_computation::setparam(double X, double Y, double Scale, double Rot)
{
	m_centroid[0]=X;
	m_centroid[1]=Y;
	m_scale=Scale;
	m_rot=Rot;
}
