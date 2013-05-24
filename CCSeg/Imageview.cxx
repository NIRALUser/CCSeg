#include "Imageview.h"
#include <iostream>

/*********************************************************************************
 * Constructor
 ********************************************************************************/
Imageview::Imageview (QString imageFileName, bool rot90, bool rot180, bool rot270, bool debug, QWidget *parent,
Qt::WFlags f) : QLabel (parent, f)
{
	/* Init */
	m_debug=debug;
	m_parent = parent;
	m_sample_pts = new Point2[256];
	m_painteventzoom=false;
	m_painteventpoint=false;
	m_visurot90 = rot90;
	m_visurot180 = rot180;
	m_visurot270 = rot270;
	m_factor=1.0;
	
	if(imageFileName.compare("")!=0)
	{
		/* Dowmload the pixmap */
		if(m_graypixmap.load (imageFileName,"PNG"))
			m_Size = m_graypixmap.size();
		else
			std::cout<<"Error loading GrayPixmap."<<std::endl;
		m_pixmap = m_graypixmap;
	}
	/* set the pixmap */
	setPixmap(m_graypixmap);
	
	/* Set some informations about m_viewimage */
	setBackgroundRole(QPalette::Base);
	setSizePolicy(QSizePolicy::Ignored, QSizePolicy::Ignored);
	setScaledContents(true);
	
	setFixedSize(m_Size);
}

/*********************************************************************************
 * Function of QLabel paintEvent for drawing the curve with the pts
 ********************************************************************************/
void Imageview::paintEvent (QPaintEvent *)
{
	if(m_painteventzoom==true)
	{
		/* Zoompainter on the Widget */
		QPainter zoompainter(this);
		/* Rotate 90/180/270 degres the image for the visualization if the user ask for that */
		if(m_visurot90 || m_visurot180 || m_visurot270)
		{
			if(m_visurot90)
			{
				zoompainter.translate(m_zoompixmap.height(),0);
				zoompainter.rotate(90);
				setFixedSize(QSize(m_zoompixmap.height(),m_zoompixmap.width()));
			}
			if(m_visurot180)
			{
				zoompainter.translate(m_zoompixmap.width(),m_zoompixmap.height());
				zoompainter.rotate(180);
			}
			if(m_visurot270)
			{
				zoompainter.translate(0,m_zoompixmap.width());
				zoompainter.rotate(270);
				setFixedSize(QSize(m_zoompixmap.height(),m_zoompixmap.width()));
			}
		}
		zoompainter.drawPixmap(0,0,m_zoompixmap);
		zoompainter.end();
	}
	else if( m_painteventpoint==true)
	{
		/* First painter on the Pixmap */
		QPainter firstpaint(&m_repulpixmap);
		firstpaint.setPen(QPen (Qt::yellow, 2));
		for(unsigned int i=0;i<m_RepulPoints.size();i=i+2)
		{
			/* Draw the point */
			firstpaint.drawPoint(static_cast<int>(m_RepulPoints[i]), static_cast<int>(m_RepulPoints[i+1]));
		}
		firstpaint.end();
		
		/* Second painter on the Pixmap */
		QPainter secpaint(&m_repulpixmap);
		secpaint.setPen(QPen (Qt::yellow, 1));
		for(unsigned int j=0;j<m_RepulContourPoints.size();j=j+2)
		{
			/* Draw the point */
			secpaint.drawPoint(static_cast<int>(m_RepulContourPoints[j]),
					   static_cast<int>(m_RepulContourPoints[j+1]));
		}
		secpaint.end();
		
		Zoom(m_factor);
		/* Third painter on the Widget */
		QPainter thirdpainter(this);
		thirdpainter.drawPixmap(0,0,m_pixmap);
		/* Rotate 90/180/270 degres the image for the visualization if the user ask for that */
		if(m_visurot90 || m_visurot180 || m_visurot270)
		{
			if(m_visurot90)
			{
				thirdpainter.translate(m_pixmap.height(),0);
				thirdpainter.rotate(90);
				setFixedSize(QSize(m_pixmap.height(),m_pixmap.width()));
			}
			if(m_visurot180)
			{
				thirdpainter.translate(m_pixmap.width(),0);
				thirdpainter.rotate(180);
			}
			if(m_visurot270)
			{
				thirdpainter.translate(0,m_pixmap.width());
				thirdpainter.rotate(270);
				setFixedSize(QSize(m_pixmap.height(),m_pixmap.width()));
			}
		}
		thirdpainter.end();
	}
	else
	{
		/* First painter on the Pixmap */
		QPainter paint(&m_pixmap);
		paint.setPen(QPen (Qt::red, 1));
		/* Draw the curve */
		if(m_num_samples==1)
			paint.drawPoint(static_cast<int>(m_sample_pts[0][0]),static_cast<int>(m_sample_pts[0][1]));
		else
		{
			// loop to join all of the points
			for(int i=0;i<m_num_samples-1;i++)
			{
				paint.drawLine(static_cast<int>(m_sample_pts[i][0]),static_cast<int>(m_sample_pts[i][1]),
					static_cast<int>(m_sample_pts[i+1][0]),static_cast<int>(m_sample_pts[i+1][1]));
			}
			// Line between the last and the first point
			paint.drawLine(static_cast<int>(m_sample_pts[m_num_samples-1][0]),
				static_cast<int>(m_sample_pts[m_num_samples-1][1]),
				static_cast<int>(m_sample_pts[0][0]),
				static_cast<int>(m_sample_pts[0][1]));
		}
		/* Draw the center */
		paint.setPen(QPen (Qt::red, 2));
		paint.drawPoint(static_cast<int>(m_X), static_cast<int>(m_Y));
		paint.end();
		
		/*second painter if vector of pts is not empty */
		if(!m_RepulPoints.empty())
		{
			QPainter ptspaint(&m_pixmap);
			ptspaint.setPen(QPen (Qt::yellow, 2));
			for(unsigned int i=0;i<m_RepulPoints.size();i=i+2)
			{
				/* Draw the point */
				ptspaint.drawPoint(static_cast<int>(m_RepulPoints[i]),
						static_cast<int>(m_RepulPoints[i+1]));
			}
			ptspaint.end();
		}
		
		
		/* Second painter on the Widget */
		QPainter painter(this);
		/* Rotate 90/180/270 degres the image for the visualization if the user ask for that */
		if(m_visurot90 || m_visurot180 || m_visurot270)
		{
			if(m_visurot90)
			{
				painter.translate(m_pixmap.height(),0);
				painter.rotate(90);
				setFixedSize(QSize(m_pixmap.height(),m_pixmap.width()));
			}
			if(m_visurot180)
			{
				painter.translate(m_pixmap.width(),m_pixmap.height());
				painter.rotate(180);
			}
			if(m_visurot270)
			{
				painter.translate(0,m_pixmap.width());
				painter.rotate(270);
				setFixedSize(QSize(m_pixmap.height(),m_pixmap.width()));
			}
		}
		painter.drawPixmap(0,0,m_pixmap);
		painter.end();
		
		m_zoompixmap = m_pixmap;
		m_repulpixmap = m_pixmap;
	}
}


/*********************************************************************************
 *  Function of QLabel paintEvent for MouseEvent
 ********************************************************************************/
int Imageview::CalculPixelvalue(int  X, int Y)
{
	QRgb pix;
	QImage Image;
	Image = m_graypixmap.toImage();
	pix = Image.pixel( static_cast<int>(X/m_factor), static_cast<int>(Y/m_factor));
	QColor color(pix);
	return color.value();
}


void Imageview::setGrayPixmap(QString imageFileName)
{
	/* Dowmload the pixmap */
	if(m_graypixmap.load (imageFileName,"PNG"))
		m_Size = m_graypixmap.size();
	else
		std::cout<<"Error loading GrayPixmap."<<std::endl;
	m_pixmap = m_graypixmap;
	
	setPixmap(m_graypixmap);
	
	setFixedSize(m_Size);
}

QPixmap* Imageview::getPixmap()
{
	return &m_pixmap;
}

QPixmap* Imageview::getsavePixmap()
{
	return &m_savepixmap;
}

QSize Imageview::Getsize()
{
	return m_Size;
}

std::vector<double>* Imageview::GetRepulPoints()
{
	return &m_RepulPoints;
}

void Imageview::setVisuFlip(bool rot90, bool rot180, bool rot270)
{
	m_visurot90 = rot90;
	m_visurot180 = rot180;
	m_visurot270 = rot270;
}

/*********************************************************************************
 * Copy the number of pts and the pts
 ********************************************************************************/
void Imageview::GetPts(Point2 *sample_pts, int num_samples, double X, double Y)
{
	m_painteventzoom=false;
	m_painteventpoint=false;
	/* Copy the graypixmap in pixmap */
	m_pixmap = m_graypixmap;
	/* Center */
	m_X = X;
	m_Y = Y;
	/* number of pts */
	m_num_samples=num_samples;
	/* pts */
	for(int i=0; i<m_num_samples;i++)
	{
		m_sample_pts[i][0]=sample_pts[i][0];
		m_sample_pts[i][1]=sample_pts[i][1];
	}
	update();
}

/*********************************************************************************
 * Zoom
 ********************************************************************************/
void Imageview::Zoom(double factor)
{
	m_painteventzoom=true;
	m_painteventpoint=false;
	m_factor = factor;
	m_zoompixmap = m_repulpixmap.scaled(static_cast<int>(m_factor*m_repulpixmap.width()),
		static_cast<int>(m_factor*m_repulpixmap.height()),Qt::KeepAspectRatioByExpanding);
	setPixmap(m_zoompixmap);
	setFixedSize(m_zoompixmap.size());
}

void Imageview::normalsize()
{
	/* reset the variables and set the size of the window */
	setFixedSize(m_pixmap.size());
	m_factor=1.0;
	m_painteventzoom=false;
	m_painteventpoint=false;
	
	m_zoompixmap = m_repulpixmap;
}

/*********************************************************************************
 * Draw the repulsive or attractive point in yellow and the zone of the effect 
 * when it's more than isoline's value
 ********************************************************************************/
void Imageview::drawrepul_attrac_point(double X,double Y, int isoline, std::string A, double r,double bestfitmean)
{
	/* Pixmap for drawing repulsive point and the area */
	m_repulpixmap = m_pixmap;
	/* Point X,Y */
	/* transform the coordinate in function of the visualization rotation */
	if(m_visurot90)
	{
		m_RepulPoints.push_back(Y/m_factor); //X
		m_RepulPoints.push_back((m_pixmap.height()-X)/m_factor); //Y
	}
	else if(m_visurot180)
	{
		m_RepulPoints.push_back((m_pixmap.width()-X)/m_factor); //X
		m_RepulPoints.push_back((m_pixmap.height()-Y)/m_factor); //Y
	}
	else if(m_visurot270)
	{
		m_RepulPoints.push_back((m_pixmap.width()-Y)/m_factor); //X
		m_RepulPoints.push_back(X/m_factor); //Y
	}
	else
	{
		m_RepulPoints.push_back(X/m_factor); //X
		m_RepulPoints.push_back(Y/m_factor); //Y
	}
	
	/* Calcul the Zone of repulsion */
	RepulsionZone(isoline, A, r, bestfitmean);
	
	if(m_debug)
		std::cout<<"Repulsives Point ("<<m_RepulPoints[0]<<","<<m_RepulPoints[1]<<")"<<std::endl;
}

void Imageview::resetR_Apts()
{
	m_painteventzoom=false;
	m_painteventpoint=false;
	
	m_RepulPoints.clear();
	m_RepulContourPoints.clear();
	
	/* Copy the graypixmap in pixmap */
	m_pixmap = m_graypixmap;
	
	update();
}

void Imageview::RepulsionZone(int isoline,std::string A, double r,double bestfitmean)
{
	m_painteventzoom=false;
	m_painteventpoint=true;
	
	double pixelpenalty;
	if(!m_RepulPoints.empty())
	{
		m_RepulContourPoints.clear();
		m_repulpixmap = m_pixmap;
		
		if(m_debug)
			std::cout<<"Repulsion Zone"<<std::endl;
		
		/*coef A*/
		if(A.compare("default")==0)
			m_A = 100.0*bestfitmean;
		else if(A.compare("")==0)
			m_A = 0.0;
		else
			m_A = atof(A.c_str());
		m_r = r;
		
		if(m_debug)
		{
			std::cout<<"Bestfitmean : "<<bestfitmean<<std::endl;
			std::cout<<"A : "<<m_A<<std::endl;
			std::cout<<"r : "<<m_r<<std::endl;
			std::cout<<"Isoline : "<<isoline<<std::endl;
			std::cout<<"Isoline*bestfitmean : "<<isoline*bestfitmean<<std::endl;
		}
		
		/* Calcul of the penalty and points with a penalty > isoline */
		for(int i=1;i<(m_pixmap.size()).width()-1;i++)
		{
			for(int j=1;j<(m_pixmap.size()).height()-1;j++)
			{
				pixelpenalty = penaltyvalue(i, j);
				/********************************************************
				* calcul of the 8 neighboor of the pixel (i,j) :
				* if the penalty of one of them is < isoline*bestfitmean
				* this pixel is one of the contour
				*********************************************************/
				if( pixelpenalty > isoline*bestfitmean) {
					if(penaltyvalue(i-1, j) < isoline*bestfitmean || 
									      penaltyvalue(i-1, j+1)< isoline*bestfitmean ||
									      penaltyvalue(i, j+1)< isoline*bestfitmean || 
									      penaltyvalue(i+1, j+1)< isoline*bestfitmean ||
									      penaltyvalue(i+1, j)< isoline*bestfitmean ||
									      penaltyvalue(i+1, j-1)< isoline*bestfitmean ||
									      penaltyvalue(i, j-1)< isoline*bestfitmean ||
									      penaltyvalue(i-1, j-1)< isoline*bestfitmean)
					{
						if(m_debug)
						{
							std::cout<<"Point ("<<i<<","<<j
									<<") add to the contour of the repulsion zone"
									<<std::endl;
						}
						m_RepulContourPoints.push_back(i); //X
						m_RepulContourPoints.push_back(j); //Y
					}
				}
			}
		}
		update();
	}
}

/*********************************************************************************
 * Function calculating the penalty
 ********************************************************************************/
double Imageview::penaltyvalue(int X, int Y)
{
	double distance;
	double penalty = 0.0;
	if(!m_RepulPoints.empty())
	{
		for(unsigned int i=0;i<m_RepulPoints.size();i=i+2)
		{
			/* Calcul of the distance */
			distance = sqrt((m_RepulPoints[i]-X)*(m_RepulPoints[i]-X) 
					+(m_RepulPoints[i+1]-Y)*(m_RepulPoints[i+1]-Y));
			//std::cout<<"Distance : "<<distance<<std::endl;
			penalty += m_A*exp(-distance*m_r);
		}
	}
	return penalty;
}

/*********************************************************************************
 * Save the image with the segmentation
 ********************************************************************************/
void Imageview::SaveImage(QString file, Point2 *sample_pts, int num_samples, double X, double Y)
{
	/* Copy the graypixmap in pixmap */
	m_savepixmap = m_graypixmap;
	/* Center */
	m_X = X;
	m_Y = Y;
	/* number of pts */
	m_num_samples=num_samples;
	/* pts */
	for(int i=0; i<m_num_samples;i++)
	{
		m_sample_pts[i][0]=sample_pts[i][0];
		m_sample_pts[i][1]=sample_pts[i][1];
	}
	
	/* First painter on the Pixmap */
	QPainter paint(&m_savepixmap);
	paint.setPen(QPen (Qt::red, 1));
	/* Draw the curve */
	if(m_num_samples==1)
		paint.drawPoint(static_cast<int>(m_sample_pts[0][0]),static_cast<int>(m_sample_pts[0][1]));
	else
	{
			// loop to join all of the points
		for(int i=0;i<m_num_samples-1;i++)
		{
			paint.drawLine(static_cast<int>(m_sample_pts[i][0]),static_cast<int>(m_sample_pts[i][1]),
				       static_cast<int>(m_sample_pts[i+1][0]),static_cast<int>(m_sample_pts[i+1][1]));
		}
			// Line between the last and the first point
		paint.drawLine(static_cast<int>(m_sample_pts[m_num_samples-1][0]),
			       static_cast<int>(m_sample_pts[m_num_samples-1][1]),
			       static_cast<int>(m_sample_pts[0][0]),
			       static_cast<int>(m_sample_pts[0][1]));
	}
	
	/* Draw the center */
	paint.setPen(QPen (Qt::red, 2));
	paint.drawPoint(static_cast<int>(m_X), static_cast<int>(m_Y));
	paint.end();
	
}

