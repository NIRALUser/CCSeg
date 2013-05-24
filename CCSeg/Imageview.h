#ifndef _IMAGEVIEW_H_
#define _IMAGEVIEW_H_

#include <QWidget>
#include <QLabel>
#include <QPixmap>
#include <QImage>
#include <QPainter>
#include <QPoint>
#include <QColor>
#include <string>
#include <ostream>
#include <fstream>
#include <vector>
#include <math.h>

typedef double Point2[2];
/***********************************************************
* Imageview
* Display an image.
************************************************************/

class Imageview : public QLabel
{
	public:
		Imageview (QString imageFileName="", bool rot90 =false, bool rot180= false, bool rot270= false,
			   bool debug=false, QWidget *parent=NULL,Qt::WFlags f = 0);
		virtual ~Imageview(){};
		void setGrayPixmap(QString imageFileName);
		QPixmap *getPixmap();
		QSize Getsize();
		QPixmap* getsavePixmap();
		std::vector<double>* GetRepulPoints();
		void setVisuFlip(bool rot90, bool rot180, bool rot270);
		/* Zoom functions */
		void Zoom(double factor);
		void normalsize();
		/* PaintEvent */
		void paintEvent (QPaintEvent*);
		/* Drawing function */
		void GetPts(Point2 *sample_pts, int num_samples, double X, double Y);
		void SaveImage(QString file, Point2 *sample_pts, int num_samples, double X, double Y);
		int CalculPixelvalue(int  X, int Y);
		void drawrepul_attrac_point(double X ,double Y, int isoline, std::string A, double r,double bestfitmean);
		void resetR_Apts();
		void RepulsionZone(int isoline, std::string A, double r,double bestfitmean);
		double penaltyvalue(int X, int Y);
		
	private:
		/* Size of pixmap */
		QSize m_Size;
		
		/* Visualization */
		bool m_visurot90;
		bool m_visurot180;
		bool m_visurot270;
		
		/* Zoom */
		double m_factor;
		bool   m_painteventzoom;
		
		/* Draw point */
		bool    m_painteventpoint;
		std::vector<double>  m_RepulPoints;//i for X, i+1 for Y
		std::vector<double>  m_RepulContourPoints;//i for X, i+1 for Y
		
		/* Points for drawing */
		double  m_X;
		double  m_Y;
		int     m_num_samples;
		double  m_A, m_r;
		Point2* m_sample_pts;
		
		/* Widget */
		QWidget* m_parent;
		QPixmap  m_pixmap;
		QPixmap  m_savepixmap;
		QPixmap  m_graypixmap;
		QPixmap  m_zoompixmap;
		QPixmap  m_repulpixmap;
		QPixmap  m_visupixmap;
		
		bool m_debug;
};

#endif
