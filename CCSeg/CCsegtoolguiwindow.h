#ifndef _CCSEGTOOLGUIWINDOW_H_
#define _CCSEGTOOLGUIWINDOW_H_

#include "ui_CCsegtoolguiwindow.h"
#include "Globalfunc.h"
#include <QtGui>
#include <QObject>
#include <QString>
#include <QChar>
#include <QMainWindow>
#include <QFileDialog>
#include <QLabel>
#include <QImage>
#include <QMessageBox>
#include <QAction>
#include <QScrollArea>
#include <QScrollBar>
#include <QMouseEvent>
#include <QSlider>

class CCsegtoolguiwindow : public QMainWindow, public Ui::MainWindow 
{
	Q_OBJECT
	public:
		CCsegtoolguiwindow(float scalefactor=1.35, bool Debug = false, QWidget * parent = 0, Qt::WFlags f = 0 );
		void setvalue();
		void showparameters(double lastX, double lastY, double lastScale, double lastRot);
		void guicompute(int run_preview_continue);
		
	private slots:
		void InitPreview();
		void ConPreview();
		void Preview(int resetVisu);
		void running();
		void ComputeProbabilityModel();
		void setDebug();
		void defaultparaminput();
		void defaultparamcomputation();
		void browserImage();
		void browserMask();
		void browserCCAtlasDirectory();
		void browserOutput();
		void help();
		void newparam();
		void saveparam();
		void openparam();
		void saveimage();
		void update();
		void updateVisuImage();
		void resetpts();
		void setisoline();
		void changeisoline();
		void changecoef();
		
	private:
		/* Filter for event just on m_viewimage */
		bool eventFilter(QObject *object, QEvent *event);
		void mouseMoveEvent(QMouseEvent* event);
		void mouseWheelEvent(QWheelEvent *event);
		
		/* Zoom functions */
		void zoomOut();
		void zoomIn();
		void scaleImage(double factor);
		void adjustScrollBar(QScrollBar *scrollBar, double factor);
		void normalSize();
		double m_scaleFactor;
		void fitToWindow();
		
		/*Parameters for input*/
		std::string m_InputImage;
		std::string m_OldInputImage;
		std::string m_MaskImage;
		std::string m_CCAtlasDirectory;
		std::string m_outputfolder;
		std::string m_nameofproject;
		QString     m_oldfilename;	
		
		/*Parameters for computation*/
		bool        m_interpolationlinear;
		bool        m_vesselremove;
		bool        m_permute_x_y;
		bool        m_old_permute_x_y;
		bool        m_reflectX;
		bool        m_oldreflectX;
		bool        m_reflectY;
		bool        m_oldreflectY;
		bool        m_opening;
		bool        m_double;
		bool        m_unconstrained;
		bool        m_seglabel;
		int         m_slicedir;
		std::string m_MidPlaneSliceNumber;
		int         m_averagenum;
		bool        m_FSXForm;
		float       m_PSDistance;
		int         m_number_pts;
		bool        m_Debug;
		float       m_scalefactor;
		
		/*Parameters useful*/
		std::string m_MPSDisplacement;
		std::string m_loop;
		std::string m_WMvalue;
		int         m_Lambdamax;
		std::string m_Coefofoptim;
		std::string m_A;
		double      m_r;
		
		/* Repulsive Zone variable */
		int    m_isoline;
		
		/* Voxelsize */
		float  m_VoxelSizeX, m_VoxelSizeY;
		
		/* X/Y/Scale/Rotation */
		double m_lastX, m_lastY, m_lastScale, m_lastRot;
		
		/* Widget help */
		QMessageBox *m_information;
		
		/* Output image */
		QPixmap *m_outputimage;
		
		/* Drawing Widget */
		Imageview *m_viewimage;
		bool       m_ViewImageExist;
		
		/* Widget for zoom */
		bool         m_Scrollexist;
		QScrollArea* m_scrollArea;
		
		/* Parameters used for continue*/
		CCsegtool_parameters *m_last_parameters;
		
		/* Correct componant or angle */
		int  m_angle; //0 for 0, 1 for 90, 2 for 180, 3 for 270 degrees
		bool m_othercomponant; //change the composant during the initialization
		
		/* Visualisation of the image */
		bool m_rot90;
		bool m_rot180;
		bool m_rot270;
		bool m_noRotation;
		
		/* Run called */
		bool m_runcalled;
};

#endif
