#ifndef _CCSEGTOOL_PARAMETERS_H_
#define _CCSEGTOOL_PARAMETERS_H_

//Specific librairies
#include "CCsegtool_initialization.h"

//Standard librairies
#include <math.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <vector>

//VTK librairies 
#include <vtkSmartPointer.h>
#include <vtkCellArray.h>
#include <vtkCellData.h>
#include <vtkDoubleArray.h>
#include <vtkPoints.h>
#include <vtkPolyLine.h>
#include <vtkPolyData.h>
#include <vtkPolyDataMapper.h>
#include <vtkActor.h>
#include <vtkRenderWindow.h>
#include <vtkRenderer.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkPolyDataWriter.h>

//Defines
#define MAX_STEPS 256
#define MAX_COEFS 150
#define MAX_PROFS MAX_STEPS
#define MAX_PROFLEN 51
#define MAX_MODES MAX_COEFS*4

#define COLMAJ 1
#define ROWMAJ 0

class CCsegtool_parameters
{
	//Methods
	public :
		//constructor
		CCsegtool_parameters(ExtractImageConstPointer _Image, std::string way, InputImageType::Pointer Image3D, 
				     bool permute_x_y, bool reflectXOn,bool reflectYOn, bool doubleOn, int sliceDir, 
				     int numSlice, int X, int Y,float Rotation, float Scale, bool FSXForm, 
				     double PSDistance, int WMIntensity,int MPSDisplacement, int OSteps, 
				     SizeType ImageSize, bool Debug=false);
		//to accede to the private attributs
		//Methods for Parameter Specifications
		void readf(std::string way, int nbfile);
		void readf4(std::string way, int nbfile);
		void getsizeof(int *i, int nbfile);
		void writeoutput(std::string path_output, std::string _nameofproject);
		void setparam(double X, double Y, double Scale, double Rot);
		void transform2DIndex_3Dindex(Point2* pts, int numpts, std::string path_output, std::string nameofproject,
					      double Xsize, double Ysize);
		
		//Methods for parameters
		void SetX(int X);
		int GetX();
		void SetY(int Y);
		int GetY();
		void SetRotation(float _Rotation);
		float GetRotation();
		void SetScale(float Scale);
		float GetScale();
		void SetFSXForm(bool FSXForm);
		bool GetFSXForm();
		void SetUnconstrained(bool Unconstrained);
		bool GetUnconstrained();
		void SetPSDistance(double Distance);
		double GetPSDistance();
		void SetWMIntensity(int Wmintensity);
		int GetWMIntensity();
		void SetMPSDisplacement(int Displacement);
		int GetMPSDisplacement();
		void SetOiter(int iter);
		int GetOiter();
		void SetResetMean(bool Resetmean);
		bool GetResetMean();
		
		//Methods for Input Port
		void SetImage(ExtractImageConstPointer _Image);
		ExtractImageConstPointer GetImage();
		std::vector<float> GetSMMean();
		float GetSMMeanval(int x,int y);
		std::vector<float> GetEigenvectors();
		float GetEigenvectorsval(int x,int y, int z);
		std::vector<float> GetBound();
		std::vector<float> GetPMMean();
		float GetPMMeanval(int x,int y);
		std::vector<float> GetSigmaInv();
		float GetSigmaInvval(int x,int y, int z);
		//Methods for Output Port
		void SetSSCoefs(std::vector<float> _SSCoefs);
		std::vector<float> GetSSCoefs();
		void SetIProfiles(std::vector<float> _IProfiles);
		std::vector<float> GetIProfiles();
		void SetGoodnessoffit(std::vector<double> _Goodnessoffit);
		std::vector<double> GetGoodnessoffit();
		void SetPPShifts(std::vector<float> _PPShifts);
		std::vector<float> GetPPShifts();
		void Seteigloads(std::vector<float> _eigloads);
		void Seteigloadsval(int i, float value);
		std::vector<float> Geteigloads();
		/*size in/out vector method*/
		int GetSMMeansize(int number);
		int GetEigenvectorssize(int number);
		int GetBoundsize();
		int GetPMMeansize(int number);
		int GetSigmaInvsize(int number);
		int GetNumPointsProba();
		int GetNumSlice();
		
		/* for coefs */
		void SetCoefs(Point4* coefs);
		Point4* GetCoefs();
		void SetNumcoefs(unsigned int num_coefs);
		unsigned int GetNumcoefs();
		bool Getcoefsgiven();
		vtkSmartPointer<vtkPoints> GetPointsvtk();
		ImagePointer GetImage3D();
		
		/* from initialization */
		void  SetGlobalvalue(short WMvalue, float voxelsizeX, float voxelsizeY);
		short GetWM();
		float GetVoxelSizeX();
		float GetVoxelSizeY();
		
	private :
		/*Parameters: */
		/* Pointer of Image */
		InputImageType::Pointer          m_3DImage;
		//OrientedInputImageType::Pointer  m_Oriented3DImage;
		ExtractImageConstPointer         m_Image;
		
		/* Input Port Specifications */
		std::vector<float> m_SMMean;
		std::vector<float> m_Eigenvectors;
		std::vector<float> m_Bound;
		std::vector<float> m_PMMean;
		std::vector<float> m_SigmaInv;
		Point4 *      m_coefs;
	
		/* Output Port Specifications */
		std::vector<float>  m_SSCoefs;
		std::vector<float>  m_IProfiles;
		std::vector<double> m_Goodnessoffit;
		std::vector<float>  m_PPShifts;
		std::vector<float>  m_eigloads;
	
	
		/* Parameter Specifications */
		int m_X;
		int m_Y;
		float m_Rotation;
		float m_Scale;
		bool m_FSXForm;
		bool m_Unconstrained;
		double m_PSDistance;
		int m_WMIntensity;
		int m_MPSDisplacement;
		int m_Oiter;
		bool m_ResetMean;
		int m_numpts;
		
		/* size of vector */
		int m_SMMean_xsize, m_SMMean_ysize;
		int m_Eigenvectors_xsize, m_Eigenvectors_ysize, m_Eigenvectors_zsize;
		int m_Bound_xsize;
		int m_PMMean_xsize, m_PMMean_ysize;
		int m_SigmaInv_xsize, m_SigmaInv_ysize, m_SigmaInv_zsize;
		
		/* Debug */
		bool         debug;
		bool         coefsgiven;
		unsigned int m_num_coefs;
		
		/* Parameters for the convertion 2Dindex to 3Dindex */
		bool        m_permute_x_y;
		bool        m_reflectXOn;
		bool        m_reflectYOn;
		bool        m_doubleOn;
		int         m_sliceDir;
		int         m_numSlice;
		SizeType    m_ImageSize;
		
		/* vtk/itk */
		vtkSmartPointer<vtkPoints> m_points;
		ImagePointer m_Image3D;
		
		/* From initialization */
		short m_WMvalue;
		float m_voxelsizeX;
		float m_voxelsizeY;
};

#endif
