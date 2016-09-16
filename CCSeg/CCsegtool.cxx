/**
  Segmentation of Corpus Callosum, main file, parsing input etc, B. YVERNAULT

test command line :
    With default WMintensity/MPSDisplacement/Number of iteration and compute probability model
bin/CCsegtool -I /biomed-resimg/Autism2/IBIS2/people/yben/5002-003-01_10_T1.hdr -M /biomed-resimg/Autism2/IBIS2/people/yben/5002-003-01_10_t2_fit_hard.hdr -A /tools/CCsegmenter_autism/Model/ -O /home/yben/test_output/ -N test --sliceDir 0 --ReflectX --double --ComputeProbaModel

    With other WMintensity/MPSDisplacement/iteration values choose in the line command
bin/CCsegtool -I /biomed-resimg/Autism2/IBIS2/people/yben/5002-003-01_10_T1.hdr -M /biomed-resimg/Autism2/IBIS2/people/yben/5002-003-01_10_t2_fit_hard.hdr -A /tools/CCsegmenter_autism/Model/ -O /home/yben/test_output/ -N test --sliceDir 0 --ReflectX --double -w 110,90,90 -i 50,5,3 -d 10,3,2 

    With the Graphical Interface
bin/CCsegtool --gui
 **/

#include "CCsegtoolCLP.h"
#include "CCsegtoolguiwindow.h"
#include <QApplication>

int main(int argc, char *argv[])
{
	PARSE_ARGS;
	QApplication app(argc, argv);
	
	/* With the GUI */
	if(gui)
	{
		if(debug)
			std::cout<<"CCSEGTOOL WITH GUI"<<std::endl;
		/* Set and show the window */
		CCsegtoolguiwindow CCsegtoolwindow(scalefactor,debug);
		CCsegtoolwindow.show();
		return app.exec();
	}
	/* Without the otpion GUI */
	else
	{
		if(debug)
			std::cout<<"CCSEGTOOL COMMAND LINE"<<std::endl;
		
		/*Name of the project */
		if(nameofproject.compare("")==0)
		{
			nameofproject = Image_filename.substr(Image_filename.find_last_of("/")+1,
					Image_filename.find_last_of(".")-Image_filename.find_last_of("/")-1);
		}
		
		/* Create the parameters for the command line */
		CCsegtool_parameters *parameters;
		Imageview * ViewImage=NULL;
		/* X/Y/Scale/Rotation */
		double lastX, lastY, lastScale, lastRot;
		
		/* Execute CCsegtool */
		/* Preview */
		if(debug)
			std::cout<<"INITIALIZE"<<std::endl;
		if (!initialize(Image_filename, Mask_filename, CCAtlasDirectory, OuputFolder, nameofproject,
				interpolationlinear, vesselRemoveOn, segLabel, averageNum, permute_x_y, reflectXOn,
				reflectYOn, openOn, doubleOn, sliceDir, MidPlaneSliceNumber, Number_Pts,
				FSXForm, PSDistance, Unconstrained, WMintensity, MPSDisplacement, Number_iteration,
				Lambdamax, Coefofoptim, debug, false, ViewImage, parameters, &lastX, &lastY, &lastScale,
				&lastRot,scalefactor, othercompo,addangle))
		  return 0;
		
		/* Run */
		if(debug)
			std::cout<<"RUN/CONTINUE"<<std::endl;
		compute(OuputFolder, nameofproject,Number_Pts,Unconstrained, WMintensity, MPSDisplacement, Number_iteration,
			Lambdamax, Coefofoptim, debug, false, ViewImage, parameters, &lastX, &lastY, &lastScale,
			&lastRot);
		
		if(loop != 0)
		{
			/* Continue */
			for(int i=0;i<loop;i++)
			{
				parameters->setparam(lastX, lastY, lastScale, lastRot);
				compute(OuputFolder, nameofproject,Number_Pts, Unconstrained, WMintensity, MPSDisplacement,
					Number_iteration,Lambdamax, Coefofoptim, debug, false, ViewImage, parameters, 
					&lastX, &lastY, &lastScale, &lastRot);
			}
		}
		
		/* Compute probability model if call */
		if(ComputeProbaModel)
		{
			// Create the Probacurve and set the parameters
			CCcurveProba *Probacurve;
			Probacurve = new CCcurveProba(parameters->GetImage3D(), parameters->GetPointsvtk(),
					CCAtlasDirectory, OuputFolder, nameofproject, vesselRemoveOn, segLabel,
					averageNum, permute_x_y, reflectXOn, reflectYOn, openOn, doubleOn, sliceDir,
					parameters->GetNumSlice(), parameters->GetNumPointsProba(), debug);
			
			/* Compute Proba Model */
			Probacurve->compute_proba();
		}
		return 0;
	}
	
	if(debug)
		std::cout<<"END OF CCSEGTOOL"<<std::endl;
}

