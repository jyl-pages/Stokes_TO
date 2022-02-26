//////////////////////////////////////////////////////////////////////////
// Opengl viewer
// Copyright (c) (2018-), Bo Zhu
// This file is part of SimpleX, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#include "ParseArgs.h"
#include "AuxFunc.h"
#include "OpenGLViewerDriver.h"
#include "OpenGLViewerConfigable.h"
#include "OpenGLViewerTestCase.h"
#include "OpenGLViewerInteractive.h"

int main(int argc,char* argv[])
{
    ////parse arguments
    ParseArgs parse_args;
	parse_args.Add_String_Argument("-o", "D:/Codes/particle-ale-exams/data/heat11-fine", "data path");
	parse_args.Add_String_Argument("-oo","","offscreen rendering output");
	parse_args.Add_String_Argument("-m","particle_ale");
    parse_args.Add_Integer_Argument("-d",3);
    parse_args.Add_Double_Argument("-vmin",0);
    parse_args.Add_Double_Argument("-vmax",1);
	parse_args.Add_Option_Argument("-anim");
	parse_args.Add_Double_Argument("-scale",1.);
	parse_args.Add_String_Argument("-cf","opengl_viewer_config");	////in the default path
	parse_args.Add_String_Argument("-cff","");						////in a customized path
	parse_args.Add_String_Argument("-cmd","");
	parse_args.Add_Option_Argument("-min","minimal mode");
	parse_args.Add_Integer_Argument("-test",1,"test");

    parse_args.Parse(argc,argv);

	std::string mode=parse_args.Get_String_Value("-m");
	std::string offscreen_output=parse_args.Get_String_Value("-oo");
	if(offscreen_output=="")offscreen_output=parse_args.Get_String_Value("-o")+"/_images";

	std::string cmd=parse_args.Get_String_Value("-cmd");

	std::shared_ptr<OpenGLViewer> viewer=nullptr;
	if(mode=="config"){viewer.reset(new OpenGLViewerConfiguable());}
	else if(mode=="test"){viewer.reset(new OpenGLViewerTestCase());}
	else if(mode=="int"){viewer.reset(new OpenGLViewerInteractive());}

	////fluid viewers
	else if(mode=="fluid"){viewer.reset(new OpenGLViewerFluidEuler());}
	else if(mode=="fluid_hr"){viewer.reset(new OpenGLViewerFluidEulerHighRes());}
	else if(mode=="fluid_lag"){viewer.reset(new OpenGLViewerFluidLagrangian());}
	else if(mode=="particles") {viewer.reset(new OpenGLViewerOrientedParticle());}
	else if(mode=="fluid_mesh"){viewer.reset(new OpenGLViewerFluidMesh());}
	else if(mode=="fluid_vor"){viewer.reset(new OpenGLViewerVortex());}
	else if(mode=="poisson"){viewer.reset(new OpenGLViewerPoisson());}
	else if(mode=="levelset"){viewer.reset(new OpenGLViewerLevelSet());}
	else if(mode=="ms"){viewer.reset(new OpenGLViewerMatchStick());}
	else if(mode=="imp"){viewer.reset(new OpenGLViewerImpulse());}
	else if(mode=="imp_demo"){viewer.reset(new OpenGLViewerImpulseDemo());}
	else if(mode=="sph_bubble"){ viewer.reset(new OpenGLViewerFluidSPHBubble()); }
	else if (mode == "particle_ale") { viewer.reset(new OpenGLViewerParticleALEFilm()); }
	else if (mode == "fluid_dec") { viewer.reset(new OpenGLViewerFluidDEC()); }
	else if (mode == "swe") { viewer.reset(new OpenGLViewerSWE()); }
	else if (mode == "microfluidic") { viewer.reset(new OpenGLViewerMicroFluidic()); }

	////solid viewers
	else if(mode=="solid"){viewer.reset(new OpenGLViewerSolid());}
	else if(mode=="drone"){viewer.reset(new OpenGLViewerDrone());}
	else if(mode=="fem"){viewer.reset(new OpenGLViewerFem());}
	else if(mode=="topo"){viewer.reset(new OpenGLViewerTopo());}
	else if(mode=="cat"){viewer.reset(new OpenGLViewerCat());}
	else if (mode == "elasticity_dec") { viewer.reset(new OpenGLViewerElasticityDEC()); }

	////geometry viewers
	else if(mode=="mesh"){viewer.reset(new OpenGLViewerMesh());}
	else if(mode=="geo"){viewer.reset(new OpenGLViewerGeometry());}
	else if(mode=="point"){viewer.reset(new OpenGLViewerPoints());}
	else if(mode=="curve"){viewer.reset(new OpenGLViewerCurve());}  
	else if(mode=="voronoi"){viewer.reset(new OpenGLViewerVoronoi());}

	else{std::cout<<"Invalid viewer mode"<<std::endl;return 0;}

	//////////////////////////////////////////////////////////////////////////
	////These parameters need to be set before initialization
	viewer->output_dir=parse_args.Get_String_Value("-o");
	int dim=parse_args.Get_Integer_Value("-d");
	viewer->use_2d_display=(dim==2);
	
	std::string config_file_local=parse_args.Get_String_Value("-cff");
	if(config_file_local!="")viewer->config_file_name=config_file_local;
	else viewer->config_file_name=Path::Script()+"/"+parse_args.Get_String_Value("-cf")+".json";;

	viewer->default_scale=parse_args.Get_Double_Value("-scale");
	viewer->minimal_mode=parse_args.Get_Option_Value("-min");

	viewer->test=parse_args.Get_Integer_Value("-test");

	//////////////////////////////////////////////////////////////////////////

	viewer->Initialize();
	viewer->Set_Offscreen_Output_Dir(offscreen_output); // Just set where to output rendered images
	
	if(cmd!="")viewer->Toggle_Command(cmd);

	viewer->Run(); // It equals to opengl_window->Run()
	viewer->Finish();

	return 0;
}