//////////////////////////////////////////////////////////////////////////
// Opengl viewer driver
// Copyright (c) (2018-), Bo Zhu
// This file is part of SimpleX, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////
////This file contains customized viewers for different applications
//////////////////////////////////////////////////////////////////////////

#ifndef __OpenGLViewerDriver_h__
#define __OpenGLViewerDriver_h__
#include "OpenGLViewer.h"
class OpenGLTextArray3D;

//////////////////////////////////////////////////////////////////////////
////fluid viewers

class OpenGLViewerFluidEuler : public OpenGLViewer
{typedef OpenGLViewer Base;
public:
	virtual void Initialize();
	virtual void Initialize_Data();
};

class OpenGLViewerFluidEulerHighRes : public OpenGLViewer
{typedef OpenGLViewer Base;
public:
	virtual void Initialize_Data();
};

class OpenGLViewerFluidLagrangian : public OpenGLViewer
{typedef OpenGLViewer Base;
public:
	virtual void Initialize();
	virtual void Initialize_Data();
};

class OpenGLViewerFluidSPHBubble : public OpenGLViewer
{
	typedef OpenGLViewer Base;
public:
	virtual void Initialize();
	virtual void Initialize_Data();
};

class OpenGLViewerParticleALEFilm : public OpenGLViewer
{
	typedef OpenGLViewer Base;
public:
	virtual void Initialize();
	virtual void Initialize_Data();
};

class OpenGLViewerOrientedParticle : public OpenGLViewer
{
	typedef OpenGLViewer Base;
public:
	virtual void Initialize();
	virtual void Initialize_Data();
};

class OpenGLViewerPoisson : public OpenGLViewer
{typedef OpenGLViewer Base;
public:
	virtual void Initialize_Data();
};

class OpenGLViewerLevelSet : public OpenGLViewer
{typedef OpenGLViewer Base;
public:
	virtual void Initialize_Data();
};

class OpenGLViewerVortex : public OpenGLViewer
{typedef OpenGLViewer Base;
public:
	virtual void Initialize();
	virtual void Initialize_Data();
};

class OpenGLViewerFluidMesh : public OpenGLViewer
{typedef OpenGLViewer Base;
public:
	virtual void Initialize();
	virtual void Initialize_Data();
};

class OpenGLViewerFluidDEC : public OpenGLViewer
{
	typedef OpenGLViewer Base;
public:
	virtual void Initialize();
	virtual void Initialize_Data();
};

class OpenGLViewerElasticityDEC : public OpenGLViewer
{
	typedef OpenGLViewer Base;
public:
	virtual void Initialize();
	virtual void Initialize_Data();
};

class OpenGLViewerSWE : public OpenGLViewer
{typedef OpenGLViewer Base;
public:
	virtual void Initialize();
	virtual void Initialize_Data();
};

class OpenGLViewerMicroFluidic : public OpenGLViewer
{
	typedef OpenGLViewer Base;
public:
	virtual void Initialize();
	virtual void Initialize_Data();
};

//////////////////////////////////////////////////////////////////////////
////solid viewers

class OpenGLViewerSolid : public OpenGLViewer
{typedef OpenGLViewer Base;
public:
	virtual void Initialize_Data();
};

class OpenGLViewerDrone : public OpenGLViewer
{typedef OpenGLViewer Base;
public:
	virtual void Initialize_Data();
};

class OpenGLViewerFem : public OpenGLViewer
{typedef OpenGLViewer Base;
public:
	virtual void Initialize_Data();
};

class OpenGLViewerTopo : public OpenGLViewer
{typedef OpenGLViewer Base;
public:
	virtual void Initialize();
	virtual void Initialize_Data();
};

class OpenGLViewerCat : public OpenGLViewer
{typedef OpenGLViewer Base;
public:
	virtual void Initialize();
	virtual void Initialize_Data();
};

//////////////////////////////////////////////////////////////////////////
////geometry viewers

class OpenGLViewerGeometry : public OpenGLViewer
{typedef OpenGLViewer Base;
public:
	virtual void Initialize();
	virtual void Initialize_Data();
};

class OpenGLViewerMesh : public OpenGLViewer
{typedef OpenGLViewer Base;
public:
	OpenGLTextArray3D* opengl_text_vertex=nullptr;
	virtual void Toggle_Next_Frame();
	virtual void Initialize();
	virtual void Initialize_Data();
};

class OpenGLViewerPoints : public OpenGLViewer
{typedef OpenGLViewer Base;
public:
	virtual void Initialize();
	virtual void Initialize_Data();
};

class OpenGLViewerCurve : public OpenGLViewer
{typedef OpenGLViewer Base;
public:
	virtual void Initialize();
	virtual void Initialize_Data();
};

class OpenGLViewerVoronoi : public OpenGLViewer
{typedef OpenGLViewer Base;
public:
	virtual void Initialize();
	virtual void Initialize_Data();
};

class OpenGLViewerMatchStick : public OpenGLViewer
{typedef OpenGLViewer Base;
public:
	virtual void Initialize();
	virtual void Initialize_Data();
};

class OpenGLViewerImpulse : public OpenGLViewer
{typedef OpenGLViewer Base;
public:
	virtual void Initialize();
	virtual void Initialize_Data();
};

class OpenGLViewerImpulseDemo : public OpenGLViewer
{typedef OpenGLViewer Base;
public:
	virtual void Initialize();
	virtual void Initialize_Data();
};

#endif