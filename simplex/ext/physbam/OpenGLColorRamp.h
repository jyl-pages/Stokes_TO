//#####################################################################
// Copyright 2004, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license 
// contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// class OPENGL_COLOR_RAMP
//##################################################################### 
#ifndef __OPENGL_COLOR_RAMP__
#define __OPENGL_COLOR_RAMP__
#include "Common.h"
#include "OpenGLColor.h"

class OPENGL_COLOR_RAMP
{
public:
    Array<real> color_x;
    Array<OPENGL_COLOR> less_colors,equal_colors,greater_colors;

    OPENGL_COLOR_RAMP(){}
    ~OPENGL_COLOR_RAMP(){}

    void Add_Color(real x,const OPENGL_COLOR& color)
    {Add_Color(x,color,color,color);}

    void Add_Color(real x,const OPENGL_COLOR& color,const OPENGL_COLOR& exact_color)
    {Add_Color(x,color,exact_color,color);}

    OPENGL_COLOR Lookup(real x) const
    {
        int left_index=0,right_index=0;
        for(int i=0;i<(int)color_x.size();i++){
            if(x>color_x[i])left_index=i;
            else if(x<color_x[i]){right_index=i;break;}
            else return equal_colors[i];}
        if(left_index&&right_index){real alpha=(x-color_x[left_index])/(color_x[right_index]-color_x[left_index]);
            return real(alpha)*less_colors[right_index]+real(1.0-alpha)*greater_colors[left_index];}
        else if(left_index)return greater_colors[left_index];
        else if(right_index)return less_colors[right_index];
        return OPENGL_COLOR(1,0,0);
    }

    void Add_Color(real x,const OPENGL_COLOR& less_color,const OPENGL_COLOR& exact_color,const OPENGL_COLOR& greater_color)
    {
        assert(color_x.size()==0||x>color_x[(int)color_x.size()-1]);
        color_x.push_back(x);
        less_colors.push_back(less_color);
        equal_colors.push_back(exact_color);
        greater_colors.push_back(greater_color);
    }

    static OPENGL_COLOR_RAMP *Matlab_Jet(real value_min,real value_max)
    {
        OPENGL_COLOR_RAMP *jet=new OPENGL_COLOR_RAMP;
        real interval_width=value_max-value_min;
        jet->Add_Color(interval_width*0+value_min,OPENGL_COLOR(0,0,0.5608f));
        jet->Add_Color(interval_width*(real)0.1406+value_min,OPENGL_COLOR(0,0,1));
        jet->Add_Color(interval_width*(real)0.3594+value_min,OPENGL_COLOR(0,1,1));
        jet->Add_Color(interval_width*(real)0.6094+value_min,OPENGL_COLOR(1,1,0));
        jet->Add_Color(interval_width*(real)0.8594+value_min,OPENGL_COLOR(1,0,0));
        jet->Add_Color(interval_width*1+value_min,OPENGL_COLOR(0.5f,0,0));
        return jet;
    }

    static OPENGL_COLOR_RAMP *Matlab_Hot(real value_min,real value_max)
    {
        OPENGL_COLOR_RAMP *hot=new OPENGL_COLOR_RAMP;
        real interval_width=value_max-value_min;
        hot->Add_Color(interval_width*0+value_min,OPENGL_COLOR(0,0,0,0));
        hot->Add_Color(interval_width*(real)0.3750+value_min,OPENGL_COLOR(1,0,0));
        hot->Add_Color(interval_width*(real)0.7656+value_min,OPENGL_COLOR(1,1,0));
        hot->Add_Color(interval_width*1+value_min,OPENGL_COLOR(1,1,1));
        return hot;
    }

    static OPENGL_COLOR_RAMP *Two_Color_Ramp(real value_min,real value_max,const OPENGL_COLOR& color_min,const OPENGL_COLOR& color_min_exact,
        const OPENGL_COLOR& color_max,const OPENGL_COLOR& color_max_exact)
    {
        OPENGL_COLOR_RAMP *ramp=new OPENGL_COLOR_RAMP;
        ramp->Add_Color(value_min,color_min,color_min_exact);ramp->Add_Color(value_max,color_max,color_max_exact);
        return ramp;
    }

    static OPENGL_COLOR_RAMP *Two_Color_Ramp(real value_min,real value_max,const OPENGL_COLOR& color_min,const OPENGL_COLOR& color_max)
    {
        OPENGL_COLOR_RAMP *ramp=new OPENGL_COLOR_RAMP;
        ramp->Add_Color(value_min,color_min,color_min);ramp->Add_Color(value_max,color_max,color_max);
        return ramp;
    }
};   

typedef OPENGL_COLOR_RAMP OpenGLColorRamp;
#endif
