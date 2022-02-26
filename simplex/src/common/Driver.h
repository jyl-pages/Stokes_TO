//////////////////////////////////////////////////////////////////////////
// Driver template, general temporal evolution
// Copyright (c) (2018-), Bo Zhu
// This file is part of SimpleX, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#ifndef __Driver_h__
#define __Driver_h__
#include "Common.h"
#include "File.h"
#include "Timer.h"

class Driver
{
public:
	int test=1;
	std::string output_dir="output";
	std::string frame_dir;
	int first_frame=0,last_frame=200,current_frame=0;
	int current_step = 0;
	real frame_rate=50;
	real time=(real)0,current_time=(real)0;
	real cfl=(real)1;
	int max_iter_per_frame = 2000;//set to -1 to disable max iteration limit
	int scale=1;
	bool verbose=false;
	int snapshot_stride = 2;
	Timer<real> timer;

	real Time_At_Frame(const int frame){return (real)frame/frame_rate;}
	int Frame_At_Time(const real time){return (int)((real)time*frame_rate);}
	virtual real CFL() const {return cfl;}

	virtual void Initialize(){}

	virtual void Run()
	{
		if (current_frame == 0) Write_Output_Files(current_frame);
		while (current_frame < last_frame) {
			current_frame++;
			Advance_To_Target_Time(Time_At_Frame(current_frame));
			Write_Output_Files(current_frame);
		}
	}

	virtual void Advance_To_Target_Time(const real target_time)
	{
		bool done=false;for(int substep=1;!done;substep++){
			real dt=CFL();
			if(time+dt>=target_time){dt=target_time-time;done=true;}
			else if(time+2*dt>=target_time){dt=(real).5*(target_time-time);}
			Advance_One_Time_Step(dt,time);
			time+=dt;}
	}

	virtual void Advance_One_Time_Step(const real dt,const real time){}

	virtual void Write_Output_Files(const int frame)
	{	
		if(frame==0){
			if(!File::Directory_Exists(output_dir.c_str()))
				File::Create_Directory(output_dir);}

		frame_dir=output_dir+"/"+std::to_string(frame);
		if(!File::Directory_Exists(frame_dir.c_str()))File::Create_Directory(frame_dir);
		
		{std::string file_name=output_dir+"/0/last_frame.txt";
		File::Write_Text_To_File(file_name,std::to_string(frame));}

		if(verbose)std::cout<<"Write output files for frame "<<frame<<std::endl;
	}

	//////////////////////////////////////////////////////////////////////////
	////streaming the file IO
	virtual void Run_Stream()
	{
		while(current_frame<last_frame){
			#pragma omp parallel sections
			{
				#pragma omp section
				{Write_Output_Files(current_frame);}

				#pragma omp section
				{Advance_To_Target_Time(Time_At_Frame(current_frame+1));}
			}
			current_frame++;}	
		Write_Output_Files(current_frame);
	}

	virtual void Update_IO_Buffer_Stream(const int frame){}

	virtual void Advance_One_Time_Step_Stream(const real dt,const real time){}

	virtual void Write_Output_Files_Stream(const int frame){Write_Output_Files(frame);}
};

#endif