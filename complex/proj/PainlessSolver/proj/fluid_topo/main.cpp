//#####################################################################
// Main
// Copyright (c) (2018-), Bo Zhu, boolzhu@gmail.com
// This file is part of SLAX, whose distribution is governed by the LICENSE file.
//#####################################################################
#include <iostream>
#include "ParseArgs.h"
#include "FluidTopoOptDriver.h"
#include "FluidTopoOptDriver3D.h"

#ifndef __Main_cpp__
#define __Main_cpp__

int main(int argc,char* argv[])
{
    const int d=2;

    ParseArgs parse_args;
    parse_args.Add_String_Argument("-o","output","output path");
    parse_args.Add_Integer_Argument("-s",64,"resolution");
    parse_args.Add_Integer_Argument("-test",2,"test");
	parse_args.Add_Integer_Argument("-driver",1,"driver");
	parse_args.Add_Integer_Argument("-lf",200,"last frame");
	parse_args.Add_Double_Argument("-frac", 200, "frac");
    parse_args.Parse(argc,argv);

    std::string output_dir=parse_args.Get_String_Value("-o");
    const int scale=parse_args.Get_Integer_Value("-s");
	const int driver=parse_args.Get_Integer_Value("-driver");
	const int test=parse_args.Get_Integer_Value("-test");
	const int last_frame=parse_args.Get_Integer_Value("-lf");
	const double frac= parse_args.Get_Double_Value("-frac");

	Eigen::setNbThreads(4);

	if (driver == 0)
	{
		FluidTopoOptDriver testDriver;
		testDriver.grid_size = scale;
		testDriver.hole_size = scale / 8;
		testDriver.test = test;
		testDriver.last_frame = last_frame;
		testDriver.frac = frac;
		testDriver.Initialize();
		testDriver.Run();
	}
	if (driver == 1)
	{
		FluidTopoOptDriver3D testDriver;
		testDriver.grid_size = scale;
		testDriver.hole_size = scale / 8;
		testDriver.test = test;
		testDriver.last_frame = last_frame;
		testDriver.frac = frac;
		testDriver.Initialize();
		testDriver.Run();
	}
	return 0;
}

#endif