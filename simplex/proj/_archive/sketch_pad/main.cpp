#include <iostream>
#include "SketchPadDriver.h"

#ifndef __Main_cpp__
#define __Main_cpp__

int main(int argc,char* argv[])
{
	int driver=1;

	switch(driver){
	case 1:{
		SketchPadDriver<2> driver;
		driver.Initialize();
		driver.Run();	
	}break;
	}
}

#endif
