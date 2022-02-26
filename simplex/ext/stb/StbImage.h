//#####################################################################
// Stb image interface
// Interface functions for SLAX
//#####################################################################
#ifndef __StbImage_h__
#define __StbImage_h__
#include <string>

namespace Stb
{
int Write_Png(char const *filename,int w,int h,int comp,const void *data,int stride_in_bytes);
unsigned char* Read_Image_8(char const *filename,int *x,int *y,int *channels_in_file,int desired_channels);
unsigned short* Read_Image_16(char const *filename,int *x,int *y,int *channels_in_file,int desired_channels);
void Set_Flip_Image_Rows(const int flip);

template<class T_VAL> void Read_Image(const std::string& name,int& width,int& height,int& channels,T_VAL* & image);
};

#endif