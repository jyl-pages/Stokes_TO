//#####################################################################
// Stb image interface
// Bo Zhu
//#####################################################################
#include "StbImage.h"

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

namespace Stb{
	void Set_Flip_Image_Rows(const int flip){stbi_set_flip_vertically_on_load(flip);}

	int Write_Png(char const *filename,int w,int h,int comp,const void *data,int stride_in_bytes)
	{return stbi_write_png(filename,w,h,comp,data,stride_in_bytes);}

	unsigned char* Read_Image_8(char const *filename, int *x, int *y, int *channels_in_file, int desired_channels)
	{return stbi_load(filename,x,y,channels_in_file,desired_channels);}

	unsigned short* Read_Image_16(char const *filename, int *x, int *y, int *channels_in_file, int desired_channels)
	{
		return stbi_load_16(filename, x, y, channels_in_file, desired_channels);
	}

	template<class T_VAL> void Read_Image(const std::string& name,int& width,int& height,int& channels,T_VAL* & image){}
	template void Read_Image(const std::string&,int& width,int& height,int& channels,float* &);

	template<> void Read_Image<unsigned short>(const std::string& file_name,int& width,int& height,int& channels,unsigned short* & image)
	{image=Read_Image_16(file_name.c_str(),&width,&height,&channels,0);}
	template void Read_Image(const std::string&,int& width,int& height,int& channels,unsigned short* &);
	template<> void Read_Image<unsigned char>(const std::string& file_name,int& width,int& height,int& channels,unsigned char* & image)
	{image=Read_Image_8(file_name.c_str(),&width,&height,&channels,0);}
	template void Read_Image(const std::string&,int& width,int& height,int& channels,unsigned char* &);
};

