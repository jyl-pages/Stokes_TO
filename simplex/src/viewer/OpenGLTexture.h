//////////////////////////////////////////////////////////////////////////
// Opengl texture
// Copyright (c) (2018-), Bo Zhu
// This file is part of SimpleX, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#ifndef __OpenGLTexture_h__
#define __OpenGLTexture_h__
#include <string>
#include <fstream>
#include <iostream>
#include <type_traits>
#include "glm.hpp"
#include "StbImage.h"
#include "Common.h"
#include "Hashtable.h"
#include "AuxFunc.h"
#include "Field.h"
#include "OpenGLCommon.h"

namespace OpenGLTextures{
//template<class T_STORAGE> GLenum GL_Storage_Format(T_STORAGE* p=nullptr){return GL_BYTE;}
inline GLenum GL_Storage_Format(ushort* p=nullptr){return GL_UNSIGNED_SHORT;}
inline GLenum GL_Storage_Format(uchar* p=nullptr){return GL_UNSIGNED_BYTE;}
inline GLenum GL_Storage_Format(GLfloat* p=nullptr){return GL_FLOAT;}
inline GLenum GL_Element_Format_Helper(std::integral_constant<int,1>* p=nullptr){return GL_RED;}
inline GLenum GL_Element_Format_Helper(std::integral_constant<int,2>* p=nullptr){return GL_RG;}
inline GLenum GL_Element_Format_Helper(std::integral_constant<int,3>* p=nullptr){return GL_RGB;}
inline GLenum GL_Element_Format_Helper(std::integral_constant<int,4>* p=nullptr){return GL_RGBA;}

class OpenGLTexture
{public:
	GLuint buffer_index=0;
	std::string name="";
	GLsizei width=0;
	GLsizei height=0;
	int channels=0;
	OpenGLTexture(){}
	virtual void Bind(GLuint idx=0){glActiveTexture(GL_TEXTURE0+idx);glBindTexture(GL_TEXTURE_2D,buffer_index);}
};

template<class T_STORAGE> class OpenGLTextureInstance : public OpenGLTexture
{typedef OpenGLTexture Base;
public:
	using Base::name;using Base::buffer_index;
	T_STORAGE* image=nullptr;

	OpenGLTextureInstance():Base(){}
	virtual void Initialize(const std::string& _prefix,const std::string& _type)
	{
		name=_prefix+"."+_type;std::string file_name=Path::Data()+"/textures/"+name;
		Stb::Set_Flip_Image_Rows(1);	////flip image rows
		Stb::Read_Image(file_name,width,height,channels,image);
		glGenTextures(1,&buffer_index);
		glBindTexture(GL_TEXTURE_2D,buffer_index);
		glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_WRAP_S,GL_REPEAT);
		glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_WRAP_T,GL_REPEAT);
		glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MIN_FILTER,GL_LINEAR);
		glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MAG_FILTER,GL_LINEAR);
		glTexImage2D(GL_TEXTURE_2D,0,GL_RGB,width,height,0,GL_RGB,GL_Storage_Format(image),image);
		delete image;
		glBindTexture(GL_TEXTURE_2D,0);
		std::cout<<"[OpenGLTextureInstance] Initialize texture 2d: "<<file_name<<", w: "<<width<<", h: "<<height<<", c: "<<channels<<", buf_idx: "<<buffer_index<<std::endl;
	}
};

template<class T_STORAGE> class OpenGLCubeMap : public OpenGLTextureInstance<T_STORAGE>
{typedef OpenGLTextureInstance<T_STORAGE> Base;
public:
	using Base::name;using Base::buffer_index;using Base::width;using Base::height;using Base::channels;using Base::image;
	OpenGLCubeMap():Base(){}
	virtual void Initialize(const std::string& _prefix,const std::string& _type)
	{
		name=_prefix+"."+_type;
		glGenTextures(1,&buffer_index);
		glBindTexture(GL_TEXTURE_CUBE_MAP,buffer_index);
		for(int i=0;i<6;i++){	////0-right,1-left,2-top,3-bottom,4-back,5-front
			std::string file_name=Path::Data()+"/textures/"+_prefix+"_"+std::to_string(i)+"."+_type;
			Stb::Read_Image(file_name,width,height,channels,image);
			glTexImage2D(GL_TEXTURE_CUBE_MAP_POSITIVE_X+i,0,GL_RGB,width,height,0,GL_RGB,GL_Storage_Format(image),image);
			std::cout<<"[OpenGLCubeMap] Initialize texture cubemap: "<<file_name<<", w: "<<width<<", h: "<<height<<", c: "<<channels<<", buf_idx: "<<buffer_index<<std::endl;
			delete image;}
		glTexParameteri(GL_TEXTURE_CUBE_MAP,GL_TEXTURE_MAG_FILTER,GL_LINEAR);
		glTexParameteri(GL_TEXTURE_CUBE_MAP,GL_TEXTURE_MIN_FILTER,GL_LINEAR);
		glTexParameteri(GL_TEXTURE_CUBE_MAP,GL_TEXTURE_WRAP_S,GL_CLAMP_TO_EDGE);
		glTexParameteri(GL_TEXTURE_CUBE_MAP,GL_TEXTURE_WRAP_T,GL_CLAMP_TO_EDGE);
		glTexParameteri(GL_TEXTURE_CUBE_MAP,GL_TEXTURE_WRAP_R,GL_CLAMP_TO_EDGE);
		glBindTexture(GL_TEXTURE_CUBE_MAP,0);
	}

//	virtual void Bind(GLuint idx=0){glActiveTexture(GL_TEXTURE0+idx);glBindTexture(GL_TEXTURE_CUBE_MAP,buffer_index);}	////TOFLX
};

template<class T_STORAGE=GLfloat> class OpenGLTextureVolume : public OpenGLTextureInstance<T_STORAGE>
{typedef OpenGLTextureInstance<T_STORAGE> Base;
public:
	using Base::name;using Base::buffer_index;
	OpenGLTextureVolume():Base(){}

	virtual void Initialize(const std::string& _name,const Field<T_STORAGE,3>* field=nullptr)
	{
		name=_name;
		glGenTextures(1,&buffer_index);
		Update_Data_To_Render(field);
	}

	virtual void Bind(GLuint idx=0){glActiveTexture(GL_TEXTURE0+idx);glBindTexture(GL_TEXTURE_3D,buffer_index);}

	void Update_Data_To_Render(const Field<T_STORAGE,3>* field=nullptr)
	{
		T_STORAGE* p=nullptr;
		if(field==nullptr)return;
		glBindTexture(GL_TEXTURE_3D,buffer_index);
		GLsizei w=(GLsizei)(*field).counts[0];GLsizei h=(GLsizei)(*field).counts[1];GLsizei dp=(GLsizei)(*field).counts[2];
		glTexImage3D(GL_TEXTURE_3D,0,GL_RED,w,h,dp,0,GL_RED,GL_Storage_Format(p),&field->array[0]);
		glTexParameteri(GL_TEXTURE_3D,GL_TEXTURE_MAG_FILTER,GL_LINEAR);
		glTexParameteri(GL_TEXTURE_3D,GL_TEXTURE_MIN_FILTER,GL_LINEAR);
		//glTexParameteri(GL_TEXTURE_3D,GL_TEXTURE_WRAP_S,GL_CLAMP_TO_EDGE);
		//glTexParameteri(GL_TEXTURE_3D,GL_TEXTURE_WRAP_T,GL_CLAMP_TO_EDGE);
		//glTexParameteri(GL_TEXTURE_3D,GL_TEXTURE_WRAP_R,GL_CLAMP_TO_EDGE);	

		float border_color[]={.0f};
		glTexParameterfv(GL_TEXTURE_3D,GL_TEXTURE_BORDER_COLOR,border_color);
		glTexParameteri(GL_TEXTURE_3D,GL_TEXTURE_WRAP_S,GL_CLAMP_TO_BORDER);
		glTexParameteri(GL_TEXTURE_3D,GL_TEXTURE_WRAP_T,GL_CLAMP_TO_BORDER);
		glTexParameteri(GL_TEXTURE_3D,GL_TEXTURE_WRAP_R,GL_CLAMP_TO_BORDER);

		glBindTexture(GL_TEXTURE_3D,0);
	}
};

template<class T_STORAGE=ushort,int d=4> class OpenGLTextureVectorVolume : public OpenGLTextureInstance<T_STORAGE>
{typedef OpenGLTextureInstance<T_STORAGE> Base;
public:
	using Base::name;using Base::buffer_index;
	OpenGLTextureVectorVolume():Base(){}

	virtual void Initialize(const std::string& _name,const Field<Vector<T_STORAGE,d>,3>* field=nullptr)
	{
		name=_name;
		glGenTextures(1,&buffer_index);
		Update_Data_To_Render(field);
	}

	virtual void Bind(GLuint idx=0){glActiveTexture(GL_TEXTURE0+idx);glBindTexture(GL_TEXTURE_3D,buffer_index);}

	void Update_Data_To_Render(const Field<Vector<T_STORAGE,d>,3>* field=nullptr)
	{
		T_STORAGE* p=nullptr;
		if(field==nullptr)return;
		glBindTexture(GL_TEXTURE_3D,buffer_index);
		GLsizei w=(GLsizei)(*field).counts[0];GLsizei h=(GLsizei)(*field).counts[1];GLsizei dp=(GLsizei)(*field).counts[2];
		Array<T_STORAGE> img(w*h*dp*d);AuxFunc::Fill(img,(T_STORAGE)0);for(auto i=0;i<field->array.size();i++)for(int j=0;j<d;j++)img[i*d+j]=field->array[i][j];
		glTexImage3D(GL_TEXTURE_3D,0,GL_Element_Format<d>(),w,h,dp,0,GL_Element_Format<d>(),GL_Storage_Format(p),&img[0]);
		glTexParameteri(GL_TEXTURE_3D,GL_TEXTURE_MAG_FILTER,GL_LINEAR);
		glTexParameteri(GL_TEXTURE_3D,GL_TEXTURE_MIN_FILTER,GL_LINEAR);
		glTexParameteri(GL_TEXTURE_3D,GL_TEXTURE_WRAP_S,GL_CLAMP_TO_EDGE);
		glTexParameteri(GL_TEXTURE_3D,GL_TEXTURE_WRAP_T,GL_CLAMP_TO_EDGE);
		glTexParameteri(GL_TEXTURE_3D,GL_TEXTURE_WRAP_R,GL_CLAMP_TO_EDGE);	
		glBindTexture(GL_TEXTURE_3D,0);
	}

	template<int dim> GLenum GL_Element_Format() const {std::integral_constant<int,dim>* p=nullptr;return GL_Element_Format_Helper(p);}	////d=4
};

template<class T_STORAGE> class OpenGLTexture1D : public OpenGLTextureInstance<T_STORAGE>
{typedef OpenGLTextureInstance<T_STORAGE> Base;
public:
	using Base::name;using Base::buffer_index;
	OpenGLTexture1D():Base(){}

	virtual void Initialize(const std::string& _name,const Array<T_STORAGE>* array=nullptr)
	{
		name=_name;
		glGenTextures(1,&buffer_index);
		Update_Data_To_Render(array);
	}

	virtual void Bind(GLuint idx=0){glActiveTexture(GL_TEXTURE0+idx);glBindTexture(GL_TEXTURE_1D,buffer_index);}

	void Update_Data_To_Render(const Array<T_STORAGE>* array=nullptr)
	{
		T_STORAGE* p=nullptr;
		if(array==nullptr)return;
		glBindTexture(GL_TEXTURE_1D,buffer_index);
		GLsizei w=(GLsizei)(*array).size();
		glTexImage1D(GL_TEXTURE_1D,0,GL_RGBA,w/4,0,GL_RGBA,GL_Storage_Format(p),&(*array)[0]);
		glTexParameteri(GL_TEXTURE_1D,GL_TEXTURE_MAG_FILTER,GL_LINEAR);
		glTexParameteri(GL_TEXTURE_1D,GL_TEXTURE_MIN_FILTER,GL_LINEAR);
		glTexParameteri(GL_TEXTURE_1D,GL_TEXTURE_WRAP_S,GL_CLAMP_TO_EDGE);
		glBindTexture(GL_TEXTURE_1D,0);
	}
};

class Texture_Library
{public:
	static Texture_Library* Instance(){static Texture_Library instance;return &instance;}
	std::shared_ptr<OpenGLTexture> Get(const std::string& name)
	{
		auto search=texture_hashtable.find(name);
		if(search!=texture_hashtable.end())return search->second;
		else return Lazy_Initialize_Texture(name);
	}

	std::shared_ptr<OpenGLTexture> Get(const std::string& name,const TextureType type)
	{
		auto search=texture_hashtable.find(name);
		if(search!=texture_hashtable.end())return search->second;
		else return Lazy_Initialize_Texture(name,type);
	}

protected:
	Hashtable<std::string,std::shared_ptr<OpenGLTexture> > texture_hashtable;

	Texture_Library(){}

	template<class T_TEX> std::shared_ptr<OpenGLTexture> Initialize_Texture(const std::string& name,T_TEX* p=nullptr)
	{
		T_TEX* tex=new T_TEX();tex->Initialize(name);
		std::shared_ptr<OpenGLTexture> ptr=std::make_shared<OpenGLTexture>();ptr.reset(tex);
		texture_hashtable.insert(std::make_pair(tex->name,ptr));return ptr;
	}

	bool Parse_Name(const std::string& file_name,std::string& prefix,std::string& ext)
	{size_t p=file_name.rfind('.');if(p==std::string::npos)return false;
	prefix=file_name.substr(0,p);ext=file_name.substr(p+1);return true;}

	std::shared_ptr<OpenGLTexture> Initialize_Texture(const std::string& name,OpenGLTextureInstance<ushort>* p=nullptr)
	{
		std::string prefix,ext;bool has_ext=Parse_Name(name,prefix,ext);if(!has_ext)return nullptr;
		OpenGLTextureInstance<ushort>* tex=new OpenGLTextureInstance<ushort>();tex->Initialize(prefix,ext);
		std::shared_ptr<OpenGLTexture> ptr=std::make_shared<OpenGLTexture>();ptr.reset(tex);
		texture_hashtable.insert(std::make_pair(tex->name,ptr));return ptr;
	}

	std::shared_ptr<OpenGLTexture> Initialize_Texture(const std::string& name,OpenGLCubeMap<ushort>* p=nullptr)
	{
		std::string prefix,ext;bool has_ext=Parse_Name(name,prefix,ext);if(!has_ext)return nullptr;
		OpenGLCubeMap<ushort>* tex=new OpenGLCubeMap<ushort>();tex->Initialize(prefix,ext);
		std::shared_ptr<OpenGLTexture> ptr=std::make_shared<OpenGLTexture>();ptr.reset(tex);
		texture_hashtable.insert(std::make_pair(tex->name,ptr));return ptr;	
	}

	std::shared_ptr<OpenGLTexture> Lazy_Initialize_Texture(const std::string& name,const TextureType type=TextureType::Tx2d)
	{
		switch(type){
		case TextureType::Tx1d:{OpenGLTexture1D<ushort>*p =nullptr;return Initialize_Texture(name,p);}
		case TextureType::Tx2d:{OpenGLTextureInstance<ushort>* p=nullptr;return Initialize_Texture(name,p);}
		case TextureType::Tx3d:{OpenGLTextureVolume<>* p=nullptr;return Initialize_Texture< >(name,p);}
		case TextureType::Tx3d2:{OpenGLTextureVectorVolume<ushort,2>* p=nullptr;return Initialize_Texture(name,p);}
		case TextureType::Tx3d3:{OpenGLTextureVectorVolume<ushort,3>* p=nullptr;return Initialize_Texture(name,p);}
		case TextureType::Tx3d4:{OpenGLTextureVectorVolume<ushort,4>* p=nullptr;return Initialize_Texture(name,p);}
		case TextureType::TxCube:{OpenGLCubeMap<ushort>* p=nullptr;return Initialize_Texture(name,p);}}
		return nullptr;

		//////image_names, image_types, tex_types
		//////tex_types: 0-tex 2d, 1-cubemap, 2-tex 3d, 3-tex 1d, 4-vec2 tex 3d, 5-vec3 tex 3d, 6-vec4 tex 3d
		//Array<std::tuple<std::string,std::string,int> > tex_info=
		//	{{"tiles","bmp",0},{"skybox","tga",1},{"moon","bmp",1},{"volume","",2},{"1d","",3},{"1d2","",3},{"1d3","",3},{"1d4","",3},
		//	{"vec2_vol","",4},{"vec3_vol","",5},{"vec4_vol","",6}};
		//std::shared_ptr<OpenGLTexture> ptr=nullptr;
		//for(auto i=0;i<tex_info.size();i++){
		//	std::string image_name=std::get<0>(tex_info[i]);if(image_name!=name)continue;
		//	std::string image_type=std::get<1>(tex_info[i]);int tex_type=std::get<2>(tex_info[i]);
		//	switch(tex_type){
		//	case 0:{OpenGLTextureInstance<ushort>* tex=new OpenGLTextureInstance<ushort>();
		//		tex->Initialize(image_name,image_type);ptr.reset(tex);
		//		texture_hashtable.insert(std::make_pair(tex->name,ptr));}return ptr;
		//	case 1:{OpenGLCubeMap<ushort>* tex=new OpenGLCubeMap<ushort>();
		//		tex->Initialize(image_name,image_type);ptr.reset(tex);
		//		texture_hashtable.insert(std::make_pair(tex->name,ptr));}return ptr;
		//	case 2:{OpenGLTextureVolume<>* tex=new OpenGLTextureVolume<>();
		//		tex->Initialize(image_name);ptr.reset(tex);
		//		texture_hashtable.insert(std::make_pair(tex->name,ptr));}return ptr;
		//	case 3:{OpenGLTexture1D<ushort>* tex=new OpenGLTexture1D<ushort>();
		//		tex->Initialize(image_name);ptr.reset(tex);
		//		texture_hashtable.insert(std::make_pair(tex->name,ptr));}return ptr;
		//	case 4:{OpenGLTextureVectorVolume<ushort,2>* tex=new OpenGLTextureVectorVolume<ushort,2>();
		//		tex->Initialize(image_name);ptr.reset(tex);
		//		texture_hashtable.insert(std::make_pair(tex->name,ptr));}return ptr;
		//	case 5:{OpenGLTextureVectorVolume<ushort,3>* tex=new OpenGLTextureVectorVolume<ushort,3>();
		//		tex->Initialize(image_name);ptr.reset(tex);
		//		texture_hashtable.insert(std::make_pair(tex->name,ptr));}return ptr;
		//	case 6:{OpenGLTextureVectorVolume<ushort,4>* tex=new OpenGLTextureVectorVolume<ushort,4>();
		//		tex->Initialize(image_name);ptr.reset(tex);
		//		texture_hashtable.insert(std::make_pair(tex->name,ptr));}return ptr;}}

		//return ptr;
	}
};

inline std::shared_ptr<OpenGLTexture> Get_Texture(const std::string& name,const TextureType type)
{return Texture_Library::Instance()->Get(name,type);}

inline bool Bind_Texture(const std::string& name,const TextureType type,const GLuint idx=0)
{std::shared_ptr<OpenGLTexture> tex=Get_Texture(name,type);if(tex==nullptr)return false;tex->Bind(idx);return true;}

inline void Unbind_Texture(){glBindTexture(GL_TEXTURE_2D,0);}
};
#endif