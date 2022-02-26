#ifndef __OpenGLMaterial_h__
#define __OpenGLMaterial_h__
#include "glm.hpp"

class Material
{public:
	glm::vec4 mat_amb=glm::vec4(1.f);
	glm::vec4 mat_dif=glm::vec4(1.f,1.f,1.f,1.f);
	glm::vec4 mat_spec=glm::vec4(1.f);
	glm::vec4 mat_shinness=glm::vec4(32.f,0.f,0.f,0.f);
};

#endif
