#pragma once

#define GLEW_STATIC true

#include <GL/glew.h>
#include <GLFW/glfw3.h>
#include <IL/il.h>
#include <IL/ilu.h>
#include <IL/ilut.h>
#include <stdlib.h>
#include <time.h>

#include "glm/glm.hpp"
#include "glm/gtc/matrix_transform.hpp"

#include <fstream>
#include <sstream>
#include <iostream>
#include <string>
#include <vector>
#include <list>
#include <unordered_map>
#include <algorithm>

#define ASSERT(x) if (!(x)) __debugbreak();
#define GLCall(x) GLClearError();\
    x;\
    ASSERT(GLLogCall(#x, __FILE__, __LINE__)) // # turns func into string 



inline void GLClearError() {                // inline 안해주면 LNK1169
    while (glGetError() != GL_NO_ERROR);
}

inline bool GLLogCall(const char* function, const char* file, int line) {  // inline 안해주면 LNK1169
    while (GLenum error = glGetError()) {
        std::cout << "[OpenGL Error (" << error << ")" << function << file << line << std::endl;
        return false;
    }
    return true;
}

inline float sqr(const glm::vec3& vec) {
    return vec.x * vec.x + vec.y * vec.y + vec.z * vec.z;
}

inline float length(const glm::vec3& vec) {
    return sqrt(sqr(vec));
}

inline void printmat3(const glm::mat3& mat) {
    std::cout << "[ " << mat[0][0] << " " << mat[0][1] << " " << mat[0][2] << " ]" << std::endl;
    std::cout << "[ " << mat[1][0] << " " << mat[1][1] << " " << mat[1][2] << " ]" << std::endl;
    std::cout << "[ " << mat[2][0] << " " << mat[2][1] << " " << mat[2][2] << " ]" << std::endl;
}