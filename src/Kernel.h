#pragma once
#include "Global.h"

#define PI 3.141592f

// 1. Replace length
// 2. what if r_diff is zero vector -> should deal in particle.cpp ?

float W_poly6(const glm::vec3& r_diff, const float& support) {
	float dist = glm::sqrt(glm::pow(r_diff.x, 2.0f) + glm::pow(r_diff.y, 2.0f) + glm::pow(r_diff.z, 2.0f));
	if (dist <= support) {
		return 315.0f / (float)(64.0f * PI * pow(support, 9)) * pow((pow(support, 2) - pow(dist,2)), 3);
	}
	else 
		return 0.0f;
}

glm::vec3 Gradient_W_poly6(const glm::vec3& r_diff, const float& support) {
	glm::vec3 normalized = glm::normalize(r_diff);
	if (length(r_diff) <= support) {
		return (float)(-945.0f / 32.0f * PI * pow(support, 9)) * (float)(pow(support, 2) - pow(length(r_diff), 2), 2) * normalized;
	}
	else
		return glm::vec3(0.0f, 0.0f, 0.f);
}

float Laplacian_W_poly6(const glm::vec3& r_diff, const float& support) {
	float dist = glm::sqrt(glm::pow(r_diff.x, 2.0f) + glm::pow(r_diff.y, 2.0f) + glm::pow(r_diff.z, 2.0f));
	if (dist <= support) {
		return (float)(945.0f / 8.0f * PI * pow(support, 9) * (pow(support, 2) - dist) * (float)(dist - 3.0f * (pow(support, 2) - dist) / 4.0f));
	}
	else
		return 0.0f;
}

glm::vec3 Gradient_W_spiky(const glm::vec3& r_diff, const float& support) {
	glm::vec3 normalized = glm::normalize(r_diff);
	if (length(normalized) != 1)
		return glm::vec3(0.0f, 0.0f, 0.0f);
	if (length(r_diff) <= support) {
		return (float)(-45.0f * pow((support - length(r_diff)), 2) / (PI * pow(support, 6))) *  normalized;
	}
	else
		return glm::vec3(0.0f, 0.0f, 0.0f);
}

float Laplacian_W_viscosity(const glm::vec3& r_diff, const float& support) {
	/* W_viscosity 
	15.0f / (2 * PI * pow(support, 3)) * (-pow(dist, 3) / (2.0 * pow(support, 3)) + (pow(dist, 2) / pow(support, 2)) + (support / 2.0 * dist) - 1.0);
	*/
	
	if (length(r_diff) <= support) {
		return 45.0f / (PI * pow(support, 6)) * (support - length(r_diff));
	}
	else
		return 0.0f;
}
