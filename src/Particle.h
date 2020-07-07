#pragma once
#include "Global.h"
#include "Parameter.h"


class Particle {
private:
	const unsigned int m_particle_id;
	float m_mass;
	float m_density;
	float m_pressure;
	bool isStatic;
	
	glm::vec3 m_position;
	glm::vec3 m_velocity;

	glm::vec3 m_force;
	glm::vec3 m_pressureforce;
	glm::vec3 m_accerlation;
public:
	std::list<const Particle*> neighbors;
	glm::vec3 predicted_velocity;
	glm::vec3 predicted_position;

public:
	Particle(const unsigned int& id, const glm::vec3& position,const glm::vec3& velocity, bool isStatic);
	~Particle();

	void Init(const unsigned int& m_particle_id, float& mass, float& density);

	void ComputeDensity_SPH(const float& support);
	void ComputePressure_SPH();
	glm::vec3 ComputePressureForce_SPH(const float& support);
	glm::vec3 ComputeViscosity_SPH(const float& support);
	glm::vec3 ComputeSurfaceTension_SPH(const float& support);

	//void Move(glm::vec3& dir);

	
	float GetDensity() const;
	float GetMass() const;
	float GetPressure() const;

	int GetID() const;
	glm::vec3 GetPos() const;
	glm::vec3 GetVelocity() const;
	glm::vec3 GetForce() const;

	void SetPos(glm::vec3& pos);
	void SetPressure(float& pressure);
	void SetVelocity(glm::vec3& vel);

	void SetForce(glm::vec3& force);
	void SetPressureForce(glm::vec3 & force);
};