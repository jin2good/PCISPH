#pragma once
#include "Global.h"
#include "Parameter.h"


class Particle {
public:
	const unsigned int m_particle_id;
	float m_mass;
	float m_density;
	float m_pressure;
	glm::mat3 stress;
	bool isStatic;
	
	glm::vec3 m_position;
	glm::vec3 m_velocity;
	glm::vec3 m_accerlation;

	glm::vec3 m_force;
	glm::vec3 m_pressureforce;
	glm::vec3 m_viscosityforce;
	glm::vec3 m_extforce;


	std::list<const Particle*> neighbors;

	glm::vec3 last_velocity;
	glm::vec3 last_position;

public:
	Particle(const unsigned int& id, const glm::vec3& position,const glm::vec3& velocity, bool isStatic);
	~Particle();

	void ComputeDensity_SPH(const float& support);
	void ComputePressure_SPH();
	glm::vec3 ComputePressureForce_SPH(const float& support);
	glm::vec3 ComputeViscosity_SPH(const float& support);
	glm::vec3 ComputeSurfaceTension_SPH(const float& support);

	
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