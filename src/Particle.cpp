#include "Particle.h"
#include "Kernel.h"

inline bool isvecnan(glm::vec3& vec) {
	return (std::isnan(vec.x) || std::isnan(vec.y) || std::isnan(vec.z));
}



Particle::Particle(const unsigned int& id,const glm::vec3& position,const glm::vec3& velocity, bool isStatic)
	:m_particle_id(id),
	 m_mass(MASS),
	 m_density(REST_DENSITY),
	 m_position(position),
	 m_velocity(velocity),
	 isStatic(isStatic)
{
	neighbors.clear();
}

Particle::~Particle()
{
}



void Particle::ComputeDensity_SPH(const float& support)
{
	float new_density = 0.0f;
	for (const auto& particle : neighbors) {
		float neighbor_mass = particle->GetMass();
		glm::vec3 neighbor_pos = particle->GetPos();

		glm::vec3 r_diff = this->GetPos() - neighbor_pos;
		new_density += neighbor_mass * W_poly6(r_diff, support);

	}
	
	m_density = new_density;
	if (ENABLE_DEBUG_MODE && SHOW_DENSITY && (this->GetID() == WATCH_PARTICLE)) {
		std::cout << "Computed Density for Particle[" << this->GetID() << "] is " << new_density << std::endl;
	}
}


void Particle::ComputePressure_SPH()
{
	float pressure = GAS_CONSTANT * (this->m_density - REST_DENSITY);
	if (ENABLE_DEBUG_MODE && SHOW_PUESSURE && (this->GetID() == WATCH_PARTICLE)) {
		std::cout << "Pressure for Particle[" << this->GetID() << "] is " << pressure << std::endl;
	}
	m_pressure = pressure;
}

glm::vec3 Particle::ComputePressureForce_SPH(const float& support)
{
	glm::vec3 pf = glm::vec3(0 ,0 ,0);

	for (const auto& particle : neighbors){
		if (particle->GetID() == this->GetID())
			continue;
		float neighbor_mass = particle->GetMass();
		glm::vec3 neighbor_pos = particle->GetPos();
		float neighbor_pressure = particle->GetPressure();
		float neighbor_density = particle->GetDensity();

		glm::vec3 r_diff = this->GetPos() - neighbor_pos;

		pf -= neighbor_mass * (this->GetPressure() + neighbor_pressure) / (2.0f * neighbor_density) * Gradient_W_spiky(r_diff, support);
	}

	if (ENABLE_DEBUG_MODE && SHOW_PRESSUREFORCE && (this->GetID() == WATCH_PARTICLE)) {
		std::cout << "PressureForce for Particle[" << this->GetID() << "] is " << pf[0] << "," << pf[1] << "," << pf[2] << std::endl;
	}
	return pf;
}

glm::vec3 Particle::ComputeViscosity_SPH(const float& support)
{
	glm::vec3 vf = glm::vec3(0, 0, 0);
	for (const auto& particle : neighbors) {
		if (particle->GetID() == this->GetID()) // 자기자신은 neighbor 아님
			continue;

		glm::vec3 neighbor_pos = particle->GetPos();
		float neighbor_mass = particle->GetMass();
		float neighbor_density = particle->GetDensity();
		glm::vec3 neighbor_vel = particle->GetVelocity();

		glm::vec3 r_diff = this->GetPos() - neighbor_pos;
		glm::vec3 per_neighbor = (float)VISCOSITY_CONSTANT * neighbor_mass * (neighbor_vel - this->GetVelocity()) / neighbor_density * Laplacian_W_viscosity(r_diff, support);
		
		if(!isvecnan(per_neighbor))
			vf += per_neighbor;
	}
	if (ENABLE_DEBUG_MODE && SHOW_VISCOSITY && (this->GetID() == WATCH_PARTICLE)) {
		std::cout << "ViscosityForce for Particle[" << this->GetID() << "] is " << vf[0] << "," << vf[1] << "," << vf[2] << std::endl;
	}
	return vf;
}

glm::vec3 Particle::ComputeSurfaceTension_SPH(const float& support)
{
	glm::vec3 sf = glm::vec3(0, 0, 0);

	float colorField = 0.0f;
	glm::vec3 surfaceNormal = glm::vec3(0, 0, 0);
	float curvature = 0.0f;
	float Laplacian = 0.0f;
	glm::vec3 surfaceTraction = glm::vec3(0, 0, 0);

	for (const auto& particle : neighbors) {
		if (particle->GetID() == this->GetID())
			continue;

		glm::vec3 neighbor_pos = particle->GetPos();
		float neighbor_mass = particle->GetMass();
		float neighbor_density = particle->GetDensity();

		glm::vec3 r_diff = this->GetPos() - neighbor_pos;

		colorField += (neighbor_mass / neighbor_density) * W_poly6(r_diff, support);
		surfaceNormal += (neighbor_mass / neighbor_density) * Gradient_W_poly6(r_diff, support);
		Laplacian += (neighbor_mass / neighbor_density) * Laplacian_W_poly6(r_diff, support);
	}
	if (length(surfaceNormal) > ST_THRESHOLD) {
		curvature = -Laplacian / length(surfaceNormal);
		//surfaceTraction = (VISCOSITY_CONSTANT * curvature / length(surfaceNormal)) * surfaceNormal;
		sf = (float)TENSION_COEFF * curvature * surfaceNormal;
		if (isvecnan(sf))
			sf = glm::vec3(0, 0, 0);
	}
	if (ENABLE_DEBUG_MODE && SHOW_SURFACETENSION && (this->GetID() == WATCH_PARTICLE)) {
		std::cout << "Surface Tension for Particle[" << this->GetID() << "] is " << sf[0] << "," << sf[1] << "," << sf[2] << std::endl;
	}
	return sf;
}


float Particle::GetDensity() const
{
	return m_density;
}

float Particle::GetMass() const
{
	return m_mass;
}

float Particle::GetPressure() const
{
	return m_pressure;
}

int Particle::GetID() const
{
	return m_particle_id;
}

glm::vec3 Particle::GetPos() const
{
	return m_position;
}

glm::vec3 Particle::GetVelocity() const
{
	return m_velocity;
}

glm::vec3 Particle::GetForce() const
{
	return this->m_force;
}

void Particle::SetPos(glm::vec3& pos)
{
	if (isStatic)
	{
		std::cout << "Tried to change Position of Static Particle[" << m_particle_id << "]" << std::endl;
	}
	else {
		this->last_position = this->m_position;
		this->m_position = pos;
	}
}

void Particle::SetPressure(float& pressure)
{
	this->m_pressure = pressure;
}

void Particle::SetVelocity(glm::vec3& vel)
{
	if (isStatic)
	{
		std::cout << "Tried to change Veloctiy of Static Particle" << m_particle_id << "]" << std::endl;
	}
	else {
		this->last_velocity = this->m_velocity;
		this->m_velocity = vel;
	}
		
}


void Particle::SetForce(glm::vec3& force)
{
	this->m_force = force;
}

void Particle::SetPressureForce(glm::vec3& force)
{
	this->m_pressureforce = force;
}

