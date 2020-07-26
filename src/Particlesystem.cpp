#include "Particlesystem.h"

#include "VertexBufferLayout.h"

ParticleSystem::ParticleSystem(int mode)
	:SPHmode(mode)
{
	this->Cell = std::vector<std::list<const Particle*>>(NUM_CELL);
	Init();
}

ParticleSystem::~ParticleSystem()
{
}

void ParticleSystem::Init()
{
	this->sphStart = clock();

	SetInitialParticlePosition();

	if (ENABLE_DEBUG_MODE) {
		std::cout << "Particle System" << std::endl;
		std::cout << "Number of Fluid Particles Created:  " << PARTICLE_COUNT << std::endl;
		std::cout << "Number of Static Particles Created: " << particle_list.size() - PARTICLE_COUNT << std::endl;
		std::cout << "Number of Total Particles Created:  " << particle_list.size() << std::endl;
	}
	/* Load particle position into particle_locaiton */
	LoadParticleVectorPosition();

	#ifdef USE_PRECOMPUTED_VALUE
	// Calculate Precompued Values
	#endif

	this->sphEnd = clock();
	sphtime += (double)(sphEnd - sphStart);
}

void ParticleSystem::Update(double deltaTime)
{
	if( ENABLE_DEBUG_MODE && SHOW_TIME)
		this->sphStart = clock();
	
	CellInsert(); // Cell
	if (SPHmode == SPH)
		Update_SPH(deltaTime);
	if (SPHmode == PCISPH)
		Update_PCISPH(deltaTime);
	if (SPHmode == GRANULAR)
		Update_GRANULAR(deltaTime);

	if (ENABLE_DEBUG_MODE && SHOW_TIME) {
		this->sphEnd = clock();
		sphtime += (double)(sphEnd - sphStart);
	}
}

void ParticleSystem::Update_SPH(double deltaTime)	
{
	timestep = dt;
	
	/* Find Neighbors */
	FindNeighbors();

	if (GET_NEIGHBOR_INFO) {
		int max = 0;
		for (int i = 0; i < particle_list.size(); i++) {
			auto particle = particle_list[i];
			int neighbor_count = particle.neighbors.size();
			max = std::max(neighbor_count, max);
		}
		std::cout << max << std::endl;
	}

	/* Compute Density */
	ComputeDensity();
	/* Compute Pressure */
	ComputePressure();

	/* Compute Forces */
	ComputeForce();

	/* Update Velocity and Position */
	ComputeVelocityandPosition(timestep);

	LoadParticleVectorPosition();
	if (ENABLE_DEBUG_MODE && SHOW_TIME) {
		std::cout << "Elapsed Time: " << sphtime << std::endl;
	}
}

void ParticleSystem::Update_PCISPH(double deltaTime)
{
	timestep = dt;

	FindNeighbors();

	for (unsigned int particle_id = 0; particle_id < particle_list.size(); particle_id++)
	{
		auto& particle = particle_list[particle_id];
		
		/* Compute Force v,g,ext */
		particle.m_viscosityforce = particle.ComputeViscosity_SPH(KERNEL) / particle.GetDensity();
		particle.m_extforce = glm::vec3(0, -GRAVITY, 0);
		
		/* Initialied Pressure to 0.0 */
		float initial_pressure = 0.0f;
		particle.SetPressure(initial_pressure);

		/* Initialized PressureForce to zero */
		glm::vec3 zerovector = glm::vec3(0.0f, 0.0f, 0.f);
		particle.SetPressureForce(zerovector);
		particle.m_force = particle.m_extforce + particle.m_viscosityforce + particle.m_pressureforce;
		
		if (ENABLE_DEBUG_MODE && SHOW_NONPUESSUREFORCE && particle_id == WATCH_PARTICLE) {
			std::cout << "Computed nonPressureForce for Particle[" << particle_id << "] is " << particle.m_force[0] << "," << particle.m_force[1] << "," << particle.m_force[2] << std::endl;
		}
	}

	int iter = 0;
	float average_density_error = 0.0f;
	bool error_tolerable = false;

	if (ENABLE_DEBUG_MODE && SHOW_TIME) {
		EstimatePressureStart = clock();
	}
	while ((!error_tolerable || (iter <= MINITERATIONS)) && (iter <= MAXITERATIONS)) // density_error
	{
		if (ENABLE_DEBUG_MODE && SHOW_TIME) {
			this->PredictNewPosVelStart = clock();
		}
		//Predict New Velocity and Position
		ComputeVelocityandPosition(timestep);

		if (ENABLE_DEBUG_MODE && SHOW_TIME) {
			this->PredictNewPosVelEnd = clock();
		}
	#pragma omp parallel for
		for (int particle_id = 0; particle_id < particle_list.size(); particle_id++)
		{
			auto& particle = particle_list[particle_id];
			/* predict density */
			particle.ComputeDensity_SPH(KERNEL);
			/* predict density variation */
			float density_err = (std::max(particle.m_density, REST_DENSITY) - REST_DENSITY); //when density is estimated to be big enough(aka. compressed)

			average_density_error += density_err;

			//float delta = PRECOMPUTED_VALUE * (REST_DENSITY * REST_DENSITY)  / (2.0f * timestep * timestep * MASS * MASS);
			if (ENABLE_DEBUG_MODE && SHOW_PUESSURE && particle_id == WATCH_PARTICLE) {
				std::cout << "Density Error : " << density_err << std::endl;
			}

			/* Update Pressure */
			particle.m_pressure += density_err * PRECOMPUTED_VALUE;
			if (ENABLE_DEBUG_MODE && SHOW_PUESSURE && particle_id == WATCH_PARTICLE) {
				std::cout << "Computed Pressure for Particle[" << particle.GetID() << "] is " << particle.m_pressure << std::endl;
			}
		}
		average_density_error /= PARTICLE_COUNT;
		
		if (average_density_error < DENSITY_FLUCTUATION_THRESHOLD)
			error_tolerable = true;
		else false;

	#pragma omp parallel for
		for (int particle_id = 0; particle_id < PARTICLE_COUNT; particle_id++)
		{
			/* compute pressure force */
			auto& particle = particle_list[particle_id];
			particle.m_pressureforce = particle.ComputePressureForce_SPH(KERNEL) / particle.GetDensity();
			
			particle.m_force = particle.m_extforce + particle.m_viscosityforce + particle.m_pressureforce;
			particle.m_velocity = particle.last_velocity;
			particle.m_position = particle.last_position;
		}
		if (ENABLE_DEBUG_MODE)
			std::cout << std::endl;
		
		iter++;
	}

	if (ENABLE_DEBUG_MODE && SHOW_TIME) {
		this->EstimatePressureEnd = clock();
		this->estimatepressuretime += (double)(this->EstimatePressureEnd - this->EstimatePressureStart);
	}

//#pragma omp parallel for
	for (int particle_id = 0; particle_id < PARTICLE_COUNT; particle_id++)
	{
		auto& particle = particle_list[particle_id];
		particle.m_force = particle.m_extforce + particle.m_viscosityforce + particle.m_pressureforce;
		glm::vec3 force = particle.m_force;
		if (ENABLE_DEBUG_MODE && SHOW_TOTAL_FORCE && (particle_id == WATCH_PARTICLE)) {
			std::cout << "Total Force for Particle[" << particle_id << "] is " << force[0] << "," << force[1] << "," << force[2] << std::endl;
		}
		particle.m_position = particle.last_position;
		particle.m_velocity = particle.last_velocity;
		ComputeVelocityandPosition_SPH(particle_id, timestep);
	}

	LoadParticleVectorPosition();

	if (ENABLE_DEBUG_MODE && SHOW_TIME) {
		std::cout << "Elapsed Time:      " << this->sphtime << std::endl;
		std::cout << "NeighborSearch:    " << this->neighborsearchtime << std::endl;
		std::cout << "Estimate Pressure: " << this->estimatepressuretime << std::endl;
	}
}

void ParticleSystem::Update_GRANULAR(double deltaTime) {
	timestep = dt;

	FindNeighbors();

	for (unsigned int particle_id = 0; particle_id < particle_list.size(); particle_id++)
	{
		auto& particle = particle_list[particle_id];

		/* Compute Force v,g,ext */
		particle.m_viscosityforce = particle.ComputeViscosity_SPH(KERNEL) / particle.GetDensity();
		particle.m_extforce = glm::vec3(0, -GRAVITY, 0);

		/* Initialied Pressure to 0.0 */
		float initial_pressure = 0.0f;
		particle.SetPressure(initial_pressure);

		/* Initialized PressureForce to zero */
		glm::vec3 zerovector = glm::vec3(0.0f, 0.0f, 0.f);
		particle.SetPressureForce(zerovector);

		/* Initialize stress*/
		particle.stress = glm::mat3(zerovector, zerovector, zerovector);

		if (ENABLE_DEBUG_MODE && SHOW_NONPUESSUREFORCE && particle_id == WATCH_PARTICLE) {
			std::cout << "Computed nonPressureForce for Particle[" << particle_id << "] is " << particle.m_force[0] << "," << particle.m_force[1] << "," << particle.m_force[2] << std::endl;
		}
	}

	int iter = 0;
	float average_density_error = 0.0f;
	float average_friction = 0.0f;
	float average_cohesion = 0.0f;
	bool error_tolerable = false;
	bool stress_tolerable = false;


	if (ENABLE_DEBUG_MODE && SHOW_TIME) {
		EstimatePressureStart = clock();
	}
	while ( (!error_tolerable || !stress_tolerable || (iter <= MINITERATIONS) ) && (iter <= MAXITERATIONS)) // density_error
	{

#pragma omp parallel for /* Update Force*/
		for (int particle_id = 0; particle_id < particle_list.size(); particle_id++) { 
			auto& particle = particle_list[particle_id];

			particle.m_pressureforce = particle.ComputePressureForce_SPH(KERNEL) / particle.GetDensity();
			particle.m_friction_cohesion = particle.ComputeFrictionCohesion(KERNEL);
			particle.m_force = particle.m_extforce + particle.m_viscosityforce + particle.m_pressureforce + particle.m_friction_cohesion;
		}

		if (ENABLE_DEBUG_MODE && SHOW_TIME) {
			this->PredictNewPosVelStart = clock();
		}

		//Predict New Velocity and Position
		ComputeVelocityandPosition(timestep);

		if (ENABLE_DEBUG_MODE && SHOW_TIME) {
			this->PredictNewPosVelEnd = clock();
		}

#pragma omp parallel for
		for (int particle_id = 0; particle_id < particle_list.size(); particle_id++)
		{
			auto& particle = particle_list[particle_id];
			/* predict density */
			particle.ComputeDensity_SPH(KERNEL);
			
			/* predict density variation */
			float density_err = (std::max(particle.m_density, REST_DENSITY) - REST_DENSITY); //when density is estimated to be big enough(aka. compressed)

			average_density_error += density_err;

			if (ENABLE_DEBUG_MODE && SHOW_PUESSURE && particle_id == WATCH_PARTICLE) {
				std::cout << "Density Error : " << density_err << std::endl;
			}

			/* update pressure */
			particle.m_pressure += density_err * PRECOMPUTED_VALUE;
			if (ENABLE_DEBUG_MODE && SHOW_PUESSURE && particle_id == WATCH_PARTICLE) {
				std::cout << "Computed Pressure for Particle[" << particle.GetID() << "] is " << particle.m_pressure << std::endl;
			}

			if (ENABLE_DEBUG_MODE && SHOW_TIME) {
				this->EstimatePressureEnd = clock();
				this->estimatepressuretime += (double)(this->EstimatePressureEnd - this->EstimatePressureStart);
			}
			/*-------------------------------------------------------------------------------------------------------------------------------------------------*/
			/* Predict Strain Rate */
			glm::vec3 zerovector = glm::vec3(0.0f, 0.0f, 0.f);
			glm::mat3 velocity_gradient = particle.VelocityGradient(KERNEL);
			glm::mat3 strain_rate = glm::mat3(0.5, 0, 0, 0, 0.5, 0, 0, 0, 0.5) * (velocity_gradient + glm::transpose(velocity_gradient));
			
			std::cout << "VG: " << velocity_gradient[0][0] << velocity_gradient[0][1] << velocity_gradient[0][2] << std::endl << velocity_gradient[1][0] << velocity_gradient[1][1] << velocity_gradient[1][2] << std::endl << velocity_gradient[2][0] << velocity_gradient[2][1] << velocity_gradient[2][2] << std::endl;
			/* Compute Dissipative Stress */
			#ifdef USE_PRECOMPUTED_VALUE
			particle.stress += inv_precomputed_stress * strain_rate;
			#endif

			/* Test yield and cohesion */
			glm::mat3 deviatoric = particle.stress - (0.3f * trace(particle.stress)) * glm::mat3(1.0f);

			//average_friction += Fnorm(deviatoric) < (2.0 / 3.0)* std::sin(30 * 3.14159 / 180.0) * std::sin(30 * 3.14159 / 180.0) * particle.GetPressure() * particle.GetPressure();
			//average_cohesion += Fnorm(particle.stress) < COHESION;
		}

		average_density_error /= PARTICLE_COUNT;

		if (average_density_error < DENSITY_FLUCTUATION_THRESHOLD)
			error_tolerable = true;
		else false;


#pragma omp parallel for /* Rollback Velocity and Position */
		for (int particle_id = 0; particle_id < PARTICLE_COUNT; particle_id++)
		{
			
			auto& particle = particle_list[particle_id];
			
			particle.m_velocity = particle.last_velocity;
			particle.m_position = particle.last_position;
		}
		
		if (ENABLE_DEBUG_MODE)
			std::cout << std::endl;

		iter++;
	}



	//#pragma omp parallel for
	for (int particle_id = 0; particle_id < PARTICLE_COUNT; particle_id++)
	{
		auto& particle = particle_list[particle_id];
		particle.m_force = particle.m_extforce + particle.m_viscosityforce + particle.m_pressureforce;
		glm::vec3 force = particle.m_force;
		if (ENABLE_DEBUG_MODE && SHOW_TOTAL_FORCE && (particle_id == WATCH_PARTICLE)) {
			std::cout << "Total Force for Particle[" << particle_id << "] is " << force[0] << "," << force[1] << "," << force[2] << std::endl;
		}
		ComputeVelocityandPosition_SPH(particle_id, timestep);
	}

	LoadParticleVectorPosition();

	if (ENABLE_DEBUG_MODE && SHOW_TIME) {
		std::cout << "Elapsed Time:      " << this->sphtime << std::endl;
		std::cout << "NeighborSearch:    " << this->neighborsearchtime << std::endl;
		std::cout << "Estimate Pressure: " << this->estimatepressuretime << std::endl;
	}
}

float* ParticleSystem::GetParticlePositionArray()
{
	if (float* particle_loc = (float*)std::malloc((SHOW_PARTICLE_BOX ? this->particle_count : PARTICLE_COUNT) * 3 * sizeof(float))) {
		particle_loc = particle_location.data();
		return particle_loc;
	}
	else {
		std::cout << "malloc for function:GetParticleLocation failed" << std::endl;
		return nullptr;
	}
}


void ParticleSystem::SetInitialParticlePosition()
{
	/* Initialize Position */
	float particle_pos_X =  INITIAL_DIST;
	float particle_pos_Y =  INITIAL_DIST;
	float particle_pos_Z =  INITIAL_DIST;
	glm::vec3 vel = glm::vec3(0, 0, 0);

	for (unsigned int i = 0; i < PARTICLE_COUNT; i++) {

		if (particle_pos_X >= (DAM_BREAKING_MODE ? PARTICLE_INITIAL_BOUNDARY_X/2 - INITIAL_DIST : PARTICLE_INITIAL_BOUNDARY_X - INITIAL_DIST)) {
			particle_pos_X = 2 * INITIAL_DIST;
			particle_pos_Z += INITIAL_DIST;
		}

		if (particle_pos_Z >= PARTICLE_INITIAL_BOUNDARY_Z - 2 * INITIAL_DIST) {
			particle_pos_Z = 2 * INITIAL_DIST;
			particle_pos_Y += INITIAL_DIST;
		}

		if (particle_pos_Y >= PARTICLE_INITIAL_BOUNDARY_Y) {
			std::cout << " Too Many Particles to be initialized" << std::endl;
		}

		glm::vec3 pos = glm::vec3(particle_pos_X, particle_pos_Y, particle_pos_Z);
		particle_list.push_back(Particle(i, pos, vel, false));

		particle_pos_X += INITIAL_DIST;
	}
	
	if (BOUNDARY_MODE == WALL_PARTICLE)
	{
		/* Create Box */
		unsigned int i = PARTICLE_COUNT;

		/* Create Bottom */
		for (particle_pos_Y = - INITIAL_DIST; particle_pos_Y < INITIAL_DIST; particle_pos_Y += INITIAL_DIST / WALL_PARTICLE_DENSITY / 2) {
			for (particle_pos_X = -INITIAL_DIST; particle_pos_X < PARTICLE_INITIAL_BOUNDARY_X + 3.0f * INITIAL_DIST; particle_pos_X += INITIAL_DIST / WALL_PARTICLE_DENSITY) {
				for (particle_pos_Z = -INITIAL_DIST; particle_pos_Z < PARTICLE_INITIAL_BOUNDARY_Z + 3.0f * INITIAL_DIST; particle_pos_Z += INITIAL_DIST / WALL_PARTICLE_DENSITY) {
					glm::vec3 pos = glm::vec3(particle_pos_X, particle_pos_Y, particle_pos_Z);
					particle_list.push_back(Particle(i, pos, vel, true));
					i++;
				}
			}
		}
		/* Create xy Plane */
		for (particle_pos_Z = -INITIAL_DIST; particle_pos_Z < 1.0f * INITIAL_DIST; particle_pos_Z += INITIAL_DIST / WALL_PARTICLE_DENSITY) {
			for (particle_pos_X = 0; particle_pos_X < PARTICLE_INITIAL_BOUNDARY_X + 2.0f * INITIAL_DIST; particle_pos_X += INITIAL_DIST / WALL_PARTICLE_DENSITY) {
				for (particle_pos_Y = 0; particle_pos_Y < PARTICLE_INITIAL_BOUNDARY_Y + 2.0f * INITIAL_DIST; particle_pos_Y += INITIAL_DIST / WALL_PARTICLE_DENSITY) {
					glm::vec3 pos = glm::vec3(particle_pos_X, particle_pos_Y, particle_pos_Z);
					particle_list.push_back(Particle(i, pos, vel, true));
					i++;
				}
			}
		}

		for (particle_pos_Z = PARTICLE_INITIAL_BOUNDARY_Z + INITIAL_DIST; particle_pos_Z < PARTICLE_INITIAL_BOUNDARY_Z + 3.0f * INITIAL_DIST; particle_pos_Z += INITIAL_DIST / WALL_PARTICLE_DENSITY) {
			for (particle_pos_X = 0; particle_pos_X < PARTICLE_INITIAL_BOUNDARY_X + 2.0f * INITIAL_DIST; particle_pos_X += INITIAL_DIST / WALL_PARTICLE_DENSITY) {
				for (particle_pos_Y = 0; particle_pos_Y < PARTICLE_INITIAL_BOUNDARY_Y + 2.0f * INITIAL_DIST; particle_pos_Y += INITIAL_DIST / WALL_PARTICLE_DENSITY) {
					glm::vec3 pos = glm::vec3(particle_pos_X, particle_pos_Y, particle_pos_Z);
					particle_list.push_back(Particle(i, pos, vel, true));
					i++;
				}
			}
		}
		/* Create zy Plane */
		for (particle_pos_X = -INITIAL_DIST; particle_pos_X < 1.0f * INITIAL_DIST; particle_pos_X += INITIAL_DIST / WALL_PARTICLE_DENSITY) {
			for (particle_pos_Z = 0; particle_pos_Z < PARTICLE_INITIAL_BOUNDARY_Z + 2.0f * INITIAL_DIST; particle_pos_Z += INITIAL_DIST / WALL_PARTICLE_DENSITY) {
				for (particle_pos_Y = 0; particle_pos_Y < PARTICLE_INITIAL_BOUNDARY_Y + 2.0f * INITIAL_DIST; particle_pos_Y += INITIAL_DIST / WALL_PARTICLE_DENSITY) {
					glm::vec3 pos = glm::vec3(particle_pos_X, particle_pos_Y, particle_pos_Z);
					particle_list.push_back(Particle(i, pos, vel, true));
					i++;
				}
			}
		}

		for (particle_pos_X = PARTICLE_INITIAL_BOUNDARY_X + INITIAL_DIST; particle_pos_X < PARTICLE_INITIAL_BOUNDARY_X + 3.0f * INITIAL_DIST; particle_pos_X += INITIAL_DIST / WALL_PARTICLE_DENSITY) {
			for (particle_pos_Z = 0; particle_pos_Z < PARTICLE_INITIAL_BOUNDARY_Z + 2.0f * INITIAL_DIST; particle_pos_Z += INITIAL_DIST / WALL_PARTICLE_DENSITY) {
				for (particle_pos_Y = 0; particle_pos_Y < PARTICLE_INITIAL_BOUNDARY_Y + 2.0f * INITIAL_DIST; particle_pos_Y += INITIAL_DIST / WALL_PARTICLE_DENSITY) {
					glm::vec3 pos = glm::vec3(particle_pos_X, particle_pos_Y, particle_pos_Z);
					particle_list.push_back(Particle(i, pos, vel, true));
					i++;
				}
			}
		}
	}
	particle_count = particle_list.size();
	std::cout << "fluid particle count " << PARTICLE_COUNT << std::endl;
	std::cout << "wall particle count " << particle_count - PARTICLE_COUNT << std::endl;
	std::cout << "total particle count " << particle_count << std::endl;
}

void ParticleSystem::LoadParticleVectorPosition()
{
	particle_location.clear();
	/* Create VertexBuffer with vector<float>(particle_location) for position */
	for (unsigned int particle = 0; particle < (SHOW_PARTICLE_BOX ? this->particle_count : PARTICLE_COUNT); particle++)
	{
		if (ENABLE_DEBUG_MODE && SHOW_LOAD_PARTICLE_VECTOR && particle == WATCH_PARTICLE) {
			std::cout << particle_list[particle].GetPos().x << ",";
			std::cout << particle_list[particle].GetPos().y << ",";
			std::cout << particle_list[particle].GetPos().z << std::endl;
		}
		
		particle_location.push_back(particle_list[particle].GetPos().x); //push in x
		particle_location.push_back(particle_list[particle].GetPos().y); //push in y
		particle_location.push_back(particle_list[particle].GetPos().z); //push in z
	}
}

void ParticleSystem::FindNeighbors()
{
	if (ENABLE_DEBUG_MODE && SHOW_TIME) {
		this->NeighborsearchStart = clock();
	}
	
	if (GRID_BASED_NEIGHBORHOOD_SEARCH) {
#pragma omp parallel for
		for (int i = 0; i < PARTICLE_COUNT; i++) {
			auto& particle = particle_list[i];

			glm::vec3 ind = GetCellIndex(particle.GetPos());
			particle.neighbors.clear();
			for (int ix = -1; ix < 2; ix++) {
				for (int iy = -1; iy < 2; iy++) {
					for (int iz = -1; iz < 2; iz++) {

						int nx = (int)ind.x + ix;
						int ny = (int)ind.y + iy;
						int nz = (int)ind.z + iz;

						if (nx > GRID_X || nx < 0)
							continue;
						if (ny > GRID_Y || ny < 0)
							continue;
						if (nz > GRID_Z || nz < 0)
							continue;

						for (auto& elem : Cell[CellCoord(nx, ny, nz)]) {
							glm::vec3 rdiff = particle.GetPos() - elem->GetPos();
							if (sqr(rdiff) < KERNEL * KERNEL)
							{
								particle.neighbors.push_back(elem);
							}
						}
					}
				}
			}
		}
	}
	else {
		for (unsigned int i = 0; i < PARTICLE_COUNT; i++)
		{
			auto& particle = particle_list[i];
			NaiveNeighborSearch(i);
		}
	}
	if (ENABLE_DEBUG_MODE && SHOW_TIME) {
		this->NeighborsearchEnd = clock();
		this->neighborsearchtime += (double)(NeighborsearchEnd - NeighborsearchStart);
	}
}

void ParticleSystem::ComputeDensity()
{
#pragma omp parallel for
	for(int i=0; i<particle_list.size(); i++) {
		auto& particle = particle_list[i];
		if (particle.isStatic) {
			particle.m_density = REST_DENSITY;
		}
		else {
			particle.ComputeDensity_SPH(KERNEL);
		}
	}
}

void ParticleSystem::ComputePressure()
{
#pragma omp parallel for
	for (int i = 0; i < particle_list.size(); i++) {
		auto& particle = particle_list[i];
		particle.ComputePressure_SPH();
	}
}

void ParticleSystem::ComputeForce()
{
#pragma omp parallel for
	for (int particle = 0; particle < particle_list.size(); particle++)
	{
		ComputeForce_SPH(particle);
	}
}

void ParticleSystem::ComputeVelocityandPosition(double timestep)
{
#pragma omp parallel for
	for (int particle = 0; particle < PARTICLE_COUNT; particle++)
	{
		ComputeVelocityandPosition_SPH(particle, timestep);
	}
}


void ParticleSystem::NaiveNeighborSearch(const unsigned int& particle_id)
{
	auto& particle = particle_list[particle_id];
	particle.neighbors.clear();
	glm::vec3 pos = particle.GetPos();
	for (const auto& p: particle_list)
	{
		glm::vec3 another_pos = p.GetPos();
		if (glm::distance(pos, another_pos) < CELL_SIZE) {
			particle.neighbors.push_back(&p);
			if (ENABLE_DEBUG_MODE && SHOW_NEIGHBORSEARCH_RESULT) {
				std::cout << "For Particle[" << particle_id << "] Added Particle [" << p.GetID() << "] as Neighbors" << std::endl;
			}
		}
	}
}

void ParticleSystem::CellInsert()
{	
	/* Initialize Grid and GridIndices */
	for (auto& list : Cell) {
		list.clear();
	}

	for (const auto& p : particle_list) {
		auto i = &p - &particle_list[0];  // i is an index
		
		glm::vec3 index = GetCellIndex(p.GetPos());
		Cell[CellCoord((int)index.x, (int)index.y, (int)index.z)].push_back(&p);
	}
}

glm::vec3 ParticleSystem::GetCellIndex(const glm::vec3& pos)
{
	int ind_x = (int)(pos.x / CELL_SIZE);
	int ind_y = (int)(pos.y / CELL_SIZE);
	int ind_z = (int)(pos.z / CELL_SIZE);
	ind_x = std::max(0, std::min(GRID_X, ind_x));
	ind_y = std::max(0, std::min(GRID_Y, ind_y));
	ind_z = std::max(0, std::min(GRID_Z, ind_z));

	return glm::vec3(ind_x, ind_y, ind_z);
}

int ParticleSystem::CellCoord(const int& x, const int& y, const int& z)
{
	return x + y * GRID_X + z * GRID_X * GRID_Y;
}


void ParticleSystem::ComputeForce_SPH(const unsigned int& particle_id)
{
	auto& particle = particle_list[particle_id];
	glm::vec3 force = glm::vec3(0, 0, 0);
	
	glm::vec3 pressureforce = particle.ComputePressureForce_SPH(KERNEL) / particle.GetDensity();
	
	glm::vec3 extforce = glm::vec3(0, 0, 0); // Considered After New Position is Evaluated
	
	glm::vec3 gravity = glm::vec3(0, -GRAVITY, 0);

	glm::vec3 viscosity = particle.ComputeViscosity_SPH(KERNEL) / particle.GetDensity();
	//viscosity = glm::vec3(0.0f, 0.0f, 0.0f);
	glm::vec3 surfacetension = particle.ComputeSurfaceTension_SPH(KERNEL) / particle.GetDensity();
	//surfacetension = glm::vec3(0.0f, 0.0f, 0.0f);
	force = (pressureforce + extforce + gravity + viscosity + surfacetension);
	if (ENABLE_DEBUG_MODE && DISABLE_FORCE) {
		force = glm::vec3(0, 0, 0);
	}
	particle.SetForce(force);
	if (ENABLE_DEBUG_MODE && SHOW_TOTAL_FORCE && (particle_id == WATCH_PARTICLE)) {
		std::cout << "Total Force for Particle[" << particle_id << "] is " << force[0] << "," << force[1] << "," << force[2] << std::endl;
	}
}


void ParticleSystem::ComputeVelocityandPosition_SPH(const unsigned int& particle_id, double timestep)
{
	float ts = (float)timestep;
	auto& particle = particle_list[particle_id];
	
	glm::vec3 force = particle.GetForce();
	glm::vec3 lastVelocity = particle.GetVelocity();
	glm::vec3 lastPosition = particle.GetPos();

	glm::vec3 accerlation = force ;

	/* Leap Frog Scheme */
	glm::vec3 newVelocity = lastVelocity + accerlation * ts; // velocity at half timestep
	glm::vec3 posdiff = newVelocity * ts;
	if ((double)(newVelocity.x / posdiff.x) + (double)(newVelocity.y / posdiff.y) + (double)(newVelocity.z / posdiff.z) * ts > 1.0f) { // CFL Condition
		posdiff.x = newVelocity.x * ts;
		posdiff.y = newVelocity.y * ts;
		posdiff.z = newVelocity.z * ts;
	}
	glm::vec3 newPos = lastPosition + posdiff;
	if (BOUNDARY_MODE == WALL_DAMPING) {
		if (newPos.x >= PARTICLE_INITIAL_BOUNDARY_X) {
			newPos.x = PARTICLE_INITIAL_BOUNDARY_X;
			newVelocity.x = newVelocity.x * WALL_DAMPING_VALUE;
		}
		if (newPos.x < 0.0f) {
			newPos.x = 0.0f;
			newVelocity.x = newVelocity.x * WALL_DAMPING_VALUE;
		}
		if (newPos.y >= PARTICLE_INITIAL_BOUNDARY_Y) {
			newPos.y = PARTICLE_INITIAL_BOUNDARY_Y;
			newVelocity.y = newVelocity.y * WALL_DAMPING_VALUE;
		}
		if (newPos.y < 0.0f) {
			newPos.y = 0.0f;
			newVelocity.y = newVelocity.y * WALL_DAMPING_VALUE;
		}
		if (newPos.z >= PARTICLE_INITIAL_BOUNDARY_Z) {
			newPos.z = PARTICLE_INITIAL_BOUNDARY_Z;
			newVelocity.z = newVelocity.z * WALL_DAMPING_VALUE;
		}
		if (newPos.z < 0.0f) {
			newPos.z = 0.0f;
			newVelocity.z = newVelocity.z * WALL_DAMPING_VALUE;
		}
	}
	

	particle.SetPos(newPos);
	particle.SetVelocity(newVelocity);
	if (ENABLE_DEBUG_MODE && SHOW_POS && (particle_id == WATCH_PARTICLE)) {
		std::cout << "New Position for Particle[" << particle_id << "] is " << newPos[0] << "," << newPos[1] << "," << newPos[2] << std::endl;
	}
	if (ENABLE_DEBUG_MODE && SHOW_VEL && (particle_id == WATCH_PARTICLE)) {
		std::cout << "New Velocity for Particle[" << particle_id << "] is " << newVelocity[0] << "," << newVelocity[1] << "," << newVelocity[2] << std::endl;
	}
}

float ParticleSystem::trace(const glm::mat3& mat)
{
	return mat[0][0] + mat[1][1] + mat[2][2];
}

float ParticleSystem::Fnorm(const glm::mat3& mat)
{
	float fnorm = 0.0f;
	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			fnorm += mat[i][j] * mat[i][j];
		}
	}
	return (float)fnorm;
}

