#include "Global.h"

#include "Particle.h"
#include "VertexArray.h"
#include "VertexBuffer.h"
#include "Parameter.h"



class ParticleSystem {
public:
	int particle_count = 0;
	clock_t sphStart, sphEnd;
	double sphtime = 0;

	clock_t NeighborsearchStart, NeighborsearchEnd;
	double neighborsearchtime = 0;

	clock_t EstimatePressureStart, EstimatePressureEnd;
	double estimatepressuretime = 0;

	clock_t PredictNewPosVelStart, PredictNewPosVelEnd;
	double estimatePosVeltime = 0;
private:
	int SPHmode;
	std::vector<Particle> particle_list;
	std::vector<float> particle_location;

	
	std::vector<std::list<const Particle*>> Cell;
	double timestep = 0;
public:
	ParticleSystem(int mode);
	~ParticleSystem();

	void Init();
	void Update(double deltaTime);
    void Update_SPH(double deltaTime);
	void Update_PCISPH(double deltaTime);
	
	float* GetParticlePositionArray();
private:
	void SetInitialParticlePosition();
	void LoadParticleVectorPosition();

	void FindNeighbors();
	void ComputeDensity();
	void ComputePressure();

	void ComputeForce();

	void ComputeVelocityandPosition(double timestep);
	
	/* For Grid Based Search */
	void CellInsert();
	glm::vec3 GetCellIndex(const glm::vec3& pos);
	int CellCoord(const int& x, const int& y, const int& z);

	/* Naive Neighbor Search*/
	void NaiveNeighborSearch(const unsigned int& particle_id);
	
	/* For SPH */
	void ComputeForce_SPH(const unsigned int& particle_id);

	void ComputeVelocityandPosition_SPH(const unsigned int& particle_id, double timestep);

	/* For PCISPH */

	

	
	
};