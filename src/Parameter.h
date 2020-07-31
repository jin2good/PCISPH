#pragma once
/* PARAMETERS */
#define PARTICLE_INITIAL_BOUNDARY_X 0.5f
#define PARTICLE_INITIAL_BOUNDARY_Y 1.0f
#define PARTICLE_INITIAL_BOUNDARY_Z 0.5f

#define KERNEL 0.05f
#define CELL_SIZE KERNEL 
#define INITIAL_DIST KERNEL * 0.5f

#define GRID_X ((int)((PARTICLE_INITIAL_BOUNDARY_X + 4 * INITIAL_DIST) / CELL_SIZE) + 1)
#define GRID_Y ((int)((PARTICLE_INITIAL_BOUNDARY_Y + 4 * INITIAL_DIST) / CELL_SIZE) + 1)
#define GRID_Z ((int)((PARTICLE_INITIAL_BOUNDARY_Z + 4 * INITIAL_DIST) / CELL_SIZE) + 1)

#define NUM_CELL (GRID_X + 1) * (GRID_Y + 1) * (GRID_Z + 1)
#define PARTICLE_COUNT 4000

//Physics Setting
#define dt 0.016f // timestep
#define GRAVITY 0.98f

#define MASS 0.02f         
#define REST_DENSITY 1000.0f
#define GAS_CONSTANT 2.5f
#define VISCOSITY_CONSTANT 3.0f
#define VISCOSITY_THRESHOLD 0.1f
#define ST_THRESHOLD 1.0f
#define TENSION_COEFF 1.0f

//PCISPH Setting
#define MINITERATIONS 3
#define MAXITERATIONS 3
#define DENSITY_FLUCTUATION_THRESHOLD 0.05f
#define PRECOMPUTED_VALUE 0.25f

//GRANULAR Setting
#define ACTIVATE_COHESION true;
#define COHESION 1.0f

//Simulation Settings
#define SPH 1
#define PCISPH 2
#define GRANULAR 3
#define DAM_BREAKING_MODE false
#define ENABLE_DEBUG_MODE true
#define BOUNDARY_MODE WALL_PARTICLE
#define SHOW_TIME true
#define USE_PRECOMPUTED_VALUE true

//Boundary Box
#define WALL_DAMPING 1
#define WALL_PARTICLE 2
#define SHOW_BOX true
#define WALL_DAMPING_VALUE -0.2f
#define WALL_PARTICLE_DENSITY 1.0f
#define SHOW_PARTICLE_BOX false
#define WALL_PARTICLE_MASS 4.0f

//Neighborhood
#define GRID_BASED_NEIGHBORHOOD_SEARCH true
#define SHOW_NEIGHBORSEARCH_RESULT false
#define GET_NEIGHBOR_INFO true

//Watch Particle
#define WATCH_PARTICLE 1000

#define SHOW_DENSITY false
#define SHOW_PUESSURE false

#define DISABLE_FORCE false
#define SHOW_TOTAL_FORCE true
#define SHOW_PRESSUREFORCE true
#define SHOW_NONPUESSUREFORCE false
#define SHOW_VISCOSITY true
#define SHOW_SURFACETENSION false
#define SHOW_FRICTIONCOHESION true
#define SHOW_STRESS true

#define SHOW_VEL false
#define SHOW_POS false
#define SHOW_LOAD_PARTICLE_VECTOR false

// Input Setting
#define ENABLE_KEY_INPUT true
#define SHOW_CAMERA_POS false

