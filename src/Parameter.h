#pragma once
/* PARAMETERS */
#define PARTICLE_INITIAL_BOUNDARY_X 0.64f
#define PARTICLE_INITIAL_BOUNDARY_Y 0.32f
#define PARTICLE_INITIAL_BOUNDARY_Z 0.64f

#define KERNEL 0.05f
#define CELL_SIZE KERNEL 
#define INITIAL_DIST KERNEL * 0.5f

#define GRID_X ((int)((PARTICLE_INITIAL_BOUNDARY_X + 4 * INITIAL_DIST) / CELL_SIZE) + 1)
#define GRID_Y ((int)((PARTICLE_INITIAL_BOUNDARY_Y + 4 * INITIAL_DIST) / CELL_SIZE) + 1)
#define GRID_Z ((int)((PARTICLE_INITIAL_BOUNDARY_Z + 4 * INITIAL_DIST) / CELL_SIZE) + 1)

#define NUM_CELL (GRID_X + 1) * (GRID_Y + 1) * (GRID_Z + 1)
#define PARTICLE_COUNT 5000

#define GRAVITY 0.98f

#define MASS 0.02f         
#define REST_DENSITY 4000.0f
#define GAS_CONSTANT 0.3f
#define VISCOSITY_CONSTANT 1.5f
#define VISCOSITY_THRESHOLD 0.1f
#define ST_THRESHOLD 6.0f
#define TENSION_COEFF 1.0f

//PCISPH Setting
#define MINITERATIONS 3
#define DENSITY_FLUCTUATION_THRESHOLD 0.05f
#define PRECOMPUTED_VALUE 1.0f

// Input Setting
#define ENABLE_KEY_INPUT true
#define SHOW_CAMERA_POS false

//Simulation Settings
#define SPH 1
#define PCISPH 2

//Initial Setting
#define DAM_BREAKING_MODE false

#define dt 0.05f // timestep

#define ENABLE_DEBUG_MODE true

//Boundary Box
#define WALL_DAMPING 1
#define WALL_PARTICLE 2
#define SHOW_BOX true
#define BOUNDARY_MODE WALL_PARTICLE
#define WALL_DAMPING_VALUE -0.5f
#define WALL_PARTICLE_DENSITY 2.0f
#define SHOW_PARTICLE_BOX false

//Neighborhood
#define GRID_BASED_NEIGHBORHOOD_SEARCH true
#define SHOW_NEIGHBORSEARCH_RESULT false
#define GET_NEIGHBOR_INFO false

//Watch Particle
#define WATCH_PARTICLE 1000

#define SHOW_DENSITY true
#define SHOW_PUESSURE true

#define DISABLE_FORCE false
#define SHOW_TOTAL_FORCE true
#define SHOW_PRESSUREFORCE true
#define SHOW_VISCOSITY false
#define SHOW_SURFACETENSION true

#define SHOW_VEL false
#define SHOW_POS true
#define SHOW_LOAD_PARTICLE_VECTOR false


