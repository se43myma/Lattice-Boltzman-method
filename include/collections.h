#ifndef LBM_COLLECTIONS_H
#define LBM_COLLECTIONS_H

#include "array"

// Enum for cell flags
enum Flags {Fluid, NoSlipBoundary, VelocityBoundary, DensityBoundary};

// Enum for directions
enum Directions {C, N, NE, E, SE, S, SW, W, NW};

// Extern variables
extern std::array<unsigned int, 9> reflections;
extern std::array<double, 9> weights;
extern std::array<std::array<int, 2> , 9> directionVectors;

#endif //LBM_COLLECTIONS_H
