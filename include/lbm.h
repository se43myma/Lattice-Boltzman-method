#ifndef LBM_LBM_H
#define LBM_LBM_H

#include "lattice.h"

void Set_edgeFlags(Flag_Field & flagField);

void Set_initialConditions(Density_Field & densityField, Velocity_Field & velocityField, PDF_Field & pdfField);

// Function to get PDF from neighbour
// neighbourDir : direction vector pointing towards the neighbour
// pdfDir : direction of the required pdf in neighbouring cell
inline double Get_pdfFromNeighbour(PDF_Field const & pdfField, std::array<size_t, 2> const & currentPos,
                                   std::array<int, 2> const & neighbourDir, uint const & pdfDir);

inline std::array<double, 2> Get_velocityFromNeighbour(Velocity_Field const & velocityField,
                                                       std::array<size_t, 2> const & currentPos,
                                                       std::array<int, 2> const & neighbourDir);

inline std::array<double, 9> Get_cellPDFs(PDF_Field const & pdfField, std::array<size_t, 2> const & currentPos);

// Function for setting up no-slip boundary cells before pull streaming
void Setup_noSlipCells(PDF_Field & pdfField, Flag_Field const & flagField,
                       size_t obstaclePosX, size_t obstaclePosY, double obstacleDiameter);

// Function for setting up velocity boundary cells before pull streaming
void Setup_velocityCells(PDF_Field & pdfField, std::array<double, 2> const & inletVelocity);

// Function for setting up density boundary cells before pull streaming
void Setup_densityCells(PDF_Field & pdfField, Velocity_Field const & velocityField);

// Function for updating velocity lattice
void Update_velocityLattice(Velocity_Field & velocityField, PDF_Field const & pdfField, Flag_Field const & flagField);

// Function for updating density lattice
void Update_densityLattice(Density_Field & densityField, PDF_Field const & pdfField, Flag_Field const & flagField);

// Function for implementing stream step (pull)
void Perform_streamStep(PDF_Field & source, PDF_Field & destination, Flag_Field & flagField);

// Function for implementing collision step and computing density, velocity lattices
void Perform_collisionStep(PDF_Field & pdfField, Flag_Field const & flagField, Density_Field & densityField,
                           Velocity_Field & velocityField, double const tau);

#endif //LBM_LBM_H
