#include "../include/parameters.h"
#include "../include/lattice.h"
#include "../include/lbm.h"
#include "../include/obstacles.h"
#include "../include/collections.h"

// Array for reflected directions
// Reflections.at(Directions::NE) = 6 i.e. SW
std::array<unsigned int, 9> reflections = {0, 5, 6, 7, 8, 1, 2, 3, 4};
// Array for weighting constants used in the collide step
double const centralWeight = 4.0/9;
double const cardinalWeight = 1.0/9;
double const ordinalWeight = 1.0/36;
// directions order : C N NE E SE S SW W NW
std::array<double, 9> weights = {centralWeight, cardinalWeight, ordinalWeight, cardinalWeight, ordinalWeight,
           cardinalWeight, ordinalWeight, cardinalWeight, ordinalWeight};
// Array for direction vectors
// directions order : C N NE E SE S SW W NW
std::array<std::array<int, 2> , 9> directionVectors{{{0,0}, {0,1}, {1,1}, {1,0}, {1,-1}, {0,-1}, {-1,-1}, {-1,0}, {-1,1}}};

#include "iostream"
inline std::array<double, 9> Get_cellPDFs(PDF_Field const & pdfField, std::array<size_t, 2> const & currentPos){
    std::array<double, 9> pdfs;
    for(uint i=0; i<9; ++i){
        pdfs[i] = pdfField(currentPos.at(0), currentPos.at(1), i);
    }
    return pdfs;
}

int main(int argc, char *argv[]) {

    std::string const parameterFilePath = argv[1];

    // Read simulation parameters
    SimulationParameters const parameters(parameterFilePath);
    parameters.Print_simulationParameters();

    // Set up required lattices
    Flag_Field flagField(parameters.sizeX_+2, parameters.sizeY_+2);
    Density_Field densityField(parameters.sizeX_+2, parameters.sizeY_+2);
    Velocity_Field velocityField(parameters.sizeX_+2, parameters.sizeY_+2);
    PDF_Field source(parameters.sizeX_+2, parameters.sizeY_+2);
    PDF_Field destination(parameters.sizeX_+2, parameters.sizeY_+2);
    Set_edgeFlags(flagField);
    Introduce_circularObstacle(flagField, parameters.sphereX_, parameters.sphereY_, parameters.diameter_);
    Set_initialConditions(densityField, velocityField, source);

    // Simulate flow
    size_t vtkCtr = 0;
    for(size_t step=0; step<parameters.timeSteps_; ++ step){
        // Set up helper cells before pull streaming
        Setup_noSlipCells(source, flagField,
                          parameters.sphereX_, parameters.sphereY_, parameters.diameter_);
        Setup_velocityCells(source, {parameters.uIn_, 0});
        Setup_densityCells(source, velocityField);
        // Stream step
        Perform_streamStep(source, destination, flagField);
        //Collision step
        Perform_collisionStep(source, flagField, densityField, velocityField, parameters.tau_);
        if(parameters.vtkStep_ != 0 && step%parameters.vtkStep_ == 0){
            Write_vtkFile(flagField, densityField, velocityField,
                          parameters.vtkFileName_ + std::to_string(vtkCtr));
            ++vtkCtr;
        }
    }
}