#ifndef LBM_PARAMETERS_H
#define LBM_PARAMETERS_H

#include "string"

struct SimulationParameters{

    // discretization details
    size_t sizeX_;
    size_t sizeY_;
    size_t timeSteps_;
    // flow details
    double uIn_;
    double Re_;
    double nu_; // kinematic viscosity
    double tau_; // relaxation time
    // obstacle details (circle centered at a grid point)
    size_t sphereX_;
    size_t sphereY_;
    double diameter_;
    // visualization details
    std::string vtkFileName_;
    size_t vtkStep_;

    SimulationParameters(std::string pathToParameterFile);

    void Print_simulationParameters() const;
};

#endif //LBM_PARAMETERS_H
