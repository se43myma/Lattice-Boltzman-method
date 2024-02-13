#include "../include/parameters.h"

#include "fstream"
#include "iostream"

SimulationParameters::SimulationParameters(std::string pathToParameterFile) {
    std::ifstream parameterFile (pathToParameterFile);
    std::string key, value[10];

    for (int i=0; i<10; ++i){
        parameterFile>>key;
        parameterFile>>value[i];
    }

    // discretization details
    sizeX_ = std::stoull(value[0]);
    sizeY_ = std::stoull(value[1]);
    timeSteps_ = std::stoull(value[2]);

    // flow details
    uIn_ = std::stod(value[3]);
    Re_ = std::stod(value[4]);
    nu_ =  double(uIn_*sizeY_)/Re_;
    tau_ = 3 * nu_ + 0.5;

    // obstacle details
    sphereX_ = std::stoull(value[5]);
    sphereY_ = std::stoull(value[6]);
    diameter_ = std::stod(value[7]);

    // visualization details
    vtkFileName_ = value[8];
    vtkStep_ = std::stoull(value[9]);

    parameterFile.close();
}

void SimulationParameters::Print_simulationParameters() const {
    std::cout<<"\n\n*************************\n";
    std::cout<<"* Simulation Parameters *\n";
    std::cout<<"*************************\n\n";
    std::cout<<"sizex: "<<sizeX_<<"\n";
    std::cout<<"sizey: "<<sizeY_<<"\n";
    std::cout<<"timesteps: "<<timeSteps_<<"\n";
    std::cout<<"uin: "<<uIn_<<"\n";
    std::cout<<"Re: "<<Re_<<"\n";
    std::cout<<"nu (kinematic viscosity): "<<nu_<<"\n";
    std::cout<<"tau (relaxation time): "<<tau_<<"\n";
    std::cout<<"spherex: "<<sphereX_<<"\n";
    std::cout<<"spherey: "<<sphereY_<<"\n";
    std::cout<<"diameter: "<<diameter_<<"\n";
    std::cout<<"vtk_file: "<<vtkFileName_<<"\n";
    std::cout<<"vtk_step: "<<vtkStep_<<"\n\n";
}