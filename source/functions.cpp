#include "../include/functions.h"
#include "../include/collections.h"

inline std::array<double, 2> scaleVector(std::array<int, 2> const a, double const scaleFactor){
    return {scaleFactor*a.at(0), scaleFactor*a.at(1)};
}

inline double Compute_cellDensity(std::array<double, 9> & cellPDFs){
    double density = 0;
    for(double pdf : cellPDFs){
        density += pdf;
    }
    return density;
}

inline std::array<double, 2> Compute_cellVelocity(std::array<double, 9> & cellPDFs){
    std::array<double, 2> velocity={0, 0};
    for(int i=0; i<9; ++i){
        velocity = addVectors<double, double>(velocity, scaleVector(directionVectors.at(i), cellPDFs.at(i)));
    }
    return velocity;
}

std::array<double, 9> Compute_equilibriumPDFs(double & cellDensity, std::array<double, 2> & cellVelocity){
    // directions order : C N NE E SE S SW W NW
    std::array<double, 9> equilibriumPDFs {0, 0, 0, 0, 0, 0, 0, 0, 0};
    double velocitySq = dotProduct<double, double>(cellVelocity, cellVelocity);
    for(int i=0; i<9; ++i){
        double velocityComp = dotProduct<int, double>(directionVectors.at(i), cellVelocity);
        equilibriumPDFs[i] = weights.at(i)*(cellDensity + 3*velocityComp + 4.5*velocityComp*velocityComp
                             - 1.5*velocitySq);
    }
    return  equilibriumPDFs;
}