#ifndef LBM_FUNCTIONS_H
#define LBM_FUNCTIONS_H

#include "array"

// Function for computing the sum of 2D vectors
template <typename Type1, typename Type2>
inline std::array<Type1, 2> addVectors(std::array<Type1, 2> const a, std::array<Type2, 2> const b){
    return {a.at(0)+b.at(0), a.at(1)+b.at(1)};
}
// Function for computing dot product of 2D vectors
template <typename Type1, typename Type2>
inline double dotProduct(std::array<Type1, 2> const a, std::array<Type2, 2> const b){
    return a.at(0)*b.at(0) + a.at(1)*b.at(1);
}

// Function for scaling a 2D vector
inline std::array<double, 2> scaleVector(std::array<int, 2> const a, double const scaleFactor);

// Function for computing cell density
inline double Compute_cellDensity(std::array<double, 9> & cellPDFs);

// Function for computing cell velocity
inline std::array<double, 2> Compute_cellVelocity(std::array<double, 9> & cellPDFs);

// Function for computing equilibrium PDFs of a cell
std::array<double, 9> Compute_equilibriumPDFs(double & cellDensity, std::array<double, 2> & cellVelocity);

#endif //LBM_FUNCTIONS_H
