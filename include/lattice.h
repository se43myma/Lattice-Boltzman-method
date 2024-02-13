#ifndef LBM_LATTICE_H
#define LBM_LATTICE_H

#include "array"
#include "assert.h"
#include "fstream"
#include "string"
#include "vector"

using uint = unsigned int;

template <typename Type, uint cellSize>
struct Lattice{
    size_t sizeX_;
    size_t sizeY_;
    std :: vector <Type> data_; // Contiguous 1D representation of the 2D data grid

    inline Lattice (size_t sizeX , size_t sizeY);

    // Access operator
    // fth entry of the data for a cell located at x, y
    inline Type & operator ()(size_t x, size_t y, uint f); // read and write
    inline Type operator ()(size_t x, size_t y, uint f) const; // read only

    // Swapping member function
    void swap(Lattice & lattice);

};

// Lattice type definitions
using PDF_Field = Lattice<double, 9>;
using Velocity_Field = Lattice<double, 2>;
using Density_Field = Lattice<double, 1>;
using Flag_Field = Lattice<uint, 1>;

template <typename Type, uint cellSize>
Lattice<Type, cellSize>::Lattice(size_t sizeX , size_t sizeY)
        : sizeX_ (sizeX)
        , sizeY_ (sizeY)
        , data_ (std::vector<Type>(cellSize * sizeX * sizeY)){}

// Implementation of the access operator
template <typename Type, uint cellSize>
inline Type & Lattice<Type, cellSize>::operator()(size_t x, size_t y, uint f){
    assert (x < sizeX_ && y < sizeY_ && f < cellSize);
    return data_[y * sizeX_ * cellSize + x * cellSize + f];
}

template <typename Type, uint cellSize>
inline Type Lattice<Type, cellSize>::operator ()(size_t x, size_t y, uint f) const{
    assert (x < sizeX_ && y < sizeY_ && f < cellSize);
    return data_[y * sizeX_ * cellSize + x * cellSize + f];
}

// Implementation of the swapping member function
template <typename Type, uint cellSize>
void Lattice<Type, cellSize>::swap(Lattice & lattice){
    std::swap(sizeX_, lattice.sizeX_);
    std::swap(sizeY_, lattice.sizeY_);
    std::swap(data_, lattice.data_);
}

// Global swap function
template < typename Type , uint cellSize >
inline void swap(Lattice<Type, cellSize> & a, Lattice<Type, cellSize> & b){
    a.swap(b);
}

// Function to write VTK files
void Write_vtkFile(Flag_Field const & flagField, Density_Field const & densityField,
                   Velocity_Field const & velocityField, std::string const & vtkFileName);

#endif //LBM_LATTICE_H
