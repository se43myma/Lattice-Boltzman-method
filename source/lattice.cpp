#include "../include/lattice.h"
#include "../include/collections.h"

void Write_vtkFile(Flag_Field const & flagField, Density_Field const & densityField,
                   Velocity_Field const & velocityField, std::string const & vtkFileName) {

    std::ofstream vtkFile("../outputs/" + vtkFileName + ".vtk");

    // write the header section
    vtkFile << "# vtk DataFile Version 4.0\n";
    vtkFile << "SiWiRVisFile\n";
    vtkFile << "ASCII\n";
    vtkFile << "DATASET STRUCTURED_POINTS\n";
    vtkFile << "DIMENSIONS "<<flagField.sizeX_-2<<" "<<flagField.sizeY_-2<<" 1\n";
    vtkFile << "ORIGIN 0 0 0\n";
    vtkFile << "Spacing 1 1 1\n";
    vtkFile << "Spacing 1 1 1\n";
    vtkFile << "POINT_DATA "<<(flagField.sizeX_-2)*(flagField.sizeY_-2)<<"\n";

    // write the flag field data
    vtkFile << "\nSCALARS flags unsigned_int 1\n"
               "LOOKUP_TABLE default\n";
    for(size_t row=1; row<flagField.sizeY_-1; ++row){
        for(size_t col=1; col<flagField.sizeX_-1; ++col){
            vtkFile << flagField(col, row, 0) << "\n";
        }
    }

    // write the cell density data
    vtkFile << "\nSCALARS density double 1\n"
               "LOOKUP_TABLE default\n";
    for(size_t row=1; row<densityField.sizeY_-1; ++row){
        for(size_t col=1; col<densityField.sizeX_-1; ++col){
            vtkFile << densityField(col, row, 0) << "\n";
        }
    }

    // write the cell velocity data
    vtkFile << "\nVECTORS velocity double\n";
    for(size_t row=1; row<velocityField.sizeY_-1; ++row){
        for(size_t col=1; col<velocityField.sizeX_-1; ++col){
            vtkFile << velocityField(col, row, 0) << " " << velocityField(col, row, 1) << " 0\n";
        }
    }

    vtkFile.close();
}
