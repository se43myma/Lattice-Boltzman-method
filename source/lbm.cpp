#include "../include/lbm.h"
#include "../include/functions.h"
#include "../include/collections.h"

#include "cmath"
#include "iostream"

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

void Set_edgeFlags(Flag_Field & flagField){
    // mark the no-slip boundaries
    for(size_t col=0; col<flagField.sizeX_; ++col){
        flagField(col, 0, 0) = Flags::NoSlipBoundary; // south wall
        flagField(col, flagField.sizeY_-1, 0) = Flags::NoSlipBoundary; // north wall
    }

    // mark the velocity and density boundaries
    for(size_t row=1; row<flagField.sizeY_-1; ++row){
        flagField(0, row, 0) = Flags::VelocityBoundary; // west wall
        flagField(flagField.sizeX_-1, row, 0) = Flags::DensityBoundary; // east wall
    }
}

void Set_initialConditions(Density_Field & densityField, Velocity_Field & velocityField, PDF_Field & pdfField){
    // Set density of cells to 1.0
    std::fill(densityField.data_.begin(), densityField.data_.end(), 1.0);
    // Set velocity of cells to {0,0}
    std::fill(velocityField.data_.begin(), velocityField.data_.end(), 0.0);
    // Set PDFs to equilibrium PDFs
    double density = 1.0;
    std::array<double, 2> initialVelocity {0, 0};
    std::array<double, 9> equilibriumPDFs = Compute_equilibriumPDFs(density, initialVelocity);
    for(size_t row=1; row<pdfField.sizeY_-1; ++row){
        for(size_t col=1; col<pdfField.sizeX_-1; ++col){
            for(uint i=0; i<9; ++i){
                pdfField(col, row, i) = equilibriumPDFs.at(i);
            }
        }
    }
}


inline double Get_pdfFromNeighbour(PDF_Field const & pdfField, std::array<size_t, 2> const & currentPos,
                                   std::array<int, 2> const & neighbourDir, uint const & pdfDir){
    std::array<size_t, 2> const neighbourPos = addVectors<size_t, int>(currentPos, neighbourDir);
    return pdfField(neighbourPos.at(0), neighbourPos.at(1), pdfDir);
}

inline std::array<double, 2> Get_velocityFromNeighbour(Velocity_Field const & velocityField,
                                                       std::array<size_t, 2> const & currentPos,
                                                       std::array<int, 2> const & neighbourDir){
    std::array<size_t, 2> const neighbourPos = addVectors<size_t, int>(currentPos, neighbourDir);
    return {velocityField(neighbourPos.at(0), neighbourPos.at(1), 0),
            velocityField(neighbourPos.at(0), neighbourPos.at(1), 1)};
}

inline std::array<double, 9> Get_cellPDFs(PDF_Field const & pdfField, std::array<size_t, 2> const & currentPos){
    std::array<double, 9> pdfs;
    for(uint i=0; i<9; ++i){
        pdfs[i] = pdfField(currentPos.at(0), currentPos.at(1), i);
    }
    return pdfs;
}

void Setup_noSlipCells(PDF_Field & pdfField, Flag_Field const & flagField,
                       size_t obstaclePosX, size_t obstaclePosY, double obstacleDiameter){
    // Set up the in cells north and south walls
    // For no-slip cells at the corners (SW, NW, NE, SE):
    size_t yMax = pdfField.sizeY_ - 1;
    size_t xMax = pdfField.sizeX_ - 1;
    pdfField(0, 0, Directions::NE) =
            Get_pdfFromNeighbour(pdfField, {0, 0},
                                 directionVectors.at(Directions::NE), Directions::SW);
    pdfField(0, yMax, Directions::SE) =
            Get_pdfFromNeighbour(pdfField, {0, yMax},
                                 directionVectors.at(Directions::SE), Directions::NW);
    pdfField(xMax, yMax, Directions::SW) =
            Get_pdfFromNeighbour(pdfField, {xMax, yMax},
                                 directionVectors.at(Directions::SW), Directions::NE);
    pdfField(xMax, 0, Directions::NW) =
            Get_pdfFromNeighbour(pdfField, {xMax, 0},
                                 directionVectors.at(Directions::NW), Directions::SE);

    // For remaining no-slip cells:
    for(size_t col=1; col<xMax; ++col){
        // south wall ( reflecting SE, S and SW PDFs of the fluid cells )
        pdfField(col, 0, Directions::NW) =
                Get_pdfFromNeighbour(pdfField, {col, 0},
                                     directionVectors.at(Directions::NW), Directions::SE);
        pdfField(col, 0, Directions::N) =
                Get_pdfFromNeighbour(pdfField, {col, 0},
                                     directionVectors.at(Directions::N), Directions::S);
        pdfField(col, 0, Directions::NE) =
                Get_pdfFromNeighbour(pdfField, {col, 0},
                                     directionVectors.at(Directions::NE), Directions::SW);

        // north wall ( reflecting NE, N and NW PDFs of the fluid cells )
        pdfField(col, yMax, Directions::SW) =
                Get_pdfFromNeighbour(pdfField, {col, yMax},
                                     directionVectors.at(Directions::SW), Directions::NE);
        pdfField(col, yMax, Directions::S) =
                Get_pdfFromNeighbour(pdfField, {col, yMax},
                                     directionVectors.at(Directions::S), Directions::N);
        pdfField(col, yMax, Directions::SE) =
                Get_pdfFromNeighbour(pdfField, {col, yMax},
                                     directionVectors.at(Directions::SE), Directions::NW);
    }


    // Set up the cells in obstacle
    size_t limit = std::ceil(0.5*obstacleDiameter);
    // Offset
    size_t posX = obstaclePosX + 1;
    size_t posY = obstaclePosY + 1;
    // For cells in the bounding box :
    for(size_t row=(posY)-limit; row<=(posY)+limit; ++row){
        for(size_t col=(posX)-limit; col<=(posX)+limit; ++col){
            if(flagField(col, row, 0) == Flags::NoSlipBoundary){
                for(uint i=1; i<9; ++i){
                    // directions order : N NE E SE S SW W NW
                    pdfField(col, row, i) = Get_pdfFromNeighbour(pdfField, {col, row},
                                                                 directionVectors.at(i), reflections.at(i));
                }

            }
        }
    }
}

void Setup_velocityCells(PDF_Field & pdfField, std::array<double, 2> const & inletVelocity){
    double contributionFromSW =
            6*weights.at(Directions::SW)*dotProduct(directionVectors.at(Directions::SW), inletVelocity);
    double contributionFromW =
            6*weights.at(Directions::W)*dotProduct(directionVectors.at(Directions::W), inletVelocity);
    double contributionFromNW =
            6*weights.at(Directions::NW)*dotProduct(directionVectors.at(Directions::NW), inletVelocity);

    for(size_t row=1; row<pdfField.sizeY_-1; ++row){
        // Set the NE, E and SE PDFs of the cells
        pdfField(0, row, Directions::NE) =
                Get_pdfFromNeighbour(pdfField, {0, row},
                                     directionVectors.at(Directions::NE), Directions::SW)
                                     - contributionFromSW;
        pdfField(0, row, Directions::E) =
                Get_pdfFromNeighbour(pdfField, {0, row},
                                     directionVectors.at(Directions::E), Directions::W)
                                     - contributionFromW;
        pdfField(0, row, Directions::SE) =
                Get_pdfFromNeighbour(pdfField, {0, row},
                                     directionVectors.at(Directions::SE), Directions::NW)
                                     - contributionFromNW;
    }
}

void Setup_densityCells(PDF_Field & pdfField, Velocity_Field const & velocityField){
    size_t yMax = pdfField.sizeY_ - 1;
    size_t xMax = pdfField.sizeX_ - 1;
    std::array<double,2> cellVelocity;
    for(size_t row=1; row<yMax; ++row){

        // Set the NW, W and SW PDFs of the cells
        cellVelocity = Get_velocityFromNeighbour(velocityField, {xMax, row},
                                                 directionVectors.at(Directions::NW));
        pdfField(xMax, row, Directions::NW) =
                -1*Get_pdfFromNeighbour(pdfField, {xMax, row},
                                        directionVectors.at(Directions::NW), Directions::SE)
                + 2*weights.at(Directions::SE)
                * (1 + 4.5*std::pow(dotProduct(directionVectors.at(Directions::SE), cellVelocity), 2)
                - 1.5* dotProduct(cellVelocity,cellVelocity));

        cellVelocity = Get_velocityFromNeighbour(velocityField, {xMax, row},
                                                 directionVectors.at(Directions::W));
        pdfField(xMax, row, Directions::W) =
                -1*Get_pdfFromNeighbour(pdfField, {xMax, row},
                                        directionVectors.at(Directions::W), Directions::E)
                + 2*weights.at(Directions::E)
                * (1 + 4.5*std::pow(dotProduct(directionVectors.at(Directions::E), cellVelocity), 2)
                - 1.5* dotProduct(cellVelocity,cellVelocity));

        cellVelocity = Get_velocityFromNeighbour(velocityField, {xMax, row},
                                                 directionVectors.at(Directions::SW));
        pdfField(xMax, row, Directions::SW) =
                -1*Get_pdfFromNeighbour(pdfField, {xMax, row},
                                        directionVectors.at(Directions::SW), Directions::NE)
                + 2*weights.at(Directions::NE)
                * (1 + 4.5*std::pow(dotProduct(directionVectors.at(Directions::NE), cellVelocity), 2)
                - 1.5* dotProduct(cellVelocity,cellVelocity));
    }
}

void Perform_streamStep(PDF_Field & source, PDF_Field & destination, Flag_Field & flagField){
    for(size_t row=1; row<source.sizeY_-1; ++row){
        for(size_t col=1; col<source.sizeX_-1; ++col){
            if (flagField(col, row, 0) == Flags::Fluid){
                // perform a pull from neighbours
                for(uint i=0; i<9; ++i){
                    destination(col, row, i) = Get_pdfFromNeighbour(source, {col, row},
                                                                    directionVectors.at(reflections.at(i)),
                                                                    i);
                }
            }
        }
    }
    swap(source, destination);
}

void Perform_collisionStep(PDF_Field & pdfField, Flag_Field const & flagField, Density_Field & densityField,
                           Velocity_Field & velocityField, double const tau){
    double factor = 1.0/tau;
    // update velocity and density lattices after streaming
    Update_densityLattice(densityField, pdfField, flagField);
    Update_velocityLattice(velocityField, pdfField, flagField);
    for(size_t row=1; row<pdfField.sizeY_-1; ++row){
        for(size_t col=1; col<pdfField.sizeX_-1; ++col){
            if (flagField(col, row, 0) == Flags::Fluid){
                std::array<double,2> cellVelocity = {velocityField(col, row, 0), velocityField(col, row, 1)};
                double cellDensity = densityField(col, row, 0);
                std::array<double, 9> equilibriumPDFs =
                        Compute_equilibriumPDFs(cellDensity, cellVelocity);
                // perform collision step
                for(uint i=0; i<9; ++i){
                    pdfField(col, row, i) = pdfField(col, row, i) -
                            factor*(pdfField(col, row, i) - equilibriumPDFs.at(i));
                }
            }
        }
    }
    Update_velocityLattice(velocityField, pdfField, flagField);
}

void Update_velocityLattice(Velocity_Field & velocityField, PDF_Field const & pdfField, Flag_Field const & flagField){
    for(size_t row=1; row<flagField.sizeY_-1; ++row){
        for(size_t col=1; col<flagField.sizeX_-1; ++col){
            if (flagField(col, row, 0) == Flags::Fluid){
                std::array<double,9> cellPDFs = Get_cellPDFs(pdfField, {col, row});
                std::array<double,2> cellVelocity = Compute_cellVelocity(cellPDFs);
                velocityField(col, row, 0) = cellVelocity.at(0);
                velocityField(col, row, 1) = cellVelocity.at(1);
            }
            else{ // no-slip cell
                velocityField(col, row, 0) = 0;
                velocityField(col, row, 1) = 0;
            }
        }
    }
}

void Update_densityLattice(Density_Field & densityField, PDF_Field const & pdfField, Flag_Field const & flagField){
    for(size_t row=1; row<flagField.sizeY_-1; ++row){
        for(size_t col=1; col<flagField.sizeX_-1; ++col){
            if (flagField(col, row, 0) == Flags::Fluid){
                std::array<double,9> cellPDFs = Get_cellPDFs(pdfField, {col, row});
                densityField(col, row, 0) = Compute_cellDensity(cellPDFs);
            }
            else{ // no-slip cell
                densityField(col, row, 0) = 0;
            }
        }
    }
}
