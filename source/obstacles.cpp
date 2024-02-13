#include "../include/obstacles.h"
#include "../include/collections.h"

#include "assert.h"
#include "cmath"

void Introduce_circularObstacle(Flag_Field & flagField, size_t posX, size_t posY, double diameter){
    // check if obstacle fits in the flow channel
    assert(posX + 0.5*diameter < flagField.sizeX_ && posX - 0.5*diameter > 0 &&
           posY + 0.5*diameter < flagField.sizeY_ && posY - 0.5*diameter > 0);

    // find all the lattice points in bounding box (defined by the radius of sphere) that lie in the circle.
    size_t limit = std::ceil(0.5*diameter);
    double radiusSq = 0.25*diameter*diameter;
    // offset
    posX = posX + 1;
    posY = posY + 1;
    // For cells in the bounding box :
    for(size_t row=(posY)-limit; row<=(posY)+limit; ++row){
        for(size_t col=(posX)-limit; col<=(posX)+limit; ++col){
            double distanceSq = (posY - row)*(posY - row) + (posX - col)*(posX - col);
            if(distanceSq <= radiusSq){
                // set the flags of lattice points to NoSlipBoundary
                flagField(col,row,0) = Flags::NoSlipBoundary;
            }
        }
    }
}