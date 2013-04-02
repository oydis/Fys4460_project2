#include "Cylinder.h"

Cylinder::Cylinder()
{
}

void Cylinder::testCylinder(Atom* atom, double R, double L){

    vec v_new (3);
    vec r = atom->getPosition();

    double R_test = sqrt(pow (r(0)-L/2,2) + pow (r(1)-L/2,2));  // cylinder

    if(R_test < R){
        atom->setMove(true);
    }
    else{
        atom->setMove(false);
        v_new(0) = 0.0;
        v_new(1) = 0.0;
        v_new(2) = 0.0;
        atom->setVelocity(v_new);
    }
}
