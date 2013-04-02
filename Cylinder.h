#ifndef CYLINDER_H
#define CYLINDER_H
#include <armadillo>
#include "Atom.h"

using namespace arma;

class Cylinder
{
public:
    Cylinder();
    void testCylinder(Atom* atom, double R, double L);
};

#endif // CYLINDER_H
