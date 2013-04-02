#ifndef INTEGRATOR_H
#define INTEGRATOR_H

#include <iostream>
#include <cstdlib>
//#include <list>
#include <cmath>

#include "Atom.h"
#include "Cell_container.h"
#include "Cylinder.h"
#include "Spheres.h"

using namespace std;

class Integrator
{
protected:
    Atom* atoms;
    //Atom* atoms_dt;
    double dt;
    double m;
    int numpart;
    int dim;
    double b;
    int Nc;
    double sig;
    double eps;
    double T_bath;

public:
    Integrator(double dt, double m, int numpart, int dim, Atom *atoms, double b, int Nc, double sig, double eps, double T_bath);
    void integrate();


};

#endif // INTEGRATOR_H
