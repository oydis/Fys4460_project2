#ifndef LATTICE_H
#define LATTICE_H

#include <iostream>
#include <cstdlib>
#include <iostream>

#include "Atom.h"

using namespace std;

class Lattice
{
protected:
    Atom* atoms;
public:
    Lattice();
    void createLattice(int Nc, double b, int dim, double sigma);
    Atom* getAtoms();

};

#endif // LATTICE_H
