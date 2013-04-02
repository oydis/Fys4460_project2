#ifndef CELL_H
#define CELL_H
#include <armadillo>

#include "Atom.h"

//using namespace arma;

class Cell
{
private:
    std::vector<Atom*> myAtoms;
public:
    Cell();
    void setAtoms(Atom *atom);
    std::vector<Atom*> getAtoms();
};

#endif // CELL_H
