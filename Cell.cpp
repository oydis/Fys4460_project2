#include "Cell.h"

Cell::Cell()
{
}

void Cell::setAtoms(Atom* atom){

    myAtoms.push_back(atom);

}

std::vector<Atom *> Cell::getAtoms(){
    return myAtoms;
}
