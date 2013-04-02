#include <iostream>
#include <armadillo>
#include "Lattice.h"
#include "Atom.h"
#include "Integrator.h"

using namespace std;
using namespace arma;
int main()
{
    int Nc = 20;
    int dim = 3;
//    double b = 1.545; // = 5.260;
    double b = 1.68; // 5.72
    int numpart = 4*Nc*Nc*Nc;

    double T_bath = 300/119.74; // K
    double eps = 1.0318e-2; // eV
    double m = 39.948; // u
    double v0 = sqrt(T_bath);

    Lattice Lat;
    Lat.createLattice(Nc,b,dim,v0);
    Atom* atoms = Lat.getAtoms();

    double dt = 0.006; //0.006;
    double sig = 3.405; //Å


    Integrator Int(dt,m,numpart,dim, atoms,b,Nc,sig,eps,T_bath);
    Int.integrate();

}
