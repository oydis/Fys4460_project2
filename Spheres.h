#ifndef SPHERES_H
#define SPHERES_H
#include <armadillo>
#include "Atom.h"

using namespace arma;

class Spheres
{
private:
    int spheres;
    double L;
//    vec sphere_pos_x;
//    vec sphere_pos_y;
//    vec sphere_pos_z;
//    vec sphere_radius;

public:
    Spheres(int spheres, double L);
    void findSpheres(Atom* atom, vec sphere_pos_x, vec sphere_pos_y, vec sphere_pos_z, vec sphere_radius);
    double findSphereVolume(double pi, vec sphere_radius);
};

#endif // SPHERES_H
