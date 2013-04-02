#include "Spheres.h"

Spheres::Spheres(int spheres, double L) :
    spheres(spheres),
    L(L)
{
}

void Spheres::findSpheres(Atom* atom, vec sphere_pos_x, vec sphere_pos_y, vec sphere_pos_z, vec sphere_radius){
    int i,q;
    int dim = 3;
    double R_test;
    int test = 0;
    vec d (dim);
    vec r_min (dim);
    vec v_new (dim);
    vec r = atom->getPosition();

    for(i=0;i<spheres;i++){  // spheres

        d(0) = r(0) - sphere_pos_x(i);
        d(1) = r(1) - sphere_pos_y(i);
        d(2) = r(2) - sphere_pos_z(i);

        for(q=0;q<dim;q++){   // minimum image convention
            if(fabs(d(q)-L) < fabs(d(q))){
               if(fabs(d(q)-L) < fabs(d(q)+L)){
                  r_min(q) = d(q)-L;
               }
            }
            if(fabs(d(q)) < fabs(d(q)+L)){
               if(fabs(d(q)) < fabs(d(q)-L)){
                  r_min(q) = d(q);
               }
            }
            if(fabs(d(q)+L) < fabs(d(q)-L)){
               if(fabs(d(q)+L) < fabs(d(q))){
                  r_min(q) = d(q)+L;
               }
            }
        }
        R_test = sqrt(r_min(0)*r_min(0) + r_min(1)*r_min(1) + r_min(2)*r_min(2));

        if(R_test < sphere_radius(i)){
            test = 1;
        }
    }

    if(test == 1){
        atom->setMove(false);
        v_new(0) = 0.0;
        v_new(1) = 0.0;
        v_new(2) = 0.0;
        atom->setVelocity(v_new);
    }
    else{
        atom->setMove(true);
    }
}

double Spheres::findSphereVolume(double pi, vec sphere_radius){
    double sphere_volume = 0;
    for(int i=0;i<spheres;i++){
        //cout << "rad" << sphere_radius(i) << endl;
        sphere_volume = sphere_volume + 4.0/3.0*pi*pow (sphere_radius(i),3);
    }
    return sphere_volume;
}
