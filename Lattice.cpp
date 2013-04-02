#include "Lattice.h"

Lattice::Lattice()
{
    cout << "A lattice is created" << endl;
}



void Lattice::createLattice(int Nc, double b, int dim, double sigma){

    int i,j,k;
    double x,y,z;
    int counter=0;

    vec r0 (dim);
    vec r1 (dim);
    vec r2 (dim);
    vec r3 (dim);
    vec v0 (dim);
    vec v1 (dim);
    vec v2 (dim);
    vec v3 (dim);

    vec sum (dim);
    sum[0] = 0;
    sum[1] = 0;
    sum[2] = 0;

    int numpart = 4*Nc*Nc*Nc;
    atoms =new Atom[numpart];

    for(i=0; i<Nc; i++){
        for(j=0; j<Nc; j++){
            for(k=0; k<Nc; k++){
                x = i*b;
                y = j*b;
                z = k*b;
                r0 << x << y << z;
                atoms[counter].setPosition(r0);
                atoms[counter].setPosInit(r0);
//                v0 = randn(3)*sigma;
                v0 = randu(3)*2*sigma - sigma;
                sum = sum + v0;
                atoms[counter].setVelocity(v0);
                counter++;

                r1 << x+b/2 << y+b/2 << z;
                atoms[counter].setPosition(r1);
                atoms[counter].setPosInit(r1);
//                v1 = randn(3)*sigma;
                v1 = randu(3)*2*sigma - sigma;
                sum = sum + v1;
                atoms[counter].setVelocity(v1);
                counter++;

                r2 << x << y+b/2 << z+b/2;
                atoms[counter].setPosition(r2);
                atoms[counter].setPosInit(r2);
//                v2 = randn(3)*sigma;
                v2 = randu(3)*2*sigma - sigma;
                sum = sum + v2;
                atoms[counter].setVelocity(v2);
                counter++;

                r3 << x+b/2 << y << z+b/2;
                atoms[counter].setPosition(r3);
                atoms[counter].setPosInit(r3);
//                v3 = randn(3)*sigma;
                v3 = randu(3)*2*sigma - sigma;
                sum = sum + v3;
                atoms[counter].setVelocity(v3);
                counter++;


            }
        }
    }

    vec vel;
    vec vel_new (dim);

    for(i=0;i<numpart;i++){

        vel = atoms[i].getVelocity();

        for(k=0;k<dim;k++){

            vel_new(k) = vel(k) - sum(k)/numpart;

        }

        atoms[i].setVelocity(vel_new);
        atoms[i].setMove(true);
        atoms[i].remove(false);
    }
}

Atom *Lattice::getAtoms()
{
    return atoms;
}

