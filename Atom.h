#ifndef ATOM_H
#define ATOM_H
#include <armadillo>

using namespace arma;

class Atom
{
private:
    vec position;
    vec velocity;
    vec force;
//    vec potential;
    vec position_0;
    bool canMove;
    bool removed;

public:
    Atom();
    void setPosition(const vec& position);
    vec getPosition();
    void setVelocity(vec velocity);
    vec getVelocity();
    void setForce(vec force);
    vec getForce();
    void setPosInit(vec position_0);
    vec getPosInit();
    void setMove(bool canMove);
    bool getMove();
    void remove(bool removed);
    bool testRemoved();
};

#endif // ATOM_H
