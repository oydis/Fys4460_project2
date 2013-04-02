#include "Atom.h"

// test
Atom::Atom()
{

}

vec Atom::getPosition()
{
    return position;
}

void Atom::setPosition(const vec& pos){

    position = pos;

}

vec Atom::getVelocity(){
    return velocity;
}

void Atom::setVelocity(vec vel){

    velocity = vel;
}

vec Atom::getForce(){
    return force;
}

void Atom::setForce(vec f){
    force = f;
}

vec Atom::getPosInit()
{
    return position_0;
}

void Atom::setPosInit(vec pos0){

    position_0 = pos0;

}

bool Atom::getMove(){
    return canMove;
}

void Atom::setMove(bool cM){
    canMove = cM;
}

bool Atom::testRemoved(){
    return removed;
}

void Atom::remove(bool rm){
    removed = rm;
}


