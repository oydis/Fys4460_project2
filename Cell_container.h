#ifndef CELL_CONTAINER_H
#define CELL_CONTAINER_H


#include "Cell.h"


class Cell_container
{
private:
    std::vector<std::vector<std::vector<Cell*> > >  myCells;
public:
    Cell_container(int cells_x);
    void container(int nx, int ny, int nz, Atom* atom);
    std::vector<Atom *> getContainer(int nx, int ny, int nz);
    ~Cell_container();

protected:
    int m_cells_x;
};

#endif // CELL_CONTAINER_H
