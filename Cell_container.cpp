#include "Cell_container.h"

Cell_container::Cell_container(int cells_x) :
    m_cells_x(cells_x)
{


    myCells.resize(cells_x);
    for (int i = 0; i < cells_x; i++) {
      myCells[i].resize(cells_x);

      for (int j = 0; j < cells_x; j++){
        myCells[i][j].resize(cells_x);

      }

    }
    for(int i=0;i<cells_x;i++){
        for(int j=0;j<cells_x;j++){
            for(int k=0;k<cells_x;k++){
                myCells[i][j][k] =(new Cell());
            }
        }
    }
}

Cell_container::~Cell_container() {

    for(int i=0;i<m_cells_x;i++){
        for(int j=0;j<m_cells_x;j++){
            for(int k=0;k<m_cells_x;k++){
                delete myCells[i][j][k];
            }
        }
    }

}

void Cell_container::container(int nx, int ny, int nz, Atom *atom){

    myCells[nx][ny][nz]->setAtoms(atom);

}

std::vector<Atom *> Cell_container::getContainer(int nx, int ny, int nz){
    return myCells[nx][ny][nz]->getAtoms();
}
