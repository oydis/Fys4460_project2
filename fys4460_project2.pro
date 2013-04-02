TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += main.cpp \
    Integrator.cpp \
    Lattice.cpp \
    Cell_container.cpp \
    Cell.cpp \
    Atom.cpp \
    Cylinder.cpp \
    Spheres.cpp

HEADERS += \
    Integrator.h \
    Lattice.h \
    Cell_container.h \
    Cell.h \
    Atom.h \
    Cylinder.h \
    Spheres.h

release {
    DEFINES += ARMA_NO_DEBUG
    QMAKE_CXXFLAGS_RELEASE -= -O2
    QMAKE_CXXFLAGS_RELEASE += -O3
}
