#include "aux.h"


#include <iostream>         // standard library for reading inputs
#include <fstream>          // standard library for showing outputs
#include <string>           // standard library for manipultaing strings
#include <cassert>          // error handling library, function assert to terminate the program

#include <vector>
#include <complex> // implements the complex class to contain complex numbers in cartesian form
#include <fftw3.h>          // FFTW library
#include <cmath>            // declares some common mathematical operations and transformation
#include <cstdlib>          // several general purpose functions, including dynamic memory management,
                            // random number generation, communication with the environment,
                            // integer arithmetics, searching, sorting and converting




aux::aux()
{
    initialize_aux(0, 0, 0);
}

void aux::initialize_aux(int n, double xMax, double kT)
{
    aux1.initializeAWT(n, xMax, kT);
    aux2.initializeAWT(n, xMax, kT);
    aux3.initializeAWT(n, xMax, kT);
    aux4.initializeAWT(n, xMax, kT);
    aux5.initializeAWT(n, xMax, kT);
    aux6.initializeAWT(n, xMax, kT);

    aux1.set_zero();
    aux2.set_zero();
    aux3.set_zero();
    aux4.set_zero();
    aux5.set_zero();
    aux6.set_zero();
}

aux::~aux()
{

}

