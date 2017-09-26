#ifndef AUX_H
#define AUX_H

#include "raw.h"
#include "awt.h"

#include <iostream>         // standard library for reading inputs
#include <fstream>          // standard library for showing outputs
#include <string>           // standard library for manipultaing strings
#include <cassert>          // error handling library, function assert to terminate the program

#include <vector>
#include <complex>          // implements the complex class to contain complex numbers in cartesian form
#include <fftw3.h>          // FFTW library
#include <cmath>            // declares some common mathematical operations and transformation
#include <cstdlib>          // several general purpose functions, including dynamic memory management,
                            // random number generation, communication with the environment,
                            // integer arithmetics, searching, sorting and converting

using namespace std;

class aux
{
    public:
        aux();
        ~aux();

        AWT aux1;
        AWT aux2;
        AWT aux3;
        AWT aux4;
        AWT aux5;
        AWT aux6;

        void initialize_aux(int, double, double);


};

#endif // AUX_H
