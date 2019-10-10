////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//          Calculates mean-field theory of frequency dependent vertices according to Janis may 2018
//                                  Vertices are momentum independent
////////////////////////////////////////////////////////////////////////////////////////////////////////////////


#ifndef EFFECTIVE_DYNSPINMAT_H
#define EFFECTIVE_DYNSPINMAT_H


#include "../Calculator.h"
#include "../global.h"


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




// PARQUETS EFFECTIVE INTERACTION APPROX
void static_spinpolarized_SSN(string);
void check_transformation(string);


#endif // EFFECTIVE_DYNSPINMAT_H
