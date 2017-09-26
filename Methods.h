#ifndef METHODS_H
#define METHODS_H

#include "Calculator.h"

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
void Effective_nonZeroT(string);
void Lambda_nonZeroT(double, double, AWT &, AWT &, AWT &, AWT &, double &, double &, AWT &, AWT &, AWT &);

pair<double, double> Brent_nT_nonZeroT(double, double, double, AWT &, AWT &, AWT &);
pair<double, double> Brent_nS_nonZeroT(double, double, double, AWT &, AWT &, AWT &, AWT &, AWT &);

#endif // METHODS_H
