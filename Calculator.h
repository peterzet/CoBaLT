#ifndef CALCULATOR_H
#define CALCULATOR_H

//#include "ArrayWithTails.h"
#include "raw.h"
#include "awt.h"
#include "aux.h"

#include <boost/math/tools/minima.hpp>

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

//#include "gnuplot-iostream.h"

/////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////           THE SOLVER FOR SPECIFIC PROBLEMS IS DEFINED HERE                      /////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////


// assesment functions
double NormL2(AWT &, AWT &);
void asymmetryTest(rawF &, AWT &, int, int);

// integrating functions
double integrateImP(AWT &);
double integrateReP(AWT &);

double integrateReLeft (AWT &);
double integrateReRight(AWT &);
double integrateImLeft (AWT &);
double integrateImRight(AWT &);

void convolution(AWT &, AWT &, AWT &, aux &, int, int);



#endif // CALCULATOR_H
