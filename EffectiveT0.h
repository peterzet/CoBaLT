#ifndef EFFECTIVET0_H
#define EFFECTIVET0_H

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
void Effective_zeroT(string);



void Lambda_zeroT(double, AWT &, AWT &, double &, AWT &, double &, AWT &, AWT &, AWT &, AWT &, AWT &, AWT &, AWT &, AWT &, AWT &);

double Brent_Lambda_zeroT(double, AWT &, AWT &, AWT &, AWT &, AWT &, AWT &, AWT &, AWT &, AWT &, AWT &, AWT &);



pair<double, double> Brent_nT_zeroT(double, double, double, AWT &, AWT &, AWT &);
pair<double, double> Brent_nS_zeroT(double, double, double, AWT &, AWT &, AWT &, AWT &, AWT &);


void drawLambda(double, rawF &, AWT &, AWT &, AWT &, AWT &, AWT &);


// PARQUETS LOW FREQUENCY APPROX
void SolverParquetLowFreq(double, AWT &, AWT &, AWT &, AWT &, AWT &);
void ParquetLowFreq(int, double, int, double, double, double, AWT &, AWT &, AWT &, AWT &, AWT &);

// FLEX
void SolverFLEX(double, double, AWT &, AWT &, AWT&, AWT &, AWT &, AWT &);
void Flex(double, double, int, int, double, double, double, double, AWT &, AWT &, AWT &);

#endif // METHODS_H
