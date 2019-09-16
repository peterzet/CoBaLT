#ifndef RAWFUNCTION_H
#define RAWFUNCTION_H


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

class rawF
{
    public:
        rawF();
        ~rawF();
        void initRawF(int, double, double);

        int n;
        double xMax;
        double kT;
        vector<complex <double> >  f;      // array of function values

        void exportFUNasFUN(string, int, double, double);
        complex<double>& operator[](int);

};

#endif // RAWFUNCTION_H
