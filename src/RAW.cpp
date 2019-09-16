#include "RAW.h"


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




rawF::rawF()
{
    initRawF(0, 0, 0);
}


void rawF::initRawF(int _n, double _xMax, double _kT)
{
    n    = _n;                 // 2_n + 1 is the number of the mesh points
    xMax = _xMax;
    kT = _kT;             // <-xMax, xMax> is the interval on which the function is defined
    f.resize(2*n+1);
    complex<double> u(0,0);
    int i;
    for(i=0; i < 2*n+1; i++)    f[i] = u;

}

rawF::~rawF()
{
    //dtor
}

void rawF::exportFUNasFUN(string name, int div, double range, double mult)
{
    ofstream output(name.c_str());
    int limit=n*range/xMax;

    int i;
    // first the negative values are printed
    for(i = n -limit; i < n + limit; i=i+div)
    {
        output << (-n + i)*xMax/n << "  " << mult*real(f[i]) << "  " << mult*imag(f[i]) << "  " << endl;
    }
}


complex<double> & rawF::operator[](int i)
{
    assert( (i<=n) && (i>=-n) );
    return(f[i >= 0 ? i : 2*n + i]);
}


