#include "Calculator.h"
//#include "physics.h"

#include <boost/math/tools/minima.hpp>

#include <vector>
#include <iostream>         // standard library for reading inputs
#include <fstream>          // standard library for showing outputs
#include <string>           // standard library for manipultaing strings
#include <cassert>          // error handling library, function assert to terminate the program

#include <complex> // implements the complex class to contain complex numbers in cartesian form
#include <fftw3.h>          // FFTW library
#include <cmath>            // declares some common mathematical operations and transformation
#include <cstdlib>          // several general purpose functions, including dynamic memory management,
                            // random number generation, communication with the environment,
                            // integer arithmetics, searching, sorting and converting

//#include "gnuplot-iostream.h"

using namespace std;

// general functions
double NormL2(AWT & G, AWT & Test)
{
    if(G.n != Test.n)   {cerr << "NormL2: AWTs are not compatible!" << endl;   exit(1); }

    double L2norm = 0.0;
    double L2diff = 0.0;

    for(int i=0; i<=G.n+1; i++)
    {
        L2norm += abs(G.y[i] * Test.y[i]);
        L2diff += abs(( G.y[i] - Test.y[i] ) * ( G.y[i] - Test.y[i] )  );
    }
    for(int i=3*G.n+4; i<=4*G.n+4; i++)
    {
        L2norm += abs(G.y[i] * Test.y[i]);
        L2diff += abs(( G.y[i] - Test.y[i] ) * ( G.y[i] - Test.y[i] )  );
    }

    return(sqrt(L2diff/L2norm));

}

// simple integration
complex<double> Trapezoid(AWT & X)
{
    // trapezoidal rule is now applied, later changing to Simpson
    double real = Trapezoid_Re(X);
    double imag = Trapezoid_Im(X);
    complex<double> u(0,1);
    complex<double> integral = real + u*imag;

    return integral;
}

// simple integration
double Trapezoid_Re(AWT & X)
{
    // trapezoidal rule is now applied, later changing to Simpson
    double dx = X.xMax/X.n;
    double integral = 0.5*real(X.y[X.n]) * dx + 0.5*real(X.y[4*X.n+3])*dx ;

    for(int i=0; i<X.n; i++) integral = integral + real(X.y[i]) * dx;
    for(int i=3*X.n+4; i<4*X.n+4; i++) integral = integral + real(X.y[i]) * dx;

    return integral;
}

double Trapezoid_Im(AWT & X)
{
    // trapezoidal rule is now applied, later changing to Simpson
    double dx = X.xMax/X.n;
    double integral = 0.5*imag(X.y[X.n]) * dx + 0.5*imag(X.y[4*X.n+3])*dx ;

    for(int i=0; i<X.n; i++) integral = integral + imag(X.y[i]) * dx;
    for(int i=3*X.n+4; i<4*X.n+4; i++) integral = integral + imag(X.y[i]) * dx;

    return integral;
}


// simple integration
complex<double> Simpson(AWT & X)
{
    // trapezoidal rule is now applied, later changing to Simpson
    double real = Simpson_Re(X);
    double imag = Simpson_Im(X);
    complex<double> u(0,1);
    complex<double> integral = real + u*imag;

    return integral;
}


// Simpson's rule integrals
double Simpson_Re(AWT & X)
{
    double dx = X.xMax/X.n;

    // correction for edge points
    double integral = 1.0*real(X.y[X.n]) + 1.0*real(X.y[4*X.n+3]);


    int point_iter = 0;
    int simpson_coeff=2.0;

    // negative part of the interval
    for(int i=3*X.n+4; i<4*X.n+4; i++)
    {
        if( (point_iter %2) == 1)   simpson_coeff = 2.0;
        else                        simpson_coeff = 4.0;

        integral = integral + simpson_coeff * real(X.y[i]);
        point_iter = point_iter + 1;
    }

    // positive part
    for(int i=0; i<X.n; i++)
    {
        if( (point_iter %2) == 1)   simpson_coeff = 2.0;
        else                        simpson_coeff = 4.0;

        integral = integral + simpson_coeff * real(X.y[i]);
        point_iter = point_iter + 1;
    }

    integral = integral * dx / 3.0;

    return integral;

}

// Simpson's rule integrals
double Simpson_Im(AWT & X)
{
    double dx = X.xMax/X.n;

    // correction for edge points
    double integral = +1.0*imag(X.y[X.n]) + 1.0*imag(X.y[4*X.n+3]);

    int point_iter = 0;
    int simpson_coeff=2.0;

    // negative part of the interval
    for(int i=3*X.n+4; i<4*X.n+4; i++)
    {
        if( (point_iter %2) == 1)   simpson_coeff = 2.0;
        else                        simpson_coeff = 4.0;

        integral = integral + simpson_coeff * imag(X.y[i]);
        point_iter = point_iter + 1;
    }

    // positive part
    for(int i=0; i<X.n; i++)
    {
        if( (point_iter %2) == 1)   simpson_coeff = 2.0;
        else                        simpson_coeff = 4.0;

        integral = integral + simpson_coeff * imag(X.y[i]);
        point_iter = point_iter + 1;
    }

    integral = integral * dx / 3.0;

    return integral;

}
