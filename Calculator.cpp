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
double integrateReP(AWT & X)
{
    // trapezoidal rule is now applied, later changing to Simpson
    double dx = X.xMax/X.n;
    double integral = 0.5*real(X.y[X.n]) * dx + 0.5*real(X.y[4*X.n+3])*dx ;

    for(int i=0; i<X.n; i++) integral = integral + real(X.y[i]) * dx;
    for(int i=3*X.n+4; i<4*X.n+4; i++) integral = integral + real(X.y[i]) * dx;

    return integral;
}

double integrateImP(AWT & X)
{
    // trapezoidal rule is now applied, later changing to Simpson
    double dx = X.xMax/X.n;
    double integral = 0.5*imag(X.y[X.n]) * dx + 0.5*imag(X.y[4*X.n+3])*dx ;

    for(int i=0; i<X.n; i++) integral = integral + imag(X.y[i]) * dx;
    for(int i=3*X.n+4; i<4*X.n+4; i++) integral = integral + imag(X.y[i]) * dx;

    return integral;
}

// simple integration
double integrateReLeft(AWT & X)
{
    // trapezoidal rule is now applied, later changing to Simpson
    double dx = X.xMax/X.n;
    double integral = 0.5*real(X.y[3*X.n+4]) * dx + 0.5*real(X.y[4*X.n+3])*dx;

    // adding left part
    for(int i=3*X.n+4; i<4*X.n+4; i++) integral = integral + real(X.y[i]) * dx;

    return integral;
}

double integrateReRight(AWT & X)
{
    // trapezoidal rule is now applied, later changing to Simpson
    double dx = X.xMax/X.n;
    double integral = 0.5*real(X.y[0]) * dx + 0.5*real(X.y[X.n])*dx;

    // adding right part
    for(int i=0; i<X.n+1; i++) integral = integral + real(X.y[i]) * dx;

    return integral;
}


double integrateImLeft(AWT & X)
{
    // trapezoidal rule is now applied, later changing to Simpson
    double dx = X.xMax/X.n;
    double integral = 0.5*imag(X.y[3*X.n+4]) * dx + 0.5*imag(X.y[4*X.n+3])*dx;

    // adding left part
    for(int i=3*X.n+4; i<4*X.n+4; i++) integral = integral + imag(X.y[i]) * dx;

    return integral;
}

double integrateImRight(AWT & X)
{
    // trapezoidal rule is now applied, later changing to Simpson
    double dx = X.xMax/X.n;
    double integral = 0.5*imag(X.y[0]) * dx + 0.5*imag(X.y[X.n])*dx;

    // adding right part
    for(int i=0; i<X.n+1; i++) integral = integral + imag(X.y[i]) * dx;

    return integral;
}

void convolution(AWT & out, AWT & X, AWT & Y, aux & help, int s1, int s2)
{

    if(s1*s2==-1)
    {
        help.aux1.conjugateY(X);
        help.aux1.forwardDFT();
        help.aux1.conjugateDFT(help.aux1);
    }
    else
    {
        help.aux1.loadAWTtoAWT(X);
        if(X.dftKnown == 0 )     {  X.forwardDFT();  }        // DFT of X if not available
    }

    if(s1==-1)
    {
        help.aux2.conjugateY(Y);
        help.aux2.forwardDFT();
        help.aux2.conjugateDFT(help.aux2);
    }
    else
    {
        help.aux2.loadAWTtoAWT(Y);
        if(Y.dftKnown == 0 )     {  Y.forwardDFT();  }        // DFT of X if not available
    }

    double dx = X.xMax/X.n;
    for(int i=0; i<X.nn; i++)
    {
        out.yDFT[i] = help.aux1.yDFT[i] * help.aux2.yDFT[i] * dx;
    }

    // the output AWT is always set closed for next convolutionAWT instance
    out.dftKnown = true;
    out.backwardDFT();

    // auxiliary AWTs are always set open for next convolutionAWT instance
    help.aux1.dftKnown = false;
    help.aux2.dftKnown = false;
}
