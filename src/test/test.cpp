#include "../physics.h"

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
#include <boost/lexical_cast.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/filesystem.hpp>
//#include "gnuplot-iostream.h"

using namespace std;

#define KNOWN_TEST cout << aux1.dftKnown << " " << aux2.dftKnown << " " << aux3.dftKnown << " " << aux4.dftKnown << " " << aux5.dftKnown << " " << aux6.dftKnown << endl;

/////////////////////////////////////////////////////////////////////////////////////////////////////

void test_convolutions(int n, double xMax, double range, string export_mode)
{
    cout << "Test of implementation of the convolution theorem in the program" << endl;
    cout << endl;
    cout << "Output is written as separate txt files in directory /Test" << endl;

    string name;

        // these four AWTs are used as help and storage place, are once allocated in the stack
    AWT aux1;
    aux1.initializeAWT(n, xMax, 0);
    AWT aux2;
    aux2.initializeAWT(n, xMax, 0);
    AWT aux3;
    aux3.initializeAWT(n, xMax, 0);
    AWT aux4;
    aux4.initializeAWT(n, xMax, 0);
    AWT aux5;
    aux5.initializeAWT(n, xMax, 0);
    AWT aux6;
    aux6.initializeAWT(n, xMax, 0);

    AWT test;
    test.initializeAWT(n, xMax, 0);

    for(int i=0;       i < n + 1;    i++)     test.y[i] = ( (i*xMax/n)  )*exp( -1.0*(i*xMax/n)*(i*xMax/n) )  ;
    for(int i=n+1;     i < 3*n + 4;  i++)     test.y[i] = 0;
    for(int i=3*n+4;   i < 4*n + 4;  i++)     test.y[i] = ( -(xMax/n)*(4*n + 4 - i)  )*exp( -1.0*( (4*n + 4 - i)*xMax/n)*( (4*n + 4 -i)*xMax/n) )  ;

    AWT test1;
    test1.initializeAWT(n, xMax, 0);

    for(int i=0;       i < n + 1;    i++)     test1.y[i] = ( (i*xMax/n) - 1.0 )*exp( -1.0*(i*xMax/n)*(i*xMax/n) )  ;
    for(int i=n+1;     i < 3*n + 4;  i++)     test1.y[i] = 0;
    for(int i=3*n+4;   i < 4*n + 4;  i++)     test1.y[i] = ( -(xMax/n)*(4*n + 4 - i) - 1.0 )*exp( -1.0*( (4*n + 4 - i)*xMax/n)*( (4*n + 4 -i)*xMax/n) )  ;


    name ="Test/testre";
    test.exportAWTasFUN(name, 10, range, 1.0, export_mode);

    name ="Test/testre1";
    test1.exportAWTasFUN(name, 10, range, 1.0, export_mode);

    AWT real1;
    real1.initializeAWT(n, xMax, 0);
    real1.convolutionAWT(test, test1, 1, 1, aux1, aux2);
    name ="Test/real1";
    real1.exportAWTasFUN(name, 10, range, 1.0, export_mode);


    AWT real2;
    real2.initializeAWT(n, xMax, 0);
    real2.convolutionAWT(test, test1, 1, -1, aux1, aux2);
    name ="Test/real2";
    real2.exportAWTasFUN(name, 10, range, 1.0, export_mode);

    AWT real3;
    real3.initializeAWT(n, xMax, 0);
    real3.convolutionAWT(test, test1, -1, 1, aux1, aux2);
    name ="Test/real3";
    real3.exportAWTasFUN(name, 10, range, 1.0, export_mode);

    AWT real4;
    real4.initializeAWT(n, xMax, 0);
    real4.convolutionAWT(test, test1, -1, -1, aux1, aux2);
    name ="Test/real4";
    real4.exportAWTasFUN(name, 10, range, 1.0, export_mode);

    cout << endl;
    cout << "Purely real valued functions have been convolved" << endl;
    cout << "Convolved functions are stored as testre and testre1" << endl;
    cout << "Results according to the convolution theorem are stored in real[1-4]" << endl;
    cout << "Results according to Newton integrals are stored in math_real[1-4]" << endl;
    cout << endl;

    AWT math_real1;
    math_real1.initializeAWT(n, xMax, 0);

    for(int i=0;       i < n + 1;    i++)
    {
        double xx = (i*xMax/n);
        math_real1.y[i] = 0.31332534*(-1.0 -2.0*xx +xx*xx )*exp(-0.5*xx*xx);
    }
    for(int i=n+1;     i < 3*n + 4;  i++)
    {
        math_real1.y[i] = 0;
    }
    for(int i=3*n+4;   i < 4*n + 4;  i++)
    {
        double xx =-(xMax/n)*(4*n + 4 - i);
        math_real1.y[i] = 0.31332534*(-1.0 -2.0*xx +xx*xx )*exp(-0.5*xx*xx);
    }

    name ="Test/math_real1";
    math_real1.exportAWTasFUN(name, 10, range, 1.0, export_mode);

    AWT math_real2;
	math_real2.initializeAWT(n, xMax, 0);

    for(int i=0;       i < n + 1;    i++)
    {
        double xx = (i*xMax/n);
        math_real2.y[i] = 0.31332534*(1.0 +2.0*xx -xx*xx )*exp(-0.5*xx*xx);
    }
    for(int i=n+1;     i < 3*n + 4;  i++)
    {
        math_real2.y[i] = 0;
    }
    for(int i=3*n+4;   i < 4*n + 4;  i++)
    {
        double xx =-(xMax/n)*(4*n + 4 - i);
        math_real2.y[i] = 0.31332534*(1.0 +2.0*xx -xx*xx )*exp(-0.5*xx*xx);
    }

    name ="Test/math_real2";
    math_real2.exportAWTasFUN(name, 10, range, 1.0, export_mode);

    AWT math_real3;
    math_real3.initializeAWT(n, xMax, 0);

    for(int i=0;       i < n + 1;    i++)
    {
        double xx = (i*xMax/n);
        math_real3.y[i] = 0.31332534*(-1.0 +2.0*xx +xx*xx )*exp(-0.5*xx*xx);
    }
    for(int i=n+1;     i < 3*n + 4;  i++)
    {
        math_real3.y[i] = 0;
    }
    for(int i=3*n+4;   i < 4*n + 4;  i++)
    {
        double xx =-(xMax/n)*(4*n + 4 - i);
        math_real3.y[i] = 0.31332534*(-1.0 +2.0*xx +xx*xx )*exp(-0.5*xx*xx);
    }

    name ="Test/math_real3";
    math_real3.exportAWTasFUN(name, 10, range, 1.0, export_mode);

    AWT math_real4;
    math_real4.initializeAWT(n, xMax, 0);

    for(int i=0;       i < n + 1;    i++)
    {
        double xx = (i*xMax/n);
        math_real4.y[i] = 0.31332534*(1.0 -2.0*xx -xx*xx )*exp(-0.5*xx*xx);
    }
    for(int i=n+1;     i < 3*n + 4;  i++)
    {
        math_real4.y[i] = 0;
    }
    for(int i=3*n+4;   i < 4*n + 4;  i++)
    {
        double xx =-(xMax/n)*(4*n + 4 - i);
        math_real4.y[i] = 0.31332534*(1.0 -2.0*xx -xx*xx )*exp(-0.5*xx*xx);
    }

    name ="Test/math_real4";
    math_real4.exportAWTasFUN(name, 10, range, 1.0, export_mode);

    AWT testim;
    testim.initializeAWT(n, xMax, 0);

    for(int i=0;       i < n + 1;    i++)
    {
        double a = ( (i*xMax/n)  )*exp( -1.0*(i*xMax/n)*(i*xMax/n) )  ;
        complex<double> u(0,a);
        testim.y[i] = u;
    }
    for(int i=n+1;     i < 3*n + 4;  i++)     testim.y[i] = 0;
    for(int i=3*n+4;   i < 4*n + 4;  i++)
    {
        double a = ( -(xMax/n)*(4*n + 4 - i)  )*exp( -1.0*( (4*n + 4 - i)*xMax/n)*( (4*n + 4 -i)*xMax/n) )  ;
        complex<double> u(0,a);
        testim.y[i] = u;
    }

    AWT testim1;
    testim1.initializeAWT(n, xMax, 0);

    for(int i=0;       i < n + 1;    i++)
    {
        double a = ( (i*xMax/n) - 1.0 )*exp( -1.0*(i*xMax/n)*(i*xMax/n) )  ;
        complex<double> u(0,a);
        testim1.y[i] = u;
    }
    for(int i=n+1;     i < 3*n + 4;  i++)     testim1.y[i] = 0;
    for(int i=3*n+4;   i < 4*n + 4;  i++)
    {
        double a = ( -(xMax/n)*(4*n + 4 - i) - 1.0 )*exp( -1.0*( (4*n + 4 - i)*xMax/n)*( (4*n + 4 -i)*xMax/n) )  ;
        complex<double> u(0,a);
        testim1.y[i] = u;
    }


    name ="Test/testim";
    testim.exportAWTasFUN(name, 10, range, 1.0, export_mode);

    name ="Test/testim1";
    testim1.exportAWTasFUN(name, 10, range, 1.0, export_mode);

    AWT imag1;
    imag1.initializeAWT(n, xMax, 0);
    imag1.convolutionAWT(testim, testim1, 1, 1, aux1, aux2);
    name ="Test/imag1";
    imag1.exportAWTasFUN(name, 10, range, 1.0, export_mode);


    AWT imag2;
    imag2.initializeAWT(n, xMax, 0);
    imag2.convolutionAWT(testim, testim1, 1, -1, aux1, aux2);
    name ="Test/imag2";
    imag2.exportAWTasFUN(name, 10, range, 1.0, export_mode);

    AWT imag3;
    imag3.initializeAWT(n, xMax, 0);
    imag3.convolutionAWT(testim, testim1, -1, 1, aux1, aux2);
    name ="Test/imag3";
    imag3.exportAWTasFUN(name, 10, range, 1.0, export_mode);

    AWT imag4;
    imag4.initializeAWT(n, xMax, 0);
    imag4.convolutionAWT(testim, testim1, -1, -1, aux1, aux2);
    name ="Test/imag4";
    imag4.exportAWTasFUN(name, 10, range, 1.0, export_mode);

    cout << endl;
    cout << "Purely imaginary valued functions have been convolved" << endl;
    cout << "Convolved functions are stored as testim and testim1" << endl;
    cout << "Results according to the convolution theorem are stored in imag[1-4]" << endl;
    cout << "Results according to Newton integrals are stored in math_imag[1-4]" << endl;
    cout << endl;


    AWT math_imag1;
    math_imag1.initializeAWT(n, xMax, 0);

    for(int i=0;       i < n + 1;    i++)
    {
        double xx = (i*xMax/n);
        double a = -0.31332534*(-1.0 -2.0*xx +xx*xx )*exp(-0.5*xx*xx);
        complex<double> u(a,0);
        math_imag1.y[i] = u;
    }
    for(int i=n+1;     i < 3*n + 4;  i++)
    {
        math_imag1.y[i] = 0;
    }
    for(int i=3*n+4;   i < 4*n + 4;  i++)
    {
        double xx =-(xMax/n)*(4*n + 4 - i);
        double a = -0.31332534*(-1.0 -2.0*xx +xx*xx )*exp(-0.5*xx*xx);
        complex<double> u(a,0);
        math_imag1.y[i] = u;
    }

    name ="Test/math_imag1";
    math_imag1.exportAWTasFUN(name, 10, range, 1.0, export_mode);

    AWT math_imag2;
	math_imag2.initializeAWT(n, xMax, 0);

    for(int i=0;       i < n + 1;    i++)
    {
        double xx = (i*xMax/n);
        double a = -0.31332534*(1.0 +2.0*xx -xx*xx )*exp(-0.5*xx*xx);
        complex<double> u(a,0);
        math_imag2.y[i] = u;
    }
    for(int i=n+1;     i < 3*n + 4;  i++)
    {
        math_imag2.y[i] = 0;
    }
    for(int i=3*n+4;   i < 4*n + 4;  i++)
    {
        double xx =-(xMax/n)*(4*n + 4 - i);
        double a = -0.31332534*(1.0 +2.0*xx -xx*xx )*exp(-0.5*xx*xx);
        complex<double> u(a,0);
        math_imag2.y[i] = u;
    }

    name ="Test/math_imag2";
    math_imag2.exportAWTasFUN(name, 10, range, 1.0, export_mode);

    AWT math_imag3;
    math_imag3.initializeAWT(n, xMax, 0);

    for(int i=0;       i < n + 1;    i++)
    {
        double xx = (i*xMax/n);
        double a = -0.31332534*(-1.0 +2.0*xx +xx*xx )*exp(-0.5*xx*xx);
        complex<double> u(a,0);
        math_imag3.y[i] = u;
    }
    for(int i=n+1;     i < 3*n + 4;  i++)
    {
        math_imag3.y[i] = 0;
    }
    for(int i=3*n+4;   i < 4*n + 4;  i++)
    {
        double xx =-(xMax/n)*(4*n + 4 - i);
        double a = -0.31332534*(-1.0 +2.0*xx +xx*xx )*exp(-0.5*xx*xx);
        complex<double> u(a,0);
        math_imag3.y[i] = u;
    }

    name ="Test/math_imag3";
    math_imag3.exportAWTasFUN(name, 10, range, 1.0, export_mode);

    AWT math_imag4;
    math_imag4.initializeAWT(n, xMax, 0);

    for(int i=0;       i < n + 1;    i++)
    {
        double xx = (i*xMax/n);
        double a = -0.31332534*(1.0 -2.0*xx -xx*xx )*exp(-0.5*xx*xx);
        complex<double> u(a,0);
        math_imag4.y[i] = u;
    }
    for(int i=n+1;     i < 3*n + 4;  i++)
    {
        math_imag4.y[i] = 0;
    }
    for(int i=3*n+4;   i < 4*n + 4;  i++)
    {
        double xx =-(xMax/n)*(4*n + 4 - i);
        double a = -0.31332534*(1.0 -2.0*xx -xx*xx )*exp(-0.5*xx*xx);
        complex<double> u(a,0);
        math_imag4.y[i] = u;
    }

    name ="Test/math_imag4";
    math_imag4.exportAWTasFUN(name, 10, range, 1.0, export_mode);

    AWT testco;
    testco.initializeAWT(n, xMax, 0);

    for(int i=0;       i < n + 1;    i++)
    {
        double xx = (i*xMax/n);
        double a = xx*exp(-xx*xx);
        double b = (xx+1.0)*exp(-xx*xx);
        complex<double> u(a,b);
        testco.y[i] = u;
    }
    for(int i=n+1;     i < 3*n + 4;  i++)     testco.y[i] = 0;
    for(int i=3*n+4;   i < 4*n + 4;  i++)
    {
        double xx = -(xMax/n)*(4*n + 4 - i);
        double a = xx*exp(-xx*xx);
        double b = (xx+1.0)*exp(-xx*xx);
        complex<double> u(a,b);
        testco.y[i] = u;
    }

    AWT testco1;
    testco1.initializeAWT(n, xMax, 0);

    for(int i=0;       i < n + 1;    i++)
    {
        double xx = (i*xMax/n);
        double a = (xx-1.0)*exp(-xx*xx);
        double b = (xx-2.0)*exp(-xx*xx);
        complex<double> u(a,b);
        testco1.y[i] = u;
    }
    for(int i=n+1;     i < 3*n + 4;  i++)     testco1.y[i] = 0;
    for(int i=3*n+4;   i < 4*n + 4;  i++)
    {
        double xx = -(xMax/n)*(4*n + 4 - i);
        double a = (xx-1.0)*exp(-xx*xx);
        double b = (xx-2.0)*exp(-xx*xx);
        complex<double> u(a,b);
        testco1.y[i] = u;
    }


    name ="Test/testco";
    testco.exportAWTasFUN(name, 10, range, 1.0, export_mode);

    name ="Test/testco1";
    testco1.exportAWTasFUN(name, 10, range, 1.0, export_mode);

    AWT comp1;
    comp1.initializeAWT(n, xMax, 0);
    comp1.convolutionAWT(testco, testco1, 1, 1, aux1, aux2);
    name ="Test/comp1";
    comp1.exportAWTasFUN(name, 10, range, 1.0, export_mode);


    AWT comp2;
    comp2.initializeAWT(n, xMax, 0);
    comp2.convolutionAWT(testco, testco1, 1, -1, aux1, aux2);
    name ="Test/comp2";
    comp2.exportAWTasFUN(name, 10, range, 1.0, export_mode);

    AWT comp3;
    comp3.initializeAWT(n, xMax, 0);
    comp3.convolutionAWT(testco, testco1, -1, 1, aux1, aux2);
    name ="Test/comp3";
    comp3.exportAWTasFUN(name, 10, range, 1.0, export_mode);

    AWT comp4;
    comp4.initializeAWT(n, xMax, 0);
    comp4.convolutionAWT(testco, testco1, -1, -1, aux1, aux2);
    name ="Test/comp4";
    comp4.exportAWTasFUN(name, 10, range, 1.0, export_mode);

    cout << endl;
    cout << "Complex valued functions have been convolved" << endl;
    cout << "Convolved functions are stored as testco and testco1" << endl;
    cout << "Results according to the convolution theorem are stored in comp[1-4]" << endl;
    cout << "Results according to Newton integrals are stored in math_comp[1-4]" << endl;
    cout << endl;



    AWT math_comp1;
    math_comp1.initializeAWT(n, xMax, 0);

    for(int i=0;       i < n + 1;    i++)
    {
        double xx = (i*xMax/n);
        double a = 1.253314137*(2.0 )*exp(-0.5*xx*xx);
        double b = 1.253314137*(-1.5 -xx +0.5*xx*xx )*exp(-0.5*xx*xx);
        complex<double> u(a,b);
        math_comp1.y[i] = u;
    }
    for(int i=n+1;     i < 3*n + 4;  i++)
    {
        math_comp1.y[i] = 0;
    }
    for(int i=3*n+4;   i < 4*n + 4;  i++)
    {
        double xx =-(xMax/n)*(4*n + 4 - i);
        double a = 1.253314137*(2.0 )*exp(-0.5*xx*xx);
        double b = 1.253314137*(-1.5 -xx +0.5*xx*xx )*exp(-0.5*xx*xx);
        complex<double> u(a,b);
        math_comp1.y[i] = u;
    }

    name ="Test/math_comp1";
    math_comp1.exportAWTasFUN(name, 10, range, 1.0, export_mode);

    AWT math_comp2;
	math_comp2.initializeAWT(n, xMax, 0);

    for(int i=0;       i < n + 1;    i++)
    {
        double xx = (i*xMax/n);
        double a = 1.253314137*(2.0 -xx )*exp(-0.5*xx*xx);
        double b = 1.253314137*(-0.5 + 2.0*xx -0.5*xx*xx )*exp(-0.5*xx*xx);
        complex<double> u(a,b);
        math_comp2.y[i] = u;
    }
    for(int i=n+1;     i < 3*n + 4;  i++)
    {
        math_comp2.y[i] = 0;
    }
    for(int i=3*n+4;   i < 4*n + 4;  i++)
    {
        double xx =-(xMax/n)*(4*n + 4 - i);
        double a = 1.253314137*(2.0 -xx )*exp(-0.5*xx*xx);
        double b = 1.253314137*(-0.5 + 2.0*xx -0.5*xx*xx )*exp(-0.5*xx*xx);
        complex<double> u(a,b);
        math_comp2.y[i] = u;
    }

    name ="Test/math_comp2";
    math_comp2.exportAWTasFUN(name, 10, range, 1.0, export_mode);

    AWT math_comp3;
    math_comp3.initializeAWT(n, xMax, 0);

    for(int i=0;       i < n + 1;    i++)
    {
        double xx = (i*xMax/n);
        double a = 1.253314137*(2.0 )*exp(-0.5*xx*xx);
        double b = 1.253314137*(-1.5 + 1.0*xx +0.5*xx*xx )*exp(-0.5*xx*xx);
        complex<double> u(a,b);
        math_comp3.y[i] = u;
    }
    for(int i=n+1;     i < 3*n + 4;  i++)
    {
        math_comp3.y[i] = 0;
    }
    for(int i=3*n+4;   i < 4*n + 4;  i++)
    {
        double xx =-(xMax/n)*(4*n + 4 - i);
        double a = 1.253314137*(2.0 )*exp(-0.5*xx*xx);
        double b = 1.253314137*(-1.5 + 1.0*xx +0.5*xx*xx )*exp(-0.5*xx*xx);
        complex<double> u(a,b);
        math_comp3.y[i] = u;
    }

    name ="Test/math_comp3";
    math_comp3.exportAWTasFUN(name, 10, range, 1.0, export_mode);

    AWT math_comp4;
    math_comp4.initializeAWT(n, xMax, 0);

    for(int i=0;       i < n + 1;    i++)
    {
        double xx = (i*xMax/n);
        double a = 1.253314137*(2.0 +xx)*exp(-0.5*xx*xx);
        double b = 1.253314137*(-0.5 - 2.0*xx -0.5*xx*xx )*exp(-0.5*xx*xx);
        complex<double> u(a,b);
        math_comp4.y[i] = u;
    }
    for(int i=n+1;     i < 3*n + 4;  i++)
    {
        math_comp4.y[i] = 0;
    }
    for(int i=3*n+4;   i < 4*n + 4;  i++)
    {
        double xx =-(xMax/n)*(4*n + 4 - i);
        double a = 1.253314137*(2.0 +xx)*exp(-0.5*xx*xx);
        double b = 1.253314137*(-0.5 - 2.0*xx -0.5*xx*xx )*exp(-0.5*xx*xx);
        complex<double> u(a,b);
        math_comp4.y[i] = u;
    }

    name ="Test/math_comp4";
    math_comp4.exportAWTasFUN(name, 10, range, 1.0, export_mode);

	cout << "convolution test finished" << endl;
}
