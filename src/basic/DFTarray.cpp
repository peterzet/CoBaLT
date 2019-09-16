#include "DFTarray.h"

#include <iostream>         // standard library for reading inputs
#include <fstream>          // standard library for showing outputs
#include <string>           // standard library for manipultaing strings
#include <cassert>          // error handling library, function assert to terminate the program


#include <fftw3.h>          // FFTW library
#include <cmath>            // declares some common mathematical operations and transformation
#include <cstdlib>          // several general purpose functions, including dynamic memory management,
                            // random number generation, communication with the environment,
                            // integer arithmetics, searching, sorting and converting

#include <iomanip>          // awkward, maybe for Set decimal precision?

using std::cout;
using std::endl;
using std::cin;


////////////////    SELECTING FFTW algorithm     /////////////////

#define FFTW_METHOD FFTW_ESTIMATE           // instead of actual measurements of different algorithms,
                                            // a simple heuristic is used to pick a (probably sub-optimal)
                                            // plan quickly.

// #define FFTW_METHOD FFTW_MEASURE
// #define FFTW_METHOD FFTW_PATIENT
// #define FFTW_METHOD FFTW_EXHAUSTIVE


DFTarray::DFTarray()
{
    n = 0;
    xMax = 0;
    dftKnown = false;
    cout << "DFTarray constructed with n = " << n << " ."<< endl;
}


DFTarray::DFTarray(int _n, double _xMax)
{
    n = _n;                 // 2n + 1 is the number of mesh points for the desired function
    xMax = _xMax;           // <-xMax, xMax> is the interval on which the function is defined
    dx = xMax / n;          // dx is the spacing between next stored values of the function
    nn = 4 * n + 4;         // nn is the size of the array used here to store the function
    nDFT = nn;              // I supose nDFT is the size for the Fourier transformed array

    // data is allocated for the input array y by fftw_malloc which is a native function of FFTW
    // (complex<double> *) menas that our input array is complex valued
    // argument of fftw_malloc tells the computer the size of the input array is nn = 4*n +4
    y = (complex<double> *) fftw_malloc(nn * sizeof(complex<double>));
    // this informs you whenever there is no space to store y
    if(y == NULL)
        {
        std::cerr << "insufficient memory!" << std::endl;  // std::cerr
        exit(1);
        }
    // data is allocated for the input array y by fftw_malloc which is a native function of FFTW
    // (complex<double> *) menas that our input array is complex valued
    // argument of fftw_malloc tells the computer the size of the input array is nn = 4*n +4
    yINTER = (fftw_complex *) fftw_malloc(nDFT * sizeof(fftw_complex));
    // this informs you whenever there is no space to store yDFT
    if(yINTER == NULL)
        {
        std::cerr << "insufficient memory!" << std::endl;
        exit(1);
        }
    // ASK!!!!, it is probably required because of variable type mismatch
    yDFT = (complex<double> *) yINTER;
    // A plan forwardFFT is created for forward fast Fourier transform
    forwardFFT = fftw_plan_dft_1d(nn,
                                  (fftw_complex *) y,
                                  (fftw_complex *) yDFT,
                                  FFTW_FORWARD,                         // FFTW_FORWARD = -1, it is integer
                                  FFTW_METHOD | FFTW_PRESERVE_INPUT);
    // A plan fbackwardFFT is created forbackward fast Fourier transform
    backwardFFT = fftw_plan_dft_1d(nn,
                                  (fftw_complex *) yDFT,
                                  (fftw_complex *) y,
                                  FFTW_BACKWARD,
                                  FFTW_METHOD | FFTW_PRESERVE_INPUT);

    // the array is ZERO PADDED, initial filling of the data
    for(int i=0; i<nn; i++)
        {
        y[i] = 0.0;
        yDFT[i] = 0.0;
        }
    dftKnown = false;
    cout << "DFTarray constructed with n = " << n << "." << endl;
}


DFTarray::~DFTarray()
{
    fftw_free(y);
    fftw_free(yINTER);
    //fftw_destroy_plan(forwardFFT);
    //fftw_destroy_plan(backwardFFT);
    cout << "DFTarray deconstructed." << endl;
}

/////////////              SET FUNCTIONS             //////////////////////

void DFTarray::setZERO()
{
    // the array is ZERO PADDED, initial filling of the data
    for(int i=0; i<nn; i++)
        {
        y[i] = 0.0;
        yDFT[i] = 0.0;
        }
}

/////////////            RETURN FUNCTIONS             //////////////////////

int DFTarray::readN()
{
    return n;
}

double DFTarray::readXMAX()
{
    return xMax;
}


int DFTarray::readNN()
{
    return nn;
}

double DFTarray::readDX()
{
    return dx;
}

complex<double> DFTarray::readPOSITION(int iii)
{
    // the array is ZERO PADDED, initial filling of the data
    return y[iii];
}


//////////////////        the PLANS are EXECUTED  via appropriate functions     ////////////

/////////////////      forward fast Fourier transform        //////////////////////////////
void DFTarray::forwardDFT()
{
    if(dftKnown)
    return;

    fftw_execute(forwardFFT);

    dftKnown = true;
}

////////////////      backward fast Fourier transform        //////////////////////////
void DFTarray::backwardDFT()
{
    assert(dftKnown);

    fftw_execute(backwardFFT);

    // normalize with the size of the transformed array and erase garbage at the place of the zero padding
    for(int i=0; i<=n; i++)
        y[i] /= nn;
    for(int i=n+1; i<nn-n; ++i)
        y[i] = 0;
    for(int i=nn-n; i<nn; i++)
        y[i] /= nn;
}





