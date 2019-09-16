#ifndef DFTARRAY_H
#define DFTARRAY_H

#include <complex> // implements the complex class to contain complex numbers in cartesian form
#include <fftw3.h> // subroutine library for computing the discrete Fourier transform (DFT)
#include <iostream>

using namespace std;

class DFTarray
{
    public:
        DFTarray();
        DFTarray(int, double);
        ~DFTarray();


        //////////////////////////      PRIMARY VARIABLES           /////////////////////////////////////////
        int n;                  // valid range is 0..n (positive x) and nn-n..nn-1 (negative x)
        double xMax;            // xMax is the value of the interval <-xMax,xMax>, where the function is stored
        fftw_plan forwardFFT;
        fftw_plan backwardFFT;


        //////////////////////////      DERIVED VARIABLES           /////////////////////////////////////////
        complex<double> * y;     // array of function values
        complex<double> * yDFT;   // DFT of y, DFT means discrete Fourier transform
        fftw_complex * yINTER;    // DFT of y

        int nn;
        int nnDFT;
        int nDFT;                 // number of mesh points of the yDFT array
        double dx;                // grid spacing, x = index * dx;
        bool dftKnown;

        /////////////////////////       FOURIER TRANSFORMS          /////////////////////////////////////////
        void forwardDFT();       // constructor of an empty function forwardDFT
        void backwardDFT();      // constructor of an empty function backwardDFT
        void conjugate();

        /////////////////////////          SET TRANSFORMS          /////////////////////////////////////////
        void setZERO();
        void setIN(int,double);

        void readOUT(int,double);


        /////////////////////////        READ TRANSFORMS          /////////////////////////////////////////
        int readN();
        double readXMAX();
        int readNN();
        double readDX();
        complex<double> readPOSITION(int);

};


#endif // DFTARRAY_H
