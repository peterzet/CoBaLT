#ifndef ARRAYWITHTAILS_H
#define ARRAYWITHTAILS_H

#include "raw.h"

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
#include <stdlib.h>



//// Acronym AWT stands for ArrayWithTails

class AWT: public rawF
{
    public:
        AWT();
        void initializeAWT(int,double,double);     // function to actually construct the object class

        ~AWT();                                    // deconstructor of the class object

        /////////////   STATISTICS  //////////////////////////////////////////////////////
        int nn;                                    // primary variable, division of the mesh points
        // n, xMax, kT are inherited from rawF

        /////////////   MEMBERS   //////////////////////////////////////////////////////
        complex<double> * y;      // array of function values
        complex<double> * yDFT;   // DFT of y, DFT means discrete Fourier transform
        fftw_complex * yINTER;    // DFT of y

        ////////////////  FOURIER TRANSFORMs  ////////////////////////////////////////////
        bool dftKnown;
        fftw_plan forwardFFT;
        fftw_plan backwardFFT;
        void forwardDFT();        // constructor of an empty function forwardDFT
        void backwardDFT();       // constructor of an empty function backwardDFT


        ////////////////  SET AWTs  ///////////////////////////

        void set_zero();
        void set_real(double);
        void set_imag(double);
        void set_FD();
        void set_BE();
        void set_K3();

        ////////////////   IMPORT/EXPORT    ///////////////////////////
        void importAWTasAWT(string &, double, double);
        void importPOKasAWT(string &, double, double, int, int);
        void importFUNasAWT(string &, double, double);


        void exportAWTasAWT(string &, int, string &);
        void exportAWTasFUN(string &, int, double, double, string &);

        void loadAWTtoAWT(AWT &);
        void loadFUNtoAWT(rawF &);

        void asymmetryTest(string &, int, int, int, double, string &);

        void cleanOUT(int);

        ////////////////             OPERATIONS ON AWTs       ///////////////////////////
        ///////////////              which overwrite it!      //////////////////////////

        void AWTForwardDFT();
        void AWTBackDFT();


        void normalFDtimesIM(AWT &, AWT &);
        void inverseFDtimesIM(AWT &, AWT &);
        void normalBEtimesIM(AWT &, AWT &);
        void inverseBEtimesIM(AWT &, AWT &);

        void conjugateY(AWT &);
        void conjugateDFT(AWT &);
        void multiplyAWT(AWT &, complex<double> );

        void deleteReal(AWT &);
        void deleteImag(AWT &);
        void KrammersKronig(AWT &, AWT &, AWT &, AWT &, AWT &, AWT &);
        void KrammersKronigDown(AWT &, AWT &, AWT &, AWT &, AWT &, AWT &);

        ////////////////           OPERATIONS BETWEEN AWTs       ///////////////////////

        void doubleNormal(AWT &, AWT &);
        void doubleInverse(AWT &);



        void convolutionAWT(AWT &, AWT &, int, int, AWT &, AWT &);

        void fermMatsubaraAWT(AWT &, AWT &, int, int, AWT &, AWT &, AWT &, AWT &, AWT &, AWT &, AWT &);
        void fermMatsubaraImP(AWT &, AWT &, int, int, AWT &, AWT &, AWT &, AWT &, AWT &, AWT &, AWT &);
        void fermMatsubaraReP(AWT &, AWT &, int, int, AWT &, AWT &, AWT &, AWT &, AWT &, AWT &, AWT &);

        void boseMatsubaraAWT(AWT &, AWT &, int, int, AWT &, AWT &, AWT &, AWT &, AWT &, AWT &, AWT &, AWT &);
        void boseMatsubaraImP(AWT &, AWT &, int, int, AWT &, AWT &, AWT &, AWT &, AWT &, AWT &, AWT &, AWT &);
        void boseMatsubaraReP(AWT &, AWT &, int, int, AWT &, AWT &, AWT &, AWT &, AWT &, AWT &, AWT &, AWT &);


        complex<double> coth(complex<double>);
        complex<double> tanh(complex<double>);




};

#endif // AWT_H
