#include "AWT.h"

#include <vector>
#include <iostream>         // standard library for reading inputs
#include <fstream>          // standard library for showing outputs
#include <string>           // standard library for manipultaing strings
#include <cassert>          // error handling library, function assert to terminate the program

#include <complex>          // implements the complex class to contain complex numbers in cartesian form
#include <fftw3.h>          // FFTW library
#include <cmath>            // declares some common mathematical operations and transformation
#include <cstdlib>          // several general purpose functions, including dynamic memory management,
                            // random number generation, communication with the environment,
                            // integer arithmetics, searching, sorting and converting
#include <stdlib.h>


////////////////    SELECTING FFTW algorithm     /////////////////

#define CX11
// #define OLD

#define FFTW_METHOD FFTW_ESTIMATE           // instead of actual measurements of different algorithms,
                                            // a simple heuristic is used to pick a (probably sub-optimal)
                                            // plan quickly.

// #define FFTW_METHOD FFTW_MEASURE
// #define FFTW_METHOD FFTW_PATIENT
// #define FFTW_METHOD FFTW_EXHAUSTIVE

// Formal constructor of AWT
AWT::AWT(){   initializeAWT(0, 0, 0);  }

// The initializer of AWT
void AWT::initializeAWT(int _n, double _x, double _kT)
{
    n    = _n;              // 2_n + 1 is the number of the mesh points
    xMax = _x;              // <-xMax, xMax> is the interval on which the function is defined
    nn   = 4 * n  + 4;      // nnHelp is the size of the array to store the function
    kT = _kT;

    // data is allocated for y in the heap, the argument of fftw_malloc defines the size of the input array as nnHelp = 4*n +4
    y = (complex<double> *) fftw_malloc(nn * sizeof(complex<double>));

    //cout << "taking space " << nn * sizeof(complex<double>) << endl;
    // this informs you whenever there is no space to store y
    if(y == NULL)
    {
        std::cerr << "insufficient memory!" << std::endl;  // std::cerr
        exit(1);
    }

    // ASK!!!!, it is probably required because of variable type mismatch
    yDFT = (complex<double> *) fftw_malloc(nn * sizeof(fftw_complex));
    if(yDFT == NULL)
    {
        std::cerr << "insufficient memory!" << std::endl;
        exit(1);
    }
    // A plan forwardFFT is created for forward fast Fourier transform
    forwardFFT  = fftw_plan_dft_1d(nn, (fftw_complex *) y, (fftw_complex *) yDFT, FFTW_FORWARD, FFTW_METHOD | FFTW_PRESERVE_INPUT);
    // A plan fbackwardFFT is created forbackward fast Fourier transform
    backwardFFT = fftw_plan_dft_1d(nn, (fftw_complex *) yDFT, (fftw_complex *) y, FFTW_BACKWARD, FFTW_METHOD | FFTW_PRESERVE_INPUT);

}


// The destructor of AWT
AWT::~AWT()
{
    fftw_free(y);
    fftw_free(yDFT);
    fftw_destroy_plan(forwardFFT);
    fftw_destroy_plan(backwardFFT);
}

//////////////////////////////////////////////////////////////////////////
/////////////          FOURIER TRANSFORMs          //////////////////////
////////////////////////////////////////////////////////////////////////

// Forward and Backward discrete Fourier transform
void AWT::forwardDFT()
{
    fftw_execute(forwardFFT);
}

void AWT::backwardDFT()
{
    fftw_execute(backwardFFT);

    // normalize with the size of the transformed array and erase garbage at the place of the zero padding
    for(int i=0;     i<=n; i++)      y[i] /= nn;
    for(int i=n+1; i<nn-n; i++)      y[i] = 0;
    for(int i=nn-n;  i<nn; i++)      y[i] /= nn;
}

//////////////////////////////////////////////////////////////////////
/////////////              SET AWTs            //////////////////////
////////////////////////////////////////////////////////////////////

// the array is filled with complex zero
void AWT::setZero()
{
    for(int i=0; i<nn; i++)
    {
        complex<double> u(0,0);
        y[i] = u;
    }
}

void AWT::set_real(double numb)
{
    for(int i=0; i<nn; i++)
    {
        complex<double> u(numb,0);
        y[i] = u;
    }
}

// The initializer of FERMI-DIRAC and BOSE-EINSTEIN distribution
void AWT::setFD(int _n, double _x, double _kT)
{
    if( _kT == 0 )
    {
                                                 y[0] = 0.5;
        for(int i=1;      i < 3*_n + 4; i++)     y[i] = 0;
        for(int i=3*_n+4; i < 4*_n + 4; i++)     y[i] = 1.0;
    }
    else
    {
                                                 y[0] = 1.0/( exp( +1e-12) + 1 );
        for(int i=1;      i < _n + 1; i++)       y[i] = 1.0/( exp(i*_x/ (_n * _kT)  ) + 1 );
        for(int i=_n + 1; i < 3*_n + 4; i++)     y[i] = 0;
        for(int i=3*_n+4; i < 4*_n + 4; i++)     y[i] = 1.0/( exp(  (i - 4*_n -4)*_x / (_n *_kT)  ) + 1 );
    }
}

// The initializer of THE FIRST DERIVATIVE OF FERMI-DIRAC function
void AWT::setFDder(int _n, double _x, double _kT)
{
    // dirac delta function
    if( _kT == 0 )
    {
                                                 y[0] = 1;
        for(int i=1;      i < 3*_n + 4; i++)     y[i] = 0;
        for(int i=3*_n+4; i < 4*_n + 4; i++)     y[i] = 0;
    }
    else
    {
                                                 y[0] = -1.0*exp( +1e-12)/( exp( +1e-12) + 1 )/( exp( +1e-12) + 1 );
        for(int i=1;      i < _n + 1; i++)       y[i] = -1.0*exp(i*_x/ (_n * _kT) )/( exp(i*_x/ (_n * _kT)  ) + 1 )/( exp(i*_x/ (_n * _kT)  ) + 1 )/_kT;
        for(int i=_n + 1; i < 3*_n + 4; i++)     y[i] = 0;
        for(int i=3*_n+4; i < 4*_n + 4; i++)     y[i] = -1.0*exp(  (i - 4*_n -4)*_x / (_n *_kT)  ) /( exp(  (i - 4*_n -4)*_x / (_n *_kT)  ) + 1 )/( exp(  (i - 4*_n -4)*_x / (_n *_kT)  ) + 1 )/_kT;
    }
}

void AWT::setBE(int _n, double _x, double _kT)
{
    if( _kT == 0 )
    {
                                                 y[0] = -0.5;
        for(int i=1;      i < 3*_n + 4; i++)     y[i] = 0;
        for(int i=3*_n+4; i < 4*_n + 4; i++)     y[i] = -1.0;
    }
    else
    {
                                                 y[0] = -0.5; // 1.0/( exp( +1e-12/kT ) - 1 );
        for(int i=1;      i <   _n + 1; i++)     y[i] = 1.0/( exp(i*_x/ (_n * kT)  ) - 1 );
        for(int i=_n + 1; i < 3*_n + 4; i++)     y[i] = 0;
        for(int i=3*_n+4; i < 4*_n + 4; i++)     y[i] = 1.0/( exp( (i - 4*_n -4)*_x / (_n *_kT) ) - 1 );
    }
}




////////////////////////////////////////////////////////////////////////
/////////                  Read AWT in                       //////////
//////////////////////////////////////////////////////////////////////
void AWT::importAWTasAWT(string & name, double x, double kT)
{

    ifstream ListFile;
    ListFile.open (name.c_str());   // list contains all the file names to be processed
    string item;                    // used for temporary store the data
    int counter = 0;                // counter is initialized
    double re;                      // temporary storage for real part of complex<double>
    double im;                      // temporary storage for imaginery part of complex<double>

    // 1. the size of the stored AWT is determined
    if (ListFile.fail()) {cerr<< "Error opening the file" << endl; exit(1);}

    while (!ListFile.eof())  // if not end of file
    {
        ListFile >> item;
        counter = counter+1;
    }
    ListFile.close();

    // 2. AWT with the right size is initialized
    initializeAWT((counter-21)/20, x, kT);

    // 3. the initialized AWT is filled with data
    counter =0;
    ListFile.open (name.c_str());
    if (ListFile.fail()) {cerr<< "Error opening the file" << endl; exit(1);}
    while (!ListFile.eof())  // if not end of file
    {
        ListFile >> item;
        if(counter % 5 == 1)    re = atof(item.c_str());      // scans for entries in the second column in the file


        if(counter % 5 == 2)                // scans for entries in the third column in the file
        {
            im = atof(item.c_str());
            complex<double> u(re, im);
            y[(int) counter/5] = u;
        }
        if(counter % 5 == 3)    re = atof(item.c_str());     // scans for entries in the 4th column in the file

        if(counter % 5 == 4)                // scans for entries in the 5th column in the file
        {
            im = atof(item.c_str());
            complex<double> u(re, im);
            yDFT[(int) counter/5] = u;
        }
        counter = counter+1;
    }

    ListFile.close();
}


////////////////////////////////////////////////////////////////////////
/////////                  import POK in                       //////////
//////////////////////////////////////////////////////////////////////
void AWT::importPOKasAWT(string & name, double x, double kT, int columns, int position)
{

    ifstream ListFile;
    ListFile.open (name.c_str());   // list contains all the file names to be processed
    string item;                    // used for temporary store the data
    int counter = 0;                // counter is initialized
    double re;                      // temporary storage for real part of complex<double>
    double im;                      // temporary storage for imaginery part of complex<double>

    // 1. the size of the stored AWT is determined
    if (ListFile.fail()) {cerr<< "Error opening the file" << endl; exit(1);}

    while (!ListFile.eof())  // if not end of file
    {
        ListFile >> item;
        counter = counter+1;
    }
    ListFile.close();

    // 2. AWT with the right size is initialized
    initializeAWT((counter-1-4*columns)/(4*columns), x, kT);

    // 3. the initialized AWT is filled with data
    counter =0;
    ListFile.open (name.c_str());
    if (ListFile.fail()) {cerr<< "Error opening the file" << endl; exit(1);}
    while (!ListFile.eof())  // if not end of file
    {
        ListFile >> item;
        if(counter % columns == position)    re = atof(item.c_str());      // scans for entries in the second column in the file


        if(counter % columns == position + 1)                // scans for entries in the third column in the file
        {
            im = atof(item.c_str());
            complex<double> u(re, im);
            y[(int) counter/5] = u;
        }
        if(counter % 5 == 3)    re = atof(item.c_str());     // scans for entries in the 4th column in the file

        if(counter % 5 == 4)                // scans for entries in the 5th column in the file
        {
            im = atof(item.c_str());
            complex<double> u(re, im);
            yDFT[(int) counter/5] = u;
        }
        counter = counter+1;
    }

    ListFile.close();
}

////////////////////////////////////////////////////////////////////////
/////////             Read function in                       //////////
//////////////////////////////////////////////////////////////////////
void AWT::importFUNasAWT(string & name, double x, double kT)
{
    ifstream ListFile;
    ListFile.open (name.c_str());   // list contains all the file names to be processed
    string item;                    // used for temporary store the data
    int counter = 0;                // counter is initialized
    double re;                      // temporary storage for real part of complex<double>
    double im;                      // temporary storage for imaginery part of complex<double>

    int j=0;

    // 1. the size of the needed AWT is determined from the number of lines
    if (ListFile.fail()) {cerr<< "Error opening the file" << endl; exit(1);}
    while (!ListFile.eof())  // if not end of file
    {
        ListFile >> item;
        counter = counter+1;
    }
    int n = (counter-4)/6+1;
    ListFile.close();

    // 2. AWT with the right size is initialized
    initializeAWT(n, x, kT);

    // 3. the initialized AWT is filled with data
    counter =0;     // reset the counter
    ListFile.open (name.c_str());
    if (ListFile.fail()) {cerr<< "Error opening the file" << endl; exit(1);}

    while (counter < 6*n+3 )  // while not end of file
    {
        ListFile >> item;
        if(counter % 3 == 1)    re = atof(item.c_str());       // scans for entries in the second column in the file

        if(counter % 3 == 2)                                   // scans for entries in the third column in the file
        {
            im = atof(item.c_str());
            if(counter < 3*n)
            {
                j = (counter-2)/3;
                complex<double> u(re, im);
                y[3*n + 4 + j] = u;
            }
            else
            {
                j = (counter - 3*n - 2)/3;
                complex<double> u(re, im);
                y[j] = u;
            }
        }

        counter = counter+1;
    }
    // zero padding of the centre of the array
    for(j=n+1; j < 3*n +4; j++)
    {
        complex<double> u(0, 0);
        y[j] = u;
    }
    // initializing the yDFT part of the body to zero
    for(j=0; j < 4*n +4; j++)
    {
        complex<double> u(0, 0);
        yDFT[j] = u;
    }

    ListFile.close();
}

void AWT::loadAWTtoAWT(AWT & in)
{
    if(n != in.n)
    {
        cerr<< "Error: AWT and AWT are not compatible" << endl; exit(1);
    }
    int i;
    for(i=0; i < 4*n+4; i++)
    {
        y[i] = in.y[i];
        yDFT[i] = in.yDFT[i];
    }
}

// load just the DFT part
void AWT::loadDFTtoAWT(AWT & in)
{
    if(n != in.n)   cerr<< "Error: DFT and AWT are not compatible" << endl; exit(1);
    for(int i=0; i < 4*n+4; i++)
    {
        yDFT[i] = in.yDFT[i];
    }
}

void AWT::loadFUNtoAWT(rawF & in)
{
    if(in.n != n)   cerr<< "Error: function and AWT are not compatible" << endl; exit(1);

    for(int i=0;     i < n+1;   i++)  y[i] = in.f[i+n];
    for(int i=n+1;   i < 3*n+4; i++)  y[i] = 0;
    for(int i=3*n+4; i < 4*n+4; i++)  y[i] = in.f[i-3*n-4];
}

///////////////////////////////////////////////////////////////////////////
////     Read out the data from the program to a predefined table     ////
/////////////////////////////////////////////////////////////////////////
void AWT::exportAWTasAWT(string & name, int div, string & output_mode)
{
    if(output_mode == "old")
    {
        ofstream output(name.c_str());
        for(  int i=0;   i < 4*n+4;  i=i+div)   output << i << "  " << real(y[i]) << "  " << imag(y[i]) << "  " << real(yDFT[i]) << "  " << imag(yDFT[i]) << "" << endl;
    }
    if(output_mode == "cx11")
    {
        ofstream output;
        output.open(name);
        for(  int i=0;   i < 4*n+4;  i=i+div)   output << i << "  " << real(y[i]) << "  " << imag(y[i]) << "  " << real(yDFT[i]) << "  " << imag(yDFT[i]) << "" << endl;
    }

}

///////////////////////////////////////////////////////////////////////////
////     Read out the data from the program to a predefined table     ////
/////////////////////////////////////////////////////////////////////////
void AWT::exportAWTasFUN(string & name, int div, double range, double mult, string & output_mode)
{
    int limit=n*range/xMax;
    //cout << name << " " << n << " " << range << " " << xMax << endl;
    #ifdef OLD

        ofstream output(name.c_str());
        // first the negative values are printed, then, the positive values are printed
        for(int i=4*n+4-limit; i<4*n+4; i=i+div)    output << (i-4*n -4)*xMax/n << "  " << mult*real(y[i]) << "  " << mult*imag(y[i]) << "  " << endl;
        for(int i=0;           i<limit; i=i+div)    output << i*xMax/n << "  " << mult*real(y[i]) << "  " << mult*imag(y[i]) << "  " << endl;

    #endif // OLD
    #ifdef CX11

        ofstream output;
        output.open(name);
        // first the negative values are printed, then, the positive values are printed
        for(int i=4*n+4-limit; i<4*n+4; i=i+div)    output << (i-4*n -4)*xMax/n << "  " << mult*real(y[i]) << "  " << mult*imag(y[i]) << "  " << endl;
        for(int i=0;           i<limit; i=i+div)    output << i*xMax/n << "  " << mult*real(y[i]) << "  " << mult*imag(y[i]) << "  " << endl;

    #endif // CX11
}

void AWT::exportDFTasFUN(string & name, int div, double range, string & output_mode)
{
    int limit=n*range/xMax;

    #ifdef OLD

        ofstream output(name.c_str());
        // first the negative values are printed, then, the positive values are printed
        for(int i=4*n+4-limit; i<4*n+4; i=i+div)    output << (i-4*n -4)*xMax/n << "  " << real(yDFT[i]) << "  " << imag(yDFT[i]) << "  " << endl;
        for(int i=0;           i<limit; i=i+div)    output << i*xMax/n << "  " << real(yDFT[i]) << "  " << imag(yDFT[i]) << "  " << endl;

    #endif // OLD
    #ifdef CX11

        ofstream output;
        output.open(name);
        // first the negative values are printed, then, the positive values are printed
        for(int i=4*n+4-limit; i<4*n+4; i=i+div)    output << (i-4*n -4)*xMax/n << "  " << real(yDFT[i]) << "  " << imag(yDFT[i]) << "  " << endl;
        for(int i=0;           i<limit; i=i+div)    output << i*xMax/n << "  " << real(yDFT[i]) << "  " << imag(yDFT[i]) << "  " << endl;

    #endif // CX11
}


////////////////////////////////////////////////////////////////////////
////              To clean up very small entries                   ////
///////////////////////////////////////////////////////////////////////
void AWT::cleanOUT(int exp)
{
    for(int ii=0; ii < nn; ii++)
    {
        if( abs(imag(yDFT[ii])) <= pow(10,exp))
        {
            complex<double> u(real(yDFT[ii]), 0);
            yDFT[ii]= u;
        }
        if( abs(real(yDFT[ii])) <= pow(10,exp))
        {
            complex<double> u(0, imag(yDFT[ii]));
            yDFT[ii]= u;
        }
        if(    abs(imag(y[ii])) <= pow(10,exp))
        {
            complex<double> u(real(y[ii]), 0);
            y[ii]= u;
        }
        if(    abs(real(y[ii])) <= pow(10,exp))
        {
            complex<double> u(0, imag(y[ii]));
            y[ii]= u;
        }
    }
}


//////////////////////////////////////////////////////////////////////////
/////////////            OPERATIONS ON AWTs        //////////////////////
////////////////////////////////////////////////////////////////////////


// Kernel3 for KRAMMERS KRONIG
void AWT::Kernel3()
{
    y[1] =  0;
    y[1] =  4*log(3.0/2.0);
    y[2] = 10*log(4.0/3.0) -6*log(3.0/2.0);

    int i;
    for(i=3; i<n+1; i++)
    {
        if(i<1000)
        {
            y[i] =     (1-i*i)*(i-2)         *atanh(1.0/(1 - 2*i))
                   +   (1-i*i)*(i+2)         *atanh(1.0/(1 + 2*i))
                   +   (i*i*i-6*i*i+11*i-6)  *atanh(1.0/(3 - 2*i))/3.0
                   +   (i+3)*(i*i+3*i+2)     *atanh(1.0/(3 + 2*i))/3.0;
        }
        else    y[i] = 1.0/i;

    }
    for(i=n+1;   i<3*n+4; i++)         y[i] =  0;
    for(i=3*n+4; i<4*n+4; i++)         y[i] =  -y[4*n+4-i];
    for(i=0;     i<4*n+4; i++)      yDFT[i] =  0;

}

void AWT::cosine(double constant)
{
    int i;
                                          y[0] = cos(0);
    for( i = 1;      i < n+1;      i++ )  y[i] = cos(constant*i*xMax/n);
    for( i = n+1;    i < 3*n + 4;  i++ )  y[i] = 0;
    for( i = 3*n+4;  i < 4*n + 4;  i++ )  y[i] = cos(-constant*(4*n+4-i)*xMax/n);

}
// multiply AWT with its Fermi-Dirac function
void AWT::normalFDtimesIM(AWT & X, AWT & f_FD)
{
    for(int i=0; i<X.nn; i++){y[i] = imag(X.y[i])*f_FD.y[i];  }
}

// multiply AWT  with its Fermi-Dirac function in the inverse order of i
void AWT::inverseFDtimesIM(AWT & X, AWT & f_FD)
{

    int i;
                                              y[0] = imag(X.y[0]) * f_FD.y[0];
    for( i = 1;        i < X.n+1;      i++ ){ y[i] = imag(X.y[4*X.n+4-i]) * f_FD.y[i];}
    for( i = X.n+1;    i < 3*X.n + 4;  i++ ){ y[i] = 0;}
    for( i = 3*X.n+4;  i < 4*X.n + 4;  i++ ){ y[i] = imag(X.y[4*X.n+4-i]) * f_FD.y[i];}
}

// multiply AWT with its Bose-Einstein function
void AWT::normalBEtimesIM(AWT & X, AWT & f_BE)
{
    for(int i=0; i<X.nn; i++){  y[i] = imag(X.y[i])*f_BE.y[i];  }
}

// multiply AWT  with its Bose-Einstein function in the inverse order of i
void AWT::inverseBEtimesIM(AWT & X, AWT & f_BE)
{
    int i;
                                              y[0] = imag(X.y[0]) * f_BE.y[0];
    for( i = 1;        i < X.n+1;      i++ ){ y[i] = imag(X.y[4*X.n+4-i]) * f_BE.y[i];}
    for( i = X.n+1;    i < 3*X.n + 4;  i++ ){ y[i] = 0;}
    for( i = 3*X.n+4;  i < 4*X.n + 4;  i++ ){ y[i] = imag(X.y[4*X.n+4-i]) * f_BE.y[i];}
}

// conjugate AWT body with specifying the output
void AWT::conjugateY(AWT & inX)
{
    for(int i=0; i < inX.nn; i++)  y[i]  =  conj(inX.y[i]);
}

// conjugate just the DFT part
void AWT::conjugateDFT(AWT & inX)
{
    for(int i=0; i < inX.nn; i++)  yDFT[i]  =  conj(inX.yDFT[i]);
}

void AWT::multiplyAWT(AWT & timesX, AWT & timesY)
{
    for(int i=0; i<timesX.nn; i++){  y[i] = timesX.y[i] * timesY.y[i]; }
}
// real part is deleted
void AWT::deleteReal(AWT & Gin)
{
    for(int i=0; i < nn; i++)
    {
        complex<double> u(0, imag(Gin.y[i]) );
        y[i]= u;
    }
}

// imag part is deleted
void AWT::deleteImag(AWT & Gin)
{
    for(int i=0; i < nn; i++)
    {
        complex<double> u(real(Gin.y[i]), 0);
        y[i]= u;
    }
}


// Krammers Kronig for the upper complex plane
void AWT::KrammersKronig(AWT & G, AWT & Kernel, AWT & ReG, AWT & ImG, AWT & aux1, AWT & aux2)
{
    double Pi=3.14159265359;

    int i;
    for(i=0; i<4*G.n+4; i++)
    {
        complex<double> u( imag(G.y[i]), 0 );
        ImG.y[i] = u;
    }

    ReG.convolutionAWT(ImG, Kernel, 1, 1, aux1, aux2);

    for(i=0; i<4*G.n+4; i++)
    {
        complex<double> u( -G.n*real(ReG.y[i])/Pi/G.xMax, imag(G.y[i]) );
        y[i] = u;
    }

}

// Krammers Kronig for the lower complex plane
void AWT::KrammersKronigDown(AWT & G, AWT & Kernel, AWT & ReG, AWT & ImG, AWT & aux1, AWT & aux2)
{
    double Pi=3.14159265359;

    int i;
    for(i=0; i<4*G.n+4; i++)
    {
        complex<double> u( imag(G.y[i]), 0 );
        ImG.y[i] = u;
    }

    ReG.convolutionAWT(ImG, Kernel, 1, 1, aux1, aux2);

    for(i=0; i<4*G.n+4; i++)
    {
        complex<double> u( G.n*real(ReG.y[i])/Pi/G.xMax, imag(G.y[i]) );
        y[i] = u;
    }
}

void AWT::multiplyAWT(AWT & X, complex<double> u)
{
    int i;
    for(i=0; i<X.nn; i++)   y[i] = u*X.y[i];
}

void AWT::derivativeAWT(AWT & X)
{
    // numerical derivation of given AWT

    // positive frequencies
                                                   y[0]       = ( X.y[1]       - X.y[4*X.n+3]   ) / ( 2 * X.xMax / X.n);
    for(int i = 1;        i < X.n;        i++ )    y[i]       = ( X.y[i+1]     - X.y[i-1]       ) / ( 2 * X.xMax / X.n);
                                                   y[X.n]     = ( X.y[X.n]     - X.y[X.n-1]     ) / (     X.xMax / X.n);

    // zero padding
    for(int i = X.n+1;    i < 3*X.n + 4;  i++ )    y[i]       = 0;

    // negative frequencies
                                                   y[3*X.n+4] = ( X.y[3*X.n+5] - X.y[3*X.n+4]   ) / (     X.xMax / X.n);
    for(int i = 3*X.n+5;  i < 4*X.n + 3;  i++ )    y[i]       = ( X.y[i+1]     - X.y[i-1]       ) / ( 2 * X.xMax / X.n);
                                                   y[4*X.n+3] = ( X.y[0]       - X.y[4*X.n+2]   ) / ( 2 * X.xMax / X.n);
}

///////////////////////////////////////////////////////////////////////////
//////////             OPERATIONS BETWEEN AWTs                  //////////
/////////////////////////////////////////////////////////////////////////

// Multiply two AWTs
void AWT::doubleNormal(AWT & X, AWT & Y)
{
    int i;
    for( i=0;  i<4*X.n + 4;  i++)     y[i] = X.y[i] * Y.y[i];


}

// Double GG
void AWT::doubleInverse(AWT & Gin)
{
    int i;
                                            y[0] = Gin.y[0] * conj(Gin.y[0]);
    for( i=1;          i<Gin.n+1;      i++){y[i] = Gin.y[i] * conj(Gin.y[4*Gin.n+4-i]);}
    for( i=Gin.n+1;    i<3*Gin.n + 4;  i++){y[i] = 0;}
    for( i=3*Gin.n+4;  i<4*Gin.n + 4;  i++){y[i] = Gin.y[i] * conj(Gin.y[4*Gin.n+4-i]);}
}



// convolve two AWTs
void AWT::convolutionAWT(AWT & inX, AWT & inY, int sig1, int sig2, AWT & helpX, AWT & helpY)
{
    inX.forwardDFT();    // DFT of inX if not available
    inY.forwardDFT();    // DFT of inY if not available

    helpX.loadAWTtoAWT(inX);
    helpY.loadAWTtoAWT(inY);

    if(sig1*sig2==-1)
    {
        helpX.conjugateY(helpX);
        helpX.forwardDFT();
        helpX.conjugateDFT(helpX);
    }

    if(sig1==-1)
    {
        helpY.conjugateY(helpY);
        helpY.forwardDFT();
        helpY.conjugateDFT(helpY);
    }

    double dx = inX.xMax/inX.n;
    for(int i=0; i<inX.nn; i++)
    {
        yDFT[i] = helpX.yDFT[i] * helpY.yDFT[i] * dx;
    }

    // the output AWT is always set closed for next convolutionAWT instance
    backwardDFT();

};



// perform Fermi-Dirac summation over two AWTs
void AWT::fermMatsubaraImP(AWT & X, AWT & Y, int sig1, int sig2, AWT & f_FD, AWT & x1, AWT & x2, AWT & y1, AWT & y2, AWT & aux1, AWT & aux2 )
{
    double Pi=3.14159265359;

    // FIRST PART of the fermi sum of two AWTs is calculated
    x1.normalFDtimesIM(X, f_FD);
    y1.deleteReal(Y);
    x1.convolutionAWT(x1, y1, sig1, -sig2, aux1, aux2);           // integral over real arguments
    x1.multiplyAWT(x1, -sig1/Pi);                        // multiply with proper factor

    // SECOND PART of the bose sum of two AWTs is calculated
    if(sig2 == 1)         x2.normalFDtimesIM(Y, f_FD);
    if(sig2 ==-1)         x2.inverseFDtimesIM(Y, f_FD);
    y2.deleteReal(X);
    x2.convolutionAWT(x2, y2, -sig1*sig2, -1, aux1, aux2);       // integral over real arguments
    x2.multiplyAWT(x2, sig1/Pi);                      // multiply with proper factor

    // first and second part are added tgether and stored in the out AWT
    for(int i=0; i < X.nn; i++)         y[i] = x1.y[i] + x2.y[i];
}



// perform Bose-Einstein summation over two AWTs
void AWT::boseMatsubaraImP(AWT & Xb, AWT & Y, int sig1, int sig2, AWT & f_FD, AWT & f_BE, AWT & x1, AWT & x2, AWT & y1, AWT & y2, AWT & aux1, AWT & aux2 )
{
    double Pi=3.14159265359;

    // FIRST PART of the fermi sum of two AWTs is calculated
    x1.normalBEtimesIM(Xb, f_BE);
    y1.deleteReal(Y);
    if(sig1 == -1)  y1.conjugateY(y1);                // conjugate the 1st part if X must be taken below the real axis

    x1.convolutionAWT(x1, y1, sig1, -sig2, aux1, aux2);           // integral over real arguments
    x1.multiplyAWT(x1, 1/Pi);                        // multiply with proper factor

    // SECOND PART of the bose sum of two AWTs is calculated
    if(sig2 == 1)         x2.normalFDtimesIM(Y, f_FD);
    else                  x2.inverseFDtimesIM(Y, f_FD);
    if(sig2 == -1)        x2.multiplyAWT(x2, -1);

    y2.deleteReal(Xb);
    if(sig1*sig2==1)    y2.conjugateY(y2);               // conjugate the 2nd part if X must be taken below the real axis

    x2.convolutionAWT(x2, y2, -sig1*sig2, -1, aux1, aux2);       // integral over real arguments
    x2.multiplyAWT(x2, -1/Pi);                      // multiply with proper factor

    // first and second part are added tgether and stored in the out AWT
    for(int i=0; i < Xb.nn; i++)         y[i] = x1.y[i] + x2.y[i];
}



