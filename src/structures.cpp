#include "AWT.h"
#include "structures.h"

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

#include <boost/lexical_cast.hpp>
#include <boost/filesystem.hpp>
#include <iomanip> // setprecision


auxiliary::auxiliary(){    aux_initialize(0, 0, 0);}

void auxiliary::aux_initialize(int n, double x, double kT)
{
    aux1.initializeAWT(n, x, kT);
    aux2.initializeAWT(n, x, kT);
    aux3.initializeAWT(n, x, kT);
    aux4.initializeAWT(n, x, kT);
    aux5.initializeAWT(n, x, kT);
    aux6.initializeAWT(n, x, kT);
}

auxiliary::~auxiliary()
{
        fftw_free(aux1.y);
    fftw_free(aux1.yINTER);
    fftw_destroy_plan(aux1.forwardFFT);
    fftw_destroy_plan(aux1.backwardFFT);

        fftw_free(aux2.y);
    fftw_free(aux2.yINTER);
    fftw_destroy_plan(aux2.forwardFFT);
    fftw_destroy_plan(aux2.backwardFFT);

        fftw_free(aux3.y);
    fftw_free(aux3.yINTER);
    fftw_destroy_plan(aux3.forwardFFT);
    fftw_destroy_plan(aux3.backwardFFT);

        fftw_free(aux4.y);
    fftw_free(aux4.yINTER);
    fftw_destroy_plan(aux4.forwardFFT);
    fftw_destroy_plan(aux4.backwardFFT);

        fftw_free(aux5.y);
    fftw_free(aux5.yINTER);
    fftw_destroy_plan(aux5.forwardFFT);
    fftw_destroy_plan(aux5.backwardFFT);

        fftw_free(aux6.y);
    fftw_free(aux6.yINTER);
    fftw_destroy_plan(aux6.forwardFFT);
    fftw_destroy_plan(aux6.backwardFFT);
}

void info::importing(string name)
{

    ifstream ListFile;
    ListFile.open (name.c_str());   // list contains all the file names to be processed
    string item;                    // used for temporary store the data

    if (ListFile.fail()) {cerr<< "Error opening the file" << endl; exit(1);}

    for(int i=0; i<58; i++)
    {
        ListFile >> item;

        if(i == 3)   delta        = boost::lexical_cast<double>(item);
        if(i == 5)   mu           = boost::lexical_cast<double>(item);
        if(i == 7)   gap          = boost::lexical_cast<double>(item);
        if(i == 9)   gammaS       = boost::lexical_cast<double>(item);
        if(i == 11)  phi          = boost::lexical_cast<double>(item);


        if(i == 13)  model        = boost::lexical_cast<string>(item);

        if(i == 18)  kT_min       = boost::lexical_cast<double>(item);
        if(i == 21)  kT_increment = boost::lexical_cast<double>(item);
        if(i == 24)  kT_max       = boost::lexical_cast<double>(item);

        if(i == 27)  U_min        = boost::lexical_cast<double>(item);
        if(i == 30)  U_increment  = boost::lexical_cast<double>(item);
        if(i == 33)  U_max        = boost::lexical_cast<double>(item);

        if(i == 36)  x_min        = boost::lexical_cast<double>(item);
        if(i == 39)  x_increment  = boost::lexical_cast<double>(item);
        if(i == 42)  x_max        = boost::lexical_cast<double>(item);

        if(i == 48)  h_min        = boost::lexical_cast<double>(item);
        if(i == 48)  h_increment  = boost::lexical_cast<double>(item);
        if(i == 51)  h_max        = boost::lexical_cast<double>(item);


        if(i == 55)  xMax         = boost::lexical_cast<double>(item);
        if(i == 57)  n            = boost::lexical_cast<int>(item);

        if(i == 62)  display      = boost::lexical_cast<int>(item);
        if(i == 64)  range        = boost::lexical_cast<int>(item);
        if(i == 67)  print_mode   = boost::lexical_cast<int>(item);
        if(i == 70)  output_mode  = boost::lexical_cast<string>(item);
        if(i == 72)  physics      = boost::lexical_cast<string>(item);

    }

    ListFile.close();
}



void info::precisions()
{
    // the number of decimals in kT_min is determined
    kT_precision = 0;
    U_precision = 2;
    x_precision = 2;

    /*double kT_ = kT_min;
    while(kT_ < 1)
    {
        kT_ = 10*kT_;
        kT_precision = kT_precision +1;
    }

    // the number of decimals in U_min is determined
    U_precision = 0;
    U_precision = 1;

    // the number of decimals in x_increment is determined
    x_precision = 0;
    double x_ = x_increment;
    while(x_ < 1)
    {
        x_ = 10*kT_;
        x_precision = x_precision +1;
    }
    */
}

void  info::iterations()
{
    if(kT_increment == 0)   kT_iterations = 1;
    else                    kT_iterations = ( kT_max - kT_min) / kT_increment + 1;

    // variable for maximal number of U iteration is initialized and set
    if(U_increment == 0)    U_iterations = 1;
    else                    U_iterations  = ( U_max -  U_min) /  U_increment + 1;

    // variable for maximal number of kT iteration is initialized and set
    if(x_increment == 0)    x_iterations = 1;
    else                    x_iterations  = ( x_max - x_min  ) /  x_increment + 1;

    if(h_increment == 0)    h_iterations = 1;
    else                    h_iterations = (h_max - h_min)/h_increment;

    /////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////                            SYSTEM PROTECTION                                            ///////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////

    if( (kT_iterations > 10000) || (U_iterations > 10000) || (x_iterations > 10000)  )
    {
        cout << endl;
        cout << "Program stopped! Too many iterations. Possible bug?" << endl;
        exit(1);
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void cout_vect(vect X)
{
    cout << X.p << " " << X.m << endl;
}

void output_vect(vect X, string name)
{
    ofstream output;
    output.open(name, std::ios_base::app);
    output << X.p << " " << X.m << endl;
    output << endl;
}

void set_matrix_real(matrix & X, complex<double> val)
{
    X.pp = val;
    X.mp = val;
    X.pm = val;
    X.mm = val;
}

void set_matrix(matrix & X, matrix & Y)
{
    X.pp = Y.pp;
    X.mp = Y.mp;
    X.pm = Y.pm;
    X.mm = Y.mm;
}

void aver_matrix(matrix & X, matrix & oldX)
{
    X.pp = ( X.pp + oldX.pp ) / 2.0;
    X.mp = ( X.mp + oldX.mp ) / 2.0;
    X.pm = ( X.pm + oldX.pm ) / 2.0;
    X.mm = ( X.mm + oldX.mm ) / 2.0;
}

void cout_matrix(matrix X)
{
    cout << X.pp << " " << X.mp << " " << X.pm << " " << X.mm << endl;
}

void output_matrix(matrix X, string name)
{
    ofstream output;
    output.open(name, std::ios_base::app);
    output << X.pp << " " << X.mp << " " << X.pm << " " << X.mm << endl;
    output << endl;
}



void set_matrixAWT(matrixAWT & X, double U)
{
    X.pp.set_real(U);
    X.mp.set_real(U);
    X.pm.set_real(U);
    X.mm.set_real(U);
}

void copy_matrixAWT(matrixAWT & in, matrixAWT & out)
{
    for(int i =0; i<in.pp.nn; i++)
    {
        out.pp.y[i] = in.pp.y[i];
        out.mp.y[i] = in.mp.y[i];
        out.pm.y[i] = in.pm.y[i];
        out.mm.y[i] = in.mm.y[i];
    }
}


void aver_matrixAWT(matrixAWT & in, matrixAWT & out)
{
    for(int i =0; i<in.pp.nn; i++)
    {
        out.pp.y[i] = ( in.pp.y[i] + out.pp.y[i] ) / 2.0;
        out.mp.y[i] = ( in.mp.y[i] + out.mp.y[i] ) / 2.0;
        out.pm.y[i] = ( in.pm.y[i] + out.pm.y[i] ) / 2.0;
        out.mm.y[i] = ( in.mm.y[i] + out.mm.y[i] ) / 2.0;
    }
}

void export_matrixAWT(matrixAWT & X, double range, string name, info in)
{

    int div = in.display;
    double mult = 1.0;
    int limit=in.n*range/in.xMax;

    ofstream output;
    output.open(name);


    // first the negative values are printed, then, the positive values are printed
    for(int i=4*in.n+4-limit; i<4*in.n+4; i=i+div)
    {
        output << (i-4*in.n -4)*in.xMax/in.n << "  ";
        output << mult*real(X.pp.y[i]) << "  " << mult*imag(X.pp.y[i]) << "  ";
        output << mult*real(X.mp.y[i]) << "  " << mult*imag(X.mp.y[i]) << "  ";
        output << mult*real(X.pm.y[i]) << "  " << mult*imag(X.pm.y[i]) << "  ";
        output << mult*real(X.mm.y[i]) << "  " << mult*imag(X.mm.y[i]) << "  ";
        output << endl;
    }
    for(int i=0;           i<limit; i=i+div)
    {
        output << i*in.xMax/in.n << "  ";
        output << mult*real(X.pp.y[i]) << "  " << mult*imag(X.pp.y[i]) << "  ";
        output << mult*real(X.mp.y[i]) << "  " << mult*imag(X.mp.y[i]) << "  ";
        output << mult*real(X.pm.y[i]) << "  " << mult*imag(X.pm.y[i]) << "  ";
        output << mult*real(X.mm.y[i]) << "  " << mult*imag(X.mm.y[i]) << "  ";
        output << endl;
    }

}



 void cout_tensor(tensor X)
{
    cout << X.ppp << endl;
    cout << X.mpp << " " << X.pmp << " " << X.ppm << endl;
    cout << X.pmm << " " << X.mpm << " " << X.mmp << endl;
    cout << X.mmm << endl;
}

void output_tensor(tensor & X, string name)
{
    ofstream output;
    output.open(name, std::ios_base::app);
    output << X.ppp << endl;
    output << X.mpp << " " << X.pmp << " " << X.ppm << endl;
    output << X.pmm << " " << X.mpm << " " << X.mmp << endl;
    output << X.mmm << endl;
    output << endl;

    output.close();
}

void copy_tensor(tensor & in, tensor & out)
{
    out.ppp = in.ppp;

    out.ppm = in.ppm;
    out.pmp = in.pmp;
    out.mpp = in.mpp;

    out.mmp = in.mmp;
    out.mpm = in.mpm;
    out.pmm = in.pmm;

    out.mmm = in.mmm;
}


void aver_tensor(tensor & in, tensor & out)
{
    out.ppp = ( out.ppp + in.ppp )/2.0;

    out.ppm = ( out.ppm + in.ppm )/2.0;
    out.pmp = ( out.pmp + in.pmp )/2.0;
    out.mpp = ( out.mpp + in.mpp )/2.0;

    out.mmp = ( out.mmp + in.mmp )/2.0;
    out.mpm = ( out.mpm + in.mpm )/2.0;
    out.pmm = ( out.pmm + in.pmm )/2.0;

    out.mmm = ( out.mmm + in.mmm )/2.0;
}

void multiply_tensor(tensor & in, complex<double> val)
{
    in.ppp = ( in.ppp ) * val;

    in.ppm = ( in.ppm ) * val;
    in.pmp = ( in.pmp ) * val;
    in.mpp = ( in.mpp ) * val;

    in.mmp = ( in.mmp ) * val;
    in.mpm = ( in.mpm ) * val;
    in.pmm = ( in.pmm ) * val;

    in.mmm = ( in.mmm ) * val;
}

void integrate_tensor(AWT & in, tensor & out, int ind1, int ind2, int ind3)
{
    if( ind1 ==  1 && ind2 ==  1 && ind3 ==  1)  { out.ppp = Simpson(in); }

    if( ind1 ==  1 && ind2 ==  1 && ind3 == -1)  { out.ppm = Simpson(in);  }
    if( ind1 ==  1 && ind2 == -1 && ind3 ==  1)  { out.pmp = Simpson(in);  }
    if( ind1 == -1 && ind2 ==  1 && ind3 ==  1)  { out.mpp = Simpson(in);  }

    if( ind1 == -1 && ind2 == -1 && ind3 ==  1)  { out.mmp = Simpson(in);  }
    if( ind1 == -1 && ind2 ==  1 && ind3 == -1)  { out.mpm = Simpson(in);  }
    if( ind1 ==  1 && ind2 == -1 && ind3 == -1)  { out.pmm = Simpson(in);  }

    if( ind1 == -1 && ind2 == -1 && ind3 == -1)  { out.mmm = Simpson(in);  }
}

void subtract_tensors(tensor & out, tensor & in1, tensor & in2)
{
    out.ppp = in1.ppp - in2.ppp;

    out.ppm = in1.ppm - in2.ppm;
    out.pmp = in1.pmp - in2.pmp;
    out.mpp = in1.mpp - in2.mpp;

    out.mmp = in1.mmp - in2.mmp;
    out.mpm = in1.mpm - in2.mpm;
    out.pmm = in1.pmm - in2.pmm;

    out.mmm = in1.mmm - in2.mmm;
}
