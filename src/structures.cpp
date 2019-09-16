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


AUX::AUX()
{
    AWT aux1;
    AWT aux2;
    AWT aux3;
    AWT aux4;
    AWT aux5;
    AWT aux6;
}

void AUX::initialize(int n, double x, double kT)
{
    aux1.initializeAWT(n, x, kT);
    aux2.initializeAWT(n, x, kT);
    aux3.initializeAWT(n, x, kT);
    aux4.initializeAWT(n, x, kT);
    aux5.initializeAWT(n, x, kT);
    aux6.initializeAWT(n, x, kT);
}

AUX::~AUX()
{

}

void info::importing(string name)
{

    ifstream ListFile;

    // 1. determine number of strings in the input file
    ListFile.open (name.c_str());
    string item;                    // temporary storege for strings of input file

    if (ListFile.fail()) {cerr<< "Error opening the file" << endl; exit(1);}
    int entries = 0;                 // entries stores number of strings
    while (!ListFile.eof())
    {
        ListFile >> item;
        entries = entries + 1;
    }
    ListFile.close();

    string input_entries[entries];

    // 2. store all strings from the input
    ListFile.open (name.c_str());   // list contains all the file names to be processed
    for(int i=0; i<entries; i++)    ListFile >> input_entries[i];
    ListFile.close();

    // 3. assign strings from the input to info variables
    for(int i=0; i<entries; i++)
    {

        if (input_entries[i].compare("delta:")         == 0)         delta = boost::lexical_cast<double>(input_entries[i+1]);
        if (input_entries[i].compare("mu:")            == 0)            mu = boost::lexical_cast<double>(input_entries[i+1]);
        if (input_entries[i].compare("gap:")           == 0)           gap = boost::lexical_cast<double>(input_entries[i+1]);
        if (input_entries[i].compare("gammaS:")        == 0)        gammaS = boost::lexical_cast<double>(input_entries[i+1]);
        if (input_entries[i].compare("phi:")           == 0)           phi = boost::lexical_cast<double>(input_entries[i+1]);
        if (input_entries[i].compare("model:")         == 0)         model = boost::lexical_cast<string>(input_entries[i+1]);


        if ( (input_entries[i].compare("kT") == 0) &&  (input_entries[i+1].compare("min:") == 0)   )
                                                                    kT_min = boost::lexical_cast<double>(input_entries[i+2]);

        if ( (input_entries[i].compare("kT") == 0) &&  (input_entries[i+1].compare("increment:") == 0)   )
                                                              kT_increment = boost::lexical_cast<double>(input_entries[i+2]);

        if ( (input_entries[i].compare("kT") == 0) &&  (input_entries[i+1].compare("max:") == 0)   )
                                                                    kT_max = boost::lexical_cast<double>(input_entries[i+2]);



        if ( (input_entries[i].compare("U") == 0) &&  (input_entries[i+1].compare("min:") == 0)   )
                                                                    U_min = boost::lexical_cast<double>(input_entries[i+2]);

        if ( (input_entries[i].compare("U") == 0) &&  (input_entries[i+1].compare("increment:") == 0)   )
                                                              U_increment = boost::lexical_cast<double>(input_entries[i+2]);

        if ( (input_entries[i].compare("U") == 0) &&  (input_entries[i+1].compare("max:") == 0)   )
                                                                    U_max = boost::lexical_cast<double>(input_entries[i+2]);



        if ( (input_entries[i].compare("x") == 0) &&  (input_entries[i+1].compare("min:") == 0)   )
                                                                    x_min = boost::lexical_cast<double>(input_entries[i+2]);

        if ( (input_entries[i].compare("x") == 0) &&  (input_entries[i+1].compare("increment:") == 0)   )
                                                              x_increment = boost::lexical_cast<double>(input_entries[i+2]);

        if ( (input_entries[i].compare("x") == 0) &&  (input_entries[i+1].compare("max:") == 0)   )
                                                                    x_max = boost::lexical_cast<double>(input_entries[i+2]);



        if ( (input_entries[i].compare("h") == 0) &&  (input_entries[i+1].compare("min:") == 0)   )
                                                                    h_min = boost::lexical_cast<double>(input_entries[i+2]);

        if ( (input_entries[i].compare("h") == 0) &&  (input_entries[i+1].compare("increment:") == 0)   )
                                                              h_increment = boost::lexical_cast<double>(input_entries[i+2]);

        if ( (input_entries[i].compare("h") == 0) &&  (input_entries[i+1].compare("max:") == 0)   )
                                                                    h_max = boost::lexical_cast<double>(input_entries[i+2]);


        if (input_entries[i].compare("xMax:")    == 0)              xMax = boost::lexical_cast<double>(input_entries[i+1]);
        if (input_entries[i].compare("n:")       == 0)                 n = boost::lexical_cast<int>(input_entries[i+1]);
        if (input_entries[i].compare("display:") == 0)           display = boost::lexical_cast<int>(input_entries[i+1]);
        if (input_entries[i].compare("range:")   == 0)             range = boost::lexical_cast<int>(input_entries[i+1]);


        if ( (input_entries[i].compare("print") == 0) &&  (input_entries[i+1].compare("mode:") == 0)   )
                                                              print_mode = boost::lexical_cast<int>(input_entries[i+2]);

        if ( (input_entries[i].compare("export") == 0) &&  (input_entries[i+1].compare("mode:") == 0)   )
                                                              output_mode = boost::lexical_cast<string>(input_entries[i+2]);

        if (input_entries[i].compare("physics:") == 0)            physics = boost::lexical_cast<string>(item);


    }


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
