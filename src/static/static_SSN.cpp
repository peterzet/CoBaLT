#include "static_SSN.h"
#include "../physics.h"
#include "../convert.h"
#include "../global.h"
#include "../structures.h"


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
#include <boost/filesystem.hpp>
#include <iomanip> // setprecision

using namespace std;


void welcome_screen_SSN(info & in, string dirname)
{
    // on screen output of starting TIME and basic model parameters
    cout << __TIME__ << ",   " << __DATE__ << endl;
    cout << "---------------------------------------"   << endl;
    cout << "static vertex with four sign combinations" << endl;
    cout << "equations according to Janis from end of November 2018  " << endl;
    cout << "---------------------------------------"   << endl;
    cout << "Coulomb interaction set to U = " << in.U   << endl;
    cout << "magnetic field set to h = "      << in.h   << endl;
    cout << "mu is set to "                   << in.mu  << endl;
    cout << "kT is set to "                   << in.kT  << endl;
    cout << "---------------------------------------"   << endl;
    cout << "this is the three-terminal set up problem with two BCS and one normal lead " << endl;
    cout << "gap is set to gap = "            << in.gap << endl;
    cout << "angle is set to phi = "          << in.phi << endl;
    cout << "---------------------------------------"   << endl;
    cout << "creating directory "            << dirname << endl;
    cout << "---------------------------------------"   << endl;
}

void lorentz_propagator(info in, AWT & G)
{
    double n = G.n;
    double x = G.xMax;

    for(int i=0;       i<n+1; i++ )     {complex<double> u(     i*x/n        + in.mu, 0 ); G.y[i] = 1.0/u; }
    for(int i=n+1;   i<3*n+3; i++ )     G.y[i] = 0;
    for(int i=3*n+4; i<4*n+4; i++ )     {complex<double> u( (i - 4*n -4)*x/n + in.mu, 0 ); G.y[i] = 1.0/u; }

    for(int i=0; i<4*n+4; i++ )         if(i*x/n + in.mu == 0)   G.y[i] = 1e24;
}

double hybridization_SIAM_const(info in, int i)
{
    //cout << "hybridization function of SIAM with constant DOS is inserted" << endl;
    return in.delta;
}

double hybridization_SSN_transformed(info in, int i)
{
    //cout << "hybridization function of the SSN structure with rotation trick inserted" << endl;
    return in.delta;
}

void insert_hybridization(info in, AWT & G)
{
    complex<double> unity(0,1);

    // default hybridization function: SIAM with constant DOS
    double (*hybridPtr)(info, int) = hybridization_SIAM_const;

    // hybridization function of the SSN structure with rotation trick
    if(in.model.compare("SSNtransformed")==0)  hybridPtr = hybridization_SSN_transformed;


    for(int i=0;        i<in.n+1;   i++ )   {complex<double> u = 1.0/G.y[i] + unity*hybridPtr(in,i);  G.y[i] = 1.0/u; }
    for(int i=in.n+1;   i<3*in.n+3; i++ )   G.y[i] = 0;
    for(int i=3*in.n+4; i<4*in.n+4; i++ )   {complex<double> u = 1.0/G.y[i] + unity*hybridPtr(in,i);  G.y[i] = 1.0/u; }
}

void insert_dynamic_hybridization(info in, AWT & G, AWT & hybrid)
{
    complex<double> unity(0,1);

    for(int i=0;        i<in.n+1;   i++ )   {complex<double> u = 1.0/G.y[i] + hybrid.y[i];  G.y[i] = 1.0/u; }
    for(int i=in.n+1;   i<3*in.n+3; i++ )   G.y[i] = 0;
    for(int i=3*in.n+4; i<4*in.n+4; i++ )   {complex<double> u = 1.0/G.y[i] + hybrid.y[i];  G.y[i] = 1.0/u; }
}

void insert_thermal_selfenergy(info in, AWT & G, complex<double> self)
{
    for(int i=0;        i<in.n+1;   i++ )   G.y[i] = 1.0 / (1.0/G.y[i] - self );
    for(int i=in.n+1;   i<3*in.n+3; i++ )   G.y[i] = 0;
    for(int i=3*in.n+4; i<4*in.n+4; i++ )   G.y[i] = 1.0 / (1.0/G.y[i] - self );
}

void insert_dynamic_selfenergy(info in, AWT & G, AWT & self)
{
    for(int i=0;        i<in.n+1; i++ )     G.y[i] = 1.0 / (1.0/G.y[i] - self.y[i] );
    for(int i=in.n+1;   i<3*in.n+3; i++ )   G.y[i] = ( 0 , 0 );
    for(int i=3*in.n+4; i<4*in.n+4; i++ )   G.y[i] = 1.0 / (1.0/G.y[i] - self.y[i] );
}




////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// BOTH ud and du variants for mag field
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void calculateGGud_SSN(AWT & GdGu1pp, AWT & GdGu1pm, AWT & GdGu2pm, AWT & GdGu2mm, AWT & GTd, AWT & GTu, auxiliary & help, int out)
{
    if(out == 1)    cout << "GG";

    complex<double> unity(0,1);

    // BUBBLE GdGu1pp
    for(int j = 0;          j < GTu.n+1;    j++)   help.aux1.y[j] = help.FD.y[j] * GTu.y[j] / ( -2.0*Pi*unity );
    for(int j = GTu.n+1;    j < 3*GTu.n+4;  j++)   help.aux1.y[j] = 0;
    for(int j = 3*GTu.n+4;  j < 4*GTu.n+4;  j++)   help.aux1.y[j] = help.FD.y[j] * GTu.y[j] / ( -2.0*Pi*unity );


    for(int j = 0;          j < GTu.n+1;    j++)   help.aux2.y[j] = GTd.y[j];
    for(int j = GTu.n+1;    j < 3*GTu.n+4;  j++)   help.aux2.y[j] = 0;
    for(int j = 3*GTu.n+4;  j < 4*GTu.n+4;  j++)   help.aux2.y[j] = GTd.y[j];

    GdGu1pp.convolutionAWT(help.aux1, help.aux2, 1, -1, help.aux3, help.aux4);  // f[x] * g[x+a]
    //GdGu1pp.KrammersKronig(GdGu1pp, K3, aux1, aux2, aux3, aux4);


    // BUBBLE GdGu1pm
    for(int j = 0;          j < GTu.n+1;    j++)   help.aux1.y[j] = help.FD.y[j] * conj( GTu.y[j] ) / ( -2.0*Pi*unity );
    for(int j = GTu.n+1;    j < 3*GTu.n+4;  j++)   help.aux1.y[j] = 0;
    for(int j = 3*GTu.n+4;  j < 4*GTu.n+4;  j++)   help.aux1.y[j] = help.FD.y[j] * conj( GTu.y[j] ) / ( -2.0*Pi*unity );

    for(int j = 0;          j < GTu.n+1;    j++)   help.aux2.y[j] = GTd.y[j];
    for(int j = GTu.n+1;    j < 3*GTu.n+4;  j++)   help.aux2.y[j] = 0;
    for(int j = 3*GTu.n+4;  j < 4*GTu.n+4;  j++)   help.aux2.y[j] = GTd.y[j];

    GdGu1pm.convolutionAWT(help.aux1, help.aux2, 1, -1, help.aux3, help.aux4);  // f[x] * g[x+a]
    //GdGu1pm.KrammersKronig(GdGu1pm, K3, aux1, aux2, aux3, aux4);


    // BUBBLE GdGu2pm
    for(int j = 0;          j < GTu.n+1;    j++)   help.aux1.y[j] = help.FD.y[j] * GTd.y[j]  / ( -2.0*Pi*unity );
    for(int j = GTu.n+1;    j < 3*GTu.n+4;  j++)   help.aux1.y[j] = 0;
    for(int j = 3*GTu.n+4;  j < 4*GTu.n+4;  j++)   help.aux1.y[j] = help.FD.y[j] * GTd.y[j]  / ( -2.0*Pi*unity );

    for(int j = 0;          j < GTu.n+1;    j++)   help.aux2.y[j] = conj( GTu.y[j] );
    for(int j = GTu.n+1;    j < 3*GTu.n+4;  j++)   help.aux2.y[j] = 0;
    for(int j = 3*GTu.n+4;  j < 4*GTu.n+4;  j++)   help.aux2.y[j] = conj( GTu.y[j] );


    GdGu2pm.convolutionAWT(help.aux1, help.aux2, -1, -1, help.aux3, help.aux4);  // f[x] * g[-a+x]
    //GdGu2pm.KrammersKronig(GdGu2pm, K3, aux1, aux2, aux3, aux4);


    // BUBBLE GdGu2mm
    for(int j = 0;          j < GTu.n+1;    j++)   help.aux1.y[j] = help.FD.y[j] * conj( GTd.y[j] )  / ( -2.0*Pi*unity );
    for(int j = GTu.n+1;    j < 3*GTu.n+4;  j++)   help.aux1.y[j] = 0;
    for(int j = 3*GTu.n+4;  j < 4*GTu.n+4;  j++)   help.aux1.y[j] = help.FD.y[j] * conj( GTd.y[j] )  / ( -2.0*Pi*unity );

    for(int j = 0;          j < GTu.n+1;    j++)   help.aux2.y[j] = conj( GTu.y[j] );
    for(int j = GTu.n+1;    j < 3*GTu.n+4;  j++)   help.aux2.y[j] = 0;
    for(int j = 3*GTu.n+4;  j < 4*GTu.n+4;  j++)   help.aux2.y[j] = conj( GTu.y[j] );


    GdGu2mm.convolutionAWT(help.aux1, help.aux2, -1, -1, help.aux3, help.aux4);  // f[x] * g[-a+x]
    //GdGu2mm.KrammersKronig(GdGu2mm, K3, aux1, aux2, aux3, aux4);

    if(out == 1)    cout << ", ";
}

void calculateGGdu_SSN(AWT & GuGd1pp, AWT & GuGd1pm, AWT & GuGd2pm, AWT & GuGd2mm, AWT & GTu, AWT & GTd, auxiliary & help, int out)
{
    if(out == 1)    cout << "GG";

    // this function is not checked for correct up down symmetries!!!
    complex<double> unity(0,1);

    // BUBBLE GuGd1pp
    for(int j = 0;          j < GTd.n+1;    j++)   help.aux1.y[j] = help.FD.y[j] * GTd.y[j] / ( -2.0*Pi*unity );
    for(int j = GTd.n+1;    j < 3*GTd.n+4;  j++)   help.aux1.y[j] = 0;
    for(int j = 3*GTd.n+4;  j < 4*GTd.n+4;  j++)   help.aux1.y[j] = help.FD.y[j] * GTd.y[j] / ( -2.0*Pi*unity );


    for(int j = 0;          j < GTd.n+1;    j++)   help.aux2.y[j] = GTu.y[j];
    for(int j = GTd.n+1;    j < 3*GTd.n+4;  j++)   help.aux2.y[j] = 0;
    for(int j = 3*GTd.n+4;  j < 4*GTd.n+4;  j++)   help.aux2.y[j] = GTu.y[j];

    GuGd1pp.convolutionAWT(help.aux1, help.aux2, 1, -1, help.aux3, help.aux4);  // f[x] * g[x+a]
    //GuGd1pp.KrammersKronig(GuGd1pp, K3, aux1, aux2, aux3, aux4);

    // BUBBLE GuGd1pm
    for(int j = 0;          j < GTd.n+1;    j++)   help.aux1.y[j] = help.FD.y[j] * conj( GTd.y[j] ) / ( -2.0*Pi*unity );
    for(int j = GTd.n+1;    j < 3*GTd.n+4;  j++)   help.aux1.y[j] = 0;
    for(int j = 3*GTd.n+4;  j < 4*GTd.n+4;  j++)   help.aux1.y[j] = help.FD.y[j] * conj( GTd.y[j] ) / ( -2.0*Pi*unity );

    for(int j = 0;          j < GTd.n+1;    j++)   help.aux2.y[j] = GTu.y[j];
    for(int j = GTd.n+1;    j < 3*GTd.n+4;  j++)   help.aux2.y[j] = 0;
    for(int j = 3*GTd.n+4;  j < 4*GTd.n+4;  j++)   help.aux2.y[j] = GTu.y[j];

    GuGd1pm.convolutionAWT(help.aux1, help.aux2, 1, -1, help.aux3, help.aux4);  // f[x] * g[x+a]
    //GuGd1pm.KrammersKronig(GuGd1pm, K3, aux1, aux2, aux3, aux4);


    // BUBBLE GuGd2pm
    for(int j = 0;          j < GTd.n+1;    j++)   help.aux1.y[j] = help.FD.y[j] * GTu.y[j]  / ( -2.0*Pi*unity );
    for(int j = GTd.n+1;    j < 3*GTd.n+4;  j++)   help.aux1.y[j] = 0;
    for(int j = 3*GTd.n+4;  j < 4*GTd.n+4;  j++)   help.aux1.y[j] = help.FD.y[j] * GTu.y[j]  / ( -2.0*Pi*unity );

    for(int j = 0;          j < GTd.n+1;    j++)   help.aux2.y[j] = conj( GTd.y[j] );
    for(int j = GTd.n+1;    j < 3*GTd.n+4;  j++)   help.aux2.y[j] = 0;
    for(int j = 3*GTd.n+4;  j < 4*GTd.n+4;  j++)   help.aux2.y[j] = conj( GTd.y[j] );


    GuGd2pm.convolutionAWT(help.aux1, help.aux2, -1, -1, help.aux3, help.aux4);  // f[x] * g[-a+x]
    //GuGd2pm.KrammersKronig(GuGd2pm, K3, aux1, aux2, aux3, aux4);


    // BUBBLE GuGd2mm
    for(int j = 0;          j < GTd.n+1;    j++)   help.aux1.y[j] = help.FD.y[j] * conj( GTu.y[j] )  / ( -2.0*Pi*unity );
    for(int j = GTd.n+1;    j < 3*GTd.n+4;  j++)   help.aux1.y[j] = 0;
    for(int j = 3*GTd.n+4;  j < 4*GTd.n+4;  j++)   help.aux1.y[j] = help.FD.y[j] * conj( GTu.y[j] )  / ( -2.0*Pi*unity );

    for(int j = 0;          j < GTd.n+1;    j++)   help.aux2.y[j] = conj( GTd.y[j] );
    for(int j = GTd.n+1;    j < 3*GTd.n+4;  j++)   help.aux2.y[j] = 0;
    for(int j = 3*GTd.n+4;  j < 4*GTd.n+4;  j++)   help.aux2.y[j] = conj( GTd.y[j] );


    GuGd2mm.convolutionAWT(help.aux1, help.aux2, -1, -1, help.aux3, help.aux4);  // f[x] * g[-a+x]
    //GuGd2mm.KrammersKronig(GuGd2mm, K3, aux1, aux2, aux3, aux4);


    if(out == 1)    cout << ", ";
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// BOTH ud and du variants for mag field
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void calculate_Dud_static_SSN(AWT & Dud, AWT & GdGu1pp, AWT & GdGu1pm, AWT & GdGu2pm, AWT & GdGu2mm, matrix & Lud, int out)
{
    if(out == 1)    cout << "Dud";

    // positive Dud frequencies
    for(int j = 0;          j < Dud.n+1;    j++)
    {
        Dud.y[j] = 1.0 + GdGu1pp.y[j] * Lud.pp
                       - GdGu1pm.y[j] * Lud.mp
                       + GdGu2pm.y[j] * Lud.mp
                       - GdGu2mm.y[j] * conj(Lud.pp)
                       - GdGu1pp.y[j] * GdGu2mm.y[j] * ( Lud.mp * conj(Lud.mp) - Lud.pp * conj(Lud.pp)  );
    }
    // zero padding of Dud
    for(int j = Dud.n+1;    j < 3*Dud.n+4;  j++)
    {
        Dud.y[j] = 0;
    }
    // negative Dud frequencies
    for(int j = 3*Dud.n+4;  j < 4*Dud.n+4;  j++)
    {
        Dud.y[j] = 1.0 + GdGu1pp.y[j] * Lud.pp
                       - GdGu1pm.y[j] * Lud.mp
                       + GdGu2pm.y[j] * Lud.mp
                       - GdGu2mm.y[j] * conj(Lud.pp)
                       - GdGu1pp.y[j] * GdGu2mm.y[j] * ( Lud.mp * conj(Lud.mp) - Lud.pp * conj(Lud.pp)  );
    }

    if(out == 1)    cout << ", ";
}

void calculate_Ddu_static_SSN(AWT & Ddu, AWT & GuGd1pp, AWT & GuGd1pm, AWT & GuGd2pm, AWT & GuGd2mm, matrix & Ldu, int out)
{
     // this function is not checked for correct up down symmetries!!!
     if(out == 1)    cout << "Ddu";

    // positive Dud frequencies
    for(int j = 0;          j < Ddu.n+1;    j++)
    {
        Ddu.y[j] = 1.0 + GuGd1pp.y[j] * Ldu.pp
                       - GuGd1pm.y[j] * Ldu.mp
                       + GuGd2pm.y[j] * Ldu.mp
                       - GuGd2mm.y[j] * conj(Ldu.pp)
                       - GuGd1pp.y[j] * GuGd2mm.y[j] * ( Ldu.mp * conj(Ldu.mp) - Ldu.pp * conj(Ldu.pp)  );
    }
    // zero padding of Dud
    for(int j = Ddu.n+1;    j < 3*Ddu.n+4;  j++)
    {
        Ddu.y[j] = 0;
    }
    // negative Dud frequencies
    for(int j = 3*Ddu.n+4;  j < 4*Ddu.n+4;  j++)
    {
        Ddu.y[j] = 1.0 + GuGd1pp.y[j] * Ldu.pp
                       - GuGd1pm.y[j] * Ldu.mp
                       + GuGd2pm.y[j] * Ldu.mp
                       - GuGd2mm.y[j] * conj(Ldu.pp)
                       - GuGd1pp.y[j] * GuGd2mm.y[j] * ( Ldu.mp * conj(Ldu.mp) - Ldu.pp * conj(Ldu.pp)  );
    }

    if(out == 1)    cout << ", ";
}

void calculate_Ddu_sym_SSN(AWT & Ddu, AWT & Dud, int out)
{
    if(out == 1)    cout << "Ddu";

                                                    Ddu.y[0] = Dud.y[0];
    for(int j = 1;          j < Ddu.n+1;    j++)    Ddu.y[j] = conj(   Dud.y[Dud.nn - j]   );
    for(int j = Ddu.n+1;    j < 3*Ddu.n+4;  j++)    Ddu.y[j] = 0;
    for(int j = 3*Ddu.n+4;  j < 4*Ddu.n+4;  j++)    Ddu.y[j] = conj(   Dud.y[Dud.nn - j]   );

}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// BOTH ud and du variants for mag field
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void calculate_Kud_static_SSN(matrix & Kud, AWT & GdGu1pp, AWT & GdGu1pm, AWT & GdGu2pm, AWT & GdGu2mm, matrix & Lud, int out)
{
    if(out == 1)    cout << "Kud";

    Kud.pp = - Lud.pp * Lud.pp * GdGu1pp.y[0] - Lud.mp * conj(Lud.mp) * conj(GdGu1pp.y[0])
             - Lud.pp * ( Lud.pp * conj(Lud.pp) - Lud.mp * conj(Lud.mp)  ) * GdGu1pp.y[0] * conj(GdGu1pp.y[0]);

    Kud.mp = - Lud.mp * Lud.pp * GdGu1pp.y[0] - Lud.mp * conj(Lud.pp) * conj(GdGu1pp.y[0])
             - Lud.mp * ( Lud.pp * conj(Lud.pp) - Lud.mp * conj(Lud.mp)  ) * GdGu1pp.y[0] * conj(GdGu1pp.y[0]);

    Kud.pm = conj(Kud.mp);
    Kud.mm = conj(Kud.pp);

    if(out == 1)    cout << ", ";
}

void calculate_Kdu_static_SSN(matrix & Kdu, AWT & GuGd1pp, AWT & GuGd1pm, AWT & GuGd2pm, AWT & GuGd2mm, matrix & Ldu, int out)
{
    // this function is not checked for correct up down symmetries!!!

    if(out == 1)    cout << "Kdu";

    Kdu.pp = - Ldu.pp * Ldu.pp * GuGd1pp.y[0] - Ldu.mp * conj(Ldu.mp) * conj(GuGd1pp.y[0])
             - Ldu.pp * ( Ldu.pp * conj(Ldu.pp) - Ldu.mp * conj(Ldu.mp)  ) * GuGd1pp.y[0] * conj(GuGd1pp.y[0]);

    Kdu.mp = - Ldu.mp * Ldu.pp * GuGd1pp.y[0] - Ldu.mp * conj(Ldu.pp) * conj(GuGd1pp.y[0])
             - Ldu.mp * ( Ldu.pp * conj(Ldu.pp) - Ldu.mp * conj(Ldu.mp)  ) * GuGd1pp.y[0] * conj(GuGd1pp.y[0]);

    Kdu.pm = conj(Kdu.mp);    Kdu.mm = conj(Kdu.pp);

    if(out == 1)    cout << ", ";
}

void calculate_Kdu_sym_SSN(matrix & Kdu, matrix & Kud, int out)
{
    // this function is not checked for correct up down symmetries!!!

    if(out == 1)    cout << "Kdu";

    Kdu.pp = Kud.pp;
    Kdu.mp = Kud.pm;
    Kdu.pm = Kud.mp;
    Kdu.mm = Kud.mm;

    if(out == 1)    cout << ", ";
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// BOTH ud and du variants for mag field
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void calculate_FDud_SSN(matrix & FDud, AWT & GTd, AWT & GTu, AWT & Dud, auxiliary & help, int out)
{

    if(out == 1)    cout << "FDud";
    int nn = help.aux1.nn;

    // FDud(+,+)
    // zero aux frequency
        help.aux1.y[0] = help.FD.y[0] * GTu.y[0] * (-1.0) * imag( GTd.y[0] ) / conj( Dud.y[0] )
                       + help.FD.y[0] * GTd.y[0] * (-1.0) * imag( GTu.y[0] ) /       Dud.y[0] ;
    // positive aux frequencies
    for(int j = 1;          j < GTu.n+1;    j++)
    {
        help.aux1.y[j] = help.FD.y[j] * GTu.y[ j  ] * (-1.0) * imag( GTd.y[nn-j] ) / conj( Dud.y[nn - j] )
                       + help.FD.y[j] * GTd.y[nn-j] * (+1.0) * imag( GTu.y[  j ] ) /     ( Dud.y[nn - j] );
    }
    // zero padding of aux
    for(int j = GTu.n+1;    j < 3*GTu.n+4;  j++)
    {
        help.aux1.y[j] = 0;
    }
    // negative aux frequencies
    for(int j = 3*GTu.n+4;  j < 4*GTu.n+4;  j++)
    {
        help.aux1.y[j] = help.FD.y[j] * GTu.y[ j  ] * (-1.0) * imag( GTd.y[nn-j] ) / conj( Dud.y[nn-j] )
                       + help.FD.y[j] * GTd.y[nn-j] * (+1.0) * imag( GTu.y[  j ] ) /     ( Dud.y[nn-j] );
    }

    FDud.pp = Simpson(help.aux1) / ( - Pi);

    // FDud(-,+)                                                    help.aux1.y[0] = help.FD.y[0] * imag( GTu.y[0] * conj( GTd.y[ 0  ] ) ) / Dud.y[ 0  ];
    for(int j = 1;          j < GTu.n+1;    j++)    help.aux1.y[j] = help.FD.y[j] * imag( GTu.y[j] * conj( GTd.y[nn-j] ) ) / Dud.y[nn-j];
    for(int j = GTu.n+1;    j < 3*GTu.n+4;  j++)    help.aux1.y[j] = 0;
    for(int j = 3*GTu.n+4;  j < 4*GTu.n+4;  j++)    help.aux1.y[j] = help.FD.y[j] * imag( GTu.y[j] * conj( GTd.y[nn-j] ) ) / Dud.y[nn-j];

    FDud.mp = Simpson(help.aux1) / ( - Pi);


    FDud.pm = conj(FDud.mp);    // FDud(-,+)
    FDud.mm = conj(FDud.pp);    // FDud(-,+)

    if(out == 1)    cout << ", ";
}

void calculate_FDdu_SSN(matrix & FDdu, AWT & GTd, AWT & GTu, AWT & Ddu, auxiliary & help, int out)
{
    // this function is not checked for correct up down symmetries!!!

    if(out == 1)    cout << "FDdu";
    int nn = help.aux1.nn;

    // FDud(+,+)
    // zero aux frequency
        help.aux1.y[0] = help.FD.y[0] * GTd.y[0] * (-1.0) * imag( GTu.y[0] ) / conj( Ddu.y[0] )
                       + help.FD.y[0] * GTu.y[0] * (-1.0) * imag( GTd.y[0] ) /       Ddu.y[0] ;
    // positive aux frequencies
    for(int j = 1;          j < GTu.n+1;    j++)
    {
        help.aux1.y[j] = help.FD.y[j] * GTd.y[ j  ] * (-1.0) * imag( GTu.y[nn-j] ) / conj( Ddu.y[nn-j] )
                       + help.FD.y[j] * GTu.y[nn-j] * (+1.0) * imag( GTd.y[  j ] ) /     ( Ddu.y[nn-j] );
    }
    // zero padding of aux
    for(int j = GTu.n+1;    j < 3*GTu.n+4;  j++)
    {
        help.aux1.y[j] = 0;
    }
    // negative aux frequencies
    for(int j = 3*GTu.n+4;  j < 4*GTu.n+4;  j++)
    {
        help.aux1.y[j] = help.FD.y[j] * GTd.y[ j  ] * (-1.0) * imag( GTu.y[nn-j] ) / conj( Ddu.y[nn-j] )
                       + help.FD.y[j] * GTu.y[nn-j] * (+1.0) * imag( GTd.y[ j  ] ) /     ( Ddu.y[nn-j] );
    }

    FDdu.pp = Simpson(help.aux1) / ( - Pi);

    // FDud(-,+)                                                    help.aux1.y[0] = help.FD.y[0] * imag( GTd.y[0] * conj( GTu.y[     0     ] ) ) / Ddu.y[    0     ];
    for(int j = 1;          j < GTd.n+1;    j++)    help.aux1.y[j] = help.FD.y[j] * imag( GTd.y[j] * conj( GTu.y[nn - j] ) ) / Ddu.y[nn-j];
    for(int j = GTu.n+1;    j < 3*GTu.n+4;  j++)    help.aux1.y[j] = 0;
    for(int j = 3*GTu.n+4;  j < 4*GTu.n+4;  j++)    help.aux1.y[j] = help.FD.y[j] * imag( GTd.y[j] * conj( GTu.y[nn - j] ) ) / Ddu.y[nn-j];

    FDdu.mp = Simpson(help.aux1) / ( - Pi);

    // FDud(-,+)
    FDdu.pm = conj(FDdu.mp);

    // FDud(-,+)
    FDdu.mm = conj(FDdu.pp);

    if(out == 1)    cout << ", ";
}

void calculate_FDdu_sym_SSN(matrix & FDdu, matrix & FDud, int out)
{
    if(out == 1)    cout << "FDdu";

    FDdu.pp = FDud.pp;
    FDdu.mp = FDud.pm;
    FDdu.pm = FDud.mp;
    FDdu.mm = FDud.mm;

    if(out == 1)    cout << ", ";
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// BOTH ud and du variants for mag field
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void calculate_BDud_SSN(matrix & BDud, AWT & GTd, AWT & GTu, AWT & Dud, auxiliary & help, int out)
{
    if(out == 1)    cout << "BDud";
    int nn = help.aux1.nn;

    // BDud(+,+)
    for(int j = 1;          j < GTu.n+1;    j++)    help.aux1.y[j] = help.BE.y[j] * GTu.y[j] * GTd.y[nn - j] * imag( 1.0 / conj( Dud.y[nn - j] )  );
    for(int j = GTu.n+1;    j < 3*GTu.n+4;  j++)    help.aux1.y[j] = 0;
    for(int j = 3*GTu.n+4;  j < 4*GTu.n+4;  j++)    help.aux1.y[j] = help.BE.y[j] * GTu.y[j] * GTd.y[nn - j] * imag( 1.0 / conj( Dud.y[nn - j] )  );
    // zero aux frequency to avoid discontinuity trap
                                                    help.aux1.y[0] = 0.5 * (help.aux1.y[4*GTu.n+3] + help.aux1.y[1] );

    BDud.pp = Simpson(help.aux1) / ( Pi);


    // BDud(-,+)

    for(int j = 1;          j < GTu.n+1;    j++)    help.aux1.y[j] = help.BE.y[j] * GTu.y[j] * conj( GTd.y[nn-j] ) * imag( 1.0 / conj( Dud.y[nn-j] )  );
    for(int j = GTu.n+1;    j < 3*GTu.n+4;  j++)    help.aux1.y[j] = 0;
    for(int j = 3*GTu.n+4;  j < 4*GTu.n+4;  j++)    help.aux1.y[j] = help.BE.y[j] * GTu.y[j] * conj( GTd.y[nn-j] ) * imag( 1.0 / conj( Dud.y[nn-j] )  );
    // zero aux frequency to avoid discontinuity trap
                                                    help.aux1.y[0] = 0.5 * (help.aux1.y[4*GTu.n+3] + help.aux1.y[1] );


    BDud.mp = Simpson(help.aux1) / ( Pi);

    // BDud(-,+)
    BDud.pm = conj(BDud.mp);

    // BDud(-,+)
    BDud.mm = conj(BDud.pp);
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


    if(out == 1)    cout << ", ";
}


void calculate_BDdu_SSN(matrix & BDdu, AWT & GTd, AWT & GTu, AWT & Ddu, auxiliary & help, int out)
{
    // this function is not checked for correct up down symmetries!!!
    if(out == 1)    cout << "BDud";
    int nn = help.aux1.nn;

    // BDud(+,+)
    for(int j = 1;          j < GTu.n+1;    j++)    help.aux1.y[j] = help.BE.y[j] * GTd.y[j] * GTu.y[nn-j] * imag( 1.0 / conj( Ddu.y[nn-j] )  );
    for(int j = GTu.n+1;    j < 3*GTu.n+4;  j++)    help.aux1.y[j] = 0;
    for(int j = 3*GTu.n+4;  j < 4*GTu.n+4;  j++)    help.aux1.y[j] = help.BE.y[j] * GTd.y[j] * GTu.y[nn-j] * imag( 1.0 / conj( Ddu.y[nn-j] )  );
    // zero aux frequency to avoid discontinuity trap
                                                    help.aux1.y[0] = 0.5 * (help.aux1.y[4*GTu.n+3] + help.aux1.y[1] );


    BDdu.pp = Simpson(help.aux1) / ( Pi);

    // BDud(-,+)
    for(int j = 1;          j < GTu.n+1;    j++)    help.aux1.y[j] = help.BE.y[j] * GTd.y[j] * conj( GTu.y[nn-j] ) * imag( 1.0 / conj( Ddu.y[nn-j] )  );
    for(int j = GTu.n+1;    j < 3*GTu.n+4;  j++)    help.aux1.y[j] = 0;
    for(int j = 3*GTu.n+4;  j < 4*GTu.n+4;  j++)    help.aux1.y[j] = help.BE.y[j] * GTd.y[j] * conj( GTu.y[nn-j] ) * imag( 1.0 / conj( Ddu.y[nn-j] )  );
    // zero aux frequency
                                                    help.aux1.y[0] = 0.5 * (help.aux1.y[4*GTu.n+3] + help.aux1.y[1] );


    BDdu.mp = Simpson(help.aux1) / ( Pi);
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // BDud(-,+)
    BDdu.pm = conj(BDdu.mp);

    // BDud(-,+)
    BDdu.mm = conj(BDdu.pp);
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


    if(out == 1)    cout << ", ";
}

void calculate_BDdu_sym_SSN(matrix & BDdu, matrix & BDud, int out)
{

    if(out == 1)    cout << "BDdu";

    BDdu.pp = BDud.pp;
    BDdu.mp = BDud.pm;
    BDdu.pm = BDud.mp;
    BDdu.mm = BDud.mm;

    if(out == 1)    cout << ", ";
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// BOTH ud and du variants for mag field
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void calculate_Lud_SSN(matrix & Lud, matrix & Kud, matrix & FDud, matrix & BDud, double U, int out)
{
    if(out == 1)    cout << "Lud";

    // xx, attention real is a cheater stuff here!
    Lud.pp = U / ( 1.0 + Kud.pp * real( BDud.pp + FDud.pp)  );
    Lud.mp = U / ( 1.0 + Kud.pm * real( BDud.pm + FDud.pm)  );

    Lud.pm = conj(Lud.mp);
    Lud.mm = conj(Lud.pp);

    if(out == 1)    cout << ", ";
}



void calculate_Ldu_SSN(matrix & Ldu, matrix & Kdu, matrix & FDdu, matrix & BDdu, double U, int out)
{
    if(out == 1)    cout << "Lud";

    // xx, attention real is a cheater stuff here!
    Ldu.pp = U / ( 1.0 + Kdu.pp * real( BDdu.pp + FDdu.pp)  );
    Ldu.mp = U / ( 1.0 + Kdu.pm * real( BDdu.pm + FDdu.pm)  );

    Ldu.pm = conj(Ldu.mp);
    Ldu.mm = conj(Ldu.pp);

    if(out == 1)    cout << ", ";
}


void calculate_Ldu_sym_SSN(matrix & Ldu, matrix & Lud, int out)
{
    if(out == 1)    cout << "Ldu";

    Ldu.pp = Lud.pp;
    Ldu.mp = Lud.pm;
    Ldu.pm = Lud.mp;
    Ldu.mm = Lud.mm;

    if(out == 1)    cout << ", ";
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// ???? for mag field
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void calculate_fGG_SSN(complex<double> & fGG, AWT & GTd, AWT & GTu, auxiliary & help)
{
    complex<double> unity(0,1);

    // positive aux frequencies
    for(int j = 0;          j < GTu.n+1;    j++)    help.aux1.y[j] = help.FD.y[j] * GTd.y[j] * GTu.y[j] ;
    // zero padding of aux
    for(int j = GTu.n+1;    j < 3*GTu.n+4;  j++)    help.aux1.y[j] = 0;
    // negative aux frequencies
    for(int j = 3*GTu.n+4;  j < 4*GTu.n+4;  j++)    help.aux1.y[j] = help.FD.y[j] * GTd.y[j] * GTu.y[j] ;

    fGG = Simpson(help.aux1) / ( - 2.0 * Pi * unity);

}


void set_lambda_iter_SSN(string mode, string dirname, double U, info in,
                     matrix & Lud,  matrix & Kud,  matrix & BDud,  matrix & FDud)
{
    complex<double> unity(0,1);

    // initial Lud set to just the bare U
    if(mode == "bare_U")    set_matrix_real(   Lud, in.U_min);

    // initial Lud set to just the bare U
    if(mode == "bare_kT")
    {
        Lud.pp = U/2.0;
        Lud.mp = U;
        Lud.pm = conj(Lud.mp);
        Lud.mm = conj(Lud.pp);
    }

    // initial Lud set to just the bare U
    if(mode == "import")
    {
        ifstream ListFile;
        ListFile.open (dirname.c_str());   // list contains all the file names to be processed
        string item;                    // used for temporary store the data

        if (ListFile.fail()) {cerr<< "Error opening the import file final" << endl; exit(1);}

        for(int i=0; i<58; i++)
        {
            ListFile >> item;
            cout << i << " " << item << endl;

            if(i == 3)   Lud.pp       = boost::lexical_cast<double>(item);
            if(i == 5)   Lud.pp       = Lud.pp + unity * boost::lexical_cast<double>(item);
            if(i == 3)   Lud.mp       = boost::lexical_cast<double>(item);
            if(i == 5)   Lud.mp       = Lud.mp + unity * boost::lexical_cast<double>(item);

            if(i == 3)   Kud.pp       = boost::lexical_cast<double>(item);
            if(i == 5)   Kud.pp       = Kud.pp + unity * boost::lexical_cast<double>(item);
            if(i == 3)   Kud.mp       = boost::lexical_cast<double>(item);
            if(i == 5)   Kud.mp       = Kud.mp + unity * boost::lexical_cast<double>(item);

            if(i == 3)   BDud.pp      = boost::lexical_cast<double>(item);
            if(i == 5)   BDud.pp      = BDud.pp + unity * boost::lexical_cast<double>(item);
            if(i == 3)   BDud.mp      = boost::lexical_cast<double>(item);
            if(i == 5)   BDud.mp      = BDud.mp + unity * boost::lexical_cast<double>(item);


            if(i == 3)   FDud.pp      = boost::lexical_cast<double>(item);
            if(i == 5)   FDud.pp      = FDud.pp + unity * boost::lexical_cast<double>(item);
            if(i == 3)   FDud.mp      = boost::lexical_cast<double>(item);
            if(i == 5)   FDud.mp      = FDud.mp + unity * boost::lexical_cast<double>(item);

        }

        ListFile.close();
    }

}

void lambda_iter_output_SSN(ofstream & log_output, int iter_lambda, int iter_sigma, complex<double> & aSTu, complex<double> & aSTd,
                        matrix Lud, matrix Ldu, matrix Kud, matrix Kdu, matrix BDud, matrix BDdu, matrix FDud, matrix FDdu)
{
    log_output << "lambda iteration: " << iter_lambda << ", sigma iteration: " << iter_sigma;
    log_output << endl << "          ";
    log_output << " Lud: " <<  real(Lud.pp) << " " <<  imag(Lud.pp) << " " <<  real(Lud.pm) << " " <<  imag(Lud.pm) << endl << "          ";
    log_output << " Kud: " <<  real(Kud.pp) << " " <<  imag(Kud.pp) << " " <<  real(Kud.pm) << " " <<  imag(Kud.pm) << endl << "          ";
    log_output << "BDud: " << real(BDud.pp) << " " << imag(BDud.pp) << " " << real(BDud.pm) << " " << imag(BDud.pm) << endl << "          ";
    log_output << "FDud: " << real(FDud.pp) << " " << imag(FDud.pp) << " " << real(FDud.pm) << " " << imag(FDud.pm) << endl << "          ";
    log_output << endl;

    log_output << "          ";
    log_output << " Ldu: " <<  real(Ldu.pp) << " " <<  imag(Ldu.pp) << " " <<  real(Ldu.pm) << " " <<  imag(Ldu.pm) << endl << "          ";
    log_output << " Kdu: " <<  real(Kdu.pp) << " " <<  imag(Kdu.pp) << " " <<  real(Kdu.pm) << " " <<  imag(Kdu.pm) << endl << "          ";
    log_output << "BDdu: " << real(BDdu.pp) << " " << imag(BDdu.pp) << " " << real(BDdu.pm) << " " << imag(BDdu.pm) << endl << "          ";
    log_output << "FDdu: " << real(FDdu.pp) << " " << imag(FDdu.pp) << " " << real(FDdu.pm) << " " << imag(FDdu.pm) << endl << "          ";
    log_output << endl << "          " << "aSTu = " << aSTu << "aSTd = " << aSTd << endl;
}

void lambda_iter_console_SSN(int out, matrix Lud, matrix Ldu, matrix Kud, matrix Kdu, matrix BDud, matrix BDdu, matrix FDud, matrix FDdu)
{
     // cout ud versions
     if(out == 1)    cout << endl << "          ";
     cout << " Lud: " << Lud.pp << " " << Lud.pm << endl << "          ";
     if(out == 1)    cout << " Kud: " <<  Kud.pp << " " <<  Kud.pm << endl << "          ";
     if(out == 1)    cout << "BDud: " << BDud.pp << " " << BDud.pm << endl << "          ";
     if(out == 1)    cout << "FDud: " << FDud.pp << " " << FDud.pm << endl << "          ";
     if(out == 1)    cout << endl;

     // cout du versions
     if(out == 1)    cout << endl << "          ";
     cout << " Ldu: " << Ldu.pp << " " << Ldu.pm << endl << "          ";
     if(out == 1)    cout << " Kdu: " <<  Kdu.pp << " " <<  Kdu.pm << endl << "          ";
     if(out == 1)    cout << "BDdu: " << BDdu.pp << " " << BDdu.pm << endl << "          ";
     if(out == 1)    cout << "FDdu: " << FDdu.pp << " " << FDdu.pm << endl << "          ";
     if(out == 1)    cout << endl;
}

void lambda_log_output_SSN(ofstream & log_output, int iter_total, matrix Lud, matrix Ldu, matrix Kud, matrix Kdu, matrix BDud, matrix BDdu, matrix FDud, matrix FDdu)
{
            log_output << iter_total
            << " " << real( Lud.pp) << " " << imag( Lud.pp) << " " << real( Lud.pm) << " " << imag( Lud.pm)
            << " " << real( Ldu.pp) << " " << imag( Ldu.pp) << " " << real( Ldu.pm) << " " << imag( Ldu.pm)

            << " " << real( Kud.pp) << " " << imag( Kud.pp) << " " << real( Kud.pm) << " " << imag( Kud.pm)
            << " " << real( Kdu.pp) << " " << imag( Kdu.pp) << " " << real( Kdu.pm) << " " << imag( Kdu.pm)

            << " " << real(BDud.pp) << " " << imag(BDud.pp) << " " << real(BDud.pm) << " " << imag(BDud.pm)
            << " " << real(BDdu.pp) << " " << imag(BDdu.pp) << " " << real(BDdu.pm) << " " << imag(BDdu.pm)

            << " " << real(FDud.pp) << " " << imag(FDud.pp) << " " << real(FDud.pm) << " " << imag(FDud.pm)
            << " " << real(FDdu.pp) << " " << imag(FDdu.pp) << " " << real(FDdu.pm) << " " << imag(FDdu.pm) << endl;
}

void calculate_aST_SSN(complex<double> & aSTu, complex<double> & aSTd,
                    matrix & Lud, matrix & Ldu,
                    AWT & GTd, AWT & GTu, auxiliary & help, int out)
{


    complex<double> unity(0,1);

    complex<double> Gup, Gdp, Gum, Gdm;

    // GT up plus integral
    for(int j = 0;          j < GTu.n+1;    j++)   help.aux1.y[j] = help.FD.y[j] * GTu.y[j];
    for(int j = GTu.n+1;    j < 3*GTu.n+4;  j++)   help.aux1.y[j] = 0;
    for(int j = 3*GTu.n+4;  j < 4*GTu.n+4;  j++)   help.aux1.y[j] = help.FD.y[j] * GTu.y[j];

    Gup = Simpson(help.aux1) / ( -2.0 * unity * Pi );

    // GT down plus integral
    for(int j = 0;          j < GTu.n+1;    j++)   help.aux1.y[j] = help.FD.y[j] * GTd.y[j];
    for(int j = GTu.n+1;    j < 3*GTu.n+4;  j++)   help.aux1.y[j] = 0;
    for(int j = 3*GTu.n+4;  j < 4*GTu.n+4;  j++)   help.aux1.y[j] = help.FD.y[j] * GTd.y[j];

    Gdp = Simpson(help.aux1) / ( -2.0 * unity * Pi );

    // GT up minus integral
    for(int j = 0;          j < GTu.n+1;    j++)   help.aux1.y[j] = help.FD.y[j] * conj( GTu.y[j] );
    for(int j = GTu.n+1;    j < 3*GTu.n+4;  j++)   help.aux1.y[j] = 0;
    for(int j = 3*GTu.n+4;  j < 4*GTu.n+4;  j++)   help.aux1.y[j] = help.FD.y[j] * conj( GTu.y[j] );

    Gum = Simpson(help.aux1) / ( -2.0 * unity * Pi );

    // GT down minus integral
    for(int j = 0;          j < GTu.n+1;    j++)   help.aux1.y[j] = help.FD.y[j] * conj( GTd.y[j] );
    for(int j = GTu.n+1;    j < 3*GTu.n+4;  j++)   help.aux1.y[j] = 0;
    for(int j = 3*GTu.n+4;  j < 4*GTu.n+4;  j++)   help.aux1.y[j] = help.FD.y[j] * conj( GTd.y[j] );

    Gdm = Simpson(help.aux1) / ( -2.0 * unity * Pi );


    aSTu = - 0.25 * (Lud.pp + Ldu.pp) * (Gup - Gdp)
           - 0.25 * (Lud.pm + Ldu.pm) * (Gup - Gdp);
    aSTd = - 0.25 * (Lud.pp + Ldu.pp) * (Gdm - Gum)
           - 0.25 * (Lud.pm + Ldu.pm) * (Gdm - Gum);

    aSTu = ( aSTu + aSTd ) / 2.0;
    aSTd = ( aSTu + aSTd ) / 2.0;

    cout << "anomalous = " << aSTu << " " << aSTd << endl;

}

void calculate_nT_mT_SSN(double & nT, double & mT, double mu, double U, AWT & GTu, AWT & GTd, auxiliary & help)
{
    for(int i=0;          i<GTu.n+1;   i++ )       help.aux1.y[i] = help.FD.y[i] * imag( GTu.y[i] );
    for(int i=GTu.n+1;    i<3*GTu.n+3; i++ )       help.aux1.y[i] = 0;
    for(int i=3*GTu.n+4;  i<4*GTu.n+4; i++ )       help.aux1.y[i] = help.FD.y[i] * imag( GTu.y[i] );



    double nTup = - Simpson_Re(help.aux1) / Pi;


    for(int i=0;          i<GTu.n+1;   i++ )       help.aux1.y[i] = help.FD.y[i] * imag( GTd.y[i] );
    for(int i=GTu.n+1;    i<3*GTu.n+3; i++ )       help.aux1.y[i] = 0;
    for(int i=3*GTu.n+4;  i<4*GTu.n+4; i++ )       help.aux1.y[i] = help.FD.y[i] * imag( GTd.y[i] );

    double nTdown = - Simpson_Re(help.aux1) / Pi;

    nT = nTup + nTdown;
    mT = nTup - nTdown;

    // add tail correction
    nT = nT - GTd.xMax * imag( GTd.y[3*GTu.n+4] ) / (Pi);
    cout << "nT up: "  << nTup << "nT down: "  << nTdown << " nT corr: "
         << - GTd.xMax * imag( GTd.y[3*GTu.n+4] ) / (Pi)  << endl;

    // if mu = 0 manually reestablish nT
    if( abs(mu) < 0.00001)
    {
        cout << "nT set to 1.0" << endl;
        nT = 1;
    }
}

void calculate_nS_mS_SSN(double & nS, double & mS, AWT & GSu, AWT & GSd, auxiliary & help)
{
    for(int i=0;          i<GSu.n+1;   i++ )       help.aux1.y[i] = help.FD.y[i] * imag( GSu.y[i] );
    for(int i=GSu.n+1;    i<3*GSu.n+3; i++ )       help.aux1.y[i] = 0;
    for(int i=3*GSu.n+4;  i<4*GSu.n+4; i++ )       help.aux1.y[i] = help.FD.y[i] * imag( GSu.y[i] );

    nS = - Simpson_Re(help.aux1) / Pi;
    mS = - Simpson_Re(help.aux1) / Pi;
    cout << "nS up = " << nS;

    for(int i=0;          i<GSu.n+1;   i++ )       help.aux1.y[i] = help.FD.y[i] * imag( GSd.y[i] );
    for(int i=GSu.n+1;    i<3*GSu.n+3; i++ )       help.aux1.y[i] = 0;
    for(int i=3*GSu.n+4;  i<4*GSu.n+4; i++ )       help.aux1.y[i] = help.FD.y[i] * imag( GSd.y[i] );

    nS = nS - Simpson_Re(help.aux1) / Pi;
    mS = mS + Simpson_Re(help.aux1) / Pi;
    cout << "  nS up + down = " << nS << endl;
    cout << "  mS = " << mS << endl;

    // add tail
    nS = nS - GSu.xMax * imag( GSu.y[3*GSu.n+4] ) / (Pi);
    //cout << "nS corr: " << - GSu.xMax * imag( GSu.y[3*GSu.n+4] ) / (Pi)  << endl;

}


void calculate_dynamic_selfenergy_SSN(AWT & SigmaUp, AWT & SigmaDown, info in, AWT & Dud, AWT & Ddu, matrix & Lud, matrix & Ldu,
                              AWT & GTd, AWT & GTu, string dirname, auxiliary & help)
{
    complex<double> unity(0,1);

    AWT ThetaPP, ThetaPM;
    ThetaPP.initializeAWT(GTu.n, GTu.xMax, in.kT);
    ThetaPM.initializeAWT(GTu.n, GTu.xMax, in.kT);

    double dx = GTu.xMax / GTu.n;
    cout << "   Self-energy Up calculation" << endl;

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // GAMMA MATRIX
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    matrix gamma;

    // positive aux frequencies
    for(int j = 0;          j < GTu.n+1;    j++)    help.aux1.y[j] = help.FD.y[j] * GTu.y[j] * GTu.y[j];
    for(int j = GTu.n+1;    j < 3*GTu.n+4;  j++)    help.aux1.y[j] = 0;
    for(int j = 3*GTu.n+4;  j < 4*GTu.n+4;  j++)    help.aux1.y[j] = help.FD.y[j] * GTu.y[j] * GTu.y[j];

    gamma.pp = Lud.pp + ( Lud.pp * conj( Lud.pp) - Lud.mp * conj( Lud.mp) ) * conj( Simpson(help.aux1) / ( -2.0*Pi*unity));
    gamma.mp = Lud.mp;
    gamma.pm = conj( gamma.mp );
    gamma.mm = conj( gamma.pp );
    cout << "   gamma-s: " << gamma.pp << " " << gamma.mp << " " << gamma.pm << " " << gamma.mm << endl;
    //////////////////////////////////////////////////////////////////////////////////////////////////////V////////


    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //  THETA CALCULATIONS Sigma UP
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // FIRST PART of THETA PLUS
    // first integrand
    for(int j = 0;          j < GTu.n+1;    j++)    help.aux4.y[j] = help.FD.y[j] * imag( GTu.y[j] ) / Pi;
    for(int j = GTu.n+1;    j < 3*GTu.n+4;  j++)    help.aux4.y[j] = 0;
    for(int j = 3*GTu.n+4;  j < 4*GTu.n+4;  j++)    help.aux4.y[j] = help.FD.y[j] * imag( GTu.y[j] ) / Pi;
    // second integrand
    for(int j = 0;          j < GTu.n+1;    j++)    help.aux5.y[j] = GTd.y[j];
    for(int j = GTu.n+1;    j < 3*GTu.n+4;  j++)    help.aux5.y[j] = 0;
    for(int j = 3*GTu.n+4;  j < 4*GTu.n+4;  j++)    help.aux5.y[j] = GTd.y[j];

    help.aux1.convolutionAWT(help.aux4, help.aux5, 1, -1, help.aux6, help.aux7);
    help.aux1.KrammersKronig(help.aux1, help.K3, help.aux4, help.aux5, help.aux6, help.aux7);

    string name;
    if(in.print_mode == 1)
    {
        cout << "BEintUp" << endl;
        name = dirname + "/ThetaPP1up";
        help.aux1.exportAWTasFUN(name, in.display, in.xMax, 1, in.output_mode);
    }


    // SECOND PART of THETA PLUS
    // first integrand
    for(int j = 0;          j < GTu.n+1;    j++)    help.aux4.y[j] = help.FD.y[j] * GTd.y[j] / (2 * Pi * unity );
    for(int j = GTu.n+1;    j < 3*GTu.n+4;  j++)    help.aux4.y[j] = 0;
    for(int j = 3*GTu.n+4;  j < 4*GTu.n+4;  j++)    help.aux4.y[j] = help.FD.y[j] * GTd.y[j] / (2 * Pi * unity );
    // second integrand
    for(int j = 0;          j < GTu.n+1;    j++)    help.aux5.y[j] = conj( GTu.y[j] );
    for(int j = GTu.n+1;    j < 3*GTu.n+4;  j++)    help.aux5.y[j] = 0;
    for(int j = 3*GTu.n+4;  j < 4*GTu.n+4;  j++)    help.aux5.y[j] = conj( GTu.y[j] );

    help.aux2.convolutionAWT(help.aux4, help.aux5, -1, -1, help.aux6, help.aux7);
    help.aux2.KrammersKronig(help.aux2, help.K3, help.aux4, help.aux5, help.aux6, help.aux7);

    if(in.print_mode == 1)
    {
        cout << "BEintUp" << endl;
        name = dirname + "/ThetaPP2up";
        help.aux2.exportAWTasFUN(name, in.display, in.xMax, 1, in.output_mode);
    }

    // THIRD PART of THETA PLUS
    // first integrand
    for(int j = 0;          j < GTu.n+1;    j++)    help.aux4.y[j] = help.FD.y[j] * conj( GTd.y[j] ) / (2 * Pi * unity );
    for(int j = GTu.n+1;    j < 3*GTu.n+4;  j++)    help.aux4.y[j] = 0;
    for(int j = 3*GTu.n+4;  j < 4*GTu.n+4;  j++)    help.aux4.y[j] = help.FD.y[j] * conj( GTd.y[j] ) / (2 * Pi * unity );
    // second integrand
    for(int j = 0;          j < GTu.n+1;    j++)    help.aux5.y[j] = conj( GTu.y[j] );
    for(int j = GTu.n+1;    j < 3*GTu.n+4;  j++)    help.aux5.y[j] = 0;
    for(int j = 3*GTu.n+4;  j < 4*GTu.n+4;  j++)    help.aux5.y[j] = conj( GTu.y[j] );

    help.aux3.convolutionAWT(help.aux4, help.aux5, -1, -1, help.aux6, help.aux7);
    help.aux3.KrammersKronig(help.aux3, help.K3, help.aux4, help.aux5, help.aux6, help.aux7);

    if(in.print_mode == 1)
    {
        cout << "BEintUp" << endl;
        name = dirname + "/ThetaPP3up";
        help.aux3.exportAWTasFUN(name, in.display, in.xMax, 1, in.output_mode);
    }


    // THETA PLUS
    for(int j = 0;          j < GTu.n+1;    j++)    ThetaPP.y[j] = -gamma.pp * help.aux1.y[j] - gamma.pp * help.aux2.y[j] + gamma.mp * help.aux3.y[j];
    for(int j = GTu.n+1;    j < 3*GTu.n+4;  j++)    ThetaPP.y[j] = 0;
    for(int j = 3*GTu.n+4;  j < 4*GTu.n+4;  j++)    ThetaPP.y[j] = -gamma.pp * help.aux1.y[j] - gamma.pp * help.aux2.y[j] + gamma.mp * help.aux3.y[j];


    if(in.print_mode == 1)
    {
        cout << "BEintUp" << endl;
        name = dirname + "/ThetaPPup";
        ThetaPP.exportAWTasFUN(name, in.display, in.xMax, 1, in.output_mode);
    }

    // FIRST PART of THETA MINUS
    // first integrand
    for(int j = 0;          j < GTu.n+1;    j++)    help.aux4.y[j] = help.FD.y[j] * imag( GTu.y[j] ) / Pi;
    for(int j = GTu.n+1;    j < 3*GTu.n+4;  j++)    help.aux4.y[j] = 0;
    for(int j = 3*GTu.n+4;  j < 4*GTu.n+4;  j++)    help.aux4.y[j] = help.FD.y[j] * imag( GTu.y[j] ) / Pi;
    // second integrand
    for(int j = 0;          j < GTu.n+1;    j++)    help.aux5.y[j] = conj( GTd.y[j] );
    for(int j = GTu.n+1;    j < 3*GTu.n+4;  j++)    help.aux5.y[j] = 0;
    for(int j = 3*GTu.n+4;  j < 4*GTu.n+4;  j++)    help.aux5.y[j] = conj( GTd.y[j] );

    help.aux1.convolutionAWT(help.aux4, help.aux5, 1, -1, help.aux6, help.aux7);
    help.aux1.KrammersKronigDown(help.aux1, help.K3, help.aux4, help.aux5, help.aux6, help.aux7);


    if(in.print_mode == 1)
    {
        cout << "BEintUp" << endl;
        name = dirname + "/ThetaPM1up";
        help.aux1.exportAWTasFUN(name, in.display, in.xMax, 1, in.output_mode);
    }


    // SECOND PART of THETA MINUS
    // first integrand
    for(int j = 0;          j < GTu.n+1;    j++)    help.aux4.y[j] = help.FD.y[j] * GTd.y[j] / (2 * Pi * unity );
    for(int j = GTu.n+1;    j < 3*GTu.n+4;  j++)    help.aux4.y[j] = 0;
    for(int j = 3*GTu.n+4;  j < 4*GTu.n+4;  j++)    help.aux4.y[j] = help.FD.y[j] * GTd.y[j] / (2 * Pi * unity );
    // second integrand is just the Green's function up

    help.aux2.convolutionAWT(help.aux4, GTu, -1, -1, help.aux6, help.aux7);
    help.aux2.KrammersKronigDown(help.aux2, help.K3, help.aux4, help.aux5, help.aux6, help.aux7);


    if(in.print_mode == 1)
    {
        cout << "BEintUp" << endl;
        name = dirname + "/ThetaPM2up";
        help.aux2.exportAWTasFUN(name, in.display, in.xMax, 1, in.output_mode);
    }


    // THIRD PART of THETA MINS
    // first integrand
    for(int j = 0;          j < GTu.n+1;    j++)    help.aux4.y[j] = help.FD.y[j] * conj( GTd.y[j] ) / (2 * Pi * unity );
    for(int j = GTu.n+1;    j < 3*GTu.n+4;  j++)    help.aux4.y[j] = 0;
    for(int j = 3*GTu.n+4;  j < 4*GTu.n+4;  j++)    help.aux4.y[j] = help.FD.y[j] * conj( GTd.y[j] ) / (2 * Pi * unity );
    // second integrand is just the Green's function up

    help.aux3.convolutionAWT(help.aux4, GTu, -1, -1, help.aux6, help.aux7);
    help.aux3.KrammersKronigDown(help.aux3, help.K3, help.aux4, help.aux5, help.aux6, help.aux7);


    if(in.print_mode == 1)
    {
        cout << "BEintUp" << endl;
        name = dirname + "/ThetaPM3up";
        help.aux3.exportAWTasFUN(name, in.display, in.xMax, 1, in.output_mode);
    }

    // THETA PLUS
    for(int j = 0;          j < GTu.n+1;    j++)    ThetaPM.y[j] = -gamma.mp * help.aux1.y[j] - gamma.pp * help.aux2.y[j] + gamma.mp * help.aux3.y[j];
    for(int j = GTu.n+1;    j < 3*GTu.n+4;  j++)    ThetaPM.y[j] = 0;
    for(int j = 3*GTu.n+4;  j < 4*GTu.n+4;  j++)    ThetaPM.y[j] = -gamma.mp * help.aux1.y[j] - gamma.pp * help.aux2.y[j] + gamma.mp * help.aux3.y[j];


    if(in.print_mode == 1)
    {
        cout << "BEintUp" << endl;
        name = dirname + "/ThetaPMup";
        ThetaPM.exportAWTasFUN(name, in.display, in.xMax, 1, in.output_mode);
    }

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //  SIGMA
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // FD integrals part
    // first integrand
    for(int j = 0;          j < GTu.n+1;    j++)    help.aux4.y[j] = in.U * help.FD.y[j] * imag( GTd.y[j] ) / Pi;
    for(int j = GTu.n+1;    j < 3*GTu.n+4;  j++)    help.aux4.y[j] = 0;
    for(int j = 3*GTu.n+4;  j < 4*GTu.n+4;  j++)    help.aux4.y[j] = in.U * help.FD.y[j] * imag( GTd.y[j] ) / Pi;
    // second integrand
    for(int j = 0;          j < GTu.n+1;    j++)    help.aux5.y[j] = ThetaPM.y[j] / ( conj( Dud.y[j] ) );
    for(int j = GTu.n+1;    j < 3*GTu.n+4;  j++)    help.aux5.y[j] = 0;
    for(int j = 3*GTu.n+4;  j < 4*GTu.n+4;  j++)    help.aux5.y[j] = ThetaPM.y[j] / ( conj( Dud.y[j] ) );

    help.aux5.deleteReal(help.aux5);

    help.aux1.convolutionAWT(help.aux4, help.aux5, -1, -1, help.aux6, help.aux7);
    help.aux1.KrammersKronig(help.aux1, help.K3, help.aux4, help.aux5, help.aux6, help.aux7);

    if(in.print_mode == 1)
    {
        cout << "BEintUp" << endl;
        name = dirname + "/FDintUp";
        help.aux1.exportAWTasFUN(name, in.display, in.xMax, 1, in.output_mode);
    }

    for(int j = 0;          j < GTu.n+1;    j++)    help.aux6.y[j] = help.aux4.y[j] * help.aux5.y[j];
    for(int j = GTu.n+1;    j < 3*GTu.n+4;  j++)    help.aux6.y[j] = 0;
    for(int j = 3*GTu.n+4;  j < 4*GTu.n+4;  j++)    help.aux6.y[j] = help.aux4.y[j] * help.aux5.y[j];
    SigmaUp.y[0] = Simpson_Re(help.aux6) + unity * Simpson_Im(help.aux6);


    for(int j = 0;          j < GTu.n+1;    j++)    help.aux4.y[j] = in.U * help.BE.y[j] * imag( ThetaPP.y[j]/Dud.y[j]  ) / Pi;
    for(int j = GTu.n+1;    j < 3*GTu.n+4;  j++)    help.aux4.y[j] = 0;
    for(int j = 3*GTu.n+4;  j < 4*GTu.n+4;  j++)    help.aux4.y[j] = in.U * help.BE.y[j] * imag( ThetaPP.y[j]/Dud.y[j]  ) / Pi;

    // second integrand
    for(int j = 0;          j < GTu.n+1;    j++)    help.aux5.y[j] = GTd.y[j];
    for(int j = GTu.n+1;    j < 3*GTu.n+4;  j++)    help.aux5.y[j] = 0;
    for(int j = 3*GTu.n+4;  j < 4*GTu.n+4;  j++)    help.aux5.y[j] = GTd.y[j];

    help.aux5.deleteReal(help.aux5);

    help.aux2.convolutionAWT(help.aux4, help.aux5, 1, -1, help.aux6, help.aux7);
    help.aux2.KrammersKronig(help.aux2, help.K3, help.aux4, help.aux5, help.aux6, help.aux7);


    if(in.print_mode == 1)
    {
        cout << "BEintUp" << endl;
        name = dirname + "/BEintUp";
        help.aux2.exportAWTasFUN(name, in.display, in.xMax, 1, in.output_mode);
    }

    for(int j = 0;          j < GTu.n+1;    j++)    help.aux6.y[j] = help.aux4.y[j] * help.aux5.y[j];
    for(int j = GTu.n+1;    j < 3*GTu.n+4;  j++)    help.aux6.y[j] = 0;
    for(int j = 3*GTu.n+4;  j < 4*GTu.n+4;  j++)    help.aux6.y[j] = help.aux4.y[j] * help.aux5.y[j];

    SigmaUp.y[0] = SigmaUp.y[0] - Simpson_Re(help.aux6) - unity * Simpson_Im(help.aux6);


    for(int j=0; j < 4*GTu.n + 4; j++ ) SigmaUp.y[j] = help.aux1.y[j] - help.aux2.y[j];


    // export
    if(in.print_mode == 1)
    {
        cout << "SSu" << endl;
        name = dirname + "/SSu";
        SigmaUp.exportAWTasFUN(name, in.display, in.xMax, 1, in.output_mode);
    }

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //  THETA CALCULATIONS Sigma DOWN
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // FIRST PART of THETA PLUS down
    // first integrand
    for(int j = 0;          j < GTu.n+1;    j++)    help.aux4.y[j] = help.FD.y[j] * imag( GTd.y[j] ) / Pi;
    for(int j = GTu.n+1;    j < 3*GTu.n+4;  j++)    help.aux4.y[j] = 0;
    for(int j = 3*GTu.n+4;  j < 4*GTu.n+4;  j++)    help.aux4.y[j] = help.FD.y[j] * imag( GTd.y[j] ) / Pi;
    // second integrand
    for(int j = 0;          j < GTu.n+1;    j++)    help.aux5.y[j] = GTu.y[j];
    for(int j = GTu.n+1;    j < 3*GTu.n+4;  j++)    help.aux5.y[j] = 0;
    for(int j = 3*GTu.n+4;  j < 4*GTu.n+4;  j++)    help.aux5.y[j] = GTu.y[j];

    help.aux1.convolutionAWT(help.aux4, help.aux5, 1, -1, help.aux6, help.aux7);
    help.aux1.KrammersKronig(help.aux1, help.K3, help.aux4, help.aux5, help.aux6, help.aux7);


    if(in.print_mode == 1)
    {
        cout << "ThetaPP1down" << endl;
        name = dirname + "/ThetaPP1down";
        help.aux1.exportAWTasFUN(name, in.display, in.xMax, 1, in.output_mode);
    }

    // SECOND PART of THETA PLUS down
    // first integrand
    for(int j = 0;          j < GTu.n+1;    j++)    help.aux4.y[j] = help.FD.y[j] * GTu.y[j] / (2 * Pi * unity );
    for(int j = GTu.n+1;    j < 3*GTu.n+4;  j++)    help.aux4.y[j] = 0;
    for(int j = 3*GTu.n+4;  j < 4*GTu.n+4;  j++)    help.aux4.y[j] = help.FD.y[j] * GTu.y[j] / (2 * Pi * unity );
    // second integrand
    for(int j = 0;          j < GTu.n+1;    j++)    help.aux5.y[j] = conj( GTd.y[j] );
    for(int j = GTu.n+1;    j < 3*GTu.n+4;  j++)    help.aux5.y[j] = 0;
    for(int j = 3*GTu.n+4;  j < 4*GTu.n+4;  j++)    help.aux5.y[j] = conj( GTd.y[j] );

    help.aux2.convolutionAWT(help.aux4, help.aux5, -1, -1, help.aux6, help.aux7);
    help.aux2.KrammersKronig(help.aux2, help.K3, help.aux4, help.aux5, help.aux6, help.aux7);


    if(in.print_mode == 1)
    {
        cout << "ThetaPP2down" << endl;
        name = dirname + "/ThetaPP2down";
        help.aux2.exportAWTasFUN(name, in.display, in.xMax, 1, in.output_mode);
    }

    // THIRD PART of THETA PLUS down
    // first integrand
    for(int j = 0;          j < GTu.n+1;    j++)    help.aux4.y[j] = help.FD.y[j] * conj( GTu.y[j] ) / (2 * Pi * unity );
    for(int j = GTu.n+1;    j < 3*GTu.n+4;  j++)    help.aux4.y[j] = 0;
    for(int j = 3*GTu.n+4;  j < 4*GTu.n+4;  j++)    help.aux4.y[j] = help.FD.y[j] * conj( GTu.y[j] ) / (2 * Pi * unity );
    // second integrand
    for(int j = 0;          j < GTu.n+1;    j++)    help.aux5.y[j] = conj( GTd.y[j] );
    for(int j = GTu.n+1;    j < 3*GTu.n+4;  j++)    help.aux5.y[j] = 0;
    for(int j = 3*GTu.n+4;  j < 4*GTu.n+4;  j++)    help.aux5.y[j] = conj( GTd.y[j] );

    help.aux3.convolutionAWT(help.aux4, help.aux5, -1, -1, help.aux6, help.aux7);
    help.aux3.KrammersKronig(help.aux3, help.K3, help.aux4, help.aux5, help.aux6, help.aux7);

    if(in.print_mode == 1)
    {
        cout << "ThetaPP3down" << endl;
        name = dirname + "/ThetaPP3down";
        help.aux3.exportAWTasFUN(name, in.display, in.xMax, 1, in.output_mode);
    }

    // THETA PLUS
    for(int j = 0;          j < GTu.n+1;    j++)    ThetaPP.y[j] = -gamma.pp * help.aux1.y[j] - gamma.pp * help.aux2.y[j] + gamma.pm * help.aux3.y[j];
    for(int j = GTu.n+1;    j < 3*GTu.n+4;  j++)    ThetaPP.y[j] = 0;
    for(int j = 3*GTu.n+4;  j < 4*GTu.n+4;  j++)    ThetaPP.y[j] = -gamma.pp * help.aux1.y[j] - gamma.pp * help.aux2.y[j] + gamma.pm * help.aux3.y[j];


    if(in.print_mode == 1)
    {
        cout << "ThetaPPdown" << endl;
        name = dirname + "/ThetaPPdown";
        ThetaPP.exportAWTasFUN(name, in.display, in.xMax, 1, in.output_mode);
    }

    // FIRST PART of THETA MINUS
    // first integrand
    for(int j = 0;          j < GTu.n+1;    j++)    help.aux4.y[j] = help.FD.y[j] * imag( GTd.y[j] ) / Pi;
    for(int j = GTu.n+1;    j < 3*GTu.n+4;  j++)    help.aux4.y[j] = 0;
    for(int j = 3*GTu.n+4;  j < 4*GTu.n+4;  j++)    help.aux4.y[j] = help.FD.y[j] * imag( GTd.y[j] ) / Pi;
    // second integrand
    for(int j = 0;          j < GTu.n+1;    j++)    help.aux5.y[j] = conj( GTu.y[j] );
    for(int j = GTu.n+1;    j < 3*GTu.n+4;  j++)    help.aux5.y[j] = 0;
    for(int j = 3*GTu.n+4;  j < 4*GTu.n+4;  j++)    help.aux5.y[j] = conj( GTu.y[j] );

    help.aux1.convolutionAWT(help.aux4, help.aux5, 1, -1, help.aux6, help.aux7);
    help.aux1.KrammersKronigDown(help.aux1, help.K3, help.aux4, help.aux5, help.aux6, help.aux7);


    if(in.print_mode == 1)
    {
        cout << "ThetaPM1down" << endl;
        name = dirname + "/ThetaPM1down";
        help.aux1.exportAWTasFUN(name, in.display, in.xMax, 1, in.output_mode);
    }

    // SECOND PART of THETA MINUS
    // first integrand
    for(int j = 0;          j < GTu.n+1;    j++)    help.aux4.y[j] = help.FD.y[j] * GTu.y[j] / (2 * Pi * unity );
    for(int j = GTu.n+1;    j < 3*GTu.n+4;  j++)    help.aux4.y[j] = 0;
    for(int j = 3*GTu.n+4;  j < 4*GTu.n+4;  j++)    help.aux4.y[j] = help.FD.y[j] * GTu.y[j] / (2 * Pi * unity );
    // second integrand is just the Green's function up

    help.aux2.convolutionAWT(help.aux4, GTd, -1, -1, help.aux6, help.aux7);
    help.aux2.KrammersKronigDown(help.aux2, help.K3, help.aux4, help.aux5, help.aux6, help.aux7);


    if(in.print_mode == 1)
    {
        cout << "ThetaPM2down" << endl;
        name = dirname + "/ThetaPM2down";
        help.aux2.exportAWTasFUN(name, in.display, in.xMax, 1, in.output_mode);
    }

    // THIRD PART of THETA MINS
    // first integrand
    for(int j = 0;          j < GTu.n+1;    j++)    help.aux4.y[j] = help.FD.y[j] * conj( GTu.y[j] ) / (2 * Pi * unity );
    for(int j = GTu.n+1;    j < 3*GTu.n+4;  j++)    help.aux4.y[j] = 0;
    for(int j = 3*GTu.n+4;  j < 4*GTu.n+4;  j++)    help.aux4.y[j] = help.FD.y[j] * conj( GTu.y[j] ) / (2 * Pi * unity );
    // second integrand is just the Green's function up

    help.aux3.convolutionAWT(help.aux4, GTd, -1, -1, help.aux6, help.aux7);
    help.aux3.KrammersKronigDown(help.aux3, help.K3, help.aux4, help.aux5, help.aux6, help.aux7);


    if(in.print_mode == 1)
    {
        cout << "ThetaPM3down" << endl;
        name = dirname + "/ThetaPM3down";
        help.aux3.exportAWTasFUN(name, in.display, in.xMax, 1, in.output_mode);
    }

    // THETA PLUS
    for(int j = 0;          j < GTu.n+1;    j++)    ThetaPM.y[j] = -gamma.mp * help.aux1.y[j] - gamma.pp * help.aux2.y[j] + gamma.mp * help.aux3.y[j];
    for(int j = GTu.n+1;    j < 3*GTu.n+4;  j++)    ThetaPM.y[j] = 0;
    for(int j = 3*GTu.n+4;  j < 4*GTu.n+4;  j++)    ThetaPM.y[j] = -gamma.mp * help.aux1.y[j] - gamma.pp * help.aux2.y[j] + gamma.mp * help.aux3.y[j];


    if(in.print_mode == 1)
    {
        cout << "ThetaPMdown" << endl;
        name = dirname + "/ThetaPMdown";
        ThetaPM.exportAWTasFUN(name, in.display, in.xMax, 1, in.output_mode);
    }

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //  SIGMA Down
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // FD integrals part
    // first integrand
    for(int j = 0;          j < GTu.n+1;    j++)    help.aux4.y[j] = in.U * help.FD.y[j] * imag( GTu.y[j] ) / Pi;
    for(int j = GTu.n+1;    j < 3*GTu.n+4;  j++)    help.aux4.y[j] = 0;
    for(int j = 3*GTu.n+4;  j < 4*GTu.n+4;  j++)    help.aux4.y[j] = in.U * help.FD.y[j] * imag( GTu.y[j] ) / Pi;
    // second integrand
    for(int j = 0;          j < GTu.n+1;    j++)    help.aux5.y[j] = ThetaPM.y[j] / ( conj( Ddu.y[j] ) );
    for(int j = GTu.n+1;    j < 3*GTu.n+4;  j++)    help.aux5.y[j] = 0;
    for(int j = 3*GTu.n+4;  j < 4*GTu.n+4;  j++)    help.aux5.y[j] = ThetaPM.y[j] / ( conj( Ddu.y[j] ) );

    help.aux5.deleteReal(help.aux5);

    help.aux1.convolutionAWT(help.aux4, help.aux5, -1, -1, help.aux6, help.aux7);
    help.aux1.KrammersKronig(help.aux1, help.K3, help.aux4, help.aux5, help.aux6, help.aux7);

    if(in.print_mode == 1)
    {
        cout << "FDintDown" << endl;
        name = dirname + "/FDintDown";
        help.aux1.exportAWTasFUN(name, in.display, in.xMax, 1, in.output_mode);
    }

    for(int j = 0;          j < GTu.n+1;    j++)    help.aux6.y[j] = help.aux4.y[j] * help.aux5.y[j];
    for(int j = GTu.n+1;    j < 3*GTu.n+4;  j++)    help.aux6.y[j] = 0;
    for(int j = 3*GTu.n+4;  j < 4*GTu.n+4;  j++)    help.aux6.y[j] = help.aux4.y[j] * help.aux5.y[j];

    SigmaDown.y[0] = Simpson_Re(help.aux6) + unity * Simpson_Im(help.aux6);


    for(int j = 0;          j < GTu.n+1;    j++)    help.aux4.y[j] = in.U * help.BE.y[j] * imag( ThetaPP.y[j]/Ddu.y[j]  ) / Pi;
    for(int j = GTu.n+1;    j < 3*GTu.n+4;  j++)    help.aux4.y[j] = 0;
    for(int j = 3*GTu.n+4;  j < 4*GTu.n+4;  j++)    help.aux4.y[j] = in.U * help.BE.y[j] * imag( ThetaPP.y[j]/Ddu.y[j]  ) / Pi;

    // second integrand
    for(int j = 0;          j < GTu.n+1;    j++)    help.aux5.y[j] = GTu.y[j];
    for(int j = GTu.n+1;    j < 3*GTu.n+4;  j++)    help.aux5.y[j] = 0;
    for(int j = 3*GTu.n+4;  j < 4*GTu.n+4;  j++)    help.aux5.y[j] = GTu.y[j];

    help.aux5.deleteReal(help.aux5);

    help.aux2.convolutionAWT(help.aux4, help.aux5, 1, -1, help.aux6, help.aux7);
    help.aux2.KrammersKronig(help.aux2, help.K3, help.aux4, help.aux5, help.aux6, help.aux7);


    if(in.print_mode == 1)
    {
        cout << "BEintDown" << endl;
        name = dirname + "/BEintDown";
        help.aux2.exportAWTasFUN(name, in.display, in.xMax, 1, in.output_mode);
    }

    for(int j = 0;          j < GTu.n+1;    j++)    help.aux6.y[j] = help.aux4.y[j] * help.aux5.y[j];
    for(int j = GTu.n+1;    j < 3*GTu.n+4;  j++)    help.aux6.y[j] = 0;
    for(int j = 3*GTu.n+4;  j < 4*GTu.n+4;  j++)    help.aux6.y[j] = help.aux4.y[j] * help.aux5.y[j];

    SigmaDown.y[0] = SigmaUp.y[0] - Simpson_Re(help.aux6) - unity * Simpson_Im(help.aux6);


    for(int j=0; j < 4*GTu.n + 4; j++ ) SigmaDown.y[j] = help.aux1.y[j] - help.aux2.y[j];


    // export
    if(in.print_mode == 1)
    {
        name = dirname + "/SSd";
        SigmaDown.exportAWTasFUN(name, in.display, in.xMax, 1, in.output_mode);
    }

}

// FIX rewrite it using auxiliary and info structs
void calculate_sigma_stat_directly_SSN(AWT & Sigma, info in, double U, double kT, AWT & Dud, AWT & Ddu, matrix & Lud, matrix & Ldu,
                                    AWT & GTd, AWT & GTu, AWT & FD, AWT & BE, AWT & K3, string dirname, int print,
                                    AWT & aux1, AWT & aux2, AWT & aux3, AWT & aux4, AWT & aux5, AWT & aux6, AWT & aux7)
{
    complex<double> unity(0,1);
    double Pi = 3.14159265359;

    AWT ThetaPP, ThetaPM;
    ThetaPP.initializeAWT(GTu.n, GTu.xMax, kT);
    ThetaPM.initializeAWT(GTu.n, GTu.xMax, kT);

    // DIRECT INTEGRALS
    aux4.setZero();
    aux5.setZero();
    aux6.setZero();
    aux7.setZero();

    int divisor = 40;

    for(int alpha = -GTu.n/divisor; alpha < GTu.n/divisor + 1; alpha++)
    {
        for(int j = -GTu.n + GTu.n/divisor;          j < GTu.n - GTu.n/divisor + 1;          j++)
        {
            int dif;
            if(-alpha + j < 0)	dif = -alpha + j + 4*GTu.n + 4;
            else			    dif = -alpha + j;

            int x;
            if(j < 0)		x = j + 4*GTu.n + 4;
            else			x = j;

            aux6.y[x] = FD.y[x] * imag( GTd.y[x] ) * ThetaPM.y[dif] / ( conj( Dud.y[dif] ) ) ;
        }
        // zero padding
        for(int j=GTu.n - GTu.n/divisor +1; j < 3*GTu.n + GTu.n/divisor + 4; j++ )
        {
            aux6.y[j] = 0;
        }

        int omega;
        if(alpha < 0)		    omega = alpha + 4*GTu.n + 4;
        else			        omega = alpha;

        aux4.y[omega] = Simpson_Re(aux6) + unity * Simpson_Im(aux6);

        if( alpha % (100) == 0 )    cout << "calculating " << alpha << endl;
    }

    for(int alpha = -GTu.n/divisor; alpha < GTu.n/divisor + 1; alpha++)
    {
        for(int j = -GTu.n + GTu.n/divisor;          j < GTu.n - GTu.n/divisor + 1;          j++)
        {
            int sum;
            if(alpha + j < 0)	sum = alpha + j + 4*GTu.n + 4;
            else			    sum = alpha + j;

            int x;
            if(j < 0)		x = j + 4*GTu.n + 4;
            else			x = j;

            aux7.y[x] = BE.y[x] * imag( ThetaPP.y[x]/Dud.y[x] ) * GTd.y[sum];
        }
        // zero padding
        for(int j=GTu.n - GTu.n/divisor +1; j < 3*GTu.n + GTu.n/divisor + 4; j++ )
        {
            aux6.y[j] = 0;
        }

        int omega;
        if(alpha < 0)		omega = alpha + 4*GTu.n + 4;
        else			    omega = alpha;

        aux5.y[omega] = Simpson_Re(aux7) + unity * Simpson_Im(aux7);

        if( alpha % (divisor*100) == 0 )    cout << "calculating " << alpha << endl;
    }


    // sum together both parts
    for(int j = 0; j < 4*GTu.n+4; j++)      aux7.y[j] = U * ( aux4.y[j] - aux5.y[j] ) / Pi;

    // export
    string name = dirname + "/SSdir";
    aux7.exportAWTasFUN(name, in.display, in.xMax, 1, in.output_mode);

}

void detailed_output_SSN(info in, string dirname, AWT & GdGu1pp, AWT & GdGu1pm, AWT & GdGu2pm, AWT & GdGu2mm,
                     AWT & Dud, AWT & Ddu, AWT & GTd, AWT & GTu)
{
    string name;

    name = dirname + "/GdGu1pp_";
    GdGu1pp.exportAWTasFUN(name, in.display, in.xMax, 1, in.output_mode);

    name = dirname + "/GdGu1pm_";
    GdGu1pm.exportAWTasFUN(name, in.display, in.xMax, 1, in.output_mode);

    name = dirname + "/GdGu2pm_";
    GdGu2pm.exportAWTasFUN(name, in.display, in.xMax, 1, in.output_mode);

    name = dirname + "/GdGu2mm_";
    GdGu2mm.exportAWTasFUN(name, in.display, in.xMax, 1, in.output_mode);

    name = dirname + "/Dud";
        Dud.exportAWTasFUN(name, in.display, in.xMax, 1, in.output_mode);

    name = dirname + "/Ddu";
        Dud.exportAWTasFUN(name, in.display, in.xMax, 1, in.output_mode);

}



void normal_output_SSN(info in, string dirname, AWT & GTd, AWT & GTu, AWT & SSd, AWT & SSu, AWT & GSd, AWT & GSu, ofstream & log_output, ofstream & data_output)
{
    string name;
    name = dirname +   "/GTd";
    GTd.exportAWTasFUN(name, in.display, in.xMax, 1, in.output_mode);

    name = dirname +   "/GTu";
    GTu.exportAWTasFUN(name, in.display, in.xMax, 1, in.output_mode);

    name = dirname + "/SSu"; // + boost::lexical_cast<string>(iter_lambda);
    SSu.exportAWTasFUN(name, in.display, in.xMax, 1, in.output_mode);

    name = dirname + "/GSu"; // + boost::lexical_cast<string>(iter_lambda);
    GSu.exportAWTasFUN(name, in.display, in.xMax, 1, in.output_mode);

    name = dirname + "/SSd"; // + boost::lexical_cast<string>(iter_lambda);
    SSd.exportAWTasFUN(name, in.display, in.xMax, 1, in.output_mode);

    name = dirname + "/GSd"; // + boost::lexical_cast<string>(iter_lambda);
    GSd.exportAWTasFUN(name, in.display, in.xMax, 1, in.output_mode);
}


void static_spinpolarized_SSN(string input)
{
    // constants
    complex<double> unity(0,1.0);


    // set iteration precision
    int oscil_control = 1;      // controls the oscillatory behaviour of the problem

    // input of values defining given model
    info in;
    in.import(input);
    in.precisions();

    // thermodynamic occupation nT and magnetization mT defined and preset
    double nT, mT;
    nT = 1.0;   // half-filling value is preset
    mT = 0.0;   // half-filling value is preset
    double nT_s, mT_s;

    // thermodynamic (static) self-energy
    complex<double>     aSTu = 0;
    complex<double>     aSTd = 0;

    // spectral (dynamic) self-energies, loaded with zeros
    AWT SSu, SSd;
    SSu.initializeAWT(in.n, in.xMax, in.kT);
    SSd.initializeAWT(in.n, in.xMax, in.kT);
    SSu.setZero();
    SSd.setZero();


    // thermodynamic and spectral GREEN FUNCTIONS arrays
    AWT GTu, GTd, GSu, GSd;
    GTu.initializeAWT(in.n, in.xMax, in.kT);
    GTd.initializeAWT(in.n, in.xMax, in.kT);
    GSu.initializeAWT(in.n, in.xMax, in.kT);
    GSd.initializeAWT(in.n, in.xMax, in.kT);


    // determinant Dud, Ddu
    AWT Dud, Ddu;
    Dud.initializeAWT(in.n, in.xMax, in.kT);
    Dud.setZero();
    Ddu.initializeAWT(in.n, in.xMax, in.kT);
    Ddu.setZero();

    // Divergent vertices ud, du
    matrix Kud, Lud, FDud, BDud;
    matrix Kdu, Ldu, FDdu, BDdu;

    // suffix _l means value of given variable from previous lambda iteration
    matrix Kud_l, Lud_l, FDud_l, BDud_l;
    matrix Kdu_l, Ldu_l, FDdu_l, BDdu_l;


    // bubbles for ud
    AWT GdGu1pp, GdGu1pm, GdGu2pm, GdGu2mm;
      GdGu1pp.initializeAWT(in.n, in.xMax, in.kT);
      GdGu1pm.initializeAWT(in.n, in.xMax, in.kT);
      GdGu2pm.initializeAWT(in.n, in.xMax, in.kT);
      GdGu2mm.initializeAWT(in.n, in.xMax, in.kT);

      GdGu1pp.setZero();
      GdGu1pm.setZero();
      GdGu2pm.setZero();
      GdGu2mm.setZero();


    // bubbles for du
    AWT GuGd1pp, GuGd1pm, GuGd2pm, GuGd2mm;
      GuGd1pp.initializeAWT(in.n, in.xMax, in.kT);
      GuGd1pm.initializeAWT(in.n, in.xMax, in.kT);
      GuGd2pm.initializeAWT(in.n, in.xMax, in.kT);
      GuGd2mm.initializeAWT(in.n, in.xMax, in.kT);

      GuGd1pp.setZero();
      GuGd1pm.setZero();
      GuGd2pm.setZero();
      GuGd2mm.setZero();

    // the auxiliary structure
    auxiliary aux;
    aux.initialize(in.n,in.xMax, in.kT);


    if(in.print_mode == 1)  cout << "Sigma-s, G-s, Lambda, K initialized" << endl;


    ofstream log_output;
    ofstream data_output;

    log_output << endl << __TIME__ << ",   " << __DATE__ << endl;


    // subdirectories for different calculations are created
    string     U_string  =  double_to_string(in.U, in.U_precision);
    string     h_string  =  double_to_string(in.h,  3);
    string     x_string  =  double_to_string(in.mu, 2);
    string     kT_string =  double_to_string(in.kT, 3);

    string input_string;

    cout << input_string << endl;

    string dirname = "Txt/U" + U_string + "_h" + h_string + "_x" + x_string + "_kT" + kT_string +  input_string;
    createDirectories(dirname);

    welcome_screen_SSN(in, dirname);

    // we start with half filling values which we then corrected iteratively
    aSTu = 0;
    aSTd = 0;


    string name;

    // logs
    name = dirname + "/log";
    log_output.open(name, std::ios_base::app);
    name = dirname + "/data";
    data_output.open(name, std::ios_base::app);

    set_matrix_real(   Lud, in.U/2);
    set_matrix_real( Lud_l, 0);
    set_matrix_real(   Kud, 0);
    set_matrix_real(  FDud, 0);
    set_matrix_real(  BDud, 0);


    set_matrix_real(   Ldu, in.U/2);
    set_matrix_real( Ldu_l,  0);
    set_matrix_real(    Kdu, 0);
    set_matrix_real(   FDdu, 0);
    set_matrix_real(   BDdu, 0);

    // output screen talks and talks
    int out = 1;

    int iter_total = 0;
    int iter_lambda;
    int iter_sigma = 0;


    //for(int iter_sigma=0; iter_sigma<sigma_iter_max; iter_sigma++)
    for(iter_sigma = 0; iter_sigma<1; iter_sigma++)
    //while(abs(real(totalSTd) - real(totalSTd_s) ) > 0.00001)
    {

        cout << endl << "   SIGMA " << iter_sigma << "-th iter: " << endl << endl;

        // iterating vertices
        iter_lambda = 0;
        while( (iter_lambda == 0) || (abs(real(Lud.mp) - real(Lud_l.mp) ) > 0.00001)  )
        //for(int iter_lambda=0; iter_lambda<1; iter_lambda++)
        {
            cout << "       " << iter_lambda << "-th iter: " ;
            if(out == 1)   cout << endl << "          ";

            // thermodynamic propagators with given hybridizations are prepared
            if(in.mu == 0)  nT = 1.0;

            // 0. step: setting the initial GTd, GTu
            lorentz_propagator(in, GTd);
            insert_hybridization(in, GTd);
            insert_thermal_selfenergy(in, GTd, in.U * (1.0 - nT) / 2.0 + in.h - aSTd);

            lorentz_propagator(in, GTu);
            insert_hybridization(in, GTu);
            insert_thermal_selfenergy(in, GTu, in.U * (1.0 - nT) / 2.0 - in.h + aSTu);

            //          calculate bubbles just once since they do not change
            if(iter_lambda==0)    calculateGGud_SSN(GdGu1pp, GdGu1pm, GdGu2pm, GdGu2mm, GTd, GTu, aux, out);


            // 1. step: calculate determinant
            calculate_Dud_static_SSN(Dud, GdGu1pp, GdGu1pm, GdGu2pm, GdGu2mm, Lud, out);
            calculate_Ddu_sym_SSN(Ddu, Dud, out);

            // 2. and 3. step: calculate Kud
            if(iter_sigma>0 && oscil_control == 1)
            {
                set_matrix( Kud_l, Kud);
                set_matrix(FDud_l, FDud);
                set_matrix(BDud_l, BDud);

                set_matrix( Kdu_l, Kdu);
                set_matrix(FDdu_l, FDdu);
                set_matrix(BDdu_l, BDdu);
            }

            calculate_Kud_static_SSN(Kud, GdGu1pp, GdGu1pm, GdGu2pm, GdGu2mm, Lud, out);
            calculate_Kdu_sym_SSN(Kdu, Kud, out);

            calculate_FDud_SSN(FDud, GTd, GTu, Dud, aux, out);
            calculate_FDdu_sym_SSN(FDdu, FDud, out);

            calculate_BDud_SSN(BDud, GTd, GTu, Dud, aux, out);
            calculate_BDdu_sym_SSN(BDdu, BDud, out);

            if(iter_sigma>0 && oscil_control == 1)
            {
                cout << "aver, ";
                aver_matrix( Kud,  Kud_l);
                aver_matrix(FDud, FDud_l);
                aver_matrix(BDud, BDud_l);

                aver_matrix( Kdu,  Kdu_l);
                aver_matrix(FDdu, FDdu_l);
                aver_matrix(BDdu, BDdu_l);
            }

            // 4. step: calculate Lud
            set_matrix(Lud_l, Lud);
            set_matrix(Ldu_l, Ldu);

            calculate_Lud_SSN(Lud, Kud, FDud, BDud, in.U, out);
            calculate_Ldu_sym_SSN(Ldu, Lud, out);


            if(oscil_control == 1)
            {
                aver_matrix(Lud, Lud_l);
                aver_matrix(Ldu, Ldu_l);
            }


            lambda_iter_console_SSN(out, Lud, Ldu, Kud, Kdu, BDud, BDdu, FDud, FDdu);

            if((iter_lambda == 0) && (iter_sigma == 0) ) log_output << endl << __TIME__ << ",   " << __DATE__ << endl << endl;

            lambda_iter_output_SSN(log_output, iter_lambda, iter_sigma, aSTu, aSTd, Lud, Ldu, Kud, Kdu, BDud, BDdu, FDud, FDdu);


            // calculate new anomalous self-energy
            calculate_aST_SSN(aSTu, aSTd, Lud, Ldu, GTd, GTu, aux, out);


            lambda_log_output_SSN(data_output, iter_total, Lud, Ldu, Kud, Kdu, BDud, BDdu, FDud, FDdu);


            iter_lambda = iter_lambda + 1;
            iter_total  = iter_total  + 1;

            calculate_nT_mT_SSN(nT, mT, in.mu, in.U, GTu, GTd, aux);

        } // end of Lambda iterations


        // store old values
        if(iter_sigma>0 && oscil_control ==0)
        {
            nT_s = nT;
            mT_s = mT;
        }

        // anomalous contribution to thermodynamic self-energy
        calculate_nT_mT_SSN(nT, mT, in.mu, in.U, GTu, GTd, aux);


        if(iter_sigma>0 && oscil_control ==0)
        {
            cout << "aver nT, sigma" << endl;
            nT = (nT + nT_s) / 2.0;
            mT = (mT + mT_s) / 2.0;
        }


              cout << "          nT, mT: " << nT << " " << mT  << endl << "          ";
        log_output << "          nT, mT: " << nT << " " << mT  << endl << endl;

        // shift the leading thermodynamic sigma iterator
        iter_sigma = iter_sigma + 1;

    } // end of Sigma iterations


    // output of final lambda iteration results
    ofstream final_output;
    name = dirname + "/final";
    final_output.open(name, std::ios_base::app);

    lambda_iter_output_SSN(final_output, iter_lambda, iter_sigma, aSTu, aSTd, Lud, Ldu, Kud, Kdu, BDud, BDdu, FDud, FDdu);



    // iterate nS
    int iter_nS = 0;
    double nS = 1;
    double mS = 0;

    double nS_old = -1.0;

    AWT SS0;
    SS0.initializeAWT(in.n, in.xMax, in.kT);

    while( ( abs(nS - nS_old) > 0.000001 ) && iter_nS < 1000 )
    {

        calculate_dynamic_selfenergy_SSN(SSu, SSd, in, Dud, Ddu, Lud, Ldu, GTd, GTu, dirname, aux);
        for(int j = 0;          j < GTu.n+1;    j++)    SS0.y[j] = ( SSu.y[j] + SSd.y[j] ) / 2.0;
        for(int j = GTu.n+1;    j < 3*GTu.n+4;  j++)    SS0.y[j] = 0;
        for(int j = 3*GTu.n+4;  j < 4*GTu.n+4;  j++)    SS0.y[j] = ( SSu.y[j] + SSd.y[j] ) / 2.0;

        //propagatorLorentzShift_SSN(in.delta, in.mu + in.U / 2.0 - in.U * nS / 2.0, - in.h + aSTu, GSu, SS0);

          lorentz_propagator(in, GSu);
        insert_hybridization(in, GSu);
        insert_thermal_selfenergy(in, GSu, in.mu + in.U / 2.0 - in.U * nS / 2.0 - in.h + aSTu);
        insert_dynamic_selfenergy(in, GSu, SS0);

        //propagatorLorentzShift_SSN(in.delta, in.mu + in.U / 2.0 - in.U * nS / 2.0, + in.h - aSTd, GSd, SS0);

          lorentz_propagator(in, GSd);
        insert_hybridization(in, GSd);
        insert_thermal_selfenergy(in, GSd, in.mu + in.U / 2.0 - in.U * nS / 2.0 + in.h - aSTd);
        insert_dynamic_selfenergy(in, GSd, SS0);

        calculate_nS_mS_SSN(nS, mS, GSu, GSd, aux);
        nS_old = nS;

        nS = ( nS + nS_old ) / 2.0;

              cout << "          nS: " << nS <<  endl << "          ";
        log_output << "          nS: " << nS << endl << endl;

        iter_nS = iter_nS + 1;
    }


    // OUTPUTS
    if(in.print_mode == 1)
    {
        cout << "print details";
        detailed_output_SSN(in, dirname, GdGu1pp, GdGu1pm, GdGu2pm, GdGu2mm, Dud, Ddu, GTd, GTu);
    }

    normal_output_SSN(in, dirname, GTd, GTu, SSd, SSu, GSd, GSu, log_output, data_output);



} // end of the whole function

void check_transformation(string input)
{
    // constants
    complex<double> unity(0.0,1.0);

        // input of values defining given model
    info in;
    in.import(input);
    in.precisions();


    // thermodynamic occupation nT and magnetization mT defined and preset
    double nT, mT;
    nT = 1.0;   // half-filling value is preset
    mT = 0.0;   // half-filling value is preset
    double nT_s, mT_s;

    // thermodynamic (static) self-energy
    complex<double>     aSTu = 0;
    complex<double>     aSTd = 0;

    // thermodynamic and spectral GREEN FUNCTIONS arrays
    AWT GTu, GTd;
    GTu.initializeAWT(in.n, in.xMax, in.kT);
    GTd.initializeAWT(in.n, in.xMax, in.kT);

    auxiliary help;
    help.initialize(in.n, in.xMax, in.kT);


    AWT hybrid;
    hybrid.initializeAWT(in.n, in.xMax, in.kT);

    for(int j = 0;         j < in.n+1;       j++)    hybrid.y[j] = 1.0;
    for(int j = in.n+1;    j < 3*in.n+4;     j++)    hybrid.y[j] = 0;
    for(int j = 3*in.n+4;  j < 4*in.n+4;     j++)    hybrid.y[j] = 1.0;


    in.gap =0.00999999;

    int i_gap = in.n*in.gap/in.xMax;

    hybrid.KrammersKronig(hybrid, help.K3, help.aux4, help.aux5, help.aux6, help.aux7);
    in.phi = 0.5;
/*

    for(int j = i_gap;     j < in.n+1;         j++)
    {
        double x = j*in.xMax/in.n;
        hybrid.y[j] = hybrid.y[j] + unity*2.0*(1-(in.gap/x)*cos(in.phi*Pi/2.0))*in.gammaS*abs(x)/sqrt(x*x - in.gap*in.gap);
    }
    for(int j = in.n+1;    j < 3*in.n+4;       j++)
    {
        hybrid.y[j] = 0;
    }
    for(int j = 3*in.n+4;  j < 4*in.n+4-i_gap; j++)
    {
        double x = (4*in.n+4-j)*in.xMax/in.n;
        hybrid.y[j] = hybrid.y[j] + unity*2.0*(1-(in.gap/x)*cos(in.phi*Pi/2.0))*in.gammaS*abs(x)/sqrt(x*x - in.gap*in.gap);
    }

*/
/*
    for(int j = i_gap;     j < in.n+1;         j++)
    {
        double x = j*in.xMax/in.n;
        hybrid.y[j] = hybrid.y[j] + unity*in.gammaS;
    }
    for(int j = in.n+1;    j < 3*in.n+4;       j++)
    {
        hybrid.y[j] = 0;
    }
    for(int j = 3*in.n+4;  j < 4*in.n+4-i_gap; j++)
    {
        double x = (4*in.n+4-j)*in.xMax/in.n;
        hybrid.y[j] = hybrid.y[j] + unity*2.0*in.gammaS;
    }
*/




    string name = "hybrid";
    hybrid.exportAWTasFUN(name, 10, 20, 1.0, in.output_mode);


    // 0. step: setting the initial GTd, GTu
    lorentz_propagator(in, GTd);
    insert_dynamic_hybridization(in, GTd, hybrid);

    lorentz_propagator(in, GTu);
    insert_dynamic_hybridization(in, GTu, hybrid);

    name = "GTd";
    GTd.exportAWTasFUN(name, 10, 10, 1.0, in.output_mode);

    name = "GTu";
    GTu.exportAWTasFUN(name, 10, 10, 1.0, in.output_mode);



}
