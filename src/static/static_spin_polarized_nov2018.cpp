#include "../physics.h"
#include "../convert.h"
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



void welcome_screen(double U, double h, double mu, double kT, string dirname)
{
    // on screen output of starting TIME and basic model parameters
    cout << __TIME__ << ",   " << __DATE__ << endl;
    cout << "---------------------------------------" << endl;
    cout << "CALCULATIONS FOR U = " << U <<  endl;
    cout << "static vertex with four sign combinations" << endl;
    cout << "equations according to Janis from end of November 2018  " << endl;
    cout << "---------------------------------------" << endl;
    cout << "---------------------------------------" << endl;
    cout << "magnetic field set to h = " << h << endl;
    cout << "mu is set to "         << mu << endl;
    cout << "kT is set to "         << kT << endl;
    cout << "creating directory "        << dirname << endl;
    cout << "---------------------------------------" << endl;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// BOTH ud and du variants for mag field
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void calculateGGud(AWT & GdGu1pp, AWT & GdGu1pm, AWT & GdGu2pm, AWT & GdGu2mm,
                    AWT & GTd, AWT & GTu, AWT & FD, AWT & K3, AUX & help, int out)
{
    if(out == 1)    cout << "GG";


    complex<double> unity(0,1);
    double Pi = 3.14159265359;

    // BUBBLE GdGu1pp
    for(int j = 0;          j < GTu.n+1;    j++)   help.aux1.y[j] = FD.y[j] * GTu.y[j] / ( -2.0*Pi*unity );
    for(int j = GTu.n+1;    j < 3*GTu.n+4;  j++)   help.aux1.y[j] = 0;
    for(int j = 3*GTu.n+4;  j < 4*GTu.n+4;  j++)   help.aux1.y[j] = FD.y[j] * GTu.y[j] / ( -2.0*Pi*unity );


    for(int j = 0;          j < GTu.n+1;    j++)   help.aux2.y[j] = GTd.y[j];
    for(int j = GTu.n+1;    j < 3*GTu.n+4;  j++)   help.aux2.y[j] = 0;
    for(int j = 3*GTu.n+4;  j < 4*GTu.n+4;  j++)   help.aux2.y[j] = GTd.y[j];

    GdGu1pp.convolutionAWT(help.aux1, help.aux2, 1, -1, help.aux3, help.aux4);  // f[x] * g[x+a]
    //GdGu1pp.KrammersKronig(GdGu1pp, K3, help.aux1, help.aux2, help.aux3, help.aux4);

    // BUBBLE GdGu1pm
    for(int j = 0;          j < GTu.n+1;    j++)   help.aux1.y[j] = FD.y[j] * conj( GTu.y[j] ) / ( -2.0*Pi*unity );
    for(int j = GTu.n+1;    j < 3*GTu.n+4;  j++)   help.aux1.y[j] = 0;
    for(int j = 3*GTu.n+4;  j < 4*GTu.n+4;  j++)   help.aux1.y[j] = FD.y[j] * conj( GTu.y[j] ) / ( -2.0*Pi*unity );

    for(int j = 0;          j < GTu.n+1;    j++)   help.aux2.y[j] = GTd.y[j];
    for(int j = GTu.n+1;    j < 3*GTu.n+4;  j++)   help.aux2.y[j] = 0;
    for(int j = 3*GTu.n+4;  j < 4*GTu.n+4;  j++)   help.aux2.y[j] = GTd.y[j];

    GdGu1pm.convolutionAWT(help.aux1, help.aux2, 1, -1, help.aux3, help.aux4);  // f[x] * g[x+a]
    //GdGu1pm.KrammersKronig(GdGu1pm, K3, help.aux1, help.aux2, help.aux3, help.aux4);


    // BUBBLE GdGu2pm
    for(int j = 0;          j < GTu.n+1;    j++)   help.aux1.y[j] = FD.y[j] * GTd.y[j]  / ( -2.0*Pi*unity );
    for(int j = GTu.n+1;    j < 3*GTu.n+4;  j++)   help.aux1.y[j] = 0;
    for(int j = 3*GTu.n+4;  j < 4*GTu.n+4;  j++)   help.aux1.y[j] = FD.y[j] * GTd.y[j]  / ( -2.0*Pi*unity );

    for(int j = 0;          j < GTu.n+1;    j++)   help.aux2.y[j] = conj( GTu.y[j] );
    for(int j = GTu.n+1;    j < 3*GTu.n+4;  j++)   help.aux2.y[j] = 0;
    for(int j = 3*GTu.n+4;  j < 4*GTu.n+4;  j++)   help.aux2.y[j] = conj( GTu.y[j] );


    GdGu2pm.convolutionAWT(help.aux1, help.aux2, -1, -1, help.aux3, help.aux4);  // f[x] * g[-a+x]
    //GdGu2pm.KrammersKronig(GdGu2pm, K3, help.aux1, help.aux2, help.aux3, help.aux4);


    // BUBBLE GdGu2mm
    for(int j = 0;          j < GTu.n+1;    j++)   help.aux1.y[j] = FD.y[j] * conj( GTd.y[j] )  / ( -2.0*Pi*unity );
    for(int j = GTu.n+1;    j < 3*GTu.n+4;  j++)   help.aux1.y[j] = 0;
    for(int j = 3*GTu.n+4;  j < 4*GTu.n+4;  j++)   help.aux1.y[j] = FD.y[j] * conj( GTd.y[j] )  / ( -2.0*Pi*unity );

    for(int j = 0;          j < GTu.n+1;    j++)   help.aux2.y[j] = conj( GTu.y[j] );
    for(int j = GTu.n+1;    j < 3*GTu.n+4;  j++)   help.aux2.y[j] = 0;
    for(int j = 3*GTu.n+4;  j < 4*GTu.n+4;  j++)   help.aux2.y[j] = conj( GTu.y[j] );


    GdGu2mm.convolutionAWT(help.aux1, help.aux2, -1, -1, help.aux3, help.aux4);  // f[x] * g[-a+x]
    //GdGu2mm.KrammersKronig(GdGu2mm, K3, help.aux1, help.aux2, help.aux3, help.aux4);


    if(out == 1)    cout << ", ";
}


void calculateGGdu(AWT & GuGd1pp, AWT & GuGd1pm, AWT & GuGd2pm, AWT & GuGd2mm, AWT & GTu, AWT & GTd,
                 AWT & FD, AWT & K3, AWT & aux1, AWT & aux2, AWT & aux3, AWT & aux4, int out)
{
    if(out == 1)    cout << "GG";

    // this function is not checked for correct up down symmetries!!!
    complex<double> unity(0,1);
    double Pi = 3.14159265359;

    // BUBBLE GuGd1pp
    for(int j = 0;          j < GTd.n+1;    j++)   aux1.y[j] = FD.y[j] * GTd.y[j] / ( -2.0*Pi*unity );
    for(int j = GTd.n+1;    j < 3*GTd.n+4;  j++)   aux1.y[j] = 0;
    for(int j = 3*GTd.n+4;  j < 4*GTd.n+4;  j++)   aux1.y[j] = FD.y[j] * GTd.y[j] / ( -2.0*Pi*unity );


    for(int j = 0;          j < GTd.n+1;    j++)   aux2.y[j] = GTu.y[j];
    for(int j = GTd.n+1;    j < 3*GTd.n+4;  j++)   aux2.y[j] = 0;
    for(int j = 3*GTd.n+4;  j < 4*GTd.n+4;  j++)   aux2.y[j] = GTu.y[j];

    GuGd1pp.convolutionAWT(aux1, aux2, 1, -1, aux3, aux4);  // f[x] * g[x+a]
    //GuGd1pp.KrammersKronig(GuGd1pp, K3, aux1, aux2, aux3, aux4);

    // BUBBLE GuGd1pm
    for(int j = 0;          j < GTd.n+1;    j++)   aux1.y[j] = FD.y[j] * conj( GTd.y[j] ) / ( -2.0*Pi*unity );
    for(int j = GTd.n+1;    j < 3*GTd.n+4;  j++)   aux1.y[j] = 0;
    for(int j = 3*GTd.n+4;  j < 4*GTd.n+4;  j++)   aux1.y[j] = FD.y[j] * conj( GTd.y[j] ) / ( -2.0*Pi*unity );

    for(int j = 0;          j < GTd.n+1;    j++)   aux2.y[j] = GTu.y[j];
    for(int j = GTd.n+1;    j < 3*GTd.n+4;  j++)   aux2.y[j] = 0;
    for(int j = 3*GTd.n+4;  j < 4*GTd.n+4;  j++)   aux2.y[j] = GTu.y[j];

    GuGd1pm.convolutionAWT(aux1, aux2, 1, -1, aux3, aux4);  // f[x] * g[x+a]
    //GuGd1pm.KrammersKronig(GuGd1pm, K3, aux1, aux2, aux3, aux4);


    // BUBBLE GuGd2pm
    for(int j = 0;          j < GTd.n+1;    j++)   aux1.y[j] = FD.y[j] * GTu.y[j]  / ( -2.0*Pi*unity );
    for(int j = GTd.n+1;    j < 3*GTd.n+4;  j++)   aux1.y[j] = 0;
    for(int j = 3*GTd.n+4;  j < 4*GTd.n+4;  j++)   aux1.y[j] = FD.y[j] * GTu.y[j]  / ( -2.0*Pi*unity );

    for(int j = 0;          j < GTd.n+1;    j++)   aux2.y[j] = conj( GTd.y[j] );
    for(int j = GTd.n+1;    j < 3*GTd.n+4;  j++)   aux2.y[j] = 0;
    for(int j = 3*GTd.n+4;  j < 4*GTd.n+4;  j++)   aux2.y[j] = conj( GTd.y[j] );


    GuGd2pm.convolutionAWT(aux1, aux2, -1, -1, aux3, aux4);  // f[x] * g[-a+x]
    //GuGd2pm.KrammersKronig(GuGd2pm, K3, aux1, aux2, aux3, aux4);


    // BUBBLE GuGd2mm
    for(int j = 0;          j < GTd.n+1;    j++)   aux1.y[j] = FD.y[j] * conj( GTu.y[j] )  / ( -2.0*Pi*unity );
    for(int j = GTd.n+1;    j < 3*GTd.n+4;  j++)   aux1.y[j] = 0;
    for(int j = 3*GTd.n+4;  j < 4*GTd.n+4;  j++)   aux1.y[j] = FD.y[j] * conj( GTu.y[j] )  / ( -2.0*Pi*unity );

    for(int j = 0;          j < GTd.n+1;    j++)   aux2.y[j] = conj( GTd.y[j] );
    for(int j = GTd.n+1;    j < 3*GTd.n+4;  j++)   aux2.y[j] = 0;
    for(int j = 3*GTd.n+4;  j < 4*GTd.n+4;  j++)   aux2.y[j] = conj( GTd.y[j] );


    GuGd2mm.convolutionAWT(aux1, aux2, -1, -1, aux3, aux4);  // f[x] * g[-a+x]
    //GuGd2mm.KrammersKronig(GuGd2mm, K3, aux1, aux2, aux3, aux4);


    if(out == 1)    cout << ", ";
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// BOTH ud and du variants for mag field
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void calculate_Dud_static(AWT & Dud, AWT & GdGu1pp, AWT & GdGu1pm, AWT & GdGu2pm, AWT & GdGu2mm, matrix & Lud, int out)
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

void calculate_Ddu_static(AWT & Ddu, AWT & GuGd1pp, AWT & GuGd1pm, AWT & GuGd2pm, AWT & GuGd2mm, matrix & Ldu, int out)
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

void calculate_Ddu_sym(AWT & Ddu, AWT & Dud, int out)
{
    if(out == 1)    cout << "Ddu";

    // zero aux frequency
    Ddu.y[0] = Dud.y[0] ;
    // positive aux frequencies
    for(int j = 1;          j < Ddu.n+1;    j++)
    {
        Ddu.y[j] =  conj(   Dud.y[Dud.nn - j]   );
    }
    // zero padding of aux
    for(int j = Ddu.n+1;    j < 3*Ddu.n+4;  j++)
    {
        Ddu.y[j] = 0;
    }
    // negative aux frequencies
    for(int j = 3*Ddu.n+4;  j < 4*Ddu.n+4;  j++)
    {
        Ddu.y[j] = conj(   Dud.y[Dud.nn - j]   );
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// BOTH ud and du variants for mag field
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void calculate_Kud_static(matrix & Kud, AWT & GdGu1pp, AWT & GdGu1pm, AWT & GdGu2pm, AWT & GdGu2mm, matrix &Lud, int out)
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

void calculate_Kdu_static(matrix & Kdu, AWT & GuGd1pp, AWT & GuGd1pm, AWT & GuGd2pm, AWT & GuGd2mm, matrix &Ldu, int out)
{
    // this function is not checked for correct up down symmetries!!!

    if(out == 1)    cout << "Kdu";

    Kdu.pp = - Ldu.pp * Ldu.pp * GuGd1pp.y[0] - Ldu.mp * conj(Ldu.mp) * conj(GuGd1pp.y[0])
             - Ldu.pp * ( Ldu.pp * conj(Ldu.pp) - Ldu.mp * conj(Ldu.mp)  ) * GuGd1pp.y[0] * conj(GuGd1pp.y[0]);

    Kdu.mp = - Ldu.mp * Ldu.pp * GuGd1pp.y[0] - Ldu.mp * conj(Ldu.pp) * conj(GuGd1pp.y[0])
             - Ldu.mp * ( Ldu.pp * conj(Ldu.pp) - Ldu.mp * conj(Ldu.mp)  ) * GuGd1pp.y[0] * conj(GuGd1pp.y[0]);

    Kdu.pm = conj(Kdu.mp);

    Kdu.mm = conj(Kdu.pp);

    if(out == 1)    cout << ", ";
}

void calculate_Kdu_sym(matrix & Kdu, matrix & Kud, int out)
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
void calculate_FDud(matrix & FDud, AWT & GTd, AWT & GTu, AWT & Dud, AWT & FD, AWT & aux, int out)
{
    double Pi = 3.14159265359;

    if(out == 1)    cout << "FDud";

    // FDud(+,+)
    // zero aux frequency
        aux.y[0] = FD.y[0] * GTu.y[0] * (-1.0) * imag( GTd.y[0] ) / conj( Dud.y[0] )
                 + FD.y[0] * GTd.y[0] * (-1.0) * imag( GTu.y[0] ) / Dud.y[0] ;
    // positive aux frequencies
    for(int j = 1;          j < GTu.n+1;    j++)
    {
        aux.y[j] = FD.y[j] * GTu.y[    j   ] * (-1.0) * imag( GTd.y[ aux.nn - j] ) / conj( Dud.y[aux.nn - j] )
                 + FD.y[j] * GTd.y[aux.nn-j] * (+1.0) * imag( GTu.y[     j     ] ) / (   Dud.y[aux.nn - j]   );
    }
    // zero padding of aux
    for(int j = GTu.n+1;    j < 3*GTu.n+4;  j++)
    {
        aux.y[j] = 0;
    }
    // negative aux frequencies
    for(int j = 3*GTu.n+4;  j < 4*GTu.n+4;  j++)
    {
        aux.y[j] = FD.y[j] * GTu.y[    j   ] * (-1.0) * imag( GTd.y[ aux.nn - j] ) / conj( Dud.y[aux.nn - j] )
                 + FD.y[j] * GTd.y[aux.nn-j] * (+1.0) * imag( GTu.y[     j     ] ) / (   Dud.y[aux.nn - j]   );
    }

    FDud.pp = Simpson(aux) / ( - Pi);

    // FDud(-,+)
    // zero aux frequency
        aux.y[0] = FD.y[0] * imag( GTu.y[0] * conj( GTd.y[     0     ] ) ) / Dud.y[    0     ];
    // positive aux frequencies
    for(int j = 1;          j < GTu.n+1;    j++)    aux.y[j] = FD.y[j] * imag( GTu.y[j] * conj( GTd.y[ aux.nn - j] ) ) / Dud.y[aux.nn - j];
    // zero padding of aux
    for(int j = GTu.n+1;    j < 3*GTu.n+4;  j++)    aux.y[j] = 0;
    // negative aux frequencies
    for(int j = 3*GTu.n+4;  j < 4*GTu.n+4;  j++)    aux.y[j] = FD.y[j] * imag( GTu.y[j] * conj( GTd.y[ aux.nn - j] ) ) / Dud.y[aux.nn - j];

    FDud.mp = Simpson(aux) / ( - Pi);

    // FDud(-,+)
    FDud.pm = conj(FDud.mp);

    // FDud(-,+)
    FDud.mm = conj(FDud.pp);

    if(out == 1)    cout << ", ";
}

void calculate_FDdu(matrix & FDdu, AWT & GTd, AWT & GTu, AWT & Ddu, AWT & FD, AWT & aux, int out)
{
    // this function is not checked for correct up down symmetries!!!
    double Pi = 3.14159265359;

    if(out == 1)    cout << "FDdu";

    // FDud(+,+)
    // zero aux frequency
        aux.y[0] = FD.y[0] * GTd.y[0] * (-1.0) * imag( GTu.y[0] ) / conj( Ddu.y[0] )
                 + FD.y[0] * GTu.y[0] * (-1.0) * imag( GTd.y[0] ) / Ddu.y[0] ;
    // positive aux frequencies
    for(int j = 1;          j < GTu.n+1;    j++)
    {
        aux.y[j] = FD.y[j] * GTd.y[    j   ] * (-1.0) * imag( GTu.y[ aux.nn - j] ) / conj( Ddu.y[aux.nn - j] )
                 + FD.y[j] * GTu.y[aux.nn-j] * (+1.0) * imag( GTd.y[     j     ] ) / (   Ddu.y[aux.nn - j]   );
    }
    // zero padding of aux
    for(int j = GTu.n+1;    j < 3*GTu.n+4;  j++)
    {
        aux.y[j] = 0;
    }
    // negative aux frequencies
    for(int j = 3*GTu.n+4;  j < 4*GTu.n+4;  j++)
    {
        aux.y[j] = FD.y[j] * GTd.y[    j   ] * (-1.0) * imag( GTu.y[ aux.nn - j] ) / conj( Ddu.y[aux.nn - j] )
                 + FD.y[j] * GTu.y[aux.nn-j] * (+1.0) * imag( GTd.y[     j     ] ) / (   Ddu.y[aux.nn - j]   );
    }

    FDdu.pp = Simpson(aux) / ( - Pi);

    // FDud(-,+)
    // zero aux frequency
        aux.y[0] = FD.y[0] * imag( GTd.y[0] * conj( GTu.y[     0     ] ) ) / Ddu.y[    0     ];
    // positive aux frequencies
    for(int j = 1;          j < GTd.n+1;    j++)    aux.y[j] = FD.y[j] * imag( GTd.y[j] * conj( GTu.y[ aux.nn - j] ) ) / Ddu.y[aux.nn - j];
    // zero padding of aux
    for(int j = GTu.n+1;    j < 3*GTu.n+4;  j++)    aux.y[j] = 0;
    // negative aux frequencies
    for(int j = 3*GTu.n+4;  j < 4*GTu.n+4;  j++)    aux.y[j] = FD.y[j] * imag( GTd.y[j] * conj( GTu.y[ aux.nn - j] ) ) / Ddu.y[aux.nn - j];

    FDdu.mp = Simpson(aux) / ( - Pi);

    // FDud(-,+)
    FDdu.pm = conj(FDdu.mp);

    // FDud(-,+)
    FDdu.mm = conj(FDdu.pp);

    if(out == 1)    cout << ", ";
}

void calculate_FDdu_sym(matrix & FDdu, matrix & FDud, int out)
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
void calculate_BDud(matrix & BDud, AWT & GTd, AWT & GTu, AWT & Dud, AWT & BE, AWT & aux, int out)
{
    if(out == 1)    cout << "BDud";
    double Pi = 3.14159265359;

    // BDud(+,+)
    // positive aux frequencies
    for(int j = 1;          j < GTu.n+1;    j++)    aux.y[j] = BE.y[j] * GTu.y[j] * GTd.y[aux.nn - j] * imag( 1.0 / conj( Dud.y[aux.nn - j] )  );
    // zero padding of aux
    for(int j = GTu.n+1;    j < 3*GTu.n+4;  j++)    aux.y[j] = 0;
    // negative aux frequencies
    for(int j = 3*GTu.n+4;  j < 4*GTu.n+4;  j++)    aux.y[j] = BE.y[j] * GTu.y[j] * GTd.y[aux.nn - j] * imag( 1.0 / conj( Dud.y[aux.nn - j] )  );
    // zero aux frequency
        aux.y[0] = 0.5 * (aux.y[4*GTu.n+3] + aux.y[1] );


    BDud.pp = Simpson(aux) / ( Pi);
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // BDud(-,+)
    // zero aux frequency
        aux.y[0] =  BE.y[0] * GTu.y[0] * conj( GTd.y[0] ) * imag( 1.0 / conj( Dud.y[0] )  );
    // positive aux frequencies
    for(int j = 1;          j < GTu.n+1;    j++)    aux.y[j] = BE.y[j] * GTu.y[j] * conj( GTd.y[aux.nn - j] ) * imag( 1.0 / conj( Dud.y[aux.nn - j] )  );
    // zero padding of aux
    for(int j = GTu.n+1;    j < 3*GTu.n+4;  j++)    aux.y[j] = 0;
    // negative aux frequencies
    for(int j = 3*GTu.n+4;  j < 4*GTu.n+4;  j++)    aux.y[j] = BE.y[j] * GTu.y[j] * conj( GTd.y[aux.nn - j] ) * imag( 1.0 / conj( Dud.y[aux.nn - j] )  );
    // zero aux frequency
        aux.y[0] = 0.5 * (aux.y[4*GTu.n+3] + aux.y[1] );


    BDud.mp = Simpson(aux) / ( Pi);
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // BDud(-,+)
    BDud.pm = conj(BDud.mp);

    // BDud(-,+)
    BDud.mm = conj(BDud.pp);
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


    if(out == 1)    cout << ", ";
}


void calculate_BDdu(matrix & BDdu, AWT & GTd, AWT & GTu, AWT & Ddu, AWT & BE, AWT & aux, int out)
{
    // this function is not checked for correct up down symmetries!!!
    if(out == 1)    cout << "BDud";
    double Pi = 3.14159265359;

    // BDud(+,+)
    // positive aux frequencies
    for(int j = 1;          j < GTu.n+1;    j++)    aux.y[j] = BE.y[j] * GTd.y[j] * GTu.y[aux.nn - j] * imag( 1.0 / conj( Ddu.y[aux.nn - j] )  );
    // zero padding of aux
    for(int j = GTu.n+1;    j < 3*GTu.n+4;  j++)    aux.y[j] = 0;
    // negative aux frequencies
    for(int j = 3*GTu.n+4;  j < 4*GTu.n+4;  j++)    aux.y[j] = BE.y[j] * GTd.y[j] * GTu.y[aux.nn - j] * imag( 1.0 / conj( Ddu.y[aux.nn - j] )  );
    // zero aux frequency
        aux.y[0] = 0.5 * (aux.y[4*GTu.n+3] + aux.y[1] );


    BDdu.pp = Simpson(aux) / ( Pi);
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // BDud(-,+)
    // zero aux frequency
        aux.y[0] =  BE.y[0] * GTu.y[0] * conj( GTu.y[0] ) * imag( 1.0 / conj( Ddu.y[0] )  );
    // positive aux frequencies
    for(int j = 1;          j < GTu.n+1;    j++)    aux.y[j] = BE.y[j] * GTd.y[j] * conj( GTu.y[aux.nn - j] ) * imag( 1.0 / conj( Ddu.y[aux.nn - j] )  );
    // zero padding of aux
    for(int j = GTu.n+1;    j < 3*GTu.n+4;  j++)    aux.y[j] = 0;
    // negative aux frequencies
    for(int j = 3*GTu.n+4;  j < 4*GTu.n+4;  j++)    aux.y[j] = BE.y[j] * GTd.y[j] * conj( GTu.y[aux.nn - j] ) * imag( 1.0 / conj( Ddu.y[aux.nn - j] )  );
    // zero aux frequency
        aux.y[0] = 0.5 * (aux.y[4*GTu.n+3] + aux.y[1] );


    BDdu.mp = Simpson(aux) / ( Pi);
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // BDud(-,+)
    BDdu.pm = conj(BDdu.mp);

    // BDud(-,+)
    BDdu.mm = conj(BDdu.pp);
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


    if(out == 1)    cout << ", ";
}

void calculate_BDdu_sym(matrix & BDdu, matrix & BDud, int out)
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
void calculate_Lud(matrix & Lud, matrix & Kud, matrix & FDud, matrix & BDud, double U, int out)
{
    if(out == 1)    cout << "Lud";

    // xx, attention real is a cheater stuff here!
    Lud.pp = U / ( 1.0 + Kud.pp * real( BDud.pp + FDud.pp)  );
    Lud.mp = U / ( 1.0 + Kud.pm * real( BDud.pm + FDud.pm)  );

    Lud.pm = conj(Lud.mp);
    Lud.mm = conj(Lud.pp);

    if(out == 1)    cout << ", ";
}



void calculate_Ldu(matrix & Ldu, matrix & Kdu, matrix & FDdu, matrix & BDdu, double U, int out)
{
    if(out == 1)    cout << "Lud";

    // xx, attention real is a cheater stuff here!
    Ldu.pp = U / ( 1.0 + Kdu.pp * real( BDdu.pp + FDdu.pp)  );
    Ldu.mp = U / ( 1.0 + Kdu.pm * real( BDdu.pm + FDdu.pm)  );

    Ldu.pm = conj(Ldu.mp);
    Ldu.mm = conj(Ldu.pp);

    if(out == 1)    cout << ", ";
}


void calculate_Ldu_sym(matrix & Ldu, matrix & Lud, int out)
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
void calculate_fGG(complex<double> & fGG, AWT & GTd, AWT & GTu, AWT & FD, AWT & aux)
{
    double Pi = 3.14159265359;
    complex<double> unity(0,1);

    // positive aux frequencies
    for(int j = 0;          j < GTu.n+1;    j++)    aux.y[j] = FD.y[j] * GTd.y[j] * GTu.y[j] ;
    // zero padding of aux
    for(int j = GTu.n+1;    j < 3*GTu.n+4;  j++)    aux.y[j] = 0;
    // negative aux frequencies
    for(int j = 3*GTu.n+4;  j < 4*GTu.n+4;  j++)    aux.y[j] = FD.y[j] * GTd.y[j] * GTu.y[j] ;

    fGG = Simpson(aux) / ( - 2.0 * Pi * unity);

}


void set_lambda_iter(string mode, string dirname, double U, info in,
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

void lambda_iter_output(ofstream & log_output, int iter_lambda, int iter_sigma, complex<double> & aSTu, complex<double> & aSTd,
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

void lambda_iter_console(int out, matrix Lud, matrix Ldu, matrix Kud, matrix Kdu, matrix BDud, matrix BDdu, matrix FDud, matrix FDdu)
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

void calculate_aST(complex<double> & aSTu, complex<double> & aSTd,
                    matrix & Lud, matrix & Ldu,
                    AWT & GTd, AWT & GTu, AWT & FD, AWT & aux, int out)
{


    complex<double> unity(0,1);
    double Pi = 3.14159265359;

    complex<double> Gup, Gdp, Gum, Gdm;

    // GT up plus integral
    for(int j = 0;          j < GTu.n+1;    j++)   aux.y[j] = FD.y[j] * GTu.y[j];
    for(int j = GTu.n+1;    j < 3*GTu.n+4;  j++)   aux.y[j] = 0;
    for(int j = 3*GTu.n+4;  j < 4*GTu.n+4;  j++)   aux.y[j] = FD.y[j] * GTu.y[j];

    Gup = Simpson(aux) / ( -2.0 * unity * Pi );

    // GT down plus integral
    for(int j = 0;          j < GTu.n+1;    j++)   aux.y[j] = FD.y[j] * GTd.y[j];
    for(int j = GTu.n+1;    j < 3*GTu.n+4;  j++)   aux.y[j] = 0;
    for(int j = 3*GTu.n+4;  j < 4*GTu.n+4;  j++)   aux.y[j] = FD.y[j] * GTd.y[j];

    Gdp = Simpson(aux) / ( -2.0 * unity * Pi );

    // GT up minus integral
    for(int j = 0;          j < GTu.n+1;    j++)   aux.y[j] = FD.y[j] * conj( GTu.y[j] );
    for(int j = GTu.n+1;    j < 3*GTu.n+4;  j++)   aux.y[j] = 0;
    for(int j = 3*GTu.n+4;  j < 4*GTu.n+4;  j++)   aux.y[j] = FD.y[j] * conj( GTu.y[j] );

    Gum = Simpson(aux) / ( -2.0 * unity * Pi );

    // GT down minus integral
    for(int j = 0;          j < GTu.n+1;    j++)   aux.y[j] = FD.y[j] * conj( GTd.y[j] );
    for(int j = GTu.n+1;    j < 3*GTu.n+4;  j++)   aux.y[j] = 0;
    for(int j = 3*GTu.n+4;  j < 4*GTu.n+4;  j++)   aux.y[j] = FD.y[j] * conj( GTd.y[j] );

    Gdm = Simpson(aux) / ( -2.0 * unity * Pi );


    aSTu = - 0.25 * (Lud.pp + Ldu.pp) * (Gup - Gdp)
           - 0.25 * (Lud.pm + Ldu.pm) * (Gup - Gdp);
    aSTd = - 0.25 * (Lud.pp + Ldu.pp) * (Gdm - Gum)
           - 0.25 * (Lud.pm + Ldu.pm) * (Gdm - Gum);

    cout << "anomalous = " << aSTu << " " << aSTd << endl;

}

void calculate_nT_mT(double & nT, double & mT, double mu, double U, AWT & GTu, AWT & GTd, AWT & FD, AWT & aux)
{
    double Pi = 3.14159265359;

    for(int i=0;          i<GTu.n+1;   i++ )       aux.y[i] = FD.y[i] * imag( GTu.y[i] );
    for(int i=GTu.n+1;    i<3*GTu.n+3; i++ )       aux.y[i] = 0;
    for(int i=3*GTu.n+4;  i<4*GTu.n+4; i++ )       aux.y[i] = FD.y[i] * imag( GTu.y[i] );



    double nTup = - Simpson_Re(aux) / Pi;


    for(int i=0;          i<GTu.n+1;   i++ )       aux.y[i] = FD.y[i] * imag( GTd.y[i] );
    for(int i=GTu.n+1;    i<3*GTu.n+3; i++ )       aux.y[i] = 0;
    for(int i=3*GTu.n+4;  i<4*GTu.n+4; i++ )       aux.y[i] = FD.y[i] * imag( GTd.y[i] );

    double nTdown = - Simpson_Re(aux) / Pi;

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

void calculate_nS_mS(double & nS, double & mS, AWT & GSu, AWT & GSd, AWT & FD, AWT & aux)
{
    double Pi = 3.14159265359;


    for(int i=0;          i<GSu.n+1;   i++ )       aux.y[i] = FD.y[i] * imag( GSu.y[i] );
    for(int i=GSu.n+1;    i<3*GSu.n+3; i++ )       aux.y[i] = 0;
    for(int i=3*GSu.n+4;  i<4*GSu.n+4; i++ )       aux.y[i] = FD.y[i] * imag( GSu.y[i] );

    nS = - Simpson_Re(aux) / Pi;
    mS = - Simpson_Re(aux) / Pi;
    cout << "nS up = " << nS;

    for(int i=0;          i<GSu.n+1;   i++ )       aux.y[i] = FD.y[i] * imag( GSd.y[i] );
    for(int i=GSu.n+1;    i<3*GSu.n+3; i++ )       aux.y[i] = 0;
    for(int i=3*GSu.n+4;  i<4*GSu.n+4; i++ )       aux.y[i] = FD.y[i] * imag( GSd.y[i] );

    nS = nS - Simpson_Re(aux) / Pi;
    mS = mS + Simpson_Re(aux) / Pi;
    cout << "  nS up + down = " << nS << endl;
    cout << "  mS = " << mS << endl;

    // add tail
    nS = nS - GSu.xMax * imag( GSu.y[3*GSu.n+4] ) / (Pi);
    //cout << "nS corr: " << - GSu.xMax * imag( GSu.y[3*GSu.n+4] ) / (Pi)  << endl;

}


void calculate_sigma_stat(AWT & SigmaUp, AWT & SigmaDown, info in, double U, double kT,
                          AWT & Dud, AWT & Ddu, matrix & Lud, matrix & Ldu,
                          AWT & GTd, AWT & GTu, AWT & FD, AWT & BE, AWT & K3, string dirname, int print,
                          AWT & aux1, AWT & aux2, AWT & aux3, AWT & aux4, AWT & aux5, AWT & aux6, AWT & aux7)
{
    complex<double> unity(0,1);
    double Pi = 3.14159265359;

    AWT ThetaPP, ThetaPM;
    ThetaPP.initializeAWT(GTu.n, GTu.xMax, kT);
    ThetaPM.initializeAWT(GTu.n, GTu.xMax, kT);

    double dx = GTu.xMax / GTu.n;
    cout << "   Self-energy Up calculation" << endl;

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // GAMMA MATRIX
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    matrix gamma;

    // positive aux frequencies
    for(int j = 0;          j < GTu.n+1;    j++)    aux1.y[j] = FD.y[j] * GTu.y[j] * GTu.y[j];
    for(int j = GTu.n+1;    j < 3*GTu.n+4;  j++)    aux1.y[j] = 0;
    for(int j = 3*GTu.n+4;  j < 4*GTu.n+4;  j++)    aux1.y[j] = FD.y[j] * GTu.y[j] * GTu.y[j];

    gamma.pp = Lud.pp + ( Lud.pp * conj( Lud.pp) - Lud.mp * conj( Lud.mp) ) * conj( Simpson(aux1) / ( -2.0*Pi*unity));
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
    for(int j = 0;          j < GTu.n+1;    j++)    aux4.y[j] = FD.y[j] * imag( GTu.y[j] ) / Pi;
    for(int j = GTu.n+1;    j < 3*GTu.n+4;  j++)    aux4.y[j] = 0;
    for(int j = 3*GTu.n+4;  j < 4*GTu.n+4;  j++)    aux4.y[j] = FD.y[j] * imag( GTu.y[j] ) / Pi;
    // second integrand
    for(int j = 0;          j < GTu.n+1;    j++)    aux5.y[j] = GTd.y[j];
    for(int j = GTu.n+1;    j < 3*GTu.n+4;  j++)    aux5.y[j] = 0;
    for(int j = 3*GTu.n+4;  j < 4*GTu.n+4;  j++)    aux5.y[j] = GTd.y[j];

    aux1.convolutionAWT(aux4, aux5, 1, -1, aux6, aux7);
    aux1.KrammersKronig(aux1, K3, aux4, aux5, aux6, aux7);

    string name;
    if(print == 1)
    {
        name = dirname + "/ThetaPP1up";
        aux1.exportAWTasFUN(name, in.display, in.xMax, 1, in.output_mode);
    }


    // SECOND PART of THETA PLUS
    // first integrand
    for(int j = 0;          j < GTu.n+1;    j++)    aux4.y[j] = FD.y[j] * GTd.y[j] / (2 * Pi * unity );
    for(int j = GTu.n+1;    j < 3*GTu.n+4;  j++)    aux4.y[j] = 0;
    for(int j = 3*GTu.n+4;  j < 4*GTu.n+4;  j++)    aux4.y[j] = FD.y[j] * GTd.y[j] / (2 * Pi * unity );
    // second integrand
    for(int j = 0;          j < GTu.n+1;    j++)    aux5.y[j] = conj( GTu.y[j] );
    for(int j = GTu.n+1;    j < 3*GTu.n+4;  j++)    aux5.y[j] = 0;
    for(int j = 3*GTu.n+4;  j < 4*GTu.n+4;  j++)    aux5.y[j] = conj( GTu.y[j] );

    aux2.convolutionAWT(aux4, aux5, -1, -1, aux6, aux7);
    aux2.KrammersKronig(aux2, K3, aux4, aux5, aux6, aux7);

    if(print == 1)
    {
        name = dirname + "/ThetaPP2up";
        aux2.exportAWTasFUN(name, in.display, in.xMax, 1, in.output_mode);
    }

    // THIRD PART of THETA PLUS
    // first integrand
    for(int j = 0;          j < GTu.n+1;    j++)    aux4.y[j] = FD.y[j] * conj( GTd.y[j] ) / (2 * Pi * unity );
    for(int j = GTu.n+1;    j < 3*GTu.n+4;  j++)    aux4.y[j] = 0;
    for(int j = 3*GTu.n+4;  j < 4*GTu.n+4;  j++)    aux4.y[j] = FD.y[j] * conj( GTd.y[j] ) / (2 * Pi * unity );
    // second integrand
    for(int j = 0;          j < GTu.n+1;    j++)    aux5.y[j] = conj( GTu.y[j] );
    for(int j = GTu.n+1;    j < 3*GTu.n+4;  j++)    aux5.y[j] = 0;
    for(int j = 3*GTu.n+4;  j < 4*GTu.n+4;  j++)    aux5.y[j] = conj( GTu.y[j] );

    aux3.convolutionAWT(aux4, aux5, -1, -1, aux6, aux7);
    aux3.KrammersKronig(aux3, K3, aux4, aux5, aux6, aux7);

    if(print == 1)
    {
        name = dirname + "/ThetaPP3up";
        aux3.exportAWTasFUN(name, in.display, in.xMax, 1, in.output_mode);
    }


    // THETA PLUS
    for(int j = 0;          j < GTu.n+1;    j++)    ThetaPP.y[j] = -gamma.pp * aux1.y[j] - gamma.pp * aux2.y[j] + gamma.mp * aux3.y[j];
    for(int j = GTu.n+1;    j < 3*GTu.n+4;  j++)    ThetaPP.y[j] = 0;
    for(int j = 3*GTu.n+4;  j < 4*GTu.n+4;  j++)    ThetaPP.y[j] = -gamma.pp * aux1.y[j] - gamma.pp * aux2.y[j] + gamma.mp * aux3.y[j];


    if(print == 1)
    {
        name = dirname + "/ThetaPPup";
        ThetaPP.exportAWTasFUN(name, in.display, in.xMax, 1, in.output_mode);
    }

    // FIRST PART of THETA MINUS
    // first integrand
    for(int j = 0;          j < GTu.n+1;    j++)    aux4.y[j] = FD.y[j] * imag( GTu.y[j] ) / Pi;
    for(int j = GTu.n+1;    j < 3*GTu.n+4;  j++)    aux4.y[j] = 0;
    for(int j = 3*GTu.n+4;  j < 4*GTu.n+4;  j++)    aux4.y[j] = FD.y[j] * imag( GTu.y[j] ) / Pi;
    // second integrand
    for(int j = 0;          j < GTu.n+1;    j++)    aux5.y[j] = conj( GTd.y[j] );
    for(int j = GTu.n+1;    j < 3*GTu.n+4;  j++)    aux5.y[j] = 0;
    for(int j = 3*GTu.n+4;  j < 4*GTu.n+4;  j++)    aux5.y[j] = conj( GTd.y[j] );

    aux1.convolutionAWT(aux4, aux5, 1, -1, aux6, aux7);
    aux1.KrammersKronigDown(aux1, K3, aux4, aux5, aux6, aux7);


    if(print == 1)
    {
        name = dirname + "/ThetaPM1up";
        aux1.exportAWTasFUN(name, in.display, in.xMax, 1, in.output_mode);
    }


    // SECOND PART of THETA MINUS
    // first integrand
    for(int j = 0;          j < GTu.n+1;    j++)    aux4.y[j] = FD.y[j] * GTd.y[j] / (2 * Pi * unity );
    for(int j = GTu.n+1;    j < 3*GTu.n+4;  j++)    aux4.y[j] = 0;
    for(int j = 3*GTu.n+4;  j < 4*GTu.n+4;  j++)    aux4.y[j] = FD.y[j] * GTd.y[j] / (2 * Pi * unity );
    // second integrand is just the Green's function up

    aux2.convolutionAWT(aux4, GTu, -1, -1, aux6, aux7);
    aux2.KrammersKronigDown(aux2, K3, aux4, aux5, aux6, aux7);


    if(print == 1)
    {
        name = dirname + "/ThetaPM2up";
        aux2.exportAWTasFUN(name, in.display, in.xMax, 1, in.output_mode);
    }


    // THIRD PART of THETA MINS
    // first integrand
    for(int j = 0;          j < GTu.n+1;    j++)    aux4.y[j] = FD.y[j] * conj( GTd.y[j] ) / (2 * Pi * unity );
    for(int j = GTu.n+1;    j < 3*GTu.n+4;  j++)    aux4.y[j] = 0;
    for(int j = 3*GTu.n+4;  j < 4*GTu.n+4;  j++)    aux4.y[j] = FD.y[j] * conj( GTd.y[j] ) / (2 * Pi * unity );
    // second integrand is just the Green's function up

    aux3.convolutionAWT(aux4, GTu, -1, -1, aux6, aux7);
    aux3.KrammersKronigDown(aux3, K3, aux4, aux5, aux6, aux7);


    if(print == 1)
    {
        name = dirname + "/ThetaPM3up";
        aux3.exportAWTasFUN(name, in.display, in.xMax, 1, in.output_mode);
    }

    // THETA PLUS
    for(int j = 0;          j < GTu.n+1;    j++)    ThetaPM.y[j] = -gamma.mp * aux1.y[j] - gamma.pp * aux2.y[j] + gamma.mp * aux3.y[j];
    for(int j = GTu.n+1;    j < 3*GTu.n+4;  j++)    ThetaPM.y[j] = 0;
    for(int j = 3*GTu.n+4;  j < 4*GTu.n+4;  j++)    ThetaPM.y[j] = -gamma.mp * aux1.y[j] - gamma.pp * aux2.y[j] + gamma.mp * aux3.y[j];


    if(print == 1)
    {
        name = dirname + "/ThetaPMup";
        ThetaPM.exportAWTasFUN(name, in.display, in.xMax, 1, in.output_mode);
    }

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //  SIGMA
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // FD integrals part
    // first integrand
    for(int j = 0;          j < GTu.n+1;    j++)    aux4.y[j] = U * FD.y[j] * imag( GTd.y[j] ) / Pi;
    for(int j = GTu.n+1;    j < 3*GTu.n+4;  j++)    aux4.y[j] = 0;
    for(int j = 3*GTu.n+4;  j < 4*GTu.n+4;  j++)    aux4.y[j] = U * FD.y[j] * imag( GTd.y[j] ) / Pi;
    // second integrand
    for(int j = 0;          j < GTu.n+1;    j++)    aux5.y[j] = ThetaPM.y[j] / ( conj( Dud.y[j] ) );
    for(int j = GTu.n+1;    j < 3*GTu.n+4;  j++)    aux5.y[j] = 0;
    for(int j = 3*GTu.n+4;  j < 4*GTu.n+4;  j++)    aux5.y[j] = ThetaPM.y[j] / ( conj( Dud.y[j] ) );

    aux5.deleteReal(aux5);

    aux1.convolutionAWT(aux4, aux5, -1, -1, aux6, aux7);
    aux1.KrammersKronig(aux1, K3, aux4, aux5, aux6, aux7);

    if(print == 1)
    {
        name = dirname + "/FDintUp";
        aux1.exportAWTasFUN(name, in.display, in.xMax, 1, in.output_mode);
    }

    for(int j = 0;          j < GTu.n+1;    j++)    aux6.y[j] = aux4.y[j] * aux5.y[j];
    for(int j = GTu.n+1;    j < 3*GTu.n+4;  j++)    aux6.y[j] = 0;
    for(int j = 3*GTu.n+4;  j < 4*GTu.n+4;  j++)    aux6.y[j] = aux4.y[j] * aux5.y[j];
    SigmaUp.y[0] = Simpson_Re(aux6) + unity * Simpson_Im(aux6);


    for(int j = 0;          j < GTu.n+1;    j++)    aux4.y[j] = U * BE.y[j] * imag( ThetaPP.y[j]/Dud.y[j]  ) / Pi;
    for(int j = GTu.n+1;    j < 3*GTu.n+4;  j++)    aux4.y[j] = 0;
    for(int j = 3*GTu.n+4;  j < 4*GTu.n+4;  j++)    aux4.y[j] = U * BE.y[j] * imag( ThetaPP.y[j]/Dud.y[j]  ) / Pi;

    // second integrand
    for(int j = 0;          j < GTu.n+1;    j++)    aux5.y[j] = GTd.y[j];
    for(int j = GTu.n+1;    j < 3*GTu.n+4;  j++)    aux5.y[j] = 0;
    for(int j = 3*GTu.n+4;  j < 4*GTu.n+4;  j++)    aux5.y[j] = GTd.y[j];

    aux5.deleteReal(aux5);

    aux2.convolutionAWT(aux4, aux5, 1, -1, aux6, aux7);
    aux2.KrammersKronig(aux2, K3, aux4, aux5, aux6, aux7);


    if(print == 1)
    {
        name = dirname + "/BEintUp";
        aux2.exportAWTasFUN(name, in.display, in.xMax, 1, in.output_mode);
    }

    for(int j = 0;          j < GTu.n+1;    j++)    aux6.y[j] = aux4.y[j] * aux5.y[j];
    for(int j = GTu.n+1;    j < 3*GTu.n+4;  j++)    aux6.y[j] = 0;
    for(int j = 3*GTu.n+4;  j < 4*GTu.n+4;  j++)    aux6.y[j] = aux4.y[j] * aux5.y[j];

    SigmaUp.y[0] = SigmaUp.y[0] - Simpson_Re(aux6) - unity * Simpson_Im(aux6);


    for(int j=0; j < 4*GTu.n + 4; j++ ) SigmaUp.y[j] = aux1.y[j] - aux2.y[j];


    // export
    if(print == 1)
    {
        name = dirname + "/SSu";
        SigmaUp.exportAWTasFUN(name, in.display, in.xMax, 1, in.output_mode);
    }

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //  THETA CALCULATIONS Sigma DOWN
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // FIRST PART of THETA PLUS down
    // first integrand
    for(int j = 0;          j < GTu.n+1;    j++)    aux4.y[j] = FD.y[j] * imag( GTd.y[j] ) / Pi;
    for(int j = GTu.n+1;    j < 3*GTu.n+4;  j++)    aux4.y[j] = 0;
    for(int j = 3*GTu.n+4;  j < 4*GTu.n+4;  j++)    aux4.y[j] = FD.y[j] * imag( GTd.y[j] ) / Pi;
    // second integrand
    for(int j = 0;          j < GTu.n+1;    j++)    aux5.y[j] = GTu.y[j];
    for(int j = GTu.n+1;    j < 3*GTu.n+4;  j++)    aux5.y[j] = 0;
    for(int j = 3*GTu.n+4;  j < 4*GTu.n+4;  j++)    aux5.y[j] = GTu.y[j];

    aux1.convolutionAWT(aux4, aux5, 1, -1, aux6, aux7);
    aux1.KrammersKronig(aux1, K3, aux4, aux5, aux6, aux7);


    if(print == 1)
    {
        name = dirname + "/ThetaPP1down";
        aux1.exportAWTasFUN(name, in.display, in.xMax, 1, in.output_mode);
    }

    // SECOND PART of THETA PLUS down
    // first integrand
    for(int j = 0;          j < GTu.n+1;    j++)    aux4.y[j] = FD.y[j] * GTu.y[j] / (2 * Pi * unity );
    for(int j = GTu.n+1;    j < 3*GTu.n+4;  j++)    aux4.y[j] = 0;
    for(int j = 3*GTu.n+4;  j < 4*GTu.n+4;  j++)    aux4.y[j] = FD.y[j] * GTu.y[j] / (2 * Pi * unity );
    // second integrand
    for(int j = 0;          j < GTu.n+1;    j++)    aux5.y[j] = conj( GTd.y[j] );
    for(int j = GTu.n+1;    j < 3*GTu.n+4;  j++)    aux5.y[j] = 0;
    for(int j = 3*GTu.n+4;  j < 4*GTu.n+4;  j++)    aux5.y[j] = conj( GTd.y[j] );

    aux2.convolutionAWT(aux4, aux5, -1, -1, aux6, aux7);
    aux2.KrammersKronig(aux2, K3, aux4, aux5, aux6, aux7);


    if(print == 1)
    {
        name = dirname + "/ThetaPP2down";
        aux2.exportAWTasFUN(name, in.display, in.xMax, 1, in.output_mode);
    }

    // THIRD PART of THETA PLUS down
    // first integrand
    for(int j = 0;          j < GTu.n+1;    j++)    aux4.y[j] = FD.y[j] * conj( GTu.y[j] ) / (2 * Pi * unity );
    for(int j = GTu.n+1;    j < 3*GTu.n+4;  j++)    aux4.y[j] = 0;
    for(int j = 3*GTu.n+4;  j < 4*GTu.n+4;  j++)    aux4.y[j] = FD.y[j] * conj( GTu.y[j] ) / (2 * Pi * unity );
    // second integrand
    for(int j = 0;          j < GTu.n+1;    j++)    aux5.y[j] = conj( GTd.y[j] );
    for(int j = GTu.n+1;    j < 3*GTu.n+4;  j++)    aux5.y[j] = 0;
    for(int j = 3*GTu.n+4;  j < 4*GTu.n+4;  j++)    aux5.y[j] = conj( GTd.y[j] );

    aux3.convolutionAWT(aux4, aux5, -1, -1, aux6, aux7);
    aux3.KrammersKronig(aux3, K3, aux4, aux5, aux6, aux7);

    if(print == 1)
    {
        name = dirname + "/ThetaPP3down";
        aux3.exportAWTasFUN(name, in.display, in.xMax, 1, in.output_mode);
    }

    // THETA PLUS
    for(int j = 0;          j < GTu.n+1;    j++)    ThetaPP.y[j] = -gamma.pp * aux1.y[j] - gamma.pp * aux2.y[j] + gamma.pm * aux3.y[j];
    for(int j = GTu.n+1;    j < 3*GTu.n+4;  j++)    ThetaPP.y[j] = 0;
    for(int j = 3*GTu.n+4;  j < 4*GTu.n+4;  j++)    ThetaPP.y[j] = -gamma.pp * aux1.y[j] - gamma.pp * aux2.y[j] + gamma.pm * aux3.y[j];


    if(print == 1)
    {
        name = dirname + "/ThetaPPdown";
        ThetaPP.exportAWTasFUN(name, in.display, in.xMax, 1, in.output_mode);
    }

    // FIRST PART of THETA MINUS
    // first integrand
    for(int j = 0;          j < GTu.n+1;    j++)    aux4.y[j] = FD.y[j] * imag( GTd.y[j] ) / Pi;
    for(int j = GTu.n+1;    j < 3*GTu.n+4;  j++)    aux4.y[j] = 0;
    for(int j = 3*GTu.n+4;  j < 4*GTu.n+4;  j++)    aux4.y[j] = FD.y[j] * imag( GTd.y[j] ) / Pi;
    // second integrand
    for(int j = 0;          j < GTu.n+1;    j++)    aux5.y[j] = conj( GTu.y[j] );
    for(int j = GTu.n+1;    j < 3*GTu.n+4;  j++)    aux5.y[j] = 0;
    for(int j = 3*GTu.n+4;  j < 4*GTu.n+4;  j++)    aux5.y[j] = conj( GTu.y[j] );

    aux1.convolutionAWT(aux4, aux5, 1, -1, aux6, aux7);
    aux1.KrammersKronigDown(aux1, K3, aux4, aux5, aux6, aux7);


    if(print == 1)
    {
        name = dirname + "/ThetaPM1down";
        aux1.exportAWTasFUN(name, in.display, in.xMax, 1, in.output_mode);
    }

    // SECOND PART of THETA MINUS
    // first integrand
    for(int j = 0;          j < GTu.n+1;    j++)    aux4.y[j] = FD.y[j] * GTu.y[j] / (2 * Pi * unity );
    for(int j = GTu.n+1;    j < 3*GTu.n+4;  j++)    aux4.y[j] = 0;
    for(int j = 3*GTu.n+4;  j < 4*GTu.n+4;  j++)    aux4.y[j] = FD.y[j] * GTu.y[j] / (2 * Pi * unity );
    // second integrand is just the Green's function up

    aux2.convolutionAWT(aux4, GTd, -1, -1, aux6, aux7);
    aux2.KrammersKronigDown(aux2, K3, aux4, aux5, aux6, aux7);


    if(print == 1)
    {
        name = dirname + "/ThetaPM2down";
        aux2.exportAWTasFUN(name, in.display, in.xMax, 1, in.output_mode);
    }

    // THIRD PART of THETA MINS
    // first integrand
    for(int j = 0;          j < GTu.n+1;    j++)    aux4.y[j] = FD.y[j] * conj( GTu.y[j] ) / (2 * Pi * unity );
    for(int j = GTu.n+1;    j < 3*GTu.n+4;  j++)    aux4.y[j] = 0;
    for(int j = 3*GTu.n+4;  j < 4*GTu.n+4;  j++)    aux4.y[j] = FD.y[j] * conj( GTu.y[j] ) / (2 * Pi * unity );
    // second integrand is just the Green's function up

    aux3.convolutionAWT(aux4, GTd, -1, -1, aux6, aux7);
    aux3.KrammersKronigDown(aux3, K3, aux4, aux5, aux6, aux7);


    if(print == 1)
    {
        name = dirname + "/ThetaPM3down";
        aux3.exportAWTasFUN(name, in.display, in.xMax, 1, in.output_mode);
    }

    // THETA PLUS
    for(int j = 0;          j < GTu.n+1;    j++)    ThetaPM.y[j] = -gamma.mp * aux1.y[j] - gamma.pp * aux2.y[j] + gamma.mp * aux3.y[j];
    for(int j = GTu.n+1;    j < 3*GTu.n+4;  j++)    ThetaPM.y[j] = 0;
    for(int j = 3*GTu.n+4;  j < 4*GTu.n+4;  j++)    ThetaPM.y[j] = -gamma.mp * aux1.y[j] - gamma.pp * aux2.y[j] + gamma.mp * aux3.y[j];


    if(print == 1)
    {
        name = dirname + "/ThetaPMdown";
        ThetaPM.exportAWTasFUN(name, in.display, in.xMax, 1, in.output_mode);
    }

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //  SIGMA Down
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // FD integrals part
    // first integrand
    for(int j = 0;          j < GTu.n+1;    j++)    aux4.y[j] = U * FD.y[j] * imag( GTu.y[j] ) / Pi;
    for(int j = GTu.n+1;    j < 3*GTu.n+4;  j++)    aux4.y[j] = 0;
    for(int j = 3*GTu.n+4;  j < 4*GTu.n+4;  j++)    aux4.y[j] = U * FD.y[j] * imag( GTu.y[j] ) / Pi;
    // second integrand
    for(int j = 0;          j < GTu.n+1;    j++)    aux5.y[j] = ThetaPM.y[j] / ( conj( Ddu.y[j] ) );
    for(int j = GTu.n+1;    j < 3*GTu.n+4;  j++)    aux5.y[j] = 0;
    for(int j = 3*GTu.n+4;  j < 4*GTu.n+4;  j++)    aux5.y[j] = ThetaPM.y[j] / ( conj( Ddu.y[j] ) );

    aux5.deleteReal(aux5);

    aux1.convolutionAWT(aux4, aux5, -1, -1, aux6, aux7);
    aux1.KrammersKronig(aux1, K3, aux4, aux5, aux6, aux7);

    if(print == 1)
    {
        name = dirname + "/FDintDown";
        aux1.exportAWTasFUN(name, in.display, in.xMax, 1, in.output_mode);
    }

    for(int j = 0;          j < GTu.n+1;    j++)    aux6.y[j] = aux4.y[j] * aux5.y[j];
    for(int j = GTu.n+1;    j < 3*GTu.n+4;  j++)    aux6.y[j] = 0;
    for(int j = 3*GTu.n+4;  j < 4*GTu.n+4;  j++)    aux6.y[j] = aux4.y[j] * aux5.y[j];

    SigmaDown.y[0] = Simpson_Re(aux6) + unity * Simpson_Im(aux6);


    for(int j = 0;          j < GTu.n+1;    j++)    aux4.y[j] = U * BE.y[j] * imag( ThetaPP.y[j]/Ddu.y[j]  ) / Pi;
    for(int j = GTu.n+1;    j < 3*GTu.n+4;  j++)    aux4.y[j] = 0;
    for(int j = 3*GTu.n+4;  j < 4*GTu.n+4;  j++)    aux4.y[j] = U * BE.y[j] * imag( ThetaPP.y[j]/Ddu.y[j]  ) / Pi;

    // second integrand
    for(int j = 0;          j < GTu.n+1;    j++)    aux5.y[j] = GTu.y[j];
    for(int j = GTu.n+1;    j < 3*GTu.n+4;  j++)    aux5.y[j] = 0;
    for(int j = 3*GTu.n+4;  j < 4*GTu.n+4;  j++)    aux5.y[j] = GTu.y[j];

    aux5.deleteReal(aux5);

    aux2.convolutionAWT(aux4, aux5, 1, -1, aux6, aux7);
    aux2.KrammersKronig(aux2, K3, aux4, aux5, aux6, aux7);


    if(print == 1)
    {
        name = dirname + "/BEintDown";
        aux2.exportAWTasFUN(name, in.display, in.xMax, 1, in.output_mode);
    }

    for(int j = 0;          j < GTu.n+1;    j++)    aux6.y[j] = aux4.y[j] * aux5.y[j];
    for(int j = GTu.n+1;    j < 3*GTu.n+4;  j++)    aux6.y[j] = 0;
    for(int j = 3*GTu.n+4;  j < 4*GTu.n+4;  j++)    aux6.y[j] = aux4.y[j] * aux5.y[j];

    SigmaDown.y[0] = SigmaUp.y[0] - Simpson_Re(aux6) - unity * Simpson_Im(aux6);


    for(int j=0; j < 4*GTu.n + 4; j++ ) SigmaDown.y[j] = aux1.y[j] - aux2.y[j];


    // export
    if(print == 1)
    {
        name = dirname + "/SSd";
        SigmaDown.exportAWTasFUN(name, in.display, in.xMax, 1, in.output_mode);
    }

}

void calculate_sigma_stat_directly(AWT & Sigma, info in, double U, double kT, AWT & Dud, AWT & Ddu, matrix & Lud, matrix & Ldu,
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

void detailed_output(info in, string dirname, AWT & GdGu1pp, AWT & GdGu1pm, AWT & GdGu2pm, AWT & GdGu2mm,
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



void normal_output(info in, string dirname, AWT & GTd, AWT & GTu, AWT & SSd, AWT & SSu, AWT & GSd, AWT & GSu, ofstream & log_output, ofstream & data_output)
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


void static_spins_matrix(string input)
{
    // constants
    double Pi = 3.14159265359;
    complex<double> unity(0,1.0);

    // set iteration precision
    double      MasterNorm = 0.0000001;
    double       SigmaNorm = MasterNorm;
    int oscil_control = 1;      // controls the oscillatory behaviour of the problem

    // intialization of variables which are then imputed by a txt file
    info in;
    in.importing(input);
    in.precisions();
    in.iterations();


    double nT, mT;
    double nT_s, mT_s;
    double kT = 0;


    // creating FERMI-DIRAC, BOSE-EINSTEIN AND K3 arrays
    AWT K3, zero;
     K3.initializeAWT(in.n, in.xMax, kT);
   zero.initializeAWT(in.n, in.xMax, kT);

     K3.Kernel3();
   zero.setZero();

    // thermodynamic and spectral self-energies, loaded with zeros
    AWT SSu, SSd;
      SSu.initializeAWT(in.n, in.xMax, kT);
      SSd.initializeAWT(in.n, in.xMax, kT);

      SSu.setZero();
      SSd.setZero();



    complex<double>     aSTu = 0;
    complex<double>     aSTd = 0;



    // GREEN FUNCTIONS arrays
    AWT GTu, GTd, GSu, GSd;
    GTu.initializeAWT(in.n, in.xMax, kT);
    GTd.initializeAWT(in.n, in.xMax, kT);
    GSu.initializeAWT(in.n, in.xMax, kT);
    GSd.initializeAWT(in.n, in.xMax, kT);


    // determinant Dud, Ddu
    AWT Dud, Ddu;
    Dud.initializeAWT(in.n, in.xMax, kT);
    Dud.setZero();

    Ddu.initializeAWT(in.n, in.xMax, kT);
    Ddu.setZero();

    // Divergent vertices ud, du
    matrix Kud, Lud, FDud, BDud;
    matrix Kdu, Ldu, FDdu, BDdu;

    // suffix _l means value of given variable from previous lambda iteration
    matrix Kud_l, Lud_l, FDud_l, BDud_l;
    matrix Kdu_l, Ldu_l, FDdu_l, BDdu_l;


    // bubbles for ud
    AWT GdGu1pp, GdGu1pm, GdGu2pm, GdGu2mm;
      GdGu1pp.initializeAWT(in.n, in.xMax, kT);
      GdGu1pm.initializeAWT(in.n, in.xMax, kT);
      GdGu2pm.initializeAWT(in.n, in.xMax, kT);
      GdGu2mm.initializeAWT(in.n, in.xMax, kT);

      GdGu1pp.setZero();
      GdGu1pm.setZero();
      GdGu2pm.setZero();
      GdGu2mm.setZero();


    // bubbles for du
    AWT GuGd1pp, GuGd1pm, GuGd2pm, GuGd2mm;
      GuGd1pp.initializeAWT(in.n, in.xMax, kT);
      GuGd1pm.initializeAWT(in.n, in.xMax, kT);
      GuGd2pm.initializeAWT(in.n, in.xMax, kT);
      GuGd2mm.initializeAWT(in.n, in.xMax, kT);

      GuGd1pp.setZero();
      GuGd1pm.setZero();
      GuGd2pm.setZero();
      GuGd2mm.setZero();

    AWT aux1, aux2, aux3, aux4, aux5, aux6, aux7;
     aux1.initializeAWT(in.n,in.xMax, kT);
     aux2.initializeAWT(in.n,in.xMax, kT);
     aux3.initializeAWT(in.n,in.xMax, kT);
     aux4.initializeAWT(in.n,in.xMax, kT);
     aux5.initializeAWT(in.n,in.xMax, kT);
     aux6.initializeAWT(in.n,in.xMax, kT);
     aux7.initializeAWT(in.n,in.xMax, kT);

    AUX help;
    help.initialize(in.n,in.xMax, kT);




    if(in.print_mode == 1)  cout << "Sigma-s, G-s, Lambda, K initialized" << endl;

    // cycle for different temperatures
    for(int iter_kT = 0; iter_kT < in.kT_iterations; iter_kT++)
    {

        kT = in.kT_min + iter_kT * in.kT_increment;

        AWT FD, BE;
        FD.initializeAWT(in.n, in.xMax, kT);
        BE.initializeAWT(in.n, in.xMax, kT);
        FD.setFD(in.n, in.xMax, kT);
        BE.setBE(in.n, in.xMax, kT);


        // cycle for different values of U
        for(int iter_U = 0; iter_U < in.U_iterations; iter_U++)
        {
            // setting the value of U in given cycle
            double U = in.U_min + iter_U * in.U_increment;

            // setting the value of h for different values of hfield
            for(int h = 0; h < in.h_increment+1; h++)
            {
                double hfield = in.h_min + h*in.h_increment;


                for(int mu_iter =0; mu_iter< in.x_iterations; mu_iter++)
                {
                    ofstream log_output;
                    ofstream data_output;

                    log_output << endl << __TIME__ << ",   " << __DATE__ << endl;

                    double mu = in.x_min + mu_iter * in.x_increment;

                    // subdirectories for different U calculations are created
                    string     U_string  =  double_to_string(U,      in.U_precision);
                    string     h_string  =  double_to_string(hfield, 3);
                    string     x_string  =  double_to_string(  mu , 2);
                    string     kT_string =  double_to_string(  kT, 3);

                    string input_string;
                    // input_string = "_L" + double_to_string(real(Lud.pp),0) + double_to_string(real(Lud.mp),0);

                    cout << input_string << endl;

                    string dirname = "Txt/U" + U_string + "_h" + h_string + "_x" + x_string + "_kT" + kT_string +  input_string;
                    createDirectories(dirname);


                    welcome_screen(U, hfield, mu, kT, dirname);

                    // we start with half filling values which we then corrected iteratively
                        aSTu = 0;
                        aSTd = 0;


                    string name;

                    // logs
                    name = dirname + "/log";
                    log_output.open(name, std::ios_base::app);
                    name = dirname + "/data";
                    data_output.open(name, std::ios_base::app);

                    set_matrix_real(   Lud, in.U_min/2);
                    set_matrix_real( Lud_l, 0);
                    set_matrix_real(   Kud, 0);
                    set_matrix_real(  FDud, 0);
                    set_matrix_real(  BDud, 0);


                    set_matrix_real(   Ldu, in.U_min/2);
                    set_matrix_real( Ldu_l,  0);
                    set_matrix_real(    Kdu, 0);
                    set_matrix_real(   FDdu, 0);
                    set_matrix_real(   BDdu, 0);

                    // output screen talks and talks
                    int out = 1;

                    // iteration enumerators
                    int iter_sigma = 0;
                    int iter_lambda = 0;
                    int iter_total = 0;



                    nT = 1;


                    //for(int iter_sigma=0; iter_sigma<sigma_iter_max; iter_sigma++)
                    for(iter_sigma = 0; iter_sigma<1; iter_sigma++)
                    //while(abs(real(totalSTd) - real(totalSTd_s) ) > 0.00001)
                    {

                        cout << endl << "   SIGMA " << iter_sigma << "-th iter: " << endl << endl;

                        // iterating vertices
                        iter_lambda = 0;
                        //while( (iter_lambda == 0)  )
                        while( (iter_lambda == 0) || (abs(real(Lud.mp) - real(Lud_l.mp) ) > 0.00001)  )
                        //for(int iter_lambda=0; iter_lambda<1; iter_lambda++)
                        {
                            cout << "       " << iter_lambda << "-th iter: " ;
                            if(out == 1)   cout << endl << "          ";


                            if(iter_lambda == 0)
                            {
                                // 0. step: setting the initial GTd, GTu
                                //          calculate bubbles

                                freePropagatorLorentzShift(in.delta, mu, U * (1.0 - nT) / 2.0 + hfield - aSTd, GTd);
                                freePropagatorLorentzShift(in.delta, mu, U * (1.0 - nT) / 2.0 - hfield + aSTu, GTu);

                                calculateGGud(GdGu1pp, GdGu1pm, GdGu2pm, GdGu2mm, GTd, GTu, FD, K3, help, out);
                            }


                            // 1. step: calculate determinant
                            calculate_Dud_static(Dud, GdGu1pp, GdGu1pm, GdGu2pm, GdGu2mm, Lud, out);
                            calculate_Ddu_sym(Ddu, Dud, out);


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


                            calculate_Kud_static(Kud, GdGu1pp, GdGu1pm, GdGu2pm, GdGu2mm, Lud, out);
                            calculate_Kdu_sym(Kdu, Kud, out);

                            calculate_FDud(FDud, GTd, GTu, Dud, FD, aux1, out);
                            calculate_FDdu_sym(FDdu, FDud, out);

                            calculate_BDud(BDud, GTd, GTu, Dud, BE, aux1, out);
                            calculate_BDdu_sym(BDdu, BDud, out);

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

                            calculate_Lud(Lud, Kud, FDud, BDud, U, out);
                            calculate_Ldu_sym(Ldu, Lud, out);


                            if(oscil_control == 1)
                            {
                                aver_matrix(Lud, Lud_l);
                                aver_matrix(Ldu, Ldu_l);
                            }


                            lambda_iter_console(out, Lud, Ldu, Kud, Kdu, BDud, BDdu, FDud, FDdu);

                            if((iter_lambda == 0) && (iter_sigma == 0) ) log_output << endl << __TIME__ << ",   " << __DATE__ << endl << endl;

                            lambda_iter_output(log_output, iter_lambda, iter_sigma, aSTu, aSTd, Lud, Ldu, Kud, Kdu, BDud, BDdu, FDud, FDdu);


                            // calculate new anomalous self-energy
                            calculate_aST(aSTu, aSTd, Lud, Ldu, GTd, GTu, FD, aux1, out);
                            aSTu = ( aSTu + aSTd ) / 2.0;
                            aSTd = ( aSTu + aSTd ) / 2.0;

                            data_output << iter_total
                            << " " << real( Lud.pp) << " " << imag( Lud.pp) << " " << real( Lud.pm) << " " << imag( Lud.pm)
                            << " " << real( Ldu.pp) << " " << imag( Ldu.pp) << " " << real( Ldu.pm) << " " << imag( Ldu.pm)

                            << " " << real( Kud.pp) << " " << imag( Kud.pp) << " " << real( Kud.pm) << " " << imag( Kud.pm)
                            << " " << real( Kdu.pp) << " " << imag( Kdu.pp) << " " << real( Kdu.pm) << " " << imag( Kdu.pm)

                            << " " << real(BDud.pp) << " " << imag(BDud.pp) << " " << real(BDud.pm) << " " << imag(BDud.pm)
                            << " " << real(BDdu.pp) << " " << imag(BDdu.pp) << " " << real(BDdu.pm) << " " << imag(BDdu.pm)

                            << " " << real(FDud.pp) << " " << imag(FDud.pp) << " " << real(FDud.pm) << " " << imag(FDud.pm)
                            << " " << real(FDdu.pp) << " " << imag(FDdu.pp) << " " << real(FDdu.pm) << " " << imag(FDdu.pm) << endl;


                            iter_lambda = iter_lambda + 1;
                            iter_total  = iter_total  + 1;

                            calculate_nT_mT(nT, mT, mu, U, GTu, GTd, FD, aux1);

                        } // end of Lambda iterations


                        // store old values
                        if(iter_sigma>0 && oscil_control ==0)
                        {
                            nT_s = nT;
                            mT_s = mT;
                        }

                        // anomalouc contribution to thermodynamic self-energy
                        //calculate_ST_anomalous(aSTd, aSTu, Lud, GTd, GTu, FD, aux1, out);

                        calculate_nT_mT(nT, mT, mu, U, GTu, GTd, FD, aux1);



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

                    lambda_iter_output(final_output, iter_lambda, iter_sigma, aSTu, aSTd, Lud, Ldu, Kud, Kdu, BDud, BDdu, FDud, FDdu);


                    calculate_sigma_stat(SSu, SSd, in, U, kT, Dud, Ddu, Lud, Ldu, GTd, GTu, FD, BE, K3, dirname, 0, aux1, aux2, aux3, aux4, aux5, aux6, aux7);


                    // iterate nS
                    int iter_nS = 0;
                    double nS = 1;
                    double mS = 0;
                    nT = 1;
                    double nS_old = -1.0;

                    while(  iter_nS < 1 )
                    //while( ( abs(nS - nS_old) > 0.000001 ) && iter_nS < 1000 )
                    {

                        AWT SS0;
                        SS0.initializeAWT(in.n, in.xMax, kT);
                        for(int j = 0;          j < GTu.n+1;    j++)    SS0.y[j] = ( SSu.y[j] + SSd.y[j] ) / 2.0;
                        for(int j = GTu.n+1;    j < 3*GTu.n+4;  j++)    SS0.y[j] = 0;
                        for(int j = 3*GTu.n+4;  j < 4*GTu.n+4;  j++)    SS0.y[j] = ( SSu.y[j] + SSd.y[j] ) / 2.0;


                        propagatorLorentzShift(in.delta, mu + U / 2.0 - U * nS / 2.0, - hfield + aSTu, GSu, SS0);
                        propagatorLorentzShift(in.delta, mu + U / 2.0 - U * nS / 2.0, + hfield - aSTd, GSd, SS0);

                        calculate_nS_mS(nS, mS, GSu, GSd, FD, aux1);
                        nS_old = nS;

                        //nS = ( nS + nS_old ) / 2.0;

                        cout << mu + U * nT / 2.0 - U * nS / 2.0 << endl;
                        cout << - hfield + aSTu << endl;
                        cout << + hfield - aSTd << endl;

                              cout << "          nS: " << nS <<  endl << "          ";
                        log_output << "          nS: " << nS << endl << endl;
                        iter_nS = iter_nS + 1;
                    }

                    // OUTPUTS
                    int print = 1;
                    if(print == 1)
                    {
                        cout << "print details";
                        detailed_output(in, dirname, GdGu1pp, GdGu1pm, GdGu2pm, GdGu2mm, Dud, Ddu, GTd, GTu);
                    }

                          normal_output(in, dirname, GTd, GTu, SSd, SSu, GSd, GSu, log_output, data_output);


                } // end of mu iteration instance

            } // end of h iteration instance

        }  // end of U iteration instance

    }  // end of kT iteration instance


} // end of the whole function

