#include "physics.h"

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


void createDirectories(string & dir_path)
{
	boost::filesystem::path dir(dir_path);
	if(boost::filesystem::create_directory(dir))
	{
		//cout << "Success" << "\n";
	}
}


void precision(double kT_min, double U_min, double x_increment, int & kT_precision, int & U_precision, int & x_precision)
{
    // the number of decimals in kT_min is determined
    kT_precision = 0;
    U_precision = 2;
    x_precision = 1;

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

void  iterations(double kT_max, double kT_min, double kT_increment,
                 double U_max, double U_min, double U_increment,
                 double x_range, double x_increment,
                 int & kT_iterations, int & U_iterations, int & x_iterations)
{
    if(kT_increment == 0)   kT_iterations = 1;
    else                    kT_iterations = ( kT_max - kT_min) / kT_increment + 1;

    // variable for maximal number of U iteration is initialized and set
    if(U_increment == 0)    U_iterations = 1;
    else                    U_iterations  = ( U_max -  U_min) /  U_increment + 1;

    // variable for maximal number of kT iteration is initialized and set
    if(x_increment == 0)    x_iterations = 1;
    else                    x_iterations  = (     x_range    ) /  x_increment + 1;

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

////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////                             Creating kT files                                 ////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////

void hartree_energies(double & Ehf,
                      int n, double xMax, double delta, int display, double range, string export_mode,
                      int x_iterations, double x_increment,
                      AWT & FD_0, AWT & BE_0)
{
    double Pi = 3.14159265359;

    double   n_0[2*x_iterations-1];
    double  E_hf[2*x_iterations-1];


    for(int j = 0; j < 2*x_iterations -1; j++)
    {

        // setting the right value of x
        double x;
        if(j > x_iterations - 1)        x = (x_iterations - 1 -j)*x_increment;
        else                            x = j*x_increment;

        // propagator is calculated and outputed
        AWT G_0;
        G_0.initializeAWT(n, xMax, 0);
        freePropagatorLorentz(delta, x, G_0);

        string fileName = "Txt/Hartree/Ghart_" + boost::lexical_cast<std::string>(x);
        G_0.exportAWTasFUN(fileName, display, range, 1, export_mode);

        // occupancy is calculated
        AWT Int;
        Int.initializeAWT(n, xMax, 0);
        for( int iii = 0;          iii<Int.n+1;     iii++ )
        Int.y[iii] = FD_0.y[iii]*imag(G_0.y[iii]);


        for( int iii = 3*Int.n+4; iii<4*Int.n+4; iii++ )
        Int.y[iii] = FD_0.y[iii]*imag(G_0.y[iii]);

        n_0[j] = -2*Trapezoid_Re(Int)/Pi;

        // total energy is calculated

        double factor_1 = xMax/n;
        double factor_2 = delta*delta;

        for( int iii = 0;               iii<Int.n+1;    iii++ )
        Int.y[iii] = FD_0.y[iii] *  (      iii*factor_1 + x     )  / (        (iii*factor_1 + x)*(iii*factor_1 + x) + factor_2          ) ;

        for( int iii = 3*Int.n+4; iii<4*Int.n+4; iii++ )
        Int.y[iii] = FD_0.y[iii]  * (   (iii-4*n-4)*factor_1 + x  ) / (   ((iii-4*n-4)*factor_1 + x)*((iii-4*n-4)*factor_1 + x)  + factor_2   ) ;

        E_hf[j] = 2.0*Trapezoid_Re(Int)/Pi;

        // results are outputed
        string hartree = "Txt/Hartree/hartree";
        ofstream hartree_output(hartree.c_str());
        // x < 0 values are outputed first
        for(int ii = 2*x_iterations-2; ii > x_iterations -1 ; ii --)
        hartree_output <<  (x_iterations - 1 -ii)*x_increment << "   " <<   E_hf[ii]    << endl;
        // x >=0 values are outputed
        for(int ii = 0; ii < x_iterations; ii++)
        hartree_output <<           ii*x_increment            << "   " <<   E_hf[ii]    << endl;

        Ehf = E_hf[0];
    }
}


void preset(int j, double & x, int x_iterations, double x_increment, double & Lambda, double & LambdaHalf, double & nT)
{


    // setting the right value of x
    if(j > x_iterations - 1)        x = (x_iterations - 1 -j)*x_increment;
    else                            x = j*x_increment;

    // setting initial value of Lambda for negative x
    if(j == x_iterations)
    {
        Lambda = LambdaHalf;
        nT = 0;
    }

    // setting initial value of Lambda for x=0
    if(j == 0)
    {
        Lambda = 0;
        nT = 1;
    }
}


void spectral_sigma_RPA(double U, AWT & PHI, AWT & GT, AWT & FD, AWT & BE, AWT & K3, AWT & SigmaSpec, AWT & Int,
                        AWT & aux1, AWT & aux2, AWT & aux3, AWT & aux4, AWT & aux5, AWT & aux6)
{

    for( int iii = 0;               iii<Int.n+1;     iii++ )
    Int.y[iii] = -U * U * PHI.y[iii] / ( 1.0 + U * PHI.y[iii]);

    for( int iii = Int.n+1;   iii<3*Int.n+4; iii++ )
    Int.y[iii] = 0;

    for( int iii = 3*Int.n+4; iii<4*Int.n+4; iii++ )
    Int.y[iii] = -U * U * PHI.y[iii] / ( 1.0 + U * PHI.y[iii]);

    SigmaSpec.boseMatsubaraImP(Int, GT, 1, 1, FD, BE, aux1, aux2, aux3, aux4, aux5, aux6);
    SigmaSpec.KrammersKronig(SigmaSpec,K3, aux1, aux2, aux3, aux4);
}



void SigmaTherm_frequency(AWT & SigmaTherm, AWT & Lambda, double U, double nT, AWT & GT, AWT & FD, AWT & BE, AWT & K3,
                        AWT & aux1, AWT & aux2, AWT & aux3, AWT & aux4, AWT & aux5, AWT & aux6)
{
    // theory from the end of august 2017
	SigmaTherm.boseMatsubaraImP(GT, Lambda, 1, 1, FD, BE, aux1, aux2, aux3, aux4, aux5, aux6);
    SigmaTherm.KrammersKronig(SigmaTherm, K3, aux1, aux2, aux3, aux4);

	for(int i=0;                i<SigmaTherm.n+1;   i++ )       SigmaTherm.y[i] = nT*U/2.0 + SigmaTherm.y[i] ;
	for(int i=SigmaTherm.n+1;   i<3*SigmaTherm.n+3; i++ )       SigmaTherm.y[i] = 0;
	for(int i=3*SigmaTherm.n+4; i<4*SigmaTherm.n+4; i++ )       SigmaTherm.y[i] = nT*U/2.0 + SigmaTherm.y[i] ;
}

void SigmaTherm_frequency_2(AWT & SigmaTherm, AWT & Lambda, double U, double nT)
{
    // theory from the end of december 2017
	for(int i=0;                i<SigmaTherm.n+1;   i++ )       SigmaTherm.y[i] = U/2 - (0.5 - nT)*Lambda.y[i] ;
	for(int i=SigmaTherm.n+1;   i<3*SigmaTherm.n+3; i++ )       SigmaTherm.y[i] = 0;
	for(int i=3*SigmaTherm.n+4; i<4*SigmaTherm.n+4; i++ )       SigmaTherm.y[i] = U/2 - (0.5 - nT)*Lambda.y[i] ;
}

void SigmaSpec_frequency(AWT & SigmaSpec, double U, AWT & LPp, AWT & GT, AWT & FD, AWT & BE, AWT & K3,
                        AWT & aux1, AWT & aux2, AWT & aux3, AWT & aux4, AWT & aux5, AWT & aux6)
{
	// theory from the end of september 2017
	AWT help;
	help.initializeAWT(GT.n, GT.xMax, GT.kT);
	for(int i=0;          i<help.n+1;   i++ )       help.y[i] = -U*LPp.y[i] / ( 1.0 + LPp.y[i] ) ;
	for(int i=help.n+1;   i<3*help.n+3; i++ )       help.y[i] = 0;
	for(int i=3*help.n+4; i<4*help.n+4; i++ )       help.y[i] = -U*LPp.y[i] / ( 1.0 + LPp.y[i] ) ;

	SigmaSpec.boseMatsubaraImP(help, GT, 1, 1, FD, BE, aux1, aux2, aux3, aux4, aux5, aux6);
    SigmaSpec.KrammersKronig(SigmaSpec, K3, aux1, aux2, aux3, aux4);
}

////////////////////////////////////////////////////////////////////////////////////////////////
//////////                      Physical quantities                                     ///////
//////////////////////////////////////////////////////////////////////////////////////////////


void LPp_function(AWT & LPp, AWT & Lambda, AWT & GT, AWT & FD, AWT & BE, AWT & K3,
                        AWT & aux1, AWT & aux2, AWT & aux3, AWT & aux4, AWT & aux5, AWT & aux6)
{

	double Pi=3.14159265359;

    AWT y1;
    y1.initializeAWT(GT.n,   GT.xMax,  GT.kT);
    AWT y2;
    y2.initializeAWT(GT.n,   GT.xMax,  GT.kT);

    // FIRST PART is calculated
    for(int i=0;         i<GT.n+1;   i++ )       aux1.y[i] = FD.y[i] * imag( GT.y[i] );
    for(int i=GT.n+1;    i<3*GT.n+3; i++ )       aux1.y[i] = 0;
    for(int i=3*GT.n+4;  i<4*GT.n+4; i++ )       aux1.y[i] = FD.y[i] * imag( GT.y[i] );

    for(int i=0;         i<GT.n+1;   i++ )       aux2.y[i] = Lambda.y[i] * GT.y[i] ;
    for(int i=GT.n+1;    i<3*GT.n+3; i++ )       aux2.y[i] = 0;
    for(int i=3*GT.n+4;  i<4*GT.n+4; i++ )       aux2.y[i] = Lambda.y[i] * GT.y[i] ;


    aux2.deleteReal(aux2);
    aux1.convolutionAWT(aux1, aux2, 1, -1, aux5, aux6);     // integral over real arguments
    aux1.multiplyAWT(aux1, -1/Pi);                          // multiply with proper factor


    // SECOND PART is calculated
    for(int i=0;         i<GT.n+1;   i++ )       aux3.y[i] = FD.y[i] * imag( Lambda.y[i] * GT.y[i]);
    for(int i=GT.n+1;    i<3*GT.n+3; i++ )       aux3.y[i] = 0;
    for(int i=3*GT.n+4;  i<4*GT.n+4; i++ )       aux3.y[i] = FD.y[i] * imag( Lambda.y[i] * GT.y[i]);

    for(int i=0;         i<GT.n+1;   i++ )       aux4.y[i] = GT.y[i] ;
    for(int i=GT.n+1;    i<3*GT.n+3; i++ )       aux4.y[i] = 0;
    for(int i=3*GT.n+4;  i<4*GT.n+4; i++ )       aux4.y[i] = GT.y[i] ;

    aux4.deleteReal(aux4);
    aux2.convolutionAWT(aux3, aux4, -1, -1, aux5, aux6);    // integral over real arguments
    aux2.multiplyAWT(aux2, 1/Pi);                           // multiply with proper factor


    // first and second part are added tgether and stored in the out AWT
    for(int i=0; i < GT.nn; i++)         LPp.y[i] = aux1.y[i] + aux2.y[i];

    LPp.KrammersKronig(LPp, K3, aux1, aux2, aux3, aux4);
}




void LLPp_function(AWT & LLPp, AWT & Lambda, AWT & GT, AWT & FD, AWT & BE, AWT & K3,
                        AWT & aux1, AWT & aux2, AWT & aux3, AWT & aux4, AWT & aux5, AWT & aux6)
{
	AWT help;
	help.initializeAWT(GT.n, GT.xMax, GT.kT);
	for(int i=0;          i<help.n+1;   i++ )       help.y[i] = Lambda.y[i] * Lambda.y[i] * GT.y[i] ;
	for(int i=help.n+1;   i<3*help.n+3; i++ )       help.y[i] = 0;
	for(int i=3*help.n+4; i<4*help.n+4; i++ )       help.y[i] = Lambda.y[i] * Lambda.y[i] * GT.y[i] ;

	LLPp.fermMatsubaraImP(GT, help, 1, 1, FD, aux1, aux2, aux3, aux4, aux5, aux6);
    LLPp.KrammersKronig(LLPp, K3, aux1, aux2, aux3, aux4);
}

void KPm_function(AWT & KPm, AWT & K, AWT & GT, AWT & FD, AWT & BE, AWT & K3,
                        AWT & aux1, AWT & aux2, AWT & aux3, AWT & aux4, AWT & aux5, AWT & aux6)
{
    double Pi=3.14159265359;

    // FIRST PART is calculated
    for(int i=0;        i<K.n+1;   i++ )       aux1.y[i] = FD.y[i] * imag( GT.y[i] );
    for(int i=K.n+1;    i<3*K.n+3; i++ )       aux1.y[i] = 0;
    for(int i=3*K.n+4;  i<4*K.n+4; i++ )       aux1.y[i] = FD.y[i] * imag( GT.y[i] );

    for(int i=0;        i<K.n+1;   i++ )       aux2.y[i] = K.y[i] * GT.y[i] ;
    for(int i=K.n+1;    i<3*K.n+3; i++ )       aux2.y[i] = 0;
    for(int i=3*K.n+4;  i<4*K.n+4; i++ )       aux2.y[i] = K.y[i] * GT.y[i] ;

    aux2.deleteReal(aux2);
    aux1.convolutionAWT(aux1, aux2, 1, 1, aux5, aux6);      // integral over real arguments
    aux1.multiplyAWT(aux1, -1/Pi);                          // multiply with proper factor


    // SECOND PART is calculated
    for(int i=0;        i<K.n+1;   i++ )       aux3.y[i] = -1.0 * FD.y[i] * imag( K.y[4*K.n +4 -i ] * GT.y[4*K.n +4 -i ] );
    for(int i=K.n+1;    i<3*K.n+3; i++ )       aux3.y[i] = 0;
    for(int i=3*K.n+4;  i<4*K.n+4; i++ )       aux3.y[i] = -1.0 * FD.y[i] * imag( K.y[4*K.n +4 -i ] * GT.y[4*K.n +4 -i ] );

    for(int i=0;        i<K.n+1;   i++ )       aux4.y[i] = GT.y[i] ;
    for(int i=K.n+1;    i<3*K.n+3; i++ )       aux4.y[i] = 0;
    for(int i=3*K.n+4;  i<4*K.n+4; i++ )       aux4.y[i] = GT.y[i] ;

    aux4.deleteReal(aux4);
    aux2.convolutionAWT(aux3, aux4, 1, -1, aux5, aux6);     // integral over real arguments
    aux2.multiplyAWT(aux2, -1/Pi);                          // multiply with proper factor


    // first and second part are added tgether and stored in the out AWT
    for(int i=0; i < K.nn; i++)         KPm.y[i] = aux1.y[i] + aux2.y[i];
    KPm.KrammersKronig(KPm, K3, aux1, aux2, aux3, aux4);


/* OLD VERSION OF CODE
	AWT help;
	help.initializeAWT(GT.n, GT.xMax, GT.kT);
	for(int i=0;          i<help.n+1;   i++ )       help.y[i] = K.y[i] * GT.y[i] ;
	for(int i=help.n+1;   i<3*help.n+3; i++ )       help.y[i] = 0;
	for(int i=3*help.n+4; i<4*help.n+4; i++ )       help.y[i] = K.y[i] * GT.y[i] ;

	KPm.fermMatsubaraImP(GT, help, 1, -1, FD);//, BE);
    KPm.KrammersKronig(KPm, K3);
*/

}




// Kvertex is obtained
void KvertexFunction(double Lambda, AWT & Phi, AWT & Kvertex)
{
    for(int i=0; i<=4*Kvertex.n+4; i++)
    Kvertex.y[i] = -1.0 * Lambda * Lambda * Phi.y[i] / ( 1.0 + Lambda * Phi.y[i] +1e-12 );
}

// Phi bubble at arbitrary T is obtained
void PhiFunction(AWT & Gin, AWT & Phi, AWT & f_FD, AWT & K3,
                        AWT & aux1, AWT & aux2, AWT & aux3, AWT & aux4, AWT & aux5, AWT & aux6)
{
    // setting Phi
    Phi.fermMatsubaraImP(Gin, Gin, 1, 1, f_FD, aux1, aux2, aux3, aux4, aux5, aux6);
    Phi.KrammersKronig(Phi,K3, aux1, aux2, aux3, aux4);
}

// thermal Phi bubble shifted into the lower complex half-plane is obtained
void PhiFunctionShift(AWT & GT, AWT & GTshift, AWT & PHIshift, AWT & FD, AWT & K3,
                        AWT & aux1, AWT & aux2, AWT & aux3, AWT & aux4, AWT & aux5, AWT & aux6)
{
    double Pi = 3.14159265359;
    AWT Int1;
    Int1.initializeAWT(GT.n, GT.xMax, GT.kT);
    Int1.normalFDtimesIM(GT, FD);

    AWT Int2;
    Int2.initializeAWT(GT.n, GT.xMax, GT.kT);
    Int2.deleteReal(GTshift);

    // convolution
    Int2.convolutionAWT(Int1, Int2, 1, -1, aux1, aux2);
    Int2.multiplyAWT(Int2, -1/Pi);

    AWT Int3;
    Int3.initializeAWT(GT.n, GT.xMax, GT.kT);
    Int3.deleteReal(GTshift);
    Int3.conjugateY(Int3);

    Int3.convolutionAWT(Int1, Int3, -1, -1, aux1, aux2);
    Int3.multiplyAWT(Int3, -1/Pi);

    for(int i=0; i < 4*GT.n+4; i++)
    PHIshift.y[i] = Int2.y[i] + Int3.y[i];
    PHIshift.KrammersKronigDown(PHIshift,K3, aux1, aux2, aux3, aux4);
}

// Psi function is obtained
void PsiFunction(AWT & Gtherm, AWT & Kvertex, AWT & BE, double & Psi)
{
    double Pi = 3.14159265359;
    AWT Integrand;
    Integrand.initializeAWT(Gtherm.n, Gtherm.xMax, Gtherm.kT);
    Integrand.y[0] = BE.y[0] * imag(   Gtherm.y[0] * conj(Gtherm.y[0]) * conj(Kvertex.y[0])  )/Pi;
    for( int i = 1;         i < Gtherm.n+1;      i++ )
    Integrand.y[i] = BE.y[i] * imag(   Gtherm.y[i] * conj(Gtherm.y[4*Gtherm.n+4-i]) * conj(Kvertex.y[4*Gtherm.n+4-i])  )/Pi;
    for( int i = Gtherm.n+1;   i < 3*Gtherm.n + 4;  i++ )
    Integrand.y[i] = 0;
    for( int i = 3*Gtherm.n+4; i < 4*Gtherm.n + 4;  i++ )
    Integrand.y[i] = BE.y[i] * imag(   Gtherm.y[i] * conj(Gtherm.y[4*Gtherm.n+4-i]) * conj(Kvertex.y[4*Gtherm.n+4-i])  )/Pi;

    Psi = Trapezoid_Re(Integrand);
}




void DensityLorentz(double del, double mu, AWT & DensityOut, AWT & Sigma)
{
    double n = DensityOut.n;
    double x = DensityOut.xMax;
    double Pi=3.14159265359;
    double SigmaT = real(Sigma.y[0]);

    for(int i=0; i<n+1; i++ )
    {
        double xx = i*x/n;
        double u =  (  1.0 - 2.0 * atan((xx-mu+SigmaT)/del)/Pi );
        DensityOut.y[i] = u;
    }

    for(int i=n+1; i<3*n+3; i++ )
    {
        double u = 0;
        DensityOut.y[i] = u;
    }


    for(int i=3*n+4; i<4*n+4; i++ )
    {
        double xx = (i - 4*n -4) *x/n;
        double u =  (  1.0 - 2.0 * atan((xx-mu+SigmaT)/del)/Pi );
        DensityOut.y[i] = u;
    }

}

////////////////////////////////////////////////////////////////////////////////////////////
/////////////////            Lorentz Propagator shifted            /////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////
// propagator is shifted from the real axis into the complex plane by the argument ShiftIm
// positive ShiftIm means propagator is G(x + i*ShiftIm)

////////////////////////////////////////////////////////////////////////////////////////////
/////////////////                 Lorentz Propagator              /////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////

void propagatorLorentz(double del, double x_, AWT & Gout, AWT & Sigma)
{
    double n = Gout.n;
    double x = Gout.xMax;

    for(int i=0; i<n+1; i++ )
    {
        //double del = 1e-7;
        double xx = i*x/n;

        double norm =  (  xx + x_ - Sigma.y[i].real() )*(  xx + x_ - Sigma.y[i].real()  )
                    +  (    del - Sigma.y[i].imag()   )*(  del - Sigma.y[i].imag()      );

        complex<double> u(
                             ( xx + x_ - Sigma.y[i].real()  )  /norm ,
                            -(   del - Sigma.y[i].imag()    )  /norm
                          );
        Gout.y[i] = u;
    }

    for(int i=n+1; i<3*n+3; i++ )
    {
        complex<double> u( 0 , 0 );
        Gout.y[i] = u;
    }


    for(int i=3*n+4; i<4*n+4; i++ )
    {
        //double del = 1e-7;
        double xx = (i - 4*n -4) *x/n;

        double norm =    (   xx + x_ - Sigma.y[i].real() ) * (   xx + x_ - Sigma.y[i].real() )
                    +    (  del - Sigma.y[i].imag() ) * (  del - Sigma.y[i].imag() );

        complex<double> u(
                           (  xx + x_  - Sigma.y[i].real()  )/norm,
                          -(  del - Sigma.y[i].imag()  )/norm
                         );
        Gout.y[i] = u;
    }

}





void propagatorLorentzShift(double del, double x_, complex<double> Shift, AWT & Gout, AWT & Sigma)
{
    double n = Gout.n;
    double x = Gout.xMax;

    for(int i=0; i<n+1; i++ )
    {
        //double del = 1e-7;
        double xx = i*x/n;

        double norm =  (  xx + x_ - real(Shift) - Sigma.y[i].real() )*(  xx + x_ - real(Shift) - Sigma.y[i].real()  )
                    +  (    del - imag(Shift) - Sigma.y[i].imag()   )*(  del - imag(Shift) - Sigma.y[i].imag()      );

        complex<double> u(
                             ( xx + x_ - real(Shift) - Sigma.y[i].real()  )  /norm ,
                            -(   del - imag(Shift) - Sigma.y[i].imag()    )  /norm
                          );
        Gout.y[i] = u;
    }

    for(int i=n+1; i<3*n+3; i++ )
    {
        complex<double> u( 0 , 0 );
        Gout.y[i] = u;
    }


    for(int i=3*n+4; i<4*n+4; i++ )
    {
        //double del = 1e-7;
        double xx = (i - 4*n -4) *x/n;

        double norm =    (   xx + x_ - real(Shift) - Sigma.y[i].real() ) * (   xx + x_ - real(Shift) - Sigma.y[i].real() )
                    +    (  del - imag(Shift) - Sigma.y[i].imag() ) * (  del - imag(Shift) - Sigma.y[i].imag() );

        complex<double> u(
                           (  xx + x_ - real(Shift) - Sigma.y[i].real()  )/norm,
                          -(  del - imag(Shift) - Sigma.y[i].imag()  )/norm
                         );
        Gout.y[i] = u;
    }

}



void propagatorNssLorentzShift(double del, double x_, double gap, double gammaS, double phi, complex<double> Shift, AWT & Gout, AWT & Sigma)
{
    double n = Gout.n;
    double x = Gout.xMax;

    double Pi = 3.14159265;

    // SMALL SHIFT TO AVOID SINGULARITIES
    gap = gap - 10e-12;


    // POSITIVE FREQUENCIES
    for(int i=0; i<n+1; i++ )
    {
        //double del = 1e-7;
        double xx = i*x/n;


        if( xx < gap)       // IN THE GAP
        {
            double norm =  (  xx + x_ - real(Shift) - Sigma.y[i].real() )*(  xx + x_ - real(Shift) - Sigma.y[i].real()  )
                        +  (    del - imag(Shift) - Sigma.y[i].imag()   )*(  del - imag(Shift) - Sigma.y[i].imag()      );

            complex<double> u(
                             ( xx + x_ - real(Shift) - Sigma.y[i].real()  )  /norm ,
                            -(   del - imag(Shift) - Sigma.y[i].imag()    )  /norm
                          );
            Gout.y[i] = u;
        }
        else                // OUT OF THE GAP
        {
            double norm =  (  xx + gammaS * xx/( sqrt(xx*xx - gap*gap) ) * ( 1.0 - (gap/xx)  * cos(Pi*phi/2) ) + x_ - real(Shift) - Sigma.y[i].real() )
                          *(  xx + gammaS * xx/( sqrt(xx*xx - gap*gap) ) * ( 1.0 - (gap/xx)  * cos(Pi*phi/2) ) + x_ - real(Shift) - Sigma.y[i].real() )
                        +  (    del - imag(Shift) - Sigma.y[i].imag()   )*(  del - imag(Shift) - Sigma.y[i].imag()      );

            complex<double> u(
                             ( xx + gammaS * xx/( sqrt(xx*xx - gap*gap) ) * ( 1.0 - (gap/xx)  * cos(Pi*phi/2) ) + x_ - real(Shift) - Sigma.y[i].real() ) /norm,
                            -(   del - imag(Shift) - Sigma.y[i].imag()    )  /norm
                             );
            Gout.y[i] = u;
        }


    }

    // ZERO PADDING
    for(int i=n+1; i<3*n+3; i++ )
    {
        complex<double> u( 0 , 0 );
        Gout.y[i] = u;
    }


    // NEGATIVE FREQUENCIES
    for(int i=3*n+4; i<4*n+4; i++ )
    {
        //double del = 1e-7;
        double xx = (i - 4*n -4) *x/n;

        if( xx > -1.0*gap)       // IN THE GAP
        {
            double norm =  (  xx + x_ - real(Shift) - Sigma.y[i].real() )*(  xx + x_ - real(Shift) - Sigma.y[i].real()  )
                        +  (    del - imag(Shift) - Sigma.y[i].imag()   )*(  del - imag(Shift) - Sigma.y[i].imag()      );

            complex<double> u(
                             ( xx + x_ - real(Shift) - Sigma.y[i].real()  )  /norm ,
                            -(   del - imag(Shift) - Sigma.y[i].imag()    )  /norm
                          );
            Gout.y[i] = u;
        }
        else                // OUT OF THE GAP
        {
            double norm =  (  xx + gammaS * xx/( sqrt(xx*xx - gap*gap) ) * ( 1.0 - (gap/xx)  * cos(Pi*phi/2) ) + x_ - real(Shift) - Sigma.y[i].real() )
                          *(  xx + gammaS * xx/( sqrt(xx*xx - gap*gap) ) * ( 1.0 - (gap/xx)  * cos(Pi*phi/2) ) + x_ - real(Shift) - Sigma.y[i].real() )
                        +  (    del - imag(Shift) - Sigma.y[i].imag()   )*(  del - imag(Shift) - Sigma.y[i].imag()      );

            complex<double> u(
                             ( xx + gammaS * xx/( sqrt(xx*xx - gap*gap) ) * ( 1.0 - (gap/xx)  * cos(Pi*phi/2) ) + x_ - real(Shift) - Sigma.y[i].real() ) /norm,
                            -(   del - imag(Shift) - Sigma.y[i].imag()    )  /norm
                             );
            Gout.y[i] = u;
        }
    }

}


void freePropagatorLorentz(double del, double mu, AWT & Gout)
{
    double n = Gout.n;
    double x = Gout.xMax;

    for(int i=0; i<n+1; i++ )
    {
        //double del = 1e-7;
        double xx = i*x/n;

        double norm =  (  xx + mu  )*(  xx + mu  )
                    +  (   del     )*(   del     );

        complex<double> u(
                             ( xx + mu   )  /norm ,
                            -(   del     )  /norm
                          );
        Gout.y[i] = u;
    }

    for(int i=n+1; i<3*n+3; i++ )
    {
        complex<double> u( 0 , 0 );
        Gout.y[i] = u;
    }


    for(int i=3*n+4; i<4*n+4; i++ )
    {
        //double del = 1e-7;
        double xx = (i - 4*n -4) *x/n;

        double norm =    (   xx + mu  )*(   xx + mu  )
                    +     (    del    )*(     del    );

        complex<double> u(
                           (  xx + mu  )/norm,
                          -(  del      )/norm
                         );
        Gout.y[i] = u;
    }

}


void freeNssPropagatorLorentz(double del, double mu, double gap, double gammaS, double phi, AWT & Gout)
{
    double n = Gout.n;
    double x = Gout.xMax;
    double Pi = 3.14159265;

    // SMALL SHIFT TO AVOID SINGULARITIES
    gap = gap - 10e-12;


    // POSITIVE FREQUENCIES
    for(int i=0; i<n+1; i++ )
    {
        //double del = 1e-7;
        double xx = i*x/n;

        if(xx < gap)        // IN THE GAP
        {
            double norm =  (  xx + mu  )
                          *(  xx + mu  )
                        +  (   del     )*(   del     );

            complex<double> u(
                             ( xx + mu   )  /norm ,
                            -(   del     )  /norm
                          );
            Gout.y[i] = u;
        }
        else             // OUT OF THE GAP
        {
            double norm =  (  xx + gammaS * xx/( sqrt(xx*xx - gap*gap) ) * ( 1.0 - (gap/xx)  * cos(Pi*phi/2) ) + mu  )
                          *(  xx + gammaS * xx/( sqrt(xx*xx - gap*gap) ) * ( 1.0 - (gap/xx)  * cos(Pi*phi/2) ) + mu  )
                        +  (   del     )*(   del     );

            complex<double> u(
                             ( xx + gammaS * xx/( sqrt(xx*xx - gap*gap) ) * ( 1.0 - (gap/xx)  * cos(Pi*phi/2) ) + mu   )  /norm ,
                            -(   del     )  /norm
                          );
            Gout.y[i] = u;
        }
    }

    // ZERO PADDING
    for(int i=n+1; i<3*n+3; i++ )
    {
        complex<double> u( 0 , 0 );
        Gout.y[i] = u;
    }

    // NEGATIVE FREQUENCIES
    for(int i=3*n+4; i<4*n+4; i++ )
    {
        //double del = 1e-7;
        double xx = (i - 4*n -4) *x/n;

        if(xx < gap)        // IN THE GAP
        {
            double norm =  (  xx + mu  )
                          *(  xx + mu  )
                        +  (   del     )*(   del     );

            complex<double> u(
                             ( xx + mu   )  /norm ,
                            -(   del     )  /norm
                          );
            Gout.y[i] = u;
        }
        else             // OUT OF THE GAP
        {
            double norm =  (  xx + gammaS * xx/( sqrt(xx*xx - gap*gap) ) * ( 1.0 - (gap/xx)  * cos(Pi*phi/2) ) + mu  )
                          *(  xx + gammaS * xx/( sqrt(xx*xx - gap*gap) ) * ( 1.0 - (gap/xx)  * cos(Pi*phi/2) ) + mu  )
                        +  (   del     )*(   del     );

            complex<double> u(
                             ( xx + gammaS * xx/( sqrt(xx*xx - gap*gap) ) * ( 1.0 - (gap/xx)  * cos(Pi*phi/2) ) + mu   )  /norm ,
                            -(   del     )  /norm
                          );
            Gout.y[i] = u;
        }
    }

}


void freePropagatorLorentzShift(double del, double mu, complex<double> Shift, AWT & Gout)
{
    double n = Gout.n;
    double x = Gout.xMax;

    for(int i=0; i<n+1; i++ )
    {
        //double del = 1e-7;
        double xx = i*x/n;

        double norm =  ( xx + mu - real(Shift) ) * ( xx  + mu - real(Shift)  )
                    +  (   del - imag(Shift)   ) * (    del - imag(Shift)    );

        complex<double> u(
                             (   xx + mu  - real(Shift)   )  /norm,
                            -(      del - imag(Shift)     )  /norm
                          );
        Gout.y[i] = u;
    }

    for(int i=n+1; i<3*n+3; i++ )
    {
        complex<double> u( 0 , 0 );
        Gout.y[i] = u;
    }


    for(int i=3*n+4; i<4*n+4; i++ )
    {
        //double del = 1e-7;
        double xx = (i - 4*n -4) *x/n;

        double norm =  ( xx + mu - real(Shift) ) * ( xx  + mu - real(Shift)  )
                    +  (   del + imag(Shift)   ) * (    del + imag(Shift)    );

        complex<double> u(
                             (   xx + mu  - real(Shift)   )  /norm,
                            -(      del - imag(Shift)     )  /norm
                          );
        Gout.y[i] = u;
    }
}



void freeNssPropagatorLorentzShift(double del, double mu, double gap, double gammaS, double phi, complex<double> Shift, AWT & Gout)
{
    double n = Gout.n;
    double x = Gout.xMax;

    double Pi = 3.14159265;

    // SMALL SHIFT TO AVOID SINGULARITIES
    gap = gap - 10e-12;


    // POSITIVE FREQUENCIES
    for(int i=0; i<n+1; i++ )
    {
        //double del = 1e-7;
        double xx = i*x/n;

        if(xx < gap)        // IN THE GAP
        {
            double norm =  ( xx + mu - real(Shift) )
                         * ( xx + mu - real(Shift)  )
                         + (   del - imag(Shift)   ) * (    del - imag(Shift)    );

            complex<double> u(
                                 (   xx + mu  - real(Shift)   )  /norm,
                                -(      del - imag(Shift)     )  /norm
                              );
            Gout.y[i] = u;
        }
        else                // OUT OF THE GAP
        {
            double norm =  ( xx + gammaS * xx/( sqrt(xx*xx - gap*gap) ) * ( 1.0 - (gap/xx)  * cos(Pi*phi/2) ) + mu - real(Shift) )
                         * ( xx + gammaS * xx/( sqrt(xx*xx - gap*gap) ) * ( 1.0 - (gap/xx)  * cos(Pi*phi/2) ) + mu - real(Shift)  )
                         + (   del - imag(Shift)   ) * (    del - imag(Shift)    );

            complex<double> u(
                                 (   xx + gammaS * xx/( sqrt(xx*xx - gap*gap) ) * ( 1.0 - (gap/xx)  * cos(Pi*phi/2) ) + mu  - real(Shift)   )  /norm,
                                -(      del - imag(Shift)     )  /norm
                              );
            Gout.y[i] = u;
        }
    }

    // ZERO PADDING
    for(int i=n+1; i<3*n+3; i++ )
    {
        complex<double> u( 0 , 0 );
        Gout.y[i] = u;
    }

    // NEGATIVE FREQUENCIES
    for(int i=3*n+4; i<4*n+4; i++ )
    {
        //double del = 1e-7;
        double xx = (i - 4*n -4) *x/n;

        if(xx < gap)        // IN THE GAP
        {
            double norm =  ( xx + mu - real(Shift) )
                         * ( xx + mu - real(Shift)  )
                         + (   del - imag(Shift)   ) * (    del - imag(Shift)    );

            complex<double> u(
                                 (   xx + mu  - real(Shift)   )  /norm,
                                -(      del - imag(Shift)     )  /norm
                              );
            Gout.y[i] = u;
        }
        else                // OUT OF THE GAP
        {
            double norm =  ( xx + gammaS * xx/( sqrt(xx*xx - gap*gap) ) * ( 1.0 - (gap/xx)  * cos(Pi*phi/2) ) + mu - real(Shift) )
                         * ( xx + gammaS * xx/( sqrt(xx*xx - gap*gap) ) * ( 1.0 - (gap/xx)  * cos(Pi*phi/2) ) + mu - real(Shift)  )
                         + (   del - imag(Shift)   ) * (    del - imag(Shift)    );

            complex<double> u(
                                 (   xx + gammaS * xx/( sqrt(xx*xx - gap*gap) ) * ( 1.0 - (gap/xx)  * cos(Pi*phi/2) ) + mu  - real(Shift)   )  /norm,
                                -(      del - imag(Shift)     )  /norm
                              );
            Gout.y[i] = u;
        }
    }
}

//////////////////////////////////////////////////
// X function
void Xfunction(double lambda, AWT & Gin, AWT & Phi, AWT & X, AWT & FD, AWT & BE, AWT & K3,
                        AWT & aux1, AWT & aux2, AWT & aux3, AWT & aux4, AWT & aux5, AWT & aux6)
{
    double Pi = 3.14159265359;
    // array of square Gin
    AWT GG;
    GG.initializeAWT(Gin.n, Gin.xMax, Gin.kT);
    for( int iii = 0;         iii<GG.n+1;     iii++ )
    GG.y[iii] = Gin.y[iii] * Gin.y[iii];

    for( int iii = Gin.n+1;   iii<3*GG.n+4;   iii++ )
    GG.y[iii] = 0;

    for( int iii = 3*Gin.n+4; iii<4*GG.n+4;   iii++ )
    GG.y[iii] = Gin.y[iii] * Gin.y[iii];

    // creating and calculating the Kappa function
    AWT Kappa;
    Kappa.initializeAWT(Gin.n, Gin.xMax, Gin.kT);
	Kappa.fermMatsubaraImP(Gin, GG, 1, 1, FD, aux1, aux2, aux3, aux4, aux5, aux6);
    Kappa.KrammersKronig(Kappa,K3, aux1, aux2, aux3, aux4);

    // imaginary part of V
    AWT V;
    V.initializeAWT(Gin.n, Gin.xMax, Gin.kT);
    for( int iii = 0;         iii<Gin.n+1;     iii++ )
    V.y[iii] =  1.0 / ( 1.0 + lambda*Phi.y[iii] )   ;

    for( int iii = Gin.n+1;   iii<3*Gin.n+4;   iii++ )
    V.y[iii] = 0;

    for( int iii = 3*Gin.n+4; iii<4*Gin.n+4;   iii++ )
    V.y[iii] =   1.0 / ( 1.0 + lambda*Phi.y[iii] )   ;

    // imaginary part of K
	AWT K;
	K.initializeAWT(Gin.n, Gin.xMax, Gin.kT);
    for( int iii = 0;         iii<Gin.n+1;     iii++ )
    K.y[iii] = lambda * Kappa.y[iii] * V.y[iii] * V.y[iii];

    for( int iii = Gin.n+1;   iii<3*Gin.n+4;   iii++ )
    K.y[iii] = 0;

    for( int iii = 3*Gin.n+4; iii<4*Gin.n+4;   iii++ )
    K.y[iii] = lambda * Kappa.y[iii] * V.y[iii] * V.y[iii];

    // imaginary part of KV
	AWT KV;
	KV.initializeAWT(Gin.n, Gin.xMax, Gin.kT);

	KV.y[0] = lambda * conj( Kappa.y[0] ) * V.y[0] * V.y[0];

    for( int iii = 1;         iii<Gin.n+1;     iii++ )
    KV.y[iii] = lambda * conj( Kappa.y[4*KV.n+4-iii] ) * V.y[iii] * V.y[iii];

    for( int iii = Gin.n+1;   iii<3*Gin.n+4;   iii++ )
    KV.y[iii] = 0;

    for( int iii = 3*Gin.n+4; iii<4*Gin.n+4;   iii++ )
    KV.y[iii] = lambda * conj( Kappa.y[4*KV.n+4-iii] ) * V.y[iii] * V.y[iii];

    // imaginary part of KV
	AWT KKV;
	KKV.initializeAWT(Gin.n, Gin.xMax, Gin.kT);

	KKV.y[0] = lambda * Kappa.y[0]  * V.y[0] * V.y[0];

    for( int iii = 1;         iii<Gin.n+1;     iii++ )
    KKV.y[iii] = lambda * Kappa.y[4*KV.n+4-iii]  * V.y[iii] * V.y[iii];

    for( int iii = Gin.n+1;   iii<3*Gin.n+4;   iii++ )
    KKV.y[iii] = 0;

    for( int iii = 3*Gin.n+4; iii<4*Gin.n+4;   iii++ )
    KKV.y[iii] = lambda * Kappa.y[4*KV.n+4-iii]  * V.y[iii] * V.y[iii];


    // performing the convolutions
    AWT u1;
    u1.initializeAWT(Gin.n, Gin.xMax, Gin.kT);

    AWT im;
    im.initializeAWT(Gin.n, Gin.xMax, Gin.kT);

    //calculate X11
	AWT X11;
	X11.initializeAWT(Gin.n, Gin.xMax, Gin.kT);
    u1.normalBEtimesIM(V, BE);
    im.deleteReal(GG);
    X11.convolutionAWT(u1, im, 1, -1, aux1, aux2);
    //calculate X21
	AWT X21;
	X21.initializeAWT(Gin.n, Gin.xMax, Gin.kT);
    u1.normalBEtimesIM(K, BE);
    im.deleteReal(Gin);
    X21.convolutionAWT(u1, im, 1, -1, aux1, aux2);

    //calculate X31
    AWT X31;
	X31.initializeAWT(Gin.n, Gin.xMax, Gin.kT);
    u1.normalBEtimesIM(KV, BE);
    im.deleteReal(Gin);
    X31.convolutionAWT(u1, im, 1, -1, aux1, aux2);

    //calculate X12
    AWT X12;
	X12.initializeAWT(Gin.n, Gin.xMax, Gin.kT);
    u1.normalBEtimesIM(GG, BE);
    im.deleteReal(V);
    X12.convolutionAWT(u1, im, -1, -1, aux1, aux2);

    //calculate X22
    AWT X22;
	X22.initializeAWT(Gin.n, Gin.xMax, Gin.kT);
    u1.normalBEtimesIM(Gin, BE);
    im.deleteReal(K);
    X22.convolutionAWT(u1, im, -1, -1, aux1, aux2);

    //calculate X32
    AWT X32;
	X32.initializeAWT(Gin.n, Gin.xMax, Gin.kT);
    u1.normalBEtimesIM(KV, BE);
    im.deleteReal(Gin);
    X32.convolutionAWT(u1, im, -1, -1, aux1, aux2);

    // calculate the whole function X
    for(int i=0;          i < Gin.n+1;   i++)
    {
        X.y[i] =  -X11.y[i] - X12.y[i] + X21.y[i] + X22.y[i] - X31.y[i] - X32.y[i];
        X.y[i] = X.y[i]/Pi;
    }
    for( int i=Gin.n+1;   i<3*Gin.n+4;   i++ )
    X.y[i] = 0;


    for( int i=3*Gin.n+4; i<4*Gin.n+4;   i++ )
    {
        X.y[i] =  -X11.y[i] - X12.y[i] + X21.y[i] + X22.y[i] - X31.y[i] - X32.y[i];
        X.y[i] = X.y[i]/Pi;
    }


    X.KrammersKronig(X, K3, aux1, aux2, aux3, aux4);
}


void kinetic_energy(double & Ekin, double Ehf, double U, double x, double nS, double delta, AWT & SigmaSpec, AWT & FD, AWT & Int)
{
    double Pi = 3.14159265359;


    double factor_1 = FD.xMax/FD.n;
    double factor_2 = U/2.0 - U*nS/2.0;
    double factor_3;
    double factor_4;

    for( int iii = 0;               iii<Int.n+1;    iii++ )
    {
        factor_3 = iii*factor_1 + x + factor_2 - real(SigmaSpec.y[iii]);
        factor_4 = delta - imag(SigmaSpec.y[iii]);
        Int.y[iii] = FD.y[iii] *     ( factor_3 ) / (  factor_3*factor_3  +  factor_4*factor_4 )    ;
    }

    for( int iii = Int.n+1;         iii<3*Int.n+4; iii++ )
    {
        Int.y[iii] = 0    ;
    }

    for( int iii = 3*Int.n+4; iii<4*Int.n+4; iii++ )
    {
        factor_3 = (iii-4*Int.n-4)*factor_1 + x + factor_2 - real(SigmaSpec.y[iii]);
        factor_4 = delta - imag(SigmaSpec.y[iii]);
        Int.y[iii] = FD.y[iii] *     ( factor_3 ) / (  factor_3*factor_3  +  factor_4*factor_4 )    ;
    }

    Ekin = 2.0*delta*Trapezoid_Re(Int)/Pi - Ehf;
}

void correlation_energy(double & Eint, double U, double nT, AWT & GT, AWT & SigmaSpec, AWT & FD, AWT & Int)
{
    double Pi = 3.14159265359;

    for( int iii = 0;            iii<Int.n+1;     iii++ )
    Int.y[iii] =  FD.y[iii] * imag(  GT.y[iii] * SigmaSpec.y[iii]  );

    for( int iii = Int.n+1;         iii<3*Int.n+4; iii++ )
    {
        Int.y[iii] = 0    ;
    }

    for( int iii = 3*Int.n+4;    iii<4*Int.n+4;   iii++ )
    Int.y[iii] =  FD.y[iii] * imag(  GT.y[iii] * SigmaSpec.y[iii]  );

    Eint = -Trapezoid_Re(Int)/Pi + U*(0.5*nT)*(0.5*nT) ;
}

void correlation_energy_alt(double & Eint, double U, double Lambda, double nT, AWT & PHI, AWT & BE, AWT & Int)
{
    double Pi = 3.14159265359;

    Int.y[0] = BE.y[0] * imag(   PHI.y[0] * conj(PHI.y[0]) /  ( 1.0 + Lambda * PHI.y[0] ) );

    for( int i = 1;         i < PHI.n+1;      i++ )
    Int.y[i] = BE.y[4*PHI.n+4-i] * imag(   PHI.y[i] * conj(PHI.y[4*PHI.n+4-i]) / ( 1.0 + Lambda * PHI.y[i]) );

    for( int i = PHI.n+1;   i < 3*PHI.n + 4;  i++ )
    Int.y[i] = 0;

    for( int i = 3*PHI.n+4; i < 4*PHI.n + 4;  i++ )
    Int.y[i] = BE.y[4*PHI.n+4-i] * imag(   PHI.y[i] * conj(PHI.y[4*PHI.n+4-i]) / ( 1.0 + Lambda * PHI.y[i]) );

    Eint = U*(0.5*nT)*(0.5*nT) - U*Lambda*Trapezoid_Re(Int)/Pi ;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////                               Peak analysis                                   ////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////


// Number of peaks in the imaginary part of AWT is determined
int peak_number(AWT & peak, double range)
{
    int peak_number = 0;

    int n = peak.n*ceil(range/peak.xMax);

    int sign[2*n + 1];
    sign[0] = 1;

    double slope;

    for(int i = 1; i < 2*n+1; i++)
    {

        if(i < n)  slope = imag(peak.y[4*peak.n + 4 - n + i ]) - imag(peak.y[4*peak.n + 4 - n + i - 1]);
        if(i == n)  slope = imag(peak.y[          0          ]) - imag(peak.y[      4*peak.n + 3      ]);
        if(i > n)  slope = imag(peak.y[        i - n        ]) - imag(peak.y[        i - n -1        ]);

        if(slope >= 0)  sign[i] = 1;
        if(slope <  0)  sign[i] = -1;
        if(sign[i]*sign[i-1] == -1)     peak_number = peak_number +1;

    }

    return peak_number/2;
}

// position of the left satelite maximum is obtained in the notation of AWT
int left_satelite_n(AWT & peak, double range)
{
    int n = peak.n*ceil(range/peak.xMax);
    int i = 1;
    double slope = -1;

    while(slope < 0)
    {
        slope = imag(peak.y[4*peak.n + 4 - n + i ]) - imag(peak.y[4*peak.n + 4 - n + i - 1]);
        i = i + 1;
        if(i == n -2)
        {
            //cout << "left satrelite missed" << endl;
            slope = 100;
        }
    }

    return 4*peak.n + 4 - n + i;
}

// Energy of the left satelite maximum is obtained
double left_satelite_x(AWT & peak, double range)
{
    int n = peak.n*ceil(range/peak.xMax);
    int i = 1;
    double slope = -1;

    while(slope < 0)
    {
        slope = imag(peak.y[4*peak.n + 4 - n + i ]) - imag(peak.y[4*peak.n + 4 - n + i - 1]);
        i = i + 1;
        if(i == n -2)
        {
            //cout << "left satrelite missed" << endl;
            slope = 100;
        }
    }

    return ( - n + i  )*peak.xMax/peak.n;
}

// position of the edge between the Kondo peak and left satelite is obtained in AWT notation
int left_low_n(AWT & peak, double range)
{
    int n = peak.n*ceil(range/peak.xMax);
    int i = 1;
    double slope = -1;

    while(slope < 0)
    {
        slope = imag(peak.y[4*peak.n + 4 - n + i ]) - imag(peak.y[4*peak.n + 4 - n + i - 1]);
        i = i + 1;
        if(i == n -2)
        {
            //cout << "left satrelite missed" << endl;
            slope = 100;
        }
    }

    while(slope >= 0)
    {
        slope = imag(peak.y[4*peak.n + 4 - n + i ]) - imag(peak.y[4*peak.n + 4 - n + i - 1]);
        i = i + 1;
        if(i == n -2)
        {
            //cout << "left satrelite missed" << endl;
            slope = -100;
        }
    }

    return (  4*peak.n + 4 - n + i   );
}


// Energy of the edge between the Kondo peak and the left satelite is obtained
double left_low_x(AWT & peak, double range)
{
    int n = peak.n*ceil(range/peak.xMax);
    int i = 1;
    double slope = -1;

    while(slope < 0)
    {
        slope = imag(peak.y[4*peak.n + 4 - n + i ]) - imag(peak.y[4*peak.n + 4 - n + i - 1]);
        i = i + 1;
        if(i == n -2)
        {
            //cout << "left satrelite missed" << endl;
            slope = 100;
        }
    }

    while(slope >= 0)
    {
        slope = imag(peak.y[4*peak.n + 4 - n + i ]) - imag(peak.y[4*peak.n + 4 - n + i - 1]);
        i = i + 1;
        if(i == n -2)
        {
            //cout << "left satrelite missed" << endl;
            slope = -100;
        }
    }

    return ( - n + i  )*peak.xMax/peak.n;
}

// the position of the right satelite maximum is obtained
int right_satelite_n(AWT & peak, double range)
{
    int n = peak.n*ceil(range/peak.xMax);
    int i = 1;
    double slope = -1;

    while(slope < 0)
    {
        slope = imag(peak.y[n - i ]) - imag(peak.y[n - i - 1]);
        i = i + 1;
        if(i == n -2)
        {
            //cout << "right satrelite missed" << endl;
            slope = 100;
        }
    }

    return n - i;
}

// Energy of the right satelite maximum is obtained
double right_satelite_x(AWT & peak, double range)
{
    int n = peak.n*ceil(range/peak.xMax);
    int i = 1;
    double slope = -1;

    while(slope < 0)
    {
        slope = imag(peak.y[ n - i - 1]) - imag(peak.y[ n - i]);
        i = i + 1;
        if(i == n -2)
        {
            //cout << "right satrelite missed" << endl;
            slope = 100;

        }
    }

    return (  n - i   )*peak.xMax/peak.n;
}

// position of the edge between the Kondo peak and the right satelite is obtained in AWT notation
int right_low_n(AWT & peak, double range)
{
    int n = peak.n*ceil(range/peak.xMax);
    int i = 1;
    double slope = -1;

    while(slope < 0)
    {
        slope = imag(peak.y[ n - i - 1]) - imag(peak.y[ n - i]);
        i = i + 1;
        if(i == n -2)
        {
            //cout << "right satrelite missed" << endl;
            slope = 100;

        }
    }

    while(slope >= 0)
    {
        slope = imag(peak.y[ n - i - 1]) - imag(peak.y[ n - i]);
        i = i + 1;
        if(i == n -2)
        {
            //cout << "right satrelite missed" << endl;
            slope = -100;

        }
    }
    return (  n - i   );
}

// Energy of the edge between the Kondo peak and the right satelite is obtained
double right_low_x(AWT & peak, double range)
{
    int n = peak.n*ceil(range/peak.xMax);
    int i = 1;
    double slope = -1;

    while(slope < 0)
    {
        slope = imag(peak.y[ n - i - 1]) - imag(peak.y[ n - i]);
        i = i + 1;
        if(i == n -2)
        {
            //cout << "right satrelite missed" << endl;
            slope = 100;

        }
    }

    while(slope >= 0)
    {
        slope = imag(peak.y[ n - i - 1]) - imag(peak.y[ n - i]);
        i = i + 1;
        if(i == n -2)
        {
            //cout << "right satrelite missed" << endl;
            slope = -100;

        }
    }
    return (  n - i   )*peak.xMax/peak.n;
}

// position of the Kondo peak is obtained in AWT notation
int  Kondo_peak_n(AWT & peak, int n_left, int n_right)
{
    int peak_n;
    int i = 1;
    int j = 0;
    double slope = -1;
    while(slope < 0)
    {
        slope = imag(peak.y[ n_left + i ]) - imag(peak.y[n_left + i - 1]);
        i = i + 1;

        if(i == 4*peak.n + 3 - n_left )
        {
            i = 0;
            j = 1;
            while(slope < 0)
            {
                slope = imag(peak.y[ n_right - j - 1]) - imag(peak.y[ n_right - j]);
                j = j + 1;
                if(j == n_right -1)     slope = 100;

            }
        }
    }

    if(j == 0)      peak_n = i + n_left;
    if(i == 0)      peak_n = n_right - j +1;

    return peak_n;
}

// energy of the Kondo peak is obtained
double Kondo_peak_x(AWT & peak, int n_left, int n_right)
{
    double peak_x;
    int i = 1;
    int j = 0;
    double slope = -1;
    while(slope < 0)
    {
        slope = imag(peak.y[ n_left + i ]) - imag(peak.y[n_left + i - 1]);
        i = i + 1;
        if(i == 4*peak.n + 3 - n_left )
        {
            i = 0;
            j = 1;
            while(slope < 0)
            {
                slope = imag(peak.y[ n_right - j - 1]) - imag(peak.y[ n_right - j]);
                j = j + 1;
                if(j == n_right -1)     slope = 100;
            }
        }
    }

    if(j == 0)      peak_x = (-i + 4*peak.n + 4 - n_left)*peak.xMax/peak.n;
    if(i == 0)      peak_x = (n_right - j +1)*peak.xMax/peak.n;

    return peak_x;
}

// Kondo scale analysis
struct half_f
{
    half_f(double y, AWT & z) : half(y), G(z) { }

    int operator()(int position)
    {
        return   abs( half - imag(G.y[position]) );
    }

    double half;
    AWT & G;
};

double hmhw_a(AWT & peak, double range)
{
    double hmhw;

    if(peak_number(peak, range) == 3)
    {

        // position of edges in AWT notation
        int  left_n =  left_low_n(peak, range);
        int right_n = right_low_n(peak, range);

        // position of Kondo peak
        int kondo_n = Kondo_peak_n(peak, left_n, right_n);

        double half_hight = imag( peak.y[kondo_n] )/2;

        int l = 0;
        int p = 0;
        int check_left = 0;
        int check_right = 0;

        // left half hight position is obtained
        if( imag(peak.y[left_n]) > half_hight )
        {
            double slope = -1;
            while(slope < 0)
            {
                slope= imag(peak.y[4*peak.n + 3 - l]) - half_hight;
                l = l + 1;
                if( l > 4*peak.n + 3 - left_n )     slope = 100;
            }
        }
        else
        {
            check_left = 1;
            //cout << "Kondo from left    ";
        }

        // right half hight position is obtained
        if( imag(peak.y[right_n]) > half_hight )
        {
            double slope = -1;
            while(slope < 0)
            {
                slope= imag(peak.y[p]) - half_hight;
                p = p + 1;
                if( p > right_n )     slope = 100;
            }
        }
        else
        {
            check_right = 1;
            //cout << "Kondo from right   ";
        }



        if( (check_left == 0) && (check_right == 0) )
        {
            //cout << "separated Kondo    ";
            hmhw = (l + p -1)*peak.xMax/peak.n;
        }
        else
        {
            hmhw = -1;
        }



    }
    else
    {
        if(peak_number(peak, range) == 2)
        {
            //cout << "Kondo in satelites ";
            hmhw = -1;
        }
        if(peak_number(peak, range) == 1)
        {
            //cout << "Kondo in satelites ";
            hmhw= -1;
        }
    }

    return hmhw;
}

// factor a is just 1 + Lambda* Real of Phi[0]
double factor_a(double Lambda, complex<double> Phi_zero)
{
    return 1 + Lambda*real(Phi_zero);
}

// bethe factor a is from Bethe ansatz
double bethe_a(double U, double Lambda, double x)
{
    double Pi = 3.14159265359;
    return exp(  (U*U/4 - x*x)  );
}

double dosf_a(double U, complex<double> G0)
{
    double Pi = 3.14159265359;
    return exp( U*imag( G0 )/Pi);
}

// calculating the Fermi-liquid quasiparticle weight (residue) Z
pair <double, double> QuasiParticleWeight(AWT & Sigma)
{
    double dE = Sigma.xMax/Sigma.n;
	double LeftSlope = real(Sigma.y[1] - Sigma.y[0])/dE;
	double RightSlope = real(Sigma.y[0] - Sigma.y[4*Sigma.n + 3])/dE;

	pair <double,double> result;
	result.first = 1.0 / ( 1.0 - 0.5*LeftSlope - 0.5*RightSlope);
	result.second = - 0.5*LeftSlope + 0.5*RightSlope;
	return result;
}

////////////////////////////////////////////////////////////////////////////////////////////
/////////////////              Semi-elliptic Propagator              //////////////////////
//////////////////////////////////////////////////////////////////////////////////////////

complex<double> propagatorSemiEl(complex<double> z)
{
    if(z.real()>=0.0)    return(2.0 * (z - sqrt(z*z-1.0)));
    else                 return(2.0 * (z + sqrt(z*z-1.0)));
}

void DensitySemiEl(double del, double mu, AWT & DensityOut, AWT & Sigma)
{
    double n = DensityOut.n;
    double x = DensityOut.xMax;
    //double Pi=3.14159265359;

}

double real_norm(AWT & X, AWT & Y)
{
	double norm = 0;

	if(X.n != Y.n)            cout << "norm of X and Y: not matching size of the meshes!";
	if(X.xMax != Y.xMax)      cout << "norm of X and Y: not matching maximal intervals!";

	for(int i=0;       i<X.n+1;   i++ )
	{
        norm = norm + abs( real(X.y[i]) - real(Y.y[i]) ) / abs( real(Y.y[i]) );
        if(abs( real(Y.y[i]) ) == 0)
        {
            //cout << "norm problem: "  << i << "  " << norm << "  " << real(Y.y[i]) << endl;
        }
    }
    for(int i=3*X.n+4; i<4*X.n+4; i++ )
    {

        norm = norm + abs( real(X.y[i]) - real(Y.y[i]) ) / abs( real(Y.y[i]) );
        if(abs( real(Y.y[i]) ) == 0)
        {
            //cout << "norm problem: "  << i << "  " << norm << "  " << real(Y.y[i]) << endl;
        }
    }

	norm = norm /( 2*X.n + 1);

	return norm;
}

complex<double> shiftedG(int i, int j, AWT & G)
{
    int n;
    int zero_condition = 0;
    complex<double> value;
    complex<double> unity(0,1);

    // in the case one needs to obtain values of Green functions outside AWTs range, we add tails
    if( i > 0        &&  i < G.n+1     &&  j > 0 &&        j < G.n+1    )  n = i+j;
    if( i > 0        &&  i < G.n+1     &&  j > 3*G.n+3 &&  j < 4*G.n+4  )  n = i+j -(4*G.n+4);
    if( i > 3*G.n+3  &&  i < 4*G.n+4   &&  j > 0 &&        j < G.n+1    )  n = i+j -(4*G.n+4);
    if( i > 3*G.n+3  &&  i < 4*G.n+4   &&  j > 3*G.n+3 &&  j < 4*G.n+4  )  n = i+j -(8*G.n+8);


    if( i < 3*G.n+4  &&  i > G.n) zero_condition = 1;
    if( j < 3*G.n+4  &&  j > G.n) zero_condition = 1;

    double pomer =  ( n / G.n );

    if(zero_condition == 0)
    {
        if(n > G.n+1                )    value = real(G.y[G.n]) / pomer / pomer + unity * imag(G.y[G.n]) / pomer ;
        if(n >-1      && n < G.n+1  )    value = G.y[n];
        if(n >-G.n-1  && n < 0      )    value = G.y[4*G.n+4-n];
        if(              n < -2*G.n )    value = real(G.y[G.n]) / pomer / pomer + unity * imag(G.y[G.n]) / pomer ;
    }
    else
    {
        value = 0;
    }

    return value;
}
