#include "EffectiveT0.h"
#include "physics.h"

//#include "ArrayWithTails.h"


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
//#include <boost/lexical_cast.hpp>


using namespace std;


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////              ZERO T PARQUETS EFFECTIVE INTERACTION APPROX                         ///////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void Effective_zeroT(string input)
{
    // intialization of model variables
    double delta;
    double mu;
    string model;

    double kT_in;
    double kT_increment;
    double kT_range;

    double U_in;
    double U_increment;
    double U_range;

    double x_increment;
    double x_range;

    // intialization of variables for mesh properties
    double xMax;
    int n;
    // intialization of variables for printing out txt files
    int display;
    double range;
    int print;

    string export_mode, physics;

    // previously initialized variables are set to values from the input file
    importInitial(input, delta, mu, model,
                  kT_in, kT_increment, kT_range,
                  U_in, U_increment, U_range,
                  x_increment, x_range,
                  xMax, n, display, range,
                  print, export_mode, physics);

    AWT aux1;
    aux1.initializeAWT(n, xMax, 0);
    AWT aux2;
    aux2.initializeAWT(n, xMax, 0);
    AWT aux3;
    aux3.initializeAWT(n, xMax, 0);
    AWT aux4;
    aux4.initializeAWT(n, xMax, 0);
    AWT aux5;
    aux3.initializeAWT(n, xMax, 0);
    AWT aux6;
    aux4.initializeAWT(n, xMax, 0);


    // this program ignores any non-zero T input!
    double kT = 0;
    double U = U_in;

    // variable for maximal number of x iteration is initialized and set
    int x_iterations = 2*ceil(x_range)/x_increment;

    /////////////////////////////////////////////////////////////////
    ///////                 Basic arrays               /////////////
    ///////////////////////////////////////////////////////////////

    // creating array with Kernel3
    AWT K3;
    K3.initializeAWT(n, xMax, kT);
    K3.set_K3();

    // creating array with Fermi-Dirac statistics
    AWT FD;
    FD.initializeAWT(n, xMax, kT);
    FD.set_FD();

    // creating array with Bose-Einstein statistics
    AWT BE;
    BE.initializeAWT(n, xMax, kT);
    BE.set_BE();

    cout << "K3, FD, BE done" << endl;

    /////////////////////////////////////////////////////////////////
    /////      Defs of quantities used in the program        ///////
    ///////////////////////////////////////////////////////////////

    // thermodynamical self-energy
    AWT SigmaTherm;
    SigmaTherm.initializeAWT(n, xMax, kT);
    SigmaTherm.set_zero();

    // spectral self-energy
    AWT SigmaSpec;
    SigmaSpec.initializeAWT(n, xMax, kT);
    SigmaSpec.set_zero();

    // thermodynamic propagator
    AWT Gtherm;
    Gtherm.initializeAWT(n, xMax, kT);
    Gtherm.set_zero();

    // spectral propagator
    AWT Gspec;
    Gspec.initializeAWT(n, xMax, kT);
    Gspec.set_zero();

    // thermodynamic DOS
    AWT DOS;
    DOS.initializeAWT(n, xMax, kT);
    DOS.set_zero();

    // bubble
    AWT Phi;
    Phi.initializeAWT(Gtherm.n, Gtherm.xMax, Gtherm.kT);
    Phi.set_zero();

    // Kvertex
    AWT Kvertex;
    Kvertex.initializeAWT(n, xMax, kT);
    Kvertex.set_zero();

    // static bubble Psi
    double Psi;

    // lambda
    double Lambda;
    double LambdaIter[x_iterations];
    // thermodynamic occupation number
    double nT;
    double nTiter[x_iterations];
    // spectral occupation number
    double nS;
    double nSiter[x_iterations];

    cout << "Sigmas, Gs, DOS, Phi, Kverex, Lambda, Psi, nT, nS initialized" << endl;

    /////////////////////////////////////////////////////////////////
    //   non-interacting half-filling Propagators and bubbles    ///
    ///////////////////////////////////////////////////////////////

    freePropagatorLorentz(delta, 0.0, Gtherm);

    Phi.fermMatsubaraImP(Gtherm, Gtherm, 1, 1, FD, aux1, aux2, aux3, aux4, aux5, aux6);
    Phi.KrammersKronig(Phi,K3, aux1, aux2, aux3, aux4);

    //Phi.exportAWTasFUN(boost::lexical_cast<string>("Txt/Phi"), display, range, 1);

    /////////////////////////////////////////////////////////////////
    ///////              Half filling case         /////////////////
    ///////////////////////////////////////////////////////////////

    // half filling value of Lambda, Psi, Kvertex, nT
    cout << "Psi, Lambda, Kvertex in half-filling is calculated" << endl;
    Lambda_zeroT(U, Gtherm, Phi, Lambda, Kvertex, Psi, FD, BE, K3, aux1, aux2, aux3, aux4, aux5, aux6);

    cout << "non-interacting Lambda at half-filling is:  " << Lambda << endl;
    cout << "non-interacting Psi at half-filling is:     " << Psi    << endl;

    double LambdaHalf = Lambda;

    for(int i = 0; i < x_iterations +1; i++)
    {

        // the parameter x is set for the out of the half filling case
        double x;
        // chemical potential for the lattice model
        if(model == "lattice")       x =  mu - U/2;       // !!!!!!!!!
        // chemical potential for the lattice model
        if(model == "lorentz")       x = ceil(x_range) - i*x_increment;

        // half filling values of nT, Lambda and SigmaTherm
        nT = 1.0;
        Lambda = LambdaHalf;
        SigmaTherm.set_zero();

        /////////////////////////////////////////////////////////////////
        ////              SigmaTherm and Gtherm                 ////////
        ////          outside of half filling case             ////////
        //////////////////////////////////////////////////////////////

        // nT_ is used to store values of nT from previous iteration
        double nT_ = 1e3;
        // Lambda_ is used to store values of nT from previous iteration
        double Lambda_ = 1e3;

        // set the precision
        double nTprec = 0.0001;
        double lambPrec = 0.0001;

        // iteration for SigmaTherm
        while( (abs( nT - nT_) > nTprec) && (abs( Lambda - Lambda_) > lambPrec) )
        {
            // values from previous iteration are stored in nT_, Lambda_
            nT_ = nT;
            Lambda_ = Lambda;

            // nT is given by an implicit equation 24 and solved with Brent
            //DensityLorentz(delta, mu, DOS, SigmaTherm);
            pair<double, double> nTbrent = Brent_nT_zeroT(delta, x, Lambda, FD, BE, K3);
            nT = nTbrent.first;

            // Gtherm in given iter calculated from previous SigmaTherm
            freePropagatorLorentz(delta, x - 0.5 * Lambda * ( nT - 1 ), Gtherm);

            // Phi in given iter calculated from previous SigmaTherm
            PhiFunction(Gtherm, Phi, FD, K3, aux1, aux2, aux3, aux4, aux5, aux6);

            // Lambda, Kvertex, Psi are recalculated
            Lambda_zeroT(U, Gtherm, Phi, Lambda, Kvertex, Psi, FD, BE, K3, aux1, aux2, aux3, aux4, aux5, aux6);

        }

        string str;
        string fileName;

        freePropagatorLorentz(delta, x - 0.5 * Lambda * ( nT - 1 ), Gtherm);

        //str = boost::lexical_cast<std::string>(x);
        fileName = "Txt/Gtherm_" + str;
        //Gtherm.exportAWTasFUN(fileName, display, range, 1);

        // using Eq. 30, SigmaSpec is calculated from Phi, Lambda and Gtherm
        Phi.fermMatsubaraImP(Gtherm, Gtherm, 1, 1, FD, aux1, aux2, aux3, aux4, aux5, aux6);
        Phi.KrammersKronig(Phi,K3, aux1, aux2, aux3, aux4);

        /////////////////////////////////////////////////////////////////
        ////      SigmaSpec outside of half filling case       /////////
        ///////////////////////////////////////////////////////////////

        AWT Integrand;
        Integrand.initializeAWT(n, xMax, kT);

        for( int iii = 0;               iii<Integrand.n+1;     iii++ )
        Integrand.y[iii] = -U * Lambda * Phi.y[iii] / ( 1.0 + Lambda * Phi.y[iii]);

        for( int iii = Integrand.n+1;   iii<3*Integrand.n+4; iii++ )
        Integrand.y[iii] = 0;

        for( int iii = 3*Integrand.n+4; iii<4*Integrand.n+4; iii++ )
        Integrand.y[iii] = -U * Lambda * Phi.y[iii] / ( 1.0 + Lambda * Phi.y[iii]);


        SigmaSpec.boseMatsubaraImP(Integrand, Gtherm, 1, 1, FD, BE, aux1, aux2, aux3, aux4, aux5, aux6);
        SigmaSpec.KrammersKronig(SigmaSpec,K3, aux1, aux2, aux3, aux4);

        //str = boost::lexical_cast<std::string>(x);
        fileName = "Txt/SigmaSpec_" + str;
        //SigmaSpec.exportAWTasFUN(fileName, display, range, 1);


        /////////////////////////////////////////////////////////////////
        ////        Gspec outside of half filling case         /////////
        ///////////////////////////////////////////////////////////////

        // using Brent nS is determined
        pair<double, double> nSbrent = Brent_nS_zeroT(delta, x, U, SigmaSpec, Gspec, FD, BE, K3);
        nS = nSbrent.first;

        propagatorLorentz(delta, x -0.5*U*nS + 0.5*U, Gspec, SigmaSpec);
        fileName = "Txt/Gspec_" + str;
        //Gspec.exportAWTasFUN(fileName, display, range, 1);

        // resultes are written or prapered to be written to the outputs
        cout << "x is: " << x << "  nT is: " << nT << "   nS is: " << nS <<  "  Lambda is:  " << Lambda  << endl;
        nTiter[i] = nT;
        nSiter[i] = nS;
        LambdaIter[i] = Lambda;
   }

    string out = "Txt/output";
    ofstream output(out.c_str());
    for(int i = 0; i < x_iterations +1; i++)
    output << ceil(x_range) - i*x_increment << "   " <<   nTiter[i]   <<  "  " << nSiter[i] << "  " << LambdaIter[i] << endl;

}


void Lambda_zeroT(double U, AWT & Gin, AWT & Phi, double & Lambda, AWT & Kvertex, double & Psi, AWT & f_FD, AWT & f_BE, AWT & K3,
                                AWT & aux1, AWT & aux2, AWT & aux3, AWT & aux4, AWT & aux5, AWT & aux6 )
{
    // the calculated lambda is stored
    Lambda = Brent_Lambda_zeroT(U, Gin, Phi, f_FD, f_BE, K3, aux1, aux2, aux3, aux4, aux5, aux6);
    // the calculated Kvertex is stored
    KvertexFunction(Lambda, Phi, Kvertex);
    // the calculated Psi is stored
    PsiFunction(Gin, Kvertex, f_BE, Psi);
}


struct lambda_implicit
{
    lambda_implicit(double y, AWT & z1, AWT & z2, AWT & z3, AWT & z4, AWT & z5, AWT & z6, AWT & z7, AWT & z8, AWT & z9, AWT & z10) :
    U(y), Gin(z1), f_FD(z2), f_BE(z3), K3(z4), aux1(z5), aux2(z6), aux3(z7), aux4(z8), aux5(z9), aux6(z10) { }

    double operator()(double lambda)
    {
        double Pi=3.14159265359;

        AWT Phi;
        Phi.initializeAWT(Gin.n, Gin.xMax, Gin.kT);
        Phi.fermMatsubaraImP(Gin, Gin, 1, 1, f_FD, aux1, aux2, aux3, aux4, aux5, aux6);
        Phi.KrammersKronig(Phi,K3, aux1, aux2, aux3, aux4);

        AWT Kvertex;
        Kvertex.initializeAWT(Gin.n, Gin.xMax, Gin.kT);
        for(int i=0; i<4*Gin.n+4; i++)   Kvertex.y[i] = -1.0 * lambda * lambda * Phi.y[i] / ( 1.0 + lambda * Phi.y[i] );

        // the static bubble psi
        AWT Integrand;
        Integrand.initializeAWT(Gin.n, Gin.xMax, Gin.kT);
        Integrand.y[0] = f_BE.y[0] * imag(   Gin.y[0] * conj(Gin.y[0]) * conj(Kvertex.y[0])  )/ Pi;
        for( int i = 1;         i < Gin.n+1;      i++ )
        Integrand.y[i] = f_BE.y[i] * imag(   Gin.y[i] * conj(Gin.y[4*Gin.n+4-i]) * conj(Kvertex.y[4*Gin.n+4-i])  )/Pi;
        for( int i = Gin.n+1;   i < 3*Gin.n + 4;  i++ )
        Integrand.y[i] = 0;
        for( int i = 3*Gin.n+4; i < 4*Gin.n + 4;  i++ )
        Integrand.y[i] = f_BE.y[i] * imag(   Gin.y[i] * conj(Gin.y[4*Gin.n+4-i]) * conj(Kvertex.y[4*Gin.n+4-i])  )/Pi;
        double Psi = integrateReP(Integrand);

        return  abs( lambda - U / ( 1.0 + Psi ) );
    }

    double U;
    AWT & Gin;
    AWT & f_FD;
    AWT & f_BE;
    AWT & K3;
    AWT & aux1, aux2, aux3, aux4, aux5, aux6;
};

double Brent_Lambda_zeroT(double U, AWT & G, AWT & Phi, AWT & FD, AWT & BE, AWT & K3, AWT & aux1, AWT & aux2, AWT & aux3, AWT & aux4, AWT & aux5, AWT & aux6 )
{
    using boost::math::tools::brent_find_minima;
    int bits = std::numeric_limits<float>::digits;

    //cout << "Lambda numeric precision set as " << bits << endl;

    double Lmin = 0.0;
    //cout << "making Lmin:    " << Lmin << endl;
    double Lmax = -1.0/real( Phi.y[0]) - 1e-10;
    //cout << "making Lmax:    " << Lmax << endl;

    pair<double, double> Lbrent = brent_find_minima( lambda_implicit(U, G, FD, BE, K3, aux1, aux2, aux3, aux4, aux5, aux6), Lmin, Lmax, bits);
    return Lbrent.first;
}

void drawLambda(double U, rawF & lambdaFile, AWT & Gin, AWT & Phi, AWT & f_FD, AWT & f_BE, AWT & K3)
{
        double Pi=3.14159265359;

        AWT Kvertex;
        Kvertex.initializeAWT(Gin.n, Gin.xMax, Gin.kT);

        double Psi;
        double lambda;

        for(int j = 0; j < lambdaFile.n; j++)
        {
            // actual lambda value
            lambda = j * Pi/lambdaFile.n;
            // Kvertex
            KvertexFunction(lambda, Phi,Kvertex);
            // the static bubble psi
            PsiFunction(Gin, Kvertex,f_BE, Psi);

            lambdaFile.f[lambdaFile.n+j] = abs( lambda - U / (1 + Psi) );
        }

}



struct nT_implicit
{
    nT_implicit(double y1, double y2, double y3, AWT & z3, AWT & z4, AWT & z5) : delta(y1), x(y2), Lambda(y3), f_FD(z3), f_BE(z4), K3(z5) { }

    double operator()(double nT)
    {
        double Pi=3.14159265359;
        //AWT DOS;
        //DOS.initializeAWT(DOS.n, DOS.xMax, DOS.kT);
        //DensityLorentz(delta, mu, DOS, St);

        // the static bubble psi
        //AWT Integrand;
        //Integrand.initializeAWT(DOS.n, DOS.xMax, DOS.kT);
        //for( int i = 0; i < 4*DOS.n + 4;  i++ )  Integrand.y[i] = DOS.y[i];

        return  abs( - nT + 1.0 +  (2.0/Pi)*atan( (x + 0.5*Lambda - 0.5*Lambda*nT )/delta )  );
    }

    double delta;
    double x;
    double Lambda;
    AWT & f_FD;
    AWT & f_BE;
    AWT & K3;
};

pair<double, double> Brent_nT_zeroT(double delta, double x, double Lambda, AWT & FD, AWT & BE, AWT & K3)
{
    using boost::math::tools::brent_find_minima;
    int bits = std::numeric_limits<float>::digits;
    //cout << "numeric precision set as " << bits << endl;

    double nMin = 0.0;
    double nMax = 2.0;

    return brent_find_minima( nT_implicit(delta, x, Lambda, FD, BE, K3), nMin, nMax, bits);
}

struct nS_implicit
{
    nS_implicit(double y1, double y2, double y3, AWT & z1, AWT & z2, AWT & z3, AWT & z4, AWT & z5) : delta(y1), x(y2), U(y3), SigmaSpec(z1), Gspec(z2), f_FD(z3), f_BE(z4), K3(z5) { }

    double operator()(double nS)
    {
        double Pi=3.14159265359;
        // the static bubble psi

        AWT help;
        help.initializeAWT(Gspec.n, Gspec.xMax, Gspec.kT);
        propagatorLorentz(delta, x-0.5*U*nS+0.5*U , help, SigmaSpec);

        AWT Integrand;
        Integrand.initializeAWT(Gspec.n, Gspec.xMax, Gspec.kT);
        for( int i=0;  i<4*Gspec.n+4; i++ )   Integrand.y[i] = -(2/Pi) * help.y[i] * f_FD.y[i];

        return  abs( - nS + integrateImP(Integrand) );
    }

    double delta;
    double x;
    double U;
    AWT & SigmaSpec;
    AWT & Gspec;
    AWT & f_FD;
    AWT & f_BE;
    AWT & K3;
};

pair<double, double> Brent_nS_zeroT(double delta, double x, double U, AWT & SigmaSpec, AWT & Gspec, AWT & FD, AWT & BE, AWT & K3)
{
    using boost::math::tools::brent_find_minima;
    int bits = std::numeric_limits<float>::digits;
    //cout << "numeric precision set as " << bits << endl;

    double nMin = 0.0;
    double nMax = 2.0;

    return brent_find_minima( nS_implicit(delta, x, U, SigmaSpec, Gspec, FD, BE, K3), nMin, nMax, bits);
}


///////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////                  PARQUETS LOW FREQ                          /////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////

void SolverParquetLowFreq(double U, AWT & G, AWT & Sigma, AWT & f_FD, AWT & f_BE, AWT & K3)
{
    // G times G
    AWT GG;
    GG.initializeAWT(G.n, G.xMax, G.kT);
    // electron-hole bubble
    AWT Chi;
    Chi.initializeAWT(G.n, G.xMax, G.kT);
    // electron-electron bubble
    AWT Psi;
    Psi.initializeAWT(G.n, G.xMax, G.kT);

    // singular part of ee vertex, Eq. (3.17)
    AWT Lsing;
    Lsing.initializeAWT(G.n, G.xMax, G.kT);
    // Equation 3.21
    AWT L;
    L.initializeAWT(G.n, G.xMax, G.kT);
    // temporary array
    AWT aux;
    aux.initializeAWT(G.n, G.xMax, G.kT);

    AWT aux1;
    aux1.initializeAWT(G.n, G.xMax, 0);
    AWT aux2;
    aux2.initializeAWT(G.n, G.xMax, 0);
    AWT aux3;
    aux3.initializeAWT(G.n, G.xMax, 0);
    AWT aux4;
    aux4.initializeAWT(G.n, G.xMax, 0);
    AWT aux5;
    aux5.initializeAWT(G.n, G.xMax, 0);
    AWT aux6;
    aux6.initializeAWT(G.n, G.xMax, 0);

    // ansatz of effec, Ubar
    double Ubar;

    // 1. Calculate GG
    GG.doubleInverse(G);

    // 2. Calculate Chi, Psi
    Chi.fermMatsubaraImP(G, G, 1, 1, f_FD, aux1, aux2, aux3, aux4, aux5, aux6);
    Chi.KrammersKronig(Chi, K3, aux1, aux2, aux3, aux4);

    Psi.fermMatsubaraImP(G, G, 1, -1, f_FD, aux1, aux2, aux3, aux4, aux5, aux6);
    Psi.KrammersKronig(Psi, K3, aux1, aux2, aux3, aux4);

    // 3. finding of Ubar
    double min = 0.0;                                   // zero temperature
    double max = - 1.0 /real(Chi.y[0]) + 1e-7;;         // crtical temperature Uc, Equation 3.15a
    double tempU;                                       // temporary U

    do
    {
        // preliminary Ubar is set
        Ubar = (min + max) / 2.0;

        // singular part of ee vertex is calculated according to Eq. (3.17)
        for(int i=0; i<=Lsing.nn; i++)      Lsing.y[i] = - Ubar * Ubar * Chi.y[i] / (1.0 + Ubar * Chi.y[i]);

        // Equation 3.21 is evaluated
        L.boseMatsubaraImP(Lsing, GG, -1, -1, f_FD, f_BE, aux1, aux2, aux3, aux4, aux5, aux6);
        L.KrammersKronig(L, K3, aux1, aux2, aux3, aux4);

        tempU = Ubar + 1.0 / real(Psi.y[0]) * real(L.y[0])*real(L.y[0])/ ( 1.0 + real(L.y[0]) );

        if(tempU > U)   max = Ubar;
        else            min = Ubar;
    }
    while(abs(U-tempU) > 0.0001);

    //cout << "a, Ubar = " << - 1.0 / real(Chi.y[0])  - Ubar << "  " << Ubar << endl;

    for(int i=0; i<=Lsing.nn; i++)     aux.y[i] = Lsing.y[i];

    // 4. Evaluate Sigma
    Sigma.boseMatsubaraImP(aux, G, 1, 1, f_FD, f_BE, aux1, aux2, aux3, aux4, aux5, aux6);
    Sigma.KrammersKronig(Sigma, K3, aux1, aux2, aux3, aux4);
}

void ParquetLowFreq(int display, double range, int iterParquets, double U, double delta, double mu, AWT & G, AWT & Sigma, AWT & FD, AWT & BE, AWT & K3)
{
        double Pi=3.14159265359;
        for(int iter=0; iter< iterParquets; iter++)
        {
            SolverParquetLowFreq(U, G, Sigma, FD, BE, K3);
            //propagatorLorentz(delta, mu, G, Sigma);

            string name;
            //name = "Txt/R" + boost::lexical_cast<string>(iter);
            //G.exportAWTasFUN(name, display, range, -1/Pi);

            //name = "Txt/S" + boost::lexical_cast<string>(iter);
            //Sigma.exportAWTasFUN(name, display, range, 1);
        }

}

//////////////////////////////////////////////////////////////////////////////////////
//////////////                       FLEX                            ////////////////
////////////////////////////////////////////////////////////////////////////////////

/*
void Flex(double xMax, double kT, int n, int display, double range, double incrementU, double delta, double seterror, AWT & FD, AWT & BE, AWT & K3)
{
    double Pi=3.14159265359;
    cout << endl;
    cout << "FLEX RUNNING" <<endl;

    // initialize values for zeroth iteration
    double U=0;
    int iter = 4/incrementU;

    // Sigma stores self-energy for given itteration of the calculation
    AWT Sigma;
    Sigma.initializeAWT(n, xMax, kT);
    Sigma.setZero();

    // Gprev stores G computed in the instance of SolverFlex called one itteration back
    AWT Gprev;
    Gprev.initializeAWT(n, xMax, kT);
    //propagatorLorentz(delta, 0, Gprev, Sigma);

    // Gprev stores G computed in the instance of SolverFlex called one itteration back
    AWT Gnow;
    Gnow.initializeAWT(n, xMax, kT);
    //propagatorLorentz(delta, 0, Gnow, Sigma);

    // Chi
    AWT Chi;
    Chi.initializeAWT(n, xMax, kT);
    Chi.fermMatsubaraImP(Gnow, Gnow, 1, 1, FD);
    Chi.KrammersKronig(Chi,K3);

    // exporting data for U = 0
    string name = "Txt/R0";
    //Gprev.exportAWTasFUN(name, display, range, -1/Pi);
    name = "Txt/S0";
    //Sigma.exportAWTasFUN(name, display, range, 1);
    name = "Txt/X0";
    //Chi.exportAWTasFUN(name, display, range, 1);

    cout << endl;
    cout << "U=0 finished" << endl;
    cout << endl;

    double error;

    // setting a variable to count iterations corresponding to U = 1, 2, 3, 4
    int counter=1;

    // iterations for ever higher U
    for(int j=0; j<iter; j++)
    {
        int jj =0;
        U = U +incrementU;

        do
        {
            // cout << "U: " << U << endl;
            SolverFLEX(U, delta, Sigma, Chi, Gnow, FD, BE, K3);
            //cout << "sigma: " << Sigma.y[0] << endl;
            error = NormL2(Gprev, Gnow);
            cout  << "U = " << U << "  L2 diference = " << error << endl;
            Gprev.loadAWTtoAWT(Gnow);

            string name;

            //name = "Txt/S" + DoubleToStr(U) + "_" + IntToStr(jj);
            //Sigma.asymmetryTest(name, -1, 1, display, range);
            //Sigma.exportAWTasFUN(name, display, range, 1);
            //name = "Txt/X" + DoubleToStr(U) + "_" + IntToStr(jj);
            //Chi.asymmetryTest(name, 1, -1, display, range);
            //Chi.exportAWTasFUN(name, display, range, 1);


            jj = jj+1;
            if(jj==1000) exit(0);
        }
        while(error > seterror);

        // module to store results for U = 1, 2, 3, 4
        if( (j +1 == 1/incrementU)  || (j +1 == 2/incrementU)  || (j +1 == 3/incrementU)  || (j +1 == 4/incrementU)   )
        {
            cout << "writting down R, S, X" << endl;
            string name;
            //name = "Txt/R" + DoubleToStr(counter);
            //Gnow.exportAWTasFUN(name, display, range, -1/Pi);
            //name = "Txt/S" + DoubleToStr(counter);
            //Sigma.exportAWTasFUN(name, display, range, 1);
            //name = "Txt/X" + DoubleToStr(counter);
            //Chi.exportAWTasFUN(name, display, range, 1);
            cout << endl;
            counter = counter +1;
        }

        Gprev.loadAWTtoAWT(Gnow);
    }
}

void SolverFLEX(double U, double delta, AWT & Sigma, AWT & Chi, AWT & Gin, AWT & f_FD, AWT & f_BE, AWT & K3)
{
    // 1. Calculate Chi
    Chi.fermMatsubaraImP(Gin, Gin, 1, 1, f_FD);
    Chi.KrammersKronig(Chi, K3);                            // real part of Chi is calculated from Kramers-Kronig

    // 2. determine the C function
    AWT C;
    C.initializeAWT(Gin.n, Gin.xMax, Gin.kT);
    for(int ii=0; ii< 4*Gin.n + 4; ii++)  C.y[ii] = -0.5*U * U * Chi.y[ii]/(1.0 + U * Chi.y[ii]);

    // 3. determine the self-energy
    Sigma.boseMatsubaraImP(C, Gin, 1, 1, f_FD, f_BE);
    Sigma.KrammersKronig(Sigma, K3);                            // real part of Sigma is calculated from Kramers-Kronig
    //cout << "SolverFlex, sigma: " << Sigma.y[0] << endl;

    //propagatorLorentz(delta, 0, Gin, Sigma);
    Gin.KrammersKronig(Gin, K3);
}

*/
