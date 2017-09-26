#include "Methods.h"
#include "physics.h"

#include "Constants.h"

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
#include <boost/lexical_cast.hpp>
#include <boost/filesystem.hpp>
#include <iomanip> // setprecision

using namespace std;




struct lambda_implicit_non_zero_T
{
    lambda_implicit_non_zero_T(double y, AWT & z1, AWT & z2, AWT & z3, AWT & z4, AWT & z5, AWT & z6, AWT & z7) :
                                        U(y), GT(z1), GTshift(z2), PHI(z3), PHIshift(z4), FD(z5), BE(z6), K3(z7) { }

    double operator()(double lambda)
    {
        AWT Integrand;
        Integrand.initializeAWT(GT.n, GT.xMax, GT.kT);

        // first integral of Eq. 19 is performed
        Integrand.y[0] = FD.y[0] * real( PHIshift.y[0] / (1.0 + lambda * PHIshift.y[0]  ) ) * imag( GT.y[0] * conj( GT.y[0] ) );

        for( int i = 1;         i < GT.n+1;      i++ )
        Integrand.y[i] = 0;

        for( int i = 3*GT.n+4; i < 4*GT.n + 4;  i++ )
        Integrand.y[i] = FD.y[i] * real( PHIshift.y[4*GT.n + 4 - i] / (1.0 + lambda * PHIshift.y[4*GT.n + 4 - i] ) ) * imag( GT.y[i] * conj( GT.y[4*GT.n + 4 - i] ) );

        // first part is stored into Psi
        double Psi = (lambda*lambda/Pi)*integrateReP(Integrand);

        // second integral of Eq. 19 is performed
        Integrand.y[0] = BE.y[0] * imag( conj(PHI.y[0]) / (1.0 + lambda * conj(PHI.y[0]) ) ) * real( conj(GTshift.y[0]) * GTshift.y[0]  );

        for( int i = 1;         i < GT.n+1;      i++ )
        Integrand.y[i] = 0;

        for( int i = 3*GT.n+4; i < 4*GT.n + 4;  i++ )
        Integrand.y[i] = BE.y[i] * imag( conj(PHI.y[4*GT.n + 4 - i]) / (1.0 + lambda * conj(PHI.y[4*GT.n + 4 - i]) ) ) * real( GTshift.y[i] * conj( GTshift.y[4*GT.n + 4 - i] ) );

        // Psi is calculated from Eq. 19
        Psi = Psi - (lambda*lambda/Pi)*integrateReP(Integrand);

        return  abs( lambda - U / ( 1.0 + Psi ) );
    }

    double U;
    AWT & GT;
    AWT & GTshift;
    AWT & PHI;
    AWT & PHIshift;
    AWT & FD;
    AWT & BE;
    AWT & K3;
};


struct lambda_implicit_non_zero_T_simple
{
    lambda_implicit_non_zero_T_simple(double y, AWT & z1, AWT & z3, AWT & z5, AWT & z6, AWT & z7) :
                                        U(y), GT(z1), PHI(z3), FD(z5), BE(z6), K3(z7) { }

    double operator()(double lambda)
    {
        AWT Integrand;
        Integrand.initializeAWT(GT.n, GT.xMax, GT.kT);

        // first integral of Eq. 19 is performed
        Integrand.y[0] = FD.y[0] * real( conj(PHI.y[0]) / (1.0 + lambda * conj(PHI.y[0]) ) ) * imag( GT.y[0] * conj( GT.y[0] ) );

        for( int i = 1;         i < GT.n+1;      i++ )
        Integrand.y[i] = 0;

        for( int i = 3*GT.n+4; i < 4*GT.n + 4;  i++ )
        Integrand.y[i] = FD.y[i] * real( PHI.y[4*GT.n + 4 - i] / (1.0 + lambda * PHI.y[4*GT.n + 4 - i] ) ) * imag( GT.y[i] * conj( GT.y[4*GT.n + 4 - i] ) );

        // first part is stored into Psi
        double Psi = (lambda*lambda/Pi)*integrateReP(Integrand);

        // second integral of Eq. 19 is performed
        Integrand.y[0] = BE.y[0] * imag( conj(PHI.y[0]) / (1.0 + lambda * conj(PHI.y[0]) ) ) * real( GT.y[0] * conj(GT.y[0])   );

        for( int i = 1;         i < GT.n+1;      i++ )
        Integrand.y[i] = 0;

        for( int i = 3*GT.n+4; i < 4*GT.n + 4;  i++ )
        Integrand.y[i] = BE.y[i] * imag( conj(PHI.y[4*GT.n + 4 - i]) / (1.0 + lambda * conj(PHI.y[4*GT.n + 4 - i]) ) ) * real( GT.y[i] * conj( GT.y[4*GT.n + 4 - i] ) );

        // Psi is calculated from Eq. 19
        Psi = Psi - (lambda*lambda/Pi)*integrateReP(Integrand);

        return  abs( lambda - U / ( 1.0 + Psi ) );
    }

    double U;
    AWT & GT;
    AWT & PHI;
    AWT & FD;
    AWT & BE;
    AWT & K3;
};

void Lambda_nonZeroT(double U, double Lmin, string physics,
                     AWT & GT, AWT & GTshift, AWT & PHI, AWT & PHIshift,
                     double & Lambda, double & Psi,
                     AWT & f_FD, AWT & f_BE, AWT & K3)
{
    // the minimal value of Lambda is passed in, the maximal value Lmax is determined from PHI
    //double Lmax = -1/real( PHI.y[0]) - 1e-10;;
    //cout << "making Lmax:    " << Lmax << endl;


    using boost::math::tools::brent_find_minima;
    int bits = std::numeric_limits<float>::digits;

    double Lmax = -1/real( PHI.y[0]) - 1e-10;

    AWT Integrand;
    Integrand.initializeAWT(GT.n, GT.xMax, GT.kT);

    if(physics == "EFF" || physics == "EFF1" || physics == "EFF2" )
    {
        pair<double, double> Lbrent = brent_find_minima( lambda_implicit_non_zero_T(U, GT, GTshift, PHI, PHIshift, f_FD, f_BE, K3), Lmin, Lmax, bits);
        Lambda = Lbrent.first;

        // first integral of Eq. 19 is performed
        Integrand.y[0] = f_FD.y[0] * real( PHIshift.y[0] / (1.0 + Lambda * PHIshift.y[0] ) ) * imag( GT.y[0] * conj( GT.y[0] ) );

        for( int i = 1; i < GT.n+1; i++ )
        Integrand.y[i] = f_FD.y[i] * real( PHIshift.y[4*GT.n + 4 - i] / (1.0 + Lambda * PHIshift.y[4*GT.n + 4 - i] ) ) * imag( GT.y[i] * conj( GT.y[4*GT.n + 4 - i] ) );

        for( int i = 3*GT.n+4; i < 4*GT.n + 4;  i++ )
        Integrand.y[i] = f_FD.y[i] * real( PHIshift.y[4*GT.n + 4 - i] / (1.0 + Lambda * PHIshift.y[4*GT.n + 4 - i] ) ) * imag( GT.y[i] * conj( GT.y[4*GT.n + 4 - i] ) );

        // first part is stored into Psi
        Psi = (Lambda*Lambda/Pi)*integrateReP(Integrand);

        // second integral of Eq. 19 is performed
        Integrand.y[0] = f_BE.y[0] * imag( conj(PHI.y[0]) / (1.0 + Lambda * conj(PHI.y[0]) ) ) * real( conj(GTshift.y[0]) * GTshift.y[0]  );

        for( int i = 1;         i < GT.n+1;      i++ )
        Integrand.y[i] = f_BE.y[i] * imag( conj(PHI.y[4*GT.n + 4 - i]) / (1.0 + Lambda * conj(PHI.y[4*GT.n + 4 - i]) ) ) * real( GTshift.y[i] * conj( GTshift.y[4*GT.n + 4 - i] ) );

        for( int i = 3*GT.n+4; i < 4*GT.n + 4;  i++ )
        Integrand.y[i] = f_BE.y[i] * imag( conj(PHI.y[4*GT.n + 4 - i]) / (1.0 + Lambda * conj(PHI.y[4*GT.n + 4 - i]) ) ) * real( GTshift.y[i] * conj( GTshift.y[4*GT.n + 4 - i] ) );
    }

    if(physics == "EFF3" )
    {
        pair<double, double> Lbrent = brent_find_minima( lambda_implicit_non_zero_T_simple(U, GT, PHI, f_FD, f_BE, K3), Lmin, Lmax, bits);
        Lambda = Lbrent.first;

        // first integral of Eq. 19 is performed
        Integrand.y[0] = f_FD.y[0] * real( PHI.y[0] / (1.0 + Lambda * PHI.y[0] ) ) * imag( GT.y[0] * conj( GT.y[0] ) );

        for( int i = 1; i < GT.n+1; i++ )
        Integrand.y[i] = f_FD.y[i] * real( PHI.y[4*GT.n + 4 - i] / (1.0 + Lambda * PHI.y[4*GT.n + 4 - i] ) ) * imag( GT.y[i] * conj( GT.y[4*GT.n + 4 - i] ) );

        for( int i = 3*GT.n+4; i < 4*GT.n + 4;  i++ )
        Integrand.y[i] = f_FD.y[i] * real( PHI.y[4*GT.n + 4 - i] / (1.0 + Lambda * PHI.y[4*GT.n + 4 - i] ) ) * imag( GT.y[i] * conj( GT.y[4*GT.n + 4 - i] ) );

        // first part is stored into Psi
        Psi = (Lambda*Lambda/Pi)*integrateReP(Integrand);

        // second integral of Eq. 19 is performed
        Integrand.y[0] = f_BE.y[0] * imag( conj(PHI.y[0]) / (1.0 + Lambda * conj(PHI.y[0]) ) ) * real( conj(GT.y[0]) * GT.y[0]  );

        for( int i = 1;         i < GT.n+1;      i++ )
        Integrand.y[i] = f_BE.y[i] * imag( conj(PHI.y[4*GT.n + 4 - i]) / (1.0 + Lambda * conj(PHI.y[4*GT.n + 4 - i]) ) ) * real( GT.y[i] * conj( GT.y[4*GT.n + 4 - i] ) );

        for( int i = 3*GT.n+4; i < 4*GT.n + 4;  i++ )
        Integrand.y[i] = f_BE.y[i] * imag( conj(PHI.y[4*GT.n + 4 - i]) / (1.0 + Lambda * conj(PHI.y[4*GT.n + 4 - i]) ) ) * real( GT.y[i] * conj( GT.y[4*GT.n + 4 - i] ) );

    }

    // Psi is calculated from Eq. 19
    Psi = Psi - (Lambda*Lambda/Pi)*integrateReP(Integrand);
}



struct nT_implicit_non_zero_T
{
    nT_implicit_non_zero_T(double y1, double y2, double y3, AWT & z3, AWT & z4, AWT & z5) : delta(y1), x(y2), Lambda(y3), f_FD(z3), f_BE(z4), K3(z5) { }

    double operator()(double nT)
    {
        // first DOS is loaded with thermal propagator
        AWT DOS;
        DOS.initializeAWT(f_FD.n, f_FD.xMax, f_FD.kT);
        freePropagatorLorentz(delta, x  + 0.5*Lambda - 0.5*Lambda*nT, DOS);

        // then it is multiplied with -2.0*FD/Pi
        for(int i=0; i< f_FD.n + 1; i++)
        DOS.y[i] = f_FD.y[i] * imag(DOS.y[i]);

        for(int i = 3*f_FD.n+4; i < 4*f_FD.n + 4;  i++ )
        DOS.y[i] = f_FD.y[i] * imag(DOS.y[i]);

        return  abs( - nT - 2.0 * integrateReP(DOS) / Pi - 2.0 * DOS.y[3*f_FD.n+4] / ( DOS.xMax * Pi ) );
    }

    double delta;
    double x;
    double Lambda;
    AWT & f_FD;
    AWT & f_BE;
    AWT & K3;
};

pair<double, double> Brent_nT_nonZeroT(double delta, double x, double Lambda, AWT & FD, AWT & BE, AWT & K3)
{
    using boost::math::tools::brent_find_minima;
    int bits = std::numeric_limits<float>::digits;
    //cout << "numeric precision set as " << bits << endl;

    double nMin = 0.0;
    double nMax = 2.0;

    return brent_find_minima( nT_implicit_non_zero_T(delta, x, Lambda, FD, BE, K3), nMin, nMax, bits);
}

struct nS_implicit_non_zero_T
{
    nS_implicit_non_zero_T(double y1, double y2, double y3, AWT & z1, AWT & z2, AWT & z3, AWT & z4, AWT & z5) : delta(y1), x(y2), U(y3), SigmaSpec(z1), Gspec(z2), f_FD(z3), f_BE(z4), K3(z5) { }

    double operator()(double nS)
    {
        // the static bubble psi

        // first DOS is loaded with spectral propagator
        AWT DOS;
        DOS.initializeAWT(Gspec.n, Gspec.xMax, Gspec.kT);
        propagatorLorentz(delta, x-0.5*U*nS+0.5*U , DOS, SigmaSpec);

        // then it is multiplied with -2.0*FD/Pi
        for(int i=0; i< f_FD.n + 1; i++)
        DOS.y[i] = f_FD.y[i] * imag(DOS.y[i]);

        for(int i = 3*f_FD.n+4; i < 4*f_FD.n + 4;  i++ )
        DOS.y[i] = f_FD.y[i] * imag(DOS.y[i]);

        return  abs( - nS - 2.0 *  integrateReP(DOS) / Pi - 2.0 *  DOS.y[3*f_FD.n+4]/ ( DOS.xMax * Pi ) );
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




pair<double, double> Brent_nS_nonZeroT(double delta, double x, double U, AWT & SigmaSpec, AWT & Gspec, AWT & FD, AWT & BE, AWT & K3)
{
    using boost::math::tools::brent_find_minima;
    int bits = std::numeric_limits<float>::digits;
    //cout << "numeric precision set as " << bits << endl;

    double nMin = 0.0;
    double nMax = 2.0;

    return brent_find_minima( nS_implicit_non_zero_T(delta, x, U, SigmaSpec, Gspec, FD, BE, K3), nMin, nMax, bits);
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////             NON ZERO T PARQUETS EFFECTVE INTERACTION APPROX                      //////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////

void Effective_nonZeroT(string input)
{
    // intialization of variables which are then imputed by a txt file
    double delta, mu;
    string model;

    // intialization of iteration parameters
    double kT_min, kT_increment, kT_max, U_min, U_increment, U_max, x_increment, x_range;

    // intialization of variables for mesh properties
    double xMax;
    int n;

    // intialization of variables for printing out txt files
    int display;
    double range;
    int print_mode;
    string export_mode, physics;

    // previously initialized variables are set to values from the input file
    importInitial(input, delta, mu, model,
                  kT_min, kT_increment, kT_max,
                  U_min, U_increment, U_max,
                  x_increment, x_range,
                  xMax, n, display, range,
                  print_mode, export_mode, physics);


    // the number of decimals in kT_min is determined
    int kT_precision = 0;
    int U_precision = 1;
    int x_precision = 0;
    //precision(kT_min, U_min, x_increment, kT_precision, U_precision, x_precision);


    // variables for maximal number of terations are initialized and set
    int kT_iterations, U_iterations, x_iterations;
    iterations(kT_max, kT_min, kT_increment,
               U_max, U_min, U_increment,
               x_range, x_increment,
               kT_iterations, U_iterations, x_iterations);


    //////////////////////////////////////////////////////////////////////////////////////////////////////
    /////////             HARTREE-FOCK AS REFERENCE SYSTEM FOR ENERGY CALCULATIONS              /////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////

    string name = "Txt/Hartree";
    createDirectories(name);

        // creating arrays with Fermi-Dirac, Bose-Einstein statistics and array with hartrre propagator G0
    AWT FD_0, BE_0, G0;
    FD_0.initializeAWT(n, xMax, 0);
    BE_0.initializeAWT(n, xMax, 0);
    G0.initializeAWT(n, xMax, 0);

    FD_0.set_FD();
    string filename = "Txt/Hartree/FD_0";
    FD_0.exportAWTasFUN(filename, display, range, 1, export_mode);

    BE_0.set_BE();
    filename = "Txt/Hartree/BE_0";
    BE_0.exportAWTasFUN(filename, display, range, 1, export_mode);

    freePropagatorLorentz(delta, 0, G0);

    // calculating the Hartree reference energy
    double Ehf;
    hartree_energies(Ehf, n, xMax, delta, display, range, export_mode, x_iterations, x_increment, FD_0, BE_0);


    /////////////////////////////////////////////////////////////////////////////////////////////////////
    /////////                  CALCULATIONS AT DIFFERENT kT, U, x                              /////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////

    // static bubble Psi
    double      Psi_iter[kT_iterations][U_iterations][2*x_iterations-1];
    double   Lambda_iter[kT_iterations][U_iterations][2*x_iterations-1];
    double       nT_iter[kT_iterations][U_iterations][2*x_iterations-1];
    double       nS_iter[kT_iterations][U_iterations][2*x_iterations-1];
    double     Ekin_iter[kT_iterations][U_iterations][2*x_iterations-1];
    double     Ecor_iter[kT_iterations][U_iterations][2*x_iterations-1];
    double    EcorA_iter[kT_iterations][U_iterations][2*x_iterations-1];
    double     Ealt_iter[kT_iterations][U_iterations][2*x_iterations-1];

    double factor_a_iter[kT_iterations][U_iterations][2*x_iterations-1];
    double   hmhw_a_iter[kT_iterations][U_iterations][2*x_iterations-1];
    double   dosf_a_iter[kT_iterations][U_iterations][2*x_iterations-1];
    double        Z_iter[kT_iterations][U_iterations][2*x_iterations-1];
    double         Z_err[kT_iterations][U_iterations][2*x_iterations-1];
    double      chi_iter[kT_iterations][U_iterations][2*x_iterations-1];
    double    gamma_iter[kT_iterations][U_iterations][2*x_iterations-1];
    double    gsgs_iter[kT_iterations][U_iterations][2*x_iterations-1];
    double    gtgt_iter[kT_iterations][U_iterations][2*x_iterations-1];

    // these four AWTs are used as help and storage place, are once allocated in the stack
    AWT aux1;
    aux1.initializeAWT(n, xMax, 0);
    AWT aux2;
    aux2.initializeAWT(n, xMax, 0);
    AWT aux3;
    aux3.initializeAWT(n, xMax, 0);
    AWT aux4;
    aux4.initializeAWT(n, xMax, 0);
    AWT aux5;
    aux5.initializeAWT(n, xMax, 0);
    AWT aux6;
    aux6.initializeAWT(n, xMax, 0);


    // iterations at different temperatures
    for(int k=0; k<1; k++)
    {
        double kT = kT_min + k*kT_increment;
        double shift = Pi*kT;

        // upper kT directory is created for storage
        stringstream kT_stream;
        kT_stream << fixed << setprecision(kT_precision) << kT;
        string kT_string = kT_stream.str();
        string name = "Txt/kT_" +  kT_string;
        createDirectories(name);

        // creating array with Kernel3
        AWT K3, FD, BE;
        K3.initializeAWT(n, xMax, kT);
        FD.initializeAWT(n, xMax, kT);
        BE.initializeAWT(n, xMax, kT);

        K3.set_K3();
        FD.set_FD();
        BE.set_BE();

        if(print_mode == 1)
        {
            name = "Txt/kT_" +  kT_string + "/FD";
            FD.exportAWTasFUN(name, display, range, 1, export_mode);

            name = "Txt/kT_" +  kT_string + "/BE";
            BE.exportAWTasFUN(name, display, range, 1, export_mode);

            cout << "K3, FD, BE, zero, done" << endl;
        }

        // thermodynamic and spectral self-energies
        AWT SigmaTherm, SigmaSpec;

        SigmaTherm.initializeAWT(n, xMax, kT);
        SigmaSpec.initializeAWT(n, xMax, kT);

        SigmaTherm.set_zero();
        SigmaSpec.set_zero();

        //
        AWT GT, GTshift, GS, PHI,PHIshift;

        GT.initializeAWT(n, xMax, kT);
        GTshift.initializeAWT(n, xMax, kT);
        GS.initializeAWT(n, xMax, kT);
        PHI.initializeAWT(n, xMax, kT);
        PHIshift.initializeAWT(n, xMax, kT);

        GT.set_zero();
        GTshift.set_zero();
        GS.set_zero();
        PHI.set_zero();
        PHIshift.set_zero();



        // static bubble Psi, Lambda, Lambda at half-filling, thermodynamic and spectral occupation numbers, energies
        double Psi, Lambda, LambdaHalf, nT, nS, Ekin, Ecor, Ecor_alt;


        if(print_mode == 1)  cout << "Sigma-s, G-s, Phi-s, Lambda, Psi, nT, nS initialized" << endl;


        // cycle for different values of U
        for(int i=0; i<U_iterations; i++)
        {

            // setting the value of U in given cycle
            double U = U_min + i*U_increment;

            // subdirectories for different U calculations are created
            stringstream U_stream;
            U_stream  << fixed << setprecision(U_precision) << U;
            string U_string = U_stream.str();
            string name = "Txt/kT_" +  kT_string + "/U" + U_string;
            createDirectories(name);


            // OUTPUT OF TIME AND CALCULATION PARAMETERS
            time_t now = time(0);
            tm* loc  = localtime(&now);
            if(k==0)
            {
                cout << endl;
                cout << asctime(loc);
                cout << "-----------------------------------------" << endl;
                cout << "---------------- U = " << U << " ----------------" << endl;
                cout << "-----------------------------------------" << endl;
                cout << "model is: " << model << " in " << physics << endl;
                cout << "reference energy is:  " << Ehf << endl;
                cout << "kT iters: " << kT_iterations << "  U iters: " << U_iterations << "  x iters: " << x_iterations << endl;
                cout << "-----------------------------------------" << endl;
            }
            //cout << "--------------------------------------------------------------------------------" << endl;
            //cout << "CALCULATIONS FOR U = " << U << " and kT = " << kT <<  endl;
            //cout << "--------------------------------------------------------------------------------" << endl;


            // setting the value of x in given cycle
            for(int j = 0; j < 2*x_iterations -1; j++)
            {
                double x;
                preset(j, x, x_iterations, x_increment, Lambda, LambdaHalf, nT);


                // SigmaTherm iteration
                if(physics == "RPA")
                {
                    nT = 1;
                    Lambda = 0;
                    LambdaHalf = 0;
                }

                if(physics == "EFF" || physics == "EFF1" || physics == "EFF2" || physics == "EFF3" )
                {

                    //  nT_, Lambda_ store values of nT, Lambda from one preceeding iterations,
                    double      nT_ = 1e3;
                    double  Lambda_ = 1e3;
                    double   nTprec = 0.00001;
                    double lambPrec = 0.00001;

                    // thermal nT and Lambda are determined
                    while( (abs( nT - nT_) > nTprec) && (abs( Lambda - Lambda_) > lambPrec) )
                    {
                            Lambda_ = Lambda;
                            nT_ = nT;
                            // Gtherm in given iter calculated from previous SigmaTherm
                            freePropagatorLorentz     (delta, x - 0.5 * Lambda * ( nT - 1 ), GT);

                            if(physics == "EFF" || physics == "EFF1" || physics == "EFF2" )
                            freePropagatorLorentzShift(delta, x - 0.5 * Lambda * ( nT - 1 ), shift, GTshift);

                            // Phi in given iter calculated from previous SigmaTherm
                            PhiFunction(GT, PHI, FD, K3, aux1, aux2, aux3, aux4, aux5, aux6);

                            if(physics == "EFF" || physics == "EFF1" || physics == "EFF2" )
                            PhiFunctionShift(GT, GTshift, PHIshift, FD, K3, aux1, aux2, aux3, aux4, aux5, aux6);

                            // setting minimum Lambda value for given x
                            double Lmin;
                            if( j == 0 )    Lmin = 0;
                            else
                            {
                                if( j == x_iterations)     Lmin = LambdaHalf;
                                else                       Lmin = Lambda_iter[k][i][j-1];
                            }

                            // Lambda and Psi are calculated and stored in the variables
                            if(physics == "EFF")         Lambda_nonZeroT(U, Lmin, physics, GT, GTshift, PHI, PHIshift, Lambda, Psi, FD, BE, K3);
                            if(physics == "EFF1")        Lambda_nonZeroT(U, Lmin, physics, GT, GTshift, PHI, PHIshift, Lambda, Psi, FD_0, BE_0, K3);
                            if(physics == "EFF2")        Lambda_nonZeroT(U, Lmin, physics, GT, GTshift, PHI, PHIshift, Lambda, Psi, FD_0, BE_0, K3);
                            if(physics == "EFF3")        Lambda_nonZeroT(U, Lmin, physics, GT, GTshift, PHI, PHIshift, Lambda, Psi, FD_0, BE_0, K3);

                            // half-filling value of Lambda is stored for later use
                            if( j == 0 )  LambdaHalf = Lambda;

                            // nT is given by an implicit equation 24 and solved with Brent
                            pair<double, double> nTbrent = Brent_nT_nonZeroT(delta, x, Lambda, FD, BE, K3);
                            nT = nTbrent.first;
                    }
                }

                if(x == 0)  nT = 1;

                // thermal G and shifted thermal G are recalculated and exported
                freePropagatorLorentz(delta, x - 0.5 * U * ( nT - 1 ), GT);
                name = "Txt/kT_" +  kT_string + "/U" + U_string +  "/Gtherm_" + boost::lexical_cast<std::string>(x);
                GT.exportAWTasFUN(name, display, range, 1, export_mode);

                if(physics == "EFF" || physics == "EFF1" || physics == "EFF2")
                {
                    freePropagatorLorentzShift(delta, x - 0.5 * Lambda * ( nT - 1 ), shift, GTshift);
                    name = "Txt/kT_" +  kT_string + "/U" + U_string + "/GTs_" + boost::lexical_cast<std::string>(x);
                    //GTshift.exportAWTasFUN(name, display, range, 1, export_mode);
                }


                // Phi bubble and Phi bubble shifted are recalculated and exported
                PhiFunction(GT, PHI, FD, K3, aux1, aux2, aux3, aux4, aux5, aux6);
                name = "Txt/kT_" +  kT_string + "/U" + U_string + "/PHI_" + boost::lexical_cast<std::string>(x);
                PHI.exportAWTasFUN(name, display, range, 1, export_mode);

                if(physics == "EFF" || physics == "EFF1" || physics == "EFF2")
                {
                    PhiFunctionShift(GT, GTshift, PHIshift, FD, K3, aux1, aux2, aux3, aux4, aux5, aux6);
                    name = "Txt/kT_" +  kT_string + "/U" + U_string + "/PHIs_" + boost::lexical_cast<std::string>(x);
                    //PHIshift.exportAWTasFUN(name, display, range, 1, export_mode);
                }

                // help function
                AWT Int;
                Int.initializeAWT(FD.n, FD.xMax, FD.kT);

                // spectral self-energy is calculated and stored
                if(physics == "RPA" )    spectral_sigma_RPA(U, PHI, GT, FD, BE, K3, SigmaSpec, Int, aux1, aux2, aux3, aux4, aux5, aux6);
                if(physics == "EFF" )    spectral_sigma_EFF(U, Lambda, PHI, GT, FD, BE, K3, SigmaSpec, Int, aux1, aux2, aux3, aux4, aux5, aux6);
                if(physics == "EFF1")    spectral_sigma_EFF(U, Lambda, PHI, GT, FD_0, BE_0, K3, SigmaSpec, Int, aux1, aux2, aux3, aux4, aux5, aux6);
                if(physics == "EFF2")    spectral_sigma_EFF(U, Lambda, PHI, GT, FD_0, BE_0, K3, SigmaSpec, Int, aux1, aux2, aux3, aux4, aux5, aux6);
                if(physics == "EFF3")    spectral_sigma_EFF(U, Lambda, PHI, GT, FD_0, BE_0, K3, SigmaSpec, Int, aux1, aux2, aux3, aux4, aux5, aux6);

                name = "Txt/kT_" +  kT_string + "/U" + U_string + "/SigmaSpec_" + boost::lexical_cast<std::string>(x);
                SigmaSpec.exportAWTasFUN(name, display, range, 1, export_mode);

                // spectral occupancy number is determined by Brent algorithm
                pair<double, double> nSbrent = Brent_nS_nonZeroT(delta, x, U, SigmaSpec, GS, FD, BE, K3);
                nS = nSbrent.first;
                if(x == 0) nS =1;

                // spectral G is recalculated and exported
                propagatorLorentz(delta, x -0.5*U*nS + 0.5*U, GS, SigmaSpec);
                name = "Txt/kT_" +  kT_string + "/U" + U_string +  "/Gspec_" + boost::lexical_cast<std::string>(x);
                //GS.exportAWTasFUN(name, display, range, 1, export_mode);



                // kinetic E is determined
                kinetic_energy(Ekin, Ehf, U, x, nS, delta, SigmaSpec, FD, Int);

                // correlation E is determined
                if(physics == "RPA" || physics == "EFF" || physics == "EFF1" || physics == "EFF3" )
                {
                    correlation_energy(Ecor, U, nT, GT, SigmaSpec, FD, Int);
                    correlation_energy_alt(Ecor_alt, U, Lambda, nT, PHI, BE_0, Int);
                }


                if(physics == "EFF2")
                correlation_energy(Ecor, U, nT, GT, SigmaSpec, FD_0, Int);

                AWT X;
                X.initializeAWT(n, xMax, kT);



                AWT GG;
                GG.initializeAWT(n, xMax, kT);
                for(int ii=0;  ii<4*GG.n + 4;  ii++)     GG.y[ii] = GT.y[ii] * GT.y[ii];
                name = "Txt/kT_" +  kT_string + "/U" + U_string +  "/GG_" + boost::lexical_cast<std::string>(x);
                GG.exportAWTasFUN(name, display, range, 1, export_mode);

                AWT f;
                f.initializeAWT(n, xMax, kT);
                for(int ii=0;  ii<4*X.n + 4;  ii++)     f.y[ii] = -PHI.y[ii]/( 1.0 + Lambda*PHI.y[ii] );

                X.boseMatsubaraImP(f, GG, 1, 1, FD, BE, aux1, aux2, aux3, aux4, aux5, aux6);
                X.KrammersKronig(X,K3, aux1, aux2, aux3, aux4);

                for(int ii=0;  ii<4*X.n + 4;  ii++)     X.y[ii] = X.y[ii];

                name = "Txt/kT_" +  kT_string + "/U" + U_string +  "/X_" + boost::lexical_cast<std::string>(x);
                X.exportAWTasFUN(name, display, range, 1, export_mode);


                // wilson ratio
                double chi;
                for(int ii=0;  ii<4*Int.n + 4;  ii++)     Int.y[ii] =  FD.y[ii]*imag(GS.y[ii] * GS.y[ii] * ( 1.0 - U*X.y[ii]/ ( 1.0 + Lambda*real(PHI.y[0])  ) ) );


                chi = (2.0/Pi)* integrateReP(Int)   ;

                for(int ii=0;  ii<4*Int.n + 4;  ii++)     Int.y[ii] = FD.y[ii]*imag(GT.y[ii] * GT.y[ii]);


                double chi_0;
                if(i == 0) chi_0 = chi;

                double gamma;
                double der;
                der = SigmaSpec.n *    real( SigmaSpec.y[1] - SigmaSpec.y[4*n+3] ) / ( 2.0 * SigmaSpec.xMax );
                gamma = Pi*( 2.0 - der) / 3.0;

                double gamma_0;
                if(i == 0) gamma_0 = gamma;

                // physical parameters are stored for later output
                      nT_iter[k][i][j] = nT;
                      nS_iter[k][i][j] = nS;
                  Lambda_iter[k][i][j] = Lambda;
                     Psi_iter[k][i][j] = Psi;

                    Ekin_iter[k][i][j] = Ekin;
                    Ecor_iter[k][i][j] = Ecor;
                   EcorA_iter[k][i][j] = Ecor_alt;
                    Ealt_iter[k][i][j] = 0;

                    chi_iter[k][i][j] = chi;
                    gamma_iter[k][i][j] = gamma;


                // Kondo scales are calculated and stored for later output
                factor_a_iter[k][i][j] = factor_a(Lambda, PHI.y[0]);
                  dosf_a_iter[k][i][j] = dosf_a(U, GS.y[0]);
                  hmhw_a_iter[k][i][j] = hmhw_a(GS, 2*U);


                 pair<double, double> Z;
                 Z = QuasiParticleWeight(SigmaSpec);
                       Z_iter[k][i][j] = Z.first;
                        Z_err[k][i][j] = Z.second;

                // resultes are written on screen
                cout //<< "x = "      << x      <<  ", nT = " << nT          << ", nS = " << nS
                         //<< "kT = " << kT
                         //<<  ",  E = "    << Ekin + Ecor
                         //<< " int_gtgt " << int_gtgt
                         //<< " chi = " << chi
                         << " g/g_0 = " << gamma/gamma_0
                         //<< " R = " << ( chi/chi_0) / ( gamma/gamma_0)
                         //<<  "  E = " << Ekin + Ecor
                         //<<  "  E_k = "    <<  Ekin
                         //<<  "  Ecor " <<  Ecor
                         //<<  "  E_c = " <<  Ecor_alt
                         <<  "  L =  " << Lambda
                         << " der = " << der
                         <<  endl;



            } // end of x iteration instance


            /////////////////////////////////////////////////////////////////////////////////////////////////////////
            ///////////                     Txt file output for each x                                 /////////////
            ///////////////////////////////////////////////////////////////////////////////////////////////////////

            // txt file for physical parameters
            string out = "Txt/kT_" +  kT_string + "/U" + U_string + ".txt";
            ofstream output(out.c_str());

            // x < 0 values are outputed first
            for(int ii = 2*x_iterations-2; ii > x_iterations -1 ; ii --)
            output << kT << "  " << (x_iterations - 1 -ii)*x_increment
                         << "  " <<    nT_iter[k][i][ii]     << "  " <<    nS_iter[k][i][ii]  << "  " << Lambda_iter[k][i][ii]
                         << "  " << chi_iter[k][i][ii]  << "  " << gamma_iter[k][i][ii]
                         << "  " << gsgs_iter[k][i][ii]  << "  " << gtgt_iter[k][i][ii]
                         << "  " <<   Ekin_iter[k][i][ii] + Ecor_iter[k][i][ii]
                         << "  " <<  Ekin_iter[k][i][ii] +  EcorA_iter[k][i][ii]
                         << "  " <<  Ekin_iter[k][i][ii]     << "  " <<  Ecor_iter[k][i][ii]  << "  " <<  EcorA_iter[k][i][ii]
                         << "  " <<   Ealt_iter[k][i][ii]    << "  " <<   Psi_iter[k][i][ii]
                         << "  " << factor_a_iter[k][i][ii]  << "  " << hmhw_a_iter[k][i][ii] << "  " << dosf_a_iter[k][i][ii]
                         << "  " <<    Z_iter[k][i][ii]      << "  " << Z_err[k][i][ii]
                         << endl;


            // x >=0 values are outputed
            for(int ii = 0; ii < x_iterations; ii++)
            output << kT << "  " << (x_iterations - 1 -ii)*x_increment
                         << "  " <<    nT_iter[k][i][ii]     << "  " <<    nS_iter[k][i][ii]  << "  " << Lambda_iter[k][i][ii]
                         << "  " <<   Ekin_iter[k][i][ii] + Ecor_iter[k][i][ii]
                         << "  " << chi_iter[k][i][ii]  << "  " << gamma_iter[k][i][ii]
                         << "  " << gsgs_iter[k][i][ii]  << "  " << gtgt_iter[k][i][ii]
                         << "  " <<  Ekin_iter[k][i][ii] +  EcorA_iter[k][i][ii]
                         << "  " <<  Ekin_iter[k][i][ii]     << "  " <<  Ecor_iter[k][i][ii]  << "  " <<  EcorA_iter[k][i][ii]
                         << "  " <<   Ealt_iter[k][i][ii]    << "  " <<   Psi_iter[k][i][ii]
                         << "  " << factor_a_iter[k][i][ii]  << "  " << hmhw_a_iter[k][i][ii] << "  " << dosf_a_iter[k][i][ii]
                         << "  " <<    Z_iter[k][i][ii]      << "  " << Z_err[k][i][ii]
                         << endl;


        }  // end of U iteration instance
    }  // end of kT iteration instance


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////                     TXT FILE OUTPUT                                 ////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

        //string out = "final_R";
        //ofstream output(out.c_str());

        string out = "final_g";
        ofstream output(out.c_str());
        // cycle for different values of U
        for(int i=0; i<U_iterations; i++)
        {
            // setting the value of U in given cycle
            double U = U_min + i*U_increment;



            for(int j = 0; j < 2*x_iterations -1; j++)
            {
                // setting the right value of x
                double x;
                if(j > x_iterations - 1)        x = (x_iterations - 1 -j)*x_increment;
                else                            x = j*x_increment;





                for(int k=0; k<kT_iterations; k++)
                {
                    double kT = kT_min + k*kT_increment;

                    output <<  U  << "  " //<<            x
                                   //<< "  " <<    nT_iter[k][i][j]     << "  " <<    nS_iter[k][i][j]   << "  " << Lambda_iter[k][i][j]
                                   //<< "  " << chi_iter[k][i][j]
                                   << "  " << gamma_iter[k][i][j]/gamma_iter[k][0][j]
                                   //<< "  " << ( chi_iter[k][i][j] /chi_iter[k][0][j] ) / ( gamma_iter[k][i][j] / gamma_iter[k][0][j] )
                                   // << "  " << gsgs_iter[k][i][j]  << "  " << gtgt_iter[k][i][j]
                                   //<< "  " <<   Ekin_iter[k][i][j] + Ecor_iter[k][i][j]
                                   //<< "  " <<  Ekin_iter[k][i][j] +  EcorA_iter[k][i][j]
                                   //<< "  " <<  Ekin_iter[k][i][j]     << "  " <<   Ecor_iter[k][i][j]  << "  " <<  EcorA_iter[k][i][j]
                                   //<< "  " <<   Ealt_iter[k][i][j]    << "  " <<    Psi_iter[k][i][j]
                                   //<< "  " << factor_a_iter[k][i][j]  << "  " << hmhw_a_iter[k][i][j]  << "  " << dosf_a_iter[k][i][j]
                                   //<< "  " <<    Z_iter[k][i][j]      << "  " <<       Z_err[k][i][j]
                                   << endl;

                }  // end of kT iteration
            }  // end of x iteration
        }  // end of U iteration




}






