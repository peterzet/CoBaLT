#include "Effective_freq.h"
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
#include <boost/lexical_cast.hpp>
#include <boost/filesystem.hpp>
#include <iomanip> // setprecision

using namespace std;





/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////             ZERO T FREQUENCY DEPENDENT PARQUETS EFFECTVE INTERACTION APPROX                      //////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////

void Effective_frequency_dep(string input)
{
    double Pi = 3.14159265359;

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

    physics = "EFF_freq";

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



    // the number of decimals in kT_min is determined
    int kT_precision, U_precision, x_precision;
    precision(kT_min, U_min, x_increment, kT_precision, U_precision, x_precision);


    // variables for maximal number of terations are initialized and set
    int kT_iterations = 0;
    int U_iterations, x_iterations;
    iterations(kT_max, kT_min, kT_increment,
               U_max, U_min, U_increment,
               x_range, x_increment,
               kT_iterations, U_iterations, x_iterations);


    /////////////////////////////////////////////////////////////////////////////////////////////////////
    /////////                  CALCULATIONS AT DIFFERENT U, x                              /////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////

    double kT = 0;
    // static bubble Psi
    double      Psi_iter[U_iterations][2*x_iterations-1];
    double   Lambda_iter[U_iterations][2*x_iterations-1];
    double       nT_iter[U_iterations][2*x_iterations-1];
    double       nS_iter[U_iterations][2*x_iterations-1];
    double     Ekin_iter[U_iterations][2*x_iterations-1];
    double     Ecor_iter[U_iterations][2*x_iterations-1];
    double    EcorA_iter[U_iterations][2*x_iterations-1];
    double     Ealt_iter[U_iterations][2*x_iterations-1];

    double factor_a_iter[U_iterations][2*x_iterations-1];
    double   hmhw_a_iter[U_iterations][2*x_iterations-1];
    double   dosf_a_iter[U_iterations][2*x_iterations-1];
    double        Z_iter[U_iterations][2*x_iterations-1];
    double         Z_err[U_iterations][2*x_iterations-1];


    // creating arrays
    AWT K3, FD, BE;
    K3.initializeAWT(n, xMax, kT);
    FD.initializeAWT(n, xMax, kT);
    BE.initializeAWT(n, xMax, kT);

    K3.set_K3();
    FD.set_FD();
    BE.set_BE();


    // thermodynamic and spectral self-energies, loaded with zeros
    AWT SigmaTherm, SigmaThermOld, SigmaSpec;

    SigmaTherm.initializeAWT(n, xMax, kT);
    SigmaThermOld.initializeAWT(n, xMax, kT);
    SigmaSpec.initializeAWT(n, xMax, kT);

    SigmaTherm.set_zero();
    SigmaThermOld.set_zero();
    SigmaSpec.set_zero();

    // different frequency dependent arrays
    AWT GT, GS;
    GT.initializeAWT(n, xMax, kT);
    GS.initializeAWT(n, xMax, kT);

    // specific functions
    AWT LPp, LLPp, KPm, Lambda, LambdaOld, K;

    LPp.initializeAWT(n, xMax, kT);
    LLPp.initializeAWT(n, xMax, kT);
    KPm.initializeAWT(n, xMax, kT);
    Lambda.initializeAWT(n, xMax, kT);
    LambdaOld.initializeAWT(n, xMax, kT);
    K.initializeAWT(n, xMax, kT);


    LPp.set_zero();
    LLPp.set_zero();
    KPm.set_zero();
    Lambda.set_zero();
    K.set_zero();

    // set thermodynamic and spectral occupation numbers, energies
    double nT, nS;
    nT =1;
    nS = 1;

    if(print_mode == 1)  cout << "Sigma-s, G-s, Lambda, K initialized" << endl;


    // cycle for different values of U
    for(int i=0; i<U_iterations; i++)
    {
        // setting the value of U in given cycle
        double U = U_min + i*U_increment;

        // subdirectories for different U calculations are created
        stringstream U_stream;
        U_stream  << fixed << setprecision(U_precision) << U;
        string U_string = U_stream.str();
        string name = "Txt/U" + U_string;
        createDirectories(name);


        // OUTPUT OF starting TIME AND CALCULATION PARAMETERS
        cout << __TIME__ << ",   " << __DATE__ << endl;
        cout << "---------------------------------------" << endl;
        cout << "CALCULATIONS FOR U = " << U <<  endl;
        cout << "frequency dependent vertex   " << endl;
        cout << "---------------------------------------" << endl;


        // setting the value of x in given cycle
        for(int j = 0; j < 2*x_iterations -1; j++)
        {

            double mu = U/2.0; //- ( x_iterations-1 ) * x_increment + j*x_increment;


            // Sigma loop strats here

            // set the SigmaTherm for first iteration
            SigmaTherm.set_real(U/2.0);
            name = "Txt/U" + U_string + "/Stherm_init";
            SigmaTherm.exportAWTasFUN(name, display, range, 1, export_mode);


            // set Sigma iteration precision
	        double SigmaPrecision = 0.5;
            int iter_sigma = 0;


            //while( real_norm(SigmaThermOld, SigmaTherm) > SigmaPrecision)
            for(iter_sigma = 0; iter_sigma <1; iter_sigma = iter_sigma +1)
            {
                //iter_sigma = iter_sigma + 1;
                cout << " entering Sigma iteration:  " << iter_sigma << endl;

                // SigmaThermOld = SigmaTherm;
                for(int i=0;                i<SigmaTherm.n+1;   i++ )  SigmaThermOld.y[i] = SigmaTherm.y[i];
                for(int i=3*SigmaTherm.n+4; i<4*SigmaTherm.n+4; i++ )  SigmaThermOld.y[i] = SigmaTherm.y[i];

                // propagators are recalculated with previously iterated Sigmas
                propagatorLorentz(delta, mu, GT, SigmaTherm);
                name = "Txt/U" + U_string + "/GT_" + boost::lexical_cast<std::string>(iter_sigma);
                GT.exportAWTasFUN(name, display, range, 1, export_mode);


                // set values for the first lambda loop iteration
                for(int i=0;                i<Lambda.n+1;     i++ )  Lambda.y[i] = U;
                for(int i=Lambda.n+1;       i<3*Lambda.n+4;   i++ )  Lambda.y[i] = 0;
                for(int i=3*Lambda.n+4;     i<4*Lambda.n+4;   i++ )  Lambda.y[i] = U;

                name = "Txt/U" + U_string + "/Lambda_init";
                Lambda.exportAWTasFUN(name, display, range, 1, export_mode);


                // set Lambda iteration precision
                double LambdaPrecision = 0.001;


                // Lambda iterations
                // while( real_norm(LambdaOld, Lambda) > LambdaPrecision)
                for(int iter_lambda = 0; iter_lambda < 1; iter_lambda++)
                {
                    // Lambda from previous iteration is stored: LambdaOld = Lambda;
                    for(int i=0;                i<Lambda.n+1;     i++ )  LambdaOld.y[i] = Lambda.y[i];
                    for(int i=Lambda.n+1;       i<3*Lambda.n+4;   i++ )  LambdaOld.y[i] = 0;
                    for(int i=3*Lambda.n+4;     i<4*Lambda.n+4;   i++ )  LambdaOld.y[i] = Lambda.y[i];

                    // LPp and LLPp functions are calculated from Lambda of the previous iteration
                    LPp_function(LPp, Lambda, GT, FD, BE, K3, aux1, aux2, aux3, aux4, aux5, aux6);

                    // export LPp
                    name = "Txt/U" + U_string + "/LPp_"
                            + boost::lexical_cast<std::string>(iter_sigma) + "_"
                            + boost::lexical_cast<std::string>(iter_lambda);
                    LPp.exportAWTasFUN(name, display, range, 1, export_mode);


                    // K is calculated from (17b);
                    for(int i=0;     i<K.n+1;   i++ )       K.y[i] = -Lambda.y[i] * (1.0 - 1.0 / ( 1.0 + LPp.y[i] ) );
                    for(int i=n+1;   i<3*K.n+3; i++ )       K.y[i] = 0;
                    for(int i=3*n+4; i<4*K.n+4; i++ )       K.y[i] = -Lambda.y[i] * (1.0 - 1.0 / ( 1.0 + LPp.y[i] ) );

                    // calculate KPm
                    KPm_function(KPm, K, GT, FD, BE, K3, aux1, aux2, aux3, aux4, aux5, aux6);

                    // export KPm
                    name = "Txt/U" + U_string + "/KPm_"
                            + boost::lexical_cast<std::string>(iter_sigma) + "_"
                            + boost::lexical_cast<std::string>(iter_lambda);
                    KPm.exportAWTasFUN(name, display, range, 1, export_mode);

                    // export K
                    name = "Txt/U" + U_string + "/K_"
                            + boost::lexical_cast<std::string>(iter_sigma) + "_"
                            + boost::lexical_cast<std::string>(iter_lambda);
                    K.exportAWTasFUN(name, display, range, 1, export_mode);


                    // Lambda is calculated from (17a);
                    for(int i=0;              i<Lambda.n+1;      i++ )  Lambda.y[i] = U / ( 1.0 + KPm.y[i] ) ;
                    // for(int i=n+1;   i<3*Lambda.n+3; i++ )  Lambda.y[i] = 0;
                    for(int i=3*Lambda.n+4;   i<4*Lambda.n+4;    i++ )  Lambda.y[i] = U / ( 1.0 + KPm.y[i] ) ;


                    // export Lambda
                    name = "Txt/U" + U_string + "/Lambda_"
                            + boost::lexical_cast<std::string>(iter_sigma) + "_"
                            + boost::lexical_cast<std::string>(iter_lambda);
                    Lambda.exportAWTasFUN(name, display, range, 1, export_mode);


                    cout << " finishing Lambda iteration:  " << real_norm(LambdaOld, Lambda) << endl;
                } // end of Lambda loop



                // SPECTRAL PROPERTIES OF THE SYSTEM
                SigmaSpec_frequency(SigmaSpec, U, LPp, GT, FD, BE, K3, aux1, aux2, aux3, aux4, aux5, aux6);
                name = "Txt/U" + U_string + "/Sspec_"
                        + boost::lexical_cast<std::string>(iter_sigma);
                SigmaSpec.exportAWTasFUN(name, display, range, 1, export_mode);

                propagatorLorentz(delta, mu - U/2, GS, SigmaSpec);
                name = "Txt/U" + U_string + "/GS_"
                        + boost::lexical_cast<std::string>(iter_sigma);
                GS.exportAWTasFUN(name, display, range, 1, export_mode);

                // THERMODYNAMIC SELF-ENERGY
                SigmaTherm_frequency(SigmaTherm, Lambda, U, nT, GT, FD, BE, K3, aux1, aux2, aux3, aux4, aux5, aux6);
                name = "Txt/U" + U_string + "/Stherm_"
                        + boost::lexical_cast<std::string>(iter_sigma);
                SigmaTherm.exportAWTasFUN(name, display, range, 1, export_mode);

                cout << " finishing Sigma iteration:  " << real_norm(SigmaThermOld, SigmaTherm) << endl;


                propagatorLorentz(delta, mu, GT, SigmaTherm);
                name = "Txt/U" + U_string + "/GT_" + boost::lexical_cast<std::string>(mu);
                GT.exportAWTasFUN(name, display, range, 1, export_mode);


            } // end of sigma loop

            // calculate nT
            nT = 1;




		} // end of x iteration instance





    }  // end of U iteration instance

} // end of the whole function
