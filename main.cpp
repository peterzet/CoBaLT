#include "Constants.h"
#include "awt.h"
#include "aux.h"
#include "Methods.h"
#include "physics.h"
#include "EffectiveT0.h"
#include "Effective_freq.h"

#include <boost/math/tools/minima.hpp>
#include <cassert>          // error handling library, function assert to terminate the program
#include <cmath>            // declares some common mathematical operations and transformation
#include <cstdlib>          // several general purpose functions, including dynamic memory management,
                            // random number generation, communication with the environment,
                            // integer arithmetics, searching, sorting and converting
#include <fftw3.h>          // FFTW library
#include <fstream>          // standard library for showing outputs

#include <iostream>         // standard library for reading inputs

#include <string>           // standard library for manipultaing strings
#include <sstream>
#include <vector>
#include <iomanip> // setprecision

#include <boost/lexical_cast.hpp>
#include <boost/filesystem.hpp>

//#define BOOST_MATH_INSTRUMENT

using namespace std;



int main()
{


    double delta;
    double mu;
    string model;

    double kT_min;
    double kT_increment;
    double kT_max;

    double U_min;
    double U_increment;
    double U_max;

    double x_increment;
    double x_max;

    // intialization of variables for mesh properties
    double xMax;
    int n;
    // intialization of variables for printing out txt files
    int display;
    double range;
    int print_mode;
    string export_mode, physics;






    string name = "input";

    importInitial(name, delta, mu, model,
                  kT_min, kT_increment, kT_max,
                  U_min, U_increment, U_max,
                  x_increment, x_max,
                  xMax, n, display, range,
                  print_mode, export_mode, physics);


    name = "Txt/def";
    createDocument(name, delta, mu, model,
                  kT_min, kT_increment, kT_max,
                  U_min, U_increment, U_max,
                  x_increment, x_max,
                  xMax, n, display, range,
                  print_mode, export_mode, physics);



    name = "input";


    // intialization of variables which are then imputed by a txt file

    // previously initialized variables are set to values from the input file
    importInitial(name, delta, mu, model,
                  kT_min, kT_increment, kT_max,
                  U_min, U_increment, U_max,
                  x_increment, x_max,
                  xMax, n, display, range,
                  print_mode, export_mode, physics);


    //test_convolutions(n, xMax, range, export_mode);
    //Effective_frequency_dep(name);

    string callme;
    AWT test;
    test.initializeAWT(1000,10,0);
    test.set_FD();
    cout << test.nn << endl;
    callme = "FD";
    test.exportAWTasAWT(callme, 1, export_mode);


    return 0;
}



