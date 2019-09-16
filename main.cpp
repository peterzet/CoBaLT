#include "src/physics.h"
#include "src/AWT.h"
#include "src/Methods.h"
#include "src/structures.h"
#include "src/test/test.h"


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

//double PII = 3.14;



int main()
{
    // TO DO 1: initialize appropriate structure

    double delta, mu;
    string model;

    double kT_min, kT_increment, kT_max;
    double U_min, U_increment, U_max;
    double x_increment, x_max;

    // intialization of variables for mesh properties
    double xMax;
    int n;
    // intialization of variables for printing out txt files
    int display;
    double range;
    int print_mode;
    string export_mode, physics;

    // TO do 1: end

   string name = "input";

/*
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


    // intialization of variables which are then imputed by a txt file

    // previously initialized variables are set to values from the input file
    importInitial(name, delta, mu, model,
                  kT_min, kT_increment, kT_max,
                  U_min, U_increment, U_max,
                  x_increment, x_max,
                  xMax, n, display, range,
                  print_mode, export_mode, physics);

*/


    double Pi = 3.14159265359;

    //test_convolutions(n, xMax, range, export_mode);
    //dynamic_spins_half(name);
    //Effective_dynamic_2(name);
    //dynamic_spins_matrix(name);
    //static_spins_matrix(name);
    //static_spins_specs(name);

    //static_Lp(name);
    //static_T_zero(name);


    return 0;
}



