#include "src/AWT.h"
#include "src/structures.h"
#include "src/physics.h"

// files for static problems
//#include "src/static/static_T0.h"
#include "src/static/static_spin_polarized_nov2018.h"

// files for dynamic problems
//#include "src/dynamic/dynamic_spinles_experimental.h"
//#include "src/dynamic/dynamic_spinles.h"
//#include "src/dynamic/dynamic_spins.h"
//#include "src/dynamic/dynamic_spins_half.h"


// files for convolution test
#include "src/test/test.h"



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

// boost libraries
#include <boost/math/tools/minima.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/filesystem.hpp>

//#define BOOST_MATH_INSTRUMENT

using namespace std;



int main()
{
    info input;
	string name = "input";
	input.importing(name);


    //test_convolutions(input.n, input.xMax, input.range, input.output_mode);
    //dynamic_spins_half(name);
    //Effective_dynamic_2(name);
    //dynamic_spins_matrix(name);
    static_spins_matrix(name);
    //static_spins_specs(name);

    //static_Lp(name);
    //static_T_zero(name);


    return 0;
}
