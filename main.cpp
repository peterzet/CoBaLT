#include "src/AWT.h"
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



int main()
{
	info input;
	string name = "input";
	input.importing(name);

	test_convolutions(input.n, input.xMax, input.range, input.output_mode);

	// various physical models are selcted here




    return 0;
}



