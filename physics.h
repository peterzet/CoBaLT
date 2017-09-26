#ifndef PHYSICS_H
#define PHYSICS_H


#include "Calculator.h"

#include <iostream>         // standard library for reading inputs
#include <fstream>          // standard library for showing outputs
#include <string>           // standard library for manipultaing strings
#include <cassert>          // error handling library, function assert to terminate the program

#include <vector>
#include <complex> // implements the complex class to contain complex numbers in cartesian form
#include <fftw3.h>          // FFTW library
#include <cmath>            // declares some common mathematical operations and transformation
#include <cstdlib>          // several general purpose functions, including dynamic memory management,
                            // random number generation, communication with the environment,
                            // integer arithmetics, searching, sorting and converting

void createDirectories(string & );

// import/export of model parameters
void createDocument(string &, double &, double &, string &,
                    double &, double &, double &,
                    double &, double &, double &,
                    double &, double &,
                    double &, int &,
                    int &, double &, int &, string &, string &);

void importInitial(string &, double &, double &, string &,
                    double &, double &, double &,
                    double &, double &, double &,
                    double &, double &,
                    double &, int &,
                    int &, double &, int &, string &, string &);

void precision(double, double, double, int &, int &, int &);
void iterations(double, double, double, double, double, double, double, double, int &, int &, int &);
void hartree_energies(double &, int, double, double, int, double, string, int, double, AWT &, AWT &);
void preset(int, double &, int, double, double &, double &, double &);

void spectral_sigma_RPA(double, AWT &, AWT &, AWT &, AWT &, AWT &, AWT &, AWT &, AWT &, AWT &, AWT &, AWT &, AWT &, AWT &);
void spectral_sigma_EFF(double, double, AWT &, AWT &, AWT &, AWT &, AWT &, AWT &, AWT &, AWT &, AWT &, AWT &, AWT &, AWT &, AWT &);

void kinetic_energy(double &, double, double, double, double, double, AWT &, AWT &, AWT &);
void correlation_energy(double &, double, double, AWT &, AWT &, AWT &, AWT &);
void correlation_energy_alt(double &, double, double, double, AWT &, AWT &, AWT &);

// Peak analysis
int         peak_number(AWT &, double range);
int     left_satelite_n(AWT &, double range);
int          left_low_n(AWT &, double range);
double  left_satelite_x(AWT &, double range);
double       left_low_x(AWT &, double range);

int    right_satelite_n(AWT &, double range);
int         right_low_n(AWT &, double range);
double right_satelite_x(AWT &, double range);
double      right_low_x(AWT &, double range);

int    Kondo_peak_n(AWT &, int n_left, int n_right);
double Kondo_peak_x(AWT &, int n_left, int n_right);

// Kondo scale analysis
double   hmhw_a(AWT &, double);
double factor_a(double, complex<double>);
double  bethe_a(double, double, double);
double   dosf_a(double, complex<double>);

pair <double,double> QuasiParticleWeight(AWT &);

// physical quantities
void PhiFunction(AWT &, AWT &, AWT &, AWT &, AWT &, AWT &, AWT &, AWT &, AWT &, AWT &);
void PhiFunctionShift(AWT &, AWT &, AWT &, AWT &, AWT &, AWT &, AWT &, AWT &, AWT &, AWT &, AWT &);

void PsiFunction(AWT &, AWT &, AWT &, double &);
void KvertexFunction(double, AWT &, AWT &);
void Xfunction(double, AWT &, AWT &, AWT &, AWT &, AWT &, AWT &, AWT &, AWT &, AWT &, AWT &, AWT &, AWT &);

void LPp_function(AWT &, AWT &, AWT &, AWT &, AWT &, AWT &, AWT &, AWT &, AWT &, AWT &, AWT &, AWT &);
void LLPp_function(AWT &, AWT &, AWT &, AWT &, AWT &, AWT &, AWT &, AWT &, AWT &, AWT &, AWT &, AWT &);
void KPm_function(AWT &, AWT &, AWT &, AWT &, AWT &, AWT &, AWT &, AWT &, AWT &, AWT &, AWT &, AWT &);

void SigmaTherm_frequency(AWT &, AWT &, double, double, AWT &, AWT &, AWT &, AWT &, AWT &, AWT &, AWT &, AWT &, AWT &, AWT &);
void SigmaSpec_frequency(AWT &, double, AWT &, AWT &, AWT &, AWT &, AWT &, AWT &, AWT &, AWT &, AWT &, AWT &, AWT &);


// LORENTZ propagator and DOS
void propagatorLorentz(double, double, AWT &, AWT &);
void freePropagatorLorentz(double, double, AWT &);
void propagatorLorentzShift(double, double, double, AWT &, AWT &);
void freePropagatorLorentzShift(double, double, double, AWT &);

void    DensityLorentz(double, double, AWT &, AWT &);

// semi eliptic propagator and DOS
complex<double> propagatorSemiEl(complex<double>);
void    densitySemiEl(double, double, AWT &, AWT &);

double real_norm(AWT &, AWT &);

// TESTING OF MATHEMATICAL CORRECTNESS OF IMPLEMENTED CONVOLUTIONS AND MATSUBARA SUMS
void test_matsubara_sums(int, double);
void test_convolutions  (int, double, double, string);




#endif // PHYSICS_H
