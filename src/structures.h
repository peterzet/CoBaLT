#ifndef STRUCTURES_H
#define STRUCTURES_H

#include "AWT.h"
#include "Calculator.h"

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
#include <stdlib.h>



//// Acronym AWT stands for ArrayWithTails

struct auxiliary
{
    public:
        auxiliary();
        ~auxiliary();

        AWT FD;
        AWT BE;
        AWT K3;

        AWT aux1;
        AWT aux2;
        AWT aux3;
        AWT aux4;
        AWT aux5;
        AWT aux6;
        AWT aux7;

    void initialize(int,double,double);


};

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
struct info
{
    public:
        double delta, gammaS, phi, gap;
        double kT, U, h, mu;
        string model;
        // intialization of iteration parameters
        double kT_min, kT_increment, kT_max;
        double U_min, U_increment, U_max;
        double h_min, h_increment, h_max;
        double x_min, x_increment, x_max;

        // precisons
        int kT_precision, U_precision, h_precision, x_precision;

        // iterations
        int kT_iterations, U_iterations, h_iterations, x_iterations;

        // intialization of variables for mesh properties
        double xMax;
        int n;
        // intialization of variables for printing out txt files
        int display, print_mode;
        double range;
        string output_mode, physics;
        // previously initialized variables are set to values from the input file


        void default_input();
        void import(string input);
        void precisions();
        void iterations();

};

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
struct vect
{
    public:
    complex<double> p;
    complex<double> m;
};

void cout_vect(vect );
void output_vect(vect, string);

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
struct matrix
{
    public:
    complex<double> pp;
    complex<double> pm;
    complex<double> mp;
    complex<double> mm;
};

void set_matrix_real(matrix &, complex<double> );
void set_matrix(matrix &,matrix & );
void aver_matrix(matrix &,matrix & );
void cout_matrix(matrix);
void output_matrix(matrix, string);

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
struct matrixAWT
{
    public:
    AWT pp;
    AWT pm;
    AWT mp;
    AWT mm;
};

void set_matrixAWT(matrixAWT &, double);
void copy_matrixAWT(matrixAWT &, matrixAWT &);
void aver_matrixAWT(matrixAWT &, matrixAWT &);
void export_matrixAWT(matrixAWT &, double, string, info);

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
struct tensor
{
    public:
    complex<double> ppp;

    complex<double> ppm;
    complex<double> pmp;
    complex<double> mpp;

    complex<double> mmp;
    complex<double> mpm;
    complex<double> pmm;

    complex<double> mmm;
 };

void cout_tensor(tensor);
void output_tensor(tensor &, string);
void copy_tensor(tensor &, tensor &);
void aver_tensor(tensor &, tensor &);
void multiply_tensor(tensor &, complex<double>);
void integrate_tensor(AWT &, tensor &, int, int, int);
void subtract_tensors(tensor &, tensor &, tensor &);


#endif // STRUCTURES_H
