## C++ compiler
CXX = g++

## C++ compiler flags
CXXFLAGS = -Wall -g

## Path for header files
INCLUDE_PATH = /usr/local/boost_1_71_0

## Path for libraries
LIBRARY_PATH =

## Name of libraries
LIBS = fftw3

## 
VPATH = src/basic:src/test

## Remove command
RM = -rm

## Target name
GOAL = a.out

## Object files
OBJS = src/RAW.o \
       src/AWT.o \
       src/structures.o \
       src/test/test.o \
       main.o
