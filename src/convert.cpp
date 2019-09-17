#ifndef CONVERTING_H_INCLUDED
#define CONVERTING_H_INCLUDED


#include "convert.h"
#include <sstream>

#include <boost/lexical_cast.hpp>
#include <boost/filesystem.hpp>
#include <iomanip>

using namespace std;

/*

// PREVOD CHAR NA INT
int CharToInt(char vloz)
{
    int numero;
    stringstream ss;
    ss << vloz;
    ss >> numero;
    return numero;
}

char StringToChar(string vloz)
{
    char numero;
    stringstream ss;
    ss << vloz;
    ss >> numero;
    return numero;
}

// PREVOD INT NA STR
string IntToStr(int vloz)
{
    string numero;
    stringstream ss;
    ss << vloz;
    ss >> numero;
    return numero;
}

// PREVOD INT NA STR
string DoubleToStr(double vloz)
{
    string numero;
    stringstream ss;
    ss << vloz;
    ss >> numero;
    return numero;
}

string SignedIntToStr(signed int vloz)
{
    string numero;
    stringstream ss;
    ss << vloz;
    ss >> numero;
    return numero;
}
*/

string double_to_string(double U, int precision)
{
        stringstream stream;
        stream  << fixed << setprecision(precision) << U;
        return stream.str();
}

string int_to_string(int U, int precision)
{
        stringstream stream;
        stream  << fixed << setprecision(precision) << U;
        return stream.str();
}

#endif // CONVERTING_H_INCLUDED

