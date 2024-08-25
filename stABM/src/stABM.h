#ifndef STABM_H
#define STABM_H

#include <string.h>
#include <iostream>
#include <complex>
#include <math.h>
#include <chrono>

#include <helper_string.h>

using namespace std::complex_literals;

void readArgument(std::string name, int *var, int argc, char *argv[]);
void readArgument(std::string name, double *var, int argc, char *argv[]);
void readArgument(std::string name, std::string *var, int argc, char *argv[]);
std::complex<double> diffEqSource(std::complex<double> diffEqFactor, std::complex<double> field);
void oneStepRK4(std::complex<double> field0, std::complex<double> *fieldFinal, std::complex<double> diffEqFactor, double stepSize);
void solveDiffEq(std::complex<double> *field0, std::complex<double> *fieldFinal, std::complex<double> diffEqFactor,
                 int notp, int nopp, double thickness);
std::complex<double> electricField(double t, double pulseDuration, double intensity);
int main(int argc, char *argv[]);

const double speedOfLight = 299.792458f; // Speed of light in vacuum in nm/fs
const double pi = 3.14159265359f;

#endif