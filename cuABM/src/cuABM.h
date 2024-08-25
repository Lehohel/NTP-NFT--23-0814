#ifndef CURK4_H
#define CURK4_H

#include <string.h>
#include <iostream>
#include <complex>
#include <math.h>
#include <chrono>

#include <cuda_runtime.h>
#include <helper_cuda.h>
#include <cuComplex.h>

// typedef float2 Complex;
// using namespace std::complex_literals;

__global__ void oneStepRK4(cuDoubleComplex *d_field, int numberOfTemporalPoints, float stepSize, int offsetInput, int offsetOutput);
__device__ void diffEqEval(cuDoubleComplex field, cuDoubleComplex *res);
__device__ void RK4Helper(cuDoubleComplex fields[6], float stepSize);
__global__ void corrector(cuDoubleComplex *d_field, int numberOfTemporalPoints, float stepSize);
__device__ void correctorHelper(cuDoubleComplex *d_field, cuDoubleComplex *solver, double stepSize);
__host__ bool printfCUDAinfo(int argc, char *argv[]);
__host__ void readArgument(std::string name, int *var, int argc, char *argv[]);
__host__ void readArgument(std::string name, float *var, int argc, char *argv[]);
__host__ void readArgument(std::string name, std::string *var, int argc, char *argv[]);
__host__ cuDoubleComplex electricField(float t, float pulseDuration, float intensity);
__host__ cuDoubleComplex *allocateDeviceMemory(int numberOfTemporalPoints);
__host__ void copyFromHostToDevice(cuDoubleComplex *h_a, cuDoubleComplex *d_a, int numberOfTemporalPoints);
__host__ void executeKernel(cuDoubleComplex *d_field0, cuDoubleComplex *d_fieldFinal, int threadsPerBlock, int notp,
                            int nopp, float thickness, float centralWavelength, float timeStep);
__host__ void copyFromDeviceToHost(cuDoubleComplex *d_a, cuDoubleComplex *h_a, int numberOfTemporalPoints);
__host__ void deallocateMemory(cuDoubleComplex *d_a);
__host__ void cleanUpDevice();
__host__ cuDoubleComplex* analyticalSolution(cuDoubleComplex *field, int notp, float nonLinearRefIndex,
                                           float centralWavelength, float linearRefIndex, float thickness);
__host__ void saveData(std::string fileName, float *time, cuDoubleComplex *field0, cuDoubleComplex *fieldFinal,
                       cuDoubleComplex *analytical, int notp);
__host__ float calculateDifference(cuDoubleComplex *numeric, cuDoubleComplex *analytic, int notp);
__host__ int main(int argc, char *argv[]);

const float speedOfLight = 299.792458f; // Speed of light in vacuum in nm/fs
const float pi = 3.14159265359f;

__constant__ cuDoubleComplex d_diffEqFactor;

#endif