#include "cuABM.h"

__global__ void oneStepRK4(cuDoubleComplex *d_field, int numberOfTemporalPoints, float stepSize, int offsetInput, int offsetOutput)
{
  int i = blockDim.x * blockIdx.x + threadIdx.x;
  if (i < numberOfTemporalPoints)
  {
    // int tid = threadIdx.x;
    cuDoubleComplex dataPoints[6];
    dataPoints[0] = d_field[numberOfTemporalPoints*offsetInput+i];
    diffEqEval(dataPoints[0], &dataPoints[1]);
    diffEqEval(cuCadd(dataPoints[0], cuCmul(dataPoints[1], make_cuDoubleComplex(stepSize / 2, 0.0))), &dataPoints[2]);
    diffEqEval(cuCadd(dataPoints[0], cuCmul(dataPoints[2], make_cuDoubleComplex(stepSize / 2, 0.0))), &dataPoints[3]);
    diffEqEval(cuCadd(dataPoints[0], cuCmul(dataPoints[3], make_cuDoubleComplex(stepSize, 0.0))), &dataPoints[4]);
    RK4Helper(dataPoints, stepSize);
    d_field[numberOfTemporalPoints*offsetOutput+i] = dataPoints[5];
  }
}

__device__ void diffEqEval(cuDoubleComplex field, cuDoubleComplex *res){
  float intens = field.x * field.x + field.y * field.y;
  *res = cuCmul(make_cuDoubleComplex(intens, 0.0), cuCmul(field, d_diffEqFactor));
}

__device__ void RK4Helper(cuDoubleComplex dataPoints[6], float stepSize)
{
  dataPoints[5].x = dataPoints[0].x + (stepSize / 6) * (dataPoints[1].x + 2 * dataPoints[2].x + 2 * dataPoints[3].x + dataPoints[4].x);
  dataPoints[5].y = dataPoints[0].y + (stepSize / 6) * (dataPoints[1].y + 2 * dataPoints[2].y + 2 * dataPoints[3].y + dataPoints[4].y);
}

__global__ void corrector(cuDoubleComplex *d_field, int numberOfTemporalPoints, float stepSize){
  int i = blockDim.x * blockIdx.x + threadIdx.x;
  if (i < numberOfTemporalPoints)
  {
    cuDoubleComplex solver[4];
    // Corrector
        for (int k=0; k<4; k++){
            diffEqEval(d_field[k*numberOfTemporalPoints+i], &solver[k]);
        }
        // Shift history
        for (int k=0; k<3; k++){
            d_field[k*numberOfTemporalPoints+i] = d_field[(k+1)*numberOfTemporalPoints+i];
        }
        // Corrector
        correctorHelper(&d_field[3*numberOfTemporalPoints+i], solver, stepSize);
  }
}

__device__ void correctorHelper(cuDoubleComplex *d_field, cuDoubleComplex *solver, double stepSize){
  (*d_field).x += (9.0 * solver[3].x + 19.0 * solver[2].x - 5.0 * solver[1].x + solver[0].x) * (stepSize / 24.0f);
  (*d_field).y += (9.0 * solver[3].y + 19.0 * solver[2].y - 5.0 * solver[1].y + solver[0].y) * (stepSize / 24.0f);
}

__host__ bool printfCUDAinfo(int argc, char *argv[])
{
  int driverVersion, runtimeVersion;
  cudaDriverGetVersion(&driverVersion);
  cudaRuntimeGetVersion(&runtimeVersion);

  printf("  CUDA Driver  Version: %d.%d\n", driverVersion / 1000,
         (driverVersion % 100) / 10);
  printf("  CUDA Runtime Version: %d.%d\n", runtimeVersion / 1000,
         (runtimeVersion % 100) / 10);

  // Min spec is SM 1.0 devices
  bool bVal = checkCudaCapabilities(1, 0);
  printf("\n");
  return bVal;
}

__host__ void readArgument(std::string name, int *var, int argc, char *argv[])
{
  char *nameChar;
  if (checkCmdLineFlag(argc, (const char **)argv, name.c_str()))
  {
    getCmdLineArgumentString(argc, (const char **)argv, name.c_str(), &nameChar);
    *var = std::stoi(nameChar);
  }
}

__host__ void readArgument(std::string name, float *var, int argc, char *argv[])
{
  char *nameChar;
  if (checkCmdLineFlag(argc, (const char **)argv, name.c_str()))
  {
    getCmdLineArgumentString(argc, (const char **)argv, name.c_str(), &nameChar);
    *var = std::stof(nameChar);
  }
}

__host__ void readArgument(std::string name, std::string *var, int argc, char *argv[])
{
  char *nameChar;
  if (checkCmdLineFlag(argc, (const char **)argv, name.c_str()))
  {
    getCmdLineArgumentString(argc, (const char **)argv, name.c_str(), &nameChar);
    *var = (std::string)(nameChar);
  }
}

__host__ cuDoubleComplex electricField(float t, float pulseDuration, float intensity)
{
  cuDoubleComplex f;
  f.x = sqrt(intensity) * exp(-1 * pow(t / pulseDuration, 2) * 2 * log(2));
  f.y = 0;
  return f;
}

__host__ cuDoubleComplex * allocateDeviceMemory(int numberOfTemporalPoints)
{
  // Allocate the device input vector A
  cuDoubleComplex *d_a = NULL;
  size_t size = numberOfTemporalPoints * sizeof(cuDoubleComplex);

  cudaError_t err = cudaMalloc(&d_a, size);
  if (err != cudaSuccess)
  {
    fprintf(stderr, "Failed to allocate device vector (error code %s)!\n", cudaGetErrorString(err));
    exit(EXIT_FAILURE);
  }

  return d_a;
}

__host__ void copyFromHostToDevice(cuDoubleComplex *h_a, cuDoubleComplex *d_a, int numberOfTemporalPoints)
{
  size_t size = numberOfTemporalPoints * sizeof(cuDoubleComplex);

  cudaError_t err = cudaMemcpy(d_a, h_a, size, cudaMemcpyHostToDevice);

  if (err != cudaSuccess)
  {
    fprintf(stderr, "Failed to copy vector from host to device (error code %s)!\n", cudaGetErrorString(err));
    exit(EXIT_FAILURE);
  }
}

__host__ void executeKernel(cuDoubleComplex *d_field, int threadsPerBlock,
                            int notp, int nopp, float thickness, float centralWavelength, float timeStep)
{
  // Launch the search CUDA Kernel
  int blocksPerGrid = (notp + threadsPerBlock - 1) / threadsPerBlock;
  float stepSize = thickness / (nopp);
  // First step
  oneStepRK4<<<blocksPerGrid, threadsPerBlock>>>(d_field, notp, stepSize, 0, 1);
  // Second step
  oneStepRK4<<<blocksPerGrid, threadsPerBlock>>>(d_field, notp, stepSize, 1, 2);
  // Third step and more
  for (int i = 0; i < nopp-3; i++){
    oneStepRK4<<<blocksPerGrid, threadsPerBlock>>>(d_field, notp, stepSize, 2, 3);
    corrector<<<blocksPerGrid, threadsPerBlock>>>(d_field, notp, stepSize);
  }

  cudaError_t err = cudaGetLastError();

  if (err != cudaSuccess)
  {
    fprintf(stderr, "Failed to launch oneStepRK4 kernel (error code %s)!\n", cudaGetErrorString(err));
    exit(EXIT_FAILURE);
  }
}

__host__ void copyFromDeviceToHost(cuDoubleComplex *d_a, cuDoubleComplex *h_a, int numberOfTemporalPoints)
{
  size_t size = numberOfTemporalPoints * sizeof(cuDoubleComplex);

  cudaError_t err = cudaMemcpy(h_a, d_a, size, cudaMemcpyDeviceToHost);

  if (err != cudaSuccess)
  {
    fprintf(stderr, "Failed to copy vector from device to host (error code %s)!\n", cudaGetErrorString(err));
    exit(EXIT_FAILURE);
  }
}

__host__ void deallocateMemory(cuDoubleComplex *d_a)
{

  cudaError_t err = cudaFree(d_a);
  if (err != cudaSuccess)
  {
    fprintf(stderr, "Failed to free device vector (error code %s)!\n", cudaGetErrorString(err));
    exit(EXIT_FAILURE);
  }
}

// Reset the device and exit
__host__ void cleanUpDevice()
{
  // cudaDeviceReset causes the driver to clean up all state. While
  // not mandatory in normal operation, it is good practice.  It is also
  // needed to ensure correct operation when the application is being
  // profiled. Calling cudaDeviceReset causes all profile data to be
  // flushed before the application exits
  cudaError_t err = cudaDeviceReset();

  if (err != cudaSuccess)
  {
    fprintf(stderr, "Failed to deinitialize the device! error=%s\n", cudaGetErrorString(err));
    exit(EXIT_FAILURE);
  }
}

__host__ cuDoubleComplex* analyticalSolution(cuDoubleComplex *field, int notp, float nonLinearRefIndex,
                                           float centralWavelength, float linearRefIndex, float thickness)
{
  cuDoubleComplex *result = (cuDoubleComplex *)malloc(sizeof(cuDoubleComplex) * notp);
  float diffEqFactor = nonLinearRefIndex * (2 * pi / centralWavelength * 1e6) / linearRefIndex;
  for (int i = 0; i < notp; i++)
  {
    float intensity = field[i].y * field[i].y + field[i].x * field[i].x;
    float phi0 = atan2(field[i].y, field[i].x);
    float phi = phi0 - diffEqFactor * thickness * intensity;
    result[i].x = sqrt(intensity) * cos(phi);
    result[i].y = sqrt(intensity) * sin(phi);
  }
  return result;
}

__host__ void saveData(std::string fileName, float *time, cuDoubleComplex *field0, cuDoubleComplex *fieldFinal, cuDoubleComplex *analytical, int notp){
  std::ofstream myfile;
  myfile.open("./data/"+fileName+".csv");
  myfile << "Time [fs],Original field real,Original field imaginary,RK4 field real,RK4 field imaginary,Analytical real,Analytical imaginary\n";
  for (int i = 0; i < notp; i += 1)
  {
    myfile << time[i] << ',' << field0[i].x << ',' << field0[i].y << ',' << fieldFinal[i].x << ',' << fieldFinal[i].y << ',' << analytical[i].x << ',' << analytical[i].y << '\n';
  }
  myfile.close();
}

__host__ float calculateDifference(cuDoubleComplex *numeric, cuDoubleComplex *analytic, int notp)
{
  printf("Comparing analytical and numerical solutions\n");
  int totalsquaredifference = 0.0f;

  for (int i = 0; i < notp; i++)
  {
    cuDoubleComplex difference;
    difference.x = analytic[i].x - numeric[i].x;
    difference.y = analytic[i].y - numeric[i].y;
    totalsquaredifference += difference.x * difference.x + difference.y * difference.y;
  }

  float meanDifference = sqrt(totalsquaredifference) / notp;
  printf("meanDifference: %f\n", meanDifference);
  return meanDifference;
}

__host__ int main(int argc, char *argv[])
{  try
  {
    int threadsPerBlock = 256;
    int numberOfTemporalPoints = 1024*2;
    int numberOfPropagationPoints = 100;
    float thickness = 1.0f; // mm
    float nonLinearRefIndex = 2.1e-4f; // cm2/TW
    // float nonLinearRefIndex = 2.1e6f; // cm2/TW
    float linearRefIndex = 1.45;
    float intensity = 10.0f;      // TW/cm2
    float pulseDuration = 30.0f; // fs
    float centralWavelength = 800.0f; // mm
    float timeStep = 0.1f; // fs
    std::string fileName = "output";

    // Read CL arguments
    readArgument("threadsPerBlock", &threadsPerBlock, argc, argv);
    readArgument("notp", &numberOfTemporalPoints, argc, argv);
    readArgument("nopp", &numberOfPropagationPoints, argc, argv);
    readArgument("thickness", &thickness, argc, argv);
    readArgument("intensity", &intensity, argc, argv);
    readArgument("tau", &pulseDuration, argc, argv);
    readArgument("wl0", &centralWavelength, argc, argv);
    readArgument("dt", &timeStep, argc, argv);
    readArgument("output", &fileName, argc, argv);
    printf("\n\n");

    // Start the clock
    auto start = std::chrono::high_resolution_clock::now();

    // allocate memory fo the arrays
    float *time = (float *) malloc(sizeof(float)*numberOfTemporalPoints);
    cuDoubleComplex *field0 = (cuDoubleComplex *) malloc(sizeof(cuDoubleComplex)*numberOfTemporalPoints*4);
    cuDoubleComplex *fieldFinal = (cuDoubleComplex *) malloc(sizeof(cuDoubleComplex)*numberOfTemporalPoints);

    // Define grid and electric field
    for (int i = 0; i < numberOfTemporalPoints; i++)
    {
      time[i] = timeStep * i - numberOfTemporalPoints / 2 * timeStep;
      field0[i] = electricField(time[i], pulseDuration, intensity);
      // Filling up the result array with zeros for easier debugging
      fieldFinal[i].x = 0;
      fieldFinal[i].y = 0;
    }

    //Allocate device memory
    cuDoubleComplex *d_field0 = allocateDeviceMemory(numberOfTemporalPoints*4);

    // Calculate and copy the differential equation's factor into constant memory
    cuDoubleComplex diffEqFactor;
    diffEqFactor.x = 0;
    diffEqFactor.y = -1*nonLinearRefIndex * (2 * pi / centralWavelength * 1e6) / linearRefIndex; // in 1/mm unit
    cudaMemcpyToSymbol(d_diffEqFactor, &diffEqFactor, sizeof(cuDoubleComplex), 0, cudaMemcpyHostToDevice);

    // Copy arrray from host to device
    copyFromHostToDevice(field0, d_field0, numberOfTemporalPoints);

    // Execute the kernel
    executeKernel(d_field0, threadsPerBlock, numberOfTemporalPoints, numberOfPropagationPoints, thickness, centralWavelength, timeStep);

    // Copy the result back to the host
    copyFromDeviceToHost(d_field0+(3*numberOfTemporalPoints), fieldFinal, numberOfTemporalPoints);

    // Calculate the analytical solution
    cuDoubleComplex *analytical = analyticalSolution(field0, numberOfTemporalPoints, nonLinearRefIndex, centralWavelength, linearRefIndex, thickness);

    // Save the data to a CSV file
    saveData(fileName, time, field0, fieldFinal, analytical, numberOfTemporalPoints);

    // Deallocate memory
    deallocateMemory(d_field0);
    cleanUpDevice();
    free(time);
    free(field0);
    free(fieldFinal);

    // Stop the clock
    auto stop = std::chrono::high_resolution_clock::now();
    
    // Calculate the duration in seconds
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);

    // Output the time
    std::cout << "Time taken by function: "
         << duration.count() << " milliseconds" << std::endl;
  }
  catch (...)
  {
    std::cerr << "Program error! An unknow type of exception occurred. \n";
    std::cerr << "Aborting." << std::endl;

    exit(EXIT_FAILURE);
    return -1;
  }

  return 0;
}