#include "stRk4.h"

void readArgument(std::string name, int *var, int argc, char *argv[])
{
  char *nameChar;
  if (checkCmdLineFlag(argc, (const char **)argv, name.c_str()))
  {
    getCmdLineArgumentString(argc, (const char **)argv, name.c_str(), &nameChar);
    *var = std::stoi(nameChar);
  }
}

void readArgument(std::string name, double *var, int argc, char *argv[])
{
  char *nameChar;
  if (checkCmdLineFlag(argc, (const char **)argv, name.c_str()))
  {
    getCmdLineArgumentString(argc, (const char **)argv, name.c_str(), &nameChar);
    *var = std::stof(nameChar);
  }
}

void readArgument(std::string name, std::string *var, int argc, char *argv[])
{
  char *nameChar;
  if (checkCmdLineFlag(argc, (const char **)argv, name.c_str()))
  {
    getCmdLineArgumentString(argc, (const char **)argv, name.c_str(), &nameChar);
    *var = (std::string)(nameChar);
  }
}

std::complex<double> electricField(double t, double pulseDuration, double intensity)
{
  std::complex<double> f;
  f = sqrt(intensity) * exp(-1 * pow(t / pulseDuration, 2) * 2 * log(2));
  return f;
}

std::complex<double> diffEqSource(std::complex<double> diffEqFactor, std::complex<double> field){
  return diffEqFactor *field * std::pow(std::abs(field), 2);
}

void solveDiffEq(std::complex<double> *field0, std::complex<double> *fieldFinal, std::complex<double> diffEqFactor, int notp, int nopp, double thickness){
  double stepSize = thickness / nopp;
  std::complex<double> *solver = (std::complex<double> *) malloc(sizeof(std::complex<double>)*5);
  for (int i = 0; i <notp; i++){
    solver[4] = field0[i];
    for (int j = 0; j<nopp; j++){
      solver[0] = diffEqSource(diffEqFactor, solver[4]);
      solver[1] = diffEqSource(diffEqFactor, solver[4] + solver[0] * (stepSize/2.0f));
      solver[2] = diffEqSource(diffEqFactor, solver[4] + solver[1] * (stepSize/2.0f));
      solver[3] = diffEqSource(diffEqFactor, solver[4] + solver[2] * stepSize);
      solver[4] += (solver[0] + 2.0 * solver[1] + 2.0 * solver[2] + solver[3]) * (stepSize / 6.0f);
    }
    fieldFinal[i] = solver[4];
  }
}

void saveData(std::string fileName, double *time, std::complex<double> *field, int notp){
  std::ofstream myfile;
  myfile.open("./data/"+fileName+".csv");
  myfile << "Time [fs],RK4 field real,RK4 field imaginary\n";
  for (int i = 0; i < notp; i += 1)
  {
    myfile << time[i] << ',' << field[i].real() << ',' << field[i].imag() << '\n';
  }
  myfile.close();
}

int main(int argc, char *argv[])
{  try
  {
    int numberOfTemporalPoints = 1024*2;
    int numberOfPropagationPoints = 100;
    double thickness = 1.0f; // mm
    double nonLinearRefIndex = 2.1e-4f; // cm2/TW
    double linearRefIndex = 1.45;
    double intensity = 10.0f;      // TW/cm2
    double pulseDuration = 30.0f; // fs
    double centralWavelength = 800.0f; // mm
    double timeStep = 0.1f; // fs
    std::string fileName = "output";

    // Read CL arguments
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
    double *time = (double *) malloc(sizeof(double)*numberOfTemporalPoints);
    std::complex<double> *field0 = (std::complex<double> *) malloc(sizeof(std::complex<double>)*numberOfTemporalPoints);
    std::complex<double> *fieldFinal = (std::complex<double> *) malloc(sizeof(std::complex<double>)*numberOfTemporalPoints);

    // Define grid and electric field
    for (int i = 0; i < numberOfTemporalPoints; i++)
    {
      time[i] = timeStep * i - numberOfTemporalPoints / 2 * timeStep;
      field0[i] = electricField(time[i], pulseDuration, intensity);
      // Filling up the result array with zeros for easier debugging
      fieldFinal[i] = 0;
    }

    // Calculate and copy the differential equation's factor into constant memory
    std::complex<double> diffEqFactor;
    diffEqFactor = -1i * nonLinearRefIndex * (2 * pi / centralWavelength * 1e6) / linearRefIndex; // in 1/mm unit
    std::cout<<"DiffEqFactor: " << diffEqFactor << std::endl;

    // Solve the differential equation
    solveDiffEq(field0, fieldFinal, diffEqFactor, numberOfTemporalPoints, numberOfPropagationPoints, thickness);

    // Save the data to a CSV file
    saveData(fileName, time, fieldFinal, numberOfTemporalPoints);

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