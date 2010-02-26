/******************************************************************************/
/*                                                                            */
/*   Flood: An Open Source Neural Networks C++ Library                        */
/*   www.cimne.com/flood                                                      */
/*                                                                            */
/*   F E D   B A T C H   F E R M E N T E R   P R O B L E M   C L A S S        */
/*                                                                            */
/*   Roberto Lopez                                                            */
/*   International Center for Numerical Methods in Engineering (CIMNE)        */
/*   Technical University of Catalonia (UPC)                                  */
/*   Barcelona, Spain                                                         */
/*   E-mail: rlopez@cimne.upc.es                                              */
/*                                                                            */
/******************************************************************************/

#include <iostream>     
#include <fstream>     
#include <math.h>   

#include "FedBatchFermenterProblem.h"

#include "../Utilities/OrdinaryDifferentialEquations.h"

namespace Flood
{

// GENERAL CONSTRUCTOR

/// General constructor. It creates a fed batch fermenter problem 
/// objective functional associated to a multilayer perceptron.
/// It also initializes all the problem parameters to their default values:
///
/// <ul>
/// <li> Final time: 54.
/// <li> Initial cell mass concentration: 1.
/// <li> Initial substrate concentration: 150.
/// <li> Initial product concentration: 0.
/// <li> Initial broth volume: 10.
/// <li> Maximum volume: 200.
/// <li> Minimum feed rate: 0.
/// <li> Maximum feed rate: 12.
/// <li> Tolerance: 1.0e-6.
/// <li> Initial size: 1000.
/// </ul>
/// 
/// @param newMultilayerPerceptron Pointer to a multilayer perceptron 
/// object.
///
/// @see ObjectiveFunctional.

FedBatchFermenterProblem
::FedBatchFermenterProblem(MultilayerPerceptron* newMultilayerPerceptron)
: ObjectiveFunctional(newMultilayerPerceptron)
{
   // Control sentence

   int numberOfInputs = multilayerPerceptron->getNumberOfInputs();
   int numberOfOutputs = multilayerPerceptron->getNumberOfOutputs();

   if(numberOfInputs != 1 || numberOfOutputs != 1)
   {
      std::cerr << std::endl
               << "Flood Error: FedBatchFermenterProblem class." << std::endl
             << "FedBatchFermenterProblem(MultilayerPerceptron*) constructor." << std::endl
             << "Number of inputs and outputs in multilayer perceptron must be 1." << std::endl
             << std::endl;

      exit(1);
   }

   finalTime = 54.0;

   initialCellMassConcentration = 1.0;
   initialSubstrateConcentration = 150.0;
   initialProductConcentration = 0.0;
   initialBrothVolume = 10.0;

   fermenterVolume = 200.0;
   
   minimumFeedRate = 0.0;
   maximumFeedRate = 12.0;

   volumeErrorWeight = 1.0e-2;
   yieldWeight = 1.0e-9;

    yieldCoefficient = 0.1; // no units
    feedSubstrateConcentration = 150.0; // g/l

    kineticConstantMu0 =0.408; // 1/h
    kineticConstantEta0 = 1.0; // 1/h
    kineticConstantKp = 16.0; // g/l
    kineticConstantKpDash = 71.5; // g/l
    kineticConstantKs = 0.22; // g/l
    kineticConstantKsDash = 0.44; // g/l

   tolerance = 1.0e-15;
   initialSize = 5000;

   ordinaryDifferentialEquations.setTolerance(tolerance);
   ordinaryDifferentialEquations.setInitialSize(initialSize);

}


// DEFAULT CONSTRUCTOR

/// Default constructor. It creates a fed batch fermenter problem 
/// objective functional not associated to any multilayer perceptron. It
/// also initializes all the problem parameters to their default values:
///
/// <ul>
/// <li> Final time: 54.
/// <li> Initial cell mass concentration: 1.
/// <li> Initial substrate concentration: 150.
/// <li> Initial product concentration: 0.
/// <li> Initial broth volume: 10.
/// <li> Maximum volume: 200.
/// <li> Minimum feed rate: 0.
/// <li> Maximum feed rate: 12.
/// <li> Tolerance: 1.0e-6.
/// <li> Initial size: 1000.
/// </ul>
///
/// @see ObjectiveFunctional.

FedBatchFermenterProblem::FedBatchFermenterProblem(void) : ObjectiveFunctional()
{
   finalTime = 54.0;

   initialCellMassConcentration = 1.0;
   initialSubstrateConcentration = 150.0;
   initialProductConcentration = 0.0;
   initialBrothVolume = 10.0;

   fermenterVolume = 200.0;
   
   minimumFeedRate = 0.0;
   maximumFeedRate = 12.0;

   volumeErrorWeight = 1.0e-2;
   yieldWeight = 1.0e-9;

   tolerance = 1.0e-15;
   initialSize = 5000;
}


// DESTRUCTOR

/// Destructor.

FedBatchFermenterProblem::~FedBatchFermenterProblem(void)
{

}


// METHODS

// double getFinalTime(void) method

/// This method returns the set time for the fermentation process (hour).

double FedBatchFermenterProblem::getFinalTime(void)
{
   return(finalTime);
}


// double getInitialCellMassConcentration(void) method

/// This method returns the initial concentration of cell mass in the fermenter (gram/litre). 

double FedBatchFermenterProblem::getInitialCellMassConcentration(void)
{
   return(initialCellMassConcentration);
}


// double getInitialSubstrateConcentration(void) method

/// This method returns the initial concentration of substrate in the fermenter (gram/litre). 

double FedBatchFermenterProblem::getInitialSubstrateConcentration(void)
{
   return(initialSubstrateConcentration);
}


// double getInitialProductConcentration(void) method

/// This method returns the initial concentration of product in the fermenter (gram/litre). 

double FedBatchFermenterProblem::getInitialProductConcentration(void)
{
   return(initialProductConcentration);
}


// double getInitialBrothVolume(void) method

/// This method returns the initial volume of culture medium in the fermenter (litre). 

double FedBatchFermenterProblem::getInitialBrothVolume(void)
{
   return(initialBrothVolume);
}


// double getFermenterVolume(void) method

/// This method returns the volume of the fermenter (litre).

double FedBatchFermenterProblem::getFermenterVolume(void)
{
   return(fermenterVolume);
}


// double getMinimumFeedRate(void) method

/// This method returns the minimum feed rate (control action) to the fermenter (litre/hour).

double FedBatchFermenterProblem::getMinimumFeedRate(void)
{
   return(minimumFeedRate);
}


// double getMaximumFeedRate(void) method

/// This method returns the maximum feed rate (control action) to the fermenter (litre/hour).

double FedBatchFermenterProblem::getMaximumFeedRate(void)
{
   return(maximumFeedRate);
}


// double getVolumeErrorWeight(void) method

/// This method returns the weight for the final volume error term in the objective functional.

double FedBatchFermenterProblem::getVolumeErrorWeight(void)
{
   return(volumeErrorWeight);
}


// double getYieldWeight(void) method

/// This method returns the weight for the final yield term in the objective functional.

double FedBatchFermenterProblem::getYieldWeight(void)
{
   return(yieldWeight);
}


// double getTolerance(void) method

/// This method returns the tolerance value to be used in the Runge-Kutta-Fehlberg
/// method for evaluating the final volume error and yield.

double FedBatchFermenterProblem::getTolerance(void)
{
   return(tolerance);
}


// int getInitialSize(void) method

/// This method returns the number of points to be reserved when using the 
/// Runge-Kutta-Fehlberg method for evaluating the final volume error and yield.

int FedBatchFermenterProblem::getInitialSize(void)
{
   return(initialSize);
}


// void setFinalTime(double) method

/// This method sets a new final time for the fermentation process. 
///
/// @param newFinalTime Final fermentation time. 

void FedBatchFermenterProblem::setFinalTime(double newFinalTime)
{
   // Set final time

   finalTime = newFinalTime;
}


// void setInitialCellMassConcentration(double) method

/// This method sets a new initial value for the cell mass concentration (state variable).
///
/// @param newInitialCellMassConcentration Initial cell mass concentration. 

void FedBatchFermenterProblem::setInitialCellMassConcentration(double newInitialCellMassConcentration)
{
   initialCellMassConcentration = newInitialCellMassConcentration;
}


// void setInitialSubstrateConcentration(double) method

/// This method sets a new initial value for the substrate concentration (state variable).
///
/// @param newInitialSubstrateConcentration Initial substrate concentration. 

void FedBatchFermenterProblem::setInitialSubstrateConcentration(double newInitialSubstrateConcentration)
{
   initialSubstrateConcentration = newInitialSubstrateConcentration;
}


// void setInitialProductConcentration(double) method

/// This method sets a new initial value for the ethanol concentration (state variable).
///
/// @param newInitialProductConcentration Initial product concentration. 

void FedBatchFermenterProblem::setInitialProductConcentration(double newInitialProductConcentration)
{
   initialProductConcentration = newInitialProductConcentration;
}


// void setInitialBrothVolume(double) method

/// This method sets a new initial value for the broth volume (state variable).
///
/// @param newInitialBrothVolume Initial broth volume. 

void FedBatchFermenterProblem::setInitialBrothVolume(double newInitialBrothVolume)
{
   initialBrothVolume = newInitialBrothVolume;
}


// void setFermenterVolume(double) method

/// This method sets a new value for the volume of the fermenter.
///
/// @param newFermenterVolume Fermenter volume. 

void FedBatchFermenterProblem::setFermenterVolume(double newFermenterVolume)
{
   fermenterVolume = newFermenterVolume;
}


// void setMinimumFeedRate(double) method

/// This method sets a new minimum value for the feed rate to the fermenter. 
/// This is used as a lower bound for the control variable. 
///
/// @param newMinimumFeedRate Minimum feed rate value. 

void FedBatchFermenterProblem::setMinimumFeedRate(double newMinimumFeedRate)
{
   // Set minimum feed rate

   minimumFeedRate = newMinimumFeedRate;
}


// void setMaximumFeedRate(double) method

/// This method sets a new maximum value for the feed rate to the fermenter. 
/// This is used as an upper bound for the control variable. 
///
/// @param newMaximumFeedRate Maximum feed rate value. 

void FedBatchFermenterProblem::setMaximumFeedRate(double newMaximumFeedRate)
{
   // Set maximum feed rate

   maximumFeedRate = newMaximumFeedRate;
}


// void setVolumeErrorWeight(double) method

/// This method sets a new value for the weight of the volume error term (constraint)
/// in the objective functional expression.
///
/// @param newVolumeErrorWeight Volume error term weight value. 

void FedBatchFermenterProblem::setVolumeErrorWeight(double newVolumeErrorWeight)
{
   volumeErrorWeight = newVolumeErrorWeight;
}


// void setYieldWeight(double) method

/// This method sets a new value for the weight of the yield term (objective)
/// in the objective functional expression.
///
/// @param newYieldWeight Yield term weight value. 

void FedBatchFermenterProblem::setYieldWeight(double newYieldWeight)
{
   yieldWeight = newYieldWeight;
}


// void setTolerance(double) method

/// This method sets a new tolerance value to be used in the Runge-Kutta-Fehlberg
/// method for evaluating the final volume error and yield.
/// The tolerance of integration must be a small value greater than zero. 
///
/// @param newTolerance Tolerance in Runge-Kutta-Fehlberg method.

void FedBatchFermenterProblem::setTolerance(double newTolerance)
{
   // Control sentence

   if(newTolerance <= 0.0)
   {
      std::cerr << std::endl
                << "Flood Error: FedBatchFermenterProblem class." << std::endl
                << "void setTolerance(double) method." << std::endl
                << "The tolerance of integration must be greater than zero." << std::endl
                << std::endl;

      exit(1);
   }

   // Set tolerance 

   tolerance = newTolerance;

   ordinaryDifferentialEquations.setTolerance(tolerance);

}


// void setInitialSize(int) method

/// This method sets a new number of points to be reserved when using the 
/// Runge-Kutta-Fehlberg method for evaluating the final volume error and yield. 
/// The initial size must be a big value greater than zero,
/// otherwise resizing will be necessary. 
///
/// @param newInitialSize Number of points to reserve in Runge-Kutta-Fehlberg method.

void FedBatchFermenterProblem::setInitialSize(int newInitialSize)
{
   // Set initial size

   initialSize = newInitialSize;

   ordinaryDifferentialEquations.setInitialSize(initialSize);
}


// double getSpecificGrowthRate(double, double) method

/// This method returns the specific growth rate in the ferementer as a function of the 
/// product and substrate concentrations. 
///
/// @param productConcentration Actual concentration of ethanol in the ferementer. 
/// @param substrateConcentration Actual concentration of substrate in the ferementer. 

double FedBatchFermenterProblem
::getSpecificGrowthRate(double productConcentration, double substrateConcentration)
{
   double specificGrowthRate = 0.0;
   
   specificGrowthRate 
   = (kineticConstantMu0/(1.0+productConcentration/kineticConstantKp))
   *(substrateConcentration/(kineticConstantKs+substrateConcentration));

   return(specificGrowthRate);
}


// double getSpecificProductivity(double, double) method

/// This method returns the specific productivity in the ferementer as a function of the 
/// product and substrate concentrations. 
///
/// @param productConcentration Actual concentration of ethanol in the ferementer. 
/// @param substrateConcentration Actual concentration of substrate in the ferementer. 

double FedBatchFermenterProblem
::getSpecificProductivity(double productConcentration, double substrateConcentration)
{
   double specificProductivity = 0.0;

   specificProductivity 
   = (kineticConstantEta0/(1.0+productConcentration/kineticConstantKpDash))
   *(substrateConcentration/(kineticConstantKsDash+substrateConcentration));

   return(specificProductivity);
}


// void getCellMassConcentrationDot(double, double, double, double, double) method

/// This method depicts the state equation for the cell mass concentration in the 
/// fed batch fermenter problem. It returns the derivative of the 
/// cell mass concentration with respect to the time as a function of actual time, 
/// cell mass concentration, substrateConcentration, product concentration and
/// broth volume.
///
/// @param time Actual time.
/// @param cellMassConcentration Actual cell mass concentration.
/// @param substrateConcentration Actual substrate concentraion.
/// @param productConcentration Actual product concentration.
/// @param brothVolume Actual broth volume.

double FedBatchFermenterProblem::getCellMassConcentrationDot(double time, 
double cellMassConcentration, double substrateConcentration, double productConcentration, double brothVolume)
{
   double cellMassConcentrationDot = 0.0;

   // Obtain control

   double feedRate = 0.0;

   Vector<double> input(1, time);

   Vector<double> output = multilayerPerceptron->calculateOutput(input);

   // Bound control

   if(output[0] < minimumFeedRate)
   {
      feedRate = minimumFeedRate;
   }
   else if(output[0] > maximumFeedRate)
   {
      feedRate = maximumFeedRate;
   }
   else
   {
       feedRate = output[0];
   }

   // Calculate cell mass concentration derivative

    double specificGrowthRate = getSpecificGrowthRate(productConcentration, substrateConcentration);

   cellMassConcentrationDot = 
   specificGrowthRate*cellMassConcentration - feedRate*(cellMassConcentration/brothVolume);
    
   return(cellMassConcentrationDot);
}


// void getSubstrateConcentrationDot(double, double, double, double, double) method

/// This method depicts the state equation for the variable y2 in the 
/// fed batch fermenter problem. It returns the derivative of the 
/// variable y2 with respect to the time as a function of time, y1, y2, y3
/// and y4.
///
/// @param time Actual time.
/// @param cellMassConcentration Actual cell mass concentration.
/// @param substrateConcentration Actual substrate concentraion.
/// @param productConcentration Actual product concentration.
/// @param brothVolume Actual broth volume.

double FedBatchFermenterProblem::getSubstrateConcentrationDot(double time, 
double cellMassConcentration, double substrateConcentration, double productConcentration, double brothVolume)
{
   double substrateConcentrationDot = 0.0;

   // Obtain control 

   double feedRate = 0.0;

   Vector<double> input(1, time);

   Vector<double> output = multilayerPerceptron->calculateOutput(input);

   if(output[0] < minimumFeedRate)
   {
      feedRate = minimumFeedRate;
   }
   else if(output[0] > maximumFeedRate)
   {
      feedRate = maximumFeedRate;
   }
   else
   {
       feedRate = output[0];
   }

   // Calculate substrate concentration derivative

    double specificGrowthRate = getSpecificGrowthRate(productConcentration, substrateConcentration);

   substrateConcentrationDot = 
   -1.0*specificGrowthRate*cellMassConcentration/yieldCoefficient 
   + feedRate*(feedSubstrateConcentration-substrateConcentration)/brothVolume;

   return(substrateConcentrationDot);
}


// void getProductConcentrationDot(double, double, double, double, double) method

/// This method depicts the state equation for the variable y3 in the 
/// fed batch fermenter problem. It returns the derivative of the 
/// variable y3 with respect to the time as a function of time, y1, y2, y3
/// and y4.
///
/// @param time Actual time.
/// @param cellMassConcentration Actual cell mass concentration.
/// @param substrateConcentration Actual substrate concentraion.
/// @param productConcentration Actual product concentration.
/// @param brothVolume Actual broth volume.

double FedBatchFermenterProblem::getProductConcentrationDot(double time, 
double cellMassConcentration, double substrateConcentration, double productConcentration, double brothVolume)
{
   double productConcentrationDot = 0.0;

   // Obtain control

   double feedRate = 0.0;

   Vector<double> input(1, time);

   Vector<double> output = multilayerPerceptron->calculateOutput(input);

   // Bound control

   if(output[0] < minimumFeedRate)
   {
      feedRate = minimumFeedRate;
   }
   else if(output[0] > maximumFeedRate)
   {
      feedRate = maximumFeedRate;
   }
   else
   {
       feedRate = output[0];
   }

   // calculate product concentration derivative

    double specificProductivity = getSpecificProductivity(productConcentration, substrateConcentration);

    productConcentrationDot = specificProductivity*cellMassConcentration
   - feedRate*productConcentration/brothVolume;

   return(productConcentrationDot);
}


// void getBrothVolumeDot(double, double, double, double, double) method

/// This method depicts the state equation for the variable y4 in the 
/// fed batch fermenter problem. It returns the derivative of the 
/// variable y4 with respect to the time as a function of time, y1, y2, y3
/// and y4.
///
/// @param time Actual time.
/// @param cellMassConcentration Actual cell mass concentration.
/// @param substrateConcentration Actual substrate concentraion.
/// @param productConcentration Actual product concentration.
/// @param brothVolume Actual broth volume.

double FedBatchFermenterProblem::getBrothVolumeDot(double time, 
double cellMassConcentration, double substrateConcentration, double productConcentration, double brothVolume)
{
   double brothVolumeDot = 0.0;

   // Obtain control

   double feedRate = 0.0;
   
   Vector<double> input(1, time);

   Vector<double> output = multilayerPerceptron->calculateOutput(input);

   // Bound control

   if(output[0] < minimumFeedRate)
   {
      feedRate = minimumFeedRate;
   }
   else if(output[0] > maximumFeedRate)
   {
      feedRate = maximumFeedRate;
   }
   else
   {
       feedRate = output[0];
   }

   // Calculate broth volume derivative

   brothVolumeDot = feedRate;

   return(brothVolumeDot);
}


// double calculateEvaluation(void) method

/// This method returns the objective functional evaluation of a multilayer perceptron
/// for the fed batch fermenter problem. 

double FedBatchFermenterProblem::calculateEvaluation()
{
    // Control sentence 

   if(multilayerPerceptron == NULL)
   {
      std::cerr << std::endl
               << "Flood Error: FedBatchFermenterProblem class." << std::endl
              << "double calculateEvaluation(void) method." << std::endl
              << "Pointer to multilayer perceptron object cannot be NULL." << std::endl
              << std::endl;

        exit(1);
   }


   // Increment number of evaluations

   numberOfEvaluations++;

   // Evaluate objective functional

   double evaluation = 0.0;

   // Solve state equations 

   Vector<double> time;
   Vector<double> cellMassConcentration;
   Vector<double> substrateConcentration;
   Vector<double> productConcentration;
   Vector<double> brothVolume;

   double initialTime = 0.0;

    int numberOfPoints = 
    ordinaryDifferentialEquations.calculateRungeKuttaFehlbergIntegral(*this,
    time, cellMassConcentration, substrateConcentration, productConcentration, brothVolume,
    &FedBatchFermenterProblem::getCellMassConcentrationDot,
    &FedBatchFermenterProblem::getSubstrateConcentrationDot,
    &FedBatchFermenterProblem::getProductConcentrationDot,
    &FedBatchFermenterProblem::getBrothVolumeDot,
    initialTime, finalTime, 
    initialCellMassConcentration, initialSubstrateConcentration, initialProductConcentration, initialBrothVolume);

   // Obtain control vector

   Vector<double> feedRate(numberOfPoints);

   Vector<double> input(1), output(1);

   for(int i = 0; i < numberOfPoints; i++)
   {
      input[0] = time[i];
      output = multilayerPerceptron->calculateOutput(input);
      feedRate[i] = output[0];
   }   

   // Calculate volume error

   double volumeError = 0.0;

   if(brothVolume[numberOfPoints-1] <= fermenterVolume)
   {
      volumeError = 0.0;
   }
   else
   {
      volumeError = brothVolume[numberOfPoints-1] - fermenterVolume;
   }

   // Calculate yield

   double yield = productConcentration[numberOfPoints-1]*brothVolume[numberOfPoints-1];

   // Calculate evaluation

   evaluation = volumeErrorWeight*pow(volumeError,2)
             -1.0*yieldWeight*pow(yield,2);
      
   return(evaluation);
}


// void saveResults(char*) method

/// This method saves the values of the independent and dependent 
/// variables for the brachistochrone problem to a data file. 
///
/// <ol>
/// <li> Time.
/// <li> Cell mass concentration value (state variable 1).
/// <li> Substrate concentration value (state variable 2).
/// <li> Product concentration value (state variable 3).
/// <li> Broth volume value (state variabl 4).
/// <li> Specific growth rate.
/// <li> Specific productivity.
/// <li> Feed rate (control variable).
/// </ol>
///
/// @param filename Filename.

void FedBatchFermenterProblem::saveResults(char* filename)
{
   Vector<double> input(1), output(1);
   
   // Solve state equations 

   Vector<double> time;
   Vector<double> cellMassConcentration;
   Vector<double> substrateConcentration;
   Vector<double> productConcentration;
   Vector<double> brothVolume;

   double initialTime = 0.0;

   int numberOfPoints = 
    ordinaryDifferentialEquations.calculateRungeKuttaFehlbergIntegral(*this,
   time, cellMassConcentration, substrateConcentration, productConcentration, brothVolume,
   &FedBatchFermenterProblem::getCellMassConcentrationDot,
   &FedBatchFermenterProblem::getSubstrateConcentrationDot,
   &FedBatchFermenterProblem::getProductConcentrationDot,
   &FedBatchFermenterProblem::getBrothVolumeDot,
   initialTime, finalTime, 
   initialCellMassConcentration, initialSubstrateConcentration, initialProductConcentration, initialBrothVolume);

   // Control vector

   Vector<double> feedRate(numberOfPoints);

   for(int i = 0; i < numberOfPoints; i++)
   {
      // Obtain control 

      input[0] = time[i];
      output = multilayerPerceptron->calculateOutput(input);

      // Bound control

      if(output[0] < minimumFeedRate)
      {
         feedRate[i] = minimumFeedRate;
      }
      else if(output[0] > maximumFeedRate)
      {
         feedRate[i] = maximumFeedRate;
      }
      else
      {
            feedRate[i] = output[0];
      }
   }   

   // Obtain yield

   double yield = productConcentration[numberOfPoints-1]*brothVolume[numberOfPoints-1];

   // Write file

   std::fstream file; 

   file.open(filename, std::ios::out);

   if(!file.is_open())
   {
      std::cout << std::endl 
               << "Cannot open fed batch fermenter problem results data file."
                << std::endl;
      
      exit(1);
   }
   else
   {
      std::cout << std::endl
               << "Saving fed batch fermenter problem results to data file..."
                << std::endl;
   }

   file << "% Yield (g): " << yield << std::endl
       << "% Batch time (h): " << finalTime << std::endl
          << "%" << std::endl;

   file << "% 1 - Time" << std::endl
        << "% 2 - Cell mass concentration (g/l)" << std::endl
       << "% 3 - Substrate concentration (g/l)" << std::endl
       << "% 4 - Product concentration (g/l)" << std::endl
       << "% 5 - Broth volume (l)" << std::endl
       << "% 6 - Specific growth rate" << std::endl
       << "% 7 - Specific productivity" << std::endl
       << "% 8 - Feed rate -control- (g/l)" << std::endl
       << std::endl;


   double specificGrowthRate = 0.0;
   double specificProductivity = 0.0;

   for(int i = 0; i < numberOfPoints; i++)
   {
      specificGrowthRate = getSpecificGrowthRate(productConcentration[i], substrateConcentration[i]);
      specificProductivity = getSpecificProductivity(productConcentration[i], substrateConcentration[i]);

      file << time[i] << " " 
           << cellMassConcentration[i] << " "
           << substrateConcentration[i] << " "
           << productConcentration[i] << " "
           << brothVolume[i] << " "
           << specificGrowthRate << " "
           << specificProductivity << " "
           << feedRate[i] << std::endl;   
   }   

   file.close();
}


// void print(void) method

/// This method prints to the screen useful information about the 
/// fed batch fermenter problem objective functional of a multilayer
/// perceptron. 

void FedBatchFermenterProblem::print()
{      
   // Solve state equations 

   Vector<double> time;
   Vector<double> cellMassConcentration;
   Vector<double> substrateConcentration;
   Vector<double> productConcentration;
   Vector<double> brothVolume;

    double initialTime = 0.0;

   int numberOfPoints = 
   ordinaryDifferentialEquations.calculateRungeKuttaFehlbergIntegral(*this,
   time, cellMassConcentration, substrateConcentration, productConcentration, brothVolume,
   &FedBatchFermenterProblem::getCellMassConcentrationDot,
   &FedBatchFermenterProblem::getSubstrateConcentrationDot,
   &FedBatchFermenterProblem::getProductConcentrationDot,
   &FedBatchFermenterProblem::getBrothVolumeDot,
   initialTime, finalTime, 
   initialCellMassConcentration, initialSubstrateConcentration, initialProductConcentration, initialBrothVolume);

   double volumeError = 0.0;
   
   if(brothVolume[numberOfPoints-1] <= fermenterVolume)
   {
      volumeError = 0.0;
   }
   else
   {
      volumeError = brothVolume[numberOfPoints-1] - fermenterVolume;
   }

   // Calculate yield

   double yield = productConcentration[numberOfPoints-1]*brothVolume[numberOfPoints-1];

   std::cout << "Volume error: " << volumeError << std::endl
             << "Yield: " << yield << std::endl;
}

}


// Flood: An Open Source Neural Networks C++ Library.
// Copyright (C) 2005-2008 Roberto Lopez 
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or any later version.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.

// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
