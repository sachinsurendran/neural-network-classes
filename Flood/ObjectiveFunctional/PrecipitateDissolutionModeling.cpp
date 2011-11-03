/****************************************************************************************************************/
/*                                                                                                              */
/*   Flood: An Open Source Neural Networks C++ Library                                                          */
/*   www.cimne.com/flood                                                                                        */
/*                                                                                                              */
/*   P R E C I P I T A T E   D I S S O L U T I O N   M O D E L I N G   C L A S S                                */
/*                                                                                                              */
/*   Roberto Lopez                                                                                              */ 
/*   International Center for Numerical Methods in Engineering (CIMNE)                                          */
/*   Technical University of Catalonia (UPC)                                                                    */
/*   Barcelona, Spain                                                                                           */
/*   E-mail: rlopez@cimne.upc.es                                                                                */ 
/*                                                                                                              */
/****************************************************************************************************************/

#include <iostream>
#include <fstream>
#include <algorithm>
#include <functional>
#include <math.h>

#include "PrecipitateDissolutionModeling.h"     

namespace Flood
{

// GENERAL CONSTRUCTOR

/// General constructor. It creates a dissolution modeling objective functional associated to a multilayer 
/// perceptron. 
/// It also initializes all the problem parameters to their default values:
///
/// <ul>
/// <li> Number of samples = 0;
/// <li> Reference time = 16.0; // s
/// <li> Universal gas constant = 8.314e-3; // kJ / mol K
/// <li> Reference temperature = 623.0; // K
/// <li> Minimum Vickers hardness = 73.8;
/// <li> Maximum Vickers Hardness = 203.1;
/// <li> Minkowski parameter = 2.0;
/// <li> Minkowski error weight = 1.0;
/// <li> Regularization weight = 0.0;
/// </ul>
///
/// @param newMultilayerPerceptron Pointer to a multilayer perceptron 
/// object.

PrecipitateDissolutionModeling::PrecipitateDissolutionModeling(MultilayerPerceptron* newMultilayerPerceptron)       
: ObjectiveFunctional(newMultilayerPerceptron)
{
   // Control sentence

   int numberOfInputs = multilayerPerceptron->getNumberOfInputs();
   int numberOfOutputs = multilayerPerceptron->getNumberOfOutputs();

   if(numberOfInputs != 1 || numberOfOutputs != 1)
   {
      std::cerr << std::endl
                << "Flood Error: PrecipitateDissolutionModeling class." << std::endl
                << "PrecipitateDissolutionModeling(MultilayerPerceptron*) constructor." << std::endl
                << "Number of inputs and outputs in multilayer perceptron must be 0." << std::endl
                << std::endl;

      exit(1);
   }

   numberOfSamples = 0;

   universalGasConstant = 8.314e-3; // kJ / mol K


   referenceTemperature = 0.0; // K
   referenceTime = 0.0; // s

   minimumVickersHardness = 0.0;
   maximumVickersHardness = 0.0;

   minkowskiParameter = 2.0;

   minkowskiErrorWeight = 1.0;
   regularizationWeight = 0.0;
}


// DEFAULT CONSTRUCTOR

/// Default constructor. It creates a thermal properties of metals estimation problem objective functional not 
/// associated to any multilayer perceptron. 
/// It also initializes all the problem parameters to their default values:
///
/// <ul>
/// <li> Number of samples = 0;
/// <li> Reference time = 16.0; // s
/// <li> Universal gas constant = 8.314e-3; // kJ / mol K
/// <li> Reference temperature = 623.0; // K
/// <li> Minimum Vickers hardness = 73.8;
/// <li> Maximum Vickers Hardness = 203.1;
/// <li> Minkowski parameter = 2.0;
/// <li> Minkowski error weight = 1.0;
/// <li> Regularization weight = 0.0;
/// </ul>

PrecipitateDissolutionModeling::PrecipitateDissolutionModeling(void) : ObjectiveFunctional()
{
   numberOfSamples = 0;

   referenceTime = 16.0; // s
   universalGasConstant = 8.314e-3; // kJ / mol K
   referenceTemperature = 623.0; // K

   minimumVickersHardness = 73.8;
   maximumVickersHardness = 203.1;

   minkowskiParameter = 2.0;

   minkowskiErrorWeight = 1.0;
   regularizationWeight = 0.0;
}


// DESTRUCTOR

/// Destructor. 

PrecipitateDissolutionModeling::~PrecipitateDissolutionModeling(void) 
{

}


// METHODS


// double getMinimumVickersHardness(void) method

/// This method returns the final hardness after complete dissolution of precipitates of an aluminium alloy. 
///
/// @see getMaximumVickersHardness.

double PrecipitateDissolutionModeling::getMinimumVickersHardness(void)
{
   return(minimumVickersHardness);
}


// double getMaximumVickersHardness(void) method

/// This method returns the initial in the hardened state of an aluminium alloy. 
///
/// @see getMaximumVickersHardness.

double PrecipitateDissolutionModeling::getMaximumVickersHardness(void)
{
   return(maximumVickersHardness);
}


// double getReferenceTemperature(void)

/// This method returns the temperature of reference associated to the time of reference for complete dissolution.
///
/// @see getReferenceTime.

double PrecipitateDissolutionModeling::getReferenceTemperature(void)
{
   return(referenceTemperature);
}


// double getReferenceTime(void)

/// This method returns the reference time for complete dissolution of the precipitates at the reference 
/// temperature.
///
/// @see getReferenceTemperature.

double PrecipitateDissolutionModeling::getReferenceTime(void)
{
   return(referenceTime);
}


// double getMinkowskiParameter(void) method

/// This method returns the Minkowski parameter (R-exponent) in the objective functional expression.

double PrecipitateDissolutionModeling::getMinkowskiParameter(void)
{
   return(minkowskiParameter);
}


// double getMinkowskiErrorWeight(void) method

/// This method returns the weight of the Minkowski error term in the objective functional expression.

double PrecipitateDissolutionModeling::getMinkowskiErrorWeight(void)
{
   return(minkowskiErrorWeight);
}


// double getRegularizationWeight(void) method

/// This method returns the weight of the regularization term (biases and synaptic weights decay) in the 
/// objective functional expression.

double PrecipitateDissolutionModeling::getRegularizationWeight(void)
{
   return(regularizationWeight);
}


// void setMinimumVickersHardness(double) method

/// This method sets a new final hardness after complete dissolution of precipitates of an aluminium alloy. 
///
/// @param newMinimumVickersHardness Minimum Vickers hardness value. 
///
/// @see setMaximumVickersHardness.

void PrecipitateDissolutionModeling::setMinimumVickersHardness(double newMinimumVickersHardness)
{
   minimumVickersHardness = newMinimumVickersHardness;
}


// void setMaximumVickersHardness(double) method

/// This method sets a new initial in the hardened state of an aluminium alloy. 
///
/// @param newMaximumVickersHardness Maximum Vickers hardness value. 
///
/// @see getMaximumVickersHardness.

void PrecipitateDissolutionModeling::setMaximumVickersHardness(double newMaximumVickersHardness)
{
   maximumVickersHardness = newMaximumVickersHardness;
}


// void setReferenceTemperature(double) method

/// This method sets a new temperature of reference associated to the time of reference for complete dissolution.
///
/// @param newReferenceTemperature Reference temperature value. 
///
/// @see setReferenceTime.

void PrecipitateDissolutionModeling::setReferenceTemperature(double newReferenceTemperature)
{
   referenceTemperature = newReferenceTemperature;
}


// void setReferenceTime(double) method

/// This method sets a new reference time for complete dissolution of the precipitates at the reference 
/// temperature.
///
/// @param newReferenceTime Reference time value. 
///
/// @see setReferenceTemperature.

void PrecipitateDissolutionModeling::setReferenceTime(double newReferenceTime)
{
   referenceTime = newReferenceTime;
}


// void setMinkowskiParameter(double) method

/// This method sets a new Minkowski parameter (R-exponent) in the objective functional expression.
/// This value must be comprised between 1 and 2. 
///
/// @param newMinkowskiParameter Minkowski R-exponent value. 

void PrecipitateDissolutionModeling::setMinkowskiParameter(double newMinkowskiParameter)
{
   // Control sentence

   if(newMinkowskiParameter < 1.0 || newMinkowskiParameter > 2.0)
   {
      std::cerr << std::endl
                << "Flood Error. PrecipitateDissolutionModeling class." << std::endl
                << "void setMinkowskiParameter(double) method." << std::endl
                << "The Minkowski parameter must be comprised between 1 and 2" << std::endl
                << std::endl;
    
      exit(1);
   }

   // Set Minkowski parameter

   minkowskiParameter = newMinkowskiParameter;
}


// void setMinkowskiErrorWeight(double) method

/// This method sets a new weight value for the Minkowski error term in the objective functional expression.
///
/// @param newMinkowskiErrorWeight Minkowski R-value. 

void PrecipitateDissolutionModeling::setMinkowskiErrorWeight(double newMinkowskiErrorWeight)
{
   // Control sentence 

   if(newMinkowskiErrorWeight <= 0.0)
   {
      std::cout << std::endl
                << "Flood Error PrecipitateDissolutionModeling class." << std::endl
                << "void setMinkowskiErrorWeight(double) method." << std::endl
                << "The Minkowski error term weight must be greater than zero." 
                << std::endl << std::endl;

      exit(1);
   }

   // Set Minkowski error term weight

   minkowskiErrorWeight = newMinkowskiErrorWeight;
}


// void setRegularizationWeight(double) method

/// This method sets a new weight value for the regularization term (biases and synaptic weights decay)
/// in the objective functional expression.
///
/// @param newRegularizationWeight Regularization term weight value. 

void PrecipitateDissolutionModeling::setRegularizationWeight(double newRegularizationWeight)
{
   // Control sentence 

   if(newRegularizationWeight < 0.0)
   {
      std::cout << std::endl
                << "Flood Error PrecipitateDissolutionModeling class." << std::endl
                << "void setRegularizationWeight(double) method." << std::endl
                << "The regularization term weight must be equal or greater than zero." 
                << std::endl << std::endl;

      exit(1);
   }

   // Set Minkowski error term weight

   regularizationWeight = newRegularizationWeight;
}


// void loadVickersHardnessTest(char*)

/// This method loads the material properties and the Vickers hardness test 
/// of an aluminium alloy from a data file:
///
/// <ul>
/// <li> Minimum Vickers hardness.
/// <li> Maximum Vickers hardness.
/// <li> Reference time (s). 
/// <li> Reference temperature (ºC). 
/// <li> Experimental data (Temperature (ºC) - Time (s) - Vickers hardness).
/// </ul> 
///
/// Please mind about the file format. 
///
/// @param filename Filename.
///
/// @see save(char*).

void PrecipitateDissolutionModeling::loadVickersHardnessTest(char* filename)
{
   std::fstream file;

   file.open(filename, std::ios::in);

   if(!file.is_open())
   {
      std::cout << std::endl
                << "void loadVickersHardnessTest(char*) method." 
                << std::endl
                << "Cannot open Vickers hardness test data file."  << std::endl
                << std::endl;

      exit(1);
   }
   else
   {
      if(display)
      {
         std::cout << std::endl
                   << "Loading Vickers harndess test from file..." << std::endl;
      }
   }


   std::string word;

   // Minimum Vickers hardness

   while(word != "MinimumVickersHardness:")
   {
      file >> word;
   }
 
   file >> minimumVickersHardness;

   // Maximum Vickers hardness

   file >> word;

   file >> maximumVickersHardness;

   // Reference time

   file >> word;

   file >> referenceTime;

   // Reference temperature

   file >> word;

   file >> referenceTemperature;

   referenceTemperature += 273.16;

   // Number of samples

   file >> word;

   file >> numberOfSamples;

   // Experimental data

   file >> word;

   Vector<double> newTemperatureData(numberOfSamples);
   Vector<double> newTimeData(numberOfSamples);
   Vector<double> newVickersHardnessData(numberOfSamples);

   for(int i = 0; i < numberOfSamples; i++)
   {
      file >> newTemperatureData[i];
      file >> newTimeData[i];
      file >> newVickersHardnessData[i];

     newTemperatureData[i] += 273.16;
   }

   temperatureData = newTemperatureData;
   timeData = newTimeData;
   vickersHardnessData = newVickersHardnessData;

   file.close();
}


// double getFullDissolutionTime(double) method

/// This method returns the time for complete dissolution of an aluminium alloy at a given temperature. 
///
/// @param temperature Temperature value at which the full dissolution time is to be obtained. 

double PrecipitateDissolutionModeling::getFullDissolutionTime(double temperature)
{
   double fullDissolutionTime = 0.0;

   double effectiveActivationEnergy = multilayerPerceptron->getSingleIndependentParameter(0);

   fullDissolutionTime = referenceTime
   *exp((effectiveActivationEnergy/universalGasConstant)*(1.0/temperature-1.0/referenceTemperature));

   return(fullDissolutionTime);
}


// double getVolumetricFraction(double) method

/// This method returns the volume fraction of precipitates as a function of the actual Vickers hardness. 
///
/// @param vickersHardness Vickers hardness value. 

double PrecipitateDissolutionModeling::getVolumetricFraction(double vickersHardness)
{
   double volumetricFraction 
   = (vickersHardness - minimumVickersHardness)/(maximumVickersHardness - minimumVickersHardness);

   return(volumetricFraction);
}


// double getVickersHardness(double) method

/// This method returns the Vickers hardness as a function of the actual volume fraction of precipitates. 
///
/// @param volumetricFraction Volumetric fraction value. 

double PrecipitateDissolutionModeling::getVickersHardness(double volumetricFraction)
{
   double vickersHardness = minimumVickersHardness 
   + (maximumVickersHardness-minimumVickersHardness)*(1.0 - volumetricFraction);

   return(vickersHardness);
}


// Vector<double> getNormalizedTimeData(void) method

/// This method returns a vector containing the normalized time (log(t/t*)) for the data samples in the Vickers 
/// hardness test. 

Vector<double> PrecipitateDissolutionModeling::getNormalizedTimeData(void)
{
   Vector<double> normalizedTimeData(numberOfSamples);

   double fullDissolutionTime = 0.0;

   for(int i = 0; i < numberOfSamples; i++)
   {
      fullDissolutionTime = getFullDissolutionTime(temperatureData[i]); 

      normalizedTimeData[i] = log10(timeData[i]/fullDissolutionTime);
   }

   return(normalizedTimeData);
}


// Vector<double> getVolumetricFractionData(void) method

/// This method returns a vector containing the volume fraction of hardening precipitates (1-f/f0) for the data 
/// samples in the Vickers hardness test. 

Vector<double> PrecipitateDissolutionModeling::getVolumetricFractionData(void)
{
   Vector<double> volumetricFractionData(numberOfSamples);

   Vector<double> normalizedTimeData = getNormalizedTimeData();

   double volumetricFraction = 0.0;

   for(int i = 0; i < numberOfSamples; i++)
   {
      volumetricFraction = getVolumetricFraction(vickersHardnessData[i]);
    
      volumetricFractionData[i] = 1.0 - volumetricFraction;
   }

   return(volumetricFractionData);
}


// double calculateEvaluation(void) method

/// This method returns the evaluation of a multilayer perceptron for modeling the hardening precipitates 
/// dissolution in fully hardened aluminium alloys. 
/// The objective functional here is composed by a Minkowski error term and a regularization term consisting on 
/// the norm of the bias and synaptic weights vector. 

double PrecipitateDissolutionModeling::calculateEvaluation()
{
    // Control sentence 

   if(multilayerPerceptron == NULL)
   {
      std::cerr << std::endl
                << "Flood Error: PrecipitateDissolutionModeling class." << std::endl
                << "double calculateEvaluation(void) method." << std::endl
                << "Pointer to multilayer perceptron object cannot be NULL." << std::endl
                << std::endl;

        exit(1);
   }


    numberOfEvaluations++;

   double evaluation = 0.0;

   Vector<double> normalizedTimeData = getNormalizedTimeData();
   Vector<double> volumetricFractionData = getVolumetricFractionData();

   Vector<double> input(1);
   Vector<double> output(1);

   Vector<double> volumetricFractionModel(numberOfSamples);

   double minkowskiError = 0.0;

   for(int i = 0; i < numberOfSamples; i++)
   {
      input[0] = normalizedTimeData[i];
      output = multilayerPerceptron->calculateOutput(input);
      
      volumetricFractionModel[i] = output[0];

      minkowskiError += pow(fabs(volumetricFractionData[i] - volumetricFractionModel[i]), minkowskiParameter);
   }

   minkowskiError = pow(minkowskiError/(double)numberOfSamples, 1.0/minkowskiParameter);
  
   Vector<double> neuralParameters = multilayerPerceptron->getNeuralParameters();

   double neuralParametersNorm = neuralParameters.calculateNorm();

   evaluation = minkowskiErrorWeight*minkowskiError + regularizationWeight*neuralParametersNorm;

   return(evaluation);
}


// void printVickersHardnessTest(void) method

/// This method prints the material properties and the Vickers hardness test of an aluminium alloy to the screen:
///
/// <ul>
/// <li> Minimum Vickers hardness.
/// <li> Maximum Vickers hardness.
/// <li> Reference time (s). 
/// <li> Reference temperature (ºC). 
/// <li> Experimental data (Temperature (ºC) - Time (s) - Vickers hardness).
/// </ul> 
///
/// @see load(char*).
/// @see save(char*).

void PrecipitateDissolutionModeling::printVickersHardnessTest(void)
{
   std::cout << std::endl;

   std::cout << "Dissolution modeling. Vickers hardness test." << std::endl;

   std::cout << "Minimum Vickers hardness:" << std::endl 
             << minimumVickersHardness << std::endl;

   std::cout << "Maximum Vickers hardness:" << std::endl 
             << maximumVickersHardness << std::endl;

   std::cout << "Reference time [s]:" << std::endl 
             << referenceTime << std::endl;

   std::cout << "Reference temperature [K]:" << std::endl 
             << referenceTemperature << std::endl;

   std::cout << "Number of samples:" << std::endl 
             << numberOfSamples << std::endl;

   std::cout << "Experimental data:" << std::endl;
 
   for(int i = 0; i < numberOfSamples; i++)
   {
      std::cout << temperatureData[i] << " " << timeData[i] << " " << vickersHardnessData[i] << std::endl;
   }
}


// void print(void) method

/// This method prints to the screen the actual effective activation energy of 
/// an aluminium alloy when training. 

void PrecipitateDissolutionModeling::print(void)
{
   double effectiveActivationEnergy = multilayerPerceptron->getSingleIndependentParameter(0);

   std::cout << "Effective activation energy: " << effectiveActivationEnergy << std::endl;
}




// void savePrecipitateDissolutionModel(char*) method

/// This method saves the volumetric fraction against the normalized time model
/// to a data file. It also saves the related data. 
///
/// @param filename Precipitate dissolution model filename.
///
/// @see saveVickersHardnessModel. 

void PrecipitateDissolutionModeling::savePrecipitateDissolutionModel(char* filename)
{
   std::fstream file; 

   file.open(filename, std::ios::out);

   // Control sentence

   if(!file.is_open())
   {
      std::cout << std::endl
                << "Cannot open precipitate dissolution model data file."
                << std::endl;
      
      exit(1);
   }
   else
   {
      std::cout << std::endl
                << "Saving precipitate dissolution model to data file..."
                << std::endl;
   }


   int numberOfInputs = multilayerPerceptron->getNumberOfInputs();
   Vector<int> numbersOfHiddenNeurons = multilayerPerceptron->getNumbersOfHiddenNeurons();
   int numberOfOutputs = multilayerPerceptron->getNumberOfOutputs();

   double effectiveActivationEnergy = multilayerPerceptron->getSingleIndependentParameter(0);

   Vector<double> normalizedTimeData = getNormalizedTimeData();
   Vector<double> volumetricFractionData = getVolumetricFractionData();

   Vector<double> input(1);
   Vector<double> output(1);

   Vector<double> volumetricFractionModel(numberOfSamples);

   for(int i = 0; i < numberOfSamples; i++)
   {
      input[0] = normalizedTimeData[i];
      output = multilayerPerceptron->calculateOutput(input);

      volumetricFractionModel[i] = output[0];
   }

   // Sort data

   Vector<double> sortedNormalizedTimeData = normalizedTimeData;

   std::sort(sortedNormalizedTimeData.begin(), sortedNormalizedTimeData.end(), std::less<double>());

   Vector<double> sortedVolumetricFractionData(numberOfSamples);
   Vector<double> sortedVolumetricFractionModel(numberOfSamples);

   for(int i = 0; i < numberOfSamples; i++)
   {
      for(int j = 0; j < numberOfSamples; j++)
      {
         if(normalizedTimeData[j] == sortedNormalizedTimeData[i])
         {
            sortedVolumetricFractionData[i] = volumetricFractionData[j];
            sortedVolumetricFractionModel[i] = volumetricFractionModel[j];
         }
      }
   }

   file << "% Flood Neural Network. Precipitate Dissolution Model Data File." << std::endl;

   file << "% Network architecture: " << numberOfInputs << ":" << numbersOfHiddenNeurons << ":" << numberOfOutputs << std::endl;
   
   file << "% Effective activation energy: " << effectiveActivationEnergy << std::endl;

   file << std::endl
       << "% Column data: " << std::endl
       << "% 1. Experimental X" << std::endl
       << "% 2. Experimental Y" << std::endl;

   
   for(int i = 0; i < numberOfSamples; i++)
   {
     file << normalizedTimeData[i] << " " << volumetricFractionData[i] << std::endl;
   }

   file << std::endl
       << "% Column data: " << std::endl
       << "% 1. Normalized time (log(t/t*))" << std::endl
       << "% 2. Volumetric fraction (1-f/f0)" << std::endl;

   for(int i = 0; i < numberOfSamples; i++)
   {
     file << sortedNormalizedTimeData[i] << " " << sortedVolumetricFractionData[i] << std::endl;
   }

   file.close();
}


// void saveVickersHardnessModel(char*) method

/// This method saves the Vickers hardness against the normalized time model
/// to a data file. It also saves the related data. 
///
/// @param filename Vickers hardness model filename.
///
/// @see saveVickersHardnessModel. 

void PrecipitateDissolutionModeling::saveVickersHardnessModel(char* filename)
{
   std::fstream file; 

   file.open(filename, std::ios::out);

   // Control sentence

   if(!file.is_open())
   {
      std::cout << std::endl
                << "Cannot open Vickers hardness model data file."
                << std::endl;
      
      exit(1);
   }
   else
   {
      std::cout << std::endl
                << "Saving Vickers hardness model to data file..."
                << std::endl;
   }

   int numberOfInputs = multilayerPerceptron->getNumberOfInputs();
   Vector<int> numbersOfHiddenNeurons = multilayerPerceptron->getNumbersOfHiddenNeurons();
   int numberOfOutputs = multilayerPerceptron->getNumberOfOutputs();

   double effectiveActivationEnergy = multilayerPerceptron->getSingleIndependentParameter(0);

   Vector<double> normalizedTimeData = getNormalizedTimeData();
   Vector<double> volumetricFractionData = getVolumetricFractionData();

   Vector<double> input(1);
   Vector<double> output(1);

   Vector<double> volumetricFractionModel(numberOfSamples);

   for(int i = 0; i < numberOfSamples; i++)
   {
      input[0] = normalizedTimeData[i];
      output = multilayerPerceptron->calculateOutput(input);

      volumetricFractionModel[i] = output[0];
   }

   // Sort data

   Vector<double> sortedNormalizedTimeData = normalizedTimeData;

   std::sort(sortedNormalizedTimeData.begin(), sortedNormalizedTimeData.end(), std::less<double>());

   Vector<double> sortedVolumetricFractionData(numberOfSamples);
   Vector<double> sortedVolumetricFractionModel(numberOfSamples);

   for(int i = 0; i < numberOfSamples; i++)
   {
      for(int j = 0; j < numberOfSamples; j++)
      {
         if(normalizedTimeData[j] == sortedNormalizedTimeData[i])
         {
            sortedVolumetricFractionData[i] = volumetricFractionData[j];
            sortedVolumetricFractionModel[i] = volumetricFractionModel[j];
         }
      }
   }

   file << "% Flood Neural Network. Precipitate Dissolution Model Data File." << std::endl;

   file << "% Network architecture: " << numberOfInputs << ":" << numbersOfHiddenNeurons << ":" << numberOfOutputs << std::endl;
   
   file << "% Effective activation energy: " << effectiveActivationEnergy << std::endl;

   file << std::endl
        << "% Column data: " << std::endl
        << "% 1. Normalized time data (log(t/t*))" << std::endl
        << "% 2. Vickers hardness data" << std::endl;

   
   for(int i = 0; i < numberOfSamples; i++)
   {
     file << normalizedTimeData[i] << " " << vickersHardnessData[i] << std::endl;
   }

   file << std::endl
        << "% Column data: " << std::endl
        << "% 1. Sorted normalized time data (log(t/t*))" << std::endl
        << "% 2. Sorted Vickers hardness model" << std::endl;

   for(int i = 0; i < numberOfSamples; i++)
   {
     file << sortedNormalizedTimeData[i] << " " << getVickersHardness(sortedVolumetricFractionModel[i]) << std::endl;
   }

   file.close();
}


// void saveInverseVickersHardnessTest(char*) method

/// This method saves to a file a Vickers hardness test generated by the 
/// precipitate dissolution model in order to compare it to the experimental data. 
///
/// @param filename Inverse Vickers hardness test filename.

void PrecipitateDissolutionModeling::saveInverseVickersHardnessTest(char* filename)
{
   std::fstream file; 

   file.open(filename, std::ios::out);

   // Control sentence

   if(!file.is_open())
   {
      std::cout << std::endl
                << "Cannot open reverse Vickers hardness test data file."
                << std::endl;
      
      exit(1);
   }
   else
   {
      std::cout << std::endl
                << "Saving reverse Vickers hardness test to data file..."
                << std::endl;
   }

   int numberOfTemperatures = 7;

   Vector<double> temperature(numberOfTemperatures);

   temperature[0] = 200.0 + 273.16;
   temperature[1] = 250.0 + 273.16;
   temperature[2] = 300.0 + 273.16;
   temperature[3] = 350.0 + 273.16;
   temperature[4] = 400.0 + 273.16;
   temperature[5] = 450.0 + 273.16;
   temperature[6] = 500.0 + 273.16;

   int numberOfPoints = 11;

   Vector<double> time(numberOfPoints);

   time[0] = 1.0;
   time[1] = 5.0;
   time[2] = 10.0;
   time[3] = 50.0;
   time[4] = 100.0;
   time[5] = 500.0;
   time[6] = 1000.0;
   time[7] = 5000.0;
   time[8] = 10000.0;
   time[9] = 50000.0;
   time[10] = 100000.0;

   // Calculate and write inverse Vickers hardness test

   for(int i  = 0; i < numberOfPoints; i++)
   {
      file << time[i] << " ";

      for(int j = 0; j < numberOfTemperatures; j++)
      {
         double fullDissolutionTime = getFullDissolutionTime(temperature[j]);
          
         Vector<double> input(1, log10(time[i]/fullDissolutionTime));
         Vector<double> output = multilayerPerceptron->calculateOutput(input);
             
         double volumetricFraction = output[0];
         double vickersHardness = getVickersHardness(volumetricFraction);

         file << vickersHardness << " ";
      }
 
      file << std::endl;
   }

   // Close file

   file.close();
}

}


// Flood: An Open Source Neural Networks C++ Library.
// Copyright (C) 2005-2008 Roberto Lopez 
//
// This library is free software; you can redistribute it and/or
// modify it under the s of the GNU Lesser General Public
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
