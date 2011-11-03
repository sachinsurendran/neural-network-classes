/****************************************************************************************************************/
/*                                                                                                              */
/*   Flood: An Open Source Neural Networks C++ Library                                                          */
/*   www.cimne.com/flood                                                                                        */
/*                                                                                                              */
/*   C O R R E C T   P R E D I C T I O N S   A N A L Y S I S   C L A S S                                        */
/*                                                                                                              */ 
/*   Roberto Lopez                                                                                              */ 
/*   International Center for Numerical Methods in Engineering (CIMNE)                                          */
/*   Technical University of Catalonia (UPC)                                                                    */
/*   Barcelona, Spain                                                                                           */
/*   E-mail: rlopez@cimne.upc.edu                                                                               */ 
/*                                                                                                              */
/****************************************************************************************************************/

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <math.h>

#include "CorrectPredictionsAnalysis.h"

namespace Flood
{

// GENERAL CONSTRUCTOR

/// General constructor. It creates a correct predictions analysis object 
/// associated to a multilayer perceptron and an input-target data set objects.
/// That multilayer perceptron is going to be validated on that input-target
/// data set.
/// This constructor also initializes the rest of class members to their default values:
///
/// <ul>
/// <li> Display: True.
/// </ul> 
///
/// @param newMultilayerPerceptron Pointer to a multilayer perceptron object.
/// @param newInputTargetDataSet Pointer to an input-target data set object.

CorrectPredictionsAnalysis::CorrectPredictionsAnalysis
(MultilayerPerceptron* newMultilayerPerceptron, InputTargetDataSet* newInputTargetDataSet)
{
   multilayerPerceptron = newMultilayerPerceptron;
   
   inputTargetDataSet = newInputTargetDataSet;

   display = true;
}


// DEFAULT CONSTRUCTOR

/// Default constructor. It creates a correct predictions analysis object 
/// not associated to any multilayer perceptron neither any input-target data set objects.
/// This constructor also initializes the rest of class members to their default values:
///
/// <ul>
/// <li> Display: True.
/// </ul> 

CorrectPredictionsAnalysis::CorrectPredictionsAnalysis(void)
{
   multilayerPerceptron = NULL;

   inputTargetDataSet = NULL;

   display = true;
}


// DESTRUCTOR

/// Destructor. 

CorrectPredictionsAnalysis::~CorrectPredictionsAnalysis()
{

}


// METHODS

// MultilayerPerceptron* getMultilayerPerceptron(void) method

/// This method returns a pointer to the multilayer perceptron object which is 
/// to be validated.

MultilayerPerceptron* CorrectPredictionsAnalysis::getMultilayerPerceptron(void)
{
   return(multilayerPerceptron);   
}


// InputTargetDataSet* getInputTargetDataSet(double) method

/// This method returns a pointer to the input-target data set object used for 
/// validating the performance of a trained multilayer perceptron.

InputTargetDataSet* CorrectPredictionsAnalysis::getInputTargetDataSet(void)
{
   return(inputTargetDataSet);   
}


// bool getDisplay(void) method

/// This method returns true if messages from this class can be displayed on the screen,
/// or false if messages from this class can't be displayed on the screen.

bool CorrectPredictionsAnalysis::getDisplay(void)
{
   return(display);     
}


// void setMultilayerPerceptron(MultilayerPerceptron*) method

/// This method sets a new multilayer perceptron which is to be validated.
///
/// @param newMultilayerPerceptron Pointer to a multilayer perceptron object.

void CorrectPredictionsAnalysis
::setMultilayerPerceptron(MultilayerPerceptron* newMultilayerPerceptron)
{
   multilayerPerceptron = newMultilayerPerceptron;   
}


// void setInputTargetDataSet(InputTargetDataSet*) method

/// This method sets a new input-target data set to be used for validating the
/// quality of a trained multilayer perceptron.
///
/// @param newInputTargetDataSet Pointer to an input-target data set object.

void CorrectPredictionsAnalysis
::setInputTargetDataSet(InputTargetDataSet* newInputTargetDataSet)
{
   inputTargetDataSet = newInputTargetDataSet;   
}


// void setDisplay(bool) method

/// This method sets a new display value. 
/// If it is set to true messages from this class are to be displayed on the screen;
/// if it is set to false messages from this class are not to be displayed on the screen.
///
/// @param newDisplay Display value.

void CorrectPredictionsAnalysis::setDisplay(bool newDisplay)
{
   display = newDisplay;
}


// Matrix<double> getOutputData(void)

/// This method retuns a Matrix containing the neural network outputs 
/// for the input data set.

Matrix<double> CorrectPredictionsAnalysis::getOutputData(void)
{
   // Control sentence (if debug)

   #ifndef NDEBUG 

   if(inputTargetDataSet == NULL)
   {
      std::cerr << std::endl
                << "Flood Error: CorrectPredictionsAnalysis class." << std::endl 
                << "Matrix<double> getOutputData(void) method." << std::endl
                << "Input-target data set object cannot be null." << std::endl
                << std::endl;

      exit(1);   
   }

   #endif
             
   int numberOfInputs = multilayerPerceptron->getNumberOfInputs();
   int numberOfOutputs = multilayerPerceptron->getNumberOfOutputs();

   #ifndef NDEBUG 

   int numberOfInputVariables = inputTargetDataSet->getNumberOfInputVariables();
   int numberOfTargetVariables = inputTargetDataSet->getNumberOfTargetVariables();

   // Control sentence
   
   if(numberOfInputs != numberOfInputVariables)
   {
      std::cerr << std::endl
                << "Flood Error: CorrectPredictionsAnalysis class." << std::endl 
                << "Matrix<double> getOutputData(void) method." << std::endl
                << "Number of inputs in multilayer perceptron " << std::endl
                << "must be equal to number of input variables in training data set." 
                << std::endl << std::endl;

      exit(1);
   }
   else if(numberOfOutputs != numberOfTargetVariables)
   {
      std::cerr << std::endl
                << "Flood Error: CorrectPredictionsAnalysis class." << std::endl
                << "Matrix<double> getOutputData(void) method." << std::endl
                << "Number of outputs in multilayer perceptron " << std::endl
                << "must be equal to number of target variables in training data set." 
                << std::endl << std::endl;

      exit(1);
   }

   #endif

   // Get output data

   int numberOfSamples = inputTargetDataSet->getNumberOfSamples();

   Matrix<double> outputData(numberOfSamples, numberOfOutputs, 0.0);

   Matrix<double>& inputData = inputTargetDataSet->getInputData(); 

   Vector<double> input(numberOfInputs, 0.0);
   Vector<double> output(numberOfOutputs, 0.0);

   for(int i = 0; i < numberOfSamples; i++)
   {
      input = inputData.getRow(i);        
     
      output = multilayerPerceptron->calculateOutput(input);

      outputData.setRow(i, output);
   }               
   
   return(outputData);
}


// Vector<double> calculateCorrectPredictionsRatio(void)

/// This method performs a correct predictions analysis between the neural network outputs and the corresponding 
/// data set targets and returns all the provided parameters in a single vector. 

double CorrectPredictionsAnalysis::calculateCorrectPredictionsRatio(void)
{
   // Control sentence (if debug)

   #ifndef NDEBUG 

   if(multilayerPerceptron == NULL)
   {
      std::cerr << std::endl
                << "Flood Error: CorrectPredictionsAnalysis class." << std::endl 
                << "Vector<double> calculateCorrectPredictionsRatio(void)." << std::endl
                << "Multilayer perceptron object cannot be null." << std::endl
                << std::endl;

      exit(1);   
   }
   else if(inputTargetDataSet == NULL)
   {
      std::cerr << std::endl
                << "Flood Error: CorrectPredictionsAnalysis class." << std::endl 
                << "Vector<double> calculateCorrectPredictionsRatio(void)." << std::endl
                << "Input-target data set object cannot be null." << std::endl
                << std::endl;

      exit(1);   
   }

   #endif

   // Multilayer perceptron stuff
                            
   int numberOfOutputs = multilayerPerceptron->getNumberOfOutputs();

   // Input target data set stuff

   Matrix<double>& inputData = inputTargetDataSet->getInputData();  
   Matrix<double>& targetData = inputTargetDataSet->getTargetData();  

   int numberOfTargetVariables = inputTargetDataSet->getNumberOfTargetVariables();
   int numberOfSamples = inputTargetDataSet->getNumberOfSamples();

   // Get output data

   Matrix<double> outputData = getOutputData();            

   // Get percentage of correct predictions in validation data set

   Vector<bool> observedPattern(numberOfTargetVariables);
   Vector<bool> predictedPattern(numberOfOutputs);

   int count = 0;

   for(int i = 0; i < numberOfSamples; i++)
   {
      // Obtain observed value
   
      for(int j = 0; j < numberOfOutputs; j++)
      {
         if(targetData[i][j] < 0.5)
         {
            observedPattern[j] = false;
         }
         else
         {
            observedPattern[j] = true; 
         }
      }

      // Obtain predicted value

      Vector<double> input = inputData.getRow(i);

      Vector<double> output = multilayerPerceptron->calculateOutput(input);

      for(int j = 0; j < numberOfOutputs; j++)
      {  
         if(output[j] < 0.5)
         {
            predictedPattern[j] = false;
         }
         else
         {
            predictedPattern[j] = true; 
         }
      }

      bool correctPrediction = true;

      for(int j = 0; j < numberOfTargetVariables; j++)
      { 
         if(observedPattern[j] != predictedPattern[j])
         {
            correctPrediction = false;
         }
      }

      if(correctPrediction == true)
	  {
	     count++;
	  }
   }
   double correctPredictionsRatio = (double)count/(double)numberOfSamples;

   return(correctPredictionsRatio);          
}


// void printResults() method

/// This method performs a correct predictions analysis between the neural network outputs and the corresponding 
/// data set targets and then prints all results in the screen:
///
/// <ul>
/// <li> Correct predictions ratio.
/// <li> Target and output data for each output variable.
/// </ul> 
///
/// @see getSaveResults(char*).

void CorrectPredictionsAnalysis::printResults(void)
{         
   std::cout << std::endl
             << "Flood Neural Network. Correct Prediction Analysis Results."
             << std::endl;

   // Multilayer perceptron stuff

   int numberOfOutputs = multilayerPerceptron->getNumberOfOutputs();

   Vector<std::string> nameOfOutputVariables = multilayerPerceptron->getNameOfOutputVariables();

   // Input target data set stuff

   int numberOfSamples = inputTargetDataSet->getNumberOfSamples();

   Matrix<double>& targetData = inputTargetDataSet->getTargetData(); 

   // Get output data

   Matrix<double> outputData = getOutputData();

   // Get regression parameters

   double correctPredictionsRatio = calculateCorrectPredictionsRatio();

   // Print results

   std::cout << std::endl 
             << "Correct predictions ratio: " << correctPredictionsRatio << std::endl;


   for(int i = 0; i < numberOfOutputs; i++)
   {    
      // Print target and output data
      
      std::cout << std::endl
	     	    << "Target and output data:" << std::endl;
    
      for(int j = 0; j < numberOfSamples; j++)
      {
         std::cout << targetData[j][i] << " " << outputData[j][i] << std::endl;      
      }                          
   }
}


// void saveResults(char*) method

/// This method performs a correct predictions analysis between the neural network outputs and the corresponding 
/// data set targets and then saves all results to a data file:
///
/// <ul>
/// <li> Correct predictions value.
/// </ul> 
///
/// @param filename Filename.
///
/// @see printResults(char*).

void CorrectPredictionsAnalysis::saveResults(char* filename)
{ 
   std::fstream file; 

   file.open(filename, std::ios::out);

   if(!file.is_open())
   {
      std::cerr << std::endl
                << "Flood Error: CorrectPredictionsAnalysis class." << std::endl
                << "void saveResults(char*) method." << std::endl
                << "Cannot open correct predictions analysis results data file."
                << std::endl << std::endl;

      exit(1);
   }
   else
   {
      if(display)
      {
         std::cout << std::endl
                   << "Saving results of correct predictions analysis to data file..."
                   << std::endl;
      }
   }
     
   // Multilayer perceptron stuff

   int numberOfOutputs = multilayerPerceptron->getNumberOfOutputs();

   Vector<std::string> nameOfOutputVariables 
   = multilayerPerceptron->getNameOfOutputVariables();

   // Input target data set stuff

   int numberOfSamples = inputTargetDataSet->getNumberOfSamples();

   Matrix<double>& targetData = inputTargetDataSet->getTargetData(); 

   // Get output data

   Matrix<double> outputData = getOutputData();

   // Get regression parameters

   double correctPredictionsRatio = calculateCorrectPredictionsRatio();

   // Save results to file

   // Write file header

   file << "% Flood Neural Network. Correct Predictions Analysis Results." << std::endl;

   file << "% Correct predictions ratio: " << correctPredictionsRatio << std::endl;


   for(int i = 0; i < numberOfOutputs; i++)
   {     

      // Write regression parameters     
           
      file << "% " << std::endl
           << "% " << nameOfOutputVariables[i] << std::endl;          
      
      // Write target and output data
      
      file << "% Target and output data:" << std::endl;
    
      for(int j = 0; j < numberOfSamples; j++)
      {
         file << targetData[j][i] << " " << outputData[j][i] << std::endl;      
      }                          
   }

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
