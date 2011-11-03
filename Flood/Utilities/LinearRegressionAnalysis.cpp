/****************************************************************************************************************/
/*                                                                                                              */
/*   Flood: An Open Source Neural Networks C++ Library                                                          */
/*   www.cimne.com/flood                                                                                        */
/*                                                                                                              */
/*   L I N E A R   R E G R E S S I O N   A N A L Y S I S   C L A S S                                            */
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

#include "LinearRegressionAnalysis.h"

namespace Flood
{

// GENERAL CONSTRUCTOR

/// General constructor. It creates a linear regression analysis object 
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

LinearRegressionAnalysis::LinearRegressionAnalysis
(MultilayerPerceptron* newMultilayerPerceptron, InputTargetDataSet* newInputTargetDataSet)
{
   multilayerPerceptron = newMultilayerPerceptron;
   
   inputTargetDataSet = newInputTargetDataSet;

   display = true;
}


// DEFAULT CONSTRUCTOR

/// Default constructor. It creates a linear regression analysis object 
/// not associated to any multilayer perceptron neither any input-target data set 
///objects.
/// This constructor also initializes the rest of class members to their default values:
///
/// <ul>
/// <li> Display: True.
/// </ul> 

LinearRegressionAnalysis::LinearRegressionAnalysis(void)
{
   multilayerPerceptron = NULL;

   inputTargetDataSet = NULL;

   display = true;
}


// DESTRUCTOR

/// Destructor. 

LinearRegressionAnalysis::~LinearRegressionAnalysis()
{

}


// METHODS

// MultilayerPerceptron* getMultilayerPerceptron(void) method

/// This method returns a pointer to the multilayer perceptron object which is 
/// to be validated.

MultilayerPerceptron* LinearRegressionAnalysis::getMultilayerPerceptron(void)
{
   return(multilayerPerceptron);   
}


// InputTargetDataSet* getInputTargetDataSet(double) method

/// This method returns a pointer to the input-target data set object used for 
/// validating the performance of a trained multilayer perceptron.

InputTargetDataSet* LinearRegressionAnalysis::getInputTargetDataSet(void)
{
   return(inputTargetDataSet);   
}


// bool getDisplay(void) method

/// This method returns true if messages from this class can be displayed on the screen,
/// or false if messages from this class can't be displayed on the screen.

bool LinearRegressionAnalysis::getDisplay(void)
{
   return(display);     
}


// void setMultilayerPerceptron(MultilayerPerceptron*) method

/// This method sets a new multilayer perceptron which is to be validated.
///
/// @param newMultilayerPerceptron Pointer to a multilayer perceptron object.

void LinearRegressionAnalysis
::setMultilayerPerceptron(MultilayerPerceptron* newMultilayerPerceptron)
{
   multilayerPerceptron = newMultilayerPerceptron;   
}


// void setInputTargetDataSet(InputTargetDataSet*) method

/// This method sets a new input-target data set to be used for validating the
/// quality of a trained multilayer perceptron.
///
/// @param newInputTargetDataSet Pointer to an input-target data set object.

void LinearRegressionAnalysis
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

void LinearRegressionAnalysis::setDisplay(bool newDisplay)
{
   display = newDisplay;
}


// Matrix<double> getOutputData(void)

/// This method retuns a Matrix containing the neural network outputs 
/// for the input data set.

Matrix<double> LinearRegressionAnalysis::getOutputData(void)
{
   // Control sentence (if debug)

   #ifndef NDEBUG 

   if(inputTargetDataSet == NULL)
   {
      std::cerr << std::endl
                << "Flood Error: LinearRegressionAnalysis class." << std::endl 
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
                << "Flood Error: LinearRegressionAnalysis class." << std::endl 
                << "Matrix<double> getOutputData(void) method." << std::endl
                << "Number of inputs in multilayer perceptron " << std::endl
                << "must be equal to number of input variables in training data set." 
                << std::endl << std::endl;

      exit(1);
   }
   else if(numberOfOutputs != numberOfTargetVariables)
   {
      std::cerr << std::endl
                << "Flood Error: LinearRegressionAnalysis class." << std::endl
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


// Matrix<double> calculateRegressionParameters(void)

/// This method performs a linear regression analysis between the neural network outputs and the corresponding 
/// data set targets and returns all the provided parameters in a single matrix. 
/// The number of rows in the matrix is equal to the number of output variables. 
/// The number of columns in the matrix is equal to the number of regression parameters (3). In this way, each 
/// row contains the regression parameters a, b and R of an output variable.

Matrix<double> LinearRegressionAnalysis::calculateRegressionParameters(void)
{
   // Control sentence (if debug)

   #ifndef NDEBUG 

   if(multilayerPerceptron == NULL)
   {
      std::cerr << std::endl
                << "Flood Error: LinearRegressionAnalysis class." << std::endl 
                << "Matrix<double> calculateRegressionParameters(void)." << std::endl
                << "Multilayer perceptron object cannot be null." << std::endl
            << std::endl;

      exit(1);   
   }
   else if(inputTargetDataSet == NULL)
   {
      std::cerr << std::endl
                << "Flood Error: LinearRegressionAnalysis class." << std::endl 
                << "Matrix<double> calculateRegressionParameters(void)." << std::endl
                << "Input-target data set object cannot be null." << std::endl
            << std::endl;

      exit(1);   
   }

   #endif

   // Multilayer perceptron stuff
                            
   int numberOfOutputs = multilayerPerceptron->getNumberOfOutputs();

   // Input target data set stuff

   Matrix<double>& targetData = inputTargetDataSet->getTargetData();  
   int numberOfSamples = inputTargetDataSet->getNumberOfSamples();

   // Get output data

   Matrix<double> outputData = getOutputData();            

   // Get regression parameters

   Vector<double> sumOfTargets(numberOfOutputs, 0.0);
   Vector<double> sumOfOutputs(numberOfOutputs, 0.0);
   Vector<double> sumOfSquaredTargets(numberOfOutputs, 0.0);
   Vector<double> sumOfSquaredOutputs(numberOfOutputs, 0.0);
   Vector<double> sumOfTargetsOutputsProduct(numberOfOutputs, 0.0);

   for(int i = 0; i < numberOfSamples; i++)
   {
      for(int j = 0; j < numberOfOutputs; j++)
      {
         sumOfTargets[j] += targetData[i][j];
         sumOfOutputs[j] += outputData[i][j];
         sumOfSquaredTargets[j] += pow(targetData[i][j],2);
         sumOfSquaredOutputs[j] += pow(targetData[i][j],2);
         sumOfTargetsOutputsProduct[j] += targetData[i][j]*outputData[i][j];
      }
   }
   
   Vector<double> a(numberOfOutputs, 0.0);
   Vector<double> b(numberOfOutputs, 0.0);
   Vector<double> R(numberOfOutputs, 0.0);

   double numerator = 0.0;
   double denominator = 0.0;
   
   // Targets = X
   // Outputs = Y
   
   for(int i = 0; i < numberOfOutputs; i++)
   {
      // a
      
      numerator = sumOfOutputs[i]*sumOfSquaredTargets[i] 
                - sumOfTargets[i]*sumOfTargetsOutputsProduct[i];
      
      denominator = numberOfSamples*sumOfSquaredTargets[i] 
                  - pow(sumOfTargets[i],2);
      
      a[i] = numerator/denominator;
      
      // b
      
      numerator = numberOfSamples*sumOfTargetsOutputsProduct[i]
                - sumOfTargets[i]*sumOfOutputs[i];
                
      denominator = numberOfSamples*sumOfSquaredTargets[i] 
                  - pow(sumOfTargets[i],2);
       
      b[i] = numerator/denominator; 
      
      // Correlation coefficient
      
      numerator = numberOfSamples*sumOfTargetsOutputsProduct[i]
                - sumOfTargets[i]*sumOfOutputs[i];    
      
      denominator = 
      sqrt((numberOfSamples*sumOfSquaredTargets[i] - pow(sumOfOutputs[i],2))
          *(numberOfSamples*sumOfSquaredOutputs[i] - pow(sumOfOutputs[i],2))); 
          
      R[i] = numerator/denominator;
   }   

   Matrix<double> regressionParameters(numberOfOutputs, 3, 0.0);   

   for(int i = 0; i < numberOfOutputs; i++)
   {
      regressionParameters[i][0] = a[i];
      regressionParameters[i][1] = b[i];
      regressionParameters[i][2] = R[i];         
   }
 
   return(regressionParameters);          
}


// void printResults() method

/// This method performs a linear regression analysis between the neural network outputs and the corresponding 
/// data set targets and then prints all results in the screen:
///
/// <ul>
/// <li> Regression parameters (a, b and R) for each output variable.
/// <li> Target and output data for each output variable.
/// </ul> 
///
/// @see getSaveResults(char*).

void LinearRegressionAnalysis::printResults(void)
{         
   std::cout << std::endl
             << "Flood Neural Network. Linear Regression Analysis Results."
             << std::endl;

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

   Matrix<double> regressionParameters = calculateRegressionParameters();

   // Print results
   
   for(int i = 0; i < numberOfOutputs; i++)
   {
      // Print regression parameters     
           
      std::cout << std::endl 
                << nameOfOutputVariables[i] << std::endl
                << "Regression parameters:" 
                << std::endl;
      
      std::cout << "a = " << regressionParameters[i][0] << ";" << std::endl;
      std::cout << "b = " << regressionParameters[i][1] << ";" << std::endl;
      std::cout << "R = " << regressionParameters[i][2] << ";" << std::endl;     
      
      // Print target and output data
      
      std::cout << "Target and output data:" << std::endl;
    
      for(int j = 0; j < numberOfSamples; j++)
      {
         std::cout << targetData[j][i] << " " << outputData[j][i] << std::endl;      
      }                          
   }
}


// void saveResults(char*) method

/// This method performs a linear regression analysis between the neural network outputs and the corresponding 
/// data set targets and then saves all results to a data file:
///
/// <ul>
/// <li> Regression parameters (a, b and R) for each output variable.
/// <li> Target and output data for each output variable.
/// </ul> 
///
/// @param filename Filename.
///
/// @see printResults(char*).

void LinearRegressionAnalysis::saveResults(char* filename)
{ 
   std::fstream file; 

   file.open(filename, std::ios::out);

   if(!file.is_open())
   {
      std::cerr << std::endl
                << "Flood Error: LinearRegressionAnalysis class." << std::endl
                << "void saveResults(char*) method." << std::endl
                << "Cannot open regression analysis results data file."
                << std::endl << std::endl;

      exit(1);
   }
   else
   {
      if(display)
      {
         std::cout << std::endl
                   << "Saving linear regression analysis results to data file..."
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

   Matrix<double> regressionParameters = calculateRegressionParameters();

   // Save results to file

   // Write file header

   file << "% Flood Neural Network. Linear Regression Analysis Results." << std::endl;

   for(int i = 0; i < numberOfOutputs; i++)
   {     

      // Write regression parameters     
           
      file << "% " << std::endl
           << "% " << nameOfOutputVariables[i] << std::endl;
           
      file << "% Linear regression parameters (y = a + bx):" << std::endl
           << "% a = " << regressionParameters[i][0] << ";" << std::endl
           << "% b = " << regressionParameters[i][1] << ";" << std::endl
           << "% R = " << regressionParameters[i][2] << ";" << std::endl;     
      
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
