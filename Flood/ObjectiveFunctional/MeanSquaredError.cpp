/****************************************************************************************************************/
/*                                                                                                              */
/*   Flood: An Open Source Neural Networks C++ Library                                                          */
/*   www.cimne.com/flood                                                                                        */
/*                                                                                                              */
/*   M E A N   S Q U A R E D   E R R O R   C L A S S                                                            */
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
#include <limits>
#include <math.h>

#include "MeanSquaredError.h"

namespace Flood
{

// GENERAL CONSTRUCTOR

/// General constructor. It creates a mean squared error objective functional object associated to a multilayer 
/// perceptron and measured on an input-target data set.
/// It also initializes all the rest of class members to their default values:
///
/// <ul>
/// <li> Number of evaluations = 0. 
/// <li> Epsilon: 1.0e-6.
/// <li> Display = true. 
/// </ul> 
///
/// @param newMultilayerPerceptron Pointer to a multilayer perceptron object.
/// @param newInputTargetDataSet Pointer to an input-target data set object.

MeanSquaredError::MeanSquaredError(MultilayerPerceptron* newMultilayerPerceptron, 
InputTargetDataSet* newInputTargetDataSet)
: ObjectiveFunctional(newMultilayerPerceptron)
{
   inputTargetDataSet = newInputTargetDataSet;
}


// DEFAULT CONSTRUCTOR

/// Default constructor. It creates a mean squared error objective functional object not associated to any 
/// multilayer perceptron and not measured on any input-target data set.
/// It also initializes all the rest of class members to their default values:
///
/// <ul>
/// <li> Number of evaluations = 0. 
/// <li> Epsilon: 1.0e-6.
/// <li> Display = true. 
/// </ul> 

MeanSquaredError::MeanSquaredError(void) : ObjectiveFunctional()
{
   inputTargetDataSet = NULL;
}


// DESTRUCTOR

/// Destructor.

MeanSquaredError::~MeanSquaredError(void)
{

}


// METHODS

// InputTargetDataSet* getInputTargetDataSet(void) method

/// This method returns a pointer to the input-target object on which the objectiveal is measured.

InputTargetDataSet* MeanSquaredError::getInputTargetDataSet(void)
{
   return(inputTargetDataSet);
}


// void setInputTargetDataSet(InputTargetDataSet*) method

/// This method sets a pointer to an input-data set object on which the objective functional is to be measured.
///
/// @param newInputTargetDataSet Pointer to an input-target data set object.

void MeanSquaredError::setInputTargetDataSet(InputTargetDataSet* newInputTargetDataSet)
{
   inputTargetDataSet = newInputTargetDataSet;
}


// double calculateEvaluation(void) method

/// This method returns the evaluation value of a multilayer perceptron according to the mean squared error on 
/// an input-target data set.
///
/// @see calculateGradient(void).

double MeanSquaredError::calculateEvaluation(void)
{
    // Control sentence 

   if(multilayerPerceptron == NULL)
   {
      std::cerr << std::endl
                << "Flood Error: MeanSquaredError class." << std::endl
                << "double calculateEvaluation(void) method." << std::endl
                << "Pointer to multilayer perceptron object cannot be NULL." << std::endl
                << std::endl;

        exit(1);
   }
   else if(inputTargetDataSet == NULL)
   {
      std::cerr << std::endl
                << "Flood Error: MeanSquaredError class." << std::endl
                << "double calculateEvaluation(void) method." << std::endl
                << "Pointer to input-target data set object cannot be NULL." << std::endl
                << std::endl;

        exit(1);
   }

   // Calculate evaluation

   double meanSquaredError = 0.0;       

   // Multilayer perceptron

   int numberOfInputs = multilayerPerceptron->getNumberOfInputs();
   int numberOfOutputs = multilayerPerceptron->getNumberOfOutputs();

   // Input-target data set

   int numberOfSamples = inputTargetDataSet->getNumberOfSamples();
   int numberOfInputVariables = inputTargetDataSet->getNumberOfInputVariables();
   int numberOfTargetVariables = inputTargetDataSet->getNumberOfTargetVariables();

   Matrix<double>& inputData = inputTargetDataSet->getInputData();
   Matrix<double>& targetData = inputTargetDataSet->getTargetData();

   if(numberOfInputs != numberOfInputVariables || numberOfOutputs != numberOfTargetVariables)
   {
      std::cout << std::endl
                << "Flood Error MeanSquaredError class." << std::endl 
                << "double calculateEvaluation(void) method." << std::endl
                << "Number of inputs and outputs in multilayer perceptron must be equal to number of input and "
        	<< "output variables in input-target data set." << std::endl 
                << std::endl;

      exit(1);
   }

   // Increment number of evaluations

   numberOfEvaluations++;

   double sumSquaredError = 0.0;

   Vector<double> input(numberOfInputs, 0.0);
   Vector<double> output(numberOfOutputs, 0.0);
   Vector<double> target(numberOfOutputs, 0.0);

   for(int sample = 0; sample < numberOfSamples; sample++)
   {
      // Input vector

      input = inputData.getRow(sample);

      // Output vector

      output = multilayerPerceptron->calculateOutput(input);

      // Target vector

     target = targetData.getRow(sample);

      // Sum of squares error

      sumSquaredError += (output-target).dot(output-target); 
   }

   meanSquaredError = sumSquaredError/(double)numberOfSamples;

   return(meanSquaredError);
}


// Vector<double> calculateGradient(void)

/// This method returns the the mean squared error function gradient of a multilayer perceptron on an 
/// input-target data set. It uses the error back-propagation method.
///
/// @see calculateEvaluation(void).

Vector<double> MeanSquaredError::calculateGradient(void)
{
   // Multilayer perceptron 

   int numberOfInputs = multilayerPerceptron->getNumberOfInputs();

   int numberOfHiddenLayers = multilayerPerceptron->getNumberOfHiddenLayers();
   Vector<int> numbersOfHiddenNeurons = multilayerPerceptron->getNumbersOfHiddenNeurons();

   int numberOfOutputs = multilayerPerceptron->getNumberOfOutputs();

   MultilayerPerceptron::PreAndPostProcessingMethod preAndPostProcessingMethod
   = multilayerPerceptron->getPreAndPostProcessingMethod();

   Vector<double> meanOfOutputVariables = multilayerPerceptron->getMeanOfOutputVariables();
   Vector<double> standardDeviationOfOutputVariables = multilayerPerceptron->getStandardDeviationOfOutputVariables();

   Vector<double> minimumOfOutputVariables = multilayerPerceptron->getMinimumOfOutputVariables();
   Vector<double> maximumOfOutputVariables = multilayerPerceptron->getMaximumOfOutputVariables();

   Vector<Perceptron>& outputLayer = multilayerPerceptron->getOutputLayer();
   Vector< Vector<Perceptron> >& hiddenLayers = multilayerPerceptron->getHiddenLayers();

   int numberOfNeuralParameters = multilayerPerceptron->getNumberOfNeuralParameters();
   Vector<double> synapticWeights;

   // Input-target data set

   int numberOfSamples = inputTargetDataSet->getNumberOfSamples();
   int numberOfTargetVariables = inputTargetDataSet->getNumberOfTargetVariables();

   Matrix<double>& inputData = inputTargetDataSet->getInputData();
   Matrix<double>& targetData = inputTargetDataSet->getTargetData();

   Vector<double> target(numberOfTargetVariables);

   // Input, outputs, gradient

   Vector<double> input(numberOfInputs);
   Vector<double> inputSignal(numberOfInputs);
   Vector< Vector<double> > netInputSignalToHiddenLayer(numberOfHiddenLayers);
   Vector< Vector<double> > outputSignalFromHiddenLayer(numberOfHiddenLayers);
   Vector< Vector<double> > outputSignalDerivativeFromHiddenLayer(numberOfHiddenLayers);
   Vector<double> netInputSignalToOutputLayer(numberOfOutputs);
   Vector<double> outputSignal(numberOfOutputs);
   Vector<double> output(numberOfOutputs);
   Vector<double> gradient(numberOfNeuralParameters, 0.0);

   // Output and hidden errors

   Vector<double> outputError(numberOfOutputs);
   Vector< Vector<double> > hiddenError(numberOfHiddenLayers);

   for(int h = 0; h < numberOfHiddenLayers; h++)
   {
      hiddenError[h].setSize(numbersOfHiddenNeurons[h]);
   }

   // Main loop   

   double sum;
   int index;

   for(int sample = 0; sample < numberOfSamples; sample++)
   {
      // Actual input

      input = inputData.getRow(sample);

      // Preprocess input to obtain input signal to the neural network 

      inputSignal = multilayerPerceptron->preprocessInput(input);

      // Calculate net input output and output signal derivative from all hidden layers

      netInputSignalToHiddenLayer[0] 
      = multilayerPerceptron->calculateNetInputSignalToHiddenLayer(0, inputSignal);

      outputSignalFromHiddenLayer[0] 
      = multilayerPerceptron->calculateOutputSignalFromHiddenLayer(0, netInputSignalToHiddenLayer[0]);

      outputSignalDerivativeFromHiddenLayer[0] 
      = multilayerPerceptron->calculateOutputSignalDerivativeFromHiddenLayer(0, netInputSignalToHiddenLayer[0]);
  
      for(int h = 1; h < numberOfHiddenLayers; h++)
      {
         netInputSignalToHiddenLayer[h] = multilayerPerceptron
         ->calculateNetInputSignalToHiddenLayer(h, outputSignalFromHiddenLayer[h-1]);
   
         outputSignalFromHiddenLayer[h] = multilayerPerceptron
         ->calculateOutputSignalFromHiddenLayer(h, netInputSignalToHiddenLayer[h]);

         outputSignalDerivativeFromHiddenLayer[h] = multilayerPerceptron
         ->calculateOutputSignalDerivativeFromHiddenLayer(h, netInputSignalToHiddenLayer[h]);         
      }

      // Get net input signal to output layer

      netInputSignalToOutputLayer = multilayerPerceptron
      ->calculateNetInputSignalToOutputLayer(outputSignalFromHiddenLayer[numberOfHiddenLayers-1]);

      // Get actual output signal from output layer

      Vector<double> outputSignal = multilayerPerceptron->calculateOutputSignal(netInputSignalToOutputLayer);

      // Get actual output signal derivative from output layer

      Vector<double> outputSignalDerivative
      = multilayerPerceptron->calculateOutputSignalDerivative(netInputSignalToOutputLayer);

      // Postprocess output signal from the neural network to get output

      Vector<double> output = multilayerPerceptron->postprocessOutputSignal(outputSignal);

      // Actual target

      Vector<double> target = targetData.getRow(sample);

      // Evaluate the error for all the output units
      
      switch(preAndPostProcessingMethod)
      {
         case MultilayerPerceptron::None:
         { 
            outputError = outputSignalDerivative*2.0*(output - target)/(double)numberOfSamples;
	 }
         break;

         case MultilayerPerceptron::MeanAndStandardDeviation:
         {
               outputError = outputSignalDerivative*2.0*(output-target)
                            *standardDeviationOfOutputVariables/(double)numberOfSamples;
         }			
         break;

         case MultilayerPerceptron::MinimumAndMaximum:
         {
               outputError = outputSignalDerivative*2.0*(output-target)*0.5
               *(maximumOfOutputVariables-minimumOfOutputVariables)/(double)numberOfSamples;
         }
         break;
      }           

      // Backpropagate the output errors to obtain the hidden errors

      // Last hidden layer

      for(int j = 0; j < numbersOfHiddenNeurons[numberOfHiddenLayers-1]; j++)
      {
         sum = 0.0;

         for(int k = 0; k < numberOfOutputs; k++)
         {
            synapticWeights = outputLayer[k].getSynapticWeights();
            sum += (synapticWeights[j])*outputError[k];
         }

         hiddenError[numberOfHiddenLayers-1][j] = outputSignalDerivativeFromHiddenLayer[numberOfHiddenLayers-1][j]*sum;
      }

      // Rest of hidden layers

      for(int h = numberOfHiddenLayers-2; h >= 0; h--) 
      {   
         for(int j = 0; j < numbersOfHiddenNeurons[h]; j++)
         {
            sum = 0.0;

            for(int k = 0; k < numbersOfHiddenNeurons[h+1]; k++)
            {
               synapticWeights = hiddenLayers[h+1][k].getSynapticWeights();
               sum += (synapticWeights[j])*hiddenError[h+1][k];
            }		 
	  
            hiddenError[h][j] = outputSignalDerivativeFromHiddenLayer[h][j]*sum;
         }
      }

      // Calculate gradient elements

      index = 0;

      // Calculate gradient elements of hidden neurons

      // First hidden layer

      for(int j = 0; j < numbersOfHiddenNeurons[0]; j++)
      {
         // Bias

         gradient[index] += hiddenError[0][j];
         index++;

         // Synaptic weights

         synapticWeights = hiddenLayers[0][j].getSynapticWeights();

         for(int k = 0; k < numberOfInputs; k++)
         {
            gradient[index] += hiddenError[0][j]*inputSignal[k];
            index++;   
         }
      }

      // Rest of hidden layers	
    
      for(int h = 1; h < numberOfHiddenLayers; h++)
      {      
         for(int j = 0; j < numbersOfHiddenNeurons[h]; j++)
         {
            // Bias

            gradient[index] += hiddenError[h][j];
            index++;

            // Synaptic weights

            synapticWeights = hiddenLayers[h][j].getSynapticWeights();

            for(int k = 0; k < numbersOfHiddenNeurons[h-1]; k++)
            {
               gradient[index] += hiddenError[h][j]*outputSignalFromHiddenLayer[h-1][k];
               index++;   
            }
	 }
      }

      // Calculate gradient elements of output neurons

      for(int j = 0; j < numberOfOutputs; j++)
      {
         // Bias

         gradient[index] += outputError[j];
         index++;

         // Synaptic weights

	 for(int k = 0; k < numbersOfHiddenNeurons[numberOfHiddenLayers-1]; k++)
         {
            gradient[index] += outputSignalFromHiddenLayer[numberOfHiddenLayers-1][k]*outputError[j];
            index++;
         }
      }
   }

   return(gradient);
}


// void saveInputTargetAndOutput(char*) method

/// This method saves to a file the inputs and the targets in the input-target data set, together with the 
/// outputs from the multilayer perceptron for that inputs:
///
/// <ul>
/// <li> Inputs.
/// <li> Targets.
/// <li> Outputs.
/// </ul> 
///
/// @param filename Filename.

void MeanSquaredError::saveInputTargetAndOutput(char* filename)
{
   int numberOfOutputs = multilayerPerceptron->getNumberOfOutputs();

   int numberOfSamples = inputTargetDataSet->getNumberOfSamples();
   int numberOfInputVariables = inputTargetDataSet->getNumberOfInputVariables();
   int numberOfTargetVariables = inputTargetDataSet->getNumberOfTargetVariables();

   Matrix<double> inputData = inputTargetDataSet->getInputData();
   Matrix<double> targetData = inputTargetDataSet->getTargetData();

   Matrix<double> outputData(numberOfSamples, numberOfOutputs, 0.0);

   Vector<double> input(numberOfInputVariables, 0.0);
   Vector<double> output(numberOfTargetVariables, 0.0);

   for(int sample = 0; sample < numberOfSamples; sample++)
   {
      input = inputData.getRow(sample);

      output = multilayerPerceptron->calculateOutput(input);

      outputData.setRow(sample, output);
   }

   std::fstream file;

   file.open(filename, std::ios::out);

   // Control sentence

   if(!file.is_open())
   {
      std::cout << std::endl
                << "Flood Error MeanSquaredError class." << std::endl
                << "void saveInputTargetAndOutput(char*) method." << std::endl
                << "Cannot open input-target-output data file." << std::endl
                << std::endl;

      exit(1);
   }
   else
   {
      if(display)
     {
         std::cout << std::endl
                   << "Saving input-target-output to data file..."
                   << std::endl;
     }
   }

   // Write file header

   file << "% Flood Neural Network. Input-target-output data file." << std::endl
        << "% 1 - Input data" << std::endl
        << "% 2 - Target data" << std::endl
        << "% 3 - Output data" << std::endl
        << std::endl;

   // Write file data

   for(int sample = 0; sample < numberOfSamples; sample++)
   {
      // Write sample input data

      for(int i = 0; i < numberOfInputVariables; i++)
      {
         file << inputData[sample][i] << " ";
      }

      // Write sample target data

      for(int i = 0; i < numberOfTargetVariables; i++)
      {
         file << targetData[sample][i] << " ";
      }

      // Write sample output data

      for(int i = 0; i < numberOfOutputs; i++)
      {
         file << outputData[sample][i] << " ";
      }

      file << std::endl;
   }

   file << std::endl;

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
