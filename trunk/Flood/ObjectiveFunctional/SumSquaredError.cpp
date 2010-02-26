/****************************************************************************************************************/
/*                                                                                                              */
/*   Flood: An Open Source Neural Networks C++ Library                                                          */
/*   www.cimne.com/flood                                                                                        */
/*                                                                                                              */
/*   S U M   S Q U A R E D   E R R O R   C L A S S                                                              */
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
#include <math.h>

#include "SumSquaredError.h"

namespace Flood
{

// GENERAL CONSTRUCTOR

/// General constructor. It creates a sum squared error objective functional object associated to a multilayer 
/// perceptron and measured on an input-target data set.
/// It also initializes all the rest of class members to their default values:
///
/// <ul>
/// <li> Numerical epsilon method: Relative.
/// <li> Numerical differentiation method: Central differences.
/// <li> Number of evaluations = 0. 
/// <li> Numerical epsilon: 1.0e-6.
/// <li> Display = true. 
/// </ul> 
///
/// @param newMultilayerPerceptron Pointer to a multilayer perceptron object.
/// @param newInputTargetDataSet Pointer to an input-target data set object.

SumSquaredError::SumSquaredError(MultilayerPerceptron* newMultilayerPerceptron, 
InputTargetDataSet* newInputTargetDataSet): ObjectiveFunctional(newMultilayerPerceptron)
{
   inputTargetDataSet = newInputTargetDataSet;
}


// DEFAULT CONSTRUCTOR

/// Default constructor. It creates a sum squared error objective functional object not associated to any 
/// multilayer perceptron and not measured on any input-target data set.
/// It also initializes all the rest of class members to their default values:
///
/// <ul>
/// <li> Numerical epsilon method: Relative.
/// <li> Numerical differentiation method: Central differences.
/// <li> Number of evaluations = 0. 
/// <li> Numerical epsilon: 1.0e-6.
/// <li> Display = true. 
/// </ul> 

SumSquaredError::SumSquaredError(void) : ObjectiveFunctional()
{
   inputTargetDataSet = NULL;
}


// DESTRUCTOR
//
/// Destructor.

SumSquaredError::~SumSquaredError(void) 
{

}


// METHODS

// InputTargetDataSet* getInputTargetDataSet(void) method

/// This method returns a pointer to the input-target data set object on which the objective functional is 
/// measured.

InputTargetDataSet* SumSquaredError::getInputTargetDataSet(void)
{
   return(inputTargetDataSet);
}


// void setInputTargetDataSet(InputTargetDataSet*) method

/// This method sets a pointer to an input-data set object on which the objective functional is to be measured.
///
/// @param newInputTargetDataSet Pointer to an input-target data set object.

void SumSquaredError::setInputTargetDataSet(InputTargetDataSet* newInputTargetDataSet)
{
   inputTargetDataSet = newInputTargetDataSet;
}


// double calculateEvaluation(void) method

/// This method returns the objective functional value of a multilayer perceptron according to the sum squared 
/// error on an input-target data set.
///
/// @see calculateGradient(void).

double SumSquaredError::calculateEvaluation(void)
{
   // Control sentence (if debug)

   #ifndef NDEBUG 

   if(multilayerPerceptron == NULL)
   {
      std::cerr << std::endl
                << "Flood Error: SumSquaredError class." << std::endl
                << "double calculateEvaluation(void) method." << std::endl
                << "Pointer to multilayer perceptron object cannot be NULL." << std::endl
                << std::endl;

        exit(1);
   }
   else if(inputTargetDataSet == NULL)
   {
      std::cerr << std::endl
                << "Flood Error: SumSquaredError class." << std::endl
                << "double calculateEvaluation(void) method." << std::endl
                << "Pointer to input-target data set object cannot be NULL." << std::endl
                << std::endl;

      exit(1);
   }

   #endif

   int numberOfInputs = multilayerPerceptron->getNumberOfInputs();
   int numberOfOutputs = multilayerPerceptron->getNumberOfOutputs();

   int numberOfSamples = inputTargetDataSet->getNumberOfSamples();

   #ifndef NDEBUG 

   int numberOfInputVariables = inputTargetDataSet->getNumberOfInputVariables();
   int numberOfTargetVariables = inputTargetDataSet->getNumberOfTargetVariables();

   if(numberOfInputs != numberOfInputVariables || numberOfOutputs != numberOfTargetVariables)
   {
      std::cout << std::endl
                << "Flood Error: SumSquaredError class." << std::endl
                << "double calculateEvaluation(void) method." << std::endl
                << "Number of inputs and outputs in multilayer perceptron must be equal to " 
                << "number of input and output variables in input-target data set." 
                << std::endl << std::endl;

      exit(1);
   }

   #endif

   Matrix<double>& inputData = inputTargetDataSet->getInputData();
   Matrix<double>& targetData = inputTargetDataSet->getTargetData();

   // Increment number of evaluations

   numberOfEvaluations++;

   double sumSquaredError = 0;   

   Vector<double> input(numberOfInputs);
   Vector<double> output(numberOfOutputs);
   Vector<double> target(numberOfOutputs);

   for(int i = 0; i < numberOfSamples; i++)
   {
      // Input vector

      input = inputData.getRow(i);

      // Output vector

      output = multilayerPerceptron->calculateOutput(input);

      // Target vector

     target = targetData.getRow(i);

      // Sum of squares error

      sumSquaredError += (output - target).dot(output - target);           
   }

   return(sumSquaredError);
}


// Vector<double> calculateGradient(void)

/// This method returns the the sum squared error function gradient of a multilayer perceptron on an input-target 
/// data set. It uses the error back-propagation method.
///
/// @see calculateEvaluation(void).

Vector<double> SumSquaredError::calculateGradient(void)
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
            outputError = outputSignalDerivative*(output-target)*2.0;
         }//end none   

         break;

         case MultilayerPerceptron::MeanAndStandardDeviation:
         {
            outputError = outputSignalDerivative*(output-target)*standardDeviationOfOutputVariables*2.0;
         }//end mean and standard deviation

         break;

         case MultilayerPerceptron::MinimumAndMaximum:
         {
            outputError = outputSignalDerivative*(output-target)*(maximumOfOutputVariables-minimumOfOutputVariables);
         }//end minimum and maximum		 

         break;
      }//end switch

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


// Matrix<double> calculateJacobian(void) method

/// This method returns the Jacobian matrix of the sum squared error function, whose elements are given by the 
/// derivatives of the squared errors data set with respect to the neural parameters.
/// The Jacobian matrix here is computed using numerical differentiation.
///
/// @see calculateJacobian(void).
///
/// @todo This method is not fully implemented, and it should be finished. 

Matrix<double> SumSquaredError::calculateJacobian(void)
{
   // Multilayer perceptron stuff

   int numberOfInputs = multilayerPerceptron->getNumberOfInputs();

   int numberOfHiddenLayers = multilayerPerceptron->getNumberOfHiddenLayers();
   Vector<int> numbersOfHiddenNeurons = multilayerPerceptron->getNumbersOfHiddenNeurons();

   int numberOfOutputs = multilayerPerceptron->getNumberOfOutputs();

   MultilayerPerceptron::PreAndPostProcessingMethod preAndPostProcessingMethod
   = multilayerPerceptron->getPreAndPostProcessingMethod();

   Vector<double> meanOfInputVariables = multilayerPerceptron->getMeanOfInputVariables();
   Vector<double> standardDeviationOfInputVariables = multilayerPerceptron->getStandardDeviationOfInputVariables();

   Vector<double> meanOfOutputVariables = multilayerPerceptron->getMeanOfOutputVariables();
   Vector<double> standardDeviationOfOutputVariables = multilayerPerceptron->getStandardDeviationOfOutputVariables();

   Vector<double> minimumOfInputVariables = multilayerPerceptron->getMinimumOfInputVariables();
   Vector<double> maximumOfInputVariables = multilayerPerceptron->getMaximumOfInputVariables();

   Vector<double> minimumOfOutputVariables = multilayerPerceptron->getMinimumOfOutputVariables();
   Vector<double> maximumOfOutputVariables = multilayerPerceptron->getMaximumOfOutputVariables();

   Vector<Perceptron>& outputLayer = multilayerPerceptron->getOutputLayer();
   Vector< Vector<Perceptron> >& hiddenLayers = multilayerPerceptron->getHiddenLayers();

   int numberOfNeuralParameters = multilayerPerceptron->getNumberOfNeuralParameters();

   // Input-target data set stuff

   int numberOfSamples = inputTargetDataSet->getNumberOfSamples();
   int numberOfTargetVariables = inputTargetDataSet->getNumberOfTargetVariables();

   Matrix<double>& inputData = inputTargetDataSet->getInputData();
   Matrix<double>& targetData = inputTargetDataSet->getTargetData();

   Vector<double> target(numberOfTargetVariables);

   // Back-popagation algorithm stuff

   Vector<double> input(numberOfInputs);
   Vector<double> inputSignal(numberOfInputs);
   Vector< Vector<double> > netInputSignalToHiddenLayer(numberOfHiddenLayers);
   Vector< Vector<double> > outputSignalFromHiddenLayer(numberOfHiddenLayers);
   Vector< Vector<double> > outputSignalDerivativeFromHiddenLayer(numberOfHiddenLayers);
   Vector<double> netInputSignalToOutputLayer(numberOfOutputs);
   Vector<double> outputSignal(numberOfOutputs);
   Vector<double> output(numberOfOutputs);

   Vector<double> gradient(numberOfNeuralParameters, 0.0);
   Matrix<double> jacobian(numberOfSamples, numberOfNeuralParameters);
/*

   Vector<double> outputError(numberOfOutputs);
   Vector< Vector<double> > hiddenError(numberOfHiddenLayers);

   // Main loop

   int index = 0;

   for(int i = 0; i < numberOfSamples; i++)
   {
      index = 0;

      // Actual input

      input = inputData.getRow(i);

      // Preprocess input to obtain input signal to the neural network 

      inputSignal = multilayerPerceptron->preprocessInput(input);      

      // Calculate net input output and output signal derivative from all hidden layers

      netInputSignalToHiddenLayer[0] 
      = multilayerPerceptron->calculateNetInputSignalToHiddenLayer(0, inputSignal);

      outputSignalFromHiddenLayer[0] 
      = multilayerPerceptron->calculateOutputSignalFromHiddenLayer(0, netInputSignalToHiddenLayer[0]);

      outputSignalDerivativeFromHiddenLayer[0] 
      = multilayerPerceptron->calculateOutputSignalDerivativeFromHiddenLayer(0, netInputSignalToHiddenLayer[0]);
  
      for(int hiddenLayer = 1; hiddenLayer < numberOfHiddenLayers; hiddenLayer++)
      {
         netInputSignalToHiddenLayer[hiddenLayer] = multilayerPerceptron
         ->calculateNetInputSignalToHiddenLayer(hiddenLayer, outputSignalFromHiddenLayer[hiddenLayer-1]);
   
         outputSignalFromHiddenLayer[hiddenLayer] = multilayerPerceptron
         ->calculateOutputSignalFromHiddenLayer(hiddenLayer, netInputSignalToHiddenLayer[hiddenLayer]);

         outputSignalDerivativeFromHiddenLayer[hiddenLayer] = multilayerPerceptron
         ->calculateOutputSignalDerivativeFromHiddenLayer(hiddenLayer, netInputSignalToHiddenLayer[hiddenLayer]);         
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

      target = targetData.getRow(i);
         
      switch(preAndPostProcessingMethod)
      {
         case MultilayerPerceptron::None:
         {
            outputError = outputSignalDerivative*2.0*(output - target);
         }//end none
         break;

         case MultilayerPerceptron::MeanAndStandardDeviation:
         {
            outputError = outputSignalDerivative*2.0*(output - target)*standardDeviationOfOutputVariables;
         }//end mean and standard deviation
         break;

         case MultilayerPerceptron::MinimumAndMaximum:
         {
            outputError = outputSignalDerivative*2.0*(output-target)
                         *0.5*(maximumOfOutputVariables-minimumOfOutputVariables);
         }//end minimum and maximum
         break;
      }

      // Backpropagate the errors of the output units to obtain the error for each hidden unit 

      double sum = 0.0;
     
      for(int j = 0; j < numbersOfHiddenNeurons[0]; j++)
      {
         sum = 0.0;

         for(int k = 0; k < numberOfOutputs; k++)
         {
            Vector<double> synapticWeights = outputLayer[k].getSynapticWeights();

            sum += (synapticWeights[j])*outputError[k];
         }

         netInputSignalToHiddenLayer[0][j] = hiddenLayers[0][j].calculateNetInputSignal(inputSignal);
         
         outputSignalDerivativeFromHiddenLayer[0][j] 
         = hiddenLayers[0][j].calculateOutputSignalDerivative(netInputSignalToHiddenLayer[j]);

         hiddenError[j] = outputSignalDerivativeFromHiddenLayer[j]*sum;
      }

      // Hidden layer

      for(int j = 0; j < numbersOfHiddenNeurons[0]; j++)
      {
         // Bias

         //gradient[index] += hiddenError[j];
         jacobian[i][index] += hiddenError[0][j];
         index++;

         // Synaptic weights

         Vector<double> synapticWeights = hiddenLayers[0][j].getSynapticWeights();

         for(int k = 0; k < numberOfInputs; k++)
         {
            //gradient[index] += hiddenError[j]*inputSignal[k];
            jacobian[i][index] += hiddenError[0][j]*inputSignal[k];

            index++;   
         }
      }

      // Output layer

      for(int j = 0; j < numberOfOutputs; j++)
      {
         // Bias

         //gradient[index] += outputError[j];
         jacobian[i][index] += outputError[j];
         index++;

         // Synaptic weights

         for(int k = 0; k < numbersOfHiddenNeurons[0]; k++)
         {
            //gradient[index] += outputSignalFromHiddenLayer[k]*outputError[j];
            jacobian[i][index] += outputSignalFromHiddenLayer[0][k]*outputError[j];
            index++;
         }
      }

  }
*/
   return(jacobian);
}


// Matrix<double> calculateJacobian(void) method

/// This method returns the Jacobian matrix of the sum squared error function, whose elements are given by the 
/// derivatives of the squared errors data set with respect to the neural parameters.
/// The Jacobian matrix here is computed using a back-propagation algorithm.
///
/// @see calculateJacobianTest(void).

Matrix<double> SumSquaredError::calculateJacobianTest(void)
{
   // Multilayer perceptron stuff

   Vector<double> neuralParameters = multilayerPerceptron->getNeuralParameters();
   int numberOfNeuralParameters = multilayerPerceptron->getNumberOfNeuralParameters();

    // Input target data set stuff

   int numberOfSamples = inputTargetDataSet->getNumberOfSamples();

   Matrix<double> jacobian(numberOfSamples, numberOfNeuralParameters);

   Vector<double> column(numberOfNeuralParameters);

   switch(numericalDifferentiationMethod)   
   {
      case ForwardDifferences:
      {
         Vector<double> squaredErrors = calculateSquaredErrors(); 
         Vector<double> squaredErrorsForward(numberOfSamples);

         for(int j = 0; j < numberOfNeuralParameters; j++)
         {
            // Perturbate neural parameters

            neuralParameters[j] += numericalEpsilon;
            multilayerPerceptron->setNeuralParameters(neuralParameters);

            // Calculate squared errors forward

            squaredErrorsForward = calculateSquaredErrors();

            // Restart biases and synaptic weights

            neuralParameters[j] -= numericalEpsilon;
            multilayerPerceptron->setNeuralParameters(neuralParameters);

            // Calculate Jacobian elements
                        
            for(int i = 0; i < numberOfSamples; i++)
            {
               jacobian[i][j] = (squaredErrorsForward[i] - squaredErrors[i])/numericalEpsilon;
            }
         }

      }// end forward differences
      break;

      case CentralDifferences:
      {
         Vector<double> squaredErrorsForward(numberOfSamples);
         Vector<double> squaredErrorsBackward(numberOfSamples);

         for(int j = 0; j < numberOfNeuralParameters; j++)
         {
            // Perturbate neural parameters

            neuralParameters[j] += numericalEpsilon;
            multilayerPerceptron->setNeuralParameters(neuralParameters);

			// Calculate squared errors     
 
            squaredErrorsForward = calculateSquaredErrors();
 
            // Restart neural parameters

            neuralParameters[j] -= numericalEpsilon;
            multilayerPerceptron->setNeuralParameters(neuralParameters);

			// Perturbate neural parameters

            neuralParameters[j] -= numericalEpsilon;
            multilayerPerceptron->setNeuralParameters(neuralParameters);

            // Calculate squared errors     
 
            squaredErrorsBackward = calculateSquaredErrors();
 
            // Restart neural parameters

            neuralParameters[j] += numericalEpsilon;
            multilayerPerceptron->setNeuralParameters(neuralParameters);
            // Calculate Jacobian elements

            for(int i = 0; i < numberOfSamples; i++)
            {
               jacobian[i][j] = (squaredErrorsForward[i] - squaredErrorsBackward[i])/(2.0*numericalEpsilon);
            }
         }
      }// end central differences
      break;
   }


   return(jacobian);
}


// Vector<double> calculateSquaredErrors(void) method

/// This method returns a vector containing the squared error values between the outputs from the neural network 
/// and the targets in the data set. 
/// The size of this vector is the number of samples in the input-target data set.

Vector<double> SumSquaredError::calculateSquaredErrors(void)
{
   int numberOfInputs = multilayerPerceptron->getNumberOfInputs();
   int numberOfOutputs = multilayerPerceptron->getNumberOfOutputs();

   // Input-target data set

   int numberOfSamples = inputTargetDataSet->getNumberOfSamples();
   int numberOfInputVariables = inputTargetDataSet->getNumberOfInputVariables();
   int numberOfTargetVariables = inputTargetDataSet->getNumberOfTargetVariables();

   if(numberOfInputs != numberOfInputVariables || numberOfOutputs != numberOfTargetVariables)
   {
      std::cout << std::endl
                << "Flood Error SumSquaredError class." << std::endl
                << "double calculateSquaredErrors(void) method." << std::endl
                << "Number of inputs and outputs in multilayer perceptron must be equal to number of input and "
                << "output variables in input-target data set." 
                << std::endl << std::endl;

      exit(1);
   }

   Matrix<double>& inputData = inputTargetDataSet->getInputData();
   Matrix<double>& targetData = inputTargetDataSet->getTargetData();

   Vector<double> squaredErrors(numberOfSamples);

   Vector<double> input(numberOfInputs);
   Vector<double> output(numberOfOutputs);
   Vector<double> target(numberOfTargetVariables);

   for(int i = 0; i < numberOfSamples; i++)
   {
      // Input vector

      input = inputData.getRow(i);

      // Output vector

      output = multilayerPerceptron->calculateOutput(input);

      // Target vector

      target = targetData.getRow(i);

      // Squared error
      
      squaredErrors[i] = (output - target).dot(output - target);
   }

   return(squaredErrors);
}


// void saveInputTargetAndOutput(char*) method

/// This method saves to a data file the inputs and the targets in the input-target data set, together with the 
/// outputs from the multilayer perceptron for that inputs. 
///
/// @param filename Filename.

void SumSquaredError::saveInputTargetAndOutput(char* filename)
{
   int numberOfInputs = multilayerPerceptron->getNumberOfInputs();
   int numberOfOutputs = multilayerPerceptron->getNumberOfOutputs();

   int numberOfSamples = inputTargetDataSet->getNumberOfSamples();
   int numberOfInputVariables = inputTargetDataSet->getNumberOfInputVariables();
   int numberOfTargetVariables = inputTargetDataSet->getNumberOfTargetVariables();

   Matrix<double> inputData = inputTargetDataSet->getInputData();
   Matrix<double> targetData = inputTargetDataSet->getTargetData();

   Matrix<double> outputData(numberOfSamples, numberOfOutputs);

   Vector<double> input(numberOfInputs);
   Vector<double> output(numberOfOutputs);

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
                << "Flood Error: SumSquaredError class." << std::endl
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
                   << "Saving input-target-output to data file..." << std::endl;
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
