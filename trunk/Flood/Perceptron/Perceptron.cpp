/****************************************************************************************************************/
/*                                                                                                              */
/*   Flood: An Open Source Neural Networks C++ Library                                                          */
/*   www.cimne.com/flood                                                                                        */
/*                                                                                                              */
/*   P E R C E P T R O N   C L A S S                                                                            */
/*                                                                                                              */
/*   Roberto Lopez                                                                                              */
/*   International Center for Numerical Methods in Engineering (CIMNE)                                          */
/*   Technical University of Catalonia (UPC)                                                                    */
/*   Barcelona, Spain                                                                                           */
/*   E-mail: rlopez@cimne.upc.edu                                                                               */
/*                                                                                                              */
/****************************************************************************************************************/


#include <iostream>
#include <stdlib.h>
#include <time.h>

#include "Perceptron.h"

namespace Flood
{

/// General constructor. It creates a perceptron object with a given number of inputs. 
/// The neuron's free paramameters (bias and synaptic weights) are initialized at random with a normal 
/// distribution of mean 0 and standard deviation 1.
/// This constructor also initializes the rest of class members to their default values:
///
/// <ul>
/// <li> Activation function: Hyperbolic tangent.
/// <li> Display: True.
/// </ul> 
///
/// @param newNumberOfInputs Number of inputs in the neuron.

Perceptron::Perceptron(int newNumberOfInputs)
{
   // Set synaptic weights

   activationFunction = HyperbolicTangent;

   numberOfInputs = newNumberOfInputs;

   synapticWeights.setSize(numberOfInputs);

   // Init bias and synaptic weights at random 
   
   initBiasNormal(0.0, 1.0);
   initSynapticWeightsNormal(0.0, 1.0);
   
   display = true;
}


/// Default constructor. It creates a perceptron object with zero inputs.
/// The neuron's bias is initialized to zero. 
/// This constructor also initializes the rest of class members to their default values:
///
/// <ul>
/// <li> Activation function: Hyperbolic tangent.
/// <li> Display: True.
/// </ul> 
///

Perceptron::Perceptron(void)
{
   activationFunction = HyperbolicTangent;

   numberOfInputs = 0;

   bias = 0.0;
   
   display = true;
}


/// Destructor.

Perceptron::~Perceptron(void)
{

}


// METHODS

// ActivationFunction getActivationFunction(void) method

/// This method returns the activation function of the neuron. 

Perceptron::ActivationFunction Perceptron::getActivationFunction(void)
{
   return(activationFunction);                           
}


// int getNumberOfInputs(void) method

/// This method returns the number of inputs in the neuron. 

int Perceptron::getNumberOfInputs(void)
{
   return(numberOfInputs);
}


// double getBias(void) method

/// This method returns the bias value of the neuron.
///
/// @see getSynapticWeights(void).
/// @see getSingleSynapticWeight(int).

double Perceptron::getBias(void)
{
   return(bias);
}


// Vector<double> getSynapticWeights(void)

/// This method returns the synaptic weight values of the neuron.
///
/// @see getBias(void).
/// @see getSingleSynapticWeight(int).

Vector<double> Perceptron::getSynapticWeights(void)
{
   return(synapticWeights);
}


// double getSingleSynapticWeight(int) method

/// This method returns the synaptic weight value with index i of the neuron.
///
/// @param i Synaptic weight index.
///
/// @see getBias(void).
/// @see getSynapticWeights(void).
 
double Perceptron::getSingleSynapticWeight(int i)
{
   // Control sentence (if debug)

   #ifndef NDEBUG 

   if(i >= numberOfInputs)
   {
      std::cerr << std::endl 
                << "Flood Error: Perceptron class." << std::endl
                << "double getSingleSynapticWeight(int) method." << std::endl
                << "Index must be less than number of inputs." << std::endl
                << std::endl;

      exit(1);
   }

   #endif

   // Get single synaptic weights

   return(synapticWeights[i]);   
}


// bool getDisplay(void) method

/// This method returns true if messages from this class are to be displayed on the screen, or false if messages 
/// from this class are not to be displayed on the screen.

bool Perceptron::getDisplay(void)
{
   return(display);
}


// void setActivationFunction(ActivationFunction) method

/// This method sets a new activation function in the neuron. 
///
/// @param newActivationFunction Activation function.

void Perceptron::setActivationFunction(Perceptron::ActivationFunction newActivationFunction)
{
   activationFunction = newActivationFunction;
}


// void setBias(double) method

/// This method sets a new bias value for the perceptron.
///
/// @param newBias Bias value.
///
/// @see setSynapticWeights(Vector<double>).
/// @see setSingleSynapticWeight(i, double).
 
void Perceptron::setBias(double newBias)
{
   bias = newBias;   
}


// void setSynapticWeights(Vector<double>) method

/// This method a new set of synaptic weights for the perceptron.
///
/// @param newSynapticWeights Synaptic weight values.
///
/// @see setBias(double).
/// @see setSingleSinapticWeight(int, double).
 
void Perceptron::setSynapticWeights(Vector<double> newSynapticWeights)
{
   // Control sentence (if debug)

   #ifndef NDEBUG 

   if(newSynapticWeights.getSize() != numberOfInputs)
   {
      std::cerr << std::endl 
                << "Flood Error: Perceptron class." << std::endl
                << "void setSynapticWeights(Vector<double>) method." << std::endl
                << "Size of synaptic weights vector must be equal to number of inputs." << std::endl
                << std::endl;

      exit(1);
   }

   #endif

   // Set synaptic weights
   
   synapticWeights = newSynapticWeights;
}


// void setSingleSynapticWeight(int, double) method

/// This method sets the synaptic weight value with index i for the neuron.
///
/// @param i Synaptic weight index.
/// @param newSynapticWeight Synaptic weight value.
///
/// @see setBias(double).
/// @see setSingleSinapticWeight(int, double).

void Perceptron::setSingleSynapticWeight(int i, double newSynapticWeight)
{
   // Control sentence (if debug)

   #ifndef NDEBUG 

   if(i >= numberOfInputs)
   {
      std::cerr << std::endl 
                << "Flood Error: Perceptron class." << std::endl
                << "void setSingleSynapticWeight(int, double) method." << std::endl
                << "Index must be less than number of inputs." << std::endl
                << std::endl;

      exit(1);
   }

   #endif

   // Set single synaptic weight

   synapticWeights[i] = newSynapticWeight;
}


// void setDisplay(bool) method

/// This method sets a new display value. 
/// If it is set to true messages from this class are to be displayed on the screen;
/// if it is set to false messages from this class are not to be displayed on the screen.
///
/// @param newDisplay Display value.
 
void Perceptron::setDisplay(bool newDisplay)
{
   display = newDisplay;   
}


// void setNumberOfInputs(int) method

/// This method sets a new number of inputs in the neuron.
/// The new free paramameters (bias and synaptic weights) are initialized at random with a normal distribution of 
/// mean 0 and standard deviation 1.
///
/// @param newNumberOfInputs Number of inputs in the neuton.
 
void Perceptron::setNumberOfInputs(int newNumberOfInputs)
{
   numberOfInputs = newNumberOfInputs;

   synapticWeights.setSize(numberOfInputs);

   // Init bias and synaptic weights at random 

   initBiasNormal(0.0,1.0);
   initSynapticWeightsNormal(0.0,1.0);
}


// int getNumberOfNeuronParameters(void) method

int Perceptron::getNumberOfNeuronParameters(void)
{
   int numberOfNeuronParameters = 1 + numberOfInputs;

   return(numberOfNeuronParameters);
}


// Vector<double> getNeuronParameters(void) method

Vector<double> Perceptron::getNeuronParameters(void)
{
   int numberOfNeuronParameters = getNumberOfNeuronParameters();

   Vector<double> neuronParameters(numberOfNeuronParameters);

   neuronParameters[0] = bias;

   for(int i = 0; i < numberOfInputs; i++)
   {
      neuronParameters[1+i] = synapticWeights[i];
   }

   return(neuronParameters); 
}


// void setNeuronParameters(Vector<double>) method

void Perceptron::setNeuronParameters(Vector<double> newNeuronParameters)
{
   bias = newNeuronParameters[0];

   for(int i = 0; i < numberOfInputs; i++)
   {
      synapticWeights[i] = newNeuronParameters[i+1];
   }
}


// void initBiasUniform(double, double) method

/// This method initializes the neuron's bias with a random value chosen from a uniform distribution.
///
/// @param minimumValue Minimum initialization value.
/// @param maximumValue Maximum initialization value.
///
/// @see initSynapticWeightsUniform(double, double).
/// @see initBiasNormal(double, double).
/// @see initSynapticWeightsNormal(double, double).

void Perceptron::initBiasUniform(double minimumValue, double maximumValue)
{
   // Control sentence (if debug)

   #ifndef NDEBUG 

   if(minimumValue > maximumValue)
   {
      std::cerr << std::endl
                << "Flood Error: Perceptron class." << std::endl 
                << "initBiasAtRandom(double, double) method." << std::endl
                << "Minimum value must be less than maximum value." << std::endl
                << std::endl;

      exit(1);
   }

   #endif

   double random = (double)rand()/(RAND_MAX+1.0);

   bias = minimumValue + (maximumValue-minimumValue)*random;
}


// void initSynapticWeightsUniform(double, double) method

/// This method initializes the neuron's synaptic weights with random values chosen from an uniform distribution.
///
/// @param minimumValue Minimum initialization value.
/// @param maximumValue Maximum initialization value.
///
/// @see initBiasUniform(double, double).
/// @see initBiasNormal(double, double).
/// @see initSynapticWeightsNormal(double, double).

void Perceptron::initSynapticWeightsUniform(double minimumValue, double maximumValue)
{
   // Control sentence (if debug)

   #ifndef NDEBUG 

   if(minimumValue > maximumValue)
   {
      std::cerr << std::endl
                << "Flood Error: Perceptron class." << std::endl 
                << "initSynapticWeightsAtRandom(double, double) method." << std::endl
                << "Minimum value must be less than maximum value." << std::endl
                << std::endl;

      exit(1);
   }

   #endif

   // Init synaptic weights

   for(int i = 0; i < numberOfInputs; i++)
   {
      double random = (double)rand()/(RAND_MAX+1.0);

      synapticWeights[i] = minimumValue + (maximumValue-minimumValue)*random;
   }
}


// void initBiasNormal(double, double) method

/// This method initializes the neuron's bias with random values chosen from a normal distribution.
///
/// @param mean Mean of normal distribution.
/// @param standardDeviation Standard deviation of normal distribution.
///
/// @see initBiasUniform(double, double).
/// @see initSynapticWeightsUniform(double, double).
/// @see initSynapticWeightsNormal(double, double).

void Perceptron::initBiasNormal(double mean, double standardDeviation)
{
   // Control sentence (if debug)

   #ifndef NDEBUG 

   if(standardDeviation < 0.0)
   {
      std::cerr << std::endl
                << "Flood Error: Perceptron class." << std::endl 
                << "initBiasNormal(double, double) method." << std::endl
                << "Standard deviation must be equal or greater than zero." << std::endl
                << std::endl;

      exit(1);
   }

   #endif

   double pi = 4.0*atan(1.0);

   double random1 = (double)rand()/(RAND_MAX+1.0);
   double random2 = (double)rand()/(RAND_MAX+1.0);

   // Box-Muller transformation

   bias = mean + sqrt(-2.0*log(random1))*sin(2.0*pi*random2)*standardDeviation;
}


// void initSynapticWeightsNormal(double, double) method

/// This method initializes the neuron's synaptic weights with random values chosen from a normal distribution.
///
/// @param mean Mean of normal distribution.
/// @param standardDeviation Standard deviation of normal distribution.
///
/// @see initBiasUniform(double, double).
/// @see initSynapticWeightsUniform(double, double).
/// @see initBiasNormal(double, double).

void Perceptron::initSynapticWeightsNormal(double mean, double standardDeviation)
{
   // Control sentence (if debug)

   #ifndef NDEBUG 

   if(standardDeviation < 0.0)
   {
      std::cerr << std::endl
                << "Flood Error: Perceptron class." << std::endl 
                << "initSynapticWeightsNormal(double, double) method." << std::endl
                << "Standard deviation must be equal or greater than zero." << std::endl
                << std::endl;

      exit(1);
   }

   #endif

   double pi = 4.0*atan(1.0);

   for(int i = 0; i < numberOfInputs; i++)
   {
      double random1 = (double)rand()/(RAND_MAX+1.0);
      double random2 = (double)rand()/(RAND_MAX+1.0);

      // Box-Muller transformation

      double normallyDistributedRandomNumber 
      = mean + sqrt(-2.0*log(random1))*sin(2.0*pi*random2)*standardDeviation;

      synapticWeights[i] = normallyDistributedRandomNumber;
   }
}


// double calculateNetInputSignal(Vector<double>) method

/// This method returns the net input signal to the neuron for a set of input signals, using the dot product 
/// combination function. 
///
/// @param inputSignal Set of input signals to the neuron.

double Perceptron::calculateNetInputSignal(Vector<double> inputSignal)
{
   // Control sentence (if debug)

   #ifndef NDEBUG 

   if(numberOfInputs == 0)
   {
      std::cerr << std::endl
                << "Flood Error: Perceptron class." << std::endl 
                << "calculateNetInputSignal(Vector<double>) method." << std::endl
                << "Number of inputs in the neuron must be greater than zero." << std::endl
                << std::endl;

      exit(1);
   }
   else if(inputSignal.getSize() != numberOfInputs)
   {
      std::cerr << std::endl
                << "Flood Error: Perceptron class." << std::endl
                << "double calculateNetInputSignal(Vector<double>) method." << std::endl
                << "Size of input signal must be equal to number of inputs." << std::endl
                << std::endl;

      exit(1);
   }

   #endif

   // Calculate net input signal

   double netInputSignal = bias + inputSignal.dot(synapticWeights);

   return(netInputSignal);   
}


// double calculateOutputSignal(double) method

/// This method returns the output signal from the neuron for a net input signal.
/// The output signal depends on the activation function used.
/// 
/// @param netInputSignal Net input signal to the neuron.
///
/// @see calculateOutputSignal(Vector<double>).
/// @see calculateOutputSignalDerivative(double).
/// @see calculateOutputSignalDerivative(Vector<double>).
/// @see getOutputSignaSecondDerivative(double).
/// @see getOutputSignaSecondDerivative(Vector<double>).

double Perceptron::calculateOutputSignal(double netInputSignal)
{
   double outputSignal = 0.0;

   switch(activationFunction)   
   {
      case Logistic:
      {
         outputSignal = 1.0/(1.0 + exp(-netInputSignal));
      }
      break;
                                     
      case HyperbolicTangent:
      {
         outputSignal = tanh(netInputSignal);
      }
      break;

      case Linear:
      {
         outputSignal = netInputSignal;        
      }
      break;
      
      default:
         std::cerr << std::endl
                   << "Flood Error: Perceptron class." << std::endl 
                   << "calculateOutputSignal(Vector<double>) method." << std::endl
                   << "Unknown activation function." << std::endl
                   << std::endl;

         exit(1);
      break;
   }

   return(outputSignal);
}


// double calculateOutputSignal(Vector<double>) method

/// This method returns the output signal from the neuron for a set of input signals.
/// The output signal depends on the activation function used.
/// 
/// @param inputSignal Set of input signals to the neuron.
///
/// @see calculateOutputSignal(double).
/// @see calculateOutputSignalDerivative(double).
/// @see calculateOutputSignalDerivative(Vector<double>).
/// @see getOutputSignaSecondDerivative(double).
/// @see getOutputSignaSecondDerivative(Vector<double>).

double Perceptron::calculateOutputSignal(Vector<double> inputSignal)
{
   // Control sentence (if debug)

   #ifndef NDEBUG 

   int size = inputSignal.getSize();

   if(size != numberOfInputs)
   {
      std::cerr << std::endl
                << "Flood Error: Perceptron class." << std::endl 
                << "double calculateOutputSignal(Vector<double>) method." << std::endl
                << "Size must be equal to number of inputs." << std::endl
                << std::endl;

      exit(1);
   }

   #endif

   // Calculate output signal

   double netInputSignal = calculateNetInputSignal(inputSignal);

   double outputSignal = calculateOutputSignal(netInputSignal);

   return(outputSignal);  
}


// double calculateOutputSignalDerivative(double) method

/// This method returns the activation derivative of the neuron for a net input signal.
/// The output signal derivative depends on the activation function used.
/// 
/// @param netInputSignal Net input signal to the neuron.
/// 
/// @see calculateOutputSignal(double).
/// @see calculateOutputSignal(Vector<double>).
/// @see calculateOutputSignalDerivative(Vector<double>).
/// @see calculateOutputSignalSecondDerivative(double).
/// @see calculateOutputSignalSecondDerivative(Vector<double>).

double Perceptron::calculateOutputSignalDerivative(double netInputSignal)
{
   double outputSignalDerivative = 0.0;

   switch(activationFunction)   
   {
      case Logistic:
      {
         outputSignalDerivative = exp(netInputSignal)/pow(1.0 + exp(netInputSignal), 2);
      }
      break;
                                     
      case HyperbolicTangent:
      {
         outputSignalDerivative = 1.0 - pow(tanh(netInputSignal), 2);
      }
      break;

      case Linear:
      {
         outputSignalDerivative = 1.0;
      }
      break;

      default:
         std::cerr << std::endl
                   << "Flood Error: Perceptron class." << std::endl 
                   << "calculateOutputSignalDerivative(Vector<double>) method." << std::endl
                   << "Unknown activation function." << std::endl
                   << std::endl;

         exit(1);
      break;
   }

   return(outputSignalDerivative);
}


// double calculateOutputSignalDerivative(Vector<double>) method

/// This method returns the activation derivative of the neuron for a set of input signals.
/// The output signal derivative depends on the activation function used.
/// 
/// @param inputSignal Set of input signals to the neuron.
/// 
/// @see calculateOutputSignal(double).
/// @see calculateOutputSignal(Vector<double>).
/// @see calculateOutputSignalDerivative(double).
/// @see calculateOutputSignalSecondDerivative(double).
/// @see calculateOutputSignalSecondDerivative(Vector<double>).

double Perceptron::calculateOutputSignalDerivative(Vector<double> inputSignal)
{
   // Control sentence (if debug)

   #ifndef NDEBUG 

   int size = inputSignal.getSize();

   if(size != numberOfInputs)
   {
      std::cerr << std::endl
                << "Flood Error: Perceptron class." << std::endl 
                << "double calculateOutputSignalDerivative(Vector<double>) method." << std::endl
                << "Size must be equal to number of inputs." << std::endl
                << std::endl;

      exit(1);
   }

   #endif

   // Calculate output signal derivative 

   double netInputSignal = calculateNetInputSignal(inputSignal);

   double outputSignalDerivative = calculateOutputSignalDerivative(netInputSignal);

   return(outputSignalDerivative);
}


// double calculateOutputSignalSecondDerivative(double) method

/// This method returns the activation second derivative of the neuron for a net input signal.
/// The second derivative of the output signal depends on the activation function used.
/// 
/// @param netInputSignal Net input signal to the neuron.
/// 
/// @see calculateOutputSignal(double).
/// @see calculateOutputSignal(Vector<double>).
/// @see calculateOutputSignalDerivative(double).
/// @see calculateOutputSignalDerivative(Vector<double>).
/// @see calculateOutputSignalSecondDerivative(Vector<double>).

double Perceptron::calculateOutputSignalSecondDerivative(double netInputSignal)
{
   double outputSignalSecondDerivative;

   switch(activationFunction)   
   {
      case Logistic:
      {
         outputSignalSecondDerivative 
         = -exp(netInputSignal)*(exp(netInputSignal) - 1.0)/pow(exp(netInputSignal + 1), 3);
      }
      break;
                                     
      case HyperbolicTangent:
      {
         outputSignalSecondDerivative = -2.0*tanh(netInputSignal)*(1.0 - pow(tanh(netInputSignal),2));
      }
      break;

      case Linear:
      {
         outputSignalSecondDerivative = 0.0;        
      }
      break;

      default:
         std::cerr << std::endl
                   << "Flood Error: Perceptron class." << std::endl 
                   << "calculateOutputSignalSecondDerivative(Vector<double>) method." << std::endl
                   << "Unknown activation function." << std::endl
                   << std::endl;

         exit(1);
      break;
   }

   return(outputSignalSecondDerivative);
}


// double calculateOutputSignalSecondDerivative(Vector<double>) method

/// This method returns the activation second derivative of the neuron for a set of input signals.
/// The second derivative of the output signal depends on the activation function used.
/// 
/// @param inputSignal Set of input signals to the neuron.
/// 
/// @see calculateOutputSignal(double).
/// @see calculateOutputSignal(Vector<double>).
/// @see calculateOutputSignalDerivative(double).
/// @see calculateOutputSignalDerivative(Vector<double>).
/// @see calculateOutputSignalSecondDerivative(double).

double Perceptron::calculateOutputSignalSecondDerivative(Vector<double> inputSignal)
{
   // Control sentence (if debug)

   #ifndef NDEBUG 

   int size = inputSignal.getSize();

   if(size != numberOfInputs)
   {
      std::cerr << std::endl
                << "Flood Error: Perceptron class." << std::endl 
                << "double calculateOutputSignalSecondDerivative(Vector<double>) method." << std::endl
                << "Size must be equal to number of inputs." << std::endl
                << std::endl;

      exit(1);
   }

   #endif

   // Calculate output signal second derivative

   double netInputSignal = calculateNetInputSignal(inputSignal);

   double outputSignalSecondDerivative = calculateOutputSignalSecondDerivative(netInputSignal);

   return(outputSignalSecondDerivative);
}


// void print(void) method

/// This method prints to the screen the number of inputs, the bias and the synaptic weight values of the 
/// perceptron.

void Perceptron::print(void)
{
    // Print header

   std::cout << std::endl
             << "Flood Neural Network. Perceptron Object." << std::endl;

   // Print number of inputs

   std::cout << "Number of inputs: " << std::endl
             << numberOfInputs << std::endl;


   // Print activation function

   std::cout << "Activation function:" << std::endl;

   switch(activationFunction)
   {
      case Logistic:   
      {
         std::cout << "Logistic" << std::endl;
      }
      break;

      case HyperbolicTangent:   
      {
         std::cout << "Hyperbolic tangent" << std::endl;
      }
      break;

      case Linear:   
      {
         std::cout << "Linear" << std::endl;
      }
      break;

      default:
         std::cerr << std::endl
                   << "Flood Error: Perceptron class." << std::endl 
                   << "void print(void) method." << std::endl
                   << "Unknown activation function." << std::endl
                   << std::endl;

         exit(1);
      break;

   }

   // Print bias

   std::cout << "Bias: " << std::endl
             << bias << std::endl;

   // Print synaptic weights

   std::cout << "Synaptic weights: " << std::endl
             << synapticWeights << std::endl;   

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
