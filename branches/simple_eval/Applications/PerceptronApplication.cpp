/****************************************************************************************************************/
/*                                                                                                              */ 
/*   Flood: An Open Source Neural Networks C++ Library                                                          */
/*   www.cimne.com/flood                                                                                        */
/*                                                                                                              */
/*   P E R C E P T R O N   A P P L I C A T I O N                                                                */
/*                                                                                                              */
/*   Roberto Lopez                                                                                              */ 
/*   International Center for Numerical Methods in Engineering (CIMNE)                                          */
/*   Technical University of Catalonia (UPC)                                                                    */
/*   Barcelona, Spain                                                                                           */
/*   E-mail: rlopez@cimne.upc.edu                                                                               */ 
/*                                                                                                              */  
/****************************************************************************************************************/

/// This application is an usage example of the Perceptron class in Flood.

// System includes

#include <iostream>
#include <time.h>

// Utilities includes

#include "../Flood/Utilities/Vector.h"
#include "../Flood/Utilities/Matrix.h"

// Neuron model includes

#include "../Flood/Perceptron/Perceptron.h"

using namespace Flood;

int main(void)
{
   std::cout << std::endl
             << "Flood NeuralNetwork. Perceptron Application." 
             << std::endl;	

//   srand((unsigned)time(NULL));

   // Construct a perceptron object with 3 inputs

   Perceptron perceptron(3);

   // Init neuron's bias and synaptic weights with random values 
   // chosen from a normal distribution with mean 0 and standard deviation 1

   perceptron.initBiasNormal(0.0,1.0);
   perceptron.initSynapticWeightsNormal(0.0,1.0); 

   // Set up a vector of input signals to the neuron 

   Vector<double> inputSignals(3);
   inputSignals[0] = -0.1; 
   inputSignals[1] = 0.0; 
   inputSignals[2] = 0.2; 

   std::cout << std::endl
             << "Input signals:" << std::endl
			 << inputSignals << std::endl;


   // Get the net input signal to the neuron

   double netInputSignal = perceptron.calculateNetInputSignal(inputSignals);

   std::cout << std::endl
             << "Net input signal:" << std::endl
             << netInputSignal << std::endl;

   // Set up a logistic activation function

   std::cout << std::endl
	         << "Logistic activation function" << std::endl;

   perceptron.setActivationFunction(Perceptron::Logistic);

   // Get the output signal from the neuron

   double outputSignalLogistic = perceptron.calculateOutputSignal(inputSignals);

   std::cout << std::endl
             << "Output signal logistic:" << std::endl
             << outputSignalLogistic << std::endl;

   // Get the activation derivative of the neuron

   double outputSignalDerivativeLogistic = perceptron.calculateOutputSignalDerivative(netInputSignal);

   std::cout << std::endl
             << "Output signal derivative logistic:" << std::endl
             << outputSignalDerivativeLogistic << std::endl;

   // Get the activation second derivative of the neuron

   double outputSignalSecondDerivativeLogistic = perceptron.calculateOutputSignalSecondDerivative(netInputSignal);

   std::cout << std::endl
             << "Output signal second derivative logistic:" << std::endl
             << outputSignalSecondDerivativeLogistic << std::endl;

   // Set up an hyperbolic tangent activation function

   std::cout << std::endl
	         << "Hyperbolic tangent activation function" << std::endl;

   perceptron.setActivationFunction(Perceptron::HyperbolicTangent);

   // Get the output signal from the neuron

   double outputSignalHyperbolicTangent = perceptron.calculateOutputSignal(inputSignals);

   std::cout << std::endl
             << "Output signal hyperbolic tangent:" << std::endl
             << outputSignalHyperbolicTangent << std::endl;

   // Get the activation derivative of the neuron

   double outputSignalDerivativeHyperbolicTangent = perceptron.calculateOutputSignalDerivative(netInputSignal);

   std::cout << std::endl
             << "Output signal derivative hyperbolic tangent:" << std::endl
             << outputSignalDerivativeHyperbolicTangent << std::endl;

   // Get the activation second derivative of the neuron

   double outputSignalSecondDerivativeHyperbolicTangent = perceptron.calculateOutputSignalSecondDerivative(netInputSignal);

   std::cout << std::endl
             << "Output signal second derivative hyperbolic tangent:" << std::endl
             << outputSignalSecondDerivativeHyperbolicTangent << std::endl;

   // Set up a linear activation function

   std::cout << std::endl
	         << "Linear activation function" << std::endl;

   perceptron.setActivationFunction(Perceptron::Linear);

   // Get the output signal from the neuron

   double outputSignalLinear = perceptron.calculateOutputSignal(inputSignals);

   std::cout << std::endl
             << "Output signal linear:" << std::endl
             << outputSignalLinear << std::endl;

   // Get the activation derivative of the neuron

   double outputSignalDerivativeLinear = perceptron.calculateOutputSignalDerivative(netInputSignal);

   std::cout << std::endl
             << "Output signal derivative linear:" << std::endl
             << outputSignalDerivativeLinear << std::endl;

   // Get the activation second derivative of the neuron

   double outputSignalSecondDerivativeLinear = perceptron.calculateOutputSignalSecondDerivative(netInputSignal);

   std::cout << std::endl
             << "Output signal second derivative linear:" << std::endl
             << outputSignalSecondDerivativeLinear << std::endl;

   std::cout << std::endl;

   return(0);
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
