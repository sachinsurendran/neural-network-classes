/****************************************************************************************************************/
/*                                                                                                              */ 
/*   Flood: An Open Source Neural Networks C++ Library                                                          */
/*   www.cimne.com/flood                                                                                        */
/*                                                                                                              */
/*   M U L T I L A Y E R   P E R C E P T R O N   A P P L I C A T I O N                                          */
/*                                                                                                              */
/*   Roberto Lopez                                                                                              */ 
/*   International Center for Numerical Methods in Engineering (CIMNE)                                          */
/*   Technical University of Catalonia (UPC)                                                                    */
/*   Barcelona, Spain                                                                                           */
/*   E-mail: rlopez@cimne.upc.edu                                                                               */ 
/*                                                                                                              */  
/****************************************************************************************************************/

/// This application is an usage example of the MultilayerPerceptron class in Flood.

// System includes

#include <iostream>
#include <time.h>

// Utilities includes

#include "../Flood/Utilities/Vector.h"
#include "../Flood/Utilities/Matrix.h"

// Network architecture includes

#include "../Flood/MultilayerPerceptron/MultilayerPerceptron.h"

using namespace Flood;

int main(void)
{
   std::cout << std::endl
             << "Flood Neural Network. MultilayerPerceptron Application." 
             << std::endl;	

//   srand((unsigned)time(NULL));

   // Construct a multilayer perceptron object with 3 inputs, two hidden layers with 6 and 5 
   // neurons and 1 output neuron

   Vector<int> numbersOfHiddenNeurons(2,1);

   MultilayerPerceptron multilayerPerceptron(3, numbersOfHiddenNeurons, 2);

   // Set new hidden and output layer activation functions

   int numberOfHiddenLayers = numbersOfHiddenNeurons.getSize();

   Vector<Perceptron::ActivationFunction> hiddenLayersActivationFunction(numberOfHiddenLayers, Perceptron::Logistic);   

   multilayerPerceptron.setHiddenLayersActivationFunction(hiddenLayersActivationFunction);

   multilayerPerceptron.setOutputLayerActivationFunction(Perceptron::Linear);

   // Use the mean and standard deviation pre and post-processing method

   multilayerPerceptron.setPreAndPostProcessingMethod(MultilayerPerceptron::None);

   // Init biases and synaptic weights with values comprised between -0.5 and 0.5

   multilayerPerceptron.initNeuralParametersUniform(-0.5, 0.5);

   // Print multilayer perceptron to screen

   multilayerPerceptron.print();

   // Save neural network to data file

   multilayerPerceptron.save("../Data/MultilayerPerceptron/MultilayerPerceptron.dat");

   // Set an input to the network
   
   Vector<double> input(3, 1.0);
   input[0] = -1.0;
   input[1] = 0.0;
   input[2] = 1.0;


   std::cout << std::endl
             << "Input: " << std::endl
             << input << std::endl;

   // Calculate output from the network

   Vector<double> output = multilayerPerceptron.calculateOutput(input);

   std::cout << std::endl
             << "Output: " << std::endl
             << output << std::endl;

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
