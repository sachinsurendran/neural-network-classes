/****************************************************************************************************************/
/*                                                                                                              */ 
/*   Flood: An Open Source Neural Networks C++ Library                                                          */
/*   www.cimne.com/flood                                                                                        */
/*                                                                                                              */
/*   E V O L U T I O N A R Y   A L G O R I T H M   A P P L I C A T I O N                                        */
/*                                                                                                              */
/*   Roberto Lopez                                                                                              */ 
/*   International Center for Numerical Methods in Engineering (CIMNE)                                          */
/*   Technical University of Catalonia (UPC)                                                                    */
/*   Barcelona, Spain                                                                                           */
/*   E-mail: rlopez@cimne.upc.edu                                                                               */ 
/*                                                                                                              */  
/****************************************************************************************************************/

/// This application is an usage example of the EvolutionaryAlgorithm class in Flood.

// System includes

#include <iostream>
#include <time.h>

// Utilities includes

#include "../Flood/Utilities/Vector.h"
#include "../Flood/Utilities/Matrix.h"

#include "../Flood/Utilities/InputTargetDataSet.h"

// Multilayer perceptron includes

#include "../Flood/MultilayerPerceptron/MultilayerPerceptron.h"

// Objective functional includes

#include "../Flood/ObjectiveFunctional/TennixTrainer.h"

// Training algorithm includes

#include "../Flood/TrainingAlgorithm/EvolutionaryAlgorithm.h"

using namespace Flood;

int main(void)
{
   std::cout << std::endl
             << "Flood Neural Network. EvolutionaryAlgorithm Application." 
             << std::endl;	

   srand((unsigned)time(NULL));

   // Input-target data set object

   //InputTargetDataSet inputTargetDataSet;

//   inputTargetDataSet.load("../Data/EvolutionaryAlgorithm/InputTargetDataSet.dat");

   //inputTargetDataSet.load("../Data/EvolutionaryAlgorithm/XOR.dat");

   // Multilayer perceptron object
   //
#define WIDTH_OF_HIDDEN_LAYER 10
   Vector<int> numbersOfHiddenNeurons (6);
   numbersOfHiddenNeurons[0] = WIDTH_OF_HIDDEN_LAYER;
   numbersOfHiddenNeurons[1] = WIDTH_OF_HIDDEN_LAYER;
   numbersOfHiddenNeurons[2] = WIDTH_OF_HIDDEN_LAYER;
   numbersOfHiddenNeurons[3] = WIDTH_OF_HIDDEN_LAYER;
   numbersOfHiddenNeurons[4] = WIDTH_OF_HIDDEN_LAYER;
   numbersOfHiddenNeurons[5] = WIDTH_OF_HIDDEN_LAYER;
//   numbersOfHiddenNeurons[6] = WIDTH_OF_HIDDEN_LAYER;

   MultilayerPerceptron multilayerPerceptron(6, numbersOfHiddenNeurons, 3);

   multilayerPerceptron.setOutputLayerActivationFunction(Perceptron::HyperbolicTangent);

   // Mean squared error object

   TennixTrainer
   tennixTrainer(&multilayerPerceptron/*, &inputTargetDataSet*/);// Now just need a NN, no inputs, coz it comes from tennix server

   // Evolutionary algorithm object

   EvolutionaryAlgorithm evolutionaryAlgorithm(&tennixTrainer);

   evolutionaryAlgorithm.setReservePopulationHistory(true);
   evolutionaryAlgorithm.setReserveMeanNormHistory(true);
   evolutionaryAlgorithm.setReserveStandardDeviationNormHistory(true);
   evolutionaryAlgorithm.setReserveBestNormHistory(true);
   evolutionaryAlgorithm.setReserveMeanEvaluationHistory(true);
   evolutionaryAlgorithm.setReserveStandardDeviationEvaluationHistory(true);
   evolutionaryAlgorithm.setReserveBestEvaluationHistory(true);

   evolutionaryAlgorithm.setPopulationSize(500);
   evolutionaryAlgorithm.initPopulationNormal(0.0,1.0);
   //evolutionaryAlgorithm.initPopulationUniform(-3,3);

   // Set stopping criteria

   evolutionaryAlgorithm.setEvaluationGoal(0.1);
   evolutionaryAlgorithm.setMaximumTime(100000.0);
   evolutionaryAlgorithm.setMaximumNumberOfGenerations(1000000);

   // Set user stuff

   evolutionaryAlgorithm.setDisplayPeriod(1);

   // Train neural network

   evolutionaryAlgorithm.train();

   // Save all training history

   //evolutionaryAlgorithm.saveTrainingHistory("../Data/EvolutionaryAlgorithm/TrainingHistory.dat");



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
