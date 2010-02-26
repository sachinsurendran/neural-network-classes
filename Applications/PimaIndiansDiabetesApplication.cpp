/****************************************************************************************************************/
/*                                                                                                              */ 
/*   Flood: An Open Source Neural Networks C++ Library                                                          */
/*   www.cimne.com/flood                                                                                        */
/*                                                                                                              */
/*   P I M A   I N D I A N S   D I A B E T E S   A P P L I C A T I O N                                          */
/*                                                                                                              */
/*   Roberto Lopez                                                                                              */ 
/*   International Center for Numerical Methods in Engineering (CIMNE)                                          */
/*   Technical University of Catalonia (UPC)                                                                    */
/*   Barcelona, Spain                                                                                           */
/*   E-mail: rlopez@cimne.upc.edu                                                                               */ 
/*                                                                                                              */  
/****************************************************************************************************************/

/// In this main function a neural network is trained to predict whether an individual of Pima Indian heritage 
/// has diabetes from personal characteristics and physical measurements.

// System includes

#include <iostream>
#include <time.h>

// Utilities includes

#include "../Flood/Utilities/Vector.h"
#include "../Flood/Utilities/Matrix.h"

#include "../Flood/Utilities/InputTargetDataSet.h"
#include "../Flood/Utilities/CorrectPredictionsAnalysis.h"

// Network architecture includes

#include "../Flood/MultilayerPerceptron/MultilayerPerceptron.h"

// Objective functional includes

#include "../Flood/ObjectiveFunctional/MeanSquaredError.h"

// Training algorithm includes

#include "../Flood/TrainingAlgorithm/QuasiNewtonMethod.h"

using namespace Flood;

int main(void)
{
   std::cout << std::endl
             << "Flood Neural Network. Pima Indians Diabetes Application." 
             << std::endl;	

   srand((unsigned)time(NULL));

   // Construct an input-target data set object for training

   InputTargetDataSet trainingDataSet;

   // Load training data set object from file

   trainingDataSet.load("../Data/PimaIndiansDiabetes/TrainingDataSet.dat");

   // Get information

   Vector< Vector<std::string> > information = trainingDataSet.getAllInformation();

   // Calculate all statistics

   Vector< Vector<double> > statistics = trainingDataSet.calculateAllStatistics();

   // Preprocess training data for zero mean and unity standard deviation 

   trainingDataSet.preprocessMeanAndStandardDeviation();

   // Construct a multilayer perceptron object with 8 inputs, 
   // 6 sigmoid neurons in the hidden layer and 1 linear output neuron

   MultilayerPerceptron multilayerPerceptron(8, 6, 1);

   // Construct a mean squared error object

   MeanSquaredError meanSquaredError(&multilayerPerceptron, &trainingDataSet);

   // Construct a quasi-Newton method object

   QuasiNewtonMethod quasiNewtonMethod(&meanSquaredError);

   quasiNewtonMethod.setReserveEvaluationHistory(true);
   quasiNewtonMethod.setReserveGradientNormHistory(true);

   // Set the stopping criteria for the quasi-Newton method training algorithm

   quasiNewtonMethod.setEvaluationGoal(0.0);
   quasiNewtonMethod.setGradientNormGoal(0.01);
   quasiNewtonMethod.setMaximumTime(100.0);
   quasiNewtonMethod.setMaximumNumberOfEpochs(100);

   quasiNewtonMethod.setDisplayPeriod(10);

   // Train neural network

   quasiNewtonMethod.train();

   // Set all information in multilayer perceptron

   multilayerPerceptron.setAllInformation(information); 


   // Set all statistics

   multilayerPerceptron.setAllStatistics(statistics); 


   // Save training history

   quasiNewtonMethod.saveTrainingHistory("../Data/PimaIndiansDiabetes/TrainingHistory.dat");

   multilayerPerceptron
   .setPreAndPostProcessingMethod(MultilayerPerceptron::MeanAndStandardDeviation);

   // Save multilayer perceptron

   multilayerPerceptron.save("../Data/PimaIndiansDiabetes/MultilayerPerceptron.dat");

   // Construct and input-target data set for validation

   InputTargetDataSet validationDataSet;

   validationDataSet.load("../Data/PimaIndiansDiabetes/ValidationDataSet.dat");

   CorrectPredictionsAnalysis correctPredictionsAnalysis(&multilayerPerceptron, &validationDataSet);

   double correctPredictionsRatio = correctPredictionsAnalysis.calculateCorrectPredictionsRatio();

   // Print results

   std::cout << std::endl 
             << "Correct predictions ratio: " << correctPredictionsRatio << std::endl;

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
