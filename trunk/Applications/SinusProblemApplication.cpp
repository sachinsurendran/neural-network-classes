/****************************************************************************************************************/
/*                                                                                                              */ 
/*   Flood: An Open Source Neural Networks C++ Library                                                          */
/*   www.cimne.com/flood                                                                                        */
/*                                                                                                              */
/*   S I N U S   P R O B L E M   A P P L I C A T I O N                                                          */
/*                                                                                                              */
/*   Roberto Lopez                                                                                              */ 
/*   International Center for Numerical Methods in Engineering (CIMNE)                                          */
/*   Technical University of Catalonia (UPC)                                                                    */
/*   Barcelona, Spain                                                                                           */
/*   E-mail: rlopez@cimne.upc.edu                                                                               */ 
/*                                                                                                              */  
/****************************************************************************************************************/

/// In this main function a neural network is trained to solve a simple function regression problem with one input
/// and one output variables. 
/// The applied training algorithm is the quasi-Newton method. 


// System includes

#include <iostream>
#include <fstream>
#include <time.h>

// Utilities includes

#include "../Flood/Utilities/Vector.h"
#include "../Flood/Utilities/Matrix.h"

#include "../Flood/Utilities/InputTargetDataSet.h"

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
             << "Flood Neural Network. Sinus Problem Application." 
             << std::endl;	

//   srand((unsigned)time(NULL));

   // Construct an input-target data set object for training

   InputTargetDataSet itds;

   // Load training data set object from file

   itds.load("../Data/SinusProblem/InputTargetDataSet.dat");

   // Get all information 

   Vector< Vector<std::string> > information = itds.getAllInformation();

   // Calculate all statistics

   Vector< Vector<double> > statistics = itds.calculateAllStatistics();
 
   Matrix<double> meanAndStandardDeviationOfInputData = itds.calculateMeanAndStandardDeviationOfInputData();
   Matrix<double> meanAndStandardDeviationOfTargetData = itds.calculateMeanAndStandardDeviationOfTargetData();

   // Preprocess training data for zero mean and unity standard deviation 

   itds.preprocessMeanAndStandardDeviation();

   // Construct a multilayer perceptron object with 8 inputs, 6 sigmoid neurons in the hidden layer and 1 linear 
   // output neuron

   MultilayerPerceptron mlp(1, 3, 1);

   mlp.initFreeParametersNormal(0.0,1.0);

   // Construct a mean squared error object

   MeanSquaredError mse(&mlp, &itds);

   // Construct a quasi-Newton method object

   QuasiNewtonMethod qnm(&mse);

   // Set the stopping criteria for the quasi-Newton method training algorithm

   qnm.setEvaluationGoal(0.0);
   qnm.setGradientNormGoal(0.0);
   qnm.setMaximumTime(1000.0);
   qnm.setMaximumNumberOfEpochs(1000);

   qnm.setDisplayPeriod(1000);

   // Train neural network

   qnm.train();

   // Set all information

   mlp.setAllInformation(information); 

   // Set all statistics

   mlp.setAllStatistics(statistics); 

   // Save training history

   itds.postprocessMeanAndStandardDeviation(meanAndStandardDeviationOfInputData, meanAndStandardDeviationOfTargetData);

   mlp.setPreAndPostProcessingMethod(MultilayerPerceptron::MeanAndStandardDeviation);

   // Save multilayer perceptron

   mlp.save("../Data/SinusProblem/MultilayerPerceptron.dat");
   mlp.saveExpression("../Data/SinusProblem/Expression.dat");

   mse.saveInputTargetAndOutput("../Data/SinusProblem/InputTargetOutput.dat");

   // Save results to file

   char* filename = "../Data/SinusProblem/Results.dat";

   std::fstream file; 
   file.open(filename, std::ios::out);

   int n = 101;

   Vector<double> sinus(n, 0.0);
   Vector<double> model(n, 0.0);

   Vector<double> input(1, 0.0);
   Vector<double> output(1, 0.0);

   double pi = 4.0*atan(1.0);

   for(int i = 0; i < n; i++)
   {      
      input[0] = 1.0*i/(n-1.0);
	  output = mlp.calculateOutput(input);

      sinus[i] = 0.5 + 0.4*sin(2.0*pi*input[0]);
      model[i] = output[0];  

	  file << input[0] << " " << sinus[i] << " " << model[i] << std::endl;
   }

   file.close();

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
