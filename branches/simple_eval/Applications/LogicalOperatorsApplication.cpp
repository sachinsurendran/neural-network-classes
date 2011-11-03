/****************************************************************************************************************/
/*                                                                                                              */ 
/*   Flood: An Open Source Neural Networks C++ Library                                                          */
/*   www.cimne.com/flood                                                                                        */
/*                                                                                                              */
/*   L O G I C A L   O P E R A T O R S   A P P L I C A T I O N                                                  */
/*                                                                                                              */
/*   Roberto Lopez                                                                                              */ 
/*   International Center for Numerical Methods in Engineering (CIMNE)                                          */
/*   Technical University of Catalonia (UPC)                                                                    */
/*   Barcelona, Spain                                                                                           */
/*   E-mail: rlopez@cimne.upc.edu                                                                               */ 
/*                                                                                                              */  
/****************************************************************************************************************/

/// In this main function a neural network is trained to learn the AND, OR, NAND, NOR, XOR and XNOR logical 
/// operations. 

// System includes

#include <iostream>
#include <math.h>
#include <time.h>

// Utilities includes

#include "../Flood/Utilities/Vector.h"
#include "../Flood/Utilities/Matrix.h"

#include "../Flood/Utilities/InputTargetDataSet.h"
#include "../Flood/Utilities/CorrectPredictionsAnalysis.h"

// Multilayer perceptron includes

#include "../Flood/MultilayerPerceptron/MultilayerPerceptron.h"

// Objective functional includes

#include "../Flood/ObjectiveFunctional/MeanSquaredError.h"

// Training algorithm includes

#include "../Flood/TrainingAlgorithm/QuasiNewtonMethod.h"

using namespace Flood;

int main(void)
{
   std::cout << std::endl
             << "Flood Neural Network. Logical Operations Application." 
			 << std::endl;	
  
   srand((unsigned)time(NULL));

   // Input-target data set object

   InputTargetDataSet inputTargetDataSet;
   
   inputTargetDataSet.load("../Data/LogicalOperations/InputTargetDataSet.dat");

   // Get all information from input-target data set

   Vector< Vector<std::string> > information = inputTargetDataSet.getAllInformation();

   // Calculate all statistics of input and target data

   Vector< Vector<double> > statistics = inputTargetDataSet.calculateAllStatistics();

   // Preprocess input and target data for minimum -1 and maximum +1

   inputTargetDataSet.preprocessMinimumAndMaximum();

   // Construct a multilayer perceptron object

   MultilayerPerceptron multilayerPerceptron(2,6,6);

   // Mean squared error object

   MeanSquaredError meanSquaredError(&multilayerPerceptron, &inputTargetDataSet);

   // Quasi-Newton method object

   QuasiNewtonMethod quasiNewtonMethod(&meanSquaredError);

   quasiNewtonMethod.setGradientNormGoal(1.0e-6);

   quasiNewtonMethod.train();

   // Set all information from input-target data set in multilayer perceptron
   
   multilayerPerceptron.setAllInformation(information); 

   // Set all statistics of input-target data set in multilayer perceptron
   
   multilayerPerceptron.setAllStatistics(statistics); 

   multilayerPerceptron.setPreAndPostProcessingMethod(MultilayerPerceptron::MinimumAndMaximum);

   // Print results to screen

   std::cout << std::endl 
             << "X " 
             << "Y " 
             << "AND "  
             << "OR " 
             << "NAND " 
             << "NOR " 
             << "XOR " 
             << "XNOR " 
             << std::endl;

   Vector<double> input(2, 0.0);
   Vector<double> output(6, 0.0);

   input[0] = 1.0;
   input[1] = 1.0;

   output = multilayerPerceptron.calculateOutput(input);

   std::cout << (int)input[0] << " " 
             << (int)input[1] << " " 
             << (int)(output[0] + 0.5) << "   "  
             << (int)(output[1] + 0.5) << "  " 
             << (int)(output[2] + 0.5) << "    " 
             << (int)(output[3] + 0.5) << "   " 
             << (int)(output[4] + 0.5) << "   " 
             << (int)(output[5] + 0.5) << std::endl;

   input[0] = 1.0;
   input[1] = 0.0;

   output = multilayerPerceptron.calculateOutput(input);

   std::cout << (int)input[0] << " " 
             << (int)input[1] << " " 
             << (int)(output[0] + 0.5) << "   "  
             << (int)(output[1] + 0.5) << "  " 
             << (int)(output[2] + 0.5) << "    " 
             << (int)(output[3] + 0.5) << "   " 
             << (int)(output[4] + 0.5) << "   " 
             << (int)(output[5] + 0.5) << std::endl;

   input[0] = 0.0;
   input[1] = 1.0;

   output = multilayerPerceptron.calculateOutput(input);

   std::cout << (int)input[0] << " " 
             << (int)input[1] << " " 
             << (int)(output[0] + 0.5) << "   "  
             << (int)(output[1] + 0.5) << "  " 
             << (int)(output[2] + 0.5) << "    " 
             << (int)(output[3] + 0.5) << "   " 
             << (int)(output[4] + 0.5) << "   " 
             << (int)(output[5] + 0.5) << std::endl;

   input[0] = 0.0;
   input[1] = 0.0;

   output = multilayerPerceptron.calculateOutput(input);

   std::cout << (int)input[0] << " " 
             << (int)input[1] << " " 
             << (int)(output[0] + 0.5) << "   "  
             << (int)(output[1] + 0.5) << "  " 
             << (int)(output[2] + 0.5) << "    " 
             << (int)(output[3] + 0.5) << "   " 
             << (int)(output[4] + 0.5) << "   " 
             << (int)(output[5] + 0.5) << std::endl;


   // Load original input-target data set

   inputTargetDataSet.load("../Data/LogicalOperations/InputTargetDataSet.dat");

   CorrectPredictionsAnalysis correctPredictionsAnalysis(&multilayerPerceptron, &inputTargetDataSet);

   double correctPredictionsRatio = correctPredictionsAnalysis.calculateCorrectPredictionsRatio();

   std::cout << std::endl
    	     << "Ratio of correct predictions: " << correctPredictionsRatio << std::endl;


   // Save multilayer perceptron object

   multilayerPerceptron.save("../Data/LogicalOperations/MultilayerPerceptron.dat");

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
