/****************************************************************************************************************/
/*                                                                                                              */ 
/*   Flood: An Open Source Neural Networks C++ Library                                                          */
/*   www.cimne.com/flood                                                                                        */
/*                                                                                                              */
/*   I S O P E R I M E T R I C   P R O B L E M   A P P L I C A T I O N                                          */
/*                                                                                                              */
/*   Roberto Lopez                                                                                              */ 
/*   International Center for Numerical Methods in Engineering (CIMNE)                                          */
/*   Technical University of Catalonia (UPC)                                                                    */
/*   Barcelona, Spain                                                                                           */
/*   E-mail: rlopez@cimne.upc.edu                                                                               */ 
/*                                                                                                              */  
/****************************************************************************************************************/

/// This main function can be used as a template for solving the isoperimetric problem by means of a multilayer
/// perceptron. 
/// It uses the quasi-Newton method training algorithm.

// System includes

#include <iostream>
#include <time.h>

// Mutilayer perceptron includes

#include "../Flood/MultilayerPerceptron/MultilayerPerceptron.h"

// Objective functional includes

#include "../Flood/ObjectiveFunctional/IsoperimetricProblem.h"

// Training algorithm includes

#include "../Flood/TrainingAlgorithm/QuasiNewtonMethod.h"


using namespace Flood;

int main(void)
{
   std::cout << std::endl 
             << "Flood Neural Network. Isoperimetric Problem Application." 
             << std::endl;	

//   srand((unsigned)time(NULL));

   // Multilayer perceptron object

   MultilayerPerceptron mlp(1, 3, 2);

   Vector<std::string> nameOfInputVariables(1, "t");
   Vector<std::string> nameOfOutputVariables(2);
   nameOfOutputVariables[0] = "x";
   nameOfOutputVariables[1] = "y";

   mlp.setNameOfInputVariables(nameOfInputVariables);
   mlp.setNameOfOutputVariables(nameOfOutputVariables);

   Vector<double> minimumOfInputVariables(1, 0.0);
   Vector<double> maximumOfInputVariables(1, 1.0);
   Vector<double> minimumOfOutputVariables(2,-1.0);
   Vector<double> maximumOfOutputVariables(2,1.0);

   mlp.setMinimumOfInputVariables(minimumOfInputVariables);
   mlp.setMaximumOfInputVariables(maximumOfInputVariables);
   mlp.setMinimumOfOutputVariables(minimumOfOutputVariables);
   mlp.setMaximumOfOutputVariables(maximumOfOutputVariables);

   // Isoperimetric problem object

   IsoperimetricProblem ip(&mlp);

   ip.setPerimeterErrorWeight(100.0);
   ip.setAreaWeight(1.0);

   ip.saveResults("../Data/IsoperimetricProblem/InitialGuessIsoperimetricProblem.dat");

   // Quasi-Newton method object

   QuasiNewtonMethod qnm(&ip);

   qnm.setReserveEvaluationHistory(true);
   qnm.setReserveGradientNormHistory(true);

   qnm.setMaximumNumberOfEpochs(5000); 
   qnm.setDisplayPeriod(10); 

   // Train neural network

   int numberOfEvaluations = 0;

   ip.setNumberOfEvaluations(numberOfEvaluations);

   qnm.train(); 

   numberOfEvaluations = ip.getNumberOfEvaluations();

   std::cout << std::endl
             << "Number of evaluations: " 
             << numberOfEvaluations <<  std::endl;	

   // Save neural network to file

   mlp.save("../Data/IsoperimetricProblem/MultilayerPerceptronIsoperimetricProblem.dat");

   mlp.saveExpression("../Data/IsoperimetricProblem/ExpressionIsoperimetricProblem.dat");

   // Save results to file

   ip.saveResults("../Data/IsoperimetricProblem/ResultsIsoperimetricProblem.dat");

   // Save training history to file

   qnm.saveTrainingHistory("../Data/IsoperimetricProblem/TrainingHistoryIsoperimetricProblem.dat");

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
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
