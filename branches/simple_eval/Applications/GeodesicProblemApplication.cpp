/****************************************************************************************************************/
/*                                                                                                              */
/*   Flood: An Open Source Neural Networks C++ Library                                                          */
/*   www.cimne.com/flood                                                                                        */
/*                                                                                                              */
/*   G E O D E S I C   P R O B L E M   A P P L I C A T I O N                                                    */
/*                                                                                                              */
/*   Roberto Lopez                                                                                              */ 
/*   International Center for Numerical Methods in Engineering (CIMNE)                                          */
/*   Technical University of Catalonia (UPC)                                                                    */
/*   Barcelona, Spain                                                                                           */
/*   E-mail: rlopez@cimne.upc.edu                                                                               */ 
/*                                                                                                              */
/****************************************************************************************************************/

/// This main function can be used as a template for solving the geodesic problem by means of a multilayer 
/// perceptron. 
/// It uses the quasi-Newton method for training. 

// System includes

#include <iostream>
#include <time.h>

// Network architecture includes

#include "../Flood/MultilayerPerceptron/MultilayerPerceptron.h"

// Objective functional includes

#include "../Flood/ObjectiveFunctional/GeodesicProblem.h"

// Training algorithm includes

#include "../Flood/TrainingAlgorithm/QuasiNewtonMethod.h"


using namespace Flood;

int main(void)
{
   std::cout << std::endl
             << "Flood Neural Network. Geodesic Problem Application." 
             << std::endl;	

   srand((unsigned)time(NULL));

   // Multilayer perceptron object

   MultilayerPerceptron mlp(1, 3, 1);

   mlp.setNameOfSingleInputVariable(0, "x");
   mlp.setNameOfSingleOutputVariable(0, "y");

   // Geodesic problem object

   GeodesicProblem gp(&mlp);

   // Save initial guess

   gp.saveResults("../Data/GeodesicProblem/InitialGuessGeodesicProblem.dat");

   double evaluation = gp.calculateEvaluation();

   std::cout << "Evaluation:" << std::endl
             << evaluation << std::endl;

   Vector<double> gradient = gp.calculateGradient();

   std::cout << "Gradient:" << std::endl
             << gradient << std::endl;

   // Quasi-Newton method object

   QuasiNewtonMethod qnm(&gp);	
    
   qnm.setReserveEvaluationHistory(true);
   qnm.setReserveGradientNormHistory(true);

   qnm.setEvaluationGoal(0.0); 

   qnm.setMinimumImprovement(0.0);
   qnm.setMaximumNumberOfEpochs(100); 
   qnm.setDisplayPeriod(1); 

   // Train neural network

   gp.setNumberOfEvaluations(0);

   qnm.train(); 

   int numberOfEvaluations = gp.getNumberOfEvaluations();

   std::cout << std::endl
             << "Number of evaluations: " 
             << numberOfEvaluations <<  std::endl;	

   // Save neural network to file

   mlp.save("../Data/GeodesicProblem/MultilayerPerceptronGeodesicProblem.dat");

   mlp.saveExpression("../Data/GeodesicProblem/ExpressionGeodesicProblem.dat");

   // Save results to file

   gp.saveResults("../Data/GeodesicProblem/ResultsGeodesicProblem.dat");

   // Save training history to file

   qnm.saveTrainingHistory("../Data/GeodesicProblem/TrainingHistoryGeodesicProblem.dat");

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
