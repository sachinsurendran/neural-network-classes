/******************************************************************************/
/*                                                                            */ 
/*   Flood: An Open Source Neural Networks C++ Library                        */
/*   www.cimne.com/flood                                                      */
/*                                                                            */
/*   B R A C H I S T O C H R O N E   P R O B L E M   A P P L I C A T I O N    */
/*                                                                            */
/*   Roberto Lopez                                                            */
/*   International Center for Numerical Methods in Engineering (CIMNE)        */
/*   Technical University of Catalonia (UPC)                                  */
/*   Barcelona, Spain                                                         */
/*   E-mail: rlopez@cimne.upc.edu                                             */ 
/*                                                                            */  
/******************************************************************************/

/// This main function can be used as a template for solving the 
/// brachistochrone problem by means of a multilayer perceptron. 
/// It uses the quasi-Newton method for training.

// System includes

#include <iostream>
#include <time.h>

// Network architecture includes

#include "../Flood/MultilayerPerceptron/MultilayerPerceptron.h"

// Objective functional includes

#include "../Flood/ObjectiveFunctional/BrachistochroneProblem.h"

// Training algorithm includes

#include "../Flood/TrainingAlgorithm/RandomSearch.h"
#include "../Flood/TrainingAlgorithm/QuasiNewtonMethod.h"

using namespace Flood;

int main(void)
{
   std::cout << std::endl
             << "Flood Neural Network. Brachistochrone Problem Application." 
             << std::endl;	

   double xA = 0.0;
   double yA = 1.0;
   double xB = 1.0;
   double yB = 0.0;

   srand((unsigned)time(NULL));

   // Multilayer perceptron object

   MultilayerPerceptron mlp(1,3,1);

   mlp.setNameOfSingleInputVariable(0, "x");
   mlp.setNameOfSingleOutputVariable(0, "y");

   mlp.setMinimumOfSingleInputVariable(0, xA);
   mlp.setMaximumOfSingleInputVariable(0, xB);

   // Brachistochrone problem object

   BrachistochroneProblem bp(&mlp);
   bp.setProblem(xA, yA, xB, yB);

   // Save initial guess to file

   bp.saveResults("../Data/BrachistochroneProblem/InitialGuessBrachistochroneProblem.dat");

   // Random search object

   RandomSearch rs(&bp);
   rs.setMaximumNumberOfEpochs(100);

   rs.train();

   // Quasi-Newton method object

   QuasiNewtonMethod qnm(&bp);	

   qnm.setReserveEvaluationHistory(true);
   qnm.setReserveGradientNormHistory(true);

   qnm.setMaximumNumberOfEpochs(1000); 
   qnm.setDisplayPeriod(10); 

   // Train neural network

   int numberOfEvaluations = 0;

   bp.setNumberOfEvaluations(numberOfEvaluations);

   qnm.train();

   numberOfEvaluations = bp.getNumberOfEvaluations();

   std::cout << std::endl
             << "Number of evaluations: " 
             << numberOfEvaluations <<  std::endl;	

   // Save trained neural network to file

   mlp.save("../Data/BrachistochroneProblem/MultilayerPerceptron.dat");

   mlp.saveExpression("../Data/BrachistochroneProblem/ExpressionBrachistochroneProblem.dat");

   // Save results to file

   bp.saveResults("../Data/BrachistochroneProblem/ResultsBrachistochroneProblem.dat");

   // Save training history to file

   qnm.saveTrainingHistory("../Data/BrachistochroneProblem/TrainingHistoryBrachistochroneProblem.dat");
   
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
