/****************************************************************************************************************/
/*                                                                                                              */
/*   Flood: An Open Source Neural Networks C++ Library                                                          */
/*   www.cimne.com/flood                                                                                        */
/*                                                                                                              */
/*   C A T E N A R Y   P R O B L E M   A P P L I C A T I O N                                                    */
/*                                                                                                              */
/*   Roberto Lopez                                                                                              */ 
/*   International Center for Numerical Methods in Engineering (CIMNE)                                          */
/*   Technical University of Catalonia (UPC)                                                                    */
/*   Barcelona, Spain                                                                                           */
/*   E-mail: rlopez@cimne.upc.edu                                                                               */ 
/*                                                                                                              */
/****************************************************************************************************************/

/// This main function can be used as a template for solving the catenary problem by means of a multilayer 
/// perceptron. 
/// It uses the quasi-Newton method for training. 

// System includes

#include <iostream>
#include <time.h>

// Multilayer perceptron includes

#include "../Flood/MultilayerPerceptron/MultilayerPerceptron.h"

// Objective functional includes

#include "../Flood/ObjectiveFunctional/CatenaryProblem.h"

// Training algorithm includes

#include "../Flood/TrainingAlgorithm/QuasiNewtonMethod.h"


using namespace Flood;

int main(void)
{
   std::cout << std::endl
             << "Flood Neural Network. Catenary Problem Application." << std::endl;	

   srand( (unsigned)time( NULL ) );

   // Multilayer perceptron object

   MultilayerPerceptron mlp(1,3,1);

   mlp.setNameOfSingleInputVariable(0, "x");
   mlp.setNameOfSingleOutputVariable(0, "y");

   // Catenary problem object

   CatenaryProblem cp(&mlp);
   cp.setPotentialEnergyWeight(1.0);
   cp.setLengthErrorWeight(100.0);

   cp.saveResults("../Data/CatenaryProblem/InitialGuessCatenaryProblem.dat");

   // Quasi-Newton method object

   QuasiNewtonMethod qnm(&cp);	

   qnm.setMaximumNumberOfEpochs(1000); 
   qnm.setDisplayPeriod(10); 

   qnm.setReserveEvaluationHistory(true);
   qnm.setReserveGradientNormHistory(true);

   // Train neural network

   cp.setNumberOfEvaluations(0);

   qnm.train(); 

   int numberOfEvaluations = cp.getNumberOfEvaluations();

   std::cout << std::endl
             << "Number of evaluations: " 
             << numberOfEvaluations <<  std::endl;	


   // Save neural network to file

   mlp.save("../Data/CatenaryProblem/MultilayerPerceptron.dat");
   mlp.saveExpression("../Data/CatenaryProblem/ExpressionCatenaryProblem.dat");

   // Save results to file

   cp.saveResults("../Data/CatenaryProblem/ResultsCatenaryProblem.dat");

   // Save training history to file

   qnm.saveTrainingHistory("../Data/CatenaryProblem/TrainingHistoryCatenaryProblem.dat");

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
