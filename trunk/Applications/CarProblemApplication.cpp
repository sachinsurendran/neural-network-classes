/****************************************************************************************************************/
/*                                                                                                              */ 
/*   Flood: An Open Source Neural Networks C++ Library                                                          */
/*   www.cimne.com/flood                                                                                        */
/*                                                                                                              */
/*   C A R   P R O B L E M    A P P L I C A T I O N                                                             */
/*                                                                                                              */
/*   Roberto Lopez                                                                                              */ 
/*   International Center for Numerical Methods in Engineering (CIMNE)                                          */
/*   Technical University of Catalonia (UPC)                                                                    */
/*   Barcelona, Spain                                                                                           */
/*   E-mail: rlopez@cimne.upc.edu                                                                               */ 
/*                                                                                                              */  
/****************************************************************************************************************/

/// This main function can be used as a template for solving the car problem by means of a multilayer perceptron. 
/// It uses the quasi-Newton method for training. 

// System includes

#include <iostream>
#include <time.h>

// Utilities includes

#include "../Flood/Utilities/Vector.h"
#include "../Flood/Utilities/Matrix.h"

// Multilayer perceptron includes

#include "../Flood/MultilayerPerceptron/MultilayerPerceptron.h"

// Objective functional includes

#include "../Flood/ObjectiveFunctional/CarProblem.h"

// Training algorithm includes

#include "../Flood/TrainingAlgorithm/RandomSearch.h"
#include "../Flood/TrainingAlgorithm/QuasiNewtonMethod.h"

using namespace Flood;

int main(void)
{
   std::cout << "Flood Neural Network. Car Problem Application." << std::endl;

   srand((unsigned)time(NULL));

   // Multilayer perceptron object

   MultilayerPerceptron mlp(1,3,2);

   mlp.setNameOfSingleInputVariable(0, "Time");

   mlp.setNameOfSingleOutputVariable(0, "ThrottleAcceleration");
   mlp.setNameOfSingleOutputVariable(1, "BrackingDeceleration");
    
   mlp.setNumberOfIndependentParameters(1);
   
   mlp.setNameOfSingleIndependentParameter(0, "FinalTime");
   
   mlp.setLowerBoundOfSingleIndependentParameter(0, 0.0);

   mlp.setPreAndPostProcessingMethod(MultilayerPerceptron::None);

   mlp.setDisplay(false);

  
   // Car problem object

   CarProblem carProblem(&mlp);
   
   carProblem.setTolerance(1.0e-12);
   carProblem.setInitialSize(1000);

   carProblem.setInitialPosition(0.0);
   carProblem.setInitialVelocity(0.0);

   carProblem.setFinalPositionGoal(1.0);
   carProblem.setFinalVelocityGoal(0.0);

   carProblem.setFinalPositionErrorWeight(1.0);
   carProblem.setFinalVelocityErrorWeight(1.0);
   carProblem.setFinalTimeWeight(1.0e-3);


   // Random search object

   RandomSearch rs(&carProblem);
   rs.train();

   carProblem.saveResults("../Data/CarProblem/InitialGuessCarProblem.dat");

   // Quasi-Newton method object

   QuasiNewtonMethod quasiNewtonMethod(&carProblem);

   quasiNewtonMethod.setReserveEvaluationHistory(true);
   quasiNewtonMethod.setReserveGradientNormHistory(true);

   quasiNewtonMethod.setMaximumNumberOfEpochs(1000);
   quasiNewtonMethod.setDisplayPeriod(10);

   int numberOfEvaluations = 0;

   carProblem.setNumberOfEvaluations(numberOfEvaluations);

   quasiNewtonMethod.train();

   numberOfEvaluations = carProblem.getNumberOfEvaluations();

   std::cout << std::endl
             << "Number of evaluations: " 
             << numberOfEvaluations <<  std::endl;	

   // Save neural network to file

   mlp.save("../Data/CarProblem/MultilayerPerceptronCarProblem.dat");
   mlp.saveExpression("../Data/CarProblem/ExpressionCarProblem.dat");

   // Save results to file

   carProblem.saveResults("../Data/CarProblem/ResultsCarProblem.dat");

   // Save training history to file

   quasiNewtonMethod.saveTrainingHistory("../Data/CarProblem/TrainingHistoryCarProblem.dat");

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
