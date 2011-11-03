/****************************************************************************************************************/
/*                                                                                                              */ 
/*   Flood: An Open Source Neural Networks C++ Library                                                          */
/*   www.cimne.com/flood                                                                                        */
/*                                                                                                              */
/*   R A N D O M   S E A R C H   A P P L I C A T I O N                                                          */
/*                                                                                                              */
/*   Roberto Lopez                                                                                              */ 
/*   International Center for Numerical Methods in Engineering (CIMNE)                                          */
/*   Technical University of Catalonia (UPC)                                                                    */
/*   Barcelona, Spain                                                                                           */
/*   E-mail: rlopez@cimne.upc.edu                                                                               */ 
/*                                                                                                              */  
/****************************************************************************************************************/

/// This application is an usage example of the RandomSearch class in Flood.

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

#include "../Flood/ObjectiveFunctional/SumSquaredError.h"

// Learning algorithm includes

#include "../Flood/TrainingAlgorithm/RandomSearch.h"

using namespace Flood;

int main(void)
{
   std::cout << std::endl
             << "Flood Neural Network. Random Search Application." 
             << std::endl;	

   // Input-target data set object

   InputTargetDataSet inputTargetDataSet;

   inputTargetDataSet.load("../Data/RandomSearch/InputTargetDataSet.dat");

   // Multilayer perceptron object 

   MultilayerPerceptron multilayerPerceptron(1,3,1);

   // Sum squared error object

   SumSquaredError 
   sumSquaredError(&multilayerPerceptron, &inputTargetDataSet);

   // Random search object

   RandomSearch randomSearch(&sumSquaredError);

   // Set reserve training history
   
   randomSearch.setReserveFreeParametersHistory(true);
   randomSearch.setReserveFreeParametersNormHistory(true);
   randomSearch.setReserveEvaluationHistory(true);
   randomSearch.setReserveElapsedTimeHistory(true);
   
   // Set stopping criteria

   randomSearch.setEvaluationGoal(0.01);
   randomSearch.setMaximumTime(100.0);
   randomSearch.setMaximumNumberOfEpochs(100);

   // Set user stuff

   randomSearch.setDisplayPeriod(10);

   // Train neural network

   randomSearch.train();

   std::cout << std::endl
             << "Number of evaluations: " << sumSquaredError.getNumberOfEvaluations() << std::endl;

   // Save training history

   randomSearch.saveTrainingHistory("../Data/RandomSearch/RandomSearchTrainingHistory.dat");

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
