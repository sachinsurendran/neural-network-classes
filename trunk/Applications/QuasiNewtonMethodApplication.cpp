/****************************************************************************************************************/
/*                                                                                                              */ 
/*   Flood: An Open Source Neural Networks C++ Library                                                          */
/*   www.cimne.com/flood                                                                                        */
/*                                                                                                              */
/*   Q U A S I - N E W T O N   M E T H O D   A P P L I C A T I O N                                              */
/*                                                                                                              */
/*   Roberto Lopez                                                                                              */ 
/*   International Center for Numerical Methods in Engineering (CIMNE)                                          */
/*   Technical University of Catalonia (UPC)                                                                    */
/*   Barcelona, Spain                                                                                           */
/*   E-mail: rlopez@cimne.upc.edu                                                                               */ 
/*                                                                                                              */  
/****************************************************************************************************************/

/// This application is an usage example of the QuasiNewtonMethod class in Flood.

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

#include "../Flood/ObjectiveFunctional/MeanSquaredError.h"

// Training algorithm includes

#include "../Flood/TrainingAlgorithm/QuasiNewtonMethod.h"

using namespace Flood;

int main(void)
{
   std::cout << std::endl
             << "Flood Neural Network. QuasiNewtonMethod Application." 
             << std::endl;	


   srand((unsigned)time(NULL));

   // Input-target data set object

   InputTargetDataSet inputTargetDataSet;

   inputTargetDataSet.load("../Data/QuasiNewtonMethod/InputTargetDataSet.dat");

   // Multilayer perceptron object 

   MultilayerPerceptron multilayerPerceptron(1,6,1);

   // Mean squared error object

   MeanSquaredError 
   meanSquaredError(&multilayerPerceptron, &inputTargetDataSet);

   // Quasi-Newton method object

   QuasiNewtonMethod quasiNewtonMethod(&meanSquaredError);
   quasiNewtonMethod.setReserveFreeParametersHistory(true);
   quasiNewtonMethod.setReserveFreeParametersNormHistory(true);
   quasiNewtonMethod.setReserveEvaluationHistory(true);
   quasiNewtonMethod.setReserveGradientHistory(true);
   quasiNewtonMethod.setReserveGradientNormHistory(true);
   quasiNewtonMethod.setReserveElapsedTimeHistory(true);
   quasiNewtonMethod.setReserveTrainingDirectionHistory(true);
   quasiNewtonMethod.setReserveTrainingRateHistory(true);
   
   // Set stopping criteria

   quasiNewtonMethod.setEvaluationGoal(0.001);
   quasiNewtonMethod.setGradientNormGoal(0.001);
   quasiNewtonMethod.setMaximumTime(1000.0);
   quasiNewtonMethod.setMaximumNumberOfEpochs(100);

   // Set user stuff

   quasiNewtonMethod.setDisplayPeriod(10);

   // Train neural network

   quasiNewtonMethod.train();

   // Save quasi-Newton method object

   quasiNewtonMethod.save("../Data/QuasiNewtonMethod/QuasiNewtonMethod.dat");

   // Save training history

   quasiNewtonMethod.saveTrainingHistory("../Data/QuasiNewtonMethod/TrainingHistory.dat");

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
