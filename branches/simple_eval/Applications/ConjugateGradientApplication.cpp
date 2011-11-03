/****************************************************************************************************************/
/*                                                                                                              */ 
/*   Flood: An Open Source Neural Networks C++ Library                                                          */
/*   www.cimne.com/flood                                                                                        */
/*                                                                                                              */
/*   C O N J U G A T E   G R A D I E N T   A P P L I C A T I O N                                                */
/*                                                                                                              */
/*   Roberto Lopez                                                                                              */ 
/*   International Center for Numerical Methods in Engineering (CIMNE)                                          */
/*   Technical University of Catalonia (UPC)                                                                    */
/*   Barcelona, Spain                                                                                           */
/*   E-mail: rlopez@cimne.upc.edu                                                                               */ 
/*                                                                                                              */  
/****************************************************************************************************************/

/// This application is an usage example of the ConjugateGradient class in Flood.

// System includes

#include <iostream>
#include <time.h>

// Utilities includes

#include "../Flood/Utilities/Vector.h"
#include "../Flood/Utilities/Matrix.h"

#include "../Flood/Utilities/InputTargetDataSet.h"

// Multilayer perceptron includes

#include "../Flood/MultilayerPerceptron/MultilayerPerceptron.h"

// Performance functional includes

#include "../Flood/ObjectiveFunctional/SumSquaredError.h"

// Learning algorithm includes

#include "../Flood/TrainingAlgorithm/ConjugateGradient.h"

using namespace Flood;

int main(void)
{
   std::cout << std::endl
             << "Flood Neural Network. Conjugate Gradient Application." << std::endl;	

   // Input-target data set object

   InputTargetDataSet inputTargetDataSet;

   inputTargetDataSet.load("../Data/ConjugateGradient/InputTargetDataSet.dat");

   // Multilayer perceptron object 

   MultilayerPerceptron multilayerPerceptron(1,3,1);

   // Sum of squares error object

   SumSquaredError sumSquaredError(&multilayerPerceptron, &inputTargetDataSet);

   // Conjugate gradient object

   ConjugateGradient conjugateGradient(&sumSquaredError);

   conjugateGradient.setReserveElapsedTimeHistory(true);
   conjugateGradient.setReserveFreeParametersHistory(true);
   conjugateGradient.setReserveFreeParametersNormHistory(true);
   conjugateGradient.setReserveEvaluationHistory(true);
   conjugateGradient.setReserveGradientHistory(true);
   conjugateGradient.setReserveGradientNormHistory(true);
   conjugateGradient.setReserveTrainingDirectionHistory(true);
   conjugateGradient.setReserveTrainingDirectionNormHistory(true);
   conjugateGradient.setReserveTrainingRateHistory(true);

   // Set stopping criteria

   conjugateGradient.setEvaluationGoal(0.01);
   conjugateGradient.setGradientNormGoal(0.0);
   conjugateGradient.setMaximumTime(100.0);
   conjugateGradient.setMaximumNumberOfEpochs(100);

   // Set user stuff

   conjugateGradient.setDisplayPeriod(10);

   // Train neural network

   conjugateGradient.train();

   conjugateGradient.saveTrainingHistory("../Data/ConjugateGradient/TrainingHistory.dat");

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
