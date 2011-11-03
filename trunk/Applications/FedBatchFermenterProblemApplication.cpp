/****************************************************************************************************************/
/*                                                                                                              */
/*   Flood: An Open Source Neural Networks C++ Library                                                          */
/*   www.cimne.com/flood                                                                                        */
/*                                                                                                              */
/*   F E D   B A T C H   F E R M E N T E R   P R O B L E M   A P P L I C A T I O N                              */
/*                                                                                                              */
/*   Roberto Lopez                                                                                              */
/*   International Center for Numerical Methods in Engineering (CIMNE)                                          */
/*   Technical University of Catalonia (UPC)                                                                    */
/*   Barcelona, Spain                                                                                           */
/*   E-mail: rlopez@cimne.upc.edu                                                                               */
/*                                                                                                              */
/****************************************************************************************************************/

/// This main function can be used as a template for solving the 
/// fed batch fermenter problem by means of a multilayer perceptron. It uses
/// the quasi-Newton method for training.

// System includes

#include <iostream>
#include <time.h>

#include "../Flood/Utilities/InputTargetDataSet.h"


// Multilayer perceptron includes

#include "../Flood/MultilayerPerceptron/MultilayerPerceptron.h"

// Objective functional includes

#include "../Flood/ObjectiveFunctional/FedBatchFermenterProblem.h"

// Training algorithm includes

#include "../Flood/TrainingAlgorithm/RandomSearch.h"
#include "../Flood/TrainingAlgorithm/EvolutionaryAlgorithm.h"
#include "../Flood/TrainingAlgorithm/QuasiNewtonMethod.h"

using namespace Flood;

int main(void)
{
   std::cout << std::endl
             << "Flood Neural Network. Fed Batch Fermenter Problem."
             << std::endl;

//   srand((unsigned)time(NULL));

   // Problem parameters

   double maximumFeedRate = 12.0;
   double finalTime = 54.0;

   // Multilayer perceptron object

   MultilayerPerceptron mlp(1,3,1);

   mlp.setNameOfSingleInputVariable(0, "Time");
   mlp.setNameOfSingleOutputVariable(0, "FeedRate");

   mlp.setMinimumOfSingleInputVariable(0, 0.0);
   mlp.setMaximumOfSingleInputVariable(0, finalTime);

   mlp.setMinimumOfSingleOutputVariable(0, 0.0);
   mlp.setMaximumOfSingleOutputVariable(0, maximumFeedRate);

   mlp.setLowerBoundOfSingleOutputVariable(0, 0.0);
   mlp.setUpperBoundOfSingleOutputVariable(0, maximumFeedRate);

   mlp.setPreAndPostProcessingMethod(MultilayerPerceptron::MinimumAndMaximum);


   // Fed batch fermenter problem object

   FedBatchFermenterProblem fbfp(&mlp);

   fbfp.setVolumeErrorWeight(1.0e-3);
   fbfp.setYieldWeight(1.0e-10);

   fbfp.setTolerance(1.0e-9);
   fbfp.setInitialSize(5000);


   // Random search object

   RandomSearch rs(&fbfp);
   rs.setMaximumNumberOfEpochs(100);
   rs.train();

   fbfp.saveResults("../Data/FedBatchFermenterProblem/InitialGuessFedBatchFermenter.dat");


   // Quasi-Newton method object

   QuasiNewtonMethod qnm(&fbfp);

   qnm.setReserveEvaluationHistory(true);
   qnm.setReserveGradientNormHistory(true);

   qnm.setMaximumNumberOfEpochs(1000);
   qnm.setDisplayPeriod(10);
  
   qnm.train();

   // Number of evaluations

   int numberOfEvaluations = fbfp.getNumberOfEvaluations();

   std::cout << std::endl
       	     << "Number of evaluations = " << numberOfEvaluations << std::endl;

   // Save all 

   mlp.save("../Data/FedBatchFermenterProblem/MultilayerPerceptronFedBatchFermenter.dat");
   mlp.saveExpression("../Data/FedBatchFermenterProblem/ExpressionFedBatchFermenter.dat");

   fbfp.saveResults("../Data/FedBatchFermenterProblem/ResultsFedBatchFermenter.dat");

   qnm.saveTrainingHistory("../Data/FedBatchFermenterProblem/TrainingHistoryFedBatchFermenter.dat");

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
