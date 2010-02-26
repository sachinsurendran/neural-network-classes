/****************************************************************************************************************/
/*                                                                                                              */ 
/*   Flood: An Open Source Neural Networks C++ Library                                                          */
/*   www.cimne.com/flood                                                                                        */
/*                                                                                                              */
/*   M I N I M U M   D R A G   P R O B L E M   A P P L I C A T I O N                                            */
/*                                                                                                              */
/*   Roberto Lopez                                                                                              */ 
/*   International Center for Numerical Methods in Engineering (CIMNE)                                          */
/*   Technical University of Catalonia (UPC)                                                                    */
/*   Barcelona, Spain                                                                                           */
/*   E-mail: rlopez@cimne.upc.edu                                                                               */ 
/*                                                                                                              */  
/****************************************************************************************************************/

/// This main function can be used as a template for solving the minimum drag problem by means of a multilayer 
/// perceptron. 
/// It uses the quasi-Newton training algorithm. 

// System includes

#include <iostream>
#include <time.h>

// Utilities includes

#include "../Flood/Utilities/Vector.h"
#include "../Flood/Utilities/Matrix.h"

// Network architecture includes

#include "../Flood/MultilayerPerceptron/MultilayerPerceptron.h"

// Objective functional includes

#include "../Flood/ObjectiveFunctional/MinimumDragProblem.h"

// Learning algorithm includes

#include "../Flood/TrainingAlgorithm/RandomSearch.h"
#include "../Flood/TrainingAlgorithm/QuasiNewtonMethod.h"

using namespace Flood;

int main(void)
{
   std::cout << std::endl
		     << "Flood Neural Network. Minimum Drag Problem." << std::endl;	

   srand((unsigned)time(NULL));

   // Multilayer perceptron object

   MultilayerPerceptron mlp(1,3,1);   

   mlp.setNameOfSingleInputVariable(0, "x");
   mlp.setNameOfSingleOutputVariable(0, "y");


   // Minimum drag problem object

   MinimumDragProblem mdp(&mlp);  

   mdp.saveResults("../Data/MinimumDragProblem/InitialGuessMinimumDragProblem.dat");


   // Random search object

   RandomSearch rs(&mdp);
   rs.setMaximumNumberOfEpochs(100);
   rs.setDisplayPeriod(100);

   rs.train();


   // Quasi-Newton method object

   QuasiNewtonMethod qnm(&mdp);

   qnm.setReserveEvaluationHistory(true);
   qnm.setReserveGradientNormHistory(true);

   qnm.setMaximumNumberOfEpochs(1000);
   qnm.setDisplayPeriod(100);

   qnm.train();
   
   int numberOfEvaluations = mdp.getNumberOfEvaluations();

   std::cout << std::endl
             << "Number of evaluations: " << numberOfEvaluations << std::endl;
   
   // Save 
   
   mlp.save("../Data/MinimumDragProblem/MultilayerPerceptronMinimumDragProblem.dat");
   mlp.saveExpression("../Data/MinimumDragProblem/ExpressionMinimumDragProblem.dat");   
   mdp.saveResults("../Data/MinimumDragProblem/Results.dat");

   qnm.saveTrainingHistory("../Data/MinimumDragProblem/TrainingHistoryMinimumDragProblem.dat");

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
