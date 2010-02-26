/****************************************************************************************************************/
/*                                                                                                              */
/*   Flood: An Open Source Neural Networks C++ Library                                                          */
/*   www.cimne.com/flood                                                                                        */
/*                                                                                                              */
/*   D E   J O N G   F U N C T I O N   A P P L I C A T I O N                                                    */
/*                                                                                                              */
/*   Roberto Lopez                                                                                              */
/*   International Center for Numerical Methods in Engineering (CIMNE)                                          */
/*   Technical University of Catalonia (UPC)                                                                    */
/*   Barcelona, Spain                                                                                           */
/*   E-mail: rlopez@cimne.upc.edu                                                                               */
/*                                                                                                              */
/****************************************************************************************************************/

/// This main function can be used as a template for optimizing the De Jong function. 
/// It uses the gradient descent method.

// System includes

#include <iostream>
#include <sstream>
#include <time.h>

// Utilities includes

#include "../Flood/Utilities/Vector.h"
#include "../Flood/Utilities/Matrix.h"

// Multilayer perceptron includes

#include "../Flood/MultilayerPerceptron/MultilayerPerceptron.h"

// Objective functional includes

#include "../Flood/ObjectiveFunctional/DeJongFunction.h"

// Training algorithm includes

#include "../Flood/TrainingAlgorithm/GradientDescent.h"

using namespace Flood;

int main(void)
{
   std::cout << std::endl
             << "Flood Neural Network. De Jong Function Application." << std::endl;	

   srand((unsigned)time(NULL));

   // Problem parameters 

   unsigned int numberOfVariables = 2;

   // Multilayer perceptron object (independent parameters)

   MultilayerPerceptron mlp;

   mlp.setNumberOfIndependentParameters(numberOfVariables);

   for(unsigned int i = 0; i < numberOfVariables; i++)
   {
      std::stringstream buffer;
      buffer << "x" << i+1;
      mlp.setNameOfSingleIndependentParameter(i, buffer.str());

      mlp.setMinimumOfSingleIndependentParameter(i, -5.12);
      mlp.setMaximumOfSingleIndependentParameter(i, 5.12); 
   }

   mlp.initIndependentParametersUniform(-5.12, 5.12);

   mlp.setPreAndPostProcessingMethod(MultilayerPerceptron::MinimumAndMaximum);

   // De Jong function object

   DeJongFunction djf(&mlp);

   djf.setNumberOfVariables(numberOfVariables);
   
   Vector<double> initialGuess(numberOfVariables, 1.0);

   mlp.setIndependentParameters(initialGuess);

   // Evaluation
  
   double evaluation = djf.calculateEvaluation();
   
   std::cout << std::endl
             << "Evaluation:" << std::endl
             << evaluation << std::endl;
             
   // Gradient vector

   Vector<double> gradient = djf.calculateGradient();
   
   std::cout << std::endl
             << "Gradient:" << std::endl
             << gradient << std::endl;
             
   // Hessian matrix

   Matrix<double> hessian = djf.calculateHessian();
   
   std::cout << std::endl
             << "Hessian:" << std::endl
             << hessian << std::endl;
             
   // Inverse Hessian matrix

   Matrix<double> inverseHessian = djf.calculateInverseHessian();
   
   std::cout << "Inverse Hessian:" << std::endl
             << inverseHessian << std::endl;

   // Gradient descent object

   GradientDescent gd(&djf);

   gd.setReserveEvaluationHistory(true);
   gd.setReserveGradientNormHistory(true);

   gd.setMaximumNumberOfEpochs(1);
   gd.setTrainRateTolerance(1.0e-12);

   gd.train();

   // Print minimal argument 

   Vector<double> minimalArgument = mlp.getIndependentParameters();

   std::cout << std::endl
	         << "Minimal argument:" << std::endl
	         << minimalArgument << std::endl;

   // Save all stuff 

   mlp.save("../Data/DeJongFunction/MultilayerPerceptron.dat");

   gd.saveTrainingHistory("../Data/DeJongFunction/TrainingHistory.dat");

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
