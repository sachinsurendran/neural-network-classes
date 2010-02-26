/****************************************************************************************************************/
/*                                                                                                              */
/*   Flood: An Open Source Neural Networks C++ Library                                                          */
/*   www.cimne.com/flood                                                                                        */
/*                                                                                                              */
/*   P L A N E - C Y L I N D E R   A P P L I C A T I O N                                                        */
/*                                                                                                              */
/*   Roberto Lopez                                                                                              */
/*   International Center for Numerical Methods in Engineering (CIMNE)                                          */
/*   Technical University of Catalonia (UPC)                                                                    */
/*   Barcelona, Spain                                                                                           */
/*   E-mail: rlopez@cimne.upc.edu                                                                               */
/*                                                                                                              */
/****************************************************************************************************************/

/// This main function is used to solve the plane-cylinder problem, which is a constrained function optimization 
/// problem. 
/// The quasi-Newton method is used here. 

// System includes

#include <iostream>
#include <time.h>
#include <stdexcept>

// Utilities includes

#include "../Flood/Utilities/Vector.h"
#include "../Flood/Utilities/Matrix.h"

// Multilayer perceptron includes

#include "../Flood/MultilayerPerceptron/MultilayerPerceptron.h"

// Objective functional includes

#include "../Flood/ObjectiveFunctional/PlaneCylinder.h"

// Training algorithm includes

#include "../Flood/TrainingAlgorithm/QuasiNewtonMethod.h"
#include "../Flood/TrainingAlgorithm/RandomSearch.h"


using namespace Flood;

int main(void)
{
   std::cout << std::endl
             << "Flood Neural Network: Plane Cylinder Application." << std::endl;

   srand( (unsigned)time( NULL ) );


   // Multilayer perceptron object

   MultilayerPerceptron multilayerPerceptron;

   multilayerPerceptron.setNumberOfIndependentParameters(2);

   multilayerPerceptron.setNameOfSingleIndependentParameter(0, "x");
   multilayerPerceptron.setNameOfSingleIndependentParameter(1, "y");
   
   multilayerPerceptron.setMinimumOfSingleIndependentParameter(0, -5.12);
   multilayerPerceptron.setMinimumOfSingleIndependentParameter(1, -5.12);
   
   multilayerPerceptron.setMaximumOfSingleIndependentParameter(0, 5.12);
   multilayerPerceptron.setMaximumOfSingleIndependentParameter(1, 5.12);

   multilayerPerceptron.initIndependentParametersUniform(-5.12, 5.12);

   multilayerPerceptron.setPreAndPostProcessingMethod(MultilayerPerceptron::MinimumAndMaximum);

   // Plane cylinder object

   PlaneCylinder planeCylinder(&multilayerPerceptron);
   planeCylinder.setPenalty(100.0);

   int numberOfVariables = multilayerPerceptron.getNumberOfIndependentParameters();
   
   Vector<double> initialGuess(numberOfVariables, 1.2);

   multilayerPerceptron.setIndependentParameters(initialGuess);

   // Print initial guess

   std::cout << std::endl
	         << "Initial guess: " << std::endl
			 << initialGuess << std::endl;
	          
   // Evaluation
  
   double evaluation = planeCylinder.calculateEvaluation();
   
   std::cout << std::endl
             << "Evaluation:" << std::endl
             << evaluation << std::endl;
             
   // Gradient vector

   Vector<double> gradient = planeCylinder.calculateGradient();
   
   std::cout << std::endl
             << "Gradient:" << std::endl
			 << gradient << std::endl;            

   // Hessian matrix

   Matrix<double> hessian = planeCylinder.calculateHessian();
   
   std::cout << std::endl
             << "Hessian:" << std::endl
			 << hessian;

   // Inverse Hessian matrix

   Matrix<double> inverseHessian = planeCylinder.calculateInverseHessian();
   
   std::cout << std::endl
             << "Inverse Hessian:" << std::endl
			 << inverseHessian;
             
   // Random search object

   RandomSearch randomSearch(&planeCylinder);
   randomSearch.train();

   // Quasi-Newton method object

   QuasiNewtonMethod quasiNewtonMethod(&planeCylinder);
   quasiNewtonMethod.train();

   // Print minimal argument

   Vector<double> minimalArgument = multilayerPerceptron.getIndependentParameters();

   std::cout << std::endl
	         << "Minimal argument:" << std::endl 
			 << minimalArgument << std::endl;

   std::cout << std::endl;         

   return(0);
}  
