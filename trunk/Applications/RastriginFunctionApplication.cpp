/****************************************************************************************************************/
/*                                                                                                              */
/*   Flood: An Open Source Neural Networks C++ Library                                                          */
/*   www.cimne.com/flood                                                                                        */
/*                                                                                                              */
/*   R A S T R I G I N   F U N C T I O N   A P P L I C A T I O N                                                */
/*                                                                                                              */
/*   Roberto Lopez                                                                                              */
/*   International Center for Numerical Methods in Engineering (CIMNE)                                          */
/*   Technical University of Catalonia (UPC)                                                                    */
/*   Barcelona, Spain                                                                                           */
/*   E-mail: rlopez@cimne.upc.edu                                                                               */
/*                                                                                                              */
/****************************************************************************************************************/

/// This main function can be used as a template for optimizing the Rastrigin's function. 
/// It uses the evolutionary algorithm to locate the optimal argument and a quasi-Newton method to refine it. 

// System includes

#include <iostream>
#include <sstream>
#include <time.h>
#include <stdexcept>

// Utilities includes

#include "../Flood/Utilities/Vector.h"
#include "../Flood/Utilities/Matrix.h"

// Multilayer perceptron includes

#include "../Flood/MultilayerPerceptron/MultilayerPerceptron.h"

// Objective functional includes

#include "../Flood/ObjectiveFunctional/RastriginFunction.h"

// Training algorithm includes

#include "../Flood/TrainingAlgorithm/QuasiNewtonMethod.h"
#include "../Flood/TrainingAlgorithm/EvolutionaryAlgorithm.h"

using namespace Flood;

int main(void)
{
   std::cout << std::endl
             << "Flood Neural Network. Rastrigin Function Application." << std::endl;

   srand((unsigned)time(NULL));

   // Problem parameters 

   unsigned int numberOfVariables = 3;

   // Multilayer perceptron object

   MultilayerPerceptron mlp;

   // Set multilayer perceptron stuff

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

   // Rastrigin function object

   RastriginFunction rf(&mlp);

   rf.setNumberOfVariables(numberOfVariables);
   
   // Evolutionary algorithm object

   EvolutionaryAlgorithm ea(&rf);

   ea.setPopulationSize(10*numberOfVariables);
   ea.setMaximumNumberOfGenerations(1000);

   ea.setDisplayPeriod(100);

   ea.initPopulationUniform(-5.12, 5.12);

   ea.train();

   // Quasi-Newton method object

   QuasiNewtonMethod qnm(&rf);

   qnm.train();

   // Print minimal argument 

   Vector<double> minimalArgument = mlp.getIndependentParameters();

   std::cout << std::endl
	         << "Minimal argument:" << std::endl
	         << minimalArgument << std::endl;

   // Save solution

   mlp.save("..//Data//RastriginFunction//MultilayerPerceptron.dat");

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
