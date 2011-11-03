/****************************************************************************************************************/
/*                                                                                                              */
/*   Flood: An Open Source Neural Networks C++ Library                                                          */
/*   www.cimne.com/flood                                                                                        */
/*                                                                                                              */
/*   T E M P L A T E   A P P L I C A T I O N                                                                    */
/*                                                                                                              */
/*   Roberto Lopez                                                                                              */
/*   International Center for Numerical Methods in Engineering (CIMNE)                                          */
/*   Technical University of Catalonia (UPC)                                                                    */
/*   Barcelona, Spain                                                                                           */
/*   E-mail: rlopez@cimne.upc.edu                                                                               */
/*                                                                                                              */
/****************************************************************************************************************/

/// This main function can be used as a template for constructing any application with Flood.  

// System includes

#include <iostream>
#include <time.h>
#include <vector>

// Utilities includes

#include "../Flood/Utilities/Vector.h"
#include "../Flood/Utilities/Matrix.h"

#include "../Flood/Utilities/InputTargetDataSet.h"
#include "../Flood/Utilities/LinearRegressionAnalysis.h"

// Perceptron includes

#include "../Flood/Perceptron/Perceptron.h"

// Multilayer perceptron includes

#include "../Flood/MultilayerPerceptron/MultilayerPerceptron.h"

// Objective functional includes

// Data modeling problems

#include "../Flood/ObjectiveFunctional/SumSquaredError.h"
#include "../Flood/ObjectiveFunctional/MeanSquaredError.h"
#include "../Flood/ObjectiveFunctional/RootMeanSquaredError.h"
#include "../Flood/ObjectiveFunctional/NormalizedSquaredError.h"
#include "../Flood/ObjectiveFunctional/MinkowskiError.h"
#include "../Flood/ObjectiveFunctional/RegularizedMinkowskiError.h"

// Clasical problems

#include "../Flood/ObjectiveFunctional/GeodesicProblem.h"
#include "../Flood/ObjectiveFunctional/BrachistochroneProblem.h"
#include "../Flood/ObjectiveFunctional/CatenaryProblem.h"
#include "../Flood/ObjectiveFunctional/IsoperimetricProblem.h"

// Optimal control problems

#include "../Flood/ObjectiveFunctional/CarProblem.h"

// Optimal shape design problems

#include "../Flood/ObjectiveFunctional/MinimumDragProblem.h"

// Function optimization problems 

#include "../Flood/ObjectiveFunctional/DeJongFunction.h"
#include "../Flood/ObjectiveFunctional/RosenbrockFunction.h"
#include "../Flood/ObjectiveFunctional/RastriginFunction.h"
#include "../Flood/ObjectiveFunctional/PlaneCylinder.h"

// Training algorithm includes

#include "../Flood/TrainingAlgorithm/GradientDescent.h"
#include "../Flood/TrainingAlgorithm/ConjugateGradient.h"
#include "../Flood/TrainingAlgorithm/QuasiNewtonMethod.h"
#include "../Flood/TrainingAlgorithm/RandomSearch.h"
#include "../Flood/TrainingAlgorithm/EvolutionaryAlgorithm.h"

using namespace Flood;

int main(void)
{
   std::cout << std::endl
             << "Flood Neural Network. Template Application." << std::endl;   

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
