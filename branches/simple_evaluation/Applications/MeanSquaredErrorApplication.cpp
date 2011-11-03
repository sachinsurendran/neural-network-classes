/****************************************************************************************************************/
/*                                                                                                              */ 
/*   Flood: An Open Source Neural Networks C++ Library                                                          */
/*   www.cimne.com/flood                                                                                        */
/*                                                                                                              */
/*   M E A N   S Q U A R E D   E R R O R   A P P L I C A T I O N                                                */
/*                                                                                                              */
/*   Roberto Lopez                                                                                              */ 
/*   International Center for Numerical Methods in Engineering (CIMNE)                                          */
/*   Technical University of Catalonia (UPC)                                                                    */
/*   Barcelona, Spain                                                                                           */
/*   E-mail: rlopez@cimne.upc.edu                                                                               */ 
/*                                                                                                              */  
/****************************************************************************************************************/

/// This application is an usage example of the MeanSquaredError class in Flood.

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

using namespace Flood;

int main(void)
{
   std::cout << std::endl
             << "Flood Neural Network. Mean Squared Error Application." 
             << std::endl;	

   srand((unsigned)time(NULL));

   // Input-target data set object

   InputTargetDataSet inputTargetDataSet;

   inputTargetDataSet.load("../Data/MeanSquaredError/InputTargetDataSet.dat");

   // Multilayer perceptron object 

   MultilayerPerceptron multilayerPerceptron(1,3,1);

   // Construct a mean squared error object associated to the multilayer perceptron and the input-target data set 
   /// objects 

   MeanSquaredError meanSquaredError(&multilayerPerceptron, &inputTargetDataSet);

   // Calculate the mean squared error between the multilayer perceptron and the input-target data set

   double evaluation = meanSquaredError.calculateEvaluation();

   std::cout << std::endl
             << "Evaluation:" << std::endl
             << evaluation << std::endl;

   // Obtain the mean squared error function gradient 

   Vector<double> gradient = meanSquaredError.calculateGradient();

   std::cout << std::endl 
             << "Gradient:" << std::endl
             << gradient << std::endl;
   
   // Get the gradient norm

   double gradientNorm = gradient.calculateNorm();	
	
   std::cout << std::endl
             << "Gradient norm:" << std::endl
             << gradientNorm << std::endl;

   std::cout << std::endl;

   return(0);
}  


// Flood: An Open Source Neural Networks C++ Library.
// Copyright (C) 2005-2007 Roberto Lopez 
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
