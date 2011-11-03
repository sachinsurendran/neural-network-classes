/****************************************************************************************************************/
/*                                                                                                              */ 
/*   Flood: An Open Source Neural Networks C++ Library                                                          */
/*   www.cimne.com/flood                                                                                        */
/*                                                                                                              */
/*   L I N E A R   R E G R E S S I O N   A N A L Y S I S   A P P L I C A T I O N                                */ 
/*                                                                                                              */
/*   Roberto Lopez                                                                                              */ 
/*   International Center for Numerical Methods in Engineering (CIMNE)                                          */
/*   Technical University of Catalonia (UPC)                                                                    */
/*   Barcelona, Spain                                                                                           */
/*   E-mail: rlopez@cimne.upc.edu                                                                               */ 
/*                                                                                                              */  
/****************************************************************************************************************/

/// This application is an usage example of the LinearRegressionAnalysis class in Flood.

// System includes

#include <iostream>
#include <fstream>
#include <time.h>

// Utilities includes

#include "../Flood/Utilities/Vector.h"
#include "../Flood/Utilities/Matrix.h"

#include "../Flood/Utilities/InputTargetDataSet.h"
#include "../Flood/Utilities/LinearRegressionAnalysis.h"

// Multilayer perceptron includes

#include "../Flood/MultilayerPerceptron/MultilayerPerceptron.h"

using namespace Flood;

int main(void)
{
   std::cout << std::endl
             << "Flood Neural Network. Linear Regression Analysis Application." 
             << std::endl;	

   srand((unsigned)time(NULL));

   // Multilayer perceptron object

   MultilayerPerceptron multilayerPerceptron(1,3,1);
   
   // Input-target data set object
 
   InputTargetDataSet inputTargetDataSet;

   inputTargetDataSet.load("../Data/LinearRegressionAnalysis/InputTargetDataSet.dat");

   // Linear regression analysis object
 
   LinearRegressionAnalysis linearRegressionAnalysis(&multilayerPerceptron, &inputTargetDataSet);
   
   linearRegressionAnalysis.saveResults("../Data/LinearRegressionAnalysis/LinearRegressionAnalysis.dat");   

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
