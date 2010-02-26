/****************************************************************************************************************/
/*                                                                                                              */ 
/*   Flood: An Open Source Neural Networks C++ Library                                                          */
/*   www.cimne.com/flood                                                                                        */
/*                                                                                                              */
/*   I N P U T   T A R G E T   D A T A   S E T   A P P L I C A T I O N                                          */
/*                                                                                                              */
/*   Roberto Lopez                                                                                              */ 
/*   International Center for Numerical Methods in Engineering (CIMNE)                                          */
/*   Technical University of Catalonia (UPC)                                                                    */
/*   Barcelona, Spain                                                                                           */
/*   E-mail: rlopez@cimne.upc.edu                                                                               */ 
/*                                                                                                              */  
/****************************************************************************************************************/

/// This application is an usage example of the InputTargetDataSet class in Flood.

// System includes

#include <iostream>

// Utilities includes

#include "../Flood/Utilities/Vector.h"
#include "../Flood/Utilities/Matrix.h"

#include "../Flood/Utilities/InputTargetDataSet.h"

using namespace Flood;

int main(void)
{
   std::cout << std::endl
             << "Flood Neural Network. Input Target Data Set Application." 
             << std::endl;	

   // Construct empty input-target data set object

   InputTargetDataSet inputTargetDataSet;

   // Load input-target data set object from file	

   inputTargetDataSet.load("../Data/InputTargetDataSet/InputTargetDataSet.dat");

   // Print input-target data set to screen

   std::cout << std::endl
             << "Input-target data set:" << std::endl;

   inputTargetDataSet.print();

   // Print input-target data set statistics

   inputTargetDataSet.printAllStatistics();

   // Preprocess the input and the target data for mean zero and standard deviation one

   inputTargetDataSet.preprocessMeanAndStandardDeviation();

   // Print the new input target data set to the screen

   std::cout << std::endl
             << "Preprocessed input-target data set:" << std::endl;
   
   inputTargetDataSet.print();	

   // Print new input-target data set statistics

   inputTargetDataSet.printAllStatistics(); 

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
