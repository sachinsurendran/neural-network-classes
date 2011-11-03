/****************************************************************************************************************/
/*                                                                                                              */ 
/*   Flood: An Open Source Neural Networks C++ Library                                                          */
/*   www.cimne.com/flood                                                                                        */
/*                                                                                                              */
/*   V E C T O R   A N D   M A T R I X   A P P L I C A T I O N                                                  */
/*                                                                                                              */
/*   Roberto Lopez                                                                                              */ 
/*   International Center for Numerical Methods in Engineering (CIMNE)                                          */
/*   Technical University of Catalonia (UPC)                                                                    */
/*   Barcelona, Spain                                                                                           */
/*   E-mail: rlopez@cimne.upc.edu                                                                               */ 
/*                                                                                                              */  
/****************************************************************************************************************/

/// This application is an usage example of the Vector and Matrix classes in Flood.

#include <iostream>

// Utilities includes

#include "../Flood/Utilities/Vector.h"
#include "../Flood/Utilities/Matrix.h"

using namespace Flood;

int main(void)
{
   std::cout << std::endl
             << "Flood Neural Network. Vector and Matrix Application." 
             << std::endl;

   // Construct an empty vector of integers
     
   Vector<int> vector1;
  
   // Construct a vector of integers with 3 elements
     
   Vector<int> vector2(3);

    // Initialize all the elements of vector2 to 1
   
   vector2[0] = 1;
   vector2[1] = 1;
   vector2[2] = 1; 
   
   // Construct a vector of integers with 5 elements and initialize them to 0
     
   Vector<int> vector3(5, 1);

   // Construct a vector which is a copy of vector2
      
   Vector<int> vector4 = vector2;

   // Get size of vector4
   
   int size4 = vector4.getSize();       

   std::cout << std::endl 
             << "Size of vector 4:" << std::endl
             << size4 << std::endl;

   // Print the elements of vector4 to the screen
   
   std::cout << std::endl 
             << "Vector 4:" << std::endl
             << vector4 << std::endl;
     
   // Construct an empty matrix of double precision numbers
     
   Matrix<double> matrix1;

   // Construct a matrix of double precision numbers with 2 rows and 3 columns
     
   Matrix<double> matrix2(2,3);

   // Initialize all the elements of matrix2 to 1
   
   matrix2[0][0] = 1.0;
   matrix2[0][1] = 1.0;
   matrix2[0][2] = 1.0;
   matrix2[1][0] = 1.0; 
   matrix2[1][1] = 1.0;
   matrix2[1][2] = 1.0;
    
   // Construct a matrix of double precision numbers with 4 rows and 2 columns, and initialize all the elements 
   // to 1 
      
   Matrix<double> matrix3(4,2,1.0);

   // Construct a matrix which is a copy of matrix2
      
   Matrix<double> matrix4 = matrix2;
   
   // Get number of rows and columns of matrix4
    
   int numberOfRows4 = matrix4.getNumberOfRows();
   int numberOfColumns4 = matrix4.getNumberOfColumns();

   std::cout << std::endl 
             << "Number of rows of matrix 4:" << std::endl
             << numberOfRows4 << std::endl
             << "Number of columns of matrix 4:" << std::endl
             << numberOfColumns4 << std::endl;

   // Print the elements of matrix4 to the screen
   
   std::cout << std::endl 
             << "Matrix 4:" << std::endl
             << matrix4 << std::endl;
                  
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
