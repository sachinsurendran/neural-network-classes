/****************************************************************************************************************/
/*                                                                                                              */
/*   Flood: An Open Source Neural Networks C++ Library                                                          */
/*   www.cimne.com/flood                                                                                        */
/*                                                                                                              */
/*   R A S T R I G I N   F U N C T I O N   C L A S S   H E A D E R                                              */
/*                                                                                                              */
/*   Roberto Lopez                                                                                              */
/*   International Center for Numerical Methods in Engineering (CIMNE)                                          */
/*   Technical University of Catalonia (UPC)                                                                    */
/*   Barcelona, Spain                                                                                           */
/*   E-mail: rlopez@cimne.upc.edu                                                                               */
/*                                                                                                              */
/****************************************************************************************************************/


#ifndef __RASTRIGINFUNCTION_H__
#define __RASTRIGINFUNCTION_H__

#include "ObjectiveFunctional.h"

namespace Flood
{

/// This class represents the Rastrigin's objective function.
/// Rastrigin's function is based on De Jong's function with the addition of cosine modulation to produce many 
/// local minima. 
/// Thus, the test function is highly multimodal. However, the location of the minima are regularly distributed.
/// Analytical expressions for the gradient vector and the Hessian matrix are provided.
/// For more information please visit
/// GEATbx - The Genetic and Evolutionary Algorithm Toolbox for Matlab, www.geatbx.com
///
/// @see ObjectiveFunctional.

class RastriginFunction : public ObjectiveFunctional
{

private:

   /// Number of variables in the Rastrigin's function. 

   int numberOfVariables;

public:

   // GENERAL CONSTRUCTOR

   RastriginFunction(MultilayerPerceptron*);

   // DEFAULT CONSTRUCTOR

   RastriginFunction(void);

   // DESTRUCTOR

   virtual ~RastriginFunction(void);


   // METHODS

   // Get methods

   int getNumberOfVariables(void);

   // Set methods

   void setNumberOfVariables(int);

   // Objective function methods

   double calculateEvaluation(void);

   // Objective function gradient vector methods

   Vector<double> calculateGradient(void);

   // Objective function Hessian matrix methods

   Matrix<double> calculateHessian(void);
   
   Matrix<double> calculateInverseHessian(void);
};

}

#endif


// Flood: An Open Source Neural Networks C++ Library.
// Copyright (C) 2005-2008 Roberto Lopez 
//
// This library is free software; you can redistribute it and/or
// modify it under the s of the GNU Lesser General Public
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
