/****************************************************************************************************************/
/*                                                                                                              */
/*   Flood: An Open Source Neural Networks C++ Library                                                          */
/*   www.cimne.com/flood                                                                                        */
/*                                                                                                              */
/*   R O S E N B R O C K   F U N C T I O N   C L A S S   H E A D E R                                            */
/*                                                                                                              */
/*   Roberto Lopez                                                                                              */
/*   International Center for Numerical Methods in Engineering (CIMNE)                                          */
/*   Technical University of Catalonia (UPC)                                                                    */
/*   Barcelona, Spain                                                                                           */
/*   E-mail: rlopez@cimne.upc.edu                                                                               */
/*                                                                                                              */
/****************************************************************************************************************/


#ifndef __ROSENBROCKFUNCTION_H__
#define __ROSENBROCKFUNCTION_H__

#include "ObjectiveFunctional.h"

namespace Flood
{

/// This class represents the Rosenbrock's objective function.
/// Rosenbrock's valley is a classic optimization problem, also known as Banana function. 
/// The global optimum is inside a long, narrow, parabolic shaped flat valley. 
/// To find the valley is trivial, however convergence to the global optimum is difficult and hence this problem 
/// has been repeatedly used in assess the performance of optimization algorithms.
/// No analytical expressions for the gradient and the Hessian are provided, so therefore they are computed by 
/// means of numerical differentiation. 
/// For more information please visit
/// GEATbx - The Genetic and Evolutionary Algorithm Toolbox for Matlab, www.geatbx.com
/// 
/// @see ObjectiveFunctional.

class RosenbrockFunction : public ObjectiveFunctional
{

private:

   /// Number of variables in the Rosenbrock's function.

   int numberOfVariables;

public:

   // GENERAL CONSTRUCTOR

   RosenbrockFunction(MultilayerPerceptron*);

   // DEFAULT CONSTRUCTOR

   RosenbrockFunction(void);

   // DESTRUCTOR

   virtual ~RosenbrockFunction(void);


   // METHODS

   // Get methods

   int getNumberOfVariables(void);

   // Set methods

   void setNumberOfVariables(int);

   // Evaluation methods

   double calculateEvaluation(void);
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
