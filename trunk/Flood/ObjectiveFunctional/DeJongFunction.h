/****************************************************************************************************************/
/*                                                                                                              */
/*   Flood: An Open Source Neural Networks C++ Library                                                          */
/*   www.cimne.com/flood                                                                                        */
/*                                                                                                              */
/*   D E   J O N G   F U N C T I O N   C L A S S   H E A D E R                                                  */
/*                                                                                                              */
/*   Roberto Lopez                                                                                              */
/*   International Center for Numerical Methods in Engineering (CIMNE)                                          */
/*   Technical University of Catalonia (UPC)                                                                    */
/*   Barcelona, Spain                                                                                           */
/*   E-mail: rlopez@cimne.upc.edu                                                                               */
/*                                                                                                              */
/****************************************************************************************************************/


#ifndef __DEJONGFUNCTION_H__
#define __DEJONGFUNCTION_H__

#include "ObjectiveFunctional.h"

namespace Flood
{

/// This class represents the De Jong's objective function.
/// The simplest function for optimization is De Jong's function. It is continuos, convex and unimodal.
/// Explicit expressions for the gradient and the Hessian are provided here, and numerical differentiation is not 
/// needed. 
/// For more information please visit
/// GEATbx - The Genetic and Evolutionary Algorithm Toolbox for Matlab, www.geatbx.com
///
/// @see ObjectiveFunctional.

class DeJongFunction : public ObjectiveFunctional
{

private:

   /// Number of variables in the De Jong's Function. 

   int numberOfVariables;

public:

   // GENERAL CONSTRUCTOR

   DeJongFunction(MultilayerPerceptron*);


   // DEFAULT CONSTRUCTOR

   DeJongFunction(void);


   // DESTRUCTOR

   virtual ~DeJongFunction(void);


   // METHODS

   // Get methods

   int getNumberOfVariables(void);

   // Set methods

   void setNumberOfVariables(int);

   // Objective function methods

   double calculateEvaluation(void);

   // Objective function gradient methods

   Vector<double> calculateGradient(void);

   // Objective function Hessian methods

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
