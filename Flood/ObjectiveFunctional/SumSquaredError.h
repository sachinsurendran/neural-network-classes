/****************************************************************************************************************/
/*                                                                                                              */
/*   Flood: An Open Source Neural Networks C++ Library                                                          */
/*   www.cimne.com/flood                                                                                        */
/*                                                                                                              */
/*   S U M   S Q U A R E D   E R R O R   C L A S S   H E A D E R                                                */
/*                                                                                                              */
/*   Roberto Lopez                                                                                              */
/*   International Center for Numerical Methods in Engineering (CIMNE)                                          */
/*   Technical University of Catalonia (UPC)                                                                    */
/*   Barcelona, Spain                                                                                           */
/*   E-mail: rlopez@cimne.upc.edu                                                                               */
/*                                                                                                              */
/****************************************************************************************************************/

#ifndef __SUMSQUAREDERROR_H__
#define __SUMSQUAREDERROR_H__

#include "ObjectiveFunctional.h"
#include "../Utilities/InputTargetDataSet.h"

namespace Flood
{

/// This class represents the sum squared error objective functional of a multilayer perceptron. This objective 
/// functional is used in data modeling problems, such as function regression, pattern recognition or time series 
/// prediction.
///
/// @see MultilayerPerceptron.
/// @see InputTargetDataSet.
/// @see MeanSquaredError.
/// @see RootMeanSquaredError. 
/// @see NormalizedSquaredError.
/// @see MinkowskiError.
/// @see RegularizedMeanSquaredError.

class SumSquaredError : public ObjectiveFunctional
{

private:

   /// Pointer to an input-target data set object.

   InputTargetDataSet* inputTargetDataSet;

public:

   // GENERAL CONSTRUCTOR

   SumSquaredError(MultilayerPerceptron*, InputTargetDataSet*);


   // DEFAULT CONSTRUCTOR

   SumSquaredError(void);


   // DESTRUCTOR

   virtual ~SumSquaredError(void);


   // METHODS

   // Get methods

   InputTargetDataSet* getInputTargetDataSet(void);

   // Set methods

   void setInputTargetDataSet(InputTargetDataSet*);

   // Objective functional evaluation methods

   double calculateEvaluation(void);

   Vector<double> calculateSquaredErrors(void);

   // Objective function gradient vector methods

   Vector<double> calculateGradient(void);

   Matrix<double> calculateJacobian(void);
   Matrix<double> calculateJacobianTest(void);

   // Utility methods

   void saveInputTargetAndOutput(char*);
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
