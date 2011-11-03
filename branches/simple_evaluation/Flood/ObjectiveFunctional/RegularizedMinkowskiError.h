/****************************************************************************************************************/
/*                                                                                                              */
/*   Flood: An Open Source Neural Networks C++ Library                                                          */
/*   www.cimne.com/flood                                                                                        */
/*                                                                                                              */
/*   R E G U L A R I Z E D   M E A N   S Q U A R E D   E R R O R   C L A S S   H E A D E R                      */
/*                                                                                                              */
/*   Roberto Lopez                                                                                              */
/*   International Center for Numerical Methods in Engineering (CIMNE)                                          */
/*   Technical University of Catalonia (UPC)                                                                    */
/*   Barcelona, Spain                                                                                           */
/*   E-mail: rlopez@cimne.upc.edu                                                                               */
/*                                                                                                              */
/****************************************************************************************************************/


#ifndef __REGULARIZEDMINKOWSKIERROR_H__
#define __REGULARIZEDMINKOWSKI_H__

#include "ObjectiveFunctional.h"
#include "../Utilities/InputTargetDataSet.h"

namespace Flood
{

/// This class represents the mean squared error with regularization objective functional of a multilayer 
/// perceptron. The regularizer form used here is free parameter decay. 
/// This objective functional is used in data modeling problems.
///
/// @see MultilayerPerceptron.
/// @see InputTargetDataSet.
///
/// @todo This class is actually implementing a regularized mean squared error. 
/// Changing the mean squared error by a Minkowski one is required. 

class RegularizedMinkowskiError : public ObjectiveFunctional
{

private:

   /// Pointer to an input-target data set object.

   InputTargetDataSet* inputTargetDataSet;

   /// Mean squared error term weight.

   double meanSquaredErrorWeight;

   /// Regularization term weight.

   double regularizationWeight;
   

public:

   // GENERAL CONSTRUCTOR

   RegularizedMinkowskiError(MultilayerPerceptron*, InputTargetDataSet*);


   // DEFAULT CONSTRUCTOR

   RegularizedMinkowskiError(void);


   // DESTRUCTOR

   virtual ~RegularizedMinkowskiError(void);


   // METHODS

   // Get methods

   InputTargetDataSet* getInputTargetDataSet(void);

   double getMeanSquaredErrorWeight(void);
   double getRegularizationWeight(void);

   // Set methods

   void setInputTargetDataSet(InputTargetDataSet*);

   void setMeanSquaredErrorWeight(double);
   void setRegularizationWeight(double);


   // Objective functional evaluation methods

   double calculateEvaluation(void);

   // Objective function gradient vector methods

   Vector<double> calculateGradient(void);

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
