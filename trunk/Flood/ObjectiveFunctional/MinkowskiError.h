/****************************************************************************************************************/
/*                                                                                                              */
/*   Flood: An Open Source Neural Networks C++ Library                                                          */
/*   www.cimne.com/flood                                                                                        */
/*                                                                                                              */
/*   M I N K O W S K I   E R R O R   C L A S S   H E A D E R                                                    */
/*                                                                                                              */
/*   Roberto Lopez                                                                                              */
/*   International Center for Numerical Methods in Engineering (CIMNE)                                          */
/*   Technical University of Catalonia (UPC)                                                                    */
/*   Barcelona, Spain                                                                                           */
/*   E-mail: rlopez@cimne.upc.edu                                                                               */
/*                                                                                                              */
/****************************************************************************************************************/


#ifndef __MINKOWSKIERROR_H__
#define __MINKOWSKIERROR_H__

#include "ObjectiveFunctional.h"
#include "../Utilities/InputTargetDataSet.h"

namespace Flood
{

/// This class represents the Minkowski error objective functional of a multilayer perceptron. This objective 
/// functional is used in data modeling problems.
///
/// @see MultilayerPerceptron.
/// @see InputTargetDataSet.

class MinkowskiError : public ObjectiveFunctional
{

private:

   /// Pointer to an input-target data set object.

   InputTargetDataSet* inputTargetDataSet;

   /// Minkowski exponent value.

   double minkowskiParameter;

public:

   // GENERAL CONSTRUCTOR

   MinkowskiError(MultilayerPerceptron*, InputTargetDataSet*);


   // DEFAULT CONSTRUCTOR

   MinkowskiError(void);


   // DESTRUCTOR

   virtual ~MinkowskiError(void);


   // METHODS

   // Get methods

   InputTargetDataSet* getInputTargetDataSet(void);

   double getMinkowskiParameter(void);

   // Set methods

   void setInputTargetDataSet(InputTargetDataSet*);

   void setMinkowskiParameter(double);

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
