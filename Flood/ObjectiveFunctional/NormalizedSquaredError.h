/****************************************************************************************************************/
/*                                                                                                              */
/*   Flood: An Open Source Neural Networks C++ Library                                                          */
/*   www.cimne.com/flood                                                                                        */
/*                                                                                                              */
/*   N O R M A L I Z E D   S Q U A R E D   E R R O R   C L A S S   H E A D E R                                  */
/*                                                                                                              */
/*   Roberto Lopez                                                                                              */ 
/*   International Center for Numerical Methods in Engineering (CIMNE)                                          */
/*   Technical University of Catalonia (UPC)                                                                    */
/*   Barcelona, Spain                                                                                           */
/*   E-mail: rlopez@cimne.upc.edu                                                                               */
/*                                                                                                              */
/****************************************************************************************************************/


#ifndef __NORMALIZEDSQUAREDERROR_H__
#define __NORMALIZEDSQUAREDERROR_H__

#include "ObjectiveFunctional.h"
#include "../Utilities/InputTargetDataSet.h"

namespace Flood
{

/// This class represents the normalized squared error objective functional of a multilayer perceptron. This 
/// objective functional is used in data modeling problems. If it has a value of unity then the neural network
/// predicting the data "in the mean" while a value of zero means perfect prediction of data.
///
/// @see MultilayerPerceptron.
/// @see InputTargetDataSet.

class NormalizedSquaredError : public ObjectiveFunctional
{

private:

   /// Pointer to an input-target data set object.

   InputTargetDataSet* inputTargetDataSet;

   /// Mean values of all the target variables. 

   Vector<double> meanOfTargetData;

public:

   // GENERAL CONSTRUCTOR

   NormalizedSquaredError(MultilayerPerceptron*, InputTargetDataSet*);


   // DEFAULT CONSTRUCTOR

   NormalizedSquaredError(void);


   // DESTRUCTOR

   virtual ~NormalizedSquaredError(void);


   // METHODS

   // Get methods

   InputTargetDataSet* getInputTargetDataSet(void);

   // Set methods

   void setInputTargetDataSet(InputTargetDataSet*);

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
