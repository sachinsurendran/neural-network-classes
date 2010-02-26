/******************************************************************************/
/*                                                                            */
/*   Flood: An Open Source Neural Networks C++ Library                        */
/*   www.cimne.com/flood                                                      */
/*                                                                            */
/*   R A N D O M   S E A R C H   C L A S S   H E A D E R                      */
/*                                                                            */ 
/*   Roberto Lopez                                                            */ 
/*   International Center for Numerical Methods in Engineering (CIMNE)        */
/*   Technical University of Catalonia (UPC)                                  */
/*   Barcelona, Spain                                                         */
/*   E-mail: rlopez@cimne.upc.edu                                             */ 
/*                                                                            */
/******************************************************************************/

#ifndef __RANDOMSEARCH_H__
#define __RANDOMSEARCH_H__

#include "TrainingAlgorithm.h"
#include "../ObjectiveFunctional/ObjectiveFunctional.h"
 
namespace Flood
{

/// This concrete class represents a random search training algorithm
/// for an objective functional of a multilayer perceptron.
///
/// @see TrainingAlgorithm.

class RandomSearch : public TrainingAlgorithm
{

private: 

   // FIELDS

   /// Minimum values for all the free parameters in the seach hypercube.
   
   Vector<double> minimumValueOfFreeParameters;
   
   /// Maximum values for all the free parameters in the seach hypercube.
   
   Vector<double> maximumValueOfFreeParameters;
   

public:

   // GENERAL CONSTRUCTOR

   RandomSearch(ObjectiveFunctional*); 


   // DEFAULT CONSTRUCTOR

   RandomSearch(void); 


   // DESTRUCTOR

   virtual ~RandomSearch(void);


   // METHODS

   // Get methods

   Vector<double> getMinimumValueOfFreeParameters(void);
   Vector<double> getMaximumValueOfFreeParameters(void);

   Matrix<double> getMinimumAndMaximumValuesOfFreeParameters(void);

   // Set methods

   void setMinimumValueOfFreeParameters(Vector<double>);
   void setMaximumValueOfFreeParameters(Vector<double>);

   void setMinimumAndMaximumValuesOfFreeParameters(Matrix<double>);

   // Train methods

   void train(void);

   // Utiltity methods

   void print(void);

   void load(char*);
   void save(char*);

   void resizeTrainingHistory(int);

   void setReserveAllTrainingHistory(bool);

   void saveTrainingHistory(char*);

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
