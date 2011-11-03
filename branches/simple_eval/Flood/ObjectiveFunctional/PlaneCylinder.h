/****************************************************************************************************************/
/*                                                                                                              */
/*   Flood: An Open Source Neural Networks C++ Library                                                          */
/*   www.cimne.com/flood                                                                                        */
/*                                                                                                              */
/*   P L A N E   C Y L I N D E R   C L A S S   H E A D E R                                                      */
/*                                                                                                              */
/*   Roberto Lopez                                                                                              */
/*   International Center for Numerical Methods in Engineering (CIMNE)                                          */
/*   Technical University of Catalonia (UPC)                                                                    */
/*   Barcelona, Spain                                                                                           */
/*   E-mail: rlopez@cimne.upc.edu                                                                               */
/*                                                                                                              */
/****************************************************************************************************************/


#ifndef __PLANECYLINDER_H__
#define __PLANECYLINDER_H__

#include "ObjectiveFunctional.h"

namespace Flood
{

/// This class represents the objective function for the plane-cylinder problem.
///
/// @see ObjectiveFunctional.

class PlaneCylinder : public ObjectiveFunctional
{

private:

   /// Penalty term ratio.

   double penalty;

public:

   // GENERAL CONSTRUCTOR

   PlaneCylinder(MultilayerPerceptron*);


   // DEFAULT CONSTRUCTOR

   PlaneCylinder(void);


   // DESTRUCTOR

   virtual ~PlaneCylinder(void);


   // METHODS

   // Get methods

   double getPenalty(void);

   // Set methods

   void setPenalty(double);

   // Objective function methods

   double calculateError(void);

   double calculateEvaluation(void);
   
   // Utility methods
   
   void print(void);
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
