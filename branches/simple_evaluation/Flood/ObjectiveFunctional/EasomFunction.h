/******************************************************************************/
/*                                                                            */
/*   Flood: An Open Source Neural Networks C++ Library                        */
/*   www.cimne.com/flood                                                      */
/*                                                                            */
/*   E A S O M   F U N C T I O N   C L A S S   H E A D E R                    */
/*                                                                            */
/*   Gilles Cadose                                                            */
/*   Carlos Vargas de la Fuente                                               */
/*   Hebert Sotelo Aedo                                                       */
/*   International Center for Numerical Methods in Engineering (CIMNE)        */
/*   Technical University of Catalonia (UPC)                                  */
/*   Barcelona, Spain                                                         */
/*   E-mail: rlopez@cimne.upc.edu                                             */
/*                                                                            */
/******************************************************************************/


#ifndef __EASOMFUNCTION_H__
#define __EASOMFUNCTION_H__

#include "ObjectiveFunctional.h"

namespace Flood
{

/// This class represents the Easom's objective function.
/// The Easom function is an unimodal test function, 
///   where the global minimum has a small area relative to the search space. 
///   The function was inverted for minimization.
/// For more information please visit
///   GEATbx - The Genetic and Evolutionary Algorithm Toolbox for Matlab, www.geatbx.com
///
/// @see ObjectiveFunctional.

class EasomFunction : public ObjectiveFunctional
{

public:

   // GENERAL CONSTRUCTOR

   EasomFunction(MultilayerPerceptron*);


   // DEFAULT CONSTRUCTOR

   EasomFunction(void);


   // DESTRUCTOR

   virtual ~EasomFunction(void);


   // METHODS

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
