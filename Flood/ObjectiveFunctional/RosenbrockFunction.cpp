/****************************************************************************************************************/
/*                                                                                                              */
/*   Flood: An Open Source Neural Networks C++ Library                                                          */
/*   www.cimne.com/flood                                                                                        */
/*                                                                                                              */
/*   R O S E N B R O C K   F U N C T I O N   C L A S S                                                          */
/*                                                                                                              */
/*   Roberto Lopez                                                                                              */
/*   International Center for Numerical Methods in Engineering (CIMNE)                                          */ 
/*   Technical University of Catalonia (UPC)                                                                    */
/*   Barcelona, Spain                                                                                           */
/*   E-mail: rlopez@cimne.upc.edu                                                                               */
/*                                                                                                              */
/****************************************************************************************************************/

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <math.h>   

#include "RosenbrockFunction.h"

namespace Flood
{

// GENERAL CONSTRUCTOR

/// General constructor. It creates a Rosenbrock's objective function object.
/// It also initializes all the rest of class members to their default values:
///
/// <ul>
/// <li> Number of variables (independent parameters):  2.
/// <li> Name of independent parameters: x1,...,xn.
/// <li> Minimum of independent parameters: -2.048,...,-2.048.
/// <li> Maximum of independent parameters: 2.048,...,2.048.
/// <li> Pre and postprocessing method: Minimum and maximum.
/// </ul> 

RosenbrockFunction::RosenbrockFunction(MultilayerPerceptron* newMultilayerPerceptron)       
: ObjectiveFunctional(newMultilayerPerceptron)
{
   numberOfVariables = 2;
}


// DESTRUCTOR

/// Destructor.

RosenbrockFunction::~RosenbrockFunction(void)
{

}


// METHODS

// int getNumberOfVariables(void) method

/// This method returns the number of variables in the Rosenbrock's function. 
/// This is the number of independent parameters in the multilayer perceptron. 

int RosenbrockFunction::getNumberOfVariables(void)
{
   return(numberOfVariables);
}


// void setNumberOfVariables(int) method

/// This method sets a new number of variables in the Rosenbrock's function. 
/// This also sets a new number of independent parameters in the multilayer perceptron. 
/// It also initializes all the rest of class members to their default values:
///
/// <ul>
/// <li> Name of independent parameters: x0,...,xn
/// <li> Minimum of independent parameters: -2.048,...,-2.048.
/// <li> Maximum of independent parameters: 2.048,...,2.048.
/// </ul> 
///
/// @param newNumberOfVariables Number of variables in the Rosenbrock's function. 

void RosenbrockFunction::setNumberOfVariables(int newNumberOfVariables)
{
   numberOfVariables = newNumberOfVariables;
}


// double calculateEvaluation(void) method

/// This method returns the Rosenbrock's function evaluation for the actual independent parameters in the 
/// multilayer perceptron.

double RosenbrockFunction::calculateEvaluation(void)
{
   // Control sentence (if debug)

   #ifndef NDEBUG 

   if(multilayerPerceptron == NULL)
   {
      std::cerr << std::endl
                << "Flood Error: RosenbrockFunction class." << std::endl
                << "double calculateEvaluation(void) method." << std::endl
                << "Pointer to multilayer perceptron object cannot be NULL." << std::endl
                << std::endl;

      exit(1);
   }

   #endif

   numberOfEvaluations++;

   double evaluation = 0.0;

   Vector<double> argument = multilayerPerceptron->getIndependentParameters();

   // Get evaluation

   for(int i = 0; i < numberOfVariables-1; i++)
   {
      evaluation += 100.0*pow(argument[i+1] - pow(argument[i],2), 2) + pow(1.0 - argument[i], 2);
   }

   return(evaluation);
}

}


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
