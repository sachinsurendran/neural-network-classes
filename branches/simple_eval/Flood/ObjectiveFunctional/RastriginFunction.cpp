/****************************************************************************************************************/
/*                                                                                                              */
/*   Flood: An Open Source Neural Networks C++ Library                                                          */
/*   www.cimne.com/flood                                                                                        */
/*                                                                                                              */
/*   R A S T R I G I N   F U N C T I O N   C L A S S                                                            */
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

#include "RastriginFunction.h"

namespace Flood
{

// GENERAL CONSTRUCTOR

/// General constructor. It creates a Rastrigin's objective function object.
/// It also initializes all the rest of class members to their default values:
///
/// <ul>
/// <li> Number of variables (independent parameters): 2.
/// <li> Minimum of independent parameters: -5.12,...,-5.12.
/// <li> Maximum of independent parameters: 5.12,...,5.12.
/// </ul> 

RastriginFunction::RastriginFunction(MultilayerPerceptron* newMultilayerPerceptron)       
: ObjectiveFunctional(newMultilayerPerceptron)
{
   numberOfVariables = 2;
}


// DESTRUCTOR

/// Destructor.

RastriginFunction::~RastriginFunction(void)
{

}


// METHODS


// int getNumberOfVariables(void) method

/// This method returns the number of variables in the Rastrigin's function. 
/// This is the number of independent parameters in the multilayer perceptron. 

int RastriginFunction::getNumberOfVariables(void)
{
   return(numberOfVariables);
}


// void setNumberOfVariables(int) method

/// This method sets a new number of variables in the Rastrigin's function. 
/// This also sets a new number of independent parameters in the multilayer perceptron. 
/// It also initializes all the rest of class members to their default values:
///
/// <ul>
/// <li> Name of independent parameters: x0,...,xn
/// <li> Minimum of independent parameters: -5.12,...,-5.12.
/// <li> Maximum of independent parameters: 5.12,...,5.12.
/// </ul> 
///
/// @param newNumberOfVariables Number of variables in the Rastrigin's function. 

void RastriginFunction::setNumberOfVariables(int newNumberOfVariables)
{
   numberOfVariables = newNumberOfVariables;
}


// double calculateEvaluation(void) method

/// This method returns the Rastrigin's function evaluation.
///
/// @see calculateGradient(void).
/// @see calculateHessian(void).
/// @see calculateInverseHessian(void).

double RastriginFunction::calculateEvaluation(void)
{
   // Control sentence (if debug)

   #ifndef NDEBUG 

   if(multilayerPerceptron == NULL)
   {
      std::cerr << std::endl
                << "Flood Error: RastriginFunction class." << std::endl
                << "double calculateEvaluation(void) method." << std::endl
                << "Pointer to multilayer perceptron object cannot be NULL." << std::endl
                << std::endl;

        exit(1);
   }

   #endif

   numberOfEvaluations++;

   double evaluation = 0.0;

   const double pi = 3.1415927;

   Vector<double> argument = multilayerPerceptron->getIndependentParameters();

   evaluation = 10.0*numberOfVariables;

   for(int i = 0; i < numberOfVariables; i++)
   {
      evaluation += pow(argument[i], 2) - 10.0*cos(2.0*pi*argument[i]);
   }

   return(evaluation);
}


// Vector<double> calculateGradient(void) method

/// This method returns the Rastrigin's analytical gradient vector. 
///
/// @see calculateEvaluation(void).
/// @see calculateHessian(void).
/// @see calculateInverseHessian(void).

Vector<double> RastriginFunction::calculateGradient(void)
{
   Vector<double> gradient(numberOfVariables);

   const double pi = 3.1415927;

   Vector<double> argument = multilayerPerceptron->getIndependentParameters();

   for(int i = 0; i < numberOfVariables; i++)
   {
      gradient[i] = 2.0*argument[i] + 10.0*sin(2.0*pi*argument[i])*2.0*pi;
   }

   return(gradient);
}


// Matrix<double> calculateHessian(void) method

/// This method returns the Rastrigin's analytical Hessian matrix. 
///
/// @see calculateEvaluation(void).
/// @see calculateGradient(void).
/// @see calculateInverseHessian(void).

Matrix<double> RastriginFunction::calculateHessian(void)
{
   Matrix<double> hessian(numberOfVariables, numberOfVariables);

   const double pi = 3.1415927;

   Vector<double> argument = multilayerPerceptron->getIndependentParameters();

   for(int i = 0; i < numberOfVariables; i++)
   {
      for(int j = 0; j < numberOfVariables; j++)
      {
         if(i == j)
         {
            hessian[i][j] = 2.0 + 10.0*cos(2.0*pi*argument[i])*4.0*pow(pi,2);
         }
         else
         {
            hessian[i][j] = 0.0;
         }
      }
   }

   return(hessian);
}


// Matrix<double> calculateInverseHessian(void) method

/// This method returns the Rastrigin's analytical inverse Hessian matrix. 
///
/// @see calculateEvaluation(void).
/// @see calculateGradient(void).
/// @see calculateHessian(void).

Matrix<double> RastriginFunction::calculateInverseHessian(void)
{
   Matrix<double> inverseHessian(numberOfVariables, numberOfVariables);

   const double pi = 3.1415927;

   Vector<double> argument = multilayerPerceptron->getIndependentParameters();

   for(int i = 0; i < numberOfVariables; i++)
   {
      for(int j = 0; j < numberOfVariables; j++)
      {
         if(i == j)
         {
            inverseHessian[i][j] = 1.0/(2.0 + 10.0*cos(2.0*pi*argument[i])*4.0*pow(pi,2));
         }
         else
         {
            inverseHessian[i][j] = 0.0;
         }
      }
   }

   return(inverseHessian);
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
