/****************************************************************************************************************/
/*                                                                                                              */
/*   Flood: An Open Source Neural Networks C++ Library                                                          */
/*   www.cimne.com/flood                                                                                        */
/*                                                                                                              */
/*   D E   J O N G   F U N C T I O N   C L A S S                                                                */
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

#include "DeJongFunction.h"

namespace Flood
{

// GENERAL CONSTRUCTOR

/// General constructor. It creates a De Jong's objective function object.
/// It also initializes all the rest of class members to their default values:
///
/// <ul>
/// <li> Number of variables (independent parameters): 2.
/// <li> Name of independent parameters: x0,...,xn
/// <li> Minimum of independent parameters: -5.12,...,-5.12.
/// <li> Maximum of independent parameters: 5.12,...,5.12.
/// <li> Pre and postprocessing method: Minimum and maximum.
/// </ul> 

DeJongFunction::DeJongFunction(MultilayerPerceptron* newMultilayerPerceptron)       
: ObjectiveFunctional(newMultilayerPerceptron)
{
   // Control sentence

   int numberOfInputs = multilayerPerceptron->getNumberOfInputs();
   int numberOfOutputs = multilayerPerceptron->getNumberOfOutputs();

   if(numberOfInputs != 0 || numberOfOutputs != 0)
   {
      std::cerr << std::endl
                << "Flood Error: DeJongFunction class." << std::endl
                << "DeJongFunction(MultilayerPerceptron*) constructor." << std::endl
                << "Number of inputs and outputs in multilayer perceptron must be 0." << std::endl
                << std::endl;

      exit(1);
   }

   numberOfVariables = 2;

}


// DESTRUCTOR

/// Destructor.

DeJongFunction::~DeJongFunction(void)
{

}


// METHODS

// int getNumberOfVariables(void) method

/// This method returns the number of variables in the De Jong's function. 
/// This is the number of independent parameters in the multilayer perceptron. 

int DeJongFunction::getNumberOfVariables(void)
{
   return(numberOfVariables);
}


// void setNumberOfVariables(int) method

/// This method sets a new number of variables in the De Jong's function. 
/// This also sets a new number of independent parameters in the multilayer perceptron. 
/// It also initializes all the rest of class members to their default values:
///
/// <ul>
/// <li> Name of independent parameters: x0,...,xn
/// <li> Minimum of independent parameters: -5.12,...,-5.12.
/// <li> Maximum of independent parameters: 5.12,...,5.12.
/// </ul> 
///
/// @param newNumberOfVariables Number of variables in the De Jong's function. 

void DeJongFunction::setNumberOfVariables(int newNumberOfVariables)
{
   numberOfVariables = newNumberOfVariables;

   // Set multilayer perceptron stuff

   multilayerPerceptron->setNumberOfIndependentParameters(numberOfVariables);
   
   for(int i = 0; i < numberOfVariables; i++)
   {
      std::stringstream buffer;
      buffer << "x" << i+1;
      multilayerPerceptron->setNameOfSingleIndependentParameter(i, buffer.str());

      multilayerPerceptron->setMinimumOfSingleIndependentParameter(i, -5.12);
      multilayerPerceptron->setMaximumOfSingleIndependentParameter(i, 5.12);  
   }

   multilayerPerceptron->initIndependentParametersUniform(-5.12, 5.12);
}


// double calculateEvaluation(void) method

/// This method returns the De Jong's function evaluation for the actual independent parameters in the multilayer 
/// perceptron.
/// 
/// @see calculateGradient(void).
/// @see calculateHessian(void).
/// @see calculateInverseHessian(void).

double DeJongFunction::calculateEvaluation(void)
{
    // Control sentence 

   if(multilayerPerceptron == NULL)
   {
      std::cerr << std::endl
                << "Flood Error: DeJongFunction class." << std::endl
                << "double calculateEvaluation(void) method." << std::endl
                << "Pointer to multilayer perceptron object cannot be NULL." << std::endl
                << std::endl;

        exit(1);
   }

   numberOfEvaluations++;

   Vector<double> argument = multilayerPerceptron->getIndependentParameters();

   double evaluation = argument.dot(argument);

   return(evaluation);
}


// Vector<double> calculateGradient(void) method

/// This method returns the De Jong's analytical gradient vector for the actual independent parameter values in 
/// the multilayer perceptron.
/// 
/// @see calculateEvaluation(void).
/// @see calculateHessian(void).
/// @see calculateInverseHessian(void).

Vector<double> DeJongFunction::calculateGradient(void)
{
   Vector<double> argument = multilayerPerceptron->getIndependentParameters();

   Vector<double> gradient = argument*2.0;

   return(gradient);
}


// Matrix<double> calculateHessian(void) method

/// This method returns the De Jong's analytical Hessian matrix for the actual independent parameter values in 
/// the multilayer perceptron.
/// 
/// @see calculateEvaluation(void).
/// @see calculateGradient(void).
/// @see calculateInverseHessian(void).

Matrix<double> DeJongFunction::calculateHessian(void)
{
   Matrix<double> hessian(numberOfVariables, numberOfVariables);

   Vector<double> argument = multilayerPerceptron->getIndependentParameters();

   for(int i = 0; i < numberOfVariables; i++)
   {
      for(int j = 0; j < numberOfVariables; j++)
      {
         if(i == j)
         {
            hessian[i][j] = 2.0;
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

/// This method returns the De Jong's analytical inverse Hessian matrix for the actual independent parameter 
/// values in the multilayer perceptron.
/// 
/// @see calculateEvaluation(void).
/// @see calculateGradient(void).
/// @see calculateHessian(void).

Matrix<double> DeJongFunction::calculateInverseHessian(void)
{
   Matrix<double> inverseHessian(numberOfVariables, numberOfVariables);

   for(int i = 0; i < numberOfVariables; i++)
   {
      for(int j = 0; j < numberOfVariables; j++)
      {
         if(i == j)
         {
            inverseHessian[i][j] = 0.5;
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
