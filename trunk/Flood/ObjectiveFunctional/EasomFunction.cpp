/****************************************************************************************************************/
/*                                                                                                              */
/*   Flood: An Open Source Neural Networks C++ Library                                                          */
/*   www.cimne.com/flood                                                                                        */
/*                                                                                                              */
/*   E A S O M   F U N C T I O N   C L A S S                                                                    */
/*                                                                                                              */
/*   Gilles Cadose                                                                                              */
/*   Carlos Vargas de la Fuente                                                                                 */
/*   Hebert Sotelo Aedo                                                                                         */
/*   International Center for Numerical Methods in Engineering (CIMNE)                                          */
/*   Technical University of Catalonia (UPC)                                                                    */
/*   Barcelona, Spain                                                                                           */
/*                                                                                                              */
/****************************************************************************************************************/


#include <iostream>
#include <fstream>
#include <math.h>

#include "EasomFunction.h"

namespace Flood
{

// GENERAL CONSTRUCTOR

/// General constructor. It creates an Easom's objective function object.
/// It also initializes all the rest of class members to their default values:
///
/// <ul>
/// <li> Number of variables = 2.
/// <li> Lower bound = 1,...,1.
/// <li> Upper bound = 10,...,10.
/// </ul> 

EasomFunction::EasomFunction(MultilayerPerceptron* newMultilayerPerceptron)       
: ObjectiveFunctional(newMultilayerPerceptron)
{
   // Control sentences

   int numberOfInputs = multilayerPerceptron->getNumberOfInputs();
   int numberOfOutputs = multilayerPerceptron->getNumberOfOutputs();

   if(numberOfInputs != 0 || numberOfOutputs != 0)
   {
      std::cerr << std::endl
                << "Flood Error: EasomFunction class." << std::endl
                << "EasomFunction(MultilayerPerceptron*) constructor." << std::endl
                << "Number of inputs and outputs in multilayer perceptron must be 0." << std::endl
                << std::endl;

      exit(1);
   }

   int numberOfIndependentParameters = multilayerPerceptron->getNumberOfIndependentParameters();

   if(numberOfIndependentParameters != 2)
   {
      std::cerr << std::endl
                << "Flood Error: EasomFunction class." << std::endl
                << "EasomFunction(MultilayerPerceptron*) constructor." << std::endl
                << "Number of independent paramters in multilayer perceptron must be 2." << std::endl
                << std::endl;

      exit(1);
   }

}


// DESTRUCTOR

/// Destructor.

EasomFunction::~EasomFunction(void)
{

}


// METHODS

// double calculateEvaluation(void) method

/// This method returns the Easom's function evaluation for a given vector of independent parameters.

double EasomFunction::calculateEvaluation(void)
{
    // Control sentence 

   if(multilayerPerceptron == NULL)
   {
      std::cerr << std::endl
                << "Flood Error: EasomFunction class." << std::endl
                << "double calculateEvaluation(void) method." << std::endl
                << "Pointer to multilayer perceptron object cannot be NULL." << std::endl
                << std::endl;

      exit(1);
   }

   double evaluation = 0.0;

   numberOfEvaluations++;

   const double pi = 3.1415927;

   Vector<double> argument = multilayerPerceptron->getIndependentParameters();
         
   double a = exp(-(pow((argument[0]-pi), 2) + pow((argument[1]-pi), 2))); 
    
   evaluation = -cos(argument[0])*cos(argument[1])*a;
    
   return(evaluation);
}


// Vector<double> calculateGradient(void) method

/// This method returns the Easom's analytical gradient vector. 

Vector<double> EasomFunction::calculateGradient(void)
{
   const int numberOfVariables = 2;

   Vector<double> gradient(numberOfVariables);

   const double pi = 3.1415927;

   Vector<double> argument = multilayerPerceptron->getIndependentParameters();
   
   double a = exp(-(pow((argument[0]-pi), 2) + pow((argument[1]-pi), 2)));
   
   gradient[0] = a*(sin(argument[0])*cos(argument[1])+2*cos(argument[0])*cos(argument[1])*(argument[0]-pi));
   gradient[1] = a*(cos(argument[0])*sin(argument[1])+2*cos(argument[0])*cos(argument[1])*(argument[1]-pi));
   
   return(gradient);
}


// Matrix<double> calculateHessian(void) method

/// This method returns the Easom's analytical Hessian matrix. 

Matrix<double> EasomFunction::calculateHessian(void)
{ 
   const int numberOfVariables = 2;

   Matrix<double> hessian(numberOfVariables, numberOfVariables);

   const double pi = 3.1415927;

   Vector<double> argument = multilayerPerceptron->getIndependentParameters();

   double a = exp(-(pow((argument[0]-pi), 2) + pow((argument[1]-pi), 2)));
   double b = cos(argument[0])*cos(argument[1]);
   double c = sin(argument[0])*cos(argument[1]);
   double d = sin(argument[0])*sin(argument[1]);
   double e = cos(argument[0])*sin(argument[1]);
   
   hessian[0][0] = a*(3.0*b - 4.0*c*(argument[0] - pi) - 4.0*b*pow((argument[0]-pi), 2));
   hessian[0][1] = a*(-d    - 2.0*c*(argument[1] - pi) - 2.0*e*(argument[0]-pi) - 4.0*b*(argument[0]-pi)*(argument[1]-pi));
   hessian[1][0] = a*(-d    - 2.0*c*(argument[1] - pi) - 2.0*e*(argument[0]-pi) - 4.0*b*(argument[0]-pi)*(argument[1]-pi));
   hessian[1][1] = a*(3.0*b - 4.0*e*(argument[1] - pi) - 4.0*b*pow((argument[1]-pi), 2));

   return(hessian);
}


// Matrix<double> calculateInverseHessian(void) method

/// This method returns the Easom's analytical inverse Hessian matrix. 

Matrix<double> EasomFunction::calculateInverseHessian(void)
{
   const int numberOfVariables = 2;

   Matrix<double> inverseHessian(numberOfVariables, numberOfVariables);
    
   Vector<double> argument = multilayerPerceptron->getIndependentParameters();

   Matrix<double> hessian = calculateHessian();
   
   double hessianDeterminant = hessian[0][0]*hessian[1][1] - hessian[0][1]*hessian[1][0];
 
   inverseHessian[0][0] = hessian[1][1]/hessianDeterminant;
   inverseHessian[0][1] = -hessian[0][1]/hessianDeterminant;
   inverseHessian[1][0] = -hessian[1][0]/hessianDeterminant;
   inverseHessian[1][1] = hessian[0][0]/hessianDeterminant;
 
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
