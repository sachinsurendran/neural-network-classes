/****************************************************************************************************************/
/*                                                                                                              */
/*   Flood: An Open Source Neural Networks C++ Library                                                          */
/*   www.cimne.com/flood                                                                                        */
/*                                                                                                              */
/*   P L A N E - C Y L I N D E R   C L A S S                                                                    */
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
#include <math.h>

#include "PlaneCylinder.h"

namespace Flood
{

// GENERAL CONSTRUCTOR

/// General constructor. It creates a plane-cylinder objective function object.
/// It also initializes all the rest of class members to their default values:
///
/// <ul>
/// <li> Minimum: -5.12,...,-5.12.
/// <li> Maximum: 5.12,...,5.12.
/// <li> Penalty: 100.
/// </ul> 

PlaneCylinder::PlaneCylinder(MultilayerPerceptron* newMultilayerPerceptron)       
: ObjectiveFunctional(newMultilayerPerceptron)
{
   // Control sentence

   int numberOfInputs = multilayerPerceptron->getNumberOfInputs();
   int numberOfOutputs = multilayerPerceptron->getNumberOfOutputs();

   if(numberOfInputs != 0 || numberOfOutputs != 0)
   {
      std::cerr << std::endl
                << "Flood Error: PlaneCylinder class." << std::endl
                << "PlaneCylinder(MultilayerPerceptron*) constructor." << std::endl
                << "Number of inputs and outputs in multilayer perceptron must be 0." << std::endl
                << std::endl;

      exit(1);
   }

   int numberOfIndependentParameters = multilayerPerceptron->getNumberOfIndependentParameters();

   if(numberOfIndependentParameters != 2)
   {
      std::cerr << std::endl
                << "Flood Error: PlaneCylinder class." << std::endl
                << "PlaneCylinder(MultilayerPerceptron*) constructor." << std::endl
                << "Number of independent parameters in multilayer perceptron must be 2." << std::endl
                << std::endl;

      exit(1);
   }

   penalty = 100.0;
}


// DESTRUCTOR

/// Destructor.

PlaneCylinder::~PlaneCylinder(void)
{

}


// METHODS

// double getPenalty(void) method

/// This method returns the penalty term ratio to be used in the plane-cylinder problem.
///
/// @see calculateError(void)
/// @see calculateEvaluation(void)

double PlaneCylinder::getPenalty(void)
{
   return(penalty);
}


// void setPenalty(double) method

/// This method sets a new penalty term ratio to be used in the plane-cylinder problem.
///
/// @param newPenalty New penalty term ratio.
///
/// @see calculateError(void)
/// @see calculateEvaluation(void)

void PlaneCylinder::setPenalty(double newPenalty)
{
   // Control sentence

   if(newPenalty <= 0.0)
   {
      std::cerr << std::endl
                << "Flood Error: PlaneCylinder class." << std::endl
                << "void setPenalty(double) method." << std::endl
                << "Penalty must be greater than 0." << std::endl
                << std::endl;

      exit(1);
   }

   // Set penalty

   penalty = newPenalty;
}


// double calculateError(void) method

/// This method returns the error made in the constraint by a given argument. 
///
/// @see calculateEvaluation(void)

double PlaneCylinder::calculateError(void)
{
   double error = 0.0;

   Vector<double> argument = multilayerPerceptron->getIndependentParameters();

   double x = argument[0];
   double y = argument[1];
   
   if(pow(x,2) + pow(y,2) <= 1.0)
   {
      error = 0.0;
   }
   else
   {
      error = pow(x,2) + pow(y,2) - 1.0;
   } 
   
   return(error);       
}


// double calculateEvaluation(void) method

/// This method returns the plane-cylinder function evaluation for a given argument.

double PlaneCylinder::calculateEvaluation(void)
{
    // Control sentence 

   if(multilayerPerceptron == NULL)
   {
      std::cerr << std::endl
                << "Flood Error: PlaneCylinder class." << std::endl
                << "double calculateEvaluation(void) method." << std::endl
                << "Pointer to multilayer perceptron object cannot be NULL." << std::endl
                << std::endl;

      exit(1);
   }

   // Calculate evaluation

   numberOfEvaluations++;

   double evaluation = 0.0;

   Vector<double> argument = multilayerPerceptron->getIndependentParameters();

   double x = argument[0];
   double y = argument[1];

   double error = calculateError();

   evaluation = x + y + penalty*pow(error, 2);

   return(evaluation);
}


// void print(void) method

/// This method prints to the screen the error made in the constraint by a given argument during the optimization 
/// process.  
///
/// @see calculateError(Vector<double>)

void PlaneCylinder::print(void)
{
   double error = calculateError();
   
   std::cout << "Flood Error " << error << std::endl;
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
