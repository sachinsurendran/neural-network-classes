/****************************************************************************************************************/
/*                                                                                                              */
/*   Flood: An Open Source Neural Networks C++ Library                                                          */
/*   www.cimne.com/flood                                                                                        */
/*                                                                                                              */
/*   M I N I M U M   D R A G   P R O B L E M   C L A S S                                                        */
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
#include <float.h>

#include "MinimumDragProblem.h"

#include "../Utilities/OrdinaryDifferentialEquations.h"
#include "../Utilities/IntegrationOfFunctions.h"

namespace Flood
{

// GENERAL CONSTRUCTOR

/// General constructor. It creates a minimum drag problem objective functional associated to a multilayer 
/// perceptron.
/// It also initializes all the problem parameters to their default values:
///
/// <ul>
/// <li> Tolerance = 1.0e-12
/// <li> Reseve = 1000000
/// </ul>
/// 
/// @param newMultilayerPerceptron Pointer to a multilayer perceptron object.

MinimumDragProblem::MinimumDragProblem(MultilayerPerceptron* newMultilayerPerceptron)
: ObjectiveFunctional(newMultilayerPerceptron)
{                     
   // Control sentence

   int numberOfInputs = multilayerPerceptron->getNumberOfInputs();
   int numberOfOutputs = multilayerPerceptron->getNumberOfOutputs();

   if(numberOfInputs != 1 || numberOfOutputs != 1)
   {
      std::cerr << std::endl
                << "Flood Error: MinimumDragProblem class." << std::endl
                << "MinimumDragProblem(MultilayerPerceptron*) constructor." << std::endl
                << "Number of inputs and outputs in multilayer perceptron must be 1." << std::endl
                << std::endl;

      exit(1);
   }

   tolerance = 1.0e-6;
   initialSize = 1000;

   ordinaryDifferentialEquations.setTolerance(tolerance);
   ordinaryDifferentialEquations.setInitialSize(initialSize);
}


// DEFAULT CONSTRUCTOR

/// Default constructor. It creates a minimum drag problem objective functional not associated to any multilayer 
/// perceptron. 
/// It also initializes all the problem parameters to their default values:
///
/// <ul>
/// <li> Tolerance = 1.0e-12
/// <li> Reseve = 1000000
/// </ul>

MinimumDragProblem::MinimumDragProblem(void) : ObjectiveFunctional()
{
   tolerance = 1.0e-12;
   initialSize = 1000;

   ordinaryDifferentialEquations.setTolerance(tolerance);
   ordinaryDifferentialEquations.setInitialSize(initialSize);
}


// DESTRUCTOR

/// Destructor.

MinimumDragProblem::~MinimumDragProblem(void) 
{

}


// METHODS

// double getTolerance(void) method

/// This method returns the tolerance value to be used in the Runge-Kutta-Fehlberg method for evaluating the drag 
/// of an axisymmetric body.

double MinimumDragProblem::getTolerance(void)
{
   return(tolerance);   
}


// int getInitialSize(void) method

/// This method returns the number of points to be reserved when using the Runge-Kutta-Fehlberg method for 
/// evaluating the drag of an axisymmetric body.

int MinimumDragProblem::getInitialSize(void)
{
   return(initialSize);   
}


// void setTolerance(double) method

/// This method sets a new tolerance value to be used in the Runge-Kutta-Fehlberg method for evaluating the drag 
/// of an axisymmetric body.
/// The tolerance of integration must be a small value greater than zero. 
///
/// @param newTolerance Tolerance in Runge-Kutta-Fehlberg method.

void MinimumDragProblem::setTolerance(double newTolerance)
{
   // Control sentence

   if(newTolerance <= 0.0)
   {
      std::cerr << std::endl
                << "Flood Error: MinimumDragProblem class." << std::endl
                << "void setTolerance(double) method." << std::endl
                << "The tolerance of integration must be greater than zero." << std::endl
                << std::endl;

      exit(1);
   }

   // Set tolerance 

   tolerance = newTolerance;

   ordinaryDifferentialEquations.setTolerance(tolerance);
}


// void setInitialSize(int) method

/// This method sets a new number of points to be reserved when using the Runge-Kutta-Fehlberg method for 
/// evaluating the drag of an axisymmetric body.
/// The initial size must be a big value greater than zero, otherwise resizing will be necessary. 
///
/// @param newInitialSize Number of points to reserve in Runge-Kutta-Fehlberg method.

void MinimumDragProblem::setInitialSize(int newInitialSize)
{
   // Set initial size

   initialSize = newInitialSize;

   ordinaryDifferentialEquations.setInitialSize(initialSize);
}


// Vector<double> calculateParticularSolution(Vector<double>) method

/// This method returns the particular solution term phi0(x) in order to satisfy the boundary conditions for the 
/// minimum drag problem. 
/// It must hold phi0(xa)=ya if there is a condition y(xa)=ya.
///
/// @param input Value of independent variable x.
///
/// @see calculateHomogeneousSolution.
/// @see calculateParticularSolutionDerivative.
/// @see calculateHomogeneousSolutionDerivative.

Vector<double> MinimumDragProblem::calculateParticularSolution(Vector<double> input)
{
   // Control sentence

   int size = input.getSize();

   if(size != 1)
   {
      std::cerr << std::endl
                << "Flood Error: MinimumDragProblem class." << std::endl
                << "Vector<double> calculateParticularSolution(Vector<double>) method." << std::endl
                << "Size of input must be 1." << std::endl
                << std::endl;

      exit(1);
   }

   // Calculate particular solution

   Vector<double> particularSolution(1);

   double x = input[0];

   double xa = 0.0;
   double ya = 0.0;

   double xb = 1.0;
   double yb = 1.0;

   particularSolution[0] = ya + (yb-ya)*(x-xa)/(xb-xa);

   return(particularSolution);
}


// Vector<double> calculateHomogeneousSolution(Vector<double>)

/// This method returns the homogeneous solution term phi1(x) in order to satisfy the boundary conditions for the 
/// minimum drag problem. 
/// It must hold phi1(xa)=0 if there is a condition y(xa)=ya.
///
/// @param input Value of independent variable x.
///
/// @see calculateParticularSolution.
/// @see calculateParticularSolutionDerivative.
/// @see calculateHomogeneousSolutionDerivative.

Vector<double> MinimumDragProblem::calculateHomogeneousSolution(Vector<double> input)
{
   // Control sentence

   int size = input.getSize();

   if(size != 1)
   {
      std::cerr << std::endl
                << "Flood Error: MinimumDragProblem class." << std::endl
                << "Vector<double> calculateHomogeneousSolution(Vector<double>) method." << std::endl
                << "Size of input must be 1." << std::endl
                << std::endl;

      exit(1);
   }

   // Calculate homogeneous solution

   Vector<double> homogeneousSolution(1);

   double x = input[0];

   double xa = 0.0;
   double xb = 1.0;

   homogeneousSolution[0] = (x-xa)*(x-xb);

   return(homogeneousSolution);
}


// Vector<double> calculateParticularSolutionDerivative(Vector<double>)

/// This method returns the derivative of the particular solution term, phi0'(x), which will be used to compute 
/// the derivative of the descent curve, y'(x), in the minimum drag problem. 
///
/// @param input Value of independent variable x.
///
/// @see getParticularSolutionSolution.
/// @see calculateHomogeneousSolution.
/// @see calculateHomogeneousSolutionDerivative.

Vector<double> MinimumDragProblem::calculateParticularSolutionDerivative(Vector<double> input)
{
   // Control sentence

   int size = input.getSize();

   if(size != 1)
   {
      std::cerr << std::endl
                << "Flood Error: MinimumDragProblem class." << std::endl
                << "Vector<double> calculateParticularSolutionDerivative(Vector<double>) method." << std::endl
                << "Size of input must be 1." << std::endl
                << std::endl;

      exit(1);
   }

   // Calculate particular solution derivative

   Vector<double> particularSolutionDerivative(1);
   
   double xa = 0.0;
   double ya = 0.0;

   double xb = 1.0;
   double yb = 1.0;

   particularSolutionDerivative[0] = (yb-ya)/(xb-xa);

   return(particularSolutionDerivative);
}


// Vector<double> calculateHomogeneousSolutionDerivative(Vector<double>)

/// This method returns the derivative of the homogeneous solution term, phi1'(x), which will be used to compute 
/// the derivative of the descent curve, y'(x), in the minimum drag problem. 
///
/// @param input Value of independent variable x.
///
/// @see getParticularSolutionSolution.
/// @see calculateHomogeneousSolution.
/// @see calculateHomogeneousSolutionDerivative.

Vector<double> MinimumDragProblem::calculateHomogeneousSolutionDerivative(Vector<double> input)
{
   // Control sentence

   int size = input.getSize();

   if(size != 1)
   {
      std::cerr << std::endl
                << "Flood Error: MinimumDragProblem class." << std::endl
                << "Vector<double> calculateHomogeneousSolutionDerivatibe(Vector<double>) method." << std::endl
                << "Size of input must be 1." << std::endl
                << std::endl;

      exit(1);
   }

   // Calculate particular solution

   Vector<double> homogeneousSolutionDerivative(1);

   double x = input[0];

   double xa = 0.0;
   double xb = 1.0;

   homogeneousSolutionDerivative[0] = (x-xa) + (x-xb);

   return(homogeneousSolutionDerivative);
}


// double getDragIntegrand(void) method

/// This method returns the drag integrand for a value of x. 
///
/// @param x Value of independent variable.
/// @param dummy Dummy variable used to make this method integrable as an ordinary differential equation.
///
/// @see getDrag.

double MinimumDragProblem::getDragIntegrand(double x, double dummy)
{
   double dragIntegrand = 0.0;

   Vector<double> input(1, x);

   // Get output

   Vector<double> output = multilayerPerceptron->calculateOutput(input, *this);

   double y = output[0];

   // Get output derivative 

   Matrix<double> jacobian = multilayerPerceptron->calculateJacobian(input, *this);
  
   double dydx = jacobian[0][0];

   if(y < 0.0)
   {
      y = 0.0;
      dydx = 0.0;
   }
   else if(y > 1.0)
   {
      y = 1.0;
      dydx = 0.0;
   }
   
   dragIntegrand = y*pow(dydx,3); 
 
   return(dragIntegrand);       
}


// double getDrag(void) method

/// This method integrates the drag integrand function to obtain the drag of the axisymmetric body represented 
/// by the neural network. 
/// Note that the integral is performed with an ordinary differential equations approach.
///
/// @see getDragIntegrand.

double MinimumDragProblem::getDrag(void)
{
   double drag = 0.0;

   Vector<double> x;
   Vector<double> y;

   int numberOfPoints = ordinaryDifferentialEquations.calculateRungeKuttaFehlbergIntegral(*this,
   x, y, 
   &MinimumDragProblem::getDragIntegrand,
   0.0, 1.0, 0.0);

   drag = y[numberOfPoints-1];
   
   return(drag);
}


// double calculateEvaluation(void) method

/// This method returns the evaluation of a multilayer perceptron for the minimum drag problem. 

double MinimumDragProblem::calculateEvaluation(void)
{
    // Control sentence 

   if(multilayerPerceptron == NULL)
   {
      std::cerr << std::endl
                << "Flood Error: MinimumDragProblem class." << std::endl
                << "double calculateEvaluation(void) method." << std::endl
                << "Pointer to multilayer perceptron object cannot be NULL." << std::endl
                << std::endl;

      exit(1);
   }

   double evaluation = 0.0;

   // Increment number of objective functional evaluations

   numberOfEvaluations++;

   double drag = getDrag();
   
   evaluation = drag; 

   return(evaluation);
}


// void saveResults(char*) method

/// This method saves the values of the independent and dependent variables for the minimum drag problem to a 
/// data file. 
///
/// <ol>
/// <li> Value of independent variariable (x).
/// <li> Value of dependent variariable (y).
/// </ol>
///
/// @param filename Filename.

void MinimumDragProblem::saveResults(char* filename)
{
   std::fstream file; 

   file.open(filename, std::ios::out);

   if(!file.is_open())
   {
      std::cout << std::endl
                << "Cannot open minimum drag problem results data file."
                << std::endl << std::endl;

      exit(1);
   }
   else
   {
      if(display)
      {
         std::cout << std::endl
                   << "Saving minimum drag problem results to data file..."
                   << std::endl;
      }
   }

   file << "% Flood Neural Network. Minimum drag problem results data file." << std::endl
        << "% Column data:" << std::endl
        << "% 1 - Independent variable (x)" << std::endl
        << "% 2 - Dependent variable (x)" << std::endl;
        
   Vector<double> input(1);
   Vector<double> output(1);

   double x = 0.0;
   double y = 0.0;
        
   int numberOfPoints = 1001;

   for(int i = 0; i < numberOfPoints; i++)
   {
      x = (double)i/(double)(numberOfPoints-1.0);
      
      input[0] = x;
      
      output = multilayerPerceptron->calculateOutput(input, *this);

      y = output[0];

      if(y < 0.0)
      {
         y = 0.0;
      } 
      else if(y > 1.0)
      {
         y = 1.0;
      }
     
      file << x << " " << y << std::endl;
   }                

   file.close();
}

}

// Flood: An Open Source Neural Networks C++ Library.
// Copyright (C) 2005-2007 Roberto Lopez 
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
