/****************************************************************************************************************/
/*                                                                                                              */
/*   Flood: An Open Source Neural Networks C++ Library                                                          */
/*   www.cimne.com/flood                                                                                        */
/*                                                                                                              */
/*   G E O D E S I C   P R O B L E M   C L A S S                                                                */
/*                                                                                                              */ 
/*   Roberto Lopez                                                                                              */ 
/*   International Center for Numerical Methods in Engineering (CIMNE)                                          */
/*   Technical University of Catalonia (UPC)                                                                    */
/*   Barcelona, Spain                                                                                           */
/*   E-mail: rlopez@cimne.upc.es                                                                                */ 
/*                                                                                                              */
/****************************************************************************************************************/


#include <iostream>     
#include <fstream>     
#include <math.h>     

#include "GeodesicProblem.h"     
#include "../Utilities/OrdinaryDifferentialEquations.h"     

namespace Flood
{

// GENERAL CONSTRUCTOR

/// General constructor. It creates a geodesic problem objective functional associated to a multilayer perceptron.
/// It also initializes all the problem parameters to their default values:
///
/// <ul>
/// <li> xA = 0.0.
/// <li> yA = 1.0.
/// <li> xB = 1.0.
/// <li> yB = 0.0.
/// <li> Tolerance = 1.0e-12.
/// <li> Initial size = 1000.
/// </ul>
///
/// @param newMultilayerPerceptron Pointer to a multilayer perceptron object.

GeodesicProblem::GeodesicProblem(MultilayerPerceptron* newMultilayerPerceptron)       
: ObjectiveFunctional(newMultilayerPerceptron)
{
   // Control sentence

   int numberOfInputs = multilayerPerceptron->getNumberOfInputs();
   int numberOfOutputs = multilayerPerceptron->getNumberOfOutputs();

   if(numberOfInputs != 1 || numberOfOutputs != 1)
   {
      std::cerr << std::endl
                << "Flood Error: GeodesicProblem class." << std::endl
                << "GeodesicProblem(MultilayerPerceptron*) constructor." << std::endl
                << "Number of inputs and outputs in multilayer perceptron must be 1." << std::endl
                << std::endl;

      exit(1);
   }

   xA = 0.0;
   yA = 1.0;
   xB = 1.0;
   yB = 0.0;

   tolerance = 1.0e-12;
   initialSize = 1000;

   ordinaryDifferentialEquations.setTolerance(tolerance);
   ordinaryDifferentialEquations.setInitialSize(initialSize);
}


// DEFAULT CONSTRUCTOR

/// Default constructor. It creates a geodesic problem objective functional not associated to any multilayer 
/// perceptron. It also initializes all the problem parameters to their default values:
///
/// <ul>
/// <li> xA = 0.0.
/// <li> yA = 1.0.
/// <li> xB = 1.0.
/// <li> yB = 0.0.
/// <li> Tolerance = 1.0e-12.
/// <li> Initial size = 1000.
///   </ul>

GeodesicProblem::GeodesicProblem(void) : ObjectiveFunctional()
{
   xA = 0.0;
   yA = 0.1;
   xB = 1.0;
   yB = 0.0;

   tolerance = 1.0e-12;
   initialSize = 1000;
}


// DESTRUCTOR

/// Destructor. 

GeodesicProblem::~GeodesicProblem(void) 
{

}


// METHODS

// double getXA(void) method

/// This method returns the x-value for the A point in the geodesic problem.

double GeodesicProblem::getXA(void)
{
   return(xA);
}


// double getYA(void) method

/// This method returns the y-value for the A point in the geodesic problem.

double GeodesicProblem::getYA(void)
{
   return(yA);
}


// double getXB(void) method

/// This method returns the x-value for the B point in the geodesic problem.

double GeodesicProblem::getXB(void)
{
   return(xB);
}

   
// double getYB(void) method

/// This method returns the y-value for the B point in the geodesic problem.

double GeodesicProblem::getYB(void)
{
   return(yB);
}
  

// double getTolerance(void) method

/// This method returns the tolerance value to be used in the Runge-Kutta-Fehlberg
/// method for evaluating the arc-length.

double GeodesicProblem::getTolerance(void)
{
   return(tolerance);   
}


// int getInitialSize(void) method

/// This method returns the number of points to be reserved when using the 
/// Runge-Kutta-Fehlberg method for evaluating the arc-length.

int GeodesicProblem::getInitialSize(void)
{
   return(initialSize);   
}


// void setXA(double) method

/// This method sets a new x-value for the A point in the geodesic problem.
///
/// @param newXA xA value.

void GeodesicProblem::setXA(double newXA)
{
   xA = newXA;
}

   
// void setYA(double) method

/// This method sets a new y-value for the A point in the geodesic problem.
///
/// @param newYA yA value.

void GeodesicProblem::setYA(double  newYA)
{
   yA = newYA;
}

   
// void setXB(double) method

/// This method sets a new x-value for the B point in the geodesic problem.
///
/// @param newXB xB value.

void GeodesicProblem::setXB(double  newXB)
{
   xB = newXB;
}


// void setYB(double) method

/// This method sets a new y-value for the B point in the geodesic problem.
///
/// @param newYB yB value.

void GeodesicProblem::setYB(double  newYB)
{
   yB = newYB;
}


// void setTolerance(double) method

/// This method sets a new tolerance value to be used in the Runge-Kutta-Fehlberg method for evaluating the 
/// arc-length.
///
/// @param newTolerance Tolerance in Runge-Kutta-Fehlberg method.

void GeodesicProblem::setTolerance(double newTolerance)
{
   // Control sentence

   if(newTolerance <= 0.0)
   {
      std::cerr << std::endl
                << "Flood Error: GeodesicProblem class." << std::endl
                << "void setTolerance(double) method." << std::endl
                << "Tolerance must be greater than 0." << std::endl
                << std::endl;

      exit(1);
   }

   // Set tolerance

   tolerance = newTolerance;

   ordinaryDifferentialEquations.setTolerance(tolerance);
}


// void setInitialSize(int) method

/// This method sets a new number of points to be reserved when using the Runge-Kutta-Fehlberg method for 
/// evaluating the arc-length.
///
/// @param newInitialSize Number of points to reserve in Runge-Kutta-Fehlberg method.

void GeodesicProblem::setInitialSize(int newInitialSize)
{
   initialSize =  newInitialSize;

   ordinaryDifferentialEquations.setInitialSize(initialSize);
}


// Vector<double> calculateParticularSolution(Vector<double>)

/// This method returns the particular solution term in order to satisfy the boundary conditions in the geodesic 
/// problem.
/// It holds phi0(a)=ya if there is a condition y(a)=ya.
///
/// @param input Value of independent variable x. 
///
/// @see calculateHomogeneousSolution.
/// @see calculateParticularSolutionDerivative.
/// @see calculateHomogeneousSolutionDerivative.

Vector<double> GeodesicProblem::calculateParticularSolution(Vector<double> input)
{
   // Control sentence

   int size = input.getSize();

   if(size != 1)
   {
      std::cerr << std::endl
                << "Flood Error: GeodesicProblem class." << std::endl
                << "Vector<double> calculateParticularSolution(Vector<double>)." << std::endl
                << "Size of input vector must be 1." << std::endl
                << std::endl;

      exit(1);
   }

   // Calculate particular solution

   Vector<double> particularSolution(1);

   double x = input[0];

   particularSolution[0] = yA + (yB-yA)*(x-xA)/(xB-xA);

   return(particularSolution);
}


// Vector<double> calculateHomogeneousSolution(Vector<double>)

/// This method returns the homogeneous solution term in order to satisfy the boundary conditions in the geodesic
/// problem.
/// It must hold phi1(a)=0 if there is a condition y(a)=ya.
///
/// @param input Value of independent variable x. 
///
/// @see calculateParticularSolution.
/// @see calculateParticularSolutionDerivative.
/// @see calculateHomogeneousSolutionDerivative.

Vector<double> GeodesicProblem::calculateHomogeneousSolution(Vector<double> input)
{
   // Control sentence

   int size = input.getSize();

   if(size != 1)
   {
      std::cerr << std::endl
                << "Flood Error: GeodesicProblem class." << std::endl
                << "Vector<double> calculateHomogeneousSolution(Vector<double>)." << std::endl
                << "Size of input vector must be 1." << std::endl
                << std::endl;

      exit(1);
   }

   // Calculate homogeneous solution

   Vector<double> homogeneousSolution(1);

   double x = input[0];

   homogeneousSolution[0] = (x-xA)*(x-xB);

   return(homogeneousSolution);
}


// Vector<double> calculateParticularSolutionDerivative(Vector<double>)

/// This method returns the derivative of the particular solution term for the boundary conditions.  
/// It is used to calculate the derivative y'(x)=dy/dx.
///
/// @param input Value of independent variable x. 
///
/// @see calculateParticularSolution.
/// @see calculateHomogeneousSolution.
/// @see calculateHomogeneousSolutionDerivative.

Vector<double> GeodesicProblem::calculateParticularSolutionDerivative(Vector<double> input)
{
   // Control sentence

   int size = input.getSize();

   if(size != 1)
   {
      std::cerr << std::endl
                << "Flood Error: GeodesicProblem class." << std::endl
                << "Vector<double> calculateParticularSolutionDerivative(Vector<double>)." << std::endl
                << "Size of input vector must be 1." << std::endl
                << std::endl;

      exit(1);
   }

   // Calculate particular solution derivative

   Vector<double> particularSolutionDerivative(1);

   particularSolutionDerivative[0] = (yB-yA)/(xB-xA);

   return(particularSolutionDerivative);
}


// Vector<double> calculateHomogeneousSolutionDerivative(Vector<double>)

/// This method returns the derivative of the homogeneous solution term for the boundary conditions.  
/// It is used to calculate the derivative y'(x)=dy/dx.
///
/// @param input Value of independent variable x.  
///
/// @see calculateParticularSolution.
/// @see calculateHomogeneousSolution.
/// @see calculateParticularSolutionDerivative.

Vector<double> GeodesicProblem::calculateHomogeneousSolutionDerivative(Vector<double> input)
{
   // Control sentence

   int size = input.getSize();

   if(size != 1)
   {
      std::cerr << std::endl
                << "Flood Error: GeodesicProblem class." << std::endl
                << "Vector<double> calculateHomogeneousSolutionDerivative(Vector<double>)." << std::endl
                << "Size of input vector must be 1." << std::endl
                << std::endl;

      exit(1);
   }

   // Calculate homogeneous solution derivative

   Vector<double> homogeneousSolutionDerivative(1);

   double x = input[0];

   homogeneousSolutionDerivative[0] = (x-xA) + (x-xB);

   return(homogeneousSolutionDerivative);
}


// double calculateArcLengthIntegrand(double, double) method

/// This method returns the arc length integrand for a value of x.
///
/// @param x Value of independent variable.
/// @param dummy Dummy variable to allow integration with Runge-Kutta-Fehlberg method.
///
/// @see calculateArcLength.

double GeodesicProblem::calculateArcLengthIntegrand(double x, double dummy)
{
   double arcLengthIntegrand = 0.0;

   Vector<double> input(1, x);

   Matrix<double> jacobian = multilayerPerceptron->calculateJacobian(input, *this);

   double dydx = jacobian[0][0];

   arcLengthIntegrand = sqrt(1.0 + pow(dydx,2));

   return(arcLengthIntegrand);   
}


// double calculateArcLength(void) method       

/// This method integrates the arc length integrand function to obtain the arc length between xa and xb of the 
/// function represented by the neural network. 
/// It uses the Runge-Kutta-Fehlberg method.
///
/// @see calculateArcLengthIntegrand.

double GeodesicProblem::calculateArcLength(void)      
{
   double arcLength = 0.0;

   Vector<double> x;
   Vector<double> y;

   int numberOfPoints 
   = ordinaryDifferentialEquations.calculateRungeKuttaFehlbergIntegral(*this,
   x, y, 
   &GeodesicProblem::calculateArcLengthIntegrand,
   xA, xB, 0.0);

   arcLength = y[numberOfPoints-1];

   return(arcLength);   
}


// double calculateEvaluation(void) method

/// This method returns the evaluation value of a multilayer perceptron for the geodesic problem. 
///
/// @see getArcLenght.

double GeodesicProblem::calculateEvaluation()
{
    // Control sentence 

   if(multilayerPerceptron == NULL)
   {
      std::cerr << std::endl
                << "Flood Error: GeodesicProblem class." << std::endl
                << "double calculateEvaluation(void) method." << std::endl
                << "Pointer to multilayer perceptron object cannot be NULL." << std::endl
                << std::endl;

        exit(1);
   }

   double evaluation = 0.0;

   // Increment number of evaluations

   numberOfEvaluations++;

   evaluation = calculateArcLength();

   return(evaluation);
}


// void saveResults(char*) method

/// This method saves the values of the independent and dependent 
/// variables for the geodesic problem to a data file. 
///
/// <ol>
/// <li> Value of independent variariable (x).
/// <li> Value of ependent variariable (y).
/// </ol>
///
/// @param filename Filename.

void GeodesicProblem::saveResults(char* filename)
{
   std::fstream file; 

   // Control sentence

   file.open(filename, std::ios::out);

   if(!file.is_open())
   {
      std::cout << "Cannot open geodesic problem results data file."  << std::endl;
      
      exit(1);
   }
   else
   {
      if(display)
      {
         std::cout << std::endl
                   << "Saving geodesic problem results to data file..." << std::endl;
      }
   }

   Vector<double> input(1);
   Vector<double> output(1);

   double x = 0.0;
   double y = 0.0;
    
   // File header

    file << "% Flood Neural Network. Geodesic problem results data file." << std::endl
         << "% Column data:" << std::endl
         << "%  1 - Independent variable (x)" << std::endl
         << "%  2 - Dependent variable (y)" << std::endl
         << std::endl;

   // File data

   int numberOfPoints = 1001;

   for(int i = 0; i < numberOfPoints; i++)
   {
      // Obtain x

      x = i/(numberOfPoints-1.0);

      // Obtain input

      input[0] = x;

      // Obtain output

      output = multilayerPerceptron->calculateOutput(input, *this);

      // Obtain y

      y = output[0];

      // Write x and y to file

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
