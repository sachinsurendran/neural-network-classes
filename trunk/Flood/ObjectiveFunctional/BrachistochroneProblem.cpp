/****************************************************************************************************************/
/*                                                                                                              */
/*   Flood: An Open Source Neural Networks C++ Library                                                          */
/*   www.cimne.com/flood                                                                                        */
/*                                                                                                              */
/*   B R A C H I S T O C H R O N E   P R O B L E M    C L A S S                                                 */
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

#include "BrachistochroneProblem.h"
#include "../Utilities/OrdinaryDifferentialEquations.h"
 
namespace Flood
{

// GENERAL CONSTRUCTOR

/// General constructor. It creates a brachistochrone problem objective functional associated to a multilayer 
/// perceptron. 
/// It also initializes all the problem parameters to their default values:
///
/// <ul>
/// <li> xa: 0.0.
/// <li> ya: 1.0.
/// <li> xb: 1.0.
/// <li> yb: 0.0.
/// <li> Tolerance: 1.0e-6.
/// <li> Initial size: 1000000.
/// </ul>
///
/// @param newMultilayerPerceptron Pointer to a multilayer perceptron object.

BrachistochroneProblem::BrachistochroneProblem(MultilayerPerceptron* newMultilayerPerceptron)
: ObjectiveFunctional(newMultilayerPerceptron)
{
   // Control sentence

   int numberOfInputs = multilayerPerceptron->getNumberOfInputs();
   int numberOfOutputs = multilayerPerceptron->getNumberOfOutputs();

   if(numberOfInputs != 1 || numberOfOutputs != 1)
   {
      std::cerr << std::endl
                << "Flood Error: BrachistochroneProblem class." << std::endl
                << "BrachistochroneProblem(MultilayerPerceptron*) constructor." << std::endl
                << "Number of inputs and outputs in multilayer perceptron must be 1." << std::endl
                << std::endl;

      exit(1);
   }

   xA = 0.0;
   yA = 1.0;
   xB = 1.0;
   yB = 0.0;

   // This problem is very sensitive to this parameter!
   tolerance = 1.0e-6;
   initialSize = 1000;

   ordinaryDifferentialEquations.setTolerance(tolerance);
   ordinaryDifferentialEquations.setInitialSize(initialSize);
}


// DEFAULT CONSTRUCTOR

/// Default constructor. It creates a brachistochrone problem objective functional not associated to any 
/// multilayer perceptron. 
/// It also initializes all the problem parameters to their default values:
///
/// <ul>
/// <li> xa: 0.0.
/// <li> ya: 1.0.
/// <li> xb: 1.0.
/// <li> yb: 0.0.
/// <li> Tolerance: 1.0e-6.
/// <li> Initial size: 1000000.
/// </ul>

BrachistochroneProblem::BrachistochroneProblem(void) : ObjectiveFunctional()
{
   xA = 0.0;
   yA = 1.0;
   xB = 1.0;
   yB = 0.0;

   // This problem is very sensitive to this parameter!!!
   tolerance = 1.0e-6;
   initialSize = 1000000;
}


// DESTRUCTOR

/// Destructor. 

BrachistochroneProblem::~BrachistochroneProblem(void) 
{

}


// METHODS

// double getXa(void) method

/// This method returns the x-value of point A in the brachistochrone problem.
///
/// @see getYa.
/// @see getXb.
/// @see getYb.

double BrachistochroneProblem::getXa(void)
{
   return(xA);
}


// double getYa(void) method

/// This method returns the y-value of point A in the brachistochrone problem.
///
/// @see getXa.
/// @see getXb.
/// @see getYb.

double BrachistochroneProblem::getYa(void)
{
   return(yA);
}


// double getXb(void) method

/// This method returns the x-value of point B in the brachistochrone problem.
///
/// @see getXa.
/// @see getYa.
/// @see getYb.

double BrachistochroneProblem::getXb(void)
{
   return(xB);
}


// double getYb(void) method

/// This method returns the y-value of point B in the brachistochrone problem.
///
/// @see getXa.
/// @see getYa.
/// @see getXb.

double BrachistochroneProblem::getYb(void)
{
   return(yB);
}


// double getTolerance(void) method

/// This method returns the tolerance value to be used in the Runge-Kutta-Fehlberg method for evaluating the 
/// descent time.

double BrachistochroneProblem::getTolerance(void)
{
   return(tolerance);   
}


// int getInitialSize(void) method

/// This method returns the number of points to be reserved when using the Runge-Kutta-Fehlberg method for 
/// evaluating the descent time.

int BrachistochroneProblem::getInitialSize(void)
{
   return(initialSize);   
}


// void setXa(double) method

/// This method sets a new x-value for point A in the brachistochrone problem.

void BrachistochroneProblem::setXa(double newXa)
{
   xA = newXa;
}


// void setYa(double) method

/// This method sets a new y-value for point A in the brachistochrone problem.

void BrachistochroneProblem::setYa(double  newYa)
{
   yA = newYa;
}


// void setXb(double) method

/// This method sets a new x-value for point B in the brachistochrone problem.

void BrachistochroneProblem::setXb(double  newXb)
{
   xB = newXb;
}


// void setYb(double) method

/// This method sets a new y-value for point B in the brachistochrone problem.

void BrachistochroneProblem::setYb(double  newYb)
{
   yB = newYb;
}


// void setTolerance(double) method

/// This method sets a new tolerance value to be used in the Runge-Kutta-Fehlberg method for evaluating the 
/// descent time.
///
/// @param newTolerance Tolerance in Runge-Kutta-Fehlberg method.

void BrachistochroneProblem::setTolerance(double newTolerance)
{
   // Control sentence

   if(newTolerance <= 0.0)
   {
      std::cerr << std::endl
                << "Flood Error: BrachistochroneProblem class." << std::endl
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
/// evaluating the descent time.
///
/// @param newInitialSize Number of points to reserve in Runge-Kutta-Fehlberg method.

void BrachistochroneProblem::setInitialSize(int newInitialSize)
{
   initialSize =  newInitialSize;

   ordinaryDifferentialEquations.setInitialSize(initialSize);
}


// void setProblem(double, double, double, double) method

/// This method sets a new pair of points A=(xa,ya) and B=(xb,yb) for the brachistochrone problem.

void BrachistochroneProblem::setProblem(double newXa, double newYa, double newXb, double newYb)
{
   xA = newXa;
   yA = newYa;
   xB = newXb;
   yB = newYb;
}


// Vector<double> calculateParticularSolution(Vector<double>)

/// This method returns the particular solution term phi0(x) in order to satisfy the boundary conditions for the 
/// brachistochrone problem. 
/// It must hold phi0(xa)=ya if there is a condition y(xa)=ya.
///
/// @param input Value of independent variable x.
///
/// @see calculateHomogeneousSolution.
/// @see calculateParticularSolutionDerivative.
/// @see calculateHomogeneousSolutionDerivative.

Vector<double> BrachistochroneProblem::calculateParticularSolution(Vector<double> input)
{
   // Control sentence

   int size = input.getSize();

   if(size != 1)
   {
      std::cerr << std::endl
                << "Flood Error: BrachistochroneProblem class." << std::endl
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

/// This method returns the homogeneous solution term phi1(x) in order to satisfy the boundary conditions for the 
/// brachistochrone problem. 
/// It must hold phi1(xa)=0 if there is a condition y(xa)=ya.
///
/// @param input Value of independent variable x.
///
/// @see calculateParticularSolution.
/// @see calculateParticularSolutionDerivative.
/// @see calculateHomogeneousSolutionDerivative.

Vector<double> BrachistochroneProblem::calculateHomogeneousSolution(Vector<double> input)
{
   // Control sentence

   int size = input.getSize();

   if(size != 1)
   {
      std::cerr << std::endl
                << "Flood Error: BrachistochroneProblem class." << std::endl
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

/// This method returns the derivative of the particular solution term, phi0'(x), which will be used to compute 
/// the derivative of the descent curve, y'(x), in the brachistochrone problem. 
///
/// @param input Value of independent variable x.
///
/// @see getParticularSolutionSolution.
/// @see calculateHomogeneousSolution.
/// @see calculateHomogeneousSolutionDerivative.

Vector<double> BrachistochroneProblem::calculateParticularSolutionDerivative(Vector<double> input)
{
   // Control sentence

   int size = input.getSize();

   if(size != 1)
   {
      std::cerr << std::endl
                << "Flood Error: BrachistochroneProblem class." << std::endl
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

/// This method returns the derivative of the homogeneous solution term, phi1'(x), which will be used to compute 
/// the derivative of the descent curve, y'(x), in the brachistochrone problem. 
///
/// @param input Value of independent variable x.
///
/// @see calculateParticularSolution.
/// @see calculateHomogeneousSolution.
/// @see calculateHomogeneousSolutionDerivative.

Vector<double> BrachistochroneProblem::calculateHomogeneousSolutionDerivative(Vector<double> input)
{
   // Control sentence

   int size = input.getSize();

   if(size != 1)
   {
      std::cerr << std::endl
                << "Flood Error: BrachistochroneProblem class." << std::endl
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


// double calculateDescentTimeIntegrand(double, double) method

/// This method returns the descent time integrand for a value of x. 
///
/// @param x Value of independent variable.
/// @param dummy Dummy variable used to make this method integrable as an ordinary differential equation.
///
/// @see calculateDescentTime.

double BrachistochroneProblem::calculateDescentTimeIntegrand(double x, double dummy)
{
   double descentTimeIntegrand = 0.0;

   Vector<double> input(1, x);


   Vector<double> output = multilayerPerceptron->calculateOutput(input, *this);

   double y = output[0];

   Matrix<double> jacobian = multilayerPerceptron->calculateJacobian(input, *this);

   double dydt = jacobian[0][0];

   // This problem is very sensitive to this parameter!

   double small = 1.0e-6;

   if(y >= yA)
   {
      y = yA - small;
      dydt = 0.0;
   }

   // Get descent time integrand

   double numerator = 1.0 + pow(dydt,2);
   double denominator = yA-y;

   descentTimeIntegrand = sqrt(numerator/denominator);
   
   return(descentTimeIntegrand);   
}


// double calculateDescentTime(void) method

/// This method integrates the descent time integrand function to obtain the time to travel from point A to point 
/// B for the chute represented by the neural network. 
/// Note that the integral is performed with an ordinary differential equations approach.
///
/// @see calculateDescentTimeIntegrand.

double BrachistochroneProblem::calculateDescentTime(void)      
{
   double descentTime = 0.0;

   double const g = 9.81;

   Vector<double> x;
   Vector<double> y;     

   // This problem is very sensitive to this parameter!
   double small = 1.0e-6;

   int numberOfPoints = 
   ordinaryDifferentialEquations.calculateRungeKuttaFehlbergIntegral(*this, x, y, 
   &BrachistochroneProblem::calculateDescentTimeIntegrand,
   xA+small, xB, 0.0);

   double integral = y[numberOfPoints-1];

   descentTime = integral/sqrt(2.0*g);

   return(descentTime);   
}


// double calculateEvaluation(void) method

/// This method returns the objective functional evaluation of a multilayer perceptron for the brachistochrone 
/// problem. 
/// This is simply the descent time of the curve represented by the neural network.

double BrachistochroneProblem::calculateEvaluation()
{
    // Control sentence 

   if(multilayerPerceptron == NULL)
   {
      std::cerr << std::endl
                << "Flood Error: BrachistochroneProblem class." << std::endl
                << "double calculateEvaluation(void) method." << std::endl
                << "Pointer to multilayer perceptron object cannot be NULL." << std::endl
                << std::endl;
      
      exit(1);
   }


   // Increment number of objective evaluations

   numberOfEvaluations++;

   // Evaluate objective functional 

   double descentTime = calculateDescentTime();   

   double evaluation = descentTime;

   return(evaluation);
}


// void save(char*) method

/// This method saves the values of the independent and dependent variables for the brachistochrone problem to a 
/// data file. 
///
/// <ol>
/// <li> Value of independent variariable (x).
/// <li> Value of dependent variariable (y).
/// </ol>
///
/// @param filename Filename.

void BrachistochroneProblem::saveResults(char* filename)
{
   std::fstream file; 

   file.open(filename, std::ios::out);

   if(!file.is_open())
   {
      std::cout << std::endl
                << "Cannot open brachistochrone problem results data file."  
                << std::endl << std::endl;

      exit(1);
   }
   else
   {
      if(display)
     {
         std::cout << std::endl
                   << "Saving brachistochrone problem results to data file..." << std::endl;
     }
   }

   // File header

   file << "% Flood Neural Network. Brachistochrone problem results data file." << std::endl
        << "% Column data:" << std::endl
        << "%  1 - Independent variable (x)" << std::endl
        << "%  2 - Dependent variable (y)" << std::endl
        << std::endl;

   // File data

   Vector<double> input(1);
   Vector<double> output(1);

   double x = 0.0;
   double y = 0.0;
   
   int numberOfPoints = 101;

   for(int i = 0; i < numberOfPoints; i++)
   {
      // Obtain x

      x = xA + (xB - xA)*i/(numberOfPoints-1.0);

      // Obtain y[i]

      input[0] = x;

      output = multilayerPerceptron->calculateOutput(input, *this);

      y = output[0];
      
      // Write x and y to file

      file << x << " " << y << std::endl;   
   }

   file.close();
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
