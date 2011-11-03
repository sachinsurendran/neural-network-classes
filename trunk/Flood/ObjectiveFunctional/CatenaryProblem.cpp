/****************************************************************************************************************/
/*                                                                                                              */
/*   Flood: An Open Source Neural Networks C++ Library                                                          */
/*   www.cimne.com/flood                                                                                        */
/*                                                                                                              */
/*   C A T E N A R Y   P R O B L E M   C L A S S                                                                */
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

#include "CatenaryProblem.h"     
#include "../Utilities/OrdinaryDifferentialEquations.h"     

namespace Flood
{

// GENERAL CONSTRUCTOR

/// General constructor. It creates a catenary problem objective functional associated to a multilayer perceptron.
/// It also initializes all the problem parameters to their default values:
///
/// <ul>
/// <li> xA = 0.0.
/// <li> yA = 1.0.
/// <li> xB = 1.0.
/// <li> yB = 1.0.
/// <li> Length goal = 1.5.
/// <li> Potential energy weight = 1.0.
/// <li> Length error weight = 1.00.
/// <li> Tolerance = 1.0e-15.
/// <li> Initial size = 1e3.
/// </ul>
///
/// @param newMultilayerPerceptron Pointer to a multilayer perceptron object.

CatenaryProblem::CatenaryProblem(MultilayerPerceptron* newMultilayerPerceptron)
: ObjectiveFunctional(newMultilayerPerceptron)
{
   // Control sentence

   int numberOfInputs = multilayerPerceptron->getNumberOfInputs();
   int numberOfOutputs = multilayerPerceptron->getNumberOfOutputs();

   if(numberOfInputs != 1 || numberOfOutputs != 1)
   {
      std::cerr << std::endl
                << "Flood Error: CatenaryProblem class." << std::endl
                << "CatenaryProblem(MultilayerPerceptron*) constructor." << std::endl
                << "Number of inputs and outputs in multilayer perceptron must be 1." << std::endl
                << std::endl;

      exit(1);
   }

   xA = 0.0;
   yA = 1.0;
   xB = 1.0;
   yB = 1.0;

   lengthGoal = 1.5;

   potentialEnergyWeight = 1.0;
   lengthErrorWeight = 10.0;

   tolerance = 1.0e-12;
   initialSize = 1000;

   ordinaryDifferentialEquations.setTolerance(tolerance);
   ordinaryDifferentialEquations.setInitialSize(initialSize);
}


// DEFAULT CONSTRUCTOR

/// Default constructor. It creates a catenary problem objective functional not associated to any multilayer 
/// perceptron. It also initializes all the problem parameters to their default values:
///
/// <ul>
/// <li> xA = 0.0
/// <li> yA = 1.0
/// <li> xB = 1.0
/// <li> yB = 1.0
/// <li> Length goal = 1.5
/// <li> Potential energy weight = 1.0
/// <li> Length error weight = 1.00
/// <li> Tolerance = 1.0e-15
/// <li> Initial size = 1e3
/// </ul>

CatenaryProblem::CatenaryProblem(void) : ObjectiveFunctional()
{
   xA = 0.0;
   yA = 1.0;
   xB = 1.0;
   yB = 1.0;

   lengthGoal = 1.5;

   potentialEnergyWeight = 1.0;
   lengthErrorWeight = 100.0;

   tolerance = 1.0e-12;
   initialSize = 1000;
}


// DESTRUCTOR

/// Destructor. 

CatenaryProblem::~CatenaryProblem(void) 
{

}


// METHODS

// double getXA(void) method

/// This method returns the x-value for the A point in the catenary problem.

double CatenaryProblem::getXA(void)
{
   return(xA);
}


// double getYA(void) method

/// This method returns the y-value for the A point in the catenary problem.

double CatenaryProblem::getYA(void)
{
   return(yA);
}


// double getXB(void) method

/// This method returns the x-value for the B point in the catenary problem.

double CatenaryProblem::getXB(void)
{
   return(xB);
}

   
// double getYB(void) method

/// This method returns the y-value for the B point in the catenary problem.

double CatenaryProblem::getYB(void)
{
   return(yB);
}
  

// double getLengthGoal(void) method

/// This method returns the chain length goal for the catenary problem. 

double CatenaryProblem::getLengthGoal(void)
{
   return(lengthGoal);
}


// double getLengthErrorWeight(void) method

/// This method returns the weight of the chain length error term in the 
/// objective functional expression.

double CatenaryProblem::getLengthErrorWeight(void)
{
   return(lengthErrorWeight);
}


// double getPotentialEnergyWeight(void) method

/// This method returns the weight of the potential energy term in the 
/// objective functional expression.

double CatenaryProblem::getPotentialEnergyWeight(void)
{
   return(potentialEnergyWeight);
}


// double getTolerance(void) method

/// This method returns the tolerance value to be used in the Runge-Kutta-Fehlberg method for evaluating the 
/// chain lenght error and the potential energy.

double CatenaryProblem::getTolerance(void)
{
   return(tolerance);   
}


// int getInitialSize(void) method

/// This method returns the number of points to be reserved when using the Runge-Kutta-Fehlberg method for 
/// evaluating the chain lenght error and the potential energy.

int CatenaryProblem::getInitialSize(void)
{
   return(initialSize);   
}


// void setXA(double) method

/// This method sets a new x-value for the A point in the catenary problem.
///
/// @param newXA xA value.

void CatenaryProblem::setXA(double newXA)
{
   // Set xa

   xA = newXA;
}

   
// void setYA(double) method

/// This method sets a new y-value for the A point in the catenary problem.
///
/// @param newYA yA value.

void CatenaryProblem::setYA(double  newYA)
{
   yA = newYA;
}

   
// void setXB(double) method

/// This method sets a new x-value for the B point in the catenary problem.
///
/// @param newXB xB value.

void CatenaryProblem::setXB(double newXB)
{
   // Set xb 

   xB = newXB;
}


// void setYB(double) method

/// This method sets a new y-value for the B point in the catenary problem.
///
/// @param newYB yB value.

void CatenaryProblem::setYB(double  newYB)
{
   yB = newYB;
}


// void setLengthGoal(double) method

/// This method sets a new chain length goal for the catenary problem. 
///
/// @param newLengthGoal Chain length goal.

void CatenaryProblem::setLengthGoal(double newLengthGoal)
{
   lengthGoal = newLengthGoal;
}


// void setLengthErrorWeight(double) method

/// This method sets a new weight value for the chain length error term in the objective functional expression.
///
/// @param newLengthErrorWeight Chain lenght error weight value.

void CatenaryProblem::setLengthErrorWeight(double newLengthErrorWeight)
{
   // Control sentence

   if(newLengthErrorWeight <= 0.0)
   {
      std::cerr << std::endl
                << "Flood Error: CatenaryProblem class." << std::endl
                << "void setLengthErrorWeight(double) method." << std::endl
                << "Length error weight must be greater than 0." << std::endl
                << std::endl;

      exit(1);
   }

   // Set length error weight

   lengthErrorWeight = newLengthErrorWeight;
}


// void setPotentialEnergyWeight(double) method

/// This method sets a new weight value for the chain length error term in the objective functional expression.
///
/// @param newPotentialEnergyWeight Potential energy weight value.

void CatenaryProblem::setPotentialEnergyWeight(double newPotentialEnergyWeight)
{
   // Control sentence

   if(newPotentialEnergyWeight <= 0.0)
   {
      std::cerr << std::endl
                << "Flood Error: CatenaryProblem class." << std::endl
                << "void setPotentialEnergyWeight(double) method." << std::endl
                << "Potential energy weight must be greater than 0." << std::endl
                << std::endl;

      exit(1);
   }

   // Set potential energy weight

   potentialEnergyWeight = newPotentialEnergyWeight;
}


// void setTolerance(double) method

/// This method sets a new tolerance value to be used in the Runge-Kutta-Fehlberg method for evaluating the chain 
/// length and the potential energy.
///
/// @param newTolerance Tolerance in Runge-Kutta-Fehlberg method.

void CatenaryProblem::setTolerance(double newTolerance)
{
   // Control sentence

   if(newTolerance <= 0.0)
   {
      std::cerr << std::endl
                << "Flood Error: CatenaryProblem class." << std::endl
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
/// evaluating the chain length and the potential energy.
///
/// @param newInitialSize Number of points to reserve in Runge-Kutta-Fehlberg method.

void CatenaryProblem::setInitialSize(int newInitialSize)
{
   initialSize =  newInitialSize;

   ordinaryDifferentialEquations.setInitialSize(initialSize);
}


// Vector<double> calculateParticularSolution(Vector<double>)

/// This method returns the particular solution term in order to satisfy the boundary conditions in the catenary 
/// problem.
/// It holds phi0(a)=ya if there is a condition y(a)=ya.
///
/// @param input Value of independent variable x. 
///
/// @see calculateHomogeneousSolution.
/// @see calculateParticularSolutionDerivative.
/// @see calculateHomogeneousSolutionDerivative.

Vector<double> CatenaryProblem::calculateParticularSolution(Vector<double> input)
{
   // Control sentence

   int size = input.getSize();

   if(size != 1)
   {
      std::cerr << std::endl
                << "Flood Error: CatenaryProblem class." << std::endl
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


// double calculateHomogeneousSolution(double)

/// This method returns the homogeneous solution term in order to satisfy the boundary conditions in the catenary 
/// problem.
/// It must hold phi1(a)=0 if there is a condition y(a)=ya.
///
/// @param input Value of independent variable x. 
///
/// @see calculateParticularSolution.
/// @see calculateParticularSolutionDerivative.
/// @see calculateHomogeneousSolutionDerivative.

Vector<double> CatenaryProblem::calculateHomogeneousSolution(Vector<double> input)
{
   // Control sentence

   int size = input.getSize();

   if(size != 1)
   {
      std::cerr << std::endl
                << "Flood Error: CatenaryProblem class." << std::endl
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

Vector<double> CatenaryProblem::calculateParticularSolutionDerivative(Vector<double> input)
{
   // Control sentence

   int size = input.getSize();

   if(size != 1)
   {
      std::cerr << std::endl
                << "Flood Error: CatenaryProblem class." << std::endl
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


// double calculateHomogeneousSolutionDerivative(double)

/// This method returns the derivative of the homogeneous solution term for the boundary conditions.  
/// It is used to calculate the derivative y'(x)=dy/dx.
///
/// @param input Value of independent variable x.  
///
/// @see calculateParticularSolution.
/// @see calculateHomogeneousSolution.
/// @see calculateParticularSolutionDerivative.

Vector<double> CatenaryProblem::calculateHomogeneousSolutionDerivative(Vector<double> input)
{
   // Control sentence

   int size = input.getSize();

   if(size != 1)
   {
      std::cerr << std::endl
                << "Flood Error: CatenaryProblem class." << std::endl
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


// double calculatePotentialEnergyIntegrand(double, double) method

/// This method returns the potential energy integrand for a value of x.
///
/// @param x Value of independent variable.
/// @param dummy Dummy variable to allow integration with Runge-Kutta-Fehlberg method.
///
/// @see calculatePotentialEnergy.

double CatenaryProblem::calculatePotentialEnergyIntegrand(double x, double dummy)
{
   double potentialEnergyIntegrand = 0.0;

   Vector<double> input(1);
   Vector<double> output(1);
   Matrix<double> jacobian(1, 1);

   input[0] = x;
   output = multilayerPerceptron->calculateOutput(input, *this);

    double y = output[0];

   jacobian = multilayerPerceptron->calculateJacobian(input, *this);

    double dydx = jacobian[0][0];

   potentialEnergyIntegrand = y*sqrt(1.0 + pow(dydx,2));

   return(potentialEnergyIntegrand);   
}


// double calculatePotentialEnergy(void) method

/// This method integrates the potential energy integrand function to obtain the potential energy for the chain 
/// represented by the neural network.
/// It uses the Runge-Kutta-Fehlberg method.
///
/// @see calculatePotentialEnergyIntegrand.

double CatenaryProblem::calculatePotentialEnergy(void)
{
   double potentialEnergy = 0.0;

   Vector<double> x;
   Vector<double> y;

   int numberOfPoints = 
   ordinaryDifferentialEquations.calculateRungeKuttaFehlbergIntegral(*this,
   x, y, 
   &CatenaryProblem::calculatePotentialEnergyIntegrand,
   xA, xB, 0.0);

   potentialEnergy = y[numberOfPoints-1];

   return(potentialEnergy);   
}


// double calculateLengthIntegrand(double) method

/// This method returns the chain length integrand for a value of x.
///
/// @param x Value of independent variable.
/// @param dummy Dummy variable to allow integration with Runge-Kutta-Fehlberg method.
///
/// @see calculateLength.

double CatenaryProblem::calculateLengthIntegrand(double x, double dummy)
{
   double lengthIntegrand = 0.0;

   Vector<double> input(1, x);

   Matrix<double> jacobian = multilayerPerceptron->calculateJacobian(input, *this);

    double dydx = jacobian[0][0];

   lengthIntegrand = sqrt(1.0 + pow(dydx,2));

   return(lengthIntegrand);   
}


// double calculateLength(void) method       

/// This method integrates the length integrand function to obtain the chain length between xA and xB of the 
/// function represented by the neural network.
///
/// @see calculateLengthIntegrand.

double CatenaryProblem::calculateLength(void)      
{
   double length = 0.0;

   Vector<double> x;
   Vector<double> y;

   int numberOfPoints = 
   ordinaryDifferentialEquations.calculateRungeKuttaFehlbergIntegral(*this,
   x, y, 
   &CatenaryProblem::calculateLengthIntegrand,
   xA, xB, 0.0);

   length = y[numberOfPoints-1];

   return(length);   
}


// double calculateLengthError(void) method

/// This method returns the error between the length of the function addressed by the neural network and the 
/// length goal.

double CatenaryProblem::calculateLengthError(void)
{
   double lengthError = 0.0;

   double length = calculateLength();

   lengthError = length - lengthGoal;

   return(lengthError);
}


// double calculateEvaluation(void) method

/// This method returns the evaluation value of a multilayer perceptron for the catenary problem. 

double CatenaryProblem::calculateEvaluation()
{
    // Control sentence 

   if(multilayerPerceptron == NULL)
   {
      std::cerr << std::endl
                << "Flood Error: CatenaryProblem class." << std::endl
                << "double calculateEvaluation(void) method." << std::endl
                << "Pointer to multilayer perceptron object cannot be NULL." << std::endl
                << std::endl;

        exit(1);
   }

   // Calculate evaluation

   double evaluation = 0.0;

   // Increment number of evaluations

   numberOfEvaluations++;

   double potentialEnergy = calculatePotentialEnergy();
   double lengthError = calculateLengthError();

   evaluation = potentialEnergyWeight*potentialEnergy + lengthErrorWeight*pow(lengthError,2);

   return(evaluation);
}


// void print(void) method

/// This method prints to the screen useful information about the catenary problem for the multilayer 
/// perceptron:
/// <ol>
/// <li> Length error.
/// <li> Potential energy.
/// </ol>

void CatenaryProblem::print(void)
{
   double potentialEnergy = calculatePotentialEnergy();
   double lengthError = calculateLengthError();

   std::cout << "Length error: " << lengthError << std::endl
             << "Potential energy: " << potentialEnergy << std::endl;
}
       

// void saveResults(char*) method

/// This method saves the values of the independent and dependent variables for the catenary problem to a data 
/// file. 
///
/// <ol>
/// <li> Value of independent variariable (x).
/// <li> Value of ependent variariable (y).
/// </ol>
///
/// @param filename Filename.

void CatenaryProblem::saveResults(char* filename)
{
   std::fstream file; 

   // Control sentence

   file.open(filename, std::ios::out);

   if(!file.is_open())
   {
	   std::cout << std::endl
                 << "Flood Error: CatenaryProblem class." << std::endl
                 << "void saveResults(char*) method." << std::endl
                 << "Cannot open catenary problem results data file."  << std::endl
				 << std::endl;
      
      exit(1);
   }
   else
   {
      if(display)
      {
         std::cout << std::endl
                   << "Saving catenary problem results to data file..." << std::endl;
      }
   }

   Vector<double> input(1);
   Vector<double> output(1);

   double x = 0.0;
   double y = 0.0;
    
   // File header

   file << "% Flood Neural Network. Catenary problem results data file." << std::endl
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
