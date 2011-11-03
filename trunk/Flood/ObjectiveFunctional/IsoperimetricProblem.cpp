/****************************************************************************************************************/
/*                                                                                                              */
/*   Flood: An Open Source Neural Networks C++ Library                                                          */
/*   www.cimne.com/flood                                                                                        */
/*                                                                                                              */
/*   I S O P E R I M E T R I C   P R O B L E M   C L A S S                                                      */
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

#include "IsoperimetricProblem.h"     
#include "../Utilities/OrdinaryDifferentialEquations.h"     

namespace Flood
{

// GENERAL CONSTRUCTOR

/// General constructor. It creates a isoperimetric problem with parametric equations objective functional 
/// associated to a multilayer perceptron. 
/// It also initializes all the problem parameters to their default values:
///
/// <ul>
/// <li> Perimeter goal: 1.0.
/// <li> Perimeter error weight: 1.0.
/// <li> Area weight: 1.0e-3.
/// <li> Tolerance: 1.0e-15.
/// <li> Initial size: 5e3.
/// </ul>
///
/// @param newMultilayerPerceptron Pointer to a multilayer perceptron object.

IsoperimetricProblem::IsoperimetricProblem(MultilayerPerceptron* newMultilayerPerceptron)       
: ObjectiveFunctional(newMultilayerPerceptron)
{
   // Control sentence

   int numberOfInputs = multilayerPerceptron->getNumberOfInputs();
   int numberOfOutputs = multilayerPerceptron->getNumberOfOutputs();

   if(numberOfInputs != 1 || numberOfOutputs != 2)
   {
      std::cerr << std::endl
                << "Flood Error: IsoperimetricProblem class." << std::endl
                << "IsoperimetricProblem(MultilayerPerceptron*) constructor." << std::endl
                << "Number of inputs and outputs in multilayer perceptron must be 1 and 2, respectively."
                << std::endl
                << std::endl;

      exit(1);
   }

   perimeterGoal = 1.0;

   perimeterErrorWeight = 1.0;
   areaWeight = 1.0e-3;

   tolerance = 1.0e-12;
   initialSize = 1000;

   ordinaryDifferentialEquations.setTolerance(tolerance);
   ordinaryDifferentialEquations.setInitialSize(initialSize);
   
}


// DEFAULT CONSTRUCTOR

/// Default constructor. It creates a isoperimetric problem with parametric equations objective functional not 
/// associated to any multilayer perceptron. 
/// It also initializes all the problem parameters to their default values:
///
/// <ul>
/// <li> Perimeter goal: 1.0.
/// <li> Perimeter error weight: 1.0.
/// <li> Area weight: 1.0e-3.
/// <li> Tolerance: 1.0e-15.
/// <li> Initial size: 5e3.
/// </ul>

IsoperimetricProblem::IsoperimetricProblem(void) : ObjectiveFunctional()
{
   perimeterGoal = 1.0;

   perimeterErrorWeight = 1.0;
   areaWeight = 1.0e-3;

   tolerance = 1.0e-12;
   initialSize = 10000;
}


// DESTRUCTOR

/// Destructor. 

IsoperimetricProblem::~IsoperimetricProblem(void) 
{

}


// METHODS

// double getPerimeterGoal(void) method

/// This method returns the closed curve perimeter goal for the isoperimetric problem. 

double IsoperimetricProblem::getPerimeterGoal(void)
{
   return(perimeterGoal);
}


// double getAreaWeight(void) method

/// This method returns the weight of the area term in the objective functional expression.

double IsoperimetricProblem::getAreaWeight(void)
{
   return(areaWeight);
}


// double getPerimeterErrorWeight(void) method

/// This method returns the weight of the perimeter error term in the objective functional expression.

double IsoperimetricProblem::getPerimeterErrorWeight(void)
{
   return(perimeterErrorWeight);
}


// double getTolerance(void) method

/// This method returns the tolerance value to be used in the Runge-Kutta-Fehlberg method for evaluating the 
/// perimeter error and the area of the closed curve.

double IsoperimetricProblem::getTolerance(void)
{
   return(tolerance);   
}


// int getInitialSize(void) method

/// This method returns the number of points to be reserved when using the Runge-Kutta-Fehlberg method for 
/// evaluating the perimeter error and the area of the closed curve.

int IsoperimetricProblem::getInitialSize(void)
{
   return(initialSize);   
}


// void setPerimeterGoal(double) method

/// This method sets a new value for the closed curve perimeter goal in the isoperimetric problem. 
///
/// @param newPerimeterGoal Closed curve perimeter goal value.

void IsoperimetricProblem::setPerimeterGoal(double newPerimeterGoal)
{
   // Control sentence    

   if(newPerimeterGoal <= 0.0)
   {
      std::cerr << std::endl
                << "Flood Error: IsoperimetricProblem class." << std::endl
                << "void setPerimeterGoal(double) method." << std::endl
                << "Perimeter goal must be greater than 0." << std::endl
                << std::endl;

      exit(1);
   }

   // Set perimeter goal

   perimeterGoal = newPerimeterGoal;
}


// void setAreaWeight(double) method

/// This method sets a new value for the weight of the area term in the objective functional expression.
///
/// @param newAreaWeight Area term weight value. 

void IsoperimetricProblem::setAreaWeight(double newAreaWeight)
{
   // Control sentence    

   if(newAreaWeight <= 0.0)
   {
      std::cerr << std::endl
                << "Flood Error: IsoperimetricProblem class." << std::endl
                << "void setAreaWeight(double) method." << std::endl
                << "Area weight must be greater than 0." << std::endl
                << std::endl;

      exit(1);
   }

   // Set area weight

   areaWeight = newAreaWeight;
}


// void setPerimeterErrorWeight(double) method

/// This method sets a new value for the weight of the perimeter error term in the objective functional 
/// expression.
///
/// @param newPerimeterErrorWeight Perimeter error term weight value. 

void IsoperimetricProblem::setPerimeterErrorWeight(double newPerimeterErrorWeight)
{
   // Control sentence    

   if(newPerimeterErrorWeight <= 0.0)
   {
      std::cerr << std::endl
                << "Flood Error: IsoperimetricProblem class." << std::endl
                << "void setPerimeterErrorWeight(double) method." << std::endl
                << "Perimeter error weight must be greater than 0." << std::endl
                << std::endl;

      exit(1);
   }

   // Set perimeter error weight

   perimeterErrorWeight = newPerimeterErrorWeight;
}


// void setTolerance(double) method

/// This method sets a new tolerance value to be used in the Runge-Kutta-Fehlberg method for evaluating the 
/// perimeter and the area of the closed curve.
///
/// @param newTolerance Tolerance in Runge-Kutta-Fehlberg method.

void IsoperimetricProblem::setTolerance(double newTolerance)
{
   // Control sentence    

   if(newTolerance <= 0.0)
   {
      std::cerr << std::endl
                << "Flood Error: IsoperimetricProblem class." << std::endl
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
/// evaluating the perimeter and the area of the closed curve.
///
/// @param newInitialSize Number of points to reserve in Runge-Kutta-Fehlberg method.

void IsoperimetricProblem::setInitialSize(int newInitialSize)
{
   initialSize =  newInitialSize;

   ordinaryDifferentialEquations.setInitialSize(initialSize);
}


// Vector<double> calculateParticularSolution(Vector<double>) method

/// This method returns the particular solution term in order to satisfy the boundary conditions in the 
/// isoperimetric problem.
/// It holds phi0(a)=ya if there is a condition y(a)=ya.
///
/// @param input Value of independent variable t. 
///
/// @see calculateHomogeneousSolution.
/// @see calculateParticularSolutionDerivative.
/// @see calculateHomogeneousSolutionDerivative.

Vector<double> IsoperimetricProblem::calculateParticularSolution(Vector<double> input)
{
   // Control sentence    
   
   int size = input.getSize();

   if(size != 1)
   {
      std::cerr << std::endl
                << "Flood Error: IsoperimetricProblem class." << std::endl
                << "Vector<double> calculateParticularSolution(Vector<double>) method." << std::endl
                << "Size of input vector must be 1." << std::endl
                << std::endl;

      exit(1);
   }

   // Calculate particular solution

   Vector<double> particularSolution(2);
   particularSolution[0] = 0.0;
   particularSolution[1] = 0.0;

   return(particularSolution);
}


// Vector<double> calculateHomogeneousSolution(Vector<double>) method

/// This method returns the homogeneous solution term in order to satisfy the boundary conditions in the 
/// isoperimetric problem.
/// It must hold phi1(a)=0 if there is a condition y(a)=ya.
///
/// @param input Value of independent variable t.  
///
/// @see calculateParticularSolution.
/// @see calculateParticularSolutionDerivative.
/// @see calculateHomogeneousSolutionDerivative.

Vector<double> IsoperimetricProblem::calculateHomogeneousSolution(Vector<double> input)
{
   // Control sentence    
   
   int size = input.getSize();

   if(size != 1)
   {
      std::cerr << std::endl
                << "Flood Error: IsoperimetricProblem class." << std::endl
                << "Vector<double> calculateHomogeneousSolution(Vector<double>) method." << std::endl
                << "Size of input vector must be 1." << std::endl
                << std::endl;

      exit(1);
   }

   // Calculate homogeneous solution

   Vector<double> homogeneousSolution(2);

   double t = input[0];

   homogeneousSolution[0] = t*(t-1.0);
   homogeneousSolution[1] = t*(t-1.0);

   return(homogeneousSolution);
}


// Vector<double> calculateParticularSolutionDerivative(double) method

/// This method returns the derivative of the particular solution term for the boundary conditions.  
/// It is used to calculate the derivatives x'(t)=dx/dt and y'(t)=dy/dt.
///
/// @param input Value of independent variable t. 
///
/// @see calculateParticularSolution.
/// @see calculateHomogeneousSolution.
/// @see calculateHomogeneousSolutionDerivative.

Vector<double> IsoperimetricProblem::calculateParticularSolutionDerivative(Vector<double> input)
{
   // Control sentence    
   
   int size = input.getSize();

   if(size != 1)
   {
      std::cerr << std::endl
                << "Flood Error: IsoperimetricProblem class." << std::endl
                << "Vector<double> calculateParticularSolutionDerivative(Vector<double>) method." << std::endl
                << "Size of input vector must be 1." << std::endl
                << std::endl;

      exit(1);
   }

   // Calculate particular solution derivative

   Vector<double> particularSolutionDerivative(2);
   particularSolutionDerivative[0] = 0.0;
   particularSolutionDerivative[1] = 0.0;

   return(particularSolutionDerivative);
}


// Vector<double> calculateHomogeneousSolutionDerivative(double) method

/// This method returns the derivative of the homogeneous solution term for the boundary conditions.  
/// It is used to calculate the derivatives x'(t)=dx/dt and y'(t)=dy/dt.
///
/// @param input Value of independent variable t. 
///
/// @see calculateParticularSolution.
/// @see calculateHomogeneousSolution.
/// @see calculateParticularSolutionDerivative.

Vector<double> IsoperimetricProblem::calculateHomogeneousSolutionDerivative(Vector<double> input)
{
   // Control sentence    
   
   int size = input.getSize();

   if(size != 1)
   {
      std::cerr << std::endl
                << "Flood Error: IsoperimetricProblem class." << std::endl
                << "Vector<double> calculateHomogeneousSolutionDerivative(Vector<double>) method." << std::endl
                << "Size of input vector must be 1." << std::endl
                << std::endl;

      exit(1);
   }

   // Calculate homogeneous solution derivative

   Vector<double> homogeneousSolutionDerivative(2);

   double t = input[0];

   homogeneousSolutionDerivative[0] = 2.0*t - 1.0;
   homogeneousSolutionDerivative[1] = 2.0*t - 1.0;

   return(homogeneousSolutionDerivative);
}


// double calculateAreaIntegrand(double) method

/// This method returns the area integrand for a value of of the independent variable t.
///
/// @param t Value of independent variable.
/// @param dummy Dummy variable to allow integration with the Runge-Kutta-Fehlberg method.
///
/// @see calculateArea.

double IsoperimetricProblem::calculateAreaIntegrand(double t, double dummy)
{
   // Control sentence    
   
   if(t < 0.0 || t > 1.0)
   {
      std::cerr << std::endl
                << "Flood Error: IsoperimetricProblem class." << std::endl
                << "double calculateAreaIntegrand(double) method." << std::endl
                << "Parameter t must lie in the range [0,1]." << std::endl
                << std::endl;

      exit(1);
   }

   // Calculate area integrand

   double areaIntegrand = 0.0;
   
   // Set input 

   Vector<double> input(1, t);

   // Get output

   Vector<double> output = multilayerPerceptron->calculateOutput(input, *this);

   double x = output[0];
   double y = output[1];

   // Get jacobian
   
   Matrix<double> jacobian = multilayerPerceptron->calculateJacobian(input, *this);

   double dxdt = jacobian[0][0];
   double dydt = jacobian[0][1];

   // Calculate area integrand

   areaIntegrand = x*dydt - y*dxdt;
   
   return(areaIntegrand);   
}


// double calculateArea(void) method       

/// This method integrates the area integrand function to obtain the area enclosed by the function represented 
/// by the neural network.
///
/// @see calculateAreaIntegrand.

double IsoperimetricProblem::calculateArea(void)      
{
   double area = 0.0;

   Vector<double> x;
   Vector<double> y;

   int numberOfPoints = 
   ordinaryDifferentialEquations.calculateRungeKuttaFehlbergIntegral(*this, x, y, 
   &IsoperimetricProblem::calculateAreaIntegrand,
   0.0, 1.0, 0.0);

   area = fabs(0.5*y[numberOfPoints-1]);

   return(area);   
}


// double calculatePerimeterIntegrand(double, double) method

/// This method returns the perimeter integrand for a value of the independent variable t.
///
/// @param t Value of independent variable.
/// @param dummy Dummy variable to allow integration with the Runge-Kutta-Fehlberg method.
///
/// @see calculatePerimeter.

double IsoperimetricProblem::calculatePerimeterIntegrand(double t, double dummy)
{
   // Control sentence    
   
   if(t < 0.0 || t > 1.0)
   {
      std::cerr << std::endl
                << "Flood Error: IsoperimetricProblem class." << std::endl
                << "double calculatePerimeterIntegrand(double) method." << std::endl
                << "Parameter t must lie in the range [0,1]." << std::endl
                << std::endl;

      exit(1);
   }

   // Calculate perimeter integrand

   double perimeterIntegrand = 0.0;

   // Set input

   Vector<double> input(1, t);

   // Get jacobian
   
   Matrix<double> jacobian = multilayerPerceptron->calculateJacobian(input, *this);

   double dxdt = jacobian[0][0];
   double dydt = jacobian[0][1];

   // Calculate perimeter integrand

   perimeterIntegrand = sqrt(pow(dxdt,2) + pow(dydt,2));

   return(perimeterIntegrand);
}


// double calculatePerimeter(void) method       
//
/// This method integrates the perimeter integrand function to obtain 
/// the arc length of the function represented by the neural network.
///
/// @see calculatePerimeterIntegrand.

double IsoperimetricProblem::calculatePerimeter(void)      
{
   double perimeter = 0.0;

   Vector<double> x;
   Vector<double> y;

   int numberOfPoints = 
   ordinaryDifferentialEquations.calculateRungeKuttaFehlbergIntegral(*this, x, y, 
   &IsoperimetricProblem::calculatePerimeterIntegrand,
   0.0, 1.0, 0.0);

   perimeter = y[numberOfPoints-1];

   return(perimeter);   
}


// double calculatePerimeterError(void) method       

/// This method returns the error between the perimeter goal and the arc
/// lenght of the function represented by the neural network.

double IsoperimetricProblem::calculatePerimeterError(void)    
{   
   double perimeter = calculatePerimeter();

   double perimeterError = perimeterGoal - perimeter;

   return(perimeterError);
}


// double calculateEvaluation(void) method

/// This method returns the evaluation of a multilayer perceptron 
/// for the isoperimetric problem with parametric equations. 

double IsoperimetricProblem::calculateEvaluation()
{
    // Control sentence 

   if(multilayerPerceptron == NULL)
   {
      std::cerr << std::endl
                << "Flood Error: IsoperimetricProblem class." << std::endl
                << "double calculateEvaluation(void) method." << std::endl
                << "Pointer to multilayer perceptron object cannot be NULL." << std::endl
                << std::endl;

        exit(1);
   }

   // Increment number of evaluations

   numberOfEvaluations++;

   double evaluation = 0.0;

   double perimeterError = calculatePerimeterError();
   double area = calculateArea();

   evaluation = -1.0*areaWeight*area + perimeterErrorWeight*pow(perimeterError,2);

   return(evaluation);
}


// void saveResults(char*) method

/// This method saves the values of the independent and dependent 
/// variables for the isoperimetric problem with parametric equations 
/// to a data file. 
///
/// <ol>
/// <li> Value of independent variable (t).
/// <li> Value of dependent variable x.
/// <li> Value of dependent variable y.
/// </ol>
///
/// @param filename Filename.

void IsoperimetricProblem::saveResults(char* filename)
{
   std::fstream file; 

   file.open(filename, std::ios::out);

   if(!file.is_open())
   {
      std::cout << std::endl
                << "Cannot open isoperimetric problem results data file." 
                << std::endl;
      
      exit(1);
   }
   else
   {
      if(display)
      {
         std::cout << std::endl
                   << "Saving isoperimetric problem results to data file..." 
                   << std::endl;
      }
   }

   // File header

   file << "% Flood Neural Network. Isoperimetric problem results data file." << std::endl
        << "% Column data:" << std::endl
        << "%  1 - Independent variable (t)" << std::endl
        << "%  2 - Dependent variable (x)" << std::endl
        << "%  3 - Dependent variable (y)" << std::endl
        << std::endl;

   // File data

   Vector<double> input(1);
   Vector<double> output(2);

   Vector<double> particularSolution(2);
   Vector<double> homogeneousSolution(2);

   double t = 0.0;
   double x = 0.0;
   double y = 0.0;

   int numberOfPoints = 1001;

   for(int i = 0; i < numberOfPoints; i++)
   {
      // Obtain t

      t = (double)i/(double)(numberOfPoints-1.0);

      // Obtain x and y

      input[0] = t;

      output = multilayerPerceptron->calculateOutput(input, *this);
      
      x = output[0];
      y = output[1];
      
      // Write t, x and y to file

      file << t << " " << x << " " << y << std::endl;   
   }   

   file.close();
}


// void print(void) method

/// This method prints to the screen useful information about the 
/// isoperimetric problem with parametric equations objective functional
/// of a multilayer perceptron:
///
/// <ul>
/// <li> Perimeter error.
/// <li> Area.
/// </ul>

void IsoperimetricProblem::print()
{
   double perimeterError = calculatePerimeterError();
   double area = calculateArea();

   std::cout << "Perimeter error: " << perimeterError << std::endl
             << "Area: " << area << std::endl;
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
