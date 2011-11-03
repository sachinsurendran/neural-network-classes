/****************************************************************************************************************/
/*                                                                                                              */
/*   Flood: An Open Source Neural Networks C++ Library                                                          */
/*   www.cimne.com/flood                                                                                        */
/*                                                                                                              */
/*   I S O P E R I M E T R I C   P R O B L E M   C L A S S   H E A D E R                                        */
/*                                                                                                              */
/*   Roberto Lopez                                                                                              */
/*   International Center for Numerical Methods in Engineering (CIMNE)                                          */
/*   Technical University of Catalonia (UPC)                                                                    */
/*   Barcelona, Spain                                                                                           */
/*   E-mail: rlopez@cimne.upc.es                                                                                */
/*                                                                                                              */
/****************************************************************************************************************/

#ifndef __ISOPERIMETRICPROBLEM_H__
#define __ISOPERIMETRICPROBLEM_H__

#include "ObjectiveFunctional.h"
#include "../Utilities/OrdinaryDifferentialEquations.h"

namespace Flood
{

/// This class represents the objective functional of a multilayer perceptron for the isoperimetric problem with 
/// parametric equations. 
/// This is a variational problem with one input and two output variables, four boundary conditions and one 
/// constraint in the objective functional, which is evaluated by integrating two functions. 
///
/// @see ObjectiveFunctional.
/// @see OrdinaryDifferentialEquations.

class IsoperimetricProblem : public ObjectiveFunctional
{

private:

   /// Closed curve perimeter goal. 

   double perimeterGoal;

   /// Perimeter error term weight in the objective functional. 

   double perimeterErrorWeight;

   /// Area term weight in the objective functional. 

   double areaWeight;

   /// Tolerance in the Runge-Kutta-Fehlberg method.

   double tolerance;

   /// Number of points to reserve for the Runge-Kutta-Fehlberg method.

   int initialSize;

   /// Ordinary differential equations object

   OrdinaryDifferentialEquations ordinaryDifferentialEquations;

public:

   // GENERAL CONSTRUCTOR

   IsoperimetricProblem(MultilayerPerceptron*);

   // DEFAULT CONSTRUCTOR

   IsoperimetricProblem(void);

   // DESTRUCTOR

   virtual ~IsoperimetricProblem(void);

   // METHODS
   
   // Get methods

   double getPerimeterGoal(void);

   double getAreaWeight(void);
   double getPerimeterErrorWeight(void);

   double getTolerance(void);
   int getInitialSize(void);

   // Set methods

   void setPerimeterGoal(double);

   void setAreaWeight(double);
   void setPerimeterErrorWeight(double);
   
   void setTolerance(double);
   void setInitialSize(int);

   // Boundary conditions methods

   Vector<double> calculateParticularSolution(Vector<double>);
   Vector<double> calculateHomogeneousSolution(Vector<double>);

   Vector<double> calculateParticularSolutionDerivative(Vector<double>);
   Vector<double> calculateHomogeneousSolutionDerivative(Vector<double>);

   // Area methods

   double calculateAreaIntegrand(double, double);
   double calculateArea(void);

   // Perimeter methods

   double calculatePerimeterIntegrand(double, double);      
   double calculatePerimeter(void);       
   double calculatePerimeterError(void);   
      
   // Evaluation methods

   double calculateEvaluation(void);

   // Utility methods

   void saveResults(char*);
   void print(void);
   
};

}
#endif


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
