/****************************************************************************************************************/
/*                                                                                                              */
/*   Flood: An Open Source Neural Networks C++ Library                                                          */
/*   www.cimne.com/flood                                                                                        */
/*                                                                                                              */
/*   C A T E N A R Y   P R O B L E M   C L A S S   H E A D E R                                                  */
/*                                                                                                              */ 
/*   Roberto Lopez                                                                                              */ 
/*   International Center for Numerical Methods in Engineering (CIMNE)                                          */
/*   Technical University of Catalonia (UPC)                                                                    */
/*   Barcelona, Spain                                                                                           */
/*   E-mail: rlopez@cimne.upc.es                                                                                */ 
/*                                                                                                              */
/****************************************************************************************************************/

#ifndef __CATENARYPROBLEM_H__
#define __CATENARYPROBLEM_H__

#include "ObjectiveFunctional.h"
#include "../Utilities/OrdinaryDifferentialEquations.h"

namespace Flood
{

/// This class represents the objective functional of a multilayer perceptron for the catenary problem.
/// This is a variational problem with one input and one output variables, two boundary conditions, and an 
/// objective functional with one constraint and requiring the integration of two functions. 
///
/// @see ObjectiveFunctional.

class CatenaryProblem : public ObjectiveFunctional
{

private:

   /// Abscissa value for the A point. 

   double xA;

   /// Ordinate value for the A point.

   double yA;

   /// Abscissa value for the B point.

   double xB;

   /// Ordinate value for the B point. 

   double yB;

   /// Chain length goal. 

   double lengthGoal;

   /// Potential energy term weight in the objective functional. 

   double potentialEnergyWeight;

   /// Chain length error term weight in the objective functional. 

   double lengthErrorWeight;

   /// Tolerance in the Runge-Kutta-Fehlberg method.

   double tolerance;

   /// Number of points to reserve for the Runge-Kutta-Fehlberg method.

   int initialSize;

   /// Ordinary differential equations object

   OrdinaryDifferentialEquations ordinaryDifferentialEquations;

public:

   // GENERAL CONSTRUCTOR

   CatenaryProblem(MultilayerPerceptron*);


   // DEFAULT CONSTRUCTOR

   CatenaryProblem(void);


   // DESTRUCTOR

   virtual ~CatenaryProblem(void);


   // METHODS

   // Get methods

   double getXA(void);
   double getYA(void);
   double getXB(void);
   double getYB(void);

   double getLengthGoal(void);

   double getPotentialEnergyWeight(void);
   double getLengthErrorWeight(void);

   double getTolerance(void);
   int getInitialSize(void);

   // Set methods

   void setXA(double);
   void setYA(double);
   void setXB(double);
   void setYB(double);

   void setLengthGoal(double);

   void setPotentialEnergyWeight(double);
   void setLengthErrorWeight(double);

   void setTolerance(double);
   void setInitialSize(int);

   // Boundary conditions methods

   Vector<double> calculateParticularSolution(Vector<double>);
   Vector<double> calculateHomogeneousSolution(Vector<double>);

   Vector<double> calculateParticularSolutionDerivative(Vector<double>);
   Vector<double> calculateHomogeneousSolutionDerivative(Vector<double>);

   // Potential energy methods

   double calculatePotentialEnergyIntegrand(double, double);
   double calculatePotentialEnergy(void);

   // Lenght methods

   double calculateLengthIntegrand(double, double);
   double calculateLength(void);
   double calculateLengthError(void);

   // Objective evaluation methods

   double calculateEvaluation(void);

   // Utility methods

   void print(void);
   void saveResults(char*);

};

}

#endif


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
