/****************************************************************************************************************/
/*                                                                                                              */
/*   Flood: An Open Source Neural Networks C++ Library                                                          */
/*   www.cimne.com/flood                                                                                        */
/*                                                                                                              */
/*   G E O D E S I C   P R O B L E M   C L A S S   H E A D E R                                                  */
/*                                                                                                              */ 
/*   Roberto Lopez                                                                                              */ 
/*   International Center for Numerical Methods in Engineering (CIMNE)                                          */
/*   Technical University of Catalonia (UPC)                                                                    */
/*   Barcelona, Spain                                                                                           */
/*   E-mail: rlopez@cimne.upc.es                                                                                */ 
/*                                                                                                              */
/****************************************************************************************************************/

#ifndef __GEODESICPROBLEM_H__
#define __GEODESICPROBLEM_H__

#include "ObjectiveFunctional.h"
#include "../Utilities/OrdinaryDifferentialEquations.h"

namespace Flood
{

/// This class represents the objective functional of a multilayer perceptron for the geodesic problem.
/// The geodesic problem for the multilayer perceptron is a variational problem with one input and one output 
/// variables, two boundary conditions and where evaluation of the objective functional is obtained by 
/// integrating a function. 
///
/// @see ObjectiveFunctional.
/// @see OrdinaryDifferentialEquations. 

class GeodesicProblem : public ObjectiveFunctional
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

   /// Tolerance in the Runge-Kutta-Fehlberg method.

   double tolerance;

   /// Number of points to reserve for the Runge-Kutta-Fehlberg method.

   int initialSize;

   /// Ordinary differential equations object.

   OrdinaryDifferentialEquations ordinaryDifferentialEquations;


public:

   // GENERAL CONSTRUCTOR

   GeodesicProblem(MultilayerPerceptron*);


   // DEFAULT CONSTRUCTOR

   GeodesicProblem(void);


   // DESTRUCTOR

   virtual ~GeodesicProblem(void);


   // METHODS

   // Get methods

   double getXA(void);
   double getYA(void);
   double getXB(void);
   double getYB(void);

   double getTolerance(void);
   int getInitialSize(void);

   // Set methods

   void setXA(double);
   void setYA(double);
   void setXB(double);
   void setYB(double);

   void setTolerance(double);
   void setInitialSize(int);

   // Arc length methods

   Vector<double> calculateParticularSolution(Vector<double>);
   Vector<double> calculateHomogeneousSolution(Vector<double>);

   Vector<double> calculateParticularSolutionDerivative(Vector<double>);
   Vector<double> calculateHomogeneousSolutionDerivative(Vector<double>);

   // Arc length methods

   double calculateArcLengthIntegrand(double, double);
   double calculateArcLength(void);

   // Objective evaluation methods

   double calculateEvaluation(void);

   // Utility methods

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
