/****************************************************************************************************************/
/*                                                                                                              */
/*   Flood: An Open Source Neural Networks C++ Library                                                          */
/*   www.cimne.com/flood                                                                                        */
/*                                                                                                              */
/*   M I N I M U M   D R A G   P R O B L E M   C L A S S   H E A D E R                                          */
/*                                                                                                              */
/*   Roberto Lopez                                                                                              */
/*   International Center for Numerical Methods in Engineering (CIMNE)                                          */
/*   Technical University of Catalonia (UPC)                                                                    */
/*   Barcelona, Spain                                                                                           */
/*   E-mail: rlopez@cimne.upc.edu                                                                               */
/*                                                                                                              */
/****************************************************************************************************************/


#ifndef __MINIMUMDRAGPROBLEM_H__
#define __MINIMUMDRAGPROBLEM_H__

#include "ObjectiveFunctional.h"
#include "../Utilities/OrdinaryDifferentialEquations.h"

namespace Flood
{

/// This class represents the objective functional of a multilayer perceptron for the minimum drag problem.
/// The minimum drag problem for the multilayer perceptron is an optimal shape design problem with one input and 
/// one output variables, besides two boundary conditions. 
/// It is defined by an unconstrained objective fucntional requiring the integration of a function. 
///
/// @see ObjectiveFunctional.

class MinimumDragProblem : public ObjectiveFunctional
{

private:

   /// Ordinary differential equations object.    

   OrdinaryDifferentialEquations ordinaryDifferentialEquations;

   /// Tolerance of integration in Runge-Kutta-Fehlberg method.

   double tolerance;

   /// Initial size of solution vectors for the Runge-Kutta-Fehlberg method.

   int initialSize;


public:

   // GENERAL CONSTRUCTOR

   MinimumDragProblem(MultilayerPerceptron*);


   // DEFAULT CONSTRUCTOR

   MinimumDragProblem(void);


   // DESTRUCTOR

   virtual ~MinimumDragProblem(void);


   // METHODS

   // Get methods

   double getTolerance(void);
   int getInitialSize(void);

   // Set methods

   void setTolerance(double);
   void setInitialSize(int);

   // Boundary conditions methods

   Vector<double> calculateParticularSolution(Vector<double>);
   Vector<double> calculateHomogeneousSolution(Vector<double>);

   Vector<double> calculateParticularSolutionDerivative(Vector<double>);
   Vector<double> calculateHomogeneousSolutionDerivative(Vector<double>);   

   // Objective functional evaluation methods

   double getDragIntegrand(double, double);  
   double getDrag(void);

   double calculateEvaluation(void);

   // Utility methods

   void saveResults(char*);
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
