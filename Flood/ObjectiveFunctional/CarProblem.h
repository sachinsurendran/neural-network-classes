/****************************************************************************************************************/
/*                                                                                                              */
/*   Flood: An Open Source Neural Networks C++ Library                                                          */
/*   www.cimne.com/flood                                                                                        */
/*                                                                                                              */
/*   C A R   P R O B L E M   C L A S S   H E A D E R                                                            */
/*                                                                                                              */
/*   Roberto Lopez                                                                                              */
/*   International Center for Numerical Methods in Engineering (CIMNE)                                          */
/*   Technical University of Catalonia (UPC)                                                                    */
/*   Barcelona, Spain                                                                                           */
/*   E-mail: rlopez@cimne.upc.es                                                                                */
/*                                                                                                              */
/****************************************************************************************************************/


#ifndef __CARPROBLEM_H__
#define __CARPROBLEM_H__

#include "ObjectiveFunctional.h"
#include "../Utilities/OrdinaryDifferentialEquations.h"

namespace Flood
{

/// This class represents the objective functional of a multilayer perceptron for the car problem.  
/// The car problem for the multilayer perceptron is an optimal control problem with two controls and two state 
/// variables.
/// It is defined by an objective functional with two constraints and requiring the integration of a system of 
/// ordinary differential equations. 
///
/// @see ObjectiveFunctional.

class CarProblem : public ObjectiveFunctional
{

private:

   /// Initial car position.

   double initialPosition;

   /// Initial car velocity.

   double initialVelocity;

   /// Desired final position of the car.

   double finalPositionGoal;

   /// Desired final velocity of the car.

   double finalVelocityGoal;

   /// Car maximum acceleration (bound).

   double maximumAcceleration;

   /// Car maximum deceleration (bound).

   double maximumDeceleration;

   /// Weight for the final position error term in the objective functional.

   double finalPositionErrorWeight;

   /// Weight for the final velocity error term in the objective functional.

   double finalVelocityErrorWeight;

   /// Weight for the final time term in the objective functional.

   double finalTimeWeight;

   /// Ordinary differential equations object.   

   OrdinaryDifferentialEquations ordinaryDifferentialEquations;

   /// Tolerance of integration in Runge-Kutta-Fehlberg method.

   double tolerance;

   /// Initial size of solution vectors for the Runge-Kutta-Fehlberg method.

   int initialSize;


public:

     // GENERAL CONSTRUCTOR

   CarProblem(MultilayerPerceptron*);

     // DEFAULT CONSTRUCTOR

   CarProblem(void);


   // DESTRUCTOR

   virtual ~CarProblem(void);


   // METHODS

   // Get methods

   double getInitialPosition(void);
   double getInitialVelocity(void);

   double getFinalPositionGoal(void);
   double getFinalVelocityGoal(void);

   double getMaximumAcceleration(void);
   double getMaximumDeceleration(void);
   
   double getFinalPositionErrorWeight(void);
   double getFinalVelocityErrorWeight(void);

   double getFinalTimeWeight(void);

   double getTolerance(void);
   int getInitialSize(void);

   // Set methods

   void setInitialPosition(double);
   void setInitialVelocity(double);

   void setFinalPositionGoal(double);
   void setFinalVelocityGoal(double);

   void setMaximumAcceleration(double);
   void setMaximumDeceleration(double);

   void setFinalPositionErrorWeight(double);
   void setFinalVelocityErrorWeight(double);

   void setFinalTimeWeight(double);

   void setTolerance(double);
   void setInitialSize(int);

   // State equation methods

   double calculatePositionDot(double, double, double);
   double calculateVelocityDot(double, double, double);

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
