/****************************************************************************************************************/
/*                                                                                                              */
/*   Flood: An Open Source Neural Networks C++ Library                                                          */
/*   www.cimne.com/flood                                                                                        */
/*                                                                                                              */
/*   B R A C H I S T O C H R O N E   P R O B L E M   C L A S S   H E A D E R                                    */
/*                                                                                                              */
/*   Roberto Lopez                                                                                              */
/*   International Center for Numerical Methods in Engineering (CIMNE)                                          */
/*   Technical University of Catalonia (UPC)                                                                    */
/*   Barcelona, Spain                                                                                           */
/*   E-mail: rlopez@cimne.upc.edu                                                                               */
/*                                                                                                              */
/****************************************************************************************************************/


#ifndef __BRACHISTOCHRONEPROBLEM_H__
#define __BRACHISTOCHRONEPROBLEM_H__

#include "ObjectiveFunctional.h"
#include "../Utilities/OrdinaryDifferentialEquations.h"

namespace Flood
{

/// This class represents the objective functional of a multilayer perceptron for the brachistochrone problem. 
/// The brachistochrone problem for the multilayer perceptron can be stated as a variational problem with one 
/// input and one output variables, two boundary conditions and an objective functional defined by an improper 
/// integral.
///
/// @see ObjectiveFunctional.
/// @see OrdinaryDifferentialEquations.

class BrachistochroneProblem : public ObjectiveFunctional
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

   /// Ordinary differential equations object used to evaluate the descent time. 

   OrdinaryDifferentialEquations ordinaryDifferentialEquations;

   /// Tolerance in the Runge-Kutta-Fehlberg method used to evaluate the descent time.  

   double tolerance;

   /// Initial size to be reserved in the Runge-Kutta-Fehlberg method used to evaluate the descent time.  

   int initialSize;


public:

   // GENERAL CONSTRUCTOR

   BrachistochroneProblem(MultilayerPerceptron*);


   // DEFAULT CONSTRUCTOR

   BrachistochroneProblem(void);


   // DESTRUCTOR

   virtual ~BrachistochroneProblem(void);


   // METHODS

   // Get methods

   double getXa(void);
   double getYa(void);
   double getXb(void);
   double getYb(void);

   double getTolerance(void);
   int getInitialSize(void);

   // Set methods

   void setXa(double);
   void setYa(double);
   void setXb(double);
   void setYb(double);
   void setProblem(double, double, double, double);

   void setTolerance(double);
   void setInitialSize(int);

   // Boundary conditions methods

   Vector<double> calculateParticularSolution(Vector<double>);
   Vector<double> calculateHomogeneousSolution(Vector<double>);

   Vector<double> calculateParticularSolutionDerivative(Vector<double>);
   Vector<double> calculateHomogeneousSolutionDerivative(Vector<double>);

   // Descent time methods

   double calculateDescentTimeIntegrand(double, double);
   double calculateDescentTime(void);

   // Objective evaluation methods

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
