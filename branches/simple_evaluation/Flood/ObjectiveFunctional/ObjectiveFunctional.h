/****************************************************************************************************************/
/*                                                                                                              */
/*   Flood: An Open Source Neural Networks C++ Library                                                          */
/*   www.cimne.com/flood                                                                                        */
/*                                                                                                              */
/*   O B J E C T I V E   F U N C T I O N A L   C L A S S   H E A D E R                                          */
/*                                                                                                              */
/*   Roberto Lopez                                                                                              */
/*   International Center for Numerical Methods in Engineering (CIMNE)                                          */
/*   Technical University of Catalonia (UPC)                                                                    */
/*   Barcelona, Spain                                                                                           */
/*   E-mail: rlopez@cimne.upc.edu                                                                               */
/*                                                                                                              */
/****************************************************************************************************************/
 

#ifndef __OBJECTIVEFUNCTIONAL_H__
#define __OBJECTIVEFUNCTIONAL_H__

#include "../MultilayerPerceptron/MultilayerPerceptron.h"

namespace Flood
{

/// This abstract class represents the concept of objective functional for the multilayer perceptron. 
/// Any derived class must implement the calculateEvaluation(void) method.
///
/// @see MultilayerPerceptron.

class ObjectiveFunctional
{

public:

   // ENUMERATIONS

   /// Enumeration of available methods for obtaining a numerical value of epsilon to be used in numerical 
   /// differentiation.

   enum NumericalEpsilonMethod{Absolute, Relative};

   /// Enumeration of available methods for calculating any derivative of the objective function numerically.

   enum NumericalDifferentiationMethod{ForwardDifferences, CentralDifferences};

protected:

   /// Pointer to a multilayer perceptron object.

   MultilayerPerceptron* multilayerPerceptron;

   /// Number of calls to the calculateEvaluation(void) method.
   ///
   /// @see calculateEvaluation(void).

   int numberOfEvaluations;

   /// Epsilon value for the calculation of the objective function gradient with numerical differentiation.

   double numericalEpsilon;

   /// Display messages to screen. 

   bool display;  

   /// Numerical epsilon methods enumeration.

   NumericalEpsilonMethod numericalEpsilonMethod;

   /// Numerical differentiation methods enumeration.

   NumericalDifferentiationMethod numericalDifferentiationMethod;

public:

   // GENERAL CONSTRUCTOR

   ObjectiveFunctional(MultilayerPerceptron*);

   // DEFAULT CONSTRUCTOR

   ObjectiveFunctional(void);


   // DESTRUCTOR

   virtual ~ObjectiveFunctional(void);


   // METHODS

   // Get methods

   MultilayerPerceptron* getMultilayerPerceptron(void);
   double getNumericalEpsilon(void);
   int getNumberOfEvaluations(void);

   bool getDisplay(void);

   NumericalEpsilonMethod getNumericalEpsilonMethod(void);
   NumericalDifferentiationMethod getNumericalDifferentiationMethod(void);

   // Set methods

   void setMultilayerPerceptron(MultilayerPerceptron*);
   void setNumericalEpsilon(double);
   void setNumberOfEvaluations(int);

   void setDisplay(bool);

   void setNumericalEpsilonMethod(NumericalEpsilonMethod);
   void setNumericalDifferentiationMethod(NumericalDifferentiationMethod);

   // Objective methods

   /// This method returns the objective value of a multilayer perceptron.
   ///
   /// @see calculateGradient(void).

   virtual double calculateEvaluation(void) = 0;

   double calculatePotentialEvaluation(Vector<double>);

   // Objective function gradient methods

   /// This method returns the objective function gradient vector for a
   /// multilayer perceptron.
   ///
   /// @see calculateEvaluation(void).

   virtual Vector<double> calculateGradient(void);

   Vector<double> calculatePotentialGradient(Vector<double>);

   // Objective function hessian methods

   /// This method returns the objective function Hessian matrix for a multilayer perceptron.
   /// Please, do not use this method. It is not yet implemented!

   virtual Matrix<double> calculateHessian(void);

   Matrix<double> calculatePotentialHessian(Vector<double>);

   // Objective function inverse Hessian methods

   virtual Matrix<double> calculateInverseHessian(void);

   Matrix<double> calculateDFPInverseHessianApproximation(
   Vector<double>, Vector<double>, Matrix<double>, 
   Vector<double>, Vector<double>);

   Matrix<double> calculateBFGSInverseHessianApproximation(
   Vector<double>, Vector<double>, Matrix<double>, 
   Vector<double>, Vector<double>);

   virtual Vector<double> calculateVectorHessianProduct(Vector<double>);

   // Utility methods
   
   virtual void print(void);   
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
