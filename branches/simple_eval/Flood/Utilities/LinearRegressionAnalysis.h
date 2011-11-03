/****************************************************************************************************************/
/*                                                                                                              */
/*   Flood: An Open Source Neural Networks C++ Library                                                          */
/*   www.cimne.com/flood                                                                                        */
/*                                                                                                              */
/*   L I N E A R   R E G R E S S I O N   A N A L Y S I S   C L A S S   H E A D E R                              */
/*                                                                                                              */
/*   Roberto Lopez                                                                                              */ 
/*   International Center for Numerical Methods in Engineering (CIMNE)                                          */
/*   Technical University of Catalonia (UPC)                                                                    */
/*   Barcelona, Spain                                                                                           */
/*   E-mail: rlopez@cimne.upc.edu                                                                               */ 
/*                                                                                                              */
/****************************************************************************************************************/

#ifndef __LINEARREGRESSIONANALYSIS_H__
#define __LINEARREGRESSIONANALYSIS_H__

#include "../Utilities/Vector.h"
#include "../Utilities/Matrix.h"

#include "../MultilayerPerceptron/MultilayerPerceptron.h"
#include "../Utilities/InputTargetDataSet.h"

namespace Flood
{

/// This class contains methods for performing linear regression analysis. 
/// This is a validation technique for data modeling problems which consists on
/// performing a linear regression between the outputs from a multilayer
/// perceptron and the corresponding targets from an input-target data set.
/// Linear regression analysis leads to 3 parameters for each output variable. 
/// The first two, a and b, correspond to the y-intercept and the slope 
/// of the best linear regression relating network outputs to targets, that is,
/// y = a + b*x. 
/// If we had a perfect fit (outputs exactly equal to targets), 
/// the y-intercept, a, would be 0 and the slope, b, would be 1. 
/// The third parameter is the correlation coefficient (R-value) between 
/// the outputs and the targets. 
/// If the correlation coefficient, R, is equal to 1, then there is 
/// perfect correlation. 
///
/// @see MultilayerPerceptron
/// @see InputTargetDataSet

class LinearRegressionAnalysis
{

private: 

   // FIELDS

   /// Pointer to a multilayer perceptron object.

   MultilayerPerceptron* multilayerPerceptron;

   /// Pointer to an input-target data set object.

   InputTargetDataSet* inputTargetDataSet;

   /// Display messages to screen.
   
   bool display;

public:  

   // GENERAL CONSTRUCTOR

   LinearRegressionAnalysis(MultilayerPerceptron*, InputTargetDataSet*);


   // DEFAULT CONSTRUCTOR

   LinearRegressionAnalysis(void);


   // DESTRUCTOR

   virtual ~LinearRegressionAnalysis();


   // METHODS

   // Get methods

   MultilayerPerceptron* getMultilayerPerceptron(void);
   InputTargetDataSet* getInputTargetDataSet(void);
   
   bool getDisplay(void);

   // Set methods

   void setMultilayerPerceptron(MultilayerPerceptron*);
   void setInputTargetDataSet(InputTargetDataSet*);

   void setDisplay(bool);

   // Output data methods
   
   Matrix<double> getOutputData(void);
   
   // Regression parameters methods
   
   Matrix<double> calculateRegressionParameters(void);

   // Utility methods

   void printResults(void);

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
