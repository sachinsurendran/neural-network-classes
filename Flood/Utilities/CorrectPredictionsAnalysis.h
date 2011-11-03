/****************************************************************************************************************/
/*                                                                                                              */
/*   Flood: An Open Source Neural Networks C++ Library                                                          */
/*   www.cimne.com/flood                                                                                        */
/*                                                                                                              */
/*   C O R R E C T   P R E D I C T I O N S   A N A L Y S I S   C L A S S   H E A D E R                          */
/*                                                                                                              */
/*   Roberto Lopez                                                                                              */ 
/*   International Center for Numerical Methods in Engineering (CIMNE)                                          */
/*   Technical University of Catalonia (UPC)                                                                    */
/*   Barcelona, Spain                                                                                           */
/*   E-mail: rlopez@cimne.upc.edu                                                                               */ 
/*                                                                                                              */
/****************************************************************************************************************/

#ifndef __CORRECTPREDICTIONSANALYSIS_H__
#define __CORRECTPREDICTIONSANALYSIS_H__

#include "../Utilities/Vector.h"
#include "../Utilities/Matrix.h"

#include "../MultilayerPerceptron/MultilayerPerceptron.h"
#include "../Utilities/InputTargetDataSet.h"

namespace Flood
{

/// This class contains methods for performing correct predictions analysis, 
/// which is a validation technique for pattern recognition problems.
/// It consists on calculating the ratio of correct predictions made by a multilayer perceptron on an independent
/// testing set. 
/// ///
/// @see MultilayerPerceptron
/// @see InputTargetDataSet

class CorrectPredictionsAnalysis
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

   CorrectPredictionsAnalysis(MultilayerPerceptron*, InputTargetDataSet*);


   // DEFAULT CONSTRUCTOR

   CorrectPredictionsAnalysis(void);


   // DESTRUCTOR

   virtual ~CorrectPredictionsAnalysis();


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
   
   // Correct predictions ratio methods
   
   double calculateCorrectPredictionsRatio(void);

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
