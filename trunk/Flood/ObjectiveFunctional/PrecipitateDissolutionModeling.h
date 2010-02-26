/****************************************************************************************************************/
/*                                                                                                              */
/*   Flood: An Open Source Neural Networks C++ Library                                                          */
/*   www.cimne.com/flood                                                                                        */
/*                                                                                                              */
/*   P R E C I P I T A T E   D I S S O L U T I O N   M O D E L I N G   C L A S S   H E A D E R                  */
/*                                                                                                              */
/*   Roberto Lopez                                                                                              */ 
/*   International Center for Numerical Methods in Engineering (CIMNE)                                          */
/*   Technical University of Catalonia (UPC)                                                                    */
/*   Barcelona, Spain                                                                                           */
/*   E-mail: rlopez@cimne.upc.es                                                                                */ 
/*                                                                                                              */
/****************************************************************************************************************/

#ifndef __PRECIPITATEDISSOLUTIONMODELING_H__
#define __PRECIPITATEDISSOLUTIONMODELING_H__

#include "ObjectiveFunctional.h"
#include "../MultilayerPerceptron/MultilayerPerceptron.h"

#include <string>

namespace Flood
{

/// This class represents the objective functional of a multilayer perceptron for modeling the dissolution rate 
/// of hardening precipitates in aluminium alloys.
///
/// @see ObjectiveFunctional.

class PrecipitateDissolutionModeling : public ObjectiveFunctional
{

private:

   /// Number of samples in Vickers hardness test. 

   int numberOfSamples;

   /// Vector containing the time measurements for each sample in the test. 

   Vector<double> timeData;

   /// Vector containing the temperature measurements for each sample in the test. 

   Vector<double> temperatureData;

   /// Vector containing the Vickers hardness measurements for each sample in the test. 

   Vector<double> vickersHardnessData;

   /// Universal gas constant R.

   double universalGasConstant;

   /// Reference temperature of aluminium alloy. 

   double referenceTemperature;

   /// Reference time of aluminium alloy. 

   double referenceTime;

   /// Minimum Vickers hardness of aluminium alloy.

   double minimumVickersHardness;

   /// Maximum Vickers hardness of aluminium alloy.

   double maximumVickersHardness;

   /// Minkowski parameter in the objective functional expression.

   double minkowskiParameter;

   /// Weight value of the Minkowski error term in the objective functional expression. 

   double minkowskiErrorWeight;

   /// Weight value of the regularization term in the objective functional expression. 

   double regularizationWeight;

public:

     // GENERAL CONSTRUCTOR

   PrecipitateDissolutionModeling(MultilayerPerceptron*);


     // DEFAULT CONSTRUCTOR

   PrecipitateDissolutionModeling(void);


   // DESTRUCTOR

   virtual ~PrecipitateDissolutionModeling(void);

   
   // METHODS

   // Evaluation methods

   // Get methods

   double getMinimumVickersHardness(void);
   double getMaximumVickersHardness(void);

   double getReferenceTemperature(void);
   double getReferenceTime(void);

   double getMinkowskiParameter(void);

   double getMinkowskiErrorWeight(void);
   double getRegularizationWeight(void);

   // Set methods

   void setMinimumVickersHardness(double);
   void setMaximumVickersHardness(double);

   void setReferenceTemperature(double);
   void setReferenceTime(double);

   void setMinkowskiParameter(double);

   void setMinkowskiErrorWeight(double);
   void setRegularizationWeight(double);

   // Evaluation methods

   void loadVickersHardnessTest(char*);

   double getFullDissolutionTime(double);

   double getVolumetricFraction(double);

   double getVickersHardness(double);

   Vector<double> getVolumetricFractionData(void);
   Vector<double> getNormalizedTimeData(void);

   double calculateEvaluation(void);
  
   // Utilities methods

   void print(void);

   void printVickersHardnessTest(void);

   void savePrecipitateDissolutionModel(char*);
   void saveVickersHardnessModel(char*);

   void saveInverseVickersHardnessTest(char*);
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
