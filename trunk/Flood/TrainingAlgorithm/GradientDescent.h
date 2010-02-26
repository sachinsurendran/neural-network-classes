/****************************************************************************************************************/
/*                                                                                                              */
/*   Flood: An Open Source Neural Networks C++ Library                                                          */
/*   www.cimne.com/flood                                                                                        */
/*                                                                                                              */
/*   G R A D I E N T   D E S C E N T   C L A S S   H E A D E R                                                  */
/*                                                                                                              */
/*   Roberto Lopez                                                                                              */
/*   International Center for Numerical Methods in Engineering (CIMNE)                                          */
/*   Technical University of Catalonia (UPC)                                                                    */
/*   Barcelona, Spain                                                                                           */
/*   E-mail: rlopez@cimne.upc.edu                                                                               */
/*                                                                                                              */
/****************************************************************************************************************/

#ifndef __GRADIENTDESCENT_H__
#define __GRADIENTDESCENT_H__

#include "TrainingAlgorithm.h"
#include "../ObjectiveFunctional/ObjectiveFunctional.h"


namespace Flood
{

/// This concrete class represents the gradient descent training algorithm for an objective functional of a 
/// multilayer perceptron.
///
/// @see TrainingAlgorithm.

class GradientDescent : public TrainingAlgorithm
{

public:

   // ENUMERATIONS

   /// Available training operators for obtaining the train rate.

   enum TrainRateMethod{Fixed, GoldenSection, BrentMethod};


private: 

   /// Initial train rate in line minimization.

   double firstTrainRate;

   /// Increase factor when bracketing a minimum.

   double bracketingFactor;

   /// Tolerance for the train rate.

   double trainRateTolerance;

   /// Train rate value at wich a warning message is written to the screen.

   double warningTrainRate;

   // Error train rate

   /// Train rate at wich the line minimization algorithm is assumed to be 
   /// unable to bracket a minimum.

   double errorTrainRate;

   /// True if the training direction history matrix is to be reserved, false otherwise.
   
   bool reserveTrainingDirectionHistory;

   /// True if the training direction norm history vector is to be reserved, false otherwise.

   bool reserveTrainingDirectionNormHistory;

   /// True if the training rate history vector is to be reserved, false otherwise.

   bool reserveTrainingRateHistory;
  
   /// Matrix containing the training direction history over the epochs.

   Matrix<double> trainingDirectionHistory;

   /// Vector containing the training direction norm history over the epochs.

   Vector<double> trainingDirectionNormHistory;

   /// Vector containing the training rate history over the epochs.

   Vector<double> trainingRateHistory;

   /// Variable containing the actual method used to obtain a suitable train rate. 

   TrainRateMethod trainRateMethod;

   // METHODS

   // Train rate methods

   double calculateGoldenSectionTrainRate(double, double, Vector<double>, Vector<double>);
   double calculateBrentMethodTrainRate(double, double, Vector<double>, Vector<double>);


public:

   // GENERAL CONSTRUCTOR

   GradientDescent(ObjectiveFunctional*);


   // DEFAULT CONSTRUCTOR

   GradientDescent(void); 


   // DESTRUCTOR

   virtual ~GradientDescent(void);


   // METHODS

   // Get methods

   TrainRateMethod getTrainRateMethod(void);

   double getFirstTrainRate(void);
   double getBracketingFactor(void);   
   double getTrainRateTolerance(void);
   double getWarningTrainRate(void);
   double getErrorTrainRate(void);

   bool getReserveTrainingDirectionHistory(void);
   bool getReserveTrainingDirectionNormHistory(void);
   bool getReserveTrainingRateHistory(void);

   Matrix<double> getTrainingDirectionHistory(void);
   Vector<double> getTrainingDirectionNormHistory(void);
   Vector<double> getTrainingRateHistory(void);

   // Set methods

   void setTrainRateMethod(TrainRateMethod);

   void setFirstTrainRate(double);
   void setBracketingFactor(double);   
   void setTrainRateTolerance(double);
   void setWarningTrainRate(double);
   void setErrorTrainRate(double);

   void setReserveTrainingDirectionHistory(bool);
   void setReserveTrainingDirectionNormHistory(bool);
   void setReserveTrainingRateHistory(bool);

   void setTrainingDirectionHistory(Matrix<double>);
   void setTrainingDirectionNormHistory(Vector<double>);
   void setTrainingRateHistory(Vector<double>);

   // Train methods

   void train(void);

   // Utility methods

   void print(void);

   void load(char*);
   void save(char*);

   void resizeTrainingHistory(int);

//   void saveTrainingDirectionHistory(char*);
//   void saveTrainingDirectionNormHistory(char*);
//   void saveTrainingRateHistory(char*);

   void setReserveAllTrainingHistory(bool);

   void saveTrainingHistory(char*);


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
