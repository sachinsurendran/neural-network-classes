/****************************************************************************************************************/
/*                                                                                                              */
/*   Flood: An Open Source Neural Networks C++ Library                                                          */
/*   www.cimne.com/flood                                                                                        */
/*                                                                                                              */
/*   Q U A S I - N E W T O N   M E T H O D    C L A S S   H E A D E R                                           */
/*                                                                                                              */
/*   Roberto Lopez                                                                                              */
/*   International Center for Numerical Methods in Engineering (CIMNE)                                          */
/*   Technical University of Catalonia (UPC)                                                                    */
/*   Barcelona, Spain                                                                                           */
/*   E-mail: rlopez@cimne.upc.edu                                                                               */
/*                                                                                                              */
/****************************************************************************************************************/

#ifndef __QUASINEWTONMETHOD_H__
#define __QUASINEWTONMETHOD_H__

#include "TrainingAlgorithm.h"
#include "../ObjectiveFunctional/ObjectiveFunctional.h"

namespace Flood
{

/// This concrete class represents a quasi-Newton training algorithm for an objective functional of a multilayer 
/// perceptron.
///
/// @see TrainingAlgorithm.

class QuasiNewtonMethod : public TrainingAlgorithm
{

public:

   // ENUMERATIONS

   /// Enumeration of the available training operators for obtaining the approximation to the inverse Hessian.

   enum InverseHessianApproximationMethod{DFP, BFGS};

   /// Enumeration of the available operators for obtaining the train rate.

   enum TrainRateMethod{Fixed, GoldenSection, BrentMethod};

private: 

   // FIELDS

   /// Initial train rate in line minimization.

   double firstTrainRate;

   /// Increase factor when bracketing a minimum.

   double bracketingFactor;

   /// Tolerance in line minimization.

   double trainRateTolerance;

   // Warning train rate

   /// Train rate at wich a warning message is written to the screen during line minimization.

   double warningTrainRate;

   // Error train rate

   /// Train rate at wich the line minimization algorithm is assumed to be unable to bracket a minimum.

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

   /// Variable containing the actual method used to obtain an approximation for the inverse Hessian matrix.

   InverseHessianApproximationMethod inverseHessianApproximationMethod;

   /// Variable containing the actual method used to obtain a suitable train rate. 

   TrainRateMethod trainRateMethod;


public:

   // GENERAL CONSTRUCTOR

   QuasiNewtonMethod(ObjectiveFunctional*);


   // DEFAULT CONSTRUCTOR

   QuasiNewtonMethod(void);


   // DESTRUCTOR

   virtual ~QuasiNewtonMethod(void);


   // METHODS

   // Get methods

   InverseHessianApproximationMethod getInverseHessianApproximationMethod(void);
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

   void setInverseHessianApproximationMethod(InverseHessianApproximationMethod);
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

   double calculateGoldenSectionTrainRate(double, double, Vector<double>, Vector<double>);
   double calculateBrentMethodTrainRate(double, double, Vector<double>, Vector<double>);

   void train(void);

   // Utility methods

   void print(void);   

   void load(char*);
   void save(char*);
   
   void resizeTrainingHistory(int);

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
