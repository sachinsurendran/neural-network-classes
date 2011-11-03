/******************************************************************************/
/*                                                                            */
/*   Flood: An Open Source Neural Networks C++ Library                        */
/*   www.cimne.com/flood                                                      */
/*                                                                            */
/*   T R A I N I N G   A L G O R I T H M   C L A S S   H E A D E R            */
/*                                                                            */
/*   Roberto Lopez                                                            */
/*   International Center for Numerical Methods in Engineering (CIMNE)        */
/*   Technical University of Catalonia (UPC)                                  */
/*   Barcelona, Spain                                                         */
/*   E-mail: rlopez@cimne.upc.edu                                             */
/*                                                                            */
/******************************************************************************/


#ifndef __TRAININGALGORITHM_H__
#define __TRAININGALGORITHM_H__

#include "../MultilayerPerceptron/MultilayerPerceptron.h"
#include "../ObjectiveFunctional/ObjectiveFunctional.h"


namespace Flood
{

/// This abstract class represents the concept of training algorithm
/// for the multilayer perceptron. Any derived class must implement the
/// train(void), print(void) and save(char*) methods.
///
/// @see ObjectiveFunctional.

class TrainingAlgorithm
{

protected:

   // FIELDS

   /// Pointer to an objective functional for a multilayer perceptron object.

   ObjectiveFunctional* objectiveFunctional;

   /// Goal value for the objective. It is used as a train stopping criterion.
   ///
   /// @see train(void).

   double evaluationGoal;

   /// Goal value for the norm of the objective function gradient.
   /// It is used as a train stopping criterion.
   ///
   /// @see train(void).

   double gradientNormGoal;

   /// Maximum training time.
   /// It is used as a train stopping criterion.
   ///
   /// @see train(void).

   double maximumTime;

   /// Minimum evaluation improvement between two successive epochs.
   /// It is used as a train stopping criterion.
   ///
   /// @see train(void).

   double minimumImprovement;

   /// Maximum number of epochs to train.
   /// It is used as a train stopping criterion.
   ///
   /// @see train(void).

   int maximumNumberOfEpochs;

   /// Value for the free parameters norm at which a warning message is written to the screen. 

   double warningFreeParametersNorm;

   /// Value for the gradient norm at which a warning message is written to the screen. 

   double warningGradientNorm;
   
   /// Display messages to screen.

   bool display;

   /// Number of epochs between the training showing progress.
   ///
   /// @see train(void).

   int displayPeriod;

   /// True if the elapsed time history vector is to be reserved, false otherwise.

   bool reserveElapsedTimeHistory;

   /// True if the free parameters history matrix is to be reserved, false otherwise.

   bool reserveFreeParametersHistory;

   /// True if the free parameters norm history vector is to be reserved, false otherwise.

   bool reserveFreeParametersNormHistory;

   /// True if the evaluation history vector is to be reserved, false otherwise.

   bool reserveEvaluationHistory;

   /// True if the gradient history matrix is to be reserved, false otherwise.

   bool reserveGradientHistory;

   /// True if the gradient norm history vector is to be reserved, false otherwise.

   bool reserveGradientNormHistory;

   /// Matrix containing the free parameters history over the training epochs.

   Matrix<double> freeParametersHistory;

   /// Vector containing the free parameters norm history over the training epochs.

   Vector<double> freeParametersNormHistory;

   /// Vector containing the evaluation history over the training epochs.

   Vector<double> evaluationHistory;

   /// Matrix containing the gradient history over the training epochs.

   Matrix<double> gradientHistory;

   /// Vector containing the gradient norm history over the training epochs.

   Vector<double> gradientNormHistory;

   /// Vector containing the elapsed time history over the training epochs.

   Vector<double> elapsedTimeHistory;

public:

   // GENERAL CONSTRUCTOR

   TrainingAlgorithm(ObjectiveFunctional*);


   // DEFAULT CONSTRUCTOR

   TrainingAlgorithm(void);


   // DESTRUCTOR

   virtual ~TrainingAlgorithm(void);


   // METHODS

   // Get methods

   ObjectiveFunctional* getObjectiveFunctional(void);

   double getEvaluationGoal(void);
   double getGradientNormGoal(void);
   double getMaximumTime(void);
   double getMinimumImprovement(void);
   int getMaximumNumberOfEpochs(void);

   double getWarningFreeParametersNorm(void);
   double getWarningGradientNorm(void);
   
   bool getDisplay(void);
   int getDisplayPeriod(void);

   bool getReserveFreeParametersHistory(void);
   bool getReserveFreeParametersNormHistory(void);
   bool getReserveEvaluationHistory(void);
   bool getReserveGradientHistory(void);
   bool getReserveGradientNormHistory(void);
   bool getReserveElapsedTimeHistory(void);

   Matrix<double> getFreeParametersHistory(void);
   Vector<double> getFreeParametersNormHistory(void);
   Vector<double> getEvaluationHistory(void);
   Matrix<double> getGradientHistory(void);
   Vector<double> getGradientNormHistory(void);
   Vector<double> getElapsedTimeHistory(void);

   // Set methods

   void setObjectiveFunctional(ObjectiveFunctional*);

   void setEvaluationGoal(double);
   void setGradientNormGoal(double);
   void setMaximumTime(double);
   void setMinimumImprovement(double);
   void setMaximumNumberOfEpochs(int);

   void setWarningFreeParametersNorm(double);
   void setWarningGradientNorm(double);

   void setDisplay(bool);
   void setDisplayPeriod(int);

   void setReserveFreeParametersHistory(bool);
   void setReserveFreeParametersNormHistory(bool);
   void setReserveEvaluationHistory(bool);
   void setReserveGradientHistory(bool);
   void setReserveGradientNormHistory(bool);
   void setReserveElapsedTimeHistory(bool);

   void setFreeParametersHistory(Matrix<double>);
   void setFreeParametersNormHistory(Vector<double>);
   void setEvaluationHistory(Vector<double>);
   void setGradientHistory(Matrix<double>);
   void setGradientNormHistory(Vector<double>);
   void setElapsedTimeHistory(Vector<double>);

   // Train methods

   /// This method trains a multilayer perceptron which has a objective
   /// functional associated. 

   virtual void train(void) = 0;

   // Utility methods

   /// This method prints to the screen the class members in a training algorithm.

   virtual void print(void) = 0;

   /// This method saves the training algorithm object to a data file.

   virtual void save(char*) = 0;

   /// This method loads the training algorithm object from a data file.

   virtual void load(char*) = 0;

   /// This method makes the training history of all variables to be reseved or not in memory.

   virtual void setReserveAllTrainingHistory(bool) = 0;

   /// This method saves the training history to a data file.

   virtual void saveTrainingHistory(char*) = 0;
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
