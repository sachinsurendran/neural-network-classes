/****************************************************************************************************************/
/*                                                                                                              */
/*   Flood: An Open Source Neural Networks C++ Library                                                          */
/*   www.cimne.com/flood                                                                                        */
/*                                                                                                              */
/*   T R A I N I N G   A L G O R I T H M   C L A S S                                                            */
/*                                                                                                              */
/*   Roberto Lopez                                                                                              */
/*   International Center for Numerical Methods in Engineering (CIMNE)                                          */
/*   Technical University of Catalonia (UPC)                                                                    */
/*   Barcelona, Spain                                                                                           */
/*   E-mail: rlopez@cimne.upc.edu                                                                               */
/*                                                                                                              */
/****************************************************************************************************************/

#include<iostream>
#include<fstream>

#include "TrainingAlgorithm.h"

namespace Flood
{

// GENERAL CONSTRUCTOR
//
/// General constructor. It creates a training algorithm object associated to an objective functional object.
///
/// @param newObjectiveFunctional Pointer to an objective functional object.

TrainingAlgorithm::TrainingAlgorithm(ObjectiveFunctional* newObjectiveFunctional)
{
   objectiveFunctional = newObjectiveFunctional;

   reserveElapsedTimeHistory = false;
   reserveFreeParametersHistory = false;
   reserveFreeParametersNormHistory = false;
   reserveEvaluationHistory = false;
   reserveGradientHistory = false;
   reserveGradientNormHistory = false;

   warningFreeParametersNorm = 10000.0;
   warningGradientNorm = 1000.0;
   
   display = true;
}


// DEFAULT CONSTRUCTOR

/// Default constructor. It creates a training algorithm object not associated to any objective functional object.  

TrainingAlgorithm::TrainingAlgorithm(void)
{
   objectiveFunctional = NULL;

   reserveFreeParametersHistory = false;
   reserveFreeParametersNormHistory = false;
   reserveEvaluationHistory = false;
   reserveGradientHistory = false;
   reserveGradientNormHistory = false;
   reserveElapsedTimeHistory = false;

   warningFreeParametersNorm = 1000.0;
   warningGradientNorm = 1000.0;

   display = true;   
}


// DESTRUCTOR 

/// Destructor.

TrainingAlgorithm::~TrainingAlgorithm(void)
{ 

}


// METHODS

// ObjectiveFunctional* getObjectiveFunctional(void) method

/// This method returns a pointer to the objective functional object to which the training algorithm is 
/// associated.

ObjectiveFunctional* TrainingAlgorithm::getObjectiveFunctional(void)
{
   return(objectiveFunctional);
}


// double getEvaluationGoal(void) method

/// This method returns the goal value for the evaluation. 
/// This is used as a stopping criterium when training a neural network.
///  
/// @see train(void).

double TrainingAlgorithm::getEvaluationGoal(void)
{
   return(evaluationGoal);
}


// double getGradientNormGoal(void) method

/// This method returns the goal value for the norm of the objective function gradient.
/// This is used as a stopping criterium when training a neural network.
///  
/// @see train(void).

double TrainingAlgorithm::getGradientNormGoal(void)
{
   return(gradientNormGoal);
}


// double getMaximumTime(void) method

/// This method returns the maximum training time.  
///  
/// @see train(void).

double TrainingAlgorithm::getMaximumTime(void)
{
   return(maximumTime);
}


// double getMinimumImprovement(void) method

/// This method returns the minimum evaluation improvement during training.  
///  
/// @see train(void).

double TrainingAlgorithm::getMinimumImprovement(void)
{
   return(minimumImprovement);
}


// int getMaximumNumberOfEpochs(void) method

/// This method returns the maximum number of epochs for training.
///  
/// @see train(void).

int TrainingAlgorithm::getMaximumNumberOfEpochs(void)
{
   return(maximumNumberOfEpochs);
}


// double getWarningFreeParametersNorm(void) method

/// This method returns the minimum value for the norm of the free parameters vector at wich a warning message is 
/// written to the screen. 

double TrainingAlgorithm::getWarningFreeParametersNorm(void)
{
   return(warningFreeParametersNorm);       
}


// double getWarningGradientNorm(void) method

/// This method returns the minimum value for the norm of the gradient vector at wich a warning message is written
/// to the screen. 

double TrainingAlgorithm::getWarningGradientNorm(void)
{
   return(warningGradientNorm);       
}


// bool getDisplay(void) method

/// This method returns true if messages from this class can be displayed on the screen, or false if messages from
/// this class can't be displayed on the screen.

bool TrainingAlgorithm::getDisplay(void)
{
   return(display);
}


// int getDisplayPeriod(void) method

/// This method returns the number of epochs between the training showing progress. 
///  
/// @see train(void).

int TrainingAlgorithm::getDisplayPeriod(void)
{
   return(displayPeriod);
}


// bool getReserveFreeParametersHistory(void) method

/// This method returns true if the free parameters history matrix is to be reserved, and false otherwise.

bool TrainingAlgorithm::getReserveFreeParametersHistory(void)
{
   return(reserveFreeParametersHistory);     
}


// bool getReserveFreeParametersNormHistory(void) method 

/// This method returns true if the free parameters norm history vector is to be reserved, and false otherwise.

bool TrainingAlgorithm::getReserveFreeParametersNormHistory(void)
{
   return(reserveFreeParametersNormHistory);     
}


// bool getReserveEvaluationHistory(void) method

/// This method returns true if the evaluation history vector is to be reserved, and false otherwise.

bool TrainingAlgorithm::getReserveEvaluationHistory(void)
{
   return(reserveEvaluationHistory);     
}


// bool getReserveGradientHistory(void) method

/// This method returns true if the gradient history matrix is to be reserved, and false otherwise.

bool TrainingAlgorithm::getReserveGradientHistory(void)
{
   return(reserveGradientHistory);     
}


// bool getReserveGradientNormHistory(void) method

/// This method returns true if the gradient norm history vector is to be reserved, and false otherwise.

bool TrainingAlgorithm::getReserveGradientNormHistory(void)
{
   return(reserveGradientNormHistory);     
}


// bool getReserveElapsedTimeHistory(void) method

/// This method returns true if the elapsed time history vector is to be reserved, and false otherwise.

bool TrainingAlgorithm::getReserveElapsedTimeHistory(void)
{
   return(reserveElapsedTimeHistory);     
}


// Matrix<double> getFreeParametersHistory(void) method

/// This method returns a matrix containing the free parameters history over the training epochs.

Matrix<double> TrainingAlgorithm::getFreeParametersHistory(void)
{
   return(freeParametersHistory);     
}


// Vector<double> getFreeParametersNormHistory(void) method

/// This method returns a vector containing the free parameters norm history over the training epochs.

Vector<double> TrainingAlgorithm::getFreeParametersNormHistory(void)
{
   return(freeParametersNormHistory);     
}


// Vector<double> getEvaluationHistory(void) method

/// This method returns a vector containing the evaluations history over the training epochs.

Vector<double> TrainingAlgorithm::getEvaluationHistory(void)
{
   return(evaluationHistory);     
}


// Matrix<double> getGradientHistory(void) method

/// This method returns a matrix containing the gradient history over the training epochs.

Matrix<double> TrainingAlgorithm::getGradientHistory(void)
{
   return(gradientHistory);     
}


// Vector<double> getGradientNormHistory(void) method

/// This method returns a vector containing the gradient norm history over the training epochs.

Vector<double> TrainingAlgorithm::getGradientNormHistory(void)
{
   return(gradientNormHistory);     
}


// Vector<double> getElapsedTimeHistory(void) method

/// This method returns a matrix containing the elapsed time history over the training epochs.

Vector<double> TrainingAlgorithm::getElapsedTimeHistory(void) 
{
   return(elapsedTimeHistory);     
}


// void setObjectiveFunctional(ObjectiveFunctional*) method

/// This method sets a pointer to an objective functional object to be associated to the training algorithm.
///
/// @param newObjectiveFunctional Pointer to an objective functional object.

void TrainingAlgorithm::setObjectiveFunctional(ObjectiveFunctional* newObjectiveFunctional)
{
   objectiveFunctional = newObjectiveFunctional;
}



// void setEvaluationGoal(double) method

/// This method sets a new goal value for the evaluation. 
/// This is used as a stopping criterium when training a neural network.
///
/// @param newEvaluationGoal Goal value for the evaluation.
/// 
/// @see train(void).

void TrainingAlgorithm::setEvaluationGoal(double newEvaluationGoal)
{
   evaluationGoal = newEvaluationGoal;
}


// void setGradientNormGoal(double) method

/// This method sets a new the goal value for the norm of the objective function gradient. 
/// This is used as a stopping criterium when training a neural network.
///
/// @param newGradientNormGoal Goal value for the norm of the objective function gradient.
/// 
/// @see train(void).

void TrainingAlgorithm::setGradientNormGoal(double newGradientNormGoal)
{
   // Control sentence (if debug)

   #ifndef NDEBUG 

   if(newGradientNormGoal < 0.0)
   {
      std::cerr << std::endl
                << "Flood Error: TrainingAlgorithm class." << std::endl
                << "void setGradientNormGoal(double) method." << std::endl
                << "Gradient norm goal must be equal or greater than 0." << std::endl 
                << std::endl;

      exit(1);
   }

   #endif

   // Set gradient norm goal

   gradientNormGoal = newGradientNormGoal;
}


// void setMaximumTime(double) method

/// This method sets a new maximum training time.  
///
/// @param newMaximumTime Maximum training time.
/// 
/// @see train(void).

void TrainingAlgorithm::setMaximumTime(double newMaximumTime)
{
   // Control sentence (if debug)

   #ifndef NDEBUG 

   if(newMaximumTime <= 0.0)
   {
      std::cerr << std::endl 
                << "Flood Error: TrainingAlgorithm class." << std::endl
                << "void setMaximumTime(double) method." << std::endl
                << "Maximum time must be greater than 0." << std::endl
                << std::endl;

      exit(1);
   }
   
   #endif

   // Set maximum time

   maximumTime = newMaximumTime;
}



// void setMinimumImprovement(double) method

/// This method sets a new minimum evaluation improvement during training.  
///
/// @param newMinimumImprovement Minimum evaluation improvement.
/// 
/// @see train(void).

void TrainingAlgorithm::setMinimumImprovement(double newMinimumImprovement)
{
   // Control sentence (if debug)

   #ifndef NDEBUG 

   if(newMinimumImprovement < 0.0)
   {
      std::cerr << std::endl 
                << "Flood Error: TrainingAlgorithm class." << std::endl
                << "void setMinimumImprovement(double) method." << std::endl
                << "Minimum improvement must be equal or greater than 0." << std::endl
                << std::endl;

      exit(1);
   }

   #endif

   // Set minimum improvement

   minimumImprovement = newMinimumImprovement;
}


// void setMaximumNumberOfEpochs(int) method

/// This method sets a maximum number of epochs for training.
///
/// @param newMaximumNumberOfEpochs Maximum number of epochs for training.
/// 
/// @see train(void).

void TrainingAlgorithm::setMaximumNumberOfEpochs(int newMaximumNumberOfEpochs)
{
   maximumNumberOfEpochs = newMaximumNumberOfEpochs;
}


// void setWarningGradientNorm(double) methodç

/// This method sets a new value for the gradient vector norm at which 
/// a warning message is written to the screen. 
///
/// @param newWarningGradientNorm Warning norm of gradient vector value. 

void TrainingAlgorithm::setWarningGradientNorm(double newWarningGradientNorm)
{
   // Control sentence (if debug)

   #ifndef NDEBUG 

   if(newWarningGradientNorm <= 0.0)
   {
      std::cerr << std::endl 
                << "Flood Error: TrainingAlgorithm class." << std::endl
                << "void setWarningGradientNorm(double) method." << std::endl
                << "Warning gradient norm must be equal or greater than 0." << std::endl
                << std::endl;

      exit(1);
   }

   #endif

   // Set warning gradient norm

   warningGradientNorm = newWarningGradientNorm;     
}


// void setWarningFreeParametersNorm(double) methodç

/// This method sets a new value for the free parameters vector norm at which a warning message is written to the 
/// screen. 
///
/// @param newWarningFreeParametersNorm Warning norm of free parameters vector value. 

void TrainingAlgorithm::setWarningFreeParametersNorm(double newWarningFreeParametersNorm)
{
   // Control sentence (if debug)

   #ifndef NDEBUG 

   if(newWarningFreeParametersNorm <= 0.0)
   {
      std::cerr << std::endl 
                << "Flood Error: TrainingAlgorithm class." << std::endl
                << "void setWarningFreeParametersNorm(double) method." << std::endl
                << "Warning free parameters norm must be equal or greater than 0." << std::endl
                << std::endl;

      exit(1);
   }

   #endif

   // Set warning free parameters norm

   warningFreeParametersNorm = newWarningFreeParametersNorm;     
}


// void setDisplay(bool) method

/// This method sets a new display value. 
/// If it is set to true messages from this class are to be displayed on the screen;
/// if it is set to false messages from this class are not to be displayed on the screen.
///
/// @param newDisplay Display value.

void TrainingAlgorithm::setDisplay(bool newDisplay)
{
   display = newDisplay;
}


// void setDisplayPeriod(int) method

/// This method sets a new number of epochs between the training showing progress. 
///
/// @param newDisplayPeriod
/// Number of epochs between the training showing progress. 
/// 
/// @see train(void).

void TrainingAlgorithm::setDisplayPeriod(int newDisplayPeriod)
{
   displayPeriod = newDisplayPeriod;
}


// void setReserveFreeParametersHistory(bool) method

/// This method makes the free parameters history matrix to be reseved or not in memory.
///
/// @param newReserveFreeParametersHistory True if the free parameters history matrix is to be reserved, false 
/// otherwise.

void TrainingAlgorithm::setReserveFreeParametersHistory(bool newReserveFreeParametersHistory)
{
   reserveFreeParametersHistory = newReserveFreeParametersHistory;     
}


// void setReserveFreeParametersNormHistory(bool) method

/// This method makes the free parameters norm history vector to be reseved or not in memory.
///
/// @param newReserveFreeParametersNormHistory True if the free parameters norm history vector is to be reserved, 
/// false otherwise.

void TrainingAlgorithm::setReserveFreeParametersNormHistory(bool newReserveFreeParametersNormHistory)
{
   reserveFreeParametersNormHistory = newReserveFreeParametersNormHistory;     
}


// void setReserveEvaluationHistory(bool) method

/// This method makes the evaluation history vector to be reseved or not in memory.
///
/// @param newReserveEvaluationHistory True if the evaluation history vector is to be reserved, false otherwise.

void TrainingAlgorithm::setReserveEvaluationHistory(bool newReserveEvaluationHistory)
{
   reserveEvaluationHistory = newReserveEvaluationHistory;     
}


// void setReserveGradientHistory(bool) method

/// This method makes the gradient history matrix to be reseved or not in memory.
///
/// @param newReserveGradientHistory True if the gradient history matrix is to be reserved, false otherwise.

void TrainingAlgorithm::setReserveGradientHistory(bool newReserveGradientHistory)
{
   reserveGradientHistory = newReserveGradientHistory;    
}


// void setReserveGradientNormHistory(bool) method

/// This method makes the gradient norm history matrix to be reseved or not in memory.
///
/// @param newReserveGradientNormHistory True if the gradient norm history matrix is to be reserved, false 
/// otherwise.

void TrainingAlgorithm::setReserveGradientNormHistory(bool newReserveGradientNormHistory)
{
   reserveGradientNormHistory = newReserveGradientNormHistory;     
}


// void setReserveElapsedTimeHistory(bool) method

/// This method makes the elapsed time over the epochs to be reseved or not in memory.
///
/// @param newReserveElapsedTimeHistory True if the elapsed time history vector is to be reserved, false 
/// otherwise.

void TrainingAlgorithm::setReserveElapsedTimeHistory(bool newReserveElapsedTimeHistory)
{
   reserveElapsedTimeHistory = newReserveElapsedTimeHistory;     
}


// void setFreeParametersHistory(Matrix<double>) method

/// This method sets a new matrix containing the free parameters history over the training epochs.
/// Each row in the matrix contains the free parameters vector of one single epoch. 
///
/// @param newFreeParametersHistory Free parameters history matrix. 

void TrainingAlgorithm::setFreeParametersHistory(Matrix<double> newFreeParametersHistory)
{
   freeParametersHistory = newFreeParametersHistory;     
}


// void setFreeParametersNormHistory(Vector<double>) method

/// This method sets a new matrix containing the free parameters norm history over the training epochs.
/// Each element in the vector contains the free parameters norm of one single epoch. 
///
/// @param newFreeParametersNormHistory Free parameters norm history vector. 

void TrainingAlgorithm::setFreeParametersNormHistory(Vector<double> newFreeParametersNormHistory)
{
   freeParametersNormHistory = newFreeParametersNormHistory;     
}


// void setEvaluationHistory(Vector<double>) method

/// This method sets a new vector containing the evaluation history over the training epochs.
/// Each row in the matrix contains the free parameters vector of one single epoch. 
///
/// @param newEvaluationHistory Evaluation history vector. 

void TrainingAlgorithm::setEvaluationHistory(Vector<double> newEvaluationHistory)
{
   evaluationHistory = newEvaluationHistory;     
}


// void setGradientHistory(Matrix<double>) method

/// This method sets a new gradient history matrix over the training epochs. 
/// The number of rows must be equal to the training size.
/// The number of columns must be equal to the number of free parameters. 
///
/// @param newGradientHistory Gradient history matrix.

void TrainingAlgorithm::setGradientHistory(Matrix<double> newGradientHistory)
{
   gradientHistory = newGradientHistory;     
}


// void setGradientNormHistory(Vector<double>) method

/// This method sets a new gradient norm history vector. 
/// The elements in the vector are the gradient norm values over the training epochs.
///
/// @param newGradientNormHistory Gradient norm history vector. 

void TrainingAlgorithm::setGradientNormHistory(Vector<double> newGradientNormHistory)
{
   gradientNormHistory = newGradientNormHistory;     
}


// void setElapsedTimeHistory(Vector<double>) method

/// This method sets a new elapsed time history vector. 
/// The elements in the vector are the measured times over the training epochs.
///
/// @param newElapsedTimeHistory Elapsed time history vector. 

void TrainingAlgorithm::setElapsedTimeHistory(Vector<double> newElapsedTimeHistory)
{
   elapsedTimeHistory = newElapsedTimeHistory;     
}

}


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
