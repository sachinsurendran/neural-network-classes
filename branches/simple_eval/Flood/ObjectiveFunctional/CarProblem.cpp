/****************************************************************************************************************/
/*                                                                                                              */
/*   Flood: An Open Source Neural Networks C++ Library                                                          */
/*   www.cimne.com/flood                                                                                        */
/*                                                                                                              */
/*   C A R   P R O B L E M   C L A S S                                                                          */
/*                                                                                                              */
/*   Roberto Lopez                                                                                              */
/*   International Center for Numerical Methods in Engineering (CIMNE)                                          */
/*   Technical University of Catalonia (UPC)                                                                    */
/*   Barcelona, Spain                                                                                           */
/*   E-mail: rlopez@cimne.upc.es                                                                                */
/*                                                                                                              */
/****************************************************************************************************************/

#include <iostream>     
#include <fstream>     
#include <math.h>     

#include "CarProblem.h"

#include "../Utilities/OrdinaryDifferentialEquations.h"     

namespace Flood
{

// GENERAL CONSTRUCTOR
//
/// General constructor. It creates a car problem objective functional associated to a multilayer perceptron.
/// It also initializes all the problem parameters to their default values:
///
/// <ul>
/// <li> Initial position: 0.
/// <li> Initial velocity: 0.
/// <li> Final position goal: 1.
/// <li> Final velocity goal: 0.
/// <li> Maximum acceleration: 1.
/// <li> Maximum decceleration: -1.
/// <li> Final position error weight: 1.
/// <li> Final velocity error weight: 1.
/// <li> Final time weight: 1.0e-3.
/// <li> Tolerance: 1.0e-12.
/// <li> Initial size: 1000.
/// </ul>
/// 
/// @param newMultilayerPerceptron Pointer to a multilayer perceptron object.

CarProblem::CarProblem(MultilayerPerceptron* newMultilayerPerceptron)
: ObjectiveFunctional(newMultilayerPerceptron)
{
   // Control sentence

   int numberOfInputs = multilayerPerceptron->getNumberOfInputs();
   int numberOfOutputs = multilayerPerceptron->getNumberOfOutputs();

   if(numberOfInputs != 1 || numberOfOutputs != 2)
   {
      std::cerr << std::endl
                << "Flood Error: CarProblem class." << std::endl
                << "CarProblem(MultilayerPerceptron*) constructor." << std::endl
                << "Number of inputs and outputs in multilayer perceptron must be 1 and 2, respectively." 
                << std::endl
                << std::endl;

      exit(1);
   }

   int numberOfIndependentParameters = multilayerPerceptron->getNumberOfIndependentParameters();

   if(numberOfIndependentParameters != 1)
   {
      std::cerr << std::endl
                << "Flood Error: CarProblem class." << std::endl
                << "CarProblem(MultilayerPerceptron*) constructor." << std::endl
                << "Number of independent parameters of multilayer perceptron must be 1." << std::endl
                << std::endl;

      exit(1);
   }
   
   initialPosition = 0.0;
   initialVelocity = 0.0;

   finalPositionGoal = 1.0;
   finalVelocityGoal = 0.0;

   maximumAcceleration = 1.0;
   maximumDeceleration = 1.0;
   
   finalPositionErrorWeight = 1.0;
   finalVelocityErrorWeight = 1.0;

   finalTimeWeight = 1.0e-9;

   tolerance = 1.0e-15;
   initialSize = 1000;

   ordinaryDifferentialEquations.setTolerance(tolerance);
   ordinaryDifferentialEquations.setInitialSize(initialSize);
}


// DEFAULT CONSTRUCTOR

/// Default constructor. It creates a car problem objective functional not associated to any multilayer 
/// perceptron. It also initializes all the problem parameters to their default values:
///
/// <ul>
/// <li> Initial position: 0.
/// <li> Initial velocity: 0.
/// <li> Final position goal: 1.
/// <li> Final velocity goal: 0.
/// <li> Maximum acceleration: 1.
/// <li> Maximum decceleration: -1.
/// <li> Final position error weight: 1.
/// <li> Final velocity error weight: 1.
/// <li> Final time weight: 1.0e-3.
/// <li> Tolerance: 1.0e-12.
/// <li> Initial size: 1000.
/// </ul>

CarProblem::CarProblem(void): ObjectiveFunctional()
{
   initialPosition = 0.0;
   initialVelocity = 0.0;

   finalPositionGoal = 1.0;
   finalVelocityGoal = 0.0;

   maximumAcceleration = 1.0;
   maximumDeceleration = 1.0;
   
   finalPositionErrorWeight = 1.0;
   finalVelocityErrorWeight = 1.0;

   finalTimeWeight = 1.0e-9;

   tolerance = 1.0e-15;
   initialSize = 10000;
}


// DESTRUCTOR

/// Destructor. 

CarProblem::~CarProblem(void)
{

}


// METHODS

// double getInitialPosition(void) method

/// This method returns the position of the car at the initial time. 

double CarProblem::getInitialPosition(void)
{
   return(initialPosition);
}


// double getInitialVelocity(void) method

/// This method returns the velocity of the car at the initial time. 

double CarProblem::getInitialVelocity(void)
{
   return(initialVelocity);
}


// double getFinalPositionGoal(void) method

/// This method returns the desired position of the car at the final time. 

double CarProblem::getFinalPositionGoal(void)
{
   return(finalPositionGoal);
}


// double getFinalVelocityGoal(void) method

/// This method returns the desired velocity of the car at the final time. 

double CarProblem::getFinalVelocityGoal(void)
{
   return(finalVelocityGoal);
}


// double getMaximumAcceleration(void) method

/// This method returns the maximum acceleration of the car, which is limited by the capability of the engine. 

double CarProblem::getMaximumAcceleration(void)
{
   return(maximumAcceleration);
}


// double getMaximumDeceleration(void) method

/// This method returns the maximum deceleration of the car, which is limited by the braking system parameters. 

double CarProblem::getMaximumDeceleration(void)
{
   return(maximumDeceleration);
}


// double getFinalPositionErrorWeight(void) method

/// This method returns the weight of the final position error term in the objective functional.

double CarProblem::getFinalPositionErrorWeight(void)
{
   return(finalPositionErrorWeight);
}


// double getFinalVelocityErrorWeight(void) method

/// This method returns the weight of the final velocity error term in the objective functional.

double CarProblem::getFinalVelocityErrorWeight(void)
{
   return(finalVelocityErrorWeight);
}


// double getFinalTimeWeight(void) method

/// This method returns the weight of the final time term in the objective functional.

double CarProblem::getFinalTimeWeight(void)
{
   return(finalTimeWeight);
}


// double getTolerance(void) method

/// This method returns the tolerance value to be used in the Runge-Kutta-Fehlberg method for evaluating the 
/// final position and velocity errors.

double CarProblem::getTolerance(void)
{
   return(tolerance);
}


// int getInitialSize(void) method

/// This method returns the number of points to be reserved when using the Runge-Kutta-Fehlberg method for 
/// evaluating the final position and velocity errors.

int CarProblem::getInitialSize(void)
{
   return(initialSize);
}


// void setInitialPosition(double) method

/// This method sets a new position of the car at the initial time. 
///
/// @param newInitialPosition Initial position value. 

void CarProblem::setInitialPosition(double newInitialPosition)
{
   initialPosition = newInitialPosition;
}

   
// void setInitialVelocity(double) method

/// This method sets a new velocity of the car at the initial time. 
///
/// @param newInitialVelocity Initial velocity value. 

void CarProblem::setInitialVelocity(double newInitialVelocity)
{
   initialVelocity = newInitialVelocity;
}


// void setFinalPositionGoal(double) method

/// This method sets a new desired position of the car at the final time. 
///
/// @param newFinalPositionGoal Desired final position.

void CarProblem::setFinalPositionGoal(double newFinalPositionGoal)
{
   finalPositionGoal = newFinalPositionGoal;
}

   
// void setFinalVelocityGoal(double) method

/// This method sets a new desired velocity of the car at the final time. 
///
/// @param newFinalVelocityGoal Desired final velocity.

void CarProblem::setFinalVelocityGoal(double newFinalVelocityGoal)
{
   finalVelocityGoal = newFinalVelocityGoal;
}


// void setMaximumAcceleration(double) method

/// This method sets a new value for the maximum acceleration of the car.
/// The maximum acceleration must be a positive number. 
///
/// @param newMaximumAcceleration Maximum acceleration. 

void CarProblem::setMaximumAcceleration(double newMaximumAcceleration)
{
   // Control sentence

   if(newMaximumAcceleration <= 0.0)
   {
      std::cerr << std::endl
                << "Flood Error: CarProblem class."
                << "void setMaximumAcceleration(double) method." << std::endl
                << "The maximum acceleration must be greater than zero." << std::endl
                << std::endl;

      exit(1);
   }

   // Set maximum acceleration

   maximumAcceleration = newMaximumAcceleration;
}


// void setMaximumDeceleration(double) method

/// This method sets a new value for the maximum deceleration of the car.
/// The maximum deceleration must be a negative number. 
///
/// @param newMaximumDeceleration Maximum deceleration. 

void CarProblem::setMaximumDeceleration(double newMaximumDeceleration)
{
   // Control sentence

   if(newMaximumDeceleration >= 0.0)
   {
      std::cerr << std::endl
                << "Flood Error: CarProblem class."
                << "void setMaximumDeceleration(double) method." << std::endl
                << "The maximum deceleration must be less than zero." << std::endl
                << std::endl;

      exit(1);
   }

   // Set maximum acceleration

   maximumDeceleration = newMaximumDeceleration;
}


// void setFinalPositionErrorWeight(double) method

/// This method sets a new weight value for the final position error term in the objective functional expression.
/// The weighs of all terms in the objective functional must be positive. 
///
/// @param newFinalPositionErrorWeight Weight of final position error term.

void CarProblem::setFinalPositionErrorWeight(double newFinalPositionErrorWeight)
{
   // Control sentence

   if(newFinalPositionErrorWeight <= 0.0)
   {
      std::cerr << std::endl
                << "Flood Error: CarProblem class."
                << "void setFinalPositionErrorWeight(double) method." << std::endl
                << "The weight of the final position error term must be greater than zero." << std::endl
                << std::endl;

      exit(1);
   }

   // Set final position error term weight

   finalPositionErrorWeight = newFinalPositionErrorWeight;
}


// void setFinalVelocityErrorWeight(double) method

/// This method sets a new weight value for the final velocity error term in the objective functional expression.
/// The weighs of all terms in the objective functional must be positive. 
///
/// @param newFinalVelocityErrorWeight Weight of final velocity error term.

void CarProblem::setFinalVelocityErrorWeight(double newFinalVelocityErrorWeight)
{
   // Control sentence

   if(newFinalVelocityErrorWeight <= 0.0)
   {
      std::cerr << std::endl
                << "Flood Error: CarProblem class."
                << "void setFinalVelocityErrorWeight(double) method." << std::endl
                << "The weight of the final velocity error term must be greater than zero." << std::endl
                << std::endl;

      exit(1);
   }

   // Set final velocity error term weight

   finalVelocityErrorWeight = newFinalVelocityErrorWeight;
}


// void setFinalTimeWeight(double) method

/// This method sets a new weight value for the final time term in the objective functional, which is to be 
/// minimized.
/// The weighs of all terms in the objective functional must be positive. 
///
/// @param newFinalTimeWeight Weight of final time term.

void CarProblem::setFinalTimeWeight(double newFinalTimeWeight)
{
   // Control sentence

   if(newFinalTimeWeight <= 0.0)
   {
      std::cerr << std::endl
                << "Flood Error: CarProblem class." << std::endl
                << "void setFinalTimeWeight(double) method." << std::endl
                << "The weight of the final time term must be greater than zero." << std::endl
                << std::endl;

      exit(1);
   }

   // Set final time term weight

   finalTimeWeight = newFinalTimeWeight;
}


// void setTolerance(double) method

/// This method sets a new tolerance value to be used in the Runge-Kutta-Fehlberg method for evaluating the final 
/// position and velocity errors.
/// The tolerance of integration must be a small value greater than zero. 
///
/// @param newTolerance Tolerance in Runge-Kutta-Fehlberg method.

void CarProblem::setTolerance(double newTolerance)
{
   // Control sentence

   if(newTolerance <= 0.0)
   {
      std::cerr << std::endl
                << "Flood Error: CarProblem class." << std::endl
                << "void setTolerance(double) method." << std::endl
                << "The tolerance of integration must be greater than zero." << std::endl
                << std::endl;

      exit(1);
   }

   // Set tolerance 

   tolerance = newTolerance;

   ordinaryDifferentialEquations.setTolerance(tolerance);
}


// void setInitialSize(int) method

/// This method sets a new number of points to be reserved when using the Runge-Kutta-Fehlberg method for 
/// evaluating the final position and velocity errors. 
/// The initial size must be a big value greater than zero, otherwise resizing will be necessary. 
///
/// @param newInitialSize Number of points to reserve in Runge-Kutta-Fehlberg method.

void CarProblem::setInitialSize(int newInitialSize)
{
   // Control sentence

   if(newInitialSize <= 0)
   {
      std::cerr << std::endl
                << "Flood Error: CarProblem class." << std::endl
                << "void setInitialsize(int) method." << std::endl
                << "The initial size must be greater than zero." << std::endl
                << std::endl;

      exit(1);
   }

   // Set initial size

   initialSize = newInitialSize;

   ordinaryDifferentialEquations.setInitialSize(initialSize);
}


// double calculatePositionDot(double, double, double) method

/// This method depicts the state equation for the position in the car problem. It returns the derivative of the 
/// position with respect to the time as a function of time, position and velocity.
///
/// @param time Actual time.
/// @param position Actual position.
/// @param velocity Actual velocity.

double CarProblem::calculatePositionDot(double time, double position, double velocity)
{
   double positionDot = velocity;

   return(positionDot);
}


// double calculateVelocityDot(double, double, double) method

/// This method depicts the state equation for the velocity in the car problem. It returns the derivative of the 
/// velocity with respect to the time as a function of time, position and velocity.
///
/// @param time Actual time.
/// @param position Actual position.
/// @param velocity Actual velocity.

double CarProblem::calculateVelocityDot(double time, double position, double velocity)
{
   Vector<double> input(1);
   input[0] = time;

   Vector<double> output = multilayerPerceptron->calculateOutput(input);

   double acceleration = output[0];

   if(acceleration < 0.0)
   {
      acceleration = 0.0;                
   }
   else if(acceleration > maximumAcceleration)
   {
      acceleration = maximumAcceleration;                
   }

   double deceleration = output[1];

   if(deceleration < -1.0*maximumDeceleration)
   {
      deceleration = -1.0*maximumDeceleration;                
   }
   else if(deceleration > 0.0)
   {
      deceleration = 0.0;                
   }

   double velocityDot = acceleration + deceleration;

   return(velocityDot);
}


// double calculateEvaluation(void) method

/// This method returns the evaluation value of a multilayer perceptron for the car problem. 

double CarProblem::calculateEvaluation(void)
{
    // Control sentence 

   if(multilayerPerceptron == NULL)
   {
      std::cerr << std::endl
                << "Flood Error: CarProblem class." << std::endl
                << "double calculateEvaluation(void) method." << std::endl
                << "Pointer to multilayer perceptron object cannot be NULL." << std::endl
                << std::endl;

      exit(1);
   }

   double evaluation = 0.0;
   
   // Increment number of evaluations

   numberOfEvaluations++;

   // Evaluate objective functional

   int numberOfIndependentParameters 
   = multilayerPerceptron->getNumberOfIndependentParameters();

   if(numberOfIndependentParameters != 1)
   {
      std::cout << std::endl
                << "Flood Error CarProblem class." << std::endl
                << "calculateEvaluation(void) method." << std::endl
                << "Number of independent parameters in multilayer perceptron must be one." << std::endl
                << std::endl;

      exit(1);                         
   }

   Vector<double> independentParameters = multilayerPerceptron->getIndependentParameters();

   double initialTime = 0.0;
   double finalTime = independentParameters[0];

   if(finalTime < 0.0)
   {
      finalTime = 0.0;             
   }

   Vector<double> time;
   Vector<double> position;
   Vector<double> velocity;

   int numberOfPoints = 
   ordinaryDifferentialEquations.calculateRungeKuttaFehlbergIntegral(*this,
   time, position, velocity, 
   &CarProblem::calculatePositionDot,
   &CarProblem::calculateVelocityDot,
   initialTime, finalTime, initialPosition, initialVelocity);

   // Obtain final position and velocity errors

   double finalPositionError = fabs(finalPositionGoal - position[numberOfPoints-1]);
   double finalVelocityError = fabs(finalVelocityGoal - velocity[numberOfPoints-1]);

   // Calculate evaluation

   evaluation
   = finalTimeWeight*finalTime
   + finalPositionErrorWeight*pow(finalPositionError, 2)
   + finalVelocityErrorWeight*pow(finalVelocityError, 2);

   return(evaluation);
}


// void saveResults(char*) method

/// This method saves the trajectory and the control signal for the car problem to a data file. 
/// The column format is as follows:
///
/// <ol>
/// <li> Time.
/// <li> Position.
/// <li> Velocity.
/// <li> Acceleration.
/// <li> Deceleration.
/// </ol>
///
/// @param filename Filename.

void CarProblem::saveResults(char* filename)
{
   std::fstream file; 

   file.open(filename, std::ios::out);
   
    if(!file.is_open())
    {
       std::cout << std::endl
                 << "Cannot open car problem results data file."  
                 << std::endl << std::endl;

       exit(1);
    }
    else
    {
       if(display)
      {
          std::cout << std::endl
                    << "Saving car problem results to data file..." << std::endl;
      }
    }

   int numberOfPoints = 1001;

   ordinaryDifferentialEquations.setNumberOfPoints(numberOfPoints);

   // Solve state equations 

   Vector<double> time(numberOfPoints);
   Vector<double> position(numberOfPoints);
   Vector<double> velocity(numberOfPoints);

   Vector<double> independentParameters 
    = multilayerPerceptron->getIndependentParameters();

   double initialTime = 0.0;
   double finalTime = independentParameters[0];

   if(finalTime < 0.0)
   {
      finalTime = 0.0;             
   }
          
   ordinaryDifferentialEquations.getRungeKuttaIntegral(*this,
   time, position, velocity, 
   &CarProblem::calculatePositionDot,
   &CarProblem::calculateVelocityDot,
   initialTime, finalTime, initialPosition, initialVelocity);

   double finalPosition = position[numberOfPoints-1];
   double finalPositionError = finalPositionGoal - finalPosition;

   double finalVelocity = velocity[numberOfPoints-1];
   double finalVelocityError = finalVelocityGoal - finalVelocity;

   // Obtain control 

   Vector<double> input(1);
   Vector<double> output(2);

   double acceleration = 0.0;
   double deceleration = 0.0;

    // Write header 

   file << "% Flood Neural Network. Car problem results data file." << std::endl
        << "% Final time: " << finalTime << std::endl
        << "% Final position error: " << finalPositionError << std::endl
        << "% Final velocity error: " << finalVelocityError << std::endl
        << "% Column data:" << std::endl
        << "%  1 - Time" << std::endl
        << "%  2 - Position" << std::endl
        << "%  3 - Velocity" << std::endl
        << "%  4 - Acceleration" << std::endl
        << "%  5 - Deceleration" << std::endl
        << std::endl;

   for(int i = 0; i < numberOfPoints; i++)
   {
      input[0] = time[i];
      output = multilayerPerceptron->calculateOutput(input);
      acceleration = output[0];
      deceleration = output[1];
      
      // Acceleration

      if(acceleration < 0.0)
      {
         acceleration = 0.0;
      }
      else if(acceleration > maximumAcceleration)
      {
         acceleration = maximumAcceleration;         
      }
      
      // Deceleration

      if(deceleration < -1.0*maximumDeceleration)
      {
         deceleration = -1.0*maximumDeceleration;
      }
      else if(deceleration > 0.0)
      {
         deceleration = 0.0;         
      }

      // Write time, trajectory and control signal to file   

      file << time[i] << " " 
           << position[i] << " "
           << velocity[i] << " "
           << acceleration << " " 
           << deceleration << std::endl;   
   }   

   file.close();
}


// void print(void) method

/// This method prints to the screen useful information about the car problem objective functional of a 
/// multilayer perceptron. 
///
/// <ul>
/// <li> Final position error.
/// <li> Final velocity error.
/// <li> Final time.
/// </ul>

void CarProblem::print()
{
   // Evaluate performance

   // Solve state equations 

   Vector<double> time;
   Vector<double> position;
   Vector<double> velocity;

   Vector<double> independentParameters 
   = multilayerPerceptron->getIndependentParameters();

   double initialTime = 0.0;
   double finalTime = independentParameters[0];

   if(finalTime < 0.0)
   {
      finalTime = 0.0;             
   }

   int numberOfPoints = 
   ordinaryDifferentialEquations.calculateRungeKuttaFehlbergIntegral(*this,
   time, position, velocity, 
   &CarProblem::calculatePositionDot,
   &CarProblem::calculateVelocityDot,
   initialTime, finalTime, initialPosition, initialVelocity);

   // Obtain final position and velocity errors

   double finalPositionError 
   = fabs(finalPositionGoal - position[numberOfPoints-1]);

   double finalVelocityError 
   = fabs(finalVelocityGoal - velocity[numberOfPoints-1]);

   std::cout << "Final position error: " << finalPositionError << std::endl      
           << "Final velocity error: " << finalVelocityError << std::endl
           << "Final time: " << finalTime << std::endl;
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
