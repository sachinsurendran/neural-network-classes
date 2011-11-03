/****************************************************************************************************************/
/*                                                                                                              */
/*   Flood: An Open Source Neural Networks C++ Library                                                          */
/*   www.cimne.com/flood                                                                                        */
/*                                                                                                              */
/*   A I R C R A F T   L A N D I N G   P R O B L E M   C L A S S                                                */
/*                                                                                                              */
/*   Roberto Lopez and Kevin Lau                                                                                */ 
/*   International Center for Numerical Methods in Engineering (CIMNE)                                          */
/*   Technical University of Catalonia (UPC)                                                                    */
/*   Barcelona, Spain                                                                                           */
/*   E-mail: rlopez@cimne.upc.es, kevin.lau@imperial.ac.uk                                                      */ 
/*                                                                                                              */
/****************************************************************************************************************/

#include <iostream>     
#include <fstream>     
#include <math.h>     

#include "AircraftLandingProblem.h"     

namespace Flood
{

// GENERAL CONSTRUCTOR

/// General constructor. It creates an aircraft landing problem objective functional associated to a multilayer 
/// perceptron. 
/// It also initializes all the problem parameters to their default values:
///
/// <ul>
/// <li> Short period gain = -0.95 (s-1)
/// <li> Path time constant = 2.5 (s)
/// <li> Short period resonant frequency = 1.0 (rad/s)
/// <li> Short period damping factor = 0.5
/// <li> Velocity = 256.0 (ft/s)
/// <li> Landing time = 20.0 (s)
/// <li> Initial pitch angle rate = 0.0 (rad/s)
/// <li> Initial pitch angle = 0.0781 (rad)
/// <li> Initial altitude = velocity*landingTime*tan(-flightPathAngle) // (m)
/// <li> Initial altitude rate = -20.0/3.28084 (m/s)
/// </ul>
///
/// @param newMultilayerPerceptron Pointer to a multilayer perceptron object. This neural network will represent 
/// the control variable for this system. 

AircraftLandingProblem::AircraftLandingProblem(MultilayerPerceptron* newMultilayerPerceptron)       
: ObjectiveFunctional(newMultilayerPerceptron)
{
   // Control sentence

   int numberOfInputs = multilayerPerceptron->getNumberOfInputs();
   int numberOfOutputs = multilayerPerceptron->getNumberOfOutputs();

   if(numberOfInputs != 1 || numberOfOutputs != 1)
   {
      std::cerr << std::endl
                << "Flood Error: AircraftLandingProblem class." << std::endl
                << "AircraftLandingProblem(MultilayerPerceptron*) constructor." << std::endl
                << "Number of inputs and outputs in multilayer perceptron must be 1." << std::endl
                << std::endl;

      exit(1);
   }

   pi = 4.0*atan(1.0);                       // pi = 3.1416

   // Aircraft properties

   shortPeriodGain = -0.95;                  // (s-1)
   pathTimeConstant = 2.5;                   // (s)
   shortPeriodResonantFrequency = 1.0;       // (rad/s)
   shortPeriodDampingFactor = 0.5;           // no units

   // Problem parameters

   velocity = 256.0;                         // (ft/s)
   landingTime = 20.0;                       // (s)

   // Initial values

   initialAltitude = 100.0;                  // (ft)
   initialAltitudeRate = -10.0*1.0;          // (ft/s)
   initialPitchAngle = (1.0/velocity)*initialAltitudeRate; 
   initialPitchAngleRate = 0.0;              // (rad/s)

   // Weight factors for the state variables in all the trajectory

   pitchAngleRateWeight = 1.0;
   pitchAngleWeight = 1.0;
   altitudeRateWeight = 1.0;
   altitudeWeight = 1.0;

   // Weight factors for the states variables at the landing time

   landingPitchAngleRateWeight = 1.0;
   landingPitchAngleWeight = 1.0;
   landingAltitudeRateWeight = 1.0;
   landingAltitudeWeight = 1.0;

   // Weight factor for the control variable

   elevatorWeight = 1.0;

   // Integration stuff

   tolerance = 1.0e-9;
   initialSize = 1000;

   ordinaryDifferentialEquations.setTolerance(tolerance);
   ordinaryDifferentialEquations.setInitialSize(initialSize);
}
 
 
// DEFAULT CONSTRUCTOR

/// Default constructor. It creates an aircraft landing problem objective functional not associated to any 
/// multilayer perceptron. 
/// It also initializes all the problem parameters to their default values:
///
/// <ul>
/// <li> Short period gain = -0.95 (s-1)
/// <li> Path time constant = 2.5 (s)
/// <li> Short period resonant frequency = 1.0 (rad/s)
/// <li> Short period damping factor = 0.5
/// <li> Velocity = 256.0 (ft/s)
/// <li> Landing time = 20.0 (s)
/// <li> Initial pitch angle rate = 0.0 (rad/s)
/// <li> Initial pitch angle = 0.0781 (rad)
/// <li> Initial altitude = velocity*landingTime*tan(-flightPathAngle) // (m)
/// <li> Initial altitude rate = -20.0/3.28084 (m/s)
/// </ul>

AircraftLandingProblem::AircraftLandingProblem(void) : ObjectiveFunctional()
{
   pi = 4.0*atan(1.0);                       // pi = 3.1415
    
    shortPeriodGain = -0.95;                  // (s-1)
   pathTimeConstant = 2.5;                   // (s)
   shortPeriodResonantFrequency = 1.0;       // (rad/s)
   shortPeriodDampingFactor = 0.5;           // no units

   velocity = 256.0;                         // (ft/s)
   landingTime = 20.0;                       // (s)

   initialPitchAngleRate = 0.0;              // (rad/s)
   initialAltitude = 100.0;                  // (ft)
   initialAltitudeRate = -20.0*1.0;              // (ft/s)
    initialPitchAngle = (1.0/velocity)*initialAltitudeRate; 

   tolerance = 1.0e-12;
   initialSize = 1000;
   
   pitchAngleRateWeight = 1.0;
   pitchAngleWeight = 1.0;
   altitudeRateWeight = 1.0;
   altitudeWeight = 1.0;

   landingPitchAngleRateWeight =1.0;
   landingPitchAngleWeight =1.0;
   landingAltitudeRateWeight = 1.0;
   landingAltitudeWeight = 1.0;
   
   elevatorWeight = 1.0;
}


// DESTRUCTOR

/// Destructor. 

AircraftLandingProblem::~AircraftLandingProblem(void) 
{

}


// METHODS

// double getShortPeriodGain(void) method

/// This method returns the short period gain value, which is an aircraf property.  

double AircraftLandingProblem::getShortPeriodGain(void)
{
   return(shortPeriodGain);
}


// double getPathTimeConstant(void) method

/// This method returns the path time constant value, which is an aircraf property.  

double AircraftLandingProblem::getPathTimeConstant(void)
{
   return(pathTimeConstant);
}


// double getShortPeriodResonantFrequency(void) method

/// This method returns the short period resonant frequency, which is an aircraf property.  

double AircraftLandingProblem::getShortPeriodResonantFrequency(void)
{
   return(shortPeriodResonantFrequency);
}


// double getShortPeriodDampingFactor(void) method

/// This method returns the short period damping factor, which is an aircraf property.  

double AircraftLandingProblem::getShortPeriodDampingFactor(void)
{
   return(shortPeriodDampingFactor);
}


// double getVelocity(void) method

/// This method returns the velocity of the aircraft, which is assumed to be constant. 

double AircraftLandingProblem::getVelocity(void)
{
   return(velocity);
}


// double getInitialPitchAngle(void) method

/// This method returns the initial value for the pitch angle (state variable).

double AircraftLandingProblem::getInitialPitchAngle(void)
{
   return(initialPitchAngle);
}


// double getInitialPitchAngleRate(void) method

/// This method returns the initial value for the pitch angle rate (state variable).

double AircraftLandingProblem::getInitialPitchAngleRate(void)
{
   return(initialPitchAngleRate);
}


// double getInitialAltitude(void) method

/// This method returns the initial value for the altitude (state variable).

double AircraftLandingProblem::getInitialAltitude(void)
{
   return(initialAltitude);
}


// double getInitialAltitudeRate(void) method

/// This method returns the initial value for the altitude rate (state variable).

double AircraftLandingProblem::getInitialAltitudeRate(void)
{
   return(initialAltitudeRate);
}


// double getLandingTime(void) method

/// This method returns the time duration of the landing process,
/// which is set to a fixed value. 

double AircraftLandingProblem::getLandingTime(void)
{
   return(landingTime);
}


// double getTolerance(void) method

/// This method returns the tolerance value to be used in the Runge-Kutta-Fehlberg method for evaluating 
/// different terms in the objective functional.

double AircraftLandingProblem::getTolerance(void)
{
   return(tolerance);
}


// int getInitialSize(void) method

/// This method returns the number of points to be reserved when using the Runge-Kutta-Fehlberg method.

int AircraftLandingProblem::getInitialSize(void)
{
   return(initialSize);
}


// void setShortPeriodGain(double) method

/// This method sets a new short period gain value. That is, it modifies the properties of the aircraft. 
///
/// @param newShortPeriodGain Short period gain value. 

void AircraftLandingProblem::setShortPeriodGain(double newShortPeriodGain)
{
   shortPeriodGain = newShortPeriodGain;
}


//void setShortPeriodResonantFrequency(double) method

/// This method sets a new short period resonant frequency. That is, it modifies the properties of the aircraft. 
///
/// @param newShortPeriodResonantFrequency Short period resonant frequency. 

void AircraftLandingProblem::setShortPeriodResonantFrequency(double newShortPeriodResonantFrequency)
{
   shortPeriodResonantFrequency = newShortPeriodResonantFrequency;
}


// void setPathTimeConstant(double) method

/// This method sets a new path time constant value. That is, it modifies the properties of the aircraft. 
///
/// @param newPathTimeConstant Path time constant value. 

void AircraftLandingProblem::setPathTimeConstant(double newPathTimeConstant)
{
   pathTimeConstant = newPathTimeConstant;
}


// void setShortPeriodDampingFactor(double) method

/// This method sets a new short period damping factor. That is, it modifies the properties of the aircraft. 
///
/// @param newShortPeriodDampingFactor Short period damping factor. 

void AircraftLandingProblem::setShortPeriodDampingFactor(double newShortPeriodDampingFactor)
{
   shortPeriodDampingFactor = newShortPeriodDampingFactor;
}


// void setVelocity(double) method

/// This method sets a new aircraft constant velocity during the landing process. 
///
/// @param newVelocity Velocity value. 

void AircraftLandingProblem::setVelocity(double newVelocity)
{
   velocity = newVelocity;
}


// void setInitialPitchAngleRate(double) method

/// This method sets a new initial value for the pitch angle rate (state variable).
///
/// @param newInitialPitchAngleRate Initial pitch angle rate. 

void AircraftLandingProblem::setInitialPitchAngleRate(double newInitialPitchAngleRate)
{
   initialPitchAngleRate = newInitialPitchAngleRate;
}


// void setInitialPitchAngle(double) method

/// This method sets a new initial value for the pitch angle (state variable).
///
/// @param newInitialPitchAngle Initial pitch angle value. 

void AircraftLandingProblem::setInitialPitchAngle(double newInitialPitchAngle)
{
   initialPitchAngle = newInitialPitchAngle;
}


// void setInitialAltitudeRate(double) method

/// This method sets a new initial value for the altitude rate (state variable).
///
/// @param newInitialAltitudeRate Initial altitude rate. 

void AircraftLandingProblem::setInitialAltitudeRate(double newInitialAltitudeRate)
{
   initialAltitudeRate = newInitialAltitudeRate;
}

// void setInitialAltitude(double) method

/// This method sets a new initial value for the altitude (state variable).
///
/// @param newInitialAltitude Initial altitude value. 

void AircraftLandingProblem::setInitialAltitude(double newInitialAltitude)
{
   initialAltitude = newInitialAltitude;

}


// void setLandingTime(double) method

/// This method sets a new fixed time for the landing process.
///
/// @param newLandingTime Landing time value. 

void AircraftLandingProblem::setLandingTime(double newLandingTime)
{
   landingTime = newLandingTime;
}


// void setTolerance(double) method

/// This method sets a new tolerance value to be used in the Runge-Kutta-Fehlberg method for evaluating the 
/// different terms in the objective functional.
///
/// @param newTolerance Tolerance in Runge-Kutta-Fehlberg method.

void AircraftLandingProblem::setTolerance(double newTolerance)
{
   tolerance = newTolerance;
}


// void setInitialSize(int) method

/// This method sets a new number of points to be reserved when using the Runge-Kutta-Fehlberg method.
///
/// @param newInitialSize Number of points to reserve in Runge-Kutta-Fehlberg method.

void AircraftLandingProblem::setInitialSize(int newInitialSize)
{
   initialSize = newInitialSize;
}


// double calculatePitchAngleRateDot(double, double, double, double, double) method

/// This method depicts the state equation for the pitch angle rate in the aircraft landing problem. 
/// It returns the derivative of the pitch angle rate with respect to the time as a function of time, pitch angle 
/// rate pitch angle, altitude rate and altitude.
///
/// @param time Actual time.
/// @param pitchAngleRate Actual pitch angle rate.
/// @param pitchAngle Actual pitch angle.
/// @param altitudeRate Actual altitude rate.
/// @param altitude Actual altitude.

double AircraftLandingProblem::calculatePitchAngleRateDot(double time, 
double pitchAngleRate, double pitchAngle, double altitudeRate, double altitude)
{
   double pitchAngleRateDot = 0.0;

   double elevatorDeflection = 0.0;

   Vector<double> input(1), output(1);

   input[0] = time;
   output = multilayerPerceptron->calculateOutput(input);

    if(output[0] < ((-35.0/360.0)*(2.0*pi)))
    {
      elevatorDeflection = ((-35.0/360.0)*(2.0*pi));
    }
    else if(output[0] > ((15.0/360.0)*(2.0*pi)))
    {
      elevatorDeflection = ((15.0/360.0)*(2.0*pi));        
    }
    else
    {
      elevatorDeflection = output[0];
    }

   double b11 = 1.0/pathTimeConstant 
              - 2.0*shortPeriodDampingFactor*shortPeriodResonantFrequency;

   double b12 = (2.0*shortPeriodDampingFactor*shortPeriodResonantFrequency) / pathTimeConstant
              - pow(shortPeriodResonantFrequency, 2.0)
              - 1.0/( pow(pathTimeConstant,2.0));

   double b13 = 1.0/(velocity*pow(pathTimeConstant,2.0))
              - ( 2.0*shortPeriodDampingFactor*shortPeriodResonantFrequency) / (velocity*pathTimeConstant)   
              + ( pow(shortPeriodResonantFrequency,2.0) / velocity );

   double c1 = pow(shortPeriodResonantFrequency,2.0)*shortPeriodGain*pathTimeConstant;

   pitchAngleRateDot = b11*pitchAngleRate + b12*pitchAngle + b13*altitudeRate
                     + c1*elevatorDeflection;

   return(pitchAngleRateDot);
}


// double calculatePitchAngleDot(double, double, double, double, double) method

/// This method depicts the state equation for the pitch angle in the aircraft landing problem. 
/// It returns the derivative of the pitch angle with respect to the time as a function of time, pitch angle 
/// rate pitch angle, altitude rate and altitude.
///
/// @param time Actual time.
/// @param pitchAngleRate Actual pitch angle rate.
/// @param pitchAngle Actual pitch angle.
/// @param altitudeRate Actual altitude rate.
/// @param altitude Actual altitude.

double AircraftLandingProblem::calculatePitchAngleDot(double time, 
double pitchAngleRate, double pitchAngle, double altitudeRate, double altitude)
{
   double pitchAngleDot = pitchAngleRate;

   return(pitchAngleDot);
}


// double calculateAltitudeRateDot(double, double, double, double, double) method

/// This method depicts the state equation for the altitude rate in the aircraft landing problem. It returns the 
/// derivative of the altitude rate with respect to the time as a function of time, pitch angle rate, pitch angle,
/// altitude rate and altitude.
///
/// @param time Actual time.
/// @param pitchAngleRate Actual pitch angle rate.
/// @param pitchAngle Actual pitch angle.
/// @param altitudeRate Actual altitude rate.
/// @param altitude Actual altitude.

double AircraftLandingProblem::calculateAltitudeRateDot(double time, 
double pitchAngleRate, double pitchAngle, double altitudeRate, double altitude)
{
   double altitudeRateDot = 0.0;

   double b32 = velocity/pathTimeConstant;

   double b33 = -(1.0/pathTimeConstant);

   altitudeRateDot = b32*pitchAngle + b33*altitudeRate;

   return(altitudeRateDot);
}


// double calculateAltitudeDot(double, double, double, double, double) method

/// This method depicts the state equation for the altitude in the aircraft landing problem. It returns the 
/// derivative of the altitude with respect to the time as a function of time, pitch angle rate, pitch angle, 
/// altitude rate and altitude.
///
/// @param time Actual time.
/// @param pitchAngleRate Actual pitch angle rate.
/// @param pitchAngle Actual pitch angle.
/// @param altitudeRate Actual altitude rate.
/// @param altitude Actual altitude.

double AircraftLandingProblem::calculateAltitudeDot(double time, 
double pitchAngleRate, double pitchAngle, double altitudeRate, double altitude)
{
   double altitudeDot = altitudeRate;

   return(altitudeDot);
}


// double getPitchAngleRateGoal(double) method

/// This method returns the desired behavior for the pith angle rate variable as a function of time in the 
/// aircraft landing problem.
///
/// @param time Time.

double AircraftLandingProblem::getPitchAngleRateGoal(double time)
{
   double pitchAngleRateGoal = 0.0;
    
   return(pitchAngleRateGoal);
}


// double getPitchAngleGoal(double) method

/// This method returns the desired behavior for the pith angle variable as a function of time in the aircraft 
/// landing problem.
///
/// @param time Time.

double AircraftLandingProblem::getPitchAngleGoal(double time)
{   
   double pitchAngleGoal = 0;
    
   if(time==20.0)
   {
      pitchAngleGoal = (1.0/180.0)*pi;  
   }
    
   return(pitchAngleGoal);
}


// double getAltitudeRateGoal(double) method

/// This method returns the desired behavior for the altitude rate variable as a function of time in the aircraft 
/// landing problem.
///
/// @param time Time.

double AircraftLandingProblem::getAltitudeRateGoal(double time)
{
   double altitudeRateGoal = 0.0;

   if( time <= 15.0 ) 
   {
      altitudeRateGoal = -20.0 * exp(-0.2 * time);   
   } 
   else if( time > 15.0 ) 
   {
      altitudeRateGoal = -1.0;   
   }

   return(altitudeRateGoal);
}


// double getAltitudeGoal(double) method

/// This method returns the desired behavior for the altitude variable as a function of time in the aircraft 
/// landing problem.
///
/// @param time Time.

double AircraftLandingProblem::getAltitudeGoal(double time)
{
   double altitudeGoal = 0.0;
   
   if( time <= 15.0  ) 
   {
      altitudeGoal = 100.0 * exp(-0.2 * time);   
   } 
   else if( time > 15.0 ) 
   {
      altitudeGoal = 20.0 - time;   
   }

   return(altitudeGoal);
}


// void setStateVariableWeights(double,double,double, double) method

/// This sets the weight factors of the state variables 
///
/// @param newPitchAngleRateWeight Pitch angle rate weight factor.
/// @param newPitchAngleWeight Pitch angle weight factor.
/// @param newAltitudeRateWeight Altitude rate weight factor.
/// @param newAltitudeWeight Altitude weight factor.

void AircraftLandingProblem::setStateVariableWeights(double newPitchAngleRateWeight,double newPitchAngleWeight,
double newAltitudeRateWeight, double newAltitudeWeight)
{
   pitchAngleRateWeight = newPitchAngleRateWeight;
   
   pitchAngleWeight = newPitchAngleWeight;
   
   altitudeRateWeight = newAltitudeRateWeight;
   
   altitudeWeight = newAltitudeWeight;
}


// void setLandingVariablesWeights(double, double, double, double) method

/// This sets the weight factors of the state variables at landing
///
/// @param newLandingPitchAngleRateWeight Landing Pitch Angle Rate Weight Factor.
/// @param newLandingPitchAngleWeight Landing Pitch Angle Weight Factor.
/// @param newLandingAltitudeRateWeight Landing Altitude Rate Weight Factor.
/// @param newLandingAltitudeWeight Landing Altitude Weight Factor.

void AircraftLandingProblem::setLandingVariablesWeights(double newLandingPitchAngleRateWeight,
double newLandingPitchAngleWeight, double newLandingAltitudeRateWeight, double newLandingAltitudeWeight)
{
   landingPitchAngleRateWeight = newLandingPitchAngleRateWeight;
    
   landingPitchAngleWeight = newLandingPitchAngleWeight;
   
   landingAltitudeRateWeight = newLandingAltitudeRateWeight;
   
   landingAltitudeWeight = newLandingAltitudeWeight;
}


// void setStateVariableWeights(double, double, double, double) method

/// This method sets the weight factors of the state variables. 
///
/// @param newElevatorWeight Elevator Control Weight Factor

void AircraftLandingProblem :: setControlWeight(double newElevatorWeight)
{
   elevatorWeight = newElevatorWeight;
}


// double getPitchAngleRateWeight(double) method

/// This method returns the weight factor for the pitch angle rate variable as a function of time in the aircraft
/// landing problem.
///
/// @param time Time.

double AircraftLandingProblem :: getPitchAngleRateWeight(double time)
{   
    double weightFactor = pitchAngleRateWeight;

    return(weightFactor);
}


// double getPitchAngleWeight(double) method

/// This method returns the weight factor for the pitch angle variable as a function of time in the aircraft 
/// landing problem.
///
/// @param time Time.

double AircraftLandingProblem :: getPitchAngleWeight(double time) 
{   
   double weightFactor = pitchAngleWeight;
     
   return(weightFactor);
 
}


// double getAltitudeRateWeight(double) method

/// This method returns the weight factor for the altitude rate variable as a function of time in the aircraft 
/// landing problem.
///
/// @param time Time.

double AircraftLandingProblem :: getAltitudeRateWeight(double time) 
{
   double weightFactor = altitudeRateWeight;
   
   return(weightFactor);   
}


// double getAltitudeWeight(double) method

/// This method returns the weight factor for the altitude variable as a function of time in the aircraft landing 
/// problem.
///
/// @param time Time.

double AircraftLandingProblem :: getAltitudeWeight(double time) 
{   
   double weightFactor = altitudeWeight;
      
   return(weightFactor);  
}


// double getLandingPitchAngleWeight(void) method

/// This method returns the weight factor at landing for the pitch angle. 

double AircraftLandingProblem :: getLandingPitchAngleRateWeight(void)
{
   
   double weightFactor = landingPitchAngleRateWeight;
      
   return(weightFactor);
 
}


// double getLandingPitchAngleWeight(void) method

/// This method returns the weight factor at landing for the pitch angle. 

double AircraftLandingProblem :: getLandingPitchAngleWeight(void)
{
   
   double weightFactor = landingPitchAngleWeight;
      
   return(weightFactor);
}


// double getLandingAltitudeWeight(void) method

/// This method returns the weight factor at landing for the altitude variable. 

double AircraftLandingProblem :: getLandingAltitudeRateWeight(void)
{
   double weightFactor = landingAltitudeRateWeight;
      
   return(weightFactor);
}


// double getLandingAltitudeWeight(void) method

/// This method returns the weight factor at landing for the altitude variable. 

double AircraftLandingProblem :: getLandingAltitudeWeight(void)
{
   double weightFactor = landingAltitudeWeight;
      
   return(weightFactor);
  
}


// double getElevatorWeight(void) method

/// This method returns the weight factor for the elevator control.

double AircraftLandingProblem :: getElevatorWeight(void)
{
   
   double weightFactor = elevatorWeight;
      
   return(weightFactor);
  
}


// void printWeights(void) method

/// This method prints to the screen all the weights used for the different terms in the objective functional 
/// expression. 

void AircraftLandingProblem :: printWeights(void)
{
    std::cout << "Altitude: " << getInitialAltitude() << std::endl;

    std::cout << "Pitch Angle Rate Weight Factor         : " 
              << getPitchAngleRateWeight(0.0) << std::endl;
              
    std::cout << "Pitch Angle Weight Factor              : " 
              << getPitchAngleWeight(0.0) << std::endl;              
    
    std::cout << "Altitude Rate Weight Factor            : " 
              << getAltitudeRateWeight(0.0) << std::endl;              
              
    std::cout << "Altitude Weight Factor                 : " 
              << getAltitudeWeight(0.0) << std::endl;              

    std::cout << "Landing Pitch Angle Rate Weight Factor : " 
              << getLandingPitchAngleRateWeight() << std::endl;                  

    std::cout << "Landing Pitch Angle Weight Factor      : " 
              << getLandingPitchAngleWeight() << std::endl;                  

    std::cout << "Landing Altitude Rate Weight Factor    : " 
              << getLandingAltitudeRateWeight() << std::endl;                  

    std::cout << "Landing Altitude Weight Factor         : " 
              << getLandingAltitudeWeight() << std::endl;                  

    std::cout << "Elevator Weight Factor                 : " 
              << getElevatorWeight() << std::endl;    
}


// double calculateEvaluation(void) method

/// This method returns the performance of a multilayer perceptron for the aircraft landing problem. 

double AircraftLandingProblem::calculateEvaluation()
{
    // Control sentence 

   if(multilayerPerceptron == NULL)
   {
      std::cerr << std::endl
                << "Flood Error: AircraftLandingProblem class." << std::endl
                << "double calculateEvaluation(void) method." << std::endl
                << "Pointer to multilayer perceptron object cannot be NULL." << std::endl
                << std::endl;
      
      exit(1);
   }

   // Increment number of performance evaluations

   numberOfEvaluations++;

   // Evaluate performance

   double evaluation = 0.0;

   Vector<double> input(1), output(1);

   // Solve state equations 

   Vector<double> time;

   Vector<double> pitchAngleRate;
   Vector<double> pitchAngle;
   Vector<double> altitudeRate;
   Vector<double> altitude;
       
   int numberOfPoints = ordinaryDifferentialEquations.calculateRungeKuttaFehlbergIntegral(*this, 
   time, pitchAngleRate, pitchAngle, altitudeRate, altitude,
   &AircraftLandingProblem::calculatePitchAngleRateDot,
   &AircraftLandingProblem::calculatePitchAngleDot,
   &AircraftLandingProblem::calculateAltitudeRateDot,
   &AircraftLandingProblem::calculateAltitudeDot,
   0.0, 20.0, 
   initialPitchAngleRate, initialPitchAngle, 
   initialAltitudeRate, initialAltitude);
       
   // Calculate performance 

   Vector<double> integrand(numberOfPoints);

   Vector<double> altitudeGoal(numberOfPoints);
   Vector<double> elevatorDeflection(numberOfPoints);

   for(int i = 0; i < numberOfPoints; i++)
   {   
      altitudeGoal[i] = getAltitudeGoal(time[i]);
            
      // Elevator deflection
      
      input[0] = time[i];
      output = multilayerPerceptron->calculateOutput(input);
      
      if(output[0] < ((-35.0/360.0)*(2.0*pi)))
      {
         elevatorDeflection[i] = ((-35.0/360.0)*(2.0*pi));
      }
      else if(output[0] > ((15.0/360.0)*(2.0*pi)))
      {
         elevatorDeflection[i] = ((15.0/360.0)*(2.0*pi));        
      }
      else
      {
         elevatorDeflection[i] = output[0];
      }   

      // Integrand

      integrand[i] = 1.0e-6*pow(altitudeGoal[i]-altitude[i], 2.0) + 1.0e-12*pow(elevatorDeflection[i], 2.0);                       
   }

   evaluation = integrationOfFunctions.calculateSimpsonIntegral(time, integrand);   
        
   evaluation += 1.0e-3*pow((1.0*pi/180.0)-pitchAngle[numberOfPoints-1], 2.0 );
   evaluation += 1.0e-6*pow(altitude[numberOfPoints-1], 2.0);
   
   return(evaluation);
}


// void saveResults(char*) method

/// This method saves the behavior goal, the trajectory and the control signal for the aircraft landing problem 
/// to a data file. 
///
/// <ol>
/// <li> Time.
/// <li> Pitch angle rate goal.
/// <li> Pitch angle goal.
/// <li> Altitude rate goal.
/// <li> Altitude goal.
/// <li> Pitch angle rate.
/// <li> Pitch angle.
/// <li> Altitude rate.
/// <li> Altitude.
/// <li> Elevator deflection.
/// </ol>
///
/// @param filename Filename.

void AircraftLandingProblem::saveResults(char* filename)
{
   Vector<double> input(1), output(1);
   
   double elevatorDeflection = 0.0;

   // Solve state equations 

   Vector<double> time;

   Vector<double> pitchAngleRate;
   Vector<double> pitchAngle;
   Vector<double> altitudeRate;
   Vector<double> altitude;

   int numberOfPoints = 
    ordinaryDifferentialEquations.calculateRungeKuttaFehlbergIntegral(*this, 
   time, pitchAngleRate, pitchAngle, altitudeRate, altitude,
   &AircraftLandingProblem::calculatePitchAngleRateDot,
   &AircraftLandingProblem::calculatePitchAngleDot,
   &AircraftLandingProblem::calculateAltitudeRateDot,
   &AircraftLandingProblem::calculateAltitudeDot,
   0.0, landingTime, 
   initialPitchAngleRate, initialPitchAngle, 
   initialAltitudeRate, initialAltitude);

   Vector<double> pitchAngleRateGoal(numberOfPoints);
   Vector<double> pitchAngleGoal(numberOfPoints);
   Vector<double> altitudeRateGoal(numberOfPoints);
   Vector<double> altitudeGoal(numberOfPoints);

   std::fstream file; 

   file.open(filename, std::ios::out);

   if(!file.is_open())
   {
      std::cout << std::endl
                << "Cannot open aircraft landing problem results data file."  << std::endl;
      
      exit(1);
   }
   else
   {
      if(display)
      { 
         std::cout << std::endl
                    << "Saving aircraft landing problem results to data file..." << std::endl;
      }
   }

   // File header

   file << "%  1 - Time" << std::endl
        << "%  2 - Pitch angle rate goal" << std::endl
        << "%  3 - Pitch angle goal" << std::endl
        << "%  4 - Altitude rate goal" << std::endl
        << "%  5 - Altitude goal" << std::endl
        << "%  6 - Pitch angle rate" << std::endl
        << "%  7 - Pitch angle" << std::endl
        << "%  8 - Altitude rate" << std::endl
        << "%  9 - Altitude" << std::endl
        << "% 10 - Elevator deflection" << std::endl << std::endl 
        << "% Weight Factors" << std::endl << std::endl
        << "% Pitch Angle Rate Weight Factor         = " << getPitchAngleRateWeight(0.0) << std::endl
        << "% Pitch Angle Weight Factor              = " << getPitchAngleWeight(0.0) << std::endl
        << "% Altitude Rate Weight Factor            = " << getAltitudeRateWeight(0.0) << std::endl
        << "% Altitude Weight Factor                 = " << getAltitudeWeight(0.0) << std::endl
        << "% Landing Pitch Angle Rate Weight Factor = " << getLandingPitchAngleRateWeight() << std::endl
        << "% Landing Pitch Angle Weight Factor      = " << getLandingPitchAngleWeight() << std::endl
        << "% Landing Altitude Rate Weight Factor    = " << getLandingAltitudeRateWeight() << std::endl
        << "% Landing Altitude Weight Factor         = " << getLandingAltitudeWeight() << std::endl
        << "% Elevator Weight Factor                 = " << getElevatorWeight() << std::endl
        << std::endl;

   // File data

   for(int i = 0; i < numberOfPoints; i++)
   {            
      pitchAngleRateGoal[i] = getPitchAngleRateGoal(time[i]);
      pitchAngleGoal[i] = getPitchAngleGoal(time[i]);
      altitudeRateGoal[i] = getAltitudeRateGoal(time[i]);
      altitudeGoal[i] = getAltitudeGoal(time[i]);

      input[0] = time[i];
      output = multilayerPerceptron->calculateOutput(input);

      if(output[0] < ((-35.0/360.0)*(2.0*pi)))
      {
         elevatorDeflection = ((-35.0/360.0)*(2.0*pi));
      }
      else if(output[0] > ((15.0/360.0)*(2.0*pi)))
      {
         elevatorDeflection = ((15.0/360.0)*(2.0*pi));        
      }
      else
      {
         elevatorDeflection = output[0];
      }

      file << time[i] << " " 
           << pitchAngleRateGoal[i] << " "
           << pitchAngleGoal[i] << " "
           << altitudeRateGoal[i] << " "
           << altitudeGoal[i] << " "
           << pitchAngleRate[i] << " "
           << pitchAngle[i] << " "
           << altitudeRate[i] << " "
           << altitude[i] << " "
           << elevatorDeflection << std::endl;   
   }  

   file.close();
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
