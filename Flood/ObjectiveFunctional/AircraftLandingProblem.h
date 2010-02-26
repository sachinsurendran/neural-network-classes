/****************************************************************************************************************/
/*                                                                                                              */
/*   Flood: An Open Source Neural Networks C++ Library                                                          */
/*   www.cimne.com/flood                                                                                        */
/*                                                                                                              */
/*   A I R C R A F T   L A N D I N G   P R O B L E M   C L A S S   H E A D E R                                  */
/*                                                                                                              */
/*   Roberto Lopez and Kevin Lau                                                                                */ 
/*   International Center for Numerical Methods in Engineering (CIMNE)                                          */
/*   Technical University of Catalonia (UPC)                                                                    */
/*   Barcelona, Spain                                                                                           */
/*   E-mail: rlopez@cimne.upc.es, kevin.lau@imperial.ac.uk                                                      */ 
/*                                                                                                              */
/****************************************************************************************************************/

#ifndef __AIRCRAFTLANDINGPROBLEM_H__
#define __AIRCRAFTLANDINGPROBLEM_H__

#include "ObjectiveFunctional.h"
#include "../Utilities/OrdinaryDifferentialEquations.h"
#include "../Utilities/IntegrationOfFunctions.h"     

namespace Flood
{


/// This class represents the objective functional of a multilayer perceptron for the aircraft landing problem. 
/// This is an optimal control problem of aeronautical engineering interest, with one control and four state 
/// variables, and where the objective functional is evaluated by integrating a system of four ordinary 
///	differential equations.
///
/// @see ObjectiveFunctional.

class AircraftLandingProblem : public ObjectiveFunctional
{

private:

   /// Pi value. 

   double pi;

   /// Short period gain aircraf property value.  

   double shortPeriodGain;

   /// Shor period resonant frequency aircraft property value. 

   double shortPeriodResonantFrequency;

   /// Short period damping factor aircraft property value. 
 
   double shortPeriodDampingFactor;

   /// Path time constant aircraft property value.

   double pathTimeConstant;

   /// Velocity of the aircraft, assumed to be constant. 

   double velocity;

   /// Initial value for the pitch angle.

   double initialPitchAngle;

   /// Initial value for the pitch angle rate.
   
   double initialPitchAngleRate;

   /// Initial value for the altitude.

   double initialAltitude;

   /// Initial value for the altitude rate.

   double initialAltitudeRate;

   /// Landing time.

   double landingTime;

   /// Weight factor for the pitch angle history term in the objective functional.

   double pitchAngleWeight;

   /// Weight factor for the pitch angle rate history term in the objective functional.

   double pitchAngleRateWeight;

   /// Weight factor for the altitude history term in the objective functional.

   double altitudeWeight;

   /// Weight factor for the altitude rate history term in the objective functional.

   double altitudeRateWeight;

   /// Weight factor for the pitch angle at the moment of landing term.

   double landingPitchAngleWeight;

   /// Weight factor for the pitch angle rate at the moment of landing term.

   double landingPitchAngleRateWeight;

   /// Weight factor for the altitude at the moment of landing term.

   double landingAltitudeWeight;

   /// Weight factor for the altitude rate at the moment of landing term.

   double landingAltitudeRateWeight;

   /// Weight factor for the elevator deflection angle term. 

   double elevatorWeight;
   
   /// Ordinary differential equations object.

   OrdinaryDifferentialEquations ordinaryDifferentialEquations;   

   /// Integration of functions object.    

   IntegrationOfFunctions integrationOfFunctions;

   /// Tolerance in the Runge-Kutta-Fehlberg method.

   double tolerance;

   /// Initial number of points to be reserved in the Runge-Kutta-Fehlberg method.

   int initialSize;


public:

     // GENERAL CONSTRUCTOR

   AircraftLandingProblem(MultilayerPerceptron*);


     // DEFAULT CONSTRUCTOR

   AircraftLandingProblem(void);


   // DESTRUCTOR

   virtual ~AircraftLandingProblem(void);

   // METHODS

   // Get methods

   double getShortPeriodGain(void);
   double getPathTimeConstant(void);
   double getShortPeriodResonantFrequency(void);
   double getShortPeriodDampingFactor(void);

   double getVelocity(void);

   double getInitialPitchAngle(void);
   double getInitialPitchAngleRate(void);

   double getInitialAltitude(void);
   double getInitialAltitudeRate(void);

   double getLandingTime(void);

   double getTolerance(void);
   int getInitialSize(void);

   // Set methods

   void setShortPeriodGain(double);
   void setShortPeriodResonantFrequency(double);
   void setShortPeriodDampingFactor(double);
   void setPathTimeConstant(double);

   void setVelocity(double);

   void setInitialPitchAngle(double);
   void setInitialPitchAngleRate(double);

   void setInitialAltitude(double);
   void setInitialAltitudeRate(double);

   void setLandingTime(double);

   void setTolerance(double);
   void setInitialSize(int);

    // Weight Factor Methods

    void setStateVariableWeights(double,double,double,double);
    void setLandingVariablesWeights(double,double,double,double);
    void setControlWeight(double);

    double getPitchAngleRateWeight(double);
    double getPitchAngleWeight(double);
    double getAltitudeRateWeight(double);    
    double getAltitudeWeight(double);
    
    double getLandingPitchAngleRateWeight(void);
    double getLandingPitchAngleWeight(void);
    double getLandingAltitudeRateWeight(void);
    double getLandingAltitudeWeight(void);

    double getElevatorWeight(void);
    
    void printWeights(void);

   // State equation methods

   double calculatePitchAngleRateDot(double, double, double, double, double);
   double calculatePitchAngleDot(double, double, double, double, double);

   double calculateAltitudeRateDot(double, double, double, double, double);
   double calculateAltitudeDot(double, double, double, double, double);

   // Behavior goal methods

   double getPitchAngleGoal(double);
   double getPitchAngleRateGoal(double);
   double getAltitudeGoal(double);
   double getAltitudeRateGoal(double);

   // Evaluation methods

   double calculateEvaluation(void);

   // Utility methods

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
