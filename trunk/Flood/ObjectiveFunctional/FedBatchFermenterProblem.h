/******************************************************************************/
/*                                                                            */
/*   Flood: An Open Source Neural Networks C++ Library                        */
/*   www.cimne.com/flood                                                      */
/*                                                                            */
/*   F E D   B A T C H   F E R M E N T E R   P R O B L E M                    */
/*                                                                            */
/*   C L A S S   H E A D E R                                                  */
/*                                                                            */
/*   Roberto Lopez                                                            */
/*   International Center for Numerical Methods in Engineering (CIMNE)        */
/*   Technical University of Catalonia (UPC)                                  */
/*   Barcelona, Spain                                                         */
/*   E-mail: rlopez@cimne.upc.es                                              */
/*                                                                            */
/******************************************************************************/

#ifndef __FEDBATCHFERMENTERPROBLEM_H__
#define __FEDBATCHFERMENTERPROBLEM_H__

#include "ObjectiveFunctional.h"
#include "../Utilities/OrdinaryDifferentialEquations.h"

namespace Flood
{

/// This class represents the objective functional of a multilayer 
/// perceptron for the fed batch fermenter problem.
/// The fed batch fermenter problem for the multilayer perceptron is 
/// an optimal control problem with one control and four state variables.
/// It is defined by an objective functional with one constraint and 
/// requiring the integration of a system of ordinary differential equations. 
///
/// @see ObjectiveFunctional.
/// @see OrdinaryDifferentialEquations. 

class FedBatchFermenterProblem : public ObjectiveFunctional
{

private:

   /// Fermentation process time (hour).

   double finalTime;

   /// Initial concentration of cell mass in the fermenter (gram/litre). 

   double initialCellMassConcentration;

   /// Initial concentration of substrate in the fermenter (gram/litre). 

   double initialSubstrateConcentration;

   /// Initial concentration of product in the fermenter (gram/litre). 

   double initialProductConcentration;

   /// Initial volume of culture medium in the fermenter (litre). 

   double initialBrothVolume;

   /// Volume of fermenter (litre).

   double fermenterVolume;

   /// Minimum feed rate (control action) to the fermenter (litre/hour).

   double minimumFeedRate;

   /// Maximum feed rate (control action) to the fermenter (litre/hour).

   double maximumFeedRate;
     
   /// Weight for the final volume error term in the objective functional.

   double volumeErrorWeight;

   /// Weight for the final yield term in the objective functional.

   double yieldWeight;
   
   /// Yield coefficient in the fermentation process (no units).

   double yieldCoefficient;

   /// Concentration of substrate in the feed to the fermenter (gram/litre).

   double feedSubstrateConcentration;

   /// Kinetic constant in the fementation process Mu0. 

   double kineticConstantMu0;

   /// Kinetic constant in the fementation process Eta0. 

   double kineticConstantEta0;

   /// Kinetic constant in the fementation process Kp. 

   double kineticConstantKp;

   /// Kinetic constant in the fementation process JpDash. 

   double kineticConstantKpDash;

   /// Kinetic constant in the fementation process Ks. 

   double kineticConstantKs;

   /// Kinetic constant in the fementation process KsDash. 

   double kineticConstantKsDash;

   /// Ordinary differential equations object.   

   OrdinaryDifferentialEquations ordinaryDifferentialEquations;

   /// Tolerance of integration in Runge-Kutta-Fehlberg method.

   double tolerance;

   /// Initial size of solution vectors for the Runge-Kutta-Fehlberg method.

   int initialSize;

public:

   // GENERAL CONSTRUCTOR  

   FedBatchFermenterProblem(MultilayerPerceptron*);


   // DEFAULT CONSTRUCTOR

   FedBatchFermenterProblem(void);


   // DESTRUCTOR

   virtual ~FedBatchFermenterProblem(void);

   // METHODS

   // Get methods

   double getFinalTime(void);

   double getInitialCellMassConcentration(void);
   double getInitialSubstrateConcentration(void);
   double getInitialProductConcentration(void);
   double getInitialBrothVolume(void);

   double getFermenterVolume(void);

   double getMinimumFeedRate(void);
   double getMaximumFeedRate(void);

   double getVolumeErrorWeight(void);
   double getYieldWeight(void);

   double getTolerance(void);
   int getInitialSize(void);

   // Set methods

   void setFinalTime(double);

   void setInitialCellMassConcentration(double);
   void setInitialSubstrateConcentration(double);
   void setInitialProductConcentration(double);
   void setInitialBrothVolume(double);

   void setFermenterVolume(double);

   void setMinimumFeedRate(double);
   void setMaximumFeedRate(double);

   void setVolumeErrorWeight(double);
   void setYieldWeight(double);

   void setTolerance(double);
   void setInitialSize(int);

   // State equation methods

   double getCellMassConcentrationDot(double, double, double, double, double);
   double getSubstrateConcentrationDot(double, double, double, double, double);
   double getProductConcentrationDot(double, double, double, double, double);
   double getBrothVolumeDot(double, double, double, double, double);

   double getSpecificGrowthRate(double, double);
   double getSpecificProductivity(double, double);

   // Objective functional evaluation methods

   double calculateEvaluation(void);

   // Utility methods

   void print(void);

   void saveResults(char*);
};

}

#endif


// Flood: An Open Source Neural Networks C++ Library.
// Copyright (C) 2005-2008 Roberto Lopez 
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
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
