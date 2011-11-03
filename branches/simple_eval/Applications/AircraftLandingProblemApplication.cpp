/****************************************************************************************************************/
/*                                                                                                              */ 
/*   Flood: An Open Source Neural Networks C++ Library                                                          */
/*   www.cimne.com/flood                                                                                        */
/*                                                                                                              */
/*   A I R C R A F T   L A N D I N G   P R O B L E M   A P P L I C A T I O N                                    */
/*                                                                                                              */
/*   Roberto Lopez and Kevin Lau                                                                                */ 
/*   International Center for Numerical Methods in Engineering (CIMNE)                                          */
/*   Technical University of Catalonia (UPC)                                                                    */
/*   Barcelona, Spain                                                                                           */
/*   E-mail: rlopez@cimne.upc.edu, kevin.lau@imperial.ac.uk                                                     */ 
/*                                                                                                              */  
/****************************************************************************************************************/

/// This main function can be used as a template for solving the aircraft landing problem by means of a multilayer 
/// perceptron. 
/// It uses the quasi-Newton method algorithm for training. 

// System includes

#include <iostream>
#include <time.h>
#include <math.h>

// Multilayer perceptron includes

#include "../Flood/MultilayerPerceptron/MultilayerPerceptron.h"

// Objective functional includes

#include "../Flood/ObjectiveFunctional/AircraftLandingProblem.h"

// Training algorithm includes

#include "../Flood/TrainingAlgorithm/RandomSearch.h"
#include "../Flood/TrainingAlgorithm/QuasiNewtonMethod.h"

using namespace Flood;

int main(void)
{
	std::cout << "Aircraft Landing Problem" << std::endl;	

	std::cout << std::endl;

	srand((unsigned)time(NULL));

	// Multilayer perceptron object

	MultilayerPerceptron mlp(1,3,1);    

    mlp.setNameOfSingleInputVariable(0, "Time");
    mlp.setNameOfSingleOutputVariable(0, "ElevatorDeflectionAngle");

    mlp.setMinimumOfSingleInputVariable(0, 0.0);
    mlp.setMaximumOfSingleInputVariable(0, 20.0);

	double pi = 3.1415927;

	mlp.setMinimumOfSingleOutputVariable(0, -2.0*pi/360.0);
    mlp.setMaximumOfSingleOutputVariable(0, 2.0*pi/360.0);

	// Aircraft landing problem object

	AircraftLandingProblem alp(&mlp);
    alp.setInitialAltitude(100.0);
    alp.setStateVariableWeights(500.0, 0.0, 1.0e-8, 5.0e-3);
    alp.setLandingVariablesWeights(0.0, 200.0, 5.0e-4, 5.0e-1);
    alp.setControlWeight(10.0);

	// Optimisation Objects

    RandomSearch rs(&alp);
	rs.setMaximumNumberOfEpochs(100);
	rs.train();

	QuasiNewtonMethod qnm(&alp);
	qnm.setMaximumNumberOfEpochs(100);
	qnm.train();

	std::cout << alp.getNumberOfEvaluations() << std::endl;

	// Save neural network to file

	mlp.save("../Data/AircraftLandingProblem/MultilayerPerceptron.dat");
    
	// Save results to file

	alp.saveResults("../Data/AircraftLandingProblem/Results.dat");
  
    return(0);  
}  


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
