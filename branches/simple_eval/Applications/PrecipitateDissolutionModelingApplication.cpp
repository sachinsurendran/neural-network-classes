/******************************************************************************/
/*                                                                            */
/*   Flood: An Open Source Neural Networks C++ Library                        */
/*   www.cimne.com/flood                                                      */
/*                                                                            */
/*   M I C R O S T R U C T U R A L   M O D E L I N G   A P P L I C A T I O N  */
/*                                                                            */
/*   Roberto Lopez                                                            */ 
/*   International Center for Numerical Methods in Engineering (CIMNE)        */
/*   Technical University of Catalonia (UPC)                                  */
/*   Barcelona, Spain                                                         */
/*   E-mail: rlopez@cimne.upc.edu                                             */ 
/*                                                                            */  
/******************************************************************************/

/// This main function can be used as a template for constructing a  
/// application with Flood.  

// System includes

#include <iostream>
#include <time.h>

// Utilities includes

#include "../Flood/Utilities/Vector.h"
#include "../Flood/Utilities/Matrix.h"

// Multilayer perceptron includes

#include "../Flood/MultilayerPerceptron/MultilayerPerceptron.h"

// Objective functional includes

#include "../Flood/ObjectiveFunctional/PrecipitateDissolutionModeling.h"

// Training algorithm includes

#include "../Flood/TrainingAlgorithm/RandomSearch.h"
#include "../Flood/TrainingAlgorithm/QuasiNewtonMethod.h"

using namespace Flood;

int main(void)
{
   std::cout << std::endl
             << "Flood Neural Network. Precipitate Dissolution Modeling Application." << std::endl;	

//   srand((unsigned)time(NULL));

   // Multilayer perceptron objects

   MultilayerPerceptron mlp(1,3,1);
   
   mlp.setNameOfSingleInputVariable(0, "log(t/t*)"); 
   mlp.setNameOfSingleOutputVariable(0, "1-f/f0"); 

   mlp.setMinimumOfSingleInputVariable(0, -6.0);
   mlp.setMaximumOfSingleInputVariable(0, 6.0);

   mlp.setMinimumOfSingleOutputVariable(0, 0.0);
   mlp.setMaximumOfSingleOutputVariable(0, 1.0);

   mlp.setLowerBoundOfSingleOutputVariable(0, 0.0);
   mlp.setUpperBoundOfSingleOutputVariable(0, 1.0);

   mlp.setNumberOfIndependentParameters(1);
   mlp.setNameOfSingleIndependentParameter(0, "EffectiveActivationEnergy");

   mlp.setMinimumOfSingleIndependentParameter(0, 100.0);
   mlp.setMaximumOfSingleIndependentParameter(0, 200.0);

   mlp.setLowerBoundOfSingleIndependentParameter(0, 0.0);
   
   mlp.setPreAndPostProcessingMethod(MultilayerPerceptron::MinimumAndMaximum);


   mlp.initFreeParametersNormal(0.0,1.0);

   mlp.setDisplay(false);

   // Activation energy object

   PrecipitateDissolutionModeling pdm(&mlp);

   pdm.setMinkowskiParameter(1.0);

   pdm.setRegularizationWeight(0.0075);

   pdm.loadVickersHardnessTest("../Data/PrecipitateDissolutionModeling/AA-2014-T6/VickersHardnessTestAA-2014-T6.dat");  

   //pdm.printVickersHardnessTest();

   // Random search object

   RandomSearch rs(&pdm);

   rs.train();

   // Quasi-Newton method object

   QuasiNewtonMethod qnm(&pdm);

   qnm.setMinimumImprovement(0.0);
   qnm.setMaximumNumberOfEpochs(1000);
   qnm.setDisplayPeriod(100);

   //ae.printExperimentalData();

   pdm.setNumberOfEvaluations(0);
   qnm.train();

   int numberOfEvaluations = pdm.getNumberOfEvaluations();

   std::cout << std::endl
	         << "Number of evaluations: " << numberOfEvaluations << std::endl; 

   mlp.save("..\\Data\\PrecipitateDissolutionModeling\\AA-7449-T79\\MultilayerPerceptronThreeNeuronsAA-7449-T79.dat");

   mlp.saveExpression("..\\Data\\PrecipitateDissolutionModeling\\AA-7449-T79\\ExpressionThreeNeuronsAA-7449-T79.dat");

   pdm.savePrecipitateDissolutionModel("..\\Data\\PrecipitateDissolutionModeling\\AA-7449-T79\\PrecipitateDissolutionModelThreeNeuronsAA-7449-T79.dat");
   pdm.saveVickersHardnessModel("..\\Data\\PrecipitateDissolutionModeling\\AA-7449-T79\\VickersHardnessModelOneNeuronAA-7449-T79.dat");
   pdm.saveInverseVickersHardnessTest("..\\Data\\PrecipitateDissolutionModeling\\AA-7449-T79\\InverseVickersHardnessTestThreeNeuronsAA-7449-T79.dat");

   std::cout << std::endl;

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
