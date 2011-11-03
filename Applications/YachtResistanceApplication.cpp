/****************************************************************************************************************/
/*                                                                                                              */ 
/*   Flood: An Open Source Neural Networks C++ Library                                                          */
/*   www.cimne.com/flood                                                                                        */
/*                                                                                                              */
/*   Y A C H T   R E S I S T A N C E   A P P L I C A T I O N                                                    */
/*                                                                                                              */
/*   Roberto Lopez                                                                                              */ 
/*   International Center for Numerical Methods in Engineering (CIMNE)                                          */
/*   Technical University of Catalonia (UPC)                                                                    */
/*   Barcelona, Spain                                                                                           */
/*   E-mail: rlopez@cimne.upc.edu                                                                               */ 
/*                                                                                                              */  
/****************************************************************************************************************/

// System includes

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cstring>
#include <time.h>

// Utilities includes

#include "../Flood/Utilities/Vector.h"
#include "../Flood/Utilities/Matrix.h"

#include "../Flood/Utilities/InputTargetDataSet.h"
#include "../Flood/Utilities/LinearRegressionAnalysis.h"

// Multilayer perceptron includes

#include "../Flood/MultilayerPerceptron/MultilayerPerceptron.h"

// Objective functional includes

#include "../Flood/ObjectiveFunctional/SumSquaredError.h"
#include "../Flood/ObjectiveFunctional/MeanSquaredError.h"
#include "../Flood/ObjectiveFunctional/RootMeanSquaredError.h"
#include "../Flood/ObjectiveFunctional/NormalizedSquaredError.h"

// Training algorithm includes

#include "../Flood/TrainingAlgorithm/RandomSearch.h"
#include "../Flood/TrainingAlgorithm/EvolutionaryAlgorithm.h"

#include "../Flood/TrainingAlgorithm/GradientDescent.h"
#include "../Flood/TrainingAlgorithm/ConjugateGradient.h"
#include "../Flood/TrainingAlgorithm/QuasiNewtonMethod.h"

#include "../Flood/TrainingAlgorithm/NewtonMethod.h"

using namespace Flood;

int main(void)
{
   std::cout << std::endl
             << "Flood Neural Network. Yacht Resistance Application." << std::endl;	

   srand((unsigned)time(NULL));

   // Input-target data set objects
  
   InputTargetDataSet trainingDataSet;
   InputTargetDataSet validationDataSet;
   InputTargetDataSet testingDataSet;
   
   trainingDataSet.load("../Data/YachtResistance/TrainingDataSet.dat");
   validationDataSet.load("../Data/YachtResistance/ValidationDataSet.dat");
   testingDataSet.load("../Data/YachtResistance/TestingDataSet.dat");

   int numberOfInputs = trainingDataSet.getNumberOfInputVariables();
   int numberOfOutputs = trainingDataSet.getNumberOfTargetVariables();

   Vector< Vector<std::string> > information = trainingDataSet.getAllInformation();
   Vector< Vector<double> > statistics = trainingDataSet.calculateAllStatistics();

   trainingDataSet.preprocessMeanAndStandardDeviation();  
      
   // Multilayer perceptron object

   int numberOfHiddenNeurons = 9;

   MultilayerPerceptron multilayerPerceptron(numberOfInputs, numberOfHiddenNeurons, numberOfOutputs);

   // Normalized squared error object

   NormalizedSquaredError trainingError(&multilayerPerceptron, &trainingDataSet);

   // Quasi-Newton method object

   QuasiNewtonMethod quasiNewtonMethod(&trainingError);
   
   quasiNewtonMethod.setMaximumNumberOfEpochs(5000);
   quasiNewtonMethod.setDisplayPeriod(50);
   quasiNewtonMethod.setMinimumImprovement(1.0e-12);
   quasiNewtonMethod.setGradientNormGoal(0.0);
   
   quasiNewtonMethod.train();

   // Set all multilayer perceptron stuff for future use

   multilayerPerceptron.setAllInformation(information);
   multilayerPerceptron.setAllStatistics(statistics);

   multilayerPerceptron.setPreAndPostProcessingMethod(MultilayerPerceptron::MeanAndStandardDeviation);

   multilayerPerceptron.save("../Data/YachtResistance/MultilayerPerceptron.dat");
   multilayerPerceptron.saveExpression("../Data/YachtResistance/Expression.dat");

   // Validation error object
   
   NormalizedSquaredError validationError(&multilayerPerceptron, &validationDataSet);
         
   std::cout << validationError.calculateEvaluation() << std::endl;
                  
   // Linear regression analysis object
   
   LinearRegressionAnalysis linearRegressionAnalysis(&multilayerPerceptron, &testingDataSet);
   
   linearRegressionAnalysis.saveResults("../Data/YachtResistance/LinearRegressionAnalysis.dat");

   std::cout << std::endl;

   return(0);
}  
