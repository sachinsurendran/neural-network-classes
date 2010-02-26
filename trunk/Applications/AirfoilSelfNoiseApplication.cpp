/****************************************************************************************************************/
/*                                                                                                              */ 
/*   Flood: An Open Source Neural Networks C++ Library                                                          */
/*   www.cimne.com/flood                                                                                        */
/*                                                                                                              */
/*   A I R F O I L   S E L F - N O I S E   A P P L I C A T I O N                                                */
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

#include "../Flood/ObjectiveFunctional/NormalizedSquaredError.h"

// Training algorithm includes

#include "../Flood/TrainingAlgorithm/QuasiNewtonMethod.h"


using namespace Flood;

int main(void)
{
   std::cout << std::endl
             << "Flood Neural Network. Airfoil Self-Noise Application." << std::endl;	

   srand((unsigned)time(NULL));

   // Input-target data set objects
  
   InputTargetDataSet trainingDataSet;
   InputTargetDataSet validationDataSet;
   InputTargetDataSet testingDataSet;
   
   trainingDataSet.load("../Data/AirfoilSelfNoise/TrainingDataSet.dat");
   validationDataSet.load("../Data/AirfoilSelfNoise/ValidationDataSet.dat");
   testingDataSet.load("../Data/AirfoilSelfNoise/TestingDataSet.dat");
   
   Vector< Vector<std::string> > information = trainingDataSet.getAllInformation();
   Vector< Vector<double> > statistics = trainingDataSet.calculateAllStatistics();

   trainingDataSet.preprocessMeanAndStandardDeviation();  
      
   // Multilayer perceptron object

   int numberOfInputs = trainingDataSet.getNumberOfInputVariables();
   int numberOfHiddenNeurons = 9;
   int numberOfOutputs = trainingDataSet.getNumberOfTargetVariables();

   MultilayerPerceptron multilayerPerceptron(numberOfInputs, numberOfHiddenNeurons, numberOfOutputs);

   // Normalized squared error object

   NormalizedSquaredError trainingError(&multilayerPerceptron, &trainingDataSet);

  // Quasi-Newton method object

   QuasiNewtonMethod quasiNewtonMethod(&trainingError);

   quasiNewtonMethod.setMaximumNumberOfEpochs(1000);
   quasiNewtonMethod.setDisplayPeriod(100);
   quasiNewtonMethod.setMinimumImprovement(1.0e-9);
   quasiNewtonMethod.setGradientNormGoal(0.0);
   
   quasiNewtonMethod.train();

   multilayerPerceptron.setAllInformation(information);
   multilayerPerceptron.setAllStatistics(statistics);

   multilayerPerceptron.setPreAndPostProcessingMethod(MultilayerPerceptron::MeanAndStandardDeviation);

   multilayerPerceptron.save("../Data/AirfoilSelfNoise/MultilayerPerceptron.dat");
   multilayerPerceptron.saveExpression("../Data/AirfoilSelfNoise/Expression.dat");

   // Validation error object
   
   NormalizedSquaredError validationError(&multilayerPerceptron, &validationDataSet);
         
   std::cout << "Validation error:" << std::endl
             << validationError.calculateEvaluation() << std::endl;

   // Linear regression analysis object
                  
   LinearRegressionAnalysis linearRegressionAnalysis(&multilayerPerceptron, &testingDataSet);
   linearRegressionAnalysis.saveResults("../Data/AirfoilSelfNoise/LinearRegressionAnalysis.dat");

   std::cout << std::endl;

   return(0);
}  
