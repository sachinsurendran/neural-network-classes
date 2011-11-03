/****************************************************************************************************************/
/*                                                                                                              */
/*   Flood: An Open Source Neural Networks C++ Library                                                          */
/*   www.cimne.com/flood                                                                                        */
/*                                                                                                              */
/*   M U L T I L A Y E R   P E R C E P T R O N   C L A S S   H E A D E R                                        */
/*                                                                                                              */
/*   Roberto Lopez                                                                                              */
/*   International Center for Numerical Methods in Engineering (CIMNE)                                          */
/*   Technical University of Catalonia (UPC)                                                                    */
/*   Barcelona, Spain                                                                                           */
/*   E-mail: rlopez@cimne.upc.edu                                                                               */
/*                                                                                                              */
/****************************************************************************************************************/

#ifndef __MULTILAYERPERCEPTRON_H__
#define __MULTILAYERPERCEPTRON_H__

#include <string>

#include "../Perceptron/Perceptron.h"

#include "../Utilities/Vector.h"
#include "../Utilities/Matrix.h"

namespace Flood
{

/// This class represents the concept of multilayer perceptron with a hidden layer 
/// and an output layer of perceptron neuron models. 
///
/// @see Perceptron.

class MultilayerPerceptron
{
public:

   // ENUMERATIONS

   /// Enumeration of available methods for pre-processing network inputs 
   /// and post-processing network outputs.
   
   enum PreAndPostProcessingMethod{None, MeanAndStandardDeviation, MinimumAndMaximum};

   /// Enumeration of available methods for calculating the Jacobian matrix of the 
   /// multilayer perceptron numerically.

   enum NumericalDifferentiationMethod{ForwardDifferences, CentralDifferences};

private:

   // MEMBERS

   /// Number of network inputs.

   int numberOfInputs;

   /// Number of neurons in the hidden layer of the neural network.

   Vector<int> numbersOfHiddenNeurons;

   /// Number of output neurons in the neural network.

   int numberOfOutputs; 

   /// Activation function of the hidden neurons.

   Vector<Perceptron::ActivationFunction> hiddenLayersActivationFunction;

   /// Activation function of the output neurons.

   Perceptron::ActivationFunction outputLayerActivationFunction;

   /// Network's hidden layer. It is built as a Vector of perceptrons.

   Vector< Vector<Perceptron> > hiddenLayers;

   /// Network's output layer. It is built as a Vector of perceptrons.

   Vector<Perceptron> outputLayer;

   /// Name of input variables.

   Vector<std::string> nameOfInputVariables;

   /// Name of output variables.

   Vector<std::string> nameOfOutputVariables;

   /// Units of input variables.

   Vector<std::string> unitsOfInputVariables;

   /// Units of output variables.

   Vector<std::string> unitsOfOutputVariables;

   /// Description of input variables.

   Vector<std::string> descriptionOfInputVariables;

   /// Description of output variables.

   Vector<std::string> descriptionOfOutputVariables;

   /// Mean of input variables.

   Vector<double> meanOfInputVariables;

   /// Standard deviation of input variables.

   Vector<double> standardDeviationOfInputVariables;

   /// Mean of output variables.

   Vector<double> meanOfOutputVariables;

   /// Standard deviation of output variables.

   Vector<double> standardDeviationOfOutputVariables;

   /// Minimum of input variables.

   Vector<double> minimumOfInputVariables;

   /// Maximum of input variables.

   Vector<double> maximumOfInputVariables;

   /// Minimum of output variables.

   Vector<double> minimumOfOutputVariables;

   /// Maximum of output variables.

   Vector<double> maximumOfOutputVariables;

   /// Lower bound of output variables

   Vector<double> lowerBoundOfOutputVariables;

   /// Upper bound of output variables

   Vector<double> upperBoundOfOutputVariables;

   /// Number of independent parameters.

   int numberOfIndependentParameters;

   /// Independent parameters.

   Vector<double> independentParameters;

   /// Name of independent parameters.

   Vector<std::string> nameOfIndependentParameters;

   /// Units of independent parameters.

   Vector<std::string> unitsOfIndependentParameters;

   /// Description of independent parameters.

   Vector<std::string> descriptionOfIndependentParameters;

   /// Mean of independent parameters.

   Vector<double> meanOfIndependentParameters;

   /// Standard deviation of independent parameters.

   Vector<double> standardDeviationOfIndependentParameters;

   /// Minimum of independent parameters.

   Vector<double> minimumOfIndependentParameters;

   /// Maximum of independent parameters.

   Vector<double> maximumOfIndependentParameters;

   /// Lower bound of independent parameters.

   Vector<double> lowerBoundOfIndependentParameters;

   /// Upper bound of independent parameters.

   Vector<double> upperBoundOfIndependentParameters;
   
   /// Numerical epsilon. It is used for computing the Jacobian matrix of the neural network for a set of input 
   /// signals by means of numerical differentiation.

   double numericalEpsilon;

   /// Display warning messages when inputs are out of the range defined by minimums and maximums.  

   bool displayOutOfRangeWarning;

   /// Display messages to screen. 

   bool display;

   /// Pre and post-processing method variable.

   PreAndPostProcessingMethod preAndPostProcessingMethod;

   /// Numerical differentiation method variable.

   NumericalDifferentiationMethod numericalDifferentiationMethod;

public:

   // GENERAL CONSTRUCTOR

   MultilayerPerceptron(int, int, int);
   MultilayerPerceptron(int, Vector<int>, int);


   // DEFAULT CONSTRUCTOR

   MultilayerPerceptron(void);


   // DESTRUCTOR

   ~MultilayerPerceptron(void);


   // METHODS

   // Get methods

   // Enumerations

   PreAndPostProcessingMethod getPreAndPostProcessingMethod(void);

   NumericalDifferentiationMethod getNumericalDifferentiationMethod(void);

   // Neural network

   int getNumberOfInputs(void);

   Vector<int> getNumbersOfHiddenNeurons(void);

   int getNumberOfOutputs(void);

   Vector<Perceptron::ActivationFunction> getHiddenLayersActivationFunction(void);
   
   Perceptron::ActivationFunction getOutputLayerActivationFunction(void);

   Vector<Perceptron>& getHiddenLayer(void);
   Vector< Vector<Perceptron> > & getHiddenLayers(void);
   Vector<Perceptron>& getOutputLayer(void);

   int getNumberOfHiddenLayers(void);
   int getNumberOfNeuralParameters(void);

   Vector<double> getNeuralParameters(void);   
   void setNeuralParameters(Vector<double>);

   // Name of input and output variables

   Vector<std::string> getNameOfInputVariables(void);
   Vector<std::string> getNameOfOutputVariables(void);

   std::string getNameOfSingleInputVariable(int);
   std::string getNameOfSingleOutputVariable(int);

   // Units of input and output variables

   Vector<std::string> getUnitsOfInputVariables(void);
   Vector<std::string> getUnitsOfOutputVariables(void);

   std::string getUnitsOfSingleInputVariable(int);
   std::string getUnitsOfSingleOutputVariable(int);

   // Description of input and output variables

   Vector<std::string> getDescriptionOfInputVariables(void);
   Vector<std::string> getDescriptionOfOutputVariables(void);

   std::string getDescriptionOfSingleInputVariable(int);
   std::string getDescriptionOfSingleOutputVariable(int);
   
   // All information

   Vector< Vector<std::string> > getAllInformation(void);

   // Mean and standard deviation of input and output variables

   Vector<double> getMeanOfInputVariables(void);
   Vector<double> getStandardDeviationOfInputVariables(void);

   Vector<double> getMeanOfOutputVariables(void);
   Vector<double> getStandardDeviationOfOutputVariables(void);

   Matrix<double> getMeanAndStandardDeviationOfInputVariables(void);
   Matrix<double> getMeanAndStandardDeviationOfOutputVariables(void);

   double getMeanOfSingleInputVariable(int);
   double getStandardDeviationOfSingleInputVariable(int);

   double getMeanOfSingleOutputVariable(int);
   double getStandardDeviationOfSingleOutputVariable(int);

   // Minimum and maximum of input and output variables

   Vector<double> getMinimumOfInputVariables(void);
   Vector<double> getMaximumOfInputVariables(void);

   Vector<double> getMinimumOfOutputVariables(void);
   Vector<double> getMaximumOfOutputVariables(void);

   Matrix<double> getMinimumAndMaximumOfInputVariables(void);
   Matrix<double> getMinimumAndMaximumOfOutputVariables(void);

   double getMinimumOfSingleInputVariable(int);
   double getMaximumOfSingleInputVariable(int);

   double getMinimumOfSingleOutputVariable(int);
   double getMaximumOfSingleOutputVariable(int);

   // All statistics

   Vector< Vector<double> > getAllStatistics(void);

   // Lower and upper bounds of input and output variables

   Vector<double> getLowerBoundOfOutputVariables(void);
   Vector<double> getUpperBoundOfOutputVariables(void);

   Matrix<double> getLowerAndUpperBoundsOfOutputVariables(void);

   double getLowerBoundOfSingleOutputVariable(int);
   double getUpperBoundOfSingleOutputVariable(int);

   // Independent parameters

   int getNumberOfIndependentParameters(void);
   
   Vector<double> scaleIndependentParameters(Vector<double>);
   Vector<double> unscaleIndependentParameters(Vector<double>);

   Vector<double> getIndependentParameters(void);   
   double getSingleIndependentParameter(int);   

   // Name of independent parameters

   Vector<std::string> getNameOfIndependentParameters(void);
   std::string getNameOfSingleIndependentParameter(int);

   // Units of independent parameters

   Vector<std::string> getUnitsOfIndependentParameters(void);
   std::string getUnitsOfSingleIndependentParameter(int);

   // Description of independent parameters

   Vector<std::string> getDescriptionOfIndependentParameters(void);
   std::string getDescriptionOfSingleIndependentParameter(int);

   // Mean and standard deviation of independent parameters

   Vector<double> getMeanOfIndependentParameters(void);
   Vector<double> getStandardDeviationOfIndependentParameters(void);
   Matrix<double> getMeanAndStandardDeviationOfIndependentParameters(void);

   double getMeanOfSingleIndependentParameter(int);
   double getStandardDeviationOfSingleIndependentParameter(int);

   // Minimum and maximum of independent parameters

   Vector<double> getMinimumOfIndependentParameters(void);
   Vector<double> getMaximumOfIndependentParameters(void);
   Matrix<double> getMinimumAndMaximumOfIndependentParameters(void);

   double getMinimumOfSingleIndependentParameter(int);
   double getMaximumOfSingleIndependentParameter(int);

   // Lower and upper bounds of independent parameters

   Vector<double> getLowerBoundOfIndependentParameters(void);
   Vector<double> getUpperBoundOfIndependentParameters(void);
   Matrix<double> getLowerAndUpperBoundsOfIndependentParameters(void);

   double getLowerBoundOfSingleIndependentParameter(int);
   double getUpperBoundOfSingleIndependentParameter(int);

   // Other

   double getNumericalEpsilon(void);
   bool getDisplayOutOfRangeWarning(void);
   bool getDisplay(void);

   // Set methods

   // Enumerations
  
   void setPreAndPostProcessingMethod(PreAndPostProcessingMethod);

   void setNumericalDifferentiationMethod(NumericalDifferentiationMethod);

   // Neural network

   void setNetworkArchitecture(int, Vector<int>, int);

   void setNumbersOfHiddenNeurons(Vector<int>);

   void setHiddenLayersActivationFunction(Vector<Perceptron::ActivationFunction>);
   void setOutputLayerActivationFunction(Perceptron::ActivationFunction);

   // Name of input and output variables

   void setNameOfInputVariables(Vector<std::string>);
   void setNameOfOutputVariables(Vector<std::string>);

   void setNameOfSingleInputVariable(int, std::string);
   void setNameOfSingleOutputVariable(int, std::string);

   // Units of input and output variables

   void setUnitsOfInputVariables(Vector<std::string>);
   void setUnitsOfOutputVariables(Vector<std::string>);

   void setUnitsOfSingleInputVariable(int, std::string);
   void setUnitsOfSingleOutputVariable(int, std::string);

   // Description of input and output variables

   void setDescriptionOfInputVariables(Vector<std::string>);
   void setDescriptionOfOutputVariables(Vector<std::string>);

   void setDescriptionOfSingleInputVariable(int, std::string);
   void setDescriptionOfSingleOutputVariable(int, std::string);

   // All information

   void setAllInformation(Vector< Vector<std::string> >);

   // Mean and standard deviation of input and output variables

   void setMeanOfInputVariables(Vector<double>);
   void setStandardDeviationOfInputVariables(Vector<double>);

   void setMeanOfOutputVariables(Vector<double>);
   void setStandardDeviationOfOutputVariables(Vector<double>);

   void setMeanAndStandardDeviationOfInputVariables(Matrix<double>);
   void setMeanAndStandardDeviationOfOutputVariables(Matrix<double>);

   void setMeanOfSingleInputVariable(int, double);
   void setStandardDeviationOfSingleInputVariable(int, double);

   void setMeanOfSingleOutputVariable(int, double);
   void setStandardDeviationOfSingleOutputVariable(int, double);

   // Minimum and maximum of input and output variables

   void setMinimumOfInputVariables(Vector<double>);
   void setMaximumOfInputVariables(Vector<double>);

   void setMinimumOfOutputVariables(Vector<double>);
   void setMaximumOfOutputVariables(Vector<double>);

   void setMinimumAndMaximumOfInputVariables(Matrix<double>);
   void setMinimumAndMaximumOfOutputVariables(Matrix<double>);

   void setMinimumOfSingleInputVariable(int, double);
   void setMaximumOfSingleInputVariable(int, double);

   void setMinimumOfSingleOutputVariable(int, double);
   void setMaximumOfSingleOutputVariable(int, double);

   // All statistics

   void setAllStatistics(Vector< Vector<double> >);

   // Lower and upper bounds of output variables

   void setLowerBoundOfOutputVariables(Vector<double>);
   void setUpperBoundOfOutputVariables(Vector<double>);

   void setLowerAndUpperBoundsOfOutputVariables(Matrix<double>);

   void setLowerBoundOfSingleOutputVariable(int, double);
   void setUpperBoundOfSingleOutputVariable(int, double);

   // Independent parameters

   void setNumberOfIndependentParameters(int);
   void setIndependentParameters(Vector<double>);

   void setSingleIndependentParameter(int, double);

   // Name of independent parameters

   void setNameOfIndependentParameters(Vector<std::string>);
   void setNameOfSingleIndependentParameter(int, std::string);

   // Units of independent parameters

   void setUnitsOfIndependentParameters(Vector<std::string>);
   void setUnitsOfSingleIndependentParameter(int, std::string);

   // Description of independent parameters

   void setDescriptionOfIndependentParameters(Vector<std::string>);
   void setDescriptionOfSingleIndependentParameter(int, std::string);

   // Mean and standard deviation of independent parameters

   void setMeanOfIndependentParameters(Vector<double>);
   void setStandardDeviationOfIndependentParameters(Vector<double>);
   
   void setMeanAndStandardDeviationOfIndependentParameters(Matrix<double>);

   void setMeanOfSingleIndependentParameter(int, double);
   void setStandardDeviationOfSingleIndependentParameter(int, double);

   // Minimum and maximum of independent parameters

   void setMinimumOfIndependentParameters(Vector<double>);
   void setMaximumOfIndependentParameters(Vector<double>);

   void setMinimumAndMaximumOfIndependentParameters(Matrix<double>);

   void setMinimumOfSingleIndependentParameter(int, double);
   void setMaximumOfSingleIndependentParameter(int, double);

   // Lower and upper bounds of independent parameters

   void setLowerBoundOfIndependentParameters(Vector<double>);
   void setUpperBoundOfIndependentParameters(Vector<double>);

   void setLowerAndUpperBoundsOfIndependentParameters(Matrix<double>);

   void setLowerBoundOfSingleIndependentParameter(int, double);
   void setUpperBoundOfSingleIndependentParameter(int, double);

   // Other

   void setDisplayOutOfRangeWarning(bool);
   void setDisplay(bool);
   void setNumericalEpsilon(double);

  // Free parameter methods

   int getNumberOfFreeParameters(void);

   Vector<double> getFreeParameters(void);   
   void setFreeParameters(Vector<double>);

   // Initialization methods

   void initNeuralParametersUniform(void);
   void initNeuralParametersUniform(double, double);
   void initNeuralParametersUniform(Vector<double>, Vector<double>);
   void initNeuralParametersUniform(Matrix<double>);

   void initNeuralParametersNormal(void);
   void initNeuralParametersNormal(double, double);
   void initNeuralParametersNormal(Vector<double>, Vector<double>);
   void initNeuralParametersNormal(Matrix<double>);

   void initIndependentParametersUniform(void);
   void initIndependentParametersUniform(double, double);
   void initIndependentParametersUniform(Vector<double>, Vector<double>);
   void initIndependentParametersUniform(Matrix<double>);

   void initIndependentParametersNormal(void);
   void initIndependentParametersNormal(double, double);
   void initIndependentParametersNormal(Vector<double>, Vector<double>);
   void initIndependentParametersNormal(Matrix<double>);

   void initFreeParametersUniform(void);
   void initFreeParametersUniform(double, double);
   void initFreeParametersUniform(Vector<double>, Vector<double>);
   void initFreeParametersUniform(Matrix<double>);

   void initFreeParametersNormal(void);
   void initFreeParametersNormal(double, double);
   void initFreeParametersNormal(Vector<double>, Vector<double>);
   void initFreeParametersNormal(Matrix<double>);

   // Output methods

   void checkInput(Vector<double>);

   Vector<double> preprocessInput(Vector<double>);

   Vector<double> calculateNetInputSignalToHiddenLayer(int, Vector<double>);
  
   Vector<double> calculateOutputSignalFromHiddenLayer(int, Vector<double>);

   Vector<double> calculateOutputSignalDerivativeFromHiddenLayer(int, Vector<double>);
   
   Vector<double> calculateNetInputSignalToOutputLayer(Vector<double>);

   Vector<double> calculateOutputSignal(Vector<double>);

   Vector<double> calculateOutputSignalDerivative(Vector<double>);

   Vector<double> postprocessOutputSignal(Vector<double>);

   Vector<double> applyLowerAndUpperBounds(Vector<double>);

   Vector<double> calculateOutput(Vector<double>);

   // Jacobian matrix methods

   Matrix<double> calculateJacobian(Vector<double>);

   Matrix<double> calculateJacobianTest(Vector<double>);

   // Sensitivity matrix methods

   Matrix<double> calculateSensitivityTest(Vector<double>);

   // Utility methods

   void print(void);

   void save(char*);
   void load(char*);

   void saveExpression(char*);


   // Vector<double> calculateOutput(Vector<double> input, BoundaryConditionsType& boundaryConditions) method

   /// This method calculates the set of outputs from the neural network in response to a set of inputs, when 
   /// boundary conditions are imposed.
   ///
   /// @param input: Set of inputs to the neural network.
   /// @param boundaryConditions: Boundary conditions given by a set of homogeneous and particular solution terms.
   ///
   /// @see checkInput(Vector<double>).
   /// @see preprocessInput(Vector<double>).
   /// @see calculateNetInputSignalToHiddenLayer(Vector<double>).
   /// @see calculateNetInputSignalToHiddenLayer(Vector<double>).
   /// @see calculateOutputSignalDerivativeFromHiddenLayer(Vector<double>).
   /// @see calculateNetInputSignalToOutputLayer(Vector<double>).
   /// @see calculateOutputSignal(Vector<double>).
   /// @see calculateOutputSignalDerivative(Vector<double>).
   /// @see postprocessOutputSignal(Vector<double>).
   /// @see calculateOutput(Vector<double>).
   ///
   /// @see GeodesicProblem. 
   /// @see BrachistochroneProblem. 
   /// @see CatenaryProblem. 
   /// @see IsoperimetricProblem. 

   template<class BoundaryConditionsType> Vector<double> calculateOutput(Vector<double> input, 
   BoundaryConditionsType& boundaryConditions)
   {
      Vector<double> output(numberOfOutputs);                     
      
      int numberOfHiddenLayers = getNumberOfHiddenLayers();
                     
      // Check input for size and maximum and minimum values
 
      checkInput(input);
 
      // Preprocess input to obtain input signal to the neural network 

      Vector<double> inputSignal = preprocessInput(input);
   
      // Get net input signal to hidden layer

      Vector< Vector<double> > netInputSignalToHiddenLayer(numberOfHiddenLayers);
      netInputSignalToHiddenLayer[0] = calculateNetInputSignalToHiddenLayer(0, inputSignal);

      // Get output signal from hidden layer

      Vector< Vector<double> > outputSignalFromHiddenLayer(numberOfHiddenLayers);
      outputSignalFromHiddenLayer[0] = calculateOutputSignalFromHiddenLayer(0, netInputSignalToHiddenLayer[0]);
   
      for(int i = 1; i < numberOfHiddenLayers; i++)
      {
         netInputSignalToHiddenLayer[i] 
         = calculateNetInputSignalToHiddenLayer(i, outputSignalFromHiddenLayer[i-1]);         

         outputSignalFromHiddenLayer[i] = calculateOutputSignalFromHiddenLayer(i, netInputSignalToHiddenLayer[i]);         
      }

      // Get net input signal to output layer

      Vector<double> netInputSignalToOutputLayer 
      = calculateNetInputSignalToOutputLayer(outputSignalFromHiddenLayer[numberOfHiddenLayers-1]);

      // Get output signal from output layer

      Vector<double> outputSignal = calculateOutputSignal(netInputSignalToOutputLayer);
  
      // Postprocess output signal from the neural network to get output

      output = postprocessOutputSignal(outputSignal);

      // Apply boundary conditions

      Vector<double> particularSolution = boundaryConditions.calculateParticularSolution(input);
      Vector<double> homogeneousSolution = boundaryConditions.calculateHomogeneousSolution(input);

      output = particularSolution + homogeneousSolution*output;

      // Apply lower and upper bounds

      output = applyLowerAndUpperBounds(output);

      return(output);
   }


   // Matrix<double> calculateJacobian(Vector<double>, BoundaryConditionsType&) method

   /// This method returns the the Jacobian matrix of the neural network for a 
   /// set of inputs, corresponding to the point in input space
   /// at which the Jacobian Matrix is to be found.
   /// This method applies when boundary conditions are imposed. 
   /// It uses the back-propagation method. 
   ///
   /// @param input: Set of inputs to the neural network.
   /// @param boundaryConditions: Boundary conditions type.
   ///
   /// @see calculateOutput(Vector<double>).
   /// @see calculateJacobianTest(Vector<double>).
   ///
   /// @todo this method is implemented only for one input and one or two output variables.
   /// An extension of the algorithm to consider any number of inputs and outputs is required. 

   template<class BoundaryConditionsType> Matrix<double> calculateJacobian(Vector<double> input,
   BoundaryConditionsType& boundaryConditions) 
   {
      Vector<double> output = calculateOutput(input);
      Matrix<double> jacobian = calculateJacobian(input);

      // Apply boundary conditions

      Vector<double> particularSolution = boundaryConditions.calculateParticularSolution(input);
      Vector<double> homogeneousSolution = boundaryConditions.calculateHomogeneousSolution(input);
      Vector<double> particularSolutionDerivative = boundaryConditions.calculateParticularSolutionDerivative(input);
      Vector<double> homogeneousSolutionDerivative = boundaryConditions.calculateHomogeneousSolutionDerivative(input);

      for(int i = 0; i < numberOfOutputs; i++)
      {
         for(int j = 0; j < numberOfInputs; j++)
         {
            jacobian[i][0] 
            = particularSolutionDerivative[i] 
            + homogeneousSolutionDerivative[i]*output[i] 
            + homogeneousSolution[i]*jacobian[i][0];           
         }
      }

      return(jacobian);
   }
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
