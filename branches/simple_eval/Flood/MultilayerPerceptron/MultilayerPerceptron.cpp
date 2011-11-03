/****************************************************************************************************************/
/*                                                                                                              */
/*   Flood: An Open Source Neural Networks C++ Library                                                          */
/*   www.cimne.com/flood                                                                                        */
/*                                                                                                              */
/*   M U L T I L A Y E R   P E R C E P T R O N   C L A S S                                                      */
/*                                                                                                              */
/*   Roberto Lopez                                                                                              */
/*   International Center for Numerical Methods in Engineering (CIMNE)                                          */
/*   Technical University of Catalonia (UPC)                                                                    */
/*   Barcelona, Spain                                                                                           */
/*   E-mail: rlopez@cimne.upc.edu                                                                               */
/*                                                                                                              */
/****************************************************************************************************************/

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <math.h>   

#include "MultilayerPerceptron.h"

namespace Flood
{

// GENERAL CONSTRUCTOR

/// General constructor. It creates a multilayer perceptron object with one hidden layer and an output layer of 
/// neurons. 
/// All the free parameters in the neural network are initialized with random values chosen from a normal 
/// distribution with mean zero and standard deviation. 
/// This constructor also initializes the rest of class members to their default values:
///
/// <ul>
/// <li> Hidden layer activation function: Hyperbolic tangent.
/// <li> Output layer activation function: Linear.
/// <li> Name of input vaviables: "InputVariable0",... 
/// <li> Name of output variables: "OutputVariable0",...
/// <li> Units of input vaviables: "None".
/// <li> Units of output variables: "None".
/// <li> Description of input vaviables: "None".
/// <li> Description of output variables: "None".
/// <li> Mean of input variables: 0.
/// <li> Standard deviation of input variables: 1.
/// <li> Mean of output variables: 0.
/// <li> Standard deviation of output variables: 1.
/// <li> Minimum of input variables: -Inf.
/// <li> Maximum of input variables: +Inf.
/// <li> Minimum of output variables: -Inf.
/// <li> Maximum of output variables: +Inf.
/// <li> Lower bound of output variables: -Inf.
/// <li> Upper bound of output variables: +Inf.
/// <li> Number of independent parameters: 0.
/// <li> Pre and post-processing method: None. 
/// <li> Numerical epsilon method: Relative. 
/// <li> Numerical differentiation method: Central differences. 
/// <li> Display out of range warning: False.
/// <li> Display: True.
/// <li> Numerical epsilon: 1.0e-6.
/// </ul> 
///
/// @param newNumberOfInputs Number of inputs in the neural network.
/// @param newNumberOfHiddenNeurons Number of neurons in the hidden layer.
/// @param newNumberOfOutputs Number of output neurons in the neural network.

MultilayerPerceptron::MultilayerPerceptron
(int newNumberOfInputs, int newNumberOfHiddenNeurons, int newNumberOfOutputs)
{
   // Set network architecture
  
   Vector<int> newNumbersOfHiddenNeurons(1, newNumberOfHiddenNeurons);

   setNetworkArchitecture(newNumberOfInputs, newNumbersOfHiddenNeurons, newNumberOfOutputs);

   // Initialize number of independent parameters to 0

   numberOfIndependentParameters = 0;

   // Pre and post processing method

   preAndPostProcessingMethod = None;

   // Numerical differentiation method

   numericalDifferentiationMethod = CentralDifferences;

   // Display out of range warning
   
   displayOutOfRangeWarning = false;

   // Display messages
   
   display = true;

   // Numerical epsilon

   numericalEpsilon = 1.0e-6;  
}


MultilayerPerceptron::MultilayerPerceptron
(int newNumberOfInputs, Vector<int> newNumbersOfHiddenNeurons, int newNumberOfOutputs)
{
   // Set network architecture

   setNetworkArchitecture(newNumberOfInputs, newNumbersOfHiddenNeurons, newNumberOfOutputs);

   // Initialize number of independent parameters to 0

   numberOfIndependentParameters = 0;

   // Pre and post processing method

   preAndPostProcessingMethod = None;

   // Numerical differentiation method

   numericalDifferentiationMethod = CentralDifferences;

   // Display out of range warning
   
   displayOutOfRangeWarning = false;

   // Display messages
   
   display = true;

   // Numerical epsilon

   numericalEpsilon = 1.0e-6;
}


// DEFAULT CONSTRUCTOR

/// Default constructor. It creates a multilayer perceptron object with zero
/// inputs, zero neurons in the hidden layer and zero output neurons.
/// It also initializes the rest of class members to their default values:
///
/// <ul>
/// <li> Hidden layer activation function: Hyperbolic tangent.
/// <li> Output layer activation function: Linear.
/// <li> Number of independent parameters: 0.
/// <li> Pre and post-processing method: None. 
/// <li> Numerical epsilon method: Relative. 
/// <li> Numerical differentiation method: Central differences. 
/// <li> Display out of range warning: False.
/// <li> Display: True.
/// <li> Numerical epsilon: 1.0e-6.
/// </ul>

MultilayerPerceptron::MultilayerPerceptron(void)
{
   // Network architecture

   numberOfInputs = 0;

   numberOfOutputs = 0;

   // Activation functions

   outputLayerActivationFunction = Perceptron::Linear;

   // Pre and post processing method

   preAndPostProcessingMethod = None;

   // Numerical differentiation method

   numericalDifferentiationMethod = CentralDifferences;

   // Initialize number of independent parameters to 0

   numberOfIndependentParameters = 0;

   // Display out of range warning
   
   displayOutOfRangeWarning = false;

   // Display messages
   
   display = true;

   // Numerical epsilon

   numericalEpsilon = 1.0e-6;   
}


// DESTRUCTOR

/// Destructor.

MultilayerPerceptron::~MultilayerPerceptron(void)
{

}


// METHODS

// PreAndPostProcessingMethod getPreAndPostProcessingMethod(void) method

/// This method returns the method used for pre-processing network inputs and post-processing network outputs.
///
/// @see preprocessInput(Vector<double>).
/// @see postprocessOutputSignal(Vector<double>).
/// @see calculateOutput(Vector<double>).
/// @see calculateJacobian(Vector<double>).

MultilayerPerceptron::PreAndPostProcessingMethod MultilayerPerceptron::getPreAndPostProcessingMethod(void)
{
   return(preAndPostProcessingMethod);                           
}


// NumericalDifferentiationMethod getNumericalDifferentiationMethod(void) method

/// This method returns the method used for numerical differentiation in order to calculate the Jacobian matrix 
/// for the multilayer perceptron.
///
/// @see calculateJacobianTest(Vector<double>).

MultilayerPerceptron::NumericalDifferentiationMethod MultilayerPerceptron::getNumericalDifferentiationMethod(void)
{
   return(numericalDifferentiationMethod);                           
}


// int getNumberOfInputs(void) method

/// This method returns the number of inputs in the neural network.
///
/// @see getNumbersOfHiddenNeurons(void).
/// @see getNumberOfOutputs(void).

int MultilayerPerceptron::getNumberOfInputs(void)
{
   return(numberOfInputs);
}


// int getNumbersOfHiddenNeurons(void) method

/// This method returns the number of neurons in the hidden layer of the neural network.
///
/// @see getNumberOfInputs(void).
/// @see getNumberOfOutputs(void).

Vector<int> MultilayerPerceptron::getNumbersOfHiddenNeurons(void)
{
   return(numbersOfHiddenNeurons);
}


// int getNumberOfOutputs(void) method

/// This method returns the number of neurons in the output layer of the neural network.
///
/// @see getNumberOfInputs(void).
/// @see getNumberOfHiddenNeurons(void).

int MultilayerPerceptron::getNumberOfOutputs(void)
{
   return(numberOfOutputs);
}


// Perceptron::ActivationFunction getHiddenLayerActivationFunction(void) method

/// This method returns the activation function of the perceptron neurons of the hidden layer. 
///
/// @see getOutputLayerActivationFunction(void).

Vector<Perceptron::ActivationFunction> MultilayerPerceptron::getHiddenLayersActivationFunction(void)
{
   return(hiddenLayersActivationFunction);
}
   

// OutputLayerActivationFunction getOutputLayerActivationFunction(void) method

/// This method returns the activation function of the perceptron neurons of the output layer. 
///
/// @see getHiddenLayerActivationFunction(void).

Perceptron::ActivationFunction MultilayerPerceptron::getOutputLayerActivationFunction(void)
{
   return(outputLayerActivationFunction);
}


// Vector<Perceptron>& getHiddenLayers(void) method

/// This method returns the Vector of neurons in the hidden layer.
///
/// @see getOutputLayer(void).

Vector <Vector<Perceptron> >& MultilayerPerceptron::getHiddenLayers(void)
{
   return(hiddenLayers);
}


// Vector<Perceptron>& getOutputLayer(void) method

/// This method returns the Vector of output neurons. 
///
/// @see getHiddenLayer(void).

Vector<Perceptron>& MultilayerPerceptron::getOutputLayer(void)
{
   return(outputLayer);
}


// Vector<std::string> getNameOfInputVariables(void) method

/// This method returns the names of the input variables of the neural network. 
/// Such names are only used to give the user basic information about the problem at hand.
///
/// @see getNameOfOutputVariables(void).
/// @see getNameOfSingleInputVariable(int).
/// @see getNameOfSingleOutputVariable(int).

Vector<std::string> MultilayerPerceptron::getNameOfInputVariables(void)
{
   return(nameOfInputVariables);
}


// Vector<std::string> getNameOfOutputVariables(void) method

/// This method returns the names of the output variables of the neural network.
/// Such names are only used to give the user basic information about the problem at hand.
///
/// @see getNameOfInputVariables(void).
/// @see getNameOfSingleInputVariable(int).
/// @see getNameOfSingleOutputVariable(int).

Vector<std::string> MultilayerPerceptron::getNameOfOutputVariables(void)
{
   return(nameOfOutputVariables);
}


// std::string getNameOfSingleInputVariable(int) method

/// This method returns the name of a single input variable of the neural network. 
/// Such a name is only used to give the user basic information about the problem at hand.
///
/// @param i Index of input variable.
///
/// @see getNameOfInputVariables(void).
/// @see getNameOfOutputVariables(void).
/// @see getNameOfSingleOutputVariable(int).

std::string MultilayerPerceptron::getNameOfSingleInputVariable(int i)
{
   // Control sentence (if debug)

   #ifndef NDEBUG 

   if(i >= numberOfInputs)
   {
      std::cerr << std::endl
                << "Flood Error: MultilayerPerceptron class." << std::endl
                << "std::string getNameOfSingleInputVariable(int) method." << std::endl
                << "Index must and less than number of inputs."
                << std::endl << std::endl;

      exit(1);
   }

   #endif

   return(nameOfInputVariables[i]);
}


// std::string getNameOfSingleOutputVariable(int) method

/// This method returns the name of a single output variable of the neural network. 
/// Such a name is only used to give the user basic information about the problem at hand.
///
/// @param i Index of output variable.
///
/// @see getNameOfInputVariables(void).
/// @see getNameOfOutputVariables(void).
/// @see getNameOfSingleInputVariable(int).

std::string MultilayerPerceptron::getNameOfSingleOutputVariable(int i)
{
   // Control sentence (if debug)

   #ifndef NDEBUG 

   if(i >= numberOfOutputs)
   {
      std::cerr << std::endl
                << "Flood Error: MultilayerPerceptron class." << std::endl  
                << "std::string getNameOfSingleOutputVariable(int) method." << std::endl
                << "Index must less than number of outputs." << std::endl 
                << std::endl;

      exit(1);
   }

   #endif

   return(nameOfOutputVariables[i]);
}


// Vector<std::string> getUnitsOfInputVariables(void) method

/// This method returns the units of the input variables of the neural network as strings. 
/// Such units are only used to give the user basic information about the problem at hand.
///
/// @see getUnitsOfOutputVariables(void).
/// @see getUnitsOfSingleInputVariable(int).
/// @see getUnitsOfSingleOutputVariable(int).

Vector<std::string> MultilayerPerceptron::getUnitsOfInputVariables(void)
{
   return(unitsOfInputVariables);
}


// Vector<std::string> getUnitsOfOutputVariables(void) method

/// This method returns the units of the output variables of the neural network as strings. 
/// Such units are only used to give the user basic information about the problem at hand.
///
/// @see getUnitsOfInputVariables(void).
/// @see getUnitsOfSingleInputVariable(int).
/// @see getUnitsOfSingleOutputVariable(int).

Vector<std::string> MultilayerPerceptron::getUnitsOfOutputVariables(void)
{
   return(unitsOfOutputVariables);
}


// std::string getUnitsOfSingleInputVariable(int) method

/// This method returns the units of a single input variable as a string. 
/// Such units are only used to give the user basic information about the problem at hand.
///
/// @param i Index of input variable.
///
/// @see getUnitsOfInputVariables(void).
/// @see getUnitsOfOutputVariables(void).
/// @see getUnitsOfSingleOutputVariable(int).

std::string MultilayerPerceptron::getUnitsOfSingleInputVariable(int i)
{
   // Control sentence (if debug)

   #ifndef NDEBUG 

   if(i >= numberOfInputs)
   {
      std::cerr << std::endl
                << "Flood Error: MultilayerPerceptron class." << std::endl 
                << "std::string getUnitsOfSingleInputVariable(int) method." << std::endl
                << "Index must be less than number of inputs." << std::endl 
                << std::endl;

      exit(1);
   }

   #endif

   return(unitsOfInputVariables[i]);
}


// std::string getUnitsOfSingleOutputVariable(int) method

/// This method returns the units of a single output variable as a string. 
/// Such units are only used to give the user basic information about the problem at hand.
///
/// @param i Index of output variable.
///
/// @see getUnitsOfInputVariables(void).
/// @see getUnitsOfOutputVariables(void).
/// @see getUnitsOfSingleInputVariable(int).

std::string MultilayerPerceptron::getUnitsOfSingleOutputVariable(int i)
{
   // Control sentence (if debug)

   #ifndef NDEBUG 

   if(i >= numberOfOutputs)
   {
      std::cerr << std::endl
                << "Flood Error: MultilayerPerceptron class." << std::endl 
                << "std::string getUnitsOfSingleOutputVariable(int) method." << std::endl
                << "Index must be less than number of outputs." << std::endl 
                << std::endl;

      exit(1);
   }

   #endif

   return(unitsOfOutputVariables[i]);
}


// Vector<std::string> getDescriptionOfInputVariables(void) method

/// This method returns the description of the input variables of the neural network as strings. 
/// Such descriptions are only used to give the user basic information about the problem at hand.
///
/// @see getDescriptionOfOutputVariables(void).
/// @see getDescriptionOfSingleInputVariable(int).
/// @see getDescriptionOfSingleOutputVariable(int).

Vector<std::string> MultilayerPerceptron::getDescriptionOfInputVariables(void)
{
   return(descriptionOfInputVariables);
}


// Vector<std::string> getDescriptionOfOutputVariables(void) method

/// This method returns the description of the output variables of the neural network as strings. 
/// Such descriptions are only used to give the user basic information about the problem at hand.
///
/// @see getDescriptionOfInputVariables(void).
/// @see getDescriptionOfSingleInputVariable(int).
/// @see getDescriptionOfSingleOutputVariable(int).

Vector<std::string> MultilayerPerceptron::getDescriptionOfOutputVariables(void)
{
   return(descriptionOfOutputVariables);
}


// std::string getDescriptionOfSingleInputVariable(int) method

/// This method returns the description of a single input variable as a string. 
/// Such a description is only used to give the user basic information about the problem at hand.
///
/// @param i Index of input variable.
///
/// @see getDescriptionOfInputVariables(void).
/// @see getDescriptionOfOutputVariables(void).
/// @see getDescriptionOfSingleOutputVariable(int).

std::string MultilayerPerceptron::getDescriptionOfSingleInputVariable(int i)
{
   // Control sentence (if debug)

   #ifndef NDEBUG 

   if(i >= numberOfInputs)
   {
      std::cerr << std::endl
                << "Flood Error: MultilayerPerceptron class." << std::endl 
                << "std::string getDescriptionOfSingleInputVariable(int) method." << std::endl
                << "Index must be less than number of inputs." << std::endl 
                << std::endl;

      exit(1);
   }

   #endif

   return(descriptionOfInputVariables[i]);
}


// std::string getDescriptionOfSingleOutputVariable(int) method

/// This method returns the description of a single input variable as a string. 
/// Such a description is only used to give the user basic information about the problem at hand.
///
/// @param i Index of output variable.
///
/// @see getDescriptionOfInputVariables(void).
/// @see getDescriptionOfOutputVariables(void).
/// @see getDescriptionOfSingleInputVariable(int).

std::string MultilayerPerceptron::getDescriptionOfSingleOutputVariable(int i)
{
   // Control sentence (if debug)

   #ifndef NDEBUG 

   if(i >= numberOfOutputs)
   {
      std::cerr << std::endl
                << "Flood Error: MultilayerPerceptron class." << std::endl 
                << "std::string getDescriptionOfSingleOutputVariable(int) method." << std::endl
                << "Index must be less than number of outputs." << std::endl 
                << std::endl;

      exit(1);
   }

   #endif

   return(descriptionOfOutputVariables[i]);
}


// Vector< Vector<std::string> > getAllInformation(void) method

Vector< Vector<std::string> > MultilayerPerceptron::getAllInformation(void)
{
   Vector< Vector<std::string> > information(6);
 
   information[0] = nameOfInputVariables;
   information[1] = nameOfOutputVariables;

   information[2] = unitsOfInputVariables;
   information[3] = unitsOfOutputVariables;

   information[4] = descriptionOfInputVariables;
   information[5] = descriptionOfOutputVariables;

   return(information);
}


// Vector<double> getMeanOfInputVariables(void) method

/// This method returns the mean values of all the input variables of the neural network.
/// Such values are to be used for preprocessing inputs with the mean and standard deviation method. 
///
/// @see getStandardDeviationOfInputVariables(void).
/// @see getMeanOfOutputVariables(void).
/// @see getStandardDeviationOfOutputVariables(void).
/// @see getMeanAndStandardDeviationOfInputVariables(void).
/// @see getMeanAndStandardDeviationOfOutputVariables(void).
/// @see getMeanOfSingleInputVariable(int).
/// @see getStandardDeviationOfSingleInputVariable(int).
/// @see getMeanOfSingleOutputVariable(int).
/// @see getStandardDeviationOfSingleOutputVariable(int).
///
/// @see preprocessInput(Vector<double>).
/// @see postprocessOutputSignal(Vector<double>).

Vector<double> MultilayerPerceptron::getMeanOfInputVariables(void)
{
   return(meanOfInputVariables);
}


// Vector<double> getStandardDeviationOfInputVariables(void) method

/// This method returns the standard deviation values of all the input variables of the neural network.
/// Such values are to be used for preprocessing inputs with the mean and standard deviation method. 
///
/// @see getMeanOfInputVariables(void).
/// @see getMeanOfOutputVariables(void).
/// @see getStandardDeviationOfOutputVariables(void).
/// @see getMeanAndStandardDeviationOfInputVariables(void).
/// @see getMeanAndStandardDeviationOfOutputVariables(void).
/// @see getMeanOfSingleInputVariable(int).
/// @see getStandardDeviationOfSingleInputVariable(int).
/// @see getMeanOfSingleOutputVariable(int).
/// @see getStandardDeviationOfSingleOutputVariable(int).
///
/// @see preprocessInput(Vector<double>).
/// @see postprocessOutputSignal(Vector<double>).

Vector<double> MultilayerPerceptron::getStandardDeviationOfInputVariables(void)
{
   return(standardDeviationOfInputVariables);
}


// Matrix<double> getMeanAndStandardDeviationOfInputVariables(void) method

/// This method returns the mean and the standard deviation values of all the input variables in a single matrix. 
/// The first row contains the mean values of the input variables.
/// The second row contains the standard deviation values of the input variables.
/// Such values are to be used for preprocessing inputs with the mean and standard deviation method. 
///
/// @see getMeanOfInputVariables(void).
/// @see getStandardDeviationOfInputVariables(void).
/// @see getMeanOfOutputVariables(void).
/// @see getStandardDeviationOfOutputVariables(void).
/// @see getMeanAndStandardDeviationOfOutputVariables(void).
/// @see getMeanOfSingleInputVariable(int).
/// @see getStandardDeviationOfSingleInputVariable(int).
/// @see getMeanOfSingleOutputVariable(int).
/// @see getStandardDeviationOfSingleOutputVariable(int).
///
/// @see preprocessInput(Vector<double>).
/// @see postprocessOutputSignal(Vector<double>).

Matrix<double> MultilayerPerceptron::getMeanAndStandardDeviationOfInputVariables(void)
{
   Matrix<double> meanAndStandardDeviationOfInputVariables(2, numberOfInputs);

   meanAndStandardDeviationOfInputVariables.setRow(0, meanOfInputVariables);
   meanAndStandardDeviationOfInputVariables.setRow(1, standardDeviationOfInputVariables);

   return(meanAndStandardDeviationOfInputVariables);
}


// Vector<double> getMeanOfOutputVariables(void)

/// This method returns the mean values of all the output variables of the neural network.
/// Such values are to be used for postprocessing output signals with the mean and standard deviation method. 
///
/// @see getMeanOfInputVariables(void).
/// @see getStandardDeviationOfInputVariables(void).
/// @see getStandardDeviationOfOutputVariables(void).
/// @see getMeanAndStandardDeviationOfInputVariables(void).
/// @see getMeanAndStandardDeviationOfOutputVariables(void).
/// @see getMeanOfSingleInputVariable(int).
/// @see getStandardDeviationOfSingleInputVariable(int).
/// @see getMeanOfSingleOutputVariable(int).
/// @see getStandardDeviationOfSingleOutputVariable(int).
///
/// @see preprocessInput(Vector<double>).
/// @see postprocessOutputSignal(Vector<double>).

Vector<double> MultilayerPerceptron::getMeanOfOutputVariables(void)
{
   return(meanOfOutputVariables);
}


// Vector<double> getStandardDeviationOfOutputVariables(void)

/// This method returns the standard deviation values of all the output variables of the neural network.
/// Such values are to be used for postprocessing output signals with the mean and standard deviation method. 
///
/// @see getMeanOfInputVariables(void).
/// @see getStandardDeviationOfInputVariables(void).
/// @see getMeanOfOutputVariables(void).
/// @see getMeanAndStandardDeviationOfInputVariables(void).
/// @see getMeanAndStandardDeviationOfOutputVariables(void).
/// @see getMeanOfSingleInputVariable(int).
/// @see getStandardDeviationOfSingleInputVariable(int).
/// @see getMeanOfSingleOutputVariable(int).
/// @see getStandardDeviationOfSingleOutputVariable(int).
///
/// @see preprocessInput(Vector<double>).
/// @see postprocessOutputSignal(Vector<double>).

Vector<double> MultilayerPerceptron::getStandardDeviationOfOutputVariables(void)
{
   return(standardDeviationOfOutputVariables);
}


// Matrix<double> getMeanAndStandardDeviationOfOutputVariables(void) method

/// This method returns the mean and the standard deviation values of all the output variables in a single matrix. 
/// The first row contains the mean values of the output variables.
/// The second row contains the standard deviation values of the output variables.
/// Such values are to be used for postprocessing output signals with the mean and standard deviation method. 
///
/// @see getMeanOfInputVariables(void).
/// @see getStandardDeviationOfInputVariables(void).
/// @see getMeanOfOutputVariables(void).
/// @see getStandardDeviationOfOutputVariables(void).
/// @see getMeanAndStandardDeviationOfInputVariables(void).
/// @see getMeanOfSingleInputVariable(int).
/// @see getStandardDeviationOfSingleInputVariable(int).
/// @see getMeanOfSingleOutputVariable(int).
/// @see getStandardDeviationOfSingleOutputVariable(int).
///
/// @see preprocessInput(Vector<double>).
/// @see postprocessOutputSignal(Vector<double>).

Matrix<double> MultilayerPerceptron::getMeanAndStandardDeviationOfOutputVariables(void)
{
   Matrix<double> meanAndStandardDeviationOfOutputVariables(2, numberOfOutputs);

   meanAndStandardDeviationOfOutputVariables.setRow(0, meanOfOutputVariables);       
   meanAndStandardDeviationOfOutputVariables.setRow(1, standardDeviationOfOutputVariables);

   return(meanAndStandardDeviationOfOutputVariables);
}


// double getMeanOfSingleInputVariable(int) method

/// This method returns the mean value of a single input variable of the neural network.
/// Such a value is to be used for preprocessing that input with the mean and standard deviation method. 
///
/// @param i  Index of input variable.
///
/// @see getMeanOfInputVariables(void).
/// @see getStandardDeviationOfInputVariables(void).
/// @see getMeanOfOutputVariables(void).
/// @see getStandardDeviationOfOutputVariables(void).
/// @see getMeanAndStandardDeviationOfInputVariables(void).
/// @see getMeanAndStandardDeviationOfOutputVariables(void).
/// @see getStandardDeviationOfSingleInputVariable(int).
/// @see getMeanOfSingleOutputVariable(int).
/// @see getStandardDeviationOfSingleOutputVariable(int).
///
/// @see preprocessInput(Vector<double>).
/// @see postprocessOutputSignal(Vector<double>).

double MultilayerPerceptron::getMeanOfSingleInputVariable(int i)
{
   // Control sentence (if debug)

   #ifndef NDEBUG 

   if(numberOfInputs == 0)
   {
      std::cerr << std::endl
                << "Flood Error: MultilayerPerceptron class." << std::endl 
                << "double getMeanOfSingleInputVariable(int) method." << std::endl
                << "Number of inputs is zero." << std::endl 
                << std::endl;

      exit(1);   
   }
   else if(i >= numberOfInputs)
   {
      std::cerr << std::endl
                << "Flood Error: MultilayerPerceptron class." << std::endl 
                << "double getMeanOfSingleInputVariable(int) method." << std::endl
                << "Index must be less than number of inputs." << std::endl 
                << std::endl;

      exit(1);
   }

   #endif

   return(meanOfInputVariables[i]);
}


// double getStandardDeviationOfSingleInputVariable(int) method

/// This method returns the standard deviation value of a single input variable of the neural network.
/// Such a value is to be used for preprocessing that input with the mean and standard deviation method. 
///
/// @param i Index of input variable.
///
/// @see getMeanOfInputVariables(void).
/// @see getStandardDeviationOfInputVariables(void).
/// @see getMeanOfOutputVariables(void).
/// @see getStandardDeviationOfOutputVariables(void).
/// @see getMeanAndStandardDeviationOfInputVariables(void).
/// @see getMeanAndStandardDeviationOfOutputVariables(void).
/// @see getMeanOfSingleInputVariable(int).
/// @see getMeanOfSingleOutputVariable(int).
/// @see getStandardDeviationOfSingleOutputVariable(int).
///
/// @see preprocessInput(Vector<double>).
/// @see postprocessOutputSignal(Vector<double>).

double MultilayerPerceptron::getStandardDeviationOfSingleInputVariable(int i)
{
   // Control sentence (if debug)

   #ifndef NDEBUG 

   if(numberOfInputs == 0)
   {
      std::cerr << std::endl
                << "Flood Error: MultilayerPerceptron class." << std::endl 
                << "double getStandardDeviationOfSingleInputVariable(int) method." << std::endl
                << "Number of inputs is zero." << std::endl 
                << std::endl;

      exit(1);   
   }
   else if(i >= numberOfInputs)
   {
      std::cerr << std::endl
                << "Flood Error: MultilayerPerceptron class." << std::endl 
                << "double getStandardDeviationOfSingleInputVariable(int) method." << std::endl
                << "Index must be less than number of inputs." << std::endl 
                << std::endl;

      exit(1);
   }

   #endif

   return(standardDeviationOfInputVariables[i]);
}


// double getMeanOfSingleOutputVariable(int) method

/// This method returns the mean values of a single input variable of the neural network.
/// Such a values is to be used for preprocessing inputs and postprocessing output signals with the mean and 
/// standard deviation method. 
///
/// @param i Index of output variable.
///
/// @see getStandardDeviationOfInputVariables(void).
/// @see getMeanOfOutputVariables(void).
/// @see getStandardDeviationOfOutputVariables(void).
/// @see getMeanAndStandardDeviationOfInputVariables(void).
/// @see getMeanAndStandardDeviationOfOutputVariables(void).
/// @see getMeanOfSingleInputVariable(int).
/// @see getStandardDeviationOfSingleInputVariable(int).
/// @see getMeanOfSingleOutputVariable(int).
/// @see getStandardDeviationOfSingleOutputVariable(int).
///
/// @see preprocessInput(Vector<double>).
/// @see postprocessOutputSignal(Vector<double>).

double MultilayerPerceptron::getMeanOfSingleOutputVariable(int i)
{
   // Control sentence (if debug)

   #ifndef NDEBUG 

   if(numberOfOutputs == 0)
   {
      std::cerr << std::endl
                << "Flood Error: MultilayerPerceptron class." << std::endl 
                << "double getMeanOfSingleOutputVariable(int) method." << std::endl
                << "Number of outputs is zero." << std::endl 
                << std::endl;

      exit(1);   
   }
   else if(i >= numberOfOutputs)
   {
      std::cerr << std::endl
                << "Flood Error: MultilayerPerceptron class." << std::endl 
                << "double getMeanOfSingleOutputVariable(int) method." << std::endl
                << "Index must be less than number of outputs." << std::endl 
                << std::endl;

      exit(1);
   }

   #endif

   return(meanOfOutputVariables[i]);
}


// double getStandardDeviationOfSingleOutputVariable(int) method

/// This method returns the standard deviation value of a single output variable of the neural network.
/// Such a value is to be used for postprocessing output signals with the mean and standard deviation method. 
///
/// @param i Index of output variable.
///
/// @see getMeanOfInputVariables(void).
/// @see getStandardDeviationOfInputVariables(void).
/// @see getMeanOfOutputVariables(void).
/// @see getStandardDeviationOfOutputVariables(void).
/// @see getMeanAndStandardDeviationOfInputVariables(void).
/// @see getMeanAndStandardDeviationOfOutputVariables(void).
/// @see getMeanOfSingleInputVariable(int).
/// @see getStandardDeviationOfSingleInputVariable(int).
/// @see getMeanOfSingleOutputVariable(int).
///
/// @see preprocessInput(Vector<double>).
/// @see postprocessOutputSignal(Vector<double>).

double MultilayerPerceptron::getStandardDeviationOfSingleOutputVariable(int i)
{
   // Control sentence (if debug)

   #ifndef NDEBUG 

   if(numberOfOutputs == 0)
   {
      std::cerr << std::endl
                << "Flood Error: MultilayerPerceptron class." << std::endl 
                << "double getStandardDeviationOfSingleOutputVariable(int) method." << std::endl
                << "Number of outputs is zero." << std::endl 
                << std::endl;

      exit(1);   
   }
   else if(i >= numberOfOutputs)
   {
      std::cerr << std::endl
                << "Flood Error: MultilayerPerceptron class." << std::endl 
                << "double getStandardDeviationOfSingleOutputVariable(int) method." << std::endl
                << "Index must be less than number of output." << std::endl 
                << std::endl;

      exit(1);
   }

   #endif

   return(standardDeviationOfOutputVariables[i]);
}


// Vector<double> getMinimumOfInputVariables(void) method

/// This method returns the minimum values of all the input variables in the neural network.
/// Such values are to be used for preprocessing inputs with the minimum and maximum method.
///
/// @see getMaximumOfInputVariables(void).
/// @see getMinimumOfOutputVariables(void).
/// @see getMaximumOfOutputVariables(void).
/// @see getMinimumAndMaximumOfInputVariables(void).
/// @see getMinimumAndMaximumOfOutputVariables(void).
/// @see getMinimumOfSingleInputVariable(int).
/// @see getMaximumOfSingleInputVariable(int).
/// @see getMinimumOfSingleOutputVariable(int).
/// @see getMaximumOfSingleOutputVariable(int).
///
/// @see preprocessInput(Vector<double>).
/// @see postprocessOutputSignal(Vector<double>).

Vector<double> MultilayerPerceptron::getMinimumOfInputVariables(void)
{
   return(minimumOfInputVariables);               
}


// Vector<double> getMaximumOfInputVariables(void) method

/// This method returns the maximum values of all the input variables in the neural network.
/// Such values are to be used for preprocessing inputs with the minimum and maximum method.
///
/// @see getMinimumOfInputVariables(void).
/// @see getMinimumOfOutputVariables(void).
/// @see getMaximumOfOutputVariables(void).
/// @see getMinimumAndMaximumOfInputVariables(void).
/// @see getMinimumAndMaximumOfOutputVariables(void).
/// @see getMinimumOfSingleInputVariable(int).
/// @see getMaximumOfSingleInputVariable(int).
/// @see getMinimumOfSingleOutputVariable(int).
/// @see getMaximumOfSingleOutputVariable(int).
///
/// @see preprocessInput(Vector<double>).
/// @see postprocessOutputSignal(Vector<double>).

Vector<double> MultilayerPerceptron::getMaximumOfInputVariables(void)
{
   return(maximumOfInputVariables);               
}


// Vector<double> getMinimumOfOutputVariables(void) method

/// This method returns the minimum values of all the output variables in the neural network.
/// Such values are to be used for postprocessing output signals with the minimum and maximum method. 
///
/// @see getMinimumOfInputVariables(void).
/// @see getMaximumOfInputVariables(void).
/// @see getMaximumOfOutputVariables(void).
/// @see getMinimumAndMaximumOfInputVariables(void).
/// @see getMinimumAndMaximumOfOutputVariables(void).
/// @see getMinimumOfSingleInputVariable(int).
/// @see getMaximumOfSingleInputVariable(int).
/// @see getMinimumOfSingleOutputVariable(int).
/// @see getMaximumOfSingleOutputVariable(int).
///
/// @see preprocessInput(Vector<double>).
/// @see postprocessOutputSignal(Vector<double>).

Vector<double> MultilayerPerceptron::getMinimumOfOutputVariables(void)
{
   return(minimumOfOutputVariables);               
}


// Vector<double> getMaximumOfOutputVariables(void) method

/// This method returns the maximum values of all the output variables in the neural network.
/// Such values are to be used for postprocessing output signals with the minimum and maximum method. 
///
/// @see getMinimumOfInputVariables(void).
/// @see getMaximumOfInputVariables(void).
/// @see getMinimumOfOutputVariables(void).
/// @see getMinimumAndMaximumOfInputVariables(void).
/// @see getMinimumAndMaximumOfOutputVariables(void).
/// @see getMinimumOfSingleInputVariable(int).
/// @see getMaximumOfSingleInputVariable(int).
/// @see getMinimumOfSingleOutputVariable(int).
/// @see getMaximumOfSingleOutputVariable(int).
///
/// @see preprocessInput(Vector<double>).
/// @see postprocessOutputSignal(Vector<double>).

Vector<double> MultilayerPerceptron::getMaximumOfOutputVariables(void)
{
   return(maximumOfOutputVariables);               
}


// Matrix<double> getMinimumAndMaximumOfInputVariables(void) method

/// This method returns the minimum and the maximum values of all the input variables in a single matrix. 
/// The first row contains the minimum values of the input variables.
/// The second row contains the maximum values of the input variables.
/// Such values are to be used for preprocessing inputs with the minimum and maximum method.
///
/// @see getMinimumOfInputVariables(void).
/// @see getMaximumOfInputVariables(void).
/// @see getMinimumOfOutputVariables(void).
/// @see getMaximumOfOutputVariables(void).
/// @see getMinimumAndMaximumOfOutputVariables(void).
/// @see getMinimumOfSingleInputVariable(int).
/// @see getMaximumOfSingleInputVariable(int).
/// @see getMinimumOfSingleOutputVariable(int).
/// @see getMaximumOfSingleOutputVariable(int).
///
/// @see preprocessInput(Vector<double>).
/// @see postprocessOutputSignal(Vector<double>).

Matrix<double> MultilayerPerceptron::getMinimumAndMaximumOfInputVariables(void)
{
   Matrix<double> minimumAndMaximumOfInputVariables(2, numberOfInputs);

   minimumAndMaximumOfInputVariables.setRow(0, minimumOfInputVariables);
   minimumAndMaximumOfInputVariables.setRow(1, maximumOfInputVariables);

   return(minimumAndMaximumOfInputVariables);
}


// Matrix<double> getMinimumAndMaximumOfOutputVariables(void) method

/// This method returns the minimum and the maximum values of all the output variables in a single matrix. 
/// The first row contains the minimum values of the output variables.
/// The second row contains the maximum values of the output variables.
/// Such values are to be used for postprocessing output signals with the minimum and maximum method. 
///
/// @see getMinimumOfInputVariables(void).
/// @see getMaximumOfInputVariables(void).
/// @see getMinimumOfOutputVariables(void).
/// @see getMaximumOfOutputVariables(void).
/// @see getMinimumAndMaximumOfInputVariables(void).
/// @see getMinimumOfSingleInputVariable(int).
/// @see getMaximumOfSingleInputVariable(int).
/// @see getMinimumOfSingleOutputVariable(int).
/// @see getMaximumOfSingleOutputVariable(int).
///
/// @see preprocessInput(Vector<double>).
/// @see postprocessOutputSignal(Vector<double>).

Matrix<double> MultilayerPerceptron::getMinimumAndMaximumOfOutputVariables(void)
{
   Matrix<double> minimumAndMaximumOfOutputVariables(2, numberOfOutputs);

   minimumAndMaximumOfOutputVariables.setRow(0, minimumOfOutputVariables);
   minimumAndMaximumOfOutputVariables.setRow(1, maximumOfOutputVariables);

   return(minimumAndMaximumOfOutputVariables);
}


// double getMinimumOfSingleInputVariable(int) method

/// This method returns the minimum value of a single input variable in the neural network.
/// Such value is to be used for preprocessing that input with the minimum and maximum method.
///
/// @see getMinimumOfInputVariables(void).
/// @see getMaximumOfInputVariables(void).
/// @see getMinimumOfOutputVariables(void).
/// @see getMaximumOfOutputVariables(void).
/// @see getMinimumAndMaximumOfInputVariables(void).
/// @see getMinimumAndMaximumOfOutputVariables(void).
/// @see getMaximumOfSingleInputVariable(int).
/// @see getMinimumOfSingleOutputVariable(int).
/// @see getMaximumOfSingleOutputVariable(int).
///
/// @see preprocessInput(Vector<double>).
/// @see postprocessOutputSignal(Vector<double>).

double MultilayerPerceptron::getMinimumOfSingleInputVariable(int i)
{
   // Control sentence (if debug)

   #ifndef NDEBUG 

   if(numberOfInputs == 0)
   {
      std::cerr << std::endl
                << "Flood Error: MultilayerPerceptron class." << std::endl 
                << "double getMinimumOfSingleInputVariable(int) method." << std::endl
                << "Number of inputs is zero." << std::endl 
                << std::endl;

      exit(1);   
   }
   else if(i >= numberOfInputs)
   {
      std::cerr << std::endl
                << "Flood Error: MultilayerPerceptron class." << std::endl 
                << "double getMinimumOfSingleInputVariable(int) method." << std::endl
                << "Index must be less than number of inputs." << std::endl 
                << std::endl;

      exit(1);
   }

   #endif

   return(minimumOfInputVariables[i]);
}


// double getMaximumOfSingleInputVariable(int) method

/// This method returns the maximum value of a single input variable in the neural network.
/// Such value is to be used for preprocessing that input with the minimum and maximum method.
///
/// @see getMinimumOfInputVariables(void).
/// @see getMaximumOfInputVariables(void).
/// @see getMinimumOfOutputVariables(void).
/// @see getMaximumOfOutputVariables(void).
/// @see getMinimumAndMaximumOfInputVariables(void).
/// @see getMinimumAndMaximumOfOutputVariables(void).
/// @see getMinimumOfSingleInputVariable(int).
/// @see getMinimumOfSingleOutputVariable(int).
/// @see getMaximumOfSingleOutputVariable(int).
///
/// @see preprocessInput(Vector<double>).
/// @see postprocessOutputSignal(Vector<double>).

double MultilayerPerceptron::getMaximumOfSingleInputVariable(int i)
{
   // Control sentence (if debug)

   #ifndef NDEBUG 

   if(numberOfInputs == 0)
   {
      std::cerr << std::endl
                << "Flood Error: MultilayerPerceptron class." << std::endl 
                << "double getMaximumOfSingleInputVariable(int) method." << std::endl
                << "Number of inputs is zero." << std::endl 
                << std::endl;

      exit(1);   
   }
   else if(i >= numberOfInputs)
   {
      std::cerr << std::endl
                << "Flood Error: MultilayerPerceptron class." << std::endl 
                << "double getMaximumOfSingleInputVariable(int) method." << std::endl
                << "Index must be less than number of inputs." << std::endl 
                << std::endl;

      exit(1);
   }

   #endif

   return(maximumOfInputVariables[i]);
}


// double getMinimumOfSingleOutputVariable(int) method

/// This method returns the minimum value of a single output variable in the neural network.
/// Such value is to be used for postprocessing that output signal with the minimum and maximum method.
///
/// @see getMinimumOfInputVariables(void).
/// @see getMaximumOfInputVariables(void).
/// @see getMinimumOfOutputVariables(void).
/// @see getMaximumOfOutputVariables(void).
/// @see getMinimumAndMaximumOfInputVariables(void).
/// @see getMinimumAndMaximumOfOutputVariables(void).
/// @see getMaximumOfSingleInputVariable(int).
/// @see getMaximumOfSingleOutputVariable(int).
///
/// @see preprocessInput(Vector<double>).
/// @see postprocessOutputSignal(Vector<double>).

double MultilayerPerceptron::getMinimumOfSingleOutputVariable(int i)
{
   // Control sentence (if debug)

   #ifndef NDEBUG 

   if(numberOfOutputs == 0)
   {
      std::cerr << std::endl
                << "Flood Error: MultilayerPerceptron class." << std::endl 
                << "double getMinimumOfSingleOutputVariable(int) method." << std::endl
                << "Number of outputs is zero." << std::endl 
                << std::endl;

      exit(1);   
   }
   else if(i >= numberOfOutputs)
   {
      std::cerr << std::endl
                << "Flood Error: MultilayerPerceptron class." << std::endl 
                << "double getMinimumOfSingleOutputVariable(int) method." << std::endl
                << "Index must be less than number of outputs." << std::endl 
                << std::endl;

      exit(1);
   }

   #endif

   return(minimumOfOutputVariables[i]);
}


// double getMaximumOfSingleOutputVariable(int) method

/// This method returns the maximum value of a single input variable in the neural network.
/// Such value is to be used for postprocessing that output signal with the minimum and maximum method.
///
/// @see getMinimumOfInputVariables(void).
/// @see getMaximumOfInputVariables(void).
/// @see getMinimumOfOutputVariables(void).
/// @see getMaximumOfOutputVariables(void).
/// @see getMinimumAndMaximumOfInputVariables(void).
/// @see getMinimumAndMaximumOfOutputVariables(void).
/// @see getMinimumOfSingleInputVariable(int).
/// @see getMaximumOfSingleInputVariable(int).
/// @see getMinimumOfSingleOutputVariable(int).
///
/// @see preprocessInput(Vector<double>).
/// @see postprocessOutputSignal(Vector<double>).

double MultilayerPerceptron::getMaximumOfSingleOutputVariable(int i)
{
   // Control sentence (if debug)

   #ifndef NDEBUG 

   if(numberOfOutputs == 0)
   {
      std::cerr << std::endl
                << "Flood Error: MultilayerPerceptron class." << std::endl 
                << "double getMaximumOfSingleOutputVariable(int) method." << std::endl
                << "Number of outputs is zero." << std::endl 
                << std::endl;

      exit(1);   
   }
   else if(i >= numberOfOutputs)
   {
      std::cerr << std::endl
                << "Flood Error: MultilayerPerceptron class." << std::endl 
                << "double getMaximumOfSingleOutputVariable(int) method." << std::endl
                << "Index must be less than number of outputs." << std::endl 
                << std::endl;

      exit(1);
   }

   #endif

   return(maximumOfOutputVariables[i]);
}


// Vector< Vector<double> > getAllStatistics(void) method

Vector< Vector<double> > MultilayerPerceptron::getAllStatistics(void)
{
   Vector< Vector<double> > statistics(8);

   statistics[0] = meanOfInputVariables;
   statistics[1] = standardDeviationOfInputVariables;

   statistics[2] = meanOfOutputVariables;
   statistics[3] = standardDeviationOfOutputVariables;

   statistics[4] = minimumOfInputVariables;
   statistics[5] = maximumOfInputVariables;

   statistics[6] = minimumOfOutputVariables;
   statistics[7] = maximumOfOutputVariables;

   return(statistics);
}

// Vector<double> getLowerBoundOfOutputVariables(void) method

/// This method returns the lower bound of all the output variables in the neural network.
/// These values are used to postprocess the outputs so that they are not less than the lower bounds. 
///
/// @see getUpperBoundOfOutputVariables(void).
/// @see getLowerAndUpperBoundsOfOutputVariables(void).
/// @see getLowerBoundOfSingleOutputVariable(void).
/// @see getUpperBoundOfSingleOutputVariable(void).
///
/// @see applyLowerAndUpperBounds(Vector<double>).

Vector<double> MultilayerPerceptron::getLowerBoundOfOutputVariables(void)
{
   return(lowerBoundOfOutputVariables);               
}


// Vector<double> getUpperBoundOfOutputVariables(void) method

/// This method returns the upper bound of all the output variables in the neural network.
/// These values are used to postprocess the outputs so that they are not greater than the upper bounds. 
///
/// @see getLowerBoundOfOutputVariables(void).
/// @see getLowerAndUpperBoundsOfOutputVariables(void).
/// @see getLowerBoundOfSingleOutputVariable(void).
/// @see getUpperBoundOfSingleOutputVariable(void).
///
/// @see applyLowerAndUpperBounds(Vector<double>).

Vector<double> MultilayerPerceptron::getUpperBoundOfOutputVariables(void)
{
   return(upperBoundOfOutputVariables);               
}


// Matrix<double> getLowerAndUpperBoundsOfOutputVariables(void) method

/// This method returns the lower bounds and the upper bounds of all the output
/// variables in a single matrix. 
/// The first row contains the lower bound values of the output variables.
/// The second row contains the upper bound values of the output variables.
/// These values are used to postprocess the outputs so that they are neither less than the lower bounds nor 
/// greater than the upper bounds.  
///
/// @see getLowerBoundOfOutputVariables(void).
/// @see getUpperBoundOfOutputVariables(void).
/// @see getLowerBoundOfSingleOutputVariable(void).
/// @see getUpperBoundOfSingleOutputVariable(void).
///
/// @see applyLowerAndUpperBounds(Vector<double>).

Matrix<double> MultilayerPerceptron::getLowerAndUpperBoundsOfOutputVariables(void)
{
   Matrix<double> lowerAndUpperBoundsOfOutputVariables(2, numberOfOutputs);

   lowerAndUpperBoundsOfOutputVariables.setRow(0, lowerBoundOfOutputVariables);
   lowerAndUpperBoundsOfOutputVariables.setRow(1, upperBoundOfOutputVariables);

   return(lowerAndUpperBoundsOfOutputVariables);
}


// double getLowerBoundOfSingleOutputVariable(int) method

/// This method returns the lower bound of a single output variable in the neural network.
/// This value is used to postprocess that output so that it is not less than the lower bound. 
///
/// @see getLowerBoundOfOutputVariables(void).
/// @see getUpperBoundOfOutputVariables(void).
/// @see getLowerAndUpperBoundsOfOutputVariables(void).
/// @see getUpperBoundOfSingleOutputVariable(void).
///
/// @see applyLowerAndUpperBounds(Vector<double>).

double MultilayerPerceptron::getLowerBoundOfSingleOutputVariable(int i)
{
   // Control sentence (if debug)

   #ifndef NDEBUG 

   if(numberOfOutputs == 0)
   {
      std::cerr << std::endl
                << "Flood Error: MultilayerPerceptron class." << std::endl 
                << "double getLowerBoundOfSingleOutputVariable(int) method." << std::endl
                << "Number of outputs is zero." << std::endl 
                << std::endl;

      exit(1);   
   }
   else if(i >= numberOfOutputs)
   {
      std::cerr << std::endl
                << "Flood Error: MultilayerPerceptron class." << std::endl 
                << "double getLowerBoundOfSingleOutputVariable(int) method." << std::endl
                << "Index must be less than number of output." << std::endl 
                << std::endl;

      exit(1);
   }

   #endif

   return(lowerBoundOfOutputVariables[i]);
}


// double getUpperBoundOfSingleOutputVariable(int) method

/// This method returns the upper bound of a single output variable in the neural network.
/// This value is used to postprocess that output so that it is not greater than the upper bound. 
///
/// @see getLowerBoundOfOutputVariables(void).
/// @see getUpperBoundOfOutputVariables(void).
/// @see getLowerAndUpperBoundsOfOutputVariables(void).
/// @see getLowerBoundOfSingleOutputVariable(void).
///
/// @see applyLowerAndUpperBounds(Vector<double>).

double MultilayerPerceptron::getUpperBoundOfSingleOutputVariable(int i)
{
   // Control sentence (if debug)

   #ifndef NDEBUG 

   if(numberOfOutputs == 0)
   {
      std::cerr << std::endl
                << "Flood Error: MultilayerPerceptron class." << std::endl 
                << "double getUpperBoundOfSingleOutputVariable(int) method." << std::endl
                << "Number of outputs is zero." << std::endl 
                << std::endl;

      exit(1);   
   }
   else if(i >= numberOfOutputs)
   {
      std::cerr << std::endl
                << "Flood Error: MultilayerPerceptron class." << std::endl 
                << "double getUpperOfSingleOutputVariable(int) method." << std::endl
                << "Index must be less than number of outputs." << std::endl 
                << std::endl;

      exit(1);
   }

   #endif

   return(upperBoundOfOutputVariables[i]);
}


// int getNumberOfIndependentParameters(void) method

/// This method returns the number of free parameters independent of the neural network.
/// Independent parameters can be used in the context of neural netwotks for many purposes.
///
/// @see CarProblem.
/// @see CarProblemNeurocomputing.

int MultilayerPerceptron::getNumberOfIndependentParameters(void)
{
   return(numberOfIndependentParameters);    
}


// Vector<std::string> getNameOfIndependentParameters(void) method

/// This method returns the name of the independent parameters. 
/// Such names are only used to give the user basic information about the problem at hand.
///
/// @see getNameOfSingleIndependentParameter(int).

Vector<std::string> MultilayerPerceptron::getNameOfIndependentParameters(void)
{
   return(nameOfIndependentParameters);    
}


// std::string getNameOfSingleIndependentParameter(int) method

/// This method returns the name of a single independent parameter. 
/// Such name is only used to give the user basic information about the problem at hand.
///
/// @param i Index of independent parameter.
///
/// @see getNameOfIndependentParameters(void).

std::string MultilayerPerceptron::getNameOfSingleIndependentParameter(int i)
{
   // Control sentence (if debug)

   #ifndef NDEBUG 

   if(numberOfIndependentParameters == 0)
   {
      std::cerr << std::endl
                << "Flood Error: MultilayerPerceptron class." << std::endl 
                << "std::string getNameOfSingleIndependentParameter(int) method." << std::endl
                << "Number of independent parameters is zero." << std::endl 
                << std::endl;

      exit(1);
   
   }
   else if(i >= numberOfIndependentParameters)
   {
      std::cerr << std::endl
                << "Flood Error: MultilayerPerceptron class." << std::endl 
                << "std::string getNameOfIndependentParameters(int) method." << std::endl
                << "Index must be less than number of independent parameters." << std::endl 
                << std::endl;

      exit(1);
   }

   #endif

   return(nameOfIndependentParameters[i]);
}


// Vector<std::string> getUnitsOfIndependentParameters(void) method

/// This method returns the units of the independent parameters. 
/// Such units are only used to give the user basic information about the problem at hand.
///
/// @see getUnitsOfSingleIndependentParameter(int).

Vector<std::string> MultilayerPerceptron::getUnitsOfIndependentParameters(void)
{
   return(unitsOfIndependentParameters);
}


// std::string getUnitsOfSingleIndependentParameter(int) method

/// This method returns the units of a single independent parameter. 
/// Such units are only used to give the user basic information about the problem at hand.
///
/// @param i Index of independent parameter.
///
/// @see getUnitsOfIndependentParameters(void).

std::string MultilayerPerceptron::getUnitsOfSingleIndependentParameter(int i)
{
   // Control sentence (if debug)

   #ifndef NDEBUG 

   if(numberOfIndependentParameters == 0)
   {
      std::cerr << std::endl
                << "Flood Error: MultilayerPerceptron class." << std::endl 
                << "std::string getUnitsOfSingleIndependentParameter(int) method." << std::endl
                << "Number of independent parameters is zero." << std::endl 
                << std::endl;

      exit(1);
   
   }
   else if(i >= numberOfIndependentParameters)
   {
      std::cerr << std::endl
                << "Flood Error: MultilayerPerceptron class." << std::endl 
                << "std::string getUnitsOfIndependentParameters(int) method." << std::endl
                << "Index must be less than number of independent parameters." << std::endl 
                << std::endl;

      exit(1);
   }

   #endif

   return(unitsOfIndependentParameters[i]);
}


// Vector<std::string> getDescriptionOfIndependentParameters(void) method

/// This method returns the description of the independent parameters. 
/// Such descriptions are only used to give the user basic information about the problem at hand.
///
/// @see getDescriptionOfSingleIndependentParameter(int).


Vector<std::string> MultilayerPerceptron::getDescriptionOfIndependentParameters(void)
{
   return(descriptionOfIndependentParameters);
}


// std::string getDescriptionOfSingleIndependentParameter(int) method

/// This method returns the description of a single independent parameter. 
/// Such description is only used to give the user basic information about the problem at hand.
///
/// @param i Index of independent parameter.
///
/// @see getDescriptionOfIndependentParameters(void).

std::string MultilayerPerceptron::getDescriptionOfSingleIndependentParameter(int i)
{
   // Control sentence (if debug)

   #ifndef NDEBUG 

   if(numberOfIndependentParameters == 0)
   {
      std::cerr << std::endl
                << "Flood Error: MultilayerPerceptron class." << std::endl 
                << "std::string getDescriptionOfSingleIndependentParameter(int) method." << std::endl
                << "Number of independent parameters is zero." << std::endl 
                << std::endl;

      exit(1);
   
   }
   else if(i >= numberOfIndependentParameters)
   {
      std::cerr << std::endl
                << "Flood Error: MultilayerPerceptron class." << std::endl 
                << "std::string getDescriptionOfIndependentParameters(int) method." << std::endl
                << "Index must be less than number of independent parameters." << std::endl 
                << std::endl;

      exit(1);
   }

   #endif

   return(descriptionOfIndependentParameters[i]);
}


// Vector<double> getMeanOfIndependentParameters(void) method

/// This method returns the mean values of all the independent parameters.
/// Such values are to be used for pre and postprocessing independent parameters with the mean and standard 
/// deviation method. 
///
/// @see getStandardDeviationOfIndependentParameters(void).
/// @see getMeanAndStandardDeviationOfIndependentParameters(void).
/// @see getMeanOfSingleIndependentParameter(int).
/// @see getStandardDeviationOfSingleIndependentParameter(int).

Vector<double> MultilayerPerceptron::getMeanOfIndependentParameters(void)
{
   return(meanOfIndependentParameters);
}


// Vector<double> getStandardDeviationOfIndependentParameters(void) method

/// This method returns the standard deviation values of all the independent parameters.
/// Such values are to be used for pre and postprocessing independent parameters with the mean and standard 
/// deviation method. 
///
/// @see getMeanOfIndependentParameters(void).
/// @see getMeanAndStandardDeviationOfIndependentParameters(void).
/// @see getMeanOfSingleIndependentParameter(int).
/// @see getStandardDeviationOfSingleIndependentParameter(int).

Vector<double> MultilayerPerceptron::getStandardDeviationOfIndependentParameters(void)
{
   return(standardDeviationOfIndependentParameters);              
}


// Matrix<double> getMeanAndStandardDeviationOfIndependentParameters(void) method

/// This method returns the mean and the standard deviation values of all the independent parameters in a single 
/// matrix. 
/// The first row contains the mean values of the independent parameters.
/// The second row contains the standard deviation values of the independent parameters.
/// Such values are to be used for pre and postprocessing independent parameters with the mean and standard 
/// deviation method. 
///
/// @see getMeanOfIndependentParameters(void).
/// @see getStandardDeviationOfIndependentParameters(void).
/// @see getMeanOfSingleIndependentParameter(int).
/// @see getStandardDeviationOfSingleIndependentParameter(int).

Matrix<double> MultilayerPerceptron::getMeanAndStandardDeviationOfIndependentParameters(void)
{
   Matrix<double> meanAndStandardDeviationOfIndependentParameters(2, numberOfIndependentParameters);

   meanAndStandardDeviationOfIndependentParameters.setRow(0, meanOfIndependentParameters);
   meanAndStandardDeviationOfIndependentParameters.setRow(1, standardDeviationOfIndependentParameters);

   return(meanAndStandardDeviationOfIndependentParameters);
}


// double getMeanOfSingleIndependentParameter(int) method

/// This method returns the mean value of a single independent parameter.
/// Such a value is to be used for pre and postprocessing that parameter with the mean and standard deviation 
/// method. 
///
/// @param i Index of independent parameter.
///
/// @see getMeanOfIndependentParameters(void).
/// @see getStandardDeviationOfIndependentParameters(void).
/// @see getMeanAndStandardDeviationOfIndependentParameters(void).
/// @see getStandardDeviationOfSingleIndependentParameter(int).

double MultilayerPerceptron::getMeanOfSingleIndependentParameter(int i)
{
   // Control sentence (if debug)

   #ifndef NDEBUG 

   if(numberOfIndependentParameters == 0)
   {
      std::cerr << std::endl
                << "Flood Error: MultilayerPerceptron class." << std::endl 
                << "double getMeanOfSingleIndependentParameter(int) method." << std::endl
                << "Number of independent parameters is zero." << std::endl 
                << std::endl;

      exit(1);
   
   }
   else if(i >= numberOfIndependentParameters)
   {
      std::cerr << std::endl
                << "Flood Error: MultilayerPerceptron class." << std::endl 
                << "double getMeanOfIndependentParameters(int) method." << std::endl
                << "Index must be less than number of independent parameters." << std::endl 
                << std::endl;

      exit(1);
   }

   #endif

   return(meanOfIndependentParameters[i]);
}


// double getStandardDeviationOfSingleIndependentParameter(int) method

/// This method returns the standard deviation value of a single independent parameter.
/// Such a value is to be used for pre and postprocessing that parameter with the mean and standard deviation 
/// method. 
///
/// @param i Index of independent parameter.
///
/// @see getMeanOfIndependentParameters(void).
/// @see getStandardDeviationOfIndependentParameters(void).
/// @see getMeanAndStandardDeviationOfIndependentParameters(void).
/// @see getMeanOfSingleIndependentParameter(int).

double MultilayerPerceptron::getStandardDeviationOfSingleIndependentParameter(int i)
{
   // Control sentence (if debug)

   #ifndef NDEBUG 

   if(numberOfIndependentParameters == 0)
   {
      std::cerr << std::endl
                << "Flood Error: MultilayerPerceptron class." << std::endl 
                << "double getStandardDeviationOfSingleIndependentParameter(int) method." << std::endl
                << "Number of independent parameters is zero." << std::endl 
                << std::endl;

      exit(1);
   
   }
   else if(i >= numberOfIndependentParameters)
   {
      std::cerr << std::endl
                << "Flood Error: MultilayerPerceptron class." << std::endl 
                << "double getStandardDeviationOfIndependentParameters(int) method." << std::endl
                << "Index must be less than number of independent parameters." << std::endl 
                << std::endl;

      exit(1);
   }

   #endif

   return(standardDeviationOfIndependentParameters[i]);
}


// Vector<double> getMinimumOfIndependentParameters(void) method

/// This method returns the minimum values of all the independent parameters.
/// Such values are to be used for pre and postprocessing independent parameters with the minimum and maximum 
/// method. 
///
/// @see getMaximumOfIndependentParameters(void).
/// @see getMinimumAndMaximumOfIndependentParameters(void).
/// @see getMinimumOfSingleIndependentParameter(int).
/// @see getMaximumOfSingleIndependentParameter(int).

Vector<double> MultilayerPerceptron::getMinimumOfIndependentParameters(void)
{
   return(minimumOfIndependentParameters);
}


// Vector<double> getMaximumOfIndependentParameters(void) method

/// This method returns the maximum values of all the independent parameters.
/// Such values are to be used for pre and postprocessing independent parameters with the minimum and maximum 
/// method. 
///
/// @see getMinimumOfIndependentParameters(void).
/// @see getMinimumAndMaximumOfIndependentParameters(void).
/// @see getMinimumOfSingleIndependentParameter(int).
/// @see getMaximumOfSingleIndependentParameter(int).

Vector<double> MultilayerPerceptron::getMaximumOfIndependentParameters(void)
{
   return(maximumOfIndependentParameters);              
}


// Matrix<double> getMinimumAndMaximumOfIndependentParameters(void) method

/// This method returns the minimum and maximum values of all the independent parameters in a single matrix. 
/// The first row contains the minimum values of the independent parameters.
/// The second row contains the maximum values of the independent parameters.
/// Such values are to be used for pre and postprocessing independent parameters with the minimum and maximum 
/// method. 
///
/// @see getMinimumOfIndependentParameters(void).
/// @see getMaximumOfIndependentParameters(void).
/// @see getMinimumOfSingleIndependentParameter(int).
/// @see getMaximumOfSingleIndependentParameter(int).

Matrix<double> MultilayerPerceptron::getMinimumAndMaximumOfIndependentParameters(void)
{
   Matrix<double> minimumAndMaximumOfIndependentParameters(2, numberOfIndependentParameters);

   minimumAndMaximumOfIndependentParameters.setRow(0, minimumOfIndependentParameters);
   minimumAndMaximumOfIndependentParameters.setRow(1, maximumOfIndependentParameters);

   return(minimumAndMaximumOfIndependentParameters);
}


// double getMinimumOfSingleIndependentParameter(int) method

/// This method returns the minimum value of a single independent parameter.
/// Such value is to be used for pre and postprocessing that independent parameter with the minimum and maximum 
/// method. 
///
/// @param i Index of independent parameter.
///
/// @see getMinimumOfIndependentParameters(void).
/// @see getMaximumOfIndependentParameters(void).
/// @see getMinimumAndMaximumOfIndependentParameters(void).
/// @see getMaximumOfSingleIndependentParameter(int).

double MultilayerPerceptron::getMinimumOfSingleIndependentParameter(int i)
{
   // Control sentence (if debug)

   #ifndef NDEBUG 

   if(numberOfIndependentParameters == 0)
   {
      std::cerr << std::endl
                << "Flood Error: MultilayerPerceptron class." << std::endl 
                << "double getMinimumOfSingleIndependentParameter(int) method." << std::endl
                << "Number of independent parameters is zero." << std::endl 
                << std::endl;

      exit(1);
   
   }
   else if(i >= numberOfIndependentParameters)
   {
      std::cerr << std::endl
                << "Flood Error: MultilayerPerceptron class." << std::endl 
                << "double getMinimumOfIndependentParameters(int) method." << std::endl
                << "Index must be less than number of independent parameters." << std::endl 
                << std::endl;

      exit(1);
   }

   #endif

   return(minimumOfIndependentParameters[i]);
}


// double getMaximumOfSingleIndependentParameter(int) method

/// This method returns the maximum value of a single independent parameter.
/// Such value is to be used for pre and postprocessing that independent parameter with the minimum and maximum 
/// method. 
///
/// @param i Index of independent parameter.
///
/// @see getMinimumOfIndependentParameters(void).
/// @see getMaximumOfIndependentParameters(void).
/// @see getMinimumAndMaximumOfIndependentParameters(void).
/// @see getMinimumOfSingleIndependentParameter(int).

double MultilayerPerceptron::getMaximumOfSingleIndependentParameter(int i)
{
   // Control sentence (if debug)

   #ifndef NDEBUG 

   if(numberOfIndependentParameters == 0)
   {
      std::cerr << std::endl
                << "Flood Error: MultilayerPerceptron class." << std::endl 
                << "double getMaximumOfSingleIndependentParameter(int) method." << std::endl
                << "Number of independent parameters is zero." << std::endl 
                << std::endl;

      exit(1);
   
   }
   else if(i >= numberOfIndependentParameters)
   {
      std::cerr << std::endl
                << "Flood Error: MultilayerPerceptron class." << std::endl 
                << "double getMaximumOfIndependentParameters(int) method." << std::endl
                << "Index must be less than number of independent parameters." << std::endl 
                << std::endl;

      exit(1);
   }

   #endif

   return(maximumOfIndependentParameters[i]);
}


// Vector<double> getLowerBoundOfIndependentParameters(void) method

/// This method returns the lower bounds of all the independent parameters.
/// These values are used to postprocess the independent parameters so that they are not less than the lower 
/// bounds. 
///
/// @see getUpperBoundOfIndependentParameters(void).
/// @see getLowerAndUpperBoundsOfIndependentParameters(void).
/// @see getLowerBoundOfSingleIndependentParameter(int).
/// @see getUpperBoundOfSingleIndependentParameter(int).

Vector<double> MultilayerPerceptron::getLowerBoundOfIndependentParameters(void)
{
   return(lowerBoundOfIndependentParameters);
}


// Vector<double> getUpperBoundOfIndependentParameters(void) method

/// This method returns the upper bounds of all the independent parameters.
/// These values are used to postprocess the independent parameters so that they are not greater than the upper 
/// bounds. 
///
/// @see getLowerBoundOfIndependentParameters(void).
/// @see getLowerAndUpperBoundsOfIndependentParameters(void).
/// @see getLowerBoundOfSingleIndependentParameter(int).
/// @see getUpperBoundOfSingleIndependentParameter(int).

Vector<double> MultilayerPerceptron::getUpperBoundOfIndependentParameters(void)
{
   return(upperBoundOfIndependentParameters);
}


// Matrix<double> getLowerAndUpperBoundsOfIndependentParameters(void) method

/// This method returns the lower and upper bounds of all the independent parameters in a single matrix. 
/// The first row contains the lower bounds of the independent parameters.
/// The second row contains the upper bounds of the independent parameters.
/// These values are used to postprocess the independent parameters so that they are neither less than the lower 
/// bounds nor greater than the upper bounds. 

/// @see getLowerBoundOfIndependentParameters(void).
/// @see getUpperBoundOfIndependentParameters(void).
/// @see getLowerBoundOfSingleIndependentParameter(int).
/// @see getUpperBoundOfSingleIndependentParameter(int).

Matrix<double> MultilayerPerceptron::getLowerAndUpperBoundsOfIndependentParameters(void)
{
   Matrix<double> lowerAndUpperBoundsOfIndependentParameters(2, numberOfIndependentParameters);

   lowerAndUpperBoundsOfIndependentParameters.setRow(0, lowerBoundOfIndependentParameters);
   lowerAndUpperBoundsOfIndependentParameters.setRow(1, upperBoundOfIndependentParameters);

   return(lowerAndUpperBoundsOfIndependentParameters);
}


// double getLowerBoundOfSingleIndependentParameter(int) method

/// This method returns the lower bound of a single independent parameter.
/// These values are used to postprocess that independent parameter so that it is not less than the lower bound. 
///
/// @param i Index of independent parameter.
///
/// @see getLowerBoundOfIndependentParameters(void).
/// @see getUpperBoundOfIndependentParameters(void).
/// @see getLowerAndUpperBoundsOfIndependentParameters(void).
/// @see getUpperBoundOfSingleIndependentParameter(int).

double MultilayerPerceptron::getLowerBoundOfSingleIndependentParameter(int i)
{
   // Control sentence (if debug)

   #ifndef NDEBUG 

   if(numberOfIndependentParameters == 0)
   {
      std::cerr << std::endl
                << "Flood Error: MultilayerPerceptron class." << std::endl 
                << "double getLowerBoundOfSingleIndependentParameter(int) method." << std::endl
                << "Number of independent parameters is zero." << std::endl 
                << std::endl;

      exit(1);   
   }
   else if(i >= numberOfIndependentParameters)
   {
      std::cerr << std::endl
                << "Flood Error: MultilayerPerceptron class." << std::endl 
                << "double getLowerBoundOfIndependentParameters(int) method." << std::endl
                << "Index must be less than number of independent parameters." << std::endl 
                << std::endl;

      exit(1);
   }

   #endif

   return(lowerBoundOfIndependentParameters[i]);
}


// double getUpperBoundOfSingleIndependentParameter(int) method

/// This method returns the upper bound of a single independent parameter.
/// These values are used to postprocess that independent parameter so that it is not greater than the upper 
/// bound. 
///
/// @param i Index of independent parameter.
///
/// @see getLowerBoundOfIndependentParameters(void).
/// @see getUpperBoundOfIndependentParameters(void).
/// @see getLowerAndUpperBoundsOfIndependentParameters(void).
/// @see getLowerBoundOfSingleIndependentParameter(int).

double MultilayerPerceptron::getUpperBoundOfSingleIndependentParameter(int i)
{
   // Control sentence (if debug)

   #ifndef NDEBUG 

   if(numberOfIndependentParameters == 0)
   {
      std::cerr << std::endl
                << "Flood Error: MultilayerPerceptron class." << std::endl 
                << "double getUpperBoundOfSingleIndependentParameter(int) method." << std::endl
                << "Number of independent parameters is zero." << std::endl 
                << std::endl;

      exit(1);  
   }
   else if(i >= numberOfIndependentParameters)
   {
      std::cerr << std::endl
                << "Flood Error: MultilayerPerceptron class." << std::endl 
                << "double getUpperBoundOfIndependentParameters(int) method." << std::endl
                << "Index must be less than number of independent parameters." << std::endl 
                << std::endl;

      exit(1);
   }

   #endif

   return(upperBoundOfIndependentParameters[i]);
}


// Vector<double> scaleIndependentParameters(Vector<double>) method

Vector<double> MultilayerPerceptron::scaleIndependentParameters(Vector<double> unscaledIndependentParameters)
{
   Vector<double> scaledIndependentParameters(numberOfIndependentParameters);

   for(int i = 0; i < numberOfIndependentParameters; i++)
   {
      switch(preAndPostProcessingMethod)   
      {
         case None:
         {
            scaledIndependentParameters[i] = unscaledIndependentParameters[i];
         }
         break;
                                     
         case MeanAndStandardDeviation:
         {
            if(standardDeviationOfIndependentParameters[i] < 1.0e-9)
            {
               if(display)
               {
                  std::cout << std::endl
                            << "Flood Warning: MultilayerPerceptron class." << std::endl
                            << "Vector<double> scaleIndependentParameters(Vector<double>) method." << std::endl
                            << "Standard deviation of independent parameter " << i << " is zero." << std::endl
                            << "That independent parameter won't be transformed." << std::endl; 
               }

               scaledIndependentParameters[i] = unscaledIndependentParameters[i];
            }
            else
            {
               scaledIndependentParameters[i] 
               = (unscaledIndependentParameters[i] - meanOfIndependentParameters[i])
               /standardDeviationOfIndependentParameters[i];
            }
         }
         break;

         case MinimumAndMaximum:
         {
            if(maximumOfIndependentParameters[i]-minimumOfIndependentParameters[i] < 1.0e-9)
            {
               if(display)
               {
                  std::cout << std::endl
                            << "Flood Warning: MultilayerPerceptron class." << std::endl
                            << "Vector<double> scaleIndependentParameters(Vector<double>) method." 
                            << std::endl
                            << "Maximum and minimum of independent parameter " << i << " are equal." << std::endl
                            << "That independent parameter won't be transformed." << std::endl; 
               }

               scaledIndependentParameters[i] = unscaledIndependentParameters[i];
            }
            else
            {
               scaledIndependentParameters[i] 
               = 2.0*(unscaledIndependentParameters[i] - minimumOfIndependentParameters[i])
               /(maximumOfIndependentParameters[i]-minimumOfIndependentParameters[i])-1.0;
            }
         }
         break;

      }// end switch
   }// end for

   return(scaledIndependentParameters);
}


// Vector<double> unscaleIndependentParameters(Vector<double>)

Vector<double> MultilayerPerceptron::unscaleIndependentParameters(Vector<double> scaledIndependentParameters)
{
   Vector<double> unscaledIndependentParameters(numberOfIndependentParameters);

   for(int i = 0; i < numberOfIndependentParameters; i++)
   {
      switch(preAndPostProcessingMethod)   
      {
         case None:
         {
            unscaledIndependentParameters[i] = scaledIndependentParameters[i]; 
         }
         break;
                                     
         case MeanAndStandardDeviation:
         {
            if(standardDeviationOfIndependentParameters[i] < 1e-9)
            {      
               if(display)
               {
                  std::cout << std::endl
                            << "Flood Warning: MultilayerPerceptron class" << std::endl
                            << "Vector<double> unscaleIndependentParameters(Vector<double>) method." 
                            << std::endl               
                            << "Standard deviation of independent parameter " << i << " is zero." << std::endl 
                            << "That independent parameter won't be transformed." << std::endl;
               }
               
               unscaledIndependentParameters[i] = scaledIndependentParameters[i];
            }      
            else
            {
               unscaledIndependentParameters[i] 
               = scaledIndependentParameters[i]*standardDeviationOfIndependentParameters[i] 
               + meanOfIndependentParameters[i];
            }
  
          }
          break;

         case MinimumAndMaximum:
         {
            if(maximumOfIndependentParameters[i]-minimumOfIndependentParameters[i] < 1e-9)
            {      
               if(display)
               {
                  std::cout << std::endl
                            << "Flood Warning: MultilayerPerceptron class" << std::endl
                            << "Vector<double> unscaleIndependentParameters(Vector<double>) method." << std::endl               
                            << "Maximum and minimum of independent parameter " << i << " are equal." << std::endl
                            << "That independent parameter won't be transformed." << std::endl; 
               }
               
               unscaledIndependentParameters[i] = scaledIndependentParameters[i];
            }      
            else
            {
               unscaledIndependentParameters[i] = 0.5*(scaledIndependentParameters[i] + 1.0)
               *(maximumOfIndependentParameters[i]-minimumOfIndependentParameters[i]) 
               + minimumOfIndependentParameters[i]; 
            }
          }
          break;

     }// end switch       
   }// end for

   return(unscaledIndependentParameters);
}


// Vector<double> getIndependentParameters(void) method

/// This method returns the values of the independent parameters.
/// These values are postprocessed so that they are neither less than the lower bounds nor greater than the upper 
/// bounds. 
///
/// @see getSingleIndependentParameter(int).

Vector<double> MultilayerPerceptron::getIndependentParameters(void)
{   
   // Apply lower and upper bounds

   for(int i = 0; i < numberOfIndependentParameters; i++)
   {
      if(independentParameters[i] < lowerBoundOfIndependentParameters[i])
     {
        independentParameters[i] = lowerBoundOfIndependentParameters[i];
     }
      else if(independentParameters[i] > upperBoundOfIndependentParameters[i])
     {
        independentParameters[i] = upperBoundOfIndependentParameters[i];
     }
   }

   return(independentParameters);    
}


// double getSingleIndependentParameter(int) method

/// This method returns the value of a single independent parameter.
/// Such a value is postprocessed so that it is neither less than the lower bound nor greater than the upper 
/// bound. 
///
/// @param i Index of independent parameter.
///
/// @see getIndependentParameters(void).

double MultilayerPerceptron::getSingleIndependentParameter(int i)
{
   // Control sentence (if debug)

   #ifndef NDEBUG 

   if(numberOfIndependentParameters == 0)
   {
      std::cerr << std::endl
                << "Flood Error: MultilayerPerceptron class." << std::endl 
                << "double getSingleIndependentParameter(int) method." << std::endl
                << "Number of independent parameters is zero." << std::endl 
            << std::endl;

      exit(1);
   
   }
   else if(i >= numberOfIndependentParameters)
   {
      std::cerr << std::endl
                << "Flood Error: MultilayerPerceptron class." << std::endl 
                << "double getSingleIndependentParameter(int) method." << std::endl
                << "Index must be less than number of independent parameters." << std::endl 
                << std::endl;

      exit(1);
   }

   #endif

   if(independentParameters[i] < lowerBoundOfIndependentParameters[i])
   {
      independentParameters[i] = lowerBoundOfIndependentParameters[i];
   }
   else if(independentParameters[i] > upperBoundOfIndependentParameters[i])
   {
      independentParameters[i] = upperBoundOfIndependentParameters[i];
   }

   return(independentParameters[i]);           
}


// bool getDisplayOutOfRangeWarning(void) method
/// This method returns true if a warning message is to be displayed on the screen when some input variables fall outside
/// the range defined by the minimum and maximum, or false otherwise.

bool MultilayerPerceptron::getDisplayOutOfRangeWarning(void)
{
   return(displayOutOfRangeWarning);
}


// bool getDisplay(void) method

/// This method returns true if messages from this class are to be displayed on the screen, or false if messages 
/// from this class are not to be displayed on the screen.

bool MultilayerPerceptron::getDisplay(void)
{
   return(display);
}


// double getNumericalEpsilon(void) method

/// This method returns the epsilon value to be used for computing the Jacobian matrix of the neural network for 
/// a set of input signals and by means of numerical differentiation.
///
/// @see calculateJacobianTest(Vector<double>).

double MultilayerPerceptron::getNumericalEpsilon(void)
{
   return(numericalEpsilon);
}


// void setPreAndPostProcessingMethod(PreAndPostProcessingMethod)

/// This method sets the method to be used for preprocessing network inputs and postprocessing network outputs.
///
/// @param newPreAndPostProcessingMethod New pre and post-processing method.

void MultilayerPerceptron::setPreAndPostProcessingMethod
(MultilayerPerceptron::PreAndPostProcessingMethod newPreAndPostProcessingMethod)
{
   preAndPostProcessingMethod = newPreAndPostProcessingMethod;
}


// void setNumericalDifferentiationMethod(NumericalDifferentiationMethod)

/// This method sets the method to be used for the numerical differentiation of the Jacobian matrix for the 
/// multilayer perceptron.
///
/// @param newNumericalDifferentiationMethod New numerical differentiation method.
///
/// @see calculateJacobianTest(Vector<double>).

void MultilayerPerceptron::setNumericalDifferentiationMethod
(MultilayerPerceptron::NumericalDifferentiationMethod newNumericalDifferentiationMethod)
{
   numericalDifferentiationMethod = newNumericalDifferentiationMethod;
}


// void setNetworkArchitecture(int, int, int) method

/// This method sets the architecture of the neural network. 
/// All the free parameters of the neural network are initialized with random values chosen from a normal 
/// distribution with mean zero and standard deviation one. 
/// It also initializes the rest of class members to their default values:
///
/// <ul>
/// <li> Name of input vaviables: "InputVariable0",... 
/// <li> Name of output variables: "OutputVariable0",...
/// <li> Units of input vaviables: "None",... 
/// <li> Units of output variables: "None",...
/// <li> Description of input vaviables: "None",... 
/// <li> Description of output variables: "None",...
/// <li> Mean of input variables: 0.
/// <li> Standard deviation of input variables: 1.
/// <li> Mean of output variables: 0.
/// <li> Standard deviation of output variables: 1.
/// <li> Minimum of input variables: -Inf.
/// <li> Maximum of input variables: +Inf.
/// <li> Minimum of output variables: -Inf.
/// <li> Maximum of output variables: +Inf.
/// <li> Lower bound of output variables: -Inf.
/// <li> Lower bound of output variables: +Inf.
/// </ul> 
///
/// @param newNumberOfInputs New number of inputs in the neural network.
/// @param newNumbersOfHiddenNeurons New numbers of neurons for the hidden layers of the neural network.
/// The size of this vector is the number of hidden layers. 
/// @param newNumberOfOutputs New number of output neurons in the neural network.

void MultilayerPerceptron::setNetworkArchitecture
(int newNumberOfInputs, Vector<int> newNumbersOfHiddenNeurons, int newNumberOfOutputs)
{
   // Set new architecture      
   
   numberOfInputs = newNumberOfInputs;
   numbersOfHiddenNeurons = newNumbersOfHiddenNeurons;
   numberOfOutputs = newNumberOfOutputs;

   int numberOfHiddenLayers = getNumberOfHiddenLayers();

   // Set activation functions

   hiddenLayersActivationFunction.setSize(numberOfHiddenLayers);

   for(int i = 0; i < numberOfHiddenLayers; i++)
   {
      hiddenLayersActivationFunction[i] = Perceptron::HyperbolicTangent;
   }

   outputLayerActivationFunction = Perceptron::Linear;

   // Hidden layer

   hiddenLayers.setSize(numberOfHiddenLayers);

   for(int i = 0; i < numberOfHiddenLayers; i++)
   {
      hiddenLayers[i].setSize(numbersOfHiddenNeurons[i]);

      for(int j = 0; j < numbersOfHiddenNeurons[i]; j++)
      {
         hiddenLayers[i][j].setActivationFunction(hiddenLayersActivationFunction[i]);

         if(i == 0)// First hidden layer
         {
            hiddenLayers[i][j].setNumberOfInputs(numberOfInputs);
         }
         else
         {
            hiddenLayers[i][j].setNumberOfInputs(numbersOfHiddenNeurons[i-1]);
         }       
      }
   }

   // Output layer

   outputLayer.setSize(numberOfOutputs);

   for(int i = 0; i < numberOfOutputs; i++)
   {
      outputLayer[i].setActivationFunction(Perceptron::Linear);

      outputLayer[i].setNumberOfInputs(numbersOfHiddenNeurons[numberOfHiddenLayers-1]);
   }

   // Name of input variables

   nameOfInputVariables.setSize(numberOfInputs);

   for(int i = 0; i < numberOfInputs; i++)
   {
      std::stringstream buffer;

      buffer << "InputVariable" << i;

      nameOfInputVariables[i] = buffer.str();
   }

   // Name of output variables

   nameOfOutputVariables.setSize(numberOfOutputs);

   for(int i = 0; i < numberOfOutputs; i++)
   {
      std::stringstream buffer;

      buffer << "OutputVariable" << i;

      nameOfOutputVariables[i] = buffer.str();
   }

   // Units of input variables

   unitsOfInputVariables.setSize(numberOfInputs);

   for(int i = 0; i < numberOfInputs; i++)
   {
      unitsOfInputVariables[i] = "None";
   }

   // Units of output variables

   unitsOfOutputVariables.setSize(numberOfOutputs);

   for(int i = 0; i < numberOfOutputs; i++)
   {
      unitsOfOutputVariables[i] = "None";
   }

   // Description of input variables

   descriptionOfInputVariables.setSize(numberOfInputs);

   for(int i = 0; i < numberOfInputs; i++)
   {
      descriptionOfInputVariables[i] = "None";
   }

   // Description of output variables

   descriptionOfOutputVariables.setSize(numberOfOutputs);

   for(int i = 0; i < numberOfOutputs; i++)
   {
      descriptionOfOutputVariables[i] = "None";
   }

   // Initialize mean of input variables to zero

   meanOfInputVariables.setSize(numberOfInputs);    

   for(int i = 0; i < numberOfInputs; i++)
   {
      meanOfInputVariables[i] = 0.0;
   }

   // Initialize standard deviation of input variables to one

   standardDeviationOfInputVariables.setSize(numberOfInputs);    

   for(int i = 0; i < numberOfInputs; i++)
   {
      standardDeviationOfInputVariables[i] = 1.0;
   }

   // Initialize mean of output variables to zero

   meanOfOutputVariables.setSize(numberOfOutputs);    

   for(int i = 0; i < numberOfOutputs; i++)
   {
      meanOfOutputVariables[i] = 0.0;
   }

   // Initialize Standard deviation of output variables to one

   standardDeviationOfOutputVariables.setSize(numberOfOutputs);    

   for(int i = 0; i < numberOfOutputs; i++)
   {
      standardDeviationOfOutputVariables[i] = 1.0;
   }

   // Initialize minimum of input variables to -1

   minimumOfInputVariables.setSize(numberOfInputs);    

   for(int i = 0; i < numberOfInputs; i++)
   {
      minimumOfInputVariables[i] = -1.0;
   }

   // Initialize maximum of input variables to +1

   maximumOfInputVariables.setSize(numberOfInputs);    

   for(int i = 0; i < numberOfInputs; i++)
   {
      maximumOfInputVariables[i] = 1.0;
   }

   // Initialize minimum of output variables to -1

   minimumOfOutputVariables.setSize(numberOfOutputs);    

   for(int i = 0; i < numberOfOutputs; i++)
   {
      minimumOfOutputVariables[i] = -1.0;
   }

   // Initialize maximum of output variables to 1

   maximumOfOutputVariables.setSize(numberOfOutputs);    

   for(int i = 0; i < numberOfOutputs; i++)
   {
      maximumOfOutputVariables[i] = 1.0;
   }

   // Initialize lower bound of output variables to -Inf

   lowerBoundOfOutputVariables.setSize(numberOfOutputs);    

   for(int i = 0; i < numberOfOutputs; i++)
   {
      lowerBoundOfOutputVariables[i] = -1.0e99;
   }

   // Initialize upper bound of output variables to Inf

   upperBoundOfOutputVariables.setSize(numberOfOutputs);    

   for(int i = 0; i < numberOfOutputs; i++)
   {
      upperBoundOfOutputVariables[i] = 1.0e99;
   }
}


void MultilayerPerceptron::setNumbersOfHiddenNeurons(Vector<int> newNumbersOfHiddenNeurons)
{
   // Set new architecture      

   numbersOfHiddenNeurons = newNumbersOfHiddenNeurons;

   int numberOfHiddenLayers = numbersOfHiddenNeurons.getSize();

   // Activation functions
   hiddenLayersActivationFunction.setSize(numberOfHiddenLayers);

   for(int i = 0; i < numberOfHiddenLayers; i++)
   {
      hiddenLayersActivationFunction[i] = Perceptron::HyperbolicTangent;
   }

   hiddenLayers.setSize(numberOfHiddenLayers);

   for(int i = 0; i < numberOfHiddenLayers; i++)
   {
      hiddenLayers[i].setSize(numbersOfHiddenNeurons[i]);

      for(int j = 0; j < numbersOfHiddenNeurons[i]; j++)
      {         hiddenLayers[i][j].setActivationFunction(hiddenLayersActivationFunction[i]);

         if(i == 0)// First hidden layer
         {
            hiddenLayers[i][j].setNumberOfInputs(numberOfInputs);
         }
         else
         {
            hiddenLayers[i][j].setNumberOfInputs(numbersOfHiddenNeurons[i-1]);
         }       
      }
   }

   // Output layer

   for(int i = 0; i < numberOfOutputs; i++)
   {
      outputLayer[i].setActivationFunction(outputLayerActivationFunction);     

      outputLayer[i].setNumberOfInputs(numbersOfHiddenNeurons[numberOfHiddenLayers-1]);
   }
}


// void setHiddenLayerActivationFunction(Perceptron::ActivationFunction) method

/// This class sets a new activation (or transfer) function in all the perceptrons composing the hidden layer. 
///
/// @param newHiddenLayersActivationFunction Activation function for the hidden layers.
/// The size of this Vector must be equal to the number of hidden layers, and each element corresponds
/// to the activation function of one hidden layer. 
///
/// @see Perceptron.
/// @see setOutputLayerActivationFunction(Perceptron::ActivationFunction).

void MultilayerPerceptron::setHiddenLayersActivationFunction(Vector<Perceptron::ActivationFunction> 
newHiddenLayersActivationFunction)
{
   hiddenLayersActivationFunction = newHiddenLayersActivationFunction;

   int numberOfHiddenLayers = numbersOfHiddenNeurons.getSize();

   for(int i = 0; i < numberOfHiddenLayers; i++)
   {
      for(int j = 0; j < numbersOfHiddenNeurons[i]; j++)
      {         hiddenLayers[i][j].setActivationFunction(hiddenLayersActivationFunction[i]);
      }
   }
}


// void setOutputLayerActivationFunction(Perceptron::ActivationFunction) method

/// This class sets a new activation (or transfer) function in all the perceptrons composing the output layer. 
///
/// @param newOutputLayerActivationFunction Activation function for the output units. 
///
/// @see Perceptron.
/// @see setHiddenLayerActivationFunction(Perceptron::ActivationFunction).

void MultilayerPerceptron::setOutputLayerActivationFunction(Perceptron::ActivationFunction 
newOutputLayerActivationFunction)
{
   outputLayerActivationFunction = newOutputLayerActivationFunction;

   for(int i = 0; i < numberOfOutputs; i++)
   {
      outputLayer[i].setActivationFunction(outputLayerActivationFunction);
   }
}


// void setNameOfInputVariables(Vector<std::string>) method

/// This method sets the names of the input variables.
/// Such values are only used to give the user basic information on the problem at hand.
///
/// @param newNameOfInputVariables New names for the input variables of the neural network.
///
/// @see setNameOfOutputVariables(Vector<std::string>).
/// @see setNameOfSingleInputVariable(int, std::string).
/// @see setNameOfSingleOutputVariable(int, std::string).

void MultilayerPerceptron::setNameOfInputVariables(Vector<std::string> newNameOfInputVariables)
{
   // Control sentence (if debug)

   #ifndef NDEBUG 

   if(newNameOfInputVariables.getSize() != numberOfInputs)
   {
      std::cerr << std::endl
                << "Flood Error: MultilayerPerceptron class." << std::endl 
                << "void setNameOfInputVariables(Vector<std::string>) method." << std::endl
                << "Size of name of input variables vector must be equal to number of inputs." << std::endl 
                << std::endl;

      exit(1);
   }

   #endif

   // Set name of input variables

   nameOfInputVariables = newNameOfInputVariables;
}


// void setNameOfOutputVariables(Vector<std::string>) method

/// This method sets the names of the output variables.
/// Such values are only used to give the user basic information on the problem at hand.
///
/// @param newNameOfOutputVariables New names for the output variables of the neural network.
///
/// @see setNameOfInputVariables(Vector<std::string>).
/// @see setNameOfSingleInputVariable(int, std::string).
/// @see setNameOfSingleOutputVariable(int, std::string).

void MultilayerPerceptron::setNameOfOutputVariables(Vector<std::string> newNameOfOutputVariables)
{
   // Control sentence (if debug)

   #ifndef NDEBUG 

   if(newNameOfOutputVariables.getSize() != numberOfOutputs)
   {
      std::cerr << std::endl
                << "Flood Error: MultilayerPerceptron class." << std::endl
                << "void setNameOfOutputVariables(Vector<double>) method." << std::endl 
                << "Size of name of output variables vector must be equal to number of outputs." << std::endl
                << std::endl;

      exit(1);
   }

   #endif

   // Set name of output variables

   nameOfOutputVariables = newNameOfOutputVariables;
}


// void setNameOfSingleInputVariable(int, std::string) method

/// This method sets the name of a single input variable.
/// Such value is only used to give the user basic information on the problem at hand.
///
/// @param i Index of input variable.
/// @param newNameOfSingleInputVariable New name for the input variable with index i.
///
/// @see setNameOfInputVariables(Vector<std::string>).
/// @see setNameOfOutputVariables(Vector<std::string>).
/// @see setNameOfSingleOutputVariable(int, std::string).

void MultilayerPerceptron::setNameOfSingleInputVariable(int i, std::string newNameOfSingleInputVariable)
{
   // Control sentence (if debug)

   #ifndef NDEBUG 

   if(i >= numberOfInputs)
   {
      std::cerr << std::endl
                << "Flood Error: MultilayerPerceptron class." << std::endl
                << "void setNameOfSingleInputVariable(int, std::string) method." << std::endl
                << "Index must be less than number of inputs." << std::endl
                << std::endl;

      exit(1);   
   }

   #endif

   // Set name of single input variable

   nameOfInputVariables[i] = newNameOfSingleInputVariable;
}


// void setNameOfSingleOutputVariable(int, std::string) method

/// This method sets the name of a single output variable.
/// Such value is only used to give the user basic information on the problem at hand.
///
/// @param i Index of output variable.
/// @param newNameOfSingleOutputVariable New name for the output variable with index i.
///
/// @see setNameOfInputVariables(Vector<std::string>).
/// @see setNameOfOutputVariables(Vector<std::string>).
/// @see setNameOfSingleInputVariable(int, std::string).

void MultilayerPerceptron::setNameOfSingleOutputVariable(int i, std::string newNameOfSingleOutputVariable)
{
   // Control sentence (if debug)

   #ifndef NDEBUG 

   if(i >= numberOfOutputs)
   {
      std::cerr << std::endl
                << "Flood Error: MultilayerPerceptron class." << std::endl
                << "void setNameOfSingleOutputVariable(int, std::string) method." << std::endl
                << "Index must be less than number of ouptuts." << std::endl
                << std::endl;

      exit(1);   
   }

   #endif

   // Set name of single output variable

   nameOfOutputVariables[i] = newNameOfSingleOutputVariable;
}


// void setUnitsOfInputVariables(Vector<std::string>) method

/// This method sets new units for all the input variables.
/// Such values are only used to give the user basic information on the problem at hand.
///
/// @param newUnitsOfInputVariables New units for the input variables of the neural network.
///
/// @see setUnitsOfOutputVariables(Vector<std::string>).
/// @see setUnitsOfSingleInputVariable(int, std::string).
/// @see setUnitsOfSingleOutputVariable(int, std::string).

void MultilayerPerceptron::setUnitsOfInputVariables(Vector<std::string> newUnitsOfInputVariables)
{
   // Control sentence (if debug)

   #ifndef NDEBUG 
 
   if(newUnitsOfInputVariables.getSize() != numberOfInputs)
   {
      std::cerr << std::endl
                << "Flood Error: MultilayerPerceptron class." << std::endl 
                << "void setUnitsOfInputVariables(Vector<double>) method." << std::endl
                << "Size must be equal to number of inputs." << std::endl
                << std::endl;

      exit(1);
   }

   #endif

   // Set units of input variables

   unitsOfInputVariables = newUnitsOfInputVariables;
}


// void setUnitsOfOutputVariables(Vector<std::string>) method

/// This method sets new units for all the output variables.
/// Such values are only used to give the user basic information on the problem at hand.
///
/// @param newUnitsOfOutputVariables New units for the input variables of the neural network.
///
/// @see setUnitsOfInputVariables(Vector<std::string>).
/// @see setUnitsOfSingleInputVariable(int, std::string).
/// @see setUnitsOfSingleOutputVariable(int, std::string).

void MultilayerPerceptron::setUnitsOfOutputVariables(Vector<std::string> newUnitsOfOutputVariables)
{
   // Control sentence (if debug)

   #ifndef NDEBUG 

   if(newUnitsOfOutputVariables.getSize() != numberOfOutputs)
   {
      std::cerr << std::endl
                << "Flood Error: MultilayerPerceptron class."  << std::endl 
                << "void setUnitsOfOutputVariables(Vector<double>) method." << std::endl
                << "Size must be equal to number of outputs." << std::endl
                << std::endl;

      exit(1);
   }

   #endif

   // Set units of output variables

   unitsOfOutputVariables = newUnitsOfOutputVariables;
}


// void setUnitsOfSingleInputVariable(int, std::string) method

/// This method sets new units for a single input variable.
/// Such value is only used to give the user basic information on the problem at hand.
///
/// @param i Index of input variable.
/// @param newUnitsOfSingleInputVariable New units for the input variable with index i.
///
/// @see setUnitsOfInputVariables(Vector<std::string>).
/// @see setUnitsOfOutputVariables(Vector<std::string>).
/// @see setUnitsOfSingleOutputVariable(int, std::string).

void MultilayerPerceptron::setUnitsOfSingleInputVariable(int i, 
std::string newUnitsOfSingleInputVariable)
{
   // Control sentence (if debug)

   #ifndef NDEBUG 

   if(i >= numberOfInputs)
   {
      std::cerr << std::endl
                << "Flood Error: MultilayerPerceptron class." << std::endl
                << "void setUnitsOfSingleInputVariable(int, std::string) method." << std::endl
                << "Index must be less than number of inputs." << std::endl 
                << std::endl;

      exit(1);
   }

   #endif

   // Set units of single input variable

   unitsOfInputVariables[i] = newUnitsOfSingleInputVariable;
}


// void setUnitsOfSingleOutputVariable(int, std::string) method

/// This method sets new units for a single output variable.
/// Such value is only used to give the user basic information on the problem at hand.
///
/// @param i Index of output variable.
/// @param newUnitsOfSingleOutputVariable New units for the output variable with index i.
///
/// @see setUnitsOfInputVariables(Vector<std::string>).
/// @see setUnitsOfOutputVariables(Vector<std::string>).
/// @see setUnitsOfSingleInputVariable(int, std::string).

void MultilayerPerceptron::setUnitsOfSingleOutputVariable(int i, 
std::string newUnitsOfSingleOutputVariable)
{
   // Control sentence (if debug)

   #ifndef NDEBUG 

   if(i >= numberOfOutputs)
   {
      std::cerr << std::endl
                << "Flood Error: MultilayerPerceptron class." << std::endl
                << "void setUnitsOfSingleOutputVariable(int, std::string) method." << std::endl 
                << "Index must be less than number of outputs."  << std::endl
                << std::endl;

      exit(1);   
   }

   #endif

   // Set units of single output variable

   unitsOfOutputVariables[i] = newUnitsOfSingleOutputVariable;
}


// void setDescriptionOfInputVariables(Vector<std::string>) method

/// This method sets new descriptions for all the input variables.
/// Such values are only used to give the user basic information on the problem at hand.
///
/// @param newDescriptionOfInputVariables New description for the input variables of the neural network.
///
/// @see setDescriptionOfOutputVariables(Vector<std::string>).
/// @see setDescriptionOfSingleInputVariable(int, std::string).
/// @see setDescriptionOfSingleOutputVariable(int, std::string).

void MultilayerPerceptron::setDescriptionOfInputVariables(Vector<std::string> newDescriptionOfInputVariables)
{
   // Control sentence (if debug)

   #ifndef NDEBUG 

   if(newDescriptionOfInputVariables.getSize() != numberOfInputs)
   {
      std::cerr << std::endl
                << "Flood Error: MultilayerPerceptron class." << std::endl  
                << "void setDescriptionOfInputVariables(Vector<double>) method." << std::endl
                << "Size must be equal to number of inputs." << std::endl
                << std::endl;

      exit(1);
   }

   #endif

   // Set description of input variables

   descriptionOfInputVariables = newDescriptionOfInputVariables;
}


// void setDescriptionOfOutputVariables(Vector<std::string>) method

/// This method sets new descriptions for all the output variables.
/// Such values are only used to give the user basic information on the problem at hand.
///
/// @param newDescriptionOfOutputVariables New description for the output variables of the neural network.
///
/// @see setDescriptionOfInputVariables(Vector<std::string>).
/// @see setDescriptionOfSingleInputVariable(int, std::string).
/// @see setDescriptionOfSingleOutputVariable(int, std::string).

void MultilayerPerceptron::setDescriptionOfOutputVariables(Vector<std::string> newDescriptionOfOutputVariables)
{
   // Control sentence (if debug)

   #ifndef NDEBUG 

   if(newDescriptionOfOutputVariables.getSize() != numberOfOutputs)
   {
      std::cerr << std::endl
                << "Flood Error: MultilayerPerceptron class." << std::endl 
                << "void setDescriptionOfOutputVariables(Vector<double>) method." << std::endl
                << "Size must be equal to number of outputs." << std::endl
                << std::endl;

      exit(1);
   }

   #endif

   // Set description of output variables

   descriptionOfOutputVariables = newDescriptionOfOutputVariables;
}


// void setDescriptionOfSingleInputVariable(int, std::string) method

/// This method sets a new description for a single input variable.
/// Such value is only used to give the user basic information on the problem at hand.
///
/// @param i Index of input variable.
/// @param newDescriptionOfSingleInputVariable New description for the input variable with index i.
///
/// @see setDescriptionOfInputVariables(Vector<std::string>).
/// @see setDescriptionOfOutputVariables(Vector<std::string>).
/// @see setDescriptionOfSingleOutputVariable(int, std::string).

void MultilayerPerceptron::setDescriptionOfSingleInputVariable(int i, 
std::string newDescriptionOfSingleInputVariable)
{
   // Control sentence (if debug)

   #ifndef NDEBUG 

   if(i >= numberOfInputs)
   {
      std::cerr << std::endl
                << "Flood Error: MultilayerPerceptron class." << std::endl
                << "void setDescriptionOfSingleInputVariable(int, std::string) method." << std::endl
                << "Index must be less than number of inputs." << std::endl 
                << std::endl;

      exit(1);   
   }

   #endif

   // Set description of single input variable

   descriptionOfInputVariables[i] = newDescriptionOfSingleInputVariable;
}


// void setDescriptionOfSingleOutputVariable(int, std::string) method

/// This method sets a new description for a single output variable.
/// Such value is only used to give the user basic information on the problem at hand.
///
/// @param i Index of output variable.
/// @param newDescriptionOfSingleOutputVariable New description for the output variable with index i.
///
/// @see setDescriptionOfInputVariables(Vector<std::string>).
/// @see setDescriptionOfOutputVariables(Vector<std::string>).
/// @see setDescriptionOfSingleInputVariable(int, std::string).

void MultilayerPerceptron::setDescriptionOfSingleOutputVariable(int i, 
std::string newDescriptionOfSingleOutputVariable)
{
   // Control sentence (if debug)

   #ifndef NDEBUG 

   if(i >= numberOfOutputs)
   {
      std::cerr << std::endl
                << "Flood Error: MultilayerPerceptron class." << std::endl
                << "void setDescriptionOfSingleOutputVariable(int, std::string) method." << std::endl 
                << "Index must be less than number of outputs." << std::endl
                << std::endl;

      exit(1);   
   }

   #endif

   // Set description of single output variable

   descriptionOfOutputVariables[i] = newDescriptionOfSingleOutputVariable;
}


// void setAllInformation(Vector< Vector<std::string> >) method

void MultilayerPerceptron::setAllInformation(Vector< Vector<std::string> > newInformation)
{
   // Control sentence (if debug)

   #ifndef NDEBUG 

   int size = newInformation.getSize();

   if(size != 6)
   {
      std::cerr << std::endl
                << "Flood Error: MultilayerPerceptron class." << std::endl
                << "void setAllInformation(Vector< Vector<std::string> >) method." << std::endl 
                << "Size of vector must be 6." << std::endl
                << std::endl;

      exit(1);   
   }

   #endif

   // Set all information

   setNameOfInputVariables(newInformation[0]);
   setNameOfOutputVariables(newInformation[1]);

   setUnitsOfInputVariables(newInformation[2]);
   setUnitsOfOutputVariables(newInformation[3]);

   setDescriptionOfInputVariables(newInformation[4]);
   setDescriptionOfOutputVariables(newInformation[5]);
}



//void setMeanOfInputVariables(Vector<double>) method

/// This method sets new mean values for all the input variables.
/// These values are used for preprocessing the inputs to the neural network with the meand and standard 
/// deviation method. 
///
/// @param newMeanOfInputVariables New set of mean values for the input variables of the neural network.
///
/// @see setStandardDeviationOfInputVariables(Vector<double>).
/// @see setMeanAndStandardDeviationOfInputVariables(Matrix<double>).
/// @see setMeanOfSingleInputVariable(int, double).
/// @see setStandardDeviationOfSingleInputVariable(int, double).

void MultilayerPerceptron::setMeanOfInputVariables(Vector<double> newMeanOfInputVariables)
{
   // Control sentence (if debug)

   #ifndef NDEBUG 

   if(newMeanOfInputVariables.getSize() != numberOfInputs)
   {
      std::cerr << std::endl
                << "Flood Error: MultilayerPerceptron class." << std::endl 
                << "void setMeanOfInputVariables(Vector<double>) method." << std::endl
                << "Size must be equal to number of inputs." << std::endl
                << std::endl;

      exit(1);
   }

   #endif

   // Set mean of input variables

   meanOfInputVariables = newMeanOfInputVariables;
}


// void setStandardDeviationOfInputVariables(Vector<double>) method

/// This method sets new standard deviation values for all the input variables.
/// These values are used for preprocessing the inputs to the neural network with the meand and standard deviation
/// method. 
///
/// @param newStandardDeviationOfInputVariables New set of standard deviation values for the input variables of 
/// the neural network.
///
/// @see setMeanOfInputVariables(Vector<double>).
/// @see setMeanAndStandardDeviationOfInputVariables(Matrix<double>).
/// @see setMeanOfSingleInputVariable(int, double).
/// @see setStandardDeviationOfSingleInputVariable(int, double).

void MultilayerPerceptron::setStandardDeviationOfInputVariables(Vector<double> newStandardDeviationOfInputVariables)
{
   // Control sentence (if debug)

   #ifndef NDEBUG 

   if(newStandardDeviationOfInputVariables.getSize() != numberOfInputs)
   {
      std::cerr << std::endl
                << "Flood Error: MultilayerPerceptron class." << std::endl
                << "void setStandardDeviationOfInputVariables(Vector<double>) method." << std::endl
                << "Size must be equal to number of inputs." << std::endl
                << std::endl;

      exit(1);
   }

   #endif

   // Set standard deviation of input variables

   standardDeviationOfInputVariables = newStandardDeviationOfInputVariables;
}


// void setMeanAndStandardDeviationOfInputVariables(Matrix<double>) method

/// This method sets both the mean and the standard deviation values of all the input variables from a single 
/// matrix.
/// The first row must contain the mean values for the input variables.
/// The second row must contain the standard deviation values for the input variables.
/// These values are used for preprocessing the inputs to the neural network with the meand and standard deviation
/// method. 
///
/// @param newMeanAndStandardDeviationOfInputVariables New set of mean and standard deviation values for the 
/// input variables of the neural network.
///
/// @see setMeanOfInputVariables(Vector<double>).
/// @see setStandardDeviationOfInputVariables(Vector<double>).
/// @see setMeanOfSingleInputVariable(int, double).
/// @see setStandardDeviationOfSingleInputVariable(int, double).

void MultilayerPerceptron::setMeanAndStandardDeviationOfInputVariables
(Matrix<double> newMeanAndStandardDeviationOfInputVariables)
{
   // Control sentence (if debug)

   #ifndef NDEBUG 

   if(newMeanAndStandardDeviationOfInputVariables.getNumberOfRows() != 2)
   {
      std::cerr << std::endl
                << "Flood Error: MultilayerPerceptron class." << std::endl
                << "void setMeanAndStandardDeviationOfInputVariables(Matrix<double>) method." << std::endl
                << "Number of rows must be 2." << std::endl
                << std::endl;

      exit(1);
   }
   else if(newMeanAndStandardDeviationOfInputVariables.getNumberOfColumns() != numberOfInputs)
   {
      std::cerr << std::endl
                << "Flood Error: MultilayerPerceptron class." << std::endl
                << "void setMeanAndStandardDeviationOfInputVariables(Matrix<double>) method." << std::endl
                << "Number of columns must be equal to number of inputs." << std::endl
                << std::endl;

      exit(1);
   }

   // Check that mean of input variables is not zero

   for(int i = 0; i < numberOfInputs; i++)
   {
      if(newMeanAndStandardDeviationOfInputVariables[1][i] < 1.0e-69)
      {
         std::cerr << std::endl
                   << "Flood Error: MultilayerPerceptron class: " << std::endl
                   << "void setMeanAndStandardDeviationOfInputVariables(Matrix<double>) method." << std::endl
                   << "Standard deviation of input variable "<< i << " is zero." << std::endl
                   << std::endl;

          exit(1);
      }
   }

   #endif

   // Set mean and standard deviation of input variables

   meanOfInputVariables = newMeanAndStandardDeviationOfInputVariables.getRow(0);

   standardDeviationOfInputVariables = newMeanAndStandardDeviationOfInputVariables.getRow(1);
}


//void setMeanOfOutputVariables(Vector< double>) method

/// This method sets new mean values for all the output variables.
/// These values are used for postprocessing the output signals from the neural network with the meand and 
/// standard deviation method. 
///
/// @param newMeanOfOutputVariables New set of mean values for the output variables of the neural network.
///
/// @see setStandardDeviationOfOutputVariables(Vector<double>).
/// @see setMeanAndStandardDeviationOfOutputVariables(Matrix<double>).
/// @see setMeanOfSingleOutputVariable(int, double).
/// @see setStandardDeviationOfSingleOutputVariable(int, double).

void MultilayerPerceptron::setMeanOfOutputVariables(Vector< double> newMeanOfOutputVariables)
{
   // Control sentence (if debug)

   #ifndef NDEBUG 

   if(newMeanOfOutputVariables.getSize() != numberOfOutputs)
   {
      std::cerr << std::endl
                << "Flood Error: MultilayerPerceptron class." << std::endl 
                << "void setMeanOfOutputVariables(Vector<double>) method." << std::endl
                << "Size must be equal to number of outputs." << std::endl
                << std::endl;

      exit(1);
   }

   #endif

   // Set mean of output variables

   meanOfOutputVariables = newMeanOfOutputVariables;
}


// void setStandardDeviationOfOutputVariables(Vector<double>) method

/// This method sets new standard deviation values for all the output variables.
/// These values are used for postprocessing the output signals from the neural network with the meand and 
/// standard deviation method. 
///
/// @param newStandardDeviationOfOutputVariables New set of standard deviation values for the output variables of 
/// the neural network.
///
/// @see setMeanOfOutputVariables(Vector<double>).
/// @see setMeanAndStandardDeviationOfOutputVariables(Matrix<double>).
/// @see setMeanOfSingleOutputVariable(int, double).
/// @see setStandardDeviationOfSingleOutputVariable(int, double).

void MultilayerPerceptron
::setStandardDeviationOfOutputVariables(Vector<double> newStandardDeviationOfOutputVariables) 
{
   // Control sentence (if debug)

   #ifndef NDEBUG 

   if(newStandardDeviationOfOutputVariables.getSize() != numberOfOutputs)
   {
      std::cerr << std::endl
                << "Flood Error: MultilayerPerceptron class." << std::endl
                << "void setStandarDeviationOfOutputVariables(Vector<double>) method." << std::endl
                << "Size must be equal to number of outputs." << std::endl
                << std::endl;

      exit(1);
   }

   #endif

   // Set standard deviation of output variables

   standardDeviationOfOutputVariables = newStandardDeviationOfOutputVariables;
}


// void setMeanAndStandardDeviationOfOutputVariables(Matrix<double>) method

/// This method sets both the mean and the standard deviation values of all the output variables from a single 
/// matrix.
/// The first row must contain the mean values for the output variables.
/// The second row must contain the standard deviation values for the output variables.
///
/// @param newMeanAndStandardDeviationOfOutputVariables New set of mean and standard deviation values for the 
/// output variables of the neural network.
///
/// @see setMeanOfOutputVariables(Vector<double>).
/// @see setStandardDeviationOfOutputVariables(Vector<double>).
/// @see setMeanOfSingleOutputVariable(int, double).
/// @see setStandardDeviationOfSingleOutputVariable(int, double).

void MultilayerPerceptron
::setMeanAndStandardDeviationOfOutputVariables(Matrix<double> newMeanAndStandardDeviationOfOutputVariables)
{
   // Control sentence (if debug)

   #ifndef NDEBUG 

   if(newMeanAndStandardDeviationOfOutputVariables.getNumberOfRows() != 2)
   {
      std::cerr << std::endl
                << "Flood Error: Multilayer Perceptron class." << std::endl
                << "void setMeanAndStandardDeviationOfOutputVariables(Matrix<double>) method." << std::endl
                << "Number of rows must be 2." << std::endl
                << std::endl;

      exit(1);
   }
   else if(newMeanAndStandardDeviationOfOutputVariables.getNumberOfColumns() != numberOfOutputs)
   {
      std::cerr << std::endl
                << "Flood Error: MultilayerPerceptron class." << std::endl
                << "void setMeanAndStandardDeviationOfOutputVariables(Matrix<double>) method." << std::endl
                << "Number of columns must be equal to number of outputs." << std::endl
                << std::endl;

      exit(1);
   }

   #endif

   // Set mean and standard deviation of output variables

   meanOfOutputVariables = newMeanAndStandardDeviationOfOutputVariables.getRow(0);
   
   standardDeviationOfOutputVariables = newMeanAndStandardDeviationOfOutputVariables.getRow(1);
}


// void setMeanOfSingleInputVariable(int, double) method

/// This method sets a new mean value for a single input variable.
/// These values are used for preprocessing the inputs to the neural network with the meand and standard 
/// deviation method. 
///
/// @param i Index of input variable.
/// @param newMeanOfSingleInputVariable New mean values for the input variable with index i.
///
/// @see setMeanOfInputVariables(Vector<double>).
/// @see setStandardDeviationOfInputVariables(Vector<double>).
/// @see setMeanAndStandardDeviationOfInputVariables(Matrix<double>).
/// @see setStandardDeviationOfSingleInputVariable(int, double).

void MultilayerPerceptron::setMeanOfSingleInputVariable(int i, double newMeanOfSingleInputVariable)
{
   // Control sentence (if debug)

   #ifndef NDEBUG 

   if(i >= numberOfInputs)
   {
      std::cerr << std::endl
                << "Flood Error: Multilayer Perceptron class." << std::endl
                << "void setMeanOfSingleInputVariable(int, double) method." << std::endl
                << "Index must be less than number of inputs." << std::endl
                << std::endl;

      exit(1);
   }

   #endif

   // Set mean of single input variable

   meanOfInputVariables[i] = newMeanOfSingleInputVariable;
}


// void setStandardDeviationOfSingleInputVariable(int, double) method

/// This method sets a new standard deviation value for a single input variable.
/// These values are used for preprocessing the inputs to the neural network with the meand and standard 
/// deviation method. 
///
/// @param i Index of input variable.
/// @param newStandardDeviationOfSingleInputVariable New standard deviation value for the input variable with 
/// index i.
///
/// @see setMeanOfInputVariables(Vector<double>).
/// @see setStandardDeviationOfInputVariables(Vector<double>).
/// @see setMeanAndStandardDeviationOfInputVariables(Matrix<double>).
/// @see setMeanOfSingleInputVariable(int, double).

void MultilayerPerceptron::setStandardDeviationOfSingleInputVariable(int i, 
double newStandardDeviationOfSingleInputVariable)
{
   // Control sentence (if debug)

   #ifndef NDEBUG 

   if(i >= numberOfInputs)
   {
      std::cerr << std::endl
                << "Flood Error: Multilayer Perceptron class." << std::endl
                << "void setStandardDeviationOfSingleInputVariable(int, double) method." << std::endl
                << "Index must be less than number of inputs." << std::endl
                << std::endl;

      exit(1);
   }

   #endif

   // Set standard deviation of single input variable

   standardDeviationOfInputVariables[i] = newStandardDeviationOfSingleInputVariable;
}


// void setMeanOfSingleOutputVariable(int, double) method

/// This method sets a new mean value for a single output variable.
/// These values are used for postprocessing the output signals form the neural network with the meand and 
/// standard deviation method. 
///
/// @param i Index of output variable.
/// @param newMeanOfSingleOutputVariable New mean value for the output variable with index i.
///
/// @see setMeanOfOutputVariables(Vector<double>).
/// @see setStandardDeviationOfOutputVariables(Vector<double>).
/// @see setMeanAndStandardDeviationOfOutputVariables(Matrix<double>).
/// @see setStandardDeviationOfSingleOutputVariable(int, double).

void MultilayerPerceptron::setMeanOfSingleOutputVariable(int i, double newMeanOfSingleOutputVariable)
{
   // Control sentence (if debug)

   #ifndef NDEBUG 

   if(i >= numberOfOutputs)
   {
      std::cerr << std::endl
                << "Flood Error: Multilayer Perceptron class." << std::endl
                << "void setMeanOfSingleOutputVariable(int, double) method." << std::endl
                << "Index must be less than number of outputs." << std::endl
                << std::endl;

      exit(1);
   }

   #endif

   // Set mean of single output variable

   meanOfOutputVariables[i] = newMeanOfSingleOutputVariable;
}


// void setStandardDeviationOfSingleOutputVariable(int, double) method

/// This method sets a new standard deviation value for a single output variable.
/// These values are used for postprocessing the output signals form the neural network with the meand and 
/// standard deviation method. 
///
/// @param i Index of output variable.
/// @param newStandardDeviationOfSingleOutputVariable New standard deviation value for the output variable with 
/// index i.
///
/// @see setMeanOfOutputVariables(Vector<double>).
/// @see setStandardDeviationOfOutputVariables(Vector<double>).
/// @see setMeanAndStandardDeviationOfOutputVariables(Matrix<double>).
/// @see setMeanOfSingleOutputVariable(int, double).

void MultilayerPerceptron::setStandardDeviationOfSingleOutputVariable(int i, 
double newStandardDeviationOfSingleOutputVariable)
{
   // Control sentence (if debug)

   #ifndef NDEBUG 

   if(i >= numberOfOutputs)
   {
      std::cerr << std::endl
                << "Flood Error: Multilayer Perceptron class." << std::endl
                << "void setStandardDeviationOfSingleOutputVariable(int, double) method." << std::endl
                << "Index must be less than number of outputs." << std::endl
                << std::endl;

      exit(1);
   }

   #endif

   // Set standard deviation of single output variable

   standardDeviationOfOutputVariables[i] = newStandardDeviationOfSingleOutputVariable;
}


// void setMinimumOfInputVariables(Vector<double>) method

/// This method sets new minimum values for all the input variables.
/// These values are used for preprocessing the inputs to the neural network with the minimum and maximum method. 
///
/// @param newMinimumOfInputVariables New set of minimum values for the input variables of the neural network.
///
/// @see setMaximumOfInputVariables(Vector<double>).
/// @see setMinimumAndMaximumOfInputVariables(Matrix<double>).
/// @see setMinimumOfSingleInputVariable(int, double).
/// @see setMaximumOfSingleInputVariable(int, double).

void MultilayerPerceptron::setMinimumOfInputVariables(Vector<double> newMinimumOfInputVariables)
{
   // Control sentence (if debug)

   #ifndef NDEBUG 

   if(newMinimumOfInputVariables.getSize() != numberOfInputs)
   {
      std::cerr << std::endl
                << "Flood Error: MultilayerPerceptron class." << std::endl
                << "void setMinimumOfInputVariables(Vector<double>) method." << std::endl
                << "Size must be equal to number of inputs." << std::endl
                << std::endl;

      exit(1);
   }

   #endif

   // Set minimum of input variables

   minimumOfInputVariables = newMinimumOfInputVariables;
}


// void setMaximumOfInputVariables(Vector<double>) method

/// This method sets new maximum values for all the input variables.
/// These values are used for preprocessing the inputs to the neural network with the minimum and maximum method. 
///
/// @param newMaximumOfInputVariables New set of maximum values for the input variables of the neural network.
///
/// @see setMinimumOfInputVariables(Vector<double>).
/// @see setMinimumAndMaximumOfInputVariables(Matrix<double>).
/// @see setMinimumOfSingleInputVariable(int, double).
/// @see setMaximumOfSingleInputVariable(int, double).

void MultilayerPerceptron::setMaximumOfInputVariables(Vector<double> newMaximumOfInputVariables)
{
   // Control sentence (if debug)

   #ifndef NDEBUG 

   if(newMaximumOfInputVariables.getSize() != numberOfInputs)
   {
      std::cerr << std::endl
                << "Flood Error: MultilayerPerceptron class." << std::endl
                << "void setMaximumOfInputVariables(Vector<double>) method." << std::endl
                << "Size must be equal to number of inputs." << std::endl
                << std::endl;

      exit(1);
   }
   
   #endif

   // Set maximum of input variables

   maximumOfInputVariables = newMaximumOfInputVariables;
}


// void setMinimumOfOutputVariables(Vector<double>) method

/// This method sets new minimum values for all the output variables.
/// These values are used for postprocessing the output signals from the neural network with the minimum and 
/// maximum method. 
///
/// @param newMinimumOfOutputVariables New set of minimum values for the output variables of the neural network.
///
/// @see setMaximumOfOutputVariables(Vector<double>).
/// @see setMinimumAndMaximumOfOutputVariables(Matrix<double>).
/// @see setMinimumOfSingleOutputVariable(int, double).
/// @see setMaximumOfSingleOutputVariable(int, double).

void MultilayerPerceptron::setMinimumOfOutputVariables(Vector<double> newMinimumOfOutputVariables)
{
   // Control sentence (if debug)

   #ifndef NDEBUG 

   if(newMinimumOfOutputVariables.getSize() != numberOfOutputs)
   {
      std::cerr << std::endl
                << "Flood Error: MultilayerPerceptron class." << std::endl
                << "void setMinimumOfOutputVariables(Vector<double>) method." << std::endl
                << "Size must be equal to number of outputs." << std::endl
                << std::endl;

      exit(1);
   }

   #endif

   // Set minimum of output variables

   minimumOfOutputVariables = newMinimumOfOutputVariables;
}


// void setMaximumOfOutputVariables(Vector<double>) method

/// This method sets new maximum values for all the output variables.
/// These values are used for postprocessing the output signals from the neural network with the minimum and 
/// maximum method. 
///
/// @param newMaximumOfOutputVariables New set of maximum values for the output variables of the neural network.
///
/// @see setMinimumOfOutputVariables(Vector<double>).
/// @see setMinimumAndMaximumOfOutputVariables(Matrix<double>).
/// @see setMinimumOfSingleOutputVariable(int, double).
/// @see setMaximumOfSingleOutputVariable(int, double).

void MultilayerPerceptron::setMaximumOfOutputVariables(Vector<double> newMaximumOfOutputVariables)
{
   // Control sentence (if debug)

   #ifndef NDEBUG 

   if(newMaximumOfOutputVariables.getSize() != numberOfOutputs)
   {
      std::cerr << std::endl
                << "Flood Error: MultilayerPerceptron class." << std::endl
                << "void setMaximumOfOutputVariables(Vector<double>) method." << std::endl
                << "Size must be equal to number of outputs." << std::endl
                << std::endl;

      exit(1);
   }

   #endif

   // Set maximum of output variables

   maximumOfOutputVariables = newMaximumOfOutputVariables;
}


// void setMinimumAndMaximumOfInputVariables(Matrix<double>) method

/// This method sets both the minimum and the maximum values of all the input variables from a single matrix.
/// The first row must contain the minimum values for the input variables.
/// The second row must contain the maximum values for the input variables.
/// These values are used for preprocessing the inputs to the neural network with the 
/// minimum and maximum method. 
///
/// @param newMinimumAndMaximumOfInputVariables 
/// New set of minimum and maximum values for the input variables of the neural network.
///
/// @see setMinimumOfInputVariables(Vector<double>).
/// @see setMaximumOfInputVariables(Matrix<double>).
/// @see setMinimumOfSingleInputVariable(int, double).
/// @see setMaximumOfSingleInputVariable(int, double).

void MultilayerPerceptron
::setMinimumAndMaximumOfInputVariables(Matrix<double> newMinimumAndMaximumOfInputVariables)
{
   // Control sentence (if debug)

   #ifndef NDEBUG 

   if(newMinimumAndMaximumOfInputVariables.getNumberOfRows() != 2)
   {
      std::cerr << std::endl
                << "Flood Error: MultilayerPerceptron class." << std::endl 
                << "void setMinimumAndMaximumOfInputVariables(Matrix<double>) method." << std::endl
                << "Number of rows must be 2." << std::endl
                << std::endl;

      exit(1);
   }
   else if(newMinimumAndMaximumOfInputVariables.getNumberOfColumns() != numberOfInputs)
   {
      std::cerr << std::endl
                << "Flood Error: MultilayerPerceptron class." << std::endl 
                << "void setMinimumAndMaximumOfInputVariables(Matrix<double>) method." << std::endl
                << "Number of columns must be equal to number of inputs." << std::endl
                << std::endl;

      exit(1);
   }

   // Check that minimum of input variables is not greater than their maximum

   for(int i = 0; i < numberOfInputs; i++)
   {
      if(newMinimumAndMaximumOfInputVariables[0][i] >= newMinimumAndMaximumOfInputVariables[1][i])
      {
         std::cerr << std::endl
                   << "Flood Error: MultilayerPerceptron class." << std::endl
                   << "void setMinimumAndMaximumOfInputVariables(Matrix<double>) method." << std::endl
                   << "Minimum of input variable " << i << " is equal or greater than maximum of that variable."
                   << std::endl << std::endl;

         exit(1);
      }
   }

   #endif

   // Set minimum and maximum of input variables

   minimumOfInputVariables = newMinimumAndMaximumOfInputVariables.getRow(0);

   maximumOfInputVariables = newMinimumAndMaximumOfInputVariables.getRow(1);
}


// void setMinimumAndMaximumOfOutputVariables(Matrix<double>) method

/// This method sets both the minimum and the maximum values of all the output variables from a single matrix.
/// The first row must contain the minimum values for the output variables.
/// The second row must contain the maximum values for the output variables.
/// These values are used for postprocessing the output signals from network with the 
/// minimum and maximum method. 
///
/// @param newMinimumAndMaximumOfOutputVariables 
/// New set of minimum and maximum values for the output variables of the neural network.
///
/// @see setMinimumOfOutputVariables(Vector<double>).
/// @see setMaximumOfOutputVariables(Matrix<double>).
/// @see setMinimumOfSingleOutputVariable(int, double).
/// @see setMaximumOfSingleOutputVariable(int, double).

void MultilayerPerceptron
::setMinimumAndMaximumOfOutputVariables(Matrix<double> newMinimumAndMaximumOfOutputVariables)
{
   // Control sentence (if debug)

   #ifndef NDEBUG 

   if(newMinimumAndMaximumOfOutputVariables.getNumberOfRows() != 2)
   {
      std::cerr << std::endl
                << "Flood Error: MultilayerPerceptron class." << std::endl 
                << "void setMinimumAndMaximumOfOutputVariables(Matrix<double>) method." << std::endl
                << "Number of rows must be 2."  << std::endl
                << std::endl;

      exit(1);
   }
   else if(newMinimumAndMaximumOfOutputVariables.getNumberOfColumns() != numberOfOutputs)
   {
      std::cerr << std::endl
                << "Flood Error: MultilayerPerceptron class." << std::endl 
                << "void setMinimumAndMaximumOfOutputVariables(Matrix<double>) method." << std::endl
                << "Number of columns must be equal to number of outputs." << std::endl
                << std::endl;

      exit(1);
   }

   // Check that minimum of output variables is not greater than their maximum

   for(int i = 0; i < numberOfOutputs; i++)
   {
      if(newMinimumAndMaximumOfOutputVariables[0][i] > newMinimumAndMaximumOfOutputVariables[1][i])
      {
         std::cerr << std::endl
                   << "Flood Error: MultilayerPerceptron class." << std::endl 
                   << "void setMinimumAndMaximumOfOutputVariables(Matrix<double>) method." << std::endl
                   << "Minimum of output variable " << i << " is greater than maximum of that variable."
                   << std::endl << std::endl;

         exit(1);
      }
   }

   #endif

   // Set minimum and maximum of output variables

   minimumOfOutputVariables = newMinimumAndMaximumOfOutputVariables.getRow(0);

   maximumOfOutputVariables = newMinimumAndMaximumOfOutputVariables.getRow(1);
}


// void setMinimumOfSingleInputVariable(int, double) method

/// This method sets a new minimum value for a single input variable.
/// This value is used for preprocessing that input to the neural network with the minimum and maximum method. 
///
/// @param i Index of input variable.
/// @param newMinimumOfSingleInputVariable New minimum value for the input variable with index i.
///
/// @see setMinimumOfInputVariables(Vector<double>).
/// @see setMaximumOfInputVariables(Vector<double>).
/// @see setMinimumAndMaximumOfInputVariables(Matrix<double>).
/// @see setMaximumOfSingleInputVariable(int, double).

void MultilayerPerceptron::setMinimumOfSingleInputVariable(int i, double newMinimumOfSingleInputVariable)
{
   // Control sentence (if debug)

   #ifndef NDEBUG 

   if(i >= numberOfInputs)
   {
      std::cerr << std::endl
                << "Flood Error: MultilayerPerceptron class." << std::endl 
                << "void setMinimumOfSingleInputVariable(int, double) method." << std::endl
                << "Index must be less than number of inputs." << std::endl
                << std::endl;

      exit(1);
   }

   #endif

   // Set minimum of single input variable

   minimumOfInputVariables[i] = newMinimumOfSingleInputVariable;
}


// void setMaximumOfSingleInputVariable(int, double) method

/// This method sets a new maximum value for a single input variable.
/// This value is used for preprocessing that input to the neural network with the minimum and maximum method. 
///
/// @param i Index of input variable.
/// @param newMaximumOfSingleInputVariable New maximum value for the input variable with index i.
///
/// @see setMinimumOfInputVariables(Vector<double>).
/// @see setMaximumOfInputVariables(Vector<double>).
/// @see setMinimumAndMaximumOfInputVariables(Matrix<double>).
/// @see setMinimumOfSingleInputVariable(int, double).

void MultilayerPerceptron::setMaximumOfSingleInputVariable(int i, double newMaximumOfSingleInputVariable)
{
   // Control sentence (if debug)

   #ifndef NDEBUG 

   if(i >= numberOfInputs)
   {
      std::cerr << std::endl
                << "Flood Error: MultilayerPerceptron class." << std::endl 
                << "void setMaximumOfSingleInputVariable(int, double) method." << std::endl
                << "Index must be less than number of inputs." << std::endl
                << std::endl;

      exit(1);
   }

   #endif

   // Set maximum of single input variable

   maximumOfInputVariables[i] = newMaximumOfSingleInputVariable;
}


// void setMinimumOfSingleOutputVariable(int, double) method

/// This method sets a new minimum value for a single output variable.
/// This value is used for postprocessing that output signals from the neural network with the minimum and 
/// maximum method. 
///
/// @param i Index of output variable.
/// @param newMinimumOfSingleOutputVariable New minimum value for the output variable with index i.
///
/// @see setMinimumOfOutputVariables(Vector<double>).
/// @see setMaximumOfOutputVariables(Vector<double>).
/// @see setMinimumAndMaximumOfOutputVariables(Matrix<double>).
/// @see setMaximumOfSingleOutputVariable(int, double).

void MultilayerPerceptron::setMinimumOfSingleOutputVariable(int i,
double newMinimumOfSingleOutputVariable)
{
   // Control sentence (if debug)

   #ifndef NDEBUG 

   if(i >= numberOfOutputs)
   {
      std::cerr << std::endl
                << "Flood Error: MultilayerPerceptron class." << std::endl 
                << "void setMinimumOfSingleOutputVariable(int, double) method." << std::endl
                << "Index must be less than number of outputs." << std::endl
                << std::endl;

      exit(1);
   }

   #endif

   // Set minimum of single output variable

   minimumOfOutputVariables[i] = newMinimumOfSingleOutputVariable;
}


// void setMaximumOfSingleOutputVariable(int, double) method

/// This method sets a new maximum value for a single output variable.
/// This value is used for postprocessing that output signals from the neural network with the minimum and 
/// maximum method. 
///
/// @param i Index of output variable.
/// @param newMaximumOfSingleOutputVariable New maximum value for the output variable with index i.
///
/// @see setMinimumOfOutputVariables(Vector<double>).
/// @see setMaximumOfOutputVariables(Vector<double>).
/// @see setMinimumAndMaximumOfOutputVariables(Matrix<double>).
/// @see setMinimumOfSingleOutputVariable(int, double).

void MultilayerPerceptron::setMaximumOfSingleOutputVariable(int i, 
double newMaximumOfSingleOutputVariable)
{
   // Control sentence (if debug)

   #ifndef NDEBUG 

   if(i >= numberOfOutputs)
   {
      std::cerr << std::endl
                << "Flood Error: MultilayerPerceptron class." << std::endl 
                << "void setMaximumOfSingleOutputVariable(int, double) method." << std::endl
                << "Index must be less than number of outputs." << std::endl
                << std::endl;

      exit(1);
   }

   #endif

   // Set maximum of single output variable

   maximumOfOutputVariables[i] = newMaximumOfSingleOutputVariable;
}


// void setAllStatistics(Vector< Vector<double> >) method

void MultilayerPerceptron::setAllStatistics(Vector< Vector<double> > newStatistics)
{
   // Control sentence (if debug)

   #ifndef NDEBUG 

   int size = newStatistics.getSize();

   if(size != 8)
   {
      std::cerr << std::endl
                << "Flood Error: MultilayerPerceptron class." << std::endl
                << "void setAllStatistics(Vector< Vector<double> >) method." << std::endl
                << "Size must be 8." << std::endl
                << std::endl;

      exit(1);
   }

   #endif

   // Set all statistics
 
   setMeanOfInputVariables(newStatistics[0]);
   setStandardDeviationOfInputVariables(newStatistics[1]);

   setMeanOfOutputVariables(newStatistics[2]);
   setStandardDeviationOfOutputVariables(newStatistics[3]);

   setMinimumOfInputVariables(newStatistics[4]);
   setMaximumOfInputVariables(newStatistics[5]);

   setMinimumOfOutputVariables(newStatistics[6]);
   setMaximumOfOutputVariables(newStatistics[7]);
}
						  

// void setLowerBoundOfOutputVariables(Vector<double>) method

/// This method sets new lower bounds for all the output variables.
/// These values are used for postprocessing the outputs so that they are not less than the lower bounds. 
///
/// @param newLowerBoundOfOutputVariables New set of lower bounds for the output variables of the neural network.
///
/// @see setUpperBoundOfOutputVariables(Vector<double>).
/// @see setLowerAndUpperBoundOfOutputVariables(Vector<double>).
/// @see setLowerBoundOfSingleOutputVariable(int, double).
/// @see setUpperBoundOfSingleOutputVariable(int, double).

void MultilayerPerceptron::setLowerBoundOfOutputVariables(Vector<double> newLowerBoundOfOutputVariables)
{
   // Control sentence (if debug)

   #ifndef NDEBUG 

   if(newLowerBoundOfOutputVariables.getSize() != numberOfOutputs)
   {
      std::cerr << std::endl
                << "Flood Error: MultilayerPerceptron class." << std::endl
                << "void setLowerBoundOfOutputVariables(Vector<double>) method." << std::endl
                << "Size must be equal to number of outputs." << std::endl
                << std::endl;

      exit(1);
   }

   #endif

   // Set lower bound of output variables

   lowerBoundOfOutputVariables = newLowerBoundOfOutputVariables;
}


// void setUpperBoundOfOutputVariables(Vector<double>) method

/// This method sets new upper bounds for all the output variables.
/// These values are used for postprocessing the outputs so that they are not greater than the upper bounds. 
///
/// @param newUpperBoundOfOutputVariables New set of upper bounds for the output variables of the neural network.
///
/// @see setLowerBoundOfOutputVariables(Vector<double>).
/// @see setLowerAndUpperBoundOfOutputVariables(Vector<double>).
/// @see setLowerBoundOfSingleOutputVariable(int, double).
/// @see setUpperBoundOfSingleOutputVariable(int, double).

void MultilayerPerceptron::setUpperBoundOfOutputVariables(Vector<double> newUpperBoundOfOutputVariables)
{
   // Control sentence (if debug)

   #ifndef NDEBUG 

   if(newUpperBoundOfOutputVariables.getSize() != numberOfOutputs)
   {
      std::cerr << std::endl
                << "Flood Error: MultilayerPerceptron class." << std::endl
                << "void setUpperBoundOfOutputVariables(Vector<double>) method." << std::endl
                << "Size must be equal to number of outputs." << std::endl
                << std::endl;

      exit(1);
   }

   #endif

   // Set upper bound of output variables

   upperBoundOfOutputVariables = newUpperBoundOfOutputVariables;
}


// void setLowerAndUpperBoundsOfOutputVariables(Matrix<double>) method

/// This method sets both the lower bounds and the upper bounds of all the output variables from a single matrix.
/// The first row must contain the lower bound values for the output variables.
/// The second row must contain the upper bound values for the output variables.
/// These values are used for postprocessing the outputs so that they are not greater than 
/// the upper bounds. 
///
/// @param newLowerAndUpperBoundsOfOutputVariables New set of lower and upper bounds for the output variables of the 
/// neural network.
///
/// @see setLowerBoundOfOutputVariables(Vector<double>).
/// @see setUpperBoundOfOutputVariables(Vector<double>).
/// @see setLowerBoundOfSingleOutputVariable(int, double).
/// @see setUpperBoundOfSingleOutputVariable(int, double).

void MultilayerPerceptron
::setLowerAndUpperBoundsOfOutputVariables(Matrix<double> newLowerAndUpperBoundsOfOutputVariables)
{
   // Control sentence (if debug)

   #ifndef NDEBUG 

   int numberOfRows = newLowerAndUpperBoundsOfOutputVariables.getNumberOfRows();
   int numberOfColumns = newLowerAndUpperBoundsOfOutputVariables.getNumberOfColumns();

   if(numberOfRows != 2)
   {
      std::cerr << std::endl
                << "Flood Error: MultilayerPerceptron class." << std::endl
                << "void setLowerAndUpperBoundsOfOutputVariables(Matrix<double>) method." << std::endl
                << "Number of rows must be 2." << std::endl
                << std::endl;

      exit(1);
   }
   else if(numberOfColumns != numberOfOutputs)
   {
      std::cerr << std::endl
                << "Flood Error: MultilayerPerceptron class." << std::endl
                << "void setLowerAndUpperBoundsOfOutputVariables(Matrix<double>) method." << std::endl
                << "Number of columns must be equal to number of outputs." << std::endl
                << std::endl;

      exit(1);
   }

   #endif

   // Set lower and upper bounds of output variables

   lowerBoundOfOutputVariables = newLowerAndUpperBoundsOfOutputVariables.getRow(0);
   upperBoundOfOutputVariables = newLowerAndUpperBoundsOfOutputVariables.getRow(1);
}


// void setLowerBoundOfSingleOutputVariable(int, double) method

/// This method sets a new lower bound for a single output variable.
/// This value is used for postprocessing that output so that it is not less than the lower bound. 
///
/// @param i Index of output variable.
/// @param newLowerBoundOfSingleOutputVariable New lower bound for the output variable with index i.
///
/// @see setLowerBoundOfOutputVariables(Vector<double>).
/// @see setUpperBoundOfOutputVariables(Vector<double>).
/// @see setLowerAndUpperBoundOfOutputVariables(Vector<double>).
/// @see setUpperBoundOfSingleOutputVariable(int, double).

void MultilayerPerceptron::setLowerBoundOfSingleOutputVariable(int i, 
double newLowerBoundOfSingleOutputVariable)
{
   // Control sentence (if debug)

   #ifndef NDEBUG 

   if(i >= numberOfOutputs)
   {
      std::cerr << std::endl
                << "Flood Error: MultilayerPerceptron class." << std::endl
                << "void setLowerBoundOfSingleOutputVariable(int, double) method." << std::endl
                << "Index must be less than number of outputs." << std::endl
                << std::endl;

      exit(1);
   }

   #endif

   // Set lower bound of single output variable

   lowerBoundOfOutputVariables[i] = newLowerBoundOfSingleOutputVariable;
}


// void setUpperBoundOfSingleOutputVariable(int, double) method

/// This method sets a new upper bound for a single output variable.
/// This value is used for postprocessing that output so that it is not greater than the upper bound. 
///
/// @param i Index of output variable.
/// @param newUpperBoundOfSingleOutputVariable New upper bound for the output variable with index i.
///
/// @see setLowerBoundOfOutputVariables(Vector<double>).
/// @see setUpperBoundOfOutputVariables(Vector<double>).
/// @see setLowerAndUpperBoundOfOutputVariables(Vector<double>).
/// @see setLowerBoundOfSingleOutputVariable(int, double).

void MultilayerPerceptron::setUpperBoundOfSingleOutputVariable(int i,
double newUpperBoundOfSingleOutputVariable)
{
   // Control sentence (if debug)

   #ifndef NDEBUG 

   if(i >= numberOfOutputs)
   {
      std::cerr << std::endl
                << "Flood Error: MultilayerPerceptron class." << std::endl
                << "void setUpperBoundOfSingleOutputVariable(int, double) method." << std::endl
                << "Index must be less than number of outputs." << std::endl
                << std::endl;

      exit(1);
   }

   #endif

   // Set upper bound of single output variable

   upperBoundOfOutputVariables[i] = newUpperBoundOfSingleOutputVariable;
}


// void setNumberOfIndependentParameters(int) method

/// This method sets a new number of free parameters independent of the neural network. 
/// It initializes the independent parameters with random values chosen from a normal distribution with mean zero
/// and standard deviation one.
/// This method also initializes all the class members related to independent with their default value:
///
/// <ul>
/// <li> Name of independent parameters: "IndependentParameter0", ... 
/// <li> Units of independent parameters: "None". 
/// <li> Description of independent parameters: "None". 
/// <li> Mean of independent parameters: 0.
/// <li> Standard deviation of independent parameters: 1.
/// <li> Minimum of independent parameters: -Inf.
/// <li> Maximum of independent parameters: +Inf.
/// <li> Lower bound of independent parameters: -Inf.
/// <li> Upper bound of independent parameters: +Inf.
/// </ul>

void MultilayerPerceptron::setNumberOfIndependentParameters(int newNumberOfIndependentParameters)
{
   numberOfIndependentParameters = newNumberOfIndependentParameters;

   // Independent parameters

   independentParameters.setSize(numberOfIndependentParameters);

   // Name of independent parameters

   nameOfIndependentParameters.setSize(numberOfIndependentParameters);

   for(int i = 0; i < numberOfIndependentParameters; i++)
   {
      std::stringstream buffer;

      buffer << "IndependentParameter" << i;

      nameOfIndependentParameters[i] = buffer.str();
   }

   // Units of independent parameters

   unitsOfIndependentParameters.setSize(numberOfIndependentParameters);

   for(int i = 0; i < numberOfIndependentParameters; i++)
   {
      unitsOfIndependentParameters[i] = "None";
   }

   // Description of independent parameters

   descriptionOfIndependentParameters.setSize(numberOfIndependentParameters);

   for(int i = 0; i < numberOfIndependentParameters; i++)
   {
      descriptionOfIndependentParameters[i] = "None";
   }

   // Mean of independent parameters

   meanOfIndependentParameters.setSize(numberOfIndependentParameters);

   for(int i = 0; i < numberOfIndependentParameters; i++)
   {
      meanOfIndependentParameters[i] = 0.0;
   }
      
   // Standard deviation of independent parameters

   standardDeviationOfIndependentParameters.setSize(numberOfIndependentParameters);

   for(int i = 0; i < numberOfIndependentParameters; i++)
   {
      standardDeviationOfIndependentParameters[i] = 1.0;
   }

   // Minimum of independent parameters

   minimumOfIndependentParameters.setSize(numberOfIndependentParameters);

   for(int i = 0; i < numberOfIndependentParameters; i++)
   {
      minimumOfIndependentParameters[i] = -1.0e-99;
   }
      
   // Maximum of independent parameters

   maximumOfIndependentParameters.setSize(numberOfIndependentParameters);

   for(int i = 0; i < numberOfIndependentParameters; i++)
   {
      maximumOfIndependentParameters[i] = 1.0e99;
   }

   // Lower bound of independent parameters

   lowerBoundOfIndependentParameters.setSize(numberOfIndependentParameters);

   for(int i = 0; i < numberOfIndependentParameters; i++)
   {
      lowerBoundOfIndependentParameters[i] = -1.0e99;
   }
      
   // Upper bound of independent parameters

   upperBoundOfIndependentParameters.setSize(numberOfIndependentParameters);

   for(int i = 0; i < numberOfIndependentParameters; i++)
   {
      upperBoundOfIndependentParameters[i] = 1.0e99;
   }

   // Initialize independent parameters 
  
   initIndependentParametersNormal(0.0,1.0);
}


// void setNameOfIndependentParameters(Vector<std::string>) method

/// This method sets new names for the independent parameters.
/// Such values are only used to give the user basic information on the problem at hand.
///
/// @param newNameOfIndependentParameters New names for the independent parameters of the neural network.
///
/// @ see setNameOfSingleIndependentParameter(int, std::string).

void MultilayerPerceptron::setNameOfIndependentParameters(Vector<std::string> newNameOfIndependentParameters)
{
   // Control sentence (if debug)

   #ifndef NDEBUG 
   
   if(newNameOfIndependentParameters.getSize() != numberOfIndependentParameters)
   {
      std::cerr << std::endl
                << "Flood Error: MultilayerPerceptron class." << std::endl
                << "void setNameOfIndependentParameters(Vector<std::string>) method."  << std::endl
                << "Name of independent parameters vector size must be equal to number of independent parameters."
                << std::endl << std::endl;

      exit(1);
   }

   #endif

   // Set name of independent parameters

   nameOfIndependentParameters = newNameOfIndependentParameters;
}


// void setNameOfSingleIndependentParameter(int, std::string) method

/// This method sets a new name for a single independent parameter.
/// Such a value is only used to give the user basic information on the problem at hand.
///
/// @param i Index of independent parameter.
/// @param newNameOfSingleIndependentParameter New name for the independent parameter of index i.
///
/// @ see setNameOfIndependentParameters(Vector<std::string>).

void MultilayerPerceptron::setNameOfSingleIndependentParameter(int i, 
std::string newNameOfSingleIndependentParameter)
{
   // Control sentence (if debug)

   #ifndef NDEBUG 

   if(i >= numberOfIndependentParameters)
   {
      std::cerr << std::endl
                << "Flood Error: MultilayerPerceptron class." << std::endl
                << "void setNameOfSingleIndependentParameter(int, std::string) method." << std::endl
                << "Index must be less than number of independent parameters." << std::endl 
                << std::endl;

      exit(1);   
   }

   #endif

   // Set name of single independent parameter

   nameOfIndependentParameters[i] = newNameOfSingleIndependentParameter;
}


// void setUnitsOfIndependentParameters(Vector<std::string>) method

/// This method sets new units for the independent parameters.
/// Such values are only used to give the user basic information on the problem at hand.
///
/// @param newUnitsOfIndependentParameters New units for the independent parameters of the neural network.
///
/// @ see setUnitsOfSingleIndependentParameter(int, std::string).

void MultilayerPerceptron::setUnitsOfIndependentParameters(Vector<std::string> newUnitsOfIndependentParameters)
{
   // Control sentence (if debug)

   #ifndef NDEBUG 

   if(newUnitsOfIndependentParameters.getSize() != numberOfIndependentParameters)
   {
      std::cerr << std::endl
                << "Flood Error: MultilayerPerceptron class." << std::endl 
                << "void setUnitsOfIndependentParameters(Vector<std::string>) method." << std::endl
                << "Size must be equal to number of independent parameters." << std::endl 
                << std::endl;

      exit(1);
   }

   #endif

   // Set units of independent parameters

   unitsOfIndependentParameters = newUnitsOfIndependentParameters;
}


// void setUnitsOfSingleIndependentParameter(int, std::string) method

/// This method sets new units for a single independent parameter.
/// Such a value is only used to give the user basic information on the problem at hand.
///
/// @param i Index of independent parameter.
/// @param newUnitsOfSingleIndependentParameter New units for the independent parameter of index i.
///
/// @ see setUnitsOfIndependentParameters(Vector<std::string>).

void MultilayerPerceptron::setUnitsOfSingleIndependentParameter(int i, 
std::string newUnitsOfSingleIndependentParameter)
{
   // Control sentence (if debug)

   #ifndef NDEBUG 

   if(i >= numberOfIndependentParameters)
   {
      std::cerr << std::endl
                << "Flood Error: MultilayerPerceptron class." << std::endl
                << "void setUnitsOfSingleIndependentParameter(int, std::string) method." << std::endl
                << "Index must be less than number of independent parameters." << std::endl 
                << std::endl;

      exit(1);   
   }

   #endif

   // Set units of single independent parameter

   unitsOfIndependentParameters[i] = newUnitsOfSingleIndependentParameter;
}


// void setDescriptionOfIndependentParameters(Vector<std::string>) method

/// This method sets new descriptions for the independent parameters.
/// Such values are only used to give the user basic information on the problem at hand.
///
/// @param newDescriptionOfIndependentParameters New description for the independent parameters of the neural 
/// network.
///
/// @ see setDescriptionOfSingleIndependentParameter(int, std::string).

void MultilayerPerceptron
::setDescriptionOfIndependentParameters(Vector<std::string> newDescriptionOfIndependentParameters)
{
   // Control sentence (if debug)

   #ifndef NDEBUG 

   if(newDescriptionOfIndependentParameters.getSize() != numberOfIndependentParameters)
   {
      std::cerr << std::endl
                << "Flood Error: MultilayerPerceptron class." << std::endl 
                << "void setDescriptionOfIndependentParameters(Vector<std::string>) method." << std::endl
                << "Size must be equal to number of independent parameters." << std::endl 
                << std::endl;

      exit(1);
   }

   #endif

   // Set description of independent parameters

   descriptionOfIndependentParameters = newDescriptionOfIndependentParameters;
}


// void setDescriptionOfSingleIndependentParameter(int, std::string) method

/// This method sets a new description for a single independent parameter.
/// Such a value is only used to give the user basic information on the problem at hand.
///
/// @param i Index of independent parameter.
/// @param newDescriptionOfSingleIndependentParameter New description for the independent parameter of index i.
///
/// @ see setDescriptionOfIndependentParameters(Vector<std::string>).

void MultilayerPerceptron::setDescriptionOfSingleIndependentParameter(int i, 
std::string newDescriptionOfSingleIndependentParameter)
{
   // Control sentence (if debug)

   #ifndef NDEBUG 

   if(i >= numberOfIndependentParameters)
   {
      std::cerr << std::endl
                << "Flood Error: MultilayerPerceptron class." << std::endl
                << "void setDescriptionOfSingleIndependentParameter(int, std::string) method." 
                << std::endl
                << "Index must be less than number of independent parameters." 
                << std::endl << std::endl;

      exit(1);   
   }

   #endif

   // Set description of single independent parameter
   
   descriptionOfIndependentParameters[i] = newDescriptionOfSingleIndependentParameter;
}


// void setIndependentParameters(Vector<double>) method

/// This method sets new values for all the independent parameters.
///
/// @param newIndependentParameters Independent parameters values.
///
/// @see setSingleIndependentParameter(int, double).

void MultilayerPerceptron::setIndependentParameters(Vector<double> newIndependentParameters)
{
   // Control sentence (if debug)

   #ifndef NDEBUG 

   if(newIndependentParameters.getSize() != numberOfIndependentParameters)
   {
      std::cerr << std::endl
                << "Flood Error: MultilayerPerceptron class." << std::endl
                << "void setIndependentParametersSet(Vector<double>) method." << std::endl
                << "Independent parameters vector size must be equal to number of independent parameters."
                << std::endl << std::endl;

      exit(1);
   }

   #endif

   // Set independent parameters

   independentParameters = newIndependentParameters;     
}


// void setSingleIndependentParameter(int, double) method

/// This method sets a new value for a single independent parameter.
///
/// @param i Independent parameter index.
/// @param newSingleIndependentParameter Independent parameter value.
///
/// @see setIndependentParameters(Vector<double>).

void MultilayerPerceptron::setSingleIndependentParameter(int i, double newSingleIndependentParameter)
{
   // Control sentence (if debug)

   #ifndef NDEBUG 

   if(i >= numberOfIndependentParameters)
   {
      std::cerr << std::endl
                << "Flood Error: MultilayerPerceptron class." << std::endl
                << "void setSingleIndependentParameter(int, double) method." << std::endl
                << "Index must be less than number of independent parameters." << std::endl 
                << std::endl;

      exit(1);   
   }

   #endif

   // Set single independent parameter

   independentParameters[i] = newSingleIndependentParameter;
}


// void setMeanOfIndependentParameters(Vector<double>) method

/// This method sets the mean values of all the independent parameters.
/// These values are used for pre and postprocessing the independent parameters with the mean and standard 
/// deviation method. 
///
/// @param newMeanOfIndependentParameters New set of mean values for the independent parameters of the neural 
/// network.
///
/// @see setStandardDeviationOfIndependentParameters(Vector<double>).
/// @see setMeanOfSingleIndependentParameter(int, double).
/// @see setStandardDeviationOfSingleIndependentParameter(int, double).

void MultilayerPerceptron::setMeanOfIndependentParameters(Vector<double> newMeanOfIndependentParameters)
{
   // Control sentence (if debug)

   #ifndef NDEBUG 

   if(newMeanOfIndependentParameters.getSize() != numberOfIndependentParameters)
   {
      std::cerr << std::endl
                << "Flood Error: MultilayerPerceptron class." << std::endl 
                << "void setMeanOfIndependentParameters(Vector<double>) method." << std::endl
                << "Size must be equal to number of independent parameters." << std::endl 
                << std::endl;

      exit(1);
   }

   #endif

   // Set mean of independent parameters

   meanOfIndependentParameters = newMeanOfIndependentParameters;                                                   
}


// void setStandardDeviationOfIndependentParameters(Vector<double>) method

/// This method sets the standard deviation values of all the independent parameters.
/// These values are used for pre and postprocessing the independent parameters with the mean and standard 
/// deviation method. 
///
/// @param newStandardDeviationOfIndependentParameters New set of standard deviation values for the 
/// independent parameters of the neural network.
///
/// @see setMeanOfIndependentParameters(Vector<double>).
/// @see setMeanOfSingleIndependentParameter(int, double).
/// @see setStandardDeviationOfSingleIndependentParameter(int, double).

void MultilayerPerceptron
::setStandardDeviationOfIndependentParameters(Vector<double> newStandardDeviationOfIndependentParameters)
{
   // Control sentence (if debug)

   #ifndef NDEBUG 

   if(newStandardDeviationOfIndependentParameters.getSize() != numberOfIndependentParameters)
   {
      std::cerr << std::endl
                << "Flood Error: MultilayerPerceptron class." << std::endl 
                << "void setStandardDeviationOfIndependentParameters(Vector<double>) method." << std::endl
                << "Size must be equal to number of independent parameters." << std::endl 
                << std::endl;

      exit(1);
   }

   #endif

   // Set standard deviation of independent parameters

   standardDeviationOfIndependentParameters = newStandardDeviationOfIndependentParameters;  
}


// void setMeanAndStandardDeviationOfIndependentParameters(Matrix<double>) method

/// This method sets both the mean and the standard deviation values of all the independent parameters from a 
/// single matrix.
/// The first row must contain the mean values values for the independent parameters.
/// The second row must contain the standard deviation values for the independent parameters.
/// These values are used for pre and postprocessing the independent parameters with the mean and standard 
/// deviation method. 
///
/// @param newMeanAndStandardDeviationOfIndependentParameters 
/// New set of mean and standard deviation values for the independent parameters of the neural network.
///
/// @see setMeanOfIndependentParameters(Vector<double>).
/// @see setStandardDeviationOfIndependentParameters(Vector<double>).

void MultilayerPerceptron::setMeanAndStandardDeviationOfIndependentParameters
(Matrix<double> newMeanAndStandardDeviationOfIndependentParameters)
{
   // Control sentence (if debug)

   #ifndef NDEBUG 

   int numberOfRows = newMeanAndStandardDeviationOfIndependentParameters.getNumberOfRows();
   int numberOfColumns = newMeanAndStandardDeviationOfIndependentParameters.getNumberOfRows();

   if(numberOfRows != 2)
   {
      std::cerr << std::endl
                << "Flood Error: MultilayerPerceptron class." << std::endl
                << "void setMeanAndStandardDeviationOfIndependentParameters(Matrix<double>) method." << std::endl
                << "Number of rows must be 2." << std::endl 
                << std::endl;

      exit(1);
   }
   else if(numberOfColumns != numberOfIndependentParameters)
   {
      std::cerr << std::endl
                << "Flood Error: MultilayerPerceptron class." << std::endl
                << "void setMeanAndStandardDeviationOfIndependentParameters(Matrix<double>) method." << std::endl
                << "Number of columns must be equal to number of independent parameters." << std::endl
                << std::endl;

      exit(1);
   }

   // Check that standard deviation of input variables is not zero

   if(display)
   {
      for(int i = 0; i < numberOfInputs; i++)
      {
         if(newMeanAndStandardDeviationOfIndependentParameters[1][i] < 1.0e-69)
         {
            std::cerr << std::endl
                      << "Flood Warning: MultilayerPerceptron class: " << std::endl
                      << "void setMeanAndStandardDeviationOfIndependentParameters(Matrix<double>) method."
                      << std::endl
                      << "Standard deviation of independent parameter " << i << " is zero." << std::endl 
                      << std::endl;
         }
      }
   }

   #endif

   // Set mean and standard deviation of independent parameters

   meanOfIndependentParameters = newMeanAndStandardDeviationOfIndependentParameters.getRow(0);
   standardDeviationOfIndependentParameters = newMeanAndStandardDeviationOfIndependentParameters.getRow(1);
}


// void setMeanOfSingleIndependentParameter(int, double) method

/// This method sets a new mean value for a single independent parameter.
/// Such a value is used for pre and postprocessing the independent parameters with the mean and standard 
/// deviation method. 
///
/// @param i Index of independent parameter.
/// @param newMeanOfSingleIndependentParameter New mean value for the independent parameter of index i.
///
/// @see setMeanOfIndependentParameters(Vector<double>).
/// @see setStandardDeviationOfIndependentParameters(Vector<double>).
/// @see setStandardDeviationOfSingleIndependentParameter(int, double).

void MultilayerPerceptron::setMeanOfSingleIndependentParameter(int i, double newMeanOfSingleIndependentParameter)
{
   // Control sentence (if debug)

   #ifndef NDEBUG 

   if(i >= numberOfIndependentParameters)
   {
      std::cerr << std::endl
                << "Flood Error: MultilayerPerceptron class." << std::endl
                << "void setMeanOfSingleIndependentParameter(int, double) method." << std::endl
                << "Index must be less than number of independent parameters." << std::endl 
                << std::endl;

      exit(1);   
   }

   #endif

   // Set mean of single independent parameter

   meanOfIndependentParameters[i] = newMeanOfSingleIndependentParameter;
}


// void setStandardDeviationOfSingleIndependentParameter(int, double) method

/// This method sets a new standard deviation value for a single independent parameter.
/// Such a value is used for pre and postprocessing the independent parameters with the mean and standard 
/// deviation method. 
///
/// @param i Index of independent parameter.
/// @param newStandardDeviationOfSingleIndependentParameter New standard deviation value for the 
/// independent parameter of index i.
///
/// @see setMeanOfIndependentParameters(Vector<double>).
/// @see setStandardDeviationOfIndependentParameters(Vector<double>).
/// @see setMeanOfSingleIndependentParameter(int, double).

void MultilayerPerceptron::setStandardDeviationOfSingleIndependentParameter(int i, 
double newStandardDeviationOfSingleIndependentParameter)
{
   // Control sentence (if debug)

   #ifndef NDEBUG 

   if(i >= numberOfIndependentParameters)
   {
      std::cerr << std::endl
                << "Flood Error: MultilayerPerceptron class." << std::endl
                << "void setStandardDeviationOfSingleIndependentParameter(int, double) method." 
                << std::endl
                << "Index must be less than number of independent parameters." 
                << std::endl << std::endl;

      exit(1);   
   }

   #endif

   // Set standard deviation of single independent parameter

   standardDeviationOfIndependentParameters[i] = newStandardDeviationOfSingleIndependentParameter;
}


// void setMinimumOfIndependentParameters(Vector<double>) method

/// This method sets the minimum values of all the independent parameters.
/// These values are used for pre and postprocessing the independent parameters with the minimum and maximum 
/// method.
///
/// @param newMinimumOfIndependentParameters New set of minimum values for the independent parameters of the 
/// neural network.
///
/// @see setMaximumOfIndependentParameters(Vector<double>).
/// @see setMinimumAndMaximumOfIndependentParameters(Matrix<double>).
/// @see setMinimumOfSingleIndependentParameter(int, double).
/// @see setMaximumOfSingleIndependentParameter(int, double).

void MultilayerPerceptron::setMinimumOfIndependentParameters(Vector<double> newMinimumOfIndependentParameters)
{
   // Control sentence (if debug)

   #ifndef NDEBUG 

   if(newMinimumOfIndependentParameters.getSize() != numberOfIndependentParameters)
   {
      std::cerr << std::endl
                << "Flood Error: MultilayerPerceptron class." << std::endl
                << "void setMinimumOfIndependentParameters(Vector<double>) method." << std::endl
                << "Size must be equal to number of independent parameters." << std::endl 
                << std::endl;

      exit(1);
   }

   #endif

   // Set minimum of independent parameters

   minimumOfIndependentParameters = newMinimumOfIndependentParameters;                                                   
}


// void setMaximumOfIndependentParameters(Vector<double>) method

/// This method sets the maximum values of all the independent parameters.
/// These values are used for pre and postprocessing the independent parameters with the minimum and maximum 
/// method.
///
/// @param newMaximumOfIndependentParameters New set of maximum values for the independent parameters of the 
/// neural network.
///
/// @see setMinimumOfIndependentParameters(Vector<double>).
/// @see setMinimumAndMaximumOfIndependentParameters(Matrix<double>).
/// @see setMinimumOfSingleIndependentParameter(int, double).
/// @see setMaximumOfSingleIndependentParameter(int, double).

void MultilayerPerceptron::setMaximumOfIndependentParameters(Vector<double> newMaximumOfIndependentParameters)
{
   // Control sentence (if debug)

   #ifndef NDEBUG 

   if(newMaximumOfIndependentParameters.getSize() != numberOfIndependentParameters)
   {
      std::cerr << std::endl
                << "Flood Error: MultilayerPerceptron class." << std::endl
                << "void setMaximumOfIndependentParameters(Vector<double>) method." << std::endl
                << "Size must be equal to number of independent parameters." << std::endl 
                << std::endl;

      exit(1);
   }

   #endif

   // Set maximum of independent parameters

   maximumOfIndependentParameters = newMaximumOfIndependentParameters;  
}


// void setMinimumAndMaximumOfIndependentParameters(Matrix<double>) method

/// This method sets both the minimum and the values of all the independent parameters from a single matrix.
/// The first row must contain the minimum values values for the independent parameters.
/// The second row must contain the maximum values for the independent parameters.
/// These values are used for pre and postprocessing the independent parameters with the
/// minimum and maximum method.
///
/// @param newMinimumAndMaximumOfIndependentParameters New set of minimum and maximum values for the independent 
/// parameters of the neural network.
///
/// @see setMinimumOfIndependentParameters(Vector<double>).
/// @see setMaximumOfIndependentParameters(Vector<double>).
/// @see setMinimumOfSingleIndependentParameter(int, double).
/// @see setMaximumOfSingleIndependentParameter(int, double).

void MultilayerPerceptron
::setMinimumAndMaximumOfIndependentParameters(Matrix<double> newMinimumAndMaximumOfIndependentParameters)
{
   // Control sentence (if debug)

   #ifndef NDEBUG 

   if(newMinimumAndMaximumOfIndependentParameters.getNumberOfRows() != 2)
   {
      std::cerr << std::endl
                << "Flood Error: MultilayerPerceptron class." << std::endl 
                << "void setMinimumAndMaximumOfIndependentParameters(Matrix<double>) method." << std::endl
                << "Number of rows must be 2." << std::endl 
                << std::endl;

      exit(1);
   }
   else if(newMinimumAndMaximumOfIndependentParameters.getNumberOfColumns() != numberOfIndependentParameters)
   {
      std::cerr << std::endl
                << "Flood Error: MultilayerPerceptron class." << std::endl 
                << "void setMinimumAndMaximumOfIndependentParameters(Matrix<double>) method." << std::endl
                << "Number of columns must be equal to number of independent parameters." << std::endl
                << std::endl;

      exit(1);
   }

   // Check that minimum of independent parameters is not greater than their maximum

   for(int i = 0; i < numberOfIndependentParameters; i++)
   {
      if(newMinimumAndMaximumOfIndependentParameters[0][i] >= newMinimumAndMaximumOfIndependentParameters[1][i])
      {
         std::cerr << std::endl
                   << "Flood Error: MultilayerPerceptron class." << std::endl
                   << "void setMinimumAndMaximumOfIndependentParameters(Matrix<double>) method." << std::endl
                   << "Minimum of independent parameter "<< i 
				   << " is equal or greater than maximum of that parameter."
                   << std::endl << std::endl;

         exit(1);
      }
   }

   #endif

   // Set minimum and maximum of independent parameters

   minimumOfIndependentParameters = newMinimumAndMaximumOfIndependentParameters.getRow(0);
   maximumOfIndependentParameters = newMinimumAndMaximumOfIndependentParameters.getRow(1);
}


// void setMinimumOfSingleIndependentParameter(int, double) method

/// This method sets a minimum value for a single independent parameter.
/// Such a value is used for pre and postprocessing that independent parameter with the minimum and maximum 
/// method.
///
/// @param i Index of independent parameter.
/// @param newMinimumOfSingleIndependentParameter New minimum value for the independent parameter of index i.
///
/// @see setMinimumOfIndependentParameters(Vector<double>).
/// @see setMaximumOfIndependentParameters(Vector<double>).
/// @see setMinimumAndMaximumOfIndependentParameters(Matrix<double>).
/// @see setMaximumOfSingleIndependentParameter(int, double).

void MultilayerPerceptron::setMinimumOfSingleIndependentParameter(int i, 
double newMinimumOfSingleIndependentParameter)
{
   // Control sentence (if debug)

   #ifndef NDEBUG 

   if(i >= numberOfIndependentParameters)
   {
      std::cerr << std::endl
                << "Flood Error: MultilayerPerceptron class." << std::endl
                << "void setMinimumOfSingleIndependentParameter(int, double) method." << std::endl
                << "Index must be less than number of independent parameters." << std::endl 
                << std::endl;

      exit(1);   
   }

   #endif

   // Set minimum of single independent parameter

   minimumOfIndependentParameters[i] = newMinimumOfSingleIndependentParameter;
}


// void setMaximumOfSingleIndependentParameter(int, double) method

/// This method sets a maximum value for a single independent parameter.
/// Such a value is used for pre and postprocessing that independent parameter with the minimum and maximum 
/// method.
///
/// @param i Index of independent parameter.
/// @param newMaximumOfSingleIndependentParameter New maximum value for the independent parameter of index i.
///
/// @see setMinimumOfIndependentParameters(Vector<double>).
/// @see setMaximumOfIndependentParameters(Vector<double>).
/// @see setMinimumAndMaximumOfIndependentParameters(Matrix<double>).
/// @see setMinimumOfSingleIndependentParameter(int, double).

void MultilayerPerceptron::setMaximumOfSingleIndependentParameter(int i, 
double newMaximumOfSingleIndependentParameter)
{
   // Control sentence (if debug)

   #ifndef NDEBUG 

   if(i >= numberOfIndependentParameters)
   {
      std::cerr << std::endl
                << "Flood Error: MultilayerPerceptron class." << std::endl
                << "void setMaximumOfSingleIndependentParameter(int, double) method." << std::endl
                << "Index must be less than number of independent parameters." << std::endl
                << std::endl;

      exit(1);   
   }

   #endif

   // Set maximum of single independent parameter

   maximumOfIndependentParameters[i] = newMaximumOfSingleIndependentParameter;
}


// void setLowerBoundOfIndependentParameters(Vector<double>) method

/// This method sets the lower bound of all the independent parameters.
/// These values are used for postprocessing the independent parameters so that they are not less than the lower 
/// bounds. 
///
/// @param newLowerBoundOfIndependentParameters New set of lower bounds for the independent parameters of the 
/// neural network.
///
/// @see setUpperBoundOfIndependentParameters(Vector<double>).
/// @see setLowerAndUpperBoundsOfIndependentParameters(Matrix<double>).
/// @see setLowerBoundOfSingleIndependentParameters(int, double).
/// @see setUpperBoundOfSingleIndependentParameters(int, double).

void MultilayerPerceptron
::setLowerBoundOfIndependentParameters(Vector<double> newLowerBoundOfIndependentParameters)
{
   // Control sentence (if debug)

   #ifndef NDEBUG 

   if(newLowerBoundOfIndependentParameters.getSize() != numberOfIndependentParameters)
   {
      std::cerr << std::endl
                << "Flood Error: MultilayerPerceptron class." << std::endl 
                << "void setLowerBoundOfIndependentParameters(Vector<double>) method." << std::endl
                << "Size must be equal to number of independent parameters." << std::endl
                << std::endl;

      exit(1);
   }

   #endif

   // Set lower bound of independent parameters

   lowerBoundOfIndependentParameters = newLowerBoundOfIndependentParameters;                                                   
}


// void setUpperBoundOfIndependentParameters(Vector<double>) method

/// This method sets the upper bound of all the independent parameters.
/// These values are used for postprocessing the independent parameters so that they are not greater than the 
/// upper bounds. 
///
/// @param newUpperBoundOfIndependentParameters New set of upper bounds for the independent parameters of the 
/// neural network.
///
/// @see setLowerBoundOfIndependentParameters(Vector<double>).
/// @see setLowerAndUpperBoundsOfIndependentParameters(Matrix<double>).
/// @see setLowerBoundOfSingleIndependentParameters(int, double).
/// @see setUpperBoundOfSingleIndependentParameters(int, double).

void MultilayerPerceptron::setUpperBoundOfIndependentParameters(Vector<double> newUpperBoundOfIndependentParameters)
{
   // Control sentence (if debug)

   #ifndef NDEBUG 

   if(newUpperBoundOfIndependentParameters.getSize() != numberOfIndependentParameters)
   {
      std::cerr << std::endl
                << "Flood Error: MultilayerPerceptron class." << std::endl 
                << "void setUpperBoundOfIndependentParameters(Vector<double>) method." << std::endl
                << "Size must be equal to number of independent parameters." << std::endl
                << std::endl;

      exit(1);
   }

   #endif

   // Set upper bound of independent parameters

   upperBoundOfIndependentParameters = newUpperBoundOfIndependentParameters;                                                   
}


// void setLowerAndUpperBoundsOfIndependentParameters(Matrix<double>) method

/// This method sets both the lower and the upper bounds of all the independent parameters from a single matrix.
/// The first row must contain the lower bound values values for the independent parameters.
/// The second row must contain the upper bound values for the independent parameters.
/// These values are used for postprocessing the independent parameters so that they are neither less than the 
/// lower bounds nor greater than the upper bounds. 
///
/// @param newLowerAndUpperBoundsOfIndependentParameters New set of lower and upper bounds for the independent 
/// parameters of the neural network.
///
/// @see setLowerBoundOfIndependentParameters(Vector<double>).
/// @see setUpperBoundOfIndependentParameters(Matrix<double>).
/// @see setLowerBoundOfSingleIndependentParameters(int, double).
/// @see setUpperBoundOfSingleIndependentParameters(int, double).

void MultilayerPerceptron
::setLowerAndUpperBoundsOfIndependentParameters(Matrix<double> newLowerAndUpperBoundsOfIndependentParameters)
{
   // Control sentence (if debug)

   #ifndef NDEBUG 

   int numberOfRows = newLowerAndUpperBoundsOfIndependentParameters.getNumberOfRows();
   int numberOfColumns = newLowerAndUpperBoundsOfIndependentParameters.getNumberOfColumns();

   if(numberOfRows != 2)
   {
      std::cerr << std::endl
                << "Flood Error: MultilayerPerceptron class." << std::endl 
                << "void setLowerAndUpperBoundsOfIndependentParameters(Matrix<double>) method." << std::endl
                << "Number of rows must be 2." << std::endl
                << std::endl;

      exit(1);
   }
   else if(numberOfColumns != numberOfIndependentParameters)
   {
      std::cerr << std::endl
                << "Flood Error: MultilayerPerceptron class." << std::endl 
                << "void setLowerAndUpperBoundsOfIndependentParameters(Matrix<double>) method." << std::endl
                << "Number of columns must be equal to number of independent parameters." << std::endl
                << std::endl;

      exit(1);
   }

   #endif

   // Set lower and upper bounds of independent parameters

   lowerBoundOfIndependentParameters = newLowerAndUpperBoundsOfIndependentParameters.getRow(0);
   upperBoundOfIndependentParameters = newLowerAndUpperBoundsOfIndependentParameters.getRow(1);
}


// void setLowerBoundOfSingleIndependentParameter(int, double) method

/// This method sets the lower bound of a single independent parameter.
/// Such a value is used for postprocessing that independent parameter so that it is not less than its lower 
/// bound. 
///
/// @param i Index of independent parameter.
/// @param newLowerBoundOfSingleIndependentParameter New lower bound for the independent parameter of index i.
///
/// @see setLowerBoundOfIndependentParameters(Vector<double>).
/// @see setUpperBoundOfIndependentParameters(Vector<double>).
/// @see setLowerAndUpperBoundsOfIndependentParameters(Matrix<double>).
/// @see setUpperBoundOfSingleIndependentParameters(int, double).

void MultilayerPerceptron::setLowerBoundOfSingleIndependentParameter(int i, 
double newLowerBoundOfSingleIndependentParameter)
{
   // Control sentence (if debug)

   #ifndef NDEBUG 

   if(i >= numberOfIndependentParameters)
   {
      std::cerr << std::endl
                << "Flood Error: MultilayerPerceptron class." << std::endl
                << "void setLowerBoundOfSingleIndependentParameter(int, double) method." << std::endl
                << "Index must be less than number of independent parameters." << std::endl 
                << std::endl;

      exit(1);   
   }

   #endif

   // Set lower bound of single independent parameter

   lowerBoundOfIndependentParameters[i] = newLowerBoundOfSingleIndependentParameter;
}


// void setUpperBoundOfSingleIndependentParameter(int, double) method

/// This method sets the upper bound of a single independent parameter.
/// Such a value is used for postprocessing that independent parameter so that it is not greater than its upper 
/// bound. 
///
/// @param i Index of independent parameter.
/// @param newUpperBoundOfSingleIndependentParameter New upper bound for the independent parameter of index i.
///
/// @see setLowerBoundOfIndependentParameters(Vector<double>).
/// @see setUpperBoundOfIndependentParameters(Vector<double>).
/// @see setLowerAndUpperBoundsOfIndependentParameters(Matrix<double>).
/// @see setLowerBoundOfSingleIndependentParameters(int, double).

void MultilayerPerceptron::setUpperBoundOfSingleIndependentParameter(int i, 
double newUpperBoundOfSingleIndependentParameter)
{
   // Control sentence (if debug)

   #ifndef NDEBUG 

   if(i >= numberOfIndependentParameters)
   {
      std::cerr << std::endl
                << "Flood Error: MultilayerPerceptron class." << std::endl
                << "void setUpperBoundOfSingleIndependentParameter(int, double) method." << std::endl
                << "Index must be less than number of independent parameters." << std::endl
                << std::endl;

      exit(1);   
   }

   #endif

   // Set upper bound of single independent parameter

   upperBoundOfIndependentParameters[i] = newUpperBoundOfSingleIndependentParameter;
}


// void setNumericalEpsilon(double) method

/// This method sets a new value of epsilon for the computation of the Jacobian matrix by means of numerical 
/// differentiation.
///
/// @param newNumericalEpsilon New value for epsilon. 
///
/// @see calculateJacobianTest(Vector<double>).

void MultilayerPerceptron::setNumericalEpsilon(double newNumericalEpsilon)
{
   // Control sentence (if debug)

   #ifndef NDEBUG 

   if(newNumericalEpsilon <= 0.0)
   {
      std::cerr << std::endl
                << "Flood Error: MultilayerPerceptron class." << std::endl
                << "void setNumericalEpsilon(double) method." << std::endl
                << "Epsilon must be greater than zero." << std::endl
                << std::endl;

      exit(1);
   }

   #endif

   // Set numerical epsilon

   numericalEpsilon = newNumericalEpsilon;
}


// void setDisplayOutOfRangeWarning(bool) method

/// This method sets a new display out of range warning value. 
/// If it is set to true a warning message is to be displayed when inputs fall outside minimum and maximums.
/// if it is set to false that warning message is not to be displayed
///
/// @param newDisplayOutOfRangeWarning Display value for the out of range warning.

void MultilayerPerceptron::setDisplayOutOfRangeWarning(bool newDisplayOutOfRangeWarning)
{
   displayOutOfRangeWarning = newDisplayOutOfRangeWarning;
}


// void setDisplay(bool) method

/// This method sets a new display value. 
/// If it is set to true messages from this class are to be displayed on the screen;
/// if it is set to false messages from this class are not to be displayed on the screen.
///
/// @param newDisplay Display value.

void MultilayerPerceptron::setDisplay(bool newDisplay)
{
   display = newDisplay;
}


int MultilayerPerceptron::getNumberOfHiddenLayers(void)
{
   int numberOfHiddenLayers = numbersOfHiddenNeurons.getSize();

   return(numberOfHiddenLayers);
}


// int getNumberOfNeuralParameters(void) method

/// This method returns the number of biases and synaptic weights in the neural network. 
///
/// @see getNumberOfIndependentParameters(void).
/// @see getNumberOfFreeParameters(void).

int MultilayerPerceptron::getNumberOfNeuralParameters(void)
{
   int numberOfHiddenLayers = getNumberOfHiddenLayers();

   int numberOfNeuralParameters = 0;

   // Hidden layers parameters

   for(int i = 0; i < numberOfHiddenLayers; i++)
   {
      for(int j = 0; j < numbersOfHiddenNeurons[i]; j++)
      {         numberOfNeuralParameters += hiddenLayers[i][j].getNumberOfNeuronParameters();
      }
   }

   // Output layer parameters

   for(int i = 0; i < numberOfOutputs; i++)
   {      numberOfNeuralParameters += outputLayer[i].getNumberOfNeuronParameters();
   }

   return(numberOfNeuralParameters);
}


// int getNumberOfFreeParameters(void) method

/// This method returns the number of free parameters in the neural network. 
/// The number of free parameters is the sum of all the biases, synaptic weights and independent parameters.
///
/// @see getNumberOfNeuralParameters(void).
/// @see getNumberOfIndependentParameters(void).

int MultilayerPerceptron::getNumberOfFreeParameters(void)
{
   int numberOfFreeParameters = 0;
  
   int numberOfNeuralParameters = getNumberOfNeuralParameters();

   numberOfFreeParameters = numberOfNeuralParameters + numberOfIndependentParameters;

   return(numberOfFreeParameters);
}


// Vector<double> getNeuralParameters(void) method

/// This method returns the values of all the biases and synaptic weights in the neural network as a single 
/// vector.
///
/// @see getIndependentParameters(void).
/// @see getFreeParameters(void).

Vector<double> MultilayerPerceptron::getNeuralParameters(void)
{
   int numberOfHiddenLayers = getNumberOfHiddenLayers();

   int numberOfNeuralParameters = getNumberOfNeuralParameters();
   int numberOfHiddenNeuronParameters;
   int numberOfOutputNeuronParameters;

   Vector<double> neuralParameters(numberOfNeuralParameters, 0.0);
   Vector<double> hiddenNeuronParameters;
   Vector<double> outputNeuronParameters;

   // Hidden layers parameters

   int position = 0;

   for(int i = 0; i < numberOfHiddenLayers; i++)
   {
      for(int j = 0; j < numbersOfHiddenNeurons[i]; j++)
      {
         hiddenNeuronParameters = hiddenLayers[i][j].getNeuronParameters();
         numberOfHiddenNeuronParameters = hiddenLayers[i][j].getNumberOfNeuronParameters();
         neuralParameters.insert(position, hiddenNeuronParameters);
         position += numberOfHiddenNeuronParameters; 
      }
   }

   // Output layer parameters

   for(int i = 0; i < numberOfOutputs; i++)
   {
      outputNeuronParameters = outputLayer[i].getNeuronParameters();
      numberOfOutputNeuronParameters = outputLayer[i].getNumberOfNeuronParameters();
      neuralParameters.insert(position, outputNeuronParameters);       
      position += numberOfOutputNeuronParameters; 
   }

   return(neuralParameters);
}


// Vector<double> getFreeParameters(void) method

/// This method returns the values of the free parameters in the neural network as a single vector.
/// This contains all the biases, synaptic weights and scaled independent parameters. 

Vector<double> MultilayerPerceptron::getFreeParameters(void)
{
   Vector<double> neuralParameters = getNeuralParameters();

   Vector<double> independentParameters = getIndependentParameters();

   Vector<double> scaledIndependentParameters = scaleIndependentParameters(independentParameters);

   Vector<double> freeParameters = neuralParameters.assemble(scaledIndependentParameters);

   return(freeParameters);
}


// void setNeuralParameters(Vector<double>) method

/// This method sets all the biases and synaptic weights in the neural network from a single vector.
///
/// @param newNeuralParameters New set of biases and synaptic weights values. 

void MultilayerPerceptron::setNeuralParameters(Vector<double> newNeuralParameters)
{
   // Control sentence (if debug)

   #ifndef NDEBUG 

   int size = newNeuralParameters.getSize();

   int numberOfNeuralParameters = getNumberOfNeuralParameters();

   if(size != numberOfNeuralParameters)
   {
      std::cerr << std::endl
                << "Flood Error: MultilayerPerceptron class." << std::endl
                << "void setNeuralParameters(Vector<double>) method." << std::endl
                << "Size must be equal to number of biases and synaptic weights." << std::endl
                << std::endl;

      exit(1);
   }

   #endif

   int numberOfHiddenLayers = getNumberOfHiddenLayers();

   // Hidden layers parameters

   int index = 0;

   for(int i = 0; i < numberOfHiddenLayers; i++)
   {
      for(int j = 0; j < numbersOfHiddenNeurons[i]; j++)
      {         hiddenLayers[i][j]
         .setNeuronParameters(newNeuralParameters.extract(index, hiddenLayers[i][j].getNumberOfNeuronParameters()));

         index += hiddenLayers[i][j].getNumberOfNeuronParameters();
      }
   }

   // Output layer parameters

   for(int i = 0; i < numberOfOutputs; i++)
   {      outputLayer[i]
     .setNeuronParameters(newNeuralParameters.extract(index, outputLayer[i].getNumberOfNeuronParameters()));

      index += outputLayer[i].getNumberOfNeuronParameters();
   }
}


// void setFreeParameters(Vector<double>) method

/// This method sets all the free parameters in the neural network from a single vector.
///
/// @param newFreeParameters New set of free parameter values. 

void MultilayerPerceptron::setFreeParameters(Vector<double> newFreeParameters)
{
   // Control sentence (if debug)

   #ifndef NDEBUG 

   int size = newFreeParameters.getSize();

   int numberOfFreeParameters = getNumberOfFreeParameters();

   if(size != numberOfFreeParameters)
   {
      std::cerr << std::endl
                << "Flood Error: MultilayerPerceptron class." << std::endl
                << "void setFreeParameters(Vector<double>) method." << std::endl
                << "Size must be equal to number of free parameters." << std::endl
                << std::endl;

      exit(1);
   }

   #endif

   // Biases and synaptic weights

   int numberOfNeuralParameters = getNumberOfNeuralParameters();

   Vector<double> newNeuralParameters(numberOfNeuralParameters);

   for(int i = 0; i < numberOfNeuralParameters; i++)
   {
      newNeuralParameters[i] = newFreeParameters[i];
   }

   setNeuralParameters(newNeuralParameters);

   // Scaled independent parameters
   
   Vector<double> scaledIndependentParameters(numberOfIndependentParameters);

   for(int i = 0; i < numberOfIndependentParameters; i++)
   {
      scaledIndependentParameters[i] = newFreeParameters[numberOfNeuralParameters+i];
   }

   Vector<double> unscaledIndependentParameters = unscaleIndependentParameters(scaledIndependentParameters);

   setIndependentParameters(unscaledIndependentParameters);
}


// void initNeuralParametersUniform(void) method

/// This method initializes all the biases and synaptic weights in the newtork at random with values comprised 
/// between -1 and +1.

void MultilayerPerceptron::initNeuralParametersUniform(void)
{
   int numberOfNeuralParameters = getNumberOfNeuralParameters();

   Vector<double> neuralParameters(numberOfNeuralParameters);

   double random;

   for(int i = 0; i < numberOfNeuralParameters; i++)
   {
      random = (double)rand()/(RAND_MAX+1.0);

      neuralParameters[i] = -1.0 + 2.0*random;
   }

   setNeuralParameters(neuralParameters);
}


// void initNeuralParametersUniform(double, double) method

/// This method initializes all the biases and synaptic weights in the neural network at random with values 
/// comprised between a minimum and a maximum values.
///
/// @param minimum Minimum initialization value.
/// @param maximum Maximum initialization value.

void MultilayerPerceptron::initNeuralParametersUniform(double minimum, double maximum)
{
   // Control sentence (if debug)

   #ifndef NDEBUG 

   if(minimum > maximum)
   {
      std::cerr << std::endl
                << "Flood Error: MultilayerPerceptron class." << std::endl 
                << "void initNeuralParametersUniform(double, double) method." << std::endl
                << "Minimum value must be less or equal than maximum value." << std::endl
                << std::endl;

      exit(1);
   }

   #endif

   int numberOfNeuralParameters = getNumberOfNeuralParameters();

   Vector<double> neuralParameters(numberOfNeuralParameters);

   double random;

   for(int i = 0; i < numberOfNeuralParameters; i++)
   {
      random = (double)rand()/(RAND_MAX+1.0);

      neuralParameters[i] = minimum + (maximum - minimum)*random;
   }

   setNeuralParameters(neuralParameters);
}


// void initNeuralParametersUniform(Vector<double>, Vector<double>) method

/// This method initializes all the biases and synaptic weights in the neural network at random, with values 
/// comprised between different minimum and maximum numbers for each parameter.
///
/// @param minimum Vector of minimum initialization values.
/// @param maximum Vector of maximum initialization values.

void MultilayerPerceptron::initNeuralParametersUniform(Vector<double> minimum, Vector<double> maximum)
{
   int numberOfNeuralParameters = getNumberOfNeuralParameters();

   // Control sentence (if debug)

   #ifndef NDEBUG 

   int minimumSize = minimum.getSize();
   int maximumSize = maximum.getSize();

   if(minimumSize != numberOfNeuralParameters || maximumSize != numberOfNeuralParameters)
   {
      std::cerr << std::endl
                << "Flood Error: MultilayerPerceptron class." << std::endl
                << "void initNeuralParametersUniform(Vector<double>, Vector<double>) method." << std::endl
                << "Sizes must be equal to number of biases and synaptic weights." << std::endl
                << std::endl;

      exit(1);
   }

   for(int i = 0; i < numberOfNeuralParameters; i++)
   {
      if(minimum[i] > maximum[i])
      {
         std::cerr << std::endl
                   << "Flood Error: MultilayerPerceptron class." << std::endl
                   << "void initNeuralParametersUniform(Vector<double>, Vector<double>) method." 
                   << std::endl
                   << "Minimum value must be less or equal than maximum value." << std::endl
                   << std::endl;

         exit(1);
      }
   }

   #endif

   Vector<double> neuralParameters(numberOfNeuralParameters);

   double random;

   for(int i = 0; i < numberOfNeuralParameters; i++)
   {
      random = (double)rand()/(RAND_MAX+1.0);

      neuralParameters[i] = minimum[i] + (maximum[i] - minimum[i])*random;
   }

   setNeuralParameters(numberOfNeuralParameters);
}


// void initNeuralParametersUniform(Matrix<double>) method

/// This method initializes all the biases and synaptic weights in the neural network at random, with values 
/// comprised between a different minimum and maximum numbers for each parameter.
/// All minimum are maximum initialization values must be given from a single matrix.
/// The first row must contain the minimum inizizalization value for each parameter.
/// The second row must contain the maximum inizizalization value for each parameter.
///
/// @param minimumAndMaximum Matrix of minimum and maximum initialization values.

void MultilayerPerceptron::initNeuralParametersUniform(Matrix<double> minimumAndMaximum)
{
   int numberOfNeuralParameters = getNumberOfNeuralParameters();

   // Control sentence (if debug)

   #ifndef NDEBUG 

   int numberOfRows = minimumAndMaximum.getNumberOfRows();
   int numberOfColumns = minimumAndMaximum.getNumberOfColumns(); 

   if(numberOfRows != 2 || numberOfColumns != numberOfNeuralParameters)
   {
      std::cerr << std::endl
                << "Flood Error: MultilayerPerceptron class." << std::endl 
                << "void initNeuralParametersUniform(Matrix<double>) method." << std::endl
                << "Numbers of rows and columns must be two and number of biases and synaptic weights, "
                << "respectively." << std::endl 
                << std::endl;

      exit(1);
   }

   #endif

   Vector<double> minimum = minimumAndMaximum.getRow(0);
   Vector<double> maximum = minimumAndMaximum.getRow(1) ;

   #ifndef NDEBUG 

   for(int i = 0; i < numberOfNeuralParameters; i++)
   {
      if(minimum[i] > maximum[i])
      {
         std::cerr << std::endl
                   << "Flood Error: MultilayerPerceptron class." << std::endl
                   << "void initNeuralParametersUniform(Vector<double>, Vector<double>) method." 
                   << std::endl
                   << "Minimum value must be less or equal than maximum value." << std::endl
                   << std::endl;

         exit(1);
      }
   }

   #endif

   Vector<double> neuralParameters(numberOfNeuralParameters);

   double random;

   for(int i = 0; i < numberOfNeuralParameters; i++)
   {
      random = (double)rand()/(RAND_MAX+1.0);

      neuralParameters[i] = minimum[i] + (maximum[i] - minimum[i])*random;
   }

   setNeuralParameters(neuralParameters);
}


// void initNeuralParametersNormal(void) method

/// This method initializes all the biases and synaptic weights in the newtork with random values chosen from a 
/// uniform distribution with mean 0 and standard deviation 1.

void MultilayerPerceptron::initNeuralParametersNormal(void)
{
   int numberOfNeuralParameters = getNumberOfNeuralParameters();

   Vector<double> neuralParameters(numberOfNeuralParameters);

   double random;

   for(int i = 0; i < numberOfNeuralParameters; i++)
   {
      random = (double)rand()/(RAND_MAX+1.0);

      neuralParameters[i] = -1.0 + 2.0*random;
   }

   setNeuralParameters(neuralParameters);
}


// void initNeuralParametersNormal(double, double) method

/// This method initializes all the biases and synaptic weights in the neural network with random random values 
/// chosen from a normal distribution with a given mean and a given standard deviation.
///
/// @param mean Mean of normal distribution.
/// @param standardDeviation Standard deviation of normal distribution.

void MultilayerPerceptron
::initNeuralParametersNormal(double mean, double standardDeviation)
{
   // Control sentence (if debug)

   #ifndef NDEBUG 

   if(standardDeviation < 0.0)
   {
      std::cerr << std::endl
                << "Flood Error: MultilayerPerceptron class."  << std::endl
                << "void initNeuralParametersNormal(double, double) method." << std::endl
                << "Standard deviation must be equal or greater than zero." << std::endl
                << std::endl;

      exit(1);
   }

   #endif

   double pi = 4.0*atan(1.0);
   double random1 = 0.0;
   double random2 = 0.0;

   int numberOfNeuralParameters = getNumberOfNeuralParameters();

   Vector<double> neuralParameters(numberOfNeuralParameters);

   for(int i = 0; i < numberOfNeuralParameters; i++)
   {
      do 
      {
         random1 = (double)rand()/(RAND_MAX+1.0);

      }while(random1 == 0.0);
      
      random2 = (double)rand()/(RAND_MAX+1.0);

      // Box-Muller transformation

      neuralParameters[i]
      = mean + sqrt(-2.0*log(random1))*sin(2.0*pi*random2)*standardDeviation;
   }

   setNeuralParameters(neuralParameters);
}


// void initNeuralParametersNormal(Vector<double>, Vector<double>) method

/// This method initializes all the biases an synaptic weights in the neural network with random values chosen 
/// from normal distributions with different mean and standard deviation for each parameter.
///
/// @param mean Vector of mean values.
/// @param standardDeviation Vector of standard deviation values.

void MultilayerPerceptron::initNeuralParametersNormal(Vector<double> mean, 
Vector<double> standardDeviation)
{
   int numberOfNeuralParameters = getNumberOfNeuralParameters();

   // Control sentence (if debug)

   #ifndef NDEBUG 

   int meanSize = mean.getSize();
   int standardDeviationSize = standardDeviation.getSize();

   if(meanSize != numberOfNeuralParameters || standardDeviationSize != numberOfNeuralParameters)
   {
      std::cerr << std::endl
                << "Flood Error: MultilayerPerceptron class." << std::endl
                << "void initNeuralParametersNormal(Vector<double>, Vector<double>) method."  << std::endl
                << "Mean and standard deviation sizes must be equal to number of biases and synaptic weights." 
                << std::endl
                << std::endl;

      exit(1);
   }

   for(int i = 0; i < numberOfNeuralParameters; i++)
   {
      if(standardDeviation[i] < 0.0)
      {
         std::cerr << std::endl
                   << "Flood Error: MultilayerPerceptron class." << std::endl
                   << "void initNeuralParametersNormal(Vector<double>, Vector<double>) method." 
                   << std::endl
                   << "Standard deviations must equal or greater than zero." << std::endl
                   << std::endl;

         exit(1);
      }
   }

   #endif

   Vector<double> neuralParameters(numberOfNeuralParameters);

   double pi = 4.0*atan(1.0);

   double random1 = 0.0;
   double random2 = 0.0;
    
   for(int i = 0; i < numberOfNeuralParameters; i++)
   {
     do
     {
      random1 = (double)rand()/(RAND_MAX+1.0);

     }while(random1 == 0.0); 

      random2 = (double)rand()/(RAND_MAX+1.0);

      // Box-Muller transformation

      neuralParameters[i]
      = mean[i] + sqrt(-2.0*log(random1))*sin(2.0*pi*random2)*standardDeviation[i];
   }

   setNeuralParameters(neuralParameters);
}


// void initNeuralParametersNormal(Matrix<double>) method

/// This method initializes all the biases and synaptic weights in the neural network with random values chosen 
/// from normal distributions with different mean and standard deviation for each parameter.
/// All mean and standard deviation values are given from a single matrix.
/// The first row must contain the mean value for each parameter.
/// The second row must contain the standard deviation value for each parameter.
///
/// @param meanAndStandardDeviation Matrix of mean and standard deviation values.

void MultilayerPerceptron::initNeuralParametersNormal(Matrix<double> meanAndStandardDeviation)
{
   int numberOfNeuralParameters = getNumberOfNeuralParameters();

   // Control sentence (if debug)

   #ifndef NDEBUG 

   int numberOfRows = meanAndStandardDeviation.getNumberOfRows();
   int numberOfColumns = meanAndStandardDeviation.getNumberOfColumns();

   if(numberOfRows != 2 || numberOfColumns != numberOfNeuralParameters)
   {
      std::cerr << std::endl
                << "Flood Error: MultilayerPerceptron class." << std::endl 
                << "void initNeuralParametersNormal(Matrix<double>) method." << std::endl
                << "Numbers of rows and columns must be two and number of biases and synaptic weights, "
				<< "respectively." << std::endl 
                << std::endl;

      exit(1);
   }

   #endif

   Vector<double> mean = meanAndStandardDeviation.getRow(0);
   Vector<double> standardDeviation = meanAndStandardDeviation.getRow(1);

   #ifndef NDEBUG 

   for(int i = 0; i < numberOfNeuralParameters; i++)
   {
      if(standardDeviation[i] < 0.0)
      {
         std::cerr << std::endl
                   << "Flood Error: MultilayerPerceptron class." << std::endl
                   << "void initNeuralParametersNormal(Vector<double>, Vector<double>) method." 
                   << std::endl
                   << "Standard deviations must be equal or greater than zero." << std::endl
                   << std::endl;

         exit(1);
      }
   }

   #endif

   Vector<double> neuralParameters(numberOfNeuralParameters);

   double pi = 4.0*atan(1.0);

   double random1 = 0.0;
   double random2 = 0.0;

   for(int i = 0; i < numberOfNeuralParameters; i++)
   {
      do
      {
         random1 = (double)rand()/(RAND_MAX+1.0);

      }while(random1 == 0.0);

      random2 = (double)rand()/(RAND_MAX+1.0);

      // Box-Muller transformation

      neuralParameters[i] = mean[i] + sqrt(-2.0*log(random1))*sin(2.0*pi*random2)*standardDeviation[i];
   }

   setNeuralParameters(neuralParameters);
}


// void initIndependentParametersUniform(void) method

/// This method initializes the independent parameters associated to the newtork at random with values comprised 
/// between -1 and +1.

void MultilayerPerceptron::initIndependentParametersUniform(void)
{
   double random = 0.0;

   for(int i = 0; i < numberOfIndependentParameters; i++)
   {
      random = (double)rand()/(RAND_MAX+1.0);

      independentParameters[i] = -1.0 + 2.0*random;
   }
}


// void initIndependentParametersUniform(double, double) method

/// This method initializes the independent parameters associated to the neural network at random with values 
/// comprised between a minimum and a maximum values.
///
/// @param minimum Minimum initialization value.
/// @param maximum Maximum initialization value.

void MultilayerPerceptron::initIndependentParametersUniform(double minimum, double maximum)
{
   // Control sentence (if debug)

   #ifndef NDEBUG 

   if(minimum > maximum)
   {
      std::cerr << std::endl
                << "Flood Error: MultilayerPerceptron class." << std::endl 
                << "void initIndependentParametersUniform(double, double) method." << std::endl
                << "Minimum value must be less or equal than maximum value." << std::endl
                << std::endl;

      exit(1);
   }

   #endif

   double random = 0.0;

   for(int i = 0; i < numberOfIndependentParameters; i++)
   {
      random = (double)rand()/(RAND_MAX+1.0);

      independentParameters[i] = minimum + (maximum - minimum)*random;
   }
}


// void initIndependentParametersUniform(Vector<double>, Vector<double>) method

/// This method initializes the independent parameters associated to the neural network at random with values 
/// comprised between different minimum and maximum numbers for each independent parameter.
///
/// @param minimum Vector of minimum initialization values.
/// @param maximum Vector of maximum initialization values.

void MultilayerPerceptron::initIndependentParametersUniform(Vector<double> minimum, Vector<double> maximum)
{
   // Control sentence (if debug)

   #ifndef NDEBUG 

   int minimumSize = minimum.getSize();
   int maximumSize = maximum.getSize();

   if(minimumSize != numberOfIndependentParameters || maximumSize != numberOfIndependentParameters)
   {
      std::cerr << std::endl
                << "Flood Error: MultilayerPerceptron class." << std::endl
                << "void initIndependentParametersUniform(Vector<double>, Vector<double>) method." << std::endl
                << "Sizes must be equal to number of independent parameters." << std::endl
                << std::endl;

      exit(1);
   }

   for(int i = 0; i < numberOfIndependentParameters; i++)
   {
      if(minimum[i] > maximum[i])
      {
         std::cerr << std::endl
                   << "Flood Error: MultilayerPerceptron class." << std::endl
                   << "void initIndependentParametersUniform(Vector<double>, Vector<double>) method." 
                   << std::endl
                   << "Minimum value must be less or equal than maximum value." << std::endl
                   << std::endl;

         exit(1);
      }
   }

   #endif

   double random = 0.0;

   for(int i = 0; i < numberOfIndependentParameters; i++)
   {
      random = (double)rand()/(RAND_MAX+1.0);

      independentParameters[i] 
     = minimum[i] + (maximum[i] - minimum[i])*random;
   }
}


// void initIndependentParametersUniform(Matrix<double>) method

/// This method initializes the independent parameters associated to the neural network at random values comprised
/// between different minimum and maximum numbers for each independent parameter.
/// All minimum and maximum values are given from a single matrix.
/// The first row must contain the minimum inizizalization value for each independent parameter.
/// The second row must contain the maximum inizizalization value for each independent parameter.
///
/// @param minimumAndMaximum Matrix of minimum and maximum initialization values.

void MultilayerPerceptron::initIndependentParametersUniform(Matrix<double> minimumAndMaximum)
{
   // Control sentence (if debug)

   #ifndef NDEBUG 

   int numberOfRows = minimumAndMaximum.getNumberOfRows();
   int numberOfColumns = minimumAndMaximum.getNumberOfColumns();
   
   if(numberOfRows != 2 || numberOfColumns != numberOfIndependentParameters)
   {
      std::cerr << std::endl
                << "Flood Error: MultilayerPerceptron class." << std::endl 
                << "void initIndependentParametersUniform(Matrix<double>) method." 
                << std::endl
                << "Numbers of rows and columns must be two and number of independent parameters, respectively" 
                << std::endl << std::endl;

      exit(1);
   }

   #endif

   Vector<double> minimum = minimumAndMaximum.getRow(0);
   Vector<double> maximum = minimumAndMaximum.getRow(1);

   #ifndef NDEBUG 

   for(int i = 0; i < numberOfIndependentParameters; i++)
   {
      if(minimum[i] > maximum[i])
      {
         std::cerr << std::endl
                   << "Flood Error: MultilayerPerceptron class." << std::endl
                   << "void initIndependentParametersUniform(Vector<double>, Vector<double>) method." 
                   << std::endl
                   << "Minimum value must be less or equal than maximum value." << std::endl
                   << std::endl;

         exit(1);
      }
   }

   #endif

   // Init independent parameters

   double random = 0.0;

   for(int i = 0; i < numberOfIndependentParameters; i++)
   {
      random = (double)rand()/(RAND_MAX+1.0);

      independentParameters[i] = minimum[i] + (maximum[i] - minimum[i])*random;
   }
}


// void initIndependentParametersNormal(void) method

/// This method initializes the independent parameters associated to the newtork with random values chosen from a 
/// normal distribution with mean 0 and standard deviation 1.

void MultilayerPerceptron::initIndependentParametersNormal(void)
{
   double random = 0.0;
  
   for(int i = 0; i < numberOfIndependentParameters; i++)
   {
      random = (double)rand()/(RAND_MAX+1.0);

      independentParameters[i] = -1.0 + 2.0*random;
   }
}


// void initIndependentParametersNormal(double, double) method

/// This method initializes the independent parameters associated to the neural network with random values chosen 
/// from a normal distribution with a given mean and a given standard deviation.
///
/// @param mean Mean of normal distribution.
/// @param standardDeviation Standard deviation of normal distribution.

void MultilayerPerceptron::initIndependentParametersNormal(double mean, double standardDeviation)
{
   // Control sentence (if debug)

   #ifndef NDEBUG 

   if(standardDeviation < 0.0)
   {
      std::cerr << std::endl
                << "Flood Error: MultilayerPerceptron class."  << std::endl
                << "void initIndependentParametersNormal(double, double) method." << std::endl
                << "Standard deviation must be equal or greater than zero." << std::endl
                << std::endl;

      exit(1);
   }

   #endif

   // Init independent parameters
   
   double pi = 4.0*atan(1.0);

   double random1 = 0.0;
   double random2 = 0.0;

   for(int i = 0; i < numberOfIndependentParameters; i++)
   {
      do
      {
         random1 = (double)rand()/(RAND_MAX+1.0);

      }while(random1 == 0.0);

      random2 = (double)rand()/(RAND_MAX+1.0);

      // Box-Muller transformation

      independentParameters[i] = mean + sqrt(-2.0*log(random1))*sin(2.0*pi*random2)*standardDeviation;
   }
}


// void initIndependentParametersNormal(Vector<double>, Vector<double>) method

/// This method initializes the independent parameters associated to the neural network with random values chosen 
/// from normal distributions with different mean and standard deviation for each independent parameter.
///
/// @param mean Vector of mean values.
/// @param standardDeviation Vector of standard deviation values.

void MultilayerPerceptron::initIndependentParametersNormal(Vector<double> mean, Vector<double> standardDeviation)
{
   // Control sentence (if debug)

   #ifndef NDEBUG 

   int meanSize = mean.getSize();
   int standardDeviationSize = standardDeviation.getSize();

   if(meanSize != numberOfIndependentParameters || standardDeviationSize != numberOfIndependentParameters)
   {
      std::cerr << std::endl
                << "Flood Error: MultilayerPerceptron class." << std::endl
                << "void initIndependentParametersNormal(Vector<double>, Vector<double>) method." << std::endl
                << "Mean and standard deviation sizes must be equal to number of independent parameters." 
                << std::endl
                << std::endl;

      exit(1);
   }

   for(int i = 0; i < numberOfIndependentParameters; i++)
   {
      if(standardDeviation[i] < 0.0)
      {
         std::cerr << std::endl
                   << "Flood Error: MultilayerPerceptron class." << std::endl
                   << "void initIndependentParametersNormal(Vector<double>, Vector<double>) method." << std::endl
                   << "Standard deviations must equal or greater than zero." << std::endl
                   << std::endl;

         exit(1);
      }
   }

   #endif

   double pi = 4.0*atan(1.0);
   
   double random1 = 0.0;
   double random2 = 0.0;

   for(int i = 0; i < numberOfIndependentParameters; i++)
   {
      do
      {
         random1 = (double)rand()/(RAND_MAX+1.0);
     
      }while(random1 == 0.0);

      random2 = (double)rand()/(RAND_MAX+1.0);

      // Box-Muller transformation

      independentParameters[i] = mean[i] + sqrt(-2.0*log(random1))*sin(2.0*pi*random2)*standardDeviation[i];
   }
}


// void initIndependentParametersNormal(Matrix<double>) method

/// This method initializes the independent parameters associated to the neural network with random values chosen 
/// from normal distributions with different mean and standard deviation for each independent parameter.
/// All mean and standard deviation values are given from a single matrix
/// The first row must contain the mean value for each independent parameter.
/// The second row must contain the standard deviation value for each independent parameter.
///
/// @param meanAndStandardDeviation Matrix of mean and standard deviation values.

void MultilayerPerceptron::initIndependentParametersNormal(Matrix<double> meanAndStandardDeviation)
{
   // Control sentence (if debug)

   #ifndef NDEBUG 

   int numberOfRows = meanAndStandardDeviation.getNumberOfRows();
   int numberOfColumns = meanAndStandardDeviation.getNumberOfColumns();

   if(numberOfRows != 2 || numberOfColumns != numberOfIndependentParameters)
   {
      std::cerr << std::endl
                << "Flood Error: MultilayerPerceptron class." << std::endl 
                << "void initIndependentParametersNormal(Matrix<double>) method." << std::endl
                << "Numbers of rows and columns must be two and number of independent parameters, respectively" 
                << std::endl << std::endl;

      exit(1);
   }

   #endif

   Vector<double> mean = meanAndStandardDeviation.getRow(0);
   Vector<double> standardDeviation = meanAndStandardDeviation.getRow(1);

   #ifndef NDEBUG 

   for(int i = 0; i < numberOfIndependentParameters; i++)
   {
      if(standardDeviation[i] < 0.0)
      {
         std::cerr << std::endl
                   << "Flood Error: MultilayerPerceptron class." << std::endl
                   << "void initIndependentParametersNormal(Vector<double>, Vector<double>) method." << std::endl
                   << "Standard deviations must be equal or greater than zero." << std::endl
                   << std::endl;

         exit(1);
      }
   }

   #endif

   double pi = 4.0*atan(1.0);

   double random1 = 0.0;
   double random2 = 0.0;

   for(int i = 0; i < numberOfIndependentParameters; i++)
   {
      do
      {
          random1 = (double)rand()/(RAND_MAX+1.0);

      }while(random1 == 0.0);

      random2 = (double)rand()/(RAND_MAX+1.0);

      // Box-Muller transformation

      independentParameters[i] = mean[i] + sqrt(-2.0*log(random1))*sin(2.0*pi*random2)*standardDeviation[i];
   }
}


// void initFreeParametersUniform(void) method

/// This method initializes all the free parameters in the newtork (biases and synaptic weiths + independent 
/// parameters) at random with values comprised between -1 and +1.

void MultilayerPerceptron::initFreeParametersUniform(void)
{
   int numberOfFreeParameters = getNumberOfFreeParameters();

   Vector<double> freeParameters(numberOfFreeParameters);

   double random;

   for(int i = 0; i < numberOfFreeParameters; i++)
   {
      random = (double)rand()/(RAND_MAX+1.0);

      freeParameters[i] = -1.0 + 2.0*random;
   }

   setFreeParameters(freeParameters);
}


// void initFreeParametersUniform(double, double) method

/// This method initializes all the free parameters in the newtork (biases and synaptic weiths + independent 
/// parameters) at random with values comprised between a given minimum and a given maximum values.
///
/// @param minimum Minimum initialization value.
/// @param maximum Maximum initialization value.

void MultilayerPerceptron::initFreeParametersUniform(double minimum, double maximum)
{
   // Control sentence (if debug)

   #ifndef NDEBUG 

   if(minimum > maximum)
   {
      std::cerr << std::endl
                << "Flood Error: MultilayerPerceptron class." << std::endl 
                << "void initFreeParametersUniform(double, double) method." << std::endl
                << "Minimum value must be less or equal than maximum value." << std::endl
                << std::endl;

      exit(1);
   }

   #endif

   int numberOfFreeParameters = getNumberOfFreeParameters();

   Vector<double> freeParameters(numberOfFreeParameters);

   double random;

   for(int i = 0; i < numberOfFreeParameters; i++)
   {
      random = (double)rand()/(RAND_MAX+1.0);

      freeParameters[i] = minimum + (maximum - minimum)*random;
   }

   setFreeParameters(freeParameters);
}


// void initFreeParametersUniform(Vector<double>, Vector<double>) method

/// This method initializes all the free parameters in the newtork (biases and synaptic weiths + independent 
/// parameters) at random with values comprised between a different minimum and maximum numbers for each free 
/// parameter.
///
/// @param minimum Vector of minimum initialization values.
/// @param maximum Vector of maximum initialization values.

void MultilayerPerceptron::initFreeParametersUniform(Vector<double> minimum, Vector<double> maximum)
{
   int numberOfFreeParameters = getNumberOfFreeParameters();

   // Control sentence (if debug)

   #ifndef NDEBUG 

   int minimumSize = minimum.getSize();
   int maximumSize = maximum.getSize();

   if(minimumSize != numberOfFreeParameters || maximumSize != numberOfFreeParameters)
   {
      std::cerr << std::endl
                << "Flood Error: MultilayerPerceptron class." << std::endl
                << "void initFreeParametersUniform(Vector<double>, Vector<double>) method." << std::endl
                << "Sizes must be equal to number of free parameters." << std::endl
                << std::endl;

      exit(1);
   }

   for(int i = 0; i < numberOfFreeParameters; i++)
   {
      if(minimum[i] > maximum[i])
      {
         std::cerr << std::endl
                   << "Flood Error: MultilayerPerceptron class." << std::endl
                   << "void initFreeParametersUniform(Vector<double>, Vector<double>) method." << std::endl
                   << "Minimum value must be less or equal than maximum value." << std::endl
                   << std::endl;

         exit(1);
      }
   }

   #endif

   Vector<double> freeParameters(numberOfFreeParameters);
  
   double random;

   for(int i = 0; i < numberOfFreeParameters; i++)
   {
      random = (double)rand()/(RAND_MAX+1.0);

      freeParameters[i] = minimum[i] + (maximum[i] - minimum[i])*random;
   }

   setFreeParameters(freeParameters);
}


// void initFreeParametersUniform(Matrix<double>) method

/// This method initializes all the free parameters in the newtork (biases and synaptic weiths + independent 
/// parameters) values comprised between a different minimum and maximum numbers for each free parameter.
/// Minimum and maximum initialization values are given from a single matrix.
/// The first row must contain the minimum inizizalization value for each free parameter.
/// The second row must contain the maximum inizizalization value for each free parameter.
///
/// @param minimumAndMaximum Matrix of minimum and maximum initialization values.

void MultilayerPerceptron::initFreeParametersUniform(Matrix<double> minimumAndMaximum)
{
   int numberOfFreeParameters = getNumberOfFreeParameters();

   // Control sentence (if debug)

   #ifndef NDEBUG 

   int numberOfRows = minimumAndMaximum.getNumberOfRows();
   int numberOfColumns = minimumAndMaximum.getNumberOfColumns();

   if(numberOfRows != 2 || numberOfColumns != numberOfFreeParameters)
   {
      std::cerr << std::endl
                << "Flood Error: MultilayerPerceptron class." << std::endl 
                << "void initFreeParametersUniform(Matrix<double>) method." << std::endl
                << "Numbers of rows and columns must be two and number of free parameters, respectively." 
                << std::endl << std::endl;

      exit(1);
   }

   #endif

   Vector<double> minimum = minimumAndMaximum.getRow(0);
   Vector<double> maximum = minimumAndMaximum.getRow(1);

   #ifndef NDEBUG 

   for(int i = 0; i < numberOfFreeParameters; i++)
   {
      if(minimum[i] > maximum[i])
      {
         std::cerr << std::endl
                   << "Flood Error: MultilayerPerceptron class." << std::endl
                   << "void initFreeParametersUniform(Vector<double>, Vector<double>) method." << std::endl
                   << "Minimum value must be less or equal than maximum value." << std::endl
                   << std::endl;

         exit(1);
      }
   }

   #endif

   Vector<double> freeParameters(numberOfFreeParameters);

   double random;

   for(int i = 0; i < numberOfFreeParameters; i++)
   {
      random = (double)rand()/(RAND_MAX+1.0);

      freeParameters[i] = minimum[i] + (maximum[i] - minimum[i])*random;
   }

   setFreeParameters(freeParameters);
}


// void initFreeParametersNormal(void) method

/// This method initializes all the free parameters in the newtork (biases and synaptic weiths + independent 
/// parameters) at random with values chosen from a normal distribution with mean 0 and standard deviation 1.

void MultilayerPerceptron::initFreeParametersNormal(void)
{
   int numberOfFreeParameters = getNumberOfFreeParameters();

   Vector<double> freeParameters(numberOfFreeParameters);

   double random;

   for(int i = 0; i < numberOfFreeParameters; i++)
   {
      random = (double)rand()/(RAND_MAX+1.0);

      freeParameters[i] = -1.0 + 2.0*random;
   }

   setFreeParameters(freeParameters);
}


// void initFreeParametersNormal(double, double) method

/// This method initializes all the free parameters in the newtork (biases and synaptic weiths + independent 
/// parameters) at random with values chosen from a normal distribution with a given mean and a given standard 
/// deviation.
///
/// @param mean Mean of normal distribution.
/// @param standardDeviation Standard deviation of normal distribution.

void MultilayerPerceptron::initFreeParametersNormal(double mean, double standardDeviation)
{
   // Control sentence (if debug)

   #ifndef NDEBUG 

   if(standardDeviation < 0.0)
   {
      std::cerr << std::endl
                << "Flood Error: MultilayerPerceptron class." << std::endl
                << "void initFreeParametersNormal(double, double) method." << std::endl
                << "Standard deviation must be equal or greater than zero." << std::endl
                << std::endl;

      exit(1);
   }

   #endif

   int numberOfFreeParameters = getNumberOfFreeParameters();

   Vector<double> freeParameters(numberOfFreeParameters);

   double pi = 4.0*atan(1.0);

   double random1 = 0.0;
   double random2 = 0.0;

   for(int i = 0; i < numberOfFreeParameters; i++)
   {
      do
      {
         random1 = (double)rand()/(RAND_MAX+1.0);

      }while(random1 == 0.0);
     
      random2 = (double)rand()/(RAND_MAX+1.0);

      // Box-Muller transformation

      freeParameters[i] = mean + sqrt(-2.0*log(random1))*sin(2.0*pi*random2)*standardDeviation;
   }

   setFreeParameters(freeParameters);
}


// void initFreeParametersNormal(Vector<double>, Vector<double>) method

/// This method initializes all the free parameters in the neural newtork (biases and synaptic weiths + 
/// independent parameters) at random with values chosen from normal distributions with a given mean and a given 
/// standard deviation for each free parameter.
///
/// @param mean Vector of minimum initialization values.
/// @param standardDeviation Vector of maximum initialization values.

void MultilayerPerceptron::initFreeParametersNormal(Vector<double> mean, Vector<double> standardDeviation)
{
   int numberOfFreeParameters = getNumberOfFreeParameters();

   // Control sentence (if debug)

   #ifndef NDEBUG 

   int meanSize = mean.getSize();
   int standardDeviationSize = standardDeviation.getSize();

   if(meanSize != numberOfFreeParameters || standardDeviationSize != numberOfFreeParameters)
   {
      std::cerr << std::endl
                << "Flood Error: MultilayerPerceptron class." << std::endl
                << "void initFreeParametersNormal(Vector<double>, Vector<double>) method." << std::endl
                << "Mean and standard deviation sizes must be equal to number of free parameters." << std::endl 
                << std::endl;

      exit(1);
   }

   for(int i = 0; i < numberOfFreeParameters; i++)
   {
      if(standardDeviation[i] < 0.0)
      {
         std::cerr << std::endl
                   << "Flood Error: MultilayerPerceptron class." << std::endl
                   << "void initFreeParametersNormal(Vector<double>, Vector<double>) method." << std::endl
                   << "Standard deviations must equal or greater than zero." << std::endl
                   << std::endl;

         exit(1);
      }
   }

   #endif

   Vector<double> freeParameters(numberOfFreeParameters);

   double pi = 4.0*atan(1.0);
   
   double random1 = 0.0;
   double random2 = 0.0;

   for(int i = 0; i < numberOfFreeParameters; i++)
   {
      do
      {
         random1 = (double)rand()/(RAND_MAX+1.0);
      
      }while(random1 == 0.0);

      random2 = (double)rand()/(RAND_MAX+1.0);

      // Box-Muller transformation

      freeParameters[i] = mean[i] + sqrt(-2.0*log(random1))*sin(2.0*pi*random2)*standardDeviation[i];
   }

   setFreeParameters(freeParameters);
}


// void initFreeParametersNormal(Matrix<double>) method

/// This method initializes all the free parameters in the newtork (biases and synaptic weiths + independent 
/// parameters) at random with values chosen from normal distributions with a given mean and a given standard 
/// deviation for each free parameter.
/// All mean and standard deviation values are given from a single matrix.
/// The first row must contain the mean value for each free parameter.
/// The second row must contain the standard deviation value for each free parameter.
///
/// @param meanAndStandardDeviation Matrix of mean and standard deviation values.

void MultilayerPerceptron::initFreeParametersNormal(Matrix<double> meanAndStandardDeviation)
{
   int numberOfFreeParameters = getNumberOfFreeParameters();

   // Control sentence (if debug)

   #ifndef NDEBUG 

   int numberOfRows = meanAndStandardDeviation.getNumberOfRows();
   int numberOfColumns = meanAndStandardDeviation.getNumberOfColumns();

   if(numberOfRows != 2 || numberOfColumns != numberOfFreeParameters)
   {
      std::cerr << std::endl
                << "Flood Error: MultilayerPerceptron class." << std::endl 
                << "void initFreeParametersNormal(Matrix<double>) method." << std::endl
                << "Numbers of rows and columns must be two and number of free parameters, respectively." 
                << std::endl << std::endl;

      exit(1);
   }

   #endif

   Vector<double> mean = meanAndStandardDeviation.getRow(0);
   Vector<double> standardDeviation = meanAndStandardDeviation.getRow(1);

   #ifndef NDEBUG 

   for(int i = 0; i < numberOfFreeParameters; i++)
   {
      if(standardDeviation[i] < 0.0)
      {
         std::cerr << std::endl
                   << "Flood Error: MultilayerPerceptron class." << std::endl
                   << "void initFreeParametersNormal(Vector<double>, Vector<double>) method." << std::endl
                   << "Standard deviations must be equal or greater than zero." << std::endl
                   << std::endl;

         exit(1);
      }
   }

   #endif

   Vector<double> freeParameters(numberOfFreeParameters);

   double pi = 4.0*atan(1.0);
  
   double random1 = 0.0;
   double random2 = 0.0;

   for(int i = 0; i < numberOfFreeParameters; i++)
   {
      do
      {
         random1 = (double)rand()/(RAND_MAX+1.0);
     
      }while(random1 == 0.0);

      random2 = (double)rand()/(RAND_MAX+1.0);

      // Box-Muller transformation

      freeParameters[i] = mean[i] + sqrt(-2.0*log(random1))*sin(2.0*pi*random2)*standardDeviation[i];
   }

   setFreeParameters(freeParameters);
}


// void checkInput(Vector<double>) method

/// This method chechs whether the inputs to the neural network have the right size. 
/// If not, it displays an error message and exits the program. 
/// It also checks whether the input values are inside the range defined by the minimum and maximum values, and 
/// displays a warning message if they are outside.
///
/// @param input Set of inputs to the neural network.
///
/// @see calculateOutput.
/// @see calculateJacobian.

void MultilayerPerceptron::checkInput(Vector<double> input)
{
   // Control sentence (if debug)

   #ifndef NDEBUG 

   int size = input.getSize();

   if(size != numberOfInputs) 
   {
      std::cerr << std::endl 
                << "Flood Error: MultilayerPerceptron class." << std::endl
                << "void checkInput(Vector<double>) method." << std::endl
                << "Size of input must be equal to number of inputs." << std::endl
                << std::endl;

      exit(1);
   }   

   #endif

   // Check input

   if(displayOutOfRangeWarning)
   {
      for(int i = 0; i < numberOfInputs; i++)
      {
         if(input[i] < minimumOfInputVariables[i])
         {
            std::cout << std::endl 
                      << "Flood Warning: MultilayerPerceptron class." << std::endl
                      << "void checkInput(Vector<double>) method."  << std::endl
                      << "Input variable " << i << " is less than minimum." << std::endl;                  
         }
         if(input[i] > maximumOfInputVariables[i])
         {
            std::cout << std::endl 
                      << "Flood Warning: MultilayerPerceptron class." << std::endl
                      << "void checkInput(Vector<double>) method." << std::endl
                      << "Input variable " << i << " is greater than maximum." << std::endl;
         }
      }
   }
}


// Vector<double> preprocessInput(Vector<double>) method

/// This method preprocesses the inputs to the neural network in order to obtain a set of scaled input signals. 
///
/// @param input Set of inputs to the neural network.
///
/// @see calculateOutput.
/// @see calculateJacobian.

Vector<double> MultilayerPerceptron::preprocessInput(Vector<double> input)
{
   // Control sentence (if debug)

   #ifndef NDEBUG 

   int size = input.getSize();

   if(size != numberOfInputs) 
   {
      std::cerr << std::endl 
                << "Flood Error: MultilayerPerceptron class." << std::endl
                << "Vector<double> preprocessInput(Vector<double>) method." << std::endl
                << "Size of input must be equal to number of inputs." << std::endl
                << std::endl;

      exit(1);
   }   

   #endif

   // Preprocess input

   Vector<double> inputSignal(numberOfInputs);

   // Preprocess input depending on the pre and postprocessing method

   switch(preAndPostProcessingMethod)
   {
      case None:
      {
         inputSignal = input;
      }
      break;
                                     
      case MeanAndStandardDeviation:
      {
         for(int i = 0; i < numberOfInputs; i++)
         {            
            if(standardDeviationOfInputVariables[i] < 1e-9)
            {      
               if(display)
               {
                  std::cout << std::endl 
                            << "Flood Warning: MultilayerPerceptron class." << std::endl
                            << "Vector<double> preprocessInput(Vector<double>) method." << std::endl
                            << "Standard deviation of input variable " << i << " is zero." << std::endl
                            << "Those inputs won't be transformed." << std::endl;
               }
               
               inputSignal[i] = input[i];
            }      
            else
            {             
               inputSignal[i] = (input[i] - meanOfInputVariables[i])/standardDeviationOfInputVariables[i];
            }
         }       
      }
      break;

      case MinimumAndMaximum:
      {
         for(int i = 0; i < numberOfInputs; i++)
         {
            if(maximumOfInputVariables[i]-minimumOfInputVariables[i] < 1e-9)
            {      
               if(display)
               {
                  std::cout << "Flood Warning: MultilayerPerceptron class" << std::endl
                            << "Vector<double> preprocessInput(Vector<double>) method." << std::endl               
                            << "Minimum and maximum values of input variable " << i << " are equal." << std::endl
                            << "Those inputs won't be transformed." << std::endl;
               }
               
               inputSignal[i] = input[i];
            }      
            else
            {             
               inputSignal[i] = 2.0*(input[i] - minimumOfInputVariables[i])
               /(maximumOfInputVariables[i]-minimumOfInputVariables[i]) - 1.0;
            }
         }
      }            
      break;

      default:
      {
         std::cerr << "Flood Warning: MultilayerPerceptron class" << std::endl
                   << "Vector<double> preprocessInput(Vector<double>) method." << std::endl               
                   << "Unknown pre and postprocessing method." << std::endl;

         exit(1); 
      }
      break;
   }

   return(inputSignal);
}  


// Vector<double> calculateNetInputSignalToHiddenLayer(int, Vector<double>) method

/// This method returns the net input signal to every perceptron in the hidden layer as a function of the set of 
/// input signals to the neural network. 
///
/// @param index Hidden layer index. 
/// @param inputSignalToHiddenLayer Set of scaled input signals to the hidden layer with the previous index. 
///
/// @see calculateOutput.
/// @see calculateJacobian.

Vector<double> MultilayerPerceptron::calculateNetInputSignalToHiddenLayer(int index, 
Vector<double> inputSignalToHiddenLayer)
{
   // Control sentence (if debug)

//   #ifndef NDEBUG 

//   int size = inputSignalToHiddenLayer.getSize();

//   if(size != numberOfInputs) 
//   {
//      std::cerr << std::endl 
//                << "Flood Error: MultilayerPerceptron class." << std::endl
//                << "Vector<double> calculateNetInputSignalToHiddenLayer(int, Vector<double>) method." 
//                << std::endl
//                << "Size of input must be equal to number of inputs." << std::endl
//                << std::endl;

//      exit(1);
//   }   

//   #endif

   // Calculate net input signal to hidden layer

   Vector<double> netInputSignalToHiddenLayer(numbersOfHiddenNeurons[index]);
   
   // Calculate net input signal to hidden layer

   for(int j = 0; j < numbersOfHiddenNeurons[index]; j++)
   {         
      netInputSignalToHiddenLayer[j] = hiddenLayers[index][j].calculateNetInputSignal(inputSignalToHiddenLayer);
   }             
   
   return(netInputSignalToHiddenLayer);
}


// Vector<double> calculateOutputSignalFromHiddenLayer(int, Vector<double>) method

/// This method returns the output signals from every perceptron in the hidden layer
/// as a function of their net input signals. 
///
/// @param index Hidden layer index. 
/// @param netInputSignalToHiddenLayer Net input signal to every neuron in the hidden layer with the previous index. 
///
/// @see calculateOutput.
/// @see calculateJacobian.

Vector<double> MultilayerPerceptron::calculateOutputSignalFromHiddenLayer(int index, 
Vector<double> netInputSignalToHiddenLayer)
{
   // Control sentence (if debug)

   #ifndef NDEBUG 

   int size = netInputSignalToHiddenLayer.getSize();

   if(size != numbersOfHiddenNeurons[index]) 
   {
      std::cerr << std::endl 
                << "Flood Error: MultilayerPerceptron class." << std::endl
                << "Vector<double> calculateOutputSignalFromHiddenLayer(int, Vector<double>) method." 
                << std::endl
                << "Size must be equal to number of hidden neurons." << std::endl
                << std::endl;

      exit(1);
   }   

   #endif

   // Calculate output signal from hidden layer

   Vector<double> outputSignalFromHiddenLayer(numbersOfHiddenNeurons[index]);

   // Calculate output signal from hidden layer

   for(int j = 0; j < numbersOfHiddenNeurons[index]; j++)
   {
      outputSignalFromHiddenLayer[j] = hiddenLayers[index][j].calculateOutputSignal(netInputSignalToHiddenLayer[j]);
   }

   return(outputSignalFromHiddenLayer);
}  


// Vector<double> calculateOutputSignalDerivativeFromHiddenLayer(int, Vector<double>) method

/// This method returns the activation derivative from every perceptron in the hidden layer
/// as a function of their net input signals. 
///
/// @param index Hidden layer index. 
/// @param netInputSignalToHiddenLayer Net input signal to every neuron in the hidden layer with the previous index. 
///
/// @see calculateOutput.
/// @see calculateJacobian.

Vector<double> MultilayerPerceptron::calculateOutputSignalDerivativeFromHiddenLayer(int index,
Vector<double> netInputSignalToHiddenLayer)
{
   // Control sentence (if debug)

//   #ifndef NDEBUG 

//   int size = netInputSignalToHiddenLayer.getSize();

//   if(size != numbersOfHiddenNeurons) 
//   {
//      std::cerr << std::endl 
//                << "Flood Error: MultilayerPerceptron class." << std::endl
//                << "Vector<double> calculateOutputSignalDerivativeFromHiddenLayer(int, Vector<double>) method." 
//                << std::endl
//                << "Size must be equal to number of hidden neurons." << std::endl
//                << std::endl;

//      exit(1);
//   }   

//   #endif

   // Calculate output signal derivative from hidden layer

   Vector<double> outputSignalDerivativeFromHiddenLayer(numbersOfHiddenNeurons[index]);

   // Calculate output signal derivative from hidden layer

   for(int j = 0; j < numbersOfHiddenNeurons[index]; j++)
   {          
      outputSignalDerivativeFromHiddenLayer[j] 
      = hiddenLayers[index][j].calculateOutputSignalDerivative(netInputSignalToHiddenLayer[j]);
      
   }
   
   return(outputSignalDerivativeFromHiddenLayer);
}


// Vector<double> calculateNetInputSignalToOutputLayer(Vector<double>) method

/// This method returns the net input signal to every perceptron in the ouptut layer
/// as a function of the output signals from the hidden neurons. 
///
/// @param outputSignalFromHiddenLayer Output signal from every neuron in the hidden layer. 
///
/// @see calculateOutput.
/// @see calculateJacobian.

Vector<double> MultilayerPerceptron::calculateNetInputSignalToOutputLayer(Vector<double> outputSignalFromHiddenLayer)
{
   // Control sentence (if debug)

//   #ifndef NDEBUG 

//   int size = outputSignalFromHiddenLayer.getSize();

//   if(size != numbersOfHiddenNeurons) 
//   {
//      std::cerr << std::endl 
//                << "Flood Error: MultilayerPerceptron class." << std::endl
//                << "Vector<double> calculateNetInputSignalToOutputLayer(Vector<double>) method." << std::endl
//                << "Size must be equal to number of hidden neurons." << std::endl
//                << std::endl;

//      exit(1);
//   }   

//   #endif

   // Calculate net input signal to hidden layer

   Vector<double> netInputSignalToOutputLayer(numberOfOutputs);            
               
   // Calculate net input signal to output layer

   for(int i = 0; i < numberOfOutputs; i++)
   {
      netInputSignalToOutputLayer[i] = outputLayer[i].calculateNetInputSignal(outputSignalFromHiddenLayer);
   }
   
   return(netInputSignalToOutputLayer);                  
}


// Vector<double> calculateOutputSignal(Vector<double>) method

/// This method returns the output signal from every perceptron in the ouptut layer as a function of their net 
/// input signals. 
///
/// @param netInputSignalToOutputLayer Net input signal to every neuron in the output layer. 
///
/// @see calculateOutput.
/// @see calculateJacobian.

Vector<double> MultilayerPerceptron::calculateOutputSignal(Vector<double> netInputSignalToOutputLayer)
{
   // Control sentence (if debug)

   #ifndef NDEBUG 

   int size = netInputSignalToOutputLayer.getSize();

   if(size != numberOfOutputs) 
   {
      std::cerr << std::endl 
                << "Flood Error: MultilayerPerceptron class." << std::endl
                << "Vector<double> calculateOutputSignal(Vector<double>) method." << std::endl
                << "Size must be equal to number of outputs." << std::endl
                << std::endl;

      exit(1);
   }   

   #endif

   // Calculate output signal

   Vector<double> outputSignal(numberOfOutputs);

   // Get output signal from output layer

   for(int i = 0; i < numberOfOutputs; i++)
   {
      outputSignal[i] = outputLayer[i].calculateOutputSignal(netInputSignalToOutputLayer[i]);
   }

   return(outputSignal);
}  


// Vector<double> calculateOutputSignalDerivative(Vector<double>) method

/// This method returns the output signal derivative from every perceptron in the ouptut layer as a function of 
/// their net input signals. 
///
/// @param netInputSignalToOutputLayer Net input signal to every neuron in the output layer. 
///
/// @see calculateOutput.
/// @see calculateJacobian.

Vector<double> MultilayerPerceptron::calculateOutputSignalDerivative(Vector<double> netInputSignalToOutputLayer)
{
   // Control sentence (if debug)

   #ifndef NDEBUG 

   int size = netInputSignalToOutputLayer.getSize();

   if(size != numberOfOutputs) 
   {
      std::cerr << std::endl 
                << "Flood Error: MultilayerPerceptron class." << std::endl
                << "Vector<double> calculateOutputSignalDerivative(Vector<double>) method." << std::endl
                << "Size must be equal to number of outputs." << std::endl
                << std::endl;

      exit(1);
   }   

   #endif

   // Calculate output signal derivative

   Vector<double> outputSignalDerivative(numberOfOutputs);

   // Get output signal derivative from output layer

   for(int i = 0; i < numberOfOutputs; i++)
   {
      outputSignalDerivative[i] = outputLayer[i].calculateOutputSignalDerivative(netInputSignalToOutputLayer[i]);
   }

   return(outputSignalDerivative);
}  


// Vector<double> postprocessOutputSignal(Vector<double>) method

/// This method postprocesses the output signals from the neural network in order to obtain a set of unscaled 
/// outputs. 
///
/// @param outputSignal Set of output signals from the neural network.
///
/// @see calculateOutput.
/// @see calculateJacobian.

Vector<double> MultilayerPerceptron::postprocessOutputSignal(Vector<double> outputSignal)
{
   // Control sentence (if debug)

   #ifndef NDEBUG 

   int size = outputSignal.getSize();

   if(size != numberOfOutputs) 
   {
      std::cerr << std::endl 
                << "Flood Error: MultilayerPerceptron class." << std::endl
                << "Vector<double> postprocessOutputSignal(Vector<double>) method." << std::endl
                << "Size must be equal to number of outputs." << std::endl
                << std::endl;

      exit(1);
   }   

   #endif

   // Postprocess output signal

   Vector<double> output(numberOfOutputs);

   switch(preAndPostProcessingMethod)
   {
      case None:
      {
         output = outputSignal;
      }// end none
      break;

      case MeanAndStandardDeviation:
      {
         for(int i = 0; i < numberOfOutputs; i++)
         {
            if(standardDeviationOfOutputVariables[i] < 1e-9)
            {      
               if(display)
               {
                  std::cout << "Flood Warning: MultilayerPerceptron class." << std::endl
                            << "Vector<double> postprocessOutputSignal(Vector<double>) method." << std::endl               
                            << "Standard deviation of output variable " << i << " is zero." << std::endl 
                            << "Those output signals won't be transformed." << std::endl;
               }
               
               output[i] = outputSignal[i];
            }      
            else
            {
               output[i] = outputSignal[i]*standardDeviationOfOutputVariables[i] 
               + meanOfOutputVariables[i];
            }
         }
      }// end mean and standard deviation
      break;

      case MinimumAndMaximum:
      {
         for(int i = 0; i < numberOfOutputs; i++)
         {
            if(maximumOfOutputVariables[i]-minimumOfOutputVariables[i] < 1e-9)
            {      
               if(display)
               {
                  std::cout << "Flood Warning: MultilayerPerceptron class." << std::endl
                            << "Vector<double> postprocessOutputSignal(Vector<double>) method." << std::endl               
                            << "Minimum and maximum values of output variable " << i << " are equal." << std::endl
                            << "Those output signals won't be transformed."
                            << std::endl;
               }

               output[i] = outputSignal[i];
            }      
            else
            {
               output[i] = 0.5*(outputSignal[i] + 1.0)*(maximumOfOutputVariables[i]-minimumOfOutputVariables[i]) 
               + minimumOfOutputVariables[i]; 
            }
         }
      }// end minimum and maximum
      break;       

      default:
      {
         std::cout << "Flood Warning: MultilayerPerceptron class." << std::endl
                   << "Vector<double> postprocessOutputSignal(Vector<double>) method." << std::endl               
                   << "Unknown pre and postprocessing method." << std::endl
                   << std::endl;
 
         exit(1);
      }// end default
      break;

   }// end swithch

   return(output);
}  


// Vector<double> applyLowerAndUpperBounds(Vector<double>) method

/// This method postprocesses the outputs from the neural network in order to obtain a new set of outputs falling 
/// in the range defined by the lower and the upper bounds. 
///
/// @param output Set of outputs from the neural network.
///
/// @see calculateOutput.
/// @see calculateJacobian.

Vector<double> MultilayerPerceptron::applyLowerAndUpperBounds(Vector<double> output)
{            
   // Control sentence (if debug)

   #ifndef NDEBUG 

   int size = output.getSize();

   if(size != numberOfOutputs) 
   {
      std::cerr << std::endl 
                << "Flood Error: MultilayerPerceptron class." << std::endl
                << "Vector<double> applyLowerAndUpperBounds(Vector<double>) method." << std::endl
                << "Size must be equal to number of outputs." << std::endl
                << std::endl;

      exit(1);
   }   

   #endif

   // Apply lower and upper bounds

   Vector<double> newOutput(numberOfOutputs);

   for(int i = 0; i < numberOfOutputs; i++)
   {
      if(output[i] < lowerBoundOfOutputVariables[i])
      {
         newOutput[i] = lowerBoundOfOutputVariables[i];
      }
      else if(output[i] > upperBoundOfOutputVariables[i])
      {
         newOutput[i] = upperBoundOfOutputVariables[i];
      }
      else
      {
         newOutput[i] =  output[i];
      }
   }

   return(newOutput);
}  


// Vector<double> calculateOutput(Vector<double>) method

/// This method calculates the set of outputs from the neural network in response to a set of inputs, when no 
/// boundary conditions ar imposed.
///
/// @param input Set of inputs to the neural network.
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
/// applyLowerAndUpperBounds(Vector<double>).
/// @see calculateJacobian(Vector<double>).
/// @see calculateJacobianTest(Vector<double>).

Vector<double> MultilayerPerceptron::calculateOutput(Vector<double> input)
{
   // Control sentence (if debug)

   #ifndef NDEBUG 

   int size = input.getSize();

   if(size != numberOfInputs) 
   {
      std::cerr << std::endl 
                << "Flood Error: MultilayerPerceptron class." << std::endl
                << "Vector<double> calculateOutput(Vector<double>) method." << std::endl
                << "Size must be equal to number of inputs." << std::endl
                << std::endl;

      exit(1);
   }   

   #endif

   int numberOfHiddenLayers = getNumberOfHiddenLayers();

   // Calculate output

   Vector<double> output(numberOfOutputs);
   
   // Check input for size and maximum and minimum values
 
   checkInput(input);
 
   // Preprocess input to obtain input signal to the neural network 

   Vector<double> inputSignal = preprocessInput(input);

   // Calculate net input and output signals from all hidden layers

   Vector< Vector<double> > netInputSignalToHiddenLayer(numberOfHiddenLayers);
   netInputSignalToHiddenLayer[0] = calculateNetInputSignalToHiddenLayer(0, inputSignal);

   Vector< Vector<double> > outputSignalFromHiddenLayer(numberOfHiddenLayers);
   outputSignalFromHiddenLayer[0] = calculateOutputSignalFromHiddenLayer(0, netInputSignalToHiddenLayer[0]);
   
   for(int i = 1; i < numberOfHiddenLayers; i++)
   {
      netInputSignalToHiddenLayer[i] = calculateNetInputSignalToHiddenLayer(i, outputSignalFromHiddenLayer[i-1]);         
      outputSignalFromHiddenLayer[i] = calculateOutputSignalFromHiddenLayer(i, netInputSignalToHiddenLayer[i]);         
   }

   // Get net input signal to output layer

   Vector<double> netInputSignalToOutputLayer 
   = calculateNetInputSignalToOutputLayer(outputSignalFromHiddenLayer[numberOfHiddenLayers-1]);

   // Get output signal from output layer

   Vector<double> outputSignal = calculateOutputSignal(netInputSignalToOutputLayer);
  
   // Postprocess output signal from the neural network to get output

   output = postprocessOutputSignal(outputSignal);

   // Apply lower and upper bounds

   output = applyLowerAndUpperBounds(output);

   return(output);
}


// Matrix<double> calculateJacobian(Vector<double>) method

/// This method returns the the Jacobian Matrix of the neural network for a set of inputs, corresponding to the 
/// point in input space at which the Jacobian Matrix is to be found. It uses the back-propagation method.
///
/// @param input Set of inputs to the neural network.
///
/// @see calculateOutput(Vector<double>).
/// @see calculateJacobianTest(Vector<double>).
///
/// @todo This method is only implemented for multilayer perceptrons with one hidden layer. 
/// An extension of this algorithm so as to consider any number of hidden layers must be implemented.

Matrix<double> MultilayerPerceptron::calculateJacobian(Vector<double> input)
{   
   // Check input for size and maximum and minimum values
 
   checkInput(input);

   int numberOfHiddenLayers = getNumberOfHiddenLayers();

   // Calculate jacobian

   Matrix<double> jacobian(numberOfOutputs, numberOfInputs, 0.0);

   // Preprocess input to obtain input signal to the neural network 

   Vector<double> inputSignal = preprocessInput(input);

   // Calculate net input output and output signal derivative from all hidden layers

   Vector< Vector<double> > netInputSignalToHiddenLayer(numberOfHiddenLayers);
   netInputSignalToHiddenLayer[0] = calculateNetInputSignalToHiddenLayer(0, inputSignal);

   Vector< Vector<double> > outputSignalFromHiddenLayer(numberOfHiddenLayers);
   outputSignalFromHiddenLayer[0] = calculateOutputSignalFromHiddenLayer(0, netInputSignalToHiddenLayer[0]);

   Vector< Vector<double> > outputSignalDerivativeFromHiddenLayer(numberOfHiddenLayers);
   outputSignalDerivativeFromHiddenLayer[0] 
   = calculateOutputSignalDerivativeFromHiddenLayer(0, netInputSignalToHiddenLayer[0]);
   
   for(int i = 1; i < numberOfHiddenLayers; i++)
   {
      netInputSignalToHiddenLayer[i] = calculateNetInputSignalToHiddenLayer(i, outputSignalFromHiddenLayer[i-1]);

      outputSignalFromHiddenLayer[i] = calculateOutputSignalFromHiddenLayer(i, netInputSignalToHiddenLayer[i]);

      outputSignalDerivativeFromHiddenLayer[i] 
      = calculateOutputSignalDerivativeFromHiddenLayer(i, netInputSignalToHiddenLayer[i]);         
   }

   // Get net input signal to output layer

   Vector<double> netInputSignalToOutputLayer 
   = calculateNetInputSignalToOutputLayer(outputSignalFromHiddenLayer[numberOfHiddenLayers-1]);

   // Calculate output signal derivative

   Vector<double> outputSignalDerivative = calculateOutputSignalDerivative(netInputSignalToOutputLayer);

   // Calculate Jacobian matrix
  
   double deltaOutput = 0.0;
   Vector<double> deltaHidden(numberOfHiddenLayers, 0.0);
   
   double sum = 0.0;

   for(int k = 0; k < numberOfOutputs; k++)
   {
      // delta output
            
      deltaOutput = outputSignalDerivative[k];
    
      for(int i = 0; i < numberOfInputs; i++)
      {      
         for(int h = numberOfHiddenLayers-1; h >= 0; h--)
         {
            for(int j = 0; j < numbersOfHiddenNeurons[h]; j++)
            {
               sum = 0.0;
                
               sum += outputLayer[k].getSingleSynapticWeight(j)*deltaOutput;
            
               deltaHidden[h] = outputSignalDerivativeFromHiddenLayer[h][j]*sum;

               // Jacobian matrix 
  
               if(h == 0)
               {  
                  jacobian[k][i] += hiddenLayers[h][j].getSingleSynapticWeight(i)*deltaHidden[h];
               }
               else
               {

               } 
            }
         }
      }
   }

   return(jacobian);
}


// Matrix<double> calculateJacobianTest(Vector<double>) method

/// This method returns the the Jacobian matrix of the neural network for a set of inputs corresponding to the 
/// point in input space at which the Jacobian Matrix is to be found. It uses numerical differentiation.
///
/// @param input: Set of inputs to the neural network.
///
/// @see numericalEpsilon.
/// @see calculateOutput(Vector<double>).
/// @see calculateJacobian(Vector<double>).

Matrix<double> MultilayerPerceptron::calculateJacobianTest(Vector<double> input)
{
   // Check input for size and maximum and minimum values
 
   checkInput(input);

   int numberOfHiddenLayers = getNumberOfHiddenLayers();

   // Calculate jacobian

   Matrix<double> jacobian(numberOfOutputs, numberOfInputs);

   // Preprocess input to obtain input signal to the neural network 

   Vector<double> inputSignal = preprocessInput(input);

   switch(numericalDifferentiationMethod)   
   {
      case ForwardDifferences:
      {
         // Calculate output signal for input signal

         Vector< Vector<double> > netInputSignalToHiddenLayer(numberOfHiddenLayers);
         netInputSignalToHiddenLayer[0] = calculateNetInputSignalToHiddenLayer(0, inputSignal);

         Vector< Vector<double> > outputSignalFromHiddenLayer(numberOfHiddenLayers);
         outputSignalFromHiddenLayer[0] = calculateOutputSignalFromHiddenLayer(0, netInputSignalToHiddenLayer[0]);
   
         for(int i = 1; i < numberOfHiddenLayers; i++)
         {
            netInputSignalToHiddenLayer[i] 
            = calculateNetInputSignalToHiddenLayer(i, outputSignalFromHiddenLayer[i-1]);         
      
            outputSignalFromHiddenLayer[i] 
            = calculateOutputSignalFromHiddenLayer(i, netInputSignalToHiddenLayer[i]);         
         }

         Vector<double> netInputSignalToOutputLayer 
         = calculateNetInputSignalToOutputLayer(outputSignalFromHiddenLayer[numberOfHiddenLayers-1]);
         Vector<double> outputSignal = calculateOutputSignal(netInputSignalToOutputLayer);
         // Calculate output signal for input signal forward

         Vector<double> inputSignalForward = inputSignal;         Vector< Vector<double> > netInputSignalToHiddenLayerForward(numberOfHiddenLayers);  
         Vector< Vector<double> > outputSignalFromHiddenLayerForward(numberOfHiddenLayers);  
         Vector<double> netInputSignalToOutputLayerForward(numberOfOutputs);  
         Vector<double> outputSignalForward(numberOfOutputs); 

         for(int i = 0; i < numberOfInputs; i++)
         {
            // Perturbate input signals

            inputSignalForward[i] += numericalEpsilon;
            // Calculate output signal             

            netInputSignalToHiddenLayerForward[0] = calculateNetInputSignalToHiddenLayer(0, inputSignalForward);
            outputSignalFromHiddenLayerForward[0] 
            = calculateOutputSignalFromHiddenLayer(0, netInputSignalToHiddenLayerForward[0]);
   
            for(int j = 1; j < numberOfHiddenLayers; j++)
            {
               netInputSignalToHiddenLayerForward[j] 
               = calculateNetInputSignalToHiddenLayer(j, outputSignalFromHiddenLayerForward[j-1]);         
      
               outputSignalFromHiddenLayerForward[j] 
               = calculateOutputSignalFromHiddenLayer(j, netInputSignalToHiddenLayerForward[j]);         
            }

            netInputSignalToOutputLayerForward 
            = calculateNetInputSignalToOutputLayer(outputSignalFromHiddenLayerForward[numberOfHiddenLayers-1]);

            outputSignalForward = calculateOutputSignal(netInputSignalToOutputLayerForward);
            // Calculate jacobian

            for(int k = 0; k < numberOfOutputs; k++)
            {
               jacobian[k][i] = (outputSignalForward[k] - outputSignal[k])/numericalEpsilon;
            }

            // Restart original input signals

            inputSignalForward[i] -= numericalEpsilon;
         }
      }                               
      break;

      case CentralDifferences:
      {
         Vector<double> inputSignalForward = inputSignal;         Vector<double> inputSignalBackward = inputSignal;
         Vector< Vector<double> > netInputSignalToHiddenLayerForward(numberOfHiddenLayers);  
         Vector< Vector<double> > netInputSignalToHiddenLayerBackward(numberOfHiddenLayers);  

         Vector< Vector<double> > outputSignalFromHiddenLayerForward(numberOfHiddenLayers);  
         Vector< Vector<double> > outputSignalFromHiddenLayerBackward(numberOfHiddenLayers);  

         Vector<double> netInputSignalToOutputLayerForward(numberOfOutputs);  
         Vector<double> netInputSignalToOutputLayerBackward(numberOfOutputs);  

         Vector<double> outputSignalForward(numberOfOutputs); 
         Vector<double> outputSignalBackward(numberOfOutputs); 

         for(int i = 0; i < numberOfInputs; i++)
         {
            // Perturbate input signals

            inputSignalForward[i] += numericalEpsilon;
            inputSignalBackward[i] -= numericalEpsilon;

            // Calculate output signals             

            netInputSignalToHiddenLayerForward[0] = calculateNetInputSignalToHiddenLayer(0, inputSignalForward);
            netInputSignalToHiddenLayerBackward[0] = calculateNetInputSignalToHiddenLayer(0, inputSignalBackward);
            outputSignalFromHiddenLayerForward[0] 
            = calculateOutputSignalFromHiddenLayer(0, netInputSignalToHiddenLayerForward[0]);

            outputSignalFromHiddenLayerBackward[0] 
            = calculateOutputSignalFromHiddenLayer(0, netInputSignalToHiddenLayerBackward[0]);
   
            for(int j = 1; j < numberOfHiddenLayers; j++)
            {
               netInputSignalToHiddenLayerForward[j] 
               = calculateNetInputSignalToHiddenLayer(j, outputSignalFromHiddenLayerForward[j-1]);         

               netInputSignalToHiddenLayerBackward[j] 
               = calculateNetInputSignalToHiddenLayer(j, outputSignalFromHiddenLayerBackward[j-1]);         
      
               outputSignalFromHiddenLayerForward[j] 
               = calculateOutputSignalFromHiddenLayer(j, netInputSignalToHiddenLayerForward[j]);         

               outputSignalFromHiddenLayerBackward[j] 
               = calculateOutputSignalFromHiddenLayer(j, netInputSignalToHiddenLayerBackward[j]);         
            }

            netInputSignalToOutputLayerForward 
            = calculateNetInputSignalToOutputLayer(outputSignalFromHiddenLayerForward[numberOfHiddenLayers-1]);
            netInputSignalToOutputLayerBackward 
            = calculateNetInputSignalToOutputLayer(outputSignalFromHiddenLayerBackward[numberOfHiddenLayers-1]);

            outputSignalForward = calculateOutputSignal(netInputSignalToOutputLayerForward);
            outputSignalBackward = calculateOutputSignal(netInputSignalToOutputLayerBackward);
            // Calculate jacobian

            for(int k = 0; k < numberOfOutputs; k++)
            {
               jacobian[k][i] = (outputSignalForward[k] - outputSignalBackward[k])/(2.0*numericalEpsilon);
            }

            // Restart original input signals

            inputSignalForward[i] -= numericalEpsilon;
            inputSignalBackward[i] += numericalEpsilon;
         }

      }// end central differences
      break;

      default:
      {
         std::cerr << std::endl 
                   << "Flood Error: MultilayerPerceptron class." << std::endl
                   << "Matrix<double> calculateJacobianTest(Vector<double>) method." << std::endl
                   << "Unknown numerical differentiation method." << std::endl
                   << std::endl;
 
         exit(1);
      }// end default
      break;

   }// end numerical differentiation method switch

   return(jacobian);
}


// Matrix<double> calculateSensitivityTest(Vector<double> input)

/// @todo The calculus of the sensitivity matrix is only implemented using numerical differentiation.
/// A back-propagation algorithm for that derivatives should be developed. 

Matrix<double> MultilayerPerceptron::calculateSensitivityTest(Vector<double> input)
{
   // Check input for size and maximum and minimum values
 
   checkInput(input);

   int numberOfHiddenLayers = getNumberOfHiddenLayers();

   int numberOfNeuralParameters = getNumberOfNeuralParameters();

   // Calculate sensitivity

   Matrix<double> sensitivity(numberOfOutputs, numberOfNeuralParameters);

   // Preprocess input to obtain input signal to the neural network 

   Vector<double> inputSignal = preprocessInput(input);

   switch(numericalDifferentiationMethod)   
   {
      case ForwardDifferences:
      {
         // Calculate output signal for input signal and neural parameters

         Vector<double> neuralParameters = getNeuralParameters();

         Vector< Vector<double> > netInputSignalToHiddenLayer(numberOfHiddenLayers);
         netInputSignalToHiddenLayer[0] = calculateNetInputSignalToHiddenLayer(0, inputSignal);

         Vector< Vector<double> > outputSignalFromHiddenLayer(numberOfHiddenLayers);
         outputSignalFromHiddenLayer[0] = calculateOutputSignalFromHiddenLayer(0, netInputSignalToHiddenLayer[0]);
   
         for(int i = 1; i < numberOfHiddenLayers; i++)
         {
            netInputSignalToHiddenLayer[i] 
            = calculateNetInputSignalToHiddenLayer(i, outputSignalFromHiddenLayer[i-1]);         
      
            outputSignalFromHiddenLayer[i] 
            = calculateOutputSignalFromHiddenLayer(i, netInputSignalToHiddenLayer[i]);         
         }

         Vector<double> netInputSignalToOutputLayer 
         = calculateNetInputSignalToOutputLayer(outputSignalFromHiddenLayer[numberOfHiddenLayers-1]);

         Vector<double> outputSignal = calculateOutputSignal(netInputSignalToOutputLayer);

         // Calculate output signal for input signal and neural parameters forward

         Vector<double> neuralParametersForward = neuralParameters;

         Vector< Vector<double> > netInputSignalToHiddenLayerForward(numberOfHiddenLayers);  
         Vector< Vector<double> > outputSignalFromHiddenLayerForward(numberOfHiddenLayers);  
         Vector<double> netInputSignalToOutputLayerForward(numberOfOutputs);  
         Vector<double> outputSignalForward(numberOfOutputs); 

         for(int i = 0; i < numberOfNeuralParameters; i++)
         {
            // Perturbate single neural parameter

            neuralParametersForward[i] += numericalEpsilon;
            setNeuralParameters(neuralParametersForward);

            // Calculate output signal             

            netInputSignalToHiddenLayerForward[0] = calculateNetInputSignalToHiddenLayer(0, inputSignal);
            outputSignalFromHiddenLayerForward[0] 
            = calculateOutputSignalFromHiddenLayer(0, netInputSignalToHiddenLayerForward[0]);
   
            for(int j = 1; j < numberOfHiddenLayers; j++)
            {
               netInputSignalToHiddenLayerForward[j] 
               = calculateNetInputSignalToHiddenLayer(j, outputSignalFromHiddenLayerForward[j-1]);         
      
               outputSignalFromHiddenLayerForward[j] 
               = calculateOutputSignalFromHiddenLayer(j, netInputSignalToHiddenLayerForward[j]);         
            }

            netInputSignalToOutputLayerForward 
            = calculateNetInputSignalToOutputLayer(outputSignalFromHiddenLayerForward[numberOfHiddenLayers-1]);

            outputSignalForward = calculateOutputSignal(netInputSignalToOutputLayerForward);
            // Calculate sensitivity

            for(int k = 0; k < numberOfOutputs; k++)
            {
               sensitivity[k][i] = (outputSignalForward[k] - outputSignal[k])/numericalEpsilon;
            }

            // Restart original neural parameters

            neuralParametersForward[i] -= numericalEpsilon;
            setNeuralParameters(neuralParametersForward);
         }
      }                               
      break;

      case CentralDifferences:
      {
         Vector<double> neuralParametersForward = getNeuralParameters();
         Vector<double> neuralParametersBackward = getNeuralParameters();

         Vector< Vector<double> > netInputSignalToHiddenLayerForward(numberOfHiddenLayers);  
         Vector< Vector<double> > netInputSignalToHiddenLayerBackward(numberOfHiddenLayers);  

         Vector< Vector<double> > outputSignalFromHiddenLayerForward(numberOfHiddenLayers);  
         Vector< Vector<double> > outputSignalFromHiddenLayerBackward(numberOfHiddenLayers);  

         Vector<double> netInputSignalToOutputLayerForward(numberOfOutputs);  
         Vector<double> netInputSignalToOutputLayerBackward(numberOfOutputs);  

         Vector<double> outputSignalForward(numberOfOutputs); 
         Vector<double> outputSignalBackward(numberOfOutputs); 

         for(int i = 0; i < numberOfNeuralParameters; i++)
         {
            // Perturbate neural parameters

            neuralParametersForward[i] += numericalEpsilon;
            setNeuralParameters(neuralParametersForward);

            // Calculate output signal             

            netInputSignalToHiddenLayerForward[0] = calculateNetInputSignalToHiddenLayer(0, inputSignal);

            outputSignalFromHiddenLayerForward[0] 
            = calculateOutputSignalFromHiddenLayer(0, netInputSignalToHiddenLayerForward[0]);
   
            for(int j = 1; j < numberOfHiddenLayers; j++)
            {
               netInputSignalToHiddenLayerForward[j] 
               = calculateNetInputSignalToHiddenLayer(j, outputSignalFromHiddenLayerForward[j-1]);         
      
               outputSignalFromHiddenLayerForward[j] 
               = calculateOutputSignalFromHiddenLayer(j, netInputSignalToHiddenLayerForward[j]);         
            }

            netInputSignalToOutputLayerForward 
            = calculateNetInputSignalToOutputLayer(outputSignalFromHiddenLayerForward[numberOfHiddenLayers-1]);

            outputSignalForward = calculateOutputSignal(netInputSignalToOutputLayerForward);

            // Restart original neural parameters

            neuralParametersForward[i] -= numericalEpsilon;
            setNeuralParameters(neuralParametersForward);

            // Perturbate neural parameters

            neuralParametersBackward[i] -= numericalEpsilon;
            setNeuralParameters(neuralParametersBackward);

            // Calculate output signal             

            netInputSignalToHiddenLayerBackward[0] = calculateNetInputSignalToHiddenLayer(0, inputSignal);
            outputSignalFromHiddenLayerBackward[0] 
            = calculateOutputSignalFromHiddenLayer(0, netInputSignalToHiddenLayerBackward[0]);
   
            for(int j = 1; j < numberOfHiddenLayers; j++)
            {
               netInputSignalToHiddenLayerBackward[j] 
               = calculateNetInputSignalToHiddenLayer(j, outputSignalFromHiddenLayerBackward[j-1]);         
      
               outputSignalFromHiddenLayerBackward[j] 
               = calculateOutputSignalFromHiddenLayer(j, netInputSignalToHiddenLayerBackward[j]);         
            }

            netInputSignalToOutputLayerBackward 
            = calculateNetInputSignalToOutputLayer(outputSignalFromHiddenLayerBackward[numberOfHiddenLayers-1]);

            outputSignalBackward = calculateOutputSignal(netInputSignalToOutputLayerBackward);
            // Restart original neural parameters

            neuralParametersBackward[i] += numericalEpsilon;
            setNeuralParameters(neuralParametersBackward);       

            // Calculate sensitivity

            for(int k = 0; k < numberOfOutputs; k++)
            {
               sensitivity[k][i] = (outputSignalForward[k] - outputSignalBackward[k])/(2.0*numericalEpsilon);
            }

         }

      }// end central differences
      break;

      default:
      {
         std::cerr << std::endl 
                   << "Flood Error: MultilayerPerceptron class." << std::endl
                   << "Matrix<double> calculateSensitivityTest(Vector<double>) method." << std::endl
                   << "Unknown numerical differentiation method." << std::endl
                   << std::endl;
 
         exit(1);
      }// end default
      break;

   }// end numerical differentiation method switch

   return(sensitivity);
}


// void print(void) method

/// This method prints to the screen the members of a multilayer perceptron object:
///
/// <ul>
/// <li> Number of inputs.
/// <li> Number of hidden neurons.
/// <li> Number of outputs.
/// <li> Hidden layer activation function.
/// <li> Output layer activation function.
/// <li> Name of input vaviables. 
/// <li> Name of output variables.
/// <li> Units of input vaviables. 
/// <li> Units of output variables.
/// <li> Description of input vaviables. 
/// <li> Description of output variables.
/// <li> Mean of input variables.
/// <li> Standard deviation of input variables.
/// <li> Mean of output variables.
/// <li> Standard deviation of output variables.
/// <li> Minimum of input variables.
/// <li> Maximum of input variables.
/// <li> Minimum of output variables.
/// <li> Maximum of output variables.
/// <li> Lower bound of output variables.
/// <li> Upper bound of output variables.
/// <li> Biases and synaptic weights.
/// <li> Pre and postprocessing method.
/// <li> Number of independent parameters.
/// <li> Name of independent parameters.
/// <li> Units of independent parameters.
/// <li> Description of independent parameters.
/// <li> Mean of independent parameters.
/// <li> Standard deviation of independent parameters.
/// <li> Minimum of independent parameters.
/// <li> Maximum of independent parameters.
/// <li> Lower bound of independent parameters.
/// <li> Upper bound of independent parameters.
/// <li> Independent parameters.
/// <li> Pre and postprocessing method.
/// <li> Numerical epsilon method.
/// <li> Numerical differentiation method.
/// <li> Display out of range warning.
/// <li> Display.
/// <li> Numerical epsilon.
/// </ul> 

void MultilayerPerceptron::print(void)
{

   // Print header

   std::cout << std::endl
             << "Flood Neural Network. Multilayer Perceptron Object." << std::endl;

   // Print network architecture

   int numberOfHiddenLayers = getNumberOfHiddenLayers();

   std::cout << "Number of inputs:" << std::endl 
             << numberOfInputs << std::endl
             << "Number of hidden layers:" << std::endl
             << numberOfHiddenLayers << std::endl
             << "Numbers of hidden neurons:" << std::endl
             << numbersOfHiddenNeurons << std::endl
             << "Number of outputs:" << std::endl
             << numberOfOutputs << std::endl;

   // Print hidden layer activation function

   std::cout << "Hidden layers activation function:" << std::endl;

   for(int i = 0; i < numberOfHiddenLayers; i++)
   {
      switch(hiddenLayersActivationFunction[i])
      {
         case Perceptron::Logistic:
         {
            std::cout << "Logistic ";
         }
         break;

         case Perceptron::HyperbolicTangent:
         {
            std::cout << "HyperbolicTangent ";
         }
         break;

         case Perceptron::Linear:
         {
            std::cout << "Linear ";
         }
         break;

         default:
         {
            std::cerr << std::endl 
                      << "Flood Error: MultilayerPerceptron class." << std::endl
                      << "void print(void) method." << std::endl
                      << "Unknown hidden layer activation function." << std::endl
                      << std::endl;
 
            exit(1);
         }
         break;
      }
   }

   std::cout << std::endl;

   // Print output layer activation function

   std::cout << "Output layer activation function:" << std::endl;

   switch(outputLayerActivationFunction)
   {
      case Perceptron::Logistic:
     {
        std::cout << "Logistic" << std::endl;
     }
     break;

     case Perceptron::HyperbolicTangent:
     {
        std::cout << "HyperbolicTangent" << std::endl;
     }
     break;

     case Perceptron::Linear:
     {
        std::cout << "Linear" << std::endl;
     }
     break;

     default:
     {
         std::cerr << std::endl 
                   << "Flood Error: MultilayerPerceptron class." << std::endl
                   << "void print(void) method." << std::endl
                   << "Unknown output layer activation function." << std::endl
                   << std::endl;
 
         exit(1);
     }
     break;
   }

   // Print name of input variables

   std::cout << "Name of input variables:" << std::endl
             << nameOfInputVariables << std::endl;

   // Print name of output variables

   std::cout << "Name of output variables:" << std::endl
             << nameOfOutputVariables << std::endl;

   // Print units of input variables

   std::cout << "Units of input variables:" << std::endl
             << unitsOfInputVariables << std::endl;

   // Print units of output variables

   std::cout << "Units of output variables:" << std::endl
             << unitsOfOutputVariables << std::endl;

   // Print description of input variables

   std::cout << "Description of input variables:" << std::endl
             << descriptionOfInputVariables << std::endl;

   // Print description of output variables

   std::cout << "Description of output variables:" << std::endl
             << descriptionOfOutputVariables << std::endl;

   // Print mean of input variables

   std::cout << "Mean of input variables:" << std::endl
             << meanOfInputVariables << std::endl;

   // Print standard deviation of input variables 

   std::cout << "Standard deviation of input variables:" << std::endl
             << standardDeviationOfInputVariables << std::endl;

   // Print mean of output variables

   std::cout << "Mean of output variables:" << std::endl
             << meanOfOutputVariables << std::endl;

   // Print standard deviation of output variables

   std::cout << "Standard deviation of output variables:" << std::endl
             << standardDeviationOfOutputVariables << std::endl;
 
   // Print minimum of input variables

   std::cout << "Minimum of input variables:" << std::endl
             << minimumOfInputVariables << std::endl;

   // Print maximum of input variables

   std::cout << "Maximum of input variables:" << std::endl
             << maximumOfInputVariables << std::endl;
   
   // Print minimum of output variables

   std::cout << "Minimum of output variables:" << std::endl
             << minimumOfOutputVariables << std::endl;

   // Print maximum of output variables

   std::cout << "Maximum of output variables:" << std::endl
             << maximumOfOutputVariables << std::endl;

   // Print lower bound of output variables

   std::cout << "Lower bound of output variables:" << std::endl
             << lowerBoundOfOutputVariables << std::endl;

   // Print upper bound of output variables

   std::cout << "Upper bound of output variables:" << std::endl
             << upperBoundOfOutputVariables << std::endl;

   // Print biases and synaptic weigths

   Vector<double> neuralParameters = getNeuralParameters();

   std::cout << "Neural parameters: " << std::endl 
             << neuralParameters << std::endl;

   // Print number of independent parameters

   std::cout << "Number of independent parameters: " << std::endl 
             << numberOfIndependentParameters << std::endl;
             
   // Print name of independent parameters

   std::cout << "Name of independent parameters:" << std::endl
             << nameOfIndependentParameters << std::endl;

   // Print units of independent parameters

   std::cout << "Units of independent parameters:" << std::endl
             << unitsOfIndependentParameters << std::endl;

   // Print description of independent parameters

   std::cout << "Description of independent parameters:" << std::endl
             << descriptionOfIndependentParameters << std::endl;

   // Print mean of independent parameters

   std::cout << "Mean of independent parameters:" << std::endl
             << meanOfIndependentParameters << std::endl;

   // Print standard deviation of independent parameters

   std::cout << "Standard deviation of independent parameters:" << std::endl
             << standardDeviationOfIndependentParameters << std::endl;

   // Print minimum of independent parameters

   std::cout << "Minimum of independent parameters:" << std::endl
             << minimumOfIndependentParameters << std::endl;

   // Print maximum of independent parameters

   std::cout << "Maximum of independent parameters:" << std::endl
             << maximumOfIndependentParameters << std::endl;

   // Print lower bound of independent parameters

   std::cout << "Lower bound of independent parameters:" << std::endl
             << lowerBoundOfIndependentParameters << std::endl;

   // Print upper bound of independent parameters

   std::cout << "Upper bound of independent parameters:" << std::endl
             << upperBoundOfIndependentParameters << std::endl;
 
   // Print independent parameters

   std::cout << "Independent parameters:" << std::endl 
             << independentParameters << std::endl;

   // Print pre and post processing method 

   std::cout << "Pre and post processing method:" << std::endl;

   if(preAndPostProcessingMethod == None)
   {
      std::cout << "None" << std::endl;
   }
   else if(preAndPostProcessingMethod == MeanAndStandardDeviation)
   {
      std::cout << "Mean and standard deviation" << std::endl;
   }
   else if(preAndPostProcessingMethod == MinimumAndMaximum)
   {
      std::cout << "Minimum and maximum" << std::endl;
   }
   else
   {
      std::cerr << std::endl 
                << "Flood Error: MultilayerPerceptron class." << std::endl
                << "void print(void) method." << std::endl
                << "Unknown pre and postprocessing method." << std::endl
                << std::endl;
 
      exit(1);
   }

   // Print numerical differentiation method 

   std::cout << "Numerical differentiation method:" << std::endl;

   if(numericalDifferentiationMethod == ForwardDifferences)
   {
      std::cout << "Forward differences" << std::endl;
   }
   else if(numericalDifferentiationMethod == CentralDifferences)
   {
      std::cout << "Central differences" << std::endl;
   }
   else
   {
      std::cerr << std::endl 
                << "Flood Error: MultilayerPerceptron class." << std::endl
                << "void print(void) method." << std::endl
                << "Unknown numerical differentiation method." << std::endl
                << std::endl;
 
      exit(1);
   }

   // Print display out of range warning

   std::cout << "Display out of range warning:" << std::endl
             << displayOutOfRangeWarning << std::endl;

   // Print display 

   std::cout << "Display:" << std::endl
             << display << std::endl;

   // Print numerical epsilon

   std::cout << "Numercial epsilon:" << std::endl
             << numericalEpsilon << std::endl;

}


// void save(char*) method 

/// This method saves to a data file the members of a multilayer perceptron object:
///
/// <ul>
/// <li> Number of inputs.
/// <li> Number of hidden neurons.
/// <li> Number of outputs.
/// <li> Hidden layer activation function.
/// <li> Output layer activation function.
/// <li> Name of input vaviables. 
/// <li> Name of output variables.
/// <li> Mean of input variables.
/// <li> Standard deviation of input variables.
/// <li> Mean of output variables.
/// <li> Standard deviation of output variables.
/// <li> Minimum of input variables.
/// <li> Maximum of input variables.
/// <li> Minimum of output variables.
/// <li> Maximum of output variables.
/// <li> Biases and synaptic weights.
/// <li> Pre and postprocessing method.
/// <li> Number of independent parameters.
/// <li> Name of independent parameters.
/// <li> Units of independent parameters.
/// <li> Description of independent parameters.
/// <li> Mean of independent parameters.
/// <li> Standard deviation of independent parameters.
/// <li> Minimum of independent parameters.
/// <li> Maximum of independent parameters.
/// <li> Independent parameters.
/// <li> Display out of range warning.
/// <li> Display.
/// <li> Numerical epsilon.
/// </ul> 
///
/// @param filename Filename.
///
/// @see load(char*).

void MultilayerPerceptron::save(char* filename)
{
   std::fstream file; 

   file.open(filename, std::ios::out);

   if(!file.is_open())
   {
      std::cerr << std::endl 
                << "Flood Error: MultilayerPerceptron class." << std::endl
                << "void save(char*) method." << std::endl
                << "Cannot open multilayer perceptron data file." << std::endl
                << std::endl;

      exit(1);
   }
   else
   {
      if(display)
      {
         std::cout << std::endl 
                   << "Saving multilayer perceptron to data file..." 
                   << std::endl;
      }
   }

   // Write file header

   file << "% Flood Neural Network. Multilayer Perceptron Object." << std::endl;

   // Write network architecture

   int numberOfHiddenLayers = getNumberOfHiddenLayers();

   file << "NumberOfInputs:" << std::endl
        << numberOfInputs << std::endl
        << "NumberOfHiddenLayers:" << std::endl
        << numberOfHiddenLayers << std::endl
        << "NumbersOfHiddenNeurons:" << std::endl
        << numbersOfHiddenNeurons << std::endl
        << "NumberOfOutputs:" << std::endl
        << numberOfOutputs << std::endl;

   // Write hidden layer activation function

   file << "HiddenLayersActivationFunction:" << std::endl;

   for(int i = 0; i < numberOfHiddenLayers; i++)
   {
      switch(hiddenLayersActivationFunction[i])   
      {
         case Perceptron::Logistic:
         {
            file << "Logistic ";
         }
         break;

         case Perceptron::HyperbolicTangent:
         {
            file << "HyperbolicTangent ";
         }
         break;

         case Perceptron::Linear:
         {
            file << "Linear ";
         }
         break;

         default:
         {
            std::cerr << std::endl 
                      << "Flood Error: MultilayerPerceptron class." << std::endl
                      << "void save(char*) method." << std::endl
                      << "Unknown hidden layer activation function." << std::endl
                      << std::endl;
    
            exit(1);
         }
         break;
      }
   }

   file << std::endl;

   // Write output layer activation function

   file << "OutputLayerActivationFunction:" << std::endl;

   switch(outputLayerActivationFunction)   
   {
      case Perceptron::Logistic:
      {
         file << "Logistic" << std::endl;
      }
      break;

      case Perceptron::HyperbolicTangent:
      {
         file << "HyperbolicTangent" << std::endl;
      }
      break;

      case Perceptron::Linear:
      {
         file << "Linear" << std::endl;
      }
      break;

      default:
      {
         std::cerr << std::endl 
                   << "Flood Error: MultilayerPerceptron class." << std::endl
                   << "void save(char*) method." << std::endl
                   << "Unknown output layer activation function." << std::endl
                   << std::endl;
 
         exit(1);
      }
      break;
   }

   // Write name of input variables

   file << "NameOfInputVariables:" << std::endl 
        << nameOfInputVariables << std::endl;

   // Write name of output variables

   file << "NameOfOutputVariables:" << std::endl
        << nameOfOutputVariables << std::endl;

   // Write units of input variables

   file << "UnitsOfInputVariables:" << std::endl 
        << unitsOfInputVariables << std::endl;

   // Write units of output variables

   file << "UnitsOfOutputVariables:" << std::endl
        << unitsOfOutputVariables << std::endl;

   // Write description of input variables

   file << "DescriptionOfInputVariables:" << std::endl 
        << descriptionOfInputVariables << std::endl;

   // Write description of output variables

   file << "DescriptionOfOutputVariables:" << std::endl
        << descriptionOfOutputVariables << std::endl;

   // Write mean of input variables

   file << "MeanOfInputVariables:" << std::endl
        << meanOfInputVariables << std::endl;

    // Write standard deviation of input variables

   file << "StandardDeviationOfInputVariables:" << std::endl
        << standardDeviationOfInputVariables << std::endl;
  
   // Write mean output variables

   file << "MeanOfOutputVariables:" << std::endl
        << meanOfOutputVariables << std::endl;

   // Write standard deviation of output variables

   file << "StandardDeviationOfOutputVariables:" << std::endl
        << standardDeviationOfOutputVariables << std::endl;

   // Write minimum of input variables

   file << "MinimumOfInputVariables:" << std::endl
        << minimumOfInputVariables << std::endl;

   // Write maximum of input variables

   file << "MaximumOfInputVariables:" << std::endl
        << maximumOfInputVariables << std::endl;

   // Write minimum of output variables

   file << "MinimumOfOutputVariables:" << std::endl 
        << minimumOfOutputVariables << std::endl;

   // Write maximum of output variables

   file << "MaximumOfOutputVariables:" << std::endl
        << maximumOfOutputVariables << std::endl;

   // Write lower bound of output variables

   file << "LowerBoundOfOutputVariables:" << std::endl 
        << lowerBoundOfOutputVariables << std::endl;

   // Write upper bound of output variables

   file << "UpperBoundOfOutputVariables:" << std::endl
        << upperBoundOfOutputVariables << std::endl;

   // Write biases and synaptic weights

   Vector<double> neuralParameters =  getNeuralParameters();

   file << "NeuralParameters:" << std::endl
        << neuralParameters << std::endl;

   // Write number of independent parameters

   file << "NumberOfIndependentParameters:" << std::endl
        << numberOfIndependentParameters << std::endl;

   // Write name of independent parameters

   file << "NameOfIndependentParameters:" << std::endl
        << nameOfIndependentParameters << std::endl;

   // Write units of independent parameters

   file << "UnitsOfIndependentParameters:" << std::endl
        << unitsOfIndependentParameters << std::endl;

   // Write description of independent parameters

   file << "DescriptionOfIndependentParameters:" << std::endl
        << descriptionOfIndependentParameters << std::endl;

   // Write mean of independent parameters

   file << "MeanOfIndependentParameters:" << std::endl
        << meanOfIndependentParameters << std::endl;

   // Write standard deviation of independent parameters

   file << "StandardDeviationOfIndependentParameters:" << std::endl
        << standardDeviationOfIndependentParameters << std::endl;

   // Write minimum of independent parameters

   file << "MinimumOfIndependentParameters:" << std::endl
        << minimumOfIndependentParameters << std::endl;

   // Write maximum of independent parameters

   file << "MaximumOfIndependentParameters:" << std::endl
        << maximumOfIndependentParameters << std::endl;

   // Write lower bound of independent parameters

   file << "LowerBoundOfIndependentParameters:" << std::endl
        << lowerBoundOfIndependentParameters << std::endl;

   // Write upper bound of independent parameters

   file << "UpperBoundOfIndependentParameters:" << std::endl
        << upperBoundOfIndependentParameters << std::endl;

   // Write independent parameters

   file << "IndependentParameters:" << std::endl
        << independentParameters << std::endl;

   // Write pre and post processing method 

   file << "PreAndPostProcessingMethod:" << std::endl;

   switch(preAndPostProcessingMethod)   
   {
      case None:
      {
         file << "None" << std::endl;
      }
      break;

      case MeanAndStandardDeviation:
      {
         file << "MeanAndStandardDeviation" << std::endl;
      }
      break;

      case MinimumAndMaximum:
      {
         file << "MinimumAndMaximum" << std::endl;
      }
      break;

      default:
      {
         std::cerr << std::endl 
                   << "Flood Error: MultilayerPerceptron class." << std::endl
                   << "void save(char*) method." << std::endl
				   << "Unknown pre and postprocesing method: " << preAndPostProcessingMethod << "." << std::endl
                   << std::endl;
 
         exit(1);
      }
      break;

   }// end switch

   // Write numerical differentiation method 

   file << "NumericalDifferentiationMethod:" << std::endl;

   switch(numericalDifferentiationMethod)   
   {
      case ForwardDifferences:
      {
         file << "ForwardDifferences" << std::endl;
      }
      break;

      case CentralDifferences:
      {
         file << "CentralDifferences" << std::endl;
      }
      break;

      default:
      {
         std::cerr << std::endl 
                   << "Flood Error: MultilayerPerceptron class." << std::endl
                   << "void save(char*) method." << std::endl
                   << "Unknown numerical differentiation method." << std::endl
                   << std::endl;
 
         exit(1);
      }
      break;

   }

   // Write display out of range warning

   file << "DisplayOutOfRangeWarning:" << std::endl
        << displayOutOfRangeWarning << std::endl;

   // Write display

   file << "Display:" << std::endl
        << display << std::endl;

   // Write numerical epsilon

   file << "NumericalEpsilon:" << std::endl
        << numericalEpsilon << std::endl;

   // Close file

   file.close();
}


// void load(char*) method

/// This method loads from a data file the members of a multilayer perceptron object: 
///
/// <ul>
/// <li> Number of inputs.
/// <li> Number of hidden neurons.
/// <li> Number of outputs.
/// <li> Hidden layer activation function.
/// <li> Output layer activation function.
/// <li> Name of input vaviables. 
/// <li> Name of output variables.
/// <li> Units of input vaviables. 
/// <li> Units of output variables.
/// <li> Description of input vaviables. 
/// <li> Description of output variables.
/// <li> Mean of input variables.
/// <li> Standard deviation of input variables.
/// <li> Mean of output variables.
/// <li> Standard deviation of output variables.
/// <li> Minimum of input variables.
/// <li> Maximum of input variables.
/// <li> Minimum of output variables.
/// <li> Maximum of output variables.
/// <li> Lower bound of output variables.
/// <li> Upper bound of output variables.
/// <li> Biases and synaptic weights.
/// <li> Number of independent parameters.
/// <li> Name of independent parameters.
/// <li> Mean of independent parameters.
/// <li> Standard deviation of independent parameters.
/// <li> Minimum of independent parameters.
/// <li> Maximum of independent parameters.
/// <li> Lower bound of independent parameters.
/// <li> Upper bound of independent parameters.
/// <li> Independent parameters.
/// <li> Pre and postprocessing method.
/// <li> Numerical epsilon method.
/// <li> Numerical differentiation method.
/// <li> Display out of range warning.
/// <li> Display.
/// <li> Numerical epsilon.
/// </ul> 
///
/// Please mind about the file format, which is specified in the User's Guide. 
///
/// @param filename Filename.
///
/// @see save(char*).

void MultilayerPerceptron::load(char* filename)
{

   std::fstream file;

   file.open(filename, std::ios::in);

   if(!file.is_open())
   {
      std::cerr << std::endl
                << "Flood Error: MultilayerPerceptron class." << std::endl
                << "void load(char*) method." << std::endl
                << "Cannot open multilayer perceptron data file." << std::endl
                << std::endl;

      exit(1);
   }
   else
   {
      if(display)
      {
         std::cout << std::endl
                   << "Loading multilayer perceptron from data file..." << std::endl;
      }
   }

   std::string word;

   // Number of inputs

   while(word != "NumberOfInputs:")
   {
      file >> word;
   }

   file >> numberOfInputs;

   // Number of hidden layers 

   file >> word;

   int numberOfHiddenLayers;

   file >> numberOfHiddenLayers;

   numbersOfHiddenNeurons.resize(numberOfHiddenLayers);

   // Numbers of hidden neurons

   file >> word;

   file >> numbersOfHiddenNeurons;

   // Number of outputs 

   file >> word;

   file >> numberOfOutputs;

   // Set architecture

   setNetworkArchitecture(numberOfInputs, numbersOfHiddenNeurons, numberOfOutputs);

   // Hidden layers activation function

   file >> word;
 
   for(int i = 0; i < numberOfHiddenLayers; i++)
   {
      file >> word;

      if(word == "Logistic")
      {
         hiddenLayersActivationFunction[i] = Perceptron::Logistic;
      }
      else if(word == "HyperbolicTangent")
      {
         hiddenLayersActivationFunction[i] = Perceptron::HyperbolicTangent;
      }
      else if(word == "Linear")
      {
         hiddenLayersActivationFunction[i] = Perceptron::Linear;
      }
      else
      {
         std::cerr << std::endl
                   << "Flood Error: MultilayerPerceptron class." << std::endl
                   << "void load(char*) method." << std::endl
                   << "Unknown file format. " << word << "???" <<std::endl
                   << std::endl;
   
         exit(1);   
      }
   }   
    
   setHiddenLayersActivationFunction(hiddenLayersActivationFunction);

   // Output layer activation function

   file >> word;

   file >> word;

   if(word == "Logistic")
   {
      outputLayerActivationFunction = Perceptron::Logistic;
   }
   else if(word == "HyperbolicTangent")
   {
      outputLayerActivationFunction = Perceptron::HyperbolicTangent;
   }
   else if(word == "Linear")
   {
      outputLayerActivationFunction = Perceptron::Linear;
   }
   else
   {
      std::cerr << std::endl
                << "Flood Error: MultilayerPerceptron class." << std::endl
                << "void load(char*) method." << std::endl
                << "Unknown file format. " << word << "???" <<std::endl
                << std::endl;

      exit(1);   
   }

   setOutputLayerActivationFunction(outputLayerActivationFunction);

   // Name of input variables

   file >> word;

   file >> nameOfInputVariables;

   // Name of output variables

   file >> word;

   file >> nameOfOutputVariables;

   // Units of input variables

   file >> word;

   file >> unitsOfInputVariables;

   // Units of output variables

   file >> word;

   file >> unitsOfOutputVariables;

   // Description of input variables

   file >> word;

   file >> descriptionOfInputVariables;

   // Description of output variables

   file >> word;

   file >> descriptionOfOutputVariables;

   // Mean of input variables

   file >> word;

   file >> meanOfInputVariables;

   // Standard deviation of input variables

   file >> word;

   file >> standardDeviationOfInputVariables;

   // Mean of output variables

   file >> word;

   file >> meanOfOutputVariables;

   // Standard deviation of output variables

   file >> word;

   file >> standardDeviationOfOutputVariables;

   // Minimum of input variables

   file >> word;

   file >> minimumOfInputVariables;

    // Maximum of input variables

   file >> word;

   file >> maximumOfInputVariables;

   // Minimum of output variables

   file >> word;

   file >> minimumOfOutputVariables;

    // Maximum of output variables

   file >> word;

   file >> maximumOfOutputVariables;

   // Lower bound of output variables

   file >> word;

   file >> lowerBoundOfOutputVariables;

    // Upper bound of output variables

   file >> word;

   file >> upperBoundOfOutputVariables;

   // Neural parameters

   file >> word;

   int numberOfNeuralParameters = getNumberOfNeuralParameters();

   Vector<double> neuralParameters(numberOfNeuralParameters);

   file >> neuralParameters;

   setNeuralParameters(neuralParameters);

   // Number of independent parameters

   file >> word;
   
   file >> numberOfIndependentParameters;

   // Set number of independent parameters
   
   setNumberOfIndependentParameters(numberOfIndependentParameters);

   // Name of independent parameters

   file >> word; 

   file >> nameOfIndependentParameters;

   // Units of independent parameters

   file >> word; 

   file >> unitsOfIndependentParameters;

   // Description of independent parameters

   file >> word; 

   file >> descriptionOfIndependentParameters;

   // Mean of independent parameters

   file >> word; 

   file >> meanOfIndependentParameters;

   // Standard deviation of independent parameters

   file >> word; 

   file >> standardDeviationOfIndependentParameters;

   // Minimum of independent parameters

   file >> word; 

   file >> minimumOfIndependentParameters;

   // Maximum of independent parameters

   file >> word; 

   file >> maximumOfIndependentParameters;

   // Lower bound of independent parameters

   file >> word; 

   file >> lowerBoundOfIndependentParameters;

   // Upper bound of independent parameters

   file >> word; 

   file >> upperBoundOfIndependentParameters;

   // Independent parameters

   file >> word; 

   file >> independentParameters;

   // Pre and postprocessing method

   file >> word;

   file >> word;

   if(word == "None")
   {
      preAndPostProcessingMethod = None;
   }
   else if(word == "MeanAndStandardDeviation")
   {
      preAndPostProcessingMethod = MeanAndStandardDeviation;
   }
   else if(word == "MinimumAndMaximum")
   {
      preAndPostProcessingMethod = MinimumAndMaximum;
   }
   else
   {
      std::cerr << std::endl
                << "Flood Error: MultilayerPerceptron class." << std::endl
                << "void load(char*) method." << std::endl
                << "Unknown file format. " << word << "???" <<std::endl
                << std::endl;

      exit(1);   
   }

   // Numerical differentiation method

   file >> word;

   file >> word;

   if(word == "ForwardDifferences")
   {
      numericalDifferentiationMethod = ForwardDifferences;
   }
   else if(word == "CentralDifferences")
   {
      numericalDifferentiationMethod = CentralDifferences;
   }
   else
   {
      std::cerr << std::endl
                << "Flood Error: MultilayerPerceptron class." << std::endl
                << "void load(char*) method." << std::endl
                << "Unknown file format. " << word << "???" <<std::endl
                << std::endl;

      exit(1);   
   }

   // Display out of range warning

   file >> word;
   
   file >> displayOutOfRangeWarning;

   // Display

   file >> word;
   
   file >> display;

   // Numerical epsilon

   file >> word;
   
   file >> numericalEpsilon;

   // Close file

   file.close();
}


// void saveExpression(char*) method

/// This method saves the explicit mathematical expression addressed by 
/// the multilayer perceptron to a file. 
///
/// @param filename Filename.
///
/// @todo This method is not fully implemented. 

void MultilayerPerceptron::saveExpression(char* filename)
{
/*
   std::fstream file;

   file.open(filename, std::ios::out);

   if(!file.is_open())
   {
      std::cerr << std::endl 
                << "Flood Error: MultilayerPerceptron class." << std::endl
                << "void saveExpression(char*) method." << std::endl
                << "Cannot open expression data file."  << std::endl
                << std::endl;

      exit(1);
   }
   else
   {
      if(display)
      {
         std::cout << std::endl
                   << "Saving expression to data file..." << std::endl;
      }
   }

   // Write file header

   file << "% Flood Neural Network. Multilayer perceptron expression." << std::endl
        << std::endl;


   file.precision(6); 

   // Preprocess inputs

   switch(preAndPostProcessingMethod)
   {
      case MultilayerPerceptron::None:
      {
         // Do nothing
      }// end none
      break;
   
      case MultilayerPerceptron::MeanAndStandardDeviation:
      {
         // Preprocess inputs

         for(int i = 0; i < numberOfInputs; i++)
         {
            file << nameOfInputVariables[i] << "=(" << nameOfInputVariables[i] << "-" << meanOfInputVariables[i] << ")/" 
                 << standardDeviationOfInputVariables[i] << ";" << std::endl;
         }
      }// end mean and standard deviation
      break;
   
      case MultilayerPerceptron::MinimumAndMaximum:
      {
         // Preprocess inputs

         for(int i = 0; i < numberOfInputs; i++)
         {
            file << nameOfInputVariables[i] << "=2*(" << nameOfInputVariables[i] << "-" << minimumOfInputVariables[i] 
                 << ")/(" << maximumOfInputVariables[i] << "-" << minimumOfInputVariables[i] << ")-1;" << std::endl;
         }
      }// end minimum and maximum
      break;

      default:
      {
         std::cerr << std::endl 
                   << "Flood Error: MultilayerPerceptron class." << std::endl
                   << "void saveExpression(char*) method." << std::endl
                   << "Unknown pre and postprocessing method." << std::endl
                   << std::endl;
 
         exit(1);
      }// end default
      break;

   }

   int numberOfHiddenLayers = getNumberOfHiddenLayers();

   // Write expression

   for(int output = 0; output < numberOfOutputs; output++)
   {
      file << nameOfOutputVariables[output] << "=" << outputLayer[output].getBias();

      for(int hiddenLayer = 0; hiddenLayer < numberOfHiddenLayers; hiddenLayer++) 
      {
         for(int hiddenNeuron = 0; hiddenNeuron < numbersOfHiddenNeurons[hiddenLayer]; hiddenNeuron++)
         {
//            file << "+" << outputLayer[output].getSingleSynapticWeight(j);
       
            switch(hiddenLayersActivationFunction[hiddenLayer])
            {
               case Perceptron::Logistic:
               {
                  file << "*logistic(";
               }
               break;

               case Perceptron::HyperbolicTangent:
               {
                  file << "*tanh(";
               }          
               break;

               case Perceptron::Linear:
               {
                  file << "*(";
               }
               break;
            }

//            file << hiddenLayer[j].getBias();

            for(int input = 0; input < numberOfInputs; input++)
            {
//               file << "+" << hiddenLayer[j].getSingleSynapticWeight(input) << "*" << nameOfInputVariables[input];
            }

            file << ")";
         }       

         file << ";" << std::endl;
      }
   }

   // Postprocess outputs

   switch(preAndPostProcessingMethod)   
   {
      case MultilayerPerceptron::None:
      {
         // Do nothing
      }
      break;
  
      case MultilayerPerceptron::MeanAndStandardDeviation:
      {
         // Postprocess outputs

         for(int i = 0; i < numberOfOutputs; i++)
         {
            file << nameOfOutputVariables[i] << "=" <<  meanOfOutputVariables[i] << "+" 
                 << standardDeviationOfOutputVariables[i] << "*" << nameOfOutputVariables[i] << ";" << std::endl;
         }
      }
      break;

      case MultilayerPerceptron::MinimumAndMaximum:
      {   
         // Postprocess outputs

         for(int i = 0; i < numberOfOutputs; i++)
         {
            file << nameOfOutputVariables[i] << "=0.5*(" << nameOfOutputVariables[i] << "+1.0)*(" 
                 << maximumOfOutputVariables[i] << "-" << minimumOfOutputVariables[i] << ")+" 
                 << minimumOfOutputVariables[i] << ";" << std::endl;       
         }
      }
      break;
   }

   file.close();
*/
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
