/****************************************************************************************************************/
/*                                                                                                              */
/*   Flood: An Open Source Neural Networks C++ Library                                                          */
/*   www.cimne.com/flood                                                                                        */
/*                                                                                                              */
/*   I N P U T   T A R G E T   D A T A   S E T   C L A S S                                                      */
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

#include "InputTargetDataSet.h"

namespace Flood
{

// GENERAL CONSTRUCTOR

/// General constructor. It creates an input-target data set object involving an input
/// Matrix and an Target matrix. All input-target data is initialized to zero.
/// It also initializes the rest of class members to their default values:
///
/// <ul>
/// <li> Name of input vaviables: InputVariable0,... 
/// <li> Name of target variables: TargetVariable0,...
/// <li> Units of input vaviables: None. 
/// <li> Units of target variables: None.
/// <li> Description of input vaviables: None. 
/// <li> Descritpion of target variables: None.
/// <li> Display messages to screen: True.
/// </ul> 
///
/// @param newNumberOfSamples Number of samples in the data set.
/// @param newNumberOfInputVariables Number of input variables.
/// @param newNumberOfTargetVariables Number of target variables.

InputTargetDataSet::InputTargetDataSet
(int newNumberOfSamples, int newNumberOfInputVariables, int newNumberOfTargetVariables)
{
   // Number of samples 

   numberOfSamples = newNumberOfSamples;

   // Number of input and target variables

   numberOfInputVariables = newNumberOfInputVariables;
   numberOfTargetVariables = newNumberOfTargetVariables;

   // Name of input variables

   nameOfInputVariables.resize(numberOfInputVariables);

   for(int i = 0; i < numberOfInputVariables; i++)
   {
      std::stringstream buffer;
  
      buffer << "InputVariable" << i;

      nameOfInputVariables[i] = buffer.str();
   }

   // Name of target variables

   nameOfTargetVariables.resize(numberOfTargetVariables);

   for(int i = 0; i < numberOfTargetVariables; i++)
   {
      std::stringstream buffer;
  
      buffer << "TargetVariable" << i;

      nameOfTargetVariables[i] = buffer.str();
   }

   // Units of input variables

   unitsOfInputVariables.resize(numberOfInputVariables);

   for(int i = 0; i < numberOfInputVariables; i++)
   {
      unitsOfInputVariables[i] = "None";
   }

   // Units of target variables

   unitsOfTargetVariables.resize(numberOfTargetVariables);

   for(int i = 0; i < numberOfTargetVariables; i++)
   {
      unitsOfTargetVariables[i] = "None";
   }

   // Description of input variables

   descriptionOfInputVariables.resize(numberOfInputVariables);

   for(int i = 0; i < numberOfInputVariables; i++)
   {
      descriptionOfInputVariables[i] = "None";
   }

   // Description of target variables

   descriptionOfTargetVariables.resize(numberOfTargetVariables);

   for(int i = 0; i < numberOfTargetVariables; i++)
   {
      descriptionOfTargetVariables[i] = "None";
   }

   // Input data 

   inputData.resize(numberOfSamples, numberOfInputVariables);

   // Target data 

   targetData.resize(numberOfSamples, numberOfTargetVariables);
    
   // Display

   display = true;
}


// DEFAULT CONSTRUCTOR

/// Default constructor. It creates an input-target data set object with zero
/// samples and zero input and target variables. 
/// It also initializes the rest of class members to their default values:
///
/// <ul>
/// <li> Display messages to screen: True.
/// </ul> 

InputTargetDataSet::InputTargetDataSet(void)
{
   numberOfSamples = 0;
   numberOfInputVariables = 0;
   numberOfTargetVariables = 0;
   
   display = true;
}


// DESTRUCTOR

/// Destructor. 

InputTargetDataSet::~InputTargetDataSet()
{

}


// METHODS


// int getNumberOfSamples(void) method

/// This method returns the number of samples in the input-target data set.

int InputTargetDataSet::getNumberOfSamples(void)
{
   return(numberOfSamples);   
}


// int getNumberOfInputVariables(void) method

/// This method returns the number of input variables of the input-target data set.

int InputTargetDataSet::getNumberOfInputVariables(void)
{
   return(numberOfInputVariables);   
}


// int getNumberOfTargetVariables(void) method

/// This method returns the number of target variables of the input-target data set.

int InputTargetDataSet::getNumberOfTargetVariables(void)
{
   return(numberOfTargetVariables);   
}


// Matrix<double>& getInputData(void) method

/// This method returns a Matrix containing the input data of the data set.

Matrix<double>& InputTargetDataSet::getInputData(void)
{
   return(inputData);   
}


// Matrix<double>& getTargetData(void) method

/// This method returns a Matrix containing the target data of the data set.

Matrix<double>& InputTargetDataSet::getTargetData(void)
{
   return(targetData);   
}


// Matrix<double> getInputAndTargetData(void) method

/// This method returns a Matrix containing both the input and the target data 
/// in the data set.

Matrix<double> InputTargetDataSet::getInputAndTargetData(void)
{
   int numberOfVariables = numberOfInputVariables + numberOfTargetVariables;

   Matrix<double> inputAndTargetData(numberOfSamples, numberOfVariables, 0.0);

   for(int i = 0; i < numberOfSamples; i++)
   {
      for(int j = 0; j < numberOfInputVariables; j++)
      {
        inputAndTargetData[i][j] = inputData[i][j];
      }

      for(int j = 0; j < numberOfTargetVariables; j++)
      {
        inputAndTargetData[i][numberOfInputVariables+j] = targetData[i][j];
      }
   }

   return(inputAndTargetData);   
}


// Vector<double> getSample(int) method

/// This method returns the input and target values of a single sample in the 
/// input-target data set. 
///
/// @param i Index of the sample. 

Vector<double> InputTargetDataSet::getSample(int i)
{
   // Control sentence (if debug)

   #ifndef NDEBUG 

   if(i >= numberOfSamples)
   {
      std::cerr << std::endl
                << "Flood Error: InputTargetDataSet class." << std::endl
                << "Vector<double> getSample(int) method." << std::endl
                << "Index of sample must be less than number of samples." << std::endl
                << std::endl;

      exit(1);   
   }

   #endif

   // Get sample

   Vector<double> sample(numberOfInputVariables+numberOfTargetVariables, 0.0);

   // Get input data

   for(int j = 0; j < numberOfInputVariables; j++)
   {
      sample[j] = inputData[i][j];
   }

   // Get target data

   for(int j = 0; j < numberOfTargetVariables; j++)
   {
      sample[numberOfInputVariables+j] = targetData[i][j];
   }

   return(sample);
}


// void setSample(int, Vector<double>)

/// This method sets new input and target values of a single sample in the input-target data set. 
///
/// @param i Index of the sample. 
/// @param sample New input and target values of the sample.

void InputTargetDataSet::setSample(int i, Vector<double> sample)
{
   // Control sentence (if debug)

   #ifndef NDEBUG 

   int size = sample.getSize();

   if(i >= numberOfSamples)
   {
      std::cerr << std::endl
                << "Flood Error: InputTargetDataSet class." << std::endl
                << "void setSample(int, Vector<double>) method." << std::endl
                << "Index of sample must be less than number of samples." << std::endl
                << std::endl;

      exit(1);   
   }
   else if(size != numberOfInputVariables+numberOfTargetVariables)
   {
      std::cerr << std::endl
                << "Flood Error: InputTargetDataSet class." << std::endl
                << "void setSample(int, Vector<double>) method." << std::endl
                << "Size of sample must be equal to number of input variables plus number of target variables." 
				<< std::endl
                << std::endl;

      exit(1);   
   }

   #endif

   // Set input data

   for(int j = 0; j < numberOfInputVariables; j++)
   {
      inputData[i][j] = sample[j];
   }

   // Set target data

   for(int j = 0; j < numberOfTargetVariables; j++)
   {
      targetData[i][j] = sample[numberOfInputVariables+j];
   }
}


// void addSample(Vector<double>) method

/// This method adds a new input-target sample to the input-target data set. 
/// Note that resizing is here necessary and therefore computationally expensive. 
///
/// @param sample Input and target values of the sample to be added. 
///
/// @see subtractSample(int).

void InputTargetDataSet::addSample(Vector<double> sample)
{
   // Control sentence (if debug)

   #ifndef NDEBUG 

   int size = sample.getSize();

   if(size != numberOfInputVariables+numberOfTargetVariables)
   {
      std::cerr << std::endl
                << "Flood Error: InputTargetDataSet class." << std::endl
                << "void addSample(Vector<double>) method." << std::endl
                << "Size of sample must be equal to number of input variables plus number of target variables." 
				<< std::endl
                << std::endl;

      exit(1);   
   }

   #endif

   // Add input values from sample

   Matrix<double> newInputData(numberOfSamples+1, numberOfInputVariables, 0.0);

   for(int i = 0; i < numberOfSamples; i++)
   {
      for(int j = 0; j < numberOfInputVariables; j++)
      {
         newInputData[i][j] = inputData[i][j];
      }
   }

   for(int j = 0; j < numberOfInputVariables; j++)
   {
      newInputData[numberOfSamples][j] = sample[j];
   }

   inputData = newInputData;

   // Add target values from sample

   Matrix<double> newTargetData(numberOfSamples+1, numberOfTargetVariables, 0.0);

   for(int i = 0; i < numberOfSamples; i++)
   {
      for(int j = 0; j < numberOfTargetVariables; j++)
      {
        newTargetData[i][j] = targetData[i][j];
      }
   }

   for(int j = 0; j < numberOfTargetVariables; j++)
   {
      newTargetData[numberOfSamples][j] = sample[numberOfInputVariables+j];
   }

   targetData = newTargetData;

   // Set new number of samples

   numberOfSamples++;
}


// void subtractSample(int) method

/// This method substract the input-target sample with a given index from the input-target data set.
/// Note that resizing is here necessary and therefore computationally expensive. 
///
/// @param sampleIndex Index of sample to be removed. 
///
/// @see addSample(Vector<double>).

void InputTargetDataSet::subtractSample(int sampleIndex)
{
   // Control sentence (if debug)

   #ifndef NDEBUG 

   if(sampleIndex >= numberOfSamples)
   {
      std::cerr << std::endl
                << "Flood Error: InputTargetDataSet class." << std::endl
                << "void subtractSample(int) method." << std::endl
                << "Index of sample must be less than number of samples." << std::endl
                << std::endl;

      exit(1);   
   }

   #endif

   // Subtract input values from sample

   Matrix<double> newInputData(numberOfSamples-1, numberOfInputVariables, 0.0);

   for(int i = 0; i < sampleIndex; i++)
   {
      for(int j = 0; j < numberOfInputVariables; j++)
      {
        newInputData[i][j] = inputData[i][j];
      }
   }

   for(int i = sampleIndex+1; i < numberOfSamples; i++)
   {
      for(int j = 0; j < numberOfInputVariables; j++)
      {
        newInputData[i-1][j] = inputData[i][j];
      }
   }

   inputData = newInputData;

   // Subtract target values from sample

   Matrix<double> newTargetData(numberOfSamples-1, numberOfTargetVariables, 0.0);

   for(int i = 0; i < sampleIndex; i++)
   {
      for(int j = 0; j < numberOfTargetVariables; j++)
      {
        newTargetData[i][j] = targetData[i][j];
      }
   }

   for(int i = sampleIndex+1; i < numberOfSamples; i++)
   {
      for(int j = 0; j < numberOfTargetVariables; j++)
      {
        newTargetData[i-1][j] = targetData[i][j];
      }
   }

   targetData = newTargetData;

   // Set new number of samples

   numberOfSamples--;
}


// Vector<std::string> getNameOfInputVariables(void) method

/// This method returns the name of the input variables in the input-target data set.
/// Such names are only used to give the user basic information about the problem at hand.
///
/// @see getNameOfTargetVariables(void).
/// @see getNameOfSingleInputVariable(void).
/// @see getNameOfSingleTargetVariable(void).

Vector<std::string> InputTargetDataSet::getNameOfInputVariables(void)
{
   return(nameOfInputVariables);
}


// Vector<std::string> getNameOfTargetVariables(void) method

/// This method returns the name of the target variables in the input-target data set.
/// Such names are only used to give the user basic information about the problem at hand.
///
/// @see getNameOfInputVariables(void).
/// @see getNameOfSingleInputVariable(void).
/// @see getNameOfSingleTargetVariable(void).

Vector<std::string> InputTargetDataSet::getNameOfTargetVariables(void)
{
   return(nameOfTargetVariables);
}


// std::string getNameOfSingleInputVariable(int) method

/// This method returns the name of a single input variable in the input-target data set.
/// Such a name is only used to give the user basic information about the problem at hand.
///
/// @param i Index of input variable.
///
/// @see getNameOfInputVariables(void).
/// @see getNameOfTargetVariables(void).
/// @see getNameOfSingleTargetVariable(void).

std::string InputTargetDataSet::getNameOfSingleInputVariable(int i)
{
   // Control sentence (if debug)

   #ifndef NDEBUG 

   if(i >= numberOfInputVariables)
   {
      std::cerr << std::endl
                << "Flood Error: InputTargetDataSet class." << std::endl
                << "std::string getNameOfSingleInputVariable(int) method." << std::endl
                << "Index must be less than number of input variables." << std::endl
                << std::endl;

      exit(1);   
   }

   #endif

   // Return name of input variable

   return(nameOfInputVariables[i]);
}


// std::string getNameOfSingleTargetVariable(int) method

/// This method returns the name of a single target variable in the input-target data set.
/// Such a name is only used to give the user basic information about the problem at hand.
///
/// @param i Index of input variable.
///
/// @see getNameOfInputVariables(void).
/// @see getNameOfTargetVariables(void).
/// @see getNameOfSingleInputVariable(void).

std::string InputTargetDataSet::getNameOfSingleTargetVariable(int i)
{
   // Control sentence (if debug)

   #ifndef NDEBUG 

   if(i >= numberOfTargetVariables)
   {
      std::cerr << std::endl
                << "Flood Error: InputTargetDataSet class." << std::endl
                << "std::string getNameOfSingleTargetVariable(int) method." << std::endl
                << "Index must be less than number of target variables." << std::endl
                << std::endl;

      exit(1);   
   }

   #endif

   // Return name of target variable

   return(nameOfTargetVariables[i]);
}


// Vector<std::string> getUnitsOfInputVariables(void) method

/// This method returns the units of the input variables in the input-target data set.
/// Such names are only used to give the user basic information about the problem at hand.
///
/// @see getUnitsOfTargetVariables(void).
/// @see getUnitsOfSingleInputVariable(void).
/// @see getUnitsOfSingleTargetVariable(void).

Vector<std::string> InputTargetDataSet::getUnitsOfInputVariables(void)
{
   return(unitsOfInputVariables);
}


// Vector<std::string> getUnitsOfTargetVariables(void) method

/// This method returns the units of the target variables in the input-target data set.
/// Such names are only used to give the user basic information about the problem at hand.
///
/// @see getUnitsOfInputVariables(void).
/// @see getUnitsOfSingleInputVariable(void).
/// @see getUnitsOfSingleTargetVariable(void).

Vector<std::string> InputTargetDataSet::getUnitsOfTargetVariables(void)
{
   return(unitsOfTargetVariables);
}


// std::string getUnitsOfSingleInputVariable(int) method

/// This method returns the units of a single input variable in the input-target data set.
/// Such units are only used to give the user basic information about the problem at hand.
///
/// @param i Index of input variable.
///
/// @see getUnitsOfInputVariables(void).
/// @see getUnitsOfTargetVariables(void).
/// @see getUnitsOfSingleTargetVariable(void).

std::string InputTargetDataSet::getUnitsOfSingleInputVariable(int i)
{
   // Control sentence (if debug)

   #ifndef NDEBUG 

   if(i >= numberOfInputVariables)
   {
      std::cerr << std::endl
                << "Flood Error: InputTargetDataSet class." << std::endl
                << "std::string getUnitsOfSingleInputVariable(int) method." << std::endl
                << "Index must be less than number of input variables." << std::endl
                << std::endl;

      exit(1);   
   }

   #endif

   // Return units of input variable

   return(unitsOfInputVariables[i]);
}


// std::string getUnitsOfSingleTargetVariable(int) method

/// This method returns the units of a single target variable in the input-target data set.
/// Such units are only used to give the user basic information about the problem at hand.
///
/// @param i Index of input variable.
///
/// @see getUnitsOfInputVariables(void).
/// @see getUnitsOfTargetVariables(void).
/// @see getUnitsOfSingleInputVariable(void).

std::string InputTargetDataSet::getUnitsOfSingleTargetVariable(int i)
{
   // Control sentence (if debug)

   #ifndef NDEBUG 

   if(i >= numberOfTargetVariables)
   {
      std::cerr << std::endl
                << "Flood Error: InputTargetDataSet class." << std::endl
                << "std::string getUnitsOfSingleTargetVariable(int) method." << std::endl
                << "Index must be less than number of target variables." << std::endl
                << std::endl;

      exit(1);   
   }

   #endif

   // Return units of target variable

   return(unitsOfTargetVariables[i]);
}



// Vector<std::string> getDescriptionOfInputVariables(void) method

/// This method returns the current description of the input variables in the input-target data set.
/// Such names are only used to give the user basic information about the problem at hand.
///
/// @see getDescriptionOfTargetVariables(void).
/// @see getDescriptionOfSingleInputVariable(void).
/// @see getDescriptionOfSingleTargetVariable(void).

Vector<std::string> InputTargetDataSet::getDescriptionOfInputVariables(void)
{
   return(descriptionOfInputVariables);
}


// Vector<std::string> getDescriptionOfTargetVariables(void) method

/// This method returns the current description of the target variables in the input-target data set.
/// Such names are only used to give the user basic information about the problem at hand.
///
/// @see getDescriptionOfInputVariables(void).
/// @see getDescriptionOfSingleInputVariable(void).
/// @see getDescriptionOfSingleTargetVariable(void).

Vector<std::string> InputTargetDataSet::getDescriptionOfTargetVariables(void)
{
   return(descriptionOfTargetVariables);
}


// std::string getDescriptionOfSingleInputVariable(int) method

/// This method returns the description of a single input variable in the input-target data set.
/// Such a description is only used to give the user basic information about the problem at hand.
///
/// @param i Index of input variable.
///
/// @see getDescriptionOfInputVariables(void).
/// @see getDescriptionOfTargetVariables(void).
/// @see getDescriptionOfSingleTargetVariable(void).

std::string InputTargetDataSet::getDescriptionOfSingleInputVariable(int i)
{
   // Control sentence (if debug)

   #ifndef NDEBUG 

   if(i >= numberOfInputVariables)
   {
      std::cerr << std::endl
                << "Flood Error: InputTargetDataSet class." << std::endl
                << "std::string getDescriptionOfSingleInputVariable(int) method." << std::endl
                << "Index must be less than number of input variables." << std::endl
                << std::endl;

      exit(1);   
   }

   #endif

   // Return description of input variable

   return(descriptionOfInputVariables[i]);
}


// std::string getDescriptionOfSingleTargetVariable(int) method

/// This method returns the description of a single target variable in the input-target data set.
/// Such a description is only used to give the user basic information about the problem at hand.
///
/// @param i Index of input variable.
///
/// @see getDescriptionOfInputVariables(void).
/// @see getDescriptionOfTargetVariables(void).
/// @see getDescriptionOfSingleInputVariable(void).

std::string InputTargetDataSet::getDescriptionOfSingleTargetVariable(int i)
{
   // Control sentence (if debug)

   #ifndef NDEBUG 

   if(i >= numberOfTargetVariables)
   {
      std::cerr << std::endl
                << "Flood Error: InputTargetDataSet class." << std::endl
                << "std::string getUnitsOfSingleTargetVariable(int) method." << std::endl
                << "Index must be less than number of target variables." << std::endl
                << std::endl;

      exit(1);   
   }

   #endif

   // Return description of target variable

   return(descriptionOfTargetVariables[i]);
}


// Vector< Vector<std::string> > getAllInformation(void) method

Vector< Vector<std::string> > InputTargetDataSet::getAllInformation(void)
{
   Vector< Vector<std::string> > information(6);

   information[0] = nameOfInputVariables;
   information[1] = nameOfTargetVariables;

   information[2] = unitsOfInputVariables;
   information[3] = unitsOfTargetVariables;

   information[4] = descriptionOfInputVariables;
   information[5] = descriptionOfTargetVariables;

   return(information);
}


// bool getDisplay(void) method

/// This method returns true if messages from this class can be displayed on the screen,
/// or false if messages from this class can't be displayed on the screen.

bool InputTargetDataSet::getDisplay(void)
{
   return(display);   
}


// void setInputTargetDataSet(int, int, int) method

/// This method sets new numbers of samples and input and target variables in the input-target data set.
/// It also initializes the rest of class members to their default values:
///
/// <ul>
/// <li> Name of input vaviables: InputVariable0,... 
/// <li> Name of target variables: TargetVariable0,...
/// <li> Units of input vaviables: None. 
/// <li> Units of target variables: None.
/// <li> Description of input vaviables: None. 
/// <li> Description of target variables: None.
/// </ul> 
///
/// @param newNumberOfSamples Number of samples.
/// @param newNumberOfInputVariables Number of input variables.
/// @param newNumberOfTargetVariables Number of target variables.

void InputTargetDataSet::setInputTargetDataSet
(int newNumberOfSamples, int newNumberOfInputVariables, int newNumberOfTargetVariables)
{
   // Number of samples

   numberOfSamples = newNumberOfSamples;

   // Number of input and target variables

   numberOfInputVariables = newNumberOfInputVariables;   
   numberOfTargetVariables = newNumberOfTargetVariables;   

   // Name of input variables

   nameOfInputVariables.resize(numberOfInputVariables);

   for(int i = 0; i < numberOfInputVariables; i++)
   {
      std::stringstream buffer;

      buffer << "InputVariable" << i;

      nameOfInputVariables[i] = buffer.str();
   }

   // Name of target variables

   nameOfTargetVariables.resize(numberOfTargetVariables);

   for(int i = 0; i < numberOfTargetVariables; i++)
   {
      std::stringstream buffer;

      buffer << "TargetVariable" << i;

      nameOfTargetVariables[i] = buffer.str();
   }

   // Units of input variables

   unitsOfInputVariables.resize(numberOfInputVariables);

   for(int i = 0; i < numberOfInputVariables; i++)
   {
      unitsOfInputVariables[i] = "None";
   }

   // Units of target variables

   unitsOfTargetVariables.resize(numberOfTargetVariables);

   for(int i = 0; i < numberOfTargetVariables; i++)
   {
      unitsOfTargetVariables[i] = "None";
   }

   // Description of input variables

   descriptionOfInputVariables.resize(numberOfInputVariables);

   for(int i = 0; i < numberOfInputVariables; i++)
   {
      descriptionOfInputVariables[i] = "None";
   }

   // Description of target variables

   descriptionOfTargetVariables.resize(numberOfTargetVariables);

   for(int i = 0; i < numberOfTargetVariables; i++)
   {
      descriptionOfTargetVariables[i] = "None";
   }

   // Input data 

   inputData.resize(numberOfSamples, numberOfInputVariables);  

   // Target data

   targetData.resize(numberOfSamples, numberOfTargetVariables);  
}


// void setInputData(Matrix<double>) method

/// This method sets a new input data matrix. 
/// The number of rows must be equal to the number of samples.
/// The number of columns must be equal to the number of input variables.
///
/// @param newInputData Input data matrix.

void InputTargetDataSet::setInputData(Matrix<double> newInputData)
{
   // Control sentence (if debug)

   #ifndef NDEBUG 

   if(newInputData.getNumberOfRows() != numberOfSamples)
   {
      std::cerr << std::endl
                << "Flood Error: InputTargetDataSet class." << std::endl
                << "void setInputData(void) method." << std::endl
                << "Number of rows must be equal to number of samples." << std::endl
                << std::endl;

      exit(1);
   }
   else if(newInputData.getNumberOfColumns() != numberOfInputVariables)
   {
      std::cerr << std::endl
                << "Flood Error: InputTargetDataSet class." << std::endl
                << "void setInputData(void) method." << std::endl
                << "Number of columns must be equal to number of input variables." << std::endl
                << std::endl;

      exit(1);
   }

   #endif

   // Set input data
   
   inputData = newInputData;   
}


// void setTargetData(Matrix<double>) method

/// This method sets a new target data matrix. 
/// The number of rows must be equal to the number of samples.
/// The number of columns must be equal to the number of target variables.
///
/// @param newTargetData Target data matrix.

void InputTargetDataSet::setTargetData(Matrix<double> newTargetData)
{
   // Control sentence (if debug)

   #ifndef NDEBUG 
 
   if(newTargetData.getNumberOfRows() != numberOfSamples)
   {
      std::cout << std::endl
                << "Flood Error: InputTargetDataSet class." << std::endl
                << "void setTargetData(void) method." << std::endl
                << "Number of rows must be equal to number of samples." << std::endl
                << std::endl;

      exit(1);
   }
   else if(newTargetData.getNumberOfColumns() != numberOfTargetVariables)
   {
      std::cout << std::endl
                << "Flood Error: InputTargetDataSet class." << std::endl
                << "void setTargetData(void) method." << std::endl
                << "Number of columns must be equal to number of target variables." << std::endl
                << std::endl;

      exit(1);
   }

   #endif

   // Set target data

   targetData = newTargetData;   
   
}


// void setInputAndTargetData(Matrix<double>) method

/// This method sets both a new input set and a new target set from a single matrix. 
/// The number of rows must be equal to the number of samples.
/// The number of columns must be equal to the total number of input and target variables.
///
/// @param newInputAndTargetData Input and target data matrix.

void InputTargetDataSet::setInputAndTargetData(Matrix<double> newInputAndTargetData)
{
   // Control sentence (if debug)

   #ifndef NDEBUG 

   int numberOfRows = newInputAndTargetData.getNumberOfRows();
   int numberOfColumns = newInputAndTargetData.getNumberOfColumns();

   if(numberOfRows != numberOfSamples)
   {
      std::cerr << std::endl
                << "Flood Error: InputTargetDataSet class." << std::endl
                << "void setInputAndTargetData(void) method." << std::endl
                << "Number of rows must be equal to number of samples." << std::endl
                << std::endl;

      exit(1);
   }
   else if(numberOfColumns != numberOfInputVariables + numberOfTargetVariables)
   {
      std::cerr << std::endl
                << "Flood Error: InputTargetDataSet class." << std::endl
                << "void setInputAndTargetData(void) method." << std::endl
                << "Number of columns must be equal to " << std::endl
                << "number of input variables  plus number of target variables." 
                << std::endl
                << std::endl;

      exit(1);
   }
   
   #endif

   for(int i = 0; i < numberOfSamples; i++)
   {
      for(int j = 0; j < numberOfInputVariables; j++)
      {
         inputData[i][j] = newInputAndTargetData[i][j];            
      }
      for(int k = numberOfInputVariables; k < numberOfInputVariables+numberOfTargetVariables; k++)
      {
         targetData[i][k-numberOfInputVariables] = newInputAndTargetData[i][k];                          
      }
   }
}


// void setNameOfInputVariables(Vector<std::string>) method

/// This method sets a new name for the input variables in the input-target data set. 
///
/// @param newNameOfInputVariables Names for the input variables.

void InputTargetDataSet::setNameOfInputVariables(Vector<std::string> newNameOfInputVariables)
{
   // Control sentence (if debug)

   #ifndef NDEBUG 
 
   if(newNameOfInputVariables.getSize() != numberOfInputVariables)
   {
      std::cerr << std::endl
                << "Flood Error: InputTargetDataSet class." << std::endl
                << "void setNameOfInputVariables(Vector<std::string>) method." << std::endl
                << "Size of vector must be equal to number of input variables." << std::endl
                << std::endl;

      exit(1);
   }

   #endif

   // Set name of input variables

   nameOfInputVariables = newNameOfInputVariables;
}


// void setNameOfTargetVariables(Vector<std::string>) method

/// This method sets a new name for the target variables in the input-target data set. 
///
/// @param newNameOfTargetVariables Names for the target variables.

void InputTargetDataSet::setNameOfTargetVariables(Vector<std::string> newNameOfTargetVariables)
{
   // Control sentence (if debug)

   #ifndef NDEBUG 

   if(newNameOfTargetVariables.getSize() != numberOfTargetVariables)
   {
      std::cerr << std::endl
                << "Flood Error: InputTargetDataSet class." << std::endl
                << "void setNameOfTargetVariables(Vector<std::string>) method." << std::endl
                << "Size of vector must be equal to number of target variables." << std::endl
                << std::endl;

      exit(1);
   }

   #endif

   // Set name of target variables

   nameOfTargetVariables = newNameOfTargetVariables;
}


// void setNameOfSingleInputVariable(int, std::string) method

/// This method sets the name of a single input variable in the data set.
/// Such value is only used to give the user basic information on the problem at hand.
///
/// @param i Index of input variable.
/// @param newNameOfSingleInputVariable New name for the input variable with index i.
///
/// @see setNameOfInputVariables(Vector<std::string>).
/// @see setNameOfTargetVariables(Vector<std::string>).
/// @see setNameOfSingleTargetVariable(int, std::string).

void InputTargetDataSet::setNameOfSingleInputVariable(int i, std::string newNameOfSingleInputVariable)
{
   // Control sentence (if debug)

   #ifndef NDEBUG 

   if(i >= numberOfInputVariables)
   {
      std::cerr << std::endl
                << "Flood Error: InputTargetDataSet class." << std::endl
                << "void setNameOfSingleInputVariable(int, std::string) method." << std::endl 
                << "Index must be less than number of input variables." << std::endl
                << std::endl;

      exit(1);   
   }

   #endif

   // Set name of single input variable

   nameOfInputVariables[i] = newNameOfSingleInputVariable;
}


// void setNameOfSingleTargetVariable(int, std::string) method

/// This method sets the name of a single target variable in the data set.
/// Such value is only used to give the user basic information on the problem at hand.
///
/// @param i Index of target variable.
/// @param newNameOfSingleTargetVariable New name for the target variable with index i.
///
/// @see setNameOfInputVariables(Vector<std::string>).
/// @see setNameOfTargetVariables(Vector<std::string>).
/// @see setNameOfSingleInputVariable(int, std::string).

void InputTargetDataSet::setNameOfSingleTargetVariable(int i, std::string newNameOfSingleTargetVariable)
{
   // Control sentence (if debug)

   #ifndef NDEBUG 

   if(i >= numberOfTargetVariables)
   {
      std::cerr << std::endl
                << "Flood Error: InputTargetDataSet class." << std::endl
                << "void setNameOfSingleTargetVariable(int, std::string) method." << std::endl 
                << "Index must be less than number of target variables." << std::endl 
                << std::endl;

      exit(1);   
   }

   #endif

   // Set name of single input variable

   nameOfTargetVariables[i] = newNameOfSingleTargetVariable;
}


// void setUnitsOfInputVariables(Vector<std::string>) method

/// This method sets new units for the input variables in the input-target data set. 
///
/// @param newUnitsOfInputVariables Units for the input variables.

void InputTargetDataSet::setUnitsOfInputVariables(Vector<std::string> newUnitsOfInputVariables)
{
   // Control sentence (if debug)

   #ifndef NDEBUG 

   if(newUnitsOfInputVariables.getSize() != numberOfInputVariables)
   {
      std::cerr << std::endl
                << "Flood Error: InputTargetDataSet class." << std::endl
                << "void setUnitsOfInputVariables(Vector<std::string>) method." << std::endl
                << "Size of vector must be equal to number of input variables." << std::endl
                << std::endl;

      exit(1);
   }

   #endif

   // Set units of input variables

   unitsOfInputVariables = newUnitsOfInputVariables;
}


// void setUnitsOfTargetVariables(Vector<std::string>) method

/// This method sets new units for the target variables in the input-target data set. 
///
/// @param newUnitsOfTargetVariables Units for the target variables.

void InputTargetDataSet::setUnitsOfTargetVariables(Vector<std::string> newUnitsOfTargetVariables)
{
   // Control sentence (if debug)

   #ifndef NDEBUG 

   if(newUnitsOfTargetVariables.getSize() != numberOfTargetVariables)
   {
      std::cout << std::endl
                << "Flood Error: InputTargetDataSet class." << std::endl
                << "void setUnitsOfTargetVariables(Vector<std::string>) method." << std::endl
                << "Size of vector must be equal to number of target variables." << std::endl
                << std::endl;

      exit(1);
   }

   #endif

   // Set units of target variables

   unitsOfTargetVariables = newUnitsOfTargetVariables;
}


// void setUnitsOfSingleInputVariable(int, std::string) method

/// This method sets new units for a single input variable in the data set.
/// Such value is only used to give the user basic information on the problem at hand.
///
/// @param i Index of input variable.
/// @param newUnitsOfSingleInputVariable New units for the input variable with index i.
///
/// @see setUnitsOfInputVariables(Vector<std::string>).
/// @see setUnitsOfTargetVariables(Vector<std::string>).
/// @see setUnitsOfSingleTargetVariable(int, std::string).

void InputTargetDataSet::setUnitsOfSingleInputVariable(int i, std::string newUnitsOfSingleInputVariable)
{
   // Control sentence (if debug)

   #ifndef NDEBUG 

   if(i >= numberOfInputVariables)
   {
      std::cerr << std::endl
                << "Flood Error: InputTargetDataSet class." << std::endl 
                << "void setUnitsOfSingleInputVariable(int, std::string) method." << std::endl
                << "Index must be less than number of input variables." << std::endl 
				<< std::endl;

      exit(1);
   }

   #endif

   // Set units of single input variable

   unitsOfInputVariables[i] = newUnitsOfSingleInputVariable;
}


// void setUnitsOfSingleTargetVariable(int, std::string) method

/// This method sets new units for a single target variable in the data set.
/// Such value is only used to give the user basic information on the problem at hand.
///
/// @param i Index of target variable.
/// @param newUnitsOfSingleTargetVariable New units for the target variable with index i.
///
/// @see setUnitsOfInputVariables(Vector<std::string>).
/// @see setUnitsOfTargetVariables(Vector<std::string>).
/// @see setUnitsOfSingleInputVariable(int, std::string).

void InputTargetDataSet::setUnitsOfSingleTargetVariable(int i, std::string newUnitsOfSingleTargetVariable)
{
   // Control sentence (if debug)

   #ifndef NDEBUG 

   if(i >= numberOfTargetVariables)
   {
      std::cerr << std::endl
                << "Flood Error: InputTargetDataSet class." << std::endl
                << "void setUnitsOfSingleTargetVariable(int, std::string) method." << std::endl
                << "Index must be less than number of target variables." << std::endl 
                << std::endl;

      exit(1);   
   }

   #endif

   // Set units of single target variable

   unitsOfTargetVariables[i] = newUnitsOfSingleTargetVariable;
}


// void setDescriptionOfInputVariables(Vector<std::string>) method

/// This method sets a new description for the input variables in the input-target data set. 
///
/// @param newDescriptionOfInputVariables Description for the input variables.

void InputTargetDataSet::setDescriptionOfInputVariables(Vector<std::string> newDescriptionOfInputVariables)
{
   // Control sentence (if debug)

   #ifndef NDEBUG 
 
   if(newDescriptionOfInputVariables.getSize() != numberOfInputVariables)
   {
      std::cerr << std::endl
                << "Flood Error: InputTargetDataSet class." << std::endl
                << "void setDescriptionOfInputVariables(Vector<std::string>) method." << std::endl
                << "Size of vector must be equal to number of input variables." << std::endl
                << std::endl;

      exit(1);
   }

   #endif

   // Set description of input variables

   descriptionOfInputVariables = newDescriptionOfInputVariables;
}


// void setDescriptionOfTargetVariables(Vector<std::string>) method

/// This method sets a new description for the target variables in the input-target data set. 
///
/// @param newDescriptionOfTargetVariables Description for the target variables.

void InputTargetDataSet::setDescriptionOfTargetVariables(Vector<std::string> newDescriptionOfTargetVariables)
{
   // Control sentence (if debug)

   #ifndef NDEBUG 

   if(newDescriptionOfTargetVariables.getSize() != numberOfTargetVariables)
   {
      std::cerr << std::endl
                << "Flood Error: InputTargetDataSet class." << std::endl
                << "void setDescriptionOfTargetVariables(Vector<std::string>) method." << std::endl
                << "Size of vector must be equal to number of target variables." << std::endl
                << std::endl;

      exit(1);
   }

   #endif

   // Set description of target variables

   descriptionOfTargetVariables = newDescriptionOfTargetVariables;
}


// void setDescriptionOfSingleInputVariable(int, std::string) method

/// This method sets a new description for a single input variable in the data set.
/// Such value is only used to give the user basic information on the problem at hand.
///
/// @param i Index of input variable.
/// @param newDescriptionOfSingleInputVariable New description for the input variable with index i.
///
/// @see setDescriptionOfInputVariables(Vector<std::string>).
/// @see setDescriptionOfTargetVariables(Vector<std::string>).
/// @see setDescriptionOfSingleTargetVariable(int, std::string).

void InputTargetDataSet::setDescriptionOfSingleInputVariable(int i, std::string newDescriptionOfSingleInputVariable)
{
   // Control sentence (if debug)

   #ifndef NDEBUG 

   if(i >= numberOfInputVariables)
   {
      std::cerr << std::endl
                << "Flood Error: InputTargetDataSet class." << std::endl
                << "void setDescriptionOfSingleInputVariable(int, std::string) method." << std::endl
                << "Index must be less than number of input variables." << std::endl 
                << std::endl;

      exit(1);   
   }

   #endif

   // Set description of single input variable

   descriptionOfInputVariables[i] = newDescriptionOfSingleInputVariable;
}


// void setDescriptionOfSingleTargetVariable(int, std::string) method

/// This method sets a new description for a single target variable in the data set.
/// Such value is only used to give the user basic information on the problem at hand.
///
/// @param i Index of target variable.
/// @param newDescriptionOfSingleTargetVariable New description for the target variable with index i.
///
/// @see setDescriptionOfInputVariables(Vector<std::string>).
/// @see setDescriptionOfTargetVariables(Vector<std::string>).
/// @see setDescriptionOfSingleInputVariable(int, std::string).

void InputTargetDataSet::setDescriptionOfSingleTargetVariable(int i, std::string newDescriptionOfSingleTargetVariable)
{
   // Control sentence (if debug)

   #ifndef NDEBUG 

   if(i >= numberOfTargetVariables)
   {
      std::cerr << std::endl
                << "Flood Error: InputTargetDataSet class." << std::endl
                << "void setDescriptionOfSingleTargetVariable(int, std::string) method." << std::endl
                << "Index must be less than number of targets." << std::endl 
                << std::endl;

      exit(1);   
   }

   #endif

   // Set description of single target variable

   descriptionOfTargetVariables[i] = newDescriptionOfSingleTargetVariable;
}


// void setDisplay(bool) method

/// This method sets a new display value. 
/// If it is set to true messages from this class are to be displayed on the screen;
/// if it is set to false messages from this class are not to be displayed on the screen.
///
/// @param newDisplay Display value.

void InputTargetDataSet::setDisplay(bool newDisplay)
{
   display = newDisplay;
}


// Matrix<double> calculateMeanAndStandardDeviationOfInputData(void) method

/// This method returns a matrix containing both the mean and the standard deviation 
/// of all the input variables.
/// The first row contains the mean of the input variables.
/// The second row contains the standard deviation of the input variables.

Matrix<double> InputTargetDataSet::calculateMeanAndStandardDeviationOfInputData(void)
{
   // Control sentence (if debug)

   #ifndef NDEBUG 

   if(numberOfSamples == 1)
   {
      std::cerr << std::endl
                << "Flood Error: InputTargetDataSet class." 
                << "Matrix<double> calculateMeanAndStandardDeviationOfInputData(void) method." << std::endl 
                << "Number of samples must be greater than one." << std::endl 
                << std::endl;

      exit(1);
   }

   #endif

   // Calculate mean and standard deviation

   Matrix<double> meanAndStandardDeviationOfInputData(2, numberOfInputVariables, 0.0);

   // Mean of input variables   

   Vector<double> meanOfInputData(numberOfInputVariables, 0.0);

   for(int j = 0; j < numberOfInputVariables; j++)
   {
      meanOfInputData[j] = 0;

      for(int i = 0; i < numberOfSamples; i++)
      {
         meanOfInputData[j] = meanOfInputData[j] + inputData[i][j];
      }

      meanOfInputData[j] = meanOfInputData[j]/(double)numberOfSamples;
   }

   // Standard deviation of input data

   Vector<double> standardDeviationOfInputData(numberOfInputVariables, 0.0);

   for(int j = 0; j < numberOfInputVariables; j++)
   {
      standardDeviationOfInputData[j] = 0;

      for(int i = 0; i < numberOfSamples; i++)
      {
         standardDeviationOfInputData[j] += pow(inputData[i][j] - meanOfInputData[j], 2);
      }

      standardDeviationOfInputData[j] = sqrt(standardDeviationOfInputData[j]/(numberOfSamples-1.0));
   }

   // Mean and standard deviation of input data

   meanAndStandardDeviationOfInputData.setRow(0, meanOfInputData);
   meanAndStandardDeviationOfInputData.setRow(1, standardDeviationOfInputData);

   return(meanAndStandardDeviationOfInputData);
}


// Matrix<double> calculateMeanAndStandardDeviationOfTargetData(void)

/// This method returns a matrix containing both the mean and the standard deviation 
/// of all the target variables.
/// The first row contains the mean of the target variables.
/// The second row contains the standard deviation of the target variables.

Matrix<double> InputTargetDataSet::calculateMeanAndStandardDeviationOfTargetData(void)
{
   Matrix<double> meanAndStandardDeviationOfTargetData(2, numberOfTargetVariables, 0.0);

   // Mean of target data

   Vector<double> meanOfTargetData(numberOfTargetVariables, 0.0);

   for(int j = 0; j < numberOfTargetVariables; j++)
   {
      meanOfTargetData[j] = 0;

      for(int i = 0; i < numberOfSamples; i++)
      {
         meanOfTargetData[j] += targetData[i][j];
      }

      meanOfTargetData[j] = meanOfTargetData[j]/(double)numberOfSamples;
    }

   // Standard deviation of target data

   Vector<double> standardDeviationOfTargetData(numberOfTargetVariables, 0.0);


   for(int j = 0; j < numberOfTargetVariables; j++)
   {
      standardDeviationOfTargetData[j] = 0;

      for(int i = 0; i < numberOfSamples; i++)
      {
         standardDeviationOfTargetData[j] += pow(targetData[i][j] - meanOfTargetData[j], 2);
      }

      standardDeviationOfTargetData[j] = sqrt(standardDeviationOfTargetData[j]/(numberOfSamples-1.0));
   }

   // Mean and standard deviation of target data

   meanAndStandardDeviationOfTargetData.setRow(0, meanOfTargetData);
   meanAndStandardDeviationOfTargetData.setRow(1, standardDeviationOfTargetData);

   return(meanAndStandardDeviationOfTargetData);
}


// Matrix<double> calculateMinimumAndMaximumOfInputData(void)

/// This method returns a matrix containing both the minimum and the maximum 
/// of all the input variables.
/// The first row contains the minimum of the input variables.
/// The second row contains the maximum of the input variables.

Matrix<double> InputTargetDataSet::calculateMinimumAndMaximumOfInputData(void)
{
   Matrix<double> minimumAndMaximumOfInputData(2, numberOfInputVariables, 0.0);

   Vector<double> minimumOfInputData(numberOfInputVariables,  1.0e69);
   Vector<double> maximumOfInputData(numberOfInputVariables, -1.0e69);

   for(int j = 0; j < numberOfInputVariables; j++)
   {
      for(int i = 0; i < numberOfSamples; i++)
      {    
         // Minimum of input data

         if(inputData[i][j] < minimumOfInputData[j])
         {
            minimumOfInputData[j] = inputData[i][j];
         }

         // Maximum of input data

         if(inputData[i][j] > maximumOfInputData[j])
         {
            maximumOfInputData[j] = inputData[i][j];
         }
      }
   }

   // Minimum and maximum of input data

   minimumAndMaximumOfInputData.setRow(0, minimumOfInputData);
   minimumAndMaximumOfInputData.setRow(1, maximumOfInputData);

   return(minimumAndMaximumOfInputData);
}


// Matrix<double> calculateMinimumAndMaximumOfTargetData(void)

/// This method returns a matrix containing both the minimum and the maximum 
/// of all the target variables.
/// The first row contains the minimum of the target variables.
/// The second row contains the maximum of the target variables.

Matrix<double> InputTargetDataSet::calculateMinimumAndMaximumOfTargetData(void)
{
   Matrix<double> minimumAndMaximumOfTargetData(2, numberOfTargetVariables, 0.0);

   Vector<double> minimumOfTargetData(numberOfTargetVariables,  1.0e69);
   Vector<double> maximumOfTargetData(numberOfTargetVariables, -1.0e69);

   for(int j = 0; j < numberOfTargetVariables; j++)
   {
      for(int i = 0; i < numberOfSamples; i++)
      {    
         // Minimum of input data

         if(targetData[i][j] < minimumOfTargetData[j])
         {
            minimumOfTargetData[j] = targetData[i][j];
         }

         // Maximum of input data

         if(targetData[i][j] > maximumOfTargetData[j])
         {
            maximumOfTargetData[j] = targetData[i][j];
         }
      }
   }

   // Minimum and maximum of target data

   minimumAndMaximumOfTargetData.setRow(0, minimumOfTargetData);
   minimumAndMaximumOfTargetData.setRow(1, maximumOfTargetData);
  
   return(minimumAndMaximumOfTargetData);
}

   
// Vector< Vector<double> > calculateAllStatistics(void) method

Vector< Vector<double> > InputTargetDataSet::calculateAllStatistics(void)
{
   Vector< Vector<double> > statistics(8);

   Matrix<double> meanAndStandardDeviationOfInputData = calculateMeanAndStandardDeviationOfInputData();
   Matrix<double> meanAndStandardDeviationOfTargetData = calculateMeanAndStandardDeviationOfTargetData();

   Matrix<double> minimumAndMaximumOfInputData = calculateMinimumAndMaximumOfInputData();
   Matrix<double> minimumAndMaximumOfTargetData = calculateMinimumAndMaximumOfTargetData();

   // Mean of input data

   statistics[0] = meanAndStandardDeviationOfInputData.getRow(0);

   // Standard deviation of input data 

   statistics[1] = meanAndStandardDeviationOfInputData.getRow(1);

   // Mean of target data

   statistics[2] = meanAndStandardDeviationOfTargetData.getRow(0);

   // Standard deviation of target data

   statistics[3] = meanAndStandardDeviationOfTargetData.getRow(1);

   // Minimum of input data 

   statistics[4] = minimumAndMaximumOfInputData.getRow(0);

   // Maximum of input data

   statistics[5] = minimumAndMaximumOfInputData.getRow(1);

   // Minimum of target data 

   statistics[6] = minimumAndMaximumOfTargetData.getRow(0);

   // Maximum of target data

   statistics[7] = minimumAndMaximumOfTargetData.getRow(1);

   return(statistics);
}


// void preprocessMeanAndStandardDeviation(void) method

/// This method preprocesses all input and target data so that the
/// mean of all input and target variables is zero and 
/// their standard deviation is one. It then updates:
///
/// <ul>
/// <li> Input data matrix. 
/// <li> Target data matrix.
/// </ul>

/// @see preprocessMeanAndStandardDeviation(Matrix<double>, Matrix<double>)

void InputTargetDataSet::preprocessMeanAndStandardDeviation(void)
{
   Matrix<double> newInputData(numberOfSamples, numberOfInputVariables, 0.0);
 
   Matrix<double> newTargetData(numberOfSamples, numberOfTargetVariables, 0.0);

   // Mean and standard deviation of input data

   Matrix<double> meanAndStandardDeviationOfInputData = calculateMeanAndStandardDeviationOfInputData();

   Vector<double> meanOfInputData = meanAndStandardDeviationOfInputData.getRow(0);
   
   Vector<double> standardDeviationOfInputData = meanAndStandardDeviationOfInputData.getRow(1);

   // Mean and standard deviation of target data

   Matrix<double> meanAndStandardDeviationOfTargetData = calculateMeanAndStandardDeviationOfTargetData();

   Vector<double> meanOfTargetData = meanAndStandardDeviationOfTargetData.getRow(0);
   
   Vector<double> standardDeviationOfTargetData = meanAndStandardDeviationOfTargetData.getRow(1);

   // Rescale input data

   for(int j = 0; j < numberOfInputVariables; j++)
   {
      if(standardDeviationOfInputData[j] < 1e-9)
      {
         if(display)
         {                                          
            std::cout << std::endl
                      << "Flood Warning: InputTargetDataSet class." << std::endl
                      << "void preprocessMeanAndStandardDeviation(void) method." << std::endl
                      << "Standard deviation of input variable " << j << " is zero." << std::endl
                      << "Those inputs won't be transformed." << std::endl;
         }

         for(int i = 0; i < numberOfSamples; i++)
         {
            newInputData[i][j] = inputData[i][j];
         }
      }
      else
      {
         for(int i = 0; i < numberOfSamples; i++)
         {
            newInputData[i][j] = (inputData[i][j] - meanOfInputData[j])/standardDeviationOfInputData[j];
         }
      }
   }

   inputData = newInputData;


   // Rescale target data

   for(int j = 0; j < numberOfTargetVariables; j++)
   {
      if(standardDeviationOfTargetData[j] < 1e-9)
      {
         if(display)
         {
            std::cout << std::endl
                      << "Flood Warning: InputTargetDataSet class." << std::endl
                      << "void preprocessMeanAndStandardDeviation(void) method." << std::endl   
                      << "Standard deviation of target variable " <<  j << " is zero." << std::endl
                      << "Those targets won't be transformed." << std::endl;
         }
         
         for(int i = 0; i < numberOfSamples; i++)
         {
            newTargetData[i][j] = targetData[i][j];
         }
      }
      else
      {
         for(int i = 0; i < numberOfSamples; i++)
         {
            newTargetData[i][j] = (targetData[i][j] - meanOfTargetData[j])/standardDeviationOfTargetData[j];
         }
      }
   }

   targetData = newTargetData;
}


// void preprocessMeanAndStandardDeviation(Matrix<double>, Matrix<double>) method

/// This method preprocesses all input and target data with given mean and standardard 
/// deviations for the input and target variables. It then updates:
///
/// <ul>
/// <li> Input data matrix. 
/// <li> Target data matrix.
/// </ul>
///
/// @param meanAndStandardDeviationOfInputData 
/// Mean and standard deviation values for the input variables to be used for preprocessing.
/// @param meanAndStandardDeviationOfTargetData 
/// Mean and standard deviation values for the target variables to be used for preprocessing.
///
/// @see preprocessMeanAndStandardDeviation(void)

void InputTargetDataSet::preprocessMeanAndStandardDeviation(Matrix<double> meanAndStandardDeviationOfInputData, Matrix<double> meanAndStandardDeviationOfTargetData)
{
   Matrix<double> newInputData(numberOfSamples, numberOfInputVariables, 0.0); 
   Matrix<double> newTargetData(numberOfSamples, numberOfTargetVariables, 0.0);

   // Mean and standard deviation of input data

   Vector<double> meanOfInputData = meanAndStandardDeviationOfInputData.getRow(0);
   
   Vector<double> standardDeviationOfInputData = meanAndStandardDeviationOfInputData.getRow(1);

   // Mean and standard deviation of target data

   Vector<double> meanOfTargetData = meanAndStandardDeviationOfTargetData.getRow(0);
   
   Vector<double> standardDeviationOfTargetData = meanAndStandardDeviationOfTargetData.getRow(1);

   // Rescale input data

   for(int j = 0; j < numberOfInputVariables; j++)
   {
      if(standardDeviationOfInputData[j] < 1e-9)
      {
         if(display)
         {                                          
            std::cout << std::endl
                      << "Flood Warning: InputTargetDataSet class." << std::endl
                      << "void preprocessMeanAndStandardDeviation(void) method." << std::endl
                      << "Standard deviation of input variable " << j << " is zero." << std::endl
                      << "Those inputs won't be transformed."
                      << std::endl;
         }

         for(int i = 0; i < numberOfSamples; i++)
         {
            newInputData[i][j] = inputData[i][j];
         }
      }
      else
      {
         for(int i = 0; i < numberOfSamples; i++)
         {
            newInputData[i][j] = (inputData[i][j] - meanOfInputData[j])/standardDeviationOfInputData[j];
         }
      }
   }

   inputData = newInputData;

   // Rescale target data

   for(int j = 0; j < numberOfTargetVariables; j++)
   {
      if(standardDeviationOfTargetData[j] < 1e-9)
      {
         if(display)
         {
            std::cout << std::endl
                      << "Flood Warning: InputTargetDataSet class." << std::endl
                      << "void preprocessMeanAndStandardDeviation(Matrix<double>, Matrix<double>) method." 
                      << std::endl   
                      << "Standard deviation of target variable " <<  j << " is zero." << std::endl
                      << "Those targets won't be transformed."
                      << std::endl;
         }
         
         for(int i = 0; i < numberOfSamples; i++)
         {
            newTargetData[i][j] = targetData[i][j];
         }
      }
      else
      {
         for(int i = 0; i < numberOfSamples; i++)
         {
            newTargetData[i][j] = (targetData[i][j] - meanOfTargetData[j])/standardDeviationOfTargetData[j];
         }
      }
   }

   targetData = newTargetData;     
}


// void preprocessMinimumAndMaximum(void) method

/// This method preprocesses all input and target data so that the
/// minimum of all input and target variables is -1 and 
/// their maximum is 1. It then updates:
///
/// <ul>
/// <li> Input data matrix. 
/// <li> Target data matrix.
/// </ul>
///
/// @see preprocessMinimumAndMaximum(Matrix<double>, Matrix<double>).

void InputTargetDataSet::preprocessMinimumAndMaximum(void)
{
   Matrix<double> newInputData(numberOfSamples, numberOfInputVariables, 0.0);
   Matrix<double> newTargetData(numberOfSamples, numberOfTargetVariables, 0.0);

   // Mean and standard deviation of input data

   Matrix<double> minimumAndMaximumOfInputData 
   = calculateMinimumAndMaximumOfInputData();

   Vector<double> minimumOfInputData = minimumAndMaximumOfInputData.getRow(0);
   Vector<double> maximumOfInputData = minimumAndMaximumOfInputData.getRow(1);

   // Minimum and maximum of target data

   Matrix<double> minimumAndMaximumOfTargetData
   = calculateMinimumAndMaximumOfTargetData();

   Vector<double> minimumOfTargetData = minimumAndMaximumOfTargetData.getRow(0);
   Vector<double> maximumOfTargetData = minimumAndMaximumOfTargetData.getRow(1);

   // Rescale input data

   for(int j = 0; j < numberOfInputVariables; j++)
   {
      if(maximumOfInputData[j] - minimumOfInputData[j] < 1e-9)
      {
         if(display)
         {
            std::cout << std::endl
                      << "Flood Warning: InputTargetDataSet class." << std::endl
                      << "void preprocessMinimumAndMaximum(void) method." << std::endl
                      << "Minimum and maximum values of input variable " << j << " are equal." << std::endl
                      << "Those inputs won't be transformed." << std::endl;
         }

         for(int i = 0; i < numberOfSamples; i++)
         {
            newInputData[i][j] = inputData[i][j];
         }
      }
      else
      {
         for(int i = 0; i < numberOfSamples; i++)
         {
            newInputData[i][j] 
            = 2.0*(inputData[i][j] - minimumOfInputData[j])/(maximumOfInputData[j]-minimumOfInputData[j])-1.0;
         }
      }
   }

   inputData = newInputData;


   // Rescale target data

   for(int j = 0; j < numberOfTargetVariables; j++)
   {
      if(maximumOfTargetData[j] - minimumOfTargetData[j] < 1e-9)
      {
         if(display)
         {
            std::cout << std::endl
                      << "Flood Warning: InputTargetDataSet class." << std::endl
                      << "void preprocessMinimumAndMaximum(void) method." << std::endl   
                      << "Minimum and maximum values of target variable " << j << " are equal." << std::endl
                      << "Those targets won't be transformed." << std::endl;
         }
         
         for(int i = 0; i < numberOfSamples; i++)
         {
            newTargetData[i][j] = targetData[i][j];
         }
      }
      else
      {
         for(int i = 0; i < numberOfSamples; i++)
         {
            newTargetData[i][j] = 2.0*(targetData[i][j] - minimumOfTargetData[j])
            /(maximumOfTargetData[j]-minimumOfTargetData[j]) - 1.0;
         }
      }
   }

   targetData = newTargetData;
   
}


// void preprocessMinimumAndMaximum(Matrix<double>, Matrix<double>) method

/// This method preprocesses all input and target data with given minimums and maximums
/// for the input and target variables. It then updates:
///
/// <ul>
/// <li> Input data matrix. 
/// <li> Target data matrix.
/// </ul>
///
/// @param minimumAndMaximumOfInputData 
/// Minimum and maximum values for the input variables to be used for preprocessing.
/// @param minimumAndMaximumOfTargetData 
/// Minimum and maximum values for the target variables to be used for preprocessing.
///
/// @see preprocessMinimumAndMaximum(void)

void InputTargetDataSet::preprocessMinimumAndMaximum(Matrix<double> minimumAndMaximumOfInputData, Matrix<double> minimumAndMaximumOfTargetData)
{
   Matrix<double> newInputData(numberOfSamples, numberOfInputVariables, 0.0);
   Matrix<double> newTargetData(numberOfSamples, numberOfTargetVariables, 0.0);

   // Mean and standard deviation of input data

   Vector<double> minimumOfInputData = minimumAndMaximumOfInputData.getRow(0);
   Vector<double> maximumOfInputData = minimumAndMaximumOfInputData.getRow(1);

   // Minimum and maximum of target data

   Vector<double> minimumOfTargetData = minimumAndMaximumOfTargetData.getRow(0);
   Vector<double> maximumOfTargetData = minimumAndMaximumOfTargetData.getRow(1);

   // Rescale input data

   for(int j = 0; j < numberOfInputVariables; j++)
   {
      if(maximumOfInputData[j] - minimumOfInputData[j] < 1e-9)
      {
         if(display)
         {
            std::cout << std::endl
                      << "Flood Warning: InputTargetDataSet class." << std::endl
                      << "void preprocessMinimumAndMaximum(Matrix<double>, Matrix<double>) method." << std::endl
                 << "Minimum and maximum values of input variable " << j << " are equal." << std::endl
                      << "Those inputs won't be transformed."
                      << std::endl;
         }

         for(int i = 0; i < numberOfSamples; i++)
         {
            newInputData[i][j] = inputData[i][j];
         }
      }
      else
      {
         for(int i = 0; i < numberOfSamples; i++)
         {
            newInputData[i][j] 
            = 2.0*(inputData[i][j] - minimumOfInputData[j])/(maximumOfInputData[j]-minimumOfInputData[j])-1.0;
         }
      }
   }

   inputData = newInputData;


   // Rescale target data

   for(int j = 0; j < numberOfTargetVariables; j++)
   {
      if(maximumOfTargetData[j] - minimumOfTargetData[j] < 1e-9)
      {
         if(display)
         {
            std::cout << std::endl
                      << "Flood Warning: InputTargetDataSet class." << std::endl
                      << "void preprocessMinimumAndMaximum(void) method." << std::endl   
                 << "Minimum and maximum values of target variable " << j << " are equal." << std::endl
                      << "Those targets won't be transformed."
                      << std::endl;
         }
         
         for(int i = 0; i < numberOfSamples; i++)
         {
            newTargetData[i][j] = targetData[i][j];
         }
      }
      else
      {
         for(int i = 0; i < numberOfSamples; i++)
         {
            newTargetData[i][j] 
            = 2.0*(targetData[i][j] - minimumOfTargetData[j])
            /(maximumOfTargetData[j]-minimumOfTargetData[j]) - 1.0;
         }
      }
   }

   targetData = newTargetData;     
}


// void postprocessMeanAndStandardDeviation(Matrix<double>, Matrix<double>) method 

/// This method postprocesses all input and target data with given mean and standardard 
/// deviations for the input and target variables. It then updates:
///
/// <ul>
/// <li> Input data matrix. 
/// <li> Target data matrix.
/// </ul>
///
/// @param meanAndStandardDeviationOfInputData 
/// Mean and standard deviation values for the input variables to be used for postprocessing.
/// @param meanAndStandardDeviationOfTargetData 
/// Mean and standard deviation values for the target variables to be used for postprocessing.

void InputTargetDataSet::postprocessMeanAndStandardDeviation(Matrix<double> meanAndStandardDeviationOfInputData, Matrix<double> meanAndStandardDeviationOfTargetData)
{
   Matrix<double> newInputData(numberOfSamples, numberOfInputVariables, 0.0);
   Matrix<double> newTargetData(numberOfSamples, numberOfTargetVariables, 0.0);

   // Mean and standard deviation of input data

   Vector<double> meanOfInputData 
   = meanAndStandardDeviationOfInputData.getRow(0);

   Vector<double> standardDeviationOfInputData 
   = meanAndStandardDeviationOfInputData.getRow(1);

   // Mean and standard deviation of target data

   Vector<double> meanOfTargetData 
   = meanAndStandardDeviationOfTargetData.getRow(0);
   
   Vector<double> standardDeviationOfTargetData 
   = meanAndStandardDeviationOfTargetData.getRow(1);

   // Postprocess input data

   for(int j = 0; j < numberOfInputVariables; j++)
   {
      if(standardDeviationOfInputData[j] < 1e-9)
      {
         if(display)
         {                                          
            std::cout << std::endl
                      << "Flood Warning: InputTargetDataSet class." << std::endl
                      << "void postprocessMeanAndStandardDeviation(Matrix<double>, Matrix<double>) method." << std::endl
                 << "Standard deviation of input variable " << j << " is zero." << std::endl
                      << "Those inputs won't be transformed."
                      << std::endl;
         }

         for(int i = 0; i < numberOfSamples; i++)
         {
            newInputData[i][j] = inputData[i][j];
         }
      }
      else
      {
         for(int i = 0; i < numberOfSamples; i++)
         {
            newInputData[i][j] = inputData[i][j]*standardDeviationOfInputData[j] 
            + meanOfInputData[j];
         }
      }
   }

   inputData = newInputData;


   // Rescale target data

   for(int j = 0; j < numberOfTargetVariables; j++)
   {
      if(standardDeviationOfTargetData[j] < 1e-9)
      {
         if(display)
         {
            std::cout << std::endl
                      << "Flood Warning: InputTargetDataSet class." << std::endl
                      << "void postprocessMeanAndStandardDeviation(Matrix<double>, Matrix<double>) method." << std::endl   
                 << "Standard deviation of target variable " <<  j << " is zero." << std::endl
                      << "Those targets won't be transformed."
                      << std::endl;
         }
         
         for(int i = 0; i < numberOfSamples; i++)
         {
            newTargetData[i][j] = targetData[i][j];
         }
      }
      else
      {
         for(int i = 0; i < numberOfSamples; i++)
         {
            newTargetData[i][j] = targetData[i][j]*standardDeviationOfTargetData[j] 
            + meanOfTargetData[j];
         }
      }
   }

   targetData = newTargetData;
}


// void postprocessMinimumAndMaximum(Matrix<double>, Matrix<double>) method

/// This method postprocesses all input and target data with given minimum and maximums
/// for the input and target variables. It then updates:
///
/// <ul>
/// <li> Input data matrix. 
/// <li> Target data matrix.
/// </ul>
///
/// @param minimumAndMaximumOfInputData 
/// Minimum and maximum values for the input variables to be used for postprocessing.
/// @param minimumAndMaximumOfTargetData 
/// Minimum and maximum values for the target variables to be used for postprocessing.

void InputTargetDataSet::postprocessMinimumAndMaximum(Matrix<double> minimumAndMaximumOfInputData, Matrix<double> minimumAndMaximumOfTargetData) 
{
   Matrix<double> newInputData(numberOfSamples, numberOfInputVariables, 0.0);
   Matrix<double> newTargetData(numberOfSamples, numberOfTargetVariables, 0.0);

   // Minimum and maximum of input data

   Vector<double> minimumOfInputData 
   = minimumAndMaximumOfInputData.getRow(0);

   Vector<double> maximumOfInputData 
   = minimumAndMaximumOfInputData.getRow(1);

   // Minimum and maximum of target data

   Vector<double> minimumOfTargetData = minimumAndMaximumOfTargetData.getRow(0);
   Vector<double> maximumOfTargetData = minimumAndMaximumOfTargetData.getRow(1);

   // Postprocess input data

   for(int j = 0; j < numberOfInputVariables; j++)
   {
      if(maximumOfInputData[j] - minimumOfInputData[j] < 1e-9)
      {
         if(display)
         {
            std::cout << std::endl
                      << "Flood Warning: InputTargetDataSet class." << std::endl
                      << "void postprocessMinimumAndMaximum(Matrix<double>, Matrix<double>) method." << std::endl
                 << "Minimum and maximum values of input variable " << j << " are equal." << std::endl
                      << "Those inputs won't be transformed."
                      << std::endl;
         }

         for(int i = 0; i < numberOfSamples; i++)
         {
            newInputData[i][j] = inputData[i][j];
         }
      }
      else
      {
         for(int i = 0; i < numberOfSamples; i++)
         {
         newInputData[i][j] = 0.5*(inputData[i][j] + 1.0)*(maximumOfInputData[j]-minimumOfInputData[j]) 
            + minimumOfInputData[j]; 
         }
      }
   }

   inputData = newInputData;


   // Postprocess target data

   for(int j = 0; j < numberOfTargetVariables; j++)
   {
      if(maximumOfTargetData[j] - minimumOfTargetData[j] < 1e-9)
      {
         if(display)
         {
            std::cout << std::endl
                      << "Flood Warning: InputTargetDataSet class." << std::endl
                      << "void postprocessMinimumAndMaximum(Matrix<double>, Matrix<double>) method." << std::endl   
                 << "Minimum and maximum values of target variable " << j << " are equal." << std::endl
                      << "Those targets won't be transformed."
                      << std::endl;
         }
         
         for(int i = 0; i < numberOfSamples; i++)
         {
            newTargetData[i][j] = targetData[i][j];
         }
      }
      else
      {
         for(int i = 0; i < numberOfSamples; i++)
         {
         newTargetData[i][j] = 0.5*(targetData[i][j] + 1.0)*(maximumOfTargetData[j]-minimumOfTargetData[j]) 
            + minimumOfTargetData[j]; 
         }
      }
   }

   targetData = newTargetData;
}


// void save(char*) method

/// This method saves the members of an input-target data set object to a data file:
///
/// <ul>
/// <li> Number of samples.
/// <li> Number of input variables.
/// <li> Number of target variables
/// <li> Units of input variables.
/// <li> Units of target variables
/// <li> Description of input variables. 
/// <li> Description of target variables.
/// <li> Input data.
/// <li> Target data.
/// </ul> 
///
/// @param filename Filename.
///
/// @see load(char*).

void InputTargetDataSet::save(char* filename)
{
   std::fstream file; 

   file.open(filename, std::ios::out);

   if(!file.is_open())
   {
      std::cerr << std::endl
                << "Flood Error: InputTargetDataSet class." << std::endl 
                << "void save(char*) method." << std::endl
                << "Cannot open input-target data set file."  << std::endl
                << std::endl;

      exit(1);
   }
   else
   {
      if(display)
      {
         std::cout << std::endl
                   << "Saving input-target data set to file..." << std::endl;
      }
   }
  
   // Write file header

   file << "% Flood Neural Network. Input-Target Data Set Object." << std::endl;

   // Write number of samples and number of input and target variables

   file << "NumberOfSamples:" << std::endl
        << numberOfSamples << std::endl
        << "NumberOfInputVariables:" << std::endl
        << numberOfInputVariables << std::endl
        << "NumberOfTargetVariables:" << std::endl
        << numberOfTargetVariables << std::endl;

   // Write name of input variables

   file << "NameOfInputVariables:" << std::endl
       << nameOfInputVariables << std::endl;
 
   // Write name of target variables

   file << "NameOfTargetVariables:" << std::endl
        << nameOfTargetVariables << std::endl;

   // Write units of input variables

   file << "UnitsOfInputVariables:" << std::endl
       << unitsOfInputVariables << std::endl;
 
   // Write units of target variables

   file << "UnitsOfTargetVariables:" << std::endl
        << unitsOfTargetVariables << std::endl;

   // Write description of input variables

   file << "DescriptionOfInputVariables:" << std::endl
       << descriptionOfInputVariables << std::endl;
 
   // Write description of target variables

   file << "DescriptionOfTargetVariables:" << std::endl
        << descriptionOfTargetVariables << std::endl;

   // Write input-target data 

   Matrix<double> inputAndTargetData = getInputAndTargetData();

   file << "InputAndTargetData:" << std::endl
        << inputAndTargetData;

   file.close();
}


// void load(char*) method

/// This method loads the members of an input-target data set object from a data file:
///
/// <ul>
/// <li> Number of samples.
/// <li> Number of input variables.
/// <li> Number of target variables
/// <li> Name of input variables. 
/// <li> Name of target variables.
/// <li> Units of input variables. 
/// <li> Units of target variables.
/// <li> Description of input variables. 
/// <li> Description of target variables.
/// <li> Input data.
/// <li> Target data.
/// </ul> 
///
/// Please mind about the file format. This is specified in the User's Guide.
///
/// @param filename Filename.
///
/// @see save(char*).

void InputTargetDataSet::load(char* filename)
{
   std::fstream file; 

   file.open(filename, std::ios::in);

   if(!file.is_open())
   {
      std::cerr << std::endl 
                << "Flood Error: InputTargetDataSet class." << std::endl
                << "void load(char*) method." << std::endl
                << "Cannot open input-target data set file."  << std::endl
                << std::endl;

      exit(1);
   }
   else
   {
      if(display)
      {
         std::cout << std::endl
                   << "Loading input-target data set from file..." << std::endl;
      }
   }

   std::string word;

   while(word != "NumberOfSamples:")
   {
      file >> word;
   }

   file >> numberOfSamples;

   // Number of input variables

   file >> word;
   
   if(word != "NumberOfInputVariables:")
   {
      std::cerr << std::endl 
                << "Flood Error: InputTargetDataSet class." << std::endl
                << "void load(char*) method." << std::endl
                << "Unknown file format."  << std::endl
                << std::endl;

      exit(1);   
   }

   file >> numberOfInputVariables;

   // Number of target variables 

   file >> word;

   if(word != "NumberOfTargetVariables:")
   {
      std::cerr << std::endl 
                << "Flood Error: InputTargetDataSet class." << std::endl
                << "void load(char*) method." << std::endl
                << "Unknown file format."  << std::endl
                << std::endl;

      exit(1);   
   }

   file >> numberOfTargetVariables;

   // Set input-target data set

   setInputTargetDataSet(numberOfSamples, numberOfInputVariables, numberOfTargetVariables);

   // Name of input variables

   file >> word;

   if(word != "NameOfInputVariables:")
   {
      std::cerr << std::endl 
                << "Flood Error: InputTargetDataSet class." << std::endl
                << "void load(char*) method." << std::endl
                << "Unknown file format."  << std::endl
                << std::endl;

      exit(1);   
   }

   file >> nameOfInputVariables;

   // Name of target variables

   file >> word;

   if(word != "NameOfTargetVariables:")
   {
      std::cerr << std::endl 
                << "Flood Error: InputTargetDataSet class." << std::endl
                << "void load(char*) method." << std::endl
                << "Unknown file format."  << std::endl
                << std::endl;

      exit(1);   
   }

   file >> nameOfTargetVariables;

   // Units of input variables

   file >> word;

   if(word != "UnitsOfInputVariables:")
   {
      std::cerr << std::endl 
                << "Flood Error: InputTargetDataSet class." << std::endl
                << "void load(char*) method." << std::endl
                << "Unknown file format."  << std::endl
                << std::endl;

      exit(1);   
   }

   file >> unitsOfInputVariables;

   // Units of target variables

   file >> word;

   if(word != "UnitsOfTargetVariables:")
   {
      std::cerr << std::endl 
                << "Flood Error: InputTargetDataSet class." << std::endl
                << "void load(char*) method." << std::endl
                << "Unknown file format."  << std::endl
                << std::endl;

      exit(1);   
   }

   file >> unitsOfTargetVariables;

   // Description of input variables

   file >> word;

   if(word != "DescriptionOfInputVariables:")
   {
      std::cerr << std::endl 
                << "Flood Error: InputTargetDataSet class." << std::endl
                << "void load(char*) method." << std::endl
                << "Unknown file format."  << std::endl
                << std::endl;

      exit(1);   
   }

   file >> descriptionOfInputVariables;

   // Description of target variables

   file >> word;

   if(word != "DescriptionOfTargetVariables:")
   {
      std::cerr << std::endl 
                << "Flood Error: InputTargetDataSet class." << std::endl
                << "void load(char*) method." << std::endl
                << "Unknown file format."  << std::endl
                << std::endl;

      exit(1);   
   }

   file >> descriptionOfTargetVariables;

   // Input-target data 

   file >> word;

   for(int i = 0; i < numberOfSamples; i++) 
   {
      for(int j = 0; j < numberOfInputVariables; j++)
      {
         file >> inputData[i][j];
      }

      for(int j = numberOfInputVariables; j < numberOfInputVariables + numberOfTargetVariables; j++)
      {
         file >> targetData[i][j-numberOfInputVariables];
      }
   }
 
   file.close(); 
}


// void print() method

/// This method prints to the screen the members of an input-target data set object:
///
/// <ul>
/// <li> Number of samples.
/// <li> Number of input variables.
/// <li> Number of target variables
/// <li> Name of input variables. 
/// <li> Name of target variables.
/// <li> Description of input variables. 
/// <li> Description of target variables.
/// <li> Input data.
/// <li> Target data.
/// </ul> 

void InputTargetDataSet::print(void)
{
   std::cout << std::endl
             << "Flood Neural Network. Input-target Data Set Object." << std::endl;

   // Print number of samples and number of input and target variables

   std::cout << "Number of samples:" << std::endl
             << numberOfSamples << std::endl
             << "Number of input variables:" << std::endl
             << numberOfInputVariables << std::endl
             << "Number of target variables:" << std::endl
             << numberOfTargetVariables << std::endl;

   // Print name of input variables

   std::cout << "Name of input variables:" << std::endl
             << nameOfInputVariables << std::endl;

   // Print name of target variables

   std::cout << "Name of target variables:" << std::endl
             << nameOfTargetVariables << std::endl;

   // Print units of input variables

   std::cout << "Units of input variables:" << std::endl
             << unitsOfInputVariables << std::endl;

   // Print units of target variables

   std::cout << "Units of target variables:" << std::endl
             << unitsOfTargetVariables << std::endl;

   // Print description of input variables

   std::cout << "Description of input variables:" << std::endl
             << descriptionOfInputVariables << std::endl;

   // Print description of target variables

   std::cout << "Description of target variables:" << std::endl
             << descriptionOfTargetVariables << std::endl;

   // Print input-target data 

   Matrix<double> inputAndTargetData = getInputAndTargetData();

   std::cout << "Input and target data:" << std::endl
             << inputAndTargetData;
}


// void saveAllStatistics(char*) method

/// This saves some basic statistics of the input-target data set to a data file.
/// This includes:
///
/// <ul>
/// <li> Mean and standard deviation of input data.
/// <li> Mean and standard deviation of target data.
/// <li> Minimum and maximum of input data.
/// <li> Minimum and maximum of target data.
/// </ul> 

void InputTargetDataSet::saveAllStatistics(char* filename)
{
   std::fstream file; 

   file.open(filename, std::ios::out);

   if(!file.is_open())
   {
      std::cerr << std::endl
                << "Flood Error: InputTargetDataSet class." << std::endl 
                << "void saveAllStatistics(char*) method." << std::endl
                << "Cannot open statistics data file."  << std::endl
                << std::endl;

      exit(1);
   }
   else
   {
      if(display)
      {
         std::cout << std::endl
                   << "Saving statistics to data file..." << std::endl;
      }
   }
  
   // Write file header

   file << "% Flood Neural Network. Input-Target Data Set Statistics." << std::endl;

   file << "Name of input variables:" << std::endl 
        << nameOfInputVariables << std::endl;

   file << "Units of input variables:" << std::endl 
        << unitsOfInputVariables << std::endl;

   file << "Description of input variables:" << std::endl 
        << descriptionOfInputVariables << std::endl;

   Vector< Vector<double> > statistics = calculateAllStatistics();

   file << "Mean of input data:" << std::endl
        << statistics[0] << std::endl
        << "Standard deviation of input data:" << std::endl
        << statistics[1] << std::endl 
        << "Mean of target data:"<< std::endl
        << statistics[2] << std::endl
        << "Standard deviation of target data:" << std::endl
        << statistics[3] << std::endl
        << "Minimum of input data:" << std::endl 
        << statistics[4] << std::endl
        << "Maximum of input data:" << std::endl
        << statistics[5] << std::endl
        << "Minimum of target data:" << std::endl
        << statistics[6]  << std::endl
        << "Maximum of target data:"  << std::endl
        << statistics[7] << std::endl;

   // Close file

   file.close();
}


// void printAllStatistics(void) method

void InputTargetDataSet::printAllStatistics(void)
{
   Vector< Vector<double> > statistics = calculateAllStatistics();

   std::cout << std::endl
             << "Flood Neural Network. Input-Target Data Set Statistics." << std::endl;

   std::cout << "Mean of input data:" << std::endl
             << statistics[0] << std::endl
             << "Standard deviation of input data:" << std::endl
             << statistics[1] << std::endl 
             << "Mean of target data:"<< std::endl
             << statistics[2] << std::endl
             << "Standard deviation of target data:" << std::endl
             << statistics[3] << std::endl
             << "Minimum of input data:" << std::endl 
             << statistics[4] << std::endl
             << "Maximum of input data:" << std::endl
             << statistics[5] << std::endl
             << "Minimum of target data:" << std::endl
             << statistics[6]  << std::endl
             << "Maximum of target data:"  << std::endl
             << statistics[7] << std::endl;
}


// Vector<InputTargetDataSet> splitTrainingValidationAndTesting(double, double, double) method

/// This method returns tree InputTargetDataSet objects, a training set with 50% 
/// of the data, a validation set with 25% of the data and a testing set with 25% 
/// of the data. 

Vector<InputTargetDataSet> InputTargetDataSet::splitTrainingValidationAndTesting(double percentageOfTrainingSamples, 
double percentageOfValidationSamples, double percentageOfTestingSamples)
{
   // Control sentence 

   if(percentageOfTrainingSamples + percentageOfValidationSamples + percentageOfTestingSamples != 100.0)
   {
      std::cerr << std::endl
                << "Flood Error: InputTargetDataSet class." << std::endl 
                << "void splitTrainingValidationAndTesting(double, double, double) method." << std::endl
                << "Sum of training, validation and testing percentages is not 100."  << std::endl
                << std::endl;

      exit(1);
   }


   Vector<InputTargetDataSet> trainingValidationAndTesting(3);

   // Get number of samples for training, validation and testing

   int numberOfTestingSamples = (int)(percentageOfTestingSamples*numberOfSamples/100.0);
   int numberOfValidationSamples = (int)(percentageOfValidationSamples*numberOfSamples/100.0);
   int numberOfTrainingSamples = numberOfSamples - numberOfValidationSamples - numberOfTestingSamples;
   
   // Construct training, validation and testing data set objects

   InputTargetDataSet testingDataSet(numberOfTestingSamples, numberOfInputVariables, numberOfTargetVariables);
   InputTargetDataSet validationDataSet(numberOfValidationSamples, numberOfInputVariables, numberOfTargetVariables);
   InputTargetDataSet trainingDataSet(numberOfTrainingSamples, numberOfInputVariables, numberOfTargetVariables);

   Vector<double> sample(numberOfInputVariables+numberOfTargetVariables, 0.0); 

   double random = 0.0;
   int subset = 0;

   int countTrainingSamples = 0;
   int countValidationSamples = 0;
   int countTestingSamples = 0;

   for(int i = 0; i < numberOfSamples; i++)
   {
      random = rand()/(RAND_MAX+1.0);
      subset = (int)(3*random);

      sample = getSample(i); 
      
      if(subset == 0)  
      {
         if(countTrainingSamples != numberOfTrainingSamples)
         {
            trainingDataSet.setSample(countTrainingSamples, sample);
            countTrainingSamples++;   
         }
         else if(countValidationSamples != numberOfValidationSamples)
         {
            validationDataSet.setSample(countValidationSamples, sample);
            countValidationSamples++;   
         }
         else if(countTestingSamples != numberOfTestingSamples)
         {
            testingDataSet.setSample(countTestingSamples, sample);
            countTestingSamples++;   
         }
      }
      else if(subset == 1) 
      {
         if(countValidationSamples != numberOfValidationSamples)
         {
            validationDataSet.setSample(countValidationSamples, sample);
            countValidationSamples++;   
         }
         else if(countTestingSamples != numberOfTestingSamples)
         {
            testingDataSet.setSample(countTestingSamples, sample);
            countTestingSamples++;   
         }
         else if(countTrainingSamples != numberOfTrainingSamples)
         {
            trainingDataSet.setSample(countTrainingSamples, sample);
            countTrainingSamples++;   
         }
      }
      else if(subset == 2) 
      {
         if(countTestingSamples != numberOfTestingSamples)
         {
            testingDataSet.setSample(countTestingSamples, sample);
            countTestingSamples++;   

         }
         else if(countTrainingSamples != numberOfTrainingSamples)
         {
            trainingDataSet.setSample(countTrainingSamples, sample);
            countTrainingSamples++;   
         }
         else if(countValidationSamples != numberOfValidationSamples)
         {
            validationDataSet.setSample(countValidationSamples, sample);
            countValidationSamples++;   
         }
      }
      else
      {
         std::cerr << std::endl
                   << "Flood Error: InputTargetDataSet class." << std::endl 
                   << "void splitTrainingValidationAndTesting(double, double, double) method." << std::endl
                   << "Subset is neiter 0, 1 or 2." << std::endl
                   << std::endl;

         exit(1);
      }
   }

   // Control sentence
 
   if(countTrainingSamples != numberOfTrainingSamples 
   || countValidationSamples != numberOfValidationSamples 
   || countTestingSamples != numberOfTestingSamples)
   {
      std::cerr << std::endl
                << "Flood Error: InputTargetDataSet class." << std::endl 
                << "void splitTrainingValidationAndTesting(double, double, double) method." << std::endl
                << "Count training validation and testing samples is wrong." << std::endl
                << std::endl;

      exit(1);
   }

   trainingValidationAndTesting[0] = trainingDataSet;
   trainingValidationAndTesting[1] = validationDataSet;
   trainingValidationAndTesting[2] = testingDataSet;

   return(trainingValidationAndTesting);
}


Vector<InputTargetDataSet> InputTargetDataSet::split(double percentage0, double percentage1)
{
   // Control sentence 

   if(percentage0 + percentage1 != 100.0)
   {
      std::cerr << std::endl
                << "Flood Error: InputTargetDataSet class." << std::endl 
                << "void split(double, double) method." << std::endl
                << "Sum of percentages is not 100."  << std::endl
                << std::endl;

      exit(1);
   }

   Vector<InputTargetDataSet> inputTargetDataSetVector(2);

   // Get number of samples for training, validation and testing

   int numberOfSamples0 = (int)(percentage0*numberOfSamples/100.0);
   int numberOfSamples1 = numberOfSamples - numberOfSamples0;
   
   // Construct training, validation and testing data set objects

   inputTargetDataSetVector[0].setInputTargetDataSet(numberOfSamples0, numberOfInputVariables, numberOfTargetVariables);
   inputTargetDataSetVector[1].setInputTargetDataSet(numberOfSamples1, numberOfInputVariables, numberOfTargetVariables);

   Vector<double> sample(numberOfInputVariables+numberOfTargetVariables, 0.0); 

   double random = 0.0;
   int subset = 0;

   int countSamples0 = 0;
   int countSamples1 = 0;

   for(int i = 0; i < numberOfSamples; i++)
   {
      random = rand()/(RAND_MAX+1.0);
      subset = (int)(2*random);

      sample = getSample(i); 
      
      if(subset == 0)  
      {
         if(countSamples0 != numberOfSamples0)
         {
            inputTargetDataSetVector[0].setSample(countSamples0, sample);
            countSamples0++;   
         }
         else if(countSamples1 != numberOfSamples1)
         {
            inputTargetDataSetVector[1].setSample(countSamples1, sample);
            countSamples1++;   
         }
      }
      else if(subset == 1) 
      {
         if(countSamples1 != numberOfSamples1)
         {
            inputTargetDataSetVector[1].setSample(countSamples1, sample);
            countSamples1++;   
         }
         else if(countSamples0 != numberOfSamples0)
         {
            inputTargetDataSetVector[0].setSample(countSamples0, sample);
            countSamples0++;   
         }
      }
      else
      {
         std::cerr << std::endl
                   << "Flood Error: InputTargetDataSet class." << std::endl 
                   << "void splitTrainingValidationAndTesting(double, double, double) method." << std::endl
                   << "Subset is neiter 0 or 1." << std::endl
                   << std::endl;

         exit(1);
      }
   }

   // Control sentence
 
   if(countSamples0 != numberOfSamples0 || countSamples1 != numberOfSamples1)
   {
      std::cerr << std::endl
                << "Flood Error: InputTargetDataSet class." << std::endl 
                << "void split(double, double) method." << std::endl
                << "Count samples is wrong." << std::endl
                << std::endl;

      exit(1);
   }

   return(inputTargetDataSetVector);
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
