/****************************************************************************************************************/
/*                                                                                                              */
/*   Flood: An Open Source Neural Networks C++ Library                                                          */
/*   www.cimne.com/flood                                                                                        */
/*                                                                                                              */
/*   I N P U T - T A R G E T   D A T A   S E T   C L A S S   H E A D E R                                        */
/*                                                                                                              */ 
/*   Roberto Lopez                                                                                              */ 
/*   International Center for Numerical Methods in Engineering (CIMNE)                                          */
/*   Technical University of Catalonia (UPC)                                                                    */
/*   Barcelona, Spain                                                                                           */
/*   E-mail: rlopez@cimne.upc.edu                                                                               */ 
/*                                                                                                              */
/****************************************************************************************************************/


#ifndef __INPUTTARGETDATASET_H__
#define __INPUTTARGETDATASET_H__

#include <string>

#include "../Utilities/Vector.h"
#include "../Utilities/Matrix.h"

namespace Flood
{

/// This class represents the concept of input-target data set for data modelling problems, such as function
/// regression, pattern recognition and time series prediction.
///
/// @see SumSquaredError
/// @see MeanSquaredError
/// @see RootMeanSquaredError
/// @see NormalizedSquaredError
/// @see MinkowskiError
/// @see RegularizedMeanSquaredError

class InputTargetDataSet 
{

private: 

   // FIELDS

   /// Number of samples in the data set. 

   int numberOfSamples;

   /// Number of input variables. 

   int numberOfInputVariables;

   /// Number of target variables. 

   int numberOfTargetVariables;

   /// Name of input variables.

   Vector<std::string> nameOfInputVariables;

   /// Name of target variables.

   Vector<std::string> nameOfTargetVariables;

   /// Units of input variables.

   Vector<std::string> unitsOfInputVariables;

   /// Units of target variables.

   Vector<std::string> unitsOfTargetVariables;

   /// Description of input variables.

   Vector<std::string> descriptionOfInputVariables;

   /// Description of target variables.

   Vector<std::string> descriptionOfTargetVariables;

   /// Input data Matrix.

   Matrix<double> inputData;

   /// Target data Matrix.

   Matrix<double> targetData;

   /// Display messages to screen.
   
   bool display;

public:  

   // GENERAL CONSTRUCTOR

   InputTargetDataSet(int, int, int);


   // DEFAULT CONSTRUCTOR

   InputTargetDataSet(void);


   // DESTRUCTOR

   virtual ~InputTargetDataSet();


   // METHODS

   // Get methods

   int getNumberOfSamples(void);

   int getNumberOfInputVariables(void);
   int getNumberOfTargetVariables(void);

   Vector<std::string> getNameOfInputVariables(void);
   Vector<std::string> getNameOfTargetVariables(void);

   std::string getNameOfSingleInputVariable(int);
   std::string getNameOfSingleTargetVariable(int);

   Vector<std::string> getUnitsOfInputVariables(void);
   Vector<std::string> getUnitsOfTargetVariables(void);

   std::string getUnitsOfSingleInputVariable(int);
   std::string getUnitsOfSingleTargetVariable(int);   

   Vector<std::string> getDescriptionOfInputVariables(void);
   Vector<std::string> getDescriptionOfTargetVariables(void);

   std::string getDescriptionOfSingleInputVariable(int);
   std::string getDescriptionOfSingleTargetVariable(int);

   Vector< Vector<std::string> > getAllInformation(void);

   Matrix<double>& getInputData(void);
   Matrix<double>& getTargetData(void);
   Matrix<double> getInputAndTargetData(void);

   Vector<double> getSample(int);
   void setSample(int, Vector<double>);

   void addSample(Vector<double>);
   void subtractSample(int);

   bool getDisplay(void);

   // Set methods

   void setInputTargetDataSet(int, int, int);

   void setNameOfInputVariables(Vector<std::string>);
   void setNameOfTargetVariables(Vector<std::string>);

   void setNameOfSingleInputVariable(int, std::string);
   void setNameOfSingleTargetVariable(int, std::string);

   void setUnitsOfInputVariables(Vector<std::string>);
   void setUnitsOfTargetVariables(Vector<std::string>);

   void setUnitsOfSingleInputVariable(int, std::string);
   void setUnitsOfSingleTargetVariable(int, std::string);

   void setDescriptionOfInputVariables(Vector<std::string>);
   void setDescriptionOfTargetVariables(Vector<std::string>);

   void setDescriptionOfSingleInputVariable(int, std::string);
   void setDescriptionOfSingleTargetVariable(int, std::string);

   void setAllInformation(Vector< Vector<std::string> >);
 
   void setInputData(Matrix<double>);
   void setTargetData(Matrix<double>);
   void setInputAndTargetData(Matrix<double>);

   void setDisplay(bool);
   
   // Mean and standard deviation of input and target data methods

   Matrix<double> calculateMeanAndStandardDeviationOfInputData(void);
   Matrix<double> calculateMeanAndStandardDeviationOfTargetData(void);

   // Minimum and maximum of input and target data methods

   Matrix<double> calculateMinimumAndMaximumOfInputData(void);
   Matrix<double> calculateMinimumAndMaximumOfTargetData(void);

   // All statistics

   Vector< Vector<double> > calculateAllStatistics(void);

   // Postprocess input and target data

   void preprocessMeanAndStandardDeviation(void);
   void preprocessMeanAndStandardDeviation(Matrix<double>, Matrix<double>);

   void preprocessMinimumAndMaximum(void);
   void preprocessMinimumAndMaximum(Matrix<double>, Matrix<double>);

   // Postprocess input and target data

   void postprocessMeanAndStandardDeviation(Matrix<double>, Matrix<double>);
   void postprocessMinimumAndMaximum(Matrix<double>, Matrix<double>);

   // Utility methods

   void save(char*);
   void load(char*);

   void print(void);

   void saveAllStatistics(char*);
   void printAllStatistics(void);

   Vector<InputTargetDataSet> splitTrainingValidationAndTesting(double, double, double);
   Vector<InputTargetDataSet> split(double, double);
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
