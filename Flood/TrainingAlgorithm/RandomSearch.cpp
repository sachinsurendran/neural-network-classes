/****************************************************************************************************************/
/*                                                                                                              */
/*   Flood: An Open Source Neural Networks C++ Library                                                          */
/*   www.cimne.com/flood                                                                                        */
/*                                                                                                              */
/*   R A N D O M   S E A R C H   C L A S S                                                                      */
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
#include <algorithm>
#include <functional>
#include <math.h>
#include <time.h>

#include "RandomSearch.h"

namespace Flood
{

// GENERAL CONSTRUCTOR 

/// General constructor. It creates a random search training algorithm object associated to an objective 
/// functional object. 
/// It also initializes the class members to their default values:
///
/// Training parameters:
/// <ul> 
/// <li> Minimum value of free parameters: -1.0.
/// <li> Maximum value of free parameters:  1.0.
/// </ul> 
///
/// Stopping criteria:
/// <ul> 
/// <li> Evaluation goal: -1.0e99.
/// <li> Maximum time: 1.0e6.
/// <li> Maximum number of epochs: 100. 
/// </ul> 
///  
/// User stuff: 
/// <ul>
/// <li> Display: true.
/// <li> Display period: 25.
/// </ul>
///
/// @param newObjectiveFunctional Pointer to an objective functional object.

RandomSearch::RandomSearch(ObjectiveFunctional* newObjectiveFunctional)
: TrainingAlgorithm(newObjectiveFunctional)
{   
   // Multilayer perceptron

   MultilayerPerceptron* multilayerPerceptron = objectiveFunctional->getMultilayerPerceptron();

   int numberOfFreeParameters = multilayerPerceptron->getNumberOfFreeParameters();

   // Init minimum value of free parameters to -1

   Vector<double> newMinimumValueOfFreeParameters(numberOfFreeParameters, -1.0);

   minimumValueOfFreeParameters = newMinimumValueOfFreeParameters;

   // Init maximum value of free parameters to +1

   Vector<double> newMaximumValueOfFreeParameters(numberOfFreeParameters, 1.0);

   maximumValueOfFreeParameters = newMaximumValueOfFreeParameters;    
   
   // Stopping criteria

   evaluationGoal = -1.0e99;
   maximumTime = 1.0e6;
   maximumNumberOfEpochs = 100; 
   
   // User stuff
   
   displayPeriod = 100;
}


// DEFAULT CONSTRUCTOR

/// Default constructor. It creates a random search training algorithm object not associated to any objective 
/// functional object. 
/// It also initializes the class members to their default values:
///
/// Stopping criteria:
/// <ul> 
/// <li> Evaluation goal: -1.0e99.
/// <li> Maximum training time: 1.0e6.
/// <li> Maximum number of epochs: 100. 
/// </ul> 
///  
/// User stuff: 
/// <ul>
/// <li> Display: true. 
/// <li> Display period: 10. 
/// </ul>

RandomSearch::RandomSearch(void) : TrainingAlgorithm()
{
   // Stopping criteria

   evaluationGoal = -1.0e99;
   maximumTime = 1.0e6;
   maximumNumberOfEpochs = 100; 
   
   // User stuff
   
   displayPeriod = 100;
}


// DESTRUCTOR

/// Destructor.

RandomSearch::~RandomSearch(void)
{

}


// Vector<double> getMinimumValueOfFreeParameters(void) method

/// This method returns the minimum values for all the free parameters in the search hypercube.

Vector<double> RandomSearch::getMinimumValueOfFreeParameters(void)
{
   return(minimumValueOfFreeParameters);
}

// Vector<double> getMaximumValueOfFreeParameters(void) method

/// This method returns the maximum values for all the free parameters in the search hypercube.

Vector<double> RandomSearch::getMaximumValueOfFreeParameters(void)
{
   return(maximumValueOfFreeParameters);
}

// Matrix<double> getMinimumAndMaximumValuesOfFreeParameters(void) method

/// This method returns a single matrix containing the maximum and minimum values for all the free parameters in 
/// the search hypercube.
/// The first row contains the minimum value of the free parameters.
/// The first row contains the maximum value of the free parameters.
///
/// @see getMinimumValueOfFreeParameters(void).
/// @see getMaximumValueOfFreeParameters(void).

Matrix<double> RandomSearch::getMinimumAndMaximumValuesOfFreeParameters(void)
{
   // Multilayer perceptron

   MultilayerPerceptron* multilayerPerceptron = objectiveFunctional->getMultilayerPerceptron();

   int numberOfFreeParameters = multilayerPerceptron->getNumberOfFreeParameters();

   Matrix<double> minimumAndMaximumValuesOfFreeParameters(2, numberOfFreeParameters);

   minimumAndMaximumValuesOfFreeParameters.setRow(0, minimumValueOfFreeParameters); 
   minimumAndMaximumValuesOfFreeParameters.setRow(1, maximumValueOfFreeParameters); 

   return(minimumAndMaximumValuesOfFreeParameters);
}


// void setMinimumValueOfFreeParameters(Vector<double>) method

/// This method sets the minimum values for all the free parameters in the search hypercube.
///
/// @param newMinimumValueOfFreeParameters 

void RandomSearch::setMinimumValueOfFreeParameters(Vector<double> newMinimumValueOfFreeParameters)
{
   // Control sentence (if debug)

   #ifndef NDEBUG 

   int size = newMinimumValueOfFreeParameters.getSize();

   MultilayerPerceptron* multilayerPerceptron = objectiveFunctional->getMultilayerPerceptron();

   int numberOfFreeParameters = multilayerPerceptron->getNumberOfFreeParameters();

   if(size != numberOfFreeParameters) 
   {
      std::cerr << std::endl
                << "Flood Error: RandomSearch class." << std::endl 
                << "void setMinimumValueOfFreeParameters(Vector<double>) method." << std::endl
                << "Size must be equal to number of free parameters." << std::endl
                << std::endl;

      exit(1);
   }

   // Check that minimum value of free parameters is not greater than their maximum value
      
   for(int i = 0; i < numberOfFreeParameters; i++)
   {
      if(newMinimumValueOfFreeParameters[i] > maximumValueOfFreeParameters[i])
      {
         std::cerr << std::endl
                   << "Flood Error: RandomSearch class." << std::endl
                   << "void setMinimumValueOfFreeParameters(Vector<double>) method." << std::endl
                   << "Minimum value of free parameter " << i 
                   << " is greater than maximum value of that parameter." << std::endl 
                   << std::endl;

         exit(1);
      }
   }

   #endif

   // Set minimum value of free parameters

   minimumValueOfFreeParameters = newMinimumValueOfFreeParameters;
}


// void setMaximumValueOfFreeParameters(Vector<double>) method

/// This method sets the maximum values for all the free parameters in the search hypercube.
///
/// @param newMaximumValueOfFreeParameters 

void RandomSearch::setMaximumValueOfFreeParameters(Vector<double> newMaximumValueOfFreeParameters)
{
   // Control sentence (if debug)

   #ifndef NDEBUG 

   int size = newMaximumValueOfFreeParameters.getSize();

   MultilayerPerceptron* multilayerPerceptron = objectiveFunctional->getMultilayerPerceptron();

   int numberOfFreeParameters = multilayerPerceptron->getNumberOfFreeParameters();

   if(size != numberOfFreeParameters) 
   {
      std::cerr << std::endl
                << "Flood Error: RandomSearch class." << std::endl 
                << "void setMaximumValueOfFreeParameters(Vector<double>) method." << std::endl
                << "Size must be equal to number of free parameters." << std::endl
                << std::endl;

      exit(1);
   }
   
   // Check that maximum value of free parameters is not less than their minimum value
      
   for(int i = 0; i < numberOfFreeParameters; i++)
   {
      if(newMaximumValueOfFreeParameters[i] < minimumValueOfFreeParameters[i])
      {
         std::cerr << std::endl
                   << "Flood Error: RandomSearch class." << std::endl
                   << "void setMaximumValueOfFreeParameters(Vector<double>) method." << std::endl
                   << "Maximum value of free parameter " << i 
                   << " is less than minimum value of that parameter." << std::endl 
                   << std::endl;

         exit(1);
      }
   }

   #endif

   // Set maximum value of free parameters

   maximumValueOfFreeParameters = newMaximumValueOfFreeParameters;
}


// void setMinimumAndMaximumValuesOfFreeParameters(Matrix<double>) method

/// This method sets both the minimum and the maximum values for all the free parameters in the search hypercube 
/// from a single matrix.
/// The first row must contain the minimum values for the free parameters.
/// The second row must contain the maximum values for the free parameters.
///
/// @param newMinimumAndMaximumValuesOfFreeParameters Set of minimum and maximum values for the free parameters 
/// of the neural network.
///
/// @see setMinimumValueOfFreeParameters(Vector<double>).
/// @see setMaximumValueOfFreeParameters(Vector<double>).

void RandomSearch
::setMinimumAndMaximumValuesOfFreeParameters(Matrix<double> newMinimumAndMaximumValuesOfFreeParameters)
{
   MultilayerPerceptron* multilayerPerceptron = objectiveFunctional->getMultilayerPerceptron();

   int numberOfFreeParameters = multilayerPerceptron->getNumberOfFreeParameters();

   // Control sentence (if debug)

   #ifndef NDEBUG 

   int numberOfRows = newMinimumAndMaximumValuesOfFreeParameters.getNumberOfRows();
   int numberOfColumns = newMinimumAndMaximumValuesOfFreeParameters.getNumberOfColumns();


   if(numberOfRows != 2 || numberOfColumns != numberOfFreeParameters) 
   {
      std::cerr << std::endl
                << "Flood Error: RandomSearch class." << std::endl
                << "void setMinimumAndMaximumValuesOfFreeParameters(Matrix<double>) method." << std::endl
                << "Number of rows must be two and number of columns must be equal to "
                << "number of free parameters." << std::endl
                << std::endl;

      exit(1);
   }

   #endif

   minimumValueOfFreeParameters = newMinimumAndMaximumValuesOfFreeParameters.getRow(0);
   maximumValueOfFreeParameters = newMinimumAndMaximumValuesOfFreeParameters.getRow(1);
   
   for(int i = 0; i < numberOfFreeParameters; i++)
   {
      if(minimumValueOfFreeParameters[i] > maximumValueOfFreeParameters[i])
      {
         std::cout << std::endl
                   << "Flood Error RandomSearch class." << std::endl
                   << "void setMinimumAndMaximumValuesOfFreeParameters(Matrix<double>) method." << std::endl
                   << "Minimum value must be less or equal than maximum value." << std::endl
                   << std::endl;

         exit(1);
      }
   }   
}


// void setReserveAllTrainingHistory(bool) method

/// This method makes the training history of all variables to reseved or not in memory.
///
/// @param newReserveAllTrainingHistory True if the training history of all variables is to be reserved, 
/// false otherwise.

void RandomSearch::setReserveAllTrainingHistory(bool newReserveAllTrainingHistory)
{
   reserveElapsedTimeHistory = newReserveAllTrainingHistory;
   reserveFreeParametersHistory = newReserveAllTrainingHistory;
   reserveFreeParametersNormHistory = newReserveAllTrainingHistory;
   reserveEvaluationHistory = newReserveAllTrainingHistory;
}


// void train(void) method

/// This method trains a multilayer perceptron with an associated 
/// objective function according to the random search training algorithm.
/// Training occurs according to the training parameters. 

void RandomSearch::train(void)
{
   // Control sentence (if debug)

   #ifndef NDEBUG 

   if(objectiveFunctional == NULL)
   {
      std::cerr << std::endl
                << "Flood Error: RandomSearch class." << std::endl
                << "void train(void) method." << std::endl
                << "Pointer to objective functional object cannot be NULL." << std::endl
                << std::endl;

      exit(1);
   }

   #endif

   time_t beginningTime, currentTime;
   double elapsedTime = 0.0;

   // Set beginning training time 

   time(&beginningTime);

   if(display)
   {
      std::cout << std::endl
                << "Training with random search..." 
                << std::endl;
   }

   // Multilayer perceptron stuff

   MultilayerPerceptron* multilayerPerceptron = objectiveFunctional->getMultilayerPerceptron();

   int numberOfFreeParameters = multilayerPerceptron->getNumberOfFreeParameters();

   // Resize training history

   resizeTrainingHistory(1+maximumNumberOfEpochs);

   Vector<double> freeParameters = multilayerPerceptron->getFreeParameters();

   if(reserveFreeParametersHistory)
   {
      freeParametersHistory.setRow(0, freeParameters);                                
   }

   double freeParametersNorm = freeParameters.calculateNorm();

   if(reserveFreeParametersNormHistory)
   {
      freeParametersNormHistory[0] = freeParametersNorm; 
   }

   if(display && (freeParametersNorm >= warningFreeParametersNorm))
   {
      std::cout << std::endl
                << "Flood Warning: Initial free parameters norm is " << freeParametersNorm << "." << std::endl;          
   }

   Vector<double> potentialFreeParameters(numberOfFreeParameters);

   double potentialFreeParametersNorm = 0.0;

   double potentialEvaluation = 0.0;

   // Initial evaluation
   
   double evaluation = objectiveFunctional->calculateEvaluation();

   if(reserveEvaluationHistory)
   {
      evaluationHistory[0] = evaluation;
   }

   // Elapsed time

   time(&currentTime);
   elapsedTime = difftime(currentTime, beginningTime);

   if(reserveElapsedTimeHistory)
   {
      elapsedTimeHistory[0] = elapsedTime;                             
   }

   if(evaluation <= evaluationGoal)
   {          
      if(display)
      {
         std::cout << std::endl
                   << "Initial evaluation is less than goal." << std::endl;

         std::cout << "Elapsed time: " << elapsedTime << ";" << std::endl;
         std::cout << "Initial parameters norm: " << freeParametersNorm << std::endl;
         std::cout << "Initial evaluation: " << evaluation << std::endl;          

         objectiveFunctional->print();
      }
      
      // Resize training history

      resizeTrainingHistory(1);

      return;
   }
   else
   {
      if(display)
      {
         std::cout << "Initial elapsed time: " << elapsedTime << ";" << std::endl;
         std::cout << "Initial parameters norm: " << freeParametersNorm << std::endl;
         std::cout << "Initial evaluation: " <<  evaluation << std::endl;

         objectiveFunctional->print();
      }
   }

   // Main loop

   for(int epoch = 1; epoch <= maximumNumberOfEpochs; epoch++)
   {
      // Free parameters           
           
      for(int i = 0; i < numberOfFreeParameters; i++)
      {
         double random = (double)rand()/(RAND_MAX+1.0);

         potentialFreeParameters[i] = minimumValueOfFreeParameters[i] 
         + (maximumValueOfFreeParameters[i]- minimumValueOfFreeParameters[i])*random;
      }

      if(reserveFreeParametersHistory)
      {
         freeParametersHistory.setRow(epoch, potentialFreeParameters);                                
      }  

      // Free parameters norm
      
      potentialFreeParametersNorm = potentialFreeParameters.calculateNorm();
      
      if(reserveFreeParametersNormHistory)
      {
         freeParametersNormHistory[epoch] = potentialFreeParametersNorm;
      }  
      
      // Evaluation

      potentialEvaluation 
      = objectiveFunctional->calculatePotentialEvaluation(potentialFreeParameters);

      if(reserveEvaluationHistory)
      {
         evaluationHistory[epoch] = potentialEvaluation;                            
      }

      // Elapsed time

      time(&currentTime);
      elapsedTime = difftime(currentTime, beginningTime);
    
      if(reserveElapsedTimeHistory)
      {
         elapsedTimeHistory[epoch] = elapsedTime;                             
      }

      // Check for evaluation improvement

      if(potentialEvaluation < evaluation)
      {
         evaluation = potentialEvaluation;

         multilayerPerceptron->setFreeParameters(potentialFreeParameters);
         
         freeParametersNorm = potentialFreeParameters.calculateNorm();
      }

      // Stopping Criteria

      // Evaluation goal 

      if(evaluation <= evaluationGoal)
      {
         if(display)
         {
            std::cout << std::endl
                      << "Epoch " << epoch << ": "
                      << "Evaluation goal reached." << std::endl;

            std::cout << "Elapsed time: " << elapsedTime << ";" << std::endl;
            std::cout << "Final parameters norm: " << freeParametersNorm << std::endl;
            std::cout << "Final evaluation: " << evaluation << std::endl;

            objectiveFunctional->print();
         }

         resizeTrainingHistory(1+epoch);
         
         break;
      }

      // Maximum optimization time

      if(elapsedTime >= maximumTime)
      {
         if(display)
         {
            std::cout << std::endl << "Epoch " << epoch << "Maximum training time reached." << std::endl;

            std::cout << "Elapsed time: " << elapsedTime << ";" << std::endl;
            std::cout << "Final parameters norm: " << freeParametersNorm << std::endl;
            std::cout << "Final evaluation: " << evaluation << std::endl;

            objectiveFunctional->print();
         }

         resizeTrainingHistory(1+epoch);
  
         break;
      }

      // Maximum number of epochs

      if(epoch == maximumNumberOfEpochs)
      {
         if(display)
         {
            std::cout << std::endl
                      << "Epoch " << epoch << ": Maximum number of epochs reached." << std::endl;

            std::cout << "Elapsed time: " << elapsedTime << ";" << std::endl;
            std::cout << "Final parameters norm: " << freeParametersNorm << std::endl;
            std::cout << "Final evaluation: " << evaluation << std::endl;

            objectiveFunctional->print();
         }

         break;
      }

      // Progress

      if(display && epoch % displayPeriod == 0)
      {
         std::cout << std::endl
                   << "Epoch " << epoch << ";" << std::endl;

         std::cout << "Elapsed time: " << elapsedTime << ";" << std::endl;
         std::cout << "Free parameters norm: " << freeParametersNorm << std::endl;
         std::cout << "Evaluation: " << evaluation << std::endl;
         
         objectiveFunctional->print();
      }
   }
}


// void print(void) method

/// This method prints to the screen the training parameters, the stopping criteria
/// and other user stuff concerning the random search object:
///
/// Training parameters:
/// <ul> 
/// <li> Minimum value of free parameters.
/// <li> Maximum value of free parameters.
/// </ul> 
///  
/// Stopping criteria:
/// <ul> 
/// <li> Evaluation goal.
/// <li> Maximum time.
/// <li> Maximum number of epochs. 
/// </ul> 
///  
/// User stuff: 
/// <ul>
/// <li> Display. 
/// <li> Display period. 
/// </ul>

void RandomSearch::print(void)
{
   // Multilayer perceptron

   std::cout << std::endl
             << "Random Search Training Algorithm Object." << std::endl;

   // Training parameters

   // Minimum value of free parameters
   
   std::cout << "Minimum value of free parameters:" << std::endl 
             << minimumValueOfFreeParameters << std::endl;

   // Maximum value of free parameters
   
   std::cout << "Maximum value of free parameters:" << std::endl
             << maximumValueOfFreeParameters << std::endl;
   
   // Stopping criteria

   std::cout << "Evaluation goal: " << std::endl
             << evaluationGoal << std::endl
             << "Maximum time: " << std::endl
             << maximumTime << std::endl 
             << "Maximum number of epochs: " << std::endl
             << maximumNumberOfEpochs << std::endl; 

   // User stuff

   std::cout << "Display: " << std::endl
             << display << std::endl
             << "Display period: " << std::endl
             << displayPeriod << std::endl;
}



// void save(char*) method

/// This method saves the random search object to a data file. 
///
/// Training parameters:
/// <ul> 
/// <li> Minimum value of free parameters.
/// <li> Maximum value of free parameters.
/// </ul> 
///  
/// Stopping criteria:
/// <ul> 
/// <li> Evaluation goal.
/// <li> Maximum time.
/// <li> Maximum number of epochs. 
/// </ul> 
///  
/// User stuff: 
/// <ul>
/// <li> Display. 
/// <li> Display period. 
/// </ul>
///
/// @param filename Filename.
///
/// @see load(char*).

void RandomSearch::save(char* filename)
{
   // File

   std::fstream file;

   file.open(filename, std::ios::out);

   if(!file.is_open())
   {
      std::cerr << std::endl 
                << "Flood Error: RandomSearch class." << std::endl
                << "void save(char*) method." << std::endl
                << "Cannot open random search object data file."  << std::endl
                << std::endl;

      exit(1);
   }
   else
   {
      if(display)
      {
         std::cout << std::endl
                   << "Saving random search object to data file..." << std::endl;
      }
   }

   // Write file header

   file << "% Flood Neural Network. Random Seach Object." << std::endl;

   // Training parameters

   // Minimum value of free parameters
   
   file << "MinimumValueOfFreeParameters: " << std::endl
        << minimumValueOfFreeParameters << std::endl;

   // Maximum value of free parameters
   
   file << "MaximumValueOfFreeParameters: " << std::endl
        << maximumValueOfFreeParameters << std::endl;

   // Stopping criteria

   file << "EvaluationGoal:" << std::endl
        << evaluationGoal << std::endl
        << "MaximumTime: " << std::endl
        << maximumTime << std::endl
        << "MaximumNumberOfEpochs: " << std::endl
        << maximumNumberOfEpochs << std::endl;

   // User stuff

   file << "Display: " << std::endl
        << display << std::endl
        << "DisplayPeriod: " << std::endl
        << displayPeriod << std::endl;

   file.close();
}


// void load(char*) method

/// This method loads a random search object from a data file. 
/// Please mind about the file format, wich is specified in the User's Guide. 
///
/// Training parameters:
/// <ul> 
/// <li> Minimum value of free parameters.
/// <li> Maximum value of free parameters.
/// </ul> 
///  
/// Stopping criteria:
/// <ul> 
/// <li> Evaluation goal.
/// <li> Maximum training time.
/// <li> Maximum number of epochs. 
/// </ul> 
///  
/// User stuff: 
/// <ul>
/// <li> Display. 
/// <li> Display period. 
/// </ul>
///
/// @param filename Filename.
///
/// @see save(char*).

void RandomSearch::load(char* filename)
{
   // File

   std::fstream file;

   file.open(filename, std::ios::in);

   if(!file.is_open())
   {
      std::cerr << std::endl
                << "Flood Error: RandomSearch class." << std::endl
                << "void load(char*) method." << std::endl
                << "Cannot open random search object data file."  << std::endl;

      exit(1);
   }
   else
   {
      if(display)
      {
         std::cout << std::endl
                   << "Loading random seach object from data file..."  << std::endl;
      }
   }

   std::string word;

   // Training parameters

   // Minimum value of free parameters

   while(word != "MinimumValueOfFreeParameters:")
   {
      file >> word;
   }

   file >> minimumValueOfFreeParameters;

   // Maximum value of free parameters

   file >> word;

   file >> maximumValueOfFreeParameters;

   // Stopping criteria: 

   // Evalaution goal

   file >> word;

   file >> evaluationGoal;

   // Maximum time

   file >> word;

   file >> maximumTime;

   // Maximum number of epochs

   file >> word;

   file >> maximumNumberOfEpochs;

   // User stuff: 

   // Display

   file >> word;

   file >> display;

   // Display period

   file >> word;

   file >> displayPeriod;

   // Close file

   file.close();
}


// void resizeTrainingHistory(int) method

/// This method resizes the vectors or matrices containing training history information 
/// to a new size:
///
/// <ul>
/// <li> Elapsed time history vector.
/// <li> Free parameters history matrix.
/// <li> Free parameters norm history vector. 
/// <li> Evaluation history vector.
/// </ul>
///
/// @param newSize Size of training history. 

void RandomSearch::resizeTrainingHistory(int newSize)
{
   // Free parameters history matrix

   if(reserveFreeParametersHistory)
   {
      MultilayerPerceptron* multilayerPerceptron = objectiveFunctional->getMultilayerPerceptron();
                                        
      int numberOfFreeParameters = multilayerPerceptron->getNumberOfFreeParameters();
                                   
      freeParametersHistory.resize(newSize, numberOfFreeParameters);
   }

   // Free parameters norm history vector

   if(reserveFreeParametersNormHistory)
   {
      freeParametersNormHistory.resize(newSize);
   }

   // Evaluation history vector

   if(reserveEvaluationHistory)
   {
      evaluationHistory.resize(newSize);
   }

   // Elapsed time history vector

   if(reserveElapsedTimeHistory)
   {
      elapsedTimeHistory.resize(newSize);
   }
}


// void saveTrainingHistory(char*) method

/// This method saves the training history to a data file. 
///
/// @param filename Training history filename. 

void RandomSearch::saveTrainingHistory(char* filename)

{
   std::fstream file; 

   file.open(filename, std::ios::out);

   // Write file header 

   if(!file.is_open())
   {
      std::cerr << std::endl 
                << "Flood Error: RandomSearch class." << std::endl
                << "void saveTrainingHistory(char*) method." << std::endl
                << "Cannot open training history data file." << std::endl
                << std::endl;

      exit(1);
   }
   else
   {
      if(display)
      {
         std::cout << std::endl 
                   << "Saving training history to data file..." << std::endl;
      }
   }

   // Write file header

   file << "% Flood Neural Network. Random Search Training History." << std::endl;

   // Write file data

   if(reserveElapsedTimeHistory)
   {
      file << "ElapsedTimeHistory:" << std::endl;      file << elapsedTimeHistory << std::endl;      
   }
   if(reserveFreeParametersHistory)
   {
      file << "FreeParametersHistory:" << std::endl;      file << freeParametersHistory << std::endl;      
   }
   if(reserveFreeParametersNormHistory)
   {
      file << "FreeParametersNormHistory:" << std::endl;      file << freeParametersHistory << std::endl;      
   }
   if(reserveEvaluationHistory)
   {
      file << "EvaluationHistory:" << std::endl;      file << evaluationHistory << std::endl;      
   }

   file.close();     
     
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
