/****************************************************************************************************************/
/*                                                                                                              */
/*   Flood: An Open Source Neural Networks C++ Library                                                          */
/*   www.cimne.com/flood                                                                                        */
/*                                                                                                              */
/*   E V O L U T I O N A R Y   A L G O R I T H M   C L A S S                                                    */
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

#include "EvolutionaryAlgorithm.h"

namespace Flood
{

// GENERAL CONSTRUCTOR 
//
/// General constructor. It creates a evolutionary training algorithm object associated to an objective functional
/// object.
/// It also initializes the class members to their default values:
///
/// Training operators:
/// <ul>
/// <li> Fitness assignment method: Linear ranking.
/// <li> Selection method: Stochastic universal sampling.
/// <li> Recombination method: Intermediate.
/// <li> Mutation method: Normal.
/// </ul>
///
/// Training parameters:
/// <ul>
/// <li> Population size: 10 * [number of free parameters].
/// <li> Selective pressure: 1.5.
/// <li> Recombination size: 0.25.
/// <li> Mutation rate: = 1.0 / [number of free parameters].
/// <li> Mutation range: = 0.1
/// </ul>
///
/// Stopping criteria:
/// <ul> 
/// <li> Evaluation goal: -1.0e99.
/// <li> Maximum training time: 1.0e6.
/// <li> Maximum number of generations: 100. 
/// </ul> 
///  
/// User stuff: 
/// <ul>
/// <li> Display: 1. 
/// <li> Display period: 1. 
/// </ul>
///
/// @param newObjectiveFunctional Pointer to an objective functional object.

EvolutionaryAlgorithm::EvolutionaryAlgorithm(ObjectiveFunctional* newObjectiveFunctional)
: TrainingAlgorithm(newObjectiveFunctional)
{
   // Multilayer perceptron

   MultilayerPerceptron* multilayerPerceptron = objectiveFunctional->getMultilayerPerceptron();

   int numberOfFreeParameters = multilayerPerceptron->getNumberOfFreeParameters();

   // Fitness assignment method

   fitnessAssignmentMethod = LinearRanking;

   // Selection method

   selectionMethod = StochasticUniversalSampling;

   // Recombination method

   recombinationMethod = Intermediate;

   // Mutation method

   mutationMethod = Normal;

   // Training parameters

   populationSize = 10*numberOfFreeParameters;

   selectivePressure = 1.5;

   recombinationSize = 0.25;

   mutationRate = 1.0/(double)numberOfFreeParameters;
   mutationRange = 0.1;

   // Stopping criteria

   evaluationGoal = -1.0e99;
   maximumTime = 1.0e12;

   maximumNumberOfGenerations = 100;

   displayPeriod = 1;

   // Population matrix

   population.resize(populationSize, numberOfFreeParameters);

   initPopulationNormal();

   // Evaluation vector

   evaluation.resize(populationSize);

   // Fitness vector

   fitness.resize(populationSize);

   // Selection vector

   selection.resize(populationSize);

   reservePopulationHistory = false;
   reserveMeanNormHistory = false;
   reserveStandardDeviationNormHistory = false;
   reserveBestNormHistory = false;
   reserveMeanEvaluationHistory = false;
   reserveStandardDeviationEvaluationHistory = false;
   reserveBestEvaluationHistory = false;
}


// DEFAULT CONSTRUCTOR

/// Default constructor. It creates a evolutionary training algorithm object not associated to any objective 
/// functional object.
/// It also initializes the class members to their default values:
///
/// Training operators:
/// <ul>
/// <li> Fitness assignment method: Linear ranking.
/// <li> Selection method: Stochastic universal sampling.
/// <li> Recombination method: Intermediate.
/// <li> Mutation method: Normal.
/// </ul>
///
/// Training parameters:
/// <ul>
/// <li> Population size: 0.
/// <li> Selective pressure: 1.5.
/// <li> Recombination size: 0.25.
/// <li> Mutation rate: = 0.1.
/// <li> Mutation range: = 0.1
/// </ul>
///
/// Stopping criteria:
/// <ul> 
/// <li> Evaluation goal: -1.0e99.
/// <li> Maximum training time: 1.0e6.
/// <li> Maximum number of generations: 100. 
/// </ul> 
///  
/// User stuff: 
/// <ul>
/// <li> Display: 1. 
/// <li> Display period: 1. 
/// </ul>

EvolutionaryAlgorithm::EvolutionaryAlgorithm(void) : TrainingAlgorithm()
{
   // Fitness assignment method

   fitnessAssignmentMethod = LinearRanking;

   // Selection method

   selectionMethod = StochasticUniversalSampling;

   // Recombination method

   recombinationMethod = Intermediate;

   // Mutation method

   mutationMethod = Normal;

   // Training parameters

   populationSize = 0;

   selectivePressure = 1.5;

   recombinationSize = 0.25;

   mutationRate = 0.1;
   mutationRange = 0.1;

   // Stopping criteria

   evaluationGoal = -1.0e99;
   maximumTime = 1.0e12;

   maximumNumberOfGenerations = 100;

   displayPeriod = 1;

   reservePopulationHistory = false;
   reserveMeanNormHistory = false;
   reserveStandardDeviationNormHistory = false;
   reserveBestNormHistory = false;
   reserveMeanEvaluationHistory = false;
   reserveStandardDeviationEvaluationHistory = false;
   reserveBestEvaluationHistory = false;
}


// DESTRUCTOR

/// Destructor.

EvolutionaryAlgorithm::~EvolutionaryAlgorithm(void)
{

}


// METHODS

// int getPopulationSize(void) method

/// This method returns the number of individuals in the population.

int EvolutionaryAlgorithm::getPopulationSize(void)
{
   return(populationSize);
}


// Vector<double> getMeanEvaluationHistory(void) method

/// This method returns a history with the mean evaluation of the population during training.
///
/// @see getBestEvaluationHistory(void).
/// @see getStandardDeviationEvaluationHistory(void).

Vector<double> EvolutionaryAlgorithm::getMeanEvaluationHistory(void)
{
   return(meanEvaluationHistory);
}


// Vector<double> getStandardDeviationEvaluationHistory(void) method

/// This method returns a history with the standard deviation of the population evaluation during training.
///
/// @see getBestEvaluationHistory(void).
/// @see getMeanEvaluationHistory(void).

Vector<double> EvolutionaryAlgorithm::getStandardDeviationEvaluationHistory(void)
{
   return(standardDeviationEvaluationHistory);
}


// Vector<double> getBestEvaluationHistory(void) method

/// This method returns a history with the evaluation value of the best individual ever during training.
///
/// @see getMeanEvaluationHistory(void).
/// @see getStandardDeviationEvaluationHistory(void).

Vector<double> EvolutionaryAlgorithm::getBestEvaluationHistory(void)
{
   return(bestEvaluationHistory);
}


// void setPopulationSize(int) method
//
/// This method sets a new population with a new number of individuals.  
/// The new population size must be an even number equal or greater than four. 
///
/// @param newPopulationSize Number of individuals in the population. This must be an even number equal or 
/// greater than four. 

void EvolutionaryAlgorithm::setPopulationSize(int newPopulationSize)
{
   MultilayerPerceptron* multilayerPerceptron = objectiveFunctional->getMultilayerPerceptron();

   int numberOfFreeParameters = multilayerPerceptron->getNumberOfFreeParameters();

   if(newPopulationSize < 4)
   {
      std::cerr << std::endl
                << "Flood Error: EvolutionaryAlgorithm class." << std::endl
                << "void setPopulationSize(int) method." << std::endl
                << "New population size must be equal or greater than 4." << std::endl
                << std::endl;

      exit(1);
   }
   else if(newPopulationSize%2 != 0)
   {
      std::cerr << std::endl
                << "Flood Error: EvolutionaryAlgorithm class." << std::endl
                << "void setPopulationSize(int) method."
                << std::endl
                << "New population size is not divisible by 2." << std::endl
                << std::endl;

      exit(1);
   }
   else
   {
      populationSize = newPopulationSize;

      // Set population matrix

      population.resize(populationSize, numberOfFreeParameters);

      initPopulationNormal();

      // Set evaluation vector

      evaluation.resize(populationSize);

      // Set fitness vector

      fitness.resize(populationSize);

      // Set selection vector

      selection.resize(populationSize);
   }
}


// Population stuff

// Matrix<double> getPopulation(void) method

/// This method returns the population Matrix.

Matrix<double> EvolutionaryAlgorithm::getPopulation(void)
{
   return(population);
}


// Vector<double> getEvaluation(void) method

/// This method returns the actual evaluation value of all individuals in the population.
///
/// @see getFitness(void).
/// @see getSelection(void).

Vector<double> EvolutionaryAlgorithm::getEvaluation(void)
{
   return(evaluation);
}


// Vector<double> getFitness(void) method

/// This method returns the actual fitness value of all individuals in the population.
///
/// @see calculateEvaluation(void).
/// @see getSelection(void).

Vector<double> EvolutionaryAlgorithm::getFitness(void)
{
   return(fitness);
}


// Vector<bool> getSelection(void) method

/// This method returns the actual selection value of all individuals in the population.
///
/// @see calculateEvaluation(void).
/// @see getFitness(void).

Vector<bool> EvolutionaryAlgorithm::getSelection(void)
{
   return(selection);
}


// bool getReservePopulationHistory(void) method

/// This method returns true if the population history vector of matrices is to be reserved, and false otherwise.

bool EvolutionaryAlgorithm::getReservePopulationHistory(void)
{
   return(reservePopulationHistory);
}


// bool getReserveMeanNormHistory(void) method

/// This method returns true if the mean population norm history vector is to be reserved, and false otherwise.

bool EvolutionaryAlgorithm::getReserveMeanNormHistory(void)
{
   return(reserveMeanNormHistory);
}


// bool getReserveStandardDeviationNormHistory(void) method

/// This method returns true if the standard deviation of the population norm history vector is to be reserved,
/// and false otherwise.

bool EvolutionaryAlgorithm::getReserveStandardDeviationNormHistory(void)
{
   return(reserveStandardDeviationNormHistory);
}


// bool getReserveBestNormHistory(void) method

/// This method returns true if the norm of the best individual in the population history vector is to be 
/// reserved, and false otherwise.

bool EvolutionaryAlgorithm::getReserveBestNormHistory(void)
{
   return(reserveBestNormHistory);
}


// bool getReserveMeanEvaluationHistory(void) method

/// This method returns true if the mean evaluation history vector is to be reserved, and false otherwise.

bool EvolutionaryAlgorithm::getReserveMeanEvaluationHistory(void)
{
   return(reserveMeanEvaluationHistory);
}


// bool getReserveStandardDeviationEvaluationHistory(void) method

/// This method returns true if the standard deviation of the evaluation history vector is to be reserved,
/// and false otherwise.

bool EvolutionaryAlgorithm::getReserveStandardDeviationEvaluationHistory(void)
{
   return(reserveStandardDeviationEvaluationHistory);
}


// bool getReserveBestEvaluationHistory(void) method

/// This method returns true if the best evaluation history vector is to be reserved, and false otherwise.

bool EvolutionaryAlgorithm::getReserveBestEvaluationHistory(void)
{
   return(reserveBestEvaluationHistory);
}


// Vector< Matrix<double> > getPopulationHistory(void) method

/// This method returns the population history over the training epochs, which is a vector of matrices. 

Vector< Matrix<double> > EvolutionaryAlgorithm::getPopulationHistory(void)
{
   return(populationHistory);
}


// Vector<double> getMeanNormHistory(void) method

/// This method returns the mean norm history.

Vector<double> EvolutionaryAlgorithm::getMeanNormHistory(void)
{
   return(meanNormHistory);
}

// Vector<double> getStandardDeviationNormHistory(void) method

/// This method returns the standard deviation norm history.

Vector<double> EvolutionaryAlgorithm::getStandardDeviationNormHistory(void)
{
   return(standardDeviationNormHistory);
}


// Vector<double> getBestNormHistory(void) method

/// This method returns the best norm history.

Vector<double> EvolutionaryAlgorithm::getBestNormHistory(void)
{
   return(bestNormHistory);
}


// void setPopulation(Matrix<double>) method

/// This method sets a new population.
///
/// @param newPopulation Population Matrix.

void EvolutionaryAlgorithm::setPopulation(Matrix<double> newPopulation)
{
   // Control sentence (if debug)

   #ifndef NDEBUG 

   MultilayerPerceptron* multilayerPerceptron = objectiveFunctional->getMultilayerPerceptron();

   int numberOfFreeParameters = multilayerPerceptron->getNumberOfFreeParameters();   

   if(newPopulation.getNumberOfRows() != populationSize)
   {
      std::cerr << std::endl
                << "Flood Error: EvolutionaryAlgorithm class." << std::endl
                << "void setPopulation(Matrix<double>) method." << std::endl
                << "New population size is not equal to population size." << std::endl
                << std::endl;

      exit(1);
   }
   else if(newPopulation.getNumberOfColumns() != numberOfFreeParameters)
   {
      std::cerr << std::endl
                << "Flood Error: EvolutionaryAlgorithm class." << std::endl
                << "void setPopulation(Matrix<double>) method." << std::endl
                << "New number of free parameters is not equal to number of free parameters." << std::endl
                << std::endl;

      exit(1);
   }

   #endif

   // Set population

   population = newPopulation;
}


// void setEvaluation(Vector<double>) method

/// This method sets a new population evaluation vector.
///
/// @param newEvaluation Population evaluation values.

void EvolutionaryAlgorithm::setEvaluation(Vector<double> newEvaluation)
{
   // Control sentence (if debug)

   #ifndef NDEBUG 

   if(newEvaluation.getSize() != populationSize)
   {
      std::cerr << std::endl
                << "Flood Error: EvolutionaryAlgorithm class." << std::endl
                << "void setEvaluation(Vector<double>) method." << std::endl
                << "Size is not equal to population size." << std::endl;

      exit(1);
   }

   #endif

   // Set evaluation

   evaluation = newEvaluation;
}


// void setFitness(Vector<double>) method

/// This method sets a new population fitness vector.
///
/// @param newFitness Population fitness values.

void EvolutionaryAlgorithm::setFitness(Vector<double> newFitness)
{
   // Control sentence (if debug)

   #ifndef NDEBUG 

   if(newFitness.getSize() != populationSize)
   {
      std::cerr << std::endl
                << "Flood Error: EvolutionaryAlgorithm class." << std::endl
                << "void setFitness(Vector<double>) method." << std::endl
                << "Size is not equal to population size." << std::endl
				<< std::endl;

      exit(1);
   }

   #endif

   // Set fitness

   fitness = newFitness;
}


// void setSelection(Vector<bool>) method

/// This method sets a new population selection vector.
///
/// @param newSelection Population selection values.

void EvolutionaryAlgorithm::setSelection(Vector<bool> newSelection)
{
   // Control sentence (if debug)

   #ifndef NDEBUG 

   if(newSelection.getSize() != populationSize)
   {
      std::cerr << std::endl
                << "Flood Error: EvolutionaryAlgorithm class." << std::endl
                << "void setSelection(Vector<double>) method." << std::endl
                << "Size is not equal to population size." << std::endl
                << std::endl;

      exit(1);
   }

   #endif

   // Set selection

   selection = newSelection;
}


// void setReservePopulationHistory(bool) method

/// This method makes the population history vector of matrices to be reseved or not in memory.
///
/// @param newReservePopulationHistory True if the population history vector of matrices is to be reserved, false 
/// otherwise.

void EvolutionaryAlgorithm::setReservePopulationHistory(bool newReservePopulationHistory)
{
   reservePopulationHistory = newReservePopulationHistory;
}


// void setReserveMeanNormHistory(bool) method

/// This method makes the mean norm history vector to be reseved or not in memory.
///
/// @param newReserveMeanNormHistory True if the mean norm history vector is to be reserved, false otherwise.

void EvolutionaryAlgorithm::setReserveMeanNormHistory(bool newReserveMeanNormHistory)
{
   reserveMeanNormHistory = newReserveMeanNormHistory;
}


// void setReserveStandardDeviationNormHistory(bool) method

/// This method makes the standard deviation norm history vector to be reseved or not in memory.
///
/// @param newReserveStandardDeviationNormHistory True if the standard deviation norm history vector is to be 
/// reserved, false otherwise.

void EvolutionaryAlgorithm::setReserveStandardDeviationNormHistory(bool newReserveStandardDeviationNormHistory)
{
   reserveStandardDeviationNormHistory = newReserveStandardDeviationNormHistory;
}


// void setReserveBestNormHistory(bool) method

/// This method makes the best norm history vector to be reseved or not in memory.
///
/// @param newReserveBestNormHistory True if the best norm history vector is to be reserved, false otherwise.

void EvolutionaryAlgorithm::setReserveBestNormHistory(bool newReserveBestNormHistory)
{
   reserveBestNormHistory = newReserveBestNormHistory;
}


// void setReserveMeanEvaluationHistory(bool) method

/// This method makes the mean evaluation history vector to be reseved or not in memory.
///
/// @param newReserveMeanEvaluationHistory True if the mean evaluation history vector is to be reserved, false 
/// otherwise.

void EvolutionaryAlgorithm::setReserveMeanEvaluationHistory(bool newReserveMeanEvaluationHistory) 
{
   reserveMeanEvaluationHistory = newReserveMeanEvaluationHistory;
}


// void setReserveStandardDeviationEvaluationHistory(bool) method

/// This method makes the standard deviation evaluation history vector to be reseved or not in memory.
///
/// @param newReserveStandardDeviationEvaluationHistory True if the standard deviation evaluation history vector 
/// is to be reserved, false otherwise.

void EvolutionaryAlgorithm
::setReserveStandardDeviationEvaluationHistory(bool newReserveStandardDeviationEvaluationHistory)
{
   reserveStandardDeviationEvaluationHistory = newReserveStandardDeviationEvaluationHistory;
}


// void setReserveBestEvaluationHistory(bool) method

/// This method makes the best evaluation history vector to be reseved or not in memory.
///
/// @param newReserveBestEvaluationHistory True if the best evaluation history vector is to be reserved, 
/// false otherwise.

void EvolutionaryAlgorithm::setReserveBestEvaluationHistory(bool newReserveBestEvaluationHistory)
{
   reserveBestEvaluationHistory = newReserveBestEvaluationHistory;
}


// void setReserveAllTrainingHistory(bool) method

/// This method makes the training history of all variables to reseved or not in memory.
///
/// @param newReserveAllTrainingHistory True if the training history of all variables is to be reserved, 
/// false otherwise.

void EvolutionaryAlgorithm::setReserveAllTrainingHistory(bool newReserveAllTrainingHistory)
{
   reserveElapsedTimeHistory = newReserveAllTrainingHistory;
   reserveFreeParametersHistory = newReserveAllTrainingHistory;
   reserveFreeParametersNormHistory = newReserveAllTrainingHistory;
   reserveEvaluationHistory = newReserveAllTrainingHistory;

   reservePopulationHistory = newReserveAllTrainingHistory;
   reserveMeanNormHistory = newReserveAllTrainingHistory;
   reserveStandardDeviationNormHistory = newReserveAllTrainingHistory;
   reserveBestNormHistory = newReserveAllTrainingHistory;
   reserveMeanEvaluationHistory = newReserveAllTrainingHistory;
   reserveStandardDeviationEvaluationHistory = newReserveAllTrainingHistory;
   reserveBestEvaluationHistory = newReserveAllTrainingHistory;
}


// void setPopulationHistory(Vector< Matrix<double> >) method

/// This method sets a new matrix containing the training direction history over the training epochs.
/// Each element in the vector must contain the population matrix of one single epoch. 
///
/// @param newPopulationHistory Population history vector of matrices. 

void EvolutionaryAlgorithm::setPopulationHistory(Vector< Matrix<double> > newPopulationHistory)
{
   populationHistory = newPopulationHistory;
}


// void setMeanNormHistory(Vector<double>) method

/// This method sets a new vector containing the mean norm history over the training epochs.
/// Each element in the vector must contain the mean norm of one single epoch. 
///
/// @param newMeanNormHistory Mean norm history vector.  

void EvolutionaryAlgorithm::setMeanNormHistory(Vector<double> newMeanNormHistory)
{
   meanNormHistory = newMeanNormHistory; 
}


// void setStandardDeviationNormHistory(Vector<double>) method

/// This method sets a new vector containing the standard deviation norm history over the training epochs.
/// Each element in the vector must contain the standard deviation norm of one single epoch. 
///
/// @param newStandardDeviationNormHistory Standard deviation norm history vector.  

void EvolutionaryAlgorithm::setStandardDeviationNormHistory(Vector<double> newStandardDeviationNormHistory)
{
   standardDeviationNormHistory = newStandardDeviationNormHistory;
}


// void setBestNormHistory(Vector<double>) method

/// This method sets a new vector containing the best norm history over the training epochs.
/// Each element in the vector must contain the best norm of one single epoch. 
///
/// @param newBestNormHistory Best norm history vector.  

void EvolutionaryAlgorithm::setBestNormHistory(Vector<double> newBestNormHistory)
{
   bestNormHistory = newBestNormHistory;
}


// void setMeanEvaluationHistory(Vector<double>) method

/// This method sets a new vector containing the mean evaluation history over the training epochs.
/// Each element in the vector must contain the mean evaluation of one single epoch. 
///
/// @param newMeanEvaluationHistory Mean evaluation history vector.  

void EvolutionaryAlgorithm::setMeanEvaluationHistory(Vector<double> newMeanEvaluationHistory)
{
   meanEvaluationHistory = newMeanEvaluationHistory;
}


// void setStandardDeviationEvaluationHistory(Vector<double>) method

/// This method sets a new vector containing the standard deviation evaluation history over the training epochs.
/// Each element in the vector must contain the standard deviation evaluation of one single epoch. 
///
/// @param newStandardDeviationEvaluationHistory Standard deviation evaluation history vector.  

void EvolutionaryAlgorithm
::setStandardDeviationEvaluationHistory(Vector<double> newStandardDeviationEvaluationHistory)
{
   standardDeviationEvaluationHistory = newStandardDeviationEvaluationHistory;
}


// void setBestEvaluationHistory(Vector<double>) method

/// This method sets a new vector containing the best evaluation history over the training epochs.
/// Each element in the vector must contain the best evaluation of one single epoch. 
///
/// @param newBestEvaluationHistory Best evaluation history vector.  

void EvolutionaryAlgorithm::setBestEvaluationHistory(Vector<double> newBestEvaluationHistory)
{
   bestEvaluationHistory = newBestEvaluationHistory;
}


// Vector<double> getIndividual(int) method

/// This method returns the Vector of free parameters corresponding to the individual i in the population.
///
/// @param i Index of individual in the population.

Vector<double> EvolutionaryAlgorithm::getIndividual(int i)
{
   // Control sentence (if debug)

   #ifndef NDEBUG 

   if(i >= populationSize)
   {
      std::cerr << std::endl
                << "Flood Error: EvolutionaryAlgorithm class." << std::endl
                << "Vector<double> getIndividual(int) method." << std::endl
                << "Index must be less than population size." << std::endl
				<< std::endl;

      exit(1);
   }
  
   #endif

   // Get individual

   Vector<double> individual = population.getRow(i);

   return(individual);
}


// setIndividual(int, Vector<double>) method

/// This method sets a new Vector of free parameters to the individual i in the population. 
///
/// @param i Index of individual in the population.
/// @param individual Vector of free parameters to be assigned to individual i.

void EvolutionaryAlgorithm::setIndividual(int i, Vector<double> individual)
{
   // Control sentence (if debug)

   #ifndef NDEBUG 
 
   int size = individual.getSize();

   MultilayerPerceptron* multilayerPerceptron = objectiveFunctional->getMultilayerPerceptron();

   int numberOfFreeParameters = multilayerPerceptron->getNumberOfFreeParameters();

   if(i >= populationSize)
   {
      std::cerr << std::endl
                << "Flood Error: EvolutionaryAlgorithm class." << std::endl
                << "setIndividual(int, Vector<double>) method." << std::endl
                << "Index must be less than population size." << std::endl
				<< std::endl;

      exit(1);
   }
   else if(size != numberOfFreeParameters)
   {
      std::cerr << std::endl
                << "Flood Error: EvolutionaryAlgorithm class." << std::endl
                << "setIndividual(int, Vector<double>) method." << std::endl
                << "Size must be equal to number of free parameters." << std::endl
				<< std::endl;

      exit(1);
   }
  
   #endif

   // Get individual

   population.setRow(i, individual);
}


// Vector<double> calculatePopulationNorm(void) method

/// This method returns a vector containing the norm of each individual in the population.

Vector<double> EvolutionaryAlgorithm::calculatePopulationNorm(void)
{
   Vector<double> populationNorm(populationSize);
      
   for(int i = 0; i < populationSize; i++)
   {
      Vector<double> individual = getIndividual(i); 
           
      populationNorm[i] = individual.calculateNorm();     
   }               
   
   return(populationNorm);            
}


// Training parameters

// double getSelectivePressure(void) method

/// This method returns the selective pressure value.
///
/// @see performLinearRankingFitnessAssignment(void).

double EvolutionaryAlgorithm::getSelectivePressure(void)
{
   return(selectivePressure);
}


// double getRecombinationSize(void) method

/// This method returns the recombination size value.
///
/// @see performIntermediateRecombination(void)
/// @see performLineRecombination(void).

double EvolutionaryAlgorithm::getRecombinationSize(void)
{
   return(recombinationSize);
}


// double getMutationRate(void) method

/// This method returns the mutation rate value.
///
///  @see performNormalMutation(void).
///  @see performUniformMutation(void).

double EvolutionaryAlgorithm::getMutationRate(void)
{
   return(mutationRate);
}


// double getMutationRange(void) method

/// This method returns the mutation range value.
///
///  @see performNormalMutation(void).
///  @see performUniformMutation(void).

double EvolutionaryAlgorithm::getMutationRange(void)
{
   return(mutationRange);
}


// double getMaximumNumberOfGenerations(void) method

/// This method returns the maximum number of generations to train.

double EvolutionaryAlgorithm::getMaximumNumberOfGenerations(void)
{
   return(maximumNumberOfGenerations);
}

// FitnessAssignmentMethod getFitnessAssignmentMethod(void) method

/// This method returns the fitness assignment method used for training.
///
/// @see performLinearRankingFitnessAssignment(void).
/// @see train(void).
 
EvolutionaryAlgorithm::FitnessAssignmentMethod EvolutionaryAlgorithm::getFitnessAssignmentMethod(void)
{
   return(fitnessAssignmentMethod);
}


// SelectionMethod getSelectionMethod(void) method

/// This method returns the selection method used for training.
///
/// @see performRouletteWheelSelection(void).
/// @see performStochasticUniversalSamplingSelection(void).
/// @see train(void).

EvolutionaryAlgorithm::SelectionMethod EvolutionaryAlgorithm::getSelectionMethod(void)
{
   return(selectionMethod);
}


// RecombinationMethod getRecombinationMethod(void) method

/// This method returns the recombination method used for training.
///
/// @see performIntermediateRecombination(void).
/// @see performLineRecombination(void).
/// @see train(void).

EvolutionaryAlgorithm::RecombinationMethod EvolutionaryAlgorithm::getRecombinationMethod(void)
{
   return(recombinationMethod);
}


// MutationMethod getMutationMethod(void) method

/// This method returns the mutation method used for training.
///
/// @see performNormalMutation(void).
/// @see performUniformMutation(void).
/// @see train(void).

EvolutionaryAlgorithm::MutationMethod EvolutionaryAlgorithm::getMutationMethod(void)
{
   return(mutationMethod);
}


// void setSelectivePressure(double) method

/// This method sets a new value for the selective pressure parameter.
/// Linear ranking allows values for the selective pressure between 1 and 2.
///
/// @param newSelectivePressure Selective pressure value. This must be between 1 and 2 for linear ranking fitness
/// assignment. 
///
/// @see performLinearRankingFitnessAssignment(void).

void EvolutionaryAlgorithm::setSelectivePressure(double newSelectivePressure)
{
   switch(fitnessAssignmentMethod)
   {
      case LinearRanking:

         if((newSelectivePressure < 1.0) || (newSelectivePressure > 2.0))
         {
            std::cerr << std::endl
                      << "Flood Error: EvolutionaryAlgorithm class." << std::endl
                      << "void setSelectivePressure(double) method. "
                      << "Case linear ranking." << std::endl
                      << "Selective pressure must be a value between 1 and 2." << std::endl 
                      << std::endl;

            exit(1);
         }
         
		 // Set selective pressure

		 selectivePressure = newSelectivePressure;

      break;
   }
}


// void setRecombinationSize(double) method

/// This method sets a new value for the recombination size parameter.
/// The recombination size value must be equal or greater than 0.
///
/// @param newRecombinationSize Recombination size value. This must be equal or greater than 0.
///
/// @see performIntermediateRecombination(void)
/// @see performLineRecombination(void).

void EvolutionaryAlgorithm::setRecombinationSize(double newRecombinationSize)
{
   if(newRecombinationSize < 0.0)
   {
      std::cerr << std::endl
                << "Flood Error: EvolutionaryAlgorithm class." << std::endl
                << "void setRecombinationSize(double) method." << std::endl
                << "Recombination size must be equal or greater than 0." << std::endl
                << std::endl;

      exit(1);
   }

   // Set recombination size

   recombinationSize = newRecombinationSize;
}


// void setMutationRate(double) method

/// This method sets a new value for the mutation rate parameter.
/// The mutation rate value must be between 0 and 1.
///
/// @param newMutationRate Mutation rate value. This value must lie in the interval [0,1]. 
///
/// @see performNormalMutation(void).
/// @see performUniformMutation(void).

void EvolutionaryAlgorithm::setMutationRate(double newMutationRate)
{
   // Control sentence

   if(newMutationRate < 0.0 || newMutationRate > 1.0)
   {
      std::cerr << std::endl
                << "Flood Error: EvolutionaryAlgorithm class." << std::endl
                << "void setMutationRate(double) method." << std::endl
                << "Mutation rate must be a value between 0 and 1. " << std::endl
                << std::endl;

      exit(1);
   }

   // Set mutation rate

   mutationRate = newMutationRate;
}


// void setMutationRange(double) method

/// This method sets a new value for the mutation range parameter.
/// The mutation range value must be 0 or a positive number. 
///
/// @param newMutationRange Mutation range value. This must be equal or greater than 0.
///
/// @see performNormalMutation(void).
/// @see performUniformMutation(void).

void EvolutionaryAlgorithm::setMutationRange(double newMutationRange)
{
   // Control sentence

   if(newMutationRange < 0.0)
   {
      std::cerr << std::endl
                << "Flood Error: EvolutionaryAlgorithm class." << std::endl
                << "void setMutationRange(double) method." << std::endl
                << "Mutation range must be equal or greater than 0. " << std::endl
                << std::endl;

      exit(1);
   }

   // Set mutation range

   mutationRange = newMutationRange;
}


// void setMaximumNumberOfGenerations(int) method

/// This method sets a new value for the maximum number of generations to train.
/// The maximum number of generations value must be a positive number. 
///
/// @param newMaximumNumberOfGenerations Maximum number of generations value.

void EvolutionaryAlgorithm::setMaximumNumberOfGenerations(int newMaximumNumberOfGenerations)
{
   // Control sentence (if debug)

   #ifndef NDEBUG 

   if(newMaximumNumberOfGenerations == 0)
   {
      std::cerr << std::endl
                << "Flood Error: EvolutionaryAlgorithm class." << std::endl
                << "void setMaximumNumberOfGenerations(int) method." << std::endl
                << "Maximum number of generations must be greater than 0. " 
                << std::endl << std::endl;

      exit(1);
   }

   #endif

   // Set maximum number of generations

   maximumNumberOfGenerations = newMaximumNumberOfGenerations;
}


// void setFitnessAssignmentMethod(FitnessAssignmentMethod) method

/// This method sets a new fitness assignment method to be used for training.
///
/// @param newFitnessAssignmentMethod Fitness assignment method chosen for training.
///
/// @see performLinearRankingFitnessAssignment(void).

void EvolutionaryAlgorithm::setFitnessAssignmentMethod
(EvolutionaryAlgorithm::FitnessAssignmentMethod newFitnessAssignmentMethod)
{
   fitnessAssignmentMethod = newFitnessAssignmentMethod;
}


// void setSelectionMethod(SelectionMethod) method

/// This method sets a new selection method to be used for training.
///
/// @param newSelectionMethod Selection method chosen for training.
///
/// @see performRouletteWheelSelection(void).
/// @see performStochasticUniversalSamplingSelection(void).

void EvolutionaryAlgorithm::setSelectionMethod(EvolutionaryAlgorithm::SelectionMethod newSelectionMethod)
{
   selectionMethod = newSelectionMethod;
}


// void setRecombinationMethod(RecombinationMethod) method

/// This method sets a new recombination method to be used for training.
///
/// @param newRecombinationMethod Recombination method chosen for training. 
///
/// @see performIntermediateRecombination(void).
/// @see performLineRecombination(void).

void EvolutionaryAlgorithm
::setRecombinationMethod(EvolutionaryAlgorithm::RecombinationMethod newRecombinationMethod)
{
   recombinationMethod = newRecombinationMethod;
}


// void setMutationMethod(MutationMethod) method

/// This method sets a new mutation method to be used for training.
///
/// @param newMutationMethod Mutation method chosen for training. 
///
/// @see performNormalMutation(void).
/// @see performUniformMutation(void).


void EvolutionaryAlgorithm::setMutationMethod(EvolutionaryAlgorithm::MutationMethod newMutationMethod)
{
   mutationMethod = newMutationMethod;
}


// void initPopulationUniform(void) method

/// This method initializes the free parameters of all the individuals in the population at random, with values 
/// comprised between -1 and 1.

void EvolutionaryAlgorithm::initPopulationUniform(void)
{
   MultilayerPerceptron* multilayerPerceptron = objectiveFunctional->getMultilayerPerceptron();

   int numberOfFreeParameters = multilayerPerceptron->getNumberOfFreeParameters();   

   Vector<double> individual(numberOfFreeParameters);

   for(int i = 0; i < populationSize; i++)
   {
      individual = getIndividual(i);

      for(int j = 0; j < numberOfFreeParameters; j++)
      {
         double random = (double)rand()/(RAND_MAX+1.0);

         individual[j] = -1.0 + 2.0*random;
      }

      setIndividual(i, individual);
   }
}


// void initPopulationUniform(double, double) method

/// This method initializes the free parameters of all the individuals in the population at random, with values 
/// comprised between a minimum and a maximum value.
///
/// @param minimum Minimum initialization value.
/// @param maximum Maximum initialization value.

void EvolutionaryAlgorithm::initPopulationUniform(double minimum, double maximum)
{
   // Control sentence (if debug)

   #ifndef NDEBUG 

   if(minimum > maximum)
   {
      std::cerr << std::endl
                << "Flood Error: EvolutionaryAlgorithm class." << std::endl
                << "void initPopulationUniform(double, double) method." << std::endl
                << "Minimum initialization value must be less than maximum value." << std::endl
				<< std::endl;

      exit(1);
   }

   #endif

   MultilayerPerceptron* multilayerPerceptron = objectiveFunctional->getMultilayerPerceptron();

   int numberOfFreeParameters = multilayerPerceptron->getNumberOfFreeParameters();   

   Vector<double> individual(numberOfFreeParameters);

   for(int i = 0; i < populationSize; i++)
   {
      individual = getIndividual(i);

      for(int j = 0; j < numberOfFreeParameters; j++)
      {
         double random = (double)rand()/(RAND_MAX+1.0);

         individual[j] = minimum + (maximum-minimum)*random;
      }

      setIndividual(i, individual);
   }
}


// void initPopulationUniform(Vector<double>, Vector<double>) method

/// This method initializes the free parameters of all the individuals in the population at random, with values 
/// comprised between different minimum and maximum values for each variable.
///
/// @param minimum Vector of minimum initialization values.
/// @param maximum Vector of maximum initialization values.

void EvolutionaryAlgorithm::initPopulationUniform(Vector<double> minimum, Vector<double> maximum)
{
   MultilayerPerceptron* multilayerPerceptron = objectiveFunctional->getMultilayerPerceptron();

   int numberOfFreeParameters = multilayerPerceptron->getNumberOfFreeParameters();   

   // Control sentence (if debug)

   #ifndef NDEBUG 

   int minimumSize = minimum.getSize();
   int maximumSize = maximum.getSize();

   if(minimumSize != numberOfFreeParameters || maximumSize != numberOfFreeParameters)   
   {
      std::cerr << std::endl
                << "Flood Error: EvolutionaryAlgorithm class." << std::endl
                << "void initPopulationUniform(Vector<double>, Vector<double>)." << std::endl
                << "Minimum value and maximum value sizes must be equal to number of free parameters." 
                << std::endl << std::endl;
 
      exit(1);
   }

   #endif

   Vector<double> individual(numberOfFreeParameters);

   for(int i = 0; i < populationSize; i++)
   {
      individual = getIndividual(i);

      for(int j = 0; j < numberOfFreeParameters; j++)
      {
         double random = (double)rand()/(RAND_MAX+1.0);

         individual[j] = minimum[j] + (maximum[j]-minimum[j])*random;
      }

      setIndividual(i, individual);
   }
}


// void initPopulationUniform(Matrix<double>) method

/// This method initializes the free parameters of all the individuals in the population at random, with values 
/// comprised between different minimum and maximum values for each variable.
/// All minimum and maximum values are given from a single matrix.
/// The first row must contain the minimum value for each free parameter.
/// The second row must contain the maximum value for each free parameter.
///
/// @param minimumAndMaximum Matrix of minimum and maximum initialization values.

void EvolutionaryAlgorithm::initPopulationUniform(Matrix<double> minimumAndMaximum)
{
   MultilayerPerceptron* multilayerPerceptron = objectiveFunctional->getMultilayerPerceptron();

   int numberOfFreeParameters = multilayerPerceptron->getNumberOfFreeParameters();   

   // Control sentence (if debug)

   #ifndef NDEBUG 

   int numberOfRows = minimumAndMaximum.getNumberOfRows();
   int numberOfColumns = minimumAndMaximum.getNumberOfColumns();

   if(numberOfRows != 2 || numberOfColumns != numberOfFreeParameters)   
   {
      std::cerr << std::endl
                << "Flood Error: EvolutionaryAlgorithm class." << std::endl
                << "void initPopulationUniform(Matrix<double>)." << std::endl
                << "Number of rows must be two and number of columns must be equal to number of free parameters."
                << std::endl << std::endl;
 
      exit(1);
   }

   #endif

   Vector<double> minimum = minimumAndMaximum.getRow(0);
   Vector<double> maximum = minimumAndMaximum.getRow(1);

   Vector<double> individual(numberOfFreeParameters);

   for(int i = 0; i < populationSize; i++)
   {
      individual = getIndividual(i);

      for(int j = 0; j < numberOfFreeParameters; j++)
      {
         double random = (double)rand()/(RAND_MAX+1.0);

         individual[j] 
         = minimum[j] + (maximum[j]-minimum[j])*random;
      }

      setIndividual(i, individual);
   }
}


// void initPopulationNormal(void) method

/// This method initializes the free parameters of all the individuals in the population with random values chosen
/// from a normal distribution with mean 0 and standard deviation 1.

void EvolutionaryAlgorithm::initPopulationNormal(void)
{
   double mean = 0.0;
   double standardDeviation = 1.0;     

   MultilayerPerceptron* multilayerPerceptron
   = objectiveFunctional->getMultilayerPerceptron();

   int numberOfFreeParameters = multilayerPerceptron->getNumberOfFreeParameters();   

   Vector<double> individual(numberOfFreeParameters);

   double pi = 4.0*atan(1.0);
   
   double random1 = 0.0;
   double random2 = 0.0;

   for(int i = 0; i < populationSize; i++)
   {
      individual = getIndividual(i);

      for(int j = 0; j < numberOfFreeParameters; j++)
      {
         do
         {
            random1 = (double)rand()/(RAND_MAX+1.0);

         }while(random1 == 0.0);

         random2 = (double)rand()/(RAND_MAX+1.0);

         // Box-Muller transformation

         double normallyDistributedRandomNumber 
         = mean + sqrt(-2.0*log(random1))*sin(2.0*pi*random2)*standardDeviation;

         individual[j] = normallyDistributedRandomNumber;
      }

      setIndividual(i, individual);
   }
}


// void initPopulationNormal(double, double) method

/// This method initializes the free parameters of all the individuals in the population with random values chosen
/// from a normal distribution with a given mean and a given standard deviation.
///
/// @param mean Mean of normal distribution.
/// @param standardDeviation Standard deviation of normal distribution.

void EvolutionaryAlgorithm::initPopulationNormal(double mean, double standardDeviation)
{
   MultilayerPerceptron* multilayerPerceptron = objectiveFunctional->getMultilayerPerceptron();

   int numberOfFreeParameters = multilayerPerceptron->getNumberOfFreeParameters();   

   Vector<double> individual(numberOfFreeParameters);

   double pi = 4.0*atan(1.0);

   double random1 = 0.0;
   double random2 = 0.0;

   for(int i = 0; i < populationSize; i++)
   {
      individual = getIndividual(i);

      for(int j = 0; j < numberOfFreeParameters; j++)
      {
         do
         {
            random1 = (double)rand()/(RAND_MAX+1.0);
       
         }while(random1 == 0.0);

         random2 = (double)rand()/(RAND_MAX+1.0);

         // Box-Muller transformation

         double normallyDistributedRandomNumber 
         = mean + sqrt(-2.0*log(random1))*sin(2.0*pi*random2)*standardDeviation;

         individual[j] = normallyDistributedRandomNumber;
      }

      setIndividual(i, individual);
   }
}


// void initPopulationNormal(Vector<double>, Vector<double>) method

/// This method initializes the free parameters of all the individuals in the population with random values chosen
/// from normal distributions with different mean and standard deviation for each free parameter.
///
/// @param mean Vector of mean values.
/// @param standardDeviation Vector of standard deviation values.

void EvolutionaryAlgorithm::initPopulationNormal(Vector<double> mean, Vector<double> standardDeviation)
{
   MultilayerPerceptron* multilayerPerceptron = objectiveFunctional->getMultilayerPerceptron();

   int numberOfFreeParameters = multilayerPerceptron->getNumberOfFreeParameters();   

   // Control sentence (if debug)

   #ifndef NDEBUG 

   int meanSize = mean.getSize();
   int standardDeviationSize = standardDeviation.getSize();

   if(meanSize != numberOfFreeParameters || standardDeviationSize != numberOfFreeParameters)   
   {
      std::cerr << std::endl
                << "Flood Error: EvolutionaryAlgorithm class." << std::endl
                << "void initPopulationNormal(Vector<double>, Vector<double>)." << std::endl
                << "Mean and standard deviation sizes must be equal to number of free parameters." 
                << std::endl << std::endl;
 
      exit(1);
   }

   #endif

   Vector<double> individual(numberOfFreeParameters);

   double pi = 4.0*atan(1.0);

   double random1 = 0.0;
   double random2 = 0.0;

   for(int i = 0; i < populationSize; i++)
   {
      individual = getIndividual(i);

      for(int j = 0; j < numberOfFreeParameters; j++)
      {
         do
         {
            random1 = (double)rand()/(RAND_MAX+1.0);

         }while(random1 == 0.0);

         random2 = (double)rand()/(RAND_MAX+1.0);

         // Box-Muller transformation

         double normallyDistributedRandomNumber 
         = mean[j] + sqrt(-2.0*log(random1))*sin(2.0*pi*random2)*standardDeviation[j];

         individual[j] = normallyDistributedRandomNumber;
      }

      setIndividual(i, individual);
   }

}

// void initPopulationNormal(Matrix<double>) method

/// This method initializes the free parameters of all the individuals in the population with random values chosen
/// from normal distributions with different mean and standard deviation for each parameter.
/// All mean and standard deviation values are given from a single matrix.
/// The first row must contain the mean value for each free parameter.
/// The second row must contain the standard deviation value for each free parameter.
///
/// @param meanAndStandardDeviation Matrix of mean and standard deviation values.

void EvolutionaryAlgorithm::initPopulationNormal(Matrix<double> meanAndStandardDeviation)
{
   MultilayerPerceptron* multilayerPerceptron = objectiveFunctional->getMultilayerPerceptron();

   int numberOfFreeParameters = multilayerPerceptron->getNumberOfFreeParameters();   

   // Control sentence (if debug)

   #ifndef NDEBUG 

   int numberOfRows = meanAndStandardDeviation.getNumberOfRows();
   int numberOfColumns = meanAndStandardDeviation.getNumberOfColumns();

   if(numberOfRows != 2 || numberOfColumns != numberOfFreeParameters)   
   {
      std::cerr << std::endl
                << "Flood Error: EvolutionaryAlgorithm class." << std::endl
                << "void initPopulationNormal(Matrix<double>)." << std::endl
                << "Number of rows must be two and number of columns must be equal to "
                << "number of free parameters." 
                << std::endl << std::endl;
 
      exit(1);
   }

   #endif

   Vector<double> mean = meanAndStandardDeviation.getRow(0);
   Vector<double> standardDeviation = meanAndStandardDeviation.getRow(1);

   Vector<double> individual(numberOfFreeParameters);

   double pi = 4.0*atan(1.0);

   double random1 = 0.0;
   double random2 = 0.0;

   for(int i = 0; i < populationSize; i++)
   {
      individual = getIndividual(i);

      for(int j = 0; j < numberOfFreeParameters; j++)
      {
         do
         {
            random1 = (double)rand()/(RAND_MAX+1.0);

         }while(random1 == 0.0);

         random2 = (double)rand()/(RAND_MAX+1.0);

         // Box-Muller transformation

         double normallyDistributedRandomNumber 
         = mean[j] + sqrt(-2.0*log(random1))*sin(2.0*pi*random2)*standardDeviation[j];

         individual[j] = normallyDistributedRandomNumber;
      }

      setIndividual(i, individual);
   }
}


// void evaluatePopulation(void) method

/// This method evaluates the objective functional of all individuals in the population. 
/// Results are stored in the evaluation vector.
///
/// @see train(void).

void EvolutionaryAlgorithm::evaluatePopulation(void)
{
   // Multilayer perceptron

   MultilayerPerceptron* multilayerPerceptron = objectiveFunctional->getMultilayerPerceptron();

   int numberOfFreeParameters = multilayerPerceptron->getNumberOfFreeParameters();

   Vector<double> individual(numberOfFreeParameters);

   // Evaluate objective functional for all individuals

   for(int i = 0; i < populationSize; i++)
   {
      individual = getIndividual(i);

      evaluation[i] = objectiveFunctional->calculatePotentialEvaluation(individual);
      
      if(!(evaluation[i] > -1.0e69 && evaluation[i] < 1.0e69))
      {
         std::cerr << std::endl
                   << "Flood Error: EvolutionaryAlgorithm class." << std::endl
                   << "void evaluatePopulation(void) method." << std::endl
                   << "Evaluation of individual " << i << " is not a real number." << std::endl
                   << std::endl;

         exit(1);
      }                
   }
}


// void performLinearRankingFitnessAssignment(void) method
//
/// This method ranks all individuals in the population by their objective evaluation, so that the least fit 
/// individual has rank 1 and the fittest individual has rank [population size].
/// It then assigns them a fitness value linearly proportional to their rank. Results are stored in the fitness 
/// vector.
///
/// @see train(void).

void EvolutionaryAlgorithm::performLinearRankingFitnessAssignment(void)
{
   // Sorted evaluation vector

   Vector<double> sortedEvaluation(populationSize);

   sortedEvaluation = evaluation;

   std::sort(sortedEvaluation.begin(), sortedEvaluation.end(), std::less<double>());

   // Rank vector

   Vector<int> rank(populationSize);

   for(int i = 0; i < populationSize; i++)
   {
      for(int j = 0; j < populationSize; j++)
      {
         if(evaluation[j] == sortedEvaluation[i])
         {
            rank[j] = populationSize - i;
         }
      }
   }

   // Perform linear ranking fitness assignment

   for(int i = 0; i < populationSize; i++)
   {
      fitness[i] = 2.0 - selectivePressure
      + 2.0*(selectivePressure - 1.0)*(rank[i] - 1.0)/(populationSize - 1.0);
      
      if(!(fitness[i] > -1.0e69 && fitness[i] < 1.0e69))
      {
         std::cerr << std::endl
                   << "Flooe Error: EvolutionaryAlgorithm class." << std::endl
                   << "void performLinearRankingFitnessAssignment(void) method." << std::endl
                   << "Fitness of individual " << i << " is not a real number." << std::endl
                   << std::endl;

         exit(1);
      }          
   }
}

// void performRouletteWheelSelection(void) method

/// This metod performs selection with roulette wheel selection. It selects half of the individuals from the 
/// population. 
/// Results are stored in the selection vector. 
///
/// @see performStochasticUniversalSamplingSelection(void) method
/// @see train(void).

void EvolutionaryAlgorithm::performRouletteWheelSelection(void)
{
   // Set selection vector to false 

   for(int i = 0; i < populationSize; i++)
   {
      selection[i] = false;
   }

   int numberOfSelectedIndividuals = populationSize/2;

   Vector<double> cumulativeFitness(populationSize);

   // Cumulative fitness vector

   cumulativeFitness[0] = fitness[0]; 

   for(int i = 1; i < populationSize; i++)
   {
      cumulativeFitness[i] = cumulativeFitness[i-1] + fitness[i];
   }

   // Select individuals until the desired number of selections is obtained

   int countNumberOfSelectedIndividuals = 0;

   do
   {
      // Random number between 0 and total cumulative fitness

      double random = (double)rand()/(RAND_MAX+1.0);

      double pointer = cumulativeFitness[populationSize-1]*random;

      // Perform selection

      if(pointer < cumulativeFitness[0])
      {
         if(selection[0] == false)
         {
            selection[0] = true;
            countNumberOfSelectedIndividuals++;
         }
      }
      
      for(int i = 1; i < populationSize; i++)
      {
         if(pointer < cumulativeFitness[i] && pointer >= cumulativeFitness[i-1])
         {
            if(selection[i] == false)
            {
               selection[i] = true;
               countNumberOfSelectedIndividuals++;
            }
         }
      }
   }while(countNumberOfSelectedIndividuals != numberOfSelectedIndividuals);

   // Control sentence

   if(countNumberOfSelectedIndividuals != numberOfSelectedIndividuals)
   {
      std::cerr << std::endl 
                << "Flood Error: EvolutionaryAlgorithm class." << std::endl
                << "void performRouletteWheelSelection(void) method." << std::endl
                << "Count number of selected individuals is not equal to number of selected individuals."
                << std::endl << std::endl;

      exit(1);
   }
}


// void performStochasticUniversalSamplingSelection(void) method
//
/// This metod performs selection with stochastic universal sampling. It selects half of the individuals from the
/// population. 
/// Results are stored in the selection vector. 
///
/// @see performRouletteWheelSelection(void) method.
/// @see train(void).


void EvolutionaryAlgorithm::performStochasticUniversalSamplingSelection(void)
{
   // Set selection vector to false

   for(int i = 0; i < populationSize; i++)
   {
      selection[i] = false;
   }
 
   int numberOfSelectedIndividuals = populationSize/2;

   Vector<double> cumulativeFitness(populationSize);

   Vector<double> pointer(numberOfSelectedIndividuals);

   // Cumulative fitness vector

   cumulativeFitness[0] = fitness[0];

   for(int i = 1; i < populationSize; i++)
   {  
      cumulativeFitness[i] = cumulativeFitness[i-1] + fitness[i];
   }


   // Pointer vector

   // Random number between 0 and totalCumulativeFitnees/(double)numberOfSelectedIndividuals 

   double random = (double)rand()/(RAND_MAX+1.0);

   pointer[0] = random
   *cumulativeFitness[populationSize-1]/(double)numberOfSelectedIndividuals;

   for(int i = 1; i < numberOfSelectedIndividuals; i++)
   {
      pointer[i] = pointer[i-1] 
      + cumulativeFitness[populationSize-1]/(double)numberOfSelectedIndividuals;
   }

   // Selection vector

   int countNumberOfSelectedIndividuals = 0;

   if(pointer[0] <= cumulativeFitness[0])
   {
      selection[0] = true;
      countNumberOfSelectedIndividuals++;
   }

   for(int i = 0; i < numberOfSelectedIndividuals; i++)
   {
      for(int j = 1; j < populationSize; j++)
      {
         if(pointer[i] <= cumulativeFitness[j] && pointer[i] > cumulativeFitness[j-1])
         {
            selection[j] = true;
            countNumberOfSelectedIndividuals++;
         }
      }
   }

   // Number of selected individuals control sentence 

   if(countNumberOfSelectedIndividuals != numberOfSelectedIndividuals)
   {
      std::cerr << std::endl
                << "Flood Error: EvolutionaryAlgorithm class." << std::endl
                << "void performStochasticUniversalSamplingSelection(void) method." << std::endl
                << "Count number of selected individuals is not equal to number of selected individuals." 
                << std::endl << std::endl;

      exit(1);
   }
}


// void performIntermediateRecombination(void) method

/// This method performs inediate recombination between pairs of selected individuals to generate a new 
/// population. 
/// Each selected individual is to be recombined with two other selected individuals chosen at random. 
/// Results are stored in the population matrix.
///
/// @see recombinationSize.
///
/// @see performLineRecombination(void).
/// @see train(void).

void EvolutionaryAlgorithm::performIntermediateRecombination(void)
{
   MultilayerPerceptron* multilayerPerceptron = objectiveFunctional->getMultilayerPerceptron();
     
   int numberOfFreeParameters = multilayerPerceptron->getNumberOfFreeParameters();

   Matrix<double> newPopulation(populationSize, numberOfFreeParameters);

   int count = 0;
   for(int i = 0; i < populationSize; i++)
   {
      if(selection[i] == true)
      count ++;        
   }

   Vector<double> parent1(numberOfFreeParameters);
   Vector<double> parent2(numberOfFreeParameters);

   Vector<double> offspring(numberOfFreeParameters);

   Matrix<int> recombination(populationSize, 2);

   // Start recombination   

   int countNewPopulationSize = 0;

   for(int i = 0; i < populationSize; i++)
   {
      if(selection[i] == true)
      {
         // Set parent 1

         parent1 = getIndividual(i);

         // Generate 2 offspring with parent 1

         for(int j = 0; j < 2; j++)
         {
            // Choose parent 2 at random among selected individuals   

            bool parent2Candidate = false;

            do{
               // Integer random number beteen 0 and population size

               double random = (double)rand()/(RAND_MAX+1.0);

               int parent2CandidateIndex = (int)(populationSize*random);

               // Check if candidate for parent 2 is ok

               if(selection[parent2CandidateIndex] == true && parent2CandidateIndex != i)
               {
                  parent2Candidate = true;

                  recombination[countNewPopulationSize][0] = i;

                  recombination[countNewPopulationSize][1] = parent2CandidateIndex;

                  parent2 = getIndividual(parent2CandidateIndex);

                  // Perform inediate recombination between parent 1 and parent 2

                  for(int j = 0; j < numberOfFreeParameters; j++)
                  {
                     // Choose the scaling factor to be a random number between
                     // -recombinationSize and 1+recombinationSize for each
                     // variable anew.

                     double random = (double)rand()/(RAND_MAX+1.0);

                     double scalingFactor = -1.0*recombinationSize + (1.0 + recombinationSize)*random;

                     offspring[j] = scalingFactor*parent1[j] + (1.0 - scalingFactor)*parent2[j];
                  }

                  // Add offspring to newPopulation matrix

                  newPopulation.setRow(countNewPopulationSize, offspring);   
                  
                  countNewPopulationSize++;
               }
            }while(parent2Candidate != true);
         }
      }
   }

   // Count number of new individuals control sentence

   if(countNewPopulationSize != populationSize)
   {
      std::cerr << std::endl 
                << "Flood Error: EvolutionaryAlgorithm class." << std::endl
                << "void performIntermediateRecombination(void) method." << std::endl
                << "Count new population size is not equal to population size." 
                << std::endl
                << std::endl;

      exit(1);
   }

   // Set new population

   population = newPopulation;
}


// void performLineRecombination(void) method

/// This method performs line recombination between pairs of selected individuals to generate a new population. 
/// Each selected individual is to be recombined with two other selected individuals chosen at random. 
/// Results are stored in the population matrix.
///
/// @see recombinationSize.
///
/// @see performIntermediateRecombination(void).
/// @see train(void).

void EvolutionaryAlgorithm::performLineRecombination(void)
{
   MultilayerPerceptron* multilayerPerceptron = objectiveFunctional->getMultilayerPerceptron();
     
   int numberOfFreeParameters = multilayerPerceptron->getNumberOfFreeParameters();

   Matrix<double> newPopulation(populationSize, numberOfFreeParameters);

   Vector<double> parent1(numberOfFreeParameters);
   Vector<double> parent2(numberOfFreeParameters);

   Vector<double> offspring(numberOfFreeParameters);

   Matrix<int> recombination(populationSize, 2);

   // Start recombination   

   int countNewPopulationSize = 0;

   for(int i = 0; i < populationSize; i++)
   {
      if(selection[i] == true)
      {
         // Set parent 1

         parent1 = getIndividual(i);

         // Generate 2 offspring with parent 1

         for(int j = 0; j < 2; j++)
         {
            // Choose parent 2 at random among selected individuals   

            bool parent2Candidate = false;

            do
            {
               // Integer random number beteen 0 and population size

               double random = (double)rand()/(RAND_MAX + 1.0);

               int parent2CandidateIndex = (int)(populationSize*random);

               // Check if candidate for parent 2 is ok

               if(selection[parent2CandidateIndex] == true && parent2CandidateIndex != i)
               {
                  parent2Candidate = true;

                  recombination[countNewPopulationSize][0] = i;
                  recombination[countNewPopulationSize][1] = parent2CandidateIndex;

                  parent2 = getIndividual(parent2CandidateIndex);

                  // Perform inediate recombination between parent 1
                  // and parent 2

                  // Choose the scaling factor to be a random number between
                  // -recombinationSize and 1+recombinationSize for all
                  // variables.

                  double random = (double)rand()/(RAND_MAX+1.0);

                  double scalingFactor = -1.0*recombinationSize 
                  + (1.0 + recombinationSize)*random;

                  offspring = parent1*scalingFactor + parent2*(1.0 - scalingFactor);

                  // Add offspring to newPopulation matrix

                  newPopulation.setRow(countNewPopulationSize, offspring);   

                  countNewPopulationSize++;
               }
            }while(parent2Candidate == false);
         }
      }
   }

   // Count new population size control sentence

   if(countNewPopulationSize != populationSize)
   {
      std::cerr << std::endl
                << "Flood Error: EvolutionaryAlgorithm class." << std::endl
                << "void performLineRecombination(void) method." << std::endl
                << "Count new population size is not equal to population size." << std::endl 
                << std::endl;

      exit(1);
   }

   // Set new population

   population = newPopulation;
}


// void performNormalMutation(void) method

/// This method performs normal mutation to all individuals in order to generate a new population. 
/// It uses the Box-Muller transformation to generate random numbers with normal distribution. 
/// Results are stored in the population matrix.
///
/// @see mutationRate.
/// @see mutationRange.
///
/// @see performUniformMutation(void).
/// @see train(void).


void EvolutionaryAlgorithm::performNormalMutation(void)
{
   const double pi = 3.141592654;

   MultilayerPerceptron* multilayerPerceptron = objectiveFunctional->getMultilayerPerceptron();
     
   int numberOfFreeParameters = multilayerPerceptron->getNumberOfFreeParameters();

   Vector<double> individual(numberOfFreeParameters);

   for(int i = 0; i < populationSize; i++)
   {
      individual = getIndividual(i);

      for(int j = 0; j < numberOfFreeParameters; j++)
      {
         // Random number between 0 and 1

         double pointer = (double)rand()/(RAND_MAX+1.0);

         if(pointer < mutationRate)
         {
            // Random numbers between 0 and 1

            double random1 = 0.0;
            
            do // random1 cannot be zero
            {
               random1 = (double)rand()/(RAND_MAX+1.0);
            
            }while(random1 == 0);

            double random2 = (double)rand()/(RAND_MAX+1.0);

            // Box-Muller transformation

            double mean = 0.0;
            double standardDeviation = mutationRange;

            double normallyDistributedRandomNumber 
            = mean + sqrt(-2.0*log(random1))*sin(2.0*pi*random2)*standardDeviation;

            individual[j] += normallyDistributedRandomNumber;
         }
      }

      setIndividual(i, individual);
   }
}  


// void performUniformMutation(void) method

/// This method performs uniform mutation to all individuals in order to generate a new population. 
/// Results are stored in the population matrix.
///
/// @see mutationRate.
/// @see mutationRange.
///
/// @see performNormalMutation(void).
/// @see train(void).

void EvolutionaryAlgorithm::performUniformMutation(void)
{
   MultilayerPerceptron* multilayerPerceptron = objectiveFunctional->getMultilayerPerceptron();
     
   int numberOfFreeParameters = multilayerPerceptron->getNumberOfFreeParameters();
   
   Vector<double> individual(numberOfFreeParameters, 0.0);

   for(int i = 0; i < populationSize; i++)
   {
      individual = getIndividual(i);

      for(int j = 0; j < numberOfFreeParameters; j++)
      {
         // random number between 0 and 1

         double pointer = (double)rand()/(RAND_MAX+1.0);

         if(pointer < mutationRate)
         {
            // random number between 0 and 1

            double random = (double)rand()/(RAND_MAX+1.0);

            double uniformlyDistributedRandomNumber
            = (-1.0 + 2.0*random)*mutationRange;

            individual[j] += uniformlyDistributedRandomNumber;
         }
      }

      setIndividual(i, individual);
   }
}


// void train(void) method

/// This method trains a multilayer perceptron with an associated
/// objective function according to the evolutionary algorithm.
/// Training occurs according to the training operators and their related
/// parameters.
///
/// @see FitnessAssignmentMethod.
/// @see SelectionMethod.
/// @see RecombinationMethod.
/// @see MutationMethod.

void EvolutionaryAlgorithm::train(void)
{
   // Control sentence (if debug)

   #ifndef NDEBUG 

   if(objectiveFunctional == NULL)
   {
      std::cout << std::endl
                << "Flood Error: EvolutionaryAlgorithm class." << std::endl
                << "void train(void) method." << std::endl
                << "Pointer to objective functional object cannot be NULL." << std::endl
                << std::endl;

      exit(1);
   }

   #endif

   // Multilayer perceptron

   MultilayerPerceptron* multilayerPerceptron = objectiveFunctional->getMultilayerPerceptron();

   int numberOfFreeParameters = multilayerPerceptron->getNumberOfFreeParameters();   

   // Training time stuff

   time_t beginningTime, currentTime;

   double elapsedTime = 0.0;

   // Start training

   time(&beginningTime);

   if(display)
   {
      std::cout << std::endl
                << "Training with the evolutionary algorithm..." 
                << std::endl;
   }

   // Evaluation of population

   evaluatePopulation();

   // Check for best individual

   Vector<double> bestIndividual(numberOfFreeParameters, 0.0);

   double bestNorm = 0.0;

   double bestEvaluation = 1.0e99;

   for(int i = 0; i < populationSize; i++)
   {
      if(evaluation[i] < bestEvaluation)
      {
         bestIndividual = getIndividual(i);

         bestNorm = bestIndividual.calculateNorm();

         bestEvaluation = evaluation[i];       
         
         // Set best individual free paramterers to multilayer perceptron

         multilayerPerceptron->setFreeParameters(bestIndividual);        
      }
   }

   resizeTrainingHistory(maximumNumberOfGenerations+1);

   // Elapsed time

   time(&currentTime);
   elapsedTime = difftime(currentTime, beginningTime);

   if(reserveElapsedTimeHistory)
   {
      elapsedTimeHistory[0] = elapsedTime;
   }

   // Population

   if(reservePopulationHistory)
   {
      populationHistory[0] = population; 
   }

   Vector<double> populationNorm = calculatePopulationNorm();

   // Mean norm 

   double meanNorm = populationNorm.calculateMean();

   if(reserveMeanNormHistory)
   {
      meanNormHistory[0] = meanNorm;
   }

   // Standard deviation of norm

   double standardDeviationNorm = populationNorm.calculateStandardDeviation();

   if(reserveMeanNormHistory)
   {
      standardDeviationNormHistory[0] = standardDeviationNorm;                                
   }

   // Best individual norm

   if(reserveBestNormHistory)
   {
      bestNormHistory[0] = bestNorm;
   }

   // Mean evaluation

   double meanEvaluation = evaluation.calculateMean();

   if(reserveMeanEvaluationHistory)
   {
      meanEvaluationHistory[0] = meanEvaluation;
   }

   // Standard deviation of evaluation

   double standardDeviationEvaluation = evaluation.calculateStandardDeviation();

   if(reserveStandardDeviationEvaluationHistory)
   {
      standardDeviationEvaluationHistory[0] = standardDeviationEvaluation;
   }

   // Best individual evaluation

   if(reserveBestEvaluationHistory)
   {
      bestEvaluationHistory[0] = bestEvaluation;
   }

   // Stopping criteria

   if(bestEvaluation <= evaluationGoal)
   {          
      if(display)
      {
         std::cout << std::endl
                   << "Initial evaluation is less than goal." << std::endl;

         std::cout << "Elapsed time: " << elapsedTime << ";" << std::endl;

         std::cout << "Initial mean norm: " << meanNorm << std::endl;
         std::cout << "Initial standard deviation of norm: " << standardDeviationNorm << std::endl;
         std::cout << "Initial best norm: " << bestNorm << std::endl;

         std::cout << "Initial mean evaluation: " << meanEvaluation << std::endl;
         std::cout << "Initial standard deviation of evaluation: " << standardDeviationEvaluation << std::endl;
         std::cout << "Initial best evaluation: " << bestEvaluation << std::endl;

         objectiveFunctional->print();   
      }

      resizeTrainingHistory(1);

      return;   
   }
   else if(elapsedTime >= maximumTime)
   {
      if(display)
      {
         std::cout << std::endl
                   << "Maximum training time reached." << std::endl;

         std::cout << "Elapsed time: " << elapsedTime << ";" << std::endl;

         std::cout << "Initial mean norm: " << meanNorm << std::endl;
         std::cout << "Initial standard deviation of norm: " << standardDeviationNorm << std::endl;
         std::cout << "Initial best norm: " << bestNorm << std::endl;

         std::cout << "Initial mean evaluation: " << meanEvaluation << std::endl;
         std::cout << "Initial standard deviation of evaluation: " << standardDeviationEvaluation << std::endl;
         std::cout << "Initial best evaluation: " << bestEvaluation << std::endl;

         objectiveFunctional->print();   
      }

      resizeTrainingHistory(1);

      return;      
   }
   else
   {
      if(display)
      {
         std::cout << "Elapsed time: " << elapsedTime << ";" << std::endl;

         // Print multilayer perceptron stuff

         std::cout << "Mean norm: " << meanNorm << std::endl;
         std::cout << "Standard deviation of norm: " << standardDeviationNorm << std::endl;
         std::cout << "Best norm: " << bestNorm << std::endl;

         std::cout << "Mean evaluation: " << meanEvaluation << std::endl;
         std::cout << "Standard deviation of evaluation: " << standardDeviationEvaluation << std::endl;
         std::cout << "Best evaluation: " << bestEvaluation << std::endl;

         objectiveFunctional->print();
      }
   }

   // Main loop

   for(int generation = 1; generation <= maximumNumberOfGenerations; generation++)
   {
      // Generate new population  

      // Fitness assignment

      switch(fitnessAssignmentMethod)
      {
         case LinearRanking:
         { 
            performLinearRankingFitnessAssignment();
         }

         break;
      }

      // Selection

      switch(selectionMethod)
      {
         case RouletteWheel:
         {
            performRouletteWheelSelection();
         }

         break;

         case StochasticUniversalSampling:
         {
            performStochasticUniversalSamplingSelection();
         }

         break;
      }

      // Recombination

      switch(recombinationMethod)
      {
         case Intermediate:
         {
            performIntermediateRecombination();
         }

         break;

         case Line:
         {
            performLineRecombination();
         } 

         break;
      }

      // Mutation

      switch(mutationMethod)
      {
         case Normal:
         {
            performNormalMutation();
         }

         break;

         case Uniform:
         {
            performUniformMutation();
         }

         break;
      }

      // Evaluation of new population

      evaluatePopulation();

      // Check for best individual

      for(int i = 0; i < populationSize; i++)
      {
         if(evaluation[i] < bestEvaluation)
         {
            bestIndividual = getIndividual(i);

            bestNorm = bestIndividual.calculateNorm();

            bestEvaluation = evaluation[i];
            
            // Set best individual free paramterers to multilayer perceptron

            multilayerPerceptron->setFreeParameters(bestIndividual);              
         }
      }

      // Elapsed time

      elapsedTime = difftime(currentTime, beginningTime);

      if(reserveElapsedTimeHistory)
      {
         elapsedTimeHistory[generation] = elapsedTime;
      }

      // Population

     if(reservePopulationHistory)
     {
        populationHistory[generation] = population;
     }

     populationNorm = calculatePopulationNorm();

      // Mean norm 

      double meanNorm = populationNorm.calculateMean();

      if(reserveMeanNormHistory)
      {
         meanNormHistory[generation] = meanNorm;
      }

      // Standard deviation of norm

      double standardDeviationNorm = populationNorm.calculateStandardDeviation();

      if(reserveMeanNormHistory)
      {
         standardDeviationNormHistory[generation] = standardDeviationNorm;                                
      }
 
      // Best individual norm

      if(reserveBestNormHistory)
      {
         bestNormHistory[generation] = bestNorm;
      }

      // Mean evaluation

      double meanEvaluation = evaluation.calculateMean();

      if(reserveMeanEvaluationHistory)
      {
         meanEvaluationHistory[generation] = meanEvaluation;
      }

      // Standard deviation of evaluation

      double standardDeviationEvaluation = evaluation.calculateStandardDeviation();

      if(reserveStandardDeviationEvaluationHistory)
      {
         standardDeviationEvaluationHistory[generation] = standardDeviationEvaluation;
      }

      // Best individual evaluation

      if(reserveBestEvaluationHistory)
      {
         bestEvaluationHistory[generation] = bestEvaluation;
      }

      // Stopping criteria

      if(bestEvaluation <= evaluationGoal)
      {
         if(display)
         {
            std::cout << std::endl
                      << "Generation " << generation << ": "
                      << "Evaluation goal reached." << std::endl;

            std::cout << "Elapsed time: " << elapsedTime << ";" << std::endl;


            std::cout << "Mean free norm: " << meanNorm << std::endl;
            std::cout << "Standard deviation of norm: " << standardDeviationNorm << std::endl;
            std::cout << "Best norm: " << bestNorm << std::endl;

            std::cout << "Mean evaluation: " << meanEvaluation << std::endl;
            std::cout << "Standard deviation of evaluation: " << standardDeviationEvaluation << std::endl;
            std::cout << "Best evaluation: " << bestEvaluation << std::endl;                   

            objectiveFunctional->print();
         }
         
         resizeTrainingHistory(generation);

         break;
      }
      else if(elapsedTime >= maximumTime)
      {
         if(display)
         {
            std::cout << std::endl
                      << "Generation " << generation << ": "
                      << "Maximum training time reached." << std::endl;

            std::cout << "Elapsed time: " << elapsedTime << ";" << std::endl;

            std::cout << "Mean norm: " << meanNorm << std::endl;
            std::cout << "Standard deviation of norm: " << standardDeviationNorm << std::endl;
            std::cout << "Best norm: " << bestNorm << std::endl;

            std::cout << "Mean evaluation: " << meanEvaluation << std::endl;
            std::cout << "Standard deviation of evaluation: " << standardDeviationEvaluation  << std::endl;
            std::cout << "Best evaluation: " << bestEvaluation << std::endl;

            objectiveFunctional->print();
         }

         resizeTrainingHistory(generation);

         break;
      }
      else if(generation == (maximumNumberOfGenerations - 1))
      {
         if(display)
         {
            std::cout << std::endl
                      << "Generation " << generation << ": "
                      << "Maximum number of generations reached." << std::endl;

            std::cout << "Elapsed time: " << elapsedTime << ";" << std::endl;

            std::cout << "Mean norm: " << meanNorm << std::endl;
            std::cout << "Standard deviation of norm: " << standardDeviationNorm << std::endl;
            std::cout << "Best norm: " << bestNorm << std::endl;

            std::cout << "Mean evaluation: " << meanEvaluation << std::endl;
            std::cout << "Standard deviation of evaluation: " << standardDeviationEvaluation << std::endl;
            std::cout << "Best evaluation: " << bestEvaluation << std::endl;                  

            objectiveFunctional->print();
         }

         break;
      }

      // Progress

      if(display && generation % displayPeriod == 0)
      {
         std::cout << std::endl
                   << "Generation " << generation << "; " << std::endl;

         std::cout << "Elapsed time: " << elapsedTime << ";" << std::endl;

         std::cout << "Mean norm: " << meanNorm << std::endl;
         std::cout << "Standard deviation of norm: " << standardDeviationNorm << std::endl;
         std::cout << "Best norm: " << bestNorm << std::endl;

         std::cout << "Mean evaluation: " << meanEvaluation << std::endl;
         std::cout << "Standard deviation of evaluation: " << standardDeviationEvaluation << std::endl;
         std::cout << "Best evaluation: " << bestEvaluation << std::endl;

         objectiveFunctional->print();
      }

      // Reset selection vector

      Vector<bool> newSelection(populationSize, false);

      selection = newSelection;
   }
}


// void print(void) method

/// This method prints to the screen the members of the evolutionary algorithm object.
///
/// Training operators:
/// <ul>
/// <li> Fitness assignment method.
/// <li> Selection method.
/// <li> Recombination method.
/// <li> Mutation method.
/// </ul>
///
/// Training parameters:
/// <ul>
/// <li> Population size.
/// <li> Selective pressure.
/// <li> Recombination size.
/// <li> Mutation rate.
/// <li> Mutation range.
/// </ul>
///
/// Stopping criteria:
/// <ul> 
/// <li> Evaluation goal.
/// <li> Maximum time.
/// <li> Maximum number of generations. 
/// </ul> 
///  
/// User stuff: 
/// <ul>
/// <li> Display. 
/// <li> Display period. 
/// <li> Reserve elapsed time.
/// <li> Reserve mean norm history.
/// <li> Reserve standard deviation of norm history.
/// <li> Reserve best norm history.
/// <li> Reserve mean evaluation history.
/// <li> Reserve standard deviation of evaluation history.
/// <li> Reserve best evaluation history.
/// </ul>
///
/// Population matrix. 

void EvolutionaryAlgorithm::print(void)
{
   int numberOfFreeParameters = 0;   
   
   if(objectiveFunctional != NULL)
   {
      MultilayerPerceptron* multilayerPerceptron = objectiveFunctional->getMultilayerPerceptron();

      numberOfFreeParameters = multilayerPerceptron->getNumberOfFreeParameters();
   }

   std::cout << std::endl
             << "Flood Neural Network. Evolutionary Algorithm Object." << std::endl;

   std::cout << "Population size: " << std::endl
             << populationSize << std::endl
             << "Number of free parameters: " << std::endl
             << numberOfFreeParameters  << std::endl;

   // Training operators

   // Fitness assingment method

   std::cout << "Fitness assignment method:" << std::endl;

   if(fitnessAssignmentMethod == LinearRanking)
   {
      std::cout << "Linear ranking" << std::endl;
   }

   // Selection method

   std::cout << "Selection method:" << std::endl;

   if(selectionMethod == RouletteWheel)
   {
      std::cout << "Roulette wheel" << std::endl;
   }
   else if(selectionMethod == StochasticUniversalSampling)
   {
      std::cout << "Stochastic universal sampling" << std::endl;
   }

   // Recombination method

   std::cout << "Recombination method:" << std::endl;

   if(recombinationMethod == Line)
   {
      std::cout << "Line" << std::endl;
   }
   else if(recombinationMethod == Intermediate)
   {
      std::cout << "Intermediate" << std::endl;
   }

   // Mutation method

   std::cout << "Mutation method:" << std::endl;

   if(mutationMethod == Normal)
   {
      std::cout << "Normal" << std::endl;
   }
   else if(mutationMethod == Uniform)
   {
      std::cout << "Uniform" << std::endl;
   }


   // Training parameters

   std::cout << "Selective pressure: " << std::endl
             << selectivePressure << std::endl
             << "Recombination size: " << std::endl
             << recombinationSize << std::endl
             << "Mutation rate: " << std::endl
             << mutationRate << std::endl
             << "Mutation range: " << std::endl
             << mutationRange << std::endl;

   // Stopping criteria

   std::cout << "Evaluation goal: " << std::endl
             << evaluationGoal << std::endl
             << "Maximum time: " << std::endl
             << maximumTime << std::endl
             << "Maximum number of generations: " << std::endl
             << maximumNumberOfGenerations << std::endl;

   // User stuff

   std::cout << "Display: " << std::endl
             << display << std::endl
             << "Display period: " << std::endl
             << displayPeriod << std::endl;
   
   std::cout << "Reserve elapsed time history:" << std::endl
             << reserveElapsedTimeHistory << std::endl
             << "Reserve mean norm history:" << std::endl
             << reserveMeanNormHistory << std::endl
             << "Reserve standard deviation norm history:" << std::endl
             << reserveStandardDeviationNormHistory << std::endl
             << "Reserve best norm history:" << std::endl
             << reserveBestNormHistory << std::endl
             << "Reserve mean evaluation history:" << std::endl
             << reserveMeanEvaluationHistory << std::endl
             << "Reserve standard deviation evaluation history:" << std::endl
             << reserveStandardDeviationEvaluationHistory << std::endl
             << "Reserve best evaluation history:" << std::endl
             << reserveBestEvaluationHistory << std::endl;

   std::cout << "Population:" << std::endl;

   Vector<double> individual(numberOfFreeParameters);

   for(int i = 0; i < populationSize; i++)
   {
      individual = getIndividual(i);

      std::cout << "Individual " << i << ":" << std::endl
                << individual << std::endl;     
   }
}


// void save(char*) method

/// This method saves the evolutionary algorithm object to a data file. 
///
/// Training operators:
/// <ul>
/// <li> Fitness assignment method.
/// <li> Selection method.
/// <li> Recombination method.
/// <li> Mutation method.
/// </ul>
///
/// Training parameters:
/// <ul>
/// <li> Population size.
/// <li> Selective pressure.
/// <li> Recombination size.
/// <li> Mutation rate.
/// <li> Mutation range.
/// </ul>
///
/// Stopping criteria:
/// <ul> 
/// <li> Evaluation goal.
/// <li> Maximum time.
/// <li> Maximum number of generations. 
/// </ul> 
///  
/// User stuff: 
/// <ul>
/// <li> Display. 
/// <li> Display period. 
/// <li> Reserve elapsed time history.
/// <li> Reserve mean norm history.
/// <li> Reserve standard deviation of norm history.
/// <li> Reserve best norm history.
/// <li> Reserve mean evaluation history.
/// <li> Reserve standard deviation of evaluation history.
/// <li> Reserve best evaluation history.
/// </ul>
///
/// @param filename Filename.
///
/// @see load(char*).

void EvolutionaryAlgorithm::save(char* filename)
{
   int numberOfFreeParameters = 0;   
   
   if(objectiveFunctional != NULL)
   {
      MultilayerPerceptron* multilayerPerceptron = objectiveFunctional->getMultilayerPerceptron();

      numberOfFreeParameters = multilayerPerceptron->getNumberOfFreeParameters();
   }


   // File

   std::fstream file;

   file.open(filename, std::ios::out);

   if(!file.is_open())
   {
      std::cerr << std::endl 
                << "Flood Error: EvolutionaryAlgorithm class." << std::endl
                << "void save(char*) method." << std::endl
                << "Cannot open evolutionary algorithm object data file."  << std::endl
                << std::endl;

      exit(1);
   }
   else
   {
      if(display)
      {   
         std::cout << std::endl
                   << "Saving evolutionary algorithm object to data file..." << std::endl;
      }
   }

   // Write file header


   file << "% Flood Neural Network. Evolutionary Algorithm Object." << std::endl;

   file << "PopulationSize:" << std::endl
        << populationSize << std::endl
        << "NumberOfFreeParameters:" << std::endl
        << numberOfFreeParameters  << std::endl;

   // Training operators

   // Fitness assingment method

   file << "FitnessAssignmentMethod:" << std::endl;

   if(fitnessAssignmentMethod == LinearRanking)
   {
      file << "LinearRanking" << std::endl;
   }

   // Selection method

   file << "SelectionMethod:" << std::endl;

   if(selectionMethod == RouletteWheel)
   {
      file << "RouletteWheel" << std::endl;
   }
   else if(selectionMethod == StochasticUniversalSampling)
   {
      file << "StochasticUniversalSampling" << std::endl;
   }

   // Recombination method

   file << "RecombinationMethod:" << std::endl;

   if(recombinationMethod == Line)
   {
      file << "Line" << std::endl;
   }
   else if(recombinationMethod == Intermediate)
   {
      file << "Intermediate" << std::endl;
   }

   // Mutation method

   file << "MutationMethod:" << std::endl;

   if(mutationMethod == Normal)
   {
      file << "Normal" << std::endl;
   }
   else if(mutationMethod == Uniform)
   {
      file << "Uniform" << std::endl;
   }

   // Training parameters

   file << "SelectivePressure:" << std::endl
        << selectivePressure << std::endl
        << "RecombinationSize:" << std::endl
        << recombinationSize << std::endl
        << "MutationRate:" << std::endl
        << mutationRate << std::endl
        << "MutationRange: " << std::endl
        << mutationRange << std::endl;

   // Stopping criteria

   file << "EvaluationGoal: " << std::endl
        << evaluationGoal << std::endl
        << "MaximumTime: " << std::endl
        << maximumTime << std::endl
        << "MaximumNumberOfGenerations: " << std::endl
        << maximumNumberOfGenerations << std::endl;

   // User stuff

   file << "Display: " << std::endl
        << display << std::endl
        << "DisplayPeriod: " << std::endl
        << displayPeriod << std::endl;

   file << "ReserveElapsedTimeHistory:" << std::endl
        << reserveElapsedTimeHistory << std::endl
        << "ReserveMeanNormHistory:" << std::endl
        << reserveMeanNormHistory << std::endl
        << "ReserveStandardDeviationNormHistory:" << std::endl
        << reserveStandardDeviationNormHistory << std::endl
        << "ReserveBestNormHistory:" << std::endl
        << reserveBestNormHistory << std::endl
        << "ReserveMeanEvaluationHistory:" << std::endl
        << reserveMeanEvaluationHistory << std::endl
        << "ReserveStandardDeviationEvaluationHistory:" << std::endl
        << reserveStandardDeviationEvaluationHistory << std::endl
        << "ReserveBestEvaluationHistory:" << std::endl
        << reserveBestEvaluationHistory << std::endl;

   // Population

   file << "Population:" << std::endl;

   for(int i = 0; i < populationSize; i++)
   {
      Vector<double> individual = getIndividual(i);

      file << "Individual" << i << ":" << std::endl; 

      file << individual << std::endl;

      file << std::endl;
   }

   file.close();
}


// void load(char*) method

/// This method loads a evolutionary algorithm object from a data file. 
/// Please mind about the file format, wich is specified in the User's Guide. 
///
/// Training operators:
/// <ul>
/// <li> Fitness assignment method.
/// <li> Selection method.
/// <li> Recombination method.
/// <li> Mutation method.
/// </ul>
///
/// Training parameters:
/// <ul>
/// <li> Population size.
/// <li> Selective pressure.
/// <li> Recombination size.
/// <li> Mutation rate.
/// <li> Mutation range.
/// </ul>
///
/// Stopping criteria:
/// <ul> 
/// <li> Evaluation goal.
/// <li> Maximum time.
/// <li> Maximum number of generations. 
/// </ul> 
///  
/// User stuff: 
/// <ul>
/// <li> Display. 
/// <li> Display period. 
/// <li> Reserve elapsed time history.
/// <li> Reserve mean norm history.
/// <li> Reserve standard deviation of norm history.
/// <li> Reserve best norm history.
/// <li> Reserve mean evaluation history.
/// <li> Reserve standard deviation of evaluation history.
/// <li> Reserve best evaluation history.
/// </ul>
///
/// @param filename Filename.
///
/// @see save(char*).

void EvolutionaryAlgorithm::load(char* filename)
{
   // Multilayer perceptron

   MultilayerPerceptron* multilayerPerceptron = objectiveFunctional->getMultilayerPerceptron();

   int numberOfFreeParameters = multilayerPerceptron->getNumberOfFreeParameters();

   // File

   std::fstream file;

   file.open(filename, std::ios::in);

   if(!file.is_open())
   {
      std::cerr << std::endl
                << "Flood Error: EvolutionaryAlgorithm class." << std::endl
                << "void load(char*) method." << std::endl
                << "Cannot open evolutionary algorithm object data file."  << std::endl
				<< std::endl;

      exit(1);
   }
   else
   {
      if(display)
      {
         std::cout << std::endl
                   << "Loading evolutionary algorithm object from data file..."  
                   << std::endl;
      }
   }

   // Load new number of individuals and new number of free parameters form file

   int newPopulationSize = 0;
   int newNumberOfFreeParameters = 0;

   std::string word;

   // Training parameters

   // Population size

   while(word != "PopulationSize:")
   {
      file >> word;
   }

   file >> newPopulationSize;

   if(newPopulationSize != populationSize)
   {
      std::cerr << std::endl
                << "Flood Error: EvolutionaryAlgorithm class." << std::endl
                << "void load(char*) method." << std::endl
                << "New population size is not equal to population size." << std::endl
                << std::endl;

      exit(1);
  }

   // Number of free parameters

   file >> word;

   file >> newNumberOfFreeParameters;

   if(newNumberOfFreeParameters != numberOfFreeParameters)
   {
      std::cerr << std::endl
                << "Flood Error: EvolutionaryAlgorithm class." << std::endl
                << "void load(char*) method." << std::endl
                << "New number of free parameters is not equal to number of free parameters."
                << std::endl << std::endl;

      exit(1);
   }


   // Training operators

   // Fitness assingment method

   file >> word;

   if(word != "FitnessAssingmentMethod:")
   {
      std::cerr << std::endl
                << "Flood Error: EvolutionaryAlgorithm class." << std::endl
                << "void load(char*) method." << std::endl
                << "Unknown file format."  << std::endl;

      exit(1);   
   }

   file >> word;

   if(word == "LinearRanking")
   {
      fitnessAssignmentMethod = LinearRanking;
   }
   else
   {
      std::cerr << std::endl
                << "Flood Error: EvolutionaryAlgorithm class." << std::endl
                << "void load(char*) method." << std::endl
                << "Unknown file format."  << std::endl;

      exit(1);   
   }

   // Selection method
 
   file >> word;

   if(word != "SelectionMethod:")
   {
      std::cerr << std::endl
                << "Flood Error: EvolutionaryAlgorithm class." << std::endl
                << "void load(char*) method." << std::endl
                << "Unknown file format."  << std::endl;

      exit(1);   
   }

   file >> word;

   if(word == "RouletteWheel")
   {
      selectionMethod = RouletteWheel;
   }
   else if(word == "StochasticUniversalSampling")
   {
      selectionMethod = StochasticUniversalSampling;
   } 
   else
   {
      std::cerr << std::endl
                << "Flood Error: EvolutionaryAlgorithm class." << std::endl
                << "void load(char*) method." << std::endl
                << "Unknown file format."  << std::endl;

      exit(1);   
   }

   // Recombination method

   file >> word;

   if(word != "RecombinationMethod:")
   {
      std::cerr << std::endl
                << "Flood Error: EvolutionaryAlgorithm class." << std::endl
                << "void load(char*) method." << std::endl
                << "Unknown file format."  << std::endl;

      exit(1);   
   }

   file >> word;

   if(word == "Line")
   {
      recombinationMethod = Line;
   }
   else if(word == "Intermediate")
   {
      recombinationMethod = Intermediate;
   } 
   else
   {
      std::cerr << std::endl
                << "Flood Error: EvolutionaryAlgorithm class." << std::endl
                << "void load(char*) method." << std::endl
                << "Unknown file format."  << std::endl;

      exit(1);   
   }

   // Mutation method

   file >> word;

   if(word != "MutationMethod:")
   {
      std::cerr << std::endl
                << "Flood Error: EvolutionaryAlgorithm class." << std::endl
                << "void load(char*) method." << std::endl
                << "Unknown file format."  << std::endl;

      exit(1);   
   }

   file >> word;

   if(word == "Normal")
   {
      mutationMethod = Normal;
   }
   else if(word == "Uniform")
   {
      mutationMethod = Uniform;
   } 
   else
   {
      std::cerr << std::endl
                << "Flood Error: EvolutionaryAlgorithm class." << std::endl
                << "void load(char*) method." << std::endl
                << "Unknown file format."  << std::endl;

      exit(1);   
   }

   // Training parameters

   // Selective pressure

   file >> word;

   if(word != "SelectivePressure:")
   {
      std::cerr << std::endl
                << "Flood Error: EvolutionaryAlgorithm class." << std::endl
                << "void load(char*) method." << std::endl
                << "Unknown file format. " << word << "???" <<std::endl
                << std::endl;

      exit(1);   
   }

   file >> selectivePressure;

   // Recombination size

   file >> word;

   if(word != "RecombinationSize")
   {
      std::cerr << std::endl
                << "Flood Error: EvolutionaryAlgorithm class." << std::endl
                << "void load(char*) method." << std::endl
                << "Unknown file format. " << word << "???" <<std::endl
                << std::endl;

      exit(1);   
   }

   file >> recombinationSize;

   // Mutation rate

   file >> word;

   if(word != "MutationRate:")
   {
      std::cerr << std::endl
                << "Flood Error: EvolutionaryAlgorithm class." << std::endl
                << "void load(char*) method." << std::endl
                << "Unknown file format. " << word << "???" <<std::endl
                << std::endl;

      exit(1);   
   }

   file >> mutationRate;

   // Mutation range
 
   file >> word;

   if(word != "MutationRange:")
   {
      std::cerr << std::endl
                << "Flood Error: EvolutionaryAlgorithm class." << std::endl
                << "void load(char*) method." << std::endl
                << "Unknown file format. " << word << "???" <<std::endl
                << std::endl;

      exit(1);   
   }

   file >> mutationRange;

   // Stopping criteria: 

   // Evaluation goal

   file >> word;

   if(word != "EvaluationGoal:")
   {
      std::cerr << std::endl
                << "Flood Error: EvolutionaryAlgorithm class." << std::endl
                << "void load(char*) method." << std::endl
                << "Unknown file format. " << word << "???" <<std::endl
                << std::endl;

      exit(1);   
   }

   file >> evaluationGoal;

   // Maximum time

   file >> word;

   if(word != "MaximumTime:")
   {
      std::cerr << std::endl
                << "Flood Error: EvolutionaryAlgorithm class." << std::endl
                << "void load(char*) method." << std::endl
                << "Unknown file format. " << word << "???" <<std::endl
                << std::endl;

      exit(1);   
   }

   file >> maximumTime;
   
   // Maximum number of generations

   file >> word;

   if(word != "MaximumNumberOfGenerations:")
   {
      std::cerr << std::endl
                << "Flood Error: EvolutionaryAlgorithm class." << std::endl
                << "void load(char*) method." << std::endl
                << "Unknown file format. " << word << "???" <<std::endl
                << std::endl;

      exit(1);   
   }

   file >> maximumNumberOfGenerations;

   // User stuff: 

   // Display

   file >> word;

   if(word != "Display:")
   {
      std::cerr << std::endl
                << "Flood Error: EvolutionaryAlgorithm class." << std::endl
                << "void load(char*) method." << std::endl
                << "Unknown file format. " << word << "???" <<std::endl
                << std::endl;

      exit(1);   
   }

   file >> display;

   // Display period

   file >> word;

   if(word != "DisplayPeriod:")
   {
      std::cerr << std::endl
                << "Flood Error: EvolutionaryAlgorithm class." << std::endl
                << "void load(char*) method." << std::endl
                << "Unknown file format. " << word << "???" <<std::endl
                << std::endl;

      exit(1);   
   }

   file >> displayPeriod;

   // Reserve elapsed time history

   file >> word;

   if(word != "ReserveElapsedTimeHistory:")
   {
      std::cerr << std::endl
                << "Flood Error: EvolutionaryAlgorithm class." << std::endl
                << "void load(char*) method." << std::endl
                << "Unknown file format. " << word << "???" <<std::endl
                << std::endl;

      exit(1);   
   }

   file >> reserveElapsedTimeHistory;

   // Reserve mean norm history

   file >> word;

   if(word != "ReserveMeanNormHistory:")
   {
      std::cerr << std::endl
                << "Flood Error: EvolutionaryAlgorithm class." << std::endl
                << "void load(char*) method." << std::endl
                << "Unknown file format. " << word << "???" <<std::endl
                << std::endl;

      exit(1);   
   }

   file >> reserveMeanNormHistory;

   // Reserve standard deviation norm history

   file >> word;

   if(word != "ReserveStandardDeviationNormHistory:")
   {
      std::cerr << std::endl
                << "Flood Error: EvolutionaryAlgorithm class." << std::endl
                << "void load(char*) method." << std::endl
                << "Unknown file format. " << word << "???" <<std::endl
                << std::endl;

      exit(1);   
   }

   file >> reserveStandardDeviationNormHistory;

   // Reserve best norm history

   file >> word;

   if(word != "ReserveBestNormHistory:")
   {
      std::cerr << std::endl
                << "Flood Error: EvolutionaryAlgorithm class." << std::endl
                << "void load(char*) method." << std::endl
                << "Unknown file format. " << word << "???" <<std::endl
                << std::endl;

      exit(1);   
   }

   file >> reserveBestNormHistory;

   // Reserve mean evaluation history

   file >> word;

   if(word != "ReserveMeanEvaluationHistory:")
   {
      std::cerr << std::endl
                << "Flood Error: EvolutionaryAlgorithm class." << std::endl
                << "void load(char*) method." << std::endl
                << "Unknown file format. " << word << "???" <<std::endl
                << std::endl;

      exit(1);   
   }

   file >> reserveMeanEvaluationHistory;

   // Reserve standard deviation evaluation history

   file >> word;

   if(word != "ReserveStandardDeviationEvaluationHistory:")
   {
      std::cerr << std::endl
                << "Flood Error: EvolutionaryAlgorithm class." << std::endl
                << "void load(char*) method." << std::endl
                << "Unknown file format. " << word << "???" <<std::endl
                << std::endl;

      exit(1);   
   }

   file >> reserveStandardDeviationEvaluationHistory;

   // Reserve best evaluation history

   file >> word;

   if(word != "ReserveBestEvaluationHistory:")
   {
      std::cerr << std::endl
                << "Flood Error: EvolutionaryAlgorithm class." << std::endl
                << "void load(char*) method." << std::endl
                << "Unknown file format. " << word << "???" <<std::endl
                << std::endl;

      exit(1);   
   }

   file >> reserveBestEvaluationHistory;

   // Population

   file >> word;

   for(int i = 0; i < populationSize; i++)
   {
      file >> word;

      for(int j = 0; j < numberOfFreeParameters; j++)
      {
         file >> population[i][j];
      }
   }

   // Close file

   file.close();
}


// void resizeTrainingHistory(int) method

/// This method resizes the vectors or matrices containing training history information 
/// to a new size:
///
/// <ul>
/// <li> Elapsed time history.
/// <li> Population history.
/// <li> Mean norm history.
/// <li> Standard deviation norm history.
/// <li> Best norm history.
/// <li> Mean evaluation history.
/// <li> Standard deviation evaluation history.
/// <li> Best evaluation history.
/// </ul>
///
/// @param newSize Size of training history. 

void EvolutionaryAlgorithm::resizeTrainingHistory(int newSize)
{
   // Elapsed time history vector

   if(reserveElapsedTimeHistory)
   {
      elapsedTimeHistory.resize(newSize);
   }

   // Population history

   if(reservePopulationHistory)
   {
      populationHistory.resize(newSize);
   }

   // Mean population norm history vector

   if(reserveMeanNormHistory)
   {
      meanNormHistory.resize(newSize);
   }

   // Standard deviation population norm history vector

   if(reserveStandardDeviationNormHistory)
   {
      standardDeviationNormHistory.resize(newSize);
   }

   // Best individual norm history vector

   if(reserveBestNormHistory)
   {
      bestNormHistory.resize(newSize);
   }

   // Mean evaluation history history vector

   if(reserveMeanEvaluationHistory)
   {
      meanEvaluationHistory.resize(newSize);
   }
 
   // Standard deviation of evaluation history vector

   if(reserveStandardDeviationEvaluationHistory)
   {
      standardDeviationEvaluationHistory.resize(newSize);
   }

   // Best evaluation history vector

   if(reserveBestEvaluationHistory)
   {
      bestEvaluationHistory.resize(newSize);
   }
}


// void saveTrainingHistory(char*) method

/// This method saves the training history to a data file. 
///
/// @param filename Training history filename. 

void EvolutionaryAlgorithm::saveTrainingHistory(char* filename)

{
   std::fstream file; 

   file.open(filename, std::ios::out);

   // Write file header 

   if(!file.is_open())
   {
      std::cerr << std::endl 
                << "Flood Error: Evolutionary Algorithm class." << std::endl
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

   file << "% Flood Neural Network. Evolutionary Training History." << std::endl;

   // Write file data

   if(reserveElapsedTimeHistory)
   {
      file << "ElapsedTimeHistory:" << std::endl;      file << elapsedTimeHistory << std::endl;      
   }
   if(reservePopulationHistory)
   {
      file << "PopulationHistory:" << std::endl;      file << populationHistory << std::endl;      
   }
   if(reserveMeanNormHistory)
   {
      file << "MeanNormHistory:" << std::endl;      file << meanNormHistory << std::endl;      
   }
   if(reserveStandardDeviationNormHistory)
   {
      file << "StandardDeviationNormHistory:" << std::endl;      file << standardDeviationNormHistory << std::endl;      
   }
   if(reserveBestNormHistory)
   {
      file << "BestNormHistory:" << std::endl;      file << bestNormHistory << std::endl;      
   }
   if(reserveMeanEvaluationHistory)
   {
      file << "MeanEvaluationHistory:" << std::endl;      file << meanEvaluationHistory << std::endl;      
   }
   if(reserveStandardDeviationEvaluationHistory)
   {
      file << "StandardDeviationEvaluationHistory:" << std::endl;      file << standardDeviationEvaluationHistory << std::endl;      
   }
   if(reserveBestEvaluationHistory)
   {
      file << "BestEvaluationHistory:" << std::endl;      file << bestEvaluationHistory << std::endl;      
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
