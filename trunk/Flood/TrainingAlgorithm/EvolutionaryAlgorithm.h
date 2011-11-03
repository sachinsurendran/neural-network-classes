/****************************************************************************************************************/
/*                                                                                                              */
/*   Flood: An Open Source Neural Networks C++ Library                                                          */
/*   www.cimne.com/flood                                                                                        */
/*                                                                                                              */
/*   E V O L U T I O N A R Y   A L G O R I T H M   C L A S S   H E A D E R                                      */
/*                                                                                                              */
/*   Roberto Lopez                                                                                              */
/*   International Center for Numerical Methods in Engineering (CIMNE)                                          */
/*   Technical University of Catalonia (UPC)                                                                    */
/*   Barcelona, Spain                                                                                           */
/*   E-mail: rlopez@cimne.upc.edu                                                                               */
/*                                                                                                              */
/****************************************************************************************************************/


#ifndef __EVOLUTIONARYALGORITHM_H__
#define __EVOLUTIONARYALGORITHM_H__


#include "TrainingAlgorithm.h"
#include "../ObjectiveFunctional/ObjectiveFunctional.h"

namespace Flood
{
 
/// This concrete class represents an evolutionary training algorithm for an objective functional of a multilayer 
/// perceptron.
///
/// @see TrainingAlgorithm.

class EvolutionaryAlgorithm : public TrainingAlgorithm
{

public:

   // ENUMERATIONS

   /// Avinash: Enumeration of the available evaluation methods.
   enum EvaluationMethod{SingleTrial, MultiTrial};

   /// Enumeration of the available training operators for fitness assignment.

   enum FitnessAssignmentMethod{LinearRanking, DecendingRanking};

   /// Enumeration of the available training operators for selection. 

   enum SelectionMethod{RouletteWheel, StochasticUniversalSampling, EliteSampling, None};

   /// Enumeration of the available training operators for recombination.

   enum RecombinationMethod{Line, Intermediate, Standard};

   /// Enumeration of the available training operators for mutation.

   enum MutationMethod{Normal, Uniform, OffspringsOnly};

private:

   // FIELDS

   // Population stuff

   /// Number of individuals in the population. 

   int populationSize;

   /// Population matrix.

   Matrix<double> population;

   /// Evaluation of population.

   /// Number of evaluations trials per individual
   int numberOfTrials;

   Vector<double> evaluation;

   /// Fitness of population.

   Vector<double> fitness;

   /// Selected individuals in population.

   Vector<bool> selection;

   /// Rank of individuals
   Vector<int> rank;
   
   // Training parameters

   /// Selective pressure. 
   /// Linear ranking allows values for the selective pressure between 1 and 2.

   double selectivePressure;

   /// Recombination size. 
   /// The recombination size value must be equal or greater than 0.

   /* Avinash: Percentage of top performers to chosen from the population for crossover.
               The crossover percentage must be value between 0 and 100 */
   double crossoverPercentage;

   double recombinationSize;

   /// Mutation rate.
   /// The mutation rate value must be between 0 and 1.

   double mutationRate;

   /// Mutation range.
   /// The mutation range value must be 0 or a positive number. 

   double mutationRange;

   /// Maximum number of generations to train.
   
   int maximumNumberOfGenerations;

   /// True if the population history, which is a vector of matrices, is to be reserved, false otherwise.
   /// Reserving the population history can be compuationally expensive if the number of free parameters, 
   /// the population size and the number of generations are big numbers. 

   bool reservePopulationHistory;

   /// True if the mean norm history vector is to be reserved, false otherwise.

   bool reserveMeanNormHistory;

   /// True if the standard deviation of norm history vector is to be reserved, false otherwise.

   bool reserveStandardDeviationNormHistory;

   /// True if the best norm history vector is to be reserved, false otherwise.

   bool reserveBestNormHistory;
   
   /// True if the mean evaluation history vector is to be reserved, false otherwise.

   bool reserveMeanEvaluationHistory;

   /// True if the standard deviation of evaluation history vector is to be reserved, false otherwise.

   bool reserveStandardDeviationEvaluationHistory;

   /// True if the best evaluation history vector is to be reserved, false otherwise.

   bool reserveBestEvaluationHistory;

   /// Population history over the epochs, which is a vector of matrices. 
   /// The element 0 of the vector is the initial population, and so on.

   Vector< Matrix<double> > populationHistory;

   /// Mean population norm training history.
   /// This vector is of the form (MeanOfPopulationNormAtEpoch0,...,MeanOfPopulationNormAtEpochN)

   Vector<double> meanNormHistory;

   /// Standard deviation of population norm training history.
   /// This vector is of the form 
   /// (StandardDeviationOfPopulationNormAtEpoch0,...,StandardDeviationOfPopulationNormAtEpochN)

   Vector<double> standardDeviationNormHistory;

   /// Best individual ever norm history.

   Vector<double> bestNormHistory;

   /// Mean evaluation training history.

   Vector<double> meanEvaluationHistory;

   /// Standard deviation of evaluation training history.

   Vector<double> standardDeviationEvaluationHistory;

   /// Best evaluation ever training history.

   Vector<double> bestEvaluationHistory;

   /// Evaluation method training operator.
   EvaluationMethod evaluationMethod;

   /// Fitness assignment training operators enumeration.

   FitnessAssignmentMethod fitnessAssignmentMethod;

   /// Selection training operators enumeration.

   SelectionMethod selectionMethod;

   /// Recombination training operators enumeration.

   RecombinationMethod recombinationMethod;

   /// Mutation training operators enumeration.

   MutationMethod mutationMethod;

   std::string save_to_filename;

public:

   // GENERAL CONSTRUCTOR

   EvolutionaryAlgorithm(ObjectiveFunctional*);


   // DEFAULT CONSTRUCTOR

   EvolutionaryAlgorithm(void);


   // DESTRUCTOR

   virtual ~EvolutionaryAlgorithm(void);


   // METHODS

   // Get methods

   int getPopulationSize(void);

   Matrix<double> getPopulation(void);

   Vector<double> getEvaluation(void);
   Vector<double> getFitness(void);
   Vector<bool> getSelection(void);

   double getSelectivePressure(void);
   double getRecombinationSize(void);
   double getMutationRate(void);
   double getMutationRange(void);

   double getMaximumNumberOfGenerations(void);

   bool getReservePopulationHistory(void);
   bool getReserveMeanNormHistory(void);
   bool getReserveStandardDeviationNormHistory(void);
   bool getReserveBestNormHistory(void);
   bool getReserveMeanEvaluationHistory(void);
   bool getReserveStandardDeviationEvaluationHistory(void);
   bool getReserveBestEvaluationHistory(void);

   Vector< Matrix<double> > getPopulationHistory(void);

   Vector<double> getMeanNormHistory(void);
   Vector<double> getStandardDeviationNormHistory(void);
   Vector<double> getBestNormHistory(void);

   Vector<double> getMeanEvaluationHistory(void);
   Vector<double> getStandardDeviationEvaluationHistory(void);
   Vector<double> getBestEvaluationHistory(void);

   FitnessAssignmentMethod getFitnessAssignmentMethod(void);
   SelectionMethod getSelectionMethod(void);
   RecombinationMethod getRecombinationMethod(void);
   MutationMethod getMutationMethod(void);

   // Set methods

   void setPopulationSize(int);

   void setPopulation(Matrix<double>);

   void setEvaluation(Vector<double>);
   void setFitness(Vector<double>);
   void setSelection(Vector<bool>);

   void setSelectivePressure(double);
   void setRecombinationSize(double);

   void setMutationRate(double);
   void setMutationRange(double);

   void setMaximumNumberOfGenerations(int);

   void setFitnessAssignmentMethod(FitnessAssignmentMethod);
   void setSelectionMethod(SelectionMethod);
   void setRecombinationMethod(RecombinationMethod);
   void setMutationMethod(MutationMethod);

   void setReservePopulationHistory(bool);
   void setReserveMeanNormHistory(bool);
   void setReserveStandardDeviationNormHistory(bool);
   void setReserveBestNormHistory(bool);
   void setReserveMeanEvaluationHistory(bool);
   void setReserveStandardDeviationEvaluationHistory(bool);
   void setReserveBestEvaluationHistory(bool);

   void setPopulationHistory(Vector< Matrix<double> >);

   void setMeanNormHistory(Vector<double>);
   void setStandardDeviationNormHistory(Vector<double>);
   void setBestNormHistory(Vector<double>);

   void setMeanEvaluationHistory(Vector<double>);
   void setStandardDeviationEvaluationHistory(Vector<double>);
   void setBestEvaluationHistory(Vector<double>);

   // Population methods

   Vector<double> getIndividual(int);
   void setIndividual(int, Vector<double>);

   void initPopulationUniform(void);
   void initPopulationUniform(double, double);
   void initPopulationUniform(Vector<double>, Vector<double>);
   void initPopulationUniform(Matrix<double>);

   void initPopulationNormal(void);
   void initPopulationNormal(double, double);
   void initPopulationNormal(Vector<double>, Vector<double>);
   void initPopulationNormal(Matrix<double>);

   Vector<double> calculatePopulationNorm(void);

   // Population evaluation methods

   void evaluatePopulation_singleTrial(void);
   void evaluatePopulation_multiTrial(void);

   // Fitness assignment methods

   void performLinearRankingFitnessAssignment(void);
   void performDecendingRanking(void);

   // Selection methods

   void performRouletteWheelSelection(void);
   void performStochasticUniversalSamplingSelection(void);
   void performEliteSelection(void);


   // Recombination methods

   void performIntermediateRecombination(void);
   void performLineRecombination(void);
   void performStandardRecombination(void);

   // Mutation methods

   void performNormalMutation(void);
   void performUniformMutation(void);
   void performOffspringsOnlyMutation(void);

   // Train methods

   void train(void);

   // Utility methods

   void resizeTrainingHistory(int);

   void print(void);

   void load(char*);
   void save(char*);

   void setReserveAllTrainingHistory(bool);

   void saveTrainingHistory(char*);

   void setSaveToFilename(char *);

};

}

#endif


// Flood: An Open Source Neural Networks C++ Flood.
// Copyright (C) 2005-2008 Roberto Lopez 
//
// This Flood is free software; you can redistribute it and/or
// modify it under the s of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or any later version.
//
// This Flood is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.

// You should have received a copy of the GNU Lesser General Public
// License along with this Flood; if not, write to the Free Software
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
