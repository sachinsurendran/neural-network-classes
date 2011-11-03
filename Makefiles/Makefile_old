#########################################################################
#                                                                       #
#   F L O O D   N E U R A L   N E T W O R K   M A K E F I L E           #
#                                                                       #
#   Roberto Lopez                                                       # 
#   International Center for Numerical Methods in Engineering (CIMNE)   #
#   Technical University of Catalonia (UPC)                             #
#   Barcelona, Spain                                                    #
#   E-mail: rlopez@cimne.upc.edu                                        #
#                                                                       #
#########################################################################

objects = TemplateApplication.o \
          Perceptron.o \
          MultilayerPerceptron.o \
          ObjectiveFunctional.o \
    	  SumSquaredError.o \
          MeanSquaredError.o \
          RootMeanSquaredError.o \
          NormalizedSquaredError.o \
          MinkowskiError.o \
          RegularizedMinkowskiError.o \
          BrachistochroneProblem.o \
          IsoperimetricProblem.o \
          CarProblem.o \
          MinimumDragProblem.o \
          TrainingAlgorithm.o \
          GradientDescent.o \
          ConjugateGradient.o \
          QuasiNewtonMethod.o \
    	  RandomSearch.o \
          EvolutionaryAlgorithm.o \
    	  InputTargetDataSet.o \
    	  LinearRegressionAnalysis.o \
    	  CorrectPredictionsAnalysis.o \
          IntegrationOfFunctions.o \
          OrdinaryDifferentialEquations.o 


Flood:  $(objects)
	g++ $(objects) -o Flood 


# Main function

TemplateApplication.o: ../Applications/TemplateApplication.cpp
	g++ -c ../Applications/TemplateApplication.cpp -o TemplateApplication.o


# Perceptron neuron model

Perceptron.o: ../Flood/Perceptron/Perceptron.cpp
	g++ -c ../Flood/Perceptron/Perceptron.cpp -o Perceptron.o


# Multilayer perceptron network architecture

MultilayerPerceptron.o: ../Flood/MultilayerPerceptron/MultilayerPerceptron.cpp
	g++ -c ../Flood/MultilayerPerceptron/MultilayerPerceptron.cpp -o MultilayerPerceptron.o


# Objective functional

ObjectiveFunctional.o: ../Flood/ObjectiveFunctional/ObjectiveFunctional.cpp
	g++ -c ../Flood/ObjectiveFunctional/ObjectiveFunctional.cpp -o ObjectiveFunctional.o

SumSquaredError.o: ../Flood/ObjectiveFunctional/SumSquaredError.cpp
	g++ -c ../Flood/ObjectiveFunctional/SumSquaredError.cpp -o SumSquaredError.o

MeanSquaredError.o: ../Flood/ObjectiveFunctional/MeanSquaredError.cpp
	g++ -c ../Flood/ObjectiveFunctional/MeanSquaredError.cpp -o MeanSquaredError.o

RootMeanSquaredError.o: ../Flood/ObjectiveFunctional/RootMeanSquaredError.cpp
	g++ -c ../Flood/ObjectiveFunctional/RootMeanSquaredError.cpp -o RootMeanSquaredError.o

NormalizedSquaredError.o: ../Flood/ObjectiveFunctional/NormalizedSquaredError.cpp
	g++ -c ../Flood/ObjectiveFunctional/NormalizedSquaredError.cpp -o NormalizedSquaredError.o

MinkowskiError.o: ../Flood/ObjectiveFunctional/MinkowskiError.cpp
	g++ -c ../Flood/ObjectiveFunctional/MinkowskiError.cpp -o MinkowskiError.o

RegularizedMinkowskiError.o: ../Flood/ObjectiveFunctional/RegularizedMinkowskiError.cpp
	g++ -c ../Flood/ObjectiveFunctional/RegularizedMinkowskiError.cpp -o RegularizedMinkowskiError.o

BrachistochroneProblem.o: ../Flood/ObjectiveFunctional/BrachistochroneProblem.cpp
	g++ -c ../Flood/ObjectiveFunctional/BrachistochroneProblem.cpp -o BrachistochroneProblem.o

IsoperimetricProblem.o: ../Flood/ObjectiveFunctional/IsoperimetricProblem.cpp
	g++ -c ../Flood/ObjectiveFunctional/IsoperimetricProblem.cpp -o IsoperimetricProblem.o

CarProblem.o: ../Flood/ObjectiveFunctional/CarProblem.cpp
	g++ -c ../Flood/ObjectiveFunctional/CarProblem.cpp -o CarProblem.o

MinimumDragProblem.o: ../Flood/ObjectiveFunctional/MinimumDragProblem.cpp
	g++ -c ../Flood/ObjectiveFunctional/MinimumDragProblem.cpp -o MinimumDragProblem.o


# Training algorithm

TrainingAlgorithm.o: ../Flood/TrainingAlgorithm/TrainingAlgorithm.cpp
	g++ -c ../Flood/TrainingAlgorithm/TrainingAlgorithm.cpp -o TrainingAlgorithm.o

GradientDescent.o: ../Flood/TrainingAlgorithm/GradientDescent.cpp
	g++ -c ../Flood/TrainingAlgorithm/GradientDescent.cpp -o GradientDescent.o

ConjugateGradient.o: ../Flood/TrainingAlgorithm/ConjugateGradient.cpp
	g++ -c ../Flood/TrainingAlgorithm/ConjugateGradient.cpp -o ConjugateGradient.o

QuasiNewtonMethod.o: ../Flood/TrainingAlgorithm/QuasiNewtonMethod.cpp
	g++ -c ../Flood/TrainingAlgorithm/QuasiNewtonMethod.cpp -o QuasiNewtonMethod.o

RandomSearch.o: ../Flood/TrainingAlgorithm/RandomSearch.cpp
	g++ -c ../Flood/TrainingAlgorithm/RandomSearch.cpp -o RandomSearch.o

EvolutionaryAlgorithm.o: ../Flood/TrainingAlgorithm/EvolutionaryAlgorithm.cpp
	g++ -c ../Flood/TrainingAlgorithm/EvolutionaryAlgorithm.cpp -o EvolutionaryAlgorithm.o


# Utilities

InputTargetDataSet.o: ../Flood/Utilities/InputTargetDataSet.cpp
	g++ -c ../Flood/Utilities/InputTargetDataSet.cpp -o InputTargetDataSet.o

LinearRegressionAnalysis.o: ../Flood/Utilities/LinearRegressionAnalysis.cpp
	g++ -c ../Flood/Utilities/LinearRegressionAnalysis.cpp -o LinearRegressionAnalysis.o

CorrectPredictionsAnalysis.o: ../Flood/Utilities/CorrectPredictionsAnalysis.cpp
	g++ -c ../Flood/Utilities/CorrectPredictionsAnalysis.cpp -o CorrectPredictionsAnalysis.o

IntegrationOfFunctions.o: ../Flood/Utilities/IntegrationOfFunctions.cpp
	g++ -c ../Flood/Utilities/IntegrationOfFunctions.cpp -o IntegrationOfFunctions.o

OrdinaryDifferentialEquations.o: ../Flood/Utilities/OrdinaryDifferentialEquations.cpp
	g++ -c ../Flood/Utilities/OrdinaryDifferentialEquations.cpp -o OrdinaryDifferentialEquations.o


clean:
	rm $(objects) Flood


# Flood: An Open Source Neural Networks C++ Library.
# Copyright (C) 2005-2008 Roberto Lopez 
#
# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public
# License as published by the Free Software Foundation; either
# version 2.1 of the License, or any later version.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Lesser General Public License for more details.

# You should have received a copy of the GNU Lesser General Public
# License along with this library; if not, write to the Free Software
# Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
