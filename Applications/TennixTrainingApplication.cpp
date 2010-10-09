/****************************************************************************************************************/
/*                                                                                                              */ 
/*   Flood: An Open Source Neural Networks C++ Library                                                          */
/*   www.cimne.com/flood                                                                                        */
/*                                                                                                              */
/*   E V O L U T I O N A R Y   A L G O R I T H M   A P P L I C A T I O N                                        */
/*                                                                                                              */
/*   Roberto Lopez                                                                                              */ 
/*   International Center for Numerical Methods in Engineering (CIMNE)                                          */
/*   Technical University of Catalonia (UPC)                                                                    */
/*   Barcelona, Spain                                                                                           */
/*   E-mail: rlopez@cimne.upc.edu                                                                               */ 
/*                                                                                                              */  
/****************************************************************************************************************/

/// This application is an usage example of the EvolutionaryAlgorithm class in Flood.

// System includes

#include <iostream>
#include <time.h>
#include <string.h>

// Utilities includes

#include "../Flood/Utilities/Vector.h"
#include "../Flood/Utilities/Matrix.h"

#include "../Flood/Utilities/InputTargetDataSet.h"

// Multilayer perceptron includes

#include "../Flood/MultilayerPerceptron/MultilayerPerceptron.h"

// Objective functional includes

#include "../Flood/ObjectiveFunctional/TennixTrainer.h"

// Training algorithm includes

#include "../Flood/TrainingAlgorithm/EvolutionaryAlgorithm.h"


using namespace Flood;

enum command_line_options {
    NO_ARGS,
    LOAD_FROM_FILE,
    SAVE_TO_FILE,
    SAVE_TO_FILE_ON_SIGNAL /* Not implemented */
};

/*
 *  Usage:
 *  Flood -l <Filename to load Genome from>
 *  Flood -s <Filename to save Genome to upon meeting required conditions>
 *  Flood
 *
 */

static int command_parser(int argc, char* argv[])
{
    int command_type;

#ifdef COMMAND_PARSER_DEBUG
   /* Command line Arguments dump */
   std::cout << "argc = " << argc << std::endl; 
   for(int i = 0; i < argc; i++) {
       std::cout << "argv[" << i << "] = " << argv[i] << std::endl; 
   }
#endif

   if (argc == 3)
   {
       if (strncmp(argv[1],"-l",2) == 0)
       {
           /* -l options for load from file */
#ifdef COMMAND_PARSER_DEBUG
           std::cout << "Loading file" << std::endl;
#endif /* COMMAND_PARSER_DEBUG */
           command_type = LOAD_FROM_FILE;
       }
       else if (strncmp(argv[1],"-s",2) == 0)
       {
           /* -s option for save to file */
#ifdef COMMAND_PARSER_DEBUG
           std::cout << "Save file" << std::endl;
#endif /* COMMAND_PARSER_DEBUG */
           command_type = SAVE_TO_FILE;
       }
       else {
#ifdef COMMAND_PARSER_DEBUG
           std::cout << "No args" << std::endl;
#endif /* COMMAND_PARSER_DEBUG */
           command_type = NO_ARGS;
       }

   } else {
       std::cout << "No args" << std::endl;
       command_type = NO_ARGS;
   }

   return command_type;

}

int main(int argc, char* argv[])
{

   int command_type;

   std::cout << std::endl
             << "Flood Neural Network. EvolutionaryAlgorithm Application." 
             << std::endl;	

   srand((unsigned)time(NULL));

   
   command_type = command_parser(argc, argv);


   // Input-target data set object

   //InputTargetDataSet inputTargetDataSet;

//   inputTargetDataSet.load("../Data/EvolutionaryAlgorithm/InputTargetDataSet.dat");

   //inputTargetDataSet.load("../Data/EvolutionaryAlgorithm/XOR.dat");

   // Multilayer perceptron object
   //
#define WIDTH_OF_HIDDEN_LAYER 10
   Vector<int> numbersOfHiddenNeurons (6);
   numbersOfHiddenNeurons[0] = WIDTH_OF_HIDDEN_LAYER;
   numbersOfHiddenNeurons[1] = WIDTH_OF_HIDDEN_LAYER;
   numbersOfHiddenNeurons[2] = WIDTH_OF_HIDDEN_LAYER;
   numbersOfHiddenNeurons[3] = WIDTH_OF_HIDDEN_LAYER;
   numbersOfHiddenNeurons[4] = WIDTH_OF_HIDDEN_LAYER;
   numbersOfHiddenNeurons[5] = WIDTH_OF_HIDDEN_LAYER;
//   numbersOfHiddenNeurons[6] = WIDTH_OF_HIDDEN_LAYER;

   MultilayerPerceptron multilayerPerceptron(6, numbersOfHiddenNeurons, 3);

   multilayerPerceptron.setOutputLayerActivationFunction(Perceptron::HyperbolicTangent);

   // Mean squared error object

   TennixTrainer
   tennixTrainer(&multilayerPerceptron/*, &inputTargetDataSet*/);// Now just need a NN, no inputs, coz it comes from tennix server

   // Evolutionary algorithm object
   //

   EvolutionaryAlgorithm evolutionaryAlgorithm(&tennixTrainer);

   evolutionaryAlgorithm.setPopulationSize(100);

   if (command_type == LOAD_FROM_FILE) {
       evolutionaryAlgorithm.load(argv[2]);
   } else {
       /* Initialise the Evolutionary algo instance */


       evolutionaryAlgorithm.setMaximumNumberOfGenerations(1000000);
       evolutionaryAlgorithm.setReservePopulationHistory(true);
       evolutionaryAlgorithm.setReserveMeanNormHistory(true);
       evolutionaryAlgorithm.setReserveStandardDeviationNormHistory(true);
       evolutionaryAlgorithm.setReserveBestNormHistory(true);
       evolutionaryAlgorithm.setReserveMeanEvaluationHistory(true);
       evolutionaryAlgorithm.setReserveStandardDeviationEvaluationHistory(true);
       evolutionaryAlgorithm.setReserveBestEvaluationHistory(true);

       evolutionaryAlgorithm.initPopulationNormal(0.0,1.0);
       //evolutionaryAlgorithm.initPopulationUniform(-3,3);

       // Set stopping criteria

       evolutionaryAlgorithm.setEvaluationGoal(0.1);
       evolutionaryAlgorithm.setMaximumTime(100000.0);

       // Set user stuff

       evolutionaryAlgorithm.setDisplayPeriod(1); 
       // Filename to save evolutionary object upon completion of objective
       if (command_type == SAVE_TO_FILE) 
       {
           /* If the command line arg specifies a file save, remember the filename to save to */
            evolutionaryAlgorithm.setSaveToFilename(argv[2]);
       }
   }

   // Train neural network

   evolutionaryAlgorithm.train();

   // Save all training history

   //evolutionaryAlgorithm.saveTrainingHistory("../Data/EvolutionaryAlgorithm/TrainingHistory.dat");



   std::cout << std::endl;

   return(0);
}


// Flood: An Open Source Neural Networks C++ Library.
// Copyright (C) 2005-2008 Roberto Lopez 
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
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
