/****************************************************************************************************************/
/*                                                                                                              */
/*   Flood: An Open Source Neural Networks C++ Library                                                          */
/*   www.cimne.com/flood                                                                                        */
/*                                                                                                              */
/*   M E A N   S Q U A R E D   E R R O R   C L A S S                                                            */
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
#include <limits>
#include <math.h>

#include "TennixTrainer.h"
#include "socket/Evaluator.h"

namespace Flood
{

// GENERAL CONSTRUCTOR

/// General constructor. It creates a mean squared error objective functional object associated to a multilayer 
/// perceptron and measured on an input-target data set.
/// It also initializes all the rest of class members to their default values:
///
/// <ul>
/// <li> Number of evaluations = 0. 
/// <li> Epsilon: 1.0e-6.
/// <li> Display = true. 
/// </ul> 
///
/// @param newMultilayerPerceptron Pointer to a multilayer perceptron object.
/// @param newInputTargetDataSet Pointer to an input-target data set object.

TennixTrainer::TennixTrainer(MultilayerPerceptron* newMultilayerPerceptron, 
InputTargetDataSet* newInputTargetDataSet)
: ObjectiveFunctional(newMultilayerPerceptron)
{
   inputTargetDataSet = newInputTargetDataSet;
}

/*
 * Constructor for taking in a Multilayer Percepteron
 */

TennixTrainer::TennixTrainer(MultilayerPerceptron* newMultilayerPerceptron)
: ObjectiveFunctional(newMultilayerPerceptron)
{
    ;
}

// DEFAULT CONSTRUCTOR

/// Default constructor. It creates a mean squared error objective functional object not associated to any 
/// multilayer perceptron and not measured on any input-target data set.
/// It also initializes all the rest of class members to their default values:
///
/// <ul>
/// <li> Number of evaluations = 0. 
/// <li> Epsilon: 1.0e-6.
/// <li> Display = true. 
/// </ul> 

TennixTrainer::TennixTrainer(void) : ObjectiveFunctional()
{
   inputTargetDataSet = NULL;
}


// DESTRUCTOR

/// Destructor.

TennixTrainer::~TennixTrainer(void)
{

}


// METHODS

// InputTargetDataSet* getInputTargetDataSet(void) method

/// This method returns a pointer to the input-target object on which the objectiveal is measured.

InputTargetDataSet* TennixTrainer::getInputTargetDataSet(void)
{
   return(inputTargetDataSet);
}


// void setInputTargetDataSet(InputTargetDataSet*) method

/// This method sets a pointer to an input-data set object on which the objective functional is to be measured.
///
/// @param newInputTargetDataSet Pointer to an input-target data set object.

void TennixTrainer::setInputTargetDataSet(InputTargetDataSet* newInputTargetDataSet)
{
   inputTargetDataSet = newInputTargetDataSet;
}


// double calculateEvaluation(void) method

/// This method returns the evaluation value of a multilayer perceptron according to the mean squared error on 
/// an input-target data set.
///
/// @see calculateGradient(void).

double TennixTrainer::calculateEvaluation(void)
{
    // Control sentence 
    //

   if(multilayerPerceptron == NULL)
   {
      std::cerr << std::endl
                << "Flood Error: TennixTrainer class." << std::endl
                << "double calculateEvaluation(void) method." << std::endl
                << "Pointer to multilayer perceptron object cannot be NULL." << std::endl
                << std::endl;

        exit(1);
   }
   else if(inputTargetDataSet == NULL)
   {
      std::cerr << std::endl
                << "Flood Error: TennixTrainer class." << std::endl
                << "double calculateEvaluation(void) method." << std::endl
                << "Pointer to input-target data set object cannot be NULL." << std::endl
                << std::endl;

        exit(1);
   }

   // Calculate evaluation

   double meanSquaredError = 0.0;       

   // Multilayer perceptron

   int numberOfInputs = multilayerPerceptron->getNumberOfInputs();
   int numberOfOutputs = multilayerPerceptron->getNumberOfOutputs();
/*
 * Dont need any input set : sachins
 */
   // Input-target data set
/*
   int numberOfSamples = inputTargetDataSet->getNumberOfSamples();
   int numberOfInputVariables = inputTargetDataSet->getNumberOfInputVariables();
   int numberOfTargetVariables = inputTargetDataSet->getNumberOfTargetVariables();

   Matrix<double>& inputData = inputTargetDataSet->getInputData();
   Matrix<double>& targetData = inputTargetDataSet->getTargetData();

   if(numberOfInputs != numberOfInputVariables || numberOfOutputs != numberOfTargetVariables)
   {
      std::cout << std::endl
                << "Flood Error TennixTrainer class." << std::endl 
                << "double calculateEvaluation(void) method." << std::endl
                << "Number of inputs and outputs in multilayer perceptron must be equal to number of input and "
                << "output variables in input-target data set." << std::endl 
                << std::endl;

      exit(1);
   }
*/
   // Increment number of evaluations

   numberOfEvaluations++;

   double sumSquaredError = 0.0;

   Vector<double> input(numberOfInputs, 0.0);
   Vector<double> output(numberOfOutputs, 0.0);
   Vector<double> target(numberOfOutputs, 0.0);

   Evaluator evaluator;
   int game_over = false;
   GameState gamestate;

   enum game_state_variables {
       DARWIN_X =0,
       DARWIN_Y,
       OPPONENT_X,
       OPPONENT_Y,
       BALL_X,
       BALL_Y
   };

   enum NN_outputs {
       UP = 0,
       DOWN = 1,
       HIT = 2
   };

   while (game_over == false)
   {
       int keys[3] = {0,0,0};
       evaluator.get_game_state(&gamestate);

       input[DARWIN_X]   = gamestate.darwin_x;
       input[DARWIN_Y]   = gamestate.darwin_y;

       //std::cout << "Darwin Y = " << gamestate.darwin_y << "Darwin X = " << gamestate.darwin_x << std::endl;

       input[OPPONENT_X] = gamestate.opponent_x;
       input[OPPONENT_Y] = gamestate.opponent_y;

       input[BALL_X]     = gamestate.ball_x;
       input[BALL_Y]     = gamestate.ball_y;

       if (gamestate.game_ended == true)
       {
           game_over = true;
           continue;
       }


      //input = inputData.getRow(sample);


      // Output vector

      output = multilayerPerceptron->calculateOutput(input);

      //std::cout << "output[UP] = " << output[UP] << " output[DOWN] = " << output[DOWN] << "output[HIT] = " << output[HIT] << std::endl;


      keys[UP] = keys[DOWN] = keys[HIT] = 0;
      
      if (output[UP] > 0) {
          //std::cout << "UP" << std::endl;
          keys[UP] = 1;
      }
      if (output[DOWN] > 0 ) {
          //std::cout << "DOWN" << std::endl;
          keys[DOWN] = 1;
      }
      if (output[HIT] > 0 ) {
          //std::cout << "HIT" << std::endl;
          keys[HIT] = 1;
      }

      if (keys[UP] && keys[DOWN])
      {
          /* Make the stronger one be sent */
          if (output[UP] >= output[DOWN])
          {
              keys[UP] = 1;
              keys[DOWN] = 0;
          } else {
              keys[UP] = 0;
              keys[DOWN] = 1;
          }
      }


      evaluator.send_NN_response(keys, 3);





      // Target vector

     //target = targetData.getRow(sample);

      // Sum of squares error

      //sumSquaredError += (output-target).dot(output-target); 
   }



#define BEST_AI_FITNESS 70000 // changes with fitness weightages, so recalculate everytime any change

#define PENALISE_NON_MOVER 1 // Comment this line , if you dont want to penalise non movers

#ifdef PENALISE_NON_MOVER

   if (gamestate.fitness < 500 && gamestate.fitness > -500)  /* Found range for inactive players */
   {
       std::cout << "Fitness = " << BEST_AI_FITNESS + 10000 << std::endl;
       return (BEST_AI_FITNESS + 10000); /* When player does not move a bit, penalise him */  
   }
#endif

   /* Best AI players fitness is BEST_AI_FITNESS (Note using current calculation method, needs update if it changes ) */

//Lower the fitness value the better
   std::cout << "Fitness = " << (BEST_AI_FITNESS - gamestate.fitness) << std::endl;

   return (BEST_AI_FITNESS - gamestate.fitness);
}




// void saveInputTargetAndOutput(char*) method

/// This method saves to a file the inputs and the targets in the input-target data set, together with the 
/// outputs from the multilayer perceptron for that inputs:
///
/// <ul>
/// <li> Inputs.
/// <li> Targets.
/// <li> Outputs.
/// </ul> 
///
/// @param filename Filename.

void TennixTrainer::saveInputTargetAndOutput(char* filename)
{
   int numberOfOutputs = multilayerPerceptron->getNumberOfOutputs();

   int numberOfSamples = inputTargetDataSet->getNumberOfSamples();
   int numberOfInputVariables = inputTargetDataSet->getNumberOfInputVariables();
   int numberOfTargetVariables = inputTargetDataSet->getNumberOfTargetVariables();

   Matrix<double> inputData = inputTargetDataSet->getInputData();
   Matrix<double> targetData = inputTargetDataSet->getTargetData();

   Matrix<double> outputData(numberOfSamples, numberOfOutputs, 0.0);

   Vector<double> input(numberOfInputVariables, 0.0);
   Vector<double> output(numberOfTargetVariables, 0.0);

   for(int sample = 0; sample < numberOfSamples; sample++)
   {
      input = inputData.getRow(sample);

      output = multilayerPerceptron->calculateOutput(input);

      outputData.setRow(sample, output);
   }

   std::fstream file;

   file.open(filename, std::ios::out);

   // Control sentence

   if(!file.is_open())
   {
      std::cout << std::endl
                << "Flood Error TennixTrainer class." << std::endl
                << "void saveInputTargetAndOutput(char*) method." << std::endl
                << "Cannot open input-target-output data file." << std::endl
                << std::endl;

      exit(1);
   }
   else
   {
      if(display)
     {
         std::cout << std::endl
                   << "Saving input-target-output to data file..."
                   << std::endl;
     }
   }

   // Write file header

   file << "% Flood Neural Network. Input-target-output data file." << std::endl
        << "% 1 - Input data" << std::endl
        << "% 2 - Target data" << std::endl
        << "% 3 - Output data" << std::endl
        << std::endl;

   // Write file data

   for(int sample = 0; sample < numberOfSamples; sample++)
   {
      // Write sample input data

      for(int i = 0; i < numberOfInputVariables; i++)
      {
         file << inputData[sample][i] << " ";
      }

      // Write sample target data

      for(int i = 0; i < numberOfTargetVariables; i++)
      {
         file << targetData[sample][i] << " ";
      }

      // Write sample output data

      for(int i = 0; i < numberOfOutputs; i++)
      {
         file << outputData[sample][i] << " ";
      }

      file << std::endl;
   }

   file << std::endl;

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
