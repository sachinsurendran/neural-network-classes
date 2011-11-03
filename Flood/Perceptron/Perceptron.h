/****************************************************************************************************************/
/*                                                                                                              */
/*   Flood: An Open Source Neural Networks C++ Library                                                          */
/*   www.cimne.com/flood                                                                                        */
/*                                                                                                              */
/*   P E R C E P T R O N   C L A S S   H E A D E R                                                              */
/*                                                                                                              */
/*   Roberto Lopez                                                                                              */
/*   International Center for Numerical Methods in Engineering (CIMNE)                                          */
/*   Technical University of Catalonia (UPC)                                                                    */
/*   Barcelona, Spain                                                                                           */
/*   E-mail: rlopez@cimne.upc.edu                                                                               */
/*                                                                                                              */
/****************************************************************************************************************/

#ifndef __PERCEPTRON_H__
#define __PERCEPTRON_H__

#include "../Utilities/Vector.h"

namespace Flood
{

/// This class represents the concept of perceptron neuron model.

class Perceptron
{
public:

   // ENUMERATIONS

   /// Enumeration of available activation functions for the perceptron neuron model. 

   enum ActivationFunction{Logistic, HyperbolicTangent, Linear};

protected: 

   // FIELDS

   /// Activation function variable. 

   ActivationFunction activationFunction;

   /// Number of inputs in the neuron.

   int numberOfInputs;

   /// Bias value.

   double bias;

   /// Synaptic weight values.

   Vector<double> synapticWeights;

   /// Display messages to screen. 

   bool display;

public:

   // GENERAL CONSTRUCTOR

   Perceptron(int);


   // DEFAULT CONSTRUCTOR

   Perceptron(void);


   // DESTRUCTOR

   virtual ~Perceptron(void);


   // METHODS

   // Get methods

   ActivationFunction getActivationFunction(void);

   int getNumberOfInputs(void);

   double getBias(void);   

   Vector<double> getSynapticWeights(void);
   double getSingleSynapticWeight(int);

   int getNumberOfNeuronParameters(void);
   Vector<double> getNeuronParameters(void);

   bool getDisplay(void);

   // Set methods

   void setActivationFunction(ActivationFunction);

   void setNumberOfInputs(int);

   void setBias(double);
   void setSynapticWeights(Vector<double>);
   void setSingleSynapticWeight(int, double);

   void setNeuronParameters(Vector<double>);

   void setDisplay(bool);

   // Bias and synaptic weights methods

   void initBiasUniform(double, double);
   void initBiasNormal(double, double);

   void initSynapticWeightsUniform(double, double);
   void initSynapticWeightsNormal(double, double);

   // Net input methods

   double calculateNetInputSignal(Vector<double>);

   // Output signal methods

   double calculateOutputSignal(double);
   double calculateOutputSignal(Vector<double>);

   // Output signal derivative methods

   double calculateOutputSignalDerivative(double);
   double calculateOutputSignalDerivative(Vector<double>);

   // Output signal second derivative methods
  
   double calculateOutputSignalSecondDerivative(double);
   double calculateOutputSignalSecondDerivative(Vector<double>);

   // Utility methods

   void print(void);
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

