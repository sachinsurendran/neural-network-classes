/****************************************************************************************************************/
/*                                                                                                              */
/*   Flood: An Open Source Neural Networks C++ Library                                                          */
/*   www.cimne.com/flood                                                                                        */
/*                                                                                                              */
/*   O B J E C T I V E   F U N C T I O N A L   C L A S S                                                        */
/*                                                                                                              */
/*   Roberto Lopez                                                                                              */
/*   International Center for Numerical Methods in Engineering (CIMNE)                                          */
/*   Technical University of Catalonia (UPC)                                                                    */
/*   Barcelona, Spain                                                                                           */
/*   E-mail: rlopez@cimne.upc.edu                                                                               */
/*                                                                                                              */
/****************************************************************************************************************/


#include "ObjectiveFunctional.h"

#include<iostream>
#include <math.h>

namespace Flood
{

// GENERAL CONSTRUCTOR

/// General constructor. It creates an objective functional object associated to a multilayer perceptron object. 
/// It also initializes all the rest of class members to their default values:
///
/// <ul>
/// <li> Numerical epsilon method: Relative.
/// <li> Numerical differentiation method: Central differences.
/// <li> Number of evaluations = 0.
/// <li> Numerical epsilon: 1.0e-6.
/// <li> Display: true.
/// </ul> 
///
/// @param newMultilayerPerceptron Pointer to a multilayer perceptron object.

ObjectiveFunctional::ObjectiveFunctional(MultilayerPerceptron* newMultilayerPerceptron)
{
   multilayerPerceptron = newMultilayerPerceptron;

   numericalEpsilonMethod = Relative;

   numericalDifferentiationMethod = CentralDifferences;

   numberOfEvaluations = 0;

   numericalEpsilon = 1.0e-6;
      
   display = true;

}


// DEFAULT CONSTRUCTOR

/// General constructor. It creates an objective functional object not associated to any multilayer perceptron 
/// object.
/// It also initializes all the rest of class members to their default values:
///
/// <ul>
/// <li> Numerical epsilon method: Relative.
/// <li> Numerical differentiation method: Central differences.
/// <li> Number of evaluations = 0.
/// <li> Numerical epsilon: 1.0e-6.
/// <li> Display: true.
/// </ul> 

ObjectiveFunctional::ObjectiveFunctional(void)
{
   multilayerPerceptron = NULL;

   numericalEpsilonMethod = Relative;

   numericalDifferentiationMethod = CentralDifferences;

   numberOfEvaluations = 0;

   numericalEpsilon = 1.0e-6;
      
   display = true;
}


// DESTRUCTOR

/// Destructor.

ObjectiveFunctional::~ObjectiveFunctional(void)
{

}


// METHODS


// MultilayerPerceptron* getMultilayerPerceptron(void) method

/// This method returns a pointer to the multilayer perceptron object associated to the objective functional.

MultilayerPerceptron* ObjectiveFunctional::getMultilayerPerceptron(void)
{
   return(multilayerPerceptron);
}


// double getNumericalEpsilon(void) method

/// This method returns the epsilon value for the calculation of the objective function gradient by means of 
/// numerical differentiation.
///
/// @see calculateGradient(void)
/// @see calculateHessian(void)
/// @see calculateVectorHessianProduct(Vector<double>)

double ObjectiveFunctional::getNumericalEpsilon(void)
{
   return(numericalEpsilon);
}


// int getNumberOfEvaluations(void) method

/// This method returns the number of calls to the calculateEvaluation(void) method.
///
/// @see calculateEvaluation(void).

int ObjectiveFunctional::getNumberOfEvaluations(void)
{
   return(numberOfEvaluations);
}


// bool getDisplay(void) method

/// This method returns true if messages from this class can be displayed on the screen, or false if messages 
/// from this class can't be displayed on the screen.

bool ObjectiveFunctional::getDisplay(void)
{
   return(display);
}


// NumericalEpsilonMethod getNumericalEpsilonMethod(void) method

/// This method returns the method used for numerical epsilon.

ObjectiveFunctional::NumericalEpsilonMethod 
ObjectiveFunctional::getNumericalEpsilonMethod(void)
{
   return(numericalEpsilonMethod);                           
}


// NumericalDifferentiationMethod getNumericalDifferentiationMethod(void) method

/// This method returns the method used for numerical differentiation.

ObjectiveFunctional::NumericalDifferentiationMethod 
ObjectiveFunctional::getNumericalDifferentiationMethod(void)
{
   return(numericalDifferentiationMethod);                           
}


// void setMultilayerPerceptron(MultilayerPerceptron*) method

/// This method sets a pointer to a multilayer perceptron object which is to be associated to the objective 
/// functional.
///
/// @param newMultilayerPerceptron Pointer to a multilayer percepron object to be associated to the objective 
/// functional.

void ObjectiveFunctional::setMultilayerPerceptron(MultilayerPerceptron* newMultilayerPerceptron)
{
   multilayerPerceptron = newMultilayerPerceptron;
}


// void setNumericalEpsilon(double) method

/// This method sets a new epsilon value for the calculation of the objective function gradient by means of 
/// numerical differentiation.
///
/// @param newNumericalEpsilon New value for epsilon.
///
/// @see calculateGradient(void)

void ObjectiveFunctional::setNumericalEpsilon(double newNumericalEpsilon)
{
   // Control sentence

   if(newNumericalEpsilon <= 0.0)
   {
      std::cerr << std::endl
                << "Flood Error: ObjectiveFunction class." << std::endl
                << "void setNumericalEpsilon(double) method." << std::endl
                << "Numerical epsilon must be greater than 0." << std::endl
                << std::endl;

      exit(1);
   }

   // Set numerical epsilon

   numericalEpsilon = newNumericalEpsilon;

}


// void setNumberOfEvaluations(int) method

/// This method sets the number of calls to the getObjective(void) method to a new value. 
///
/// @param newNumberOfEvaluations New number of calls to the calculateEvaluation(void) method.
///
/// @see calculateEvaluation(void).

void ObjectiveFunctional::setNumberOfEvaluations(int newNumberOfEvaluations)
{
   // Control sentence

   if(newNumberOfEvaluations < 0)
   {
      std::cerr << std::endl
                << "Flood Error: ObjectiveFunction class." << std::endl
                << "void setNumberOfEvaluations(double) method." << std::endl
                << "Number of evaluations must be equal or greater than 0." << std::endl
                << std::endl;

      exit(1);
   }

   // Set number of evaluations

   numberOfEvaluations = newNumberOfEvaluations;
}


// void setDisplay(bool) method

/// This method sets a new display value. 
/// If it is set to true messages from this class are to be displayed on the screen;
/// if it is set to false messages from this class are not to be displayed on the screen.
///
/// @param newDisplay Display value.

void ObjectiveFunctional::setDisplay(bool newDisplay)
{
   display = newDisplay;
}


// void setNumericalEpsilonMethod(NumericalEpsilonMethod)

/// This method sets the method to be used for obtaining a numerical epsilon value to be used in numerical 
/// differentiation. 
///
/// @param newNumericalEpsilonMethod New numerical epsilon method.
///
/// @see calculateGradient(void).
/// @see calculateHessian(void).

void ObjectiveFunctional
::setNumericalEpsilonMethod(ObjectiveFunctional::NumericalEpsilonMethod newNumericalEpsilonMethod)
{
   numericalEpsilonMethod = newNumericalEpsilonMethod;
}


// void setNumericalDifferentiationMethod(NumericalDifferentiationMethod)

/// This method sets the method to be used for the numerical differentiation of the Jacobian matrix for the 
/// multilayer perceptron.
///
/// @param newNumericalDifferentiationMethod New numerical differentiation method.
///
/// @see calculateGradient(void).
/// @see calculateHessian(void).

void ObjectiveFunctional::setNumericalDifferentiationMethod
(ObjectiveFunctional::NumericalDifferentiationMethod newNumericalDifferentiationMethod)
{
   numericalDifferentiationMethod = newNumericalDifferentiationMethod;
}


// double calculatePotentialEvaluation(Vector<double>) method

/// This method returns which would be the objective evaluation of a multilayer perceptron for an hypothetical 
/// vector of free parameters. It does not set that vector of free parameters to the multilayer perceptron. 
///
/// @param potentialFreeParameters Vector of a potential free parameters for the multilayer perceptron associated 
/// to the objective functional.
///
/// @see calculateEvaluation(void).

double ObjectiveFunctional::calculatePotentialEvaluation(Vector<double> potentialFreeParameters)
{
   // Control sentence

   int size = potentialFreeParameters.getSize();

   int numberOfFreeParameters = multilayerPerceptron->getNumberOfFreeParameters();

   if(size != numberOfFreeParameters)
   {
      std::cerr << std::endl
                << "Flood Error: ObjectiveFunction class." << std::endl
                << "double calculatePotentialEvaluation(Vector<double>) method." << std::endl
                << "Size must be equal to number of free parameters." << std::endl
                << std::endl;

      exit(1);
   }

   // Calculate potential evaluation

   double potentialEvaluation = 0.0;

   Vector<double> originalFreeParameters = multilayerPerceptron->getFreeParameters();

   // Set potential parameters

   multilayerPerceptron->setFreeParameters(potentialFreeParameters);

   // Get Objective

   potentialEvaluation = calculateEvaluation();

   // Restart original free parameters

   multilayerPerceptron->setFreeParameters(originalFreeParameters);


   return(potentialEvaluation);
}


// double calculateGradient(void) method

/// This method returns the default objective function gradient vector, which is computed using numerical 
/// differentiation.
///
/// @see NumericalEpsilonMethod.
/// @see NumericalDifferentiationMethod.

Vector<double> ObjectiveFunctional::calculateGradient(void)
{
   int numberOfFreeParameters = multilayerPerceptron->getNumberOfFreeParameters();

   Vector<double> gradient(numberOfFreeParameters);

   Vector<double> potentialFreeParameters = multilayerPerceptron->getFreeParameters();   

   double actualEpsilon = 0.0;

   // Numerical differentiation method

   switch(numericalDifferentiationMethod)   
   {
      case ForwardDifferences:
      {
         double evaluation = calculateEvaluation();

         double evaluationForward = 0.0;

         for(int i = 0; i < numberOfFreeParameters; i++)
         {
            // Numerical epsilon method

            switch(numericalEpsilonMethod)   
            {
               case Absolute:
               {
                  actualEpsilon = numericalEpsilon;
               }             
               break;
            
               case Relative:
               {
                  if(fabs(potentialFreeParameters[i]) < 1.0)
                  {
                     actualEpsilon = numericalEpsilon;              
                  }
                  else
                  {     
                     actualEpsilon = numericalEpsilon*fabs(potentialFreeParameters[i]);
                  }
               }                          
               break;           

               default:
               {
                  std::cerr << std::endl 
                            << "Flood Error: ObjectiveFunctional class." << std::endl
                            << "Vector<double> calculateGradient(void) method." << std::endl
                            << "Unknown numerical epsilon method." << std::endl
                            << std::endl;
 
                  exit(1);
               }// end default
               break;
 
            }// end numerical epsilon switch

            // Add epsilon to free parameter

            potentialFreeParameters[i] += actualEpsilon;

            // Get potential evaluation

            evaluationForward = calculatePotentialEvaluation(potentialFreeParameters);

            // Restart free parameter

            potentialFreeParameters[i] -= actualEpsilon;

            // Calculate derivative

            gradient[i] = (evaluationForward - evaluation)/actualEpsilon;
         }// end for 
      }

      break;

      case CentralDifferences:
      {
         double evaluationForward = 0.0;
         double evaluationBackward = 0.0;

         for(int i = 0; i < numberOfFreeParameters; i++)
         {
            // Numerical epsilon method

            switch(numericalEpsilonMethod)   
            {
               case Absolute:
               {
                  actualEpsilon = numericalEpsilon;
               }             
               break;
            
               case Relative:
               {
                  if(fabs(potentialFreeParameters[i]) < 1.0)
                  {
                     actualEpsilon = numericalEpsilon;              
                  }
                  else
                  {     
                     actualEpsilon = numericalEpsilon*fabs(potentialFreeParameters[i]);
                  }
               }                            
               break;           

               default:
               {
                  std::cerr << std::endl 
                            << "Flood Error: ObjectiveFunctional class." << std::endl
                            << "Vector<double> calculateGradient(void) method." << std::endl
                            << "Unknown numerical epsilon method." << std::endl
                            << std::endl;
 
                  exit(1);
               }// end default
               break;
 
            }// end numerical epsilon method switch

            // Add epsilon to free parameter

            potentialFreeParameters[i] += actualEpsilon;

            // Get potential evaluation

            evaluationForward = calculatePotentialEvaluation(potentialFreeParameters);

            // Restart free parameter

            potentialFreeParameters[i] -= actualEpsilon;

            // Substract epsilon from free parameter

            potentialFreeParameters[i] -= actualEpsilon;

            // Get potential evaluation

            evaluationBackward = calculatePotentialEvaluation(potentialFreeParameters);

            // Restart free parameter

            potentialFreeParameters[i] += actualEpsilon;

            // Calculate derivative

            gradient[i] = (evaluationForward - evaluationBackward)/(2.0*actualEpsilon);
         }// end for
      }

      break;
   }// end numerical differentiation method switch


   return(gradient);
}


// Vector<double> calculatePotentialGradient(Vector<double>) method

/// This method returns which would be the objective function gradient of a multilayer perceptron for an 
/// hypothetical vector of free parameters.
/// It does not set that vector of free parameters to the multilayer perceptron.
///
///
/// @param potentialFreeParameters Vector of a potential free parameters for the multilayer perceptron associated 
/// to the objective functional.
///
/// @see calculateGradient(void).

Vector<double> ObjectiveFunctional::calculatePotentialGradient(Vector<double> potentialFreeParameters)
{
   // Control sentence

   int size = potentialFreeParameters.getSize();

   int numberOfFreeParameters = multilayerPerceptron->getNumberOfFreeParameters();

   if(size != numberOfFreeParameters)
   {
      std::cerr << std::endl
                << "Flood Error: ObjectiveFunction class." << std::endl
                << "double calculatePotentialGradient(Vector<double>) method." << std::endl
                << "Size must be equal to number of free parameters." << std::endl
                << std::endl;

      exit(1);
   }

   // Calculate potential gradient

   Vector<double> originalFreeParameters = multilayerPerceptron->getFreeParameters();

   Vector<double> potentialGradient(numberOfFreeParameters);

   // Set potential parameters

   multilayerPerceptron->setFreeParameters(potentialFreeParameters);

   // Get Objective function gradient

   potentialGradient = calculateGradient();

   // Restart original free parameters

   multilayerPerceptron->setFreeParameters(originalFreeParameters);

   return(potentialGradient);
}


// Matrix<double> calculateHessian(void) method

/// This method returns the default objective function Hessian matrix, which is computed using numerical 
/// differentiation.
///
/// @see NumericalEpsilonMethod.
/// @see NumericalDifferentiationMethod.

Matrix<double> ObjectiveFunctional::calculateHessian(void)
{
   int numberOfFreeParameters = multilayerPerceptron->getNumberOfFreeParameters();

   Vector<double> potentialFreeParameters = multilayerPerceptron->getFreeParameters();

   Matrix<double> hessian(numberOfFreeParameters, numberOfFreeParameters);

   double actualEpsilonI = 0.0;
   double actualEpsilonJ = 0.0;

   switch(numericalDifferentiationMethod)   
   {
      case ForwardDifferences:
      {                      
         std::cerr << std::endl
                   << "Flood Error: ObjectiveFunction class." << std::endl
                   << "Matrix<double> calculateHessian(void) method." << std::endl
                   << "Numerical differentiation of Hessian matrix with forward differences is not yet implemented." 
                   << std::endl
                   << "Please use central differences." << std::endl
                   << std::endl;

         exit(1);
      }

      break;
                                     
      case CentralDifferences:
      {
         double evaluationBackwardBackward = 0.0; 
         double evaluationBackwardForward = 0.0;
         double evaluationForwardBackward = 0.0; 
         double evaluationForwardForward = 0.0;

         // Obtain the upper part of the Hessian matrix

         for(int i = 0; i < numberOfFreeParameters; i++)
         {
            for(int j = i; j < numberOfFreeParameters; j++)
            {
               // Numerical epsilon method
  
               switch(numericalEpsilonMethod)   
               {
                  case Absolute:
                  {
                     actualEpsilonI = numericalEpsilon;
                     actualEpsilonJ = numericalEpsilon;
                  } 
            
                  break;
            
                  case Relative:
                  {
                     // Obtain epsilon for free parameter i
                      
                     if(fabs(potentialFreeParameters[i]) < 1.0)
                     {
                        actualEpsilonI = numericalEpsilon;                   
                     }
                     else
                     {
                        actualEpsilonI = numericalEpsilon*fabs(potentialFreeParameters[i]);                
                     }

                     // Obtain epsilon for free parameter j

                     if(potentialFreeParameters[j] < 1.0)
                     {
                        actualEpsilonJ = numericalEpsilon;                   
                     }
                     else
                     {
                        actualEpsilonJ = numericalEpsilon*fabs(potentialFreeParameters[j]);                
                     }
                  }                
            
                  break;            
               }
             
               // Perturb potential free parameters i and j

               potentialFreeParameters[i] += actualEpsilonI;
               potentialFreeParameters[j] += actualEpsilonJ;

               // Calculate evaluation

               evaluationForwardForward = calculatePotentialEvaluation(potentialFreeParameters);

               // Restart potential free parameters i and j

               potentialFreeParameters[i] -= actualEpsilonI;
               potentialFreeParameters[j] -= actualEpsilonJ;

               // Perturb potential free parameters i and j

               potentialFreeParameters[i] += actualEpsilonI;
               potentialFreeParameters[j] -= actualEpsilonJ;

               // Calculate evaluation

               evaluationForwardBackward = calculatePotentialEvaluation(potentialFreeParameters);

               // Restart potential free parameters i and j

               potentialFreeParameters[i] -= actualEpsilonI;
               potentialFreeParameters[j] += actualEpsilonJ;

               // Perturb potential free parameters i and j

               potentialFreeParameters[i] -= actualEpsilonI;
               potentialFreeParameters[j] += actualEpsilonJ;

               // Calculate evaluation

               evaluationBackwardForward = calculatePotentialEvaluation(potentialFreeParameters);

               // Restart potential free parameters i and j

               potentialFreeParameters[i] += actualEpsilonI;
               potentialFreeParameters[j] -= actualEpsilonJ;

               // Perturb potential free parameters i and j

               potentialFreeParameters[i] -= actualEpsilonI;
               potentialFreeParameters[j] -= actualEpsilonJ;

               // Calculate evaluation

               evaluationBackwardBackward = calculatePotentialEvaluation(potentialFreeParameters);

               // Restart potential free parameters i and j

               potentialFreeParameters[i] += actualEpsilonI;
               potentialFreeParameters[j] += actualEpsilonJ;

               // Calculate second derivative

               hessian[i][j] = (evaluationForwardForward - evaluationForwardBackward 
                              - evaluationBackwardForward + evaluationBackwardBackward)
                              /(4.0*pow(numericalEpsilon,2));
            }
         }
      }

      break;
   }// end numerical differentiation method switch

   // Obtain the rest of elements by symmetry

   for(int i = 0; i < numberOfFreeParameters; i++)
   {
      for(int j = 0; j < i; j++)
      {
         hessian[i][j] = hessian[j][i];
      }
   }

   return(hessian);
}


// Vector<double> calculatePotentialHessian(Vector<double>) method
//
/// This method returns which would be the objective function Hessian of a multilayer perceptron for an 
/// hypothetical vector of free parameters.
/// It does not set that vector of free parameters to the multilayer perceptron.
///
///
/// @param potentialFreeParameters Vector of a potential free parameters for the multilayer perceptron associated 
/// to the objective functional.
///
/// @see calculateHessian(void).

Matrix<double> ObjectiveFunctional::calculatePotentialHessian(Vector<double> potentialFreeParameters)
{
   // Control sentence

   int size = potentialFreeParameters.getSize();

   int numberOfFreeParameters = multilayerPerceptron->getNumberOfFreeParameters();

   if(size != numberOfFreeParameters)
   {
      std::cerr << std::endl
                << "Flood Error: ObjectiveFunction class." << std::endl
                << "double calculatePotentialHessian(Vector<double>) method." << std::endl
                << "Size must be equal to number of free parameters." << std::endl
                << std::endl;

      exit(1);
   }

   // Calculate potential Hessian

   Vector<double> originalFreeParameters = multilayerPerceptron->getFreeParameters();

   Matrix<double> potentialHessian(numberOfFreeParameters, numberOfFreeParameters);

   // Set potential parameters

   multilayerPerceptron->setFreeParameters(potentialFreeParameters);

   // Get objective function Hessian

   potentialHessian = calculateHessian();

   // Restart original free parameters

   multilayerPerceptron->setFreeParameters(originalFreeParameters);

   return(potentialHessian);
}


// Matrix<double> calculateInverseHessian(void) method

/// This method returns inverse matrix of the Hessian.
/// It first computes the Hessian matrix and then computes its inverse. 

Matrix<double> ObjectiveFunctional::calculateInverseHessian(void)
{  
   Matrix<double> hessian = calculateHessian();
   
   Matrix<double> inverseHessian = hessian.calculateInverse();
      
   return(inverseHessian);               
}


// Matrix<double> calculateDFPInverseHessianApproximation(
// Vector<double>, Vector<double>, Matrix<double>,
// Vector<double>, Vector<double>) method

/// This method returns an approximation of the inverse Hessian matrix according to the Davidon-Fletcher-Powel 
/// (DFP) algorithm. 
///
/// @param oldFreeParameters A previous set of free parameters.
/// @param oldGradient The gradient of the objective function for that previous set of free parameters.
/// @param oldInverseHessian The Hessian of the objective function for that previous set of free parameters.
/// @param freeParameters Actual set of free parameters.
/// @param gradient The gradient of the objective function for the actual set of free parameters.

Matrix<double> ObjectiveFunctional::calculateDFPInverseHessianApproximation(
Vector<double> oldFreeParameters, Vector<double> oldGradient, Matrix<double> oldInverseHessian,
Vector<double> freeParameters, Vector<double> gradient)
{
   // Control sentence

   int oldFreeParametersSize = oldFreeParameters.getSize();
   int oldGradientSize = oldGradient.getSize();
   int freeParametersSize = freeParameters.getSize();
   int gradientSize = gradient.getSize();

   int numberOfRows = oldInverseHessian.getNumberOfRows();
   int numberOfColumns = oldInverseHessian.getNumberOfColumns();

   int numberOfFreeParameters = multilayerPerceptron->getNumberOfFreeParameters();            
    
   if(oldFreeParametersSize != numberOfFreeParameters)
   {
      std::cerr << std::endl
                << "Flood Error: ObjectiveFunction class." << std::endl
                << "calculateDFPInverseHessianApproximation method." << std::endl
                << "Size of old free parameters vector must be equal to number of free parameters." << std::endl
                << std::endl;

      exit(1);
   }
   else if(oldGradientSize != numberOfFreeParameters)
   {
      std::cerr << std::endl
                << "Flood Error: ObjectiveFunction class." << std::endl
                << "calculateDFPInverseHessianApproximation method." << std::endl
                << "Size of old gradient vector must be equal to number of free parameters." << std::endl
                << std::endl;

      exit(1);
   }
   else if(freeParametersSize != numberOfFreeParameters)
   {
      std::cerr << std::endl
                << "Flood Error: ObjectiveFunction class." << std::endl
                << "calculateDFPInverseHessianApproximation method." << std::endl
                << "Size of free parameters vector must be equal to number of free parameters." << std::endl
                << std::endl;

      exit(1);
   }
   else if(gradientSize != numberOfFreeParameters)
   {
      std::cerr << std::endl
                << "Flood Error: ObjectiveFunction class." << std::endl
                << "calculateDFPInverseHessianApproximation method." << std::endl
                << "Size of gradient vector must be equal to number of free parameters." << std::endl
                << std::endl;

      exit(1);
   }
   else if(numberOfRows != numberOfFreeParameters)
   {
      std::cerr << std::endl
                << "Flood Error: ObjectiveFunction class." << std::endl
                << "calculateDFPInverseHessianApproximation method." << std::endl
                << "Number of rows in old inverse Hessian must be equal to number of free parameters." 
                << std::endl << std::endl;

      exit(1);
   }
   else if(numberOfColumns != numberOfFreeParameters)
   {
      std::cerr << std::endl
                << "Flood Error: ObjectiveFunction class." << std::endl
                << "calculateDFPInverseHessianApproximation method." << std::endl
                << "Number of columns in old inverse Hessian must be equal to number of free parameters." 
                << std::endl << std::endl;

      exit(1);
   }

   // Calculate inverse Hessian approximation

   Matrix<double> inverseHessianApproximation(numberOfFreeParameters, numberOfFreeParameters);
   
   // Free parameters difference Vector
   
   Vector<double> freeParametersDifference = freeParameters - oldFreeParameters;
   
   // Gradient difference Vector
   
   Vector<double> gradientDifference = gradient - oldGradient;
   
   // Inverse Hessian approximation
  
   inverseHessianApproximation = oldInverseHessian 
   + (freeParametersDifference.outer(freeParametersDifference))/(freeParametersDifference.dot(gradientDifference))
   - (oldInverseHessian*gradientDifference).outer(oldInverseHessian*gradientDifference)/
   ((gradientDifference*oldInverseHessian).dot(gradientDifference));

   return(inverseHessianApproximation);               
}


// Matrix<double> calculateBFGSInverseHessianApproximation(
// Vector<double>, Vector<double>, Matrix<double>,
// Vector<double>, Vector<double>) method

/// This method returns an approximation of the inverse Hessian matrix according to the 
/// Broyden-Fletcher-Goldfarb-Shanno (BGFS) algorithm. 
///
/// @param oldFreeParameters A previous set of free parameters.
/// @param oldGradient The gradient of the objective function for that previous set of free parameters.
/// @param oldInverseHessian The Hessian of the objective function for that previous set of free parameters.
/// @param freeParameters Actual set of free parameters.
/// @param gradient The gradient of the objective function for the actual set of free parameters.

Matrix<double> ObjectiveFunctional::calculateBFGSInverseHessianApproximation(
Vector<double> oldFreeParameters, Vector<double> oldGradient, Matrix<double> oldInverseHessian,
Vector<double> freeParameters, Vector<double> gradient)
{
   // Control sentence

   int oldFreeParametersSize = oldFreeParameters.getSize();
   int oldGradientSize = oldGradient.getSize();
   int freeParametersSize = freeParameters.getSize();
   int gradientSize = gradient.getSize();

   int numberOfRows = oldInverseHessian.getNumberOfRows();
   int numberOfColumns = oldInverseHessian.getNumberOfColumns();

   int numberOfFreeParameters = multilayerPerceptron->getNumberOfFreeParameters();            
    
   if(oldFreeParametersSize != numberOfFreeParameters)
   {
      std::cerr << std::endl
                << "Flood Error: ObjectiveFunction class." << std::endl
                << "calculateBFGSInverseHessianApproximation method." << std::endl
                << "Size of old free parameters vector must be equal to number of free parameters." << std::endl
                << std::endl;

      exit(1);
   }
   else if(oldGradientSize != numberOfFreeParameters)
   {
      std::cerr << std::endl
                << "Flood Error: ObjectiveFunction class." << std::endl
                << "calculateBFGSInverseHessianApproximation method." << std::endl
                << "Size of old gradient vector must be equal to number of free parameters." << std::endl
                << std::endl;

      exit(1);
   }
   else if(freeParametersSize != numberOfFreeParameters)
   {
      std::cerr << std::endl
                << "Flood Error: ObjectiveFunction class." << std::endl
                << "calculateBFGSInverseHessianApproximation method." << std::endl
                << "Size of free parameters vector must be equal to number of free parameters." << std::endl
                << std::endl;

      exit(1);
   }
   else if(gradientSize != numberOfFreeParameters)
   {
      std::cerr << std::endl
                << "Flood Error: ObjectiveFunction class." << std::endl
                << "calculateBFGSInverseHessianApproximation method." << std::endl
                << "Size of gradient vector must be equal to number of free parameters." << std::endl
                << std::endl;

      exit(1);
   }
   else if(numberOfRows != numberOfFreeParameters)
   {
      std::cerr << std::endl
                << "Flood Error: ObjectiveFunction class." << std::endl
                << "calculateBFGSInverseHessianApproximation method." << std::endl
                << "Number of rows in old inverse Hessian must be equal to number of free parameters." << std::endl 
                << std::endl;

      exit(1);
   }
   else if(numberOfColumns != numberOfFreeParameters)
   {
      std::cerr << std::endl
                << "Flood Error: ObjectiveFunction class." << std::endl
                << "calculateBFGSInverseHessianApproximation method." << std::endl
                << "Number of columns in old inverse Hessian must be equal to number of free parameters." << std::endl 
                << std::endl;

      exit(1);
   }

   // Calculate inverse Hessian approximation

   Matrix<double> inverseHessianApproximation(numberOfFreeParameters, numberOfFreeParameters);

   // Free parameters difference Vector
   
   Vector<double> freeParametersDifference = freeParameters - oldFreeParameters;
   
   // Gradient difference Vector
   
   Vector<double> gradientDifference = gradient - oldGradient;
     
   // BGFS Vector
   
   Vector<double> bgfs = freeParametersDifference/(freeParametersDifference.dot(gradientDifference)) 
   - (oldInverseHessian*gradientDifference)/((gradientDifference*oldInverseHessian).dot(gradientDifference));
     
   // Inverse Hessian approximation
   
   inverseHessianApproximation = oldInverseHessian
   + (freeParametersDifference.outer(freeParametersDifference))/(freeParametersDifference.dot(gradientDifference)) 
   - ((oldInverseHessian*gradientDifference).outer((gradientDifference*oldInverseHessian))/
     ((gradientDifference*oldInverseHessian).dot(gradientDifference))
   + (bgfs.outer(bgfs))*((gradientDifference*oldInverseHessian).dot(gradientDifference)));   
   
   return(inverseHessianApproximation);               
}


// Vector<double> calculateVectorHessianProduct(Vector<double>) method

/// This method returns the default product of some vector with the objective function Hessian matrix, which is 
/// computed using numerical differentiation.
///
/// @see NumericalEpsilonMethod.
/// @see NumericalDifferentiationMethod.

Vector<double> ObjectiveFunctional::calculateVectorHessianProduct(Vector<double> vector)
{
   // Control sentence

   int size = vector.getSize();

   int numberOfFreeParameters = multilayerPerceptron->getNumberOfFreeParameters();

   if(size != numberOfFreeParameters)
   {
      std::cerr << std::endl
                << "Flood Error: ObjectiveFunction class." << std::endl
                << "Vector<double> calculateVectorHessianProduct(Vector<double>) method." << std::endl
                << "Size of vector must be equal to number of free parameters." << std::endl
                << std::endl;

      exit(1);
   }

   // Calculate vector Hessian product

   Vector<double> potentialFreeParameters = multilayerPerceptron->getFreeParameters();

   Vector<double> vectorHessianProduct(numberOfFreeParameters);

   switch(numericalDifferentiationMethod)   
   {
      case ForwardDifferences:
      {
         std::cerr << std::endl
                   << "Flood Error: ObjectiveFunction class." << std::endl
                   << "Vector<double> calculateVectorHessianProduct(Vector<double>) method." << std::endl
                   << "Numerical differentiation of vector-Hessian matrix with forward differences is not yet "
                   << "implemented." << std::endl
                   << "Please use central differences." << std::endl
                   << std::endl;

         exit(1);
      }

      break;
                                     
      case CentralDifferences:
      {
         Vector<double> gradientForward(numberOfFreeParameters);
         Vector<double> gradientBackward(numberOfFreeParameters);

         double actualEpsilon = 0.0;

         for(int i = 0; i < numberOfFreeParameters; i++)
         {
            // Obtain numerical epsilon value

            switch(numericalEpsilonMethod)   
            {
               case Absolute:
               {
                  actualEpsilon = numericalEpsilon;
               } 
            
               break;
            
               case Relative:
               {
                  if(fabs(potentialFreeParameters[i]) < 1.0)
                  {
                     actualEpsilon = numericalEpsilon;              
                  }
                  else
                  {     
                     actualEpsilon = numericalEpsilon*fabs(potentialFreeParameters[i]);
                  }
               }                
   	       break;
            }

            // Add actualEpsilon*vector[i] to free parameter i

            potentialFreeParameters[i] += actualEpsilon*vector[i];

            // Calculate objective function gradient

            gradientForward = calculatePotentialGradient(potentialFreeParameters);

            // Restart free parameter i

            potentialFreeParameters[i] -= actualEpsilon*vector[i];

            // Substract numericalEpsilon*vector[i] to free parameter i

            potentialFreeParameters[i] -= actualEpsilon*vector[i];

            // Calculate objective function gradient

            gradientBackward = calculatePotentialGradient(potentialFreeParameters);

            // Restart potential free parameter i

            potentialFreeParameters[i] += actualEpsilon*vector[i];

            // Calculate the product

            vectorHessianProduct[i] = (gradientForward[i] - gradientBackward[i])/(2.0*actualEpsilon);
         }
      }

      break;
   }// end numerical differentiation method switch

   return(vectorHessianProduct);
}


// virtual void print(void) method

/// This method prints any useful information about the objective function during training. By default it prints 
/// nothing.

void ObjectiveFunctional::print(void)
{
   // Do nothing by default
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
