/****************************************************************************************************************/
/*                                                                                                              */
/*   Flood: An Open Source Neural Networks C++ Library                                                          */
/*   www.cimne.com/flood                                                                                        */
/*                                                                                                              */
/*   G R A D I E N T   D E S C E N T   C L A S S                                                                */
/*                                                                                                              */
/*   Roberto Lopez                                                                                              */
/*   International Center for Numerical Methods in Engineering (CIMNE)                                          */
/*   Technical University of Catalonia (UPC)                                                                    */
/*   Barcelona, Spain                                                                                           */
/*   E-mail: rlopez@cimne.upc.edu                                                                               */
/*                                                                                                              */
/****************************************************************************************************************/

#include "GradientDescent.h"

#include <iostream>
#include <fstream>
#include <algorithm>
#include <functional>
#include <limits>
#include <math.h>
#include <time.h>

namespace Flood
{

// GENERAL CONSTRUCTOR 
//
/// General constructor. It creates a gradient descent training algorithm object associated to an objective 
/// functional object.
/// It also initializes the class members to their default values:
///
/// Training operators:
/// <ul>
/// <li> Train rate method = Brent;
/// </ul>
///
/// Training parameters:
/// <ul>
/// <li> First train rate: 1.0.
/// <li> Bracketing factor: 2.0.
/// <li> Tolerance in train rate: 1.0e-3.
/// </ul>
///
/// Stopping criteria:
/// <ul> 
/// <li> Evaluation goal: -1.0e99.
/// <li> Gradient norm goal: 0.0.
/// <li> Maximum training time: 1.0e6.
/// <li> Maximum number of epochs: 100. 
/// </ul> 
///  
/// User stuff: 
/// <ul>
/// <li> Warning train rate: 1.0e6.
/// <li> Error train rate: 1.0e12.
/// <li> Display: true.
/// <li> Display period: 25.
/// </ul>
///
/// @param newObjectiveFunctional Pointer to an objective functional object.
///
/// @see TrainingAlgorithm.

GradientDescent::GradientDescent(ObjectiveFunctional* newObjectiveFunctional)
: TrainingAlgorithm(newObjectiveFunctional)
{
   // Training operators

   trainRateMethod = BrentMethod;

   // Training parameters

   firstTrainRate = 1.0e-6;
   bracketingFactor = 2.0;
   trainRateTolerance = 1.0e-6;

   // Stopping criteria

   evaluationGoal = -1.0e99;
   gradientNormGoal = 0.0;
   maximumTime = 1.0e6;
   maximumNumberOfEpochs = 100;

   // User stuff

   warningTrainRate = 1.0e6;
   errorTrainRate = 1.0e12;
   displayPeriod = 25;
   
   reserveTrainingDirectionHistory = false;
   reserveTrainingDirectionNormHistory = false;
   reserveTrainingRateHistory = false;
   
}


// DEFAULT CONSTRUCTOR
//
/// Default constructor. It creates a gradient descent training algorithm object not associated to any objective 
/// functional object.
/// It also initializes the class members to their default values:
///
/// Training operators:
/// <ul>
/// <li> Train rate method = Brent;
/// </ul>
///
/// Training parameters:
/// <ul>
/// <li> First train rate: 1.0.
/// <li> Bracketing factor: 2.0.
/// <li> Tolerance in train rate: 1.0e-3.
/// </ul>
///
/// Stopping criteria:
/// <ul>
/// <li> Evaluation goal: -1.0e99.
/// <li> Gradient norm goal: 0.0.
/// <li> Maximum training time: 1.0e6.
/// <li> Maximum number of epochs: 100.
/// </ul>
///
/// User stuff:
/// <ul>
/// <li> Warning train rate: 1.0e6.
/// <li> Error train rate: 1.0e12.
/// <li> Display: true.
/// <li> Display period: 25.
/// </ul>
///
/// @see TrainingAlgorithm.

GradientDescent::GradientDescent(void) : TrainingAlgorithm()
{
   // Training operators

   trainRateMethod = BrentMethod;

   // Training parameters

   firstTrainRate = 1.0;
   bracketingFactor = 2.0;
   trainRateTolerance = 1.0e-3;

   // Stopping criteria

   evaluationGoal = -1.0e99;
   gradientNormGoal = 0.0;
   maximumTime = 1.0e6;
   maximumNumberOfEpochs = 100;

   // User stuff

   warningTrainRate = 1.0e6;
   errorTrainRate = 1.0e12;
   displayPeriod = 25;
   
   reserveTrainingDirectionHistory = false;
   reserveTrainingDirectionNormHistory = false;
   reserveTrainingRateHistory = false;  
}


// DESTRUCTOR

/// Destructor.

GradientDescent::~GradientDescent(void)
{

}


// METHODS

// TrainRateMethod getTrainRateMethod(void) method

/// This method returns the train rate method used for training.
///
/// @see calculateGoldenSectionTrainRate(double, double, Vector<double>, Vector<double>).
/// @see calculateBrentMethodTrainRate(double, double, Vector<double>, Vector<double>).
/// @see train(void).

GradientDescent::TrainRateMethod GradientDescent::getTrainRateMethod(void)
{
   return(trainRateMethod);
}


// double getFirstTrainRate(void) method

/// This method returns the initial train rate value in line minimization.
///
/// @see calculateGoldenSectionTrainRate(double, double, Vector<double>, Vector<double>).
/// @see calculateBrentMethodTrainRate(double, double, Vector<double>, Vector<double>).
/// @see train(void).

double GradientDescent::getFirstTrainRate(void)
{
   return(firstTrainRate);
}


// double getBracketingFactor(void) method

/// This method returns the increase factor when bracketing a minimum in line minimization. 
///
/// @see calculateGoldenSectionTrainRate(double, double, Vector<double>, Vector<double>).
/// @see calculateBrentMethodTrainRate(double, double, Vector<double>, Vector<double>).
/// @see train(void).

double GradientDescent::getBracketingFactor(void)
{
   return(bracketingFactor);       
}


// double getTrainRateTolerance(void) method

/// This method returns the tolerance value in line minimization.
///
/// @see calculateGoldenSectionTrainRate(double, double, Vector<double>, Vector<double>).
/// @see calculateBrentMethodTrainRate(double, double, Vector<double>, Vector<double>).
/// @see train(void).

double GradientDescent::getTrainRateTolerance(void)
{
   return(trainRateTolerance);
}


// double getWarningTrainRate(void) method

/// This method returns the train rate value at wich a warning message is written to the screen during line 
/// minimization.
///
/// @see calculateGoldenSectionTrainRate(double, double, Vector<double>, Vector<double>).
/// @see calculateBrentMethodTrainRate(double, double, Vector<double>, Vector<double>).
/// @see train(void).

double GradientDescent::getWarningTrainRate(void)
{
   return(warningTrainRate);
}


// double getErrorTrainRate(void) method

/// This method returns the train rate value at wich the line minimization algorithm is assumed to fail when 
/// bracketing a minimum.
///
/// @see calculateGoldenSectionTrainRate(double, double, Vector<double>, Vector<double>).
/// @see calculateBrentMethodTrainRate(double, double, Vector<double>, Vector<double>).
/// @see train(void).

double GradientDescent::getErrorTrainRate(void)
{
   return(errorTrainRate);
}


// bool getReserveTrainingDirectionHistory(void) method

/// This method returns true if the training direction history matrix is to be reserved, and false otherwise.

bool GradientDescent::getReserveTrainingDirectionHistory(void)
{
   return(reserveTrainingDirectionHistory);     
}


// bool getReserveTrainingDirectionNormHistory(void) method

/// This method returns true if the training direction norm history vector is to be reserved, and false otherwise.

bool GradientDescent::getReserveTrainingDirectionNormHistory(void)
{
   return(reserveTrainingDirectionHistory);     
}


// bool getReserveTrainingRateHistory(void) method

/// This method returns true if the training rate history vector is to be reserved, and false otherwise.

bool GradientDescent::getReserveTrainingRateHistory(void)
{
   return(reserveTrainingRateHistory);     
}


// Matrix<double> getTrainingDirectionHistory(void) method

/// This method returns the training direction history matrix. 

Matrix<double> GradientDescent::getTrainingDirectionHistory(void)
{
   return(trainingDirectionHistory);               
}


// Vector<double> getTrainingDirectionNormHistory(void) method

/// This method returns the training direction norm history vector. 

Vector<double> GradientDescent::getTrainingDirectionNormHistory(void)
{
   return(trainingDirectionNormHistory);               
}


// Vector<double> getTrainingRateHistory(void) method

/// This method returns the training rate history vector. 

Vector<double> GradientDescent::getTrainingRateHistory(void)
{
   return(trainingRateHistory);               
}


// void setTrainRateMethod(TrainRateMethod) method

/// This method sets a new train rate method to be used for training.
///
/// @param newTrainRateMethod Train rate method.
///
/// @see calculateGoldenSectionTrainRate(double, double, Vector<double>, Vector<double>).
/// @see calculateBrentMethodTrainRate(double, double, Vector<double>, Vector<double>).
/// @see train(void).

void GradientDescent::setTrainRateMethod(GradientDescent::TrainRateMethod newTrainRateMethod)
{
   trainRateMethod = newTrainRateMethod;
}


// void setFirstTrainRate(double) method

/// This method sets a new value to be used as an initial train rate in line  minimization.
///
/// @param newFirstTrainRate Initial train rate value.
///
/// @see calculateGoldenSectionTrainRate(double, double, Vector<double>, Vector<double>).
/// @see calculateBrentMethodTrainRate(double, double, Vector<double>, Vector<double>).
/// @see train(void).

void GradientDescent::setFirstTrainRate(double newFirstTrainRate)
{
   // Control sentence (if debug)

   #ifndef NDEBUG 
     
   if(newFirstTrainRate <= 0.0)
   {
      std::cerr << std::endl
                << "Flood Error: GradientDescent class." << std::endl
                << "void setFirstTrainRate(double) method." << std::endl
                << "First train rate must be greater than 0." << std::endl
                << std::endl;

      exit(1);
   }

   #endif

   // Set first train rate

   firstTrainRate = newFirstTrainRate;
}


// void setBracketingFactor(double) method

/// This method sets a new increase factor value to be used for line minimization when bracketing a minimum.
///
/// @param newBracketingFactor Bracketing factor value.
///
/// @see calculateGoldenSectionTrainRate(double, double, Vector<double>, Vector<double>).
/// @see calculateBrentMethodTrainRate(double, double, Vector<double>, Vector<double>).
/// @see train(void).

void GradientDescent::setBracketingFactor(double newBracketingFactor)
{ 
   // Control sentence (if debug)

   #ifndef NDEBUG 

   if(newBracketingFactor <= 0.0) 
   {
      std::cerr << std::endl
                << "Flood Error: GradientDescent class." << std::endl
                << "void setBracketingFactor(double) method." << std::endl
                << "Bracketing factor must be greater than 0." << std::endl
                << std::endl;

      exit(1);
   }

   #endif

   bracketingFactor = newBracketingFactor;
}


// void setTrainRateTolerance(double) method

/// This method sets a new tolerance value to be used in line minimization.
///
/// @param newTrainRateTolerance Tolerance value in line minimization.
///
/// @see calculateGoldenSectionTrainRate(double, double, Vector<double>, Vector<double>).
/// @see calculateBrentMethodTrainRate(double, double, Vector<double>, Vector<double>).
/// @see train(void).

void GradientDescent::setTrainRateTolerance(double newTrainRateTolerance)
{
   // Control sentence (if debug)

   #ifndef NDEBUG 
                                      
   if(newTrainRateTolerance <= 0.0)
   {
      std::cerr << std::endl
                << "Flood Error: GradientDescent class." << std::endl
                << "void setTrainRateTolerance(double) method." << std::endl
                << "Tolerance must be greater than 0." << std::endl
                << std::endl;

      exit(1);
   }

   #endif

   // Set train rate tolerance

   trainRateTolerance = newTrainRateTolerance;
}


// void setWarningTrainRate(double) method

/// This method sets a new train rate value at wich a warning message is written to the screen during line 
/// minimization.
///
/// @param newWarningTrainRate Warning train rate value.
///
/// @see calculateGoldenSectionTrainRate(double, double, Vector<double>, Vector<double>).
/// @see calculateBrentMethodTrainRate(double, double, Vector<double>, Vector<double>).
/// @see train(void).

void GradientDescent::setWarningTrainRate(double newWarningTrainRate)
{
   // Control sentence (if debug)

   #ifndef NDEBUG 

   if(newWarningTrainRate <= 0.0)
   {
      std::cerr << std::endl
                << "Flood Error: GradientDescent class." << std::endl 
                << "void setWarningTrainRate(double) method." << std::endl
                << "Warning train rate must be greater than 0." << std::endl
                << std::endl;

      exit(1);
   }

   #endif

   warningTrainRate = newWarningTrainRate;
}


// void setErrorTrainRate(double) method

/// This method sets a new train rate value at wich a the line minimization algorithm is assumed to fail when 
/// bracketing a minimum.
///
/// @param newErrorTrainRate Error train rate value.
///
/// @see calculateGoldenSectionTrainRate(double, double, Vector<double>, Vector<double>).
/// @see calculateBrentMethodTrainRate(double, double, Vector<double>, Vector<double>).
/// @see train(void).

void GradientDescent::setErrorTrainRate(double newErrorTrainRate)
{
   // Control sentence (if debug)

   #ifndef NDEBUG 

   if(newErrorTrainRate <= 0.0)
   {
      std::cerr << std::endl
                << "Flood Error: GradientDescent class." << std::endl
                << "void setErrorTrainRate(double) method." << std::endl
                << "Error train rate must be greater than 0." << std::endl
                << std::endl;

      exit(1);
   }

   #endif

   // Set error train rate

   errorTrainRate = newErrorTrainRate;
}


// void setReserveTrainingDirectionHistory(bool) method

/// This method makes the training direction history matrix to be reseved or not in memory.
///
/// @param newReserveTrainingDirectionHistory True if the training direction history matrix is to be reserved, 
/// false otherwise.

void GradientDescent::setReserveTrainingDirectionHistory(bool newReserveTrainingDirectionHistory)
{
   reserveTrainingDirectionHistory = newReserveTrainingDirectionHistory;          
}


// void setReserveTrainingDirectionNormHistory(bool) method

/// This method makes the history vector of training direction norms to be reseved or not in memory.
///
/// @param newReserveTrainingDirectionNormHistory True if the training direction norm history to be reserved, 
/// false otherwise.

void GradientDescent::setReserveTrainingDirectionNormHistory(bool newReserveTrainingDirectionNormHistory)
{
   reserveTrainingDirectionNormHistory = newReserveTrainingDirectionNormHistory;          
}


// void setReserveTrainingRateHistory(bool) method

/// This method makes the training rate history vector to be reseved or not in memory.
///
/// @param newReserveTrainingRateHistory True if the training rate history vector is to be reserved, false 
/// otherwise.

void GradientDescent::setReserveTrainingRateHistory(bool newReserveTrainingRateHistory)
{
   reserveTrainingRateHistory = newReserveTrainingRateHistory;          
}


// void setReserveAllTrainingHistory(bool) method

/// This method makes the training history of all variables to reseved or not in memory.
///
/// @param newReserveAllTrainingHistory True if the training history of all variables is to be reserved, 
/// false otherwise.

void GradientDescent::setReserveAllTrainingHistory(bool newReserveAllTrainingHistory)
{
   reserveElapsedTimeHistory = newReserveAllTrainingHistory;
   reserveFreeParametersHistory = newReserveAllTrainingHistory;
   reserveFreeParametersNormHistory = newReserveAllTrainingHistory;
   reserveEvaluationHistory = newReserveAllTrainingHistory;
   reserveGradientHistory = newReserveAllTrainingHistory;
   reserveGradientNormHistory = newReserveAllTrainingHistory;
   reserveTrainingDirectionHistory = newReserveAllTrainingHistory;
   reserveTrainingDirectionNormHistory = newReserveAllTrainingHistory;
   reserveTrainingRateHistory = newReserveAllTrainingHistory;
}


// void setTrainingDirectionHistory(Matrix<double>) method

/// This method sets a new matrix containing the training direction history over the training epochs.
/// Each row in the matrix contains the training direction vector of one single epoch. 
///
/// @param newTrainingDirectionHistory Training direction history matrix. 

void GradientDescent::setTrainingDirectionHistory(Matrix<double> newTrainingDirectionHistory)
{
   trainingDirectionHistory = newTrainingDirectionHistory;     
}


// void setTrainingDirectionNormHistory(Vector<double>) method

/// This method sets a new vector containing the training direction norm history over the training epochs.
/// Each element in the vector contains the training direction norm of one single epoch. 
///
/// @param newTrainingDirectionNormHistory Training direction norm history vector. 

void GradientDescent::setTrainingDirectionNormHistory(Vector<double> newTrainingDirectionNormHistory)
{
   trainingDirectionNormHistory = newTrainingDirectionNormHistory;     
}


// void setTrainingRateHistory(Vector<double>) method

/// This method sets a new vector containing the training rate history over the training epochs.
/// Each element in the vector contains the training rate of one single epoch. 
///
/// @param newTrainingRateHistory Training rate history vector. 

void GradientDescent::setTrainingRateHistory(Vector<double> newTrainingRateHistory)
{
   trainingRateHistory = newTrainingRateHistory;     
}


// double calculateGoldenSectionTrainRate(double, double, Vector<double>, Vector<double>)
// method

/// This method returns the train rate by searching in a given direction to locate the minimum of the objective 
/// function in that direction. It uses the golden section method.
///
/// @param initialTrainRate Initial train rate in line minimization.
/// @param evaluation Network's evaluation value.
/// @param freeParameters Network's free parameters vector.
/// @param trainDirection Train direction vector.
///
/// @see calculateBrentMethodTrainRate(double, double, Vector<double>, Vector<double>).
/// @see train(void).

double GradientDescent::calculateGoldenSectionTrainRate
(double initialTrainRate, double evaluation, Vector<double> freeParameters, 
Vector<double> trainDirection)
{
   double trainRate = 0.0;

   // Multilayer perceptron

   MultilayerPerceptron* multilayerPerceptron = objectiveFunctional->getMultilayerPerceptron();

   int numberOfFreeParameters = multilayerPerceptron->getNumberOfFreeParameters();

   Vector<double> potentialFreeParameters(numberOfFreeParameters);

   double a = 0.0;
   double evaluationA = 0.0;
   double b = 0.0;
   double evaluationB = 0.0;

   double c = 0.0;
   double evaluationC = 0.0;
   double d = 0.0;
   double evaluationD = 0.0;

   double tau = (3.0-sqrt(5.0))/2.0; // 0.382

   // Start golden section search

   // Set initial a point

   a = 0.0;

   // Calculate evaluation for a

   evaluationA = evaluation;

   // Set initial b point

   b = initialTrainRate;

   // Calculate evaluation for b

   potentialFreeParameters = freeParameters + trainDirection*b;

   evaluationB = objectiveFunctional->calculatePotentialEvaluation(potentialFreeParameters);

   // Find initial interval where minimum evaluation occurs

   while(evaluationA > evaluationB)
   {
      // Set new b

      b *= bracketingFactor;

      if(b >= errorTrainRate)
      {
         std::cerr << std::endl
                   << "Flood Error: GradientDescent class." << std::endl
                   << "double calculateGoldenSectionTrainRate"
                   << "(double, double, Vector<double>, Vector<double>) method." << std::endl
                   << "Unable to bracket a minimum." << std::endl;
                   
         exit(1);
      }
      else if(display && b >= warningTrainRate)
      {
         std::cout << std::endl
                   << "Flood Warning: Train rate is " << b
                   << std::endl;
      }

      // Calculate evaluation for new b

      potentialFreeParameters = freeParameters + trainDirection*b;

      evaluationB = objectiveFunctional->calculatePotentialEvaluation(potentialFreeParameters);
   }

   // Initialize c and d (interior points for line minimization)

   // Initialize c point

   c = a + tau*(b-a);

   // Calculate evaluation for c

   potentialFreeParameters = freeParameters + trainDirection*c;

   evaluationC = objectiveFunctional->calculatePotentialEvaluation(potentialFreeParameters);

   // Initialize d point

   d = b - tau*(b-a);

   // Calculate evaluation for d

   potentialFreeParameters = freeParameters + trainDirection*d;

   evaluationD = objectiveFunctional->calculatePotentialEvaluation(potentialFreeParameters);

   // Reduce the interval with the golden section algorithm

   while(b-a > trainRateTolerance)
   {
      Vector<double> evaluationVectorLeft(3);
      evaluationVectorLeft[0] = evaluationA;
      evaluationVectorLeft[1] = evaluationC;
      evaluationVectorLeft[2] = evaluationD;

      double minimumEvaluationLeft = evaluationVectorLeft.calculateMinimum();

      Vector<double> evaluationVectorRight(3);
      evaluationVectorRight[0] = evaluationB;
      evaluationVectorRight[1] = evaluationC;
      evaluationVectorRight[2] = evaluationD;

      double minimumEvaluationRight = evaluationVectorRight.calculateMinimum();

      if((evaluationC <= evaluationD && evaluationB >= minimumEvaluationLeft)
      || (evaluationA <=  minimumEvaluationRight))

      // There is a minimum between a and b
      {
         b=d;
         d=c;

         evaluationB = evaluationD;
         evaluationD = evaluationC;

         // Set new c point

         c = a + tau*(b-a);

         // Calculate evaluation for new c

         potentialFreeParameters = freeParameters + trainDirection*c;

         evaluationC = objectiveFunctional->calculatePotentialEvaluation(potentialFreeParameters);
      }
      else if((evaluationD <= evaluationC && evaluationA >= minimumEvaluationRight)
      || (evaluationB <= minimumEvaluationLeft))

      // There is a minimum between c and b
      {
         a = c;
         c = d;

         evaluationA = evaluationC;
         evaluationC = evaluationD;

         // Set new d point

         d = b - tau*(b-a);

         // Calculate evaluation for new d

         potentialFreeParameters = freeParameters + trainDirection*d;

         evaluationD = objectiveFunctional->calculatePotentialEvaluation(potentialFreeParameters);
      }
      else
      {
         std::cerr << std::endl
                   << "Flood Error: GradientDescent class." << std::endl
                   << "double calculateGoldenSectionTrainRate"
                   << "(double, double, Vector<double>, Vector<double>) method." << std::endl
                   << "Unable to find were the minimum is." << std::endl;

         exit(1);
      }
   }

   // Get minimum evaluation and train rate among points A, B, C and D

   double minimumEvaluation = evaluation;
   trainRate = 0.0;

   if(evaluationA < minimumEvaluation)
   {
      minimumEvaluation  = evaluationA;
      trainRate = a;
   }
   else if(evaluationB < minimumEvaluation)
   {
      minimumEvaluation = evaluationB;
      trainRate = b;
   }
   else if(evaluationC < minimumEvaluation)
   {
      minimumEvaluation = evaluationC;
      trainRate = c;
   }
   else if(evaluationD < minimumEvaluation)
   {
      minimumEvaluation = evaluationD;
      trainRate = d;
   }

   return(trainRate);
}


// double calculateBrentMethodTrainRate(double, double, Vector<double>, Vector<double>)
// method

/// This method returns the train rate by searching in a given direction to locate the minimum of the evaluation 
/// function in that direction. It uses the Brent's method.
///
/// @param initialTrainRate Initial train rate in line minimization.
/// @param evaluation Network's evaluation value.
/// @param freeParameters Network's free parameters vector.
/// @param trainDirection Train direction vector.
///
/// @see calculateGoldenSectionTrainRate(double, double, Vector<double>, Vector<double>).
/// @see train(void).

double GradientDescent::calculateBrentMethodTrainRate(double initialTrainRate, 
double evaluation, Vector<double> freeParameters, Vector<double> trainDirection)
{
   double trainRate = 0.0;

   MultilayerPerceptron* multilayerPerceptron = objectiveFunctional->getMultilayerPerceptron();

   int numberOfFreeParameters = multilayerPerceptron->getNumberOfFreeParameters();

   Vector<double> potentialFreeParameters(numberOfFreeParameters);

   double a = 0.0;
   double evaluationA = evaluation;
   double b = 0.0;
   double evaluationB = evaluation;

   double u = 0.0;
   double evaluationU = evaluation;
   double v = 0.0;
   double evaluationV = evaluation;
   double w = 0.0;
   double evaluationW = evaluation;
   double x = 0.0;
   double evaluationX = evaluation;

   // Golden number

   double tau = (3.0-sqrt(5.0))/2.0; // 0.382

   // Start Brent search

   // Set initial a point

   a = 0.0;

   // Calculate evaluation for a

   evaluationA = evaluation;

   // Set initial b point

   b = initialTrainRate;

   // Calculate evaluation for b

   potentialFreeParameters = freeParameters + trainDirection*b;

   evaluationB = objectiveFunctional->calculatePotentialEvaluation(potentialFreeParameters);

   // Find initial interval where minimum evaluation occurs

   while(evaluationA > evaluationB)
   {
      // Reset a point and corresponding evaluation 

      a = b;
      evaluationA = evaluationB;

      // Set new b

      b *= bracketingFactor;

      if(b >= errorTrainRate)
      {
         std::cerr << std::endl
                   << "Flood Error: GradientDescent class." << std::endl
                   << "double calculateGoldenSectionTrainRate"
				   << "(double, double, Vector<double>, Vector<double>) method." << std::endl
                   << "Unable to bracket a minimum." << std::endl;
                   
         exit(1);
      }
      else if(display && b >= warningTrainRate)
      {
         std::cout << std::endl
                   << "Flood Warning: Train rate is  " << b
                   << std::endl;
      }

      // Calculate evaluation for new b

      potentialFreeParameters = freeParameters + trainDirection*b;

      evaluationB = objectiveFunctional->calculatePotentialEvaluation(potentialFreeParameters);
   }

   // Get inediate point V

   v = a + tau*(b-a);

   // Calculate evaluation for V

   potentialFreeParameters = freeParameters + trainDirection*v;

   evaluationV = objectiveFunctional->calculatePotentialEvaluation(potentialFreeParameters);

   // Set initial W and X points

   w = v;
   evaluationW = evaluationV;

   x = v;
   evaluationX = evaluationV;

   // Maximum and minimum intervals ???

   bool goldenSection = false;

   // Reduce the interval

   while(b-a > trainRateTolerance)
   {
      // Quadratic interpolation

      if(w != x && w != v && x != v) // Can construct parabola
      {
         // zz vector

         Vector<double> trainRateVector(3);
         trainRateVector[0] = v;
         trainRateVector[1] = w;
         trainRateVector[2] = x;

         std::sort(trainRateVector.begin(), trainRateVector.end(), std::less<double>());

         // pp vector

         Vector<double> evaluationVector(3);

         for(int i = 0; i < 3; i++)
         {
            if(trainRateVector[i] == v)
            {
               evaluationVector[i] = evaluationV;
            }
            else if(trainRateVector[i] == w)
            {
               trainRateVector[i] = evaluationW;
            }
            else if(trainRateVector[i] == x)
            {
               trainRateVector[i] = evaluationX;
            }
            else
            {
               std::cerr << std::endl
                         << "Flood Error: GradientDescent class." << std::endl
                         << "double calculateBrentMethodTrainRate" << std::endl
                         << "(double, double, Vector<double>, Vector<double>) method." << std::endl
                         << "Unable to construct train rate and evaluation vectors right." << std::endl 
                         << std::endl;

               exit(1);
            }
         }

         // xStar is the minimum of the parabola through the three train rate points

         double numerator
         = (pow(trainRateVector[2],2) - pow(trainRateVector[1],2))*evaluationVector[0]
         + (pow(trainRateVector[1],2) - pow(trainRateVector[0],2))*evaluationVector[2]
         + (pow(trainRateVector[0],2) - pow(trainRateVector[2],2))*evaluationVector[1];

         double denominator
         = (trainRateVector[2] - trainRateVector[1])*evaluationVector[0]
         + (trainRateVector[1] - trainRateVector[0])*evaluationVector[2]
         + (trainRateVector[0] - trainRateVector[2])*evaluationVector[1];

         double xStar = 0.5*numerator/denominator;

         if(xStar < b && a < xStar) // xStar is in [a,b]
         {
            u = xStar;

            // Good, no need to perform golden section

            goldenSection = false;
         }
         else // xStar is not in [a,b]
         {
            // Bad, need to perform golden section

            goldenSection = true;
         }
      }
      else // Cannot construct parabola
      {
         // Bad, need to perform golden section

         goldenSection = true;
      }

      //
      // Golden section
      //

      if(goldenSection == true)
      {
         if(x >= (a+b)/2.0)
         {
            u = x-tau*(x-a);
         }
         else
         {
            u = x+tau*(b-x);
         }
      }

      // Calculate evaluation for U

      potentialFreeParameters = freeParameters + trainDirection*u;

      evaluationU = objectiveFunctional->calculatePotentialEvaluation(potentialFreeParameters);

      // Update points

      if(evaluationU <= evaluationX)
      {
         if(u < x)
         {
            b = x;
            evaluationB = evaluationX;
         }
         else
         {
            a = x;
            evaluationA = evaluationX;
         }

         v = w;
         evaluationV = evaluationW;

         w = x;
         evaluationW = evaluationX;

         x = u;
         evaluationX = evaluationU;
      }
      else
      {
         if(u < x)
         {
            a = u;
            evaluationA = evaluationU;
         }
         else
         {
            b = u;
            evaluationB = evaluationU;
         }

         if((evaluationU <= evaluationW) || (w == x))
         {
             v = w;
             evaluationV = evaluationW;

             w = u;
             evaluationW = evaluationU;
         }
         else if((evaluationU <= evaluationV) || (v == x) || (v == w))
         {
            v = u;
            evaluationV = evaluationU;
         }
      }
   } // while loop

   // Get minimum evaluation and train rate among points A, B, V, W and X

   double minimumEvaluation = evaluation;
   trainRate = 0.0;

   if(evaluationA < minimumEvaluation)
   {
      minimumEvaluation = evaluationA;
      trainRate = a;
   }
   else if(evaluationB < minimumEvaluation)
   {
      minimumEvaluation = evaluationB;
      trainRate = b;
   }
   else if(evaluationU < minimumEvaluation)
   {
      minimumEvaluation = evaluationU;
      trainRate = u;
   }
   else if(evaluationV < minimumEvaluation)
   {
      minimumEvaluation = evaluationV;
      trainRate = v;
   }
   else if(evaluationW < minimumEvaluation)
   {
      minimumEvaluation = evaluationW;
      trainRate = w;
   }
   else if(evaluationX < minimumEvaluation)
   {
      minimumEvaluation = evaluationX;
      trainRate = x;
   }

   return(trainRate);
}


// void train(void) method

/// This method trains a multilayer perceptron with an associated objective functional according to the gradient 
/// descent algorithm.
/// Training occurs according to the training parameters. 
///
/// @see TrainRateMethod.

void GradientDescent::train(void)
{
   // Control sentence (if debug)

   #ifndef NDEBUG 

   if(objectiveFunctional == NULL)
   {
      std::cerr << std::endl
                << "Flood Error: GradientDescent class." << std::endl
                << "void train(void) method." << std::endl
                << "Pointer to objective functional object cannot be NULL." << std::endl
                << std::endl;

      exit(1);
   }

   #endif

   // Multilayer perceptron stuff

   MultilayerPerceptron* multilayerPerceptron = objectiveFunctional->getMultilayerPerceptron();

   int numberOfFreeParameters = multilayerPerceptron->getNumberOfFreeParameters();

   resizeTrainingHistory(maximumNumberOfEpochs+1);

   // Free parameters

   Vector<double> freeParameters = multilayerPerceptron->getFreeParameters();

   if(reserveFreeParametersHistory)
   {
      freeParametersHistory.setRow(0, freeParameters);                                
   }

   // Free parameters norm

   double freeParametersNorm = freeParameters.calculateNorm();

   if(reserveFreeParametersNormHistory)
   {
      freeParametersNormHistory[0] = freeParametersNorm; 
   }

   if(display && (freeParametersNorm >= warningFreeParametersNorm))
   {
      std::cout << "Flood Warning: Initial free parameters norm is " << freeParametersNorm << "." << std::endl;          
   }
   
   // Training time stuff

   time_t beginningTime, currentTime;
   double elapsedTime = 0.0;

   // Start training 

   time(&beginningTime);

   if(display)
   {
      std::cout << std::endl
                << "Training with gradient descent..." 
                << std::endl;
   }
   
   // Initial evaluation
   
   double evaluation = objectiveFunctional->calculateEvaluation();

   if(reserveEvaluationHistory)
   {
      evaluationHistory[0] = evaluation;                                
   }

   // Initial objective function gradient

   Vector<double> gradient = objectiveFunctional->calculateGradient();

   if(reserveGradientHistory)
   {
      gradientHistory.setRow(0, gradient);                                
   }

   // Initial gradient norm

   double gradientNorm = gradient.calculateNorm();

   if(display && (gradientNorm >= warningGradientNorm))
   {
      std::cout << "Flood Warning: Initial gradient norm is " << gradientNorm << "." << std::endl;          
   }

   // Free parameters history 

   if(reserveFreeParametersHistory)
   {
      freeParametersHistory.setRow(0, freeParameters);   
   }

   // Evaluation history

   if(reserveEvaluationHistory)
   {
      evaluationHistory[0] = evaluation;
   }

   // Gradient norm history
   
   if(reserveGradientNormHistory)
   {
      gradientNormHistory[0] = gradientNorm;
   }

   if(display && (gradientNorm >= warningGradientNorm))
   {
      std::cout << "Flood Warning: Initial gradient norm is " << gradientNorm << "." << std::endl;          
   }

   // Elapsed time

   time(&currentTime);
   elapsedTime = difftime(currentTime, beginningTime);

   if(reserveElapsedTimeHistory)
   {
      elapsedTimeHistory[0] = elapsedTime;
   }

   // Stopping criteria

   if(evaluation <= evaluationGoal)
   {
      if(display)
      {          
         std::cout << std::endl
                   << "Initial evaluation is less than goal." << std::endl;

         std::cout << "Elapsed time: " << elapsedTime << ";" << std::endl;
         std::cout << "Initial parameters norm: " << freeParametersNorm << std::endl;
         std::cout << "Initial evaluation: " << evaluation << std::endl;          
         std::cout << "Initial gradient norm: " << gradientNorm << std::endl;

         objectiveFunctional->print();
      }
 
      resizeTrainingHistory(1);
      
      return;
   }
   else if(gradientNorm <= gradientNormGoal)
   {  
      if(display)
      {        
         std::cout << std::endl
                   << "Initial gradient norm is less than goal." << std::endl;

         std::cout << "Elapsed time: " << elapsedTime << ";" << std::endl;
         std::cout << "Initial parameters norm: " << freeParametersNorm << std::endl;
         std::cout << "Initial evaluation: " << evaluation << std::endl;
         std::cout << "Initial gradient norm: " << gradientNorm << std::endl;

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
         std::cout << "Initial parameters norm: " << freeParametersNorm << std::endl;
         std::cout << "Initial evaluation: " <<  evaluation << std::endl;
         std::cout << "Initial gradient norm: " <<  gradientNorm << std::endl;      

         objectiveFunctional->print();
      }
   }

   // Initialize tain direction vector

   Vector<double> trainDirection(numberOfFreeParameters);

   double trainDirectionNorm = 0.0;

   // Initialize train rate stuff

   double initialTrainRate = 0.0;
   double trainRate = 0.0;
   double oldTrainRate = 0.0;

   // Main loop

   for(int epoch = 1; epoch <= maximumNumberOfEpochs; epoch++)
   {
      // Get train direction

      trainDirection = gradient*(-1.0);

      if(reserveTrainingDirectionHistory)
      {
         trainingDirectionHistory.setRow(epoch-1, trainDirection);                     
      }
      
      // Train direction norm
      
      trainDirectionNorm = trainDirection.calculateNorm();
      
      if(reserveTrainingDirectionNormHistory)
      {
         trainingDirectionNormHistory[epoch-1] = trainDirectionNorm; 
      }
      
      // Get initial train rate

      if(epoch == 1)
      {
         initialTrainRate = firstTrainRate;
      }
      else
      {
         initialTrainRate = oldTrainRate;
      }
    
      // Get train rate

      switch(trainRateMethod)
      {
         case Fixed:

            trainRate = firstTrainRate;

            break;

         case GoldenSection:

            trainRate = calculateGoldenSectionTrainRate
            (initialTrainRate, evaluation, freeParameters, trainDirection);

            break;

         case BrentMethod:

            trainRate = calculateBrentMethodTrainRate
            (initialTrainRate, evaluation, freeParameters, trainDirection);

         break;
      }

      if(reserveTrainingRateHistory)
      {
         trainingRateHistory[epoch-1] = trainRate;
      }

      // Train rate stopping criteria

      if(trainRate == 0.0)
      {
         if(display)
         {
            std::cout << std::endl
                      << "Epoch " << epoch << ": "
                      << "Train rate is zero." << std::endl;

            std::cout << "Elapsed time: " << elapsedTime << ";" << std::endl;
            std::cout << "Final parameters norm: " << freeParametersNorm << std::endl;
            std::cout << "Final evaluation: " << evaluation << std::endl;
            std::cout << "Final gradient norm: " << gradientNorm << std::endl;

            objectiveFunctional->print();
         }

         resizeTrainingHistory(epoch);

         break;
      }

      // Get new free parameters

      freeParameters = freeParameters + trainDirection*trainRate;

      // Set new free parameters

      multilayerPerceptron->setFreeParameters(freeParameters);

      // Free parameters history 

     if(reserveFreeParametersHistory)
     {
      freeParametersHistory.setRow(epoch, freeParameters);
     }

      freeParametersNorm = freeParameters.calculateNorm();

      if(display && (freeParametersNorm >= warningFreeParametersNorm))
      {
         std::cout << "Flood Warning: Free parameters norm is " << freeParametersNorm << "." << std::endl;          
      }
      
     if(reserveFreeParametersNormHistory)
     {
      freeParametersNormHistory[epoch] = freeParametersNorm;
     }
      
      // Objective function evaluation

     double oldEvaluation = evaluation;
   
      evaluation = objectiveFunctional->calculateEvaluation();

     if(reserveEvaluationHistory)
     {
      evaluationHistory[epoch] = evaluation;
     }

      // Objective function gradient

      gradient = objectiveFunctional->calculateGradient();

     if(reserveGradientHistory)
     {
      gradientHistory.setRow(epoch, gradient);
     }

      // Gradient norm

      gradientNorm = gradient.calculateNorm();

      if(display && (gradientNorm >= warningGradientNorm))
      {
         std::cout << "Flood Warning: Gradient norm is " << gradientNorm << "." << std::endl;          
      }

     if(reserveGradientNormHistory)
     {
      gradientNormHistory[epoch] = gradientNorm;
     }

      // Elapsed time

      time(&currentTime);
      elapsedTime = difftime(currentTime, beginningTime);

      if(reserveElapsedTimeHistory)
      {
         elapsedTimeHistory[epoch] = elapsedTime;                             
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
            std::cout << "Final gradient norm: " << gradientNorm << std::endl;

            objectiveFunctional->print();
         }

         resizeTrainingHistory(1+epoch);

         break;
      }

      // Gradient norm goal 

      if(gradientNorm <= gradientNormGoal)
      {
         if(display)
         {
            std::cout << std::endl
                      << "Epoch " << epoch << ": "
                      << "Gradient norm goal reached."
                      << std::endl;  

            std::cout << "Elapsed time: " << elapsedTime << ";" << std::endl;
            std::cout << "Final parameters norm: " << freeParametersNorm << std::endl;
            std::cout << "Final evaluation: " << evaluation << std::endl;
            std::cout << "Final gradient norm: " << gradientNorm << std::endl;

            objectiveFunctional->print();
         }

         resizeTrainingHistory(1+epoch);
 
         break;
      }

      // Minimum evaluation improvement 

      double improvement = fabs(evaluation - oldEvaluation);

      if(improvement <= minimumImprovement)
      {
         if(display)
         {
            std::cout << std::endl
                      << "Epoch " << epoch << ": "
                      << "Minimum evaluation improvement reached."
                      << std::endl;

            std::cout << "Elapsed time: " << elapsedTime << ";" << std::endl;
            std::cout << "Final parameters norm: " << freeParametersNorm << std::endl;
            std::cout << "Final evaluation: " << evaluation << std::endl;
            std::cout << "Final gradient norm: " << gradientNorm << std::endl;

            objectiveFunctional->print();
         }

         resizeTrainingHistory(1+epoch);

       break;
      }

      // Maximum optimization time

      time(&currentTime);

      elapsedTime = difftime(currentTime, beginningTime);

      if(elapsedTime >= maximumTime)
      {
         if(display)
         {
            std::cout << std::endl
                      << "Epoch " << epoch << ": "
                      << "Maximum training time reached."
                      << std::endl;

            std::cout << "Elapsed time: " << elapsedTime << ";" << std::endl;
            std::cout << "Final parameters norm: " << freeParametersNorm << std::endl;
            std::cout << "Final evaluation: " << evaluation << std::endl;
            std::cout << "Final gradient norm: " << gradientNorm << std::endl;

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
                      << "Epoch " << epoch << ": "
                      << "Maximum number of epochs reached."
                      << std::endl;

            std::cout << "Elapsed time: " << elapsedTime << ";" << std::endl;
            std::cout << "Final parameters norm: " << freeParametersNorm << std::endl;
            std::cout << "Final evaluation: " << evaluation << std::endl;
            std::cout << "Final gradient norm: " << gradientNorm << std::endl;

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
         std::cout << "Gradient norm: " << gradientNorm << std::endl;
         
         objectiveFunctional->print();
      }

      // Update

      oldTrainRate = trainRate;
   } 
}


// void print(void) method

/// This method prints to the screen the training parameters, the stopping criteria and other user stuff 
/// concerning the gradient descent object:
///
/// Training operators:
/// <ul>
/// <li> Train rate method.
/// </ul>
///
/// Training parameters:
/// <ul>
/// <li> Initial train rate.
/// <li> Train rate tolerance.
/// </ul>
///
/// Stopping criteria:
/// <ul> 
/// <li> Evaluation goal.
/// <li> Gradient norm goal.
/// <li> Maximum time.
/// <li> Maximum number of epochs. 
/// </ul> 
///  
/// User stuff: 
/// <ul>
/// <li> Warning train rate.
/// <li> Error train rate.
/// <li> Display.
/// <li> Display period.
/// <li> Reserve elapsed time history.
/// <li> Reserve free parameters history.
/// <li> Reserve free parameters norm history.
/// <li> Reserve evaluation history. 
/// <li> Reserve gradient history. 
/// <li> Reserve gradient norm history.
/// <li> Reserve training direction history. 
/// <li> Reserve training direction norm history.
/// <li> Reserve training rate history.
/// </ul>

void GradientDescent::print(void)
{
   std::cout << std::endl
             << "Gradient Descent Object." << std::endl;

   // Training operators

   // Train rate method

   std::cout << "Train rate method:" << std::endl;

   switch(trainRateMethod)
   {
      case Fixed:

         std::cout << "Fixed" << std::endl;

      break;

      case GoldenSection:

         std::cout << "Golden section" << std::endl;

      break;

      case BrentMethod:

         std::cout << "Brent Method" << std::endl;

      break;
   }


   // Training parameters

   std::cout << "First train rate: " << std::endl
             << firstTrainRate << std::endl
             << "Train rate tolerance: " << std::endl
             << trainRateTolerance << std::endl;

   // Stopping criteria

   std::cout << "Evaluation goal:" << std::endl
             << evaluationGoal << std::endl
             << "Gradient norm goal:" << std::endl
             << gradientNormGoal << std::endl
             << "Maximum time:" << std::endl
             << maximumTime << std::endl 
             << "Maximum number of epochs:" << std::endl
             << maximumNumberOfEpochs << std::endl; 

   // User stuff

   std::cout << "Warning train rate:" << std::endl
             << warningTrainRate << std::endl
             << "Error train rate:" << std::endl
             << errorTrainRate << std::endl
             << "Display:" << std::endl
             << display << std::endl
             << "Display period:" << std::endl
             << displayPeriod << std::endl;

   std::cout << "Reserve elapsed time history:" << std::endl
             << reserveElapsedTimeHistory << std::endl 
             << "Reserve free parameters history:" << std::endl
             << reserveFreeParametersHistory << std::endl 
             << "Reserve free parameters norm history:" << std::endl
             << reserveFreeParametersNormHistory << std::endl 
             << "Reserve evaluation history:" << std::endl
             << reserveEvaluationHistory << std::endl 
             << "Reserve gradient history:" << std::endl
             << reserveGradientHistory << std::endl 
             << "Reserve gradient norm history:" << std::endl
             << reserveGradientNormHistory << std::endl 
             << "Reserve training direction history:" << std::endl
             << reserveTrainingDirectionHistory << std::endl 
             << "Reserve training direction norm history:" << std::endl
             << reserveTrainingDirectionNormHistory << std::endl 
             << "Reserve training rate history:"  << std::endl
             << reserveTrainingRateHistory << std::endl;

}


// void save(char*) method

/// This method saves the gradient descent object to a data file. 
///
/// Training operators:
/// <ul>
/// <li> Train rate method.
/// </ul>
///
/// Training parameters:
/// <ul>
/// <li> First train rate.
/// <li> Train rate tolerance.
/// </ul>
///
/// Stopping criteria:
/// <ul> 
/// <li> Evaluation goal.
/// <li> Gradient norm goal.
/// <li> Maximum time.
/// <li> Maximum number of epochs. 
/// </ul> 
///  
/// User stuff: 
/// <ul>
/// <li> Warning train rate.
/// <li> Error train rate.
/// <li> Display.
/// <li> Display period.
/// <li> Reserve elapsed time history.
/// <li> Reserve free parameters history.
/// <li> Reserve free parameters norm history.
/// <li> Reserve evaluation history. 
/// <li> Reserve gradient history. 
/// <li> Reserve gradient norm history.
/// <li> Reserve training direction history. 
/// <li> Reserve training direction norm history.
/// <li> Reserve training rate history.
/// </ul>
///
/// @param filename Filename.
///
/// @see load(char*).

void GradientDescent::save(char* filename)
{
   // File

   std::fstream file;

   file.open(filename, std::ios::out);

   if(!file.is_open())
   {
      std::cerr << std::endl 
                << "Flood Error: GradientDescent class." << std::endl
                << "void save(char*) method." << std::endl
                << "Cannot open gradient descent object data file."  << std::endl
                << std::endl;

      exit(1);
   }
   else
   {
      if(display)
      {
         std::cout << std::endl
                   << "Saving gradient descent object to data file..." << std::endl;
      }
   }

   // Write file header


   file << "% Flood Neural Network. Gradient Descent Object." << std::endl;


   // Training operators

   // Train direction method

   // Train rate method

   file << "TrainRateMethod:" << std::endl;

   switch(trainRateMethod)
   {
      case Fixed:

      file << "Fixed" << std::endl;

      break;

      case GoldenSection:

      file << "GoldenSection" << std::endl;

      break;

      case BrentMethod:

      file << "BrentMethod" << std::endl;

      break;
   }

   // Training parameters

   file << "FirstTrainRate:" << std::endl
        << firstTrainRate << std::endl
        << "TrainRateTolerance:" << std::endl
        << trainRateTolerance << std::endl;


   // Stopping criteria

   file << "EvaluationGoal:" << std::endl
        << evaluationGoal << std::endl
        << "GradientNormGoal:" << std::endl
        << gradientNormGoal << std::endl
        << "MaximumTime: " << std::endl
        << maximumTime << std::endl
        << "MaximumNumberOfEpochs: " << std::endl
        << maximumNumberOfEpochs << std::endl;

   // User stuff

   file << "WarningTrainRate:" << std::endl
        << warningTrainRate << std::endl
        << "ErrorTrainRate:" << std::endl
        << errorTrainRate << std::endl
        << "Display:" << std::endl
        << display << std::endl
        << "DisplayPeriod:" << std::endl
        << displayPeriod << std::endl;

   file << "ReserveElapsedTimeHistory:" << std::endl
        << reserveElapsedTimeHistory << std::endl 
        << "ReserveFreeParametersHistory:" << std::endl
        << reserveFreeParametersHistory << std::endl 
        << "ReserveFreeParametersNormHistory:" << std::endl
        << reserveFreeParametersNormHistory << std::endl 
        << "ReserveEvaluationHistory:" << std::endl
        << reserveEvaluationHistory << std::endl 
        << "ReserveGradientHistory:" << std::endl
        << reserveGradientHistory << std::endl 
        << "ReserveGradientNormHistory:" << std::endl
        << reserveGradientNormHistory << std::endl 
        << "ReserveTrainingDirectionHistory:" << std::endl
        << reserveTrainingDirectionHistory << std::endl 
        << "ReserveTrainingDirectionNormHistory:" << std::endl
        << reserveTrainingDirectionNormHistory << std::endl 
        << "ReserveTrainingRateHistory:"  << std::endl
        << reserveTrainingRateHistory << std::endl;

   file.close();
}


// void load(char*) method

/// This method loads a gradient descent object from a data file. 
/// Please mind about the file format, wich is specified in the User's Guide. 
///
/// Training operators:
/// <ul>
/// <li> Train rate method.
/// </ul>
///
/// Training parameters:
/// <ul>
/// <li> Initial train rate.
/// <li> Train rate tolerance.
/// </ul>
///
/// Stopping criteria:
/// <ul> 
/// <li> Evaluation goal.
/// <li> Gradient norm goal.
/// <li> Maximum time.
/// <li> Maximum number of epochs. 
/// </ul> 
///  
/// User stuff: 
/// <ul>
/// <li> Warning train rate.
/// <li> Error train rate.
/// <li> Display.
/// <li> Display period.
/// <li> Reserve elapsed time history.
/// <li> Reserve free parameters history.
/// <li> Reserve free parameters norm history.
/// <li> Reserve evaluation history. 
/// <li> Reserve gradient history. 
/// <li> Reserve gradient norm history.
/// <li> Reserve training direction history. 
/// <li> Reserve training direction norm history.
/// <li> Reserve training rate history.
/// </ul>
///
/// @param filename Filename.
///
/// @see save(char*).

void GradientDescent::load(char* filename)
{
   std::string word;

   // File

   std::fstream file;

   file.open(filename, std::ios::in);

   if(!file.is_open())
   {
      std::cerr << std::endl
                << "Flood Error: GradientDescent class." << std::endl
                << "void load(char*) method."
                << std::endl
                << "Cannot open gradient descent object data file."  << std::endl;

      exit(1);
   }
   else
   {
      if(display)
      { 
         std::cout << std::endl
                   << "Loading gradient descent object from data file..."
                   << std::endl;
      }
   }

   // Train rate method

   file >> word;

   file >> word;

   if(word == "Fixed")
   {
      trainRateMethod = Fixed;
   }
   else if(word == "GoldenSection")
   {
      trainRateMethod = GoldenSection;
   }
   else if(word == "BrentMethod")
   {
      trainRateMethod = BrentMethod;
   }
   else
   {
      std::cerr << std::endl
                << "Flood Error: GradientDescent class." << std::endl
                << "void load(char*) method." << std::endl
                << "Unknown file format."  << std::endl;

      exit(1);   
   }

   // Training parameters

   // First train rate

   file >> word;

   if(word != "FirstTrainRate:")
   {
      std::cerr << std::endl
                << "Flood Error: GradientDescent class." << std::endl
                << "void load(char*) method." << std::endl
                << "Unknown file format."  << std::endl;

      exit(1);   
   }

   file >> firstTrainRate;

   // Tolerance

   file >> word;

   if(word != "TrainRateTolerance:")
   {
      std::cerr << std::endl
                << "Flood Error: GradientDescent class." << std::endl
                << "void load(char*) method." << std::endl
                << "Unknown file format."  << std::endl;

      exit(1);   
   }

   file >> trainRateTolerance;

   // Stopping criteria: 

   // Evaluation goal

   file >> word;

   if(word != "EvaluationGoal:")
   {
      std::cerr << std::endl
                << "Flood Error: GradientDescent class." << std::endl
                << "void load(char*) method." << std::endl
                << "Unknown file format."  << std::endl;

      exit(1);   
   }

   file >> evaluationGoal;
   
   // Gradient norm goal

   file >> word;

   if(word != "GradientNormGoal:")
   {
      std::cerr << std::endl
                << "Flood Error: GradientDescent class." << std::endl
                << "void load(char*) method." << std::endl
                << "Unknown file format."  << std::endl;

      exit(1);   
   }

   file >> gradientNormGoal;

   // Maximum time

   file >> word;

   if(word != "MaximumTime:")
   {
      std::cerr << std::endl
                << "Flood Error: GradientDescent class." << std::endl
                << "void load(char*) method." << std::endl
                << "Unknown file format."  << std::endl;

      exit(1);   
   }

   file >> maximumTime;

   // Maximum number of epochs

   file >> word;

   if(word != "MaximumNumberOfEpochs:")
   {
      std::cerr << std::endl
                << "Flood Error: GradientDescent class." << std::endl
                << "void load(char*) method." << std::endl
                << "Unknown file format."  << std::endl;

      exit(1);   
   }

   file >> maximumNumberOfEpochs;

   // User stuff: 

   // Warning train rate

   file >> word;

   if(word != "WarningTrainRate:")
   {
      std::cerr << std::endl
                << "Flood Error: GradientDescent class." << std::endl
                << "void load(char*) method." << std::endl
                << "Unknown file format."  << std::endl;

      exit(1);   
   }

   file >> warningTrainRate;

   // Error train rate

   file >> word;

   if(word != "ErrorTrainRate:")
   {
      std::cerr << std::endl
                << "Flood Error: GradientDescent class." << std::endl
                << "void load(char*) method." << std::endl
                << "Unknown file format."  << std::endl;

      exit(1);   
   }

   file >> errorTrainRate;

   // Display

   file >> word;

   if(word != "Display:")
   {
      std::cerr << std::endl
                << "Flood Error: GradientDescent class." << std::endl
                << "void load(char*) method." << std::endl
                << "Unknown file format."  << std::endl;

      exit(1);   
   }

   file >> display;

   // Display period

   file >> word;

   if(word != "DisplayPeriod:")
   {
      std::cerr << std::endl
                << "Flood Error: GradientDescent class." << std::endl
                << "void load(char*) method." << std::endl
                << "Unknown file format."  << std::endl;

      exit(1);   
   }

   file >> displayPeriod;

   // Reserve elapsed time history

   file >> word;

   if(word != "ReserveElapsedTimeHistory:")
   {
      std::cerr << std::endl
                << "Flood Error: GradientDescent class." << std::endl
                << "void load(char*) method." << std::endl
                << "Unknown file format."  << std::endl;

      exit(1);   
   }

   file >> reserveElapsedTimeHistory; 

   // Reserve free parameters history

   file >> word;

   if(word != "ReserveFreeParametersHistory:")
   {
      std::cerr << std::endl
                << "Flood Error: GradientDescent class." << std::endl
                << "void load(char*) method." << std::endl
                << "Unknown file format."  << std::endl;

      exit(1);   
   }

   file >> reserveFreeParametersHistory;

   // Reserve free parameters norm history

   file >> word;

   if(word != "ReserveFreeParametersNormHistory:")
   {
      std::cerr << std::endl
                << "Flood Error: GradientDescent class." << std::endl
                << "void load(char*) method." << std::endl
                << "Unknown file format."  << std::endl;

      exit(1);   
   }

   file >> reserveFreeParametersNormHistory;

   // Reserve evaluation history

   file >> word;

   if(word != "ReserveEvaluationHistory:")
   {
      std::cerr << std::endl
                << "Flood Error: GradientDescent class." << std::endl
                << "void load(char*) method." << std::endl
                << "Unknown file format."  << std::endl;

      exit(1);   
   }

   file >> reserveEvaluationHistory;

   // Reserve gradient history

   file >> word;

   if(word != "ReserveGradientHistory:")
   {
      std::cerr << std::endl
                << "Flood Error: GradientDescent class." << std::endl
                << "void load(char*) method." << std::endl
                << "Unknown file format."  << std::endl;

      exit(1);   
   }

   file >> reserveGradientHistory;

   // Reserve gradient norm history

   file >> word;

   if(word != "ReserveGradientNormHistory:")
   {
      std::cerr << std::endl
                << "Flood Error: GradientDescent class." << std::endl
                << "void load(char*) method." << std::endl
                << "Unknown file format."  << std::endl;

      exit(1);   
   }

   file >> reserveGradientNormHistory;

   // Reserve training direction history

   file >> word;

   if(word != "ReserveTrainingDirectionHistory:")
   {
      std::cerr << std::endl
                << "Flood Error: GradientDescent class." << std::endl
                << "void load(char*) method." << std::endl
                << "Unknown file format."  << std::endl;

      exit(1);   
   }

   file >> reserveTrainingDirectionHistory;

   // Reserve training direction norm history

   file >> word;

   if(word != "ReserveTrainingDirectionNormHistory:")
   {
      std::cout << std::endl
                << "Flood Error: GradientDescent class." << std::endl
                << "void load(char*) method." << std::endl
                << "Unknown file format."  << std::endl;

      exit(1);   
   }

   file >> reserveTrainingDirectionNormHistory;

   // Reserve training rate history

   file >> word;

   if(word != "ReserveTrainingRateHistory:")
   {
      std::cerr << std::endl
                << "Flood Error: GradientDescent class." << std::endl
                << "void load(char*) method." << std::endl
                << "Unknown file format."  << std::endl;

      exit(1);   
   }

   file >> reserveTrainingRateHistory;

   // Close file

   file.close();
}


// void resizeTrainingHistory(int) method

/// This method resizes the vectors or matrices containing training history information to a new size:
///
/// <ul>
/// <li> Elapsed time history vector.
/// <li> Free parameters history matrix.
/// <li> Free parameters norm history vector. 
/// <li> Evaluation history vector.
/// <li> Gradient history matrix.
/// <li> Gradient norm history vector. 
/// <li> Training direction history matrix.
/// <li> Training direction norm history vector
/// <li> Training rate history vector. 
/// </ul>
///
/// @param newSize Size of training history. 

void GradientDescent::resizeTrainingHistory(int newSize)
{
   // Elapsed time history vector

   if(reserveElapsedTimeHistory)
   {
      elapsedTimeHistory.resize(newSize);
   }

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

   // Gradient history matrix

   if(reserveGradientHistory)
   {
      MultilayerPerceptron* multilayerPerceptron = objectiveFunctional->getMultilayerPerceptron();
                                        
      int numberOfFreeParameters = multilayerPerceptron->getNumberOfFreeParameters();

      gradientHistory.resize(newSize, numberOfFreeParameters);
   }
 
   // Gradient norm history vector

   if(reserveGradientNormHistory)
   {
      gradientNormHistory.resize(newSize);
   }

   // Training direction history matrix

   if(reserveTrainingDirectionHistory)
   {
      MultilayerPerceptron* multilayerPerceptron = objectiveFunctional->getMultilayerPerceptron();
                                        
      int numberOfFreeParameters = multilayerPerceptron->getNumberOfFreeParameters();

      trainingDirectionHistory.resize(newSize, numberOfFreeParameters);
   }

   // Training direction norm history vector

   if(reserveTrainingDirectionNormHistory)
   {
      trainingDirectionNormHistory.resize(newSize);
   }

   // Training rate history vector

   if(reserveTrainingRateHistory)
   {
      trainingRateHistory.resize(newSize);
   }
}


// void saveTrainingHistory(char*) method

/// This method saves the training history to a data file. 
///
/// @param filename Training history filename. 

void GradientDescent::saveTrainingHistory(char* filename)

{
   std::fstream file; 

   file.open(filename, std::ios::out);

   // Write file header 

   if(!file.is_open())
   {
      std::cerr << std::endl 
                << "Flood Error: GradientDescent class." << std::endl
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

   file << "% Flood Neural Network. Gradient Descent Training History." << std::endl;

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
   if(reserveGradientHistory)
   {
      file << "GradientHistory:" << std::endl;      file << gradientHistory << std::endl;      
   }
   if(reserveGradientNormHistory)
   {
      file << "GradientNormHistory:" << std::endl;      file << gradientNormHistory << std::endl;      
   }
   if(reserveTrainingDirectionHistory)
   {
      file << "TrainingDirectionHistory:" << std::endl;      file << trainingDirectionHistory << std::endl;      
   }
   if(reserveTrainingDirectionNormHistory)
   {
      file << "TrainingDirectionNormHistory:" << std::endl;      file << trainingDirectionNormHistory << std::endl;      
   }
   if(reserveTrainingRateHistory)
   {
      file << "TrainingRateHistory:" << std::endl;      file << trainingRateHistory << std::endl;      
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
