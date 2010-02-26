/****************************************************************************************************************/
/*                                                                                                              */
/*   Flood: An Open Source Neural Networks C++ Library                                                          */
/*   www.cimne.com/flood                                                                                        */
/*                                                                                                              */
/*   O R D I N A R Y   D I F F E R E N T I A L   E Q U A T I O N S   C L A S S                                  */
/*                                                                                                              */ 
/*   Roberto Lopez                                                                                              */ 
/*   International Center for Numerical Methods in Engineering (CIMNE)                                          */
/*   Technical University of Catalonia (UPC)                                                                    */
/*   Barcelona, Spain                                                                                           */
/*   E-mail: rlopez@cimne.upc.edu                                                                               */ 
/*                                                                                                              */
/****************************************************************************************************************/

#include "OrdinaryDifferentialEquations.h"

#include<iostream>
#include<fstream>
#include<limits>
#include<math.h>

namespace Flood
{

// GENERAL CONSTRUCTOR

/// General constructor. 

OrdinaryDifferentialEquations::OrdinaryDifferentialEquations(void)
{   
   // Runge-Kutta method parameters

   numberOfPoints = 1001;

   // Runge-Kutta-Fehlberg method parameters

   tolerance = 1.0e-9;

   initialSize = 1000;
   warningSize = 1000000;
   errorSize = 1000000000;

   display = true;
}


// DESTRUCTOR

/// Destructor.

OrdinaryDifferentialEquations::~OrdinaryDifferentialEquations(void)
{

}


// METHODS

// GET METHODS

// int getNumberOfPoints(void) method

/// This method returns the number of integration points to be used in the Runge-Kutta method. 

int OrdinaryDifferentialEquations::getNumberOfPoints(void)
{
   return(numberOfPoints);
}


// double getTolerance(void) method

/// This method returns the tolerance to be used in the Runge-Kutta-Fehlberg method. 

double OrdinaryDifferentialEquations::getTolerance(void)
{
   return(tolerance);
}


// int getInitialSize(void) method

/// This method returns the initial size of the solution vectors for the Runge-Kutta-Fehlberg method. 
/// If during integration the number of points required is greater, the algorithm resizes that vectors. 

int OrdinaryDifferentialEquations::getInitialSize(void)
{
   return(initialSize);
}


// int getWarningSize(void) method

/// This method returns the size for the solution vectors at wich a warning message is written to the screen 
/// during Runge-Kutta-Fehlberg integration.
///
/// @see getErrorSize(void).

int OrdinaryDifferentialEquations::getWarningSize(void)
{
   return(warningSize);
}


// int getErrorSize(void) method

/// This method returns the size of the solution vectors at wich the Runge-Kutta-Fehlberg method throws an error.
///
/// @see getWarningSize(void).

int OrdinaryDifferentialEquations::getErrorSize(void)
{
   return(errorSize);
}


// bool getDisplay(void) method

/// This method returns true if messages from this class are to be displayed on the screen, or false if messages 
/// from this class are not to be displayed on the screen.

bool OrdinaryDifferentialEquations::getDisplay(void)
{
   return(display);
}


// SET METHODS

// void setNumberOfPoints(int) method

/// This method returns sets a new number of integration points to be used in the Runge-Kutta method. 
///
/// @param newNumberOfPoints Number of integration points.

void OrdinaryDifferentialEquations::setNumberOfPoints(int newNumberOfPoints)
{
   numberOfPoints = newNumberOfPoints;
}


// void setTolerance(double) method

/// This method sets a new tolerance to be used in the Runge-Kutta-Fehlberg method. 
///
/// @param newTolerance Tolerance in the integration error. 

void OrdinaryDifferentialEquations::setTolerance(double newTolerance)
{
   tolerance = newTolerance;
}


// void setInitialSize(int) method

/// This method returns sets a new initial size of the solution vectors for the Runge-Kutta-Fehlberg method. 
/// If during integration the number of points required is greater, the algorithm resizes that vectors. 
///
/// @param newInitialSize Initial size for the solution vectors. 

void OrdinaryDifferentialEquations::setInitialSize(int newInitialSize)
{
   initialSize = newInitialSize;
}


// void setWarningSize(int) method

/// This method a new size for the solution vectors at wich a warning message is written to the screen during 
/// Runge-Kutta-Fehlberg integration.
///
/// @see setErrorSize(int).

void OrdinaryDifferentialEquations::setWarningSize(int newWarningSize)
{
   warningSize = newWarningSize;
}


// void setErrorSize(int) method

/// This method returns sets a new size for the solution vectors at wich the Runge-Kutta-Fehlberg method throws 
/// an error.
///
/// @see setWarningSize(int).

void OrdinaryDifferentialEquations::setErrorSize(int newErrorSize)
{
   errorSize = newErrorSize;
}


// void setDisplay(bool) method

/// This method sets a new display value. 
/// If it is set to true messages from this class are to be displayed on the screen;
/// if it is set to false messages from this class are not to be displayed on the screen.
///
/// @param newDisplay Display value.

void OrdinaryDifferentialEquations::setDisplay(bool newDisplay)
{
   display = newDisplay;
}


// void OrdinaryDifferentialEquations::getRungeKuttaIntegral(
// Vector<double>& x, Vector<double>& y, 
// double (*f)(double, double), 
// double a, double b, double ya, int n)

/// This method integrates the differential equation
///
/// <ul>
/// <li>  y' = f(x, y),
/// </ul>
///
/// where y' is given as a C-style function, from x=a to x=b with initial condition y(a)=ya. 
/// The method is based on an explicit fourth order Runge-Kutta formula.  
///
/// @param x Pointer to vector with x values. 
/// @param y Pointer to vector with y values.
/// @param f Pointer to C-style function f (state equation for variable y).
/// @param a Lower integration limit.
/// @param b Upper integration limit.
/// @param ya Initial condition for variable y.

void OrdinaryDifferentialEquations::getRungeKuttaIntegral(
Vector<double>& x, 
Vector<double>& y, 
double (*f)(double, double), 
double a, double b, double ya)
{
   // Integration step

   double h = (b-a)/(numberOfPoints-1.0);      

   // Fourth order Runge-Kutta coefficients

   Vector<double> k(4);

   // Independent variable

   x[0] = a;

   for(int i = 0;  i < numberOfPoints-2; i++)
   {
      x[i+1] = x[i] + h;
   }

   x[numberOfPoints-1] = b;

   // Dependent variable

   y[0] = ya;

   for(int i = 0;  i < numberOfPoints-1; i++)
   {   
      // Runge-Kutta coefficients

      k[0] = f(x[i], y[i]);
      k[1] = f(x[i]+h/2.0, y[i]+h*k[0]/2.0);
      k[2] = f(x[i]+h/2.0, y[i]+h*k[1]/2.0);
      k[3] = f(x[i+1], y[i]+h*k[2]);

      // Dependent variable

      y[i+1] = y[i] + h*(k[0] + 2.0*k[1] + 2.0*k[2] + k[3])/6.0;
   }
}


// void getRungeKuttaIntegral(
// double (*f1)(double, double, double), 
// double (*f2)(double, double,double),
// Vector<double>&, Vector<double>&, Vector<double>&, 
// double, double, double) method
//
/// This method integrates the system of two differential equations 
///
/// <ul>
/// <li> y1' = f1(x, y1, y2),
/// <li> y2' = f2(x, y1, y2),
/// </ul>
///
/// where y1' and y2' are given as C-style functions, from x=a to x=b with initial conditions y1(a)=y1a and 
/// y2(a)=y2a. The method is based on an 
/// explicit fourth order Runge-Kutta formula.  
///
/// @param x Pointer to vector with x values. 
/// @param y1 Pointer to vector with y1 values.
/// @param y2 Pointer to vector with y2 values.
/// @param f1 Pointer to C-style function f1' (state equation for variable y1).
/// @param f2 Pointer to C-style function f2' (state equation for variable y2).
/// @param a Lower integration limit.
/// @param b Upper integration limit.
/// @param y1a Initial condition for variable y1.
/// @param y2a Initial condition for variable y2.

void OrdinaryDifferentialEquations::getRungeKuttaIntegral(
Vector<double>& x, 
Vector<double>& y1, 
Vector<double>& y2, 
double (*f1)(double, double,double), 
double (*f2)(double, double,double), 
double a, double b, 
double y1a, double y2a)
{
   // Integration step

   double h = (b-a)/(numberOfPoints-1.0);      

   // Fourth order Runge-Kutta coefficients

   Vector<double> k1(4);
   Vector<double> k2(4);

   // Independent variable

   x[0] = a;

   for(int i = 0;  i < numberOfPoints-2; i++)
   {
      x[i+1] = x[i] + h;
   }

   x[numberOfPoints-1] = b;

   // Dependent variables
  
   y1[0] = y1a;
   y2[0] = y2a;

   for(int i = 0;  i < numberOfPoints-1; i++)
   {
      // Runge-Kutta coefficients

      k1[0] = f1(x[i], y1[i], y2[i]);
      k2[0] = f2(x[i], y1[i], y2[i]);

      k1[1] = f1(x[i]+h/2.0, y1[i]+h*k1[0]/2.0, y2[i]+h*k2[0]/2.0);
      k2[1] = f2(x[i]+h/2.0, y1[i]+h*k1[0]/2.0, y2[i]+h*k2[0]/2.0);

      k1[2] = f1(x[i]+h/2.0, y1[i]+h*k1[1]/2.0, y2[i]+h*k2[1]/2.0);
      k2[2] = f2(x[i]+h/2.0, y1[i]+h*k1[1]/2.0, y2[i]+h*k2[1]/2.0);

      k1[3] = f1(x[i+1], y1[i]+h*k1[2], y2[i]+h*k2[2]);
      k2[3] = f2(x[i+1], y1[i]+h*k1[2], y2[i]+h*k2[2]);

      // Dependent variables

      y1[i+1] = y1[i] + h*(k1[0] + 2.0*k1[1] + 2.0*k1[2] + k1[3])/6.0;
      y2[i+1] = y2[i] + h*(k2[0] + 2.0*k2[1] + 2.0*k2[2] + k2[3])/6.0;
   }
}


// void getRungeKuttaIntegral(
// double (*f1)(double, double, double, double), 
// double (*f2)(double, double, double, double),
// double (*f3)(double, double, double, double),   
// Vector<double>&, Vector<double>&, Vector<double>&, Vector<double>&, 
// double, double, double, double, double, int) method
//
/// This method integrates the system of three differential equations 
///
/// <ul>
/// <li> y1' = f1(x, y1, y2, y3),
/// <li> y2' = f2(x, y1, y2, y3),
/// <li> y3' = f3(x, y1, y2, y3),
/// </ul>
///
/// where y1', y2' and y3' are given as C-style functions, from x=a to x=b
/// with initial conditions y1(a)=y1a, y2(a)=y2a and 
/// y3(a)=y3a. 
///
/// The method is based on an explicit fourth order Runge-Kutta formula.  
///
/// @param x Pointer to vector with x values. 
/// @param y1 Pointer to vector with y1 values.
/// @param y2 Pointer to vector with y2 values.
/// @param y3 Pointer to vector with y3 values.
/// @param f1 Pointer to C-style function f1 (state equation for variable y1).
/// @param f2 Pointer to C-style function f2 (state equation for variable y2).
/// @param f3 Pointer to C-style function f3 (state equation for variable y3).
/// @param a Lower integration limit.
/// @param b Upper integration limit.
/// @param y1a Initial condition for variable y1.
/// @param y2a Initial condition for variable y2.
/// @param y3a Initial condition for variable y3.

void OrdinaryDifferentialEquations::getRungeKuttaIntegral(
Vector<double> &x, 
Vector<double>& y1, 
Vector<double>& y2, 
Vector<double>& y3, 
double (*f1)(double, double, double, double), 
double (*f2)(double, double, double, double), 
double (*f3)(double, double, double, double), 
double a, double b, 
double y1a, double y2a, double y3a)
{
   // Integration step

   double h = (b-a)/(numberOfPoints-1.0);      

   // Fourth order Runge-Kutta coefficients

   Vector<double> k1(4);
   Vector<double> k2(4);
   Vector<double> k3(4);

   // Independent variable

   x[0] = a;

   for(int i = 0;  i < numberOfPoints-2; i++)
   {
      x[i+1] = x[i] + h;
   }

   x[numberOfPoints-1] = b;

   // Dependent variables

   y1[0] = y1a;
   y2[0] = y2a;
   y3[0] = y3a;

   for(int i = 0;  i < numberOfPoints-1; i++)
   {
      // Runge-Kutta coefficients

      k1[0] = f1(x[i], y1[i], y2[i], y3[i]);
      k2[0] = f2(x[i], y1[i], y2[i], y3[i]);
      k3[0] = f3(x[i], y1[i], y2[i], y3[i]);

      k1[1] = f1(x[i]+h/2.0, y1[i]+h*k1[0]/2.0, y2[i]+h*k2[0]/2.0, y3[i]+h*k3[0]/2.0);
      k2[1] = f2(x[i]+h/2.0, y1[i]+h*k1[0]/2.0, y2[i]+h*k2[0]/2.0, y3[i]+h*k3[0]/2.0);
      k3[1] = f3(x[i]+h/2.0, y1[i]+h*k1[0]/2.0, y2[i]+h*k2[0]/2.0, y3[i]+h*k3[0]/2.0);

      k1[2] = f1(x[i]+h/2.0, y1[i]+h*k1[1]/2.0, y2[i]+h*k2[1]/2.0, y3[i]+h*k3[1]/2.0);
      k2[2] = f2(x[i]+h/2.0, y1[i]+h*k1[1]/2.0, y2[i]+h*k2[1]/2.0, y3[i]+h*k3[1]/2.0);
      k3[2] = f3(x[i]+h/2.0, y1[i]+h*k1[1]/2.0, y2[i]+h*k2[1]/2.0, y3[i]+h*k3[1]/2.0);

      k1[3] = f1(x[i+1], y1[i]+h*k1[2], y2[i]+h*k2[2], y3[i]+h*k3[2]);
      k2[3] = f2(x[i+1], y1[i]+h*k1[2], y2[i]+h*k2[2], y3[i]+h*k3[2]);
      k3[3] = f2(x[i+1], y1[i]+h*k1[2], y2[i]+h*k2[2], y3[i]+h*k3[2]);

      // Dependent variables

      y1[i+1] = y1[i] + h*(k1[0] + 2.0*k1[1] + 2.0*k1[2] + k1[3])/6.0;
      y2[i+1] = y2[i] + h*(k2[0] + 2.0*k2[1] + 2.0*k2[2] + k2[3])/6.0;
      y3[i+1] = y3[i] + h*(k3[0] + 2.0*k3[1] + 2.0*k3[2] + k3[3])/6.0;
   }
}


// void getRungeKuttaIntegral(
// double (*f1)(double, double, double, double, double), 
// double (*f2)(double, double, double, double, double),
// double (*f3)(double, double, double, double, double),   
// double (*f4)(double, double, double, double, double),   
// Vector<double>&, 
// Vector<double>&, Vector<double>&, Vector<double>&, Vector<double>&, 
// double, double, double, double, double, double, int) method

/// This method integrates the system of four differential equations 
///   
/// <ul>
/// <li> y1' = f1(x, y1, y2, y3, y4)
/// <li> y2' = f2(x, y1, y2, y3, y4)
/// <li> y3' = f3(x, y1, y2, y3, y4)
/// <li> y4' = f4(x, y1, y2, y3, y4)
/// </ul>
/// 
/// where y1', y2', y3' and y4' are given as C-style functions, from x=a 
/// to x=b with initial conditions y1(a)=y1a, y2(a)=y2a, y3(a)=y3a and
/// y4(a)=y4a. The method is based on an explicit fourth order Runge-Kutta
/// formula.  
///
/// @param x Pointer to vector with x values. 
/// @param y1 Pointer to vector with y1 values.
/// @param y2 Pointer to vector with y2 values.
/// @param y3 Pointer to vector with y3 values.
/// @param y4 Pointer to vector with y4 values.
/// @param f1 Pointer to C-style function f1 (state equation for variable y1).
/// @param f2 Pointer to C-style function f2 (state equation for variable y2).
/// @param f3 Pointer to C-style function f3 (state equation for variable y3).
/// @param f4 Pointer to C-style function f4 (state equation for variable y4).
/// @param a Lower integration limit.
/// @param b Upper integration limit.
/// @param y1a Initial condition for variable y1.
/// @param y2a Initial condition for variable y2.
/// @param y3a Initial condition for variable y3.
/// @param y4a Initial condition for variable y4.

void OrdinaryDifferentialEquations::getRungeKuttaIntegral(
Vector<double>& x, 
Vector<double>& y1, 
Vector<double>& y2, 
Vector<double>& y3, 
Vector<double>& y4, 
double (*f1)(double, double, double, double, double), 
double (*f2)(double, double, double, double, double), 
double (*f3)(double, double, double, double, double), 
double (*f4)(double, double, double, double, double), 
double a, double b, 
double y1a, double y2a, double y3a, double y4a)
{
   // Integration step

   double h = (b-a)/(numberOfPoints-1.0);

   // Fourth order Runge-Kutta coefficients

   Vector<double> k1(4);
   Vector<double> k2(4);
   Vector<double> k3(4);
   Vector<double> k4(4);

   // Independent variable

   x[0] = a;

   for(int i = 0;  i < numberOfPoints-2; i++)
   {
      x[i+1] = x[i] + h;
   }

   x[numberOfPoints-1] = b;

   // Dependent variables

   y1[0] = y1a;
   y2[0] = y2a;
   y3[0] = y3a;
   y4[0] = y4a;

   for(int i = 0;  i < numberOfPoints-1; i++)
   {
      // Runge-Kutta coefficients

      k1[0] = f1(x[i], y1[i], y2[i], y3[i], y4[i]);
      k2[0] = f2(x[i], y1[i], y2[i], y3[i], y4[i]);
      k3[0] = f3(x[i], y1[i], y2[i], y3[i], y4[i]);
      k4[0] = f4(x[i], y1[i], y2[i], y3[i], y4[i]);

      k1[1] = f1(x[i]+h/2.0, 
      y1[i]+h*k1[0]/2.0, 
      y2[i]+h*k2[0]/2.0, 
      y3[i]+h*k3[0]/2.0, 
      y4[i]+h*k4[0]/2.0);

      k2[1] = f2(x[i]+h/2.0, 
      y1[i]+h*k1[0]/2.0, 
      y2[i]+h*k2[0]/2.0, 
      y3[i]+h*k3[0]/2.0, 
      y4[i]+h*k4[0]/2.0);

      k3[1] = f3(x[i]+h/2.0, 
      y1[i]+h*k1[0]/2.0, 
      y2[i]+h*k2[0]/2.0, 
      y3[i]+h*k3[0]/2.0, 
      y4[i]+h*k4[0]/2.0);

      k4[1] = f4(x[i]+h/2.0, 
      y1[i]+h*k1[0]/2.0, 
      y2[i]+h*k2[0]/2.0, 
      y3[i]+h*k3[0]/2.0, 
      y4[i]+h*k4[0]/2.0);

      k1[2] = f1(x[i]+h/2.0, 
      y1[i]+h*k1[1]/2.0, 
      y2[i]+h*k2[1]/2.0, 
      y3[i]+h*k3[1]/2.0, 
      y4[i]+h*k4[1]/2.0);

      k2[2] = f2(x[i]+h/2.0, 
      y1[i]+h*k1[1]/2.0, 
      y2[i]+h*k2[1]/2.0, 
      y3[i]+h*k3[1]/2.0, 
      y4[i]+h*k4[1]/2.0);

      k3[2] = f3(x[i]+h/2.0, 
      y1[i]+h*k1[1]/2.0, 
      y2[i]+h*k2[1]/2.0, 
      y3[i]+h*k3[1]/2.0, 
      y4[i]+h*k4[1]/2.0);

      k4[2] = f4(x[i]+h/2.0, 
      y1[i]+h*k1[1]/2.0, 
      y2[i]+h*k2[1]/2.0, 
      y3[i]+h*k3[1]/2.0, 
      y4[i]+h*k4[1]/2.0);

      k1[3] = f1(x[i+1], y1[i]+h*k1[2], y2[i]+h*k2[2], y3[i]+h*k3[2], y4[i]+h*k4[2]);
      k2[3] = f2(x[i+1], y1[i]+h*k1[2], y2[i]+h*k2[2], y3[i]+h*k3[2], y4[i]+h*k4[2]);
      k3[3] = f2(x[i+1], y1[i]+h*k1[2], y2[i]+h*k2[2], y3[i]+h*k3[2], y4[i]+h*k4[2]);
      k4[3] = f4(x[i+1], y1[i]+h*k1[2], y2[i]+h*k2[2], y3[i]+h*k3[2], y4[i]+h*k4[2]);

      // Dependent variables

      y1[i+1] = y1[i] + h*(k1[0] + 2.0*k1[1] + 2.0*k1[2] + k1[3])/6.0;
      y2[i+1] = y2[i] + h*(k2[0] + 2.0*k2[1] + 2.0*k2[2] + k2[3])/6.0;
      y3[i+1] = y3[i] + h*(k3[0] + 2.0*k3[1] + 2.0*k3[2] + k3[3])/6.0;
      y4[i+1] = y4[i] + h*(k4[0] + 2.0*k4[1] + 2.0*k4[2] + k4[3])/6.0;
   }
}


// void getRungeKuttaIntegral(
// double (*f1)(double, double, double, double, double), 
// double (*f2)(double, double, double, double, double),
// double (*f3)(double, double, double, double, double),   
// double (*f4)(double, double, double, double, double),   
// Vector<double>&, 
// Vector<double>&, Vector<double>&, Vector<double>&, Vector<double>&, 
// double, double, double, double, double, double) method

/// This method integrates the system of five differential equations 
///   
/// <ul>
/// <li> y1' = f1(x, y1, y2, y3, y4, y5)
/// <li> y2' = f2(x, y1, y2, y3, y4, y5)
/// <li> y3' = f3(x, y1, y2, y3, y4, y5)
/// <li> y4' = f4(x, y1, y2, y3, y4, y5)
/// <li> y5' = f5(x, y1, y2, y3, y4, y5)
/// </ul>
/// 
/// where y1', y2', y3', y4' and y5' are given as C-style functions, from x=a 
/// to x=b with initial conditions y1(a)=y1a, y2(a)=y2a, y3(a)=y3a, y4(a)=y4a and
/// y5(a)=y5a. The method is based on an explicit fourth order Runge-Kutta
/// formula.  
///
/// @param x Pointer to vector with x values. 
/// @param y1 Pointer to vector with y1 values.
/// @param y2 Pointer to vector with y2 values.
/// @param y3 Pointer to vector with y3 values.
/// @param y4 Pointer to vector with y4 values.
/// @param y5 Pointer to vector with y5 values.
/// @param f1 Pointer to C-style function f1 (state equation for variable y1).
/// @param f2 Pointer to C-style function f2 (state equation for variable y2).
/// @param f3 Pointer to C-style function f3 (state equation for variable y3).
/// @param f4 Pointer to C-style function f4 (state equation for variable y4).
/// @param f5 Pointer to C-style function f5 (state equation for variable y5).
/// @param a Lower integration limit.   
/// @param b Upper integration limit.
/// @param y1a Initial condition for variable y1.
/// @param y2a Initial condition for variable y2.
/// @param y3a Initial condition for variable y3.
/// @param y4a Initial condition for variable y4.
/// @param y5a Initial condition for variable y5.

void OrdinaryDifferentialEquations::getRungeKuttaIntegral(
Vector<double>& x, 
Vector<double>& y1, 
Vector<double>& y2, 
Vector<double>& y3, 
Vector<double>& y4, 
Vector<double>& y5, 
double (*f1)(double, double, double, double, double, double), 
double (*f2)(double, double, double, double, double, double), 
double (*f3)(double, double, double, double, double, double), 
double (*f4)(double, double, double, double, double, double), 
double (*f5)(double, double, double, double, double, double), 
double a, double b, 
double y1a, double y2a, double y3a, double y4a, double y5a)
{
   // Integration step

   double h = (b-a)/(numberOfPoints-1.0);

   // Fourth order Runge-Kutta coefficients

   Vector<double> k1(4);
   Vector<double> k2(4);
   Vector<double> k3(4);
   Vector<double> k4(4);
   Vector<double> k5(4);

   // Independent variable

   x[0] = a;

   for(int i = 0;  i < numberOfPoints-2; i++)
   {
      x[i+1] = x[i] + h;
   }

   x[numberOfPoints-1] = b;

   // Dependent variables

   y1[0] = y1a;
   y2[0] = y2a;
   y3[0] = y3a;
   y4[0] = y4a;
   y5[0] = y5a;

   for(int i = 0;  i < numberOfPoints-1; i++)
   {
      // Runge-Kutta coefficients

      k1[0] = f1(x[i], y1[i], y2[i], y3[i], y4[i], y5[i]);
      k2[0] = f2(x[i], y1[i], y2[i], y3[i], y4[i], y5[i]);
      k3[0] = f3(x[i], y1[i], y2[i], y3[i], y4[i], y5[i]);
      k4[0] = f4(x[i], y1[i], y2[i], y3[i], y4[i], y5[i]);
      k5[0] = f5(x[i], y1[i], y2[i], y3[i], y4[i], y5[i]);

      k1[1] = f1(x[i]+h/2.0, y1[i]+h*k1[0]/2.0, y2[i]+h*k2[0]/2.0, y3[i]+h*k3[0]/2.0, y4[i]+h*k4[0]/2.0, y5[i]+h*k5[0]/2.0);
      k2[1] = f2(x[i]+h/2.0, y1[i]+h*k1[0]/2.0, y2[i]+h*k2[0]/2.0, y3[i]+h*k3[0]/2.0, y4[i]+h*k4[0]/2.0, y5[i]+h*k5[0]/2.0);
      k3[1] = f3(x[i]+h/2.0, y1[i]+h*k1[0]/2.0, y2[i]+h*k2[0]/2.0, y3[i]+h*k3[0]/2.0, y4[i]+h*k4[0]/2.0, y5[i]+h*k5[0]/2.0);
      k4[1] = f4(x[i]+h/2.0, y1[i]+h*k1[0]/2.0, y2[i]+h*k2[0]/2.0, y3[i]+h*k3[0]/2.0, y4[i]+h*k4[0]/2.0, y5[i]+h*k5[0]/2.0);
      k5[1] = f5(x[i]+h/2.0, y1[i]+h*k1[0]/2.0, y2[i]+h*k2[0]/2.0, y3[i]+h*k3[0]/2.0, y4[i]+h*k4[0]/2.0, y5[i]+h*k5[0]/2.0);

      k1[2] = f1(x[i]+h/2.0, y1[i]+h*k1[1]/2.0, y2[i]+h*k2[1]/2.0, y3[i]+h*k3[1]/2.0, y4[i]+h*k4[1]/2.0, y5[i]+h*k5[1]/2.0);
      k2[2] = f2(x[i]+h/2.0, y1[i]+h*k1[1]/2.0, y2[i]+h*k2[1]/2.0, y3[i]+h*k3[1]/2.0, y4[i]+h*k4[1]/2.0, y5[i]+h*k5[1]/2.0);
      k3[2] = f3(x[i]+h/2.0, y1[i]+h*k1[1]/2.0, y2[i]+h*k2[1]/2.0, y3[i]+h*k3[1]/2.0, y4[i]+h*k4[1]/2.0, y5[i]+h*k5[1]/2.0);
      k4[2] = f4(x[i]+h/2.0, y1[i]+h*k1[1]/2.0, y2[i]+h*k2[1]/2.0, y3[i]+h*k3[1]/2.0, y4[i]+h*k4[1]/2.0, y5[i]+h*k5[1]/2.0);
      k5[2] = f5(x[i]+h/2.0, y1[i]+h*k1[1]/2.0, y2[i]+h*k2[1]/2.0, y3[i]+h*k3[1]/2.0, y4[i]+h*k4[1]/2.0, y5[i]+h*k5[1]/2.0);

      k1[3] = f1(x[i+1], y1[i]+h*k1[2], y2[i]+h*k2[2], y3[i]+h*k3[2], y4[i]+h*k4[2], y5[i]+h*k5[2]);
      k2[3] = f2(x[i+1], y1[i]+h*k1[2], y2[i]+h*k2[2], y3[i]+h*k3[2], y4[i]+h*k4[2], y5[i]+h*k5[2]);
      k3[3] = f2(x[i+1], y1[i]+h*k1[2], y2[i]+h*k2[2], y3[i]+h*k3[2], y4[i]+h*k4[2], y5[i]+h*k5[2]);
      k4[3] = f4(x[i+1], y1[i]+h*k1[2], y2[i]+h*k2[2], y3[i]+h*k3[2], y4[i]+h*k4[2], y5[i]+h*k5[2]);
      k5[3] = f5(x[i+1], y1[i]+h*k1[2], y2[i]+h*k2[2], y3[i]+h*k3[2], y4[i]+h*k4[2], y5[i]+h*k5[2]);

      // Dependent variables

      y1[i+1] = y1[i] + h*(k1[0] + 2.0*k1[1] + 2.0*k1[2] + k1[3])/6.0;
      y2[i+1] = y2[i] + h*(k2[0] + 2.0*k2[1] + 2.0*k2[2] + k2[3])/6.0;
      y3[i+1] = y3[i] + h*(k3[0] + 2.0*k3[1] + 2.0*k3[2] + k3[3])/6.0;
      y4[i+1] = y4[i] + h*(k4[0] + 2.0*k4[1] + 2.0*k4[2] + k4[3])/6.0;
      y5[i+1] = y5[i] + h*(k5[0] + 2.0*k5[1] + 2.0*k5[2] + k5[3])/6.0;
   }
}


// int calculateRungeKuttaFehlbergIntegral(Vector<double>& x, Vector<double>& y, 
// double (*f)(double, double), 
// double a, double b, double ya) method

/// This method integrates the differential equation
///
/// <ul>
/// <li>  y' = f(x, y),
/// </ul>
///
/// where y' is given as a C-style function, from x=a to x=b with initial
/// condition y(a)=ya. The method is based on an explicit fourth order 
/// Runge-Kutta-Fehlberg formula.  
///
/// @param xOut Pointer to vector with x values. 
/// @param yOut Pointer to vector with y values.
/// @param f Pointer to C-style function f (state equation for variable y).
/// @param a Lower integration limit.
/// @param b Upper integration limit.
/// @param ya Initial condition for variable y.

int OrdinaryDifferentialEquations
::calculateRungeKuttaFehlbergIntegral(Vector<double>& xOut, Vector<double>& yOut, 
double (*f)(double, double), double a, double b, double ya)
{
   double error = 0.0;

   int count = 0;
   int size = initialSize;

   Vector<double> tempX(size);
   Vector<double> tempY(size);

   tempX[0] = a;
   tempY[0] = ya;

   count++;

   double epsilon = std::numeric_limits<double>::epsilon();

   double a2 = 1.0/5.0;
   double a3 = 3.0/10.0;
   double a4 = 3.0/5.0;
   double a5 = 1.0;
   double a6 = 7.0/8.0;

   double b21 = 1.0/5.0;   
   double b31 = 3.0/40.0;
   double b32 = 9.0/40.0;   
   double b41 = 3.0/10.0;
   double b42 = -9.0/10.0;   
   double b43 = 6.0/5.0; 
   double b51 = -11.0/54.0;
   double b52 = 5.0/2.0;
   double b53 = -70.0/27.0;
   double b54 = 35.0/27.0; 
   double b61 = 1631.0/55296.0;
   double b62 = 175.0/512.0;   
   double b63 = 575.0/13824.0;   
   double b64 = 44275.0/110592.0;   
   double b65 = 253.0/4096.0;
   
   double c41 = 37.0/378.0;
   double c42 = 0.0;
   double c43 = 250.0/621.0;
   double c44 = 125.0/594.0;
   double c45 = 0.0;
   double c46 = 512.0/1771.0;

   double c51 = 2825.0/27648.0;  
   double c52 = 0.0;  
   double c53 = 18575.0/48384.0;   
   double c54 = 13525.0/55296.0;
   double c55 = 277.0/14336.0;
   double c56 = 1.0/4.0;

   double d1 = c41 - c51;
   double d2 = c42 - c52;
   double d3 = c43 - c53;
   double d4 = c44 - c54;
   double d5 = c45 - c55;
   double d6 = c46 - c56;

   double k1 = 0.0;
   double k2 = 0.0;
   double k3 = 0.0;
   double k4 = 0.0;
   double k5 = 0.0;
   double k6 = 0.0;

   // Main loop

   double x = a;
   double y = ya;

   double hmin = 16.0*epsilon*fabs(x);

   double h = (b-a)*1.0e-3;

   while(x < b)
   { 
      // Set smallest allowable stepsize

      hmin = 32.0*epsilon*fabs(x);

      if(h < hmin)
      {
         std::cout << std::endl
                   << "Flood Warning: OrdinaryDifferentialEquations class." << std::endl
                   << "calculateRungeKuttaFehlbergIntegral() method." << std::endl
                   << "Step size is less than smallest allowable." << std::endl
                   << std::endl;

         h = hmin;
      }

      if(x + h > b)
      {
         h = b - x;
      }

      // Runge-Kutta-Fehlberg coefficients

      k1 = h*f(x, y);
      k2 = h*f(x+a2*h, y+b21*k1);
      k3 = h*f(x+a3*h, y+b31*k1+b32*k2);
      k4 = h*f(x+a4*h, y+b41*k1+b42*k2+b43*k3); 
      k5 = h*f(x+a5*h, y+b51*k1+b52*k2+b53*k3+b54*k4); 
      k6 = h*f(x+a6*h, y+b61*k1+b62*k2+b63*k3+b64*k4+b65*k5);

      // Error estimate

      error = fabs(d1*k1+d2*k2+d3*k3+d4*k4+d5*k5+d6*k6);

      if(error <= tolerance)
      {
         x += h;
         y += c51*k1+c52*k2+c53*k3+c54*k4+c55*k5+c56*k6;

         tempX[count] = x;
         tempY[count] = y;

         count++;

         if(error != 0.0)
         {
            h *= 0.9*pow(fabs(tolerance/error),0.2);      
         }

         if(count >= size)
         {
            std::cerr << std::endl
                      << "Flood Error: calculateRungeKuttaFehlbergIntegral() method." << std::endl
                      << "Insufficient memory reserved." << std::endl
                      << std::endl;

            exit(1);
         }
      }
      else
      {
         h *= 0.9*pow(fabs(tolerance/error),0.25);
      }
   }// end while loop   

   Vector<double> newX(count);
   Vector<double> newY(count);

   for(int i = 0; i < count; i++)
   {
      newX[i] = tempX[i];
      newY[i] = tempY[i];
   }

   xOut = newX;
   yOut = newY;

   return(count);
}


// int calculateRungeKuttaFehlbergIntegral(Vector<double>&, Vector<double>&, Vector<double>&, 
// double (*f1)(double, double, double), 
// double (*f2)(double, double, double), 
// double, double, double, double)

/// This method integrates the system of two differential equations 
///
/// <ul>
/// <li> y1' = f1(x, y1, y2),
/// <li> y2' = f2(x, y1, y2),
/// </ul>
///
/// where y1' and y2' are given as C-syle functions, from
/// x=a to x=b with initial conditions y1(a)=y1a and y2(a)=y2a.
/// The method is based on an explicit fourth order 
/// Runge-Kutta-Fehlberg formula.  
///
/// @param xOut Pointer to vector with x values. 
/// @param y1Out Pointer to vector with y1 values.
/// @param y2Out Pointer to vector with y2 values.
/// @param f1 Member method to integrate (state equation for variable y1).
/// @param f2 Member method to integrate (state equation for variable y2).
/// @param a Lower integration limit. 
/// @param b: Upper integration limit. 
/// @param y1a Initial condition for variable y1.
/// @param y2a Initial condition for variable y2.

int OrdinaryDifferentialEquations::calculateRungeKuttaFehlbergIntegral(
Vector<double>& xOut, Vector<double>& y1Out, Vector<double>& y2Out, 
double (*f1)(double, double, double), 
double (*f2)(double, double, double), 
double a, double b, double y1a, double y2a)
{
   double error = 0.0, error1 = 0.0, error2 = 0.0;

   int count = 0;
   int size = initialSize;

   Vector<double> tempX(size);
   Vector<double> tempY1(size);
   Vector<double> tempY2(size);

   tempX[0] = a;

   tempY1[0] = y1a;
   tempY2[0] = y2a;

   count++;

   double epsilon = std::numeric_limits<double>::epsilon();

   double a2 = 1.0/5.0;
   double a3 = 3.0/10.0;
   double a4 = 3.0/5.0;
   double a5 = 1.0;
   double a6 = 7.0/8.0;

   double b21 = 1.0/5.0;   
   double b31 = 3.0/40.0;
   double b32 = 9.0/40.0;   
   double b41 = 3.0/10.0;
   double b42 = -9.0/10.0;   
   double b43 = 6.0/5.0; 
   double b51 = -11.0/54.0;
   double b52 = 5.0/2.0;
   double b53 = -70.0/27.0;
   double b54 = 35.0/27.0; 
   double b61 = 1631.0/55296.0;
   double b62 = 175.0/512.0;   
   double b63 = 575.0/13824.0;   
   double b64 = 44275.0/110592.0;   
   double b65 = 253.0/4096.0;

   double c41 = 37.0/378.0;
   double c42 = 0.0;
   double c43 = 250.0/621.0;
   double c44 = 125.0/594.0;
   double c45 = 0.0;
   double c46 = 512.0/1771.0;

   double c51 = 2825.0/27648.0;  
   double c52 = 0.0;  
   double c53 = 18575.0/48384.0;   
   double c54 = 13525.0/55296.0;
   double c55 = 277.0/14336.0;
   double c56 = 1.0/4.0;

   double d1 = c41 - c51;
   double d2 = c42 - c52;
   double d3 = c43 - c53;
   double d4 = c44 - c54;
   double d5 = c45 - c55;
   double d6 = c46 - c56;

   double k11 = 0.0, k21 = 0.0;
   double k12 = 0.0, k22 = 0.0;
   double k13 = 0.0, k23 = 0.0;
   double k14 = 0.0, k24 = 0.0;
   double k15 = 0.0, k25 = 0.0;
   double k16 = 0.0, k26 = 0.0;

   // Main loop

   double x = a;
   double y1 = y1a;
   double y2 = y2a;

   double hmin = 16.0*epsilon*fabs(x);

   double h = (b-a)*1.0e-3;

   while(x < b)
   {
      // Set smallest allowable stepsize

      hmin = 32.0*epsilon*fabs(x);

      if(h < hmin)
      {
         std::cout << std::endl
                   << "Flood Warning: OrdinaryDifferentialEquations class." << std::endl
                   << "calculateRungeKuttaFehlbergIntegral() method." << std::endl
                   << "Step size is less than smallest allowable." << std::endl
                   << std::endl;

       h = hmin;
      }

      if(x + h > b)
      {
         h = b - x;
      }

      // Runge-Kutta-Fehlberg coefficients

      k11 = h*f1(x, y1, y2);
      k21 = h*f2(x, y1, y2);

      k12 = h*f1(x+a2*h, y1+b21*k11, y2+b21*k21);
      k22 = h*f2(x+a2*h, y1+b21*k11, y2+b21*k21);

      k13 = h*f1(x+a3*h, y1+b31*k11+b32*k12, y2+b31*k21+b32*k22);
      k23 = h*f2(x+a3*h, y1+b31*k11+b32*k12, y2+b31*k21+b32*k22);

      k14 = h*f1(x+a4*h, y1+b41*k11+b42*k12+b43*k13, y2+b41*k21+b42*k22+b43*k23); 
      k24 = h*f2(x+a4*h, y1+b41*k11+b42*k12+b43*k13, y2+b41*k21+b42*k22+b43*k23); 

      k15 = h*f1(x+a5*h, 
      y1+b51*k11+b52*k12+b53*k13+b54*k14, 
      y2+b51*k21+b52*k22+b53*k23+b54*k24); 

      k25 = h*f2(x+a5*h, 
      y1+b51*k11+b52*k12+b53*k13+b54*k14, 
      y2+b51*k21+b52*k22+b53*k23+b54*k24); 

      k16 = h*f1(x+a6*h, 
      y1+b61*k11+b62*k12+b63*k13+b64*k14+b65*k15, 
      y2+b61*k21+b62*k22+b63*k23+b64*k24+b65*k25);

      k26 = h*f2(x+a6*h, 
      y1+b61*k11+b62*k12+b63*k13+b64*k14+b65*k15, 
      y2+b61*k21+b62*k22+b63*k23+b64*k24+b65*k25);

      // Error estimate

      error1 = fabs(d1*k11+d2*k12+d3*k13+d4*k14+d5*k15+d6*k16);
      error2 = fabs(d1*k21+d2*k22+d3*k23+d4*k24+d5*k25+d6*k26);

      error = 0.0;
      if(error1 > error)
      {
         error = error1;
      }
      if(error2 > error)
      {
         error = error2;
      }

      if(error <= tolerance)
      {
         x += h;
         y1 += c51*k11+c52*k12+c53*k13+c54*k14+c55*k15+c56*k16;
         y2 += c51*k21+c52*k22+c53*k23+c54*k24+c55*k25+c56*k26;

         tempX[count] = x;
         tempY1[count] = y1;
         tempY2[count] = y2;
         count++;

         if(error != 0.0)
         {
            h *= 0.9*pow(fabs(tolerance/error),0.2);
         }

         if(count >= size)
         {
            std::cerr << std::endl
                      << "Flood Error: calculateRungeKuttaFehlbergIntegral() method." << std::endl
                      << "Insufficient memory reserved." << std::endl
                      << std::endl;

            exit(1);
         }
      }
      else
      {
         h *= 0.9*pow(fabs(tolerance/error),0.25);
      }
   } // end while loop   

   Vector<double> newX(count);
   Vector<double> newY1(count);
   Vector<double> newY2(count);

   for(int i = 0; i < count; i++)
   {
      newX[i] = tempX[i];
      newY1[i] = tempY1[i];
      newY2[i] = tempY2[i];
   }

   xOut = newX;
   y1Out = newY1;
   y2Out = newY2;

   return(count);
}


// int calculateRungeKuttaFehlbergIntegral(Vector<double>&,
// Vector<double>&, Vector<double>&, Vector<double>&, 
// double (*f1)(double, double, double, double), 
// double (*f2)(double, double, double, double), 
// double (*f3)(double, double, double, double), 
// double, double, double, double, double) method

/// This method integrates the system of three differential equations 
///
/// <ul>
/// <li> y1' = f1(x, y1, y2, y3),
/// <li> y2' = f2(x, y1, y2, y3),
/// <li> y3' = f3(x, y1, y2, y3),
/// </ul>
///
/// where y1', y2' and y3' are given as C-style functions, from
/// x=a to x=b with initial conditions y1(a)=y1a, y2(a)=y2a and y3(a)=y3a.
/// The method is based on an explicit fourth order 
/// Runge-Kutta-Fehlberg formula.  
///
/// @param xOut Pointer to vector with x values. 
/// @param y1Out Pointer to vector with y1 values.
/// @param y2Out Pointer to vector with y2 values.
/// @param y3Out Pointer to vector with y3 values.
/// @param f1 Member method to integrate (state equation for variable y1).
/// @param f2 Member method to integrate (state equation for variable y2).
/// @param f3 Member method to integrate (state equation for variable y3).
/// @param a Lower integration limit. 
/// @param b Upper integration limit. 
/// @param y1a Initial condition for variable y1.
/// @param y2a Initial condition for variable y2.
/// @param y3a Initial condition for variable y3.

int OrdinaryDifferentialEquations::calculateRungeKuttaFehlbergIntegral(
Vector<double>& xOut, Vector<double>& y1Out, Vector<double>& y2Out, Vector<double>& y3Out,
double (*f1)(double, double, double, double), 
double (*f2)(double, double, double, double), 
double (*f3)(double, double, double, double), 
double a, double b, double y1a, double y2a, double y3a)
{
   double error = 0.0, error1 = 0.0, error2 = 0.0, error3 = 0.0;

   int size = initialSize;

   int count = 0;

   Vector<double> tempX(size);
   Vector<double> tempY1(size);
   Vector<double> tempY2(size);
   Vector<double> tempY3(size);

   tempX[0] = a;

   tempY1[0] = y1a;
   tempY2[0] = y2a;
   tempY3[0] = y3a;

   count++;

   double epsilon = std::numeric_limits<double>::epsilon();

   double a2 = 1.0/5.0;
   double a3 = 3.0/10.0;
   double a4 = 3.0/5.0;
   double a5 = 1.0;
   double a6 = 7.0/8.0;

   double b21 = 1.0/5.0;   
   double b31 = 3.0/40.0;
   double b32 = 9.0/40.0;   
   double b41 = 3.0/10.0;
   double b42 = -9.0/10.0;   
   double b43 = 6.0/5.0; 
   double b51 = -11.0/54.0;
   double b52 = 5.0/2.0;
   double b53 = -70.0/27.0;
   double b54 = 35.0/27.0; 
   double b61 = 1631.0/55296.0;
   double b62 = 175.0/512.0;   
   double b63 = 575.0/13824.0;   
   double b64 = 44275.0/110592.0;   
   double b65 = 253.0/4096.0;
   
   double c41 = 37.0/378.0;
   double c42 = 0.0;
   double c43 = 250.0/621.0;
   double c44 = 125.0/594.0;
   double c45 = 0.0;
   double c46 = 512.0/1771.0;

   double c51 = 2825.0/27648.0;  
   double c52 = 0.0;  
   double c53 = 18575.0/48384.0;   
   double c54 = 13525.0/55296.0;
   double c55 = 277.0/14336.0;
   double c56 = 1.0/4.0;

   double d1 = c41 - c51;
   double d2 = c42 - c52;
   double d3 = c43 - c53;
   double d4 = c44 - c54;
   double d5 = c45 - c55;
   double d6 = c46 - c56;

   double k11 = 0.0, k21 = 0.0, k31 = 0.0;
   double k12 = 0.0, k22 = 0.0, k32 = 0.0;
   double k13 = 0.0, k23 = 0.0, k33 = 0.0;
   double k14 = 0.0, k24 = 0.0, k34 = 0.0;
   double k15 = 0.0, k25 = 0.0, k35 = 0.0;
   double k16 = 0.0, k26 = 0.0, k36 = 0.0;

   // Main loop

   double x = a;
   double y1 = y1a;
   double y2 = y2a;
   double y3 = y3a;

   double hmin = 16.0*epsilon*fabs(x);

   double h = (b-a)*1.0e-3;

   while(x < b)
   {   
      // Set smallest allowable stepsize

      hmin = 32.0*epsilon*fabs(x);

      if(h < hmin)
      {
         std::cout << std::endl
                   << "Flood Warning: OrdinaryDifferentialEquations class." << std::endl
                   << "calculateRungeKuttaFehlbergIntegral() method." << std::endl
                   << "Step size is less than smallest allowable." << std::endl
                   << std::endl;

       h = hmin;
      }

      if(x + h > b)
      {
         h = b - x;
      }

      // Runge-Kutta-Fehlberg coefficients

      k11 = h*f1(x, y1, y2, y3);
      k21 = h*f2(x, y1, y2, y3);
      k31 = h*f3(x, y1, y2, y3);

      k12 = h*f1(x+a2*h, y1+b21*k11, y2+b21*k21, y3+b21*k31);
      k22 = h*f2(x+a2*h, y1+b21*k11, y2+b21*k21, y3+b21*k31);
      k32 = h*f3(x+a2*h, y1+b21*k11, y2+b21*k21, y3+b21*k31);

      k13 = h*f1(x+a3*h, y1+b31*k11+b32*k12, y2+b31*k21+b32*k22, y3+b31*k31+b32*k32);
      k23 = h*f2(x+a3*h, y1+b31*k11+b32*k12, y2+b31*k21+b32*k22, y3+b31*k31+b32*k32);
      k33 = h*f3(x+a3*h, y1+b31*k11+b32*k12, y2+b31*k21+b32*k22, y3+b31*k31+b32*k32);

      k14 = h*f1(x+a4*h, y1+b41*k11+b42*k12+b43*k13, y2+b41*k21+b42*k22+b43*k23, y3+b41*k31+b42*k32+b43*k33); 
      k24 = h*f2(x+a4*h, y1+b41*k11+b42*k12+b43*k13, y2+b41*k21+b42*k22+b43*k23, y3+b41*k31+b42*k32+b43*k33); 
      k34 = h*f3(x+a4*h, y1+b41*k11+b42*k12+b43*k13, y2+b41*k21+b42*k22+b43*k23, y3+b41*k31+b42*k32+b43*k33); 

      k15 = h*f1(x+a5*h, y1+b51*k11+b52*k12+b53*k13+b54*k14, y2+b51*k21+b52*k22+b53*k23+b54*k24, y3+b51*k31+b52*k32+b53*k33+b54*k34); 
      k25 = h*f2(x+a5*h, y1+b51*k11+b52*k12+b53*k13+b54*k14, y2+b51*k21+b52*k22+b53*k23+b54*k24, y3+b51*k31+b52*k32+b53*k33+b54*k34); 
      k35 = h*f3(x+a5*h, y1+b51*k11+b52*k12+b53*k13+b54*k14, y2+b51*k21+b52*k22+b53*k23+b54*k24, y3+b51*k31+b52*k32+b53*k33+b54*k34); 

      k16 = h*f1(x+a6*h, y1+b61*k11+b62*k12+b63*k13+b64*k14+b65*k15, y2+b61*k21+b62*k22+b63*k23+b64*k24+b65*k25, y3+b61*k31+b62*k32+b63*k33+b64*k34+b65*k35);
      k26 = h*f2(x+a6*h, y1+b61*k11+b62*k12+b63*k13+b64*k14+b65*k15, y2+b61*k21+b62*k22+b63*k23+b64*k24+b65*k25, y3+b61*k31+b62*k32+b63*k33+b64*k34+b65*k35);
      k36 = h*f3(x+a6*h, y1+b61*k11+b62*k12+b63*k13+b64*k14+b65*k15, y2+b61*k21+b62*k22+b63*k23+b64*k24+b65*k25, y3+b61*k31+b62*k32+b63*k33+b64*k34+b65*k35);

      // Error estimate

      error1 = fabs(d1*k11+d2*k12+d3*k13+d4*k14+d5*k15+d6*k16);
      error2 = fabs(d1*k21+d2*k22+d3*k23+d4*k24+d5*k25+d6*k26);
      error3 = fabs(d1*k31+d2*k32+d3*k33+d4*k34+d5*k35+d6*k36);

      error = 0.0;
      if(error1 > error)
      {
         error = error1;
      }
      if(error2 > error)
      {
         error = error2;
      }
      if(error3 > error)
      {
         error = error3;
      }
      if(error <= tolerance)
      {
         x += h;

         y1 += c51*k11+c52*k12+c53*k13+c54*k14+c55*k15+c56*k16;
         y2 += c51*k21+c52*k22+c53*k23+c54*k24+c55*k25+c56*k26;
         y3 += c51*k31+c52*k32+c53*k33+c54*k34+c55*k35+c56*k36;

         tempX[count] = x;
         tempY1[count] = y1;
         tempY2[count] = y2;
         tempY3[count] = y3;

         count++;

         if(error != 0.0)
         {
            h *= 0.9*pow(fabs(tolerance/error),0.2);
         }

         if(count >= size)
         {
            std::cout << std::endl
                      << "Flood Error calculateRungeKuttaFehlbergIntegral() method." << std::endl
                      << "Insufficient memory reserved." << std::endl
                      << std::endl;

            exit(1);
         }
      }
      else
      {
         h *= 0.9*pow(fabs(tolerance/error),0.25);      
      }
   } // end while loop   

   Vector<double> newX(count);
   Vector<double> newY1(count);
   Vector<double> newY2(count);
   Vector<double> newY3(count);

   for(int i = 0; i < count; i++)
   {
      newX[i] = tempX[i];
      newY1[i] = tempY1[i];
      newY2[i] = tempY2[i];
      newY3[i] = tempY3[i];
   }

   xOut = newX;
   y1Out = newY1;
   y2Out = newY2;
   y3Out = newY3;

   return(count);
}


// int calculateRungeKuttaFehlbergIntegral(Vector<double>&,
// Vector<double>&, Vector<double>&, Vector<double>&, Vector<double>&,
// double (*f1)(double, double, double, double, double), 
// double (*f2)(double, double, double, double, double), 
// double (*f3)(double, double, double, double, double), 
// double (*f4)(double, double, double, double, double), 
// double, double, double, double, double, double) method

/// This method integrates the system of four differential equations 
///
/// <ul>
/// <li> y1' = f1(x, y1, y2, y3, y4),
/// <li> y2' = f2(x, y1, y2, y3, y4),
/// <li> y3' = f3(x, y1, y2, y3, y4),
/// <li> y4' = f4(x, y1, y2, y3, y4),
/// </ul>
///
/// where y1', y2', y3' and y4' are given as C-style functions, from
/// x=a to x=b with initial conditions y1(a)=y1a, y2(a)=y2a, y3(a)=y3a and y4(a)=y4a.
/// The method is based on an explicit fourth order 
/// Runge-Kutta-Fehlberg formula.  
///
/// @param xOut Pointer to vector with x values. 
/// @param y1Out Pointer to vector with y1 values.
/// @param y2Out Pointer to vector with y2 values.
/// @param y3Out Pointer to vector with y3 values.
/// @param y4Out Pointer to vector with y4 values.
/// @param f1 Member method to integrate (state equation for variable y1).
/// @param f2 Member method to integrate (state equation for variable y2).
/// @param f3 Member method to integrate (state equation for variable y3).
/// @param f4 Member method to integrate (state equation for variable y4).
/// @param a Lower integration limit. 
/// @param b Upper integration limit. 
/// @param y1a Initial condition for variable y1.
/// @param y2a Initial condition for variable y2.
/// @param y3a Initial condition for variable y3.
/// @param y4a Initial condition for variable y4.

int OrdinaryDifferentialEquations::calculateRungeKuttaFehlbergIntegral(
Vector<double>& xOut, Vector<double>& y1Out, Vector<double>& y2Out, Vector<double>& y3Out, Vector<double>& y4Out,
double (*f1)(double, double, double, double, double), 
double (*f2)(double, double, double, double, double), 
double (*f3)(double, double, double, double, double), 
double (*f4)(double, double, double, double, double), 
double a, double b, double y1a, double y2a, double y3a, double y4a)
{
   double error = 0.0, error1 = 0.0, error2 = 0.0, error3 = 0.0, error4 = 0.0;

   int size = initialSize;

   int count = 0;

   Vector<double> tempX(size);
   Vector<double> tempY1(size);
   Vector<double> tempY2(size);
   Vector<double> tempY3(size);
   Vector<double> tempY4(size);

   tempX[0] = a;

   tempY1[0] = y1a;
   tempY2[0] = y2a;
   tempY3[0] = y3a;
   tempY4[0] = y4a;

   count++;

   double epsilon = std::numeric_limits<double>::epsilon();

   double a2 = 1.0/5.0;
   double a3 = 3.0/10.0;
   double a4 = 3.0/5.0;
   double a5 = 1.0;
   double a6 = 7.0/8.0;

   double b21 = 1.0/5.0;   
   double b31 = 3.0/40.0;
   double b32 = 9.0/40.0;   
   double b41 = 3.0/10.0;
   double b42 = -9.0/10.0;   
   double b43 = 6.0/5.0; 
   double b51 = -11.0/54.0;
   double b52 = 5.0/2.0;
   double b53 = -70.0/27.0;
   double b54 = 35.0/27.0; 
   double b61 = 1631.0/55296.0;
   double b62 = 175.0/512.0;   
   double b63 = 575.0/13824.0;   
   double b64 = 44275.0/110592.0;   
   double b65 = 253.0/4096.0;
   
   double c41 = 37.0/378.0;
   double c42 = 0.0;
   double c43 = 250.0/621.0;
   double c44 = 125.0/594.0;
   double c45 = 0.0;
   double c46 = 512.0/1771.0;

   double c51 = 2825.0/27648.0;  
   double c52 = 0.0;  
   double c53 = 18575.0/48384.0;   
   double c54 = 13525.0/55296.0;
   double c55 = 277.0/14336.0;
   double c56 = 1.0/4.0;

   double d1 = c41 - c51;
   double d2 = c42 - c52;
   double d3 = c43 - c53;
   double d4 = c44 - c54;
   double d5 = c45 - c55;
   double d6 = c46 - c56;

   double k11 = 0.0, k21 = 0.0, k31 = 0.0, k41 = 0.0;
   double k12 = 0.0, k22 = 0.0, k32 = 0.0, k42 = 0.0;
   double k13 = 0.0, k23 = 0.0, k33 = 0.0, k43 = 0.0;
   double k14 = 0.0, k24 = 0.0, k34 = 0.0, k44 = 0.0;
   double k15 = 0.0, k25 = 0.0, k35 = 0.0, k45 = 0.0;
   double k16 = 0.0, k26 = 0.0, k36 = 0.0, k46 = 0.0;

   // Main loop

   double x = a;
   double y1 = y1a;
   double y2 = y2a;
   double y3 = y3a;
   double y4 = y4a;

   double hmin = 16.0*epsilon*fabs(x);

   double h = (b-a)*1.0e-3;

   while(x < b)
   {   
      // Set smallest allowable stepsize

      hmin = 32.0*epsilon*fabs(x);

      if(h < hmin)
      {
         std::cout << std::endl
                   << "Flood Warning: OrdinaryDifferentialEquations class." << std::endl
                   << "calculateRungeKuttaFehlbergIntegral() method." << std::endl
                   << "Step size is less than smallest allowable." << std::endl
                   << std::endl;

       h = hmin;
      }

      if(x + h > b)
      {
         h = b - x;
      }

      // Runge-Kutta-Fehlberg coefficients

      k11 = h*f1(x, y1, y2, y3, y4);
      k21 = h*f2(x, y1, y2, y3, y4);
      k31 = h*f3(x, y1, y2, y3, y4);
      k41 = h*f4(x, y1, y2, y3, y4);

      k12 = h*f1(x+a2*h, y1+b21*k11, y2+b21*k21, y3+b21*k31, y4+b21*k41);
      k22 = h*f2(x+a2*h, y1+b21*k11, y2+b21*k21, y3+b21*k31, y4+b21*k41);
      k32 = h*f3(x+a2*h, y1+b21*k11, y2+b21*k21, y3+b21*k31, y4+b21*k41);
      k42 = h*f4(x+a2*h, y1+b21*k11, y2+b21*k21, y3+b21*k31, y4+b21*k41);

      k13 = h*f1(x+a3*h, y1+b31*k11+b32*k12, y2+b31*k21+b32*k22, y3+b31*k31+b32*k32, y4+b31*k41+b32*k42);
      k23 = h*f2(x+a3*h, y1+b31*k11+b32*k12, y2+b31*k21+b32*k22, y3+b31*k31+b32*k32, y4+b31*k41+b32*k42);
      k33 = h*f3(x+a3*h, y1+b31*k11+b32*k12, y2+b31*k21+b32*k22, y3+b31*k31+b32*k32, y4+b31*k41+b32*k42);
      k43 = h*f4(x+a3*h, y1+b31*k11+b32*k12, y2+b31*k21+b32*k22, y3+b31*k31+b32*k32, y4+b31*k41+b32*k42);

      k14 = h*f1(x+a4*h, y1+b41*k11+b42*k12+b43*k13, y2+b41*k21+b42*k22+b43*k23, y3+b41*k31+b42*k32+b43*k33, y4+b41*k41+b42*k42+b43*k43); 
      k24 = h*f2(x+a4*h, y1+b41*k11+b42*k12+b43*k13, y2+b41*k21+b42*k22+b43*k23, y3+b41*k31+b42*k32+b43*k33, y4+b41*k41+b42*k42+b43*k43); 
      k34 = h*f3(x+a4*h, y1+b41*k11+b42*k12+b43*k13, y2+b41*k21+b42*k22+b43*k23, y3+b41*k31+b42*k32+b43*k33, y4+b41*k41+b42*k42+b43*k43); 
      k44 = h*f4(x+a4*h, y1+b41*k11+b42*k12+b43*k13, y2+b41*k21+b42*k22+b43*k23, y3+b41*k31+b42*k32+b43*k33, y4+b41*k41+b42*k42+b43*k43); 

      k15 = h*f1(x+a5*h, y1+b51*k11+b52*k12+b53*k13+b54*k14, y2+b51*k21+b52*k22+b53*k23+b54*k24, y3+b51*k31+b52*k32+b53*k33+b54*k34, y4+b51*k41+b52*k42+b53*k43+b54*k44); 
      k25 = h*f2(x+a5*h, y1+b51*k11+b52*k12+b53*k13+b54*k14, y2+b51*k21+b52*k22+b53*k23+b54*k24, y3+b51*k31+b52*k32+b53*k33+b54*k34, y4+b51*k41+b52*k42+b53*k43+b54*k44); 
      k35 = h*f3(x+a5*h, y1+b51*k11+b52*k12+b53*k13+b54*k14, y2+b51*k21+b52*k22+b53*k23+b54*k24, y3+b51*k31+b52*k32+b53*k33+b54*k34, y4+b51*k41+b52*k42+b53*k43+b54*k44); 
      k45 = h*f4(x+a5*h, y1+b51*k11+b52*k12+b53*k13+b54*k14, y2+b51*k21+b52*k22+b53*k23+b54*k24, y3+b51*k31+b52*k32+b53*k33+b54*k34, y4+b51*k41+b52*k42+b53*k43+b54*k44); 

      k16 = h*f1(x+a6*h, y1+b61*k11+b62*k12+b63*k13+b64*k14+b65*k15, y2+b61*k21+b62*k22+b63*k23+b64*k24+b65*k25, y3+b61*k31+b62*k32+b63*k33+b64*k34+b65*k35, y4+b61*k41+b62*k42+b63*k43+b64*k44+b65*k45);
      k26 = h*f2(x+a6*h, y1+b61*k11+b62*k12+b63*k13+b64*k14+b65*k15, y2+b61*k21+b62*k22+b63*k23+b64*k24+b65*k25, y3+b61*k31+b62*k32+b63*k33+b64*k34+b65*k35, y4+b61*k41+b62*k42+b63*k43+b64*k44+b65*k45);
      k36 = h*f3(x+a6*h, y1+b61*k11+b62*k12+b63*k13+b64*k14+b65*k15, y2+b61*k21+b62*k22+b63*k23+b64*k24+b65*k25, y3+b61*k31+b62*k32+b63*k33+b64*k34+b65*k35, y4+b61*k41+b62*k42+b63*k43+b64*k44+b65*k45);
      k46 = h*f4(x+a6*h, y1+b61*k11+b62*k12+b63*k13+b64*k14+b65*k15, y2+b61*k21+b62*k22+b63*k23+b64*k24+b65*k25, y3+b61*k31+b62*k32+b63*k33+b64*k34+b65*k35, y4+b61*k41+b62*k42+b63*k43+b64*k44+b65*k45);

      // Error estimate

      error1 = fabs(d1*k11+d2*k12+d3*k13+d4*k14+d5*k15+d6*k16);
      error2 = fabs(d1*k21+d2*k22+d3*k23+d4*k24+d5*k25+d6*k26);
      error3 = fabs(d1*k31+d2*k32+d3*k33+d4*k34+d5*k35+d6*k36);
      error4 = fabs(d1*k41+d2*k42+d3*k43+d4*k44+d5*k45+d6*k46);

      error = 0.0;
      if(error1 > error)
      {
         error = error1;
      }
      if(error2 > error)
      {
         error = error2;
      }
      if(error3 > error)
      {
         error = error3;
      }
      if(error4 > error)
      {
         error = error4;
      }
      if(error <= tolerance)
      {
         x += h;

         y1 += c51*k11+c52*k12+c53*k13+c54*k14+c55*k15+c56*k16;
         y2 += c51*k21+c52*k22+c53*k23+c54*k24+c55*k25+c56*k26;
         y3 += c51*k31+c52*k32+c53*k33+c54*k34+c55*k35+c56*k36;
         y4 += c51*k41+c52*k42+c53*k43+c54*k44+c55*k45+c56*k46;

         tempX[count] = x;
         tempY1[count] = y1;
         tempY2[count] = y2;
         tempY3[count] = y3;
         tempY4[count] = y4;

         count++;

         if(error != 0.0)
         {
            h *= 0.9*pow(fabs(tolerance/error),0.2);
         }

         if(count >= size)
         {
            std::cout << std::endl
                      << "Flood Error calculateRungeKuttaFehlbergIntegral() method." << std::endl
                      << "Insufficient memory reserved." << std::endl
                      << std::endl;

            exit(1);
         }
      }
      else
      {
         h *= 0.9*pow(fabs(tolerance/error),0.25);      
      }
   } // end while loop   

   Vector<double> newX(count);
   Vector<double> newY1(count);
   Vector<double> newY2(count);
   Vector<double> newY3(count);
   Vector<double> newY4(count);

   for(int i = 0; i < count; i++)
   {
      newX[i] = tempX[i];
      newY1[i] = tempY1[i];
      newY2[i] = tempY2[i];
      newY3[i] = tempY3[i];
      newY4[i] = tempY4[i];
   }

   xOut = newX;
   y1Out = newY1;
   y2Out = newY2;
   y3Out = newY3;
   y4Out = newY4;

   return(count);
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
