/****************************************************************************************************************/
/*                                                                                                              */
/*   Flood: An Open Source Neural Networks C++ Library                                                          */
/*   www.cimne.com/flood                                                                                        */
/*                                                                                                              */
/*   V E C T O R   C O N T A I N E R                                                                            */
/*                                                                                                              */
/*   Roberto Lopez                                                                                              */
/*   International Center for Numerical Methods in Engineering (CIMNE)                                          */
/*   Technical University of Catalonia (UPC)                                                                    */
/*   Barcelona, Spain                                                                                           */
/*   E-mail: rlopez@cimne.upc.edu                                                                               */
/*                                                                                                              */
/****************************************************************************************************************/

#ifndef __VECTOR_H__
#define __VECTOR_H__

// System includes

#include <iostream>
#include <fstream>
#include <math.h>
#include <stdlib.h>


namespace Flood
{

/// Forward declaration of Matrix template

template <class Type> class Matrix;

/// This template class defines a vector for general purpose use.
///
/// @see Matrix.

template <class Type>
class Vector
{

private:

   /// Size of vector.

   int size;

   /// Pointer to a Type.

   Type* vector;

public:


   // CONSTRUCTORS

   Vector(void);

   Vector(int);

   Vector(int, const Type&);

   Vector(int, const Type*);

   Vector(const Vector&);


   // ASSINGMENT OPERATOR

   Vector& operator=(const Vector&);


   // REFERENCE OPERATORS

   inline Type& operator[](const int);

   inline const Type& operator[](const int) const;

   // METHODS

   inline int getSize(void);
   inline void setSize(int);

   inline void resize(int);

   inline void fillAtRandom(void);
   inline void fillAtRandom(double, double);

   inline Type calculateMean(void);
   inline Type calculateStandardDeviation(void);

   inline Type calculateMinimum(void);
   inline Type calculateMaximum(void);

   inline int calculateMinimalIndex(void);
   inline int calculateMaximalIndex(void);

   inline double calculateNorm(void);

   inline Vector<Type> operator+(Type);
   inline Vector<Type> operator+(Vector<Type>);

   inline Vector<Type> operator-(Type);
   inline Vector<Type> operator-(Vector<Type>);

   inline Vector<Type> operator*(Type);
   inline Type dot(Vector<Type>);
   inline Vector<Type> operator*(Vector<Type>);
   inline Vector<Type> operator*(Matrix<Type>);
   inline Matrix<Type> outer(Vector<Type>);

   inline Vector<Type> operator/(Type);
   inline Vector<Type> operator/(Vector<Type>);

   inline Type* begin();
   inline Type* end();

   inline void insert(int, Vector<Type>);
   inline Vector<Type> extract(int, int);
   
   inline Vector<Type> assemble(Vector<Type>);

   void load(char*);
   void save(char*);

   // DESTRUCTOR

   ~Vector();
};


// CONSTRUCTORS

/// Default constructor. It creates a vector of size zero.

template <typename Type>
Vector<Type>::Vector(void) : size(0), vector(0)
{

}


/// Constructor. It creates a vector of size n, containing n copies of the default value for Type.
///
/// @param newSize Size of Vector.

template <typename Type>
Vector<Type>::Vector(int newSize) : size(newSize), vector(new Type[newSize])
{

}


/// Constructor. It creates a vector of size n, containing n copies of the type value of Type. 
///
/// @param newSize Size of Vector.
/// @param type Value of Type.

template <typename Type> Vector<Type> ::Vector(int newSize, const Type& type) 
: size(newSize), vector(new Type[newSize])
{
   for(int i = 0; i < newSize; i++)
   {
      vector[i] = type;
   }
}


/// Constructor. It creates a vector of size n, containing n copies of the type value of Type. 
///
/// @param newSize Size of Vector.
/// @param type Value of Type.

template <typename Type> Vector<Type>::Vector(int newSize, const Type* type) 
: size(newSize), vector(new Type[newSize])
{
   for(int i = 0; i < newSize; i++)
   {
      vector[i] = *type++;
   }
}


/// Copy constructor. It creates a copy of an existing Vector. 
///
/// @param oldVector Vector to be copied.

template <typename Type> Vector<Type>::Vector(const Vector<Type>& oldVector) 
: size(oldVector.size), vector(new Type[size])
{
   for(int i = 0; i < size; i++)
   {
      vector[i] = oldVector[i];
   }
}


// ASSIGNMENT OPERATORS

/// Assignment operator. It assigns to self a copy of an existing Vector.
///
/// @param oldVector Vector to be assigned.

template <typename Type>
Vector<Type>& Vector<Type>::operator=(const Vector<Type>& oldVector)
{
   if(this != &oldVector)
   {
      if(size != oldVector.size)
      {
         if(vector != 0)
         {
            delete [] (vector);
         }

         size = oldVector.size;

         vector = new Type[size];
      }

      for(int i = 0; i < size; i++)
      {
         vector[i] = oldVector[i];
      }
   }

   return(*this);
}


// REFERENCE OPERATORS

/// Reference operator. 

template <typename Type>
inline Type& Vector<Type>::operator[](const int i) 
{
   // Control sentence (if debug)

   #ifndef NDEBUG 

   if(i >= size)
   {
      std::cerr << std::endl
                << "Flood Error: Vector Template. " << std::endl
                << "Reference operator []." << std::endl
                << "Index is " << i << " and it must be less than " << size << "." << std::endl
                << std::endl;

      exit(1);
   }

   #endif

   // Return vector element

   return(vector[i]);
}


/// Reference operator. 

template <typename Type>
inline const Type& Vector<Type>::operator[](const int i) const 
{
   // Control sentence (if debug)
   
   #ifndef NDEBUG 

   if(i >= size)
   {
      std::cerr << std::endl
                << "Flood Error: Vector Template. " << std::endl
                << "Reference operator []." << std::endl
                << "Index is " << i << " and it must be less than " << size << "." << std::endl
                << std::endl;

      exit(1);
   }

   #endif

   return(vector[i]);
}


// METHODS

// int getSize(void) method

/// This method returns the number of elements in the vector. 

template <typename Type>
inline int Vector<Type>::getSize()
{
   return(size);
}


// void setSize(int) method

/// This method sets a new size in the vector.
/// It also initializes the new vector with the previous values.
///
/// @param newSize New number of elements.

template <typename Type>
inline void Vector<Type>::setSize(int newSize)
{
   vector = new Type[newSize];

   size = newSize;
}


// void resize(int) method

/// This method sets a new size in the vector.
/// It also initializes the new vector with the previous values.
///
/// @param newSize New number of elements.

template <typename Type>
inline void Vector<Type>::resize(int newSize)
{
   if(newSize > size)
   {
      Vector<Type> newVector(newSize);   

      for(int i = 0; i < size; i++)
      {
        newVector[i] = vector[i];
      }

      vector = new Type[newSize];

      for(int i = 0; i < size; i++)
      {
         vector[i] = newVector[i];      
      }

      size = newSize;
   }
   else if(newSize < size)
   {      
      Vector<Type> newVector(newSize);   

      for(int i = 0; i < newSize; i++)
      {
        newVector[i] = vector[i];
      }

      vector = new Type[newSize];

      for(int i = 0; i < newSize; i++)
      {
         vector[i] = newVector[i];      
      }

      size = newSize;
   }
   else
   {
      // Do nothing
   }
}


// void fillAtRandom(void) method

/// This method assigns a random value comprised between -1 and 1 to each element in the vector. 

template <class Type>
inline void Vector<Type>::fillAtRandom(void)
{
   double random = 0.0;

   for(int i = 0; i < size; i++)
   {
      random = (double)rand()/(RAND_MAX+1.0);

      vector[i] = -1.0 + 2.0*random;
   }
}


// void fillAtRandom(double, double) method

/// This method assigns a random value comprised between a minimum and a maximum values to each element in the 
/// vector. 
///
/// @param minimum Minimum filling value.  
/// @param maximum Maximum filling value.

template <class Type>
inline void Vector<Type>::fillAtRandom(double minimum, double maximum)
{
   // Control sentence (if debug)

   #ifndef NDEBUG 

   if(minimum > maximum)
   {
      std::cerr << std::endl
                << "Flood Error: Vector Template." << std::endl 
                << "void fillAtRandom(double, double) method." << std::endl
                << "Minimum value must be less or equal than maximum value." << std::endl
                << std::endl;

      exit(1);
   }

   #endif

   double random = 0.0;

   for(int i = 0; i < size; i++)
   {
      random = (double)rand()/(RAND_MAX+1.0);

      vector[i] =  minimum + (maximum - minimum)*random;
   }
}


// Type calculateMinimum(void) method

/// This method returns the smallest element in the vector.

template <typename Type>
inline Type Vector<Type>::calculateMinimum()
{
   Type minimum = vector[0];

   for(int i = 1; i < size; i++)
   {
      if(vector[i] < minimum)
      {
         minimum = vector[i];
      }
   }
   
   return(minimum);
}


// Type calculateMaximum(void) method

/// This method returns the largest element in the vector.

template <typename Type>
inline Type Vector<Type>::calculateMaximum()
{
   Type maximum = vector[0];

   for(int i = 1; i < size; i++)
   {
      if(vector[i] > maximum)
      {
         maximum = vector[i];
      }
   }
   
   return(maximum);
}


// int calculateMinimalIndex(void) method

/// This method returns the index of the smallest element in the vector.

template <typename Type>
inline int Vector<Type>::calculateMinimalIndex()
{
   Type minimum = vector[0];
   int minimalIndex = 0;

   for(int i = 1; i < size; i++)
   {
      if(vector[i] < minimum)
      {
         minimum = vector[i];
         minimalIndex = i;
      }
   }
   
   return(minimalIndex);
}


// int calculateMaximalIndex(void) method

/// This method returns the index of the largest element in the vector.

template <typename Type>
inline int Vector<Type>::calculateMaximalIndex()
{
   Type maximum = vector[0];
   int maximalIndex = 0;

   for(int i = 1; i < size; i++)
   {
      if(vector[i] > maximum)
      {
         maximum = vector[i];
         maximalIndex = i;
      }
   }
   
   return(maximalIndex);
}


// Type calculateMean(void) method

/// This method returns the mean of the elements in the vector.

template <typename Type>
inline Type Vector<Type>::calculateMean(void)
{
   Type mean = 0;

   double sum = 0.0;

   for(int i = 0; i < size; i++)
   {
      sum += vector[i];
   }

   mean = sum/(double)size;
   
   return(mean);
}


// Type calculateStandardDeviation(void) method

/// This method returns the standard deviation of the elements in the vector.

template <typename Type>
inline Type Vector<Type>::calculateStandardDeviation(void)
{
   Type standardDeviation = 0;

   double mean = calculateMean();

   double sum = 0.0;
   
   for(int i = 0; i < size; i++)
   {
      sum += pow(vector[i] - mean, 2);
   }

   standardDeviation = sqrt(sum/(double)size);
   
   return(standardDeviation);
}


// double calculateNorm(void) method

/// This element returns the vector norm.

template <typename Type>
inline double Vector<Type>::calculateNorm()
{
   // Control sentence (if debug)

   #ifndef NDEBUG 

   if(size == 0)
   {
      std::cerr << std::endl
                << "Flood Error: Vector Template." << std::endl 
                << "double calculateNorm(void) method." << std::endl
                << "Size of vector is zero." << std::endl
                << std::endl;

      exit(1);
   }

   #endif

   double norm = 0.0;

   for(int i = 0; i < size; i++)
   {
      norm += vector[i]*vector[i];
   }
   
   norm = sqrt(norm);
   
   return(norm);
}


// Vector<Type> operator + (Type) method 

/// Sum vector+scalar arithmetic operator. 

template <typename Type>
inline Vector<Type> Vector<Type>::operator + (Type scalar)
{       
   Vector<Type> sum(size);

   for(int i = 0; i < size; i++)
   {
      sum[i] = vector[i] + scalar;
   }
   
   return(sum);
}


// Vector<Type> operator + (Vector<Type>)

/// Sum vector+vector arithmetic operator.

template <typename Type>
inline Vector<Type> Vector<Type>::operator + (Vector<Type> otherVector)
{       
   // Control sentence (if debug)

   #ifndef NDEBUG 

   int otherSize = otherVector.getSize();

   if(otherSize != size)
   {
      std::cerr << std::endl
                << "Flood Error: Vector Template. " << std::endl
                << "Vector<Type> operator + (Vector<Type>)." << std::endl
                << "Size of vectors is " << size << " and " << otherSize << " and they must be the same." 
                << std::endl
                << std::endl;

      exit(1);          
   }

   #endif


   Vector<Type> sum(size);
  
   for(int i = 0; i < size; i++)
   {
      sum[i] = vector[i] + otherVector[i];
   }
   
   return(sum);
}


//Vector<Type> operator - (Type) method 

/// Difference vector-scalar arithmetic operator.

template <typename Type>
inline Vector<Type> Vector<Type>::operator - (Type scalar)
{       
   Vector<Type> difference(size);

   for(int i = 0; i < size; i++)
   {
      difference[i] = vector[i] - scalar;
   }
   
   return(difference);
}


// Vector<Type> operator - (Vector<Type>)

/// Difference vector-vector arithmetic operator.

template <typename Type>
inline Vector<Type> Vector<Type>::operator - (Vector<Type> otherVector)
{       
   // Control sentence (if debug)      
       
   #ifndef NDEBUG 

   int otherSize = otherVector.getSize();

   if(otherSize != size)
   {
      std::cerr << std::endl
                << "Flood Error: Vector Template." << std::endl 
                << "Vector<Type> operator - (Vector<Type>)." << std::endl
                << "Size of vectors is " << size << " and " << otherSize << " and they must be the same." 
                << std::endl
                << std::endl;

      exit(1);          
   }
      
   #endif

   Vector<Type> difference(size);
  
   for(int i = 0; i < size; i++)
   {
      difference[i] = vector[i] - otherVector[i];
   }
   
   return(difference);
}


// Vector<Type> operator * (Type) method 

/// Product vector*scalar arithmetic operator.

template <typename Type>
inline Vector<Type> Vector<Type>::operator * (Type scalar)
{       
   Vector<Type> product(size);

   for(int i = 0; i < size; i++)
   {
      product[i] = vector[i]*scalar;
   }
   
   return(product);
}


// Type operator * (Vector<Type>)  

/// Element by element product vector*vector arithmetic operator.

template <typename Type>
inline Vector<Type> Vector<Type>::operator * (Vector<Type> otherVector)
{       
   // Control sentence (if debug)      
       
   #ifndef NDEBUG 

   int otherSize = otherVector.getSize();

   if(otherSize != size)
   {
      std::cerr << std::endl
                << "Flood Error: Vector Template." << std::endl 
                << "Vector<Type> operator * (Vector<Type>)." << std::endl
                << "Both vector sizes must be the same." << std::endl
                << std::endl;

      exit(1);          
   }

   #endif

   Vector<Type> product(size);
  
   for(int i = 0; i < size; i++)
   {
      product[i] = vector[i]*otherVector[i];
   }
   
   return(product);
}


// Vector<Type> operator * (Matrix<Type>)  

/// Product vector*matrix arithmetic operator.

template <typename Type>
inline Vector<Type> Vector<Type>::operator * (Matrix<Type> matrix)
{
   int numberOfRows = matrix.getNumberOfRows();

   // Control sentence (if debug)      

   #ifndef NDEBUG 

   if(numberOfRows != size)
   {
      std::cerr << std::endl
                << "Flood Error: Vector Template." << std::endl 
                << "Type operator * (Matrix<Type>)." << std::endl
                << "Matrix number of rows must be equal to vector size." << std::endl
                << std::endl;

      exit(1);          
   }

   #endif

   int numberOfColumns = matrix.getNumberOfColumns();

   Vector<Type> product(numberOfColumns);
  
   for(int j = 0; j < numberOfColumns; j++)
   {     
      product[j] = 0;      

      for(int i = 0; i < numberOfRows; i++)
      {
         product[j] += vector[i]*matrix[i][j];
      }
   }
    
   return(product);
}


// Vector<Type> dot(Vector<Type>) method

/// Dot product vector*vector arithmetic operator.

template <typename Type>
inline Type Vector<Type>::dot(Vector<Type> otherVector)
{            
   // Control sentence (if debug)      

   #ifndef NDEBUG 

   int otherSize = otherVector.getSize();

   if(otherSize != size)
   {
      std::cerr << std::endl
                << "Flood Error: Vector Template." << std::endl 
                << "Type dot(Vector<Type>) method." << std::endl
                << "Both vector sizes must be the same." << std::endl
                << std::endl;

      exit(1);          
   }

   #endif

   Type dotProduct = 0;
  
   for(int i = 0; i < size; i++)
   {
      dotProduct += vector[i]*otherVector[i];
   }
   
   return(dotProduct);
}


// Matrix<Type> outer(Vector<Type>) method

/// Outer product vector*vector arithmetic operator.

template <typename Type>
inline Matrix<Type> Vector<Type>::outer(Vector<Type> otherVector)
{            
   int otherSize = otherVector.getSize();

   // Control sentence (if debug)      

   #ifndef NDEBUG 

   if(otherSize != size)
   {
      std::cerr << std::endl
                << "Flood Error: Vector Template." << std::endl 
                << "Matrix<Type> outer(Vector<Type>) method." << std::endl
                << "Both vector sizes must be the same." << std::endl
                << std::endl;

      exit(1);          
   }

   #endif

   int numberOfRows = size;
   int numberOfColumns = otherSize;
   
   Matrix<Type> outer(numberOfRows, numberOfColumns);
     
   for(int i = 0;  i < numberOfRows; i++)
   {
      for(int j = 0;  j < numberOfColumns; j++)
      {
         outer[i][j] = vector[i]*otherVector[j];
      }           
   }
   
   return(outer);
}


//Vector<Type> operator / (Type) method 

/// Cocient vector/scalar arithmetic operator.

template <typename Type>
inline Vector<Type> Vector<Type>::operator / (Type scalar)
{       
   Vector<Type> cocient(size);

   for(int i = 0; i < size; i++)
   {
      cocient[i] = vector[i]/scalar;
   }
   
   return(cocient);
}


// Vector<Type> operator - (Vector<Type>)

/// Cocient vector/vector arithmetic operator.

template <typename Type>
inline Vector<Type> Vector<Type>::operator / (Vector<Type> otherVector)
{       
   int otherSize = otherVector.getSize();

   // Control sentence (if debug)            

   #ifndef NDEBUG 

   if(otherSize != size)
   {
      std::cerr << std::endl
                << "Flood Error: Vector Template." << std::endl 
                << "Vector<Type> operator - (Vector<Type>)." << std::endl
                << "Both vector sizes must be the same." << std::endl
                << std::endl;

      exit(1);          
   }

   #endif

   Vector<Type> cocient(size);
  
   for(int i = 0; i < size; i++)
   {
      cocient[i] = vector[i]/otherVector[i];
   }
   
   return(cocient);
}


// Type* begin(void) method

/// This method returns a pointer to the first element in the container.

template <typename Type>
inline Type* Vector<Type>::begin()
{
   return(vector);
}


// Type* end(void) method

/// This method returns a pointer to the last element in the container. 

template <typename Type>
inline Type* Vector<Type>::end()
{
   return(vector + size);
}


// DESTRUCTOR

/// Destructor. 

template <typename Type> 
Vector<Type>::~Vector()
{
   if(vector != 0)
   {
      delete[](vector);
   }
}


// Input operator

/// This method re-writes the input operator >> for the Vector template. 

template<typename Type>
std::istream& operator>>(std::istream& is, Vector<Type>& v)
{
   int size = v.getSize();
   
   for(int i = 0; i < size; i++)
   {
      is >> v[i];
   }

   return(is);
}


// Output operator

/// This method re-writes the output operator << for the Vector template. 

template<typename Type>
std::ostream& operator<<(std::ostream& os, Vector<Type>& v)
{
   int size = v.getSize();
   
   for(int i = 0; i < size; i++)
   {
      os << v[i] << " ";
   }

   return(os);
}


// void load(char*) method

/// This method loads the elements of a vector from a data file.
/// The file format is as follows:
/// 
/// NumberOfElements
/// Element_0 Element_1 ... Element_N-1
///
/// @param filename Filename.
///
/// @see save(char*).

template <class Type>
inline void Vector<Type>::load(char* filename)
{
   std::fstream file;

   // Open file
    
   file.open(filename, std::ios::in);

   if(!file.is_open())
   {
      std::cerr << std::endl
                << "Flood Error: Vector template." << std::endl
                << "void load(char*) method." << std::endl
                << "Cannot open vector data file."  << std::endl
                << std::endl;

      exit(1);
   }
   else
   {
      std::cout << std::endl
                << "Loading vector from data file..." << std::endl;
   }

   // Read file

   file >> size;   

   resize(size);

   for(int i = 0; i < size; i++)
   {
      file >> vector[i];
   }

   // Close file

   file.close();
}


// void save(char*) method

/// This method saves to a data file the elements of the vector.
///
/// @param filename Filename.
///
/// @see load(char*).

template <class Type>
inline void Vector<Type>::save(char* filename)
{
   std::fstream file; 

   // Open file

   file.open(filename, std::ios::out);

   if(!file.is_open())
   {
      std::cerr << std::endl 
                << "Flood Error: Vector template." << std::endl
                << "void save(char*) method." << std::endl
                << "Cannot open vector data file." << std::endl
                << std::endl;

      exit(1);
   }
   else
   {
      std::cout << std::endl 
                << "Saving vector to data file..." 
                << std::endl;
   }

   // Write file

   file << size << std::endl;

   for(int i = 0; i < size; i++)
   {
      file << vector[i] << " ";
   }

   file << std::endl;

   // Close file

   file.close();
}


// void insert(int, Vector<Type>) method

/// Insert another vector starting from a given position.

template <typename Type>
inline void Vector<Type>::insert(int position, Vector<Type> otherVector)
{       
   int otherSize = otherVector.getSize();

   // Control sentence (if debug)            

   #ifndef NDEBUG 

   if(position + otherSize > size)
   {
      std::cerr << std::endl
                << "Flood Error: Vector Template." << std::endl 
                << "void insert(Vector<Type>, int) method." << std::endl
                << "Cannot insert vector." << std::endl
                << std::endl;

      exit(1);          
   }

   #endif

   for(int i = 0; i < otherSize; i++)
   {
      vector[position + i] = otherVector[i];
   }
}


// Vector<Type> extract(int, int) method

/// Extract a vector of a given size from a given position

template <typename Type>
inline Vector<Type> Vector<Type>::extract(int position, int otherSize)
{            
   // Control sentence (if debug)            

   #ifndef NDEBUG 

   if(position + otherSize > size)
   {
      std::cerr << std::endl
                << "Flood Error: Vector Template." << std::endl 
                << "Vector<Type> extract(int, int) method." << std::endl
                << "Cannot extract vector." << std::endl
                << std::endl;

      exit(1);          
   }

   #endif

   Vector<Type> otherVector(otherSize);

   for(int i = 0; i < otherSize; i++)
   {
      otherVector[i] = vector[position + i];
   }

   return(otherVector);
}


// Vector<Type> assemble(Vector<Type>) method

/// Assemble two vectors.

template <typename Type>
inline Vector<Type> Vector<Type>::assemble(Vector<Type> otherVector)
{            
   int otherSize = otherVector.getSize();

   Vector<double> assembly(size + otherSize);

   for(int i = 0; i < size; i++)
   {
      assembly[i] = vector[i];
   }

   for(int i = 0; i < otherSize; i++)
   {
      assembly[size+i] = otherVector[i];
   } 
   
   return(assembly);
}


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
