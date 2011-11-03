/****************************************************************************************************************/
/*                                                                                                              */
/*   Flood: An Open Source Neural Networks C++ Library                                                          */
/*   www.cimne.com/flood                                                                                        */
/*                                                                                                              */
/*   M A T R I X   C O N T A I N E R                                                                            */
/*                                                                                                              */
/*   Roberto Lopez                                                                                              */
/*   International Center for Numerical Methods in Engineering (CIMNE)                                          */
/*   Technical University of Catalonia (UPC)                                                                    */
/*   Barcelona, Spain                                                                                           */
/*   E-mail: rlopez@cimne.upc.edu                                                                               */
/*                                                                                                              */
/****************************************************************************************************************/


#ifndef __MATRIX_H__
#define __MATRIX_H__

// System includes

#include <iostream>
#include <fstream>
#include <math.h>

// Utilities includes

#include "Vector.h"

namespace Flood
{

/// This template class defines a matrix for general purpose use.
///
/// @see Vector.

template <class Type>
class Matrix 
{

private:

   /// Number of rows in matrix.

   int numberOfRows;

   /// Number of columns in matrix.

   int numberOfColumns;

   /// Double pointer to a Type.

   Type** matrix;

public:

   // CONSTRUCTORS

   Matrix(void);

   Matrix(int, int);

   Matrix(int, int, const Type&);

   Matrix(int, int, const Type*);

   Matrix(const Matrix&);


   // ASSIGNMENT OPERATOR

   Matrix& operator=(const Matrix&);


   // REFERENCE OPERATORS

   inline Type* operator[](const int);

   inline const Type* operator[](const int) const;


   // METHODS

   inline int getNumberOfRows(void) const;
   inline int getNumberOfColumns(void) const;  

   inline void resize(int, int);
    
   inline Vector<Type> getRow(int);
   inline Vector<Type> getColumn(int);

   inline Matrix<Type> getTranspose(void);

   inline void setRow(int, Vector<Type>);
   inline void setColumn(int, Vector<Type>);

   inline Matrix<Type> operator+(Type);
   inline Matrix<Type> operator+(Matrix<Type>);

   inline Matrix<Type> operator-(Type);
   inline Matrix<Type> operator-(Matrix<Type>);

   inline Matrix<Type> operator*(Type);
   inline Vector<Type> operator*(Vector<Type>);
   inline Matrix<Type> operator*(Matrix<Type>);

   inline Matrix<Type> operator/(Type);

   inline void fillAtRandom(void);
   inline void fillAtRandom(double, double);

   inline void setToIdentity(void);

   inline Type calculateDeterminant(void);

   inline Matrix<Type> calculateInverse(void);

   void load(char*);
   void save(char*);

   // DESTRUCTOR

   ~Matrix();
};


// CONSTRUCTORS


/// Default constructor. It creates a matrix with zero rows and zero columns.

template <class Type>
Matrix<Type>::Matrix() : numberOfRows(0), numberOfColumns(0), matrix(0) 
{

}


/// Constructor. It creates a matrix with n rows and m columns, containing n*m copies of the default value for 
/// Type. 
///
/// @param newNumberOfRows Number of rows in Matrix.
/// @param newNumberOfColumns Number of columns in Matrix.

template <class Type>
Matrix<Type>::Matrix(int newNumberOfRows, int newNumberOfColumns) 
: numberOfRows(newNumberOfRows), numberOfColumns(newNumberOfColumns), 
matrix(new Type*[newNumberOfRows])
{
   matrix[0] = new Type[numberOfColumns*numberOfRows];

   for(int i = 1; i < numberOfRows; i++)
   {
      matrix[i] = matrix[i-1] + numberOfColumns;
   }
}


/// Constructor. It creates a matrix with n rows and m columns, containing n*m copies of the type value of Type. 
///
/// @param newNumberOfRows Number of rows in Matrix.
/// @param newNumberOfColumns Number of columns in Matrix.
/// @param type Value of Type.

template <class Type>
Matrix<Type>::Matrix(int newNumberOfRows, int newNumberOfColumns, const Type& type) 
: numberOfRows(newNumberOfRows), numberOfColumns(newNumberOfColumns), 
matrix(new Type*[newNumberOfRows])
{
   matrix[0] = new Type[numberOfColumns*numberOfRows];

   for(int i = 1; i < numberOfRows; i++)
   {
      matrix[i] = matrix[i-1] + numberOfColumns;
   }

   for(int i = 0; i < numberOfRows; i++)
   {
      for(int j = 0; j < numberOfColumns; j++)
      {
         matrix[i][j] = type;
      }
   }
}



/// Constructor. It creates a matrix with n rows and m columns, containing n*m copies of the type value of Type. 
///
/// @param newNumberOfRows Number of rows in Matrix.
/// @param newNumberOfColumns Number of columns in Matrix.
/// @param type Value of Type.

template <class Type>
Matrix<Type>::Matrix(int newNumberOfRows, int newNumberOfColumns, const Type* type) 
: numberOfRows(newNumberOfRows), numberOfColumns(newNumberOfColumns), 
matrix(new Type*[newNumberOfRows])
{
   // Construct matrix

   matrix[0] = new Type[numberOfColumns*numberOfRows];

   for(int i = 1; i < numberOfRows; i++)
   {
      matrix[i] = matrix[i-1] + numberOfColumns;
   }

   // Set all elements of matrix to the type value of Type

   for(int i = 0; i < numberOfRows; i++)
   {
      for(int j = 0; j < numberOfColumns; j++)
      {
         matrix[i][j] = *type++;
     }
   }
}


/// Copy constructor. It creates a copy of an existing matrix. 
///
/// @param oldMatrix Matrix to be copied.

template <class Type>
Matrix<Type>::Matrix(const Matrix& oldMatrix) 
: numberOfRows(oldMatrix.numberOfRows), numberOfColumns(oldMatrix.numberOfColumns), 
matrix(new Type*[numberOfRows])
{
   // Construct matrix

   matrix[0] = new Type[numberOfColumns*numberOfRows];

   for(int i = 1; i < numberOfRows; i++)
   {
      matrix[i] = matrix[i-1] + numberOfColumns;
   }

   // Set all elements of matrix to the old matrix type values

   for(int i = 0; i < numberOfRows; i++)
   {
      for(int j = 0; j < numberOfColumns; j++)
      {
         matrix[i][j] = oldMatrix[i][j];
      }
   }
}


// ASSIGNMENT OPERATORS

/// Assignment operator. It assigns to self a copy of an existing Matrix.
///
/// @param oldMatrix Matrix to be assigned.

template <class Type>
Matrix<Type>& Matrix<Type>::operator=(const Matrix<Type>& oldMatrix)
{
   if(this != &oldMatrix) 
   {
      if(numberOfRows != oldMatrix.numberOfRows 
      || numberOfColumns != oldMatrix.numberOfColumns) 
      {
         if(matrix != 0) 
         {
            delete[] (matrix[0]);

            delete[] (matrix);
         }

         numberOfRows = oldMatrix.numberOfRows;

         numberOfColumns = oldMatrix.numberOfColumns;

         matrix = new Type*[numberOfRows];

         matrix[0] = new Type[numberOfColumns*numberOfRows];
      }

      for(int i = 1; i < numberOfRows; i++)
      {
         matrix[i] = matrix[i-1] + numberOfColumns;
      }

      // Set all elements of matrix to the old matrix type values

      for(int i = 0; i < numberOfRows; i++)
      {
         for(int j = 0; j < numberOfColumns; j++)
         {
            matrix[i][j] = oldMatrix[i][j];
         }
      }
   }

   return(*this);
}


// REFERENCE OPERATORS

/// Reference operator.  

template <class Type>
inline Type* Matrix<Type>::operator[](const int i) 
{
   // Control sentence (if debug)

   #ifndef NDEBUG 
   
   if(i >= numberOfRows)
   {
      std::cerr << std::endl
                << "Flood Error: Matrix Template. " << std::endl
                << "Reference operator []." << std::endl
                << "Index must be less than matrix dimensions." << std::endl
                << std::endl;

      exit(1);
   }

   #endif

   // Return matrix element

   return(matrix[i]);
}


/// Reference operator.  

template <class Type>
inline const Type* Matrix<Type>::operator[](const int i) const
{
   // Control sentence (if debug)

   #ifndef NDEBUG 
   
   if(i >= numberOfRows)
   {
      std::cerr << std::endl
                << "Flood Error: Matrix Template. " << std::endl
                << "Reference operator []." << std::endl
                << "Index must be and less than matrix dimensions." << std::endl
                << std::endl;

      exit(1);
   }

   #endif

   // Return matrix element

   return(matrix[i]);
}


// METHODS

// int getNumberOfRows(void) method

/// This method returns the number of rows in the matrix. 

template <class Type>
inline int Matrix<Type>::getNumberOfRows() const
{
   return(numberOfRows);
}


// int getNumberOfColumns(void) method

/// This method returns the number of columns in the matrix. 

template <class Type>
inline int Matrix<Type>::getNumberOfColumns() const
{
   return(numberOfColumns);
}


// void resize(int, int) method

/// This method sets new numbers of rows and columns in the vector.
/// It does not initialize the new matrix with the previous values. 
///
/// @param newNumberOfRows New number of rows.
/// @param newNumberOfColumns New number of columns.

template <class Type>
inline void Matrix<Type>::resize(int newNumberOfRows, int newNumberOfColumns) 
{
   numberOfRows = newNumberOfRows;
   numberOfColumns = newNumberOfColumns;

   matrix = new Type*[numberOfRows];

   matrix[0] = new Type[numberOfColumns*numberOfRows];

   for(int i = 1; i < numberOfRows; i++)
   {
      matrix[i] = matrix[i-1] + numberOfColumns;
   }
}


// Matrix<Type> getTranspose(void) method

/// This method returns the transpose of the matrix. 

template <class Type>
inline Matrix<Type> Matrix<Type>::getTranspose(void)
{
   Matrix<Type> transpose(numberOfColumns, numberOfRows);

   for(int i = 0; i < numberOfColumns; i++)
   {
      for(int j = 0; j < numberOfRows; j++)
      {
         transpose[i][j] = matrix[j][i];
      }     
   }

   return(transpose);
}


// Vector<Type> getRow(int) method

/// This method returns the row i of the matrix. 

template <class Type>
inline Vector<Type> Matrix<Type>::getRow(int i)
{
   Vector<Type> row(numberOfColumns, 0.0);

   for(int j = 0; j < numberOfColumns; j++)
   {
      row[j] = matrix[i][j];
   }

   return(row);
}


// Vector<Type> getColumn(int) method

/// This method returns the column j of the matrix. 

template <class Type>
inline Vector<Type> Matrix<Type>::getColumn(int j)
{
   Vector<Type> column(numberOfRows, 0.0);

   for(int i = 0; i < numberOfRows; i++)
   {
      column[i] = matrix[i][j];
   }

   return(column);
}


// void setRow(int, Vector<Type>) method

/// This method sets new values of a single row in the matrix. 
///
/// @param rowIndex Index of row. 
/// @param newRow New values of single row. 

template <class Type>
inline void Matrix<Type>::setRow(int rowIndex, Vector<Type> newRow)
{
   // Control sentence (if debug)

   #ifndef NDEBUG 

   if(rowIndex >= numberOfRows)
   {
      std::cerr << std::endl
                << "Flood Error: Matrix Template. " << std::endl
                << "setRow(int, Vector<Type>) method." << std::endl
                << "Index must be less than number of rows." << std::endl
                << std::endl;

      exit(1);
   }

   int size = newRow.getSize();

   if(size != numberOfColumns)
   {
      std::cerr << std::endl
                << "Flood Error: Matrix Template. " << std::endl
                << "setRow(int, Vector<Type>) method." << std::endl
                << "Size must be equal to number of columns." << std::endl
                << std::endl;

      exit(1);
   }

   #endif

   // Set new row

   for(int i = 0; i < numberOfColumns; i++)
   {
      matrix[rowIndex][i] = newRow[i];
   }
}


// void setColumn(int, Vector<Type>) method

/// This method sets new values of a single column in the matrix. 
///
/// @param columnIndex Index of column. 
/// @param newColumn New values of single column. 

template <class Type>
inline void Matrix<Type>::setColumn(int columnIndex, Vector<Type> newColumn)
{
   // Control sentence (if debug)

   #ifndef NDEBUG 

   if(columnIndex >= numberOfColumns)
   {
      std::cerr << std::endl
                << "Flood Error: Matrix Template. " << std::endl
                << "setColumn(int, Vector<Type>)." << std::endl
                << "Index must be less than number of columns." << std::endl
                << std::endl;

      exit(1);
   }

   int size = newColumn.getSize();

   if(size != numberOfRows)
   {
      std::cerr << std::endl
                << "Flood Error: Matrix Template. " << std::endl
                << "setColumn(int, Vector<Type>)." << std::endl
                << "Size must be equal to number of rows." << std::endl
                << std::endl;

      exit(1);
   }

   #endif

   // Set new column

   for(int i = 0; i < numberOfRows; i++)
   {
      matrix[i][columnIndex] = newColumn[i];
   }
}


// void fillAtRandom(void) method

/// This method fills all the elements in the matrix with random values comprised between -1 and 1.

template <class Type>
inline void Matrix<Type>::fillAtRandom(void)
{
   double random = 0.0;

   for(int i = 0; i < numberOfRows; i++)
   {
      for(int j = 0; j < numberOfColumns; j++)
     {
         random = (double)rand()/(RAND_MAX+1.0);

         matrix[i][j] = -1.0 + 2.0*random;
     }
   }
}


// void fillAtRandom(double, double) method

/// This method fills all the elements in the matrix with random values comprised between a minimum and a maximum
/// values.
///
/// @param minimum Minimum possible value. 
/// @param maximum Maximum possible value. 

template <class Type>
inline void Matrix<Type>::fillAtRandom(double minimum, double maximum)
{
   // Control sentence (if debug)

   #ifndef NDEBUG 

   if(minimum > maximum)
   {
      std::cerr << std::endl
                << "Flood Error: Matrix Template." << std::endl 
                << "void fillAtRandom(double, double) method." << std::endl
                << "Minimum value must be less or equal than maximum value." << std::endl
                << std::endl;

      exit(1);
   }

   #endif

   double random = 0.0;

   for(int i = 0; i < numberOfRows; i++)
   {
      for(int j = 0; j < numberOfColumns; j++)
     {
         random = (double)rand()/(RAND_MAX+1.0);

         matrix[i][j] = minimum + (maximum - minimum)*random;
     }
   }
}


// void setToIdentity(void) method

/// This method sets the diagonal elements in the matrix with ones and the rest elements with zeros. The matrix 
/// must be square. 

template <class Type>
inline void Matrix<Type>::setToIdentity(void)
{
   // Control sentence (if debug)

   #ifndef NDEBUG 
   
   if(numberOfRows != numberOfColumns)
   {
      std::cout << std::endl
                << "Flood Error: Matrix Template." << std::endl
                << "setToIdentity(void) method." << std::endl
                << "Matrix must be square." << std::endl
                << std::endl;
      
      exit(1);
   }

   #endif

   for(int i = 0; i < numberOfRows; i++)
   {
      for(int j = 0; j < numberOfColumns; j++)
      {
         if(i==j)
         {
            matrix[i][j] = 1;
         }
         else
         {
            matrix[i][j] = 0;      
         }
      }
   }
}



// Type calculateDeterminant(void) method

/// This method returns the determinant of a square matrix. 

template <class Type>
inline Type Matrix<Type>::calculateDeterminant(void)
{ 
   // Control sentence (if debug)

   #ifndef NDEBUG 

   if(numberOfRows != numberOfColumns)
   {
      std::cerr << std::endl
                << "Flood Error: Matrix Template." << std::endl
                << "calculateDeterminant(void) method." << std::endl
                << "Matrix must be square." << std::endl
                << std::endl;
      
      exit(1);
   }

   #endif

   Type determinant = 0;

   if(numberOfRows == 0)
   {
      std::cerr << std::endl
                << "Flood Error: Matrix Template." << std::endl
                << "calculateDeterminant(void) method." << std::endl
                << "Size of matrix is zero." << std::endl
                << std::endl;
      
      exit(1);                   
   }
   else if(numberOfRows == 1)
   {
      determinant = matrix[0][0];                   
   }
   else if(numberOfRows == 2)
   {
      determinant = matrix[0][0]*matrix[1][1] - matrix[1][0]*matrix[0][1];
   }
   else
   {
      for(int j1 = 0; j1 < numberOfRows; j1++) 
      {
         Matrix<double> subMatrix(numberOfRows-1, numberOfColumns-1, 0.0);     
     
         for(int i = 1; i < numberOfRows; i++) 
         {
            int j2 = 0;
      
            for(int j = 0; j < numberOfColumns; j++) 
            {
               if(j == j1)
               {
                  continue;
               }

               subMatrix[i-1][j2] = matrix[i][j];

               j2++;
            }
         }
   
         determinant += pow(-1.0, j1+2.0)*matrix[0][j1]*subMatrix.calculateDeterminant();    
      }
   }
     
   return(determinant);
}


// Matrix<Type> calculateInverse(void) method

/// This method returns the inverse of a square matrix.
/// An error message is printed if the matrix is singular.

template <class Type>
inline Matrix<Type> Matrix<Type>::calculateInverse(void)
{
   // Control sentence (if debug)

   #ifndef NDEBUG 
   
   if(numberOfRows != numberOfColumns)
   {
      std::cerr << std::endl
                << "Flood Error: Matrix Template." << std::endl
                << "calculateDeterminant(void) method." << std::endl
                << "Matrix must be square." << std::endl
                << std::endl;
      
      exit(1);
   }

   #endif

   double determinant = calculateDeterminant();

   if(determinant == 0.0)
   {
      std::cerr << std::endl
                << "Flood Error: Matrix Template." << std::endl
                << "calculateInverse(void) method." << std::endl
                << "Matrix is singular." << std::endl
                << std::endl;
      
      exit(1);
   }


   Matrix<Type> inverse(numberOfRows, numberOfColumns);
   
   // Get cofactor matrix
   
   Matrix<double> cofactor(numberOfRows, numberOfColumns, 0.0);
                  
   Matrix<double> c(numberOfRows-1, numberOfColumns-1, 0.0);

   for(int j = 0; j < numberOfRows; j++) 
   {
      for(int i = 0; i < numberOfRows; i++) 
      {
         // Form the adjoint a[i][j]

         int i1 = 0;

         for(int ii = 0; ii < numberOfRows; ii++) 
         {
            if(ii == i)
            {
               continue;
            }
            
            int j1 = 0;

            for(int jj = 0; jj < numberOfRows; jj++) 
            {
               if(jj == j)
               {
                  continue;
               }

               c[i1][j1] = matrix[ii][jj];
               j1++;
            }
            i1++;
         }

         double determinant = c.calculateDeterminant();

         cofactor[i][j] = pow(-1.0, i+j+2.0)*determinant;
      }
   }

   // Adjoint matrix is the transpose of cofactor matrix

   Matrix<double> adjoint(numberOfRows, numberOfColumns, 0.0);
   
   for(int i = 0; i < numberOfRows; i++) 
   {
      for(int j = 0; j < numberOfRows; j++) 
      {
         adjoint[i][j] = cofactor[j][i];
      }
   }

   // Inverse matrix is adjoint matrix divided by matrix determinant
   
   for(int i = 0; i < numberOfRows; i++)
   {
      for(int j = 0; j < numberOfRows; j++)
      {
         inverse[i][j] = adjoint[i][j]/determinant;
      }        
   } 

   return(inverse);
}


// Matrix<Type> operator + (Type) 

/// Sum matrix+scalar arithmetic operator. 

template <typename Type>
inline Matrix<Type> Matrix<Type>::operator + (Type scalar)
{
   Matrix<Type> sum(numberOfRows, numberOfColumns);

   for(int i = 0; i < numberOfRows; i++)
   {
      for(int j = 0; j < numberOfColumns; j++)
      {
         sum[i][j] = matrix[i][j] + scalar;    
      }     
   }

   return(sum);
}


// Matrix<Type> operator + (Matrix<Type>) 

/// Sum matrix+matrix arithmetic operator. 

template <typename Type>
inline Matrix<Type> Matrix<Type>::operator + (Matrix<Type> otherMatrix)
{
   // Control sentence (if debug)

   #ifndef NDEBUG 

   int otherNumberOfRows = otherMatrix.getNumberOfRows();    
   int otherNumberOfColumns = otherMatrix.getNumberOfColumns();    
       
   if(otherNumberOfRows != numberOfRows || otherNumberOfColumns != numberOfColumns)
   {
      std::cerr << std::endl
                << "Flood Error: Matrix Template." << std::endl 
                << "Matrix<Type> operator + (Matrix<Type>)." << std::endl
                << "Both matrix sizes must be the same." << std::endl
                << std::endl;

      exit(1);          
   }

   #endif

   Matrix<Type> sum(numberOfRows, numberOfColumns);

   for(int i = 0; i < numberOfRows; i++)
   {
      for(int j = 0; j < numberOfColumns; j++)
      {
         sum[i][j] = matrix[i][j] + otherMatrix[i][j];    
      }     
   }

   return(sum);
}


// Matrix<Type> operator - (Type) 

/// Difference matrix-scalar arithmetic operator. 

template <typename Type>
inline Matrix<Type> Matrix<Type>::operator - (Type scalar)
{
   Matrix<Type> difference(numberOfRows, numberOfColumns);

   for(int i = 0; i < numberOfRows; i++)
   {
      for(int j = 0; j < numberOfColumns; j++)
      {
         difference[i][j] = matrix[i][j] - scalar;    
      }     
   }

   return(difference);
}


// Matrix<Type> operator - (Matrix<Type>) 

/// Difference matrix-matrix arithmetic operator. 

template <typename Type>
inline Matrix<Type> Matrix<Type>::operator - (Matrix<Type> otherMatrix)
{
   // Control sentence (if debug)

   #ifndef NDEBUG 

   int otherNumberOfRows = otherMatrix.getNumberOfRows();    
   int otherNumberOfColumns = otherMatrix.getNumberOfColumns();    
       
   if(otherNumberOfRows != numberOfRows || otherNumberOfColumns != numberOfColumns)
   {
      std::cerr << std::endl
                << "Flood Error: Matrix Template." << std::endl 
                << "Matrix<Type> operator - (Type)." << std::endl
                << "Both matrix sizes must be the same." << std::endl
                << std::endl;

      exit(1);          
   }

   #endif

   Matrix<Type> difference(numberOfRows, numberOfColumns);

   for(int i = 0; i < numberOfRows; i++)
   {
      for(int j = 0; j < numberOfColumns; j++)
      {
         difference[i][j] = matrix[i][j] - otherMatrix[i][j];    
      }     
   }

   return(difference);
}


// Matrix<Type> operator * (Type) 

/// Product matrix*scalar arithmetic operator. 

template <typename Type>
inline Matrix<Type> Matrix<Type>::operator * (Type scalar)
{
   Matrix<Type> product(numberOfRows, numberOfColumns);

   for(int i = 0; i < numberOfRows; i++)
   {
      for(int j = 0; j < numberOfColumns; j++)
      {
          product[i][j] = matrix[i][j]*scalar;     
      }      
   }

   return(product);
}


// Vector<Type> operator * (Vector<Type>) 

/// Product matrix*vector arithmetic operator. 

template <typename Type>
inline Vector<Type> Matrix<Type>::operator * (Vector<Type> vector)
{
   // Control sentence (if debug)

   #ifndef NDEBUG 
       
   int size = vector.getSize();

   if(size != numberOfColumns)
   {
      std::cerr << std::endl
                << "Flood Error: Matrix Template." << std::endl 
                << "Vector<Type> operator * (Vector<Type>)." << std::endl
                << "Vector size must be equal to matrix number of columns." 
                << std::endl
                << std::endl;

      exit(1);          
   }

   #endif

   // Calculate matrix-vector poduct   
      
   Vector<Type> product(numberOfRows);

   for(int i = 0; i < numberOfRows; i++)
   {     
      product[i] = 0;      

      for(int j = 0; j < numberOfColumns; j++)
      {
         product[i] += vector[j]*matrix[i][j];
      }
   }

   return(product);
}


// Matrix<Type> operator * (Matrix<Type>) 

/// Product matrix*matrix arithmetic operator. 

template <typename Type>
inline Matrix<Type> Matrix<Type>::operator * (Matrix<Type> otherMatrix)
{
   // Control sentence

   int otherNumberOfColumns = otherMatrix.getNumberOfColumns();    

   Matrix<Type> product(numberOfRows, otherNumberOfColumns, 0.0);

   for(int i = 0; i < numberOfRows; i++)
   {
      for(int j = 0; j < otherNumberOfColumns; j++)
      { 
         for(int k = 0; k < numberOfColumns; k++)
         {
             product[i][j] += matrix[i][k]*otherMatrix[k][j];
         }
      }
   }

   return(product);
}


// Matrix<Type> operator / (Type) 

/// Cocient Matrix/scalar arithmetic operator. 

template <typename Type>
inline Matrix<Type> Matrix<Type>::operator / (Type scalar)
{
   Matrix<Type> cocient(numberOfRows, numberOfColumns);

   for(int i = 0; i < numberOfRows; i++)
   {
      for(int j = 0; j < numberOfColumns; j++)
      {
         cocient[i][j] = matrix[i][j]/scalar;     
      }      
   }

   return(cocient);
}


// DESTRUCTOR

/// Destructor. 

template <class Type>
Matrix<Type>::~Matrix()
{
   if(matrix != 0) 
   {
      delete[] (matrix[0]);

      delete[] (matrix);
   }
}


/// This method re-writes the input operator >> for the Vector template. 

template<typename Type>
std::istream& operator>>(std::istream& is, Matrix<Type>& m)
{
   int numberOfRows = m.getNumberOfRows();
   int numberOfColumns = m.getNumberOfColumns();

   for(int i = 0; i < numberOfRows; i++)
   {
      for(int j = 0; j < numberOfColumns; j++)
      {
         is >> m[i][j];
      }     
   }

   return(is);
}


// Output operator

/// This method re-writes the output operator << for the Vector template. 

template<typename Type>
std::ostream& operator<<(std::ostream& os, Matrix<Type>& m)
{
   int numberOfRows = m.getNumberOfRows();
   int numberOfColumns = m.getNumberOfColumns();
   
   for(int i = 0; i < numberOfRows; i++)
   {
      for(int j = 0; j < numberOfColumns; j++)
      {
         os << m[i][j] << " ";
      }     

      os << std::endl;
   }

   return(os);
}


// void load(char*) method

/// This method loads the numbers of rows and columns and the values of the matrix from a data file. 
///
/// @param filename Filename.

template <class Type>
inline void Matrix<Type>::load(char* filename)
{
   std::fstream file;

   // Open file
    
   file.open(filename, std::ios::in);

   if(!file.is_open())
   {
      std::cerr << std::endl
                << "Flood Error: Matrix template." << std::endl
                << "void load(char*) method." << std::endl
                << "Cannot open matrix data file."  << std::endl
                << std::endl;

      exit(1);
   }
   else
   {
      std::cout << std::endl
                << "Loading matrix from data file..." << std::endl;
   }

   // Read file

   file >> numberOfRows;   
   file >> numberOfColumns;   

   resize(numberOfRows, numberOfColumns);

   file >> matrix;

   // Close file

   file.close();
}


// void save(char*) method

/// This method saves the numbers of rows and columns and the values of the matrix to a data file. 
///
/// @param filename Filename.

template <class Type>
inline void Matrix<Type>::save(char* filename)
{
   std::fstream file; 
   
   // Open file

   file.open(filename, std::ios::out);

   if(!file.is_open())
   {
      std::cerr << std::endl 
                << "Flood Error: Matrix template." << std::endl
                << "void save(char*) method." << std::endl
                << "Cannot open matrix data file." << std::endl
                << std::endl;

      exit(1);
   }
   else
   {
      std::cout << std::endl 
                << "Saving matrix to data file..." 
                << std::endl;
   }

   // Write file

   file << numberOfRows << " " << numberOfColumns << std::endl;

   for(int i = 0; i < numberOfRows; i++)
   {
      for(int j = 0; j < numberOfColumns; j++)
     {
         file << matrix[i][j] << " ";
     }

      file << std::endl;
   }

   // Close file

   file.close();
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
