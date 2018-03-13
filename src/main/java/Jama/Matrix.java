package Jama;

import qual.Immutable;
import qual.Readonly;
import qual.Mutable;
import java.io.Serializable;
import java.text.NumberFormat;
import java.text.DecimalFormat;
import java.text.DecimalFormatSymbols;
import java.util.Locale;
import java.text.FieldPosition;
import java.io.PrintWriter;
import java.io.BufferedReader;
import java.io.StreamTokenizer;
import Jama.util.*;
import qual.ReceiverDependantMutable;

/**
   Jama = Java Matrix class.
<P>
   The Java Matrix Class provides the fundamental operations of numerical
   linear algebra.  Various constructors create Matrices from two dimensional
   arrays of double precision floating point numbers.  Various "gets" and
   "sets" provide access to submatrices and matrix elements.  Several methods
   implement basic matrix arithmetic, including matrix addition and
   multiplication, matrix norms, and element-by-element array operations.
   Methods for reading and printing matrices are also included.  All the
   operations in this version of the Matrix Class involve real matrices.
   Complex matrices may be handled in a future version.
<P>
   Five fundamental matrix decompositions, which consist of pairs or triples
   of matrices, permutation vectors, and the like, produce results in five
   decomposition classes.  These decompositions are accessed by the Matrix
   class to compute solutions of simultaneous linear equations, determinants,
   inverses and other matrix functions.  The five decompositions are:
<P><UL>
   <LI>Cholesky Decomposition of symmetric, positive definite matrices.
   <LI>LU Decomposition of rectangular matrices.
   <LI>QR Decomposition of rectangular matrices.
   <LI>Singular Value Decomposition of rectangular matrices.
   <LI>Eigenvalue Decomposition of both symmetric and nonsymmetric square matrices.
</UL>
<DL>
<DT><B>Example of use:</B></DT>
<P>
<DD>Solve a linear system A x = b and compute the residual norm, ||b - A x||.
<P><PRE>
      double[][] vals = {{1.,2.,3},{4.,5.,6.},{7.,8.,10.}};
      Matrix A = new Matrix(vals);
      Matrix b = Matrix.random(3,1);
      Matrix x = A.solve(b);
      Matrix r = A.times(x).minus(b);
      double rnorm = r.normInf();
</PRE></DD>
</DL>

@author The MathWorks, Inc. and the National Institute of Standards and Technology.
@version 5 August 1998
*/
@ReceiverDependantMutable
public class Matrix implements Cloneable, Serializable {

/* ------------------------
   Class variables
 * ------------------------ */

   /** Array for internal storage of elements.
   @serial internal array storage.
   */
   private @Immutable double @Readonly [] @Mutable [] A;

   /** Row and column dimensions.
   @serial row dimension.
   @serial column dimension.
   */
   private int m, n;

/* ------------------------
   Constructors
 * ------------------------ */

   /** Construct an m-by-n matrix of zeros.
   @param m    Number of rows.
   @param n    Number of colums.
   */

   public @ReceiverDependantMutable Matrix (int m, int n) {
      this.m = m;
      this.n = n;
      A = new @Immutable double @Immutable [m][n];
   }

   /** Construct an m-by-n constant matrix.
   @param m    Number of rows.
   @param n    Number of colums.
   @param s    Fill the matrix with this scalar value.
   */

   public @Immutable Matrix (@Immutable int m, @Immutable int n, @Immutable double s) {
      this.m = m;
      this.n = n;
      A = new @Immutable double @Immutable [m][n];
      for (@Immutable int i = 0; i < m; i++) {
         for (@Immutable int j = 0; j < n; j++) {
            A[i][j] = s;
         }
      }
   }

   /** Construct a matrix from a 2-D array.
   @param A    Two-dimensional array of doubles.
   @exception  IllegalArgumentException All rows must have the same length
   @see        #constructWithCopy
   */

   public @ReceiverDependantMutable Matrix (@Immutable double @Mutable [] @Mutable [] A) {
      m = A.length;
      n = A[0].length;
      for (@Immutable int i = 0; i < m; i++) {
         if (A[i].length != n) {
            throw new @Mutable IllegalArgumentException("All rows must have the same length.");
         }
      }
      this.A = A;
   }

   /** Construct a matrix quickly without checking arguments.
   @param A    Two-dimensional array of doubles.
   @param m    Number of rows.
   @param n    Number of colums.
   */

   public @Immutable Matrix (@Immutable double @Readonly [] @Mutable [] A, @Immutable int m, @Immutable int n) {
      this.A = A;
      this.m = m;
      this.n = n;
   }

   /** Construct a matrix from a one-dimensional packed array
   @param vals One-dimensional array of doubles, packed by columns (ala Fortran).
   @param m    Number of rows.
   @exception  IllegalArgumentException Array length must be a multiple of m.
   */

   public @Immutable Matrix (@Immutable double vals @Immutable [], @Immutable int m) {
      this.m = m;
      n = (m != 0 ? vals.length/m : 0);
      if (m*n != vals.length) {
         throw new @Mutable IllegalArgumentException("Array length must be a multiple of m.");
      }
      A = new @Immutable double @Immutable [m][n];
      for (@Immutable int i = 0; i < m; i++) {
         for (@Immutable int j = 0; j < n; j++) {
            A[i][j] = vals[i+j*m];
         }
      }
   }

/* ------------------------
   Public Methods
 * ------------------------ */

   /** Construct a matrix from a copy of a 2-D array.
   @param A    Two-dimensional array of doubles.
   @exception  IllegalArgumentException All rows must have the same length
   */

   public static @Immutable Matrix constructWithCopy(@Immutable double @Readonly [] @Readonly [] A) {
      @Immutable
      int m = A.length;
      @Immutable
      int n = A[0].length;
      @Immutable
      Matrix X = new @Immutable Matrix(m,n);
      @Immutable
      double @Readonly [] @Mutable [] C = X.getArray();
      for (@Immutable int i = 0; i < m; i++) {
         if (A[i].length != n) {
            throw new @Mutable IllegalArgumentException
               ("All rows must have the same length.");
         }
         for (@Immutable int j = 0; j < n; j++) {
            C[i][j] = A[i][j];
         }
      }
      return X;
   }

   /** Make a deep copy of a matrix
   */

   public @ReceiverDependantMutable Matrix copy (@Readonly Matrix this) {
      @ReceiverDependantMutable
      Matrix X = new @ReceiverDependantMutable Matrix(m,n);
      @Immutable
      double @Readonly [] @Mutable [] C = X.getArray();
      for (@Immutable int i = 0; i < m; i++) {
         for (@Immutable int j = 0; j < n; j++) {
            C[i][j] = A[i][j];
         }
      }
      return X;
   }

   /** Clone the Matrix object.
   */

   public @ReceiverDependantMutable Object clone (@ReceiverDependantMutable Matrix this) {
      return this.copy();
   }

   /** Access the internal two-dimensional array.
   @return     Pointer to the two-dimensional array of matrix elements.
   */

   public double @Readonly [] @Mutable [] getArray (@Readonly Matrix this) {
      return A;
   }

   /** Copy the internal two-dimensional array.
   @return     Two-dimensional array copy of matrix elements.
   */

   public double @ReceiverDependantMutable [] @Mutable [] getArrayCopy (@Readonly Matrix this) {
      @Immutable
      double @ReceiverDependantMutable [] @Mutable [] C = new @Immutable double @ReceiverDependantMutable [m][n];
      for (@Immutable int i = 0; i < m; i++) {
         for (@Immutable int j = 0; j < n; j++) {
            C[i][j] = A[i][j];
         }
      }
      return C;
   }

   /** Make a one-dimensional column packed copy of the internal array.
   @return     Matrix elements packed in a one-dimensional array by columns.
   */

   public @Immutable double @Mutable [] getColumnPackedCopy (@Readonly Matrix this) {
      @Immutable
      double @Mutable [] vals = new double @Mutable [m*n];
      for (@Immutable int i = 0; i < m; i++) {
         for (@Immutable int j = 0; j < n; j++) {
            vals[i+j*m] = A[i][j];
         }
      }
      return vals;
   }

   /** Make a one-dimensional row packed copy of the internal array.
   @return     Matrix elements packed in a one-dimensional array by rows.
   */

   public @Immutable double @Mutable [] getRowPackedCopy (@Readonly Matrix this) {
      @Immutable
      double @Mutable [] vals = new double @Mutable [m*n];
      for (@Immutable int i = 0; i < m; i++) {
         for (@Immutable int j = 0; j < n; j++) {
            vals[i*n+j] = A[i][j];
         }
      }
      return vals;
   }

   /** Get row dimension.
   @return     m, the number of rows.
   */

   public int getRowDimension (@Readonly Matrix this) {
      return m;
   }

   /** Get column dimension.
   @return     n, the number of columns.
   */

   public int getColumnDimension (@Readonly Matrix this) {
      return n;
   }

   /** Get a single element.
   @param i    Row index.
   @param j    Column index.
   @return     A(i,j)
   @exception  ArrayIndexOutOfBoundsException
   */

   public @Immutable double get (@Readonly Matrix this, @Immutable int i, @Immutable int j) {
      return A[i][j];
   }

   /** Get a submatrix.
   @param i0   Initial row index
   @param i1   Final row index
   @param j0   Initial column index
   @param j1   Final column index
   @return     A(i0:i1,j0:j1)
   @exception  ArrayIndexOutOfBoundsException Submatrix indices
   */

   public @ReceiverDependantMutable Matrix getMatrix (@Immutable Matrix this, @Immutable int i0, @Immutable int i1, @Immutable int j0, @Immutable int j1) {
      @ReceiverDependantMutable
      Matrix X = new @ReceiverDependantMutable Matrix(i1-i0+1,j1-j0+1);
      @Immutable
      double @Readonly [] @Mutable [] B = X.getArray();
      try {
         for (@Immutable int i = i0; i <= i1; i++) {
            for (@Immutable int j = j0; j <= j1; j++) {
               B[i-i0][j-j0] = A[i][j];
            }
         }
      } catch(@Mutable ArrayIndexOutOfBoundsException e) {
         throw new @Mutable ArrayIndexOutOfBoundsException("Submatrix indices");
      }
      return X;
   }

   /** Get a submatrix.
   @param r    Array of row indices.
   @param c    Array of column indices.
   @return     A(r(:),c(:))
   @exception  ArrayIndexOutOfBoundsException Submatrix indices
   */

   public @Immutable Matrix getMatrix (@Readonly Matrix this, @Immutable int @Readonly [] r, @Immutable int @Readonly [] c) {
      @Immutable
      Matrix X = new @Immutable Matrix(r.length,c.length);
      @Immutable
      double @Readonly [] @Mutable [] B = X.getArray();
      try {
         for (@Immutable int i = 0; i < r.length; i++) {
            for (@Immutable int j = 0; j < c.length; j++) {
               B[i][j] = A[r[i]][c[j]];
            }
         }
      } catch(@Mutable ArrayIndexOutOfBoundsException e) {
         throw new @Mutable ArrayIndexOutOfBoundsException("Submatrix indices");
      }
      return X;
   }

   /** Get a submatrix.
   @param i0   Initial row index
   @param i1   Final row index
   @param c    Array of column indices.
   @return     A(i0:i1,c(:))
   @exception  ArrayIndexOutOfBoundsException Submatrix indices
   */

   public @Readonly Matrix getMatrix (@Readonly Matrix this, @Immutable int i0, @Immutable int i1, @Immutable int @Readonly [] c) {
      @Readonly
      Matrix X = new @Immutable Matrix(i1-i0+1,c.length);
      @Immutable
      double @Readonly [] @Mutable [] B = X.getArray();
      try {
         for (@Immutable int i = i0; i <= i1; i++) {
            for (@Immutable int j = 0; j < c.length; j++) {
               B[i-i0][j] = A[i][c[j]];
            }
         }
      } catch(@Mutable ArrayIndexOutOfBoundsException e) {
         throw new @Mutable ArrayIndexOutOfBoundsException("Submatrix indices");
      }
      return X;
   }

   /** Get a submatrix.
   @param r    Array of row indices.
   @param j0   Initial column index
   @param j1   Final column index
   @return     A(r(:),j0:j1)
   @exception  ArrayIndexOutOfBoundsException Submatrix indices
   */

   public @Immutable Matrix getMatrix (@Readonly Matrix this, int @Readonly [] r, int j0, int j1) {
      @Immutable
      Matrix X = new @Immutable Matrix(r.length,j1-j0+1);
      @Immutable
      double @Readonly [] @Mutable [] B = X.getArray();
      try {
         for (@Immutable int i = 0; i < r.length; i++) {
            for (@Immutable int j = j0; j <= j1; j++) {
               B[i][j-j0] = A[r[i]][j];
            }
         }
      } catch(@Mutable ArrayIndexOutOfBoundsException e) {
         throw new @Mutable ArrayIndexOutOfBoundsException("Submatrix indices");
      }
      return X;
   }

   /** Set a single element.
   @param i    Row index.
   @param j    Column index.
   @param s    A(i,j).
   @exception  ArrayIndexOutOfBoundsException
   */

   public void set (@Readonly Matrix this, @Immutable int i, @Immutable int j, @Immutable double s) {
      A[i][j] = s;
   }

   /** Set a submatrix.
   @param i0   Initial row index
   @param i1   Final row index
   @param j0   Initial column index
   @param j1   Final column index
   @param X    A(i0:i1,j0:j1)
   @exception  ArrayIndexOutOfBoundsException Submatrix indices
   */

   public void setMatrix (@Readonly Matrix this, @Immutable int i0, @Immutable int i1, @Immutable int j0, @Immutable int j1, @Readonly Matrix X) {
      try {
         for (@Immutable int i = i0; i <= i1; i++) {
            for (@Immutable int j = j0; j <= j1; j++) {
               A[i][j] = X.get(i-i0,j-j0);
            }
         }
      } catch(@Mutable ArrayIndexOutOfBoundsException e) {
         throw new @Mutable ArrayIndexOutOfBoundsException("Submatrix indices");
      }
   }

   /** Set a submatrix.
   @param r    Array of row indices.
   @param c    Array of column indices.
   @param X    A(r(:),c(:))
   @exception  ArrayIndexOutOfBoundsException Submatrix indices
   */

   public void setMatrix (@Readonly Matrix this, @Immutable int @Readonly [] r, @Immutable int @Readonly [] c, @Readonly Matrix X) {
      try {
         for (@Immutable int i = 0; i < r.length; i++) {
            for (@Immutable int j = 0; j < c.length; j++) {
               A[r[i]][c[j]] = X.get(i,j);
            }
         }
      } catch(@Mutable ArrayIndexOutOfBoundsException e) {
         throw new @Mutable ArrayIndexOutOfBoundsException("Submatrix indices");
      }
   }

   /** Set a submatrix.
   @param r    Array of row indices.
   @param j0   Initial column index
   @param j1   Final column index
   @param X    A(r(:),j0:j1)
   @exception  ArrayIndexOutOfBoundsException Submatrix indices
   */

   public void setMatrix (@Readonly Matrix this, @Immutable int @Readonly [] r, @Immutable int j0, @Immutable int j1, @Readonly Matrix X) {
      try {
         for (@Immutable int i = 0; i < r.length; i++) {
            for (@Immutable int j = j0; j <= j1; j++) {
               A[r[i]][j] = X.get(i,j-j0);
            }
         }
      } catch(@Mutable ArrayIndexOutOfBoundsException e) {
         throw new @Mutable ArrayIndexOutOfBoundsException("Submatrix indices");
      }
   }

   /** Set a submatrix.
   @param i0   Initial row index
   @param i1   Final row index
   @param c    Array of column indices.
   @param X    A(i0:i1,c(:))
   @exception  ArrayIndexOutOfBoundsException Submatrix indices
   */

   public void setMatrix (@Readonly Matrix this, @Immutable int i0, @Immutable int i1, @Immutable int @Readonly [] c, @Readonly Matrix X) {
      try {
         for (@Immutable int i = i0; i <= i1; i++) {
            for (@Immutable int j = 0; j < c.length; j++) {
               A[i][c[j]] = X.get(i-i0,j);
            }
         }
      } catch(@Mutable ArrayIndexOutOfBoundsException e) {
         throw new @Mutable ArrayIndexOutOfBoundsException("Submatrix indices");
      }
   }

   /** Matrix transpose.
   @return    A'
   */

   public @Immutable Matrix transpose (@Readonly Matrix this) {
      @Immutable
      Matrix X = new @Immutable Matrix(n,m);
      @Immutable
      double @Readonly [] @Mutable [] C = X.getArray();
      for (@Immutable int i = 0; i < m; i++) {
         for (@Immutable int j = 0; j < n; j++) {
            C[j][i] = A[i][j];
         }
      }
      return X;
   }

   /** One norm
   @return    maximum column sum.
   */

   public @Immutable double norm1 (@Readonly Matrix this) {
      @Immutable
      double f = 0;
      for (@Immutable int j = 0; j < n; j++) {
         @Immutable
         double s = 0;
         for (@Immutable int i = 0; i < m; i++) {
            s += Math.abs(A[i][j]);
         }
         f = Math.max(f,s);
      }
      return f;
   }

   /** Two norm
   @return    maximum singular value.
   */

   public @Immutable double norm2 (@Readonly Matrix this) {
      return (new @Mutable SingularValueDecomposition(this).norm2());
   }

   /** Infinity norm
   @return    maximum row sum.
   */

   public @Immutable double normInf (@Readonly Matrix this) {
      @Immutable
      double f = 0;
      for (@Immutable int i = 0; i < m; i++) {
         @Immutable
         double s = 0;
         for (@Immutable int j = 0; j < n; j++) {
            s += Math.abs(A[i][j]);
         }
         f = Math.max(f,s);
      }
      return f;
   }

   /** Frobenius norm
   @return    sqrt of sum of squares of all elements.
   */

   public @Immutable double normF (@Readonly Matrix this) {
      @Immutable
      double f = 0;
      for (@Immutable int i = 0; i < m; i++) {
         for (@Immutable int j = 0; j < n; j++) {
            f = Maths.hypot(f,A[i][j]);
         }
      }
      return f;
   }

   /**  Unary minus
   @return    -A
   */

   public @Immutable Matrix uminus (@Readonly Matrix this) {
      @Immutable
      Matrix X = new @Immutable Matrix(m,n);
      @Immutable
      double @Readonly [] @Mutable [] C = X.getArray();
      for (@Immutable int i = 0; i < m; i++) {
         for (@Immutable int j = 0; j < n; j++) {
            C[i][j] = -A[i][j];
         }
      }
      return X;
   }

   /** C = A + B
   @param B    another matrix
   @return     A + B
   */

   public @Immutable Matrix plus (@Readonly Matrix this, @Readonly Matrix B) {
      checkMatrixDimensions(B);
      @Immutable
      Matrix X = new @Immutable Matrix(m,n);
      @Immutable
      double @Readonly [] @Mutable [] C = X.getArray();
      for (@Immutable int i = 0; i < m; i++) {
         for (@Immutable int j = 0; j < n; j++) {
            C[i][j] = A[i][j] + B.A[i][j];
         }
      }
      return X;
   }

   /** A = A + B
   @param B    another matrix
   @return     A + B
   */

   public @Readonly Matrix plusEquals (@Readonly Matrix this, @Readonly Matrix B) {
      checkMatrixDimensions(B);
      for (@Immutable int i = 0; i < m; i++) {
         for (@Immutable int j = 0; j < n; j++) {
            A[i][j] = A[i][j] + B.A[i][j];
         }
      }
      return this;
   }

   /** C = A - B
   @param B    another matrix
   @return     A - B
   */

   public @Immutable Matrix minus (@Readonly Matrix this, @Readonly Matrix B) {
      checkMatrixDimensions(B);
      @Immutable
      Matrix X = new @Immutable Matrix(m,n);
      @Immutable
      double @Readonly [] @Mutable [] C = X.getArray();
      for (@Immutable int i = 0; i < m; i++) {
         for (@Immutable int j = 0; j < n; j++) {
            C[i][j] = A[i][j] - B.A[i][j];
         }
      }
      return X;
   }

   /** A = A - B
   @param B    another matrix
   @return     A - B
   */

   public @Readonly Matrix minusEquals (@Readonly Matrix this, @Readonly Matrix B) {
      checkMatrixDimensions(B);
      for (@Immutable int i = 0; i < m; i++) {
         for (@Immutable int j = 0; j < n; j++) {
            A[i][j] = A[i][j] - B.A[i][j];
         }
      }
      return this;
   }

   /** Element-by-element multiplication, C = A.*B
   @param B    another matrix
   @return     A.*B
   */

   public @Immutable Matrix arrayTimes (@Readonly Matrix this, @Readonly Matrix B) {
      checkMatrixDimensions(B);
      @Immutable
      Matrix X = new @Immutable Matrix(m,n);
      @Immutable
      double @Readonly [] @Mutable [] C = X.getArray();
      for (@Immutable int i = 0; i < m; i++) {
         for (@Immutable int j = 0; j < n; j++) {
            C[i][j] = A[i][j] * B.A[i][j];
         }
      }
      return X;
   }

   /** Element-by-element multiplication in place, A = A.*B
   @param B    another matrix
   @return     A.*B
   */

   public @Readonly Matrix arrayTimesEquals (@Readonly Matrix this, @Readonly Matrix B) {
      checkMatrixDimensions(B);
      for (@Immutable int i = 0; i < m; i++) {
         for (@Immutable int j = 0; j < n; j++) {
            A[i][j] = A[i][j] * B.A[i][j];
         }
      }
      return this;
   }

   /** Element-by-element right division, C = A./B
   @param B    another matrix
   @return     A./B
   */

   public @Immutable Matrix arrayRightDivide (@Readonly Matrix this, @Readonly Matrix B) {
      checkMatrixDimensions(B);
      @Immutable
      Matrix X = new @Immutable Matrix(m,n);
      @Immutable
      double @Readonly [] @Mutable [] C = X.getArray();
      for (@Immutable int i = 0; i < m; i++) {
         for (@Immutable int j = 0; j < n; j++) {
            C[i][j] = A[i][j] / B.A[i][j];
         }
      }
      return X;
   }

   /** Element-by-element right division in place, A = A./B
   @param B    another matrix
   @return     A./B
   */

   public @Readonly Matrix arrayRightDivideEquals (@Readonly Matrix this, @Readonly Matrix B) {
      checkMatrixDimensions(B);
      for (@Immutable int i = 0; i < m; i++) {
         for (@Immutable int j = 0; j < n; j++) {
            A[i][j] = A[i][j] / B.A[i][j];
         }
      }
      return this;
   }

   /** Element-by-element left division, C = A.\B
   @param B    another matrix
   @return     A.\B
   */

   public @Immutable Matrix arrayLeftDivide (@Readonly Matrix this, @Readonly Matrix B) {
      checkMatrixDimensions(B);
      @Immutable
      Matrix X = new @Immutable Matrix(m,n);
      @Immutable
      double @Readonly [] @Mutable [] C = X.getArray();
      for (@Immutable int i = 0; i < m; i++) {
         for (@Immutable int j = 0; j < n; j++) {
            C[i][j] = B.A[i][j] / A[i][j];
         }
      }
      return X;
   }

   /** Element-by-element left division in place, A = A.\B
   @param B    another matrix
   @return     A.\B
   */

   public @Readonly Matrix arrayLeftDivideEquals (@Readonly Matrix this, @Readonly Matrix B) {
      checkMatrixDimensions(B);
      for (@Immutable int i = 0; i < m; i++) {
         for (@Immutable int j = 0; j < n; j++) {
            A[i][j] = B.A[i][j] / A[i][j];
         }
      }
      return this;
   }

   /** Multiply a matrix by a scalar, C = s*A
   @param s    scalar
   @return     s*A
   */

   public @Immutable Matrix times (@Readonly Matrix this, @Immutable double s) {
      @Immutable
      Matrix X = new @Immutable Matrix(m,n);
      @Immutable
      double @Readonly [] @Mutable [] C = X.getArray();
      for (@Immutable int i = 0; i < m; i++) {
         for (@Immutable int j = 0; j < n; j++) {
            C[i][j] = s*A[i][j];
         }
      }
      return X;
   }

   /** Multiply a matrix by a scalar in place, A = s*A
   @param s    scalar
   @return     replace A by s*A
   */

   public @Readonly Matrix timesEquals (@Readonly Matrix this, @Immutable double s) {
      for (@Immutable int i = 0; i < m; i++) {
         for (@Immutable int j = 0; j < n; j++) {
            A[i][j] = s*A[i][j];
         }
      }
      return this;
   }

   /** Linear algebraic matrix multiplication, A * B
   @param B    another matrix
   @return     Matrix product, A * B
   @exception  IllegalArgumentException Matrix inner dimensions must agree.
   */

   public @Immutable Matrix times (@Readonly Matrix this, @Readonly Matrix B) {
      if (B.m != n) {
         throw new @Mutable IllegalArgumentException("Matrix inner dimensions must agree.");
      }
      @Immutable
      Matrix X = new @Immutable Matrix(m,B.n);
      @Immutable
      double @Readonly [] @Mutable [] C = X.getArray();
      @Immutable
      double @Mutable [] Bcolj = new double @Mutable [n];
      for (@Immutable int j = 0; j < B.n; j++) {
         for (@Immutable int k = 0; k < n; k++) {
            Bcolj[k] = B.A[k][j];
         }
         for (@Immutable int i = 0; i < m; i++) {
            @Immutable
            double @Readonly [] Arowi = A[i];
            @Immutable
            double s = 0;
            for (@Immutable int k = 0; k < n; k++) {
               s += Arowi[k]*Bcolj[k];
            }
            C[i][j] = s;
         }
      }
      return X;
   }

   /** LU Decomposition
   @return     LUDecomposition
   @see LUDecomposition
   */

   public @Mutable LUDecomposition lu (@Immutable Matrix this) {
      return new @Mutable LUDecomposition(this);
   }

   /** QR Decomposition
   @return     QRDecomposition
   @see QRDecomposition
   */

   public @Immutable QRDecomposition qr (@Immutable Matrix this) {
      return new @Immutable QRDecomposition(this);
   }

   /** Cholesky Decomposition
   @return     CholeskyDecomposition
   @see CholeskyDecomposition
   */

   public @Mutable CholeskyDecomposition chol (@Readonly Matrix this) {
      return new @Mutable CholeskyDecomposition(this);
   }

   /** Singular Value Decomposition
   @return     SingularValueDecomposition
   @see SingularValueDecomposition
   */

   public @Mutable SingularValueDecomposition svd (@Readonly Matrix this) {
      return new @Mutable SingularValueDecomposition(this);
   }

   /** Eigenvalue Decomposition
   @return     EigenvalueDecomposition
   @see EigenvalueDecomposition
   */

   public @Mutable EigenvalueDecomposition eig (@Readonly Matrix this) {
      return new @Mutable EigenvalueDecomposition(this);
   }

   /** Solve A*X = B
   @param B    right hand side
   @return     solution if A is square, least squares solution otherwise
   */

   public @Immutable Matrix solve (@Immutable Matrix this, @Readonly Matrix B) {
      return (m == n ? (new @Mutable LUDecomposition(this)).solve(B) :
                       (new @Immutable QRDecomposition(this)).solve(B));
   }

   /** Solve X*A = B, which is also A'*X' = B'
   @param B    right hand side
   @return     solution if A is square, least squares solution otherwise.
   */

   public @Immutable Matrix solveTranspose (@Readonly Matrix this, @Readonly Matrix B) {
      return transpose().solve(B.transpose());
   }

   /** Matrix inverse or pseudoinverse
   @return     inverse(A) if A is square, pseudoinverse otherwise.
   */

   public @Immutable Matrix inverse (@Immutable Matrix this) {
      return solve(identity(m,m));
   }

   /** Matrix determinant
   @return     determinant
   */

   public @Immutable double det (@Immutable Matrix this) {
      return new @Mutable LUDecomposition(this).det();
   }

   /** Matrix rank
   @return     effective numerical rank, obtained from SVD.
   */

   public @Immutable int rank (@Readonly Matrix this) {
      return new @Mutable SingularValueDecomposition(this).rank();
   }

   /** Matrix condition (2 norm)
   @return     ratio of largest to smallest singular value.
   */

   public @Immutable double cond (@Readonly Matrix this) {
      return new @Mutable SingularValueDecomposition(this).cond();
   }

   /** Matrix trace.
   @return     sum of the diagonal elements.
   */

   public @Immutable double trace (@Readonly Matrix this) {
      @Immutable
      double t = 0;
      for (@Immutable int i = 0; i < Math.min(m,n); i++) {
         t += A[i][i];
      }
      return t;
   }

   /** Generate matrix with random elements
   @param m    Number of rows.
   @param n    Number of colums.
   @return     An m-by-n matrix with uniformly distributed random elements.
   */

   public static @Immutable Matrix random (@Immutable int m, @Immutable int n) {
      @Immutable
      Matrix A = new @Immutable Matrix(m,n);
      @Immutable
      double @Readonly [] @Mutable [] X = A.getArray();
      for (@Immutable int i = 0; i < m; i++) {
         for (@Immutable int j = 0; j < n; j++) {
            X[i][j] = Math.random();
         }
      }
      return A;
   }

   /** Generate identity matrix
   @param m    Number of rows.
   @param n    Number of colums.
   @return     An m-by-n matrix with ones on the diagonal and zeros elsewhere.
   */

   public static @Immutable Matrix identity (@Immutable int m, @Immutable int n) {
      @Immutable
      Matrix A = new @Immutable Matrix(m,n);
      @Immutable
      double @Readonly [] @Mutable [] X = A.getArray();
      for (@Immutable int i = 0; i < m; i++) {
         for (@Immutable int j = 0; j < n; j++) {
            X[i][j] = (i == j ? 1.0 : 0.0);
         }
      }
      return A;
   }


   /** Print the matrix to stdout.   Line the elements up in columns
     * with a Fortran-like 'Fw.d' style format.
   @param w    Column width.
   @param d    Number of digits after the decimal.
   */

   public void print (@Readonly Matrix this, @Immutable int w, @Immutable int d) {
      print(new @Mutable PrintWriter(System.out,true),w,d); }

   /** Print the matrix to the output stream.   Line the elements up in
     * columns with a Fortran-like 'Fw.d' style format.
   @param output Output stream.
   @param w      Column width.
   @param d      Number of digits after the decimal.
   */

   public void print (@Readonly Matrix this, @Mutable PrintWriter output, @Immutable int w, @Immutable int d) {
      @Mutable
      DecimalFormat format = new @Mutable DecimalFormat();
      format.setDecimalFormatSymbols(new @Mutable DecimalFormatSymbols(Locale.US));
      format.setMinimumIntegerDigits(1);
      format.setMaximumFractionDigits(d);
      format.setMinimumFractionDigits(d);
      format.setGroupingUsed(false);
      print(output,format,w+2);
   }

   /** Print the matrix to stdout.  Line the elements up in columns.
     * Use the format object, and right justify within columns of width
     * characters.
     * Note that is the matrix is to be read back in, you probably will want
     * to use a NumberFormat that is set to US Locale.
   @param format A  Formatting object for individual elements.
   @param width     Field width for each column.
   @see java.text.DecimalFormat#setDecimalFormatSymbols
   */

   public void print (@Readonly Matrix this, @Mutable NumberFormat format, @Immutable int width) {
      print(new @Mutable PrintWriter(System.out,true),format,width); }

   // DecimalFormat is a little disappointing coming from Fortran or C's printf.
   // Since it doesn't pad on the left, the elements will come out different
   // widths.  Consequently, we'll pass the desired column width in as an
   // argument and do the extra padding ourselves.

   /** Print the matrix to the output stream.  Line the elements up in columns.
     * Use the format object, and right justify within columns of width
     * characters.
     * Note that is the matrix is to be read back in, you probably will want
     * to use a NumberFormat that is set to US Locale.
   @param output the output stream.
   @param format A formatting object to format the matrix elements
   @param width  Column width.
   @see java.text.DecimalFormat#setDecimalFormatSymbols
   */

   public void print (@Readonly Matrix this, @Mutable PrintWriter output, @Mutable NumberFormat format, @Immutable int width) {
      output.println();  // start on new line.
      for (@Immutable int i = 0; i < m; i++) {
         for (@Immutable int j = 0; j < n; j++) {
            @Immutable
            String s = format.format(A[i][j]); // format the number
            @Immutable
            int padding = Math.max(1,width-s.length()); // At _least_ 1 space
            for (@Immutable int k = 0; k < padding; k++)
               output.print(' ');
            output.print(s);
         }
         output.println();
      }
      output.println();   // end with blank line.
   }

   /** Read a matrix from a stream.  The format is the same the print method,
     * so printed matrices can be read back in (provided they were printed using
     * US Locale).  Elements are separated by
     * whitespace, all the elements for each row appear on a single line,
     * the last row is followed by a blank line.
   @param input the input stream.
   */

   public static @Mutable Matrix read (@Mutable BufferedReader input) throws java.io.IOException {
      @Mutable
      StreamTokenizer tokenizer= new @Mutable StreamTokenizer(input);

      // Although StreamTokenizer will parse numbers, it doesn't recognize
      // scientific notation (E or D); however, Double.valueOf does.
      // The strategy here is to disable StreamTokenizer's number parsing.
      // We'll only get whitespace delimited words, EOL's and EOF's.
      // These words should all be numbers, for Double.valueOf to parse.

      tokenizer.resetSyntax();
      tokenizer.wordChars(0,255);
      tokenizer.whitespaceChars(0, ' ');
      tokenizer.eolIsSignificant(true);
      java.util.@Mutable Vector<@Immutable Double> vD = new java.util.@Mutable Vector<@Immutable Double>();

      // Ignore initial empty lines
      while (tokenizer.nextToken() == StreamTokenizer.TT_EOL);
      if (tokenizer.ttype == StreamTokenizer.TT_EOF)
	throw new java.io.@Mutable IOException("Unexpected EOF on matrix read.");
      do {
         vD.addElement(Double.valueOf(tokenizer.sval)); // Read & store 1st row.
      } while (tokenizer.nextToken() == StreamTokenizer.TT_WORD);

      @Immutable
      int n = vD.size();  // Now we've got the number of columns!
      @Immutable
      double row @Mutable [] = new double @Mutable [n];
      for (@Immutable int j=0; j<n; j++)  // extract the elements of the 1st row.
         row[j]=vD.elementAt(j).doubleValue();
      java.util.@Mutable Vector<@Immutable double @Mutable []> v = new java.util.@Mutable Vector<@Immutable double @Mutable []>();
      v.addElement(row);  // Start storing rows instead of columns.
      while (tokenizer.nextToken() == StreamTokenizer.TT_WORD) {
         // While non-empty lines
         v.addElement(row = new double @Mutable [n]);
         @Immutable
         int j = 0;
         do {
            if (j >= n) throw new java.io.@Mutable IOException
               ("Row " + v.size() + " is too long.");
            row[j++] = Double.valueOf(tokenizer.sval).doubleValue();
         } while (tokenizer.nextToken() == StreamTokenizer.TT_WORD);
         if (j < n) throw new java.io.@Mutable IOException
            ("Row " + v.size() + " is too short.");
      }
      @Immutable
      int m = v.size();  // Now we've got the number of rows.
      @Immutable
      double @Mutable [] @Mutable [] A = new @Immutable double @Mutable [m][];
      v.copyInto(A);  // copy the rows out of the vector
      return new @Mutable Matrix(A);
   }


/* ------------------------
   Private Methods
 * ------------------------ */

   /** Check if size(A) == size(B) **/

   private void checkMatrixDimensions (@Readonly Matrix this, @Readonly Matrix B) {
      if (B.m != m || B.n != n) {
         throw new @Mutable IllegalArgumentException("Matrix dimensions must agree.");
      }
   }

  private static final @Immutable long serialVersionUID = 1;
}
