package Jama;
import qual.Immutable;
import qual.Mutable;
import qual.Readonly;
import qual.ReceiverDependantMutable;

   /** LU Decomposition.
   <P>
   For an m-by-n matrix A with m >= n, the LU decomposition is an m-by-n
   unit lower triangular matrix L, an n-by-n upper triangular matrix U,
   and a permutation vector piv of length m so that A(piv,:) = L*U.
   If m < n, then L is m-by-m and U is m-by-n.
   <P>
   The LU decompostion with pivoting always exists, even if the matrix is
   singular, so the constructor will never fail.  The primary use of the
   LU decomposition is in the solution of square systems of simultaneous
   linear equations.  This will fail if isNonsingular() returns false.
   */

@ReceiverDependantMutable
public class LUDecomposition implements java.io.Serializable {

/* ------------------------
   Class variables
 * ------------------------ */

   /** Array for internal storage of decomposition.
   @serial internal array storage.
   */
   private @Immutable double @Immutable [] @ReceiverDependantMutable [] LU;

   /** Row and column dimensions, and pivot sign.
   @serial column dimension.
   @serial row dimension.
   @serial pivot sign.
   */
   private int m, n, pivsign; 

   /** Internal storage of pivot vector.
   @serial pivot vector.
   */
   private @Immutable int @ReceiverDependantMutable [] piv;

/* ------------------------
   Constructor
 * ------------------------ */

   /** LU Decomposition
       Structure to access L, U and piv.
   @param  A Rectangular matrix
   */

   public @Mutable LUDecomposition (@Immutable Matrix A) {

   // Use a "left-looking", dot-product, Crout/Doolittle algorithm.

      LU = A.getArrayCopy();
      m = A.getRowDimension();
      n = A.getColumnDimension();
      piv = new int @Mutable [m];
      for (@Immutable int i = 0; i < m; i++) {
         piv[i] = i;
      }
      pivsign = 1;
      @Immutable
      double @Mutable [] LUrowi;
      @Immutable
      double @Mutable [] LUcolj = new double @Mutable [m];

      // Outer loop.

      for (@Immutable int j = 0; j < n; j++) {

         // Make a copy of the j-th column to localize references.

         for (@Immutable int i = 0; i < m; i++) {
            LUcolj[i] = LU[i][j];
         }

         // Apply previous transformations.

         for (@Immutable int i = 0; i < m; i++) {
            LUrowi = LU[i];

            // Most of the time is spent in the following dot product.

            @Immutable
            int kmax = Math.min(i,j);
            @Immutable
            double s = 0.0;
            for (@Immutable int k = 0; k < kmax; k++) {
               s += LUrowi[k]*LUcolj[k];
            }

            LUrowi[j] = LUcolj[i] -= s;
         }
   
         // Find pivot and exchange if necessary.

         @Immutable
         int p = j;
         for (@Immutable int i = j+1; i < m; i++) {
            if (Math.abs(LUcolj[i]) > Math.abs(LUcolj[p])) {
               p = i;
            }
         }
         if (p != j) {
            for (@Immutable int k = 0; k < n; k++) {
               @Immutable
               double t = LU[p][k]; LU[p][k] = LU[j][k]; LU[j][k] = t;
            }
            @Immutable
            int k = piv[p]; piv[p] = piv[j]; piv[j] = k;
            pivsign = -pivsign;
         }

         // Compute multipliers.
         
         if (j < m & LU[j][j] != 0.0) {
            for (@Immutable int i = j+1; i < m; i++) {
               LU[i][j] /= LU[j][j];
            }
         }
      }
   }

/* ------------------------
   Temporary, experimental code.
   ------------------------ *\

   \** LU Decomposition, computed by Gaussian elimination.
   <P>
   This constructor computes L and U with the "daxpy"-based elimination
   algorithm used in LINPACK and MATLAB.  In Java, we suspect the dot-product,
   Crout algorithm will be faster.  We have temporarily included this
   constructor until timing experiments confirm this suspicion.
   <P>
   @param  A             Rectangular matrix
   @param  linpackflag   Use Gaussian elimination.  Actual value ignored.
   @return               Structure to access L, U and piv.
   *\

   public LUDecomposition (Matrix A, int linpackflag) {
      // Initialize.
      LU = A.getArrayCopy();
      m = A.getRowDimension();
      n = A.getColumnDimension();
      piv = new int[m];
      for (int i = 0; i < m; i++) {
         piv[i] = i;
      }
      pivsign = 1;
      // Main loop.
      for (int k = 0; k < n; k++) {
         // Find pivot.
         int p = k;
         for (int i = k+1; i < m; i++) {
            if (Math.abs(LU[i][k]) > Math.abs(LU[p][k])) {
               p = i;
            }
         }
         // Exchange if necessary.
         if (p != k) {
            for (int j = 0; j < n; j++) {
               double t = LU[p][j]; LU[p][j] = LU[k][j]; LU[k][j] = t;
            }
            int t = piv[p]; piv[p] = piv[k]; piv[k] = t;
            pivsign = -pivsign;
         }
         // Compute multipliers and eliminate k-th column.
         if (LU[k][k] != 0.0) {
            for (int i = k+1; i < m; i++) {
               LU[i][k] /= LU[k][k];
               for (int j = k+1; j < n; j++) {
                  LU[i][j] -= LU[i][k]*LU[k][j];
               }
            }
         }
      }
   }

\* ------------------------
   End of temporary code.
 * ------------------------ */

/* ------------------------
   Public Methods
 * ------------------------ */

   /** Is the matrix nonsingular?
   @return     true if U, and hence A, is nonsingular.
   */

   public @Immutable boolean isNonsingular (@Readonly LUDecomposition this) {
      for (@Immutable int j = 0; j < n; j++) {
         if (LU[j][j] == 0)
            return false;
      }
      return true;
   }

   /** Return lower triangular factor
   @return     L
   */

   public @Mutable Matrix getL (@Readonly LUDecomposition this) {
      @Mutable
      Matrix X = new @Mutable Matrix(m,n);
      @Immutable
      double @Readonly [] @Mutable [] L = X.getArray();
      for (@Immutable int i = 0; i < m; i++) {
         for (@Immutable int j = 0; j < n; j++) {
            if (i > j) {
               L[i][j] = LU[i][j];
            } else if (i == j) {
               L[i][j] = 1.0;
            } else {
               L[i][j] = 0.0;
            }
         }
      }
      return X;
   }

   /** Return upper triangular factor
   @return     U
   */

   public @Readonly Matrix getU (@Readonly LUDecomposition this) {
      @ReceiverDependantMutable
      Matrix X = new @ReceiverDependantMutable Matrix(n,n);
      @Immutable
      double @Readonly [] @Mutable [] U = X.getArray();
      for (@Immutable int i = 0; i < n; i++) {
         for (@Immutable int j = 0; j < n; j++) {
            if (i <= j) {
               U[i][j] = LU[i][j];
            } else {
               U[i][j] = 0.0;
            }
         }
      }
      return X;
   }

   /** Return pivot permutation vector
   @return     piv
   */

   public @Immutable int @Readonly [] getPivot (@Readonly LUDecomposition this) {
      @Immutable
      int @Mutable [] p = new int @Mutable [m];
      for (@Immutable int i = 0; i < m; i++) {
         p[i] = piv[i];
      }
      return p;
   }

   /** Return pivot permutation vector as a one-dimensional double array
   @return     (double) piv
   */

   public @Immutable double @Readonly [] getDoublePivot (@Readonly LUDecomposition this) {
      @Immutable
      double @Mutable [] vals = new double @Mutable [m];
      for (@Immutable int i = 0; i < m; i++) {
         vals[i] = (double) piv[i];
      }
      return vals;
   }

   /** Determinant
   @return     det(A)
   @exception  IllegalArgumentException  Matrix must be square
   */

   public @Immutable double det (@Readonly LUDecomposition this) {
      if (m != n) {
         throw new @Mutable IllegalArgumentException("Matrix must be square.");
      }
      @Immutable
      double d = (double) pivsign;
      for (@Immutable int j = 0; j < n; j++) {
         d *= LU[j][j];
      }
      return d;
   }

   /** Solve A*X = B
   @param  B   A Matrix with as many rows as A and any number of columns.
   @return     X so that L*U*X = B(piv,:)
   @exception  IllegalArgumentException Matrix row dimensions must agree.
   @exception  RuntimeException  Matrix is singular.
   */

   public @Immutable Matrix solve (@Readonly LUDecomposition this, @Readonly Matrix B) {
      if (B.getRowDimension() != m) {
         throw new @Mutable IllegalArgumentException("Matrix row dimensions must agree.");
      }
      if (!this.isNonsingular()) {
         throw new @ReceiverDependantMutable RuntimeException("Matrix is singular.");
      }

      // Copy right hand side with pivoting
      @Immutable
      int nx = B.getColumnDimension();
      @Immutable
      Matrix Xmat = B.getMatrix(piv,0,nx-1);
      @Immutable
      double @Readonly [] @Mutable [] X = Xmat.getArray();

      // Solve L*Y = B(piv,:)
      for (@Immutable int k = 0; k < n; k++) {
         for (@Immutable int i = k+1; i < n; i++) {
            for (@Immutable int j = 0; j < nx; j++) {
               X[i][j] -= X[k][j]*LU[i][k];
            }
         }
      }
      // Solve U*X = Y;
      for (@Immutable int k = n-1; k >= 0; k--) {
         for (@Immutable int j = 0; j < nx; j++) {
            X[k][j] /= LU[k][k];
         }
         for (@Immutable int i = 0; i < k; i++) {
            for (@Immutable int j = 0; j < nx; j++) {
               X[i][j] -= X[k][j]*LU[i][k];
            }
         }
      }
      return Xmat;
   }
  private static final @Immutable long serialVersionUID = 1;
}
