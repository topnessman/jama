package Jama;
import qual.Immutable;
import qual.Mutable;
import qual.Readonly;
import qual.ReceiverDependantMutable;
import Jama.util.*;

/** QR Decomposition.
<P>
   For an m-by-n matrix A with m >= n, the QR decomposition is an m-by-n
   orthogonal matrix Q and an n-by-n upper triangular matrix R so that
   A = Q*R.
<P>
   The QR decompostion always exists, even if the matrix does not have
   full rank, so the constructor will never fail.  The primary use of the
   QR decomposition is in the least squares solution of nonsquare systems
   of simultaneous linear equations.  This will fail if isFullRank()
   returns false.
*/

@ReceiverDependantMutable
public class QRDecomposition implements java.io.Serializable {

/* ------------------------
   Class variables
 * ------------------------ */

   /** Array for internal storage of decomposition.
   @serial internal array storage.
   */
   private @Immutable double @ReceiverDependantMutable [] @Mutable [] QR;

   /** Row and column dimensions.
   @serial column dimension.
   @serial row dimension.
   */
   private int m, n;

   /** Array for internal storage of diagonal of R.
   @serial diagonal of R.
   */
   private @Immutable double @Mutable [] Rdiag;

/* ------------------------
   Constructor
 * ------------------------ */

   /** QR Decomposition, computed by Householder reflections.
       Structure to access R and the Householder vectors and compute Q.
   @param A    Rectangular matrix
   */

   public @ReceiverDependantMutable QRDecomposition (@ReceiverDependantMutable Matrix A) {
      // Initialize.
      QR = A.getArrayCopy();
      m = A.getRowDimension();
      n = A.getColumnDimension();
      Rdiag = new double @Mutable [n];

      // Main loop.
      for (@Immutable int k = 0; k < n; k++) {
         // Compute 2-norm of k-th column without under/overflow.
         @Immutable
         double nrm = 0;
         for (@Immutable int i = k; i < m; i++) {
            nrm = Maths.hypot(nrm,QR[i][k]);
         }

         if (nrm != 0.0) {
            // Form k-th Householder vector.
            if (QR[k][k] < 0) {
               nrm = -nrm;
            }
            for (@Immutable int i = k; i < m; i++) {
               QR[i][k] /= nrm;
            }
            QR[k][k] += 1.0;

            // Apply transformation to remaining columns.
            for (@Immutable int j = k+1; j < n; j++) {
               @Immutable
               double s = 0.0;
               for (@Immutable int i = k; i < m; i++) {
                  s += QR[i][k]*QR[i][j];
               }
               s = -s/QR[k][k];
               for (@Immutable int i = k; i < m; i++) {
                  QR[i][j] += s*QR[i][k];
               }
            }
         }
         Rdiag[k] = -nrm;
      }
   }

/* ------------------------
   Public Methods
 * ------------------------ */

   /** Is the matrix full rank?
   @return     true if R, and hence A, has full rank.
   */

   public @Immutable boolean isFullRank (@Readonly QRDecomposition this) {
      for (@Immutable int j = 0; j < n; j++) {
         if (Rdiag[j] == 0)
            return false;
      }
      return true;
   }

   /** Return the Householder vectors
   @return     Lower trapezoidal matrix whose columns define the reflections
   */

   public @Immutable Matrix getH (@Readonly QRDecomposition this) {
      @Immutable
      Matrix X = new @Immutable Matrix(m,n);
      @Immutable
      double @Readonly [] @Mutable [] H = X.getArray();
      for (@Immutable int i = 0; i < m; i++) {
         for (@Immutable int j = 0; j < n; j++) {
            if (i >= j) {
               H[i][j] = QR[i][j];
            } else {
               H[i][j] = 0.0;
            }
         }
      }
      return X;
   }

   /** Return the upper triangular factor
   @return     R
   */

   public @Immutable Matrix getR (@Readonly QRDecomposition this) {
      @Immutable
      Matrix X = new @Immutable Matrix(n,n);
      @Immutable
      double @Readonly [] @Mutable [] R = X.getArray();
      for (@Immutable int i = 0; i < n; i++) {
         for (@Immutable int j = 0; j < n; j++) {
            if (i < j) {
               R[i][j] = QR[i][j];
            } else if (i == j) {
               R[i][j] = Rdiag[i];
            } else {
               R[i][j] = 0.0;
            }
         }
      }
      return X;
   }

   /** Generate and return the (economy-sized) orthogonal factor
   @return     Q
   */

   public @Immutable Matrix getQ (@Readonly QRDecomposition this) {
      @Immutable
      Matrix X = new @Immutable Matrix(m,n);
      @Immutable
      double @Readonly [] @Mutable [] Q = X.getArray();
      for (@Immutable int k = n-1; k >= 0; k--) {
         for (@Immutable int i = 0; i < m; i++) {
            Q[i][k] = 0.0;
         }
         Q[k][k] = 1.0;
         for (@Immutable int j = k; j < n; j++) {
            if (QR[k][k] != 0) {
               @Immutable
               double s = 0.0;
               for (@Immutable int i = k; i < m; i++) {
                  s += QR[i][k]*Q[i][j];
               }
               s = -s/QR[k][k];
               for (@Immutable int i = k; i < m; i++) {
                  Q[i][j] += s*QR[i][k];
               }
            }
         }
      }
      return X;
   }

   /** Least squares solution of A*X = B
   @param B    A Matrix with as many rows as A and any number of columns.
   @return     X that minimizes the two norm of Q*R*X-B.
   @exception  IllegalArgumentException  Matrix row dimensions must agree.
   @exception  RuntimeException  Matrix is rank deficient.
   */

   public @Immutable Matrix solve (@Readonly QRDecomposition this, @Readonly Matrix B) {
      if (B.getRowDimension() != m) {
         throw new @Mutable IllegalArgumentException("Matrix row dimensions must agree.");
      }
      if (!this.isFullRank()) {
         throw new @Immutable RuntimeException("Matrix is rank deficient.");
      }

      // Copy right hand side
      @Immutable
      int nx = B.getColumnDimension();
      @Immutable
      double @Readonly [] @Mutable [] X = B.getArrayCopy();

      // Compute Y = transpose(Q)*B
      for (@Immutable int k = 0; k < n; k++) {
         for (@Immutable int j = 0; j < nx; j++) {
            @Immutable
            double s = 0.0;
            for (@Immutable int i = k; i < m; i++) {
               s += QR[i][k]*X[i][j];
            }
            s = -s/QR[k][k];
            for (@Immutable int i = k; i < m; i++) {
               X[i][j] += s*QR[i][k];
            }
         }
      }
      // Solve R*X = Y;
      for (@Immutable int k = n-1; k >= 0; k--) {
         for (@Immutable int j = 0; j < nx; j++) {
            X[k][j] /= Rdiag[k];
         }
         for (@Immutable int i = 0; i < k; i++) {
            for (@Immutable int j = 0; j < nx; j++) {
               X[i][j] -= X[k][j]*QR[i][k];
            }
         }
      }
      return (new @Immutable Matrix(X,n,nx).getMatrix(0,n-1,0,nx-1));
   }
  private static final @Immutable long serialVersionUID = 1;
}
