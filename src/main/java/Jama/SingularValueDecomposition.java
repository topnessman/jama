package Jama;
import qual.Immutable;
import qual.Readonly;
import qual.Mutable;
import qual.ReceiverDependantMutable;
import Jama.util.*;

   /** Singular Value Decomposition.
   <P>
   For an m-by-n matrix A with m >= n, the singular value decomposition is
   an m-by-n orthogonal matrix U, an n-by-n diagonal matrix S, and
   an n-by-n orthogonal matrix V so that A = U*S*V'.
   <P>
   The singular values, sigma[k] = S[k][k], are ordered so that
   sigma[0] >= sigma[1] >= ... >= sigma[n-1].
   <P>
   The singular value decompostion always exists, so the constructor will
   never fail.  The matrix condition number and the effective numerical
   rank can be computed from this decomposition.
   */

@ReceiverDependantMutable
public class SingularValueDecomposition implements java.io.Serializable {

/* ------------------------
   Class variables
 * ------------------------ */

   /** Arrays for internal storage of U and V.
   @serial internal storage of U.
   @serial internal storage of V.
   */
   private double @Immutable[][] U, V;

   /** Array for internal storage of singular values.
   @serial internal storage of singular values.
   */
   private @Immutable double @ReceiverDependantMutable [] s;

   /** Row and column dimensions.
   @serial row dimension.
   @serial column dimension.
   */
   private int m, n;

/* ------------------------
   Constructor
 * ------------------------ */

   /** Construct the singular value decomposition
       Structure to access U, S and V.
   @param Arg    Rectangular matrix
   */

   public @Mutable SingularValueDecomposition (@Readonly Matrix Arg) {

      // Derived from LINPACK code.
      // Initialize.
      @Immutable
      double @Readonly [] @Mutable [] A = Arg.getArrayCopy();
      m = Arg.getRowDimension();
      n = Arg.getColumnDimension();

      /* Apparently the failing cases are only a proper subset of (m<n),
	 so let's not throw error.  Correct fix to come later?
      if (m<n) {
	  throw new IllegalArgumentException("Jama SVD only works for m >= n"); }
      */
      @Immutable
      int nu = Math.min(m,n);
      s = new double @Mutable [Math.min(m+1,n)];
      U = new double @Immutable [m][nu];
      V = new double @Immutable [n][n];
      @Immutable
      double @Mutable [] e = new double @Mutable [n];
      @Immutable
      double @Mutable [] work = new double @Mutable [m];
      @Immutable
      boolean wantu = true;
      @Immutable
      boolean wantv = true;

      // Reduce A to bidiagonal form, storing the diagonal elements
      // in s and the super-diagonal elements in e.

      @Immutable
      int nct = Math.min(m-1,n);
      @Immutable
      int nrt = Math.max(0,Math.min(n-2,m));
      for (@Immutable int k = 0; k < Math.max(nct,nrt); k++) {
         if (k < nct) {

            // Compute the transformation for the k-th column and
            // place the k-th diagonal in s[k].
            // Compute 2-norm of k-th column without under/overflow.
            s[k] = 0;
            for (@Immutable int i = k; i < m; i++) {
               s[k] = Maths.hypot(s[k],A[i][k]);
            }
            if (s[k] != 0.0) {
               if (A[k][k] < 0.0) {
                  s[k] = -s[k];
               }
               for (@Immutable int i = k; i < m; i++) {
                  A[i][k] /= s[k];
               }
               A[k][k] += 1.0;
            }
            s[k] = -s[k];
         }
         for (@Immutable int j = k+1; j < n; j++) {
            if ((k < nct) & (s[k] != 0.0))  {

            // Apply the transformation.

               @Immutable
               double t = 0;
               for (@Immutable int i = k; i < m; i++) {
                  t += A[i][k]*A[i][j];
               }
               t = -t/A[k][k];
               for (@Immutable int i = k; i < m; i++) {
                  A[i][j] += t*A[i][k];
               }
            }

            // Place the k-th row of A into e for the
            // subsequent calculation of the row transformation.

            e[j] = A[k][j];
         }
         if (wantu & (k < nct)) {

            // Place the transformation in U for subsequent back
            // multiplication.

            for (@Immutable int i = k; i < m; i++) {
               U[i][k] = A[i][k];
            }
         }
         if (k < nrt) {

            // Compute the k-th row transformation and place the
            // k-th super-diagonal in e[k].
            // Compute 2-norm without under/overflow.
            e[k] = 0;
            for (@Immutable int i = k+1; i < n; i++) {
               e[k] = Maths.hypot(e[k],e[i]);
            }
            if (e[k] != 0.0) {
               if (e[k+1] < 0.0) {
                  e[k] = -e[k];
               }
               for (@Immutable int i = k+1; i < n; i++) {
                  e[i] /= e[k];
               }
               e[k+1] += 1.0;
            }
            e[k] = -e[k];
            if ((k+1 < m) & (e[k] != 0.0)) {

            // Apply the transformation.

               for (@Immutable int i = k+1; i < m; i++) {
                  work[i] = 0.0;
               }
               for (@Immutable int j = k+1; j < n; j++) {
                  for (@Immutable int i = k+1; i < m; i++) {
                     work[i] += e[j]*A[i][j];
                  }
               }
               for (@Immutable int j = k+1; j < n; j++) {
                  @Immutable
                  double t = -e[j]/e[k+1];
                  for (@Immutable int i = k+1; i < m; i++) {
                     A[i][j] += t*work[i];
                  }
               }
            }
            if (wantv) {

            // Place the transformation in V for subsequent
            // back multiplication.

               for (@Immutable int i = k+1; i < n; i++) {
                  V[i][k] = e[i];
               }
            }
         }
      }

      // Set up the final bidiagonal matrix or order p.

      @Immutable
      int p = Math.min(n,m+1);
      if (nct < n) {
         s[nct] = A[nct][nct];
      }
      if (m < p) {
         s[p-1] = 0.0;
      }
      if (nrt+1 < p) {
         e[nrt] = A[nrt][p-1];
      }
      e[p-1] = 0.0;

      // If required, generate U.

      if (wantu) {
         for (@Immutable int j = nct; j < nu; j++) {
            for (@Immutable int i = 0; i < m; i++) {
               U[i][j] = 0.0;
            }
            U[j][j] = 1.0;
         }
         for (@Immutable int k = nct-1; k >= 0; k--) {
            if (s[k] != 0.0) {
               for (@Immutable int j = k+1; j < nu; j++) {
                  @Immutable
                  double t = 0;
                  for (@Immutable int i = k; i < m; i++) {
                     t += U[i][k]*U[i][j];
                  }
                  t = -t/U[k][k];
                  for (@Immutable int i = k; i < m; i++) {
                     U[i][j] += t*U[i][k];
                  }
               }
               for (@Immutable int i = k; i < m; i++ ) {
                  U[i][k] = -U[i][k];
               }
               U[k][k] = 1.0 + U[k][k];
               for (@Immutable int i = 0; i < k-1; i++) {
                  U[i][k] = 0.0;
               }
            } else {
               for (@Immutable int i = 0; i < m; i++) {
                  U[i][k] = 0.0;
               }
               U[k][k] = 1.0;
            }
         }
      }

      // If required, generate V.

      if (wantv) {
         for (@Immutable int k = n-1; k >= 0; k--) {
            if ((k < nrt) & (e[k] != 0.0)) {
               for (@Immutable int j = k+1; j < nu; j++) {
                  @Immutable
                  double t = 0;
                  for (@Immutable int i = k+1; i < n; i++) {
                     t += V[i][k]*V[i][j];
                  }
                  t = -t/V[k+1][k];
                  for (@Immutable int i = k+1; i < n; i++) {
                     V[i][j] += t*V[i][k];
                  }
               }
            }
            for (@Immutable int i = 0; i < n; i++) {
               V[i][k] = 0.0;
            }
            V[k][k] = 1.0;
         }
      }

      // Main iteration loop for the singular values.

      @Immutable
      int pp = p-1;
      @Immutable
      int iter = 0;
      @Immutable
      double eps = Math.pow(2.0,-52.0);
      @Immutable
      double tiny = Math.pow(2.0,-966.0);
      while (p > 0) {
         int k,kase;

         // Here is where a test for too many iterations would go.

         // This section of the program inspects for
         // negligible elements in the s and e arrays.  On
         // completion the variables kase and k are set as follows.

         // kase = 1     if s(p) and e[k-1] are negligible and k<p
         // kase = 2     if s(k) is negligible and k<p
         // kase = 3     if e[k-1] is negligible, k<p, and
         //              s(k), ..., s(p) are not negligible (qr step).
         // kase = 4     if e(p-1) is negligible (convergence).

         for (k = p-2; k >= -1; k--) {
            if (k == -1) {
               break;
            }
            if (Math.abs(e[k]) <=
                  tiny + eps*(Math.abs(s[k]) + Math.abs(s[k+1]))) {
               e[k] = 0.0;
               break;
            }
         }
         if (k == p-2) {
            kase = 4;
         } else {
            @Immutable
            int ks;
            for (ks = p-1; ks >= k; ks--) {
               if (ks == k) {
                  break;
               }
               @Immutable
               double t = (ks != p ? Math.abs(e[ks]) : 0.) +
                          (ks != k+1 ? Math.abs(e[ks-1]) : 0.);
               if (Math.abs(s[ks]) <= tiny + eps*t)  {
                  s[ks] = 0.0;
                  break;
               }
            }
            if (ks == k) {
               kase = 3;
            } else if (ks == p-1) {
               kase = 1;
            } else {
               kase = 2;
               k = ks;
            }
         }
         k++;

         // Perform the task indicated by kase.

         switch (kase) {

            // Deflate negligible s(p).

            case 1: {
               @Immutable
               double f = e[p-2];
               e[p-2] = 0.0;
               for (@Immutable int j = p-2; j >= k; j--) {
                  @Immutable
                  double t = Maths.hypot(s[j],f);
                  @Immutable
                  double cs = s[j]/t;
                  @Immutable
                  double sn = f/t;
                  s[j] = t;
                  if (j != k) {
                     f = -sn*e[j-1];
                     e[j-1] = cs*e[j-1];
                  }
                  if (wantv) {
                     for (@Immutable int i = 0; i < n; i++) {
                        t = cs*V[i][j] + sn*V[i][p-1];
                        V[i][p-1] = -sn*V[i][j] + cs*V[i][p-1];
                        V[i][j] = t;
                     }
                  }
               }
            }
            break;

            // Split at negligible s(k).

            case 2: {
               @Immutable
               double f = e[k-1];
               e[k-1] = 0.0;
               for (@Immutable int j = k; j < p; j++) {
                  @Immutable
                  double t = Maths.hypot(s[j],f);
                  @Immutable
                  double cs = s[j]/t;
                  @Immutable
                  double sn = f/t;
                  s[j] = t;
                  f = -sn*e[j];
                  e[j] = cs*e[j];
                  if (wantu) {
                     for (@Immutable int i = 0; i < m; i++) {
                        t = cs*U[i][j] + sn*U[i][k-1];
                        U[i][k-1] = -sn*U[i][j] + cs*U[i][k-1];
                        U[i][j] = t;
                     }
                  }
               }
            }
            break;

            // Perform one qr step.

            case 3: {

               // Calculate the shift.

               @Immutable
               double scale = Math.max(Math.max(Math.max(Math.max(
                       Math.abs(s[p-1]),Math.abs(s[p-2])),Math.abs(e[p-2])),
                       Math.abs(s[k])),Math.abs(e[k]));
               @Immutable
               double sp = s[p-1]/scale;
               @Immutable
               double spm1 = s[p-2]/scale;
               @Immutable
               double epm1 = e[p-2]/scale;
               @Immutable
               double sk = s[k]/scale;
               @Immutable
               double ek = e[k]/scale;
               @Immutable
               double b = ((spm1 + sp)*(spm1 - sp) + epm1*epm1)/2.0;
               @Immutable
               double c = (sp*epm1)*(sp*epm1);
               @Immutable
               double shift = 0.0;
               if ((b != 0.0) | (c != 0.0)) {
                  shift = Math.sqrt(b*b + c);
                  if (b < 0.0) {
                     shift = -shift;
                  }
                  shift = c/(b + shift);
               }
               @Immutable
               double f = (sk + sp)*(sk - sp) + shift;
               @Immutable
               double g = sk*ek;

               // Chase zeros.

               for (@Immutable int j = k; j < p-1; j++) {
                  @Immutable
                  double t = Maths.hypot(f,g);
                  @Immutable
                  double cs = f/t;
                  @Immutable
                  double sn = g/t;
                  if (j != k) {
                     e[j-1] = t;
                  }
                  f = cs*s[j] + sn*e[j];
                  e[j] = cs*e[j] - sn*s[j];
                  g = sn*s[j+1];
                  s[j+1] = cs*s[j+1];
                  if (wantv) {
                     for (@Immutable int i = 0; i < n; i++) {
                        t = cs*V[i][j] + sn*V[i][j+1];
                        V[i][j+1] = -sn*V[i][j] + cs*V[i][j+1];
                        V[i][j] = t;
                     }
                  }
                  t = Maths.hypot(f,g);
                  cs = f/t;
                  sn = g/t;
                  s[j] = t;
                  f = cs*e[j] + sn*s[j+1];
                  s[j+1] = -sn*e[j] + cs*s[j+1];
                  g = sn*e[j+1];
                  e[j+1] = cs*e[j+1];
                  if (wantu && (j < m-1)) {
                     for (@Immutable int i = 0; i < m; i++) {
                        t = cs*U[i][j] + sn*U[i][j+1];
                        U[i][j+1] = -sn*U[i][j] + cs*U[i][j+1];
                        U[i][j] = t;
                     }
                  }
               }
               e[p-2] = f;
               iter = iter + 1;
            }
            break;

            // Convergence.

            case 4: {

               // Make the singular values positive.

               if (s[k] <= 0.0) {
                  s[k] = (s[k] < 0.0 ? -s[k] : 0.0);
                  if (wantv) {
                     for (@Immutable int i = 0; i <= pp; i++) {
                        V[i][k] = -V[i][k];
                     }
                  }
               }

               // Order the singular values.

               while (k < pp) {
                  if (s[k] >= s[k+1]) {
                     break;
                  }
                  @Immutable
                  double t = s[k];
                  s[k] = s[k+1];
                  s[k+1] = t;
                  if (wantv && (k < n-1)) {
                     for (@Immutable int i = 0; i < n; i++) {
                        t = V[i][k+1]; V[i][k+1] = V[i][k]; V[i][k] = t;
                     }
                  }
                  if (wantu && (k < m-1)) {
                     for (@Immutable int i = 0; i < m; i++) {
                        t = U[i][k+1]; U[i][k+1] = U[i][k]; U[i][k] = t;
                     }
                  }
                  k++;
               }
               iter = 0;
               p--;
            }
            break;
         }
      }
   }

/* ------------------------
   Public Methods
 * ------------------------ */

   /** Return the left singular vectors
   @return     U
   */

   public @Immutable Matrix getU (@Mutable SingularValueDecomposition this) {
      return new @Immutable Matrix(U,m,Math.min(m+1,n));
   }

   /** Return the right singular vectors
   @return     V
   */

   public @Immutable Matrix getV (@Mutable SingularValueDecomposition this) {
      return new @Immutable Matrix(V,n,n);
   }

   /** Return the one-dimensional array of singular values
   @return     diagonal of S.
   */

   public @Immutable double @Readonly [] getSingularValues (@Readonly SingularValueDecomposition this) {
      return s;
   }

   /** Return the diagonal matrix of singular values
   @return     S
   */

   public @Immutable Matrix getS (@Readonly SingularValueDecomposition this) {
      @Immutable
      Matrix X = new @Immutable Matrix(n,n);
      @Immutable
      double @Readonly [] @Mutable [] S = X.getArray();
      for (@Immutable int i = 0; i < n; i++) {
         for (@Immutable int j = 0; j < n; j++) {
            S[i][j] = 0.0;
         }
         S[i][i] = this.s[i];
      }
      return X;
   }

   /** Two norm
   @return     max(S)
   */

   public double norm2 (@Readonly SingularValueDecomposition this) {
      return s[0];
   }

   /** Two norm condition number
   @return     max(S)/min(S)
   */

   public double cond (@Readonly SingularValueDecomposition this) {
      return s[0]/s[Math.min(m,n)-1];
   }

   /** Effective numerical matrix rank
   @return     Number of nonnegligible singular values.
   */

   public int rank (@Readonly SingularValueDecomposition this) {
      @Immutable
      double eps = Math.pow(2.0,-52.0);
      @Immutable
      double tol = Math.max(m,n)*s[0]*eps;
      @Immutable
      int r = 0;
      for (@Immutable int i = 0; i < s.length; i++) {
         if (s[i] > tol) {
            r++;
         }
      }
      return r;
   }
  private static final @Immutable long serialVersionUID = 1;
}
