
/* BiCGSTAB(L) algorithm for the n-by-n problem Ax = b */
int bicgstabL(const int           L,
              GlobalVector&       x,
              MATRIX&             A,
              void*               Adata,
              const GlobalVector& b,
              const double        tol,
              int                 maxiters,
              realnum*            work,
              const bool          quiet) {
  int n    = x.n();
  int ncomp= x.ncomp();

  vector<GlobalVector> r(L + 1, GlobalVector(ncomp, n));
  vector<GlobalVector> u(L + 1, GlobalVector(ncomp, n));
  assert(r[0].n() == n);
  assert(r[0].ncomp() == ncomp);

  double bnrm= b.norm();
  if (bnrm == 0.0)
    bnrm= 1.0;

  int    iter                 = 0;
  double last_output_wall_time= wall_time();

  vector<double> gamma(L + 1), gamma_p(L + 1), gamma_pp(L + 1), tau(L * L), sigma(L + 1);

  int          ierr    = 0;  // error code to return, if any
  const double breaktol= 1e-30;

  // rtilde = r[0] = b - Ax
  realnum* rtilde= work + (2 * L + 2) * n;
  r[0]           = b;
  vmult(A, r[0], x, -1);
  GlobalVector rtilde= r[0];

  { /* Sleipjen normalizes rtilde in his code; it seems to help slightly */
    double s= rtilde.norm();
    assert(s > 0);
    rtilde*= 1.0 / s;
  }

  u[0].zero();

  double rho= 1.0, alpha= 0, omega= 1;

  double resid;
  while ((resid= r[0].norm()) > tol * bnrm) {
    ++iter;
    cout << "residual " << d << " " << resid / bnrm << endl;

    rho= -omega * rho;
    for (int j= 0; j < L; ++j) {
      if (fabs(rho) < breaktol) {
        ierr= -1;
        goto finish;
      }
      double rho1= r[j] * rtilde;
      double beta= alpha * rho1 / rho;
      rho        = rho1;
      for (int i= 0; i <= j; ++i) {
        // for (int m = 0; m < n; ++m) u[i][m] = r[i][m] - beta * u[i][m];
        u[i]= sadd(-beta, 1.0, r[i]);
      }
      //   A(u[j], u[j+1], Adata);
      u[j + 1].zero();
      vmult(A, u[j + 1], u[j], 1.0);

      alpha= rho / (u[j + 1] * rtilde);
      for (int i= 0; i <= j; ++i)
        // xpay(n, r[i], -alpha, u[i+1]);
        r[i].add(-alpha, u[i + 1]);

      // A(r[j], r[j+1], Adata);
      r[j + 1].zero();
      vmult(A, r[j + 1], r[j], 1.0);
      // xpay(n, x, alpha, u[0]);
      x.add(alpha, u[0]);
    }

    for (int j= 1; j <= L; ++j) {
      for (int i= 1; i < j; ++i) {
        int ij = (j - 1) * L + (i - 1);
        tau[ij]= r[j] * r[i] / sigma[i];
        r[j].add(-tau[ij], r[i]);
      }
      sigma[j]  = r[j] * r[j];
      gamma_p[j]= r[0] * r[j] / sigma[j];
    }

    omega= gamma[L]= gamma_p[L];
    for (int j= L - 1; j >= 1; --j) {
      gamma[j]= gamma_p[j];
      for (int i= j + 1; i <= L; ++i)
        gamma[j]-= tau[(i - 1) * L + (j - 1)] * gamma[i];
    }
    for (int j= 1; j < L; ++j) {
      gamma_pp[j]= gamma[j + 1];
      for (int i= j + 1; i < L; ++i)
        gamma_pp[j]+= tau[(i - 1) * L + (j - 1)] * gamma[i + 1];
    }
    // xpay(n, x, gamma[1], r[0]);
    x.add(gamma[1], r[0]);
    // xpay(n, r[0], -gamma_p[L], r[L]);
    r[0].add(-gamma_p[L], r[L]);
    // xpay(n, u[0], -gamma[L], u[L]);
    u[0].add(-gamma[L], u[L]);

    for (int j= 1; j < L; ++j) { /* TODO: use blas DGEMV (for L > 2) */
      // xpay(n, x, gamma_pp[j], r[j]);
      x.add(gamma_pp[j], r[j]);
      // xpay(n, r[0], -gamma_p[j], r[j]);
      r[0].add(-gamma_p[j], r[j]);
      // xpay(n, u[0], -gamma[j], u[j]);
      u[0].add(-gamma[j], u[j]);
    }

    if (iter == maxiters) {
      ierr= 1;
      break;
    }
  }

  cout << "final residual = " << r[0].norm() / bnrm << endl;

finish:
  maxiter= iter;
  return ierr;
}
