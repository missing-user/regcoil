subroutine svd_scan

  use global_variables
  use stel_constants
  use stel_kinds

  implicit none

  integer :: iflag, tic, toc, countrate
  real(dp), dimension(:,:), allocatable :: this_current_potential
  real(dp), dimension(:), allocatable :: RHS, solution, singular_values, U_transpose_times_RHS
  real(dp), dimension(:), allocatable :: JDifference_x, JDifference_y, Jdifference_z
  real(dp), dimension(:,:), allocatable :: this_J2_times_N, svd_matrix, U, VT
  real(dp) :: factor_theta, factor_zeta, factor
  integer :: ialpha, itheta, izeta, index, n_singular_values

  ! Variables needed by LAPACK:
  integer :: INFO, LWORK, M, N, LDA, LDU, LDVT
  real(dp), dimension(:), allocatable :: WORK
  integer, dimension(:), allocatable :: IWORK
  character :: JOBZ


  deallocate(alpha)
  nalpha = num_basis_functions
  allocate(alpha(nalpha))
  alpha=0

  allocate(solution(num_basis_functions), stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'

  allocate(chi2_B(nalpha), stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'
  allocate(chi2_J(nalpha), stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'
  allocate(max_Bnormal(nalpha), stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'
  allocate(max_J(nalpha), stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'
  allocate(current_potential(ntheta_coil,nzeta_coil,nalpha), stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'
  allocate(single_valued_current_potential_thetazeta(ntheta_coil,nzeta_coil,nalpha), stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'
  allocate(this_current_potential(ntheta_coil,nzeta_coil), stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'
  allocate(single_valued_current_potential_mn(num_basis_functions,nalpha), stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'
  allocate(Bnormal_total(ntheta_plasma,nzeta_plasma,nalpha), stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'
  allocate(J2(ntheta_coil,nzeta_coil,nalpha), stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'
  allocate(JDifference_x(ntheta_coil*nzeta_coil), stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'
  allocate(JDifference_y(ntheta_coil*nzeta_coil), stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'
  allocate(JDifference_z(ntheta_coil*nzeta_coil), stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'
  allocate(this_J2_times_N(ntheta_coil,nzeta_coil), stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'

  allocate(svd_matrix(ntheta_plasma*nzeta_plasma, num_basis_functions), stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'
  allocate(RHS(ntheta_plasma*nzeta_plasma), stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'


  print *,"Beginning SVD."
  call system_clock(tic,countrate)

  if (num_basis_functions > ntheta_plasma*nzeta_plasma) then
     print *,"Error! The svd scan is only designed to work when num_basis_functions <= ntheta_plasma*nzeta_plasma"
     stop
  end if
  
  index=0
  do itheta=1,ntheta_plasma
     do izeta = 1,nzeta_plasma
        !index = index+1
        index = (izeta-1)*ntheta_plasma + itheta
        factor = sqrt(norm_normal_plasma(itheta,izeta))
        RHS(index) = (Bnormal_from_plasma_current(itheta,izeta) + Bnormal_from_net_coil_currents(itheta,izeta))*factor
        svd_matrix(index,:) = g(index,:) / factor
     end do
  end do

  !JOBZ='A'  ! Compute ALL of the singular vectors.
  JOBZ='S'  ! Compute only the first min(M,N) singular vectors, which is N in our case.
  M = ntheta_plasma*nzeta_plasma
  N = num_basis_functions
  LDA = M
  LDU = M
  LDVT = N
  ! This next formula comes from the LAPACK documentation at the end of the file.
  LWORK = max( 3*min(M,N) + max(max(M,N),7*min(M,N)), &
       3*min(M,N) + max(max(M,N),5*min(M,N)*min(M,N)+4*min(M,N)), &
       min(M,N)*(6+4*min(M,N))+max(M,N))
  
  allocate(WORK(LWORK),stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'
  allocate(IWORK(8*min(M,N)),stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'
  
  n_singular_values = min(M,N)
  allocate(singular_values(n_singular_values),stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'    
  allocate(U(M,N),stat=iflag) ! If all singular vectors were computed, U would be M*M. But here we only compute the first N singular vectors,
  ! so U is M*N.
  if (iflag .ne. 0) stop 'Allocation error!'
  allocate(VT(N,N),stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'
  
  ! Call LAPACK to do the SVD:
  ! Note that svd_matrix is destroyed by LAPACK!
  call DGESDD(JOBZ, M, N, svd_matrix, LDA, singular_values, U, LDU, &
       VT, LDVT, WORK, LWORK, IWORK, INFO)
  
  if (INFO==0) then
     print *,"SVD (DGESDD) successful."
     if (n_singular_values<5) then
        print *,"Singular values:",singular_values
     else
        print *,"First 5 singular values:",singular_values(1:5)
        print *,"Last 5 singular values:", &
             singular_values(n_singular_values-4:n_singular_values)
     end if
  else if (INFO>0) then
     print *,"Error in SVD (DGESDD): Did not converge."
  else
     print *,"Error in SVD (DGESDD): Argument",INFO," was invalid."
  end if
  
  call system_clock(toc)
  print *,"Took ",real(toc-tic)/countrate," sec."

  call system_clock(tic)
  allocate(U_transpose_times_RHS(num_basis_functions),stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'
  U_transpose_times_RHS = matmul(transpose(U),RHS)
  call system_clock(toc)
  print *,"matmul: ",real(toc-tic)/countrate," sec."

  factor_zeta  = net_poloidal_current_Amperes / twopi
  factor_theta = net_toroidal_current_Amperes / twopi

  solution = 0
  do ialpha = 1,nalpha
     print "(a,e10.3,a,i3,a,i3,a)"," Solving system for alpha=",alpha(ialpha)," (",ialpha," of ",nalpha,")"
     call system_clock(tic,countrate)

     ! Add the contribution from the ialpha-th singular vectors and singular value:
     solution = solution + VT(ialpha,:) * (1/singular_values(ialpha)) * U_transpose_times_RHS(ialpha)

     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     ! Now compute diagnostics
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

     single_valued_current_potential_mn(:,ialpha) = solution
     this_current_potential = reshape(matmul(basis_functions, solution), (/ ntheta_coil, nzeta_coil /)) ! Could use BLAS2 for this line for extra speed.
     single_valued_current_potential_thetazeta(:,:,ialpha) = this_current_potential
     do izeta = 1,nzeta_coil
        do itheta = 1,ntheta_coil
           this_current_potential(itheta,izeta) = this_current_potential(itheta,izeta) &
                + factor_zeta*zeta_coil(izeta) + factor_theta*theta_coil(itheta)
        end do
     end do
     current_potential(:,:,ialpha) = this_current_potential

     JDifference_x = d_x - matmul(f_x, solution)
     JDifference_y = d_y - matmul(f_y, solution)
     JDifference_z = d_z - matmul(f_z, solution)
     this_J2_times_N = reshape(JDifference_x*Jdifference_x + JDifference_y*JDifference_y + JDifference_z*JDifference_z, (/ ntheta_coil, nzeta_coil /)) &
          / norm_normal_coil
     chi2_J(ialpha) = nfp * dtheta_coil * dzeta_coil * sum(this_J2_times_N)
     J2(:,:,ialpha) = this_J2_times_N / norm_normal_coil

     Bnormal_total(:,:,ialpha) = (reshape(matmul(g,solution),(/ ntheta_plasma, nzeta_plasma /)) / norm_normal_plasma) &
          + Bnormal_from_plasma_current + Bnormal_from_net_coil_currents

     max_Bnormal(ialpha) = maxval(abs(Bnormal_total(:,:,ialpha)))
     max_J(ialpha) = sqrt(maxval(J2(:,:,ialpha)))

     chi2_B(ialpha) = nfp * dtheta_plasma * dzeta_plasma &
          * sum(Bnormal_total(:,:,ialpha) * Bnormal_total(:,:,ialpha) * norm_normal_plasma)

     call system_clock(toc)
     print *,"  Diagnostics: ",real(toc-tic)/countrate," sec."
     print *,"  chi2_B:",chi2_B(ialpha),"chi2_J:",chi2_J(ialpha)
  end do

end subroutine svd_scan

    ! Here is the LAPACK documentation for the relevant SVD subroutine:

!!$*       SUBROUTINE DGESDD( JOBZ, M, N, A, LDA, S, U, LDU, VT, LDVT, WORK,
!!$*                          LWORK, IWORK, INFO )
!!$* 
!!$*       .. Scalar Arguments ..
!!$*       CHARACTER          JOBZ
!!$*       INTEGER            INFO, LDA, LDU, LDVT, LWORK, M, N
!!$*       ..
!!$*       .. Array Arguments ..
!!$*       INTEGER            IWORK( * )
!!$*       DOUBLE PRECISION   A( LDA, * ), S( * ), U( LDU, * ),
!!$*      $                   VT( LDVT, * ), WORK( * )
!!$*       ..
!!$*  
!!$*
!!$*> \par Purpose:
!!$*  =============
!!$*>
!!$*> \verbatim
!!$*>
!!$*> DGESDD computes the singular value decomposition (SVD) of a real
!!$*> M-by-N matrix A, optionally computing the left and right singular
!!$*> vectors.  If singular vectors are desired, it uses a
!!$*> divide-and-conquer algorithm.
!!$*>
!!$*> The SVD is written
!!$*>
!!$*>      A = U * SIGMA * transpose(V)
!!$*>
!!$*> where SIGMA is an M-by-N matrix which is zero except for its
!!$*> min(m,n) diagonal elements, U is an M-by-M orthogonal matrix, and
!!$*> V is an N-by-N orthogonal matrix.  The diagonal elements of SIGMA
!!$*> are the singular values of A; they are real and non-negative, and
!!$*> are returned in descending order.  The first min(m,n) columns of
!!$*> U and V are the left and right singular vectors of A.
!!$*>
!!$*> Note that the routine returns VT = V**T, not V.
!!$*>
!!$*> The divide and conquer algorithm makes very mild assumptions about
!!$*> floating point arithmetic. It will work on machines with a guard
!!$*> digit in add/subtract, or on those binary machines without guard
!!$*> digits which subtract like the Cray X-MP, Cray Y-MP, Cray C-90, or
!!$*> Cray-2. It could conceivably fail on hexadecimal or decimal machines
!!$*> without guard digits, but we know of none.
!!$*> \endverbatim
!!$*
!!$*  Arguments:
!!$*  ==========
!!$*
!!$*> \param[in] JOBZ
!!$*> \verbatim
!!$*>          JOBZ is CHARACTER*1
!!$*>          Specifies options for computing all or part of the matrix U:
!!$*>          = 'A':  all M columns of U and all N rows of V**T are
!!$*>                  returned in the arrays U and VT;
!!$*>          = 'S':  the first min(M,N) columns of U and the first
!!$*>                  min(M,N) rows of V**T are returned in the arrays U
!!$*>                  and VT;
!!$*>          = 'O':  If M >= N, the first N columns of U are overwritten
!!$*>                  on the array A and all rows of V**T are returned in
!!$*>                  the array VT;
!!$*>                  otherwise, all columns of U are returned in the
!!$*>                  array U and the first M rows of V**T are overwritten
!!$*>                  in the array A;
!!$*>          = 'N':  no columns of U or rows of V**T are computed.
!!$*> \endverbatim
!!$*>
!!$*> \param[in] M
!!$*> \verbatim
!!$*>          M is INTEGER
!!$*>          The number of rows of the input matrix A.  M >= 0.
!!$*> \endverbatim
!!$*>
!!$*> \param[in] N
!!$*> \verbatim
!!$*>          N is INTEGER
!!$*>          The number of columns of the input matrix A.  N >= 0.
!!$*> \endverbatim
!!$*>
!!$*> \param[in,out] A
!!$*> \verbatim
!!$*>          A is DOUBLE PRECISION array, dimension (LDA,N)
!!$*>          On entry, the M-by-N matrix A.
!!$*>          On exit,
!!$*>          if JOBZ = 'O',  A is overwritten with the first N columns
!!$*>                          of U (the left singular vectors, stored
!!$*>                          columnwise) if M >= N;
!!$*>                          A is overwritten with the first M rows
!!$*>                          of V**T (the right singular vectors, stored
!!$*>                          rowwise) otherwise.
!!$*>          if JOBZ .ne. 'O', the contents of A are destroyed.
!!$*> \endverbatim
!!$*>
!!$*> \param[in] LDA
!!$*> \verbatim
!!$*>          LDA is INTEGER
!!$*>          The leading dimension of the array A.  LDA >= max(1,M).
!!$*> \endverbatim
!!$*>
!!$*> \param[out] S
!!$*> \verbatim
!!$*>          S is DOUBLE PRECISION array, dimension (min(M,N))
!!$*>          The singular values of A, sorted so that S(i) >= S(i+1).
!!$*> \endverbatim
!!$*>
!!$*> \param[out] U
!!$*> \verbatim
!!$*>          U is DOUBLE PRECISION array, dimension (LDU,UCOL)
!!$*>          UCOL = M if JOBZ = 'A' or JOBZ = 'O' and M < N;
!!$*>          UCOL = min(M,N) if JOBZ = 'S'.
!!$*>          If JOBZ = 'A' or JOBZ = 'O' and M < N, U contains the M-by-M
!!$*>          orthogonal matrix U;
!!$*>          if JOBZ = 'S', U contains the first min(M,N) columns of U
!!$*>          (the left singular vectors, stored columnwise);
!!$*>          if JOBZ = 'O' and M >= N, or JOBZ = 'N', U is not referenced.
!!$*> \endverbatim
!!$*>
!!$*> \param[in] LDU
!!$*> \verbatim
!!$*>          LDU is INTEGER
!!$*>          The leading dimension of the array U.  LDU >= 1; if
!!$*>          JOBZ = 'S' or 'A' or JOBZ = 'O' and M < N, LDU >= M.
!!$*> \endverbatim
!!$*>
!!$*> \param[out] VT
!!$*> \verbatim
!!$*>          VT is DOUBLE PRECISION array, dimension (LDVT,N)
!!$*>          If JOBZ = 'A' or JOBZ = 'O' and M >= N, VT contains the
!!$*>          N-by-N orthogonal matrix V**T;
!!$*>          if JOBZ = 'S', VT contains the first min(M,N) rows of
!!$*>          V**T (the right singular vectors, stored rowwise);
!!$*>          if JOBZ = 'O' and M < N, or JOBZ = 'N', VT is not referenced.
!!$*> \endverbatim
!!$*>
!!$*> \param[in] LDVT
!!$*> \verbatim
!!$*>          LDVT is INTEGER
!!$*>          The leading dimension of the array VT.  LDVT >= 1; if
!!$*>          JOBZ = 'A' or JOBZ = 'O' and M >= N, LDVT >= N;
!!$*>          if JOBZ = 'S', LDVT >= min(M,N).
!!$*> \endverbatim
!!$*>
!!$*> \param[out] WORK
!!$*> \verbatim
!!$*>          WORK is DOUBLE PRECISION array, dimension (MAX(1,LWORK))
!!$*>          On exit, if INFO = 0, WORK(1) returns the optimal LWORK;
!!$*> \endverbatim
!!$*>
!!$*> \param[in] LWORK
!!$*> \verbatim
!!$*>          LWORK is INTEGER
!!$*>          The dimension of the array WORK. LWORK >= 1.
!!$*>          If JOBZ = 'N',
!!$*>            LWORK >= 3*min(M,N) + max(max(M,N),7*min(M,N)).
!!$*>          If JOBZ = 'O',
!!$*>            LWORK >= 3*min(M,N) + 
!!$*>                     max(max(M,N),5*min(M,N)*min(M,N)+4*min(M,N)).
!!$*>          If JOBZ = 'S' or 'A'
!!$*>            LWORK >= min(M,N)*(6+4*min(M,N))+max(M,N)
!!$*>          For good performance, LWORK should generally be larger.
!!$*>          If LWORK = -1 but other input arguments are legal, WORK(1)
!!$*>          returns the optimal LWORK.
!!$*> \endverbatim
!!$*>
!!$*> \param[out] IWORK
!!$*> \verbatim
!!$*>          IWORK is INTEGER array, dimension (8*min(M,N))
!!$*> \endverbatim
!!$*>
!!$*> \param[out] INFO
!!$*> \verbatim
!!$*>          INFO is INTEGER
!!$*>          = 0:  successful exit.
!!$*>          < 0:  if INFO = -i, the i-th argument had an illegal value.
!!$*>          > 0:  DBDSDC did not converge, updating process failed.
!!$*> \endverbatim
