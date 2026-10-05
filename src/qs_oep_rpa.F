!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Dense molecular direct RPA-OEP kernel, following PyOEP methods/rpaoep.py.
!>        All spatial orbitals are retained. The correlation potential has zero source charge,
!>        independently of any HOMO constraint on the exchange potential.
! **************************************************************************************************
MODULE qs_oep_rpa
   USE kinds,                           ONLY: dp
   USE mathlib,                         ONLY: diamat_all
   USE qs_ks_types,                     ONLY: qs_oep_type
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE
   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_oep_rpa'
   PUBLIC :: oep_rpa_solve
CONTAINS

! **************************************************************************************************
!> \brief Build the charge-constrained, coupling-filtered MO integrals.
!> \param basis Coulomb integrals and normalized metric preprocessing
!> \param c complete MO coefficients
!> \param nocc number of doubly occupied spatial orbitals
!> \param threshold cutoff on the unweighted occupied-virtual Gram eigenvalues
!> \param w retained auxiliary basis transformation
!> \param mo full MO integrals in the retained basis
! **************************************************************************************************
   SUBROUTINE rpa_filtered_mo(basis, c, nocc, threshold, w, mo)
      TYPE(qs_oep_type), INTENT(IN)                      :: basis
      REAL(dp), INTENT(IN)                               :: c(:, :)
      INTEGER, INTENT(IN)                                :: nocc
      REAL(dp), INTENT(IN)                               :: threshold
      REAL(dp), ALLOCATABLE, INTENT(OUT)                 :: w(:, :), mo(:, :, :)

      INTEGER                                            :: i, naux, nmo, nret, nvir, p
      REAL(dp), ALLOCATABLE                              :: bov(:, :), d(:, :), eig(:), gram(:, :), &
                                                            matrix(:, :), projector(:, :), q(:), &
                                                            raw(:, :, :), trans(:, :)

      naux = SIZE(basis%charge)
      nmo = SIZE(c, 2)
      nvir = nmo - nocc
      IF (naux < 2) CPABORT("RPA-OEP requires at least two auxiliary functions")
      ALLOCATE (raw(naux, nmo, nmo), bov(naux, nvir*nocc), projector(naux, naux), eig(naux))
      DO p = 1, naux
         matrix = MATMUL(TRANSPOSE(c), MATMUL(basis%three_center(:, :, p), c))
         raw(p, :, :) = matrix
      END DO
      DO i = 1, nocc
         bov(:, (i - 1)*nvir + 1:i*nvir) = raw(:, nocc + 1:, i)
      END DO
      q = basis%y2/SQRT(SUM(basis%y2**2))
      projector = -SPREAD(q, 1, naux)*SPREAD(q, 2, naux)
      DO i = 1, naux
         projector(i, i) = projector(i, i) + 1.0_dp
      END DO
      CALL diamat_all(projector, eig)
      trans = MATMUL(basis%w2, projector(:, 2:))
      d = MATMUL(TRANSPOSE(trans), bov)
      gram = MATMUL(d, TRANSPOSE(d))
      DEALLOCATE (eig)
      ALLOCATE (eig(naux - 1))
      CALL diamat_all(gram, eig)
      nret = COUNT(eig >= threshold)
      IF (nret == 0) CPABORT("RPA-OEP coupling cutoff removes every auxiliary direction")
      w = MATMUL(trans, gram(:, naux - nret:))
      mo = RESHAPE(MATMUL(TRANSPOSE(w), RESHAPE(raw, [naux, nmo*nmo])), [nret, nmo, nmo])
   END SUBROUTINE rpa_filtered_mo

! **************************************************************************************************
!> \brief Solve 4 X0 vc = -rhs in the charge-only OEP subspace.
!> \param oep potential auxiliary basis
!> \param ri response auxiliary basis
!> \param c ...
!> \param eps ...
!> \param nocc ...
!> \param threshold OEP coupling cutoff
!> \param threshold_ri RI coupling cutoff
!> \param nw imaginary-frequency quadrature size
!> \param x0 imaginary-frequency quadrature scale
!> \param energy RPA correlation energy
!> \param v correlation potential AO matrix
!> \param coefficients correlation source coefficients in the original OEP basis
!> \param wc retained potential space
!> \param wr retained response space
!> \param rhs negative correlation energy derivative in the potential space
!> \param residual linear solve residual
! **************************************************************************************************
   SUBROUTINE oep_rpa_solve(oep, ri, c, eps, nocc, threshold, threshold_ri, nw, x0, &
                            energy, v, coefficients, wc, wr, rhs, residual)
      TYPE(qs_oep_type), INTENT(IN)                      :: oep, ri
      REAL(dp), INTENT(IN)                               :: c(:, :), eps(:)
      INTEGER, INTENT(IN)                                :: nocc
      REAL(dp), INTENT(IN)                               :: threshold, threshold_ri
      INTEGER, INTENT(IN)                                :: nw
      REAL(dp), INTENT(IN)                               :: x0
      REAL(dp), INTENT(OUT)                              :: energy
      REAL(dp), ALLOCATABLE, INTENT(OUT)                 :: v(:, :), coefficients(:), wc(:, :), &
                                                            wr(:, :), rhs(:)
      REAL(dp), INTENT(OUT)                              :: residual

      INTEGER                                            :: a, i, info, j, nc, nmo, nvir, p
      INTEGER, ALLOCATABLE                               :: piv(:)
      REAL(dp)                                           :: gap
      REAL(dp), ALLOCATABLE                              :: d(:, :), m(:, :, :), n(:, :, :), sol(:), &
                                                            weights(:), work(:, :), x(:, :)

      EXTERNAL dgesv

      CALL rpa_filtered_mo(oep, c, nocc, threshold, wc, n)
      CALL rpa_filtered_mo(ri, c, nocc, threshold_ri, wr, m)
      nmo = SIZE(eps)
      nvir = nmo - nocc
      nc = SIZE(wc, 2)
      ALLOCATE (rhs(nc), d(nc, nvir*nocc), weights(nvir*nocc), piv(nc))
      CALL rpa_rhs(n, m, eps, nocc, nw, x0, rhs, energy)
      DO i = 1, nocc
         DO a = 1, nvir
            j = (i - 1)*nvir + a
            gap = MAX(eps(nocc + a) - eps(i), 1.0E-6_dp)
            weights(j) = -4.0_dp/gap
            d(:, j) = n(:, nocc + a, i)
         END DO
      END DO
      x = MATMUL(d*SPREAD(weights, 1, nc), TRANSPOSE(d))
      work = x
      sol = -rhs
      CALL dgesv(nc, 1, work, nc, piv, sol, nc, info)
      IF (info /= 0) CPABORT("RPA-OEP constrained linear system is singular")
      residual = MAXVAL(ABS(MATMUL(x, sol) + rhs))
      coefficients = MATMUL(wc, sol)
      ALLOCATE (v(SIZE(c, 1), SIZE(c, 1)))
      v = 0.0_dp
      DO p = 1, SIZE(coefficients)
         v = v + coefficients(p)*oep%three_center(:, :, p)
      END DO
   END SUBROUTINE oep_rpa_solve

! **************************************************************************************************
!> \brief Direct RPA energy and negative local-potential derivative, with fixed auxiliary spaces.
!>        Includes occupied/virtual orbital rotations and orbital-energy derivatives.
!> \param n filtered potential MO integrals (potential,MO,MO)
!> \param m filtered response MO integrals (response,MO,MO)
!> \param eps orbital energies
!> \param no occupied spatial orbitals
!> \param nw frequency points
!> \param x0 frequency scale
!> \param rhs negative energy derivative
!> \param energy correlation energy
! **************************************************************************************************
   SUBROUTINE rpa_rhs(n, m, eps, no, nw, x0, rhs, energy)
      REAL(dp), INTENT(IN)                               :: n(:, :, :), m(:, :, :), eps(:)
      INTEGER, INTENT(IN)                                :: no, nw
      REAL(dp), INTENT(IN)                               :: x0
      REAL(dp), INTENT(OUT)                              :: rhs(:), energy

      INTEGER                                            :: a, i, iw, j, nc, nm, nr, nv, s
      REAL(dp)                                           :: denom, ev, omega, pi, value, wt
      REAL(dp), ALLOCATABLE :: b(:, :), eig(:), gamma(:), gap(:), inva(:, :), invi(:, :), &
         jacobi(:, :), kmat(:, :), kn(:), lam(:), mocc(:, :), mvirt(:, :), nodes(:), r(:, :, :), &
         response(:, :), rhs_omega(:), rotated(:, :), rt(:, :, :), ta(:, :), ti(:, :), weights(:), &
         y(:)

      nm = SIZE(eps)
      nv = nm - no
      nr = SIZE(m, 1)
      nc = SIZE(n, 1)
      IF (eps(no + 1) <= eps(no)) CPABORT("RPA-OEP requires a positive occupied-virtual gap")
      pi = ACOS(-1.0_dp)
      ! Golub-Welsch Gauss-Legendre rule on [-1,1], followed by PyOEP's rational map.
      ALLOCATE (jacobi(nw, nw), nodes(nw), weights(nw))
      jacobi = 0.0_dp
      DO i = 1, nw - 1
         jacobi(i, i + 1) = REAL(i, dp)/SQRT(4.0_dp*REAL(i, dp)**2 - 1.0_dp)
         jacobi(i + 1, i) = jacobi(i, i + 1)
      END DO
      CALL diamat_all(jacobi, nodes)
      weights = 2.0_dp*jacobi(1, :)**2
      ALLOCATE (gap(nv*no), lam(nv*no), gamma(nv*no), b(nr, nv*no), eig(nr), kn(nr), &
                mvirt(nm, nr*nv), mocc(nm, nr*no), invi(nm, no), inva(nm, nv), &
                r(nr, nv, no), rt(nr, no, nv), rhs_omega(nc), y(nv*no))
      DO s = 1, nm
         mvirt(s, :) = RESHAPE(m(:, no + 1:, s), [nr*nv])
         mocc(s, :) = RESHAPE(m(:, :no, s), [nr*no])
      END DO
      invi = 0.0_dp
      inva = 0.0_dp
      DO i = 1, no
         DO s = 1, nm
            denom = eps(i) - eps(s)
            IF (ABS(denom) > 1.0E-6_dp) invi(s, i) = 1.0_dp/denom
         END DO
         DO a = 1, nv
            j = (i - 1)*nv + a
            gap(j) = MAX(eps(no + a) - eps(i), 1.0E-6_dp)
            b(:, j) = m(:, no + a, i)
         END DO
      END DO
      DO a = 1, nv
         DO s = 1, nm
            denom = eps(no + a) - eps(s)
            IF (ABS(denom) > 1.0E-6_dp) inva(s, a) = 1.0_dp/denom
         END DO
      END DO
      energy = 0.0_dp
      rhs = 0.0_dp
      DO iw = 1, nw
         omega = x0*(1.0_dp + nodes(iw))/(1.0_dp - nodes(iw))
         wt = weights(iw)*2.0_dp*x0/(1.0_dp - nodes(iw))**2/(2.0_dp*pi)
         lam = 4.0_dp*gap/(gap**2 + omega**2)
         gamma = -4.0_dp/(gap**2 + omega**2) + 8.0_dp*gap**2/(gap**2 + omega**2)**2
         response = MATMUL(b*SPREAD(lam, 1, nr), TRANSPOSE(b))
         CALL diamat_all(response, eig)
         IF (MINVAL(eig) < -1.0E-8_dp) CPABORT("RPA-OEP response is not positive semidefinite")
         eig = MAX(eig, 0.0_dp)
         DO i = 1, nr
            ev = eig(i)
            ! Avoid loss of significance in log(1+ev)-ev at high frequencies.
            IF (ev < 1.0E-4_dp) THEN
               value = ev**2*(-0.5_dp + ev*(1.0_dp/3.0_dp + ev*(-0.25_dp + ev/5.0_dp)))
            ELSE
               value = LOG(1.0_dp + ev) - ev
            END IF
            energy = energy + wt*value
         END DO
         kn = -eig/(1.0_dp + eig)
         rotated = MATMUL(TRANSPOSE(response), b)
         kmat = MATMUL(response*SPREAD(kn, 1, nr), TRANSPOSE(response))
         r = RESHAPE(MATMUL(kmat, -b*SPREAD(lam, 1, nr)), [nr, nv, no])
         DO a = 1, nv
            rt(:, :, a) = r(:, a, :)
         END DO
         ti = MATMUL(mvirt, RESHAPE(r, [nr*nv, no]))*invi
         ta = MATMUL(mocc, RESHAPE(rt, [nr*no, nv]))*inva
         rhs_omega = 0.0_dp
         DO i = 1, no
            DO s = 1, nm
               rhs_omega = rhs_omega + 2.0_dp*n(:, s, i)*ti(s, i)
            END DO
         END DO
         DO a = 1, nv
            DO s = 1, nm
               rhs_omega = rhs_omega + 2.0_dp*n(:, s, no + a)*ta(s, a)
            END DO
         END DO
         y = gamma*MATMUL(kn, rotated**2)
         DO i = 1, no
            DO a = 1, nv
               j = (i - 1)*nv + a
               rhs_omega = rhs_omega + y(j)*(n(:, no + a, no + a) - n(:, i, i))
            END DO
         END DO
         rhs = rhs + wt*rhs_omega
      END DO
   END SUBROUTINE rpa_rhs
END MODULE qs_oep_rpa
