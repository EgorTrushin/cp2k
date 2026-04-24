!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Export CC4S input files (EigenEnergies + CoulombVertex) from a
!>        converged RI-HFX Gamma-point calculation.
!>
!>  Assumptions:
!>    - Gamma-point only (real MO coefficients and real 3-center integrals)
!>    - RI-HFX with Coulomb metric (same_op=.TRUE., the default)
!>    - Closed-shell reference; spin channel 1 is used
!>
!>  The Coulomb vertex written to disk is:
!>    Gamma^P_{pq} = sum_{mu,nu}  C_{mu,p}  B^P_{mu,nu}  C_{nu,q}
!>  where the Cholesky-dressed 3-center integrals are
!>    B^P_{mu,nu} = sum_Q  M_{PQ}  (mu nu | Q)
!>  and
!>    M           = t_2c_int   =  J^{-1/2}  (Coulomb metric factor)
!>    (mu nu | Q) = t_3c_int_ctr_1  in  (AO, RI, AO) block-tensor format
!>
!>  Output files (CC4S text format, Hartree atomic units):
!>    EigenEnergies.yaml / .elements
!>    CoulombVertex.yaml / .elements
!>
!> \par History
!>       2026 created
! **************************************************************************************************

MODULE cc4s_export

   USE cp_dbcsr_api,                    ONLY: dbcsr_p_type
   USE cp_dbcsr_operations,             ONLY: copy_dbcsr_to_fm
   USE cp_files,                        ONLY: get_unit_number
   USE cp_fm_diag,                      ONLY: cp_fm_geeig
   USE cp_fm_struct,                    ONLY: cp_fm_struct_create, &
                                              cp_fm_struct_release, &
                                              cp_fm_struct_type
   USE cp_fm_types,                     ONLY: cp_fm_create, &
                                              cp_fm_get_submatrix, &
                                              cp_fm_release, &
                                              cp_fm_type
   USE cp_log_handling,                 ONLY: cp_get_default_logger, &
                                              cp_logger_get_default_io_unit, &
                                              cp_logger_type
   USE dbt_api,                         ONLY: dbt_get_block, &
                                              dbt_get_info, &
                                              dbt_iterator_blocks_left, &
                                              dbt_iterator_next_block, &
                                              dbt_iterator_start, &
                                              dbt_iterator_stop, &
                                              dbt_iterator_type, &
                                              dbt_nblks_total, &
                                              dbt_type
   USE cp_control_types,                ONLY: dft_control_type
   USE hfx_ri,                          ONLY: hfx_ri_pre_scf_mo
   USE hfx_types,                       ONLY: hfx_type
   USE input_section_types,             ONLY: section_vals_type, &
                                              section_vals_val_get
   USE kinds,                           ONLY: dp
   USE message_passing,                 ONLY: mp_para_env_type
   USE qs_environment_types,            ONLY: get_qs_env, &
                                              qs_environment_type
   USE qs_mo_types,                     ONLY: get_mo_set, &
                                              mo_set_type

#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   PUBLIC :: write_cc4s_coulomb_vertex

   CHARACTER(LEN=*), PARAMETER, PRIVATE :: moduleN = 'cc4s_export'

CONTAINS

! **************************************************************************************************
!> \brief  Main entry point: gather RI-HFX tensors, build Coulomb vertex, write CC4S files.
!> \param qs_env        QS environment after converged SCF
!> \param print_section &CC4S print-key section
! **************************************************************************************************
   SUBROUTINE write_cc4s_coulomb_vertex(qs_env, print_section)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(section_vals_type), POINTER                   :: print_section

      CHARACTER(LEN=*), PARAMETER :: routineN = 'write_cc4s_coulomb_vertex'

      INTEGER                                            :: handle, homo, iMO, iNU, &
                                                           nAO, nRI, nMO, nspins, output_unit
      INTEGER, DIMENSION(3)                              :: dims3
      INTEGER, DIMENSION(2)                              :: dims2
      REAL(KIND=dp)                                      :: fermi_energy
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: eigenvalues_buf, eri_flat
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: C_buf, metric_mat, B_2d, T_Pp
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: eri_raw, Gamma_full
      REAL(KIND=dp), DIMENSION(:), POINTER               :: eigenvalues_ptr
      TYPE(cp_fm_struct_type), POINTER                   :: ao_ao_struct
      TYPE(cp_fm_type)                                   :: fm_ks, fm_s, fm_eigvec, fm_work
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_ks, matrix_s
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(hfx_type), DIMENSION(:, :), POINTER           :: x_data
      TYPE(mo_set_type), DIMENSION(:), POINTER           :: mos
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(cp_fm_type), POINTER                          :: mo_coeff_fm
      TYPE(cp_logger_type), POINTER                      :: logger
      CHARACTER(LEN=512)                                 :: outdir
      CHARACTER(LEN=32)                                  :: cc4s_method

      EXTERNAL DGEMM, DAXPY

      CALL timeset(routineN, handle)

      logger => cp_get_default_logger()
      output_unit = cp_logger_get_default_io_unit(logger)

      ! ----------------------------------------------------------------
      ! Output directory and method
      ! ----------------------------------------------------------------
      CALL section_vals_val_get(print_section, "OUTPUT_DIR", c_val=outdir)
      CALL section_vals_val_get(print_section, "METHOD", c_val=cc4s_method)
      outdir = TRIM(ADJUSTL(outdir))
      IF (LEN_TRIM(outdir) == 0) outdir = '.'

      ! ----------------------------------------------------------------
      ! Get RI data and MO data from qs_env
      ! ----------------------------------------------------------------
      NULLIFY (x_data, mos, para_env, mo_coeff_fm, eigenvalues_ptr, dft_control, &
               matrix_ks, matrix_s, ao_ao_struct)
      CALL get_qs_env(qs_env, x_data=x_data, mos=mos, para_env=para_env, dft_control=dft_control)
      nspins = dft_control%nspins

      IF (.NOT. ASSOCIATED(x_data)) THEN
         CPWARN("CC4S export: no HFX data (x_data). Enable &HF / &RI in the input.")
         CALL timestop(handle)
         RETURN
      END IF
      IF (.NOT. ASSOCIATED(x_data(1, 1)%ri_data)) THEN
         CPWARN("CC4S export: RI data not initialized. Enable &HF %RI.")
         CALL timestop(handle)
         RETURN
      END IF

      IF (output_unit > 0) WRITE (output_unit, '(/,T3,A)') &
         "CC4S_EXPORT| Starting Coulomb vertex export"

      ! ----------------------------------------------------------------
      ! Rebuild raw 3-center integrals and J^{-1/2} metric into
      ! t_3c_int_ctr_1 and t_2c_int.  These may be empty or stale if the
      ! SCF ran in ri_pmat flavor (t_2c_int is cleared after each KS
      ! update in that path).  hfx_ri_pre_scf_mo recomputes them freshly.
      ! ----------------------------------------------------------------
      CALL hfx_ri_pre_scf_mo(qs_env, x_data(1, 1)%ri_data, nspins)

      ! ----------------------------------------------------------------
      ! MO coefficients and eigenvalues (all orbitals, spin channel 1).
      ! When nMO < nAO (no ADDED_MOS), diagonalize the converged KS matrix
      ! to obtain all nAO eigenpairs automatically.
      ! ----------------------------------------------------------------
      CALL get_mo_set(mos(1), mo_coeff=mo_coeff_fm, eigenvalues=eigenvalues_ptr, &
                      homo=homo, nmo=nMO, nao=nAO)

      IF (nMO < nAO) THEN
         IF (output_unit > 0) WRITE (output_unit, '(T3,A,2(I0,A))') &
            "CC4S_EXPORT| nMO (", nMO, ") < nAO (", nAO, &
            ") — diagonalizing KS matrix for full MO set"

         CALL get_qs_env(qs_env, matrix_ks=matrix_ks, matrix_s=matrix_s)

         ! Create nAO x nAO fm struct and workspace matrices
         CALL cp_fm_struct_create(ao_ao_struct, nrow_global=nAO, ncol_global=nAO, &
                                  template_fmstruct=mo_coeff_fm%matrix_struct)
         CALL cp_fm_create(fm_ks,    matrix_struct=ao_ao_struct, name="fm_ks")
         CALL cp_fm_create(fm_s,     matrix_struct=ao_ao_struct, name="fm_s")
         CALL cp_fm_create(fm_eigvec, matrix_struct=ao_ao_struct, name="fm_eigvec")
         CALL cp_fm_create(fm_work,  matrix_struct=ao_ao_struct, name="fm_work")
         CALL cp_fm_struct_release(ao_ao_struct)

         CALL copy_dbcsr_to_fm(matrix_ks(1)%matrix, fm_ks)
         CALL copy_dbcsr_to_fm(matrix_s(1)%matrix,  fm_s)

         ALLOCATE (eigenvalues_buf(nAO), C_buf(nAO, nAO))
         CALL cp_fm_geeig(fm_ks, fm_s, fm_eigvec, eigenvalues_buf, fm_work)
         C_buf(:, :) = 0.0_dp
         CALL cp_fm_get_submatrix(fm_eigvec, C_buf)

         CALL cp_fm_release(fm_ks)
         CALL cp_fm_release(fm_s)
         CALL cp_fm_release(fm_eigvec)
         CALL cp_fm_release(fm_work)

         nMO = nAO
      ELSE
         ALLOCATE (C_buf(nAO, nMO), eigenvalues_buf(nMO))
         C_buf(:, :) = 0.0_dp
         CALL cp_fm_get_submatrix(mo_coeff_fm, C_buf)
         eigenvalues_buf(1:nMO) = eigenvalues_ptr(1:nMO)
      END IF

      fermi_energy = eigenvalues_buf(homo)
      IF (homo < nMO) &
         fermi_energy = 0.5_dp*(eigenvalues_buf(homo) + eigenvalues_buf(homo + 1))

      ! ----------------------------------------------------------------
      ! Gather 3-center integrals t_3c_int_ctr_1 -> eri_raw(nRI, nAO, nAO)
      ! Source tensor index order: (AO, RI, AO) = (mu, Q, nu).
      ! We store the result transposed to (P, mu, nu) so that the metric
      ! contraction reduces to a single DGEMM over the RI index.
      ! ----------------------------------------------------------------
      CALL dbt_get_info(x_data(1, 1)%ri_data%t_3c_int_ctr_1(1, 1), nfull_total=dims3)
      ! dims3 = (nAO, nRI, nAO) matching the (AO, RI, AO) tensor
      nRI = dims3(2)
      CPASSERT(dims3(1) == nAO .AND. dims3(3) == nAO)

      IF (output_unit > 0) WRITE (output_unit, '(T3,A,3I6)') &
         "CC4S_EXPORT| nAO, nRI, nMO =", nAO, nRI, nMO

      ALLOCATE (eri_raw(nRI, nAO, nAO))
      CALL gather_3c_pao(x_data(1, 1)%ri_data%t_3c_int_ctr_1(1, 1), &
                         para_env, nRI, nAO, eri_raw)

      ! ----------------------------------------------------------------
      ! Gather 2-center metric t_2c_int(1,1) -> metric_mat(nRI, nRI)
      ! ----------------------------------------------------------------
      CALL dbt_get_info(x_data(1, 1)%ri_data%t_2c_int(1, 1), nfull_total=dims2)
      CPASSERT(dims2(1) == nRI .AND. dims2(2) == nRI)

      ALLOCATE (metric_mat(nRI, nRI))
      CALL gather_2c_dense(x_data(1, 1)%ri_data%t_2c_int(1, 1), &
                           para_env, nRI, metric_mat)

      ! ----------------------------------------------------------------
      ! Apply metric: B_2d(P, mu + nAO*(nu-1)) = sum_Q M(P,Q) * eri_raw(Q, mu, nu)
      ! i.e.  B_2d(nRI, nAO^2) = M(nRI, nRI) @ eri_flat(nRI, nAO^2)
      ! ----------------------------------------------------------------
      ALLOCATE (eri_flat(nRI*nAO*nAO))
      eri_flat(:) = RESHAPE(eri_raw, [nRI*nAO*nAO])
      DEALLOCATE (eri_raw)

      ALLOCATE (B_2d(nRI, nAO*nAO))
      CALL DGEMM('N', 'N', nRI, nAO*nAO, nRI, 1.0_dp, &
                 metric_mat(1, 1), nRI, &
                 eri_flat(1), nRI, &
                 0.0_dp, B_2d(1, 1), nRI)
      DEALLOCATE (eri_flat, metric_mat)

      ! ----------------------------------------------------------------
      ! MO transform:
      !   Gamma(P, p, q) = sum_{mu,nu} C(mu,p) * B_2d(P, mu+nAO*(nu-1)) * C(nu,q)
      !
      ! Outer loop over nu:
      !   T_Pp(P, p) = sum_mu B_2d(P, mu + nAO*(nu-1)) * C(mu, p)   [DGEMM]
      !   Gamma(:, :, q) += T_Pp * C(nu, q)                          [DAXPY]
      ! ----------------------------------------------------------------
      ALLOCATE (Gamma_full(nRI, nMO, nMO), T_Pp(nRI, nMO))
      Gamma_full(:, :, :) = 0.0_dp

      DO iNU = 1, nAO
         ! B_slice = B_2d(:, (iNU-1)*nAO+1 : iNU*nAO), contiguous nRI x nAO block
         CALL DGEMM('N', 'N', nRI, nMO, nAO, 1.0_dp, &
                    B_2d(1, (iNU - 1)*nAO + 1), nRI, &
                    C_buf(1, 1), nAO, &
                    0.0_dp, T_Pp(1, 1), nRI)
         ! Gamma(:, :, iMO) += T_Pp * C(iNU, iMO)
         DO iMO = 1, nMO
            CALL DAXPY(nRI*nMO, C_buf(iNU, iMO), T_Pp(1, 1), 1, &
                       Gamma_full(1, 1, iMO), 1)
         END DO
      END DO
      DEALLOCATE (B_2d, T_Pp, C_buf)

      IF (output_unit > 0) WRITE (output_unit, '(T3,A,3(I0,A))') &
         "CC4S_EXPORT| Coulomb vertex shape: (", nRI, ", ", nMO, ", ", nMO, ")"

      ! ----------------------------------------------------------------
      ! Write CC4S files (ionode only)
      ! ----------------------------------------------------------------
      IF (para_env%mepos == 0) THEN
         CALL EXECUTE_COMMAND_LINE('mkdir -p '//TRIM(outdir), wait=.TRUE.)
         CALL write_eigen_energies_files(outdir, nMO, homo, eigenvalues_buf, fermi_energy)
         CALL write_coulomb_vertex_files(outdir, nRI, nMO, Gamma_full)
         CALL write_cc4s_in(outdir, TRIM(cc4s_method))
         IF (output_unit > 0) WRITE (output_unit, '(T3,A,A)') &
            "CC4S_EXPORT| Files written to: ", TRIM(outdir)
      END IF

      DEALLOCATE (Gamma_full, eigenvalues_buf)

      CALL timestop(handle)

   END SUBROUTINE write_cc4s_coulomb_vertex

! **************************************************************************************************
!> \brief  Gather the (AO, RI, AO) = (mu, Q, nu) block tensor into a dense array
!>         with index order (P, mu, nu) so P varies fastest in Fortran storage.
!>         Allreduces across all MPI ranks so every rank holds the result.
! **************************************************************************************************
   SUBROUTINE gather_3c_pao(tensor, para_env, nRI, nAO, eri_pao)
      TYPE(dbt_type), INTENT(INOUT)                      :: tensor
      TYPE(mp_para_env_type), POINTER                    :: para_env
      INTEGER, INTENT(IN)                                :: nRI, nAO
      REAL(KIND=dp), DIMENSION(nRI, nAO, nAO), INTENT(OUT) :: eri_pao

      INTEGER                                            :: imu, inu, nb1, nb2, nb3
      INTEGER, DIMENSION(3)                              :: ind3
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: off1, off2, off3
      LOGICAL                                            :: found
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: blk3
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: flat
      TYPE(dbt_iterator_type)                            :: iter

      eri_pao(:, :, :) = 0.0_dp

      ! Allocate offset arrays sized by number of blocks in each dimension
      nb1 = dbt_nblks_total(tensor, 1)
      nb2 = dbt_nblks_total(tensor, 2)
      nb3 = dbt_nblks_total(tensor, 3)
      ALLOCATE (off1(nb1), off2(nb2), off3(nb3))
      CALL dbt_get_info(tensor, blk_offset_1=off1, blk_offset_2=off2, blk_offset_3=off3)

      ! DBT iterator requires an OpenMP parallel region; each thread processes its own shard
!$OMP PARALLEL DEFAULT(NONE) &
!$OMP SHARED(tensor, off1, off2, off3, eri_pao) &
!$OMP PRIVATE(iter, ind3, blk3, found, imu, inu)
      CALL dbt_iterator_start(iter, tensor)
      DO WHILE (dbt_iterator_blocks_left(iter))
         CALL dbt_iterator_next_block(iter, ind3)
         CALL dbt_get_block(tensor, ind3, blk3, found)
         IF (.NOT. found) CYCLE

         ! Source: blk3(mu_loc, Q_loc, nu_loc), tensor index order = (AO, RI, AO)
         ! Target: eri_pao(Q, mu, nu)  — transpose first two indices
         DO inu = 1, SIZE(blk3, 3)
            DO imu = 1, SIZE(blk3, 1)
               eri_pao(off2(ind3(2)):off2(ind3(2)) + SIZE(blk3, 2) - 1, &
                        off1(ind3(1)) + imu - 1, &
                        off3(ind3(3)) + inu - 1) = blk3(imu, :, inu)
            END DO
         END DO

         DEALLOCATE (blk3)
      END DO
      CALL dbt_iterator_stop(iter)
!$OMP END PARALLEL
      DEALLOCATE (off1, off2, off3)

      ! Allreduce: collect contributions from all MPI ranks
      ALLOCATE (flat(nRI*nAO*nAO))
      flat(:) = RESHAPE(eri_pao, [nRI*nAO*nAO])
      CALL para_env%sum(flat)
      eri_pao(:, :, :) = RESHAPE(flat, [nRI, nAO, nAO])
      DEALLOCATE (flat)

   END SUBROUTINE gather_3c_pao

! **************************************************************************************************
!> \brief  Gather a (RI, RI) block tensor to a dense matrix.  Allreduces across all MPI ranks.
! **************************************************************************************************
   SUBROUTINE gather_2c_dense(tensor, para_env, nRI, dense)
      TYPE(dbt_type), INTENT(INOUT)                      :: tensor
      TYPE(mp_para_env_type), POINTER                    :: para_env
      INTEGER, INTENT(IN)                                :: nRI
      REAL(KIND=dp), DIMENSION(nRI, nRI), INTENT(OUT)    :: dense

      INTEGER                                            :: nb1, nb2
      INTEGER, DIMENSION(2)                              :: ind2
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: off1, off2
      LOGICAL                                            :: found
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: blk2
      TYPE(dbt_iterator_type)                            :: iter

      dense(:, :) = 0.0_dp

      nb1 = dbt_nblks_total(tensor, 1)
      nb2 = dbt_nblks_total(tensor, 2)
      ALLOCATE (off1(nb1), off2(nb2))
      CALL dbt_get_info(tensor, blk_offset_1=off1, blk_offset_2=off2)

!$OMP PARALLEL DEFAULT(NONE) &
!$OMP SHARED(tensor, off1, off2, dense) &
!$OMP PRIVATE(iter, ind2, blk2, found)
      CALL dbt_iterator_start(iter, tensor)
      DO WHILE (dbt_iterator_blocks_left(iter))
         CALL dbt_iterator_next_block(iter, ind2)
         CALL dbt_get_block(tensor, ind2, blk2, found)
         IF (.NOT. found) CYCLE
         dense(off1(ind2(1)):off1(ind2(1)) + SIZE(blk2, 1) - 1, &
               off2(ind2(2)):off2(ind2(2)) + SIZE(blk2, 2) - 1) = blk2
         DEALLOCATE (blk2)
      END DO
      CALL dbt_iterator_stop(iter)
!$OMP END PARALLEL
      DEALLOCATE (off1, off2)

      CALL para_env%sum(dense)

   END SUBROUTINE gather_2c_dense

! **************************************************************************************************
!> \brief  Write EigenEnergies.yaml and EigenEnergies.elements (CC4S format, Hartree units).
! **************************************************************************************************
   SUBROUTINE write_eigen_energies_files(outdir, nMO, homo, eigenvalues, fermi_energy)
      CHARACTER(LEN=*), INTENT(IN)                       :: outdir
      INTEGER, INTENT(IN)                                :: nMO, homo
      REAL(KIND=dp), INTENT(IN)                          :: eigenvalues(nMO), fermi_energy

      INTEGER                                            :: iunit, i
      CHARACTER(LEN=512)                                 :: fname

      fname = TRIM(outdir)//'/EigenEnergies.yaml'
      iunit = get_unit_number(fname)
      OPEN (UNIT=iunit, FILE=TRIM(fname), STATUS='REPLACE', ACTION='WRITE')
      WRITE (iunit, '(A)')         'version: 100'
      WRITE (iunit, '(A)')         'type: Tensor'
      WRITE (iunit, '(A)')         'scalarType: Real64'
      WRITE (iunit, '(A)')         'dimensions:'
      WRITE (iunit, '(A,I0)')      '- length: ', nMO
      WRITE (iunit, '(A)')         '  type: State'
      WRITE (iunit, '(A)')         'elements:'
      WRITE (iunit, '(A)')         '  type: TextFile'
      WRITE (iunit, '(A)')         'unit: 1.0'
      WRITE (iunit, '(A)')         'metaData:'
      WRITE (iunit, '(A,ES25.16)') '  fermiEnergy: ', fermi_energy
      WRITE (iunit, '(A,I0)')      '  homoIndex: ', homo
      WRITE (iunit, '(A)')         '  energies:'
      DO i = 1, nMO
         WRITE (iunit, '(A,ES25.16)') '  - ', eigenvalues(i)
      END DO
      CLOSE (iunit)

      fname = TRIM(outdir)//'/EigenEnergies.elements'
      iunit = get_unit_number(fname)
      OPEN (UNIT=iunit, FILE=TRIM(fname), STATUS='REPLACE', ACTION='WRITE')
      DO i = 1, nMO
         WRITE (iunit, '(ES25.16)') eigenvalues(i)
      END DO
      CLOSE (iunit)

   END SUBROUTINE write_eigen_energies_files

! **************************************************************************************************
!> \brief  Write CoulombVertex.yaml and CoulombVertex.elements (CC4S format, Hartree units).
!>         Elements are written in Fortran (column-major) order: P (AuxiliaryField) varies fastest.
! **************************************************************************************************
   SUBROUTINE write_coulomb_vertex_files(outdir, nRI, nMO, gamma)
      CHARACTER(LEN=*), INTENT(IN)                       :: outdir
      INTEGER, INTENT(IN)                                :: nRI, nMO
      REAL(KIND=dp), INTENT(IN)                          :: gamma(nRI, nMO, nMO)

      INTEGER                                            :: iunit, iP, ip_mo, iq_mo
      CHARACTER(LEN=512)                                 :: fname

      fname = TRIM(outdir)//'/CoulombVertex.yaml'
      iunit = get_unit_number(fname)
      OPEN (UNIT=iunit, FILE=TRIM(fname), STATUS='REPLACE', ACTION='WRITE')
      WRITE (iunit, '(A)')      'version: 100'
      WRITE (iunit, '(A)')      'type: Tensor'
      WRITE (iunit, '(A)')      'scalarType: Real64'
      WRITE (iunit, '(A)')      'dimensions:'
      WRITE (iunit, '(A,I0)')   '- length: ', nRI
      WRITE (iunit, '(A)')      '  type: AuxiliaryField'
      WRITE (iunit, '(A,I0)')   '- length: ', nMO
      WRITE (iunit, '(A)')      '  type: State'
      WRITE (iunit, '(A,I0)')   '- length: ', nMO
      WRITE (iunit, '(A)')      '  type: State'
      WRITE (iunit, '(A)')      'elements:'
      WRITE (iunit, '(A)')      '  type: TextFile'
      WRITE (iunit, '(A)')      'unit: 1.0'
      CLOSE (iunit)

      ! Fortran column-major (P fastest, then p, then q) — matches CC4S CTF convention
      fname = TRIM(outdir)//'/CoulombVertex.elements'
      iunit = get_unit_number(fname)
      OPEN (UNIT=iunit, FILE=TRIM(fname), STATUS='REPLACE', ACTION='WRITE')
      DO iq_mo = 1, nMO
         DO ip_mo = 1, nMO
            DO iP = 1, nRI
               WRITE (iunit, '(ES25.16)') gamma(iP, ip_mo, iq_mo)
            END DO
         END DO
      END DO
      CLOSE (iunit)

   END SUBROUTINE write_coulomb_vertex_files

! **************************************************************************************************
!> \brief  Write cc4s.in driver file.
!> \param  outdir   output directory
!> \param  method   one of: MP2, CCSD, CCSD_T, DRCCD, APXRPA
! **************************************************************************************************
   SUBROUTINE write_cc4s_in(outdir, method)
      CHARACTER(LEN=*), INTENT(IN)                       :: outdir, method

      INTEGER                                            :: iunit
      CHARACTER(LEN=512)                                 :: fname

      fname = TRIM(outdir)//'/cc4s.in'
      iunit = get_unit_number(fname)
      OPEN (UNIT=iunit, FILE=TRIM(fname), STATUS='REPLACE', ACTION='WRITE')

      ! --- common header: read tensors, slice, build integrals ---
      WRITE (iunit, '(A)') '- name: Read'
      WRITE (iunit, '(A)') '  in:'
      WRITE (iunit, '(A)') '    fileName: CoulombVertex.yaml'
      WRITE (iunit, '(A)') '  out:'
      WRITE (iunit, '(A)') '    destination: CoulombVertex'
      WRITE (iunit, '(A)') '- name: Read'
      WRITE (iunit, '(A)') '  in:'
      WRITE (iunit, '(A)') '    fileName: EigenEnergies.yaml'
      WRITE (iunit, '(A)') '  out:'
      WRITE (iunit, '(A)') '    destination: EigenEnergies'
      WRITE (iunit, '(A)') '- name: DefineHolesAndParticles'
      WRITE (iunit, '(A)') '  in:'
      WRITE (iunit, '(A)') '    eigenEnergies: EigenEnergies'
      WRITE (iunit, '(A)') '  out:'
      WRITE (iunit, '(A)') '    slicedEigenEnergies: EigenEnergies'
      WRITE (iunit, '(A)') '- name: SliceOperator'
      WRITE (iunit, '(A)') '  in:'
      WRITE (iunit, '(A)') '    slicedEigenEnergies: EigenEnergies'
      WRITE (iunit, '(A)') '    operator: CoulombVertex'
      WRITE (iunit, '(A)') '  out:'
      WRITE (iunit, '(A)') '    slicedOperator: CoulombVertex'
      WRITE (iunit, '(A)') '- name: VertexCoulombIntegrals'
      WRITE (iunit, '(A)') '  in:'
      WRITE (iunit, '(A)') '    slicedCoulombVertex: CoulombVertex'
      WRITE (iunit, '(A)') '  out:'
      WRITE (iunit, '(A)') '    coulombIntegrals: CoulombIntegrals'

      ! --- MP2 (always included) ---
      WRITE (iunit, '(A)') '- name: SecondOrderPerturbationTheory'
      WRITE (iunit, '(A)') '  in:'
      WRITE (iunit, '(A)') '    coulombIntegrals: CoulombIntegrals'
      WRITE (iunit, '(A)') '    slicedEigenEnergies: EigenEnergies'
      WRITE (iunit, '(A)') '  out:'
      WRITE (iunit, '(A)') '    energy: Mp2Energy'

      ! --- CoupledCluster (CCSD, CCSD_T, DRCCD, APXRPA) ---
      SELECT CASE (TRIM(method))
      CASE ('CCSD', 'CCSD_T')
         CALL write_coupled_cluster_step(iunit, 'Ccsd', 'Amplitudes')
      CASE ('DRCCD', 'APXRPA')
         CALL write_coupled_cluster_step(iunit, 'Drccd', 'Amplitudes')
      END SELECT

      ! --- PerturbativeTriples (CCSD_T only) ---
      IF (TRIM(method) == 'CCSD_T') THEN
         WRITE (iunit, '(A)') '- name: PerturbativeTriples'
         WRITE (iunit, '(A)') '  in:'
         WRITE (iunit, '(A)') '    slicedEigenEnergies: EigenEnergies'
         WRITE (iunit, '(A)') '    amplitudes: Amplitudes'
         WRITE (iunit, '(A)') '    coulombIntegrals: CoulombIntegrals'
         WRITE (iunit, '(A)') '  out: {}'
      END IF

      ! --- ApxRpa second pass (APXRPA only) ---
      IF (TRIM(method) == 'APXRPA') THEN
         CALL write_coupled_cluster_step(iunit, 'ApxRpa', 'ApxRpaAmplitudes')
      END IF

      CLOSE (iunit)

   END SUBROUTINE write_cc4s_in

! **************************************************************************************************
!> \brief  Write a single CoupledCluster step to the cc4s.in file.
! **************************************************************************************************
   SUBROUTINE write_coupled_cluster_step(iunit, cc_method, amp_key)
      INTEGER, INTENT(IN)                                :: iunit
      CHARACTER(LEN=*), INTENT(IN)                      :: cc_method, amp_key

      WRITE (iunit, '(A)') '- name: CoupledCluster'
      WRITE (iunit, '(A)') '  in:'
      WRITE (iunit, '(A)') '    method: '//TRIM(cc_method)
      WRITE (iunit, '(A)') '    integralsSliceSize: 100'
      WRITE (iunit, '(A)') '    slicedEigenEnergies: EigenEnergies'
      WRITE (iunit, '(A)') '    coulombIntegrals: CoulombIntegrals'
      WRITE (iunit, '(A)') '    slicedCoulombVertex: CoulombVertex'
      WRITE (iunit, '(A)') '    maxIterations: 30'
      WRITE (iunit, '(A)') '    energyConvergence: 1.0e-05'
      WRITE (iunit, '(A)') '    amplitudesConvergence: 1.0e-05'
      WRITE (iunit, '(A)') '    mixer:'
      WRITE (iunit, '(A)') '      type: DiisMixer'
      WRITE (iunit, '(A)') '      maxResidua: 5'
      WRITE (iunit, '(A)') '  out:'
      WRITE (iunit, '(A)') '    amplitudes: '//TRIM(amp_key)

   END SUBROUTINE write_coupled_cluster_step

END MODULE cc4s_export
