!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Restricted molecular EXX-OEP with Coulomb auxiliary charge functions.
!>        Implements the preprocessing of Trushin and Goerling, JCP 155, 054109 (2021),
!>        as in PyOEP: normalized Coulomb orthogonalization, charge/HOMO constraints,
!>        and filtering the UNWEIGHTED occupied-virtual coupling Gram matrix.
!>        Dense arrays are replicated; this initial implementation targets small molecules.
! **************************************************************************************************
MODULE qs_oep
   USE ai_contraction_sphi,             ONLY: ab_contract,&
                                              abc_contract
   USE atomic_kind_types,               ONLY: get_atomic_kind
   USE basis_set_types,                 ONLY: gto_basis_set_p_type,&
                                              gto_basis_set_type
   USE cell_types,                      ONLY: cell_type
   USE coulomb_integral_interface,      ONLY: compute_coulomb_2c,&
                                              compute_coulomb_3c,&
                                              coulomb_integral_cleanup,&
                                              coulomb_integral_context_type,&
                                              coulomb_integral_init,&
                                              coulomb_operator_type
   USE cp_blacs_env,                    ONLY: cp_blacs_env_type
   USE cp_control_types,                ONLY: dft_control_type
   USE cp_dbcsr_api,                    ONLY: dbcsr_p_type
   USE cp_dbcsr_operations,             ONLY: copy_dbcsr_to_fm,&
                                              copy_fm_to_dbcsr
   USE cp_files,                        ONLY: close_file,&
                                              open_file
   USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                              cp_fm_struct_release,&
                                              cp_fm_struct_type
   USE cp_fm_types,                     ONLY: cp_fm_create,&
                                              cp_fm_get_submatrix,&
                                              cp_fm_release,&
                                              cp_fm_set_submatrix,&
                                              cp_fm_type
   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_get_default_unit_nr
   USE input_constants,                 ONLY: do_potential_coulomb
   USE input_section_types,             ONLY: section_vals_get,&
                                              section_vals_get_subs_vals,&
                                              section_vals_type,&
                                              section_vals_val_get
   USE kinds,                           ONLY: default_path_length,&
                                              dp
   USE mathlib,                         ONLY: diamat_all
   USE message_passing,                 ONLY: mp_para_env_type
   USE orbital_pointers,                ONLY: indco,&
                                              init_orbital_pointers,&
                                              ncoset
   USE orbital_transformation_matrices, ONLY: init_spherical_harmonics
   USE particle_methods,                ONLY: get_particle_set
   USE particle_types,                  ONLY: particle_type
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_interactions,                 ONLY: init_interaction_radii_orb_basis
   USE qs_kind_types,                   ONLY: get_qs_kind,&
                                              qs_kind_type
   USE qs_ks_types,                     ONLY: qs_ks_env_type,&
                                              qs_oep_type
   USE qs_mo_types,                     ONLY: get_mo_set,&
                                              mo_set_type
   USE scf_control_types,               ONLY: scf_control_type
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE
   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_oep'
   PUBLIC :: oep_validate, oep_exchange_matrix

CONTAINS

! **************************************************************************************************
!> \brief Check the explicitly supported molecular OEP domain before dispatching HFX.
!> \param qs_env ...
!> \param hf_section ...
!> \param calculate_forces ...
! **************************************************************************************************
   SUBROUTINE oep_validate(qs_env, hf_section, calculate_forces)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(section_vals_type), POINTER                   :: hf_section
      LOGICAL, INTENT(IN)                                :: calculate_forces

      INTEGER                                            :: nrep, potential
      LOGICAL                                            :: ace, enabled, kpoints
      REAL(KIND=dp)                                      :: fraction, threshold
      TYPE(cell_type), POINTER                           :: cell
      TYPE(dft_control_type), POINTER                    :: dft
      TYPE(scf_control_type), POINTER                    :: scf
      TYPE(section_vals_type), POINTER                   :: section

      section => section_vals_get_subs_vals(hf_section, "OEP")
      CALL section_vals_get(section, explicit=enabled)
      IF (.NOT. enabled) RETURN
      CALL get_qs_env(qs_env, dft_control=dft, scf_control=scf, cell=cell, do_kpoints=kpoints)
      IF (calculate_forces) CPABORT("EXX-OEP currently supports ENERGY calculations only")
      IF (dft%nspins /= 1 .OR. dft%roks) CPABORT("EXX-OEP requires a closed-shell restricted calculation")
      IF (ANY(cell%perd /= 0) .OR. kpoints) CPABORT("EXX-OEP requires CELL PERIODIC NONE and no KPOINTS")
      IF (dft%do_admm .OR. qs_env%run_rtp) CPABORT("EXX-OEP does not support ADMM or real-time propagation")
      IF (scf%use_ot .OR. .NOT. scf%use_diag) CPABORT("EXX-OEP requires SCF DIAGONALIZATION and ADDED_MOS -1")
      IF (scf%level_shift /= 0.0_dp) CPABORT("EXX-OEP requires unshifted orbital energies: LEVEL_SHIFT 0")
      CALL section_vals_val_get(section, "THR_FAI_OEP", r_val=threshold)
      IF (threshold <= 0.0_dp) CPABORT("EXX-OEP THR_FAI_OEP must be positive")
      CALL section_vals_get(hf_section, n_repetition=nrep)
      IF (nrep /= 1) CPABORT("EXX-OEP requires exactly one HF section")
      CALL section_vals_val_get(hf_section, "FRACTION", r_val=fraction)
      IF (ABS(fraction - 1.0_dp) > 1.0E-14_dp) CPABORT("EXX-OEP requires HF FRACTION 1")
      CALL section_vals_val_get(hf_section, "INTERACTION_POTENTIAL%POTENTIAL_TYPE", i_val=potential)
      IF (potential /= do_potential_coulomb) CPABORT("Molecular EXX-OEP requires the COULOMB interaction")
      CALL section_vals_val_get(hf_section, "ACE%ACTIVE", l_val=ace)
      IF (ace) CPABORT("EXX-OEP requires the full exchange matrix; disable ACE")
   END SUBROUTINE oep_validate

! **************************************************************************************************
!> \brief Replace only the nonlocal exchange matrix, preserving the HFX energy.
!> \param qs_env ...
!> \param matrix_x pure nonlocal exchange on entry, local OEP matrix on exit
!> \param hf_section ...
! **************************************************************************************************
   SUBROUTINE oep_exchange_matrix(qs_env, matrix_x, hf_section)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: matrix_x
      TYPE(section_vals_type), POINTER                   :: hf_section

      CHARACTER(LEN=*), PARAMETER :: routineN = 'oep_exchange_matrix'

      CHARACTER(LEN=default_path_length)                 :: debug_file, potential_file
      INTEGER                                            :: handle, i, nao, nmo, nocc, unit_dump, &
                                                            unit_nr
      LOGICAL                                            :: enabled, homo_condition
      REAL(KIND=dp)                                      :: charge_error, residual, threshold
      REAL(KIND=dp), ALLOCATABLE                         :: c(:, :), coefficients(:), k(:, :), &
                                                            positions(:, :), reference(:), &
                                                            v(:, :), w(:, :)
      REAL(KIND=dp), DIMENSION(:), POINTER               :: eps, occ
      TYPE(cp_blacs_env_type), POINTER                   :: blacs
      TYPE(cp_fm_struct_type), POINTER                   :: layout
      TYPE(cp_fm_type)                                   :: work
      TYPE(cp_fm_type), POINTER                          :: mo_coeff
      TYPE(mo_set_type), DIMENSION(:), POINTER           :: mos
      TYPE(mp_para_env_type), POINTER                    :: para
      TYPE(particle_type), DIMENSION(:), POINTER         :: particles
      TYPE(qs_ks_env_type), POINTER                      :: ks_env
      TYPE(section_vals_type), POINTER                   :: section

      section => section_vals_get_subs_vals(hf_section, "OEP")
      CALL section_vals_get(section, explicit=enabled)
      IF (.NOT. enabled) RETURN
      CALL get_qs_env(qs_env, ks_env=ks_env, mos=mos, particle_set=particles, blacs_env=blacs, para_env=para)
      CALL get_mo_set(mos(1), nao=nao, nmo=nmo, homo=nocc, mo_coeff=mo_coeff, &
                      eigenvalues=eps, occupation_numbers=occ)
      unit_nr = -1
      IF (para%is_source()) unit_nr = cp_logger_get_default_unit_nr(cp_get_default_logger())
      IF (nmo /= nao) CPABORT("EXX-OEP requires all orbitals: set ADDED_MOS -1 and avoid overlap truncation")
      IF (nocc < 1 .OR. nocc >= nmo) CPABORT("EXX-OEP requires occupied and virtual orbitals")
      IF (ANY(ABS(occ(:nocc) - 2.0_dp) > 1.0E-10_dp) .OR. &
          ANY(ABS(occ(nocc + 1:)) > 1.0E-10_dp)) CPABORT("EXX-OEP requires occupations 2 and 0, without smearing")

      ! Atomic/random initial guesses have no orbital energies. One conventional HFX
      ! diagonalization supplies the initial complete eigensystem, as the HF start in PyOEP.
      IF (MAXVAL(ABS(eps)) < 1.0E-12_dp) THEN
         IF (unit_nr > 0) WRITE (unit_nr, '(T2,A)') "OEP| Initial HFX diagonalization to obtain orbital energies"
         RETURN
      END IF
      CALL timeset(routineN, handle)
      ALLOCATE (positions(3, SIZE(particles)))
      DO i = 1, SIZE(particles)
         positions(:, i) = particles(i)%r
      END DO
      IF (ALLOCATED(ks_env%oep)) THEN
         IF (ANY(ABS(positions - ks_env%oep%positions) > 1.0E-14_dp)) DEALLOCATE (ks_env%oep)
      END IF
      IF (.NOT. ALLOCATED(ks_env%oep)) THEN
         ALLOCATE (ks_env%oep)
         ks_env%oep%positions = positions
         CALL oep_integrals(qs_env, ks_env%oep, nao)
         CALL oep_metric_preprocess(ks_env%oep)
         IF (unit_nr > 0) WRITE (unit_nr, '(T2,A,I0)') &
            "OEP| Auxiliary charge functions: ", SIZE(ks_env%oep%charge)
      END IF
      ALLOCATE (c(nao, nmo), k(nao, nao), v(nao, nao))
      CALL cp_fm_get_submatrix(mo_coeff, c)
      NULLIFY (layout)
      CALL cp_fm_struct_create(layout, context=blacs, nrow_global=nao, ncol_global=nao, para_env=para)
      CALL cp_fm_create(work, layout, set_zero=.TRUE.)
      CALL cp_fm_struct_release(layout)
      CALL copy_dbcsr_to_fm(matrix_x(1, 1)%matrix, work)
      CALL cp_fm_get_submatrix(work, k)
      CALL section_vals_val_get(section, "THR_FAI_OEP", r_val=threshold)
      CALL section_vals_val_get(section, "HOMO_CONDITION", l_val=homo_condition)
      CALL oep_solve(ks_env%oep, c, eps, nocc, k, threshold, homo_condition, &
                     v, coefficients, reference, w, residual, charge_error)
      CALL cp_fm_set_submatrix(work, v)
      CALL copy_fm_to_dbcsr(work, matrix_x(1, 1)%matrix)
      CALL cp_fm_release(work)
      IF (unit_nr > 0) THEN
         WRITE (unit_nr, '(T2,A,I0)') "OEP| Retained auxiliary functions: ", SIZE(w, 2)
         WRITE (unit_nr, '(T2,A,ES12.4)') "OEP| Linear solve residual: ", residual
         WRITE (unit_nr, '(T2,A,ES12.4)') "OEP| Charge constraint error: ", charge_error
      END IF
      IF (para%is_source()) THEN
         CALL section_vals_val_get(section, "POTENTIAL_FILE_NAME", c_val=potential_file)
         CALL section_vals_val_get(section, "DEBUG_FILE_NAME", c_val=debug_file)
         IF (LEN_TRIM(potential_file) > 0) THEN
            IF (potential_file == debug_file) CPABORT("OEP potential and debug filenames must be different")
            CALL oep_write_potential(qs_env, TRIM(potential_file), coefficients, reference, &
                                     threshold, homo_condition, SIZE(w, 2))
         END IF
         IF (LEN_TRIM(debug_file) > 0) THEN
            OPEN (NEWUNIT=unit_dump, FILE=TRIM(debug_file), ACCESS="STREAM", FORM="UNFORMATTED", STATUS="REPLACE")
            WRITE (unit_dump) nao, SIZE(coefficients), nmo, nocc, SIZE(w, 2)
            WRITE (unit_dump) ks_env%oep%three_center, ks_env%oep%metric, ks_env%oep%charge, &
               c, eps, k, coefficients, v, w
            CLOSE (unit_dump)
         END IF
      END IF
      CALL timestop(handle)
   END SUBROUTINE oep_exchange_matrix

! **************************************************************************************************
!> \brief Export the Coulomb source as unnormalized Cartesian primitive Gaussians.
!>        Folding the spherical contractions into the coefficients makes the file independent
!>        of basis libraries, normalization conventions, and angular-function ordering.
!> \param qs_env ...
!> \param filename ...
!> \param coefficients total exchange coefficients
!> \param reference reference exchange coefficients
!> \param threshold ...
!> \param homo_condition ...
!> \param nret ...
! **************************************************************************************************
   SUBROUTINE oep_write_potential(qs_env, filename, coefficients, reference, threshold, homo_condition, nret)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      CHARACTER(LEN=*), INTENT(IN)                       :: filename
      REAL(KIND=dp), INTENT(IN)                          :: coefficients(:), reference(:), threshold
      LOGICAL, INTENT(IN)                                :: homo_condition
      INTEGER, INTENT(IN)                                :: nret

      CHARACTER(LEN=2)                                   :: element
      INTEGER                                            :: first, ia, ico, ip, iset, j, ka, nc, &
                                                            ncart, ns, offset, unit_nr, z
      INTEGER, ALLOCATABLE                               :: off_aux(:)
      LOGICAL                                            :: first_term
      REAL(KIND=dp), ALLOCATABLE                         :: cart_reference(:), cart_total(:)
      TYPE(gto_basis_set_p_type), ALLOCATABLE            :: aux(:)
      TYPE(gto_basis_set_type), POINTER                  :: basis
      TYPE(particle_type), DIMENSION(:), POINTER         :: particles
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: kinds_set
      TYPE(qs_ks_env_type), POINTER                      :: ks_env

      CALL get_qs_env(qs_env, qs_kind_set=kinds_set, particle_set=particles, ks_env=ks_env)
      ALLOCATE (aux(SIZE(kinds_set)), off_aux(SIZE(particles)))
      DO ka = 1, SIZE(kinds_set)
         CALL get_qs_kind(kinds_set(ka), basis_set=aux(ka)%gto_basis_set, basis_type="OEP")
      END DO
      CALL get_particle_set(particles, kinds_set, first_sgf=off_aux, basis=aux)
      CALL open_file(filename, unit_nr, file_status="REPLACE", file_action="WRITE")
      WRITE (unit_nr, '(A)') '{', &
         '  "format": "CP2K_OEP_POTENTIAL",', &
         '  "version": 1,', &
         '  "length_unit": "bohr",', &
         '  "potential_unit": "hartree",', &
         '  "representation": "unnormalized_cartesian_gaussians",', &
         '  "term_columns": ["lx", "ly", "lz", "exponent", "total", "reference"],'
      IF (homo_condition) THEN
         WRITE (unit_nr, '(A)') '  "homo_condition": true,'
      ELSE
         WRITE (unit_nr, '(A)') '  "homo_condition": false,'
      END IF
      WRITE (unit_nr, '(A,ES26.17E3,A)') '  "threshold": ', threshold, ','
      WRITE (unit_nr, '(A,I0,A)') '  "naux": ', SIZE(coefficients), ','
      WRITE (unit_nr, '(A,I0,A)') '  "nretained": ', nret, ','
      WRITE (unit_nr, '(A,ES26.17E3,A)') '  "total_charge": ', DOT_PRODUCT(ks_env%oep%charge, coefficients), ','
      WRITE (unit_nr, '(A,ES26.17E3,A)') '  "reference_charge": ', DOT_PRODUCT(ks_env%oep%charge, reference), ','
      WRITE (unit_nr, '(A)') '  "atoms": ['
      DO ia = 1, SIZE(particles)
         CALL get_atomic_kind(particles(ia)%atomic_kind, kind_number=ka, element_symbol=element, z=z)
         basis => aux(ka)%gto_basis_set
         WRITE (unit_nr, '(A)') '    {'
         WRITE (unit_nr, '(A)') '      "element": "'//TRIM(element)//'",'
         WRITE (unit_nr, '(A,I0,A)') '      "atomic_number": ', z, ','
         WRITE (unit_nr, '(A,2(ES26.17E3,A),ES26.17E3,A)') &
            '      "position": [', particles(ia)%r(1), ',', particles(ia)%r(2), ',', particles(ia)%r(3), '],'
         WRITE (unit_nr, '(A)') '      "terms": ['
         first_term = .TRUE.
         DO iset = 1, basis%nset
            nc = ncoset(basis%lmax(iset))
            ncart = nc*basis%npgf(iset)
            first = basis%first_sgf(1, iset)
            ns = basis%nsgf_set(iset)
            offset = off_aux(ia) + first - 1
            cart_total = MATMUL(basis%sphi(:ncart, first:first + ns - 1), coefficients(offset:offset + ns - 1))
            cart_reference = MATMUL(basis%sphi(:ncart, first:first + ns - 1), reference(offset:offset + ns - 1))
            DO ip = 1, basis%npgf(iset)
               DO ico = 1, nc
                  j = (ip - 1)*nc + ico
                  IF (cart_total(j) == 0.0_dp .AND. cart_reference(j) == 0.0_dp) CYCLE
                  IF (.NOT. first_term) WRITE (unit_nr, '(A)') ','
                  WRITE (unit_nr, '(A,3(I0,A),2(ES26.17E3,A),ES26.17E3,A)', ADVANCE="NO") &
                     '        [', indco(1, ico), ',', indco(2, ico), ',', indco(3, ico), ',', &
                     basis%zet(ip, iset), ',', cart_total(j), ',', cart_reference(j), ']'
                  first_term = .FALSE.
               END DO
            END DO
         END DO
         WRITE (unit_nr, '(A)') ''
         WRITE (unit_nr, '(A)') '      ]'
         IF (ia < SIZE(particles)) THEN
            WRITE (unit_nr, '(A)') '    },'
         ELSE
            WRITE (unit_nr, '(A)') '    }'
         END IF
      END DO
      WRITE (unit_nr, '(A)') '  ]', '}'
      CALL close_file(unit_nr)
   END SUBROUTINE oep_write_potential

! **************************************************************************************************
!> \brief Build analytic, isolated two/three-center integrals and auxiliary charges.
!> \param qs_env ...
!> \param oep ...
!> \param nao ...
! **************************************************************************************************
   SUBROUTINE oep_integrals(qs_env, oep, nao)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(qs_oep_type), INTENT(INOUT)                   :: oep
      INTEGER, INTENT(IN)                                :: nao

      INTEGER                                            :: fa, fb, fc, ga, gb, gc, ia, ib, ic, ka, &
                                                            kb, ma, max_l, mb, mc, na, natom, &
                                                            naux, nb, nc, nkind, sa, sb, sc
      INTEGER, ALLOCATABLE                               :: kind_of(:), off_ao(:), off_aux(:), &
                                                            sizes_aux(:)
      REAL(KIND=dp), ALLOCATABLE                         :: cart2(:, :), cart3(:, :, :), sph2(:, :), &
                                                            sph3(:, :, :)
      TYPE(coulomb_integral_context_type)                :: context
      TYPE(coulomb_operator_type)                        :: potential
      TYPE(gto_basis_set_p_type), ALLOCATABLE            :: ao(:), aux(:)
      TYPE(gto_basis_set_type), POINTER                  :: ba, bb, bc
      TYPE(particle_type), DIMENSION(:), POINTER         :: particles
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: kinds_set

      CALL get_qs_env(qs_env, qs_kind_set=kinds_set, particle_set=particles)
      nkind = SIZE(kinds_set)
      natom = SIZE(particles)
      ALLOCATE (ao(nkind), aux(nkind), kind_of(natom), off_ao(natom), off_aux(natom), sizes_aux(natom))
      max_l = 0
      DO ka = 1, nkind
         CALL get_qs_kind(kinds_set(ka), basis_set=ao(ka)%gto_basis_set, basis_type="ORB")
         CALL get_qs_kind(kinds_set(ka), basis_set=aux(ka)%gto_basis_set, basis_type="OEP")
         IF (.NOT. ASSOCIATED(aux(ka)%gto_basis_set)) CPABORT("EXX-OEP requires BASIS_SET OEP on every KIND")
         CALL init_interaction_radii_orb_basis(aux(ka)%gto_basis_set, 1.0E-16_dp)
         max_l = MAX(max_l, MAXVAL(ao(ka)%gto_basis_set%lmax), MAXVAL(aux(ka)%gto_basis_set%lmax))
      END DO
      CALL init_orbital_pointers(2*max_l)
      CALL init_spherical_harmonics(max_l, -1)
      CALL coulomb_integral_init(context, max_am_2c=max_l, max_am_3c=max_l)
      CALL get_particle_set(particles, kinds_set, first_sgf=off_ao, basis=ao)
      CALL get_particle_set(particles, kinds_set, first_sgf=off_aux, nsgf=sizes_aux, basis=aux)
      DO ia = 1, natom
         CALL get_atomic_kind(particles(ia)%atomic_kind, kind_number=kind_of(ia))
      END DO
      naux = SUM(sizes_aux)
      ALLOCATE (oep%three_center(nao, nao, naux), oep%metric(naux, naux), oep%charge(naux))
      oep%three_center = 0.0_dp
      oep%metric = 0.0_dp
      oep%charge = 0.0_dp
      DO ia = 1, natom
         ka = kind_of(ia)
         ba => aux(ka)%gto_basis_set
         CALL oep_auxiliary_charge(ba, oep%charge(off_aux(ia):off_aux(ia) + sizes_aux(ia) - 1))
         DO ib = 1, natom
            kb = kind_of(ib)
            bb => aux(kb)%gto_basis_set
            DO sa = 1, ba%nset
               na = ba%npgf(sa)*ncoset(ba%lmax(sa))
               fa = ba%first_sgf(1, sa)
               ma = ba%nsgf_set(sa)
               ga = off_aux(ia) + fa - 1
               DO sb = 1, bb%nset
                  nb = bb%npgf(sb)*ncoset(bb%lmax(sb))
                  fb = bb%first_sgf(1, sb)
                  mb = bb%nsgf_set(sb)
                  gb = off_aux(ib) + fb - 1
                  ALLOCATE (cart2(na, nb), sph2(ma, mb))
                  cart2 = 0.0_dp
                  CALL compute_coulomb_2c(context, ba%lmin(sa), ba%lmax(sa), bb%lmin(sb), bb%lmax(sb), &
                                          ba%npgf(sa), bb%npgf(sb), ba%zet(:, sa), bb%zet(:, sb), &
                                          ba%pgf_radius(:, sa), bb%pgf_radius(:, sb), &
                                          particles(ia)%r, particles(ib)%r, cart2, potential)
                  CALL ab_contract(sph2, cart2, ba%sphi(:, fa:), bb%sphi(:, fb:), na, nb, ma, mb)
                  oep%metric(ga:ga + ma - 1, gb:gb + mb - 1) = sph2
                  DEALLOCATE (cart2, sph2)
               END DO
            END DO
         END DO
      END DO
      DO ia = 1, natom
         ba => ao(kind_of(ia))%gto_basis_set
         DO ib = 1, natom
            bb => ao(kind_of(ib))%gto_basis_set
            DO ic = 1, natom
               bc => aux(kind_of(ic))%gto_basis_set
               DO sa = 1, ba%nset
                  na = ba%npgf(sa)*ncoset(ba%lmax(sa))
                  fa = ba%first_sgf(1, sa)
                  ma = ba%nsgf_set(sa)
                  ga = off_ao(ia) + fa - 1
                  DO sb = 1, bb%nset
                     nb = bb%npgf(sb)*ncoset(bb%lmax(sb))
                     fb = bb%first_sgf(1, sb)
                     mb = bb%nsgf_set(sb)
                     gb = off_ao(ib) + fb - 1
                     DO sc = 1, bc%nset
                        nc = bc%npgf(sc)*ncoset(bc%lmax(sc))
                        fc = bc%first_sgf(1, sc)
                        mc = bc%nsgf_set(sc)
                        gc = off_aux(ic) + fc - 1
                        ALLOCATE (cart3(na, nb, nc), sph3(ma, mb, mc))
                        cart3 = 0.0_dp
                        CALL compute_coulomb_3c(context, ba%lmin(sa), ba%lmax(sa), bb%lmin(sb), bb%lmax(sb), &
                                                bc%lmin(sc), bc%lmax(sc), ba%npgf(sa), bb%npgf(sb), bc%npgf(sc), &
                                                ba%zet(:, sa), bb%zet(:, sb), bc%zet(:, sc), &
                                                ba%pgf_radius(:, sa), bb%pgf_radius(:, sb), bc%pgf_radius(:, sc), &
                                                particles(ia)%r, particles(ib)%r, particles(ic)%r, cart3, potential)
                        CALL abc_contract(sph3, cart3, ba%sphi(:, fa:), bb%sphi(:, fb:), bc%sphi(:, fc:), &
                                          na, nb, nc, ma, mb, mc)
                        oep%three_center(ga:ga + ma - 1, gb:gb + mb - 1, gc:gc + mc - 1) = sph3
                        DEALLOCATE (cart3, sph3)
                     END DO
                  END DO
               END DO
            END DO
         END DO
      END DO
      CALL coulomb_integral_cleanup(context)
   END SUBROUTINE oep_integrals

! **************************************************************************************************
!> \brief Analytic integral of each auxiliary Gaussian (PyOEP evaluates this on a grid).
!> \param basis ...
!> \param charge ...
! **************************************************************************************************
   SUBROUTINE oep_auxiliary_charge(basis, charge)
      TYPE(gto_basis_set_type), INTENT(IN)               :: basis
      REAL(KIND=dp), INTENT(OUT)                         :: charge(:)

      INTEGER                                            :: axis, first, ico, ip, iset, j, m, nc, &
                                                            ncart, ns, power
      REAL(KIND=dp)                                      :: moment, zeta
      REAL(KIND=dp), ALLOCATABLE                         :: moments(:)

      charge = 0.0_dp
      DO iset = 1, basis%nset
         nc = ncoset(basis%lmax(iset))
         ncart = nc*basis%npgf(iset)
         first = basis%first_sgf(1, iset)
         ns = basis%nsgf_set(iset)
         ALLOCATE (moments(ncart))
         DO ip = 1, basis%npgf(iset)
            zeta = basis%zet(ip, iset)
            DO ico = 1, nc
               moment = 1.0_dp
               DO axis = 1, 3
                  power = indco(axis, ico)
                  IF (MOD(power, 2) /= 0) THEN
                     moment = 0.0_dp
                  ELSE
                     moment = moment*SQRT(ACOS(-1.0_dp)/zeta)
                     DO m = 1, power/2
                        moment = moment*REAL(2*m - 1, dp)/(2.0_dp*zeta)
                     END DO
                  END IF
               END DO
               moments((ip - 1)*nc + ico) = moment
            END DO
         END DO
         DO j = first, first + ns - 1
            charge(j) = DOT_PRODUCT(basis%sphi(:ncart, j), moments)
         END DO
         DEALLOCATE (moments)
      END DO
   END SUBROUTINE oep_auxiliary_charge

! **************************************************************************************************
!> \brief PyOEP get_WII and get_y_and_yII, with analytic auxiliary charges.
!> \param oep ...
! **************************************************************************************************
   SUBROUTINE oep_metric_preprocess(oep)
      TYPE(qs_oep_type), INTENT(INOUT)                   :: oep

      INTEGER                                            :: i, n
      REAL(KIND=dp), ALLOCATABLE                         :: eig(:), norm(:)

      n = SIZE(oep%charge)
      ALLOCATE (norm(n), eig(n), oep%w2(n, n), oep%y2(n))
      DO i = 1, n
         IF (oep%metric(i, i) <= 0.0_dp) CPABORT("EXX-OEP auxiliary metric has a nonpositive diagonal")
         norm(i) = 1.0_dp/SQRT(oep%metric(i, i))
      END DO
      oep%w2 = oep%metric*SPREAD(norm, 1, n)*SPREAD(norm, 2, n)
      CALL diamat_all(oep%w2, eig)
      IF (MINVAL(eig) <= 0.0_dp) CPABORT("EXX-OEP auxiliary Coulomb metric is not positive definite")
      DO i = 1, n
         oep%w2(:, i) = norm*oep%w2(:, i)/SQRT(eig(i))
      END DO
      oep%y2 = MATMUL(TRANSPOSE(oep%w2), oep%charge)
   END SUBROUTINE oep_metric_preprocess

! **************************************************************************************************
!> \brief Solve the constrained, filtered EXX-OEP equation in the current orbital space.
!> \param oep ...
!> \param c ...
!> \param eps ...
!> \param nocc ...
!> \param k ...
!> \param threshold ...
!> \param homo_condition ...
!> \param v ...
!> \param coefficients ...
!> \param vref ...
!> \param w ...
!> \param residual ...
!> \param charge_error ...
! **************************************************************************************************
   SUBROUTINE oep_solve(oep, c, eps, nocc, k, threshold, homo_condition, v, coefficients, vref, w, residual, charge_error)
      TYPE(qs_oep_type), INTENT(IN)                      :: oep
      REAL(KIND=dp), INTENT(IN)                          :: c(:, :), eps(:)
      INTEGER, INTENT(IN)                                :: nocc
      REAL(KIND=dp), INTENT(IN)                          :: k(:, :), threshold
      LOGICAL, INTENT(IN)                                :: homo_condition
      REAL(KIND=dp), INTENT(OUT)                         :: v(:, :)
      REAL(KIND=dp), ALLOCATABLE, INTENT(OUT)            :: coefficients(:), vref(:), w(:, :)
      REAL(KIND=dp), INTENT(OUT)                         :: residual, charge_error

      INTEGER                                            :: a, i, j, naux, nconstraint, ndeg, nfree, &
                                                            nmo, nov, nret, nvir, p
      REAL(KIND=dp)                                      :: constraints(2, 2), denom, normalization, &
                                                            TARGET(2)
      REAL(KIND=dp), ALLOCATABLE :: bov(:, :), d(:, :), delta(:), eval(:), gram(:, :), kmo(:, :), &
         mo(:, :), projector(:, :), q(:), rhs(:), solution(:), trans(:, :), u(:), vc2(:), &
         weights(:), x(:, :), z(:), z2(:)

      naux = SIZE(oep%charge)
      nmo = SIZE(c, 2)
      nvir = nmo - nocc
      nov = nocc*nvir
      ALLOCATE (bov(naux, nov), u(naux), z(naux), weights(nov), delta(nov))
      ndeg = 1
      DO i = nocc - 1, 1, -1
         IF (ABS(eps(i) - eps(nocc)) >= 1.0E-12_dp) EXIT
         ndeg = ndeg + 1
      END DO
      DO p = 1, naux
         mo = MATMUL(TRANSPOSE(c), MATMUL(oep%three_center(:, :, p), c))
         u(p) = 0.0_dp
         DO i = 1, nocc
            u(p) = u(p) + 2.0_dp*mo(i, i)
            DO a = 1, nvir
               bov(p, (i - 1)*nvir + a) = mo(nocc + a, i)
            END DO
         END DO
         z(p) = 0.0_dp
         DO i = nocc - ndeg + 1, nocc
            z(p) = z(p) + mo(i, i)/REAL(ndeg, dp)
         END DO
      END DO
      kmo = MATMUL(TRANSPOSE(c), MATMUL(k, c))
      vc2 = MATMUL(TRANSPOSE(oep%w2), u)
      z2 = MATMUL(TRANSPOSE(oep%w2), z)
      nconstraint = 1
      vref = MATMUL(oep%w2, vc2)
      normalization = DOT_PRODUCT(oep%charge, vref)
      IF (ABS(normalization) < 1.0E-12_dp) CPABORT("EXX-OEP reference charge is zero")
      vref = -vref/normalization
      IF (homo_condition) THEN
         nconstraint = 2
         constraints(1, 1) = DOT_PRODUCT(oep%y2, oep%y2)
         constraints(1, 2) = DOT_PRODUCT(oep%y2, vc2)
         constraints(2, 1) = DOT_PRODUCT(z2, oep%y2)
         constraints(2, 2) = DOT_PRODUCT(z2, vc2)
         TARGET = [-1.0_dp, kmo(nocc, nocc)]
         CALL oep_linear_solve(constraints, TARGET)
         vref = MATMUL(oep%w2, TARGET(1)*oep%y2 + TARGET(2)*vc2)
      END IF

      ! Orthogonal complement of the charge and optional HOMO constraints in WII space.
      ALLOCATE (projector(naux, naux), eval(naux))
      projector = 0.0_dp
      DO i = 1, naux
         projector(i, i) = 1.0_dp
      END DO
      q = oep%y2/SQRT(DOT_PRODUCT(oep%y2, oep%y2))
      projector = projector - SPREAD(q, 1, naux)*SPREAD(q, 2, naux)
      IF (homo_condition) THEN
         q = MATMUL(projector, z2)
         normalization = SQRT(DOT_PRODUCT(q, q))
         IF (normalization < 1.0E-12_dp) CPABORT("EXX-OEP charge and HOMO constraints are linearly dependent")
         q = q/normalization
         projector = projector - SPREAD(q, 1, naux)*SPREAD(q, 2, naux)
      END IF
      CALL diamat_all(projector, eval)
      nfree = naux - nconstraint
      IF (nfree < 1) CPABORT("EXX-OEP auxiliary basis is too small for the constraints")
      trans = MATMUL(oep%w2, projector(:, nconstraint + 1:))
      d = MATMUL(TRANSPOSE(trans), bov)
      gram = MATMUL(d, TRANSPOSE(d))
      DEALLOCATE (eval)
      ALLOCATE (eval(nfree))
      CALL diamat_all(gram, eval)
      nret = COUNT(eval >= threshold)
      IF (nret == 0) CPABORT("EXX-OEP THR_FAI_OEP removes every auxiliary direction")
      w = MATMUL(trans, gram(:, nfree - nret + 1:))
      d = MATMUL(TRANSPOSE(w), bov)
      DO i = 1, nocc
         DO a = 1, nvir
            j = (i - 1)*nvir + a
            denom = eps(i) - eps(nocc + a)
            IF (ABS(denom) < 1.0E-6_dp) denom = SIGN(1.0E-6_dp, denom + 1.0E-12_dp)
            weights(j) = 1.0_dp/denom
            delta(j) = kmo(nocc + a, i) - DOT_PRODUCT(bov(:, j), vref)
         END DO
      END DO
      x = MATMUL(d*SPREAD(weights, 1, nret), TRANSPOSE(d))
      rhs = MATMUL(d, weights*delta)
      solution = rhs
      gram = x
      CALL oep_linear_solve(gram, solution)
      residual = MAXVAL(ABS(MATMUL(x, solution) - rhs))
      coefficients = vref + MATMUL(w, solution)
      charge_error = ABS(DOT_PRODUCT(oep%charge, coefficients) + 1.0_dp)
      v = 0.0_dp
      DO p = 1, naux
         v = v + coefficients(p)*oep%three_center(:, :, p)
      END DO
   END SUBROUTINE oep_solve

! **************************************************************************************************
!> \brief Small dense linear solve, with failure reporting.
!> \param a ...
!> \param b ...
! **************************************************************************************************
   SUBROUTINE oep_linear_solve(a, b)
      REAL(KIND=dp), INTENT(INOUT)                       :: a(:, :), b(:)

      INTEGER                                            :: info, n
      INTEGER, ALLOCATABLE                               :: piv(:)

      EXTERNAL dgesv
      n = SIZE(b)
      ALLOCATE (piv(n))
      CALL dgesv(n, 1, a, n, piv, b, n, info)
      IF (info /= 0) CPABORT("EXX-OEP constrained linear system is singular")
   END SUBROUTINE oep_linear_solve
END MODULE qs_oep
