!  Solves 2D periodic Navier Stokes equation using spectral methods
!  FFT using fftw library from REAL physical space to complex spectral space
!  Written by Hossein A Kafiabad fall 2017
!  time-stepping using forward Euler, 3-step Adams-Bashforth, and Crank-Nicolson (for the moment)

! the equation that is solved:
! zz : vorticity
! Dzz/Dt = \nu.\nabla^n(zz) (n = 2 for the moment, newtonian dissipation)
! psi stream function ---> \nabla^2(psi) = -zz
!  u =   psi_y zonal velocity (along x-axis)
!  v = - psi_x meridional velocity (along y-axis)


PROGRAM MAIN
  use, intrinsic :: iso_c_binding       ! Use the intrinsic FORTRAN module "iso_c_binding," which supports C interoperability
  use param_John                        ! Use the FORTRAN file "param_John.f90," which initializes constants
  implicit none                         ! Keep FORTRAN from treating "i," "j," "k," "l," "m," and "n" as integers with every other variable treated as a real
  include "fftw3.f03"                  ! Use the FORTRAN fast Fourier transform library

! ==========================================================
! ----------------------------------------------------------
! +++++ variable declarations +++++

! tint = time interval over which numerical method runs
! tstop = final time used
! tstart = initial time

! icread = UNUSED
! icrk = initial condition in REAL SPACE

! visc = molecular viscosity in the 2D Navier-Stokes equation

! ikx = index from 1 to (size of mesh in x-direction)/2 + 1 (see param_John for size of mesh)
! iky = index from 1 to size of mesh in y-direction (see param_John for size of mesh)
! jj = variable used in calculating y-component of wave number vector
! i = index variable only used in subroutines
! j = index variable only used in subroutines
! k = index variable for output file

! kx = variable storing wave numbers created in kxa (UNNEEDED?)
! ky = variable storing wave numbers created in kya (UNNEEDED?)
! kh = square root of dissipation term k^2
! k2 = dissipation term k^2
! Lino = stand in variable for the term in Crank-Nicolson nu*k2/2
! xx = index variable only used in subroutines
! yy = index variable only used in subroutines

! plan2_ur_uk = FORTRAN plan (C pointer) that contains all the information necessary to compute the fft transform of "u" component of velocity from real to complex
! plan2_uk_ur = FORTRAN plan (C pointer) that stores information for fft transform of "u" component of velocity from complex to real
! plan2_vr_vk = FORTRAN plan (C pointer) that stores information for fft transform of "v" component of velocity from real to complex
! plan2_vk_vr = FORTRAN plan (C pointer) that stores information for fft transform of "v" component of velocity from complex to real
! plan2_zzr_zzk = FORTRAN plan (C pointer) that stores information for fft transform of vorticity from real to complex
! plan2_zzk_zzr = FORTRAN plan (C pointer) that stores information for fft transform of vorticity from complex to real
! plan2_nur_nuk = FORTRAN plan (C pointer) that stores information for fft transform of new "u" component of velocity from real to complex
! plan2_nvr_nvk = FORTRAN plan (C pointer) that stores information for fft transform of new "v" component of velocity from real to complex

! u_p = pointer referring to "u" component of velocity (used as a memory location for allocating fft array)
! v_p = pointer referring to "v" component of velocity (used as a memory location for allocating fft array)
! zz_p = pointer referring to vorticity (used as a memory location for allocating fft array)
! nu_p = pointer referring to new "u" component of velocity (used as a memory location for allocating fft array)
! nv_p = pointer referring to new "v" component of velocity (used as a memory location for allocating fft array)

! ur = real array of "u" component of velocity going into fft
! uk = complex array of "u" component of velocity coming out of fft
! vr = real array of "v" component of velocity going into fft
! vk = complex array of "v" component of velocity coming out of fft
! zzr = real array of vorticities going into fft
! zzk = complex array of vorticities coming out of fft
! nur = real array of updated velocity "u" component
! nuk = complex array of updated velocity "u" component
! nvr = real array of updated velocity "v" component
! nvk = complex array of updated velocity "v" component

! nnk = nonlinear term in numerical method
! zzknew = updated vorticity from Crank-Nicholson plus nonlinear term

! fk = UNUSED variable input to UNUSED subroutine "WriteFieldfk"
! fr = UNUSED variable input to UNUSED subroutine "WriteFieldfr"

! ----------------------------------------------------------

! time
  real                 :: tint = 0.5, tstop, tstart         ! Initialize time variables

! specify numerical method
  character (len = 20), parameter       :: method = "Crank"

! related to initial condition
  real, parameter      :: icread = 0                        ! reads the initial condition from an ext. file (UNUSED)
  real, parameter      :: icrk   = 0                        ! 0: I.C. set in real space (UNUSED)
                                                            ! 1: I.C. set in spectral space
! dissipation
  real, parameter      :: visc = 1e-5                       ! molecular viscosity


! working variables
  integer              :: ikx, iky, jj, i, j, k             ! Initialize index variables in x- and y-directions
  real                 :: kx, ky, kh, k2, Lino, xx, yy      ! Initialize variables storing wave numbers, dissipation term, square root of dissipation term, and variable storing term in Crank-Nicholson

! set the FFTW plans
  type(C_PTR), save    :: plan2_ur_uk, plan2_uk_ur, plan2_vr_vk, plan2_vk_vr                                ! Initialize fft plans for velocity components
  type(C_PTR), save    :: plan2_zzr_zzk, plan2_zzk_zzr, plan2_nur_nuk, plan2_nvr_nvk                        ! Initialize fft plans for vorticity
  type(C_PTR), save    :: u_p, v_p, zz_p, nu_p, nv_p                                                        ! Initialize pointers used as memory locations for fft arrays

! the fields whose Fourier transform is calculated
  real(C_DOUBLE), dimension(:,:), pointer          :: ur, vr, zzr, nur, nvr                                 ! Initialize pointers for storing REAL velocities and vorticities
  complex(C_DOUBLE_COMPLEX),dimension(:,:),pointer :: uk, vk, zzk, nuk, nvk                                 ! Initialize pointers for storing COMPLEX velocities and vorticities

! the fields whose Fourier transform is NEVER calculated
  complex(C_DOUBLE_COMPLEX), dimension(iktx,ikty)  :: nnk, nnk_nm1, nnk_nm2, zzknew, zzk_nm1, zzk_nm2       ! Initialize variables for Euler, Crank-Nicholson, and 3-step Adams-Bashforth

  ! temp
  complex(C_DOUBLE_COMPLEX), dimension(iktx,ikty)  :: fk                                                    ! UNUSED variable for diagnostics
  real(C_DOUBLE),    dimension(Nx,Ny)              :: fr                                                    ! UNUSED variable for diagnostics

! ==========================================================
! ----------------------------------------------------------
! +++++ Allocations and initializations +++++
! ----------------------------------------------------------

! allocate fftw arrays

  u_p = fftw_alloc_complex(int(iktx * ikty, C_SIZE_T))              ! Allocate memory to store array of complex numbers for fft of "u" component of velocity
  call c_f_pointer(u_p, ur, [2*(Nx/2+1) , Ny])                      ! Create pointer for storing 2D array "ur" with size Nx+2 by Ny
  call c_f_pointer(u_p, uk, [iktx , ikty])                          ! Create pointer for storing 2D array "uk" with size Nx/2+1 by Ny
  v_p = fftw_alloc_complex(int(iktx * ikty, C_SIZE_T))              ! Allocate memory to store array of complex numbers for fft of "v" component of velocity
  call c_f_pointer(v_p, vr, [2*(Nx/2+1) , Ny])                      ! Create pointer for storing 2D array "vr" with size Nx+2 by Ny
  call c_f_pointer(v_p, vk, [iktx , ikty])                          ! Create pointer for storing 2D array "vk" with size Nx/2+1 by Ny
  zz_p = fftw_alloc_complex(int(iktx * ikty, C_SIZE_T))             ! Allocate memory to store array of complex numbers for fft of vorticity
  call c_f_pointer(zz_p, zzr, [2*(Nx/2+1) , Ny])                    ! Create pointer for storing 2D array "zzr" with size Nx+2 by Ny
  call c_f_pointer(zz_p, zzk, [iktx , ikty])                        ! Create pointer for storing 2D array "zzk" with size Nx/2+1 by Ny
  nu_p = fftw_alloc_complex(int(iktx * ikty, C_SIZE_T))             ! Allocate memory to store array of complex numbers for fft of "nu"
  call c_f_pointer(nu_p, nur, [2*(Nx/2+1) , Ny])                    ! Create pointer for storing 2D array "nur" with size Nx+2 by Ny
  call c_f_pointer(nu_p, nuk, [iktx , ikty])                        ! Create pointer for storing 2D array "nuk" with size Nx/2+1 by Ny
  nv_p = fftw_alloc_complex(int(iktx * ikty, C_SIZE_T))             ! Allocate memory to store array of complex numbers for fft of "nv"
  call c_f_pointer(nv_p, nvr, [2*(Nx/2+1) , Ny])                    ! Create pointer for storing 2D array "nvr" with size Nx+2 by Ny
  call c_f_pointer(nv_p, nvk, [iktx , ikty])                        ! Create pointer for storing 2D array "nvk" with size Nx/2+1 by Ny

! ----------------------------------------------------------

! Initialize FFTW plans (consider dimension are passed reversely
! FFTW_MEASURE tells machine to find an optimized fft by taking multiple ffts of the problem and measuring the computational time required by each

  plan2_ur_uk = fftw_plan_dft_r2c_2d(Ny,Nx,ur,uk, FFTW_MEASURE)         ! Initialize the plan for fft of real vector of "u" to complex vector
  plan2_uk_ur = fftw_plan_dft_c2r_2d(Ny,Nx,uk,ur, FFTW_MEASURE)         ! Initialize the plan for fft of complex vector of "u" to real vector
  plan2_vr_vk = fftw_plan_dft_r2c_2d(Ny,Nx,vr,vk, FFTW_MEASURE)         ! Initialize the plan for fft of real vector of "v" to complex vector
  plan2_vk_vr = fftw_plan_dft_c2r_2d(Ny,Nx,vk,vr, FFTW_MEASURE)         ! Initialize the plan for fft of complex vector of "v" to real vector
  plan2_zzr_zzk = fftw_plan_dft_r2c_2d(Ny,Nx,zzr,zzk, FFTW_MEASURE)     ! Initialize the plan for fft of real vector of vorticity to complex vector
  plan2_zzk_zzr = fftw_plan_dft_c2r_2d(Ny,Nx,zzk,zzr, FFTW_MEASURE)     ! Initialize the plan for fft of complex vector of vorticity to real vector
  plan2_nur_nuk = fftw_plan_dft_r2c_2d(Ny,Nx,nur,nuk, FFTW_MEASURE)     ! Initialize the plan for fft of real vector of "nu" to complex vector
  plan2_nvr_nvk = fftw_plan_dft_r2c_2d(Ny,Nx,nvr,nvk, FFTW_MEASURE)     ! Initialize the plan for fft of complex vector of "nv" to real vector

! ==========================================================
! ----------------------------------------------------------
! +++++ Initialize wavenumbers and truncation mask +++++
! ----------------------------------------------------------

! Here we initialize the 2D wave number vector "k_hat" with components in x- and y-directions
! NOTE: 2*pi/Lx = 1 and 2*pi/Ly = 1

  do  ikx = 1,iktx                              ! This "do" loop runs from 1 to the total number of x-component wavenumbers, which is Nx/2 + 1, with Nx/2 + 1 being the largest index of the wavenumbers
     kxa(ikx) = float(ikx-1)*2*pi/Lx            ! Store the initial x-component of the wavenumber in the vector "kxa"
  end do

  do iky = 1,ikty                               ! This "do" loop runs from 1 to the total number of y-component wavenumbers, which is Ny
     jj = iky - 1                               ! "jj" is always one less than the current value of "iky"
     if (iky.gt.kty)   jj = jj - 2*kty          ! If iky > kty, where kty = Ny/2, then "jj" is negative
     if (iky.eq.kty+1) jj = 0                   ! If iky = kty+1, then "jj" is 0
     if (iky.gt.2*kty) jj = 0                   ! If iky > 2*kty (greater than Ny), then "jj" is 0
     kya(iky) = float(jj)*2*pi/Ly               ! Store the initial y-component of the wavenumber in the vector "kya"
  end do

  L  = 1                                        ! Initialize mask for wavenumber truncation to a 2D array of ones with dimensions (Nx/2+1) by Ny
  do iky = 1,ikty                               ! This "do" loop runs from 1 to the total number of y-component wavenumbers, which is Ny
     ky = kya(iky)                              ! Store initial y-component of wavenumber in real variable "ky"

     do  ikx = 1,iktx                           ! This "do" loop runs from 1 to total number of x-component wavenumbers
        kx = kxa(ikx)                           ! Store initial x-component of the wavenumber in real variable "kx"
        kh = sqrt(kx*kx + ky*ky)                ! Store magnitude of wavenumber in "kh"
        if (iky.eq.kty+1)                         L(ikx,iky) = 0                            ! If iky = kty+1 truncation mask value is 0
        if (kx.eq.0 .and. ky.le.0)                L(ikx,iky) = 0                            ! If kx = 0 AND ky <= 0 truncation mask value is 0
        ! truncation with 2/3 instead of 8/9 truncation
        if (abs(kx*Lx/2.0/pi).gt.int(float(Nx)*1./3. + 0.5)-0.5)    L(ikx,iky) = 0          ! If abs((kx*Lx)/(2*pi)) > int((Nx*4)/9 + 0.5) - 0.5 truncation mask value is 0
        if (abs(ky*Ly/2.0/pi).gt.int(float(Ny)*1./3. + 0.5)-0.5)    L(ikx,iky) = 0          ! If abs((ky*Ly)/(2*pi)) > int((Ny*4)/9 + 0.5) - 0.5 truncation mask value is 0
     end do

  end do

! ----------------------------------------------------------
! produce initial condition

  if (icrk.eq.0) then                               ! "If" loop always is TRUE because icrk = 0
     ! I.C. described in real space
     call initR(fr,tstart)                          ! Call the subroutine "initR" with formal arguments "fr" and "tstart"
  else
     ! I.C. described in spectral space
     ! call initK(zzk,tstart)
  end if

  zzr = fr                                          ! Set initial vorticity matrix equal to the initial condition "fr"
  call fftw_execute_dft_r2c(plan2_zzr_zzk,zzr,zzk)  ! Use fft to transform initial condition from real to spectral space: input is "zzr" and output is "zzk"
  zzk=zzk/fftnorm                                   ! Normalize vorticity

! ----------------------------------------------------------
! 3-step Adams-Bashforth starting procedure

  if (method.eq."AB") then

    call wave_numbers(uk, vk, zzk)                                  ! Return components of velocity in Fourier space based on the wavenumbers

    call fftw_execute_dft_c2r(plan2_zzk_zzr,zzk,zzr)                ! Use fft transforming from complex to real vorticities
    call fftw_execute_dft_c2r(plan2_uk_ur,uk,ur)                    ! Use fft transforming from complex to real "u" component of velocity
    call fftw_execute_dft_c2r(plan2_vk_vr,vk,vr)                    ! Use fft transforming from complex to real "v" component of velocity

    nur = ur*zzr                                                    ! New real "u" component
    nvr = vr*zzr                                                    ! New real "v" component

    call fftw_execute_dft_r2c(plan2_nur_nuk,nur,nuk)                ! fft transforming from real to complex "nu"
    nuk=nuk/fftnorm                                                 ! Normalizing "nuk"
    call fftw_execute_dft_r2c(plan2_nvr_nvk,nvr,nvk)                ! fft transforming from real to complex "nv"
    nvk=nvk/fftnorm                                                 ! Normalizing "nvk"
    call fftw_execute_dft_r2c(plan2_zzr_zzk,zzr,zzk)                ! fft transforming from real to complex "zz"
    zzk=zzk/fftnorm                                                 ! Normalizing "zzk"

    do iky = 1, ikty                                                ! This "do" loop runs from 1 to Ny
        ky = kya(iky)                                               ! Set y-component of wavenumbers

        do ikx = 1, iktx                                            ! This "do" loop runs from 1 to Nx/2+1
            kx = kxa(ikx)                                           ! Set x-component of wavenumbers
            k2 = kx*kx + ky*ky                                      ! Calculate squared magnitude of wave number k^2
            nnk(ikx,iky) = - zi * (kx * nuk(ikx,iky) + ky * nvk(ikx,iky))* L(ikx,iky)       ! Calculating "nnk" to use in numerical methods
            Lino =  0.5 * visc * k2                                                         ! Calculating "nu*k^2/2" to use in numerical methods
            zzknew(ikx, iky) = zzk(ikx, iky) + delt*(-nnk(ikx, iky) - Lino*zzk(ikx,iky))    ! Calculating updated vorticity using forward Euler
        end do

    end do

    time = tstart + delt                                            ! Add one timsestep

    zzk_nm1 = zzk                                                   ! Store vorticity values from timestep n-1 in zzk_nm1
    zzk = zzknew                                                    ! Store vorticity values from timestep n in zzk
    nnk_nm1 = nnk                                                   ! Store nonlinear terms from timestep n-1 in nnk_nm1

    call wave_numbers(uk, vk, zzk)                                  ! Calculate new velocity components based on new vorticities

    call fftw_execute_dft_c2r(plan2_zzk_zzr,zzk,zzr)                ! Use fft transforming from complex to real vorticities
    call fftw_execute_dft_c2r(plan2_uk_ur,uk,ur)                    ! Use fft transforming from complex to real "u" component of velocity
    call fftw_execute_dft_c2r(plan2_vk_vr,vk,vr)                    ! Use fft transforming from complex to real "v" component of velocity

    nur = ur*zzr                                                    ! New real "u" component
    nvr = vr*zzr                                                    ! New real "v" component

    call fftw_execute_dft_r2c(plan2_nur_nuk,nur,nuk)                ! fft transforming from real to complex "nu"
    nuk=nuk/fftnorm                                                 ! Normalizing "nuk"
    call fftw_execute_dft_r2c(plan2_nvr_nvk,nvr,nvk)                ! fft transforming from real to complex "nv"
    nvk=nvk/fftnorm                                                 ! Normalizing "nvk"
    call fftw_execute_dft_r2c(plan2_zzr_zzk,zzr,zzk)                ! fft transforming from real to complex "zz"
    zzk=zzk/fftnorm                                                 ! Normalizing "zzk"

    do iky = 1, ikty
        ky = kya(iky)

        do ikx = 1, iktx
            kx = kxa(ikx)
            k2 = kx*kx + ky*ky
            nnk(ikx,iky) = - zi * (kx * nuk(ikx,iky) + ky * nvk(ikx,iky))* L(ikx,iky)               ! Calculating "nnk" to use in numerical methods
            Lino =  0.5 * visc * k2                                                                 ! Calculating "nu*k^2/2" to use in numerical methods
            zzknew(ikx, iky) = zzk(ikx, iky) + delt*(3./2.)*(-Lino*zzk(ikx, iky) - nnk(ikx, iky)) - delt*(1./2.)*(-Lino*zzk_nm1(ikx, iky) - nnk_nm1(ikx, iky))           ! Calculating updated vorticity using 2-step Adams-Bashforth
        end do

    end do

    time = time + delt      ! Add a timestep
    zzk_nm2 = zzk_nm1       ! Store vorticities from time n-2 in zzk_nm2
    zzk_nm1 = zzk           ! Store vorticities from time n-1 in zzk_nm1
    zzk = zzknew            ! Store vorticities from tine n in zzk
    nnk_nm2 = nnk_nm1       ! Store nonlinear terms from time n-2 in nnk_nm2
    nnk_nm1 = nnk           ! Store nonlinear terms from time n-1 in nnk_nm1

  else

    time = tstart           ! Set first time to starting time

  end if
! ==========================================================
! ----------------------------------------------------------
! +++++ Executions +++++
! ----------------------------------------------------------

  tstop = tstart + tint                                         ! Initialize stopping time

  do while (time.le.tstop)                                      ! This "do-while" loop runs while time <= tstop
    time = time + delt                                          ! Increment time by one timestep

! doing one step of time-stepping
  call wave_numbers(uk, vk, zzk)

  call fftw_execute_dft_c2r(plan2_zzk_zzr,zzk,zzr)              ! Use fft transforming from complex to real vorticities
  call fftw_execute_dft_c2r(plan2_uk_ur,uk,ur)                  ! Use fft transforming from complex to real "u" component of velocity
  call fftw_execute_dft_c2r(plan2_vk_vr,vk,vr)                  ! Use fft transforming from complex to real "v" component of velocity

  nur = ur*zzr                                                  ! New real "u" component
  nvr = vr*zzr                                                  ! New real "v" component

  call fftw_execute_dft_r2c(plan2_nur_nuk,nur,nuk)              ! fft transforming from real to complex "nu"
  nuk=nuk/fftnorm                                               ! Normalizing "nuk"
  call fftw_execute_dft_r2c(plan2_nvr_nvk,nvr,nvk)              ! fft transforming from real to complex "nv"
  nvk=nvk/fftnorm                                               ! Normalizing "nvk"
  call fftw_execute_dft_r2c(plan2_zzr_zzk,zzr,zzk)              ! fft transforming from real to complex "zz"
  zzk=zzk/fftnorm                                               ! Normalizing "zzk"

  do iky = 1,ikty                                                                           ! This "do" loop runs from 1 to size of grid in y-direction
     ky = kya(iky)                                                                          ! Set "ky" equal to y-component of wavenumber

     do  ikx = 1,iktx                                                                       ! This "do" loop runs from 1 to size of grid in x-direction
        kx = kxa(ikx)                                                                       ! Set "kx" equal to x-component of wavenumber
        k2 = kx*kx + ky*ky                                                                  ! "k^2" used in numerical methods
        ! nnk is in the RHS of equation. Hence F^-1(nnk) = - ( u z_x + v z_y)
        nnk(ikx,iky) = - zi * (kx * nuk(ikx,iky) + ky * nvk(ikx,iky))* L(ikx,iky)           ! Calculating "nnk" to use in numerical methods
        Lino =  0.5 * visc * k2                                                             ! Calculating "nu*k^2/2" to use in numerical methods

        if (method.eq."Forward") then                                                       ! This "if" loop determines what equation to used based on the desired numerical method
            zzknew(ikx, iky) = zzk(ikx, iky) + delt*(-nnk(ikx, iky) - Lino*zzk(ikx,iky))    ! Calculating updated vorticity using forward Euler
        else if (method.eq."AB") then
            zzknew(ikx, iky) = zzk(ikx, iky) + delt*23./12.*(-Lino*zzk(ikx, iky) - nnk(ikx, iky)) - &                               ! Calculating updated vorticity using 3-step Adams-Bashforth
            delt*4./3.*(-Lino*zzk_nm1(ikx, iky) - nnk_nm1(ikx, iky)) + delt*5./12.*(-Lino*zzk_nm2(ikx, iky) - nnk_nm2(ikx, iky))
        else
            zzknew(ikx,iky)=((1/delt - Lino)*zzk(ikx,iky)-nnk(ikx,iky))/(1/delt + Lino) &   ! Calculating updated vorticity using Crank-Nicholson
                            * L(ikx,iky)
        end if

     end do

  end do

    if ((method.eq."Forward").or.(method.eq."Crank")) then
        zzk = zzknew
    else
        zzk_nm2 = zzk_nm1
        zzk_nm1 = zzk
        nnk_nm2 = nnk_nm1
        nnk_nm1 = nnk
        zzk = zzknew
    end if                                                                          ! Setting vorticity value to updated vorticity value to be used in next iteration of Crank-Nicholson

  end do

  call fftw_execute_dft_c2r(plan2_zzk_zzr,zzk,zzr)                                      ! Finally, perform fft from complex to real vorticity
  print*, 'zzr(1:4,1:4) = ', zzr(1:4,1:4)
  print*,  'time = ', time

  ! open file
  open (1, file='output_file.csv', status='unknown')

  do k = 1, Nx
    write (1,*) zzr(k, :)
  end do

  ! close file
  close(1)

! ==========================================================
! ----------------------------------------------------------
! +++++ Clean-up and deallocations +++++
! ----------------------------------------------------------

! Destroy fftw plans and free memory spaces
  call fftw_destroy_plan(plan2_ur_uk)
  call fftw_destroy_plan(plan2_uk_ur)
  call fftw_destroy_plan(plan2_vr_vk)
  call fftw_destroy_plan(plan2_vk_vr)
  call fftw_destroy_plan(plan2_zzr_zzk)
  call fftw_destroy_plan(plan2_zzk_zzr)
  call fftw_destroy_plan(plan2_nur_nuk)
  call fftw_destroy_plan(plan2_nvr_nvk)
  call fftw_free(u_p)
  call fftw_free(v_p)
  call fftw_free(zz_p)
  call fftw_free(nu_p)
  call fftw_free(nv_p)

  contains

    ! This subroutine calculates velocity components in Fourier space
    subroutine wave_numbers(uk, vk, zzk)

        use param_John                                                                          ! Use the FORTRAN file "param_John.f90," which initializes constants
        implicit none                                                                           ! Keep FORTRAN from treating "i," "j," "k," "l," "m," and "n" as integers with every other variable treated as a real

        integer :: ikx, iky                                                                     ! Initialize index variables
        real :: kx, ky, k2                                                                      ! Initialize variables for holding x- and y-components of wavenumbers as well as k^2
        complex(C_DOUBLE_COMPLEX), intent(in), dimension(:,:), pointer :: zzk                   ! Initialize variables to store components and magnitude squared of wavenumbers
        complex(C_DOUBLE_COMPLEX), intent(out), dimension(:,:), pointer :: uk, vk               ! Initialize pointers for storing COMPLEX velocities and vorticities


        do iky = 1,ikty                                                                         ! This "do" loop runs from 1 to total number of y-component wavenumbers
            ky = kya(iky)                                                                       ! "ky" is set equal to current y-component wavenumber

            do ikx = 1,iktx                                                                     ! This "do" loop runs from 1 to total number of x-component wavenumbers
                kx = kxa(ikx)                                                                   ! "kx" is set equal to the current x-component wavenumber
                k2 = kx*kx + ky*ky                                                              ! k^2 used in Crank-Nicholson
                if (L(ikx,iky).eq.1) then                                                       ! If truncation mask L = 1, then do stuff
                    uk(ikx,iky) =  zi * ky * zzk(ikx,iky) / k2                                  ! Only complex part of "u" component of velocity calculated here
                    vk(ikx,iky) = -zi * kx * zzk(ikx,iky) / k2                                  ! Only complex part of "v" component of velocity calculated here
                else
                    uk(ikx,iky) = 0                                                             ! Set complex part of "u" component of velocity to zero otherwise
                    vk(ikx,iky) = 0                                                             ! Set complex part of "v" component of velocity to zero otherwise
                end if
            end do

        end do

    end subroutine wave_numbers

END PROGRAM MAIN

! ==========================================================
! ----------------------------------------------------------
! +++++ SUBROUTINES  +++++
! ----------------------------------------------------------

! **********************************************************
! +++++ initialization subroutines +++++
!       initR, wavenumbers
! **********************************************************

! This subroutine derives the initial condition. The inputs are the formal arguments "zzr" and "tstart." It returns these values.
subroutine initR(zzr,tstart)

! derive the initial condition described in REAL space
  use param_John                                                    ! Use the FORTRAN file "param_John.f90," which initializes constants
  implicit none                                                     ! Keep FORTRAN from treating "i," "j," "k," "l," "m," and "n" as integers with every other variable treated as a real

  real(8), intent(out), dimension(Nx,Ny) :: zzr                     ! Initialize double precision real vorticity "zzr" with the intent that it receives its value from outside the subroutine
  real,    intent(out)                   :: tstart                  ! Initialize starting time "tstart" with the intent that it receives its value form outside the subroutine
  real                                   :: xx,yy,uu,vv,u_y,v_x     ! Initialize variables for computations
  integer                                :: i,j                     ! Initialize index variables

  do i = 1, Nx                                                      ! This "do" loop runs from 1 to size of physical grid in x-direction

     do j = 1, Ny                                                   ! This "do" loop runs from 1 to size of physical grid in y-direction
        xx = (float(i-1)/Nx)*Lx                                     ! Calculate "xx," which goes from 0 to domain size in x-direction (Lx = 2*pi = domain size in physical space)
        yy = (float(j-1)/Ny)*Ly                                     ! Calculate "yy," which goes form 0 to domain size in y-direction (Ly = 2*pi = domain size in physical space)
        uu=sin(xx)*cos(yy)                                          ! Calculate "uu," which is the initial condition for the "u" component of velocity
        vv=-cos(xx)*sin(yy)                                         ! Calculate "vv," which is the initial condition for the "v" component of velocity
        u_y = -sin(xx)*sin(yy)                                      ! Take partial with respect to "y" of "u" component of velocity
        v_x =  sin(xx)*sin(yy)                                      ! Take partial with respect to "x" of "v" component of velocity
        zzr(i,j)=v_x-u_y                                            ! Calculate initial vorticity values "zzr"
     end do

  end do
  tstart = 0.0                                                      ! "tstart" always set to 0.0

end subroutine initR

! **********************************************************
! +++++ diagnostics  +++++
!       WriteFieldk, WriteFieldr
! **********************************************************

! THIS SUBROUTINE IS NOT USED IN THE CODE
subroutine WriteFieldk(fk)

  use param_John
  implicit none

  complex(8), intent(in), dimension(iktx,ikty) :: fk
  integer              :: i,j
  write(6,'(A12)',advance="no") ' ky|  kx -->'
  write(6,5004) (kxa(i),i=1,iktx)
  write(6,'(A12)') ' v        '
  do j = 1, ikty
     write(6,'(F6.2,"   ")',advance="no") kya(j)
     write(6,'(*(F8.2,"+",F6.2,"i"))') (fk(i,j),i=1,iktx)
  end do
  print*
5004 format(*(F8.2,8x))
end subroutine WriteFieldk

! THIS SUBROURTINE IS NOT USED IN THE CODE
subroutine WriteFieldr(fr)
  use param_John
  implicit none

  real(8), intent(in), dimension(Nx,Ny) :: fr
  real                 :: xx(Nx), yy(Ny)
  integer              :: i,j

  ! grid points
  do i = 1, Nx
     xx(i) = (float(i-1)/Nx)*Lx
  end do
  do j = 1, Ny
     yy(j) = (float(j-1)/Ny)*Ly
  end do

  write(6,'(A12)',advance="no") ' yy|  xx -->'
  write(6,5004) (xx(i),i=1,Nx)
  write(6,'(A12)') ' v        '
  do j = 1, Ny
     write(6,'(F6.2,"   ")',advance="no") yy(j)
     write(6,'(*(F9.3))') (fr(i,j),i=1,Nx)
  end do
  print*
5004 format(*(F5.2,4x))
end subroutine WriteFieldr
