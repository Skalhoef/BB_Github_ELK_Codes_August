! May 2025, The Half Blood Mathematician 
! Adapted from bandstr.f90 to rhonkLOBasis.f90

!BOP
! !ROUTINE: rhonkLOBasis
! !INTERFACE:
subroutine rhonkLOBasis
! !USES:
use modmain
use modomp
! !DESCRIPTION:
!   Prints the matrix \rho_{N k} w.r.t. the [L]ocal [O]rbital Basis 
!   along the path in reciprocal space which connects
!   the vertices in the array {\tt vvlp1d}.
!EOP
!BOC
implicit none
! local variables
integer ik,ist,ispn,is,ia,ias
integer lmax,lmmax,l,m,lm,iv,nthd

! SK: Additions for printing off-diagonal elements
integer lmaxprime, lmmaxprime, lprime, mprime, lmprime
integer ispnprime

! SK: Store results for later writing
integer n_bands_in_window  ! Add this line
complex(8), allocatable :: dmat_all(:,:,:,:,:,:,:)  ! (lmmax,nspinor,lmmax,nspinor,nstsv,natmtot,nkpt)

! SK: 2025-07-14: 
! Here we have an nstsv-dimensional array of booleans. 
! At each entry it is contains either the variable .true. or .false.
! One has enbox_boolean_arr(n) is equal to .true. if band n is fully contained in the energy-window
! [E_F - 0.5 Ha, E_F + 0.5 Ha]
! and it is equal to .false. otherwise.
logical, allocatable :: enbox_boolean_arr(:)


real(8) emin,emax,sm,t1
character(256) fname
! allocatable arrays
real(8), allocatable :: evalfv(:,:),e(:,:)
! low precision for band character array saves memory
real(4), allocatable :: bc(:,:,:,:)
complex(8), allocatable :: dmat(:,:,:,:,:),apwalm(:,:,:,:,:)
complex(8), allocatable :: evecfv(:,:,:),evecsv(:,:)
! initialise universal variables
call init0
call init1

! SK: 2025-07-14:
! Allocate the boolean array after nstsv is properly set
allocate(enbox_boolean_arr(nstsv))

! allocate array for storing the eigenvalues
allocate(e(nstsv,nkpt))
! maximum angular momentum for band character
lmax=min(3,lmaxo)
lmmax=(lmax+1)**2

! SK: 2025-07-14: Starting Configuration: No band crosses the energy window
enbox_boolean_arr = .false.


if (task == 26) then
  allocate(dmat_all(lmmax,nspinor,lmmax,nspinor,nstsv,natmtot,nkpt))
  ! Maybe we should also set each entry of the array here to zero, since 
  ! Fortran might fill it up with random numbers otherwise?
end if
! read density and potentials from file
call readstate
! Fourier transform Kohn-Sham potential to G-space
call genvsig
! read Fermi energy from file
call readfermi
! find the new linearisation energies
call linengy
! generate the APW and local-orbital radial functions and integrals
call genapwlofr
! generate the spin-orbit coupling radial functions
call gensocfr
emin=1.d5
emax=-1.d5
! begin parallel loop over k-points
call holdthd(nkpt,nthd)
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(evalfv,evecfv,evecsv) &
!$OMP PRIVATE(dmat,apwalm,ist,ispn) &
!$OMP PRIVATE(ias,l,m,lm,sm) &
!$OMP NUM_THREADS(nthd)
allocate(evalfv(nstfv,nspnfv))
allocate(evecfv(nmatmax,nstfv,nspnfv),evecsv(nstsv,nstsv))
if (task >= 21) then
  allocate(dmat(lmmax,nspinor,lmmax,nspinor,nstsv))
  allocate(apwalm(ngkmax,apwordmax,lmmaxapw,natmtot,nspnfv))
end if

!$OMP DO
do ik=1,nkpt
!$OMP CRITICAL(bandstr_1)
  write(*,'("Info(bandstr): ",I6," of ",I6," k-points")') ik,nkpt
!$OMP END CRITICAL(bandstr_1)
! solve the first- and second-variational eigenvalue equations
  call eveqn(ik,evalfv,evecfv,evecsv)
  do ist=1,nstsv
! subtract the Fermi energy
    e(ist,ik)=evalsv(ist,ik)-efermi
!$OMP CRITICAL(bandstr_2)
    emin=min(emin,e(ist,ik))
    emax=max(emax,e(ist,ik))
!$OMP END CRITICAL(bandstr_2)
  end do
! compute the band characters if required
  if (task >= 21) then
! find the matching coefficients
    do ispn=1,nspnfv
      call match(ngk(ispn,ik),vgkc(:,:,ispn,ik),gkc(:,ispn,ik), &
       sfacgk(:,:,ispn,ik),apwalm(:,:,:,:,ispn))
    end do
    
    do ias=1,natmtot
      ! SK: generate the density matrix (ALSO OFF-DIAGONAL ELEMENTS!!!)
      call gendmatk(.false.,.false.,0,lmax,ias,nstsv,[0],ngk(:,ik),apwalm,evecfv,&
       evecsv,lmmax,dmat)
      dmat_all(:,:,:,:,:,ias,ik) = dmat(:,:,:,:,:)
    end do
  end if
! end loop over k-points
end do
!$OMP END DO
deallocate(evalfv,evecfv,evecsv)
if (task >= 21) deallocate(dmat,apwalm)
!$OMP END PARALLEL

! Now we determine which of the bands cross the relevant energy-window
do ist=1,nstsv
  ! Check if the band is fully contained in the energy window
  if (all(e(ist,:) >= - 0.5d0) .and. all(e(ist,:) <= + 0.5d0)) then
    enbox_boolean_arr(ist) = .true.
  end if
end do

! Now write everything sequentially outside the parallel region
if (task == 26) then
  ! Step 1/2: Print the bands
  open(50,file='BAND.OUT',form='FORMATTED',action='WRITE')
  do ist=1,nstsv
    ! Check whether the band is fully contained in the energy window
    if (.not. enbox_boolean_arr(ist)) then
      cycle  ! Skip bands not fully contained in the energy window
    end if
    do ik=1,nkpt
      write(50,'(2G18.10)') dpp1d(ik),e(ist,ik)
    end do
    write(50,*)
  end do
  close(50)
  write(*,*)
  write(*,'("Info(bandstr):")')
  write(*,'(" band structure plot written to BAND.OUT")')
  
  ! count the number of bands fully contained in the energy window
  n_bands_in_window = count(enbox_boolean_arr)

  ! Step 2/2: Print the density matrices
  ! Printing the dimensions...
  write(92, '("lmax = ", I4)') lmax
  write(92, '("lmmax = ", I4)') lmmax
  write(92, '("nspinor = ", I4)') nspinor
  write(92, '("n_bands_in_window = ", I4)') n_bands_in_window
  write(92, '("nspecies = ", I4)') nspecies  ! Changed from natmtot
  write(92, '("nkpt = ", I4)') nkpt
  
  do ik=1,nkpt
    write(92, '("k-point = ", I4)') ik
    do is=1,nspecies                    ! Loop over species instead
      ia = 1                            ! Pick first atom of this species
      ias = idxas(ia,is)               ! Get global atom index
      write(92, '("Species = ", I4, " (ias = ", I4, ")")') is, ias
      do ist=1,nstsv
        if (.not. enbox_boolean_arr(ist)) then
          cycle  ! Skip bands not fully contained in the energy window
        end if
        write(92, *) ! Somehow we need another newline here...
        write(92, '("Density matrix for state ", I4)') ist
        lm = 0
        do l = 0, lmax
          do m = -l, l
            lm = lm + 1
            do ispn = 1, nspinor
              lmprime = 0
              do lprime = 0, lmax 
                do mprime = -lprime, lprime 
                  lmprime = lmprime + 1
                  do ispnprime = 1, nspinor 
                    write(92, '(2G18.10)', advance='no') &
                      real(dmat_all(lm,ispn,lmprime,ispnprime,ist,ias,ik)), &
                      aimag(dmat_all(lm,ispn,lmprime,ispnprime,ist,ias,ik))
                  end do
                end do
              end do
            end do
          end do
        end do
        write(92, *)  ! Newline after each matrix
      end do
      write(92, *) ! Newline after each state
    end do 
    write(92, *) ! Newline after each species
  end do
  write(92, *) ! Newline after each k-point
end if

deallocate(dmat_all)

call freethd(nthd)
write(*,*)
write(*,'(" Fermi energy is at zero in plot")')
! output the vertex location lines
open(50,file='BANDLINES.OUT',form='FORMATTED',action='WRITE')
do iv=1,nvp1d
  write(50,'(2G18.10)') dvp1d(iv),emin
  write(50,'(2G18.10)') dvp1d(iv),emax
  write(50,*)
end do
close(50)
write(*,*)
write(*,'(" Vertex location lines written to BANDLINES.OUT")')
deallocate(e)
deallocate(enbox_boolean_arr)  ! Add this line
end subroutine
!EOC

