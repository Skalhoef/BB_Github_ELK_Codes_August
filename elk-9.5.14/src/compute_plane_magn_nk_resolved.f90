
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: bandstr
! !INTERFACE:
subroutine compute_plane_magn_nk_resolved
! !USES:
use modmain

! SK Begin: Insertion of additional modules
use modsebbe
! SK End

! !DESCRIPTION:
!   Produces files MAG2D_n_k.OUT in which m_{ n k } ( r ) is written in the format of task 72.
!
! !REVISION HISTORY:
!   Created June 2003 (JKD)
!EOP
!BOC
implicit none
! local variables
integer ik,ist,ispn,is,ia,ias
integer lmax,lmmax,l,m,lm,iv,nthd

! SK Begin: Insertion of additional variables

! For storing results that we need later
real(8) wppt
real(8), allocatable :: magmt_nk(:,:,:,:,:)  ! (npmtmax,natmtot,ndmag,nstsv,nkpt)
real(8), allocatable :: magir_nk(:,:,:,:)    ! (ngtot, ndmag, nstsv, nkpt)
real(8), allocatable :: occsvp(:,:)
integer idm
! SK End

real(8) emin,emax,sm,t1
character(256) fname
! allocatable arrays
real(8), allocatable :: evalfv(:,:),e(:,:)
! low precision for band character array saves memory
real(4), allocatable :: bc(:,:,:,:)
complex(8), allocatable :: apwalm(:,:,:,:,:)
complex(8), allocatable :: evecfv(:,:,:),evecsv(:,:)
! initialise universal variables
call init0
call init1
! allocate array for storing the eigenvalues
allocate(e(nstsv,nkpt))
! maximum angular momentum for band character
lmax=min(3,lmaxo)
lmmax=(lmax+1)**2


! SK Begin: Insertion of allocations (test)
allocate(magmt_nk(npmtmax,natmtot,ndmag,nstsv,nkpt))
allocate(magir_nk(ngtot,ndmag,nstsv,nkpt))
allocate(occsvp(nstsv,nkpt))
wppt = 1.0d0
! SK End

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
! SK Begin: 
occsvp(:,:)=1.d0
! SK End




! begin parallel loop over k-points
allocate(evalfv(nstfv,nspnfv))
allocate(evecfv(nmatmax,nstfv,nspnfv),evecsv(nstsv,nstsv))
if (task >= 21) then
  allocate(apwalm(ngkmax,apwordmax,lmmaxapw,natmtot,nspnfv))
end if
do ik=1,nkpt
  write(*,'("Info(compute_plane_magn_nk_resolved): ",I6," of ",I6," k-points")') ik,nkpt
! solve the first- and second-variational eigenvalue equations
  call eveqn(ik,evalfv,evecfv,evecsv)
  do ist=1,nstsv
! subtract the Fermi energy
    e(ist,ik)=evalsv(ist,ik)-efermi
    emin=min(emin,e(ist,ik))
    emax=max(emax,e(ist,ik))
  end do

  ! find the matching coefficients
  do ispn=1,nspnfv
    call match(ngk(ispn,ik),vgkc(:,:,ispn,ik),gkc(:,ispn,ik), &
     sfacgk(:,:,ispn,ik),apwalm(:,:,:,:,ispn))
  end do

  call compute_rhomag_nk_Sebbe(ngk(:,ik),igkig(:,:,ik),wppt,occsvp(:,ik),apwalm, &
     evecfv,evecsv,magmt_nk(:,:,:,:,ik),magir_nk(:,:,:,ik),ik)

  do ist=1,nstsv
    call rhomagsh_Sebbe(magmt_nk(:,:,:,ist,ik))

    ! symmetrise the magnetisations
    call symrvf(.true.,ncmag,nrcmt,nrcmti,npcmt,ngdgc,ngtc,ngvc,igfc,npmtmax, &
     magmt_nk(:,:,:,ist,ik),ngtot,magir_nk(:,:,ist,ik))

    ! convert the muffin-tin magnetisations from coarse to fine radial mesh
    ! convert the interstitial magnetisations from coarse to fine grid
    do idm=1,ndmag
      call rfmtctof(magmt_nk(:,:,idm,ist,ik))
      call rfirctof(magir_nk(:,idm,ist,ik),magir_nk(:,idm,ist,ik))
    end do
    
    call vecplot_nk_Sebbe(ist,ik,magmt_nk(:,:,:,ist,ik),magir_nk(:,:,ist,ik))
  end do 

! end loop over k-points
end do
deallocate(evalfv,evecfv,evecsv)


if (task >= 21) deallocate(apwalm)

! SK Begin: deallocate the arrays
deallocate(magmt_nk,magir_nk)
! SK End


t1=(emax-emin)*0.5d0
emin=emin-t1
emax=emax+t1



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
end subroutine
!EOC