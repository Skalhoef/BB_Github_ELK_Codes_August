! Copyright (C) 2002-2010 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: compute_rhomag_nk_Sebbe
! !INTERFACE:
subroutine compute_rhomag_nk_Sebbe(ngp,igpig,wppt,occsvp,apwalm,evecfv,evecsv,magmt_n,magir_n)
! !USES:
use modmain
! use modsebbe 
! !INPUT/OUTPUT PARAMETERS:
!   ngp    : number of G+p-vectors (in,integer(nspnfv))
!   igpig  : index from G+p-vectors to G-vectors (in,integer(ngkmax,nspnfv))
!   wppt   : weight of input p-point (in,real)
!   occsvp : occupation number for each state (in,real(nstsv))
!   apwalm : APW matching coefficients
!            (in,complex(ngkmax,apwordmax,lmmaxapw,natmtot,nspnfv))
!   evecfv : first-variational eigenvectors (in,complex(nmatmax,nstfv,nspnfv))
!   evecsv : second-variational eigenvectors (in,complex(nstsv,nstsv))
! !DESCRIPTION:
!   TBD
!
!EOP
!BOC
implicit none
! arguments
integer, intent(in) :: ngp(nspnfv),igpig(ngkmax,nspnfv)
real(8), intent(in) :: wppt,occsvp(nstsv)
complex(8), intent(in) :: apwalm(ngkmax,apwordmax,lmmaxapw,natmtot,nspnfv)
complex(8), intent(in) :: evecfv(nmatmax,nstfv,nspnfv),evecsv(nstsv,nstsv)

! SK Begin: Insertion of additional variables
! For storing results that we need later
real(8), intent(inout) :: magmt_n(npmtmax,natmtot,ndmag,nstsv)  ! (lmmaxo,nrmtmax,natmtot,ndmag,nstsv)
real(8), intent(inout) :: magir_n(ngtot,ndmag,nstsv)    ! (ngtot, ndmag, nstsv)
real(8), allocatable :: rhomt_n(:,:,:)  ! (npmtmax,natmtot,nstsv)
real(8), allocatable :: rhoir_n(:,:)    ! (ngtot, nstsv)
! real(8) :: rhomt_n(npmtmax,natmtot,nstsv)  ! (npmtmax,natmtot,nstsv)
! real(8) :: rhoir_n(ngtot,nstsv)            ! (ngtot, nstsv)
! SK End

! local variables
integer ispn,jspn,nst,ist
integer is,ias,npc,i,j,k
integer n,igp
real(8) wo,ts0,ts1
complex(8) z1
! automatic arrays
integer idx(nstsv)



! complex(8) wfir(ngtc,nspinor),wfgp(ngkmax)
! allocatable arrays
complex(8), allocatable :: wfmt(:,:,:)
complex(8), allocatable :: wfir(:,:), wfgp(:)

allocate(wfir(ngtc,nspinor), wfgp(ngkmax))
allocate(rhomt_n(npmtmax,natmtot,nstsv))
allocate(rhoir_n(ngtot,nstsv))


call timesec(ts0)
magmt_n(:,:,:,:)=0.d0
magir_n(:,:,:)=0.d0
rhomt_n(:,:,:)=0.d0
rhoir_n(:,:)=0.d0
!----------------------------------------------!
!     muffin-tin density and magnetisation     !
!----------------------------------------------!

! number of and index to states that we are considering
nst=0
do ist=1,nstsv
  ! if (abs(occsvp(ist)) < epsocc) cycle ! SK: Legacy-Code removed
  nst=nst+1
  idx(nst)=ist
end do

allocate(wfmt(npcmtmax,nspinor,nst))


do ias=1,natmtot
  is=idxis(ias)
  npc=npcmt(is)
  call wfmtsv(.false.,lradstp,is,ias,nst,idx,ngp,apwalm,evecfv,evecsv,npcmtmax,&
   wfmt)

  do j=1,nst
    k=idx(j)
    wo=occsvp(k)*wppt
    
! add to density and magnetisation
    if (spinpol) then
    ! spin-polarised
      if (ncmag) then
      ! non-collinear
        call rmk1(npc,wo,wfmt(:,1,j),wfmt(:,2,j),rhomt_n(:,ias,j),magmt_n(:,ias,1,j), &
         magmt_n(:,ias,2,j),magmt_n(:,ias,3,j))
      end if
    end if
  end do
end do
!------------------------------------------------!
!     interstitial density and magnetisation     !
!------------------------------------------------!
do j=1,nst
  k=idx(j)
  wo=occsvp(k)*wppt/omega
  if (tevecsv) then
! generate spinor wavefunction from second-variational eigenvectors
    do ispn=1,nspinor
      jspn=jspnfv(ispn)
      n=ngp(jspn)
      wfgp(1:n)=0.d0
      do ist=1,nstfv
        i=(ispn-1)*nstfv+ist
        z1=evecsv(i,k)
        if (abs(dble(z1))+abs(aimag(z1)) > epsocc) then
          wfgp(1:n)=wfgp(1:n)+z1*evecfv(1:n,ist,jspn)
        end if
      end do
      wfir(:,ispn)=0.d0
      do igp=1,n
        wfir(igfc(igpig(igp,jspn)),ispn)=wfgp(igp)
      end do
! Fourier transform wavefunction to real-space
      call zfftifc(3,ngdgc,1,wfir(:,ispn))
    end do
  else
! spin-unpolarised wavefunction
    wfir(:,1)=0.d0
    do igp=1,ngp(1)
      wfir(igfc(igpig(igp,1)),1)=evecfv(igp,k,1)
    end do
    call zfftifc(3,ngdgc,1,wfir)
  end if
! add to density and magnetisation
  if (spinpol) then
    ! spin-polarised
    if (ncmag) then
      ! non-collinear
      call rmk1(ngtc,wo,wfir,wfir(:,2),rhoir_n(:,j),magir_n(:,1,j),magir_n(:,2,j),magir_n(:,3,j))
    end if
  end if
end do
deallocate(wfmt)
call timesec(ts1)

timerho=timerho+ts1-ts0
return

contains

pure subroutine rmk1(n,wo,wf1,wf2,rho,mag1,mag2,mag3)
implicit none
! arguments
integer, intent(in) :: n
real(8), intent(in) :: wo
complex(8), intent(in) :: wf1(n),wf2(n)
real(8), intent(inout) :: rho(n),mag1(n),mag2(n),mag3(n)
! local variables
integer i
real(8) wo2,t1,t2
real(8) a1,b1,a2,b2
wo2=2.d0*wo
do i=1,n
  a1=dble(wf1(i)); b1=aimag(wf1(i))
  a2=dble(wf2(i)); b2=aimag(wf2(i))
  t1=a1**2+b1**2; t2=a2**2+b2**2
  mag1(i)=mag1(i)+wo2*(a1*a2+b1*b2)
  mag2(i)=mag2(i)+wo2*(a1*b2-b1*a2)
  mag3(i)=mag3(i)+wo*(t1-t2)
  rho(i)=rho(i)+wo*(t1+t2)
end do
end subroutine

end subroutine
!EOC

