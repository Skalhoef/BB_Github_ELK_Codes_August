
! Copyright (C) 2009 J. K. Dewhurst, S. Sharma and E. K. U. Gross
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.

!BOP
! !ROUTINE: rhomagsh_Sebbe
! !INTERFACE:
subroutine rhomagsh_Sebbe(magmt_nk)
! !USES:
use modmain
use modomp
! !DESCRIPTION:
!   Converts the muffin-tin density and magnetisation from spherical coordinates
!   to a spherical harmonic expansion. See {\tt rhomagk}.
!
! !REVISION HISTORY:
!   Created January 2009 (JKD)
!EOP
!BOC
implicit none
! local variables
integer idm,is,ias,nthd

real(8), intent(in) :: magmt_nk(npmtmax,natmtot,ndmag)

do idm=1,ndmag
  do ias=1,natmtot
    is=idxis(ias)
! convert the magnetisation to spherical harmonics
    call rfshtip(nrcmt(is),nrcmti(is),magmt_nk(:,ias,idm))
  end do
end do

end subroutine
!EOC

