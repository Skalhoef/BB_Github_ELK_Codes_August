
! 2025, The Half Blood Mathematician 
! Adapted from vecplot.f90 to vecplot_Sebbe.f90

!BOP
! !ROUTINE: vecplot
! !INTERFACE:
subroutine vecplot_nk_Sebbe(ist,ik,magmt_nk,magir_nk)
! !DESCRIPTION:
!   Outputs a 2D or 3D vector field for plotting. The vector field can be the
!   magnetisation vector field, ${\bf m}$; The magnetisation is obtained from the
!   spin density matrix, $\rho_{\alpha\beta}$, by solving
!   $$ \rho_{\alpha\beta}({\bf r})=\frac{1}{2}\left(n({\bf r})
!    \delta_{\alpha\beta}+\sigma\cdot {\bf m}({\bf r})\right), $$
!   where $n\equiv\tr\rho_{\alpha\beta}$ is the total density. In the case of 2D
!   plots, the magnetisation vectors are still 3D, but are in the coordinate
!   system of the plane.
!
! !REVISION HISTORY:
!   Created August 2004 (JKD)
!   Included electric field plots, August 2006 (JKD)
!EOP
!BOC
use modmain
use modsebbe 

implicit none
integer, intent(in) :: ist,ik
real(8), intent(in) :: magmt_nk(npmtmax,natmtot,ndmag),magir_nk(ngtot,ndmag)
character(len=50) :: filename

! allocatable arrays
real(8), allocatable :: rvfmt(:,:,:),rvfir(:,:)
allocate(rvfmt(npmtmax,natmtot,3),rvfir(ngtot,3))

! SK Begin: 
! Let's Test something.
! call readMag_Sebbe
! if (ncmag) then
!   rvfmt(:,:,:)=magmt(:,:,:)
!   rvfir(:,:)=magir(:,:)
! end if
! SK End


! magnetisation
if (ncmag) then
! non-collinear
  rvfmt(:,:,:)=magmt_nk(:,:,:)
  rvfir(:,:)=magir_nk(:,:)
end if



!$OMP CRITICAL(file_write)
write(filename,'("MAG2D_ist_",I0,"_ik_",I0,".OUT")') ist, ik
open(50,file=filename,form='FORMATTED')
call plot2d(.true.,50,3,rvfmt,rvfir)
close(50)
!$OMP END CRITICAL(file_write)

! deallocate(rvfmt,rvfir)
end subroutine
!EOC

