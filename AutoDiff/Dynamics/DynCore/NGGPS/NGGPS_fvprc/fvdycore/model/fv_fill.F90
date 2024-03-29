module fv_fill_mod
   use fv_arrays_mod,  only: REAL4, REAL8, FVPRC
   implicit none
   public fillz

!---- version number -----
   character(len=128) :: version = '$Id: fv_fill.F90,v 1.1.2.1.36.1.88.1.20.1.2.1.4.1 2017/02/16 03:47:47 aoloso Exp $'
   character(len=128) :: tagname = '$Name: Heracles-UNSTABLE_ncepdyn_Feb222017 $'

contains

 subroutine fillz(im, km, nq, q, dp)
   integer,  intent(in):: im                ! No. of longitudes
   integer,  intent(in):: km                ! No. of levels
   integer,  intent(in):: nq                ! Total number of tracers
   real(FVPRC) , intent(in)::  dp(im,km)       ! pressure thickness
   real(FVPRC) , intent(inout) :: q(im,km,nq)   ! tracer mixing ratio
! !LOCAL VARIABLES:
   logical:: zfix(im)
   real(FVPRC) ::  dm(km)
   integer i, k, ic
   real(FVPRC)  qup, qly, dup, dq, sum0, sum1, fac

   do ic=1,nq
! Top layer
      do i=1,im
         if( q(i,1,ic) < 0. ) then
             q(i,2,ic) = q(i,2,ic) + q(i,1,ic)*dp(i,1)/dp(i,2)
             q(i,1,ic) = 0.
          endif
      enddo

! Interior
      zfix(:) = .false.
      do k=2,km-1
         do i=1,im
         if( q(i,k,ic) < 0. ) then
             zfix(i) = .true.
             if ( q(i,k-1,ic) > 0. ) then
! Borrow from above
                dq = min ( q(i,k-1,ic)*dp(i,k-1), -q(i,k,ic)*dp(i,k) ) 
                q(i,k-1,ic) = q(i,k-1,ic) - dq/dp(i,k-1)
                q(i,k  ,ic) = q(i,k  ,ic) + dq/dp(i,k  )
             endif
             if ( q(i,k,ic)<0.0 .and. q(i,k+1,ic)>0. ) then
! Borrow from below:
                dq = min ( q(i,k+1,ic)*dp(i,k+1), -q(i,k,ic)*dp(i,k) ) 
                q(i,k+1,ic) = q(i,k+1,ic) - dq/dp(i,k+1)
                q(i,k  ,ic) = q(i,k  ,ic) + dq/dp(i,k  )
             endif
          endif
         enddo
      enddo
 
! Bottom layer
      k = km
      do i=1,im
         if( q(i,k,ic)<0. .and. q(i,k-1,ic)>0.) then
             zfix(i) = .true.
! Borrow from above
             qup =  q(i,k-1,ic)*dp(i,k-1)
             qly = -q(i,k  ,ic)*dp(i,k  )
             dup =  min(qly, qup)
             q(i,k-1,ic) = q(i,k-1,ic) - dup/dp(i,k-1) 
             q(i,k,  ic) = q(i,k,  ic) + dup/dp(i,k  )
          endif
      enddo

! Perform final check and non-local fix if needed
      do i=1,im
         if ( zfix(i) ) then

           sum0 = 0.
           do k=2,km
              dm(k) = q(i,k,ic)*dp(i,k)
              sum0 = sum0 + dm(k)
           enddo

           if ( sum0 > 0. ) then
             sum1 = 0.
             do k=2,km
                sum1 = sum1 + max(0., dm(k))
             enddo
             fac = sum0 / sum1
             do k=2,km
                q(i,k,ic) = max(0., fac*dm(k)/dp(i,k))
             enddo
           endif

         endif
      enddo

   enddo
 end subroutine fillz


end module fv_fill_mod
