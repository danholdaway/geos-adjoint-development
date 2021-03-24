      program compute_vcorr_length 
!
! Compute EOFs (or principal components) of vertical correlation matrices.
! Can just consider levels k=ktop,...,nlevels, where ktop<=1
!
      use netcdf
!
      implicit none
!
      logical, parameter :: l_hmean=.true. ! true if mean height used
      integer, parameter :: r4=4
      integer, parameter :: r8=8
!
      integer, parameter :: ktop=38
      integer :: krange
      integer :: unitin
      integer :: nr, k, k2, k3, k4, kf, kfx, nrx  
      integer :: imax, jmax, nlevels, k_fields, n_regions, ntimes
      real(r4) :: plev, dsig_sum, d2
      real(r4), allocatable :: pk(:), dsig(:)
      real(r4), allocatable :: xregions(:,:)
      real(r4), allocatable :: fmean(:,:,:)
      real(r4), allocatable :: vcorr(:,:,:,:)
      real(r4), allocatable :: vcross(:,:,:,:)
      real(r8), allocatable :: wcov(:,:,:,:)
      real(r8), allocatable :: evects(:,:,:,:)
      real(r8), allocatable :: evalues(:,:,:)
      real(r4), allocatable :: cprod(:,:,:,:)
      character(len=4), allocatable :: field_names(:)
      character(len=30) :: hk_header
      character(len=1) :: cdum
      data unitin/10/
!
      real(r8), ALLOCATABLE, DIMENSION(:,:) :: T_EVECTS, Q_EVECTS, U_EVECTS, V_EVECTS
      real(r8), ALLOCATABLE, DIMENSION(:)   :: T_EVALUS, Q_EVALUS, U_EVALUS, V_EVALUS

      real(r8), ALLOCATABLE, DIMENSION(:,:) :: COV_T, COV_Q, COV_U, COV_V

      !NETCDF
      integer :: rcode, ncID, dimID, varID
      integer :: w_dimid, x_dimid, y_dimid, z_dimid
      integer :: varid1, varid2, varid3, varid4, varid5, varid6, varid7, varid8, varid9, varid10, varid11, varid12, varid13
      character(len=64) :: input_file, output_file
      character(len=64) :: charbuf


      open (unitin, file='vcorr_diff.bin',form='unformatted')
      print *,'File for unit ',unitin,' opened'
      read (unitin) imax,jmax,nlevels,k_fields,n_regions,ntimes
      print *,'header=',imax,jmax,nlevels,k_fields,n_regions,ntimes
      allocate (field_names(k_fields))
      allocate (fmean(nlevels,k_fields,n_regions))
      allocate (vcorr(nlevels,nlevels,k_fields,n_regions))
      allocate (vcross(nlevels,nlevels,k_fields,n_regions))
      allocate (pk(nlevels))
      allocate (xregions(n_regions,5))
      read (unitin) ! skip file names
      read (unitin) pk  ! pressure of level k in mb
      read (unitin) ! skip regions
      do nr=1,n_regions
        do kf=1,k_fields
          read (unitin) kfx,nrx,field_names(kf),xregions(nr,:), &
                           fmean(:,kf,nr),vcorr(:,:,kf,nr),vcross(:,:,kf,nr)
        enddo
      enddo
      print *,'File read on unit=',unitin
!
!print*, field_names(1)
!print*, field_names(2)
!print*, field_names(3)
!print*, field_names(4)
!print*, xregions(1,:)
!print*, xregions(2,:)
!print*, xregions(3,:)
!print*, xregions(4,:)
!
      krange=nlevels-ktop+1
      allocate (wcov(krange,krange,k_fields,n_regions))
      allocate (evects(krange,krange,k_fields,n_regions))
      allocate (cprod(krange,krange,k_fields,n_regions))
      allocate (evalues(krange,k_fields,n_regions))
      allocate (dsig(krange))
      allocate (COV_T(krange,krange))
      allocate (COV_Q(krange,krange))
      allocate (COV_U(krange,krange))
      allocate (COV_V(krange,krange))
      allocate (T_EVECTS(krange,krange))
      allocate (Q_EVECTS(krange,krange))
      allocate (U_EVECTS(krange,krange))
      allocate (V_EVECTS(krange,krange))
      allocate (T_EVALUS(krange))
      allocate (Q_EVALUS(krange))
      allocate (U_EVALUS(krange))
      allocate (V_EVALUS(krange))
!
! Compute dsig
      dsig_sum=0.
      plev=1000.  ! surface pressure in mb
      do k=nlevels,ktop,-1
        k2=krange+k-nlevels
        dsig(k2)=2.*(plev-pk(k))
        plev=plev-dsig(k2) 
        dsig_sum=dsig_sum+dsig(k2)
      enddo  
      do k=1,krange
        dsig(k)=sqrt(dsig(k)/dsig_sum)
      enddo  
print*, 'dsig', sum(dsig)
!
! Compute weighted covariances matrix
      do nr=1,n_regions
        do kf=1,k_fields
          do k=1,krange
            k3=k+ktop-1
            do k2=1,krange
              k4=k2+ktop-1
              wcov(k2,k,kf,nr)=dsig(k2)*dsig(k)*0.5* &
                  (vcross(k4,k3,kf,nr)+vcross(k3,k4,kf,nr))
            enddo
          enddo
        enddo
      enddo
COV_T = wcov(:,:,1,3)
COV_Q = wcov(:,:,2,3)
COV_U = wcov(:,:,3,3)
COV_V = wcov(:,:,4,3)

print*, 'krange = ', krange
do k = 1,krange
print*, sqrt(wcov(k,k,4,3))
enddo
!
! Compute EOFS and eigenvalues
      do nr=1,n_regions
        do kf=1,k_fields
          call pert_eigen (krange,wcov(1,1,kf,nr),evects(1,1,kf,nr), &
                           evalues(1,kf,nr))
          call pert_normalize_evects (krange,evects(1,1,kf,nr))
          call pert_fix_evalues (krange,evalues(1,kf,nr))
          call Reorder (evects(1,1,kf,nr),evalues(1,kf,nr),krange)
          print ('(a,2i4,1p2e15.5)'),'nr, kf, max min evalues = ', &
             nr,kf,evalues(1,kf,nr),evalues(krange,kf,nr)
       enddo
     enddo
!
! Compute cross products
      do nr=1,n_regions
        do kf=1,k_fields
          do k=1,krange
            do k2=1,krange
              cprod(k2,k,kf,nr)=0.
              do k3=1,krange
                cprod(k2,k,kf,nr)=cprod(k2,k,kf,nr)+ &
                    evects(k3,k,kf,nr)*evects(k3,k2,kf,nr)
              enddo
            enddo
          enddo
       enddo
     enddo
!
! Remove weight from eigenvectors
      do nr=1,n_regions
        do kf=1,k_fields
          do k=1,krange
            do k2=1,krange
              evects(k2,k,kf,nr)=evects(k2,k,kf,nr)/dsig(k2)
            enddo
          enddo
       enddo
     enddo
print*, 'sumevec = ', (EVECTS(:,1,1,3))*dsig
!
! Print results
      do nr=1,1   ! n_regions
        do kf=3,3 ! k_fields
          print ('(2(a,i4))'),'nr=',nr,'  kf=',kf
          print ('(12x,1p10e11.1)'),evalues(1:10,kf,nr)
          do k3=1,krange
            k4=k3+ktop-1
            d2=dsig(k3)**2
            print ('(2i3,f8.2,f8.4,1p10e10.1)'), &
                    k3,k4,pk(k4),d2,evects(k3,1:10,kf,nr)
          enddo
          print ('(2(a,i4))'),'nr=',nr,'  kf=',kf
          print ('(12x,1p10e11.1)'),evalues(1:10,kf,nr)
          do k2=1,10
            print ('(i3,1p10e11.1)'),k2,cprod(k2,1:10,kf,nr)
          enddo
        enddo
      enddo
!
T_EVECTS = EVECTS(:,:,1,3)
Q_EVECTS = EVECTS(:,:,2,3)
U_EVECTS = EVECTS(:,:,3,3)
V_EVECTS = EVECTS(:,:,4,3)
T_EVALUS = EVALUES(:,1,3)
Q_EVALUS = EVALUES(:,2,3)
U_EVALUS = EVALUES(:,3,3)
V_EVALUS = EVALUES(:,4,3)

! Write Data to nc4 file
   output_file = 'EVECSEVALS.nc4'
   rcode = nf90_create(trim(adjustl(output_file)), NF90_CLOBBER, ncid)
   rcode = nf90_def_dim(ncid, "x", krange, x_dimid) 
   rcode = nf90_def_dim(ncid, "y", krange, y_dimid)
   rcode = nf90_def_var(ncid, "T_EVECTS", NF90_DOUBLE, (/ y_dimid, x_dimid /), varid1)
   rcode = nf90_def_var(ncid, "Q_EVECTS", NF90_DOUBLE, (/ y_dimid, x_dimid /), varid2)
   rcode = nf90_def_var(ncid, "U_EVECTS", NF90_DOUBLE, (/ y_dimid, x_dimid /), varid3)
   rcode = nf90_def_var(ncid, "V_EVECTS", NF90_DOUBLE, (/ y_dimid, x_dimid /), varid4)   
   rcode = nf90_def_var(ncid, "T_EVALUS", NF90_DOUBLE, (/ y_dimid /), varid5)
   rcode = nf90_def_var(ncid, "Q_EVALUS", NF90_DOUBLE, (/ y_dimid /), varid6)
   rcode = nf90_def_var(ncid, "U_EVALUS", NF90_DOUBLE, (/ y_dimid /), varid7)
   rcode = nf90_def_var(ncid, "V_EVALUS", NF90_DOUBLE, (/ y_dimid /), varid8)
   rcode = nf90_enddef(ncid)
   rcode = nf90_put_var(ncid, varid1, T_EVECTS)
   rcode = nf90_put_var(ncid, varid2, Q_EVECTS)
   rcode = nf90_put_var(ncid, varid3, U_EVECTS)
   rcode = nf90_put_var(ncid, varid4, V_EVECTS)
   rcode = nf90_put_var(ncid, varid5, T_EVALUS)
   rcode = nf90_put_var(ncid, varid6, Q_EVALUS)
   rcode = nf90_put_var(ncid, varid7, U_EVALUS)
   rcode = nf90_put_var(ncid, varid8, V_EVALUS)
   rcode = nf90_close(ncid)
!
      end program compute_vcorr_length 
!
!
!  X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X
!
      Subroutine Reorder (evectrs,evalues,nk)
!
!  Reorder eigenvalues and corresponding eigenvectors from largest 
!  eigenvalue to smallest. 
!
      integer, intent(in)  :: nk 
      real(8), intent(inout)  :: evectrs(nk,nk)
      real(8), intent(inout)  :: evalues(nk)
!
      integer :: k, kx,i
      real(8) :: copy, emax     
!     
      do k=1,nk-1
        emax=-1.e-30 
        do i=k,nk
          if (evalues(i).gt.emax) then
            emax=evalues(i)
            kx=i
          endif 
        enddo
        copy=evalues(kx)
        evalues(kx)=evalues(k)
        evalues(k)=copy
        do i=1,nk
          copy=evectrs(i,kx)
          evectrs(i,kx)=evectrs(i,k)
          evectrs(i,k)=copy
        enddo
      enddo
!
      end subroutine reorder
!
!
!  X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X
!
      subroutine pert_eigen (nk,matrix,evects,evalues)
!
! Call library routine to compute eigenvalues and vectors of a symmetric 
! matrix (stored as an array of 8-byte values).
!
      implicit none
!
      integer,  intent(in)   :: nk
      real(8), intent(in)   :: matrix(nk,nk)
      real(8), intent(out)  :: evects(nk,nk)
      real(8), intent(out)  :: evalues(nk)
!
      integer(4)     :: nwork1  
      integer(4)     :: info 
      integer(4)     :: nk4
      real(8), allocatable :: work1(:) ! must be size >= 3*nk-1
      external dsyev
!
      nk4=nk
      nwork1=3*nk-1
      allocate (work1(nwork1))
      evects=matrix
      call dsyev ( 'V','U',nk4,evects,nk4,evalues,work1,nwork1,info)
      if (info /= 0) then
        print *,'   '
        print *,' * * * * * * * * * * * * * * * * * * * * * * *  '
        print *,' info from dsyev =  ',info
        print *,'   '
      endif
!
      deallocate (work1)  
!
      end subroutine pert_eigen
!
!
!  X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X
!
      Subroutine pert_fix_evalues (nk,evalues)
!
!  Set small (relative to largest values) or negative values to 
!  a small value positive value
!
      integer,  intent(in)    :: nk
      real(8), intent(inout) :: evalues(nk)
!
      integer  :: k
      real(8) :: emax, emin
!
      emax=0.
      do k=1,nk
        if (abs(evalues(k)) > emax) emax=abs(evalues(k))
      enddo
!
      if (rkind3==8) then
        emin=emax*1.e-12
      else
        emin=emax*1.e-6
      endif
!
      do k=1,nk
        if (abs(evalues(k)) < emin) evalues(k)=emin
      enddo
!
      end subroutine pert_fix_evalues
!
!!
!  X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X
!
      Subroutine pert_normalize_evects (nk,evects)
!
!  Normalize each eigenvector such that the sum of its squared 
!  components is 1
!
      implicit none
!
      integer,  intent(in)    :: nk
      real(8), intent(inout) :: evects(nk,nk)
!
      integer  :: k, i
      real(8) :: sum
     
      do k=1,nk
        sum=0.
        do i=1,nk
          sum=sum+evects(i,k)*evects(i,k)
        enddo
        sum=1./sqrt(sum)
        evects(:,k)=evects(:,k)*sum
      enddo
!
      end subroutine pert_normalize_evects
!
