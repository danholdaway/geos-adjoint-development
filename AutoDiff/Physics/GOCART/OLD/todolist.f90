!PMAXMIN

REAL*4 :: tmp4(im,jm,km), qmin4, qmax4


      ijl = im*jm
      tmp4 = real(dz/vsettle,4)
      call pmaxmin ( 'Chem_Settling: dt', tmp4, qmin4, qmax4, ijl, km, 0. )
      qmin = dble(qmin4)
      qmax = dble(qmax4)


!CHEM_SETTLING

!Adjoint sets some variables to zero, control of this is done above, tmpu_ad should not be set to zero!


!DryDeposition also sets some variables to zero that should not be.
