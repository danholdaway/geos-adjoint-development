SUBROUTINE adv2_tracer_2d_1l(q, adq, dp0, addp0, mfx, admfx, mfy, admfy, cx, &
& adcx, cy, adcy, npx, npy, npz, nq, hord, q_split, k, q3, adq3, dt, &
& uniform_ppm, id_divg)

!==============================================
! referencing used modules
!==============================================
use adv2_tracer_2d_1l_store, only : trc_1l_nsplt_3h,trc_1l_v2_tracer_2d_1l,trc_1lnspzq_cx_2h,trc_1lnspzq_cy_3h,trc_1lnspzq_fx_13h,&
&trc_1lnspzq_fx_16h,trc_1lnspzq_fy_14h,trc_1lnspzq_fy_17h,trc_1lnspzq_mfx_10h,trc_1lnspzq_mfy_11h,trc_1lnspzq_q_12h,&
&trc_1lnspzq_q_15h,trc_1lnspzq_ra_x_8h,trc_1lnspzq_ra_y_9h,trc_1lnspzq_v2_tracer_2d_1l,trc_1lnspzq_xfx_6h,trc_1lnspzq_yfx_7h
use adv2_tp_core_mod, only : adv2_fv_tp_2d_trc
use adfv_my_mpp, only : admp_reduce_max_dummy, admpp_update_domains_dummy4

!==============================================
! all entries are defined explicitly
!==============================================
implicit none

  INTEGER, INTENT(IN) :: npx, npy, npz
  INTEGER, INTENT(IN) :: k
! number of tracers to be advected
  INTEGER, INTENT(IN) :: nq
  INTEGER, INTENT(IN) :: hord
  INTEGER, INTENT(IN) :: q_split
  logical, intent(in) :: uniform_ppm
  INTEGER, INTENT(IN) :: id_divg
  REAL, INTENT(IN) :: dt
! 2D Tracers
  REAL, INTENT(INOUT) :: q(isd:ied, jsd:jed, nq)
  REAL, INTENT(INOUT) :: adq(isd:ied, jsd:jed, nq)
! Tracers
  REAL, INTENT(INOUT) :: q3(isd:ied, jsd:jed, npz, nq)
  REAL, INTENT(INOUT) :: adq3(isd:ied, jsd:jed, npz, nq)
! DELP before dyn_core
  REAL, INTENT(INOUT) :: dp0(is:ie, js:je)
  REAL, INTENT(INOUT) :: addp0(is:ie, js:je)
! Mass Flux X-Dir
  REAL, INTENT(IN) :: mfx(is:ie+1, js:je)
  REAL :: admfx(is:ie+1, js:je)
! Mass Flux Y-Dir
  REAL, INTENT(IN) :: mfy(is:ie, js:je+1)
  REAL :: admfy(is:ie, js:je+1)
! Courant Number X-Dir
  REAL, INTENT(IN) :: cx(is:ie+1, jsd:jed)
  REAL :: adcx(is:ie+1, jsd:jed)
! Courant Number Y-Dir
  REAL, INTENT(IN) :: cy(isd:ied, js:je+1)
  REAL :: adcy(isd:ied, js:je+1)
! Local Arrays
  REAL :: mfx2(is:ie+1, js:je)
  REAL :: admfx2(is:ie+1, js:je)
  REAL :: mfy2(is:ie, js:je+1)
  REAL :: admfy2(is:ie, js:je+1)
  REAL :: cx2(is:ie+1, jsd:jed)
  REAL :: adcx2(is:ie+1, jsd:jed)
  REAL :: cy2(isd:ied, js:je+1)
  REAL :: adcy2(isd:ied, js:je+1)
  REAL :: dp1(is:ie, js:je)
  REAL :: addp1(is:ie, js:je)
  REAL :: dp2(is:ie, js:je)
  REAL :: addp2(is:ie, js:je)
  REAL :: fx(is:ie+1, js:je)
  REAL :: adfx(is:ie+1, js:je)
  REAL :: fy(is:ie, js:je+1)
  REAL :: adfy(is:ie, js:je+1)
  REAL :: ra_x(is:ie, jsd:jed)
  REAL :: adra_x(is:ie, jsd:jed)
  REAL :: ra_y(isd:ied, js:je)
  REAL :: adra_y(isd:ied, js:je)
  REAL :: xfx(is:ie+1, jsd:jed)
  REAL :: adxfx(is:ie+1, jsd:jed)
  REAL :: yfx(isd:ied, js:je+1)
  REAL :: adyfx(isd:ied, js:je+1)
  REAL :: cmax
  REAL :: frac, rdt
  INTEGER :: nsplt
  INTEGER :: i, j, it, iq
  INTRINSIC ABS
  INTRINSIC MAX
  INTRINSIC INT
  INTRINSIC REAL
  INTEGER :: branch
  REAL :: tempb4
  REAL :: tempb3
  REAL :: tempb2
  REAL :: tempb1
  REAL :: tempb0
  REAL :: x1
  REAL :: tempb
  REAL :: abs1
  REAL :: abs0
  REAL :: y1
  DO j=jsd,jed
    DO i=is,ie+1
      IF (cx(i, j) .GT. 0.) THEN
        xfx(i, j) = cx(i, j)*dxa(i-1, j)*dy(i, j)*sina_u(i, j)
        CALL PUSHCONTROL1B(1)
      ELSE
        xfx(i, j) = cx(i, j)*dxa(i, j)*dy(i, j)*sina_u(i, j)
        CALL PUSHCONTROL1B(0)
      END IF
    END DO
  END DO
  DO j=js,je+1
    DO i=isd,ied
      IF (cy(i, j) .GT. 0.) THEN
        yfx(i, j) = cy(i, j)*dya(i, j-1)*dx(i, j)*sina_v(i, j)
        CALL PUSHCONTROL1B(1)
      ELSE
        yfx(i, j) = cy(i, j)*dya(i, j)*dx(i, j)*sina_v(i, j)
        CALL PUSHCONTROL1B(0)
      END IF
    END DO
  END DO
  IF (q_split .EQ. 0) THEN
! Determine nsplt for tracer advection
    cmax = 0.
    DO j=js,je
      DO i=is,ie
        IF (cx(i, j) .GE. 0.) THEN
          abs0 = cx(i, j)
        ELSE
          abs0 = -cx(i, j)
        END IF
        x1 = abs0 + (1.-sina_u(i, j))
        IF (cy(i, j) .GE. 0.) THEN
          abs1 = cy(i, j)
        ELSE
          abs1 = -cy(i, j)
        END IF
        y1 = abs1 + (1.-sina_v(i, j))
        IF (x1 .LT. y1) THEN
          IF (y1 .LT. cmax) THEN
            cmax = cmax
          ELSE
            cmax = y1
          END IF
        ELSE IF (x1 .LT. cmax) THEN
          cmax = cmax
        ELSE
          cmax = x1
        END IF
      END DO
    END DO
    CALL MP_REDUCE_MAX(cmax)
    nsplt = INT(1.01 + cmax)
  ELSE
    nsplt = q_split
  END IF
  frac = 1./REAL(nsplt)
  DO j=jsd,jed
    DO i=is,ie+1
      cx2(i, j) = cx(i, j)*frac
      xfx(i, j) = xfx(i, j)*frac
    END DO
  END DO
  DO j=js,je
    DO i=is,ie+1
      mfx2(i, j) = mfx(i, j)*frac
    END DO
  END DO
  DO j=js,je+1
    DO i=isd,ied
      cy2(i, j) = cy(i, j)*frac
      yfx(i, j) = yfx(i, j)*frac
    END DO
  END DO
  DO j=js,je+1
    DO i=is,ie
      mfy2(i, j) = mfy(i, j)*frac
    END DO
  END DO
  DO j=jsd,jed
    DO i=is,ie
      ra_x(i, j) = area(i, j) + xfx(i, j) - xfx(i+1, j)
    END DO
  END DO
  DO j=js,je
    DO i=isd,ied
      ra_y(i, j) = area(i, j) + yfx(i, j) - yfx(i, j+1)
    END DO
  END DO
  DO j=js,je
    DO i=is,ie
      dp1(i, j) = dp0(i, j)
    END DO
  END DO
  DO it=1,nsplt
    DO j=js,je
      DO i=is,ie
        CALL PUSHREAL4(dp2(i, j))
        dp2(i, j) = dp1(i, j) + (mfx2(i, j)-mfx2(i+1, j)+mfy2(i, j)-mfy2&
&         (i, j+1))*rarea(i, j)
      END DO
    END DO
!call timing_on('COMM_TOTAL')
!call timing_on('COMM_TRAC')
    CALL PUSHREAL4ARRAY(q, (ied-isd+1)*(jed-jsd+1)*nq)
    call mpp_update_domains( q, domain, complete=.true. )
!call timing_off('COMM_TRAC')
!call timing_off('COMM_TOTAL')
    DO iq=1,nq
      CALL PUSHREAL4ARRAY(mfy2, (ie-is+1)*(je-js+2))
      CALL PUSHREAL4ARRAY(mfx2, (ie-is+2)*(je-js+1))
      CALL PUSHREAL4ARRAY(ra_y, (ied-isd+1)*(je-js+1))
      CALL PUSHREAL4ARRAY(ra_x, (ie-is+1)*(jed-jsd+1))
      CALL PUSHREAL4ARRAY(yfx, (ied-isd+1)*(je-js+2))
      CALL PUSHREAL4ARRAY(xfx, (ie-is+2)*(jed-jsd+1))
      CALL PUSHREAL4ARRAY(fy, (ie-is+1)*(je-js+2))
      CALL PUSHREAL4ARRAY(fx, (ie-is+2)*(je-js+1))
      CALL PUSHREAL4ARRAY(cy2, (ied-isd+1)*(je-js+2))
      CALL PUSHREAL4ARRAY(cx2, (ie-is+2)*(jed-jsd+1))
      CALL PUSHREAL4ARRAY(q(:, :, iq), (ied-isd+1)*(jed-jsd+1))
      CALL fv_tp_2d_trc(q(isd:, jsd:, iq), cx2, cy2, npx, npy, hord, fx, fy&
&             , xfx, yfx, area, ra_x, ra_y, mfx2, mfy2)
      IF (it .EQ. nsplt) THEN
        CALL PUSHCONTROL1B(1)
      ELSE
        DO j=js,je
          DO i=is,ie
            CALL PUSHREAL4(q(i, j, iq))
            q(i, j, iq) = (q(i, j, iq)*dp1(i, j)+(fx(i, j)-fx(i+1, j)+fy&
&             (i, j)-fy(i, j+1))*rarea(i, j))/dp2(i, j)
          END DO
        END DO
        CALL PUSHCONTROL1B(0)
      END IF
    END DO
    IF (it .NE. nsplt) THEN
      DO j=js,je
        DO i=is,ie
          CALL PUSHREAL4(dp1(i, j))
          dp1(i, j) = dp2(i, j)
        END DO
      END DO
      CALL PUSHCONTROL1B(1)
    ELSE
      CALL PUSHCONTROL1B(0)
    END IF
  END DO
! nsplt
  IF (id_divg .GT. 0) THEN
    rdt = 1./(frac*dt)
    adxfx = 0.0
    adyfx = 0.0
    DO j=je,js,-1
      DO i=ie,is,-1
        tempb4 = rarea(i, j)*rdt*addp0(i, j)
        adxfx(i+1, j) = adxfx(i+1, j) + tempb4
        adxfx(i, j) = adxfx(i, j) - tempb4
        adyfx(i, j+1) = adyfx(i, j+1) + tempb4
        adyfx(i, j) = adyfx(i, j) - tempb4
        addp0(i, j) = 0.0
      END DO
    END DO
  ELSE
    adxfx = 0.0
    adyfx = 0.0
  END IF
  adcx2 = 0.0
  admfx2 = 0.0
  addp1 = 0.0
  addp2 = 0.0
  adcy2 = 0.0
  admfy2 = 0.0
  adra_x = 0.0
  adra_y = 0.0
  DO it=nsplt,1,-1
    CALL POPCONTROL1B(branch)
    IF (branch .NE. 0) THEN
      DO j=je,js,-1
        DO i=ie,is,-1
          CALL POPREAL4(dp1(i, j))
          addp2(i, j) = addp2(i, j) + addp1(i, j)
          addp1(i, j) = 0.0
        END DO
      END DO
    END IF
    DO iq=nq,1,-1
      iq_taf = (k-1)*nq*nsplt_max+(it-1)*nq+iq
      CALL POPCONTROL1B(branch)
      IF (branch .EQ. 0) THEN
        adfx = 0.0
        adfy = 0.0
        DO j=je,js,-1
          DO i=ie,is,-1
            CALL POPREAL4(q(i, j, iq))
            tempb2 = adq(i, j, iq)/dp2(i, j)
            tempb3 = rarea(i, j)*tempb2
            addp1(i, j) = addp1(i, j) + q(i, j, iq)*tempb2
            adfx(i, j) = adfx(i, j) + tempb3
            adfx(i+1, j) = adfx(i+1, j) - tempb3
            adfy(i, j) = adfy(i, j) + tempb3
            adfy(i, j+1) = adfy(i, j+1) - tempb3
            addp2(i, j) = addp2(i, j) - (q(i, j, iq)*dp1(i, j)+rarea(i, j)&
&             *(fx(i, j)-fx(i+1, j)+fy(i, j)-fy(i, j+1)))*tempb2/dp2(i, &
&             j)
            adq(i, j, iq) = dp1(i, j)*tempb2
          END DO
        END DO
      ELSE
        adfx = 0.0
        adfy = 0.0
        DO j=je,js,-1
          DO i=ie,is,-1
            tempb0 = adq3(i, j, k, iq)/dp2(i, j)
            tempb1 = rarea(i, j)*tempb0
            adq(i, j, iq) = adq(i, j, iq) + dp1(i, j)*tempb0
            addp1(i, j) = addp1(i, j) + q(i, j, iq)*tempb0
            adfx(i, j) = adfx(i, j) + tempb1
            adfx(i+1, j) = adfx(i+1, j) - tempb1
            adfy(i, j) = adfy(i, j) + tempb1
            adfy(i, j+1) = adfy(i, j+1) - tempb1
            addp2(i, j) = addp2(i, j) - (q(i, j, iq)*dp1(i, j)+rarea(i, j)&
&             *(fx(i, j)-fx(i+1, j)+fy(i, j)-fy(i, j+1)))*tempb0/dp2(i, &
&             j)
            adq3(i, j, k, iq) = 0.0
          END DO
        END DO
      END IF
      CALL POPREAL4ARRAY(q(:, :, iq), (ied-isd+1)*(jed-jsd+1))
      CALL POPREAL4ARRAY(cx2, (ie-is+2)*(jed-jsd+1))
      CALL POPREAL4ARRAY(cy2, (ied-isd+1)*(je-js+2))
      CALL POPREAL4ARRAY(fx, (ie-is+2)*(je-js+1))
      CALL POPREAL4ARRAY(fy, (ie-is+1)*(je-js+2))
      CALL POPREAL4ARRAY(xfx, (ie-is+2)*(jed-jsd+1))
      CALL POPREAL4ARRAY(yfx, (ied-isd+1)*(je-js+2))
      CALL POPREAL4ARRAY(ra_x, (ie-is+1)*(jed-jsd+1))
      CALL POPREAL4ARRAY(ra_y, (ied-isd+1)*(je-js+1))
      CALL POPREAL4ARRAY(mfx2, (ie-is+2)*(je-js+1))
      CALL POPREAL4ARRAY(mfy2, (ie-is+1)*(je-js+2))
      CALL adv2_fv_tp_2d_trc(adq(isd:, jsd:, iq), cx2, adcx2, &
&               cy2, adcy2, npx, npy, hord, fx, adfx, fy, adfy, xfx, adxfx, &
&               yfx, adyfx, area, ra_x, adra_x, ra_y, adra_y, mfx2, &
&               admfx2, mfy2, admfy2, .true. ,iq_taf )
    END DO
    CALL POPREAL4ARRAY(q, (ied-isd+1)*(jed-jsd+1)*nq)
    call admpp_update_domains_dummy4( adq,is,ie,js,je,isd,ied,jsd,jed,1,nq )
    DO j=je,js,-1
      DO i=ie,is,-1
        CALL POPREAL4(dp2(i, j))
        tempb = rarea(i, j)*addp2(i, j)
        addp1(i, j) = addp1(i, j) + addp2(i, j)
        admfx2(i, j) = admfx2(i, j) + tempb
        admfx2(i+1, j) = admfx2(i+1, j) - tempb
        admfy2(i, j) = admfy2(i, j) + tempb
        admfy2(i, j+1) = admfy2(i, j+1) - tempb
        addp2(i, j) = 0.0
      END DO
    END DO
  END DO
  DO j=je,js,-1
    DO i=ie,is,-1
      addp0(i, j) = addp0(i, j) + addp1(i, j)
      addp1(i, j) = 0.0
    END DO
  END DO
  DO j=je,js,-1
    DO i=ied,isd,-1
      adyfx(i, j) = adyfx(i, j) + adra_y(i, j)
      adyfx(i, j+1) = adyfx(i, j+1) - adra_y(i, j)
      adra_y(i, j) = 0.0
    END DO
  END DO
  DO j=jed,jsd,-1
    DO i=ie,is,-1
      adxfx(i, j) = adxfx(i, j) + adra_x(i, j)
      adxfx(i+1, j) = adxfx(i+1, j) - adra_x(i, j)
      adra_x(i, j) = 0.0
    END DO
  END DO
  admfy = 0.0
  DO j=je+1,js,-1
    DO i=ie,is,-1
      admfy(i, j) = admfy(i, j) + frac*admfy2(i, j)
      admfy2(i, j) = 0.0
    END DO
  END DO
  adcy = 0.0
  DO j=je+1,js,-1
    DO i=ied,isd,-1
      adyfx(i, j) = frac*adyfx(i, j)
      adcy(i, j) = adcy(i, j) + frac*adcy2(i, j)
      adcy2(i, j) = 0.0
    END DO
  END DO
  admfx = 0.0
  DO j=je,js,-1
    DO i=ie+1,is,-1
      admfx(i, j) = admfx(i, j) + frac*admfx2(i, j)
      admfx2(i, j) = 0.0
    END DO
  END DO
  adcx = 0.0
  DO j=jed,jsd,-1
    DO i=ie+1,is,-1
      adxfx(i, j) = frac*adxfx(i, j)
      adcx(i, j) = adcx(i, j) + frac*adcx2(i, j)
      adcx2(i, j) = 0.0
    END DO
  END DO
  DO j=je+1,js,-1
    DO i=ied,isd,-1
      CALL POPCONTROL1B(branch)
      IF (branch .EQ. 0) THEN
        adcy(i, j) = adcy(i, j) + dya(i, j)*dx(i, j)*sina_v(i, j)*adyfx(i, &
&         j)
        adyfx(i, j) = 0.0
      ELSE
        adcy(i, j) = adcy(i, j) + dx(i, j)*sina_v(i, j)*dya(i, j-1)*adyfx(i&
&         , j)
        adyfx(i, j) = 0.0
      END IF
    END DO
  END DO
  DO j=jed,jsd,-1
    DO i=ie+1,is,-1
      CALL POPCONTROL1B(branch)
      IF (branch .EQ. 0) THEN
        adcx(i, j) = adcx(i, j) + dxa(i, j)*dy(i, j)*sina_u(i, j)*adxfx(i, &
&         j)
        adxfx(i, j) = 0.0
      ELSE
        adcx(i, j) = adcx(i, j) + dy(i, j)*sina_u(i, j)*dxa(i-1, j)*adxfx(i&
&         , j)
        adxfx(i, j) = 0.0
      END IF
    END DO
  END DO
END SUBROUTINE adv2_tracer_2d_1l

