!xtp

 elseif ( iord==333 ) then

!$AD II-LOOP
      do j=js,je+1
        do i=is,ie+1
           if( c(i,j)>0. ) then
               flux(i,j) = (2.0*u(i,j)+5.0*u(i-1,j)-u(i-2,j) )/6.0 &
                            - 0.5*c(i,j)*rdx(i-1,j)*(u(i,j)-u(i-1,j)) &
                            + (c(i,j)*rdx(i-1,j)*c(i,j)*rdx(i-1,j)/6.0)*(u(i,j)-2.0*u(i-1,j)+u(i-2,j) )
           else
               flux(i,j) = (2.0*u(i-1,j)+5.0*u(i,j)-u(i+1,j) )/6.0 &
                            - 0.5*c(i,j)*rdx(i,j)*(u(i,j)-u(i-1,j)) &
                            + (c(i,j)*rdx(i,j)*c(i,j)*rdx(i,j)/6.0)*(u(i+1,j)-2.0*u(i,j)+u(i-1,j) )
           endif
        enddo
     enddo


!ytp

 elseif ( jord==333 ) then

!$AD II-LOOP
      do j=js,je+1
         do i=is,ie+1
            if( c(i,j)>0. ) then
               flux(i,j) = (2.0*v(i,j)+5.0*v(i,j-1)-v(i,j-2) )/6.0 &
                            - 0.5*c(i,j)*rdy(i,j-1)*(v(i,j)-v(i,j-1)) &
                            + (c(i,j)*rdy(i,j-1)*c(i,j)*rdy(i,j-1)/6.0)*(v(i,j)-2.0*v(i,j-1)+v(i,j-2) ) 
            else
               flux(i,j) = (2.0*v(i,j-1)+5.0*v(i,j)-v(i,j+1) )/6.0 &
                            - 0.5*c(i,j)*rdy(i,j)*(v(i,j)-v(i,j-1)) &
                            + (c(i,j)*rdy(i,j)*c(i,j)*rdy(i,j)/6.0)*(v(i,j+1)-2.0*v(i,j)+v(i,j-1) ) 
            endif
         enddo
      enddo
