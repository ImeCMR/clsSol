      program clsSol
      parameter (nmol=50)
      integer   c
      real xg,yg,zg,xt,yt,zt,xmax,ymax,zmax,xmin,ymin,zmin,xmean,ymean
      real zmean
      dimension ax(3000),ay(3000),az(3000) !added 3000 instead of nmol
      dimension bx(nmol),by(nmol),bz(nmol)
      dimension cx(nmol),cy(nmol),cz(nmol)
      dimension ex(nmol),ey(nmol),ez(nmol)
      dimension xyz(3,50)
      dimension sx(25000,nmol),sy(25000,nmol),sz(25000,nmol)
      character ofile*60,ifile*60
      character an(3000)*4,bn(nmol)*4,cn(nmol)*4,dn(nmol)*4 !added 3000 
      character sul1*4,sul2*4,sul3*4,sul4*4 
c 
      open(unit=7,file='clsSol.dat',form='formatted',status='old')
      read(7,*) xbox
      read(7,*) ybox
      read(7,*) zbox
      read(7,*) ns1  
      read(7,*) ms1  
      read(7,*) ns2   
      read(7,*) ms2  
      read(7,*) ns3
      read(7,*) ms3 
      read(7,*) ns4
      read(7,*) ms4 
      read(7,*) rclose
      read(7,'(a)') sul1
      read(7,'(a)') sul2
      read(7,'(a)') sul3
      read(7,'(a)') sul4
      read(7,'(a)') ofile
      read(7,'(a)') ifile ! >>>>>>>>>>>>>>>>>>>>>>>>
c
      open(unit=23,file=ifile,form='formatted',status='old') !>>>>>>>
      open(unit=8,file=ofile,form='formatted',status='unknown')
      write(8,'(a)') 'TITLE   Creation of a box of molecules for MD' 
      write(8,'(a)') 'REMARK     This is MD simulaion attempt '
      write(8,73)'CRYST1',xbox,ybox,zbox,90.,90.,90.,'P ',1       
 73   format(a,3f9.3,3f7.2,1x,a2,i1)
      write(8,'(a)') 'MODEL       1'       
c
c ----- rclose: Minimum distance between two atoms of two adjacent molecules 
c
      rclose=rclose**2
c
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      !comx=0.0
      !comy=0.0
      !comz=0.0
      !do i=1,ms1   
       ! read(7,'(13x,a,13x,3f8.3)') an(i),ax(i),ay(i),az(i)
       ! comx=comx+ax(i)
       ! comy=comy+ay(i)
       ! comz=comz+az(i)
      !end do
      !ax0=comx/ms1   
      !ay0=comy/ms1  
      !az0=comz/ms1   
      !do i=1,ms1    
       ! ax(i)=ax(i)-ax0
       ! ay(i)=ay(i)-ay0
       ! az(i)=az(i)-az0
      !end do
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!
!==========================================================
! Reading the coordinate file for the cluster 
!
!      
       read(23,*)
       read(23,*)
       read(23,*)
       read(23,*)
       do i=1,(ms1*ns1)
        read(23,'(13x,a,13x,3f8.3)') an(i),ax(i),ay(i),az(i)   
        !read(23,'(25x,a,4x,3f8.3)') an(i),ax(i),ay(i),az(i)
       end do
!        write(6,*) ax
!==============================================================
! Taking the gravity point of the cluster
!
       xt=0
       yt=0
       zt=0
       !
       do i = 1,3000
         if (ax(i) /= 0) then
           xt = xt + ax(i)
         end if
         if (ay(i) /= 0) then
           yt = yt + ay(i)
         end if
         if (az(i) /= 0) then
           zt = zt + az(i)
         end if
       end do
       xg = xt/(ns1*ms1)
       yg = yt/(ns1*ms1)
       zg = zt/(ns1*ms1)

!==============================================================
! Taking the mix values and max values for x,y,z and getting the mean
! 
       xmax = maxval(ax)
       ymax = maxval(ay)
       zmax = maxval(az)
       xmin = minval(ax)
       ymin = minval(ay)
       zmin = minval(az)
       !
       !
       xmean = (xmax+xmin)/2
       ymean = (ymax+ymin)/2
       zmean = (zmax+zmin)/2
!===============================================================       

c
      if( ns2 .ne. 0 )then
      comx=0.0
      comy=0.0
      comz=0.0
      do i=1,ms2   
        read(7,'(13x,a,13x,3f8.3)') bn(i),bx(i),by(i),bz(i)
        comx=comx+bx(i)
        comy=comy+by(i)
        comz=comz+bz(i)
      end do
      bx0=comx/ms2   
      by0=comy/ms2   
      bz0=comz/ms2   
      do i=1,ms2   
        bx(i)=bx(i)-bx0
        by(i)=by(i)-by0
        bz(i)=bz(i)-bz0
      end do
      endif
c
      if( ns3 .ne. 0 )then
      comx=0.0
      comy=0.0
      comz=0.0
      do i=1,ms3   
        read(7,'(13x,a,13x,3f8.3)') cn(i),cx(i),cy(i),cz(i)
        comx=comx+cx(i)
        comy=comy+cy(i)
        comz=comz+cz(i)
      end do
      cx0=comx/ms3   
      cy0=comy/ms3   
      cz0=comz/ms3   
      do i=1,ms3   
        cx(i)=cx(i)-cx0
        cy(i)=cy(i)-cy0
        cz(i)=cz(i)-cz0
      end do
      endif
c
      if( ns4 .ne. 0 )then
      comx=0.0
      comy=0.0
      comz=0.0
      do i=1,ms4   
        read(7,'(13x,a,13x,3f8.3)') dn(i),ex(i),ey(i),ez(i)
        comx=comx+ex(i)
        comy=comy+ey(i)
        comz=comz+ez(i)
      end do
      ex0=comx/ms4   
      ey0=comy/ms4   
      ez0=comz/ms4   
      do i=1,ms4   
        ex(i)=ex(i)-ex0
        ey(i)=ey(i)-ey0
        ez(i)=ez(i)-ez0
      end do
      endif
c
c ---------------- Check --
c
      ntmol   = ns1 + ns2 + ns3 + ns4    
      ntatom  = ns1*ms1 + ns2*ms2 + ns3*ms3 + ns4*ms4   
      write(6,*) ' '
      write(6,*) ' '
      write(6,*) ' Total number of molecules           = ',ntmol
      write(6,*) ' Total number of atoms in the system = ',ntatom
      write(6,*) ' '
      write(6,*) ' '
c
      print*,' ------- NTMOL = ',ntmol 
c
      io=0
      ir=0
!==========================================================
!      do irtry=1,5000*ns1   
!        if(ir.eq.ns1) goto 50
c
c ---- If more than one "solutes", get different orientations 
c
!        if( ns1 .gt. 1 )then
!          ax0=(rand()-0.5)*xbox
!          ay0=(rand()-0.5)*ybox
!          az0=(rand()-0.5)*zbox          
!          if( ms1 .gt. 1 )call rotate( ms1,ax,ay,az )
c
!        else
c
c -- If only one "solute" then do not rotate it
c
!          ax0 = 0.0
!          ay0 = 0.0
!          az0 = 0.0
!        endif
!        do i=1,ms1   
!          xyz(1,i) = ax(i)+ax0
!          xyz(2,i) = ay(i)+ay0
!          xyz(3,i) = az(i)+az0          
!        end do
!        do jat = 1,ms1   
!        do iout=1,ir
!        do iat = 1,ms1   
!          dx=xyz(1,jat)-sx(iout,iat)
!          dy=xyz(2,jat)-sy(iout,iat)
!          dz=xyz(3,jat)-sz(iout,iat)
!          dx=dx-xbox*nint(dx/xbox)
!          dy=dy-ybox*nint(dy/ybox)
!          dz=dz-zbox*nint(dz/zbox)
!          d2=dx**2+dy**2+dz**2
!          if(d2.lt.rclose) goto 60
!        end do
!        end do
!        end do
c
      !=======================================================
      do j=1,ns1
        ir=ir+1
        do i=1,ms1     
          io=io+1
          c=(ms1*j)-(ms1-i)
          write(8,100) io,an(c),sul1,ir,ax(c),ay(c),az(c),1.,0.
          ! write(6,*) io,an(i*j),sul1,ir,ax(i*j),ay(i*j),az(i*j),1.,0.
        end do
      end do
      !write(6,*) ax
      !=================================================............
!         do iat = 1,ms1   
!          sx(ir,iat)=xyz(1,iat)
!          sy(ir,iat)=xyz(2,iat)
!          sz(ir,iat)=xyz(3,iat)
!         end do
c
!         if( mod(ir,100) .eq. 0)print*,' Counts = ',ir,' out of ',ntmol 
c
!   60   continue
!      end do
!      stop ' not all the molecules were generated'
!   50 continue
!      write(6,'(a,a,a,i6)') ' Number of ',sul1,' molecules = ',nsolu
!========================================================================== 
c -------- 2nd one 
      if( ns2 .ne. 0 )then
      do irtry=1,2000*ns2    
        if(ir.eq.ns1+ns2) goto 500
         bx0=xmean+((rand()-0.5)*xbox)
         by0=ymean+((rand()-0.5)*xbox)
         bz0=zmean+((rand()-0.5)*ybox)
         ! This lines are not giving the best fit
         !bx0=xg+((rand()-0.5)*xbox)
         !by0=yg+((rand()-0.5)*xbox)
         !bz0=zg+((rand()-0.5)*ybox)
       ! bz0=(rand()-0.5)*zbox
       ! by0=(rand()-0.5)*ybox
       ! bz0=(rand()-0.5)*zbox
        if( ms2 .gt. 1 )call rotate( ms2,bx,by,bz )
        do i=1,ms2   
          xyz(1,i) = bx(i) + bx0
          xyz(2,i) = by(i) + by0
          xyz(3,i) = bz(i) + bz0  
        end do
        do jat = 1,ms2  
        do iout=1,ns1   
        do iat = 1,ms1   
          dx=xyz(1,jat)-sx(iout,iat)
          dy=xyz(2,jat)-sy(iout,iat)
          dz=xyz(3,jat)-sz(iout,iat)
          dx=dx-xbox*nint(dx/xbox)
          dy=dy-ybox*nint(dy/ybox)
          dz=dz-zbox*nint(dz/zbox)
          d2=dx**2+dy**2+dz**2
          if(d2.lt.rclose) goto 600
        end do
        end do
        end do
        do jat  = 1,ms2   
        do iout = ns1+1,ir   
        do iat  = 1,ms2   
          dx=xyz(1,jat)-sx(iout,iat)
          dy=xyz(2,jat)-sy(iout,iat)
          dz=xyz(3,jat)-sz(iout,iat)
          dx=dx-xbox*nint(dx/xbox)
          dy=dy-ybox*nint(dy/ybox)
          dz=dz-zbox*nint(dz/zbox)
          d2=dx**2+dy**2+dz**2
          if(d2.lt.rclose) goto 600
        end do
        end do
        end do
c 
        ir=ir+1
        do i=1,ms2   
          io=io+1
          write(8,100) io,bn(i),sul2,ir,(xyz(j,i),j=1,3),1.,0.
        end do
         do iat = 1,ms2   
          sx(ir,iat)=xyz(1,iat)
          sy(ir,iat)=xyz(2,iat)
          sz(ir,iat)=xyz(3,iat)
         end do
c 
         if( mod(ir,100) .eq. 0)print*,' Counts = ',ir,' out of ',ntmol 
c
  600   continue
      end do
      stop ' not all water molecules were generated'
  500 continue
      endif 
c ---------- 3rd one 
      if( ns3 .ne. 0 )then
      do irtry=1,1000*ns3    
        if(ir.eq.ns1+ns2+ns3) goto 599
         cx0=xmean+((rand()-0.5)*xbox)
         cy0=ymean+((rand()-0.5)*xbox)
         cz0=zmean+((rand()-0.5)*ybox)
         ! this Lines are not giving the best fit
         !cx0=xg+((rand()-0.5)*xbox)
         !cy0=yg+((rand()-0.5)*xbox)
         !cz0=zg+((rand()-0.5)*ybox)
        !cx0=(103.057+10)*rand()
        !cy0=(76.355+10)*rand()
        !cz0=(81.393+10)*rand()
        !cx0=(rand()-0.5)*xbox
        !cy0=(rand()-0.5)*ybox
        !cz0=(rand()-0.5)*zbox
        if( ms3 .gt. 1 )call rotate( ms3,cx,cy,cz )
        do i=1,ms3   
          xyz(1,i) = cx(i) + cx0
          xyz(2,i) = cy(i) + cy0
          xyz(3,i) = cz(i) + cz0
        end do
        do jat = 1,ms3  
        do iout=1,ns1   
        do iat = 1,ms1   
          dx=xyz(1,jat)-sx(iout,iat)
          dy=xyz(2,jat)-sy(iout,iat)
          dz=xyz(3,jat)-sz(iout,iat)
          dx=dx-xbox*nint(dx/xbox)
          dy=dy-ybox*nint(dy/ybox)
          dz=dz-zbox*nint(dz/zbox)
          d2=dx**2+dy**2+dz**2
          if(d2.lt.rclose) goto 699
        end do
        end do
        end do
        do jat  = 1,ms3   
        do iout = ns1+1,ns1+ns2
        do iat  = 1,ms2   
          dx=xyz(1,jat)-sx(iout,iat)
          dy=xyz(2,jat)-sy(iout,iat)
          dz=xyz(3,jat)-sz(iout,iat)
          dx=dx-xbox*nint(dx/xbox)
          dy=dy-ybox*nint(dy/ybox)
          dz=dz-zbox*nint(dz/zbox)
          d2=dx**2+dy**2+dz**2
          if(d2.lt.rclose) goto 699
        end do
        end do
        end do
c 
        do jat  = 1,ms3   
        do iout = ns1+ns2+1,ir
        do iat  = 1,ms3   
          dx=xyz(1,jat)-sx(iout,iat)
          dy=xyz(2,jat)-sy(iout,iat)
          dz=xyz(3,jat)-sz(iout,iat)
          dx=dx-xbox*nint(dx/xbox)
          dy=dy-ybox*nint(dy/ybox)
          dz=dz-zbox*nint(dz/zbox)
          d2=dx**2+dy**2+dz**2
          if(d2.lt.rclose) goto 699
        end do
        end do
        end do
c 
        ir=ir+1
        do i=1,ms3   
          io=io+1
          write(8,100) io,cn(i),sul3,ir,(xyz(j,i),j=1,3),1.,0.
        end do
         do iat = 1,ms3   
          sx(ir,iat)=xyz(1,iat)
          sy(ir,iat)=xyz(2,iat)
          sz(ir,iat)=xyz(3,iat)
         end do
c 
         if( mod(ir,100) .eq. 0)print*,' Counts = ',ir,' out of ',ntmol 
c
  699   continue
      end do
      stop ' not all of the molecules were generated'
  599 continue
      endif
c --------- 4th one
      if( ns4 .ne. 0 )then 
      do irtry=1,1000*ns4    
        if(ir.eq.ns1+ns2+ns3+ns4) goto 588
        ex0=(rand()-0.5)*xbox
        ey0=(rand()-0.5)*ybox
        ez0=(rand()-0.5)*zbox
        if( ms4 .gt. 1 )call rotate( ms4,ex,ey,ez )
        do i=1,ms4   
          xyz(1,i) = ex(i) + ex0
          xyz(2,i) = ey(i) + ey0
          xyz(3,i) = ez(i) + ez0
        end do
        do jat = 1,ms4  
        do iout=1,ns1   
        do iat = 1,ms1   
          dx=xyz(1,jat)-sx(iout,iat)
          dy=xyz(2,jat)-sy(iout,iat)
          dz=xyz(3,jat)-sz(iout,iat)
          dx=dx-xbox*nint(dx/xbox)
          dy=dy-ybox*nint(dy/ybox)
          dz=dz-zbox*nint(dz/zbox)
          d2=dx**2+dy**2+dz**2
          if(d2.lt.rclose) goto 688
        end do
        end do
        end do
        do jat  = 1,ms4   
        do iout = ns1+1,ns1+ns2
        do iat  = 1,ms2   
          dx=xyz(1,jat)-sx(iout,iat)
          dy=xyz(2,jat)-sy(iout,iat)
          dz=xyz(3,jat)-sz(iout,iat)
          dx=dx-xbox*nint(dx/xbox)
          dy=dy-ybox*nint(dy/ybox)
          dz=dz-zbox*nint(dz/zbox)
          d2=dx**2+dy**2+dz**2
          if(d2.lt.rclose) goto 688
        end do
        end do
        end do
c 
        do jat  = 1,ms4   
        do iout = ns1+ns2+1,ns1+ns2+ns3
        do iat  = 1,ms3   
          dx=xyz(1,jat)-sx(iout,iat)
          dy=xyz(2,jat)-sy(iout,iat)
          dz=xyz(3,jat)-sz(iout,iat)
          dx=dx-xbox*nint(dx/xbox)
          dy=dy-ybox*nint(dy/ybox)
          dz=dz-zbox*nint(dz/zbox)
          d2=dx**2+dy**2+dz**2
          if(d2.lt.rclose) goto 688
        end do
        end do
        end do
c 
        do jat  = 1,ms4   
        do iout = ns1+ns2+ns3+1,ir 
        do iat  = 1,ms4   
          dx=xyz(1,jat)-sx(iout,iat)
          dy=xyz(2,jat)-sy(iout,iat)
          dz=xyz(3,jat)-sz(iout,iat)
          dx=dx-xbox*nint(dx/xbox)
          dy=dy-ybox*nint(dy/ybox)
          dz=dz-zbox*nint(dz/zbox)
          d2=dx**2+dy**2+dz**2
          if(d2.lt.rclose) goto 688
        end do
        end do
        end do
c 
        ir=ir+1
        do i=1,ms4   
          io=io+1
          write(8,100) io,dn(i),sul4,ir,(xyz(j,i),j=1,3),1.,0.
        end do
         do iat = 1,ms4   
          sx(ir,iat)=xyz(1,iat)
          sy(ir,iat)=xyz(2,iat)
          sz(ir,iat)=xyz(3,iat)
         end do
c 
         if( mod(ir,100) .eq. 0)print*,' Counts = ',ir,' out of ',ntmol 
c
  688   continue
      end do
      stop ' not all water molecules were generated'
  588 continue
      endif
c
  100 format('ATOM',2x,i5,2x,a4,a4,1x,i4,4x,3f8.3,2f6.2)
c
      write(8,'(a)') 'END'
c
      stop
      end
c ------------------------------------------------------------------
      subroutine rotate(n,x,y,z)
      dimension x(50),y(50),z(50),xo(50),yo(50),zo(50) 
c
      write(15,*) ' n = ', n 
      do j = 1,n
       xo(j) = x(j)
       yo(j) = y(j)
       zo(j) = z(j)
      end do
c
      q1  = (rand() - 0.5)*2
      q2  = (rand() - 0.5)*2
      q3  = (rand() - 0.5)*2
      q4  = (rand() - 0.5)*2
      fct = sqrt( q1*q1 + q2*q2 + q3*q3 + q4*q4 )
      q1  = q1/fct
      q2  = q2/fct
      q3  = q3/fct
      q4  = q4/fct
c     write(6,*) 'norm = ',sqrt(q1**2 + q2**2 + q3**2 + q4**2)
c     
      s1   = q1**2
      s2   = q2**2
      s3   = q3**2
      s4   = q4**2
c
      a11  = s1 + s2 - s3 - s4
      a22  = s1 - s2 + s3 - s4
      a33  = s1 - s2 - s3 + s4
      a12  = 2*(q2*q3 + q1*q4)
      a13  = 2*(q2*q4 - q1*q3)
      a21  = 2*(q2*q3 - q1*q4)
      a23  = 2*(q3*q4 + q1*q2)
      a31  = 2*(q2*q4 + q1*q3)
      a32  = 2*(q3*q4 - q1*q2)
c
      do j = 1,n
       x(j)  = a11*xo(j) + a12*yo(j) + a13*zo(j) 
       y(j)  = a21*xo(j) + a22*yo(j) + a23*zo(j)
       z(j)  = a31*xo(j) + a32*yo(j) + a33*zo(j)
      end do
c
      return
      end
 
             
