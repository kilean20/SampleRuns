! Simulation of the particle-core model
      implicit double precision (a-h,o-z)
      parameter(np=32)
      dimension y(2,np),yv(2)
      common/cons/eta,eta2,eta21
      common/xtra/r0
      open(unit=54,file='halo.in',status='old')
      open(unit=15,file='env.out',status='unknown')
      open(unit=19,file='halo.out',status='unknown')
      read(54,*)tfin,ns,nprnt,nskip,eta,r0
      h=tfin/ns
      t=0.
!  initial conditions:
      call ics(y,yv)
      call prnt(t,y,yv)
      twopi=asin(1.0)*4.0
!  integrate the equations:
      do 100 n=1,ns
      tt=t
      call symp4(h,tt,y,yv)
      if(mod(n-nskip,nprnt).eq.0)call prnt(t,y,yv)
  100 continue
      stop
      end

      subroutine ics(y,yv)
! routine to specify initial conditions
      implicit double precision (a-h,o-z)
      parameter(np=32)
      dimension y(2,np),yv(2)
      common/cons/eta,eta2,eta21
      common/xtra/r0
      eta2=eta*eta
      eta21=1.0-eta2
! envelope:
      yv(1)=r0
      yv(2)=0.
! trajectories:
      np2=np/2
      delx=3.0/np2
      delp=3.0/np2
      do 10 j=1,np2
      xval=j*delx
      y(1,j)=xval
      y(2,j)=0.
   10 continue
      do 20 j=1,np2
      pval=j*delp
      y(1,np2+j)=0.
      y(2,np2+j)=xval
   20 continue
      return
      end

      subroutine map1(tau,t,y,yt,yv,ytv)
      implicit double precision (a-h,o-z)
      parameter(np=32)
      dimension y(2,np),yv(2),yt(2,np),ytv(2)
      do 100 j=1,np
      yt(1,j)=y(1,j)+y(2,j)*tau
      yt(2,j)=y(2,j)
  100 continue
      ytv(1)=yv(1)+yv(2)*tau
      ytv(2)=yv(2)
      t=t+tau
      return
      end

      subroutine map2(tau,t,y,yt,yv,ytv)
      implicit double precision (a-h,o-z)
      parameter(np=32)
      dimension y(2,np),yv(2),yt(2,np),ytv(2)
      common/cons/eta,eta2,eta21
! trajectories:
      do 100 j=1,np
      if(abs(y(1,j)).lt.yv(1))then
        vcon=y(1,j)/yv(1)**2
      else
        vcon=1.0/y(1,j)
      endif
      yt(1,j)=y(1,j)
      yt(2,j)=y(2,j)+tau*eta21*vcon-tau*y(1,j)
  100 continue
! envelope:
      ytv(1)=yv(1)
      ytv(2)=yv(2)-tau*(yv(1)-eta2/yv(1)**3-eta21/yv(1))
      return
      end
!
      subroutine prnt(t,y,yv)
      implicit double precision (a-h,o-z)
      parameter(np=32)
      dimension y(2,np),yv(2)
      write(15,101)t,yv(1),yv(2)
  101 format(6(1x,e14.7))
      do 20 j=1,np
   20 write(19,101)t,y(1,j),y(2,j)
      return
      end

      subroutine symp4(h,t,y,yv)
      implicit double precision (a-h,o-z)
      parameter(np=32)
      dimension y(2,np),yv(2)
      dimension yt(2,np),ytv(2)
      alf=1.0-2.0**(1.0/3.0)
      tau2=h/(1.0+alf)
      tau1=tau2*0.50
      tau3=alf*tau1
      tau4=(alf-1.0)*tau2
! trajectories and envelope:
      tt=t
      call map1(tau1,tt,y,yt,yv,ytv)
      call map2(tau2,tt,yt,y,ytv,yv)
      call map1(tau3,tt,y,yt,yv,ytv)
      call map2(tau4,tt,yt,y,ytv,yv)
      call map1(tau3,tt,y,yt,yv,ytv)
      call map2(tau2,tt,yt,y,ytv,yv)
      call map1(tau1,tt,y,yt,yv,ytv)
! finish up:
      y(1,:)=yt(1,:)
      yv=ytv
      t=t+h
      return
      end
