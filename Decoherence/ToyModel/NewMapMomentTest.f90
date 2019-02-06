
     Program StandardMoment
     implicit none
     double precision, allocatable:: particles(:,:)
     double precision:: G,Gmean,q,p,sum,qf,pf
     double precision:: eps0,q0,delta
     integer:: n,j,npart,niterate
     open(unit=33,file='partcl.data',status='old')
     open(unit=34,file='moment.data',status='new')
     open(unit=35,file='output.data',status='new')
     open(unit=36,file='partcl_in.data',status='new')
     read(33,*) npart
     allocate(particles(npart,2))
     write(*,*) 'Enter the desired number of iterates:'
     read(*,*) niterate
     eps0 = 0.01d0 !nominal
     q0 = 0.5d0   !nominal
     Gmean = 0.0d0
!  Read data, compute initial value of Gmean:
     do j=1,npart
        read(33,*) particles(j,:)
        q = particles(j,1)
        p = particles(j,2)
!  Scale and translate input data
        q = q*dsqrt(eps0) + q0
        p = p*dsqrt(eps0)
        particles(j,1) = q
        particles(j,2) = p
        write(36,*) q,p  !nominal
!  Compute mean
        Gmean = Gmean + G(q,p)
     enddo
     Gmean = Gmean/dble(npart)
     write(34,*) 0,Gmean
!  Perform iterations
     do n=1,niterate
        Gmean = 0.0d0
        do j=1,npart
           q = particles(j,1)
           p = particles(j,2)
           call map(q,p,qf,pf)
           Gmean = Gmean + G(qf,pf)
           particles(j,1) = qf
           particles(j,2) = pf
        enddo
        Gmean = Gmean/dble(npart)
        write(34,*) n,Gmean
     enddo
     do j=1,npart
        q = particles(j,1)
        p = particles(j,2)
        write(35,*) q,p   !nominal
     enddo
     deallocate(particles)
!     close(33)
!     close(34)
!     close(35)
!     close(36)
     End

     subroutine map(q,p,qf,pf)
     implicit none
     double precision, intent(in):: q,p
     double precision, intent(out):: qf,pf
     double precision:: sum,pi,ratio,mod2pi,K
     double precision:: psi,alpha,action,beta,arg
! Map parameters
     psi = 0.3d0   !Zero for std map
     alpha = 0.1d0   !One for std map
! Unperturbed, integrable map:
     action = (q**2+p**2)/2.0d0
     arg = psi + alpha*action !nominal
     qf = q*dcos(arg) + p*dsin(arg)
     pf = -q*dsin(arg) + p*dcos(arg)
     end subroutine map

     function G(q,p)
     implicit none
     double precision, intent(in):: q,p
     double precision:: G,delta
     G = q
!     G = p
     end function G

     function mod2pi(x)
     implicit none
     double precision, intent(in):: x
     double precision:: pi,ratio,mod2pi
     pi = 4.0d0*datan(1.0d0)
     ratio = x/(2.0d0*pi)
     mod2pi = x - 2.0d0*pi*aint(ratio)
!  Be careful about negative values to ensure result in [0,2pi]
     if(mod2pi<0.0d0) then
       mod2pi = 2.0d0*pi + mod2pi
     endif
     end function mod2pi
