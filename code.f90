!******************************************
!*** Micro & Nano Flows Course Project ****
!******        A. F. Forughi         ******
!*****           (3/2011)             *****
!**** Sharif University of Technology  ****
!***      Molecular Dynamics (MD)       ***
!***      Liquid Argon Simulation       ***
!***   Reduced Units ** Periodic BC     ***
!****     Lennard-Jones potential      ****
!******************************************


implicit none

real(8) rx,ry,vx,vy,ax,ay,f,fx,fy,trefr,eps,sig,rnd1,rnd2,kb,rc,dist,distx,disty
real(8) m,h,tnewr,tstat,pot,vxave,vyave,kinet,g,dr,mvd,dv
integer i,n,mt,ii,iii,hn,nm,ndr,nmt
allocatable rx(:),ry(:),vx(:),vy(:),ax(:),ay(:),fx(:),fy(:),g(:),hn(:),mvd(:)
n=100  ! number of molecules
nmt=10000  ! number of time steps
allocate (rx(1:n),ry(1:n),vx(1:n),vy(1:n),ax(1:n),ay(1:n),fx(1:n),fy(1:n))
open(1,file="kinet.txt")
open(2,file="mvd.txt")
open(3,file="rdf.txt")
open(4,file="pot.txt")
open(5,file="tot.txt")

!********************************** Constants: ***********************************
eps=1.65e-21 !epsilon (j)
sig=3.4e-10  !sigma (m)
kb=1.3805e-23  !bultzman constant (j/k)
trefr=100.0/(eps/kb) ! reduced temperature ()
rc=2.5d0 !R cut off
m=1.0d0  ! reduced mass ()
h=0.005d0 !Reduced time step ()

!********************************* Initialize: ***********************************
!do i=1,n
!   rx(i)=(10.0/8.0)*((i-aint(real(i)/8.0)-1.0)+0.5)
!   ry(i)=(10.0/8.0)*((aint(real(i)/8.0)-1.0)+0.5)
!enddo

do i=1,n
   rx(i)=(i-1.0)-10.0*aint(0.1*(i-1.0))+0.5
   ry(i)=aint(0.1*(i-1.0))+0.5
   call RANDOM_NUMBER(rnd1)
   call RANDOM_NUMBER(rnd2)
   if ((2*rnd2-1)>0.0) then
      rnd2=1.0d0
   elseif ((2*rnd2-1)<0.0) then
      rnd2=-1.0d0
   else
      rnd2=0.0d0
   endif
   vx(i)=dsqrt(trefr*2.0/1.0)*(2.0*rnd1-1.0)!0.031838d0!
   vy(i)=dsqrt((trefr*2.0/1.0)-(vx(i))**2)*rnd2!0.0d0!
   !print*,vx(i),"  ",vy(i)
enddo

!******** Initial Temp. calculation ***********
tnewr=0
do i=1,n
   tnewr=tnewr+1.0*((vx(i))**2+(vy(i))**2)/(2.0*n)
enddo
print*,trefr," -> ",tnewr


tstat=1.0d0

vxave=0.0d0
vyave=0.0d0

do mt=1,nmt!  Time loop

   tnewr=0.0d0
   kinet=0.0d0
   !****************************** Forcing: ********************************
   fx=0.0d0
   fy=0.0d0
   pot=0.0d0

   do i=1,n
      do ii=1,n
         if (i/=ii) then
            if (dabs(rx(i)-rx(ii))<=2.5) then
               distx=rx(i)-rx(ii)
            elseif (dabs(rx(i)-rx(ii))>=7.5) then
               if (rx(i)<rx(ii)) then
                  distx=rx(i)-(rx(ii)-10.0)
               else
                  distx=rx(i)-(rx(ii)+10.0)
               endif
            else
               goto 11
            endif

            if (dabs(ry(i)-ry(ii))<=2.5) then
               disty=ry(i)-ry(ii)
            elseif (dabs(ry(i)-ry(ii))>=7.5) then
               if (ry(i)<ry(ii)) then
                  disty=ry(i)-(ry(ii)-10.0)
               else
                  disty=ry(i)-(ry(ii)+10.0)
               endif
            else
               goto 11
            endif

            dist=dsqrt(distx**2+disty**2)
            f=24.0*(2.0*((dist)**(-13.0))-((dist)**(-7.0))) !Force magnitude
            fx(i)=fx(i)+f*(distx/dist)
            fy(i)=fy(i)+f*(disty/dist)

            pot=pot+4.0*(dist**(-12.0)-dist**(-6.0))

            11 continue
         endif
      enddo
      !print*,"f:",i,fx(i),fy(i)
      !pause
   enddo
   !******************************** Forcing Done ***********************************


   !******************** A @ current time + V & R on next time step: ****************
   do i=1,n
      ax(i)=fx(i)/m
      ay(i)=fy(i)/m



      vx(i)=tstat*vx(i)+h*ax(i)-vxave !-vave= for momentum balance
      vy(i)=tstat*vy(i)+h*ay(i)-vyave


      rx(i)=rx(i)+h*vx(i)+0.5*h*h*ax(i)
      ry(i)=ry(i)+h*vy(i)+0.5*h*h*ay(i)
      !if (i==62) print*,"rr",rx(62)

      ! Periodic BC:
      if (rx(i)>=10.0) rx(i)=rx(i)-10.0
      if (rx(i)<0.0)   rx(i)=rx(i)+10.0

      if (ry(i)>=10.0) ry(i)=ry(i)-10.0
      if (ry(i)<0.0)   ry(i)=ry(i)+10.0

   enddo

   !Average velocity:
   vxave=0.0d0
   vyave=0.0d0
   do i=1,n
      vxave=vxave+vx(i)/n
      vyave=vyave+vy(i)/n
   enddo


   !New temperature for the Thermostat + kinetic energy:
   do i=1,n
      tnewr=tnewr+1.0*((vx(i)-vxave)**2+(vy(i)-vyave)**2)/(2.0*n)
      kinet=kinet+0.5*1.0*((vx(i))**2+(vy(i))**2)
   enddo
   tstat=dsqrt(trefr/tnewr) !Thermostat


   if (mod(mt-1,1000)==0) then    ! print on screen
      print*,mt,"Tnew:",tnewr," tstat:",tstat
      print*,"Vave=",vxave,vyave,"pot=",pot
      print*,"r1:",rx(44),ry(44)," v1:",vx(44),vy(44)
      !print*," "
   endif

   ! Writing:
   write(1,*) mt,kinet         !Write kinetic energy to kinet.txt
   write(4,*) mt,pot           !Write potential energy to pot.txt
   write(5,*) mt,(kinet+pot)   !Write total energy to tot.txt
   !******************************** A+V+R Done ***********************************

enddo! To the next time step



!******************** Radial distribution function (RDF): *************************
ndr=40
dr=0.1d0 !must ndr*dr=4.0
allocate (g(1:ndr),hn(0:ndr))
g=0.0d0  !RDF
hn=0     !counter

do iii=0,ndr
   do i=1,n
      do ii=1,n
         if (i/=ii) then

            if (dabs(rx(i)-rx(ii))<=(iii*dr)) then
               distx=rx(i)-rx(ii)
            elseif (dabs(rx(i)-rx(ii))>=(10.0-iii*dr)) then
               if (rx(i)<rx(ii)) then
                  distx=rx(i)-(rx(ii)-10.0)
               else
                  distx=rx(i)-(rx(ii)+10.0)
               endif
            else
               goto 12
            endif

            if (dabs(ry(i)-ry(ii))<=(iii*dr)) then
               disty=ry(i)-ry(ii)
            elseif (dabs(ry(i)-ry(ii))>=(10.0-iii*dr)) then
               if (ry(i)<ry(ii)) then
                  disty=ry(i)-(ry(ii)-10.0)
               else
                  disty=ry(i)-(ry(ii)+10.0)
               endif
            else
               goto 12
            endif

            dist=dsqrt(distx**2+disty**2)
            if ((dist<iii*dr).and.(dist>=(iii-1)*dr)) then!
               hn(iii)=hn(iii)+1
            endif

            12 continue
         endif
      enddo
   enddo
enddo

nm=0 !Count all
do iii=1,ndr
   nm=nm+hn(iii)
enddo

do iii=1,ndr
   g(iii)=g(iii)+(4.0*4.0*3.141592654)*hn(iii)/(3.141592654*(nm*2)*(iii-0.5)*dr*dr)
enddo

do iii=1,ndr !write to rdf.txt
   write(3,*) iii*dr,g(iii)
enddo
!************************************ RDF Done************************************


!*********************** molecular velocity distribution (MVD): ******************
allocate (mvd(-50:50))
mvd=0.0d0 ! Velocity Histogram
dv=0.25d0
do ii=-20,20
   do i=1,n
      if ((vy(i)<(ii*dv)).and.(vy(i)>=((ii-1)*dv))) then
         mvd(ii)=mvd(ii)+1.0/n
      endif
   enddo
enddo

do ii=-20,20
   write(2,*) ii*dv,mvd(ii)
enddo
!********************************** MVD done! ***********************************


print*,"Last kinet&pot:",kinet,pot
!pause
stop
end
