      program SHH
      implicit none
      real*8 A!surface
      real*8 Cm!capacitance in F/m**2
      real*8 pk!probability of potasium mechanism
      real*8 pna!probability of sodium mechanism
      real*8 pl!probability of other ions mechanisms
      real*8 p!total probability
      real*8 V!initial potential value
      real*8 I!inicial current value
      real*8 Vk
      real*8 Vna
      real*8 Vl
      real*8 Ik
      real*8 Ina
      real*8 Il
      real*8 gkm
      real*8 gnam
      real*8 glm
      real*8 n
      real*8 m
      real*8 h
      real*8 alphan
      real*8 alpham
      real*8 alphah
      real*8 betan
      real*8 betam
      real*8 betah
      real*8 dt!timestep
      real*8 TT!total simulation time
      real*8 Tr!lapsed time
      real*8 r!random number
      real*8 RR!ideal gas constant
      real*8 T!temperature
      real*8 F!Faraday's constant
      real*8 e!electron charge
      real*8 pt
      real*8 Iapp!Iapp
      real*8 dt2
      real*8 Tr2
      real*8 ts!time for Iapp
      real*8 tsi!initial time for Iapp
      real*8 b!1/A
      
      integer*8 ii!integer for loops
      integer*8 iw!period to write
      integer*8 is!number of iterations to spike
      integer*8 id!period to print
      integer*8 z!atomic number

      integer*8 nik,nina,nil,nni,nk,nna,nl,nn!number of events
      integer*8 spikestart,spikend! iteration at which Iapp and possible spike starts/ends



      parameter(F=9.6485d4)!C/mol
      parameter(z=1)
      parameter(TT=292.15d0)!kelvin
      parameter(RR=8.3145d0)!V*C/(mol*K)
      parameter(e=1.60217657d-19)!Coulombs

!     values

      open(20,file='Values_SHH')

      read(20,*)T
      read(20,*)Cm
      read(20,*)b
      read(20,*)iw
      read(20,*)is
      read(20,*)id
      read(20,*)I
      read(20,*)V
      read(20,*)Vk
      read(20,*)Vna
      read(20,*)Vl
      read(20,*)gkm
      read(20,*)gnam
      read(20,*)glm
      read(20,*)Iapp
      read(20,*)dt
      read(20,*)tsi
      read(20,*)ts


      spikestart = dint(tsi/dt)

      spikend = dint(ts/dt) + spikestart

c     initial values of n,m and h

      n=(alphan(V)/(alphan(V)+betan(V)))
      m=(alpham(V)/(alpham(V)+betam(V)))
      h=(alphah(V)/(alphah(V)+betah(V)))

!     initializing variables

      Tr=0.d0
      ii=0
      nk=0
      nna=0
      nl=0
      nik=0
      nina=0
      nil=0
      nn=0

      CALL init_random_seed()

!     open files to write output data

      open(10,file='Vm')
      open(17,file='Ik')
      open(18,file='Ina')
      open(19,file='Il')
      open(21,file='Iapp')
      open(22,file='n')
      open(23,file='m')
      open(24,file='h')

      do while (Tr.lt.T)

!     print progress

      if (mod(ii,id).eq.0) then

      write(*,*)dint(Tr*100/T)

      end if


      if (abs(ii-spikestart).lt.1) then

      print*,'Iapp applied'

      I=Iapp

      end if

      if (abs(ii-spikend).lt.1) then

      print*,'end Iapp'


      I=0.d0

      end if

!     new n

      n=n+dt*(alphan(V)*(1.d0-n)-betan(V)*n)

!     new m

      m=m+dt*(alpham(V)*(1.d0-m)-betam(V)*m)

!     new h

      h=h+dt*(alphah(V)*(1.d0-h)-betah(V)*h)

!     current

      Ik=gkm*n**4*(V-Vk)
      Ina=gnam*m**3*h*(V-Vna)
      Il=glm*(V-Vl)

!     Probabilities

      pk=dabs(Ik/(e*b))
      pna=dabs(Ina/(e*b))
      pl=dabs(Il/(e*b))
      p=pk+pna+pl

      pk=pk/p
      pna=pna/p
      pl=pl/p

!     write output

      if (mod(ii,iw).eq.0) then

      write(10,*)Tr,V
      write(17,*)Tr,dsign(1.d0*nik,Ik)*e*b/dt
      write(18,*)Tr,dsign(1.d0*nina,Ina)*e*b/dt
      write(19,*)Tr,dsign(1.d0*nil,Il)*e*b/dt
      write(21,*)Tr,I
      write(22,*)Tr,n
      write(23,*)Tr,m
      write(24,*)Tr,h

      end if

!     initialize some variables

      Tr2=0.d0

      nik=0
      nina=0
      nil=0
      nni=0

      do while (Tr2.lt.dt)

!     time between events

      call random_number(r)

      dt2=-(1.d0/p)*dlog(r)

!     time
 
      Tr2=Tr2+dt2

!     choose kind of event

      call random_number(r)

      if (r.lt.pk) then! potasium current

      nik=nik+1

      else if (r.lt.pk+pna) then! sodium current

      nina=nina+1

      else! leak current

      nil=nil+1

      end if

      end do

!     updates values

      ii=ii+1
      Tr=Tr+dt
      nni=dsign(1.d0*nik,Ik)+dsign(1.d0*nina,Ina)+dsign(1.d0*nil,Il)
      nk=nk+dsign(1.d0*nik,Ik)
      nna=nna+dsign(1.d0*nina,Ina)
      nl=nl+dsign(1.d0*nil,Il)



      V=V+(dt/Cm)*(I-nni*e*b/dt)



      end do

      close(10)
      close(17)
      close(18)
      close(19)
      close(20)
      close(21)
      close(22)
      close(23)
      close(24)

      stop
      end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                    FUNCTIONS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      real*8 function alphan(V)
      implicit none
      real*8 V
      alphan=(0.01d6*(-V+10.d-3))/(dexp((-V+10.d-3)/10.d-3)-1)
      return
      end

      real*8 function betan(V)
      implicit none
      real*8 V
      betan=0.125d3*dexp(-V/80.d-3)
      return
      end

      real*8 function alpham(V)
      implicit none
      real*8 V
      alpham=(0.1d6*(-V+25.d-3))/(dexp((-V+25.d-3)/10.d-3)-1)
      return
      end

      real*8 function betam(V)
      implicit none
      real*8 V
      betam=4.d3*dexp(-V/18.d-3)
      return
      end

      real*8 function alphah(V)
      implicit none
      real*8 V
      alphah=0.07d3*dexp(-V/20.d-3)
      return
      end

      real*8 function betah(V)
      implicit none
      real*8 V
      betah=1.d3/(dexp((-V+30.d-3)/10.d-3)+1)
      return
      end
      
      
      
      
      
      
      subroutine init_random_seed()
      implicit none
      integer, allocatable :: seed(:)
      integer :: i, n, un, istat, dt(8), pid, t(2), s
      integer(8) :: count, tms, getpid
      call random_seed(size = n)
      allocate(seed(n))
      open(unit=un, file="/dev/urandom", access="stream",
     +  form="unformatted", action="read", status="old", iostat=istat)
      if (istat == 0) then
        read(un) seed
        close(un)
      else
        call system_clock(count)
        if (count /= 0) then
          t = transfer(count, t)
        else
          call date_and_time(values=dt)
          tms = (dt(1) - 1970) * 365_8 * 24 * 60 * 60 * 1000
     -         + dt(2) * 31_8 * 24 * 60 * 60 * 1000
     -         + dt(3) * 24 * 60 * 60 * 60 * 1000
     -         + dt(5) * 60 * 60 * 1000
     -         + dt(6) * 60 * 1000 + dt(7) * 1000
     -         + dt(8)
          t = transfer(tms, t)
        end if
        s = ieor(t(1), t(2))
        pid = getpid() + 1099279
        s = ieor(s, pid)
        if (n.ge.3) then
          seed(1) = t(1) + 36269
          seed(2) = t(2) + 72551
          seed(3) = pid
          if (n > 3) then
            seed(4:) = s + 37 * (/ (i, i = 0, n - 4) /)
          end if
        else
          seed = s + 37 * (/ (i, i = 0, n - 1 ) /)
        end if
      end if
      call random_seed(put=seed)
      end subroutine init_random_seed


      
      
      
      
      
      
      
      
      
      
