      program SHHIC
      include 'definition.h'

!     read values
      call reader
!     initialize some variables
      call initialize1
!     random number
      CALL init_random_seed()
!     open files to write on
      call openoutputfiles

      do while (Tr.lt.T)

!     screen

      if (mod(ii,id).eq.0) then

      write(*,*)dint(Tr*1000/T)

      end if

!     write output files

      if (ii.ge.iitime1.and.ii.le.iitime2) then

      if (mod(ii-1,iw2).eq.0) then

      call write1

      end if

      else

      if (mod(ii-1,iw).eq.0) then

      call write1

      end if

      end if

!     nerst potentials

      call nerstcalc

!     write nerst potentials

      if (ii.ge.iitime1.and.ii.le.iitime2) then

      if (mod(ii-1,iw2).eq.0) then

      call write2

      end if

      else

      if (mod(ii-1,iw).eq.0) then

      call write2

      end if

      end if


!     currents and associated probabilities

      call currentandprobability

!     initialize some variables

      call initialize2

      do while (Tr2.lt.dt)

      call secondloop

      end do

!     updates values

      call update

!     check if concentrations are positive to prevent errors

      If (cko.lt.0) then

      print*,'cko concentration less than 0'

      exit

      end if

      If (cki.lt.0) then

      print*,'cki concentration less than 0'

      exit

      end if  

      If (cnao.lt.0) then

      print*,'cnao concentration less than 0'

      exit

      end if

      If (cnai.lt.0) then

      print*,'cnai concentration less than 0'

      exit

      end if

      If (clo.lt.0) then

      print*,'clo concentration less than 0'

      exit

      end if

      If (cli.lt.0) then

      print*,'cli concentration less than 0'

      exit

      end if

!     write quantities

      if (ii.ge.iitime1.and.ii.le.iitime2) then

      if (mod(ii-1,iw2).eq.0) then

      call write3

      end if

      else

      if (mod(ii-1,iw).eq.0) then

      call write3

      end if

      end if



      end do

!     close files
      call closer


      stop
      end







!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                    FUNCTIONS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      real*8 function alphan(V)
      implicit none
      real*8 V
      alphan=(0.01d6*(V+34.d-3))/(-dexp(-(V+34.d-3)/10.d-3)+1.d0)
      return
      end

      real*8 function betan(V)
      implicit none
      real*8 V
      betan=0.125d3*dexp(-(V+44.d-3)/80.d-3)
      return
      end

      real*8 function alpham(V)
      implicit none
      real*8 V
      alpham=(0.1d6*(V+30.d-3))/(-dexp(-(V+30.d-3)/10.d-3)+1.d0)
      return
      end

      real*8 function betam(V)
      implicit none
      real*8 V
      betam=4.d3*dexp(-(V+55.d-3)/18.d-3)
      return
      end

      real*8 function alphah(V)
      implicit none
      real*8 V
      alphah=0.07d3*dexp(-(V+44.d-3)/20.d-3)
      return
      end

      real*8 function betah(V)
      implicit none
      real*8 V
      betah=1.d3/(dexp(-(V+14.d-3)/10.d-3)+1)
      return
      end



      real*8 function Inaf(V,m,h,Vna,gnam,glna)
      implicit none
      real*8 V,m,h,Vna
      real*8 gnam,glna
      Inaf=(gnam*m**3*h+glna)*(V-Vna)
      return
      end

      real*8 function Ikf(V,n,Vk,gkm,glk)
      implicit none
      real*8 V,n,Vk
      real*8 gkm,glk
      Ikf=(gkm*n**4+glk)*(V-Vk)
      return
      end

      real*8 function Ilf(V,Vl,glm)
      implicit none
      real*8 V,Vl
      real*8 glm
      Ilf=glm*(V-Vl)
      return
      end

      real*8 function Ipf(cnai,cko,rho)
      implicit none
      real*8 cnai,cko
      real*8 rho
      Ipf=rho/((1.d0+dexp((25.d0-cnai)/3.d0))*(1.d0
     &+dexp(5.5d0-cko)))
      return
      end


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!                   SUBROUTINES
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


      subroutine reader
      include 'definition.h'

!     values

      open(20,file='Values_SHHIC')

      read(20,*)T
      read(20,*)Cm
      read(20,*)I
      read(20,*)gkm
      read(20,*)gnam
      read(20,*)glm
      read(20,*)glk
      read(20,*)glna
      read(20,*)lk
      read(20,*)nreg
      read(20,*)dt
      read(20,*)Am
      read(20,*)wi
      read(20,*)we
      read(20,*)TT
      read(20,*)rho
      read(20,*)iw
      read(20,*)Q
      read(20,*)iw2
      read(20,*)time1
      read(20,*)time2

!     read diferent concentrations to simulate

      open(101,file='Concentrations')

      read(101,*)cko
      read(101,*)cki
      read(101,*)cnao
      read(101,*)cnai
      read(101,*)clo
      read(101,*)cli

!     read excitation parameters      

      open(102,file='Excitation')

      read(102,*)tsi
      read(102,*)ts
      read(102,*)Iapp

      return
      end









      subroutine nerstcalc
      include 'definition.h'

      Vk=(KB*TT/e/zk)*dlog(cko/cki)
      Vna=(KB*TT/e/zna)*dlog(cnao/cnai)
      Vl=(KB*TT/e/zl)*dlog(clo/cli)

      return
      end








      subroutine openoutputfiles
      include 'definition.h'

!     open files to write output data

      open(10,file='Vm')
      open(17,file='Ik')
      open(18,file='Ina')
      open(19,file='Il')
      open(35,file='Ipump')
      open(22,file='n')
      open(23,file='m')
      open(24,file='h')
      open(25,file='cko')
      open(26,file='cki')
      open(27,file='cnao')
      open(28,file='cnai')
      open(29,file='clo')
      open(30,file='cli')
      open(31,file='Vk')
      open(32,file='Vna')
      open(33,file='Vl')
      open(41,file='Iapp')

      return
      end











      subroutine initialize1
      include 'definition.h'

!     times involved in Iapp

      spikestart = dint(tsi/dt)

      spikend = dint(ts/dt) + spikestart
      id=int(T/dt/1000.)

!     writing timing

      iitime1=dint(time1/dt)
      iitime2=dint(time2/dt)

!     Initial voltage

      V=Q/Cm
  
      cnaoeq=cnao
      cnaieq=cnai
      ckoeq=cko
      ckieq=cki
      cloeq=clo
      clieq=cli

!     initial values of n,m and h

      n=(alphan(V)/(alphan(V)+betan(V)))
      m=(alpham(V)/(alpham(V)+betam(V)))
      h=(alphah(V)/(alphah(V)+betah(V)))

!     initializing variables

      Tr=0.d0
      ii=1

      return
      end










      subroutine initialize2
      include 'definition.h'

      Tr2=0.d0

      nik=0
      nina=0
      nil=0
      nni=0
      nip=0

      return
      end











      subroutine currentandprobability
      include 'definition.h'

!     new n
      n=n+dt*(alphan(V)*(1.d0-n)-betan(V)*n)

!     new m
      m=m+dt*(alpham(V)*(1.d0-m)-betam(V)*m)

!     new h
      h=h+dt*(alphah(V)*(1.d0-h)-betah(V)*h)

!     check if it is time for Iapp

      if (ii.gt.spikestart.and.ii.lt.spikend) then

!     Iapp applied here
      		Inasp=Iapp

!     currents

   		Ik=Ikf(V,n,Vk,gkm,glk)
    		Ina=Inaf(V,m,h,Vna,gnam,glna) - Inasp
		Il=Ilf(V,Vl,glm)
    		Ip=Ipf(cnai,cko,rho)

!     if it is not time for Iapp

      else 

   		Ik=Ikf(V,n,Vk,gkm,glk)
    		Ina=Inaf(V,m,h,Vna,gnam,glna)
		Il=Ilf(V,Vl,glm)
    		Ip=Ipf(cnai,cko,rho)

      end if

      if (ii.eq.spikend) then
      		Inasp=0.d0
      end if

!     probabilities

      pk=dabs(Ik*Am/e)
      pna=dabs(Ina*Am/e)
      pl=dabs(Il*Am/e)
      pp=dabs(Ip*Am/e)
      p=pk+pna+pl+pp

      pk=pk/p
      pna=pna/p
      pl=pl/p
      pp=pp/p

      return
      end








      subroutine currentcalcwithoutexcitation
      include 'definition.h'

      Ik=(gkm*n**4+glk)*(V-Vk)
      Ina=(gnam*m**3*h+glna)*(V-Vna)
      Il=glm*(V-Vl)
      Ip=rho/((1.d0+dexp((25.d0-cnai)/3.d0))*(1.d0+dexp(5.5d0-cko)))

      return
      end













      subroutine currentcalcwithexcitation
      include 'definition.h'

      Ik=(gkm*n**4+glk)*(V-Vk)
      Ina=(gnam*m**3*h+glna)*(V-Vna) + Inasp
      Il=glm*(V-Vl)
      Ip=rho/((1.d0+dexp((25.d0-cnai)/3.d0))*(1.d0+dexp(5.5d0-cko)))

      return
      end














      subroutine secondloop
      include 'definition.h'

!     time between events

      call random_number(r)

      dt2=-(1.d0/p)*dlog(r)

!     transcurred time
 
      Tr2=Tr2+dt2

      if (Tr2.lt.dt) then 

!     kind of event

      call random_number(r)

      if (r.lt.pk) then! potasium current event

      nik=nik+1

      else if (r.lt.pk+pna) then! sodium current event

      nina=nina+1

      else if (r.lt.pk+pna+pl) then! leak current event

      nil=nil+1

      else
	
      nip=nip+1

      end if

      end if

      return
      end
      
   







      subroutine update
      include 'definition.h'

      ii=ii+1
      Tr=Tr+dt
      nik=dsign(1.d0*nik,Ik)
      nina=dsign(1.d0*nina,Ina)
      nil=dsign(1.d0*nil,Il)
      nip=dsign(1.d0*nip,Ip)
      nni=nik+nina+nil+nip
      cnao=cnao+(nina+3.d0*nip)/Na/we
      cnai=cnai-(nina+3.d0*nip)/Na/wi

!Regulation term in potasium concentration
      cko=cko+(nik-2.d0*nip)/Na/we+lk*nreg*(4.-cko)*dt

      cki=cki-(nik-2.d0*nip)/Na/wi
      clo=clo-nil/Na/we
      cli=cli+nil/Na/wi


	 V=V+(dt/Cm)*(I-nni*e/Am/(dt))

      return
      end














      subroutine write1
      include 'definition.h'

      write(10,*)Tr,V
      write(22,*)Tr,n
      write(23,*)Tr,m
      write(24,*)Tr,h
      write(25,*)Tr,cko
      write(26,*)Tr,cki
      write(27,*)Tr,cnao
      write(28,*)Tr,cnai
      write(29,*)Tr,clo
      write(30,*)Tr,cli

      return
      end


















      subroutine write2
      include 'definition.h'

      write(31,*)Tr,Vk
      write(32,*)Tr,Vna
      write(33,*)Tr,Vl

      return
      end









      subroutine write3
      include 'definition.h'

      write(17,*)Tr,nik*e/Am/dt
      write(18,*)Tr,nina*e/Am/dt
      write(19,*)Tr,nil*e/Am/dt
      write(35,*)Tr,nip*e/Am/dt
      write(41,*)Tr,Inasp


      return
      end











      subroutine closer
      include 'definition.h'

      close(10)
      close(17)
      close(18)
      close(19)
      close(20)
      close(22)
      close(23)
      close(24)
      close(25)
      close(26)
      close(27)
      close(28)
      close(29)
      close(30)
      close(31)
      close(32)
      close(33)

      close(35)
      close(41)

      close(101)
      close(102)

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


      
