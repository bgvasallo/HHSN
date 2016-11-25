      implicit none
      real*8 A!surface
      real*8 Cm!capacitance in F/m**2
      real*8 pk!probability of potasium mechanism
      real*8 pna!probability of sodium mechanism
      real*8 pl!probability of other ions mechanisms
      real*8 pp!probability of pump mechanism
      real*8 p!total probability
      real*8 V!membrane voltage
      real*8 I!membrane current
      real*8 Vk
      real*8 Vna
      real*8 Vl
      real*8 Ik
      real*8 Ina
      real*8 Il
      real*8 Ip
      real*8 gkm
      real*8 gnam
      real*8 glm
      real*8 glk
      real*8 glna
      real*8 n
      real*8 m
      real*8 h
      real*8 alphan
      real*8 alpham
      real*8 alphah
      real*8 betan
      real*8 betam
      real*8 betah
      real*8 dt
      real*8 TT!total time of simulation
      real*8 Tr!lapsed time
      real*8 r!random number
      real*8 RR!ideal gas constant
	  real*8 KB!Boltzmann constant
      real*8 T!temperature
      real*8 F!Faraday's constant
      real*8 e!electron charge
      real*8 dt2
      real*8 Tr2
      real*8 ts!
      real*8 tsi!time at which excitation is applied
      real*8 Am!!Am membrane sufrace (m**2)
      real*8 wi!wi volumen of intracelular space (m**3)
      real*8 we!we volumen of extracelular space (m**3)
      real*8 Na!avogadro number
      real*8 rho!max pump current
      real*8 cko,cnao,clo,cki,cnai,cli!concentration of each kind of ion outside and inside the membran
      real*8 ckoeq,cnaoeq,cloeq,ckieq,cnaieq,clieq ! Equilibrium concentrations
      real*8 Q
      real*8 Iapp! applied current
      real*8 Inasp
      real*8 Tr3
      real*8 Ikf,Inaf,Ilf,Ipf! currents
      real*8 time1,time2
      real*8 lk

      integer*8 ii!integers for iterations
      integer*8 iw!period to write
      integer*8 iw2! period to write 2
      integer*8 id!period to print
      integer*8 is!number of iterations to spike
      integer*8 zk,zna,zl!valencia de cada ion
      integer*8 nik,nina,nil,nip,nni!number of events
      integer*8 spikestart,spikend
      integer*8 jj! integer for loops
      integer*8 iitime1,iitime2! iteration for time1 and 2
       integer*8 nreg


      parameter(F=9.6485d4)!C/mol
      parameter(RR=8.3145d0)!V*C/(mol*K)
      parameter(e=1.60217657d-19)!Coulombs
      parameter(Na=6.022d23)!1/mol,ions/mol
	  parameter(KB=1.38064852d-23)! m2 kg s-2 K-1
      parameter(zk=1)
      parameter(zna=1)
      parameter(zl=-1)

      common/integers/ii,iw,id,is,nik,nina,nil,nreg
     &,nip,nni,spikestart,spikend,jj
     &,iw2,iitime1,iitime2
      common/reals/A,Cm,pk,pna,pl,pp,p,V,I,Vk,Vna,Vl,Ik,Ina,Il
     &,Ip,gkm,gnam,glm,glk,glna,n,m,h,dt,TT,Tr,r,T
     &,Inasp,dt2,Tr2,ts,tsi,Am,wi,we,rho,cko,cnao
     &,clo,cki,cnai,cli,Q,Iapp,Tr3
     &,time1,time2,lk,ckoeq,cnaoeq,cloeq,ckieq,cnaieq,clieq
