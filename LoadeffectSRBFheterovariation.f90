!  LoadeffectSRBFheterovariation.f90 
!
!  FUNCTIONS:
!  LoadeffectSRBFheterovariation - Entry point of console application.
!
!****************************************************************************
      program LoadeffectSRBFheterovariation
      implicit none
	character*800::dtmfl,obstmsqdfl
	character*80000::line,headmn
	integer kln,sn,astat(15),tmlth,i,j,n,minN,minN1,maxN,maxN1
	integer hepch,frow,kndrow,wghrow,itern,cntrlknd,knd,krbf,krbf1,order,order1
	real*8::rec(8000),flv(40000,3),inp(24),lnth,kwgh,dr,dr1,dpth,dpth1
	integer::status=0
!---------------------------------------------------------------------
      !read load green functions for indirect load effect 负荷格林函数(间接影响)
	!读负荷勒夫数 read the load love numbers
      flv=0.d0
      open(unit=8,file='Love_load_cm.dat',status="old",iostat=status)
      if(status/=0)goto 902 
      do i=1,6
        read(8,'(a)') line
      enddo
      n=0
	do while(.not.eof(8))
        n=n+1
	  read(8,*,end=904)i,(flv(n,j),j=1,3)
      enddo
904   close(8)
      !输入计算面高度格网文件名
      !Input the heterogeneous geodetic variation record time series file name.
      write(dtmfl,*)'dtm3m.dat'
      !输入大地测量站点记录时间序列文件名
      !Input the calculation surface height grid file name.
      write(obstmsqdfl,*)'heterobstm.txt'
      !输入站点平均间距lnth(m)和累计SRBF逼近次数
      !Input the mean distance (m) between geodetic sites and cumulative SRBF approach times
      lnth=5.d3;itern=1
      !输入主SRBF和累积SRBF1的类型krbf,krbf1及次数order,order1
      !Input type and order the main SRBF and cumulative SRBF1
      krbf=0;order=0;krbf1=1;order1=0
      !输入主SRBF中心和累积SRBF1中心的作用距离(m)Input the action distance (m) of main SRBF center and cumulative SRBF1 center
      dr=120.d3;dr1=60.d3
      !输入SRBF勒让德函数最小最大阶数minN,maxN Input minimum and maximum degree of SRBF Legendre expansion
      minN=9;maxN=900
      !输入累积SRBF1勒让德函数最小最大阶数minN1,maxN1 Input minimum and maximum degree of cumulative SRBF1 Legendre expansion
      minN1=720;maxN1=1800
      !输入SRBF和SRBF1的Bjerhammar球面埋藏深度 Input the Bjerhammar sphere burial depth for SRBF and SRBF1
      dpth=1.d3;dpth1=5.d3
      !设置可调控监测量类型及其贡献率
      !Set type (1~6) of adjustable variation and contribution rate
      cntrlknd=3;kwgh=1.d0!cntrlknd=1~6;kwgh=1.d0 意味不调控
      !输入法方程解算方法 the method of the solution of normal equation
      !knd=1 LU triangular decomposition method, =2 Cholesky decomposition, =3 least square qr decomposition,
      !   =4 Minimum norm singular value decomposition
      knd=1!knd=1-LU分解,2-Cholesky分解,3-最小二乘QR分解,4-最小范数奇异值分解
      hepch=2!记录时序头文件首个采样时刻列序号 Column ordinal number of the first epoch time in header
      frow=7!记录时序中首次采样列序号 Column ordinal number of the first variation in record
      kndrow=6!监测量类型列序号 Column ordinal number of the variation type in record
      wghrow=5!监测量权值列序号 !Column ordinal number of the variation weight in record
      !将输入参数综合成向量inp(22) merge input parameters into inp(12)
      inp(1)=hepch;inp(2)=frow;inp(3)=kndrow;inp(4)=wghrow;inp(5)=lnth
      inp(6)=itern;inp(7)=knd;inp(8)=cntrlknd;inp(9)=kwgh
      inp(10)=dr;inp(11)=dr1;inp(12)=krbf;inp(13)=krbf1;inp(14)=order;inp(15)=order1
      inp(16)=minN;inp(17)=minN1;inp(18)=maxN;inp(19)=maxN1;inp(20)=dpth;inp(21)=dpth1
      open(unit=8,file=obstmsqdfl,status="old",iostat=status)
      if(status/=0)goto 902
      read(8,'(a)') headmn
      call PickReclong(headmn,kln,rec,sn)
      if(sn<hepch)then
         close(8);goto 902!缺采样历元时刻信息 Missing epoch time information
      endif
      close(8)
      tmlth=sn-hepch+1!监测量时序长度 length of geodetic variation record time series
      write(*, *)"    Begin compulation......"
      do i=frow,frow+tmlth-1 !i=当前监测量列序号 column ordinal number of the current variations in record
        call LoadestmateSRBF(dtmfl,obstmsqdfl,flv,i,inp)
        write(*, 104)nint(rec(hepch+i-frow))
      enddo
902   continue
      pause
104   format('     Load EWH and load effects computed at Epoch time',I12)
      end
