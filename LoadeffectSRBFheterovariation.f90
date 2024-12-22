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
      !read load green functions for indirect load effect ���ɸ��ֺ���(���Ӱ��)
	!�������շ��� read the load love numbers
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
      !���������߶ȸ����ļ���
      !Input the heterogeneous geodetic variation record time series file name.
      write(dtmfl,*)'dtm3m.dat'
      !�����ز���վ���¼ʱ�������ļ���
      !Input the calculation surface height grid file name.
      write(obstmsqdfl,*)'heterobstm.txt'
      !����վ��ƽ�����lnth(m)���ۼ�SRBF�ƽ�����
      !Input the mean distance (m) between geodetic sites and cumulative SRBF approach times
      lnth=5.d3;itern=1
      !������SRBF���ۻ�SRBF1������krbf,krbf1������order,order1
      !Input type and order the main SRBF and cumulative SRBF1
      krbf=0;order=0;krbf1=1;order1=0
      !������SRBF���ĺ��ۻ�SRBF1���ĵ����þ���(m)Input the action distance (m) of main SRBF center and cumulative SRBF1 center
      dr=120.d3;dr1=60.d3
      !����SRBF���õº�����С������minN,maxN Input minimum and maximum degree of SRBF Legendre expansion
      minN=9;maxN=900
      !�����ۻ�SRBF1���õº�����С������minN1,maxN1 Input minimum and maximum degree of cumulative SRBF1 Legendre expansion
      minN1=720;maxN1=1800
      !����SRBF��SRBF1��Bjerhammar���������� Input the Bjerhammar sphere burial depth for SRBF and SRBF1
      dpth=1.d3;dpth1=5.d3
      !���ÿɵ��ؼ�������ͼ��乱����
      !Set type (1~6) of adjustable variation and contribution rate
      cntrlknd=3;kwgh=1.d0!cntrlknd=1~6;kwgh=1.d0 ��ζ������
      !���뷨���̽��㷽�� the method of the solution of normal equation
      !knd=1 LU triangular decomposition method, =2 Cholesky decomposition, =3 least square qr decomposition,
      !   =4 Minimum norm singular value decomposition
      knd=1!knd=1-LU�ֽ�,2-Cholesky�ֽ�,3-��С����QR�ֽ�,4-��С��������ֵ�ֽ�
      hepch=2!��¼ʱ��ͷ�ļ��׸�����ʱ������� Column ordinal number of the first epoch time in header
      frow=7!��¼ʱ�����״β�������� Column ordinal number of the first variation in record
      kndrow=6!�������������� Column ordinal number of the variation type in record
      wghrow=5!�����Ȩֵ����� !Column ordinal number of the variation weight in record
      !����������ۺϳ�����inp(22) merge input parameters into inp(12)
      inp(1)=hepch;inp(2)=frow;inp(3)=kndrow;inp(4)=wghrow;inp(5)=lnth
      inp(6)=itern;inp(7)=knd;inp(8)=cntrlknd;inp(9)=kwgh
      inp(10)=dr;inp(11)=dr1;inp(12)=krbf;inp(13)=krbf1;inp(14)=order;inp(15)=order1
      inp(16)=minN;inp(17)=minN1;inp(18)=maxN;inp(19)=maxN1;inp(20)=dpth;inp(21)=dpth1
      open(unit=8,file=obstmsqdfl,status="old",iostat=status)
      if(status/=0)goto 902
      read(8,'(a)') headmn
      call PickReclong(headmn,kln,rec,sn)
      if(sn<hepch)then
         close(8);goto 902!ȱ������Ԫʱ����Ϣ Missing epoch time information
      endif
      close(8)
      tmlth=sn-hepch+1!�����ʱ�򳤶� length of geodetic variation record time series
      write(*, *)"    Begin compulation......"
      do i=frow,frow+tmlth-1 !i=��ǰ���������� column ordinal number of the current variations in record
        call LoadestmateSRBF(dtmfl,obstmsqdfl,flv,i,inp)
        write(*, 104)nint(rec(hepch+i-frow))
      enddo
902   continue
      pause
104   format('     Load EWH and load effects computed at Epoch time',I12)
      end
