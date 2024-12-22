      subroutine EqHestens(BPB,XX,nn,BPL)
!�����ݶȷ�����ͷ�����BPB.XX=BPL
!2005��3��28��,�´���
!---------------------------------------------------------------
      implicit none
	integer::nn
	real*8::BPB(nn,nn),BPL(nn),XX(nn)
	real*8,allocatable::R(:),RK(:),XXK(:),P(:),tmp(:),LL(:),RK0(:)
	integer*4::i,j,ii,si,sj
	real*8 rab,alfa,bita,para,rr,st
!-----------------------------------------------------------------------
	allocate(XXK(nn),R(nn),RK(nn),P(nn),tmp(nn),RK0(nn),LL(nn))
	do i=1,nn
	  XX(i)=BPL(i)/BPB(i,i)
	  LL(i)=1.d0/dsqrt(dabs(BPB(i,i)))
	enddo
	do i=1,nn
	  tmp(i)=0.d0
	  do j=1,nn
	    tmp(i)=tmp(i)+BPB(i,j)*XX(j)
	  enddo
	enddo
	R=BPL-tmp     !�в�����
	RK=LL*R;P=LL*RK                 !������������
      rr=0.d0;
	do i=1,nn
	  rr=rr+R(i)**2
	enddo
	rr=dsqrt(rr/dble(nn))
	rab=rr*1.d-4                      !������������
	do while(rr<rab)
!	�������alfa
	  tmp(1:nn)=0.d0
	  do i=1,nn
	    tmp(i)=0.d0
	    do j=1,nn
	      tmp(i)=tmp(i)+BPB(i,j)*P(j)
	    enddo
	  enddo
	  alfa=0.d0; para=0.d0
	  do i=1,nn
	    alfa=alfa+RK(i)*RK(i)
	    para=para+tmp(i)*P(i)
	  enddo
	  alfa=alfa/para
!     -----------------
	  XXK=XX+alfa*P
	  RK=RK0-alfa*LL*tmp
	  R=BPL-tmp;rr=0.d0
	  do i=1,nn
	    rr=rr+R(i)**2
	  enddo
	  rr=dsqrt(rr/dble(nn))
!	�������bita
	  bita=0.d0; para=0.d0
	  do i=1,nn
	    bita=bita+RK(i)*RK(i)
	    para=para+RK0(i)*RK0(i)
	  enddo
	  bita=bita/para
!     ----------------
	  P=LL*RK+bita*P
	  RK0=RK
	enddo
	do i=1,nn
	  BPL(i)=dsqrt(1.d0/BPB(i,i))
	enddo
2001	deallocate(XXK,R,RK,P,tmp,RK0,LL)
	return
      end
!
!*****************************************************************
!
      subroutine EqJordon(BPB,xx,nn,BPL)
!Guass-Jordon��������BPB.XX=BPL���
!���:xx-δ֪���Ľ�,BPL-xx�ĵ�λȨ��ƽ����
!2006��12��11��,�´���
!---------------------------------------------------------------
      implicit none
	integer*4 nn,maxk
	real*8::BPB(nn,nn),BPL(nn),xx(nn),dl
	real*8,allocatable::dt(:)
	real*4,allocatable::tmp(:,:)
	integer,allocatable::chg(:)
	integer::i,j,k
!-----------------------------------------------------------------------
      xx=0.d0
	allocate(tmp(nn,nn),chg(nn),dt(nn))
	do i=1,nn !��λ��
	  tmp(i,i)=1.d0
	enddo
      !�����Ǿ���
	do k=1,nn-1
	  maxk=k
	  do i=k+1,nn !������dabs(BPB(i,i))
	    if(dabs(BPB(i,k))>dabs(BPB(k,k)))maxk=i
	  enddo
	  dt=BPB(k,:);BPB(k,:)=BPB(maxk,:);BPB(maxk,:)=dt
	  dl=BPL(k);BPL(k)=BPL(maxk);BPL(maxk)=dl
	  chg(k)=maxk
	  do i=k+1,nn
	    do j=k+1,nn
	      BPB(i,j)=BPB(i,j)-BPB(k,j)/BPB(k,k)*BPB(i,k)
	    enddo
	    do j=1,nn
	      tmp(i,j)=tmp(i,j)-tmp(k,j)/BPB(k,k)*BPB(i,k)
	    enddo
	    BPL(i)=BPL(i)-BPL(k)/BPB(k,k)*BPB(i,k)
	  enddo
	enddo
	!�Խ��߾���
	do k=nn-1,1,-1
	  do i=1,k
	    do j=1,nn
	      tmp(i,j)=tmp(i,j)-tmp(k+1,j)/BPB(k+1,k+1)*BPB(i,k+1)
	    enddo
	    BPL(i)=BPL(i)-BPL(k+1)/BPB(k+1,k+1)*BPB(i,k+1)
	  enddo
	enddo
	do i=1,nn
	  xx(i)=BPL(i)/BPB(i,i)
	  BPL(i)=dsqrt(dabs(tmp(i,i)/BPB(i,i)))
	enddo
	!�����й���
	do i=nn-1,1,-1
	  dt=tmp(i,:);tmp(i,:)=tmp(chg(i),:);tmp(chg(i),:)=dt
	  dl=xx(i);xx(i)=xx(chg(i));xx(chg(i))=dl
	  dl=BPL(i);BPL(i)=BPL(chg(i));BPL(chg(i))=dl
	enddo
	BPB=tmp
	deallocate(tmp,chg,dt)
	return
      end
!
!**************************************************************************************
!
      subroutine EqCholesky(BPB,XX,nn,BPL)
!ƽ������������BPB.XX=BPL���
!���:BPB-BPB�����,xx-δ֪���Ľ�,BPL-xx�ĵ�λȨ��ƽ����
!2006��12��11��,�´���
!---------------------------------------------------------------
      implicit none
	integer*4 nn
	real*8::BPB(nn,nn),BPL(nn),XX(nn),dl
	real*8,allocatable::yy(:)
	real*8,allocatable::tmp(:,:)
	integer::i,j,k
!-----------------------------------------------------------------------
	allocate(tmp(nn,nn),yy(nn))
	do i=1,nn !L
	  tmp(i,i)=0.d0
	enddo
      !BPB��LLT�ֽ�
	do j=1,nn
	  dl=0.d0
	  do k=1,j-1
	   dl=tmp(j,k)**2
	  enddo
	  tmp(j,j)=dsqrt(BPB(j,j)-dl)
	  do i=j+1,nn
	    dl=0.d0
	    do k=1,j-1
	      dl=dl+tmp(i,k)*tmp(j,k)
	    enddo
	    tmp(i,j)=(BPB(i,j)-dl)/tmp(j,j)
	  enddo
	enddo
	!��Ly=b��LTx=y
	do i=1,nn
	  dl=0.d0
	  do k=1,i-1
	   dl=tmp(i,k)*yy(k)
	  enddo
	  yy(i)=(BPL(i)-dl)/tmp(i,i)
	enddo
	do i=nn,1,-1
	  dl=0.d0
	  do k=i+1,nn
	   dl=tmp(k,i)*xx(k)
	  enddo
	  xx(i)=(yy(i)-dl)/tmp(i,i)
	enddo
	do i=1,nn
	  BPL(i)=dsqrt(1.d0/dabs(BPB(i,i)))
	enddo
	deallocate(tmp,yy)
	return
      end
!
!***********************************************************************
!
      subroutine EqueSVD(BPB,XX,nn,BPL)
!SVD������ͷ�����BPB.XX=BPL
!2019��7��29��,�´���
!---------------------------------------------------------------
      USE f95_precision, ONLY: WP => DP
      USE lapack95, ONLY: gesvd,gelss
      implicit none
	integer::nn
	real*8::BPB(nn,nn),BPL(nn),XX(nn)
	real*8,allocatable::s(:)
	integer*4::i,ki
!-----------------------------------------------------------------------
      allocate(s(nn))
      s=0.d0
      call gesvd(BPB,s)!f95
      do ki=2,nn
         if(dabs(s(ki))<dabs(s(1))*1.d-8)goto 9020
      enddo
9020  do i=ki,nn
        s(i)=0.d0
      enddo
      call gelss(BPB,BPL,nn,s)
      xx=BPL
2001	deallocate(s)
	return
      end
!
!*****************************************************************
!
      subroutine RidgeEstimate(BPB,XX,nn,BPL)
!����ƽ���ͷ�����BPB.XX=BPL
!2012��12��20��,�´���
!---------------------------------------------------------------
      implicit none
	integer::i,j,k,nn,kk,info,lwork,nm
	real*8::BPB(nn,nn),BPL(nn),XX(nn),work(3*nn),w(nn)
	real*8::maxp,minp,kp(80000),nta
	real*8::rx(200000),ry(80000),rx1(80000),ry1(80000),rx2(80000),ry2(80000)
	real*8,allocatable::B0(:,:),BPB0(:,:),L0(:)
!-----------------------------------------------------------------------
      nm=500
	allocate(B0(nn,nn),BPB0(nn,nn),L0(nn))
      BPB0=BPB
      call dsyev('N','U', nn, BPB0, nn, w, work, -1, info)
      do i=1,nn
         L0(i)=BPB0(i,i)
      enddo
      L0=dabs(L0)
	maxp=L0(1);minp=L0(1)
	do i=1,nn
	  if(maxp<L0(i))maxp=L0(i)
	  if(minp>L0(i))minp=L0(i)
	enddo
	if(maxp/minp<1.d3)goto 2001
	BPB0=BPB
	rx=0.d0;ry=0.d0;nta=0.d0
	do i=1,nn
	  nta=nta+L0(i)**2
	enddo
	nta=dsqrt(nta/dble(nn))*2.d-4  !nta-1.d-3�ϴ������죬����
	do k=1,nm   !1000
        B0=BPB;L0=BPL
	  do i=1,nn
	    BPB(i,i)=BPB(i,i)+dble(k)*nta
	  enddo
	  call EqHestens(BPB,xx,nn,BPL)
	  BPB=B0;BPL=L0
	  L0=matmul(BPB,xx)-BPL
	  do i=1,nn
	    rx(k)=rx(k)+dlog(xx(i)**2)
	    ry(k)=ry(k)+dlog(L0(i)**2)
	  enddo
	enddo
      do k=2,nm-1  !999
	  rx1(k)=(rx(k+1)-rx(k-1))/2.d0
	  ry1(k)=(ry(k+1)-ry(k-1))/2.d0
	enddo
      do k=3,nm-2  !998
	  rx2(k)=(rx1(k+1)-rx1(k-1))/2.d0
	  ry2(k)=(ry1(k+1)-ry1(k-1))/2.d0
	  kp(k)=dabs(rx1(k)*ry2(k)-rx2(k)*ry1(k))/(rx1(k)**2+ry1(k)**2)**1.5d0
	enddo
	kk=3;maxp=kp(k)
	do k=4,nm-2  !998
	  if(maxp<kp(k))then
	    kk=k;maxp=kp(k)
	  endif
	enddo
	BPB=BPB0
	do i=1,nn
	   BPB(i,i)=BPB(i,i)+dble(kk)*nta
	enddo
2001	call EqJordon(BPB,xx,nn,BPL)
	deallocate(B0,BPB0,L0)
	return
      end
!
!******************************************************************
!
      subroutine Equsolve(BB,xx,nn,BL,knd,bf)
!����ͷ�����BB.xx=BL
!knd=1LU�ֽ�,2Cholesky�ֽ�,3��С����QR�ֽ�,4��С��������ֵ�ֽ�
!bf(8)-�������
!---------------------------------------------------------------
      USE OMP_LIB
      implicit none
	integer::nn,knd,nm,inf,lwk,rk,astat(20)
	real*8::BB(nn,nn),BL(nn),xx(nn),bf(8),rnd
      real*8,allocatable::s(:),wk(:)
	integer,allocatable::iwk(:),ipv(:)
	integer::i,ki
!-----------------------------------------------------------------------
      nm=nn;bf=0.d0;rk=0;xx=BL
 	allocate(s(nn), stat=astat(1))
 	allocate(wk(nn*nn), stat=astat(2))
 	allocate(iwk(nn*nn), stat=astat(3))
 	allocate(ipv(nn), stat=astat(4))
	if (sum(astat(1:4)) /= 0) then
        bf(1)=1.d0;return
      endif
      lwk=-1;s=0.d0
      if(knd==3)call dgels('No transpose',nm,nm,1,BB,nm,xx,nm,wk,lwk,inf)
      if(knd==4)call dgelsd(nm,nm,1,BB,nm,xx,nm,s,-1.d0,rk,wk,lwk,iwk,inf)
      if(knd==5)call dgelss(nm,nm,1,BB,nm,xx,nm,s,-1.d0,rk,wk,lwk,inf)
      nm=nn;lwk = nm**nm
      if(knd==1)call dgesv(nm,1,BB,nm,ipv,xx,nm,inf)
      if(knd==2)call dposv('Upper',nm,1,BB,nm,xx,nm,inf)
      if(knd==3)call dgels('No transpose',nm,nm,1,BB,nm,xx,nm,wk,lwk,inf)
      if(knd==4)call dgelsd(nm,nm,1,BB,nm,xx,nm,s,-1.d0,rk,wk,lwk,iwk,inf)
      if(knd==5)call dgelss(nm,nm,1,BB,nm,xx,nm,s,-1.d0,rk,wk,lwk,inf)
      if( inf >0 ) then
        bf(2)=1.d0;goto 2001!����ʧ��-,������������¼���
      endif
      bf(3)=rk
2001	deallocate(s,wk,iwk,ipv)
	return
      end
