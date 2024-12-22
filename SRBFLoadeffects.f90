      subroutine SRBF6curve(RBF,flv,order,krbf,mpn,minN,maxN,NF,nta)
      !5种的RBF（krbf）曲线（NF+1）0~dr°
      !kobs=1高程异常,2扰动重力,3地面重力,4大地高,5正常高,6等效水高变化
      !krbf-order次0径向多极子1Possion小波
      implicit none
      integer::order,krbf,minN,maxN,NF
      real*8::RBF(NF+1,6),mpn(maxN-minN+1,NF+1),nta,flv(40000,3)
      integer::i,j,n,m,k
      real*8::onei,Bn,pi,RAD,CnmCalc,dh,dl,dk
!---------------------------------------------------------------
      pi=datan(1.d0)*4.d0;RAD=pi/180.d0
      m=order;RBF=0.d0
      onei=0.d0
      do i=1,NF+1!采用扰动位RFB系数归一化
        do n=minN,maxN
          k=n-minN+1;dh=flv(n,1);dl=flv(n,2);dk=flv(n,3)
          if(m==0)then
            if(krbf==0)Bn=1.d0/(2.d0*n+1.d0)
            if(krbf==1)Bn=1.d0
          else
            if(krbf==0)Bn=CnmCalc(n,m)*dexp(-dlog(nta)*m)/(2.d0*n+1.d0)
            if(krbf==1)Bn=(-dble(n)*dlog(nta))**m
          endif
          Bn=Bn*1.d-10
          RBF(i,1)=RBF(i,1)+(1.d0+dk)*Bn*dexp(dlog(nta)*n)*mpn(k,i)
          RBF(i,2)=RBF(i,2)+(1.d0+dk)*(n+1.d0)*Bn*dexp(dlog(nta)*(n-1.d0))*mpn(k,i)
          RBF(i,3)=RBF(i,3)+(1.d0+2.d0*dh/n-(n+1.d0)*dk/n)*(n+1.d0)*Bn*dexp(dlog(nta)*(n-1.d0))*mpn(k,i)
          RBF(i,4)=RBF(i,4)+dh*Bn*dexp(dlog(nta)*n)*mpn(k,i)
          RBF(i,5)=RBF(i,5)+(dh-1.d0-dk)*Bn*dexp(dlog(nta)*n)*mpn(k,i)
          RBF(i,6)=RBF(i,6)+(2.d0*n+1.d0)*Bn*dexp(dlog(nta)*n)*mpn(k,i)
          onei=onei+(1.d0+dk)*Bn*dexp(dlog(nta)*n)*mpn(k,i)
        enddo
      enddo
      RBF(1:NF+1,1:6)=RBF(1:NF+1,1:6)/onei
      end
!
!******************************************************************
!
      subroutine SRBF5curve(RBF,flv,order,krbf,mpn,minN,maxN,NF,nta)
      !5种的RBF（krbf）曲线（NF+1）0~dr°
      !kobs=1高程异常,2扰动重力,3地面重力,4大地高,5正常高变化
      !krbf-order次0径向多极子1Possion小波
      implicit none
      integer::order,krbf,minN,maxN,NF
      real*8::RBF(NF+1,5),mpn(maxN-minN+1,NF+1),nta,flv(40000,3)
      integer::i,j,n,m,k
      real*8::onei,Bn,pi,RAD,CnmCalc,dh,dl,dk
!---------------------------------------------------------------
      pi=datan(1.d0)*4.d0;RAD=pi/180.d0
      m=order;RBF=0.d0
      onei=0.d0
      do i=1,NF+1!采用扰动位RFB系数归一化
        do n=minN,maxN
          k=n-minN+1;dh=flv(n,1);dl=flv(n,2);dk=flv(n,3)
          if(m==0)then
            if(krbf==0)Bn=1.d0/(2.d0*n+1.d0)
            if(krbf==1)Bn=1.d0
          else
            if(krbf==0)Bn=CnmCalc(n,m)*dexp(-dlog(nta)*m)/(2.d0*n+1.d0)
            if(krbf==1)Bn=(-dble(n)*dlog(nta))**m
          endif
          Bn=Bn*1.d-10
          RBF(i,1)=RBF(i,1)+(1.d0+dk)*Bn*dexp(dlog(nta)*n)*mpn(k,i)
          RBF(i,2)=RBF(i,2)+(1.d0+dk)*(n+1.d0)*Bn*dexp(dlog(nta)*(n-1.d0))*mpn(k,i)
          RBF(i,3)=RBF(i,3)+(1.d0+2.d0*dh/n-(n+1.d0)*dk/n)*(n+1.d0)*Bn*dexp(dlog(nta)*(n-1.d0))*mpn(k,i)
          RBF(i,4)=RBF(i,4)+dh*Bn*dexp(dlog(nta)*n)*mpn(k,i)
          RBF(i,5)=RBF(i,5)+(dh-1.d0-dk)*Bn*dexp(dlog(nta)*n)*mpn(k,i)
          onei=onei+(1.d0+dk)*Bn*dexp(dlog(nta)*n)*mpn(k,i)
        enddo
      enddo
      RBF(1:NF+1,1:5)=RBF(1:NF+1,1:5)/onei
      end
!
!******************************************************************
!
      subroutine SRBF11all(RBF,flv,order,krbf,mpn,mdp,mp2,minN,maxN,NF,nta)
      !11种的RBF（krbf）曲线（NF+1）0~dr°
      !knd=1等效水高,2高程异常,3地面重力,4扰动重力,5地倾斜,6垂线偏差,7水平位移,8大地高,9正常高,10扰动梯度,11水平梯度
      !krbf-order次0径向多极子1Possion小波
      implicit none
      integer::order,krbf,minN,maxN,NF
      real*8::mpn(maxN-minN+1,NF+1),mdp(maxN-minN+1,NF+1),mp2(maxN-minN+1,NF+1)
      real*8::RBF(NF+1,11),onei,Bn,pi,RAD,CnmCalc,dh,dl,dk,nta,flv(40000,3)
      integer::i,j,n,m,k
!---------------------------------------------------------------
      pi=datan(1.d0)*4.d0;RAD=pi/180.d0
      m=order;RBF=0.d0
      onei=0.d0
      do i=1,NF+1!采用扰动位RFB系数归一化
        do n=minN,maxN
          k=n-minN+1;dh=flv(n,1);dl=flv(n,2);dk=flv(n,3)
          if(m==0)then
            if(krbf==0)Bn=1.d0/(2.d0*n+1.d0)
            if(krbf==1)Bn=1.d0
          else
            if(krbf==0)Bn=CnmCalc(n,m)*dexp(-dlog(nta)*m)/(2.d0*n+1.d0)
            if(krbf==1)Bn=(-dble(n)*dlog(nta))**m
          endif
          Bn=Bn*1.d-10
          RBF(i,1)=RBF(i,1)+(2.d0*n+1.d0)*Bn*dexp(dlog(nta)*n)*mpn(k,i)
          RBF(i,2)=RBF(i,2)+(1.d0+dk)*Bn*dexp(dlog(nta)*n)*mpn(k,i)
          RBF(i,3)=RBF(i,3)+(1.d0+2.d0*dh/n-(n+1.d0)*dk/n)*(n+1.d0)*Bn*dexp(dlog(nta)*(n-1.d0))*mpn(k,i)
          RBF(i,4)=RBF(i,4)+(1.d0+dk)*(n+1.d0)*Bn*dexp(dlog(nta)*(n-1.d0))*mpn(k,i)
          RBF(i,5)=RBF(i,5)+(1.d0+dk-dh)*Bn*dexp(dlog(nta)*n)*mdp(k,i)
          RBF(i,6)=RBF(i,6)+(1.d0+dk)*Bn*dexp(dlog(nta)*n)*mdp(k,i)
          RBF(i,7)=RBF(i,7)+dl*Bn*dexp(dlog(nta)*n)*mdp(k,i)
          RBF(i,8)=RBF(i,8)+dh*Bn*dexp(dlog(nta)*n)*mpn(k,i)
          RBF(i,9)=RBF(i,9)+(dh-1.d0-dk)*Bn*dexp(dlog(nta)*n)*mpn(k,i)
          RBF(i,10)=RBF(i,10)+(1.d0+dk)*(n+1.d0)*(n+2.d0)*Bn*dexp(dlog(nta)*(n-1.d0))*mpn(k,i)
          RBF(i,11)=RBF(i,11)+(1.d0+dk)*Bn*dexp(dlog(nta)*n)*mp2(k,i)
          onei=onei+(1.d0+dk)*Bn*dexp(dlog(nta)*n)*mpn(k,i)
        enddo
      enddo
      RBF(1:NF+1,1:11)=RBF(1:NF+1,1:11)/onei
      end
!
!******************************************************************
!
      subroutine SRBFone(RBFi,RBFn,flv,order,krbf,mpn,mdp,mp2,minN,maxN,NF,nta)
      !RBFi-11种RBF（krbf）空域曲线（NF+1）0~dr°;RBFn-RBF谱域曲线
      !knd=1等效水高,2高程异常,3地面重力,4扰动重力,5地倾斜,6垂线偏差,7水平位移,8大地高,9正常高,10扰动梯度,11水平梯度
      !krbf-order次0径向多极子1Possion小波
      implicit none
      integer::order,krbf,fknd,minN,maxN,NF
      real*8::mpn(maxN-minN+1,NF+1),mdp(maxN-minN+1,NF+1),mp2(maxN-minN+1,NF+1),nta,flv(40000,3)
      real*8::RBFi(NF+1,11),RBFn(maxN-minN+1,11),onei(11),Bn,pi,RAD,CnmCalc,dh,dl,dk,maxv,minv
      integer::i,j,n,m,k
!---------------------------------------------------------------
      pi=datan(1.d0)*4.d0;RAD=pi/180.d0
      m=order;RBFi=0.d0;RBFn=0.d0
      onei=0.d0
      do i=1,NF+1!采用扰动位RFB系数归一化
        do n=minN,maxN
          k=n-minN+1;dh=flv(n,1);dl=flv(n,2);dk=flv(n,3)
          if(m==0)then
            if(krbf==0)Bn=1.d0/(2.d0*n+1.d0)
            if(krbf==1)Bn=1.d0
          else
            if(krbf==0)Bn=CnmCalc(n,m)*dexp(-dlog(nta)*m)/(2.d0*n+1.d0)
            if(krbf==1)Bn=(-dble(n)*dlog(nta))**m
          endif
          Bn=Bn*1.d-10
          RBFi(i,1)=RBFi(i,1)+(2.d0*n+1.d0)*Bn*dexp(dlog(nta)*n)*mpn(k,i)
          RBFi(i,2)=RBFi(i,2)+(1.d0+dk)*Bn*dexp(dlog(nta)*n)*mpn(k,i)
          RBFi(i,3)=RBFi(i,3)+(1.d0+2.d0*dh/n-(n+1.d0)*dk/n)*(n+1.d0)*Bn*dexp(dlog(nta)*(n-1.d0))*mpn(k,i)
          RBFi(i,4)=RBFi(i,4)+(1.d0+dk)*(n+1.d0)*Bn*dexp(dlog(nta)*(n-1.d0))*mpn(k,i)
          RBFi(i,5)=RBFi(i,5)+(1.d0+dk-dh)*Bn*dexp(dlog(nta)*n)*mdp(k,i)
          RBFi(i,6)=RBFi(i,6)+(1.d0+dk)*Bn*dexp(dlog(nta)*n)*mdp(k,i)
          RBFi(i,7)=RBFi(i,7)+dl*Bn*dexp(dlog(nta)*n)*mdp(k,i)
          RBFi(i,8)=RBFi(i,8)+dh*Bn*dexp(dlog(nta)*n)*mpn(k,i)
          RBFi(i,9)=RBFi(i,9)+(dh-1.d0-dk)*Bn*dexp(dlog(nta)*n)*mpn(k,i)
          RBFi(i,10)=RBFi(i,10)+(1.d0+dk)*(n+1.d0)*(n+2.d0)*Bn*dexp(dlog(nta)*(n-1.d0))*mpn(k,i)
          RBFi(i,11)=RBFi(i,11)+(1.d0+dk)*Bn*dexp(dlog(nta)*n)*mp2(k,i)
          !-------------------------------------------------------------------------------
          RBFn(k,1)=RBFn(k,1)+(2.d0*n+1.d0)*Bn*dexp(dlog(nta)*n)*mpn(k,i)
          RBFn(k,2)=RBFn(k,2)+(1.d0+dk)*Bn*dexp(dlog(nta)*n)*mpn(k,i)
          RBFn(k,3)=RBFn(k,3)+(1.d0+2.d0*dh/n-(n+1.d0)*dk/n)*(n+1.d0)*Bn*dexp(dlog(nta)*n)*mpn(k,i)
          RBFn(k,4)=RBFn(k,4)+(1.d0+dk)*(n+1.d0)*Bn*dexp(dlog(nta)*n)*mpn(k,i)
          RBFn(k,5)=RBFn(k,5)+(1.d0+dk-dh)*Bn*dexp(dlog(nta)*n)*mdp(k,i)
          RBFn(k,6)=RBFn(k,6)+(1.d0+dk)*Bn*dexp(dlog(nta)*n)*mdp(k,i)
          RBFn(k,7)=RBFn(k,7)+dl*Bn*dexp(dlog(nta)*n)*mdp(k,i)
          RBFn(k,8)=RBFn(k,8)+dh*Bn*dexp(dlog(nta)*n)*mpn(k,i)
          RBFn(k,9)=RBFn(k,9)+(dh-1.d0-dk)*Bn*dexp(dlog(nta)*n)*mpn(k,i)
          RBFn(k,10)=RBFn(k,10)+(1.d0+dk)*(n+1.d0)*(n+2.d0)*Bn*dexp(dlog(nta)*n)*mpn(k,i)
          RBFn(k,11)=RBFn(k,11)+(1.d0+dk)*Bn*dexp(dlog(nta)*n)*mp2(k,i)
        enddo
      enddo
      do i=1,11
        maxv=maxval(RBFi(1:NF+1,i));minv=minval(RBFi(1:NF+1,i))
        RBFi(1:NF+1,i)=RBFi(1:NF+1,i)/(maxv-minv)!=(RBFi(1:NF+1,i)-minv)/(maxv-minv)
        maxv=maxval(RBFn(1:maxN-minN+1,i));minv=minval(RBFn(1:maxN-minN+1,i))
        RBFn(1:maxN-minN+1,i)=RBFn(1:maxN-minN+1,i)/(maxv-minv)!=(RBFn(1:maxN-minN+1,i)-minv)/(maxv-minv)
      enddo
      end
!
!******************************************************************************
!
      subroutine LegPn02(mpn,mdp,mp2,minN,maxN,NF,dr)
      !计算minN~maxN阶NF+1组勒让德函数[minN~maxN]×[0,dr]
	implicit none
	integer::minN,maxN,NF,nn,i,k,n
	real*8::dr,t,dt,pi,RAD
	real*8::mpn(maxN-minN+1,NF+1),mdp(maxN-minN+1,NF+1),mp2(maxN-minN+1,NF+1)
	real*8,allocatable::pn(:),dp1(:),dp2(:)
!---------------------------------------------------------------------------
      pi=datan(1.d0)*4.d0;RAD=pi/180.d0
      dt=dr/dble(NF)*RAD
      allocate(pn(maxN+2),dp1(maxN+2),dp2(maxN+2))
      do k=1,NF+1
        t=dcos(dble(k-1)*dt)
        call LegPn_dt2(pn,dp1,dp2,maxN+1,t)
        do n=minN,maxN
          i=n-minN+1
          mpn(i,k)=pn(n+1)
          mdp(i,k)=dp1(n+1)
          mp2(i,k)=dp2(n+1)
        enddo
      enddo
      deallocate(pn,dp1,dp2)
      end
!
!******************************************************************************
!
      subroutine RBFvalue(RBF,NF,dr,xx,val)
      !!由夹角xx°内插RBF值val
	implicit none
	integer::NF,i,k
	real*8::xx,dr,val,dt,RBF(NF+1),px,tmp
!---------------------------------------------------------------------------
      val=0.d0;px=0.d0;dt=dr/dble(NF)!间隔°
      k=nint(xx/dt)!与xx最近的数组序号
      do i=k-2,k+2
        if(i<1.or.i>NF+1)goto 1001
        tmp=(i-1.d0)*dt
        val=val+RBF(i)*tmp;px=px+tmp
1001    continue
      enddo
      if(px>1.d-12)val=val/px
      end
!
!******************************************************************
!
      subroutine drln(rln1,rln2,dln)
      !由两点球坐标计算距离与夹角°
      implicit none
      real*8::rln1(3),rln2(3),dln(2),XYZ1(3),XYZ2(3),L2
      integer::i,j,n,m,k
      real*8::onei,Bn,pi,RAD,CnmCalc
!---------------------------------------------------------------
      RAD=datan(1.d0)/45.d0
      call RLAT_XYZ(rln1,XYZ1)
      call RLAT_XYZ(rln2,XYZ2)
      L2=(XYZ1(1)-XYZ2(1))**2+(XYZ1(2)-XYZ2(2))**2+(XYZ1(3)-XYZ2(3))**2
      dln(1)=dsqrt(L2)
      dln(2)=dabs(dacos((rln1(1)**2+rln2(1)**2-L2)/2.d0/rln1(1)/rln2(1)))/RAD
      end
!
!*********************************************************
!
      real*8 function CnmCalc(n,m) 
      implicit none
      integer :: i,k,n,m 
      real*8 :: Cnm(10),tmp
      if(m==0)then
        CnmCalc=1.d0;return
      endif
      if(m==1)then
        CnmCalc=n;return
      endif
      if(m==2)then
        CnmCalc=n*(n-1)/2;return
      endif
      Cnm(1)=n;Cnm(2)=n*(n-1)/2  !n=m+1->C(m+1,m)   
      do i=3,m
        Cnm(3)=Cnm(1)+Cnm(2)
        Cnm(1)=Cnm(2);Cnm(2)=Cnm(1)
      enddo
      CnmCalc=Cnm(3)
      return
      end
