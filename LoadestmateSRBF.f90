      subroutine LoadestmateSRBF(dtmfl,obstmsqdfl,flv,dtrow,inp)
      USE f95_precision, ONLY: WP => DP
      USE lapack95
      implicit none
	character*800::dtmfl,obstmsqdfl,outfl,rntfl
	character*80000::line,headmn,nmstr(6),stm,pstr
      character(len=25)::str(8000)
	integer i,j,k,n,m,nn,mm,kk,ii,jj,kx,kw,sn,kln,astat(15),itern,mxk,sgn,kp,knd,kobs,ks,Kb,cntrlknd
	integer ki,kj,ni,nj,mk,krbf,krbf1,order,order1,minN,minN1,maxN,maxN1,Kt,mn,lvl
 	real*8 GRS(6),hd(6),wgh,RAD,rec(8000),st(6,4,5),unt(6),dlat,blat,r0,inp(24)
 	real*8 BLH(3),rln(3),ae,flv(40000,3),tsr(4),bf(8),tmp,lnth,ab(6),sta(6,4)
 	real*8 dpth,dpth1,rlnk(5),dr,dr1,factor,fk,GMr,gr,NFD(5),dln(2),rhd(4),rr,val,st0(6,4),nta,br
 	real*8 sinf,tt,cosa,sina,u15(15),mr,ff(15),fw(15),sf(15),wght(6)
      real*8 fi,fi1,la,la1,dla,dfi2,dnla2
	integer nlon,nlat,nk,hepch,frow,kndrow,wghrow,dtrow,tmlth,obsn(6),NF,gn,nd,edgn,lp(100)!
	real*8,allocatable::mpn(:,:),mdp(:,:),mp2(:,:),obs(:,:,:),obsp(:,:),chs(:,:,:),rst(:,:,:)
	real*8,allocatable::BPB(:,:,:),BPL(:,:),BB(:),xx(:),rlatlon(:,:),RBF(:,:),RBFi(:,:),RBFn(:,:)
	real*8,allocatable::APA(:,:),APL(:),tm(:),B2(:),sr(:),hgt(:,:),dl(:),lon(:,:),B15(:,:)
	integer,allocatable::nln(:),nrd(:,:),node(:),enode(:),gpnt(:,:)!格网、未知数序号,每个观测量有效节点序号
	character(len=25),allocatable::strnm(:,:)
	integer::status=0
!kobs=1高程异常mm,2扰动重力uGal,3地面重力uGal,4大地高mm,5正常高mm,6等效水高cm变化
!输出ewh,geoid,terrgrav,gravdist,grndtilt,vertdefl,horzdisp,elliphgt,orthohgt,gradient,horzgrad
!---------------------------------------------------------------------
    	RAD=datan(1.d0)/45.d0;BLH(3)=0.d0;mr=36.d5/RAD
      unt(1)=1.d3;unt(2)=1.d8;unt(3)=1.d8;unt(4)=1.d3;unt(5)=1.d3;unt(6)=1.d2
      u15(1)=1.d2;u15(2)=1.d3;u15(3)=1.d8;u15(4)=1.d8;u15(5)=mr;u15(6)=mr
      u15(7)=mr;u15(8)=mr;u15(9)=1.d3;u15(10)=1.d3;u15(11)=1.d3;u15(12)=1.d3
      u15(13)=1.d12;u15(14)=1.d9;u15(15)=1.d9!mE,E
      write(nmstr(1),'(a28,I3)')'GNSS-levelling (mm)',1
      write(nmstr(2),'(a28,I3)')'Gravity disturbance (μGal)',2
      write(nmstr(3),'(a28,I3)')'Ground gravity (μGal)',3
      write(nmstr(4),'(a28,I3)')'Ellipsoidal height (mm)',4
      write(nmstr(5),'(a28,I3)')'Levelling height (mm)',5
      write(nmstr(6),'(a28,I3)')'Hydrologic site ewh (cm)',6
      GRS(1)= 3.986004415d14; GRS(2)=6378137.d0; GRS(3)=1.0826359d-3
      GRS(4) = 7.292115d-5; GRS(5)=1.d0/298.25641153d0;ae=GRS(2)
      fk=3.d0/5.517d0*GRS(1);wght(1:6)=1.d0
      hepch=nint(inp(1));frow=nint(inp(2));kndrow=nint(inp(3));wghrow=nint(inp(4));itern=nint(inp(6))
      lnth=inp(5)/ae/RAD;knd=nint(inp(7));cntrlknd=nint(inp(8));wght(cntrlknd)=inp(9)**2;dr=inp(10)/ae/RAD
      dr1=inp(11)/ae/RAD;krbf=nint(inp(12));krbf1=nint(inp(13));order=nint(inp(14));order1=nint(inp(15))
      minN=nint(inp(16));minN1=nint(inp(17));maxN=nint(inp(18));maxN1=nint(inp(19));dpth=inp(20);dpth1=inp(21)
      open(unit=8,file=dtmfl,status="old",iostat=status)
      if(status/=0) goto 902
      read(8,'(a)') line
      call PickReclong(line,kln,rec,sn)
      if(sn<6)then
         close(8);goto 902
      endif
      hd(1:6)=rec(1:6)
	nlat=nint((hd(4)-hd(3))/hd(6))
	nlon=nint((hd(2)-hd(1))/hd(5))
	hd(5)=(hd(2)-hd(1))/dble(nlon)
	hd(6)=(hd(4)-hd(3))/dble(nlat)
 	allocate(hgt(nlat,nlon), stat=astat(1))
	if (sum(astat(1:1)) /= 0) then
         close(8);goto 902
      endif
 	do i=1,nlat
	   read(8,*,end=905)(hgt(i,j),j=1,nlon)
      enddo
905   close(8)
      kk=nint((hd(4)-hd(3))/lnth);lnth=(hd(4)-hd(3))/dble(kk)
      lvl=180.d0/lnth
      open(unit=8,file=obstmsqdfl,status="old",iostat=status)
      if(status/=0)goto 901
      read(8,'(a)') headmn
      call PickReclong(headmn,kln,rec,sn)
      mn=0!有效观测数
      if(sn<hepch)then
         close(8);goto 901!缺采样历元时刻信息
      endif
      tmlth=sn-hepch+1!采样时序长度
      do while(.not.eof(8))
        read(8,'(a)') line
        call PickReclong(line,kln,rec,sn)
        if(rec(2)<hd(1).or.rec(2)>hd(2).or.rec(3)<hd(3).or.rec(3)>hd(4))goto 505
        if(sn<kndrow.or.sn<wghrow.or.sn<frow.or.sn<dtrow)goto 505
        if(sn-frow+1<tmlth)goto 505!缺监测量采样时序
        mn=mn+1
 505    continue
	enddo
      close(8)
      if(mn<6)goto 901 !监测量太少!中止计算！
      i=nlat*nlon
 	allocate(obs(mn*3,6,6), stat=astat(1))!球坐标,监测量,权值,监测量;6种监测量(第3维)
 	allocate(obsp(mn*3+i*2,3), stat=astat(2))!全部观测点球坐标
 	allocate(strnm(mn,6), stat=astat(3))
 	allocate(chs(mn*3+i*2,6,5), stat=astat(4))
 	allocate(rst(nlat,nlon,15), stat=astat(5))
 	allocate(tm(tmlth), stat=astat(6))
	if (sum(astat(1:6)) /= 0) goto 901
      obsn=0;k=0!重新读取观测文件，转换为观测点球面坐标，计算平均地心距r0
      open(unit=8,file=obstmsqdfl,status="old",iostat=status)
      read(8,'(a)') headmn
      call PickReclong(headmn,kln,rec,sn)
      tm(1:tmlth)=rec(hepch:sn);BLH(3)=0.d0!地面
      do while(.not.eof(8))
        read(8,'(a)') line
        call PickReclong(line,kln,rec,sn)
        if(rec(2)<hd(1).or.rec(2)>hd(2).or.rec(3)<hd(3).or.rec(3)>hd(4))goto 506
        if(sn<kndrow.or.sn<wghrow.or.sn<frow.or.sn<dtrow)goto 506
        if(sn-frow+1<tmlth.or.rec(dtrow)>9900.d0)goto 506!缺监测量采样时序
        if(rec(kndrow)<0.5d0.or.rec(kndrow)>6.5d0)goto 506!监测量类型超限
        k=k+1;BLH(1)=rec(3);BLH(2)=rec(2)
        nk=nint(rec(kndrow));obsn(nk)=obsn(nk)+1
        call BLH_RLAT(GRS,BLH,rln);obsp(k,1:3)=rln(1:3)
        if(rec(wghrow)<-0.0001)rec(wghrow)=1.d0
        obs(obsn(nk),1:3,nk)=rln(1:3);obs(obsn(nk),4,nk)=rec(dtrow)!4当前监测量-迭代更新
        obs(obsn(nk),5,nk)=rec(wghrow);obs(obsn(nk),6,nk)=rec(dtrow)!6原监测量-输出
        call PickRstrlg(line,kln,str,sn)
        strnm(obsn(nk),nk)=str(1)
 506    continue
	enddo
      close(8);mn=k
      rhd(1:4)=hd(1:4);BLH(2)=(hd(1)+hd(2))/2.d0!!!!!!!目标格网范围用球坐标表示
      BLH(1)=hd(3);call BLH_RLAT(GRS,BLH,rln);rhd(3)=rln(2)
      BLH(1)=hd(4);call BLH_RLAT(GRS,BLH,rln);rhd(4)=rln(2)!!!!!!!!!
      dlat=180.d0/dble(lvl); r0=rln(1)
      nn=nint((rhd(4)-rhd(3))/dlat+0.5);mm=nint((rhd(2)-rhd(1))/dlat+0.5)!mm平行圈方向最大格网数
      allocate(nln(nn),sr(nn),dl(nn),nrd(nn,mm),gpnt(nn,mm),lon(nn,mm),rlatlon(2*(nn+mm),2),enode(2*(nn+mm)))
      call ReuterGrid(rhd,lvl,Kt,blat,nn,mm,nln,sr,dl,nrd,lon)!Kt节点数/未知数个数>>>>lvl1
 	allocate(BPB(Kt,Kt,6), stat=astat(1))!6种监测量
 	allocate(APA(Kt,Kt), stat=astat(2))
 	allocate(B15(Kt,15), stat=astat(3))
	if (sum(astat(1:3)) /= 0) goto 900
      allocate(BPL(Kt*3,6),APL(Kt),BB(Kt),B2(Kt),xx(Kt),node(Kt))
      gpnt=0!计算格网中测点数，修正Reuter格网节点数Kt,序号nrd
      call Edgnode(enode,rlatlon,lvl,edgn,lon,blat,nln,gpnt,nn,mm)
      ks=edgn+mn;BPL(1:edgn,2:3)=rlatlon(1:edgn,1:2);BPL(1:edgn,1)=r0
      BPL(edgn+1:edgn+mn,1:3)=obsp(1:mn,1:3)
      do i=1,nn!补充节点
          tmp=blat+(i-0.5d0)*dlat
          do j=1,nln(i)
            ks=ks+1;BPL(ks,2)=tmp;BPL(ks,3)=lon(i,j)
          enddo
      enddo
      BPL(edgn+mn+1:ks,1)=r0;gpnt=0;nrd=0
      call AdjReuterGrd(BPL(1:ks,1:3),ks,Kt,blat,lvl,nn,mm,nln,dl,lon,nrd,gpnt)
      chs=0.d0!计算残差并统计
      st=0.d0;rst=0.d0;kp=0!!!!!!!!!!!!!!!!!!!!!!
4444  NF=nint(dr*3600)!影响半径等分,间隔1″,NF+1→[0,dr]
      nd=nint(dr/dlat+0.5d0)!dlat格网间隔,积分半径对应的格网数>>>>lvl1
      allocate(RBF(NF+1,6),RBFi(NF+1,11),RBFn(maxN-minN+1,11))
 	allocate(mpn(maxN-minN+1,NF+1), stat=astat(1))
 	allocate(mdp(maxN-minN+1,NF+1), stat=astat(2))
 	allocate(mp2(maxN-minN+1,NF+1), stat=astat(3))
	if (sum(astat(1:3)) /= 0) goto 601
   !计算minN~maxN阶NF+1组勒让德函数[minN~maxN]×[0,dr]
      call LegPn02(mpn,mdp,mp2,minN,maxN,NF,dr)
   !由初始补偿深度dpth,计算SRFB曲线
      br=r0-dpth;nta=br/r0;rlnk(1)=br!初始补偿深度dpth和宽度参数nta
      call SRBF6curve(RBF,flv,order,krbf,mpn,minN,maxN,NF,nta)
      call SRBF11all(RBFi,flv,order,krbf,mpn,mdp,mp2,minN,maxN,NF,nta)
1011  BPB=0.d0;BPL=0.d0!构造观测方程和法方程
      do kobs=1,6
        do k=1,obsn(kobs)!obs(mn,6,kobs)-监测量
          rln(1:3)=obs(k,1:3,kobs);wgh=obs(k,5,kobs);val=obs(k,4,kobs)/unt(kobs);rr=rln(1)
          call RLAT_BLH(GRS,rln,BLH);call GNormalfd(BLH,NFD,GRS);gr=NFD(2)
          if(kobs==1.or.kobs==4.or.kobs==5)GMr=fk/gr/rr
          if(kobs==2.or.kobs==3)GMr=fk/rr/rr
          if(kobs==6)GMr=rr
          call RtGridij(rln,ki,kj,blat,lvl,nn,mm,nln,dl,lon)!计算rln在格网中的序号ki>0,kj>0
          mk=0;BB=0.d0;node=0!BB-RBF系数,node(Kt)观测量作用半径范围内有效节点序号
          do i=ki-nd,ki+nd!增加实现节点序号与格网行列号对应，即节点序号二维整数数组
            if(i<1.or.i>nn)goto 1001!lon(i,j),blat第一行平行圈格网中心地心纬度°
            rlnk(2)=blat+(i-1.d0)*dlat
            do j=kj-nd,kj+nd
              if(j<1.or.j>nln(i)) goto 1002 
              if(nrd(i,j)<1) goto 1002 !nrd(i,j)未知数序号
              rlnk(3)=lon(i,j);call drln(rln,rlnk,dln)!距离m与夹角°
              if(dln(2)>dr)goto 1002
              call RBFvalue(RBF(:,kobs),NF,dr,dln(2),tmp)!由球面角距°dln(2)内插RBF值tmp
              mk=mk+1;node(mk)=nrd(i,j)
              BB(nrd(i,j))=GMr*tmp*dexp((minN-1.0)*dlog(r0/rr))
1002          continue
            enddo
1001        continue    
          enddo
          do i=1,mk
            ki=node(i)
            BPL(ki,kobs)=BPL(ki,kobs)+BB(ki)*val*wgh
            do j=1,i
              kj=node(j)
              BPB(ki,kj,kobs)=BPB(ki,kj,kobs)+BB(ki)*BB(kj)*wgh
	      enddo
          enddo
        enddo
        if(obsn(kobs)>0.and.kp==0)call Stat1d(obs(1:obsn(kobs),6,kobs),obsn(kobs),st(kobs,1:4,kp+1))
      enddo!!!!***
      APA=0.d0;APL=0.d0;ab=0.d0
      do k=1,6!计算法方程对角线标准差
        if(obsn(k)<1)goto 3131
        do i=1,Kt
           B2(i)=BPB(i,i,k)
        enddo
        call Stat1d(B2(1:Kt),Kt,tsr)
        ab(k)=3.d0*tsr(2)+(maxval(B2(1:Kt))-minval(B2(1:Kt)))*0.1d0
3131    continue
      enddo
      do k=1,6
        if(obsn(k)<1)goto 3232
        tmp=wght(k)/ab(k)!法方程对角线最大最小值差法组合，监测量调控
        do i=1,Kt
          APL(i)=APL(i)+tmp*BPL(i,k)
          do j=1,i
            APA(i,j)=APA(i,j)+tmp*BPB(i,j,k)
	    enddo
        enddo
3232    continue
      enddo
1601  tmp=0.d0
	do i=1,Kt
	   do j=1,i-1
	      APA(j,i)=APA(i,j)
         enddo
         tmp=tmp+APA(i,i)**2/dble(Kt)
      enddo
      tmp=dsqrt(tmp)
      !以Reuter格网四周节点未知数为零组成观测方程，抑制边缘效应。
      !节点序号数组enode
      do i=1,edgn!edgn-Reuter格网四周节点数
         ki=enode(i); APA(ki,ki)=APA(ki,ki)+tmp*4.d15!/dsqrt(dble(mn))
      enddo
	do i=1,Kt
         APA(i,i)=APA(i,i)+tmp*2.d-4
      enddo
      xx=0.d0!knd=1LU分解,2Cholesky分解,3最小二乘QR分解,4最小范数奇异值分解
      call Equsolve(APA,xx,Kt,APL,knd,bf)
      do kobs=1,6!!!!!***
        do k=1,obsn(kobs)!6种观测场元-球坐标，观测量，权值obs(k,1:5,kobs)
          rln(1:3)=obs(k,1:3,kobs);rr=rln(1)
          call RLAT_BLH(GRS,rln,BLH);call GNormalfd(BLH,NFD,GRS);gr=NFD(2)
          if(kobs==1.or.kobs==4.or.kobs==5)GMr=fk/gr/rr
          if(kobs==2.or.kobs==3)GMr=fk/rr/rr
          if(kobs==6)GMr=rr
          call RtGridij(rln,ki,kj,blat,lvl,nn,mm,nln,dl,lon)!计算rln在格网中的序号ki>0,kj>0
          mk=0;BB=0.d0;node=0!s-RBF系数,node(Kt)观测量作用半径范围内有效节点序号
          do i=ki-nd,ki+nd!增加实现节点序号与格网行列号对应，即节点序号二维整数数组，0表示
            if(i<1.or.i>nn)goto 2001!lon(i,j),blat第一行平行圈格网中心地心纬度°
            rlnk(2)=blat+(i-1.d0)*dlat
            do j=kj-nd,kj+nd
              if(j<1.or.j>nln(i)) goto 2002 
              if(nrd(i,j)<1) goto 2002 !nrd(i,j)未知数序号
              rlnk(3)=lon(i,j);call drln(rln,rlnk,dln)!距离m与夹角°
              if(dln(2)>dr)goto 2002
              call RBFvalue(RBF(:,kobs),NF,dr,dln(2),tmp)!由球面角距°dln(2)内插RBF值tmp
              mk=mk+1;node(mk)=nrd(i,j)
              BB(nrd(i,j))=GMr*tmp*dexp((minN-1.0)*dlog(r0/rr))
2002          continue
            enddo
2001        continue    
          enddo
          val=0.d0
          do i=1,mk
            ki=node(i)
            val=val+BB(ki)*xx(ki)
          enddo
          chs(k,kobs,kp+1)=obs(k,4,kobs)-val*unt(kobs);obs(k,4,kobs)=chs(k,kobs,kp+1)
        enddo
        if(obsn(kobs)>0)call Stat1d(chs(1:obsn(kobs),kobs,kp+1),obsn(kobs),st(kobs,1:4,kp+2))
      enddo!!!!***
      write(stm,*)nint(tm(dtrow-frow+1))
 	write(rntfl,*) "testdata\rnt",trim(AdjustL(stm)),".txt"
      open(unit=10,file=rntfl,status="replace")
      do kobs=1,6
        if(obsn(kobs)>1)then
           write(10,'(a30)')trim(nmstr(kobs))
           do k=1,kp+2
              write(10,'(I18,80F10.4)')k-1,st(kobs,1:4,k)
           enddo
        endif
      enddo
      do kobs=1,6
        do k=1,obsn(kobs)
          write(str(1),'(a)')strnm(k,kobs);call RLAT_BLH(GRS,obs(k,1:3,kobs),BLH)
	    write(10,'(a12,2F10.4,F8.2,I3,80F12.4)')trim(str(1)),BLH(2),BLH(1),obs(k,5,kobs),kobs,obs(k,6,kobs),chs(k,kobs,1:kp+1)
        enddo
      enddo
      close(10)
	do ni=1,nlat
        BLH(1)=hd(3)+(real(ni)-0.5d0)*hd(6)
        do nj=1,nlon
          BLH(2)=hd(1)+(real(nj)-0.5d0)*hd(5);BLH(3)=hgt(ni,nj)
          call BLH_RLAT(GRS,BLH,rln);rr=rln(1)
          call GNormalfd(BLH,NFD,GRS);gr=NFD(2)
          call RtGridij(rln,ki,kj,blat,lvl,nn,mm,nln,dl,lon)!计算rln在格网中的序号ki>0,kj>0
          mk=0;B15=0.d0;node=0!s-RBF系数,node(Kt)观测量作用半径范围内有效节点序号
          ff(1)=rr;ff(2)=fk/gr/rr;ff(3)=fk/rr/rr;ff(4)=ff(3);ff(5)=ff(2)/rr;ff(6)=ff(5)
          ff(7)=-ff(5);ff(8)=ff(2);ff(9)=ff(2);ff(10)=ff(3)/rr;ff(11)=ff(5)/rr
          do i=ki-nd,ki+nd!增加实现节点序号与格网行列号对应，即节点序号二维整数数组，0表示
            if(i<1.or.i>nn)goto 3001!lon(i,j),blat第一行平行圈格网中心地心纬度°
            rlnk(2)=blat+(i-1.d0)*dlat
            do j=kj-nd,kj+nd
              if(j<1.or.j>nln(i)) goto 3002 
              if(nrd(i,j)<1) goto 3002 !nrd(i,j)未知数序号
              rlnk(3)=lon(i,j);call drln(rln,rlnk,dln)!距离m与夹角°
              if(dln(2)>dr)goto 3002
              mk=mk+1;node(mk)=nrd(i,j);sf=0.d0
              do k=1,11!!!!!!!!!!!!!!!!!!!!!!!!!11
                call RBFvalue(RBFi(:,k),NF,dr,dln(2),tmp)!由球面角距°dln(2)内插RBF值tmp
                sf(k)=ff(k)*tmp*dexp((minN-1.0)*dlog(r0/rr))
              enddo
              tt=dcos(dln(2)*RAD);sinf=dsqrt(1.d0-tt**2);if(sinf<1.d-12)sinf=1.d-12
              fi=rln(2)*RAD;fi1=rlnk(2)*RAD;la=rln(3)*RAD;la1=rln(3)*RAD;dla=(rlnk(3)-rln(3))*RAD
	        cosa=(dcos(fi)*dsin(fi1)-dsin(fi)*dcos(fi1)*dcos(dla))/sinf
	        sina=dcos(fi1)*dsin(dla)/sinf
              dfi2=-dsin(fi)*dsin(fi1)-dcos(fi)*dcos(fi1)*dcos(dla)+tt*cosa**2;dfi2=dfi2/sinf
              dnla2=-dcos(fi1)*dsin(dla)+tt*dcos(fi)*sinf**2;dnla2=dnla2/sinf
              fw(1:5)=sf(1:5);fw(6)=fw(5);fw(7)=sf(6);fw(8)=sf(6);fw(9)=sf(7);fw(10)=sf(7)
              fw(11)=sf(8);fw(12)=sf(9);fw(13)=sf(10);fw(14)=sf(11);fw(15)=-sf(11)/((dcos(fi))**2)
              fw(5)=fw(5)*cosa;fw(6)=fw(6)*sina;fw(7)=fw(7)*cosa;fw(8)=fw(8)*sina
              fw(9)=fw(9)*cosa;fw(10)=fw(10)*sina;fw(14)=fw(14)*dfi2;fw(15)=fw(15)*dnla2
              B15(nrd(i,j),1:15)=fw(1:15)
3002          continue
            enddo
3001        continue    
          enddo
          ff=0.d0
          do i=1,mk
            ki=node(i)
            do k=1,15
              ff(k)=ff(k)+B15(ki,k)*xx(ki)
            enddo
          enddo
          rst(ni,nj,1:15)=rst(ni,nj,1:15)+ff(1:15)*u15(1:15)
        enddo
      enddo
      if(kp>itern-1)goto 2222
      if(kp>0)goto 2221
      call SRBFone(RBFi,RBFn,flv,order,krbf,mpn,mdp,mp2,minN,maxN,NF,nta)
	write(outfl,*) "SRBFspc.txt"
      open(unit=10,file=outfl,status="replace")
	write(10,'(2I3,2I6,F8.2,a)')krbf,order,minN,maxN,dpth*1.d-3,
     * " ewh,hgtanomaly,terrgrav,gravdisturb,tilt,vdeflect,hdisplace,radial,ortheight,gradient,vertgrad" 
      tmp=dr/dble(NF)*RAD*ae*1.d-3
 	do i=1,NF
	   write(10,'(F12.3,11F13.5)')-(NF-i+1.0)*tmp,(RBFi(NF-i+2,j),j=1,11)
      enddo
 	do i=1,NF+1
	   write(10,'(F12.3,11F13.5)')(i-1.0)*tmp,(RBFi(i,j),j=1,11)
      enddo
      close(10)
	write(outfl,*) "SRBFdgr.txt"
      open(unit=10,file=outfl,status="replace")
	write(10,'(2I3,2I6,F8.2,a)')krbf,order,minN,maxN,dpth*1.d-3,
     * " ewh,hgtanomaly,terrgrav,gravdisturb,tilt,vdeflect,hdisplace,radial,ortheight,gradient,vertgrad" 
 	do i=minN,maxN
	   write(10,'(I8,11F13.5)')i,(RBFn(i-minN+1,j),j=1,11)
      enddo
      close(10)
2221  deallocate(mpn,mdp,mp2,RBF,RBFi,RBFn)
      krbf=krbf1;order=order1;minN=minN1;maxN=maxN1;dpth=dpth1;dr=dr1
      kp=kp+1;goto 4444
2222 	continue
      open(unit=10,file=outfl,status="replace")
      write(10,101)(hd(i),i=1,6),nint(tm(dtrow-frow+1))
	do i=1,nlat
	  write(10,'(15F12.4)')(rst(i,j,1),j=1,nlon)!等效水高m->cm
 	enddo
      close(10)
 	write(outfl,*) "testdata\SRBFgeoid",trim(AdjustL(stm)),".dat"
      open(unit=10,file=outfl,status="replace")
      write(10,101)(hd(i),i=1,6),nint(tm(dtrow-frow+1))
	do i=1,nlat
	  write(10,'(15F12.4)')(rst(i,j,2),j=1,nlon)
 	enddo
      close(10)
 	write(outfl,*) "testdata\SRBFterrgrav",trim(AdjustL(stm)),".dat"
      open(unit=10,file=outfl,status="replace")
      write(10,101)(hd(i),i=1,6),nint(tm(dtrow-frow+1))
	do i=1,nlat
	  write(10,'(15F12.4)')(rst(i,j,3),j=1,nlon)
 	enddo
      close(10)
 	write(outfl,*) "testdata\SRBFgravdist",trim(AdjustL(stm)),".dat"
      open(unit=10,file=outfl,status="replace")
      write(10,101)(hd(i),i=1,6),nint(tm(dtrow-frow+1))
	do i=1,nlat
	  write(10,'(15F12.4)')(rst(i,j,4),j=1,nlon)
 	enddo
      close(10)
 	write(outfl,*) "testdata\SRBFgrndtilt",trim(AdjustL(stm)),".dat"
      open(unit=10,file=outfl,status="replace")
      write(10,101)(hd(i),i=1,6),nint(tm(dtrow-frow+1))
	do i=1,nlat
	  write(10,'(15F12.4)')(rst(i,j,5),j=1,nlon)
 	enddo
	do i=1,nlat
	  write(10,'(15F12.4)')(rst(i,j,6),j=1,nlon)
 	enddo
      close(10)
 	write(outfl,*) "testdata\SRBFvertdefl",trim(AdjustL(stm)),".dat"
      open(unit=10,file=outfl,status="replace")
      write(10,101)(hd(i),i=1,6),nint(tm(dtrow-frow+1))
	do i=1,nlat
	  write(10,'(15F12.4)')(rst(i,j,7),j=1,nlon)
 	enddo
	do i=1,nlat
	  write(10,'(15F12.4)')(rst(i,j,8),j=1,nlon)
 	enddo
      close(10)
 	write(outfl,*) "testdata\SRBFhorzdisp",trim(AdjustL(stm)),".dat"
      open(unit=10,file=outfl,status="replace")
      write(10,101)(hd(i),i=1,6),nint(tm(dtrow-frow+1))
	do i=1,nlat
	  write(10,'(15F12.4)')(rst(i,j,9),j=1,nlon)
 	enddo
	do i=1,nlat
	  write(10,'(15F12.4)')(rst(i,j,10),j=1,nlon)
 	enddo
      close(10)
 	write(outfl,*) "testdata\SRBFelliphgt",trim(AdjustL(stm)),".dat"
      open(unit=10,file=outfl,status="replace")
      write(10,101)(hd(i),i=1,6),nint(tm(dtrow-frow+1))
	do i=1,nlat
	  write(10,'(15F12.4)')(rst(i,j,11),j=1,nlon)
 	enddo
      close(10)
 	write(outfl,*) "testdata\SRBForthohgt",trim(AdjustL(stm)),".dat"
      open(unit=10,file=outfl,status="replace")
      write(10,101)(hd(i),i=1,6),nint(tm(dtrow-frow+1))
	do i=1,nlat
	  write(10,'(15F12.4)')(rst(i,j,12),j=1,nlon)
 	enddo
      close(10)
 	write(outfl,*) "testdata\SRBFgradient",trim(AdjustL(stm)),".dat"
      open(unit=10,file=outfl,status="replace")
      write(10,101)(hd(i),i=1,6),nint(tm(dtrow-frow+1))
	do i=1,nlat
	  write(10,'(15ES14.5)')(rst(i,j,13),j=1,nlon)
 	enddo
      close(10)
 	write(outfl,*) "testdata\SRBFhorzgrad",trim(AdjustL(stm)),".dat"
      open(unit=10,file=outfl,status="replace")
      write(10,101)(hd(i),i=1,6),nint(tm(dtrow-frow+1))
	do i=1,nlat
	  write(10,'(15ES14.5)')(rst(i,j,14),j=1,nlon)
 	enddo
	do i=1,nlat
	  write(10,'(15ES14.5)')(rst(i,j,15),j=1,nlon)
 	enddo
      close(10)
1111  close(12)
      deallocate(mpn,mdp,mp2,RBF,RBFi,RBFn)
601 	deallocate(nln,sr,dl,nrd,gpnt,lon,rlatlon,enode,BPL,APL,BB,B2,xx,node,BPB,APA,B15) !调入迭代内部
900   deallocate(obs,obsp,chs,rst,strnm,tm)
901 	deallocate(hgt)
902   continue
101   format(4F12.6,2ES16.8,I14)
      end
