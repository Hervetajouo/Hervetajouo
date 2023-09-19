        SUBROUTINE n2oco(R,tht1,tht2,phee,V)
        implicit real*8 (a-h,o-z)
        parameter (ndmn=181,ndr=71,ndp=361)
        dimension xr(ndr),xt1(ndmn),xt2(ndmn),xp(ndmn)
        dimension yyR(29,13,13,13)
        dimension valu(ndmn),val(ndr)
        dimension xs(ndmn),xtp(ndr)
        dimension r(ndr), theta1(ndmn), theta2(ndmn), phi(ndp)
        logical   firstcalc
        save      xr, xt1, xt2, xp, yyR

        wn2ha=219474.63d0
        pii=dacos(-1d0)
        firstcalc=.true.
        !theta1=tht1
        !theta2=tht2
        !phi=phee
        theta1=tht1*180d0/pii
        theta2=tht2*180d0/pii
        phi=phee*180d0/pii

        if (firstcalc) then
        open (unit=8,status='old',file='total-N2O.dat')
        rewind(8)
        nr=29 !251
        nt1=13 
        nt2=13 
        np=13
        ntintm=1  !!181  !251
        nrintm=1  !!251
        npintm=1  !!91
        PI=acos(-1.0)
        do ip=1,np
        do it2=1,nt2
        do it1=1,nt1
          do ir=nr,1,-1
           read (8,*) xr(ir),xt1(it1),xt2(it2),xp(ip),yyR(ir,it1,it2,ip)
          enddo
        enddo
        enddo
        enddo
        firstcalc = .false.
        endif
c            write(*,*) '1  --------------->  (R,rheta1)'
c            write(*,*) '2  --------------->  (R,rheta2)'
c            write(*,*) '3  --------------->  (R,phi)'
c            write(*,*) '4  --------------->  (theta1,rheta2)'
c            write(*,*) '5  --------------->  (theta1,phi)'
c            write(*,*) '6  --------------->  (theta2,phi)'
c            write(*,*) '7  --------------->  1 value, 1 point'
c            read(*,*) iflag
c        select case(iflag)
c                    write(*,*) 'r, theta1, theta2, phi  ?'
c                    read(*,*)  r(1),theta1(1),theta2(1),phi(1)
                    nint1=1
                    nint2=1
                    nintp=1
                    nrint=1
c                    do i=1,nrint
c                       r(i)=xr(1)+dfloat(i-1)*0.1d0
c                    enddo
c               case default
c                   write(*,*) ' Revoir votre choix (entre 1 et 6)'
c                   stop
c        end select
        do ip=1,np
        do it2=1,nt2
        do it1=1,nt1
          do ir=5,nr-1
             yyR(ir,it1,it2,ip)=yyR(ir+1,it1,it2,ip)
          enddo
        enddo
        enddo
        enddo
          do ir=5,nr-1
             xr(ir)=xr(ir+1)
          enddo
        nr=nr-1
        do ir=1,nr
        do it1=1,nt1
        do it2=1,nt2
           do ip=1,np
              xs(ip)=yyR(ir,it1,it2,ip)
           enddo
           yp1=0.d0
           ypn=0.d0
           call spline(xp,xs,np,yp1,ypn,valu)
           do it=1,nintp
              xtp(it)=phi(it)
              xeval=xtp(it)
              ich=0
              call splint(xp,xs,valu,np,xeval,a,ich)
              val(it)=a
           enddo
           do it=1,nintp
              yyR(ir,it1,it2,it)=val(it)
           enddo
        enddo
        enddo
        enddo
           do it=1,nintp
              xp(it)=xtp(it)
           enddo
        do ir=1,nr
        do it1=1,nt1
        do ip=1,nintp
           do it2=1,nt2
              xs(it2)=yyR(ir,it1,it2,ip)
           enddo
           yp1=0.d0
           ypn=0.d0
           call spline(xt2,xs,nt2,yp1,ypn,valu)
           do it=1,nint2
              xtp(it)=theta2(it)
              xeval=xtp(it)
              ich=0
              call splint(xt2,xs,valu,nt2,xeval,a,ich)
              val(it)=a
! 100          format(f7.2,' t1',f8.2,' t2',f8.2,' p',f8.2,f20.8)
 100          format(f7.2,f8.2,f8.2,f8.2,f20.9)
           enddo
           do it=1,nint2
              yyR(ir,it1,it,ip)=val(it)
           enddo
        enddo
        enddo
        enddo
           do it=1,nint2
              xt2(it)=xtp(it)
           enddo
        do ir=1,nr
        do it2=1,nint2
        do ip=1,nintp
           do it1=1,nt1
              xs(it1)=yyR(ir,it1,it2,ip)
           enddo
           yp1=0.d0
           ypn=0.d0
           call spline(xt1,xs,nt1,yp1,ypn,valu)
           do it=1,nint1
              xtp(it)=theta1(it)
              xeval=xtp(it)
              ich=0
              call splint(xt1,xs,valu,nt1,xeval,a,ich)
              val(it)=a
           enddo
           do it1=1,nint1
              yyR(ir,it1,it2,ip)=val(it1)
           enddo
        enddo
        enddo
        enddo
           do it=1,nint1
              xt1(it)=xtp(it)
           enddo
        do it1=1,nint1
        do it2=1,nint2
        do ip=1,nintp
           do ir=1,nr
              xs(ir)=yyR(ir,it1,it2,ip)
           enddo
           yp1=(xs(2)-xs(1))/(xr(2)-xr(1))
     & +((xs(2)-xs(1))/(xr(2)-xr(1))-(xs(3)-xs(2))/(xr(3)-xr(2)))*1.d0
           ypn=0.d0
           call spline(xr,xs,nr,yp1,ypn,valu)
           do it=1,nrint
              xtp(it)=r(it)
              xeval=xtp(it)
              ich=0
              call splint(xr,xs,valu,nr,xeval,a,ich)
              V=a
c              write(12,100) xtp(it),xt1(it1),xt2(it2),xp(ip),a
c              write(6,100) xtp(it),xt1(it1),xt2(it2),xp(ip),a
           enddo
c           do it=1,nrint
c              yyR(it,it1,it2,ip)=val(it)
c           enddo
        enddo
        enddo
        enddo
c           do it=1,nrint
c              xr(it)=xtp(it)
c           enddo

        V=a/wn2ha

        END SUBROUTINE n2oco
c..............................................................................
       SUBROUTINE SPLINE(X,Y,N,YP1,YPN,Y2)
       IMPLICIT DOUBLE PRECISION (A-H,O-Z)
       PARAMETER(NMAX=2000)
c       DIMENSION X(N),Y(N),Y2(N),U(NMAX)
       DIMENSION X(N),Y(N),Y2(N),U(N)
       IF(YP1.GT..99D+30)THEN
         Y2(1)=0.0D0
         U(1)=0.0D0
       ELSE
         Y2(1)=-0.5D0
         U(1)=(3.D0/(X(2)-X(1)))*((Y(2)-Y(1))/(X(2)-X(1))-YP1)
       ENDIF
       DO 11 I=2,N-1
        SIG = (X(I)-X(I-1))/(X(I+1)-X(I-1))
        P = SIG*Y2(I-1)+2.0D0
        Y2(I)=(SIG-1.0D0)/P
c               write(*,*) x(i)
        U(I) = (6.0D0*((Y(I+1)-Y(I))/(X(I+1)-X(I))-(Y(I)-Y(I-1))
     &   /(X(I)-X(I-1)))/(X(I+1)-X(I-1)) - SIG*U(I-1))/P
 11     CONTINUE
        IF(YPN.GT..99D30)THEN
         QN=0.0D0
         UN=0.0D0
        ELSE
         QN=0.50D0
         UN=(3.0D0/(X(N)-X(N-1)))*(YPN-(Y(N)-Y(N-1))/(X(N)-X(N-1)))
        ENDIF
         Y2(N)=(UN-QN*U(N-1))/(QN*Y2(N-1)+1.)
        DO 12 K=N-1,1,-1
         Y2(K)=Y2(K)*Y2(K+1)+U(K)
 12     CONTINUE
        RETURN
        END

****************************************************************************************************
        SUBROUTINE SPLINT(XA,YA,Y2A,N,X,Y,ICH)
C Given arrays XA and YA of lengh N,which tabulate a function, and
C given the array Y2A, which is the output from SPLINE above, and given
C a value X, this routine  returns
C       for ICH=0   a cubic-spline interpolated value Y.
C       for ICH=1   the first derivative value Y.
C       for ICH=2   the second derivative value Y.
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        DIMENSION XA(N),YA(N),Y2A(N)
        KLO=1
        KHI=N
 1      IF ((KHI-KLO).GT.1)THEN
            K=(KHI+KLO)/2
            IF(XA(K).GT.X)THEN
                 KHI=K
            ELSE
                 KLO=K
            ENDIF
         GO TO 1
         ENDIF
         H = XA(KHI) - XA(KLO)
         IF(H.EQ.0) PAUSE'BAD XA INPUT'
         A = (XA(KHI)-X)/H
         B = (X-XA(KLO))/H
           IF(ICH.EQ.0)THEN
         Y = A*YA(KLO) + B*YA(KHI) +
     &          ((A**3-A)*Y2A(KLO) + (B**3-B)*Y2A(KHI))*(H**2)/6.0D0
           ELSEIF(ICH.EQ.1)THEN
         Y = (YA(KHI)-YA(KLO))/H
     &          -(3.0D0*(A**2)-1.0D0)*Y2A(KLO)*H/6.0D0
     &          +(3.0D0*(B**2)-1.0D0)*Y2A(KHI)*H/6.0D0
           ELSEIF(ICH.EQ.2)THEN
         Y = A*Y2A(KLO) + B*Y2A(KHI)
           ELSE
         Y = A*YA(KLO) + B*YA(KHI)
           ENDIF
         RETURN
         END
***************************************************************************************************************************
      SUBROUTINE SORT1(A,JMAT,NA,EKK,SIG)
C  SORTS THE FIRST NA ELEMENTS OF ARRAY A INTO ASCENDING ORDER.
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION A(1),JMAT(1)
      dimension sig(1),ekk(1)
      NAM1=NA-1
      DO 30 K=1,NAM1
      G=A(K)

       ekg = ekk(k)
       sigg = sig(k)

      JMM=JMAT(K)
      J=K
      KP1=K+1
         DO 10 I=KP1,NA
         IF(A(I).GE.G) GO TO 10
         G=A(I)

           ekg = ekk(i)
           sigg = sig(i)

         JMM=JMAT(I)
         J=I
10    CONTINUE
C  THE J-TH ELEMENT IS NOW THE LARGEST ELEMENT REMAINING TO BE SORTED
      IF(J.EQ.K) GO TO 30
      A(J)=A(K)
          ekk(j) = ekk(k)
          sig(j) = sig(k)
      JMAT(J)=JMAT(K)
      A(K)=G
          ekk(k) = ekg
          sig(k) = sigg
      JMAT(K)=JMM
30    CONTINUE
      RETURN
      END

                                                   
