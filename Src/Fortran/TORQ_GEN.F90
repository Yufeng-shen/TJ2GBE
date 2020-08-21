!***********************************************************************
  PROGRAM torque_gen
!***********************************************************************
! Generates triple junctions - data for reproduction 
! of the energy function.
! The torque term is taken into account.
! Because MINUIT crashes, the program is run from a batch file; after
! each crush the program is restarted automatically.
! The generated data are saved in 'my_file_name.dat', where 
! my_file_name is a random number.
!***********************************************************************
  IMPLICIT REAL*8 (A-H,O-Z)
  COMMON/globals/ pi, eps
  COMMON/trijun/g1(3,3),g2(3,3),g3(3,3),gx1(3,3),gx2(3,3),gx3(3,3),cn1(3),cn2(3),cn3(3),sn1(3),sn2(3),sn3(3),sl(3)
!***********************************************************************
!  'g_i' grain orientations.
!  'gx_i' grain misorientations.
!  'cn_i' normals to the boundaries (crystal coordinate system)
!  'sn_i' normals to the boundaries (sample coordinate system)
!  'sl' direction of the triple line (sample coordinate system)
!***********************************************************************
  REAL*8, DIMENSION(3) :: eulang
  CHARACTER(LEN=8) my_file_name,arg
  REAL*8 xi(3)
  CALL GET_COMMAND_ARGUMENT(1,arg)
  my_file_name=arg
  CALL GET_COMMAND_ARGUMENT(2,arg)
  READ(arg,'(i8)')npoint
  WRITE(*,*)my_file_name
  WRITE(*,*)npoint
  CALL initox(my_file_name)
!  WRITE(*,'(2x,''Number of points : '')',ADVANCE='NO')
!  READ(*,*) npoint
!  npoint=10000
  IF(npoint==0) GOTO 3
  OPEN(1,file=my_file_name//'.dat')
  OPEN(6)
!  OPEN(6,CARRIAGECONTROL='FORTRAN')
  snoraliz=0.
  WRITE(*,'(30(/),32X,''A.Morawiec, modified by Y.Shen'',5(/))')
 
  OPEN(5,FILE='groundTruth.txt',STATUS='NEW') 
!------------------------------------------------------------
  mainloop: DO i=1,npoint
    WRITE(1,'(1x,i7)') i
!    IF(MOD(i,10).EQ.0) WRITE(6,101) i,npoint

!   -------- The main subroutine -------- 
    CALL inclinations
    CALL force_vector(gx1,cn1,xi)
    WRITE(5,*) xi
    CALL force_vector(gx2,cn2,xi)
    WRITE(5,*) xi
    CALL force_vector(gx3,cn3,xi)
    WRITE(5,*) xi

!   This is auxiliary part for the determination of the normalization coefficient
!   because the original function is eq. to 1 + some cusps, i.e., the integral
!   over the whole space is less than 1.
    sumox=energy_gn(gx1,cn1)+energy_gn(gx2,cn2)+energy_gn(gx3,cn3)
    snoraliz=snoraliz+sumox
    WRITE(6,101) i,npoint
    WRITE(6,*) sumox

!   Saving results in the form of: sl vector, euler angles and sn_i vectors.
    WRITE(1,'(8x,3f10.6)')    sl(:)

    CALL MEul(g1,eulang(1),eulang(2),eulang(3))
    WRITE(1,'(8x,3f10.4,3x,3f12.6)')    (180./pi)*eulang(:)
    WRITE(1,'(8x,3f10.6)')    sn2(:)

    CALL MEul(g2,eulang(1),eulang(2),eulang(3))
    WRITE(1,'(8x,3f10.4,3x,3f12.6)')    (180./pi)*eulang(:)
    WRITE(1,'(8x,3f10.6)')    sn3(:)

    CALL MEul(g3,eulang(1),eulang(2),eulang(3))
    WRITE(1,'(8x,3f10.4,3x,3f12.6)')    (180./pi)*eulang(:)
    WRITE(1,'(8x,3f10.6)')    sn1(:)

  ENDDO mainloop
  CLOSE(1)
  CLOSE(5)
!------------------------------------------------------------

! This is auxiliary part. 
! The normalization coefficient is saved in 'normailz.res'.
  OPEN(1,file='normailz.res')
2 WRITE(1,*) 'Normailization coefficient : ',1./(snoraliz/(3.*npoint))
  CLOSE(1)

! This is auxiliary part. Creates files for drawing model function.
3 IF(npoint==0) THEN 
    CALL tescik03()
    CALL tescik07()
    CALL tescik09()
    CALL tescik11()
  ENDIF

  STOP
101 FORMAT('+',24x,'Point nr : ',i7,' -->',i7)
  END


!***********************************************************************
  SUBROUTINE inclinations
!***********************************************************************
! Determines geometry of triple junctions.
! Third part of Appendix A.
!***********************************************************************
  IMPLICIT REAL*8 (A-H,O-Z)
  COMMON/trijun/g1(3,3),g2(3,3),g3(3,3),gx1(3,3),gx2(3,3),gx3(3,3),cn1(3),cn2(3),cn3(3),sn1(3),sn2(3),sn3(3),sl(3)

2 licznik=0

! Generation of random orientations.
  CALL gener_g(g1)
  CALL gener_g(g2)
  CALL gener_g(g3)

! Calculation of misorientations.
  gx1=MATMUL(g1,TRANSPOSE(g2))
  gx2=MATMUL(g2,TRANSPOSE(g3))
  gx3=MATMUL(g3,TRANSPOSE(g1))

! If MINUIT is unable to determine inclinations satisfying the Herring relationship
! for given orientations in 5 attempts, new orientations are generated.
! Liczmik counts attempts
1 licznik=licznik+1
  IF(licznik > 5) GOTO 2

! Generation of triple line direction and normal to the first boundary. 
  CALL triple_line(sl,sn1)

! MINUIT determines normals to the other two boundaries.
! 'om2' - angle between 'sn1' and 'sn2' (clockwise).
! 'om3' - angle between 'sn1' and 'sn3' (clockwise).
  CALL start_minuit(om2,om3,iflag)
  IF(iflag.ne.0) GOTO 1

! For om2,om3 leading to deviation=0,
! the normals 'sn2' and 'sn3' to other two boundaries are calculated
  CALL triple_junction(om2,om3)
  CALL deviation(deviat)
  IF(deviat > 0.00001) GOTO 1

! Calculation of boundary normals in crystal coordinate system
  cn1=MATMUL(g1,sn1)
  cn2=MATMUL(g2,sn2)
  cn3=MATMUL(g3,sn3)
  RETURN
  END


!***********************************************************************
      SUBROUTINE start_minuit(om2,om3,iflag)
!***********************************************************************
!     Prepares parameters for MINUIT and calls MINUIT.
!***********************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
      EXTERNAL FCN
      DIMENSION argts(10)
      CHARACTER(LEN=10) pnam
      OPEN(33,STATUS='SCRATCH')
      OPEN(34,STATUS='SCRATCH')
      CALL MNINIT(5,33,34)
      pnam='parameter '
      pi=ACOS(-1.)
      step=pi/30.
      i=1
      start=(2.*pi)/3.
      bndl=pi/2.
      bndu=pi
      CALL MNPARM(i,pnam,start,step,bndl,bndu,iflag)
      IF(iflag.ne.0) THEN
          CLOSE(33)
          CLOSE(34)
          RETURN
      ENDIF
      i=2
      start=(4.*pi)/3.
      bndl=pi
      bndu=(3.*pi)/2.
      CALL MNPARM(i,pnam,start,step,bndl,bndu,iflag)
      IF(iflag.ne.0) THEN
          CLOSE(33)
          CLOSE(34)
          RETURN
      ENDIF
      zero=0.
      argts(1)=10000.
      izero=0
      ione=1
      CALL MNEXCM(FCN,'MIGRAD',argts,ione,iflag,zero)
!      CALL MNEXCM(FCN,'SIMPLEX',argts,izero,iflag,zero)
      IF(iflag.ne.0) THEN
          CLOSE(33)
          CLOSE(34)
          RETURN
      ENDIF
      i=1
      CALL MNPOUT(i,pnam,om2,error,bndl,bndu,ivabl)
      i=2
      CALL MNPOUT(i,pnam,om3,error,bndl,bndu,ivabl)
      CLOSE(33)
      CLOSE(34)
      RETURN
      END


!***********************************************************************
      SUBROUTINE FCN(NPAR,GRAD,FVAL,XVAL,IFLAG)
!***********************************************************************
!     Called by MINUIT.
!     FVAL - the value minimized by MINUIT.
!***********************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION XVAL(*), GRAD(*)
      CALL triple_junction(XVAL(1),XVAL(2))
      CALL deviation(FVAL)
      RETURN
      END

!***********************************************************************
  SUBROUTINE deviation(deviat)
!***********************************************************************
!   'deviat' is a measure of the deviation of the square of the 
!   Herring relationship's left side from zero.
!  'sn_i' normals to the boundaries (sample coordinate system)
!  'sl' direction of the triple line (sample coordinate system)
!  'g_i' grain orientations.
!  'gx_i' grain misorientations.
!  'cn_i' normals to the boundaries (crystal coordinate system)
!  'cxi_i' force vectors (crystal coordinate system)
!  'sxi_i' force vectors (sample coordinate system)
!***********************************************************************
  IMPLICIT REAL*8 (A-H,O-Z)
  COMMON/trijun/g1(3,3),g2(3,3),g3(3,3),gx1(3,3),gx2(3,3),gx3(3,3),cn1(3),cn2(3),cn3(3),sn1(3),sn2(3),sn3(3),sl(3)
  REAL*8, DIMENSION(3) :: cxi1, cxi2, cxi3, sxi1, sxi2, sxi3
  REAL*8, DIMENSION(3) :: sumsxi
  REAL*8  :: deviat

! normals in crystal coordinate system
  cn1=MATMUL(g1,sn1)
  cn2=MATMUL(g2,sn2)
  cn3=MATMUL(g3,sn3)

! force vectors (xi_i) in crystal coordinate system
  CALL force_vector(gx1,cn1, cxi1)
  CALL force_vector(gx2,cn2, cxi2)
  CALL force_vector(gx3,cn3, cxi3)

! force vectors (xi_i) in sample coordinate system
  sxi1=MATMUL(TRANSPOSE(g1),cxi1)
  sxi2=MATMUL(TRANSPOSE(g2),cxi2)
  sxi3=MATMUL(TRANSPOSE(g3),cxi3)

  sumsxi=sxi1+sxi2+sxi3

  deviat=DOT_PRODUCT(sumsxi,sumsxi)-(DOT_PRODUCT(sl,sumsxi))**2

  RETURN
  END

!***********************************************************************
  SUBROUTINE tescik03()
!***********************************************************************
! Creates file for drawing model function for the misorientation of Sigma 3.
!***********************************************************************
  IMPLICIT REAL*8 (A-H,O-Z)
  COMMON/globals/ pi, eps
  REAL*8 g(3,3),ax(3), xn(3), xi(3)
  REAL*8 eng(0:360)
  OPEN(2,file='mod_gp03.m')
! Orthogonal matrix (g) from axis (ax) angle (om) parameters.
  ax=SQRT(3.)/3.;
  om=60.*(pi/180.)
  CALL axiort(om,ax,g)
! The loops scan inclinations. 
  DO i=0,90
    DO j=0,360
      write(6,101) j
      a=i*(pi/180.)
      b=j*(pi/180.)
      xn(1)=SIN(a)*COS(b)
      xn(2)=SIN(a)*SIN(b)
      xn(3)=COS(a)
!     Capillarity vector 'xi' for misorientation 'g' and direction 'xn' is calculated. 
      CALL force_vector(g,xn, xi)
!     Energy for  direction 'xn' and capillarity vector 'xi' is calculated. 
      eng(j)=DOT_PRODUCT(xn,xi)
    ENDDO
    WRITE(2,'(2x,361f10.7)') eng
  ENDDO
  CLOSE(2)
  RETURN
101 FORMAT('+',24x,'      03 : ',i7,' --> 360')
  END 

!***********************************************************************
  SUBROUTINE tescik07()
!***********************************************************************
! Creates file for drawing model function for the misorientation of Sigma 7.
!***********************************************************************
  IMPLICIT REAL*8 (A-H,O-Z)
  COMMON/globals/ pi, eps
  REAL*8 g(3,3),ax(3), xn(3), xi(3)
  REAL*8 eng(0:360)
  OPEN(2,file='mod_gp07.m')

  ax=SQRT(3.)/3.;
  om=38.21*(pi/180.)
  CALL axiort(om,ax,g)
  DO i=0,90
    DO j=0,360
      write(6,101) j
      a=i*(pi/180.)
      b=j*(pi/180.)
      xn(1)=SIN(a)*COS(b)
      xn(2)=SIN(a)*SIN(b)
      xn(3)=COS(a)
      CALL force_vector(g,xn, xi)
      eng(j)=DOT_PRODUCT(xn,xi)
    ENDDO
    WRITE(2,'(2x,361f10.7)') eng
  ENDDO
  CLOSE(2)
  RETURN
101 FORMAT('+',24x,'      07 : ',i7,' --> 360')
  END 

!***********************************************************************
  SUBROUTINE tescik09()
!***********************************************************************
! Creates file for drawing model function for the misorientation of Sigma 9.
!***********************************************************************
  IMPLICIT REAL*8 (A-H,O-Z)
  COMMON/globals/ pi, eps
  REAL*8 g(3,3),ax(3), xn(3), xi(3)
  REAL*8 eng(0:360)
  OPEN(2,file='mod_gp09.m')

  ax(1:2)=SQRT(2.)/2.; ax(3)=0.
  om=38.94*(pi/180.)
  CALL axiort(om,ax,g)
  DO i=0,90
    DO j=0,360
      write(6,101) j
      a=i*(pi/180.)
      b=j*(pi/180.)
      xn(1)=SIN(a)*COS(b)
      xn(2)=SIN(a)*SIN(b)
      xn(3)=COS(a)
      CALL force_vector(g,xn, xi)
      eng(j)=DOT_PRODUCT(xn,xi)
    ENDDO
    WRITE(2,'(2x,361f10.7)') eng
  ENDDO
  CLOSE(2)
  RETURN
101 FORMAT('+',24x,'      09 : ',i7,' --> 360')
  END 

!***********************************************************************
  SUBROUTINE tescik11()
!***********************************************************************
! Creates file for drawing model function for the misorientation of Sigma 11.
!***********************************************************************
  IMPLICIT REAL*8 (A-H,O-Z)
  COMMON/globals/ pi, eps
  REAL*8 g(3,3),ax(3), xn(3), xi(3)
  REAL*8 eng(0:360)
  OPEN(2,file='mod_gp11.m')
  write(*,*)'11'

  ax(1:2)=SQRT(2.)/2.; ax(3)=0.
  om=50.48*(pi/180.)
  CALL axiort(om,ax,g)
  DO i=0,90
    DO j=0,360
      write(6,101) j
      a=i*(pi/180.)
      b=j*(pi/180.)
      xn(1)=SIN(a)*COS(b)
      xn(2)=SIN(a)*SIN(b)
      xn(3)=COS(a)
      CALL force_vector(g,xn, xi)
      eng(j)=DOT_PRODUCT(xn,xi)
    ENDDO
    WRITE(2,'(2x,361f10.7)') eng
  ENDDO
  CLOSE(2)
  RETURN
101 FORMAT('+',24x,'      11 : ',i7,' --> 360')
  END 


!***********************************************************************
  SUBROUTINE initox(my_file_name)
!***********************************************************************
! Initialization of constants and cusps.
! The location of a given cusp in cusp() is specified by the boundary matrix
!***********************************************************************
  IMPLICIT REAL*8 (A-H,O-Z)
  CHARACTER(LEN=8) my_file_name
  COMMON/xcusp1/number_of_cusps
  COMMON/xcusp2/depth_of_cusp(50),width_of_cusp(50),cusp(50,4,4)
  COMMON/globals/ pi, eps

  DIMENSION ax(3),gx(3,3),xn(3),b(4,4)
  INTEGER, ALLOCATABLE :: iseed(:)
  INTEGER isize
  pi=ACOS(-1.)
  eps=EPSILON(eps)
  CALL symbs_gb
  CALL RANDOM_SEED(SIZE=isize)
  ALLOCATE(iseed(isize))
  DO i=1,isize
     CALL SYSTEM_CLOCK(icount_num,icount_rate,icount_max)
     iseed(i)=icount_num
  ENDDO
!  CALL my_get_seed(iseed,isize)
!  WRITE(*,'(2x,''put filename (8 characters): '')',ADVANCE='NO')
!  READ(*,*) my_file_name
  CALL RANDOM_SEED(PUT=iseed)
!  CALL RANDOM_SEED()
!  WRITE(my_file_name,'(i8)') iseed(1)
  OPEN(9,file='SIGMA_X.DAT')
  READ(9,*) number_of_cusps
  DO i=1,number_of_cusps
     READ(9,*) sig, idummy, om, ax, xn
!    sig - sigma, om - angle, ax - axis, xn - normal to the densest plane
     om=(pi/180.)*om
     aa=ax(1)**2+ax(2)**2+ax(3)**2
     ax=ax/SQRT(aa)
     xx=xn(1)**2+xn(2)**2+xn(3)**2
     xn=xn/SQRT(xx)
!    Misorientation matrix (gx) for i-th cusp.
     CALL axiort(om,ax,gx)
!    Boundary matrix (b) for i-th cusp.
     CALL bconvert(gx,xn,b)
     cusp(i,:,:)=b
     depth_of_cusp(i)=1./SQRT(sig)
     width_of_cusp(i)=pi/(12.*SQRT(sig))
!     width_of_cusp(i)=pi/(4.*SQRT(sig))
  ENDDO
  CLOSE(9)
  RETURN
  END

!!***********************************************************************
!  SUBROUTINE my_get_seed(iseed,isize)
!!***********************************************************************
!! Seed for random number generator
!!***********************************************************************
!  INTEGER :: count_num,count_rate,count_max
!  INTEGER, ALLOCATABLE :: iseed(:)
!  INTEGER isize
!  IF(ALLOCATED(iseed))then
!      WRITE(*,*) SIZE(iseed,1)
!      DO i=1,isize
!         CALL SYSTEM_CLOCK(count_num,count_rate,count_max)
!         iseed(i)=count_num
!      ENDDO
!  ENDIF
!  END


!***********************************************************************
  FUNCTION energy_gn(g,xn)
!***********************************************************************
! Determines the value of the energy function
! for misorientation 'g' and inclination 'xn'
! First part of Appendix A.
!***********************************************************************
  IMPLICIT REAL*8 (A-H,O-Z)
  COMMON/xcusp1/number_of_cusps
  COMMON/xcusp2/depth_of_cusp(50),width_of_cusp(50),cusp(50,4,4)
  DIMENSION g(3,3), gx(3,3), xn(3), xxn(3), b(4,4), bcusp(4,4)

! The cusp at misorientation 0.
  a=depth_of_cusp(1)
  gx=cusp(1,1:3,1:3)
  x=xmdist(g,gx)/width_of_cusp(1)
  energy_gn=auxilx(x,a)
  CALL bconvert(g,xn,b)

! Remaining cusps.
  DO i=2,number_of_cusps
     a=depth_of_cusp(i)
     bcusp=cusp(i,:,:)
     x=x4dist(b,bcusp)/width_of_cusp(i)
     energy_gn=energy_gn*auxilx(x,a)
  ENDDO
  END FUNCTION energy_gn

!***********************************************************************
  FUNCTION auxilx(x,a)
!***********************************************************************
! Auxiliary function determining the shape of the energy cusp.
! Based on Read-Shockley formula.
!***********************************************************************
  IMPLICIT REAL*8 (A-H,O-Z)
  x=MAX(x,1.0e-10)
  SELECT CASE (FLOOR(x))
  CASE(0)
     auxilx=a*x*(1.-LOG(x))+(1.-a)
  CASE DEFAULT
     auxilx=1.
  END SELECT
  END FUNCTION auxilx

!***********************************************************************
  SUBROUTINE bconvert(g,xn,b)
!***********************************************************************
! Boundary matrix b from misorientation matrix g and inclination xn.
!***********************************************************************
  IMPLICIT REAL*8 (A-H,O-Z)
  DIMENSION g(3,3),xn(3),b(4,4), xn2(3)
  b(1:3,1:3)=g
  b(1:3,4)=xn
  xn2=-MATMUL(TRANSPOSE(g),xn)
  b(4,1:3)=xn2
  b(4,4)=0.
  RETURN
  END

!***********************************************************************
  SUBROUTINE gener_g(g)
!***********************************************************************
! Random orientation.
!***********************************************************************
  IMPLICIT REAL*8 (A-H,O-Z)
  COMMON/globals/ pi, eps
  DIMENSION xx(3),g(3,3)
  CALL RANDOM_NUMBER(xx)
  f1=2*pi*xx(1)
  f=acoss(2.*xx(2)-1.)
  f2=2*pi*xx(3)
  CALL EMat(f1,f,f2,g)
  RETURN
  END

!***********************************************************************
  SUBROUTINE gener_n(xn)
!***********************************************************************
! Random direction.
!***********************************************************************
  IMPLICIT REAL*8 (A-H,O-Z)
  COMMON/globals/ pi, eps
  DIMENSION xx(2),xn(3)
  CALL RANDOM_NUMBER(xx)
  p=2*pi*xx(1)
  t=acoss(2.*xx(2)-1.)
  st=SIN(t)
  xn(1)=st*COS(p)
  xn(2)=st*SIN(p)
  xn(3)=COS(t)
  RETURN
  END

      FUNCTION xmdist(g1,g2)
!***********************************************************************
!     Distance between two misorientations 'g1' and 'g2'.
!     The table gcsym(,,) contains matrices of 'nse' crystal symmetry
!     elements. 
!***********************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON/symmop/nse, nselau, gcsym(24,3,3), bsym(24,4,4)

      DIMENSION g1(3,3),g2(3,3)
      xmdist=-1.
      DO 1 ii=1,nselau
      DO 1 jj=1,nselau
         omeg=0.
         DO 2 i=1,3
         DO 2 j=1,3
         DO 2 k=1,3
         DO 2 l=1,3
  2         omeg=omeg+gcsym(ii,i,j)*g1(i,k)*gcsym(jj,k,l)*g2(j,l)
         IF(omeg.gt.xmdist) xmdist=omeg 
         omeg=0.
         DO 3 i=1,3
         DO 3 j=1,3
         DO 3 k=1,3
         DO 3 l=1,3
  3         omeg=omeg+gcsym(ii,i,j)*g1(i,k)*gcsym(jj,k,l)*g2(l,j)
         IF(omeg.gt.xmdist) xmdist=omeg 
  1   CONTINUE
!      xmdist=acoss((xmdist-1.)/2.)
      xmdist=3.-xmdist
      END

      FUNCTION x4dist(b1,b2)
!***********************************************************************
!     Distance between two 'boundaries' 'b1' and 'b2'.
!     The table bsym(,,) contains matrices of 'nse' crystal symmetry
!     elements. 
!***********************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON/symmop/nse, nselau, gcsym(24,3,3), bsym(24,4,4)
      DIMENSION b1(4,4),b2(4,4),t(4,4),tt(4,4),ttt(4,4)
      x4dist=5.

      DO 1 ii=1,nselau
      DO 1 jj=1,nselau
         ttt=MATMUL(MATMUL(bsym(ii,:,:),b2),bsym(jj,:,:))

         tt=(b1-ttt)
         t=MATMUL(tt,TRANSPOSE(tt))
         compar=t(1,1)+t(2,2)+t(3,3)+t(4,4)
         IF(compar.LT.x4dist) x4dist=compar

         tt=(TRANSPOSE(b1)-ttt)
         t=MATMUL(tt,TRANSPOSE(tt))
         compar=t(1,1)+t(2,2)+t(3,3)+t(4,4)
         IF(compar.LT.x4dist) x4dist=compar

         ttt(4,1:3)=-ttt(4,1:3)
         ttt(1:3,4)=-ttt(1:3,4)

         tt=(b1-ttt)
         t=MATMUL(tt,TRANSPOSE(tt))
         compar=t(1,1)+t(2,2)+t(3,3)+t(4,4)
         IF(compar.LT.x4dist) x4dist=compar

         tt=(TRANSPOSE(b1)-ttt)
         t=MATMUL(tt,TRANSPOSE(tt))
         compar=t(1,1)+t(2,2)+t(3,3)+t(4,4)
         IF(compar.LT.x4dist) x4dist=compar
  1   CONTINUE
      END


      SUBROUTINE spv4dist(b1,b2,v4dist,bb)
!***********************************************************************
!     Distance between two 'boundaries' 'b1' and 'b2'.
!     The table bsym(,,) contains matrices of 'nse' crystal symmetry
!     elements. 
!     It returns 'v4dist=x4dist' and 'bb' - symmetrically equivalent 
!     to 'b2' and nearest to 'b1'.
!***********************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON/symmop/nse, nselau, gcsym(24,3,3), bsym(24,4,4)
      REAL*8, DIMENSION(4,4) :: b1,b2,t,tt,ttt,bb
 
      v4dist=5.
      DO 1 ii=1,nselau
      DO 1 jj=1,nselau
         ttt=MATMUL(MATMUL(bsym(ii,:,:),b2),bsym(jj,:,:))

         tt=(b1-ttt)
         t=MATMUL(tt,TRANSPOSE(tt))
         compar=t(1,1)+t(2,2)+t(3,3)+t(4,4)
         IF(compar.LT.v4dist) THEN 
           v4dist=compar
           bb=ttt
         ENDIF

         tt=(b1-TRANSPOSE(ttt))
         t=MATMUL(tt,TRANSPOSE(tt))
         compar=t(1,1)+t(2,2)+t(3,3)+t(4,4)
         IF(compar.LT.v4dist) THEN 
           v4dist=compar
           bb=TRANSPOSE(ttt)
         ENDIF

         ttt(4,1:3)=-ttt(4,1:3)
         ttt(1:3,4)=-ttt(1:3,4)

         tt=(b1-ttt)
         t=MATMUL(tt,TRANSPOSE(tt))
         compar=t(1,1)+t(2,2)+t(3,3)+t(4,4)
         IF(compar.LT.v4dist) THEN 
           v4dist=compar
           bb=ttt
         ENDIF

         tt=(b1-TRANSPOSE(ttt))
         t=MATMUL(tt,TRANSPOSE(tt))
         compar=t(1,1)+t(2,2)+t(3,3)+t(4,4)
         IF(compar.LT.v4dist) THEN 
           v4dist=compar
           bb=TRANSPOSE(ttt)
         ENDIF

  1   CONTINUE
      END

!***********************************************************************
      SUBROUTINE symbs_gb
!***********************************************************************
!     Reads symmetry operations.
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON/symmop/nse, nselau, gcsym(24,3,3), bsym(24,4,4)
      OPEN(7,file='SYMMETRY.PRP',status='old',err=11)
      READ(7,*,err=11,END=11)
      READ(7,*,err=11,END=11)
      READ(7,*,err=11,END=11) is,nse,nselau
      bsym=0.
      DO 1 nn=1,nselau
         READ(7,*,err=11,END=11)
         DO 2 i=1,3
   2        READ(7,*,err=11,END=11) (gcsym(nn,i,j),j=1,3)
         bsym(nn,1:3,1:3)=gcsym(nn,:,:)
         bsym(nn,4,4)=1.
  1   CONTINUE
      CLOSE(7,err=11)
      RETURN
  11  WRITE(*,*) 'Cannot open (read) SYMMETRY.PRP.'
      END

!***********************************************************************
      FUNCTION acoss(ang)
!***********************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
      angl=ang
      IF(angl.gt.1) angl=1
      IF(angl.lt.-1) angl=-1
      acoss=ACOS(angl)
      END

      SUBROUTINE axiort(om,ax,g)
!***********************************************************************
!     Rotation axis and angle --> rotation matrix                      *
!***********************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION ax(3), g(3,3)
      d1=ax(1)
      d2=ax(2)
      d3=ax(3)
      co=COS(om)
      so=SIN(om)
      g(1,1)=(1.0-d1*d1)*co+d1*d1
      g(1,2)=d1*d2*(1.0-co)+d3*so
      g(1,3)=d1*d3*(1.0-co)-d2*so
!
      g(2,1)=d1*d2*(1.0-co)-d3*so
      g(2,2)=(1.0-d2*d2)*co+d2*d2
      g(2,3)=d2*d3*(1.0-co)+d1*so
!
      g(3,1)=d1*d3*(1.0-co)+d2*so
      g(3,2)=d2*d3*(1.0-co)-d1*so
      g(3,3)=(1.0-d3*d3)*co+d3*d3
      END

      SUBROUTINE EMat(f1,f,f2,g)
!***********************************************************************
!     Euler angles --> rotation matrix                                 *
!***********************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION g(3,3)
      sf=SIN(f)
      cf=COS(f)
      sf1=SIN(f1)
      cf1=COS(f1)
      sf2=SIN(f2)
      cf2=COS(f2)
      g(1,1)=cf1*cf2-sf1*sf2*cf
      g(1,2)=sf1*cf2+cf1*sf2*cf
      g(1,3)=sf2*sf
!
      g(2,1)=-cf1*sf2-sf1*cf2*cf
      g(2,2)=-sf1*sf2+cf1*cf2*cf
      g(2,3)=cf2*sf
!
      g(3,1)=sf1*sf
      g(3,2)=-cf1*sf
      g(3,3)=cf
      END

      subroutine MEul(g,a1,a,a2)
!***********************************************************************
!     Rotation matrix --> Euler angles                                 *
!***********************************************************************

      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON/globals/ pi, eps

      dimension g(3,3)
      x=g(3,3)
      sf=sqrt(1.-min(x*x,1.))
      if(sf.gt.eps) then
        a1=-g(3,2)/sf
        a1=acoss(a1)
        if(g(3,1).lt.0.) a1=2.*pi-a1
        a=acos(x)
        a2=g(2,3)/sf      
        a2=acoss(a2)
        if(g(1,3).lt.0.) a2=2.*pi-a2
      else 
        a1=g(1,1)
        a1=acoss(a1)
        if(g(1,2).lt.0.) a1=2.*pi-a1
        a=acoss(x)
        a2=0.
      endif
      end


!***********************************************************************
  SUBROUTINE force_vector(g,xn, xi)
!***********************************************************************
! Determines the force vector 'xi' for a given boundary '(g,xn)'.
! Second part of Appendix A.
!***********************************************************************
  IMPLICIT REAL*8 (A-H,O-Z)
  COMMON/xcusp1/number_of_cusps
  COMMON/xcusp2/depth_of_cusp(50),width_of_cusp(50),cusp(50,4,4)
  REAL*8, INTENT(IN) , DIMENSION(3,3) :: g(3,3), xn(3)
  REAL*8, INTENT(OUT), DIMENSION(3) :: xi(3)
  REAL*8, DIMENSION(4,4) :: b, bcusp, bb
  REAL*8, DIMENSION(3,3) :: gx, gm
  REAL*8, DIMENSION(3)   :: xnm, dgam_dn
  REAL*8    :: a, x, w, coef, coefx
  INTEGER :: i, isum

  dgam_dn=0.

  a=depth_of_cusp(1)
  gx=cusp(1,1:3,1:3)
  x=xmdist(g,gx)/width_of_cusp(1)
  coef=auxilx(x,a)
  CALL bconvert(g,xn,b)

  sumowanie: DO isum=2,number_of_cusps
     coefx=coef
     a=depth_of_cusp(isum)
     bcusp=cusp(isum,:,:)
     w=width_of_cusp(isum)
     CALL spv4dist(b,bcusp,x,bb)
     x=MAX(x/w,1.0e-10)
     IF (x >= 1.) CYCLE sumowanie
     coefx=coefx*((a*LOG(x))/(2.*w*x))
     xnm=bb(1:3,4)
     gm=bb(1:3,1:3)
     DO i=2,number_of_cusps
        IF(i.EQ.isum) CYCLE
        a=depth_of_cusp(i)
        bcusp=cusp(i,:,:)
        x=x4dist(b,bcusp)/width_of_cusp(i)
        coefx=coefx*auxilx(x,a)
     ENDDO
     dgam_dn=dgam_dn+coefx*(xnm+MATMUL(MATMUL(g,TRANSPOSE(gm)),xnm))
  END DO  sumowanie 

  xi= energy_gn(g,xn)*xn + (dgam_dn - DOT_PRODUCT(xn,dgam_dn)*xn)

  END 

!***********************************************************************
  SUBROUTINE triple_line(sl,sn)
!***********************************************************************
!  Generates the direction 'sl' of the triple line and normal 'sn'
!  to one of the boundaries.
!  Vectors 'sl' and 'sn' are perpendicular.
!  Both vectors are given in the sample coordinate system.
!***********************************************************************
  IMPLICIT REAL*8 (A-H,O-Z)
  COMMON/globals/ pi, eps
  REAL*8, INTENT(OUT), DIMENSION(3) :: sl, sn
  REAL*8 xsn

  CALL gener_n(sl) 
1 CALL gener_n(sn) 
  sn=sn-DOT_PRODUCT(sl,sn)*sl
  xsn=SQRT(DOT_PRODUCT(sn,sn))
  IF(xsn > eps) THEN
    sn=sn/xsn
  ELSE
    goto 1
  ENDIF
  RETURN
  END

!***********************************************************************
  SUBROUTINE triple_junction(om2,om3)
!***********************************************************************
!  For given direction 'sl' of the triple line and normal 'sn1'
!  to one of the boundaries, 
!  calculates the normals 'sn2' and 'sn3' to other two boundaries.
!  Vectors 'sni'are perpendicular to 'sl'.
!  Angle between 'sn1' and 'sn2' is 'om2' (clockwise).
!  Angle between 'sn1' and 'sn3' is 'om3' (clockwise).
!  All vectors are given in the sample coordinate system.
!***********************************************************************
  IMPLICIT REAL*8 (A-H,O-Z)
  COMMON/trijun/g1(3,3),g2(3,3),g3(3,3),gx1(3,3),gx2(3,3),gx3(3,3),cn1(3),cn2(3),cn3(3),sn1(3),sn2(3),sn3(3),sl(3)
  REAL*8, DIMENSION(3,3) :: g(3,3)

  CALL axiort(om2,sl,g)
  sn2=MATMUL(TRANSPOSE(g),sn1)
  CALL axiort(om3,sl,g)
  sn3=MATMUL(TRANSPOSE(g),sn1)
  RETURN
  END

!!-----------------------------------------------------------------------
!
!!***********************************************************************
!      LOGICAL FUNCTION INTRAC(DUMMY)
!!***********************************************************************
!!     Required by MINUIT.
!      IMPLICIT REAL*8 (A-H,O-Z)
!      INTRAC =.FALSE.
!      END


!-----------------------------------------------------------------------
