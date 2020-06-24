        include "DGAUSS.f"

        module minuit_mod
            implicit none
            real*8 a,b,c
            real*8 x1,x2
            integer k
            contains

            real*8 function ff(kk,x)
            real*8 x,c1,c2,eps
            integer kk
            real*8,external:: dgauss1
            k=kk
            c1=0.
            c2=1.
            eps=1d-3
            select case (k)
            case (1)
                    x2=x
            case (2)
                    x1=x
            end select
            ff=dgauss1(f1,c1,c2,eps)
            return
            end function

            real*8 function f1(x)
            real*8 x
            select case (k)
            case (1)
                    f1=f(x,x2)
            case (2)
                    f1=f(x1,x)
            end select
            return
            end function

            real*8 function f(x,y)
            real*8 x,y
            f=a*x*y+b*x+c*y
            return
            end function
        end module

      subroutine FCN(NPAR,GRAD,FVAL,XVAL,IFLAG,FUTIL)
        use minuit_mod
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        DIMENSION XVAL(3),GRAD(*)
        DIMENSION dataexp(3,2),dataerr(3,2),point(3,2)
                
        dataexp=reshape((/4.5,8.,11.5,1.5,4.,6.5/),(/3,2/))
        dataerr=reshape((/0.1,0.1,0.1,0.1,0.1,0.1/),(/3,2/))
        point=reshape((/1.,2.,3.,0.,1.,2./),(/3,2/))

        a=xval(1)
        b=xval(2)
        c=xval(3)
        Fval=0D0
        do i=1,2
                do j=1,3
                theoval=ff(i,point(j,i))
                !theoval=ff(point(j))
                fval=fval+(theoval-dataexp(j,i))**2/dataerr(j,i)**2
                end do
        end do
        return 
      end

      program main
        implicit double precision(A-H,O-Z)
        double precision VSTRT ,STP, FCN,zero
        EXTERNAL FCN
        DIMENSION NPRM(3),VSTRT(3),STP(3)
        CHARACTER*10 PNAM(3)
        DATA NPRM /1,2,3/
        DATA PNAM /'a','b','c'/
        DATA VSTRT/100d0,100d0,100d0/
        DATA STP/0.1d0,1d0,1d0/

        CALL MNINIT(5,6,7)
        
        do i=1,3
        CALL MNPARM(NPRM(i),PNAM(i),vstrt(i),stp(i),0d0,0d0,ierflg)
        end do 

        call mnseti('test1111')
        CALL MNEXCM(FCN, 'CALL FCN',1 ,1,IERFLG, 0)
        call mnexcm(fcn,'set print',0,0,ierflg,0)
        CALL MNEXCM(fcn,'MIGRAD', 0 ,0,IERFLG,0)
        call mnexcm(fcn,'minos',0,0,ierflg,0)
  
        do i=1,3
        call mnpout(i,PNAM(i),val,error,bnd1,bnd2,ivarbl)
        print *, val,error,bnd1,bnd2
        end do    
        stop
        end