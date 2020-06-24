        include "DCUHRE.f"
        
        program main
            parameter (nd=2,nf=1,nw=1500)
            real*8 am(nd),bm(nd),work(nw)
            real*8 abserr(nf),fnuest(nf),epsabs,epsrel
            external dcuhre,sub1
            am(1)=0d0
            bm(1)=1d0
            am(2)=0d0
            bm(2)=1d0
            minp=65
            maxp=2000
            key=0
            epsabs=0
            epsrel=1e-3
            call dcuhre(nd,nf,am,bm,minp,maxp,sub1,epsabs
     &  ,epsrel,key,nw,0,fnuest,abserr,neval,ifail,work)
            fnu=fnuest(1)
            print *,fnu
        end

        subroutine sub1(nd,va,nf,fun1)
            real*8 va(nd),fun1(nf)
            x1=va(1)
            x2=va(2)
            fun1(1)=x1**2+x2
            return
        end subroutine