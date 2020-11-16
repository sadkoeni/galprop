c-------------------------------------------------------------------------------
      subroutine aprtab(filename,fnlength)
c-------------------------------------------------------------------------------
c Read in the interpolation table for the apspec routine
c-------------------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      common/fit/tablef(100,7,4,2),alphap,betap
      integer fnlength
      character*(fnlength) filename
      
      open(1,file=filename,status='old')
       read(1,*)tablef,alphap,betap
      close(1)

      end
      

c-------------------------------------------------------------------------------
      double precision function apspec(e0,epbar,iap,iat) 
c-------------------------------------------------------------------------------
c pbar+nbar spectrum in lab. frame for pp, pA, Ap, & AA collisions:
c = Epbar * d[sig(E0,Epbar)]  / d[Epbar]; based on the modified QGSJET-II model
c  [Kachelriess, Moskalenko & Ostapchenko, arXiv:1502.04158]
c-------------------------------------------------------------------------------
c e0    - incident particle (total) energy (per nucleon) in GeV;
c epbar - antiproton (total) energy in GeV;
c iap   - incident particle mass number (=1 for proton);
c iat   - target mass number: 1 or 4 (only proton or helium)   
c-------------------------------------------------------------------------------
      implicit double precision(a-h,o-z)
      dimension wi(3),wk(3),ff(3,3)
      common/fit/tablef(100,7,4,2),alpha,beta

      if(iat.ne.1.and.iat.ne.4)stop'only proton or helium as a target'
      if(iap.gt.60)write(*,*)'warning: tabulation done till Fe'
      
      amp=.938d0
      apspec=0.d0
      if(e0.le.7.d0*amp)return
      
      eps0=e0/amp
      rrr=(1.d0-3.d0/eps0+4.d0/eps0**2)/(1.d0+eps0)
      xmax=.5d0*(1.d0-3.d0/eps0+dsqrt((1.d0-3.d0/eps0)**2
     *-4.d0*rrr))
      xmin=2.d0*rrr/(1.d0-3.d0/eps0+dsqrt((1.d0-3.d0/eps0)**2
     *-4.d0*rrr))

      x=epbar/e0
      if(x.le.xmin.or.x.ge.xmax)return
      
      xx=dlog(x/xmin)/dlog(xmax/xmin)*101.d0
      ix=max(1,int(xx))
      ix=min(ix,98)
      wi(2)=xx-ix
      wi(3)=wi(2)*(wi(2)-1.d0)*.5d0
      wi(1)=1.d0-wi(2)+wi(3)
      wi(2)=wi(2)-2.d0*wi(3)

      if(e0.lt.1.d2)then
       ie=1
       ee=dlog10(e0)*2.d0-1.d0
      else
       ee=dlog10(e0)+1.d0
       ie=max(3,int(ee))
       ie=min(5,ie)
      endif
      wk(2)=ee-ie
      wk(3)=wk(2)*(wk(2)-1.d0)*.5d0
      wk(1)=1.d0-wk(2)+wk(3)
      wk(2)=wk(2)-2.d0*wk(3)
      
      jj=1+(iat/4)
      if(iap.eq.1)then
       do i=1,3
       do k=1,3
        ff(i,k)=tablef(ix+i-1,ie+k-1,1,jj)
       enddo
       enddo
      elseif(iap.eq.4)then
       do i=1,3
       do k=1,3
        ff(i,k)=tablef(ix+i-1,ie+k-1,2,jj)
       enddo
       enddo
      else
       wa1=dlog(dble(iap))/dlog(4.d0)*dlog(56.d0/dble(iap))/dlog(14.d0)
       wa2=dlog(dble(iap))/dlog(56.d0)*dlog(dble(iap)/4.d0)/dlog(14.d0)
       wa0=1.d0-wa1-wa2
       do i=1,3
       do k=1,3
        ff(i,k)=tablef(ix+i-1,ie+k-1,2,jj)*wa0+tablef(ix+i-1
     *  ,ie+k-1,3,jj)*wa1+tablef(ix+i-1,ie+k-1,4,jj)*wa2
       enddo
       enddo
      endif

      do i=1,3
      do k=1,3
       apspec=apspec+ff(i,k)*wi(i)*wk(k)
      enddo
      enddo
      apspec=exp(apspec)*(1.d0-x/xmax)**alpha*(1.d0-xmin/x)**beta
      return
      end
      
