      function ran3(idum)
      implicit none
c     
c     Numerical Recipes Random Number Generator (ran3)
c
      integer idum,mbig,mseed,mz
      real*8 ran3,fac
      parameter (mbig=1000000000,mseed=161803398,mz=0,fac=1./mbig)
      integer i,iff,ii,inext,inextp,k
      integer mj,mk,ma(55)
c     
      save iff,inext,inextp,ma
      data iff /0/
c     
      if (idum.lt.0 .or. iff.eq.0) then
         iff=1
         mj=mseed-iabs(idum)
         mj=mod(mj,mbig)
         ma(55)=mj
         mk=1
         do 10 i=1,54
            ii=mod(21*i,55)
            ma(ii)=mk
            mk=mj-mk
            if(mk.lt.mz) mk=mk+mbig
            mj=ma(ii)
 10      enddo
         do 20 k=1,4
            do 30 i=1,55
               ma(i)=ma(i)-ma(1+mod(i+30,55))
               if(ma(i).lt.mz) ma(i)=ma(i)+mbig
 30         enddo
 20      enddo
         inext=0
         inextp=31
         idum=1
      end if
      inext=inext+1
      if (inext.eq.56) inext=1
      inextp=inextp+1
      if(inextp.eq.56) inextp=1
      mj=ma(inext)-ma(inextp)
      if(mj.lt.mz) mj=mj+mbig
      ma(inext)=mj
      ran3=mj*fac
      return
      end
