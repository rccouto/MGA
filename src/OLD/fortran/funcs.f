C- function "ran" that generat random namber.
C- Function ran()	
      function ran(nseed)
      implicit real*8(a-h,o-z)
      parameter(l=1029,nc=221591,m=1048576)
      nseed=mod(nseed*l+nc,m)
      ran=1d0*nseed/m
      return
      end
C-cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc	
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
      write(*,*)'inside ran3'
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
C-cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc	
      function sum_vec(n,vec)
      implicit real*8(a-h,o-z)
c
      real*8         vec(n)
c
      sum_vec=0d0
c
      do k=1,n,1
         sum_vec=sum_vec+vec(k)
      enddo
c
      return
      end
C-cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      function naleat(nseed,n,ditr_chos)
      implicit real*8(a-h,o-z)
c
      character*15   ditr_chos
c
C	x=ran(nseed) robson 05/04
	x=ran3(nseed)
	if(ditr_chos(1:10).eq.'ditr_divx1')then
           n1=2*n
           distr=1d0/(x**0.5+1d0)
           naleat=int(n1*distr)-n1/2
c     
        elseif(ditr_chos(1:10).eq.'ditr_gaus1')then
           cst=0d0
           distr=exp(-(x**2-cst)**2)
           naleat=n-int(n*distr) 
      	  
c     
      	elseif(ditr_chos(1:10).eq.'ditr_divx2')then
           n1=2*n
           distr=1d0/(x+1d0)
           naleat=int(n1*distr)-n1/2     	
c     
      	elseif(ditr_chos(1:10).eq.'ditr_ddddd')then
           n1=n
           distr=(1d0-x)**2
           naleat=int(n1*distr)
c
      	else
           write(*,*)'Warning!!! distribution error Warning!!!'
           write(*,*)'Caution!!! distribution error Caution!!! '
	endif
c
	return
	end
      

	
