	subroutine cor_const(pot_chos,n,m)
	implicit real*8(a-h,o-z)
	include "common.inc"  
c
	common/cst/    cst1,cst2,cst3,cst4,cst5,cst6,cst7
     &                ,cst8,cst9,cst10,cst11
c
	real*8         cst1(7,7),cst2(7,7),cst3(7,7),cst4(7,7)
     &                ,cst5(7,7),cst6(7,7),cst7(7,7),cst8(7,7)
     &                ,cst9(7,7),cst10(7,7),cst11(7,7)
c
	character*15    pot_chos
c
	if(n.gt.m)then
	  nax=n
	  n=m
	  m=nax
	endif
c
	if(pot_chos(1:5).eq.'tip3p')then
	   qi=cst1(n,m)
	   qj=cst2(n,m)
	   At=cst3(n,m)
	   Ct=cst4(n,m)
	   re=cst5(n,m)
c
	else
	   write(*,*)'correction constats not found'
	   write(*,*)'error: cc'
	   stop
	endif
c		
	return
	end

	   
