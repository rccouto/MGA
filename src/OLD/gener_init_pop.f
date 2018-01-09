	subroutine gener_init_pop(pot_chos,nclust,natoms,nmono,nseed
	1                        ,naux1,pop_clus)
	implicit real*8(a-h,o-z)
	include "param.inc"
c
	real*8         pop_clus(n_max_clus,n_max_atom)
	character*20   pot_chos
c	
	pi=acos(-1d0)
	namo=natoms
	if(pot_chos(1:3).eq.'tip')namo=natoms/nmono
	typ=dble(natoms)**(1d0/3d0)
c
c   	typ=typ/2d0
	typaux=typ/2d0
	do ir=1,nclust,1
	   do jr=1,namo,1
	      pop_clus(ir,jr)=typ*ran3(nseed)-typaux
	      pop_clus(ir,jr+namo)=typ*ran3(nseed)-typaux
	      pop_clus(ir,jr+2*namo)=typ*ran3(nseed)-typaux
c	write(*,*)pop_clus(ir,jr), pop_clus(ir,jr+namo)
c     &            ,pop_clus(ir,jr+2*namo)
	      if(pot_chos(1:3).eq.'tip')then
		 pop_clus(ir,jr+3*namo)=pi*ran3(nseed)
		 pop_clus(ir,jr+4*namo)=2d0*pi*ran3(nseed)
		 pop_clus(ir,jr+5*namo)=2d0*pi*ran3(nseed) 
	      endif
	   enddo
	enddo	
c
	return
	end

