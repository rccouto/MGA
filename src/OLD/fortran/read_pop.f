	subroutine read_pop_ini(pot_chos,filename,m1,nclust,natoms
     &                         ,namo,energ,v2m,v3m,pop_ini)
	implicit real*8(a-h,o-z)
	include "param.inc"
c       
	common/el/     elem	
	real*8         pop_ini(n_max_clus,n_max_atom)
	1             ,energ(n_max_clus)
	1             ,v2m(n_max_clus),v3m(n_max_clus)
	character*3    elem(n_max_atom)
	character*10   co1,co2,co3,co4	
	character*20   pot_chos,filename
c
	if(pot_chos(1:3).eq.'tip')then
	   do i=1,nclust,1
	      open(unit=9,status='old',file=filename(1:m1))
	      read(9,*)natoms
 13	      format(1x,a7,i3,2x,a8,f14.6,4x,a8,f14.6,2x,a8,f14.6)
	      read(9,13)co1,iclus,co2,energ(i)
	1	   ,co3,v2m(i),co4,v3m(i)  
	      do j=1,namo,1
 14		 format(a3,3x,f14.6,4x,f14.6,4x,f14.6
	1	      ,4x,f14.6,4x,f14.6,4x,f14.6)
		 read(9,14)elem(j),pop_ini(i,j),pop_ini(i,j+namo)
	1	      ,pop_ini(i,j+2*namo),pop_ini(i,j+3*namo)
	2	      ,pop_ini(i,j+4*namo),pop_ini(i,j+5*namo)

	      enddo
	   enddo
	else
	   do i=1,nclust,1
	      open(unit=9,status='old',file=filename(1:m1))
	      read(9,*)natoms
 23	      format(1x,a7,i3,2x,a8,f14.6,4x,a8,f14.6,2x,a8,f14.6)
	      read(9,23)co1,iclus,co2,energ(i)
	1	   ,co3,v2m(i),co4,v3m(i)  
	      do j=1,natoms
 24		 format(a3,3x,f14.6,4x,f14.6,4x,f14.6)
		 read(9,24)elem(j),pop_ini(i,j)
	1	      ,pop_ini(i,j+natoms),pop_ini(i,j+2*natoms)
	      enddo
	   enddo
	endif
	close(9)
c
	return
	end
                                                                       
