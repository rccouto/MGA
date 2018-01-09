	subroutine pimss(pot_chos,opt_ga,opt_min,nclust,natoms,nmono
	1    ,n_cl_tot,naux1,enr_min,enr_bef,enr_1,enr_nxt_m,vm2,vm3
	2    ,vm2t,vm3t,energ_bet,energ,pop_clus,pop_c)
	implicit real*8(a-h,o-z)
	include "param.inc"
c	include "common.inc"
c       
	real*8         pop_clus(n_max_clus,n_max_atom)
	1             ,pop_c(n_max_clus,n_max_atom)
	2             ,energ_bet(n_max_clus)
	3             ,energ(n_max_clus)
	4             ,vm2(n_max_clus),vm3(n_max_clus)
	5             ,vm2t(n_max_clus),vm3t(n_max_clus)
	character*20   pot_chos,opt_min,opt_ga
	logical        task1
c       
	task1=.false.
c       
	namo=natoms/nmono
c
c	write(*,*)'out read',n_cl_tot
	do l=1,n_cl_tot,1
	   vm2t(l+2*nclust)=vm2(l)
	   vm3t(l+2*nclust)=vm3(l)
	   energ(l+2*nclust)=energ_bet(l)
	enddo
c	write(*,*)'pimssin1'
	call matr1_matr2(nclust,3*natoms,pop_clus(1,1)
	1    ,pop_c(2*nclust+1,1))
c
cxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
c	   do i=1,3*nclust
c	      write(59,*)'clust',i,' energy',energ(i),vm2t(i),vm3t(i)
c	      do j=1,natoms
c		 write(59,*)j,pop_c(i,j),pop_c(i,j+natoms)
c	1	      ,pop_c(i,j+2*natoms)
c	      enddo
c	   enddo
c	   close(59)
cxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx	
c
c	write(*,*)'pimssin2'
	call sort_matr(pot_chos,3*nclust,natoms
	1    ,vm2t,vm3t,energ,pop_c)
c	write(*,*)'pimssout'
c
	enr_exs=energ(1)
c	write(*,*)'pimssout2'
	kit=3*nclust
	do 100 i=2,3*nclust,1
	   if(abs(energ(i)-enr_exs).le.enr_nxt_m  
	1	.or. energ(i).eq.0d0)then
	      kit=kit-1
	      energ(i)=0d0
	   endif
	   enr_exs=energ(i)
	   k=i
	   if(energ(k).eq.0d0)then
 111	      continue
	      k=k-1
	      if(energ(k).eq.0d0)goto 111
	      enr_exs=energ(k)
	   endif
 100	continue
c       
c	write(*,*)'predin'
	call predator(pot_chos,'next',opt_min,3*nclust,natoms
	1    ,nat,iclus,nseed,naux1,2,nmono,noff,enr_min
	2    ,enr_bef,enr_1,enr_nxt_m,vm2,vm3,energ
	3    ,pop_c,task1)
c	write(*,*)'predout'
c       
	call sort_matr(pot_chos,3*nclust,natoms
	1    ,vm2t,vm3t,energ,pop_c)
c	write(*,*)'********kit',kit
cxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
c	do i=1,2*nclust
c	   write(69,*)'clust',i,' energy',energ(i),vm2t(i),vm3t(i)
c	   do j=1,natoms
c	      write(69,*)j,pop_c(i,j),pop_c(i,j+natoms)
c	1	   ,pop_c(i,j+2*natoms)
c	   enddo
c	enddo
c	close(69)
cxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
	kt=0
	do i=kit+1,2*nclust
	   kt=kt+1
	   energ(i)=energ(kt) 
	   call matr1_matr2(1,3*natoms,pop_c(kt,1)
	1	,pop_c(i+1,1))
	enddo
cxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
c	   do i=1,2*nclust
c	      write(68,*)'clust',i,' energy',energ(i),vm2t(i),vm3t(i)
c	      do j=1,natoms
c		 write(68,*)j,pop_c(i,j),pop_c(i,j+natoms)
c	1	      ,pop_c(i,j+2*natoms)
c	      enddo
c	   enddo
c	   close(68)
cxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
c       
	if(opt_ga(1:6).eq.'noelit')then
	   do i=1,nclust,1
	      energ_bet(i)=energ(i)
	      vm2(i)=vm2t(i)
	      vm3(i)=vm3t(i)
	   enddo
	   call matr1_matr2(nclust,3*natoms,pop_c
	1	,pop_clus)
	endif
c	
	end
	






