cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine print_1(filename,m1,nclust,iclus,ngen,en_aver
	1                 ,energ_bet)
	implicit real*8(a-h,o-z)
	include "param.inc"
c	
	real*8         energ_bet(nclust)
        character*5    tet
	character*20   filename,newnam1,newnam2
c
	en_aver=sum_vec(nclust,energ_bet)/nclust
c
13	format(1x,a3,4x,a13,2x,f14.6,4x,a3)
14	format(1x,a3,4x,a13,5x,i5,10x,a3)
	write(*,*) '###########################################'
	write(*,*) '##***************************************##'
	write(*,14)'##*','generation  :',ngen,             '*##'
	write(*,*) '##***************************************##'
	write(*,13)'##*','energy_min  =',energ_bet(1),     '*##'
	write(*,13)'##*','energy_aver =',en_aver,          '*##'
	write(*,13)'##*','energy_max  =',energ_bet(nclust),'*##'
        write(*,*) '##***************************************##'
	write(*,*) '###########################################'
	write(*,*)' '
c
	write(tet,'(i3)')iclus
	if(iclus.lt.10)then
	newnam1=filename(1:m1)//'_epp.'//tet(3:3)
	newnam2=filename(1:m1)//'_eppa.'//tet(3:3)
	elseif(iclus.lt.100)then
	newnam1=filename(1:m1)//'_epp.'//tet(2:3)
	newnam2=filename(1:m1)//'_eppa.'//tet(2:3)
	else
	newnam1=filename(1:m1)//'_epp.'//tet(1:3)
	newnam2=filename(1:m1)//'_eppa.'//tet(1:3)
	endif
c	
	open(unit=2,status='unknown',file=newnam1(1:20))
23 	format(1x,i4,3x,f14.6,4x,f14.6,4x,f14.6)
	write(2,23)ngen,energ_bet(1),en_aver,energ_bet(nclust)
c	
	open(unit=3,status='unknown',file=newnam2(1:20))
	write(3,*)ngen,energ_bet
c
	return
	end	
	
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine print_2(pot_chos,filename,m1,m2,nclust,natoms,namo
	1                 ,iclus,energ_best,v2tot,v3tot,pop_clus)
	implicit real*8(a-h,o-z)
	include "param.inc"
c
	common/el/     elem
	real*8         pop_clus(n_max_clus,n_max_atom)
	1             ,clust(3*natoms)
	character*3    elem(n_max_atom)
	character*20   pot_chos,filename	
c
c	if(v2tot.gt.9999)v2tot=9999
c	if(v3tot.gt.9999)v3tot=9999
c	if(pot_chos(1:3).eq.'tip')then
c	   call matr_vec(1,6*namo,clust,pop_clus)
c	   call centre(natoms,re,clust)
c	   open(unit=m2,status='unknown',file=filename(1:m1))
c	   write(m2,'(i5)')natoms
c 13	   format(1x,a7,i3,2x,a8,f14.6,4x,a8,f14.6,2x,a8,f14.6)
c	   write(m2,13)'Frame =',iclus,'energy =',energ_best
c	1	,'v2_tot =',v2tot,'v3_tot =',v3tot   
c	   do j=1,namo,1
c 14	      format(a3,3x,f14.6,4x,f14.6,4x,f14.6
c	1	   ,4x,f14.6,4x,f14.6,4x,f14.6)
c	      write(m2,14)elem(j),clust(j),clust(j+namo),clust(j+2*namo)
c	1	   ,clust(j+3*namo),clust(j+4*namo),clust(j+5*namo)
c	   enddo
c	else
	   call matr_vec(1,3*natoms,clust,pop_clus)
	   call centre(natoms,re,clust)
	   open(unit=m2,status='unknown',file=filename(1:m1))
	   write(m2,'(i5)')natoms
 23	   format(1x,a7,i3,2x,a8,f14.6,4x,a8,f14.6,2x,a8,f14.6)
	   write(m2,23)'Frame =',iclus,'energy =',energ_best
	1	,'v2_tot =',v2tot,'v3_tot =',v3tot       
	   do j=1,natoms
 24	      format(a3,3x,f14.6,4x,f14.6,4x,f14.6)
	      write(m2,24)elem(j),clust(j),clust(j+natoms),
     &	                  clust(j+2*natoms)
	   enddo
c	endif
c       
	return
	end



