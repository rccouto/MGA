	subroutine muting(pot_chos,opt_mut,per_mut,nseed,nms,naux1
	1    ,nclust,natoms,nmono,n_count,noff,nat,pop_clus)
	implicit real*8(a-h,o-z)
	include "param.inc" 
c
	integer         nat(10)
	real*8          pop_clus(n_max_clus,n_max_atom)
	1              ,clus_mut(n_max_clus,n_max_atom)
	2              ,xyz(n_max_atom)
	character*20    pot_chos,opt_mut
c
	n_count=0
	mnt=nmono
	if(opt_mut(1:8).eq.'rot_chan')per_mut=per_mut/2d0
	if(pot_chos(1:3).eq.'tip')then
	   mnt=1
	   namo=natoms/nmono
	else
	   namo=natoms
	endif
c
	call matr1_matr2(2*nclust,3*natoms,pop_clus,clus_mut)
c       
	call sort_atoms_z(pot_chos,nclust,natoms,nmono,nat,clus_mut)
c
cxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
c	      do i=1,2*nclust,1
c		 write(41,*)'clust',i,'      energy'
c		 do j=1,namo
c		    write(41,*)j,clus_mut(i,j),clus_mut(i,j+namo)
c	1		 ,clus_mut(i,j+2*namo)
c		    if(pot_chos(1:3).eq.'tip')then
c		       write(41,*)clus_mut(i,j+3*namo),clus_mut(i,j+4*namo)
c	1		    ,clus_mut(i,j+5*namo)
c		    endif
c		 enddo
c	      enddo
c	      close(41)
cxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
c
	do n_ran=1,2*nclust,1  !robson 05/04 alterar?
	   if(per_mut/1d2.lt.ran3(nseed))goto 1111 !robson 05/04
	   n_count=n_count+1	  	
	   if(opt_mut(1:6).eq.'rotate'
	1	.or. opt_mut(1:8).eq.'rot_chan')then
	      nii=1
	      nf=0
	      call matr_vec(n_ran,3*natoms,xyz,clus_mut)
	      call centre(namo,re,xyz)
c	      nms=nat(1)
	      n_rot=nms*sqrt(ran3(nseed)) !robson 05/04 
	      do 100 k=1,mnt,1
		 if(n_rot.ge.natoms)n_rot=natoms
		 ni=nii
		 nf=nf+nat(k)
		 call rotate(namo,nseed,xyz,ni,n_rot,k)
		 n_rot=n_rot+nat(k)
		 nii=nf+1
 100	      continue
c       
	      call vec_matr(2*nclust+n_count,3*natoms,xyz,pop_clus)
c       
	      if(opt_mut(1:8).eq.'rot_chan')goto 111     
	   elseif(opt_mut(1:9).eq.'chan_part')then
 111	      continue
	      n_count=n_count+1
	      typ=natoms**(1d0/3d0)
	      typaux=typ/2d0
c       
	      do j=1,namo,1
		 pop_clus(2*nclust+n_count,j)
	1	      =clus_mut(n_ran,j)
		 pop_clus(2*nclust+n_count,j+namo)
	1	      =clus_mut(n_ran,j+namo)
		 pop_clus(2*nclust+n_count,j+2*namo)
	1	      =clus_mut(n_ran,j+2*namo)
		 if(pot_chos(1:3).eq.'tip' .and. 
	1	      20d-2.lt.ran3(nseed))then   !robson 05/04
		    pop_clus(2*nclust+n_count,j+3*namo)
	1		 =clus_mut(n_ran,j+3*namo)
		    pop_clus(2*nclust+n_count,j+4*namo)
	1		 =clus_mut(n_ran,j+4*namo)
		    pop_clus(2*nclust+n_count,j+5*namo)
	1		 =clus_mut(n_ran,j+5*namo)
		 endif
		 if(ran3(nseed).lt.20d-2)then !robson 05/04 
		    pop_clus(2*nclust+n_count,j)
	1		 =clus_mut(n_ran,j)*ran3(nseed) !robson 05/04 
		    pop_clus(2*nclust+n_count,j+namo)
	1		 =clus_mut(n_ran,j+namo)*ran3(nseed) !robson 05/04
		    pop_clus(2*nclust+n_count,j+2*namo)
	1		 =clus_mut(n_ran,j+2*namo)*ran3(nseed) !robson 05/04
		 endif
		 if(pot_chos(1:3).eq.'tip' .and. 
	1	      20d-2.lt.ran3(nseed))then    !robson 05/04
		    pop_clus(2*nclust+n_count,j+3*namo)
	1		 =clus_mut(n_ran,j+3*namo)*ran3(nseed) !robson 05/04
		    pop_clus(2*nclust+n_count,j+4*namo)
	1		 =clus_mut(n_ran,j+4*namo)*ran3(nseed) !robson 05/04
		    pop_clus(2*nclust+n_count,j+5*namo)
	1		 =clus_mut(n_ran,j+5*namo)*ran3(nseed) !robson 05/04
		 endif
	      enddo
c
	   elseif(opt_mut(1:8).eq.'new_clus')then 
	      n_count=n_count+1 	  
	      call gener_init_pop(pot_chos,1,natoms,nmono,nseed
	1	   ,naux1,pop_clus(2*nclust+n_count,1))
c       
	   else
	      write(*,*)'Muting option not found'
	      write(*,*)'error: Muting'
	      stop
	   endif	
c       
 1111	   continue
	enddo
c       
	return
	end

