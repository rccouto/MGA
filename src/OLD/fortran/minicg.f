	subroutine minicg(pot_chos,nclust,natoms,nmono,vm2,vm3
	1    ,energ,pop)
	implicit real*8(a-h,o-z)
	include "param.inc"
	external funct
c       
        integer        n,maxfn,ier,i,m       
        real*8         pop(n_max_clus,n_max_atom),w(20*natoms)
	1             ,xyz(3*natoms),energ(nclust),vm2(nclust)
	2             ,vm3(nclust),g(3*natoms)
	character*20   pot_chos
c       
c- minimisation options
        acc=1d-5
        dfpred=10
        maxfn=0
	namo=natoms/nmono
	do i=1,nclust,1
	   if(pot_chos(1:3).eq.'tip')then
	      do k=1,6*namo,1
		 xyz(k)=pop(i,k)
	      enddo
	   else
	      do k=1,3*natoms,1
		 xyz(k)=pop(i,k)
	      enddo
	   endif
c       
	   if(pot_chos(1:5).eq.'morse'
	1	.or. pot_chos(1:9).eq.'tlr_morse'  
	2	.or. pot_chos(1:5).eq.'mm_2b')then
	      call zxcgr(funct,3*natoms,acc,maxfn,dfpred,xyz,g,v,w,ier)
	      call diffv2(pot_chos,natoms,xyz,v,g)
           elseif(pot_chos(1:5).eq.'mm_cp')then
	      call zxcgr(funct,3*natoms,acc,maxfn,dfpred,xyz,g,v,w,ier)
	      call diffvmm(pot_chos,natoms,xyz,v2t,v3t,v,g)
	   elseif(pot_chos(1:5).eq.'gupta')then
	      call zxcgr(funct,3*natoms,acc,maxfn,dfpred,xyz,g,v,w,ier)
	      call diffvgp(pot_chos,natoms,xyz,v2t,v3t,v,g)
	   elseif(pot_chos(1:5).eq.'tip3p')then
	      call zxcgr(funct,6*namo,acc,maxfn,dfpred,xyz,g,v,w,ier)
	      call diffvtip3p(pot_chos,natoms,nmono,xyz,v2t,v3t,v,g)
	   elseif(pot_chos(1:5).eq.'tip4p')then
	      call zxcgr(funct,6*namo,acc,maxfn,dfpred,xyz,g,v,w,ier)
	      call diffvtip4p(pot_chos,natoms,nmono,xyz,v2t,v3t,v,g)
	   elseif(pot_chos(1:5).eq.'tip5p')then
	      call zxcgr(funct,6*namo,acc,maxfn,dfpred,xyz,g,v,w,ier)
	      call diffvtip5p(pot_chos,natoms,nmono,xyz,v2t,v3t,v,g)   
	   else
	      write(*,*)' This potential does not exist'
	      write(*,*)'error:mini'
	      stop
	   endif
c 
	   energ(i)=v
	   vm2(i)=v2t
	   vm3(i)=v3t	
c
	   if(pot_chos(1:3).eq.'tip')then
	      do k=1,6*namo,1
		 pop(i,k)=xyz(k)
	      enddo
	   else
	      do k=1,3*natoms,1
		 pop(i,k)=xyz(k)
	      enddo
	   endif
	enddo
c
	return
	end

	subroutine funct(n,xyz,v,g)
	implicit real*8(a-h,o-z)
	common/opot/pot_chos,nmono
c
	real*8         xyz(n),g(n)
	character*20   pot_chos	
c	
c-
	if(pot_chos(1:5).eq.'morse'
	1    .or. pot_chos(1:9).eq.'tlr_morse'  
	2    .or. pot_chos(1:5).eq.'mm_2b')then
	   call diffv2(pot_chos,n/3,xyz,v,g)
	elseif(pot_chos(1:5).eq.'mm_cp')then
	   call diffvmm(pot_chos,n/3,xyz,v2t,v3t,v,g)
	elseif(pot_chos(1:5).eq.'gupta')then
	   call diffvgp(pot_chos,n/3,xyz,v2t,v3t,v,g)
	elseif(pot_chos(1:5).eq.'tip3p')then
	   call diffvtip3p(pot_chos,n/2,nmono,xyz,v2t,v3t,v,g)
	elseif(pot_chos(1:5).eq.'tip4p')then
	   call diffvtip4p(pot_chos,n/2,nmono,xyz,v2t,v3t,v,g)
	elseif(pot_chos(1:5).eq.'tip5p')then
	   call diffvtip5p(pot_chos,n/2,nmono,xyz,v2t,v3t,v,g)
	else
	   write(*,*)' This potential does not exist'
	   write(*,*)'error:mini'
	   stop
	endif	
c       
	return
	end	








