	subroutine mating(pot_chos,opt_mat,ditr_chos,natoms,nclust,mmol
     &                   ,nseed,nseed1,noff,nmono,nat,fit,pop_clus
     &                   ,pop_off)
	implicit real*8(a-h,o-z)
	include "param.inc"	
c
	integer        nat(10)	
	real*8         pop_clus(n_max_clus,n_max_atom)
	1             ,pop_off(n_max_clus,n_max_atom)
	1             ,pop(n_max_clus,n_max_atom)
	2             ,pop_chd(4,n_max_atom),fit(n_max_clus)
	3             ,poff(n_max_atom)
	character*20   opt_mat,ditr_chos,pot_chos
	external       ran3
c
	if(pot_chos(1:3).eq.'tip')then
	   nmt=1
	   nmo=natoms/nmono
	else
	   nmt=nmono
	   nmo=natoms
	endif

	call matr1_matr2(nclust,3*natoms,pop_clus,pop)

	do i=1,nclust,1
	   call matr_vec(i,3*nmo,poff,pop)
	   call rotate(nmo,nseed,poff,1,nmo,1)
c	   call rotate1(nmo,poff,nseed1)
	   call vec_matr(i,3*nmo,poff,pop)
	enddo
c	
	td=ran3(nseed1)
	if(td.le.5d-1)then
	   call sort_atoms_z(pot_chos,nclust,natoms,nmono,nat,pop)
	else
	   call sort_atoms_x(pot_chos,nclust,natoms,nmono,nat,pop) 
	endif
c
cxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
c	      do i=1,nclust
c		 write(58,*)'clust',i,'      energy',energ
c		 do j=1,nmo,1
c		    write(58,*)j,pop(i,j),pop(i,j+nmo)
c	1		 ,pop(i,j+2*nmo)
c		    if(pot_chos(1:3).eq.'tip')then
c		       write(58,*)pop(i,j+3*nmo),pop(i,j+4*nmo)
c	1		    ,pop(i,j+5*nmo)
c		    endif
c		 enddo
c	      enddo
c	      close(58)
cxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
c
	do 1000 j=1,nclust,1
 11	   iclus_chos1=int((nclust-1)*ran3(nseed1))+1
	   if(ran3(nseed).gt.fit(iclus_chos1))goto 11 !robson 05/04
 21	   iclus_chos2=int((nclust-1)*ran3(nseed1))+1
	   if(iclus_chos1.eq.iclus_chos2 .or. 
	1	ran3(nseed).gt.fit(iclus_chos2))goto 21 !robson 05/04
	   if(opt_mat(1:8).eq.'one_part'
	1	.or. opt_mat(1:8).eq.'many_off')then
	      nms=mmol
	      naux=int(nms*sqrt(ran3(nseed)))  !robson 05/04
	      if(naux.lt.2)naux=2
	      if(naux.ge.nms)naux=nms-2
	      nii=1
	      nf=0
	      do 100 k=1,nmt,1
		 ni=nii
		 nf=nf+nat(k)
		 do 110 i=ni,naux,1
		    pop_chd(1,i)=pop(iclus_chos1,i)
		    pop_chd(1,i+nmo)=pop(iclus_chos1,i+nmo)
		    pop_chd(1,i+2*nmo)=pop(iclus_chos1,i+2*nmo)
		    pop_chd(2,i)=pop(iclus_chos2,i)
		    pop_chd(2,i+nmo)=pop(iclus_chos2,i+nmo)
		    pop_chd(2,i+2*nmo)=pop(iclus_chos2,i+2*nmo)
		    if(pot_chos(1:3).eq.'tip')then
		       pop_chd(1,i+3*nmo)=pop(iclus_chos1,i+3*nmo)
		       pop_chd(1,i+4*nmo)=pop(iclus_chos1,i+4*nmo)
		       pop_chd(1,i+5*nmo)=pop(iclus_chos1,i+5*nmo)
		       pop_chd(2,i+3*nmo)=pop(iclus_chos2,i+3*nmo)
		       pop_chd(2,i+4*nmo)=pop(iclus_chos2,i+4*nmo)
		       pop_chd(2,i+5*nmo)=pop(iclus_chos2,i+5*nmo)
		    endif
 110		 continue
		 do 120 i=naux+1,nf,1
		    pop_chd(1,i)=pop(iclus_chos2,i)
		    pop_chd(1,i+nmo)=pop(iclus_chos2,i+nmo)
		    pop_chd(1,i+2*nmo)=pop(iclus_chos2,i+2*nmo)
		    pop_chd(2,i)=pop(iclus_chos1,i)
		    pop_chd(2,i+nmo)=pop(iclus_chos1,i+nmo)
		    pop_chd(2,i+2*nmo)=pop(iclus_chos1,i+2*nmo)
		    if(pot_chos(1:3).eq.'tip')then
		       pop_chd(1,i+3*nmo)=pop(iclus_chos2,i+3*nmo)
		       pop_chd(1,i+4*nmo)=pop(iclus_chos2,i+4*nmo)
		       pop_chd(1,i+5*nmo)=pop(iclus_chos2,i+5*nmo)
		       pop_chd(2,i+3*nmo)=pop(iclus_chos1,i+3*nmo)
		       pop_chd(2,i+4*nmo)=pop(iclus_chos1,i+4*nmo)
		       pop_chd(2,i+5*nmo)=pop(iclus_chos1,i+5*nmo)
		    endif
 120		 continue
		 naux=naux+nat(k)
		 nii=nf+1
 100	      continue
	      n_opt=1
c
	      if(opt_mat(1:8).eq.'many_off')goto 31
	   elseif(opt_mat(1:9).eq.'two_parts')then
 31	      continue
	      nms=mmol
	      naux1=int(nms*5d-1*sqrt(ran3(nseed))) !robson 05/04
	      naux2=int(nms*(5d-1*sqrt(ran3(nseed))+5d-1)) !robson 05/04
	      if(naux1.le.1)naux1=2
	      if(naux2.le.int(nms*5d-1+1))naux1=int(nms*5d-1+2)
	      nii=1
	      nf=0
	      do 200  k=1,nmt,1
		 ni=nii
		 nf=nf+nat(k)
		 do 210 i=ni,naux1,1
		    pop_chd(3,i)=pop(iclus_chos1,i)
		    pop_chd(3,i+nmo)=pop(iclus_chos1,i+nmo)
		    pop_chd(3,i+2*nmo)=pop(iclus_chos1,i+2*nmo)
		    pop_chd(4,i)=pop(iclus_chos2,i)
		    pop_chd(4,i+nmo)=pop(iclus_chos2,i+nmo)
		    pop_chd(4,i+2*nmo)=pop(iclus_chos2,i+2*nmo)
		    if(pot_chos(1:3).eq.'tip')then
		       pop_chd(3,i+3*nmo)=pop(iclus_chos1,i+3*nmo)
		       pop_chd(3,i+4*nmo)=pop(iclus_chos1,i+4*nmo)
		       pop_chd(3,i+5*nmo)=pop(iclus_chos1,i+5*nmo)
		       pop_chd(4,i+3*nmo)=pop(iclus_chos2,i+3*nmo)
		       pop_chd(4,i+4*nmo)=pop(iclus_chos2,i+4*nmo)
		       pop_chd(4,i+5*nmo)=pop(iclus_chos2,i+5*nmo)
		    endif
 210		 continue
		 do 220 i=naux1+1,naux2,1
		    pop_chd(3,i)=pop(iclus_chos2,i)
		    pop_chd(3,i+nmo)=pop(iclus_chos2,i+nmo)
		    pop_chd(3,i+2*nmo)=pop(iclus_chos2,i+2*nmo)
		    pop_chd(4,i)=pop(iclus_chos1,i)
		    pop_chd(4,i+nmo)=pop(iclus_chos1,i+nmo)
		    pop_chd(4,i+2*nmo)=pop(iclus_chos1,i+2*nmo)
		    if(pot_chos(1:3).eq.'tip')then
		       pop_chd(3,i+3*nmo)=pop(iclus_chos2,i+3*nmo)
		       pop_chd(3,i+4*nmo)=pop(iclus_chos2,i+4*nmo)
		       pop_chd(3,i+5*nmo)=pop(iclus_chos2,i+5*nmo)
		       pop_chd(4,i+3*nmo)=pop(iclus_chos1,i+3*nmo)
		       pop_chd(4,i+4*nmo)=pop(iclus_chos1,i+4*nmo)
		       pop_chd(4,i+5*nmo)=pop(iclus_chos1,i+5*nmo)
		    endif
 220		 continue
		 do 230 i=naux2+1,nf,1
		    pop_chd(3,i)=pop(iclus_chos1,i)
		    pop_chd(3,i+nmo)=pop(iclus_chos1,i+nmo)
		    pop_chd(3,i+2*nmo)=pop(iclus_chos1,i+2*nmo)
		    pop_chd(4,i)=pop(iclus_chos2,i)
		    pop_chd(4,i+nmo)=pop(iclus_chos2,i+nmo)
		    pop_chd(4,i+2*nmo)=pop(iclus_chos2,i+2*nmo)
		    if(pot_chos(1:3).eq.'tip')then
		       pop_chd(3,i+3*nmo)=pop(iclus_chos1,i+3*nmo)
		       pop_chd(3,i+4*nmo)=pop(iclus_chos1,i+4*nmo)
		       pop_chd(3,i+5*nmo)=pop(iclus_chos1,i+5*nmo)
		       pop_chd(4,i+3*nmo)=pop(iclus_chos2,i+3*nmo)
		       pop_chd(4,i+4*nmo)=pop(iclus_chos2,i+4*nmo)
		       pop_chd(4,i+5*nmo)=pop(iclus_chos2,i+5*nmo)
		    endif
 230		 continue
		 naux1=naux1+nat(k)
		 naux2=naux2+nat(k)
		 nii=nf+1
 200	      continue
	      n_opt=3
	   else
	      write(*,*)'!! Caution mating population error !!'
	      write(*,*)'Error: Ma'
	   endif 
c
c generat 4 offsprings --------------------------------------	
	   if(opt_mat(1:11).eq.'many_off=4d')then
	      noff=4
	      do 600 k=1,4,1
		 do 610 i=1,3*natoms,1
		    pop_off((j-1)*4+k,i)=pop_chd(k,i)
 610		 continue
 600	      continue
c generat 2 offsprings of different matings -----------------
	   elseif(opt_mat(1:12).eq.'many_off=2dd')then
	      noff=2
	      ki=0
	      do 650 k=1,3,2
		 ki=ki+1
		 do 660 i=1,3*natoms,1
		    pop_off((j-1)*2+ki,i)=pop_chd(k,i)
 660		 continue
 650	      continue
c generat 2 offsprings of same matings only one peace -------
	   elseif(opt_mat(1:12).eq.'many_off=2so')then
	      noff=2 
	      do 700 k=1,2,1
		 do 710 i=1,3*natoms,1
		    pop_off((j-1)*2+k,i)=pop_chd(k,i)
 710		 continue
 700	      continue 
c generat 2 offsprings of same matings two peace ------------
	   elseif(opt_mat(1:12).eq.'many_off=2st')then
	      noff=2
	      do 750 k=3,4,1
		 do 760 i=1,3*natoms,1
		    pop_off((j-1)*2+k-2,i)=pop_chd(k,i)
 760		 continue
 750	      continue
c generat 1 offspring ---------------------------------------	
	   else 
	      noff=1
	      do 900 i=1,3*natoms,1
		 pop_off(j,i)=pop_chd(n_opt,i)
 900	      continue
	   endif
 1000	continue
c	
 3000	continue
	return
	end





