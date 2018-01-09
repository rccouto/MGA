	subroutine sort_atoms_z(pot_chos,nclust,natoms,nmono,nat,pop_clus)
	implicit real*8(a-h,o-z)
	include "param.inc"
c
	integer        nat(10)   
	real*8         pop_clus(n_max_clus,n_max_atom)
	1             ,vec(n_max_atom)
	character*20   pot_chos	
c
	namo=natoms/nmono
	if(pot_chos(1:3).eq.'tip')then
	   nm=1
	   nmo=natoms/nmono
	else
	   nm=nmono
	   nmo=natoms
	endif
c
	do k=1,nclust,1
	   nni=1 
	   nf=0
	   do i=1,nm,1
	      ni=nni
	      nf=nf+nat(i)
	      l=0
	      do j=ni,nf,1 
		 l=l+1
		 vec(l)=pop_clus(k,j)
		 vec(l+nat(i))=pop_clus(k,j+nmo)
		 vec(l+2*nat(i))=pop_clus(k,j+2*nmo)
		 if(pot_chos(1:3).eq.'tip')then
		    vec(l+3*nmo)=pop_clus(k,j+3*nmo)
		    vec(l+4*nmo)=pop_clus(k,j+4*nmo)
		    vec(l+5*nmo)=pop_clus(k,j+5*nmo)
		 endif
	      enddo
	      call sort_vecz(pot_chos,nat(i),vec)
cxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
c		 write(38,*)'clust',i,'      energy'
c		 do j=1,namo
c		    write(38,*)j,vec(j),vec(j+namo)
c	1		 ,vec(j+2*namo)
c		    if(pot_chos(1:3).eq.'tip')then
c		       write(38,*)vec(j+3*namo),vec(j+4*namo)
c	1		    ,vec(j+5*namo)
c		    endif
c		 enddo
cxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
	      l=0
	      do j=ni,nf,1
		 l=l+1
		 pop_clus(k,j)=vec(l)
		 pop_clus(k,j+nmo)=vec(l+nat(i))
		 pop_clus(k,j+2*nmo)=vec(l+2*nat(i))
		 if(pot_chos(1:3).eq.'tip')then
		    pop_clus(k,j+3*nmo)=vec(l+3*nmo)
		    pop_clus(k,j+4*nmo)=vec(l+4*nmo)
		    pop_clus(k,j+5*nmo)=vec(l+5*nmo)
		 endif
	      enddo
	      nni=nf+1
	   enddo
	enddo
c	
	return
	end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine sort_vecz(pot_chos,m,vec)
	implicit real*8(a-h,o-z)
	include "param.inc"
c	
	real*8         vecat(0:n_max_atom),vec(n_max_atom)
	character*20   pot_chos	
c	
	do i=1,3*m,1
	   vecat(i)=+9d299
	enddo
	vecat(2*m)=-9d299
c
	do k=2*m+1,3*m,1
	   do i=2*m+1,3*m,1
	      if(vec(i).lt.vecat(k) .and. 
	1	   vec(i).gt.vecat(k-1))then
		 vecat(k-2*m)=vec(i-2*m)
		 vecat(k-m)=vec(i-m)
		 vecat(k)=vec(i)
		 if(pot_chos(1:3).eq.'tip')then
		    vecat(k+m)=vec(i+m)
		    vecat(k+2*m)=vec(i+2*m)
		    vecat(k+3*m)=vec(i+3*m)
		 endif
	      endif
	   enddo
	enddo
	do i=1,m,1
	   vec(i)=vecat(i)
	   vec(i+m)=vecat(i+m)
	   vec(i+2*m)=vecat(i+2*m)
	   if(pot_chos(1:3).eq.'tip')then
	      vec(i+3*m)=vecat(i+3*m)
	      vec(i+4*m)=vecat(i+4*m)
	      vec(i+5*m)=vecat(i+5*m)
	   endif
	enddo
c
	return
	end	
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine sort_atoms_x(pot_chos,nclust,natoms,nmono,nat,pop_clus)
	implicit real*8(a-h,o-z)
	include "param.inc"
c
	integer        nat(10)   
	real*8         pop_clus(n_max_clus,n_max_atom)
	1             ,vec(n_max_atom)
	character*20   pot_chos	
c
	namo=natoms/nmono
	if(pot_chos(1:3).eq.'tip')then
	   nm=1
	   nmo=natoms/nmono
	else
	   nm=nmono
	   nmo=natoms
	endif
c
	do k=1,nclust,1
	   nni=1 
	   nf=0
	   do i=1,nm,1
	      ni=nni
	      nf=nf+nat(i)
	      l=0
	      do j=ni,nf,1 
		 l=l+1
		 vec(l)=pop_clus(k,j)
		 vec(l+nat(i))=pop_clus(k,j+nmo)
		 vec(l+2*nat(i))=pop_clus(k,j+2*nmo)
		 if(pot_chos(1:3).eq.'tip')then
		    vec(l+3*nmo)=pop_clus(k,j+3*nmo)
		    vec(l+4*nmo)=pop_clus(k,j+4*nmo)
		    vec(l+5*nmo)=pop_clus(k,j+5*nmo)
		 endif
	      enddo
	      call sort_vecx(pot_chos,nat(i),vec)
cxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
c		 write(38,*)'clust',i,'      energy'
c		 do j=1,namo
c		    write(38,*)j,vec(j),vec(j+namo)
c	1		 ,vec(j+2*namo)
c		    if(pot_chos(1:3).eq.'tip')then
c		       write(38,*)vec(j+3*namo),vec(j+4*namo)
c	1		    ,vec(j+5*namo)
c		    endif
c		 enddo
cxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
	      l=0
	      do j=ni,nf,1
		 l=l+1
		 pop_clus(k,j)=vec(l)
		 pop_clus(k,j+nmo)=vec(l+nat(i))
		 pop_clus(k,j+2*nmo)=vec(l+2*nat(i))
		 if(pot_chos(1:3).eq.'tip')then
		    pop_clus(k,j+3*nmo)=vec(l+3*nmo)
		    pop_clus(k,j+4*nmo)=vec(l+4*nmo)
		    pop_clus(k,j+5*nmo)=vec(l+5*nmo)
		 endif
	      enddo
	      nni=nf+1
	   enddo
	enddo
c	
	return
	end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine sort_vecx(pot_chos,m,vec)
	implicit real*8(a-h,o-z)
	include "param.inc"
c	
	real*8         vecat(0:n_max_atom),vec(n_max_atom)
	character*20   pot_chos	
c	
	vecat(0)=-9d299
	do i=1,m,1
	vecat(i)=+9d299
	enddo
c
	do k=1,m,1
	   do i=1,m,1
	      if(vec(i).lt.vecat(k).and.vec(i).gt.vecat(k-1))then
		 vecat(k)=vec(i)
		 vecat(k+m)=vec(i+m)
		 vecat(k+2*m)=vec(i+2*m)
		 if(pot_chos(1:3).eq.'tip')then
		    vecat(k+3*m)=vec(i+3*m)
		    vecat(k+4*m)=vec(i+4*m)
		    vecat(k+5*m)=vec(i+5*m)
		 endif
	      endif
	   enddo
	enddo
	do i=1,m,1
	   vec(i)=vecat(i)
	   vec(i+m)=vecat(i+m)
	   vec(i+2*m)=vecat(i+2*m)
	   if(pot_chos(1:3).eq.'tip')then
	      vec(i+3*m)=vecat(i+3*m)
	      vec(i+4*m)=vecat(i+4*m)
	      vec(i+5*m)=vecat(i+5*m)
	   endif
	enddo
c
	return
	end	
c	

	
	

	
	
