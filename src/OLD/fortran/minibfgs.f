	subroutine minibfgs(pot_chos,nclust,natoms,nat,nmono,vm2,vm3
	1    ,energ,pop)
	implicit real*8(a-h,o-z)
	include "param.inc"
	common/cpot/qo,qh,AA,BB
        common/pp/ec
c
        integer        nbd(n_max_atom),isave(44),iwa(3*nmax)
c       
        real*8         pop(n_max_clus,n_max_atom),l(n_max_atom)
	1             ,wa(2*mmax*3*nmax+12*nmax+12*mmax*mmax+12*mmax)
	2             ,xyz(n_max_atom),energ(n_max_clus)
	3             ,vm2(n_max_clus),vm3(n_max_clus)
	4             ,g(n_max_atom),u(n_max_atom),dsave(29)       
	character*20   pot_chos
        character*60   task,csave
        logical        lsave(4)
c       
c- minimisation options
	m=6
	factr=1d+7
	pgtol=1d-5
	iprint=-1
	namo=natoms/nmono
c
c- parallelize the structure minimisations
C$PAR DOALL
C$PAR& PRIVATE(g,xyz,v,v2t,v3t)
C$PAR& PRIVATE(nbd,u,l,wa,iwa,csave,lsave,isave,dsave)
C$PAR& SHARED(m,factr,pgtol,namo,pot_chos,nclust,natoms,nat,nmono)
C$PAR& SHARED(pop,energ,vm2,vm3)
C$PAR& MAXCPUS(nclust)
C$PAR& READONLY(m,factr,pgtol,namo,pot_chos,nclust,natoms,nat,nmono)
C$PAR& SAVELAST 
C$PAR& SCHEDTYPE(SELF[(1)]) 
	do i1p=1,nclust,1
	   i=i1p
	   task='START'
	   do k=1,6*natoms,1
	      nbd(k)=0
	      u(k)=0
	      l(k)=1000
	      g(k)=0
	      xyz(k)=pop(i,k)
	   enddo

	   write(*,*)'Indiv =', i1p

cxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
c	do i=1,noff*nclust
c	   write(60,*)'clust',i,'      energy',energ(i)
c	   do jj=1,namo
c	      write(60,*)jj,pop(i,jj),pop(i,jj+namo)
c	1	   ,pop(i,jj+2*namo)
c	      write(60,*)jj,xyz(jj),xyz(jj+namo),xyz(jj+2*namo)
c	      if(pot_chos(1:3).eq.'tip')then
c		 write(60,*)pop(i,jj+3*namo),pop(i,jj+4*namo)
c	1	      ,pop(i,jj+5*namo)
c	      endif
c	   enddo
c	enddo
c	close(60)
cxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
c       
	   if(pot_chos(1:5).eq.'morse'
	1	.or. pot_chos(1:9).eq.'tlr_morse'  
	2	.or. pot_chos(1:5).eq.'mm_2b')then
	      do while(task(1:5).eq.'START' .or. task(1:2).eq.'FG'
	1	   .or. task(1:5).eq.'NEW_X')
		 call setulb(3*natoms,m,xyz,l,u,nbd,v,g,factr
	1	      ,pgtol,wa,iwa,task,iprint,csave
	2	      ,lsave,isave,dsave)
		 if(task(1:2).eq.'FG')call diffv2(pot_chos,natoms
	1	      ,xyz,v,g) 
	      enddo  
	   elseif(pot_chos(1:5).eq.'mm_cp')then
	      do while(task(1:5).eq.'START' 
	1	   .or. task(1:2).eq.'FG'
	2	   .or. task(1:5).eq.'NEW_X')

		 write(*,*)'inside mm_cp',task

		 call setulb(3*natoms,m,xyz,l,u,nbd,v,g,factr
	1	      ,pgtol,wa,iwa,task,iprint,csave
	2	      ,lsave,isave,dsave)

		 write(*,*)'inside after setulb mm_cp',task

		 if(task(1:2).eq.'FG')call diffvmm(pot_chos,natoms
	1	      ,xyz,v2t,v3t,v,g)

		 write(*,*)'v=',v,'   g(j)=', (g(j), j=1,36,1)
		 do j=1,12,1
		    write(*,*)xyz(j), xyz(j+12), xyz(j+24)
		 enddo
c		 read(*,*)
	      enddo
	   elseif(pot_chos(1:5).eq.'gupta')then
	      do while(task(1:5).eq.'START' .or. task(1:2).eq.'FG'
	1	   .or. task(1:5).eq.'NEW_X')
		 call setulb(3*natoms,m,xyz,l,u,nbd,v,g,factr
	1	      ,pgtol,wa,iwa,task,iprint,csave
	2	      ,lsave,isave,dsave)
		 if(task(1:2).eq.'FG')call diffvgp(pot_chos,natoms
	1	      ,nat,xyz,v2t,v3t,v,g)
	      enddo
	   elseif(pot_chos(1:5).eq.'tip3p')then
	      do while(task(1:5).eq.'START' .or. task(1:2).eq.'FG'
	1	   .or. task(1:5).eq.'NEW_X')	
		 call setulb(6*namo,m,xyz,l,u,nbd,v,g,factr
	1	      ,pgtol,wa,iwa,task,iprint,csave
	1	      ,lsave,isave,dsave)
		 if(task(1:2).eq.'FG')call diffvtip3p(pot_chos,natoms
	1	      ,nmono,xyz,v2t,v3t,v,g)
	      enddo
	   elseif(pot_chos(1:5).eq.'tip4p')then
	      do while(task(1:5).eq.'START' .or. task(1:2).eq.'FG'
	1	   .or. task(1:5).eq.'NEW_X')	
		 call setulb(6*namo,m,xyz,l,u,nbd,v,g,factr
	1	      ,pgtol,wa,iwa,task,iprint,csave
	2	      ,lsave,isave,dsave)
		 if(task(1:2).eq.'FG')call diffvtip4p(pot_chos,natoms
	1	      ,nmono,xyz,v2t,v3t,v,g)
	      enddo
	   elseif(pot_chos(1:5).eq.'tip5p')then
	      do while(task(1:5).eq.'START' .or. task(1:2).eq.'FG'
	1	   .or. task(1:5).eq.'NEW_X')	
		 call setulb(6*namo,m,xyz,l,u,nbd,v,g,factr
	1	      ,pgtol,wa,iwa,task,iprint,csave
	1	      ,lsave,isave,dsave)
		 if(task(1:2).eq.'FG')call diffvtip5p(pot_chos,natoms
	1	      ,nmono,xyz,v2t,v3t,v,g)
	      enddo



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
	   do k=1,6*natoms,1
	      pop(i,k)=xyz(k)
	   enddo
	enddo
	close(60)
c
	return
	end













