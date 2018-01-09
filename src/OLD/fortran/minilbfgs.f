	subroutine minilbfgs(pot_chos,nclust,natoms,nat,nmono,vm2,vm3
	1    ,energ,dhes,pop)
	implicit real*8(a-h,o-z)
	include "param.inc"
	external lb2
	common/cpot/qo,qh,AA,BB
        common/pp/ec
c
        real*8         pop(n_max_clus,n_max_atom),hes(n_max_atom)
	1             ,w(3*n_max_atom*(2*10+1)+2*10)
	2             ,xyz(n_max_atom),energ(n_max_clus)
	3             ,vm2(n_max_clus),vm3(n_max_clus)
	4             ,dhes(n_max_clus,n_max_atom),g(n_max_atom)
	character*20   pot_chos
        logical        diagco
c       
c- minimisation options 
	m=7
	diagco='true'
	iprint=-1
	eps=1d-8
	xtol=2.2204460492503D-16
	namo=natoms/nmono
c
	do i=1,nclust,1
	   task='START'
	   do k=1,6*natoms,1
	      g(k)=0
	      xyz(k)=pop(i,k)
	   enddo
c       
	   if(pot_chos(1:5).eq.'morse'
	1	.or. pot_chos(1:9).eq.'tlr_morse'  
	2	.or. pot_chos(1:5).eq.'mm_2b')then
	      do while(iflag.ne.2)
		 write(*,*)'i,iflag'i,iflag
		 call lbfgs(3*natoms,m,xyz,v,g,diagco,diag
	1	      ,iprint,eps,xtol,w,iflag)
		 if(iflag.eq.0)call diffv2(pot_chos,natoms
	1	      ,xyz,v,g) 
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
	      dhes(i,k)=hes(k)
	   enddo
	enddo
c
	return
	end













