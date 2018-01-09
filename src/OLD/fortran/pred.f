	subroutine predator(pot_chos,opt_pred,opt_min,nclust,natoms,nat
	1                  ,iclus,nseed,naux1,ngen,nmono,noff,enr_min
	2                  ,enr_bef,enr_1,enr_nxt_m,vm2,vm3,energ
	3                  ,pop_pred,task1)
	implicit real*8(a-h,o-z)
	include "param.inc"
	include "common.inc"
c	
	integer        nat(10)
	real*8         pop_pred(n_max_clus,n_max_atom)
	1             ,energ(nclust),vm2(n_max_clus)
	2             ,vm3(n_max_clus)
	character*20   pot_chos,opt_pred,opt_min
	logical        task1
c
	task1=.false.
	n_killnext=0
	n_aux_pred=0
	nxt_clus_aux=0
	nxt_clus=0
c 
	n_kill2=0
	do i=1,nclust,1
	   n_kill=0
	   if(energ(i).lt.enr_1)then
	      enr_min=-9d299
	      enr_bef=-9d299
	      enr_1=-9d299
	      iclus=1
	      task1=.true.
c	      close(2)
c	      close(3)
	   elseif(energ(i).le.enr_min)then
	      n_aux_pred=1	
	      if(opt_pred(1:5).eq.'subst')then
		 if(ngen.eq.0)goto 111
 222		 nxt_clus=(noff+1)*nclust*ran3(nseed) !robson 05/04
		 if(energ(nxt_clus).lt.enr_min)goto 222       
		 energ(i)=energ(nxt_clus)
		 do j=1,3*natoms,1 
		    pop_pred(i,j)=pop_pred(nxt_clus,j)
		 enddo
c       
	      elseif(opt_pred(1:4).eq.'gnew')then
 111		 continue
		 n_kill=n_kill+1
		 if(n_kill.gt.1000)goto 333
		 call gener_init_pop(pot_chos,1,natoms,nmono
	1	      ,nseed,naux1,pop_pred(i,1))
		 if(opt_min(1:4).eq.'bfgs')then
		    call minibfgs(pot_chos,1,natoms,nat,nmono
	1		 ,vm2(i),vm3(i),energ(i),pop_pred(i,1))
		 elseif(opt_min(1:2).eq.'cg')then
		    call minicg(pot_chos,1,natoms,nmono
	1		 ,vm2(i),vm3(i),energ(i),pop_pred(i,1))
		 else
		    write(*,*)'minimizator not found'
		    stop
		 endif
		 if(energ(i).lt.enr_min)goto 111
	      elseif(opt_pred(1:4).eq.'next')then
		 if(ngen.eq.0)goto 111
		 n_kill2=n_kill2+1
		 if(n_kill.ge.2*nclust)goto 111
		 energ(i)=0d0 
	      else
		 write(*,*)'Warning!!! predator error Warning!!!'
		 write(*,*)'error: p1'
		 stop
 333		 continue
		 write(*,*)'predator could not find another minimum!!!'
		 stop
	      endif
	   endif
c       
	enddo
c	
	return
	end
	


