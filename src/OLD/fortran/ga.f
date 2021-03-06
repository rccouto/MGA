CcCcCc==CcCcCc==CcCcCc==CcCcCc==CcCcCc==CcCcCc==CcCcCc==CcCcCc==CcCcCc==CcCcCcc
cc/-\||/-\|/-\||/-\|/-\||/-\|/-\||/-\|/-\||/-\|/-\||/-\|/-\||/-\|/-\||/-\|/-\cc
cc\-/||\-/|\-/||\-/|\-/||\-/|\-/||\-/|\-/||\-/|\-/||\-/|\-/||\-/|\-/||\-/|\-/cc
CcCcCc==CcCcCc==CcCcCc==CcCcCc==CcCcCc==CcCcCc==CcCcCc==CcCcCc==CcCcCc==CcCcCcc
cxxxxxxxx1xxxxxxxxx2xxxxxxxxx3xxxxxxxxx4xxxxxxxxx5xxxxxxxxx6xxxxxxxxx7xxxxxxxxx
c234567890123456789012345678901234567890123456789012345678901234567890123456789
      PROGRAM mga
      implicit real*8(a-h,o-z)
      include "param.inc"
	include "common.inc"
c       
	integer        nat(10)
	real*8         pop_clus(n_max_clus,n_max_atom)
	1    ,pop_off(n_max_clus,n_max_atom)
	2    ,pop_hist(n_max_clus,n_max_atom)
	3    ,energ(n_max_clus),energ_bet(n_max_clus)
	4    ,vm2(n_max_clus),vm2_bet(n_max_clus)
	5    ,vm3(n_max_clus),vm3_bet(n_max_clus)
	6    ,enert(n_max_clus),vm2t(n_max_clus)
	7    ,vm3t(n_max_clus)
	8    ,fit(n_max_clus),vec(n_max_clus)
	character*20   pot_chos,opt_mat,opt_pred,opt_ini,opt_ga
	1    ,opt_mut,ditr_chos,inputfile,min_chos
	2    ,opt_min,opt_out
	character*30   filename
	logical        task1
c
c>>     <<<>>><<<>>><<<>>><<<>>><<<>>><<<>>><<<>>><<<>>><<<>>><<<>>><<<>>><<<>>><<<
c       - Read input parameters >><<<>>><<<>>><<<>>><<<>>><<<>>><<<>>><<<>>><<<>>><<<
	write(*,*)'read_in'
	call  read_in(nclust,nat,nmono,nseed,nbet_clust,nterm_aft_rep
	1    ,nmaxgen,per_mut,enr_nxt_m,natoms,namo,filename
	2    ,nfilename,opt_ga,ditr_chos,opt_pred,pot_chos
	3    ,opt_mat,opt_mut,opt_ini,inputfile,ifilelen
	4    ,min_chos,opt_min,opt_out)
	write(*,*)'read_in2'
c
	nseed1=nseed/2
	nm4=natoms
	if(pot_chos(1:3).eq.'tip')then
	   nm4=nat(1)
	   namo=natoms/nmono
	   nms=natoms/nmono
	else
	   nm4=natoms
	   namo=natoms
	   nms=nat(1)
	   if(nms.gt.nat(2) .and. nat(2).ge.4)nms=nat(2)
	   if(nms.gt.nat(3) .and. nat(3).ge.4)nms=nat(3) 
	   if(nms.gt.nat(4) .and. nat(4).ge.4)nms=nat(4)
	   if(nms.gt.nat(5) .and. nat(5).ge.4)nms=nat(5)
	endif
c
c       - read potential constats <<<>>><<<>>><<<>>><<<>>><<<>>><<<>>><<<>>><<<>>><<<
	write(*,*)'read_pot_const'
	call read_pot_const(pot_chos,inputfile,ifilelen,natoms,nmono,nat)
	write(*,*)'read_pot_const2'
c       
c       - set correct value for the variable ><<<>>><<<>>><<<>>><<<>>><<<>>><<<>>><<<
	task1=.false.
	naux1=0
	enr_min=-9d299
	enr_bef=-9d299
        enr_1=-9d299
c       
c       - generet struc.hlp <<<>>><<<>>><<<>>><<<>>><<<>>><<<>>><<<>>><<<>>><<<>>><<<
	do n=1,n_max_clus,1
	   energ(n)=0d0
	   energ_bet(n)=0d0
	   vm2(n)=0d0
	   vm2_bet(n)=0d0
	   vm3(n)=0d0
	   vm3_bet(n)=0d0
	   enert(n)=0d0
	   vm2t(n)=0d0
	   vm3t(n)=0d0
	   fit(n)=0d0
	   vec(n)=0d0
	   do m=1,n_max_atom,1
	      pop_clus(n,m)=0d0
	      pop_off(n,m)=0d0
	      pop_hist(n,m)=0d0
	   enddo
	enddo
	iout=0
c
c       - start search of cluster <<<>>><<<>>><<<>>><<<>>><<<>>><<<>>><<<>>><<<>>><<<
	do iclus=1,nbet_clust,1
	   iout=iout+1
c       
c       - set correct value for the variable ><<<>>><<<>>><<<>>><<<>>><<<>>><<<>>><<<
	   n_count=0
	   noff=2
	   ngen=-1
	   ncont_aux_ga=0
	   ene_avbf=0d0
c       
c       - set zero value for matrix and vectors <>>><<<>>><<<>>><<<>>><<<>>><<<>>><<<
	   do n=1,n_max_clus,1
	      energ(n)=0d0
	      energ_bet(n)=0d0
	      vm2(n)=0d0
	      vm3(n)=0d0
	      do m=1,n_max_atom,1
		 pop_clus(n,m)=0d0
		 pop_off(n,m)=0d0
	      enddo
	   enddo
c
c       - Generate initial population >><<<>>><<<>>><<<>>><<<>>><<<>>><<<>>><<<>>><<<
	   write(*,*)'opt_ini'
	   if(opt_ini(1:6).eq.'random' .and. iclus.eq.1)then
	      call gener_init_pop(pot_chos,noff*nclust,natoms,nmono
	1	   ,nseed,naux1,pop_off)
	      ncm=nclust
	   elseif(opt_ini(1:4).eq.'file' .and. iclus.eq.1)then
	      call read_pop_ini(pot_chos,filename,nfilename,2*nclust
	1	   ,natoms,namo,energ,vm2,vm3,pop_off)
	   elseif(iclus.ge.2)then
	      ngen=-1
	      noff=2
              ncm=nclust
	      n_count=5
c       
	      do ip=1,nclust+n_count,1
		 energ(ip)=enert(ip)
	      enddo
	      call matr1_matr2(nclust+n_count,3*natoms,pop_hist(1,1)
	1	   ,pop_off(1,1)) 
	      call gener_init_pop(pot_chos,nclust,natoms,nmono
	1	   ,nseed,naux1,pop_off(ncm+n_count+1,1))
c       
	   else
	      write(*,*)'Error reading intial population'
	      write(*,*)'error: GA'
	   endif
	   write(*,*)'opt_ini2'
c
 1111	   continue
	   ngen=ngen+1
	   ncm=nclust
	   n_cl_tot=noff*ncm+n_count
c       
c       - Minimiser the random clusters <<<>>><<<>>><<<>>><<<>>><<<>>><<<>>><<<>>><<<
	   write(*,*)'opt_min'
	   if(opt_min(1:4).eq.'bfgs')then
	      call minibfgs(pot_chos,n_cl_tot,natoms,nat,nmono,vm2
	1	   ,vm3,energ,pop_off)
	   elseif(opt_min(1:2).eq.'cg')then
	      call minicg(pot_chos,n_cl_tot,natoms,nmono,vm2,vm3
	1	   ,energ,pop_off)
	   else
	      write(*,*)'minimizator not found'
	      stop
	   endif
	   write(*,*)'opt_min2'
c
c       - Call routine predator >><<<>>><<<>>><<<>>><<<>>><<<>>><<<>>><<<>>><<<>>><<<
	   write(*,*)'predator'
	   call predator(pot_chos,opt_pred,opt_min,n_cl_tot,natoms
	1	,nat,iclus,nseed,naux1,ngen,nmono,noff,enr_min
	2	,enr_bef,enr_1,enr_nxt_m,vm2,vm3,energ
	3	,pop_off,task1)
	   write(*,*)'predator2'
c
c       - Kill offsprings with same energy >>><<<>>><<<>>><<<>>><<<>>><<<>>><<<>>><<<
	   if(ngen.gt.0)then
	      call sort_matr(pot_chos,n_cl_tot,natoms
	1	   ,vm2,vm3,energ,pop_off)
	      enr_exs=energ(1)
	      do i=2,n_cl_tot,1
		 if(abs(energ(i)-enr_exs).le.enr_nxt_m)then
		    energ(i)=0d0 
c       else
c       enr_exs=energ(i)
		 endif
		 enr_exs=energ(i)
		 k=i
		 if(energ(k).eq.0d0)then
 1214		    continue
		    k=k-1
		    if(energ(k).eq.0d0)goto 1214
		    enr_exs=energ(k)
		 endif
	      enddo
	      call sort_matr(pot_chos,n_cl_tot,natoms
	1	   ,vm2,vm3,energ,pop_off)
	   endif
c
c       - Rewrite matrix offsprings to matrix clusters >>><<<>>><<<>>><<<>>><<<>>><<<
	   nnc=nclust
	   ncc=0.8*nclust+nclust
	   if(ngen.eq.0)ncc=3*nclust
	   do i=1,n_cl_tot,1
	      vm2_bet(i+nnc)=vm2(i)
	      vm3_bet(i+nnc)=vm3(i)
	      energ_bet(i+nnc)=energ(i)
	   enddo
	   call matr1_matr2(n_cl_tot,3*natoms,pop_off
	1	,pop_clus(nnc+1,1))	
c       
c       - Sort matrix in decrease order of energy >><<<>>><<<>>><<<>>><<<>>><<<>>><<<
	   write(*,*)'sort_matr'
	   call sort_matr(pot_chos,ncc,natoms
	1	,vm2_bet,vm3_bet,energ_bet
	2	,pop_clus)
	   write(*,*)'sort_matr'
c       
c       - Center clusters in the center of mass <>>><<<>>><<<>>><<<>>><<<>>><<<>>><<<
	   do i=1,nclust,1
	      call matr_vec(i,3*natoms,vec,pop_clus)
	      call centre(nm4,re,vec) 
	      call vec_matr(i,3*natoms,vec,pop_clus)
	   enddo
 1213	   continue
c       
c       - Save initial population for meta-stable states ><<<>>><<<>>><<<>>><<<>>><<<
	   call pimss(pot_chos,opt_ga,opt_min,nclust,natoms,nmono
	1	,n_cl_tot,naux1,enr_min,enr_bef,enr_1
	2	,enr_nxt_m,vm2_bet,vm3_bet,vm2t,vm3t
	3	,energ_bet,enert,pop_clus,pop_hist)	      
c       
c       - Print data in screen >>><<<>>><<<>>><<<>>><<<>>><<<>>><<<>>><<<>>><<<>>><<<
	   en_aver=sum_vec(nclust,energ_bet)/nclust	      
	   if(opt_out(1:3).eq.'yes')call print_1(filename
	1	,nfilename,nclust,iout,ngen
	2	,en_aver,energ_bet)
c       
c       - Calculate the fitness for each cluster >>><<<>>><<<>>><<<>>><<<>>><<<>>><<<
	   call fitness(nclust,fit,energ_bet) 
c       
c       - Make mate of clusters >><<<>>><<<>>><<<>>><<<>>><<<>>><<<>>><<<>>><<<>>><<<
	   call mating(pot_chos,opt_mat,ditr_chos,natoms,ncm,nms
	1	,nseed,nseed1,noff,nmono,nat,fit,pop_clus
	2	,pop_off)
cxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
	   do i=1,2*nclust,1
	      write(12,*)'clust',i,'      energy',energ(i)
	      do j=1,natoms
		 write(12,*)j,pop_off(i,j),pop_off(i,j+nmono)
	1	      ,pop_off(i,j+2*nmono)
		 if(pot_chos(1:3).eq.'tip')then
		    write(12,*)pop_off(i,j+3*nmono),pop_off(i,j+4*nmono)
	1		 ,pop_off(i,j+5*nmono)
		 endif
	      enddo
	   enddo
	   close(12)
cxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx   	      
c       
c       - Rewrite matrix offsprings to matrix clusters >>><<<>>><<<>>><<<>>><<<>>><<<
	   call matr1_matr2(nclust,3*natoms,pop_off
	1	,pop_clus(nclust+1,1))
c       
c       - Mutate clusters >><<<>>><<<>>><<<>>><<<>>><<<>>><<<>>><<<>>><<<>>><<<>>><<<
	   call muting(pot_chos,opt_mut,per_mut,nseed,nms,naux1,nclust
	1	,natoms,nmono,n_count,noff,nat,pop_clus)
c       
c       - Rewrite matrix clusters to matrix offsprings >>><<<>>><<<>>><<<>>><<<>>><<<
	   call matr1_matr2(n_count,3*natoms
	1	,pop_clus(2*nclust+1,1)
	2	,pop_off(noff*nclust+1,1))
cxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
	   do i=1,3*nclust,1
	      write(13,*)'clust',i,'      energy',energ(i)
	      do j=1,natoms
		 write(13,*)j,pop_off(i,j),pop_off(i,j+nmono)
	1	      ,pop_off(i,j+2*nmono)
		 if(pot_chos(1:3).eq.'tip')then
		    write(13,*)pop_off(i,j+3*nmono),pop_off(i,j+4*nmono)
	1		 ,pop_off(i,j+5*nmono)
		 endif
	      enddo
	   enddo
	   close(12)
cxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
	   
c       
c       - Out program <<<>>><<<>>><<<>>><<<>>><<<>>><<<>>><<<>>><<<>>><<<>>><<<>>><<<
c       - select kind of ga (elit) <<>>><<<>>><<<>>><<<>>><<<>>><<<>>><<<>>><<<>>><<<
	   if(opt_ga(1:4).eq.'elit')then
c       - Out of minimization ga elit >><<<>>><<<>>><<<>>><<<>>><<<>>><<<>>><<<>>><<<
	      dif_aux=en_aver-energ_bet(1)
	      dif_aux1=en_aver-ene_avbf
	      if(abs(dif_aux).gt.1d-5)ncont_aux_ga=0 
c	1	  .and. dif_aux1.gt.1d-5)ncont_aux_ga=0
	      if(ncont_aux_ga.ge.nterm_aft_rep)goto 2222
	      if(ngen.ge.nmaxgen. and .iclus.eq.1)then
		 do i=1,nclust,1
		    energ_best=energ_bet(i)
		    v2tot=vm2_bet(1)
		    v3tot=vm3_bet(1)
		    filename=filename(1:nfilename)//'.int'
		    nfil=nfilename+4
		    call print_2(pot_chos,filename,nfil,4,nclust
	1		 ,natoms,namo,i,energ_best,v2tot,v3tot
	2		 ,pop_clus(i,1))
		 enddo
		 goto 3333
	      endif
	      if(ngen.ge.nmaxgen .and. iclus.ge.2)goto 2222
c       
c       - Loop minimization ga elit <>>><<<>>><<<>>><<<>>><<<>>><<<>>><<<>>><<<>>><<<
	      ene_avbf=en_aver
	      ncont_aux_ga=ncont_aux_ga+1
	      goto 1111
c
c       - Out of minimization ga elit >><<<>>><<<>>><<<>>><<<>>><<<>>><<<>>><<<>>><<< 
 2222	      continue
	      close(2)
	      close(3)
c       
c       - Print best clusters ga elit >><<<>>><<<>>><<<>>><<<>>><<<>>><<<>>><<<>>><<<
	      energ_best=energ_bet(1)
	      v2tot=vm2_bet(1)
	      v3tot=vm3_bet(1)
	      filename=filename(1:nfilename)//'.xyz'
	      nfil=nfilename+4
	      if(pot_chos(1:3).eq.'tip')call angcart(nclust
	1	   ,natoms,namo,pop_clus)
	      call print_2('noneed',filename,nfil,4,nclust
	1	   ,natoms,namo,i,energ_best,v2tot,v3tot
	2	   ,pop_clus)
c
c       - Set energy for predator ga elit <>>><<<>>><<<>>><<<>>><<<>>><<<>>><<<>>><<<
	      if(iclus.eq.1)enr_1=energ_bet(1)-enr_nxt_m
	      enr_bef=enr_min
	      enr_min=energ_bet(1)+enr_nxt_m
c       
c       - select kind of ga (noelit) >>><<<>>><<<>>><<<>>><<<>>><<<>>><<<>>><<<>>><<<
	   elseif(opt_ga(1:6).eq.'noelit')then
c       
c       - Out of minimization ga multno <<<>>><<<>>><<<>>><<<>>><<<>>><<<>>><<<>>><<<
	      dif_aux=en_aver-ene_avbf
	      if(abs(dif_aux).gt.1d-5)ncont_aux_ga=0
	      if(ncont_aux_ga.ge.nterm_aft_rep)goto 5555
	      if(ngen.ge.nmaxgen. and .iclus.eq.1)then
		 do i=1,nclust,1
		    energ_best=energ_bet(i)
		    v2tot=vm2_bet(1)
		    v3tot=vm3_bet(1)
		    filename=filename(1:nfilename)//'.int'
		    nfil=nfilename+4
		    call print_2(pot_chos,filename,nfil,4,nclust
	1		 ,natoms,namo,iclus,energ_best,v2tot,v3tot
	2		 ,pop_clus(i,1))
		 enddo
		 goto 3333	
	      endif	  
	      if(ngen.ge.nmaxgen. and .iclus.ge.2)goto 5555
c       
c       - Loop minimization ga multno >><<<>>><<<>>><<<>>><<<>>><<<>>><<<>>><<<>>><<<
	      ene_avbf=en_aver
	      ncont_aux_ga=ncont_aux_ga+1
	      goto 1111
c
c       - Out of minimization ga multno <<<>>><<<>>><<<>>><<<>>><<<>>><<<>>><<<>>><<< 
 5555	      continue
	      close(2)
	      close(3)
c
c       - Set energy for predator ga multno >><<<>>><<<>>><<<>>><<<>>><<<>>><<<>>><<<
	      enr_min=energ_bet(1)+enr_nxt_m
c       
c       - Print best clusters ga multno <<<>>><<<>>><<<>>><<<>>><<<>>><<<>>><<<>>><<<
	      energ_best=energ_bet(1)
	      v2tot=vm2_bet(1)
	      v3tot=vm3_bet(1)
	      filename=filename(1:nfilename)//'.xyz'
	      nfil=nfilename+4
	      if(pot_chos(1:3).eq.'tip')call angcart(nclust
	1	   ,natoms,namo,pop_clus)
	      call print_2('noneed',filename,nfil,4,nclust
	1	   ,natoms,namo,iclus,energ_best,v2tot,v3tot
	2	   ,pop_clus)
c       
	   else
	      write(*,*)'option GA type not found'
	      write(*,*)'error: ga1'
	   endif
c
c       - Loop best clusters <<>>><<<>>><<<>>><<<>>><<<>>><<<>>><<<>>><<<>>><<<>>><<<
	enddo
c
c       - Out not convergence <>>><<<>>><<<>>><<<>>><<<>>><<<>>><<<>>><<<>>><<<>>><<<
 3333	continue
	close(4)
c       
c       - Finish <<>>><<<>>><<<>>><<<>>><<<>>><<<>>><<<>>><<<>>><<<>>><<<>>><<<>>><<<
	if(opt_out(1:3).eq.'yes')then
	   write(*,*)'       DDDD   OOOOOO NN  NN  EEEEE !!'
	   write(*,*)'       DD  D  OO  OO NNN NN  EE    !!'
	   write(*,*)'       DD  DD OO  OO NN NNN  EEEE  !!'
	   write(*,*)'       DD  D  OO  OO NN  NN  EE    !!'
	   write(*,*)'       DDDD   OOOOOO NN  NN  EEEEE ..'
	   write(*,*)'       ------------------------------'
	endif
	stop
c>><<<>>><<<>>><<<>>><<<>>><<<>>><<<>>><<<>>><<<>>><<<>>><<<>>><<<>>><<<>>><<<
	end
cxxxxxxxx1xxxxxxxxx2xxxxxxxxx3xxxxxxxxx4xxxxxxxxx5xxxxxxxxx6xxxxxxxxx7xxxxxxxxx
c234567890123456789012345678901234567890123456789012345678901234567890123456789
CcCcCc==CcCcCc==CcCcCc==CcCcCc==CcCcCc==CcCcCc==CcCcCc==CcCcCc==CcCcCc==CcCcCcc
cc/-\||/-\|/-\||/-\|/-\||/-\|/-\||/-\|/-\||/-\|/-\||/-\|/-\||/-\|/-\||/-\|/-\cc
cc\-/||\-/|\-/||\-/|\-/||\-/|\-/||\-/|\-/||\-/|\-/||\-/|\-/||\-/|\-/||\-/|\-/cc
CcCcCc==CcCcCc==CcCcCc==CcCcCc==CcCcCc==CcCcCc==CcCcCc==CcCcCc==CcCcCc==CcCcCcc


