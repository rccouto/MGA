      subroutine read_pot_const(pot_chos,inputfile,ifilelen,natoms,nmono
     &                           ,nat)
	implicit real*8(a-h,o-z)
	include "param.inc"
	include "common.inc"
c
	common/el/     elem    
	common/cst/    cst1,cst2,cst3,cst4,cst5
c
	integer        nat(10)
	real*8         cst1(7,7),cst2(7,7),cst3(7,7),cst4(7,7)
     &                ,cst5(7,7),cst6(7,7),cst7(7,7),cst8(7,7)
     &                ,cst9(7,7),cst10(7,7),cst11(7,7)
	character*3    elem(n_max_atom)
	character*20   pot_chos,inputfile
	character*60   title   
c
        write(*,*)'dentro do subr read_pot_const'
        write(*,*)pot_chos
	if(pot_chos(1:5).eq.'morse')then
	  nmono=1
	  open(2,status='old',file=inputfile(1:ifilelen))
	  read(2,*)alpha
	  read(2,'(a,f5.4)')elem(501),bscale
c
	elseif(pot_chos(1:4).eq.'tlr_')then
	  nmono=1
	  open(2,status='old',file=inputfile(1:ifilelen))
	  read(2,*)alpha
	  read(2,*)c6,c8,c10
	  read(2,'(a,f5.4)')elem(501),bscale
c
	elseif(pot_chos(1:2).eq.'mm')then
           write(*,*) 'estou aqui'
	  nmono=1
	  open(2,status='old',file=inputfile(1:ifilelen))	 
	  read(2,'(a)')title
	  read(2,*)de,re
	  read(2,*)a2,a3
	  read(2,*)c0,c1,c2
	  read(2,*)c3,c4,c5	
	  read(2,*)c6,c7,c8
	  read(2,*)cutoff
c	  read(2,'(a,f5.4)')elem(501),bscale 
          read(2,*)elem(501),bscale 


          write(*,*)title
	  write(*,*)de,re
	  write(*,*)a2,a3
	  write(*,*)c0,c1,c2
	  write(*,*)c3,c4,c5	
	  write(*,*)c6,c7,c8
	  write(*,*)cutoff
	  write(*,*)elem,bscale 

c
	elseif(pot_chos(1:5).eq.'gupta')then

           write(*,*)'Inside gupta pot param'

	   open(2,status='old',file=inputfile(1:ifilelen))	 
c	   read(2,'(a)')title
	   i = 0
 10	   continue
	   i = i + 1
	     read(2,'(a)')title
	     read(2,*)cst1(i,1)              ! re
	     read(2,*)cst2(i,1),cst3(i,1)    ! A,epsilon
	     read(2,*)cst4(i,1),cst5(i,1)    ! p,q
             read(2,*)cutoff
	  if(nmono.eq.2 .and. i.lt.3)goto 10
	 
	  do i=1,nmono,1
	     read(2,'(a,f5.4)')elem(500+i),bscale 
	  enddo
c
	elseif(pot_chos(1:8).eq.'lenjones')then	
          nmono=1					
	  open(2,status='old',file=inputfile(1:ifilelen))
	  read(2,*)epsi
	  read(2,*)alpha
	  read(2,'(a,f5.4)')elem(501),bscale   
c
	elseif(pot_chos(1:8).eq.'m_nob_gas')then 
          nmono=1
	  open(2,status='old',file=inputfile(1:ifilelen))
	  read(2,*)a0,a1,a2,a3,a4
	  read(2,*)rm2
	  read(2,*)c6,c8,c10
	  read(2,'(a,f5.4)')elem(501),bscale
c
	elseif(pot_chos(1:3).eq.'tip')then
          open(2,status='old',file=inputfile(1:ifilelen))
	  read(2,'(a)')title
	  read(2,'(a)')elem(501)
	  read(2,'(a)')elem(502)
	  read(2,'(a)')elem(503)
	  read(2,'(a)')elem(504)
	  read(2,'(a)')elem(505)
	  read(2,*)dioa1,dioa2,diob1,diob2
	  read(2,*)alph1,alph2
	  read(2,*)ecag
	  read(2,*)qia1,qia2,qib1,qib2
	  read(2,*)qja1,qja2,qjb1,qjb2
	  read(2,*)At,Ct,R_L,R_U
c
	else
	  write(*,*)'potential not found const.'
	endif
c
	nf=0
c	write(*,*)'teste',nmono
	do k=501,500+nmono,1
c	   write(*,*)'teste',k
	   ni=nf+1
	   nf=nf+nat(k-500)
	   do j=ni,nf,1 
	      elem(j)=elem(k)
	   enddo
	enddo
c
	close(2)
c
	return
	end


	
	
