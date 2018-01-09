	subroutine diffvgp(pot_chos,natoms,nat,xyz,vrtot,vmtot,v,diffv)
	implicit real*8(a-h,o-z)
	include "param.inc"
	include "common.inc"
c
	common/cst/    cst1,cst2,cst3,cst4,cst5
c
	integer        nat(10)
	real*8         xyz(n_max_atom),diffv(n_max_atom)
	1             ,dv2dx(n_max_atom/3),dv2dy(n_max_atom/3)
	2             ,dv2dz(n_max_atom/3),rij(n_max_atom/3)
	3             ,cst1(7,7),cst2(7,7),cst3(7,7),cst4(7,7)
	4             ,cst5(7,7)
	character*20   pot_chos
c     
	v=0d0
	vr=0d0
	vm=0d0
	vrtot=0d0
	vmtot=0d0 
c 
	do 100 i=1,natoms,1
	rij(i)=0d0
	dv2dx(i)=0d0
	dv2dy(i)=0d0
	dv2dz(i)=0d0
100	continue
c  
	do 200 i=1,natoms,1
        rij(i)=0d0
             do 300 j=1,natoms,1
	     if(i.eq.j)goto 300
	     dxij=xyz(i)-xyz(j)
             dyij=xyz(natoms+i)-xyz(natoms+j)
             dzij=xyz(2*natoms+i)-xyz(2*natoms+j)
             rij(j)=dsqrt(dxij*dxij+dyij*dyij+dzij*dzij)
c	  write(*,*)'to'
c	  write(*,*)'dist',rij(j)
	  if(rij(j).ge.30)goto 300
	     if(i.le.nat(1) .and. j.le.nat(1))nn=1
	     if(i.gt.nat(1) .and. j.le.nat(1))nn=2
	     if(i.le.nat(1) .and. j.gt.nat(1))nn=2
 	     if(i.gt.nat(1) .and. j.gt.nat(1))nn=3
	     re=cst1(nn,1)
c	     write(*,*)'re',re,nn
c	     read(*,*)
	     rhoij=(rij(j)-re)/re
	     eta=cst3(nn,1)
	     qg=cst5(nn,1)
	     vm=vm+eta**2*exp(-2d0*qg*rhoij)
300	     continue
c
	     do 400 j=1,natoms,1	
		if(i.eq.j)goto 400
		if(rij(j).ge.30)goto 400
		dxij=xyz(i)-xyz(j)
		dyij=xyz(natoms+i)-xyz(natoms+j)
		dzij=xyz(2*natoms+i)-xyz(2*natoms+j)
c
		if(i.le.nat(1) .and. j.le.nat(1))nn=1
		if(i.gt.nat(1) .and. j.le.nat(1))nn=2
		if(i.le.nat(1) .and. j.gt.nat(1))nn=2
		if(i.gt.nat(1) .and. j.gt.nat(1))nn=3
		re=cst1(nn,1)
		rhoij=(rij(j)-re)/re
		rij1=1d0/(rij(j)+1d-99)
		eta=cst3(nn,1)
		qg=cst5(nn,1)   
		Ag=cst2(nn,1)  
		pg=cst4(nn,1)         
c
		vr=vr+Ag*exp(-pg*rhoij)
c
		dv1dr=-pg*Ag*exp(-pg*rhoij)/re 
		dvdx1=dv1dr*dxij*rij1
		dvdy1=dv1dr*dyij*rij1
		dvdz1=dv1dr*dzij*rij1
c
		dv2dr=qg*eta**2*exp(-2d0*qg*rhoij)/re
	1	     *1d0/(sqrt(vm)+1d-99)
		dvdx2=dv2dr*dxij*rij1            
		dvdy2=dv2dr*dyij*rij1
		dvdz2=dv2dr*dzij*rij1	
c
		dvdx=dvdx1+dvdx2
		dvdy=dvdy1+dvdy2
		dvdz=dvdz1+dvdz2
c
		dv2dx(i)=dv2dx(i)+dvdx
		dv2dy(i)=dv2dy(i)+dvdy
		dv2dz(i)=dv2dz(i)+dvdz
		dv2dx(j)=dv2dx(j)-dvdx
		dv2dy(j)=dv2dy(j)-dvdy
		dv2dz(j)=dv2dz(j)-dvdz	
c
 400	     continue
	     vrtot=vrtot+vr
	     vmtot=vmtot+sqrt(vm)
	     vr=0d0
	     vm=0d0
 200	  continue
c
c ... total potential energy
	v=vrtot/1d0-vmtot
c ... calculate total first derivatives
c ... diffv is a single 1-d array which stores all derivatives
c ... diffv(i)=dv/dx(i), diffv(natoms+i)=dv/dy(i), diffv(2*natoms+i)=dv/dz(i)
c
	do 500 i=1,natoms
	   diffv(i)=dv2dx(i)
	   diffv(natoms+i)=dv2dy(i)
	   diffv(2*natoms+i)=dv2dz(i)
 500	continue
c       
	return
	end
