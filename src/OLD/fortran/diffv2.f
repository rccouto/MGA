	subroutine diffv2(pot_chos,natoms,xyz,v,diffv)
	implicit real*8(a-h,o-z)
	include "param.inc"
c
	real*8         xyz(3*natoms),diffv(n_max_atom)
	1    ,dv2dx(n_max_atom/3),dv2dy(n_max_atom/3)
	2    ,dv2dz(n_max_atom/3)
c 
	character*20   pot_chos
c  
	ambda=0
	v=0d0   
	do 100 i=1,natoms
	dv2dx(i)=0d0
	dv2dy(i)=0d0
	dv2dz(i)=0d0
100	continue
c  
	do 200 i=1,natoms-1
             do 300 j=i+1,natoms
	     dxij=xyz(i)-xyz(j)
             dyij=xyz(natoms+i)-xyz(natoms+j)
             dzij=xyz(2*natoms+i)-xyz(2*natoms+j)
c       dxij=dxij+ambda
c       dyij=dyij+ambda
c       dzij=dzij+ambda
c
             rij=dsqrt(dxij*dxij+dyij*dyij+dzij*dzij) 
	     rij=(rij+ambda)/(1+ambda)
c
	     v=v+pot_2b(pot_chos,rij)
c	write(*,*)'tt',v,pot_2b(pot_chos,rij) 
c	read(*,*)
c	  
	     dvdr1=dvdl_2b(pot_chos,rij)
             dvdx=dvdr1*dxij             
             dvdy=dvdr1*dyij
             dvdz=dvdr1*dzij 
c
             dv2dx(i)=dv2dx(i)+dvdx
             dv2dy(i)=dv2dy(i)+dvdy
             dv2dz(i)=dv2dz(i)+dvdz
             dv2dx(j)=dv2dx(j)-dvdx
             dv2dy(j)=dv2dy(j)-dvdy
             dv2dz(j)=dv2dz(j)-dvdz	
c
300	     continue
200	continue
c
c ... total potential energy
c ... calculate total first derivatives
c ... diffv is a single 1-d array which stores all derivatives
c ... diffv(i)=dv/dx(i), diffv(natoms+i)=dv/dy(i), diffv(2*natoms+i)=dv/dz(i)
c
	do 400 i=1,natoms
        diffv(i)=dv2dx(i)
        diffv(natoms+i)=dv2dy(i)
        diffv(2*natoms+i)=dv2dz(i)
400	continue
c
	return
	end
