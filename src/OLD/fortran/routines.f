ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc	  
	subroutine vec_matr(n,m,vec,matr)
	implicit real*8(a-h,o-z)
	include "param.inc"
c       
	real*8         matr(n_max_clus,n_max_atom)
	1             ,vec(n_max_atom)
c       
	do j=1,m,1        
	   matr(n,j)=vec(j)
	enddo
c       
	return
	end
c	              
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc   
	subroutine matr_vec(n,m,vec,matr)
	implicit real*8(a-h,o-z)
	include "param.inc"
c       
	real*8         matr(n_max_clus,n_max_atom)
	1             ,vec(n_max_atom)	         
c	
	do j=1,m,1        
	   vec(j)=matr(n,j)
	enddo
c       
	return
	end
c		
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc	
	subroutine matr1_matr2(n,m,matr1,matr2)
	implicit real*8(a-h,o-z)
	include "param.inc"
c
	real*8         matr1(n_max_clus,n_max_atom)
     &	              ,matr2(n_max_clus,n_max_atom)         
c	         
	do i=1,n,1
	   do j=1,m,1
	      matr2(i,j)=matr1(i,j)
	   enddo
	enddo
c
	return
	end
c	
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine rotate(n,nseed,xyz,ni,n_rot,nh)
	implicit real*8(a-h,o-z)
	include "param.inc" 
c
	real*8          xyz(n_max_atom)
c
	pi=acos(-1d0)
	rot_thet=0d0
	if(nh.eq.1 .or. nh.eq.3)rot_phi=2d0*pi*ran3(nseed) !robson 05/04
	if(nh.eq.3)rot_thet=2d0*pi*ran3(nseed) !robson 05/04
	do i=ni,n_rot,1
	   rd=sqrt(xyz(2*n+i)**2+xyz(n+i)**2+xyz(i)**2) 
	   theta=acos(xyz(2*n+i)/(rd+1d-99))
	   cq=(xyz(n+i)+1d-99)/abs(xyz(n+i)+1d-99)
	   if(xyz(i).eq.0d0)then
	      phi=90d0*pi/180d0
	   else
	      phi=cq*acos(xyz(i)/(rd*sin(theta))+1d-99)
	   endif
	   phi=phi+rot_phi
	   theta=theta+rot_thet
	   xyz(i)=rd*sin(theta)*cos(phi)
	   xyz(n+i)=rd*sin(theta)*sin(phi)
	   xyz(2*n+i)=rd*cos(theta)
	enddo
c
	return
	end
c 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine rotate1(natoms,poff,idum)
	implicit real*8(a-h,o-z)
	include "param.inc" 
c	
	real*8         poff(n_max_atom)
	external       ran3
c 
	pi=acos(-1d0)
	ti=ran3(idum)*pi
	si=ran3(idum)*2d0*pi
	ppi=ran3(idum)*2d0*pi
	theta = ran3(idum)*(pi/2.d0)
	phi = ran3(idum)*pi
c	phi=0d0
c	theta=0d0
	ti=0d0
	si=0d0
	ppi=0d0
c
	rot11=cos(ti)*cos(si)*cos(ppi)-sin(si)*sin(ppi)
	rot12=cos(ti)*sin(si)*cos(ppi)+cos(si)*sin(ppi)
	rot13=-sin(ti)*cos(ppi) 
	rot21=-cos(ti)*cos(si)*sin(ppi)-sin(si)*cos(ppi) 
	rot22=-cos(ti)*sin(si)*sin(ppi)+cos(si)*cos(ppi)
	rot23=sin(ti)*sin(ppi) 
	rot31=sin(ti)*cos(ppi)
	rot32=sin(ti)*sin(ppi) 
	rot33=cos(ti)     
	
	rot11 = cos(phi)
	rot12 = -cos(theta)*sin(phi)
	rot13 = sin(theta)*sin(phi) 
	rot21 = sin(phi)
	rot22 = cos(theta)*cos(phi)
	rot23 = -sin(theta)*cos(phi)
	rot31 = 0d0
	rot32 = sin(theta)
	rot33 = cos(theta)
c
	do i=1,natoms
	   poff(i)=rot11*poff(i)+rot12*poff(i+natoms)
	1	+rot13*poff(i+2*natoms)
	   poff(i+natoms)=rot21*poff(i)+rot22*poff(i+natoms)
	1	+rot23*poff(i+2*natoms)
	   poff(i+2*natoms)=rot31*poff(i)+rot32*poff(i+natoms)
	1	+rot33*poff(i+2*natoms)
	end do
c
	return
	end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine sort_matr(pot_chos,n,m,vec2,vec3,vec,satr)
	implicit real*8(a-h,o-z)
	include "param.inc" 
c       
	real*8         satr(n_max_clus,n_max_atom)
	1             ,snewsatr(n_max_clus,n_max_atom)
	2             ,vecat(0:n_max_clus),vec(n_max_clus)
	3             ,vec2(n_max_clus),vec2cat(n_max_clus)
	4             ,vec3(n_max_clus),vec3cat(n_max_clus)
	5             ,nch(n_max_clus)
	character*20   pot_chos
c
	vecat(0)=-9d299
	do i=1,n,1
	   vecat(i)=9d299
	   vec2cat(i)=0d0
	   vec3cat(i)=0d0
	   nch(i)=0
	enddo
c
	do k=1,n,1
	   do i=1,n,1
	      if(vec(i).le.vecat(k) .and. vec(i).ge.vecat(k-1)
	1	   .and. nch(i).eq.0)then
		 vecat(k)=vec(i)
		 npos=i
	       endif
	    enddo
c
	    nch(npos)=1
	    vec2cat(k)=vec2(npos)
	    vec3cat(k)=vec3(npos)
	    do j=1,m,1	
	       snewsatr(k,j)=satr(npos,j)
	       snewsatr(k,j+m)=satr(npos,j+m)
	       snewsatr(k,j+2*m)=satr(npos,j+2*m)
	       if(pot_chos(1:3).eq.'tip')then
		 snewsatr(k,j+3*m)=satr(npos,j+3*m)
		 snewsatr(k,j+4*m)=satr(npos,j+4*m)
		 snewsatr(k,j+5*m)=satr(npos,j+5*m) 
	       endif
	    enddo
	 enddo
c
	 do i=1,n,1
	    vec2(i)=vec2cat(i)
	    vec3(i)=vec3cat(i)
	    vec(i)=vecat(i)
	      do j=1,m,1
		 satr(i,j)=snewsatr(i,j)
		 satr(i,j+m)=snewsatr(i,j+m)
		 satr(i,j+2*m)=snewsatr(i,j+2*m)
		 if(pot_chos(1:3).eq.'tip')then
		   satr(i,j+3*m)=snewsatr(i,j+3*m)
		   satr(i,j+4*m)=snewsatr(i,j+4*m)
		   satr(i,j+5*m)=snewsatr(i,j+5*m) 
	        endif
	      enddo
	   enddo
c
	return
	end
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine centre(n,re,vec)
      	implicit real*8(a-h,o-z)
	include "param.inc" 
c      	
	real*8         vec(n_max_atom)
	character*10   opt
c      
	opt='no'
	pi=acos(-1d0)	
c
	do i=1,n,1
	   xmass=xmass+vec(i)
	   ymass=ymass+vec(n+i)
	   zmass=zmass+vec(2*n+i)
	enddo
	xmass=xmass/n
	ymass=ymass/n
	zmass=zmass/n
c       
      	do i=1,n,1
	   vec(i)=vec(i)-xmass
	   vec(n+i)=vec(n+i)-ymass
	   vec(2*n+i)=vec(2*n+i)-zmass
	enddo
c
c- Check that all the atoms are in the container. If not then rescale.
c	if(opt(1:3).eq.'chk')then
c	   radius=re*(1+(3*n/(4*pi*2**0.5))**(1d0/3d0))
c	   distmax=0d0
c	   do j=1,n,1
c	      dist=vec(j)**2+vec(n+j)**2+vec(2*n+j)**2
c	      if(dist.gt.distmax)distmax=dist
c	   enddo
c       
c      	  if(distmax.gt.radius)then
c	     cort=sqrt(distmax/radius)*99d-3
c	     do j=1,n,1
c		vec(i)=vec(i)*cort
c		vec(n+i)=vec(n+i)*cort
c		vec(2*n+i)=vec(2*n+i)*cort
c	     enddo
c          endif
c	endif
c       
	return
	end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine fitness(nclust,fit,energ)
	implicit real*8(a-h,o-z)
	include "param.inc" 
c
	real*8         fit(n_max_clus),energ(n_max_clus)
c      
	best=energ(1)
	erange=(energ(nclust)-best)+1d-99
c      
	do i=1,nclust
c	   fit(i)=0.5d0*(1-dtanh(2.d0*((energ(i)-best)/erange)-1.d0))
	   fit(i)=exp(3.d0*((best-energ(i))/erange)) 
	end do
c
	return
	end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


	








