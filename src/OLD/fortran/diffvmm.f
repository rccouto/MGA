	subroutine diffvmm(pot_chos,natoms,xyz,v2tot,v3tot,v,diffv)
	implicit real*8(a-h,o-z)
	include "param.inc"
	include "common.inc"
c       
	real*8   	  xyz(n_max_atom),diffv(n_max_atom)
	1                ,dv2dx(n_max_atom/3),dv2dy(n_max_atom/3)
	2                ,dv2dz(n_max_atom/3),dv3dx(n_max_atom/3)
	3                ,dv3dy(n_max_atom/3),dv3dz(n_max_atom/3)
c       
	character*20      pot_chos 
c  
	v2tot=0d0
	v3tot=0d0
	do 100 i=1,natoms
	   dv2dx(i)=0d0
	   dv2dy(i)=0d0
	   dv2dz(i)=0d0
	   dv3dx(i)=0d0
	   dv3dy(i)=0d0
	   dv3dz(i)=0d0
 100	continue
c	
	do 200 i=1,natoms-1
	   do 300 j=i+1,natoms    
	      dxij=xyz(i)-xyz(j)
	      dyij=xyz(natoms+i)-xyz(natoms+j)
	      dzij=xyz(2*natoms+i)-xyz(2*natoms+j)
	      rij2=dxij*dxij+dyij*dyij+dzij*dzij
c
	      if(rij2.ge.cutoff)goto 300
c
	      rij=sqrt(rij2)
	      v2tot=v2tot+pot_2b(pot_chos,rij)
c 	    
	      dvdx=dvdl_2b(pot_chos,dxij,rij)
	      dvdy=dvdl_2b(pot_chos,dyij,rij)
	      dvdz=dvdl_2b(pot_chos,dzij,rij)
c	     	   	
	      dv2dx(i)=dv2dx(i)+dvdx
	      dv2dy(i)=dv2dy(i)+dvdy
	      dv2dz(i)=dv2dz(i)+dvdz
	      dv2dx(j)=dv2dx(j)-dvdx
	      dv2dy(j)=dv2dy(j)-dvdy
	      dv2dz(j)=dv2dz(j)-dvdz
c
c   	 ... 3-body loop
c	
	      do 400 k=j+1,natoms
		 dxik=xyz(i)-xyz(k)
		 dyik=xyz(natoms+i)-xyz(natoms+k)
		 dzik=xyz(2*natoms+i)-xyz(2*natoms+k)
		 dxjk=xyz(j)-xyz(k)
		 dyjk=xyz(natoms+j)-xyz(natoms+k)
		 dzjk=xyz(2*natoms+j)-xyz(2*natoms+k)
		 rki2=dxik*dxik+dyik*dyik+dzik*dzik
		 rjk2=dxjk*dxjk+dyjk*dyjk+dzjk*dzjk
		 if(rki2.ge.cutoff. or .rjk2.ge.cutoff)goto 400
		 rki=sqrt(rki2)
		 rjk=sqrt(rjk2)
c
		 sq1=0.5773502691896257645d0
		 sq2=0.7071067811865475244d0
		 sq3=0.40824829046386301635d0
		 rhoij=(rij-re)/re
		 rhojk=(rjk-re)/re
		 rhoki=(rki-re)/re
		 Q1=sq1*(rhoij+rhojk+rhoki)
		 Q2=sq2*(rhojk-rhoki)				
		 Q3=sq3*(2d0*rhoij-rhojk-rhoki)
		 PQ123=c0+c1*Q1+c2*Q1**2+c3*(Q2**2+Q3**2)+c4*Q1**3
	1	      +c5*Q1*(Q2**2+Q3**2)+c6*(Q3**3-3d0*Q3*Q2**2)
		 pot_3b=De*PQ123*damp_mm(a3,Q1)
		 v3tot=v3tot+pot_3b
c
		 dvdp=De*damp_mm(a3,Q1)
		 dvdQ1=dvdp*(c1+2d0*c2*Q1+3d0*c4*Q1**2+c5*(Q2**2+Q3**2))
	1	      +De*PQ123*ddamp_mm(a3,Q1)
		 dvdQ2=dvdp*(2d0*c3*Q2+2d0*c5*Q1*Q2-6d0*c6*Q3*Q2) 
		 dvdQ3=dvdp*(2d0*c3*Q3+2d0*c5*Q1*Q3+3d0*c6*Q3**2-3d0*c6*Q2**2)
c
		 dQ1drholf=0.5773502691896257645d0
		 dQ2drhojk=0.7071067811865475244d0
		 dQ2drhoki=-0.7071067811865475244d0
		 dQ3drhoij=.81649658092772603272d0
		 dQ3drhojk=-0.40824829046386301635
		 dQ3drhoki=-0.40824829046386301635
c
		 dvdrholf12=dvdQ1*dQ1drholf
	1	      +dvdQ3*dQ3drhoij
		 dvdrholf23=dvdQ1*dQ1drholf
	1	      +dvdQ2*dQ2drhojk
	2	      +dvdQ3*dQ3drhojk
		 dvdrholf31=dvdQ1*dQ1drholf
	1	      +dvdQ2*dQ2drhoki
	2	      +dvdQ3*dQ3drhoki
c
		 drholfdrlf=1d0/re
		 rij1=1d0/rij
		 rki1=1d0/rki
		 rjk1=1d0/rjk
c
		 dvdxij=dvdrholf12*drholfdrlf*dxij*rij1
		 dvdyij=dvdrholf12*drholfdrlf*dyij*rij1
		 dvdzij=dvdrholf12*drholfdrlf*dzij*rij1
		 dvdxik=dvdrholf31*drholfdrlf*dxij*rjk1
		 dvdyik=dvdrholf31*drholfdrlf*dyij*rjk1
		 dvdzik=dvdrholf31*drholfdrlf*dzij*rjk1
		 dvdxjk=dvdrholf23*drholfdrlf*dxij*rki1
		 dvdyjk=dvdrholf23*drholfdrlf*dyij*rki1
		 dvdzjk=dvdrholf23*drholfdrlf*dzij*rki1
c	
		 dv3dx(i)=dv3dx(i)+dvdxij+dvdxik
		 dv3dy(i)=dv3dy(i)+dvdyij+dvdyik
		 dv3dz(i)=dv3dz(i)+dvdzij+dvdzik
		 dv3dx(j)=dv3dx(j)-dvdxij+dvdxjk
		 dv3dy(j)=dv3dy(j)-dvdyij+dvdyjk
		 dv3dz(j)=dv3dz(j)-dvdzij+dvdzjk
		 dv3dx(k)=dv3dx(k)-dvdxik-dvdxjk
		 dv3dy(k)=dv3dy(k)-dvdyik-dvdyjk
		 dv3dz(k)=dv3dz(k)-dvdzik-dvdzjk
 400	      continue
 300	   continue
 200	continue   	   	   	   
c   	   	   	
c     total potential energy
	v=v2tot+v3tot      		
c     calculate total first derivatives
c     diffv is a single 1-d array which stores all derivatives
c     diffv(i)=dv/dx(i), diffv(nat+i)=dv/dy(i), diffv(2nat+i)=dv/dz(i)
c
	do 500 i=1,natoms,1
	   diffv(i)=dv2dx(i)+dv3dx(i)
	   diffv(natoms+i)=dv2dy(i)+dv3dy(i)
	   diffv(2*natoms+i)=dv2dz(i)+dv3dz(i)
 500	continue
c
   	return
   	end
c
        function damp_mm(a3,Q1)
        implicit real*8(a-h,o-z)
        arg=a3*Q1
        damp_mm=1d0/cosh(arg)
        return
        end

        function ddamp_mm(a3,q1)
        implicit real*8(a-h,o-z)
c
        arg=a3*Q1
        t=tanh(arg)
        ddamp_mm=-a3*t*damp_mm(a3,q1)
c
        return
        end

