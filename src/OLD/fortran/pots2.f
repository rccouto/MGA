	function pot_2b(pot_chos,rij)
	implicit real*8(a-h,o-z)
	include "common.inc"
c
	character*20   pot_chos
c 
	if(pot_chos(1:5).eq.'morse')then
	   afac=dexp(alpha*(1.0d0-(rij+1d-99)))
	   pot_2b=afac*(afac-2.0d0)
c	   
	elseif(pot_chos(1:9).eq.'tlr_morse')then
	   afac=dexp(alpha*(1.0d0-rij))
	   pot_2b=afac*(afac-2.0d0)+c6/rij**6
c	
	elseif(pot_chos(1:3).eq.'mm_')then
	   rhoij=(rij-re)/re
	   pot_2b=-De*((1d0+a2*rhoij)*exp(-a2*rhoij))
c
	elseif(pot_chos(1:5).eq.'gupta')then
	   rhoij=(rij-re)/re
	   Vr=Ag*exp(-pg*rhoij)
	   Vm=sqrt(eta**2*exp(-2d0*qg*rhoij))             
	   pot_2b=Vr+Vm
c	   
	elseif(pot_chos(1:8).eq.'lenjones')then   
	   pot_2b=4*epsi*((alpha/rij)**12-(alpha/rij)**6)
c
	elseif(pot_chos(1:8).eq.'m_nob_gas')then   
	   pot_2b=a0*(1+a1*rij+a2*rij**2+a3*rij**3)
     &           *exp(-a4*rij)-tanh(rij-rm2)*(c6/rij**6
     &           +c8/rij**8+c10/rij**10)
c	
	else
	   write(*,*)'Potential not Found'
	   write(*,*)'error:2b'
	   stop
	endif	
	
	return
	end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc	
	function dvdl_2b(pot_chos,rij)
	implicit real*8(a-h,o-z)
	include "common.inc"
c
	character*20   pot_chos
c
	rij=rij+1d-99
	rij1=1d0/(rij+1d-99)
	pi=acos(-1d0)
c
	if(pot_chos(1:5).eq.'morse')then
	   aa=2d0*alpha
	   afac=dexp(alpha*(1d0-(rij+1d-99))) 
	   dvdr=aa*afac*(1.0d0-afac)
	   dvdl_2b=dvdr*rij1
c	
	elseif(pot_chos(1:9).eq.'tlr_morse')then
	   aa=2d0*alpha
	   afac=dexp(alpha*(1d0-rij)) 
	   dvdr=aa*afac*(1.0d0-afac)-6d0*c6*rij1**7
	   dvdl_2b=dvdr*rij1	
c
	elseif(pot_chos(1:3).eq.'mm_')then
	   rhoij=(rij-re)/re  
	   dvdr=De*a2*dexp(-a2*rhoij)*a2*rhoij/re           
	   dvdl_2b=dvdr*rij1
c
        elseif(pot_chos(1:5).eq.'gupta')then
	   rhoij=(rij-re)/re
	   dVrdr=-pg*Ag*exp(-pg*rhoij)/re
	   dVmdr=5d-1*(eta**2*-2d0*qg*exp(-2d0*qg*rhoij))
     &          /sqrt(eta**2*exp(-2d0*qg*rhoij))*re             
	   dvdr=dVrdr+dVmdr
           dvdl_2b=dvdr*rij1
c	   
	elseif(pot_chos(1:8).eq.'lenjones')then  
	   dvdr=4*epsi*((-12d0*alpha**12/rij**13)
     &         -(-6d0*alpha**6/rij**7))
	   dvdl_2b=dvdr*rij1
c
	elseif(pot_chos(1:8).eq.'m_nob_gas')then   
	   dvdr=-a0*a4*(1+a1*rij+a2*rij**2+a3*rij**3)
     &         *exp(-a4*rij)+a0*(a1+2*a2*rij+3*a3*rij**2)
     &         *exp(-a4*rij)-tanh(rij-rm2)*(-6d0*c6/rij**7
     &         -8d0*c8/rij**9-10d0*c10/rij**11)
     &         -(1d0-tanh(rij-rm2)**2)
     &         *(c6/rij**6+c8/rij**8+c10/rij**10)
           dvdl_2b=dvdr*rij1
c
	else
	   write(*,*)'Potential derivate not Found'
	   write(*,*)'error:2b'
	   stop
	endif
c		
	return
	end
