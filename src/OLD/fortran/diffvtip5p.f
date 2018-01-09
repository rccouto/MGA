       	subroutine diffvtip5p(pot_chos,natoms,nmono,xyz,v2tot,v3tot
     &                       ,v,diffv)
	implicit real*8(a-h,o-z)
	include "param.inc"
	include "common.inc"
        common/cpot1/qa1i,qa2i,qb1i,qb2i,qa1j,qa2j,qb1j,qb2j
        common/cpot2/AA,BB,ac,as,bc,bs,ec
        common/cpot3/cst1,cst2,cst3,cst4
c
	real*8         xyz(6*natoms/nmono),diffv(6*natoms/nmono)
     &                ,dv2dxo(natoms/nmono),dv2dyo(natoms/nmono)
     &                ,dv2dzo(natoms/nmono),dv2dt(natoms/nmono)
     &                ,dv2ds(natoms/nmono),dv2dpp(natoms/nmono)
	character*20   pot_chos 
c 
	v=0d0
        v2tot=0d0
      	v3tot=0d0
	namo=natoms/nmono
	if(namo.lt.999d-3)namo=1
c     
	do i=1,namo,1
	   dv2dxo(i)=0d0
	   dv2dyo(i)=0d0
	   dv2dzo(i)=0d0
           dv2dt(i)=0d0
           dv2ds(i)=0d0
           dv2dpp(i)=0d0
	enddo
c  
        pi=acos(-1d0)
        AA=At
        BB=Ct
        ap1=alph1*pi/360d0
        ap2=alph2*pi/360d0
        ac=dioa1*cos(ap1)
        as=dioa2*sin(ap1)
        bc=-diob1*cos(ap2)
        bs=diob2*sin(ap2)
        cst1=R_L
        cst2=R_U
        ec=ecag
        qa1i=qia1
        qa2i=qia2
        qb1i=qib1
        qb2i=qib2
        qa1j=qja1
        qa2j=qja2
        qb1j=qjb1
        qb2j=qjb2
c         write(*,*)'dioa1,dioa2,diob1,diob2',dioa1,dioa2,diob1,diob2
c        write(*,*)'ap1,ap2,ac,as,ec',ap1,ap2,ac,as,ec
c        write(*,*)'bc,bs,cst1,cst2',bc,bs,cst1,cst2
c        write(*,*)'qa1i,qa2i,qb1i,qb2i',qa1i,qa2i,qb1i,qb2i
c        write(*,*)'qa1j,qa2j,qb1j,qb2j',qa1j,qa2j,qb1j,qb2j
c        read(*,*)
c        
	do i=1,namo-1,1
           Xoi=xyz(i)
           Yoi=xyz(i+namo)
           Zoi=xyz(i+2*namo)
           ti=xyz(i+3*namo)
           si=xyz(i+4*namo)
           ppi=xyz(i+5*namo)
	   do j=i+1,namo,1
              Xoj=xyz(j)
              Yoj=xyz(j+namo)
              Zoj=xyz(j+2*namo)
              tj=xyz(j+3*namo)
              sj=xyz(j+4*namo)
              ppj=xyz(j+5*namo)
              dOOx=Xoi-Xoj
              dOOy=Yoi-Yoj
              dOOz=Zoi-Zoj
              rOO=sqrt(dOOx**2+dOOy**2+dOOz**2)
c              write(*,*)'rOO',rOO
              if(rOO.le.cst1)then
c              write(*,*)'1'   
                 cst3=0d0
                 cst4=0d0
              elseif(rOO.ge.cst1 .and. rOO.le.cst2)then
c              write(*,*)'2'
                 cst3=1d0
                 cst4=0d0
              elseif(rOO.ge.cst2)then
c              write(*,*)'3'
                 cst3=0d0
                 cst4=1d0
              endif
c     
	      v=v+vpot(Xoi,Xoj,Yoi,Yoj,Zoi,Zoj,ppi,si,ti,ppj,sj,tj)
c
              dvdx=dVdXoi(Xoi,Xoj,Yoi,Yoj,Zoi,Zoj
     &             ,ppi,si,ti,ppj,sj,tj)
              dvdy=dVdYoi(Xoi,Xoj,Yoi,Yoj,Zoi,Zoj
     &             ,ppi,si,ti,ppj,sj,tj) 
              dvdz=dVdZoi(Xoi,Xoj,Yoi,Yoj,Zoi,Zoj
     &             ,ppi,si,ti,ppj,sj,tj)
c     
	      dv2dxo(i)=dv2dxo(i)+dvdx
	      dv2dyo(i)=dv2dyo(i)+dvdy
	      dv2dzo(i)=dv2dzo(i)+dvdz
	      dv2dt(i)=dv2dt(i)+dVdti(Xoi,Xoj,Yoi,Yoj,Zoi,Zoj
     &             ,ppi,si,ti,ppj,sj,tj)
	      dv2ds(i)=dv2ds(i)+dVdsi(Xoi,Xoj,Yoi,Yoj,Zoi,Zoj
     &             ,ppi,si,ti,ppj,sj,tj)
	      dv2dpp(i)=dv2dpp(i)+dVdppi(Xoi,Xoj,Yoi,Yoj,Zoi,Zoj
     &             ,ppi,si,ti,ppj,sj,tj)
              dv2dxo(j)=dv2dxo(j)-dvdx
	      dv2dyo(j)=dv2dyo(j)-dvdy
	      dv2dzo(j)=dv2dzo(j)-dvdz
	      dv2dt(j)=dv2dt(j)+dVdtj(Xoi,Xoj,Yoi,Yoj,Zoi,Zoj
     &             ,ppi,si,ti,ppj,sj,tj)
	      dv2ds(j)=dv2ds(j)+dVdsj(Xoi,Xoj,Yoi,Yoj,Zoi,Zoj
     &             ,ppi,si,ti,ppj,sj,tj)
	      dv2dpp(j)=dv2dpp(j)+dVdppj(Xoi,Xoj,Yoi,Yoj,Zoi,Zoj
     &             ,ppi,si,ti,ppj,sj,tj)
c
c
c
c              d1=dVdXoi(Xoi,Xoj,Yoi,Yoj,Zoi,Zoj
c     &             ,ppi,si,ti,ppj,sj,tj)
c              d2=dVdXoj(Xoi,Xoj,Yoi,Yoj,Zoi,Zoj
c     &             ,ppi,si,ti,ppj,sj,tj)
c              d3=dVdYoi(Xoi,Xoj,Yoi,Yoj,Zoi,Zoj
c     &             ,ppi,si,ti,ppj,sj,tj)
c              d4=dVdYoj(Xoi,Xoj,Yoi,Yoj,Zoi,Zoj
c     &             ,ppi,si,ti,ppj,sj,tj)
c              d5=dVdZoi(Xoi,Xoj,Yoi,Yoj,Zoi,Zoj
c     &             ,ppi,si,ti,ppj,sj,tj)
c              d6=dVdZoj(Xoi,Xoj,Yoi,Yoj,Zoi,Zoj
c     &             ,ppi,si,ti,ppj,sj,tj)
c              d7=dVdti(Xoi,Xoj,Yoi,Yoj,Zoi,Zoj
c     &             ,ppi,si,ti,ppj,sj,tj)
c              d8=dVdtj(Xoi,Xoj,Yoi,Yoj,Zoi,Zoj
c     &             ,ppi,si,ti,ppj,sj,tj)
c              d9=dVdsi(Xoi,Xoj,Yoi,Yoj,Zoi,Zoj
c     &             ,ppi,si,ti,ppj,sj,tj)
c              d10=dVdsj(Xoi,Xoj,Yoi,Yoj,Zoi,Zoj
c     &             ,ppi,si,ti,ppj,sj,tj)
c              d11=dVdppi(Xoi,Xoj,Yoi,Yoj,Zoi,Zoj
c     &             ,ppi,si,ti,ppj,sj,tj)
c              d12=dVdppj(Xoi,Xoj,Yoi,Yoj,Zoi,Zoj
c     &             ,ppi,si,ti,ppj,sj,tj)
c     
c              write(*,*)'i,j,d1,d2',i,j,d1,d2
c              write(*,*)'d3,d4',d3,d4
c              write(*,*)'d5,d6',d5,d6
c              write(*,*)'d7,d8',d7,d8
c              write(*,*)'d9,d10',d9,d10
c              write(*,*)'d11,d12',d11,d12
c              write(*,*)'v',v
c     read(*,*)
           enddo
        enddo
c     
	v2tot=v
	v3tot=0d0
c ... total potential energy
c ... calculate total first derivatives
c ... diffv is a single 1-d array which stores all derivatives
c ... diffv(i)=dv/dx(i), diffv(natoms+i)=dv/dy(i), diffv(2*natoms+i)=dv/dz(i)
c
	do i=1,namo,1
  	   diffv(i)=dv2dxo(i)
	   diffv(i+namo)=dv2dyo(i)
	   diffv(i+2*namo)=dv2dzo(i)
           diffv(i+3*namo)=dv2dt(i)
	   diffv(i+4*namo)=dv2ds(i)
	   diffv(i+5*namo)=dv2dpp(i)
        enddo
c
	return
	end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      function vpot(Xoi,Xoj,Yoi,Yoj,Zoi,Zoj,ppi,si,ti,ppj,sj,tj)
      implicit real*8(a-h,o-z)
      common/cpot1/qa1i,qa2i,qb1i,qb2i,qa1j,qa2j,qb1j,qb2j
      common/cpot2/AA,BB,ac,as,bc,bs,ec  
      common/cpot3/cst1,cst2,cst3,cst4
c      
      t1 = Xoi**2
      t4 = Xoj**2
      t5 = Yoi**2
      t8 = Yoj**2
      t9 = Zoi**2
      t12 = Zoj**2
      t13 = t1-2*Xoi*Xoj+t4+t5-2*Yoi*Yoj+t8+t9-2*Zoi*Zoj+t12
      t14 = sqrt(t13)
      t16 = (t14-cst1)**2
      t21 = cst2-cst1
      t22 = t21**2
      t29 = cos(ppi)
      t30 = cos(si)
      t32 = cos(ti)
      t33 = sin(si)
      t34 = t32*t33
      t35 = sin(ppi)
      t37 = t29*t30-t34*t35
      t38 = t37*ac
      t40 = t32*t30
      t43 = (t29*t33+t40*t35)*as
      t44 = cos(ppj)
      t45 = cos(sj)
      t47 = cos(tj)
      t48 = sin(sj)
      t49 = t47*t48
      t50 = sin(ppj)
      t52 = t44*t45-t49*t50
      t53 = t52*ac
      t55 = t47*t45
      t58 = (t44*t48+t55*t50)*as
      t60 = (t38+t43+Xoi-t53-t58-Xoj)**2
      t63 = -t35*t30-t34*t29
      t64 = t63*ac
      t68 = (-t35*t33+t40*t29)*as
      t71 = -t50*t45-t49*t44
      t72 = t71*ac
      t76 = (-t50*t48+t55*t44)*as
      t78 = (t64+t68+Yoi-t72-t76-Yoj)**2
      t79 = sin(ti)
      t80 = t79*t33
      t81 = t80*ac
      t83 = t79*t30*as
      t84 = sin(tj)
      t85 = t84*t48
      t86 = t85*ac
      t88 = t84*t45*as
      t90 = (t81-t83+Zoi-t86+t88-Zoj)**2
      t92 = sqrt(t60+t78+t90)
      t98 = (t38+t43+Xoi-t53+t58-Xoj)**2
      t100 = (t64+t68+Yoi-t72+t76-Yoj)**2
      t102 = (t81-t83+Zoi-t86-t88-Zoj)**2
      t104 = sqrt(t98+t100+t102)
      t109 = t52*bc
      t111 = t50*t84*bs
      t113 = (t38+t43+Xoi-t109-t111-Xoj)**2
      t114 = t71*bc
      t116 = t44*t84*bs
      t118 = (t64+t68+Yoi-t114-t116-Yoj)**2
      t119 = t85*bc
      t120 = t47*bs
      t122 = (t81-t83+Zoi-t119-t120-Zoj)**2
      t124 = sqrt(t113+t118+t122)
      t130 = (t38+t43+Xoi-t109+t111-Xoj)**2
      t132 = (t64+t68+Yoi-t114+t116-Yoj)**2
      t134 = (t81-t83+Zoi-t119+t120-Zoj)**2
      t136 = sqrt(t130+t132+t134)
      t142 = (t38-t43+Xoi-t53-t58-Xoj)**2
      t144 = (t64-t68+Yoi-t72-t76-Yoj)**2
      t146 = (t81+t83+Zoi-t86+t88-Zoj)**2
      t148 = sqrt(t142+t144+t146)
      t154 = (t38-t43+Xoi-t53+t58-Xoj)**2
      t156 = (t64-t68+Yoi-t72+t76-Yoj)**2
      t158 = (t81+t83+Zoi-t86-t88-Zoj)**2
      t160 = sqrt(t154+t156+t158)
      t166 = (t38-t43+Xoi-t109-t111-Xoj)**2
      t168 = (t64-t68+Yoi-t114-t116-Yoj)**2
      t170 = (t81+t83+Zoi-t119-t120-Zoj)**2
      t172 = sqrt(t166+t168+t170)
      t178 = (t38-t43+Xoi-t109+t111-Xoj)**2
      t180 = (t64-t68+Yoi-t114+t116-Yoj)**2
      t182 = (t81+t83+Zoi-t119+t120-Zoj)**2
      t184 = sqrt(t178+t180+t182)
      t189 = t37*bc
      t191 = t35*t79*bs
      t193 = (t189+t191+Xoi-t53-t58-Xoj)**2
      t194 = t63*bc
      t196 = t29*t79*bs
      t198 = (t194+t196+Yoi-t72-t76-Yoj)**2
      t199 = t80*bc
      t200 = t32*bs
      t202 = (t199+t200+Zoi-t86+t88-Zoj)**2
      t204 = sqrt(t193+t198+t202)
      t210 = (t189+t191+Xoi-t53+t58-Xoj)**2
      t212 = (t194+t196+Yoi-t72+t76-Yoj)**2
      t214 = (t199+t200+Zoi-t86-t88-Zoj)**2
      t216 = sqrt(t210+t212+t214)
      t222 = (t189+t191+Xoi-t109-t111-Xoj)**2
      t224 = (t194+t196+Yoi-t114-t116-Yoj)**2
      t226 = (t199+t200+Zoi-t119-t120-Zoj)**2
      t228 = sqrt(t222+t224+t226)
      t234 = (t189+t191+Xoi-t109+t111-Xoj)**2
      t236 = (t194+t196+Yoi-t114+t116-Yoj)**2
      t238 = (t199+t200+Zoi-t119+t120-Zoj)**2
      t240 = sqrt(t234+t236+t238)
      t246 = (t189-t191+Xoi-t53-t58-Xoj)**2
      t248 = (t194-t196+Yoi-t72-t76-Yoj)**2
      t250 = (t199-t200+Zoi-t86+t88-Zoj)**2
      t252 = sqrt(t246+t248+t250)
      t258 = (t189-t191+Xoi-t53+t58-Xoj)**2
      t260 = (t194-t196+Yoi-t72+t76-Yoj)**2
      t262 = (t199-t200+Zoi-t86-t88-Zoj)**2
      t264 = sqrt(t258+t260+t262)
      t270 = (t189-t191+Xoi-t109-t111-Xoj)**2
      t272 = (t194-t196+Yoi-t114-t116-Yoj)**2
      t274 = (t199-t200+Zoi-t119-t120-Zoj)**2
      t276 = sqrt(t270+t272+t274)
      t282 = (t189-t191+Xoi-t109+t111-Xoj)**2
      t284 = (t194-t196+Yoi-t114+t116-Yoj)**2
      t286 = (t199-t200+Zoi-t119+t120-Zoj)**2
      t288 = sqrt(t282+t284+t286)
      t292 = qa1i*qa1j*ec/t92+qa1i*qa2j*ec/t104+qa1i*qb1j*ec/t124+qa1i*q
     #b2j*ec/t136+qa2i*qa1j*ec/t148+qa2i*qa2j*ec/t160+qa2i*qb1j*ec/t172+
     #qa2i*qb2j*ec/t184+qb1i*qa1j*ec/t204+qb1i*qa2j*ec/t216+qb1i*qb1j*ec
     #/t228+qb1i*qb2j*ec/t240+qb2i*qa1j*ec/t252+qb2i*qa2j*ec/t264+qb2i*q
     #b1j*ec/t276+qb2i*qb2j*ec/t288
      t294 = t13**2
      t295 = t294**2
      t302 = (cst3*t16*(3*cst2-cst1-2*t14)/t22/t21+cst4)*t292+AA/t295/t2
     #94-BB/t294/t13
c
      vpot=t302
c
      return
      end
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      function dVdXoi(Xoi,Xoj,Yoi,Yoj,Zoi,Zoj,ppi,si,ti,ppj,sj,tj)
      implicit real*8(a-h,o-z)
      common/cpot1/qa1i,qa2i,qb1i,qb2i,qa1j,qa2j,qb1j,qb2j
      common/cpot2/AA,BB,ac,as,bc,bs,ec 
      common/cpot3/cst1,cst2,cst3,cst4
c
      t1 = Xoi**2
      t4 = Xoj**2
      t5 = Yoi**2
      t8 = Yoj**2
      t9 = Zoi**2
      t12 = Zoj**2
      t13 = t1-2*Xoi*Xoj+t4+t5-2*Yoi*Yoj+t8+t9-2*Zoi*Zoj+t12
      t14 = sqrt(t13)
      t15 = t14-cst1
      t19 = 3*cst2-cst1-2*t14
      t21 = cst2-cst1
      t22 = t21**2
      t24 = 1/t22/t21
      t27 = Xoi-Xoj
      t28 = 2*t24/t14*t27
      t30 = t15**2
      t31 = cst3*t30
      t34 = qa1i*qa1j
      t35 = cos(ppi)
      t36 = cos(si)
      t38 = cos(ti)
      t39 = sin(si)
      t40 = t38*t39
      t41 = sin(ppi)
      t43 = t35*t36-t40*t41
      t44 = t43*ac
      t46 = t38*t36
      t49 = (t35*t39+t46*t41)*as
      t50 = cos(ppj)
      t51 = cos(sj)
      t53 = cos(tj)
      t54 = sin(sj)
      t55 = t53*t54
      t56 = sin(ppj)
      t58 = t50*t51-t55*t56
      t59 = t58*ac
      t61 = t53*t51
      t64 = (t50*t54+t61*t56)*as
      t65 = t44+t49+Xoi-t59-t64-Xoj
      t66 = t65**2
      t69 = -t41*t36-t40*t35
      t70 = t69*ac
      t74 = (-t41*t39+t46*t35)*as
      t77 = -t56*t51-t55*t50
      t78 = t77*ac
      t82 = (-t56*t54+t61*t50)*as
      t84 = (t70+t74+Yoi-t78-t82-Yoj)**2
      t85 = sin(ti)
      t86 = t85*t39
      t87 = t86*ac
      t89 = t85*t36*as
      t90 = sin(tj)
      t91 = t90*t54
      t92 = t91*ac
      t94 = t90*t51*as
      t96 = (t87-t89+Zoi-t92+t94-Zoj)**2
      t98 = sqrt(t66+t84+t96)
      t102 = qa1i*qa2j
      t103 = t44+t49+Xoi-t59+t64-Xoj
      t104 = t103**2
      t106 = (t70+t74+Yoi-t78+t82-Yoj)**2
      t108 = (t87-t89+Zoi-t92-t94-Zoj)**2
      t110 = sqrt(t104+t106+t108)
      t114 = qa1i*qb1j
      t115 = t58*bc
      t117 = t56*t90*bs
      t118 = t44+t49+Xoi-t115-t117-Xoj
      t119 = t118**2
      t120 = t77*bc
      t122 = t50*t90*bs
      t124 = (t70+t74+Yoi-t120-t122-Yoj)**2
      t125 = t91*bc
      t126 = t53*bs
      t128 = (t87-t89+Zoi-t125-t126-Zoj)**2
      t130 = sqrt(t119+t124+t128)
      t134 = qa1i*qb2j
      t135 = t44+t49+Xoi-t115+t117-Xoj
      t136 = t135**2
      t138 = (t70+t74+Yoi-t120+t122-Yoj)**2
      t140 = (t87-t89+Zoi-t125+t126-Zoj)**2
      t142 = sqrt(t136+t138+t140)
      t146 = qa2i*qa1j
      t147 = t44-t49+Xoi-t59-t64-Xoj
      t148 = t147**2
      t150 = (t70-t74+Yoi-t78-t82-Yoj)**2
      t152 = (t87+t89+Zoi-t92+t94-Zoj)**2
      t154 = sqrt(t148+t150+t152)
      t158 = qa2i*qa2j
      t159 = t44-t49+Xoi-t59+t64-Xoj
      t160 = t159**2
      t162 = (t70-t74+Yoi-t78+t82-Yoj)**2
      t164 = (t87+t89+Zoi-t92-t94-Zoj)**2
      t166 = sqrt(t160+t162+t164)
      t170 = qa2i*qb1j
      t171 = t44-t49+Xoi-t115-t117-Xoj
      t172 = t171**2
      t174 = (t70-t74+Yoi-t120-t122-Yoj)**2
      t176 = (t87+t89+Zoi-t125-t126-Zoj)**2
      t178 = sqrt(t172+t174+t176)
      t182 = qa2i*qb2j
      t183 = t44-t49+Xoi-t115+t117-Xoj
      t184 = t183**2
      t186 = (t70-t74+Yoi-t120+t122-Yoj)**2
      t188 = (t87+t89+Zoi-t125+t126-Zoj)**2
      t190 = sqrt(t184+t186+t188)
      t194 = qb1i*qa1j
      t195 = t43*bc
      t197 = t41*t85*bs
      t198 = t195+t197+Xoi-t59-t64-Xoj
      t199 = t198**2
      t200 = t69*bc
      t202 = t35*t85*bs
      t204 = (t200+t202+Yoi-t78-t82-Yoj)**2
      t205 = t86*bc
      t206 = t38*bs
      t208 = (t205+t206+Zoi-t92+t94-Zoj)**2
      t210 = sqrt(t199+t204+t208)
      t214 = qb1i*qa2j
      t215 = t195+t197+Xoi-t59+t64-Xoj
      t216 = t215**2
      t218 = (t200+t202+Yoi-t78+t82-Yoj)**2
      t220 = (t205+t206+Zoi-t92-t94-Zoj)**2
      t222 = sqrt(t216+t218+t220)
      t226 = qb1i*qb1j
      t227 = t195+t197+Xoi-t115-t117-Xoj
      t228 = t227**2
      t230 = (t200+t202+Yoi-t120-t122-Yoj)**2
      t232 = (t205+t206+Zoi-t125-t126-Zoj)**2
      t234 = sqrt(t228+t230+t232)
      t238 = qb1i*qb2j
      t239 = t195+t197+Xoi-t115+t117-Xoj
      t240 = t239**2
      t242 = (t200+t202+Yoi-t120+t122-Yoj)**2
      t244 = (t205+t206+Zoi-t125+t126-Zoj)**2
      t246 = sqrt(t240+t242+t244)
      t250 = qb2i*qa1j
      t251 = t195-t197+Xoi-t59-t64-Xoj
      t252 = t251**2
      t254 = (t200-t202+Yoi-t78-t82-Yoj)**2
      t256 = (t205-t206+Zoi-t92+t94-Zoj)**2
      t258 = sqrt(t252+t254+t256)
      t262 = qb2i*qa2j
      t263 = t195-t197+Xoi-t59+t64-Xoj
      t264 = t263**2
      t266 = (t200-t202+Yoi-t78+t82-Yoj)**2
      t268 = (t205-t206+Zoi-t92-t94-Zoj)**2
      t270 = sqrt(t264+t266+t268)
      t274 = qb2i*qb1j
      t275 = t195-t197+Xoi-t115-t117-Xoj
      t276 = t275**2
      t278 = (t200-t202+Yoi-t120-t122-Yoj)**2
      t280 = (t205-t206+Zoi-t125-t126-Zoj)**2
      t282 = sqrt(t276+t278+t280)
      t286 = qb2i*qb2j
      t287 = t195-t197+Xoi-t115+t117-Xoj
      t288 = t287**2
      t290 = (t200-t202+Yoi-t120+t122-Yoj)**2
      t292 = (t205-t206+Zoi-t125+t126-Zoj)**2
      t294 = sqrt(t288+t290+t292)
      t298 = t34*ec/t98+t102*ec/t110+t114*ec/t130+t134*ec/t142+t146*ec/t
     #154+t158*ec/t166+t170*ec/t178+t182*ec/t190+t194*ec/t210+t214*ec/t2
     #22+t226*ec/t234+t238*ec/t246+t250*ec/t258+t262*ec/t270+t274*ec/t28
     #2+t286*ec/t294
      t303 = t98**2
      t309 = t110**2
      t315 = t130**2
      t321 = t142**2
      t327 = t154**2
      t333 = t166**2
      t339 = t178**2
      t345 = t190**2
      t351 = t210**2
      t357 = t222**2
      t363 = t234**2
      t369 = t246**2
      t375 = t258**2
      t381 = t270**2
      t387 = t282**2
      t393 = t294**2
      t399 = -2*t34*ec/t303/t98*t65-2*t102*ec/t309/t110*t103-2*t114*ec/t
     #315/t130*t118-2*t134*ec/t321/t142*t135-2*t146*ec/t327/t154*t147-2*
     #t158*ec/t333/t166*t159-2*t170*ec/t339/t178*t171-2*t182*ec/t345/t19
     #0*t183-2*t194*ec/t351/t210*t198-2*t214*ec/t357/t222*t215-2*t226*ec
     #/t363/t234*t227-2*t238*ec/t369/t246*t239-2*t250*ec/t375/t258*t251-
     #2*t262*ec/t381/t270*t263-2*t274*ec/t387/t282*t275-2*t286*ec/t393/t
     #294*t287
      t401 = t13**2
      t403 = t401**2
      t413 = (cst3*t15*t19*t28-t31*t28)*t298+(t31*t19*t24+cst4)*t399/2-1
     #2*AA/t403/t401/t13*t27+6*BB/t403*t27
c
      dVdXoi=t413
c
      return
      end
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      function dVdYoi(Xoi,Xoj,Yoi,Yoj,Zoi,Zoj,ppi,si,ti,ppj,sj,tj)
      implicit real*8(a-h,o-z)
      common/cpot1/qa1i,qa2i,qb1i,qb2i,qa1j,qa2j,qb1j,qb2j
      common/cpot2/AA,BB,ac,as,bc,bs,ec
      common/cpot3/cst1,cst2,cst3,cst4
c
      t1 = Xoi**2
      t4 = Xoj**2
      t5 = Yoi**2
      t8 = Yoj**2
      t9 = Zoi**2
      t12 = Zoj**2
      t13 = t1-2*Xoi*Xoj+t4+t5-2*Yoi*Yoj+t8+t9-2*Zoi*Zoj+t12
      t14 = sqrt(t13)
      t15 = t14-cst1
      t19 = 3*cst2-cst1-2*t14
      t21 = cst2-cst1
      t22 = t21**2
      t24 = 1/t22/t21
      t27 = Yoi-Yoj
      t28 = 2*t24/t14*t27
      t30 = t15**2
      t31 = cst3*t30
      t34 = qa1i*qa1j
      t35 = cos(ppi)
      t36 = cos(si)
      t38 = cos(ti)
      t39 = sin(si)
      t40 = t38*t39
      t41 = sin(ppi)
      t43 = t35*t36-t40*t41
      t44 = t43*ac
      t46 = t38*t36
      t49 = (t35*t39+t46*t41)*as
      t50 = cos(ppj)
      t51 = cos(sj)
      t53 = cos(tj)
      t54 = sin(sj)
      t55 = t53*t54
      t56 = sin(ppj)
      t58 = t50*t51-t55*t56
      t59 = t58*ac
      t61 = t53*t51
      t64 = (t50*t54+t61*t56)*as
      t66 = (t44+t49+Xoi-t59-t64-Xoj)**2
      t69 = -t41*t36-t40*t35
      t70 = t69*ac
      t74 = (-t41*t39+t46*t35)*as
      t77 = -t56*t51-t55*t50
      t78 = t77*ac
      t82 = (-t56*t54+t61*t50)*as
      t83 = t70+t74+Yoi-t78-t82-Yoj
      t84 = t83**2
      t85 = sin(ti)
      t86 = t85*t39
      t87 = t86*ac
      t89 = t85*t36*as
      t90 = sin(tj)
      t91 = t90*t54
      t92 = t91*ac
      t94 = t90*t51*as
      t96 = (t87-t89+Zoi-t92+t94-Zoj)**2
      t98 = sqrt(t66+t84+t96)
      t102 = qa1i*qa2j
      t104 = (t44+t49+Xoi-t59+t64-Xoj)**2
      t105 = t70+t74+Yoi-t78+t82-Yoj
      t106 = t105**2
      t108 = (t87-t89+Zoi-t92-t94-Zoj)**2
      t110 = sqrt(t104+t106+t108)
      t114 = qa1i*qb1j
      t115 = t58*bc
      t117 = t56*t90*bs
      t119 = (t44+t49+Xoi-t115-t117-Xoj)**2
      t120 = t77*bc
      t122 = t50*t90*bs
      t123 = t70+t74+Yoi-t120-t122-Yoj
      t124 = t123**2
      t125 = t91*bc
      t126 = t53*bs
      t128 = (t87-t89+Zoi-t125-t126-Zoj)**2
      t130 = sqrt(t119+t124+t128)
      t134 = qa1i*qb2j
      t136 = (t44+t49+Xoi-t115+t117-Xoj)**2
      t137 = t70+t74+Yoi-t120+t122-Yoj
      t138 = t137**2
      t140 = (t87-t89+Zoi-t125+t126-Zoj)**2
      t142 = sqrt(t136+t138+t140)
      t146 = qa2i*qa1j
      t148 = (t44-t49+Xoi-t59-t64-Xoj)**2
      t149 = t70-t74+Yoi-t78-t82-Yoj
      t150 = t149**2
      t152 = (t87+t89+Zoi-t92+t94-Zoj)**2
      t154 = sqrt(t148+t150+t152)
      t158 = qa2i*qa2j
      t160 = (t44-t49+Xoi-t59+t64-Xoj)**2
      t161 = t70-t74+Yoi-t78+t82-Yoj
      t162 = t161**2
      t164 = (t87+t89+Zoi-t92-t94-Zoj)**2
      t166 = sqrt(t160+t162+t164)
      t170 = qa2i*qb1j
      t172 = (t44-t49+Xoi-t115-t117-Xoj)**2
      t173 = t70-t74+Yoi-t120-t122-Yoj
      t174 = t173**2
      t176 = (t87+t89+Zoi-t125-t126-Zoj)**2
      t178 = sqrt(t172+t174+t176)
      t182 = qa2i*qb2j
      t184 = (t44-t49+Xoi-t115+t117-Xoj)**2
      t185 = t70-t74+Yoi-t120+t122-Yoj
      t186 = t185**2
      t188 = (t87+t89+Zoi-t125+t126-Zoj)**2
      t190 = sqrt(t184+t186+t188)
      t194 = qb1i*qa1j
      t195 = t43*bc
      t197 = t41*t85*bs
      t199 = (t195+t197+Xoi-t59-t64-Xoj)**2
      t200 = t69*bc
      t202 = t35*t85*bs
      t203 = t200+t202+Yoi-t78-t82-Yoj
      t204 = t203**2
      t205 = t86*bc
      t206 = t38*bs
      t208 = (t205+t206+Zoi-t92+t94-Zoj)**2
      t210 = sqrt(t199+t204+t208)
      t214 = qb1i*qa2j
      t216 = (t195+t197+Xoi-t59+t64-Xoj)**2
      t217 = t200+t202+Yoi-t78+t82-Yoj
      t218 = t217**2
      t220 = (t205+t206+Zoi-t92-t94-Zoj)**2
      t222 = sqrt(t216+t218+t220)
      t226 = qb1i*qb1j
      t228 = (t195+t197+Xoi-t115-t117-Xoj)**2
      t229 = t200+t202+Yoi-t120-t122-Yoj
      t230 = t229**2
      t232 = (t205+t206+Zoi-t125-t126-Zoj)**2
      t234 = sqrt(t228+t230+t232)
      t238 = qb1i*qb2j
      t240 = (t195+t197+Xoi-t115+t117-Xoj)**2
      t241 = t200+t202+Yoi-t120+t122-Yoj
      t242 = t241**2
      t244 = (t205+t206+Zoi-t125+t126-Zoj)**2
      t246 = sqrt(t240+t242+t244)
      t250 = qb2i*qa1j
      t252 = (t195-t197+Xoi-t59-t64-Xoj)**2
      t253 = t200-t202+Yoi-t78-t82-Yoj
      t254 = t253**2
      t256 = (t205-t206+Zoi-t92+t94-Zoj)**2
      t258 = sqrt(t252+t254+t256)
      t262 = qb2i*qa2j
      t264 = (t195-t197+Xoi-t59+t64-Xoj)**2
      t265 = t200-t202+Yoi-t78+t82-Yoj
      t266 = t265**2
      t268 = (t205-t206+Zoi-t92-t94-Zoj)**2
      t270 = sqrt(t264+t266+t268)
      t274 = qb2i*qb1j
      t276 = (t195-t197+Xoi-t115-t117-Xoj)**2
      t277 = t200-t202+Yoi-t120-t122-Yoj
      t278 = t277**2
      t280 = (t205-t206+Zoi-t125-t126-Zoj)**2
      t282 = sqrt(t276+t278+t280)
      t286 = qb2i*qb2j
      t288 = (t195-t197+Xoi-t115+t117-Xoj)**2
      t289 = t200-t202+Yoi-t120+t122-Yoj
      t290 = t289**2
      t292 = (t205-t206+Zoi-t125+t126-Zoj)**2
      t294 = sqrt(t288+t290+t292)
      t298 = t34*ec/t98+t102*ec/t110+t114*ec/t130+t134*ec/t142+t146*ec/t
     #154+t158*ec/t166+t170*ec/t178+t182*ec/t190+t194*ec/t210+t214*ec/t2
     #22+t226*ec/t234+t238*ec/t246+t250*ec/t258+t262*ec/t270+t274*ec/t28
     #2+t286*ec/t294
      t303 = t98**2
      t309 = t110**2
      t315 = t130**2
      t321 = t142**2
      t327 = t154**2
      t333 = t166**2
      t339 = t178**2
      t345 = t190**2
      t351 = t210**2
      t357 = t222**2
      t363 = t234**2
      t369 = t246**2
      t375 = t258**2
      t381 = t270**2
      t387 = t282**2
      t393 = t294**2
      t399 = -2*t34*ec/t303/t98*t83-2*t102*ec/t309/t110*t105-2*t114*ec/t
     #315/t130*t123-2*t134*ec/t321/t142*t137-2*t146*ec/t327/t154*t149-2*
     #t158*ec/t333/t166*t161-2*t170*ec/t339/t178*t173-2*t182*ec/t345/t19
     #0*t185-2*t194*ec/t351/t210*t203-2*t214*ec/t357/t222*t217-2*t226*ec
     #/t363/t234*t229-2*t238*ec/t369/t246*t241-2*t250*ec/t375/t258*t253-
     #2*t262*ec/t381/t270*t265-2*t274*ec/t387/t282*t277-2*t286*ec/t393/t
     #294*t289
      t401 = t13**2
      t403 = t401**2
      t413 = (cst3*t15*t19*t28-t31*t28)*t298+(t31*t19*t24+cst4)*t399/2-1
     #2*AA/t403/t401/t13*t27+6*BB/t403*t27
c
      dVdYoi=t413
c
      return
      end   
   
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      function dVdZoi(Xoi,Xoj,Yoi,Yoj,Zoi,Zoj,ppi,si,ti,ppj,sj,tj)
      implicit real*8(a-h,o-z)
      common/cpot1/qa1i,qa2i,qb1i,qb2i,qa1j,qa2j,qb1j,qb2j
      common/cpot2/AA,BB,ac,as,bc,bs,ec   
      common/cpot3/cst1,cst2,cst3,cst4   
c
      t1 = Xoi**2
      t4 = Xoj**2
      t5 = Yoi**2
      t8 = Yoj**2
      t9 = Zoi**2
      t12 = Zoj**2
      t13 = t1-2*Xoi*Xoj+t4+t5-2*Yoi*Yoj+t8+t9-2*Zoi*Zoj+t12
      t14 = sqrt(t13)
      t15 = t14-cst1
      t19 = 3*cst2-cst1-2*t14
      t21 = cst2-cst1
      t22 = t21**2
      t24 = 1/t22/t21
      t27 = Zoi-Zoj
      t28 = 2*t24/t14*t27
      t30 = t15**2
      t31 = cst3*t30
      t34 = qa1i*qa1j
      t35 = cos(ppi)
      t36 = cos(si)
      t38 = cos(ti)
      t39 = sin(si)
      t40 = t38*t39
      t41 = sin(ppi)
      t43 = t35*t36-t40*t41
      t44 = t43*ac
      t46 = t38*t36
      t49 = (t35*t39+t46*t41)*as
      t50 = cos(ppj)
      t51 = cos(sj)
      t53 = cos(tj)
      t54 = sin(sj)
      t55 = t53*t54
      t56 = sin(ppj)
      t58 = t50*t51-t55*t56
      t59 = t58*ac
      t61 = t53*t51
      t64 = (t50*t54+t61*t56)*as
      t66 = (t44+t49+Xoi-t59-t64-Xoj)**2
      t69 = -t41*t36-t40*t35
      t70 = t69*ac
      t74 = (-t41*t39+t46*t35)*as
      t77 = -t56*t51-t55*t50
      t78 = t77*ac
      t82 = (-t56*t54+t61*t50)*as
      t84 = (t70+t74+Yoi-t78-t82-Yoj)**2
      t85 = sin(ti)
      t86 = t85*t39
      t87 = t86*ac
      t89 = t85*t36*as
      t90 = sin(tj)
      t91 = t90*t54
      t92 = t91*ac
      t94 = t90*t51*as
      t95 = t87-t89+Zoi-t92+t94-Zoj
      t96 = t95**2
      t98 = sqrt(t66+t84+t96)
      t102 = qa1i*qa2j
      t104 = (t44+t49+Xoi-t59+t64-Xoj)**2
      t106 = (t70+t74+Yoi-t78+t82-Yoj)**2
      t107 = t87-t89+Zoi-t92-t94-Zoj
      t108 = t107**2
      t110 = sqrt(t104+t106+t108)
      t114 = qa1i*qb1j
      t115 = t58*bc
      t117 = t56*t90*bs
      t119 = (t44+t49+Xoi-t115-t117-Xoj)**2
      t120 = t77*bc
      t122 = t50*t90*bs
      t124 = (t70+t74+Yoi-t120-t122-Yoj)**2
      t125 = t91*bc
      t126 = t53*bs
      t127 = t87-t89+Zoi-t125-t126-Zoj
      t128 = t127**2
      t130 = sqrt(t119+t124+t128)
      t134 = qa1i*qb2j
      t136 = (t44+t49+Xoi-t115+t117-Xoj)**2
      t138 = (t70+t74+Yoi-t120+t122-Yoj)**2
      t139 = t87-t89+Zoi-t125+t126-Zoj
      t140 = t139**2
      t142 = sqrt(t136+t138+t140)
      t146 = qa2i*qa1j
      t148 = (t44-t49+Xoi-t59-t64-Xoj)**2
      t150 = (t70-t74+Yoi-t78-t82-Yoj)**2
      t151 = t87+t89+Zoi-t92+t94-Zoj
      t152 = t151**2
      t154 = sqrt(t148+t150+t152)
      t158 = qa2i*qa2j
      t160 = (t44-t49+Xoi-t59+t64-Xoj)**2
      t162 = (t70-t74+Yoi-t78+t82-Yoj)**2
      t163 = t87+t89+Zoi-t92-t94-Zoj
      t164 = t163**2
      t166 = sqrt(t160+t162+t164)
      t170 = qa2i*qb1j
      t172 = (t44-t49+Xoi-t115-t117-Xoj)**2
      t174 = (t70-t74+Yoi-t120-t122-Yoj)**2
      t175 = t87+t89+Zoi-t125-t126-Zoj
      t176 = t175**2
      t178 = sqrt(t172+t174+t176)
      t182 = qa2i*qb2j
      t184 = (t44-t49+Xoi-t115+t117-Xoj)**2
      t186 = (t70-t74+Yoi-t120+t122-Yoj)**2
      t187 = t87+t89+Zoi-t125+t126-Zoj
      t188 = t187**2
      t190 = sqrt(t184+t186+t188)
      t194 = qb1i*qa1j
      t195 = t43*bc
      t197 = t41*t85*bs
      t199 = (t195+t197+Xoi-t59-t64-Xoj)**2
      t200 = t69*bc
      t202 = t35*t85*bs
      t204 = (t200+t202+Yoi-t78-t82-Yoj)**2
      t205 = t86*bc
      t206 = t38*bs
      t207 = t205+t206+Zoi-t92+t94-Zoj
      t208 = t207**2
      t210 = sqrt(t199+t204+t208)
      t214 = qb1i*qa2j
      t216 = (t195+t197+Xoi-t59+t64-Xoj)**2
      t218 = (t200+t202+Yoi-t78+t82-Yoj)**2
      t219 = t205+t206+Zoi-t92-t94-Zoj
      t220 = t219**2
      t222 = sqrt(t216+t218+t220)
      t226 = qb1i*qb1j
      t228 = (t195+t197+Xoi-t115-t117-Xoj)**2
      t230 = (t200+t202+Yoi-t120-t122-Yoj)**2
      t231 = t205+t206+Zoi-t125-t126-Zoj
      t232 = t231**2
      t234 = sqrt(t228+t230+t232)
      t238 = qb1i*qb2j
      t240 = (t195+t197+Xoi-t115+t117-Xoj)**2
      t242 = (t200+t202+Yoi-t120+t122-Yoj)**2
      t243 = t205+t206+Zoi-t125+t126-Zoj
      t244 = t243**2
      t246 = sqrt(t240+t242+t244)
      t250 = qb2i*qa1j
      t252 = (t195-t197+Xoi-t59-t64-Xoj)**2
      t254 = (t200-t202+Yoi-t78-t82-Yoj)**2
      t255 = t205-t206+Zoi-t92+t94-Zoj
      t256 = t255**2
      t258 = sqrt(t252+t254+t256)
      t262 = qb2i*qa2j
      t264 = (t195-t197+Xoi-t59+t64-Xoj)**2
      t266 = (t200-t202+Yoi-t78+t82-Yoj)**2
      t267 = t205-t206+Zoi-t92-t94-Zoj
      t268 = t267**2
      t270 = sqrt(t264+t266+t268)
      t274 = qb2i*qb1j
      t276 = (t195-t197+Xoi-t115-t117-Xoj)**2
      t278 = (t200-t202+Yoi-t120-t122-Yoj)**2
      t279 = t205-t206+Zoi-t125-t126-Zoj
      t280 = t279**2
      t282 = sqrt(t276+t278+t280)
      t286 = qb2i*qb2j
      t288 = (t195-t197+Xoi-t115+t117-Xoj)**2
      t290 = (t200-t202+Yoi-t120+t122-Yoj)**2
      t291 = t205-t206+Zoi-t125+t126-Zoj
      t292 = t291**2
      t294 = sqrt(t288+t290+t292)
      t298 = t34*ec/t98+t102*ec/t110+t114*ec/t130+t134*ec/t142+t146*ec/t
     #154+t158*ec/t166+t170*ec/t178+t182*ec/t190+t194*ec/t210+t214*ec/t2
     #22+t226*ec/t234+t238*ec/t246+t250*ec/t258+t262*ec/t270+t274*ec/t28
     #2+t286*ec/t294
      t303 = t98**2
      t309 = t110**2
      t315 = t130**2
      t321 = t142**2
      t327 = t154**2
      t333 = t166**2
      t339 = t178**2
      t345 = t190**2
      t351 = t210**2
      t357 = t222**2
      t363 = t234**2
      t369 = t246**2
      t375 = t258**2
      t381 = t270**2
      t387 = t282**2
      t393 = t294**2
      t399 = -2*t34*ec/t303/t98*t95-2*t102*ec/t309/t110*t107-2*t114*ec/t
     #315/t130*t127-2*t134*ec/t321/t142*t139-2*t146*ec/t327/t154*t151-2*
     #t158*ec/t333/t166*t163-2*t170*ec/t339/t178*t175-2*t182*ec/t345/t19
     #0*t187-2*t194*ec/t351/t210*t207-2*t214*ec/t357/t222*t219-2*t226*ec
     #/t363/t234*t231-2*t238*ec/t369/t246*t243-2*t250*ec/t375/t258*t255-
     #2*t262*ec/t381/t270*t267-2*t274*ec/t387/t282*t279-2*t286*ec/t393/t
     #294*t291
      t401 = t13**2
      t403 = t401**2
      t413 = (cst3*t15*t19*t28-t31*t28)*t298+(t31*t19*t24+cst4)*t399/2-1
     #2*AA/t403/t401/t13*t27+6*BB/t403*t27
c
      dVdZoi=t413
c
      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      function dVdti(Xoi,Xoj,Yoi,Yoj,Zoi,Zoj,ppi,si,ti,ppj,sj,tj)
      implicit real*8(a-h,o-z)
      common/cpot1/qa1i,qa2i,qb1i,qb2i,qa1j,qa2j,qb1j,qb2j
      common/cpot2/AA,BB,ac,as,bc,bs,ec
      common/cpot3/cst1,cst2,cst3,cst4
c
      t1 = Xoi**2
      t4 = Xoj**2
      t5 = Yoi**2
      t8 = Yoj**2
      t9 = Zoi**2
      t12 = Zoj**2
      t14 = sqrt(t1-2*Xoi*Xoj+t4+t5-2*Yoi*Yoj+t8+t9-2*Zoi*Zoj+t12)
      t16 = (t14-cst1)**2
      t21 = cst2-cst1
      t22 = t21**2
      t29 = cos(ppi)
      t30 = cos(si)
      t32 = cos(ti)
      t33 = sin(si)
      t34 = t32*t33
      t35 = sin(ppi)
      t37 = t29*t30-t34*t35
      t38 = t37*ac
      t40 = t32*t30
      t43 = (t29*t33+t40*t35)*as
      t44 = cos(ppj)
      t45 = cos(sj)
      t47 = cos(tj)
      t48 = sin(sj)
      t49 = t47*t48
      t50 = sin(ppj)
      t52 = t44*t45-t49*t50
      t53 = t52*ac
      t55 = t47*t45
      t58 = (t44*t48+t55*t50)*as
      t59 = t38+t43+Xoi-t53-t58-Xoj
      t60 = t59**2
      t63 = -t35*t30-t34*t29
      t64 = t63*ac
      t68 = (-t35*t33+t40*t29)*as
      t71 = -t50*t45-t49*t44
      t72 = t71*ac
      t76 = (-t50*t48+t55*t44)*as
      t77 = t64+t68+Yoi-t72-t76-Yoj
      t78 = t77**2
      t79 = sin(ti)
      t80 = t79*t33
      t81 = t80*ac
      t82 = t79*t30
      t83 = t82*as
      t84 = sin(tj)
      t85 = t84*t48
      t86 = t85*ac
      t88 = t84*t45*as
      t89 = t81-t83+Zoi-t86+t88-Zoj
      t90 = t89**2
      t92 = sqrt(t60+t78+t90)
      t93 = t92**2
      t98 = t80*t35*ac
      t100 = t82*t35*as
      t101 = t98-t100
      t104 = t80*t29*ac
      t106 = t82*t29*as
      t107 = t104-t106
      t109 = t34*ac
      t110 = t40*as
      t111 = t109-t110
      t117 = t38+t43+Xoi-t53+t58-Xoj
      t118 = t117**2
      t119 = t64+t68+Yoi-t72+t76-Yoj
      t120 = t119**2
      t121 = t81-t83+Zoi-t86-t88-Zoj
      t122 = t121**2
      t124 = sqrt(t118+t120+t122)
      t125 = t124**2
      t136 = t52*bc
      t138 = t50*t84*bs
      t139 = t38+t43+Xoi-t136-t138-Xoj
      t140 = t139**2
      t141 = t71*bc
      t143 = t44*t84*bs
      t144 = t64+t68+Yoi-t141-t143-Yoj
      t145 = t144**2
      t146 = t85*bc
      t147 = t47*bs
      t148 = t81-t83+Zoi-t146-t147-Zoj
      t149 = t148**2
      t151 = sqrt(t140+t145+t149)
      t152 = t151**2
      t163 = t38+t43+Xoi-t136+t138-Xoj
      t164 = t163**2
      t165 = t64+t68+Yoi-t141+t143-Yoj
      t166 = t165**2
      t167 = t81-t83+Zoi-t146+t147-Zoj
      t168 = t167**2
      t170 = sqrt(t164+t166+t168)
      t171 = t170**2
      t182 = t38-t43+Xoi-t53-t58-Xoj
      t183 = t182**2
      t184 = t64-t68+Yoi-t72-t76-Yoj
      t185 = t184**2
      t186 = t81+t83+Zoi-t86+t88-Zoj
      t187 = t186**2
      t189 = sqrt(t183+t185+t187)
      t190 = t189**2
      t194 = t98+t100
      t196 = t104+t106
      t198 = t109+t110
      t204 = t38-t43+Xoi-t53+t58-Xoj
      t205 = t204**2
      t206 = t64-t68+Yoi-t72+t76-Yoj
      t207 = t206**2
      t208 = t81+t83+Zoi-t86-t88-Zoj
      t209 = t208**2
      t211 = sqrt(t205+t207+t209)
      t212 = t211**2
      t223 = t38-t43+Xoi-t136-t138-Xoj
      t224 = t223**2
      t225 = t64-t68+Yoi-t141-t143-Yoj
      t226 = t225**2
      t227 = t81+t83+Zoi-t146-t147-Zoj
      t228 = t227**2
      t230 = sqrt(t224+t226+t228)
      t231 = t230**2
      t242 = t38-t43+Xoi-t136+t138-Xoj
      t243 = t242**2
      t244 = t64-t68+Yoi-t141+t143-Yoj
      t245 = t244**2
      t246 = t81+t83+Zoi-t146+t147-Zoj
      t247 = t246**2
      t249 = sqrt(t243+t245+t247)
      t250 = t249**2
      t261 = t37*bc
      t263 = t35*t79*bs
      t264 = t261+t263+Xoi-t53-t58-Xoj
      t265 = t264**2
      t266 = t63*bc
      t268 = t29*t79*bs
      t269 = t266+t268+Yoi-t72-t76-Yoj
      t270 = t269**2
      t271 = t80*bc
      t272 = t32*bs
      t273 = t271+t272+Zoi-t86+t88-Zoj
      t274 = t273**2
      t276 = sqrt(t265+t270+t274)
      t277 = t276**2
      t282 = t80*t35*bc
      t284 = t35*t32*bs
      t285 = t282+t284
      t288 = t80*t29*bc
      t290 = t29*t32*bs
      t291 = t288+t290
      t293 = t34*bc
      t294 = t79*bs
      t295 = t293-t294
      t301 = t261+t263+Xoi-t53+t58-Xoj
      t302 = t301**2
      t303 = t266+t268+Yoi-t72+t76-Yoj
      t304 = t303**2
      t305 = t271+t272+Zoi-t86-t88-Zoj
      t306 = t305**2
      t308 = sqrt(t302+t304+t306)
      t309 = t308**2
      t320 = t261+t263+Xoi-t136-t138-Xoj
      t321 = t320**2
      t322 = t266+t268+Yoi-t141-t143-Yoj
      t323 = t322**2
      t324 = t271+t272+Zoi-t146-t147-Zoj
      t325 = t324**2
      t327 = sqrt(t321+t323+t325)
      t328 = t327**2
      t339 = t261+t263+Xoi-t136+t138-Xoj
      t340 = t339**2
      t341 = t266+t268+Yoi-t141+t143-Yoj
      t342 = t341**2
      t343 = t271+t272+Zoi-t146+t147-Zoj
      t344 = t343**2
      t346 = sqrt(t340+t342+t344)
      t347 = t346**2
      t358 = t261-t263+Xoi-t53-t58-Xoj
      t359 = t358**2
      t360 = t266-t268+Yoi-t72-t76-Yoj
      t361 = t360**2
      t362 = t271-t272+Zoi-t86+t88-Zoj
      t363 = t362**2
      t365 = sqrt(t359+t361+t363)
      t366 = t365**2
      t370 = t282-t284
      t372 = t288-t290
      t374 = t293+t294
      t380 = t261-t263+Xoi-t53+t58-Xoj
      t381 = t380**2
      t382 = t266-t268+Yoi-t72+t76-Yoj
      t383 = t382**2
      t384 = t271-t272+Zoi-t86-t88-Zoj
      t385 = t384**2
      t387 = sqrt(t381+t383+t385)
      t388 = t387**2
      t399 = t261-t263+Xoi-t136-t138-Xoj
      t400 = t399**2
      t401 = t266-t268+Yoi-t141-t143-Yoj
      t402 = t401**2
      t403 = t271-t272+Zoi-t146-t147-Zoj
      t404 = t403**2
      t406 = sqrt(t400+t402+t404)
      t407 = t406**2
      t418 = t261-t263+Xoi-t136+t138-Xoj
      t419 = t418**2
      t420 = t266-t268+Yoi-t141+t143-Yoj
      t421 = t420**2
      t422 = t271-t272+Zoi-t146+t147-Zoj
      t423 = t422**2
      t425 = sqrt(t419+t421+t423)
      t426 = t425**2
      s1 = -2*qa1i*qa1j*ec/t93/t92*(t59*t101+t77*t107+t89*t111)-2*qa1i*q
     #a2j*ec/t125/t124*(t117*t101+t119*t107+t121*t111)-2*qa1i*qb1j*ec/t1
     #52/t151*(t139*t101+t144*t107+t148*t111)-2*qa1i*qb2j*ec/t171/t170*(
     #t163*t101+t165*t107+t167*t111)-2*qa2i*qa1j*ec/t190/t189*(t182*t194
     #+t184*t196+t186*t198)-2*qa2i*qa2j*ec/t212/t211*(t204*t194+t206*t19
     #6+t208*t198)-2*qa2i*qb1j*ec/t231/t230*(t223*t194+t225*t196+t227*t1
     #98)-2*qa2i*qb2j*ec/t250/t249*(t242*t194+t244*t196+t246*t198)
      t436 = s1-2*qb1i*qa1j*ec/t277/t276*(t264*t285+t269*t291+t273*t295)
     #-2*qb1i*qa2j*ec/t309/t308*(t301*t285+t303*t291+t305*t295)-2*qb1i*q
     #b1j*ec/t328/t327*(t320*t285+t322*t291+t324*t295)-2*qb1i*qb2j*ec/t3
     #47/t346*(t339*t285+t341*t291+t343*t295)-2*qb2i*qa1j*ec/t366/t365*(
     #t358*t370+t360*t372+t362*t374)-2*qb2i*qa2j*ec/t388/t387*(t380*t370
     #+t382*t372+t384*t374)-2*qb2i*qb1j*ec/t407/t406*(t399*t370+t401*t37
     #2+t403*t374)-2*qb2i*qb2j*ec/t426/t425*(t418*t370+t420*t372+t422*t3
     #74)
      t437 = (cst3*t16*(3*cst2-cst1-2*t14)/t22/t21+cst4)*t436/2
c
      dVdti=t437
c
      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      function dVdsi(Xoi,Xoj,Yoi,Yoj,Zoi,Zoj,ppi,si,ti,ppj,sj,tj)
      implicit real*8(a-h,o-z)
      common/cpot1/qa1i,qa2i,qb1i,qb2i,qa1j,qa2j,qb1j,qb2j
      common/cpot2/AA,BB,ac,as,bc,bs,ec
      common/cpot3/cst1,cst2,cst3,cst4
c
      t1 = Xoi**2
      t4 = Xoj**2
      t5 = Yoi**2
      t8 = Yoj**2
      t9 = Zoi**2
      t12 = Zoj**2
      t14 = sqrt(t1-2*Xoi*Xoj+t4+t5-2*Yoi*Yoj+t8+t9-2*Zoi*Zoj+t12)
      t16 = (t14-cst1)**2
      t21 = cst2-cst1
      t22 = t21**2
      t29 = cos(ppi)
      t30 = cos(si)
      t32 = cos(ti)
      t33 = sin(si)
      t34 = t32*t33
      t35 = sin(ppi)
      t37 = t29*t30-t34*t35
      t38 = t37*ac
      t40 = t32*t30
      t42 = t29*t33+t40*t35
      t43 = t42*as
      t44 = cos(ppj)
      t45 = cos(sj)
      t47 = cos(tj)
      t48 = sin(sj)
      t49 = t47*t48
      t50 = sin(ppj)
      t52 = t44*t45-t49*t50
      t53 = t52*ac
      t55 = t47*t45
      t58 = (t44*t48+t55*t50)*as
      t59 = t38+t43+Xoi-t53-t58-Xoj
      t60 = t59**2
      t63 = -t35*t30-t34*t29
      t64 = t63*ac
      t67 = -t35*t33+t40*t29
      t68 = t67*as
      t71 = -t50*t45-t49*t44
      t72 = t71*ac
      t76 = (-t50*t48+t55*t44)*as
      t77 = t64+t68+Yoi-t72-t76-Yoj
      t78 = t77**2
      t79 = sin(ti)
      t80 = t79*t33
      t81 = t80*ac
      t82 = t79*t30
      t83 = t82*as
      t84 = sin(tj)
      t85 = t84*t48
      t86 = t85*ac
      t88 = t84*t45*as
      t89 = t81-t83+Zoi-t86+t88-Zoj
      t90 = t89**2
      t92 = sqrt(t60+t78+t90)
      t93 = t92**2
      t97 = -t42*ac
      t98 = t37*as
      t99 = t97+t98
      t101 = -t67*ac
      t102 = t63*as
      t103 = t101+t102
      t105 = t82*ac
      t106 = t80*as
      t107 = t105+t106
      t113 = t38+t43+Xoi-t53+t58-Xoj
      t114 = t113**2
      t115 = t64+t68+Yoi-t72+t76-Yoj
      t116 = t115**2
      t117 = t81-t83+Zoi-t86-t88-Zoj
      t118 = t117**2
      t120 = sqrt(t114+t116+t118)
      t121 = t120**2
      t132 = t52*bc
      t134 = t50*t84*bs
      t135 = t38+t43+Xoi-t132-t134-Xoj
      t136 = t135**2
      t137 = t71*bc
      t139 = t44*t84*bs
      t140 = t64+t68+Yoi-t137-t139-Yoj
      t141 = t140**2
      t142 = t85*bc
      t143 = t47*bs
      t144 = t81-t83+Zoi-t142-t143-Zoj
      t145 = t144**2
      t147 = sqrt(t136+t141+t145)
      t148 = t147**2
      t159 = t38+t43+Xoi-t132+t134-Xoj
      t160 = t159**2
      t161 = t64+t68+Yoi-t137+t139-Yoj
      t162 = t161**2
      t163 = t81-t83+Zoi-t142+t143-Zoj
      t164 = t163**2
      t166 = sqrt(t160+t162+t164)
      t167 = t166**2
      t178 = t38-t43+Xoi-t53-t58-Xoj
      t179 = t178**2
      t180 = t64-t68+Yoi-t72-t76-Yoj
      t181 = t180**2
      t182 = t81+t83+Zoi-t86+t88-Zoj
      t183 = t182**2
      t185 = sqrt(t179+t181+t183)
      t186 = t185**2
      t190 = t97-t98
      t192 = t101-t102
      t194 = t105-t106
      t200 = t38-t43+Xoi-t53+t58-Xoj
      t201 = t200**2
      t202 = t64-t68+Yoi-t72+t76-Yoj
      t203 = t202**2
      t204 = t81+t83+Zoi-t86-t88-Zoj
      t205 = t204**2
      t207 = sqrt(t201+t203+t205)
      t208 = t207**2
      t219 = t38-t43+Xoi-t132-t134-Xoj
      t220 = t219**2
      t221 = t64-t68+Yoi-t137-t139-Yoj
      t222 = t221**2
      t223 = t81+t83+Zoi-t142-t143-Zoj
      t224 = t223**2
      t226 = sqrt(t220+t222+t224)
      t227 = t226**2
      t238 = t38-t43+Xoi-t132+t134-Xoj
      t239 = t238**2
      t240 = t64-t68+Yoi-t137+t139-Yoj
      t241 = t240**2
      t242 = t81+t83+Zoi-t142+t143-Zoj
      t243 = t242**2
      t245 = sqrt(t239+t241+t243)
      t246 = t245**2
      t257 = t37*bc
      t259 = t35*t79*bs
      t260 = t257+t259+Xoi-t53-t58-Xoj
      t261 = t260**2
      t262 = t63*bc
      t264 = t29*t79*bs
      t265 = t262+t264+Yoi-t72-t76-Yoj
      t266 = t265**2
      t267 = t80*bc
      t268 = t32*bs
      t269 = t267+t268+Zoi-t86+t88-Zoj
      t270 = t269**2
      t272 = sqrt(t261+t266+t270)
      t273 = t272**2
      t282 = t30*bc
      t288 = t257+t259+Xoi-t53+t58-Xoj
      t289 = t288**2
      t290 = t262+t264+Yoi-t72+t76-Yoj
      t291 = t290**2
      t292 = t267+t268+Zoi-t86-t88-Zoj
      t293 = t292**2
      t295 = sqrt(t289+t291+t293)
      t296 = t295**2
      t310 = t257+t259+Xoi-t132-t134-Xoj
      t311 = t310**2
      t312 = t262+t264+Yoi-t137-t139-Yoj
      t313 = t312**2
      t314 = t267+t268+Zoi-t142-t143-Zoj
      t315 = t314**2
      t317 = sqrt(t311+t313+t315)
      t318 = t317**2
      t332 = t257+t259+Xoi-t132+t134-Xoj
      t333 = t332**2
      t334 = t262+t264+Yoi-t137+t139-Yoj
      t335 = t334**2
      t336 = t267+t268+Zoi-t142+t143-Zoj
      t337 = t336**2
      t339 = sqrt(t333+t335+t337)
      t340 = t339**2
      t354 = t257-t259+Xoi-t53-t58-Xoj
      t355 = t354**2
      t356 = t262-t264+Yoi-t72-t76-Yoj
      t357 = t356**2
      t358 = t267-t268+Zoi-t86+t88-Zoj
      t359 = t358**2
      t361 = sqrt(t355+t357+t359)
      t362 = t361**2
      t376 = t257-t259+Xoi-t53+t58-Xoj
      t377 = t376**2
      t378 = t262-t264+Yoi-t72+t76-Yoj
      t379 = t378**2
      t380 = t267-t268+Zoi-t86-t88-Zoj
      t381 = t380**2
      t383 = sqrt(t377+t379+t381)
      t384 = t383**2
      t398 = t257-t259+Xoi-t132-t134-Xoj
      t399 = t398**2
      t400 = t262-t264+Yoi-t137-t139-Yoj
      t401 = t400**2
      t402 = t267-t268+Zoi-t142-t143-Zoj
      t403 = t402**2
      t405 = sqrt(t399+t401+t403)
      t406 = t405**2
      t420 = t257-t259+Xoi-t132+t134-Xoj
      t421 = t420**2
      t422 = t262-t264+Yoi-t137+t139-Yoj
      t423 = t422**2
      t424 = t267-t268+Zoi-t142+t143-Zoj
      t425 = t424**2
      t427 = sqrt(t421+t423+t425)
      t428 = t427**2
      s1 = -2*qa1i*qa1j*ec/t93/t92*(t59*t99+t77*t103+t89*t107)-2*qa1i*qa
     #2j*ec/t121/t120*(t113*t99+t115*t103+t117*t107)-2*qa1i*qb1j*ec/t148
     #/t147*(t135*t99+t140*t103+t144*t107)-2*qa1i*qb2j*ec/t167/t166*(t15
     #9*t99+t161*t103+t163*t107)-2*qa2i*qa1j*ec/t186/t185*(t178*t190+t18
     #0*t192+t182*t194)-2*qa2i*qa2j*ec/t208/t207*(t200*t190+t202*t192+t2
     #04*t194)-2*qa2i*qb1j*ec/t227/t226*(t219*t190+t221*t192+t223*t194)-
     #2*qa2i*qb2j*ec/t246/t245*(t238*t190+t240*t192+t242*t194)
      t441 = s1-2*qb1i*qa1j*ec/t273/t272*(-t260*t42*bc-t265*t67*bc+t269*
     #t79*t282)-2*qb1i*qa2j*ec/t296/t295*(-t288*t42*bc-t290*t67*bc+t292*
     #t79*t282)-2*qb1i*qb1j*ec/t318/t317*(-t310*t42*bc-t312*t67*bc+t314*
     #t79*t282)-2*qb1i*qb2j*ec/t340/t339*(-t332*t42*bc-t334*t67*bc+t336*
     #t79*t282)-2*qb2i*qa1j*ec/t362/t361*(-t354*t42*bc-t356*t67*bc+t358*
     #t79*t282)-2*qb2i*qa2j*ec/t384/t383*(-t376*t42*bc-t378*t67*bc+t380*
     #t79*t282)-2*qb2i*qb1j*ec/t406/t405*(-t398*t42*bc-t400*t67*bc+t402*
     #t79*t282)-2*qb2i*qb2j*ec/t428/t427*(-t420*t42*bc-t422*t67*bc+t424*
     #t79*t282)
      t442 = (cst3*t16*(3*cst2-cst1-2*t14)/t22/t21+cst4)*t441/2
c
      dVdsi=t442
c
      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      function dVdppi(Xoi,Xoj,Yoi,Yoj,Zoi,Zoj,ppi,si,ti,ppj,sj,tj)
      implicit real*8(a-h,o-z)
      common/cpot1/qa1i,qa2i,qb1i,qb2i,qa1j,qa2j,qb1j,qb2j
      common/cpot2/AA,BB,ac,as,bc,bs,ec
      common/cpot3/cst1,cst2,cst3,cst4
c
      t1 = Xoi**2
      t4 = Xoj**2
      t5 = Yoi**2
      t8 = Yoj**2
      t9 = Zoi**2
      t12 = Zoj**2
      t14 = sqrt(t1-2*Xoi*Xoj+t4+t5-2*Yoi*Yoj+t8+t9-2*Zoi*Zoj+t12)
      t16 = (t14-cst1)**2
      t21 = cst2-cst1
      t22 = t21**2
      t29 = cos(ppi)
      t30 = cos(si)
      t32 = cos(ti)
      t33 = sin(si)
      t34 = t32*t33
      t35 = sin(ppi)
      t37 = t29*t30-t34*t35
      t38 = t37*ac
      t40 = t32*t30
      t42 = t29*t33+t40*t35
      t43 = t42*as
      t44 = cos(ppj)
      t45 = cos(sj)
      t47 = cos(tj)
      t48 = sin(sj)
      t49 = t47*t48
      t50 = sin(ppj)
      t52 = t44*t45-t49*t50
      t53 = t52*ac
      t55 = t47*t45
      t58 = (t44*t48+t55*t50)*as
      t59 = t38+t43+Xoi-t53-t58-Xoj
      t60 = t59**2
      t63 = -t35*t30-t34*t29
      t64 = t63*ac
      t68 = (-t35*t33+t40*t29)*as
      t71 = -t50*t45-t49*t44
      t72 = t71*ac
      t76 = (-t50*t48+t55*t44)*as
      t77 = t64+t68+Yoi-t72-t76-Yoj
      t78 = t77**2
      t79 = sin(ti)
      t80 = t79*t33
      t81 = t80*ac
      t83 = t79*t30*as
      t84 = sin(tj)
      t85 = t84*t48
      t86 = t85*ac
      t88 = t84*t45*as
      t90 = (t81-t83+Zoi-t86+t88-Zoj)**2
      t92 = sqrt(t60+t78+t90)
      t93 = t92**2
      t97 = t64+t68
      t99 = -t37*ac
      t100 = -t42*as
      t101 = t99+t100
      t107 = t38+t43+Xoi-t53+t58-Xoj
      t108 = t107**2
      t109 = t64+t68+Yoi-t72+t76-Yoj
      t110 = t109**2
      t112 = (t81-t83+Zoi-t86-t88-Zoj)**2
      t114 = sqrt(t108+t110+t112)
      t115 = t114**2
      t125 = t52*bc
      t127 = t50*t84*bs
      t128 = t38+t43+Xoi-t125-t127-Xoj
      t129 = t128**2
      t130 = t71*bc
      t132 = t44*t84*bs
      t133 = t64+t68+Yoi-t130-t132-Yoj
      t134 = t133**2
      t135 = t85*bc
      t136 = t47*bs
      t138 = (t81-t83+Zoi-t135-t136-Zoj)**2
      t140 = sqrt(t129+t134+t138)
      t141 = t140**2
      t151 = t38+t43+Xoi-t125+t127-Xoj
      t152 = t151**2
      t153 = t64+t68+Yoi-t130+t132-Yoj
      t154 = t153**2
      t156 = (t81-t83+Zoi-t135+t136-Zoj)**2
      t158 = sqrt(t152+t154+t156)
      t159 = t158**2
      t169 = t38-t43+Xoi-t53-t58-Xoj
      t170 = t169**2
      t171 = t64-t68+Yoi-t72-t76-Yoj
      t172 = t171**2
      t174 = (t81+t83+Zoi-t86+t88-Zoj)**2
      t176 = sqrt(t170+t172+t174)
      t177 = t176**2
      t181 = t64-t68
      t183 = t99-t100
      t189 = t38-t43+Xoi-t53+t58-Xoj
      t190 = t189**2
      t191 = t64-t68+Yoi-t72+t76-Yoj
      t192 = t191**2
      t194 = (t81+t83+Zoi-t86-t88-Zoj)**2
      t196 = sqrt(t190+t192+t194)
      t197 = t196**2
      t207 = t38-t43+Xoi-t125-t127-Xoj
      t208 = t207**2
      t209 = t64-t68+Yoi-t130-t132-Yoj
      t210 = t209**2
      t212 = (t81+t83+Zoi-t135-t136-Zoj)**2
      t214 = sqrt(t208+t210+t212)
      t215 = t214**2
      t225 = t38-t43+Xoi-t125+t127-Xoj
      t226 = t225**2
      t227 = t64-t68+Yoi-t130+t132-Yoj
      t228 = t227**2
      t230 = (t81+t83+Zoi-t135+t136-Zoj)**2
      t232 = sqrt(t226+t228+t230)
      t233 = t232**2
      t243 = t37*bc
      t245 = t35*t79*bs
      t246 = t243+t245+Xoi-t53-t58-Xoj
      t247 = t246**2
      t248 = t63*bc
      t250 = t29*t79*bs
      t251 = t248+t250+Yoi-t72-t76-Yoj
      t252 = t251**2
      t253 = t80*bc
      t254 = t32*bs
      t256 = (t253+t254+Zoi-t86+t88-Zoj)**2
      t258 = sqrt(t247+t252+t256)
      t259 = t258**2
      t263 = t248+t250
      t265 = -t37*bc
      t266 = t265-t245
      t272 = t243+t245+Xoi-t53+t58-Xoj
      t273 = t272**2
      t274 = t248+t250+Yoi-t72+t76-Yoj
      t275 = t274**2
      t277 = (t253+t254+Zoi-t86-t88-Zoj)**2
      t279 = sqrt(t273+t275+t277)
      t280 = t279**2
      t290 = t243+t245+Xoi-t125-t127-Xoj
      t291 = t290**2
      t292 = t248+t250+Yoi-t130-t132-Yoj
      t293 = t292**2
      t295 = (t253+t254+Zoi-t135-t136-Zoj)**2
      t297 = sqrt(t291+t293+t295)
      t298 = t297**2
      t308 = t243+t245+Xoi-t125+t127-Xoj
      t309 = t308**2
      t310 = t248+t250+Yoi-t130+t132-Yoj
      t311 = t310**2
      t313 = (t253+t254+Zoi-t135+t136-Zoj)**2
      t315 = sqrt(t309+t311+t313)
      t316 = t315**2
      t326 = t243-t245+Xoi-t53-t58-Xoj
      t327 = t326**2
      t328 = t248-t250+Yoi-t72-t76-Yoj
      t329 = t328**2
      t331 = (t253-t254+Zoi-t86+t88-Zoj)**2
      t333 = sqrt(t327+t329+t331)
      t334 = t333**2
      t338 = t248-t250
      t340 = t265+t245
      t346 = t243-t245+Xoi-t53+t58-Xoj
      t347 = t346**2
      t348 = t248-t250+Yoi-t72+t76-Yoj
      t349 = t348**2
      t351 = (t253-t254+Zoi-t86-t88-Zoj)**2
      t353 = sqrt(t347+t349+t351)
      t354 = t353**2
      t364 = t243-t245+Xoi-t125-t127-Xoj
      t365 = t364**2
      t366 = t248-t250+Yoi-t130-t132-Yoj
      t367 = t366**2
      t369 = (t253-t254+Zoi-t135-t136-Zoj)**2
      t371 = sqrt(t365+t367+t369)
      t372 = t371**2
      t382 = t243-t245+Xoi-t125+t127-Xoj
      t383 = t382**2
      t384 = t248-t250+Yoi-t130+t132-Yoj
      t385 = t384**2
      t387 = (t253-t254+Zoi-t135+t136-Zoj)**2
      t389 = sqrt(t383+t385+t387)
      t390 = t389**2
      s1 = -2*qa1i*qa1j*ec/t93/t92*(t59*t97+t77*t101)-2*qa1i*qa2j*ec/t11
     #5/t114*(t107*t97+t109*t101)-2*qa1i*qb1j*ec/t141/t140*(t128*t97+t13
     #3*t101)-2*qa1i*qb2j*ec/t159/t158*(t151*t97+t153*t101)-2*qa2i*qa1j*
     #ec/t177/t176*(t169*t181+t171*t183)-2*qa2i*qa2j*ec/t197/t196*(t189*
     #t181+t191*t183)-2*qa2i*qb1j*ec/t215/t214*(t207*t181+t209*t183)-2*q
     #a2i*qb2j*ec/t233/t232*(t225*t181+t227*t183)
      t399 = s1-2*qb1i*qa1j*ec/t259/t258*(t246*t263+t251*t266)-2*qb1i*qa
     #2j*ec/t280/t279*(t272*t263+t274*t266)-2*qb1i*qb1j*ec/t298/t297*(t2
     #90*t263+t292*t266)-2*qb1i*qb2j*ec/t316/t315*(t308*t263+t310*t266)-
     #2*qb2i*qa1j*ec/t334/t333*(t326*t338+t328*t340)-2*qb2i*qa2j*ec/t354
     #/t353*(t346*t338+t348*t340)-2*qb2i*qb1j*ec/t372/t371*(t364*t338+t3
     #66*t340)-2*qb2i*qb2j*ec/t390/t389*(t382*t338+t384*t340)
      t400 = (cst3*t16*(3*cst2-cst1-2*t14)/t22/t21+cst4)*t399/2
c
      dVdppi=t400
c
      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      function dVdXoj(Xoi,Xoj,Yoi,Yoj,Zoi,Zoj,ppi,si,ti,ppj,sj,tj)
      implicit real*8(a-h,o-z)
      common/cpot1/qa1i,qa2i,qb1i,qb2i,qa1j,qa2j,qb1j,qb2j
      common/cpot2/AA,BB,ac,as,bc,bs,ec  
      common/cpot3/cst1,cst2,cst3,cst4
c
      t1 = Xoi**2
      t4 = Xoj**2
      t5 = Yoi**2
      t8 = Yoj**2
      t9 = Zoi**2
      t12 = Zoj**2
      t13 = t1-2*Xoi*Xoj+t4+t5-2*Yoi*Yoj+t8+t9-2*Zoi*Zoj+t12
      t14 = sqrt(t13)
      t15 = t14-cst1
      t19 = 3*cst2-cst1-2*t14
      t21 = cst2-cst1
      t22 = t21**2
      t24 = 1/t22/t21
      t27 = -Xoi+Xoj
      t28 = 2*t24/t14*t27
      t30 = t15**2
      t31 = cst3*t30
      t34 = qa1i*qa1j
      t35 = cos(ppi)
      t36 = cos(si)
      t38 = cos(ti)
      t39 = sin(si)
      t40 = t38*t39
      t41 = sin(ppi)
      t43 = t35*t36-t40*t41
      t44 = t43*ac
      t46 = t38*t36
      t49 = (t35*t39+t46*t41)*as
      t50 = cos(ppj)
      t51 = cos(sj)
      t53 = cos(tj)
      t54 = sin(sj)
      t55 = t53*t54
      t56 = sin(ppj)
      t58 = t50*t51-t55*t56
      t59 = t58*ac
      t61 = t53*t51
      t64 = (t50*t54+t61*t56)*as
      t65 = t44+t49+Xoi-t59-t64-Xoj
      t66 = t65**2
      t69 = -t41*t36-t40*t35
      t70 = t69*ac
      t74 = (-t41*t39+t46*t35)*as
      t77 = -t56*t51-t55*t50
      t78 = t77*ac
      t82 = (-t56*t54+t61*t50)*as
      t84 = (t70+t74+Yoi-t78-t82-Yoj)**2
      t85 = sin(ti)
      t86 = t85*t39
      t87 = t86*ac
      t89 = t85*t36*as
      t90 = sin(tj)
      t91 = t90*t54
      t92 = t91*ac
      t94 = t90*t51*as
      t96 = (t87-t89+Zoi-t92+t94-Zoj)**2
      t98 = sqrt(t66+t84+t96)
      t102 = qa1i*qa2j
      t103 = t44+t49+Xoi-t59+t64-Xoj
      t104 = t103**2
      t106 = (t70+t74+Yoi-t78+t82-Yoj)**2
      t108 = (t87-t89+Zoi-t92-t94-Zoj)**2
      t110 = sqrt(t104+t106+t108)
      t114 = qa1i*qb1j
      t115 = t58*bc
      t117 = t56*t90*bs
      t118 = t44+t49+Xoi-t115-t117-Xoj
      t119 = t118**2
      t120 = t77*bc
      t122 = t50*t90*bs
      t124 = (t70+t74+Yoi-t120-t122-Yoj)**2
      t125 = t91*bc
      t126 = t53*bs
      t128 = (t87-t89+Zoi-t125-t126-Zoj)**2
      t130 = sqrt(t119+t124+t128)
      t134 = qa1i*qb2j
      t135 = t44+t49+Xoi-t115+t117-Xoj
      t136 = t135**2
      t138 = (t70+t74+Yoi-t120+t122-Yoj)**2
      t140 = (t87-t89+Zoi-t125+t126-Zoj)**2
      t142 = sqrt(t136+t138+t140)
      t146 = qa2i*qa1j
      t147 = t44-t49+Xoi-t59-t64-Xoj
      t148 = t147**2
      t150 = (t70-t74+Yoi-t78-t82-Yoj)**2
      t152 = (t87+t89+Zoi-t92+t94-Zoj)**2
      t154 = sqrt(t148+t150+t152)
      t158 = qa2i*qa2j
      t159 = t44-t49+Xoi-t59+t64-Xoj
      t160 = t159**2
      t162 = (t70-t74+Yoi-t78+t82-Yoj)**2
      t164 = (t87+t89+Zoi-t92-t94-Zoj)**2
      t166 = sqrt(t160+t162+t164)
      t170 = qa2i*qb1j
      t171 = t44-t49+Xoi-t115-t117-Xoj
      t172 = t171**2
      t174 = (t70-t74+Yoi-t120-t122-Yoj)**2
      t176 = (t87+t89+Zoi-t125-t126-Zoj)**2
      t178 = sqrt(t172+t174+t176)
      t182 = qa2i*qb2j
      t183 = t44-t49+Xoi-t115+t117-Xoj
      t184 = t183**2
      t186 = (t70-t74+Yoi-t120+t122-Yoj)**2
      t188 = (t87+t89+Zoi-t125+t126-Zoj)**2
      t190 = sqrt(t184+t186+t188)
      t194 = qb1i*qa1j
      t195 = t43*bc
      t197 = t41*t85*bs
      t198 = t195+t197+Xoi-t59-t64-Xoj
      t199 = t198**2
      t200 = t69*bc
      t202 = t35*t85*bs
      t204 = (t200+t202+Yoi-t78-t82-Yoj)**2
      t205 = t86*bc
      t206 = t38*bs
      t208 = (t205+t206+Zoi-t92+t94-Zoj)**2
      t210 = sqrt(t199+t204+t208)
      t214 = qb1i*qa2j
      t215 = t195+t197+Xoi-t59+t64-Xoj
      t216 = t215**2
      t218 = (t200+t202+Yoi-t78+t82-Yoj)**2
      t220 = (t205+t206+Zoi-t92-t94-Zoj)**2
      t222 = sqrt(t216+t218+t220)
      t226 = qb1i*qb1j
      t227 = t195+t197+Xoi-t115-t117-Xoj
      t228 = t227**2
      t230 = (t200+t202+Yoi-t120-t122-Yoj)**2
      t232 = (t205+t206+Zoi-t125-t126-Zoj)**2
      t234 = sqrt(t228+t230+t232)
      t238 = qb1i*qb2j
      t239 = t195+t197+Xoi-t115+t117-Xoj
      t240 = t239**2
      t242 = (t200+t202+Yoi-t120+t122-Yoj)**2
      t244 = (t205+t206+Zoi-t125+t126-Zoj)**2
      t246 = sqrt(t240+t242+t244)
      t250 = qb2i*qa1j
      t251 = t195-t197+Xoi-t59-t64-Xoj
      t252 = t251**2
      t254 = (t200-t202+Yoi-t78-t82-Yoj)**2
      t256 = (t205-t206+Zoi-t92+t94-Zoj)**2
      t258 = sqrt(t252+t254+t256)
      t262 = qb2i*qa2j
      t263 = t195-t197+Xoi-t59+t64-Xoj
      t264 = t263**2
      t266 = (t200-t202+Yoi-t78+t82-Yoj)**2
      t268 = (t205-t206+Zoi-t92-t94-Zoj)**2
      t270 = sqrt(t264+t266+t268)
      t274 = qb2i*qb1j
      t275 = t195-t197+Xoi-t115-t117-Xoj
      t276 = t275**2
      t278 = (t200-t202+Yoi-t120-t122-Yoj)**2
      t280 = (t205-t206+Zoi-t125-t126-Zoj)**2
      t282 = sqrt(t276+t278+t280)
      t286 = qb2i*qb2j
      t287 = t195-t197+Xoi-t115+t117-Xoj
      t288 = t287**2
      t290 = (t200-t202+Yoi-t120+t122-Yoj)**2
      t292 = (t205-t206+Zoi-t125+t126-Zoj)**2
      t294 = sqrt(t288+t290+t292)
      t298 = t34*ec/t98+t102*ec/t110+t114*ec/t130+t134*ec/t142+t146*ec/t
     #154+t158*ec/t166+t170*ec/t178+t182*ec/t190+t194*ec/t210+t214*ec/t2
     #22+t226*ec/t234+t238*ec/t246+t250*ec/t258+t262*ec/t270+t274*ec/t28
     #2+t286*ec/t294
      t303 = t98**2
      t309 = t110**2
      t315 = t130**2
      t321 = t142**2
      t327 = t154**2
      t333 = t166**2
      t339 = t178**2
      t345 = t190**2
      t351 = t210**2
      t357 = t222**2
      t363 = t234**2
      t369 = t246**2
      t375 = t258**2
      t381 = t270**2
      t387 = t282**2
      t393 = t294**2
      t399 = 2*t34*ec/t303/t98*t65+2*t102*ec/t309/t110*t103+2*t114*ec/t3
     #15/t130*t118+2*t134*ec/t321/t142*t135+2*t146*ec/t327/t154*t147+2*t
     #158*ec/t333/t166*t159+2*t170*ec/t339/t178*t171+2*t182*ec/t345/t190
     #*t183+2*t194*ec/t351/t210*t198+2*t214*ec/t357/t222*t215+2*t226*ec/
     #t363/t234*t227+2*t238*ec/t369/t246*t239+2*t250*ec/t375/t258*t251+2
     #*t262*ec/t381/t270*t263+2*t274*ec/t387/t282*t275+2*t286*ec/t393/t2
     #94*t287
      t401 = t13**2
      t403 = t401**2
      t413 = (cst3*t15*t19*t28-t31*t28)*t298+(t31*t19*t24+cst4)*t399/2-1
     #2*AA/t403/t401/t13*t27+6*BB/t403*t27
c
      dVdXoj=t413
c
      return
      end
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      function dVdYoj(Xoi,Xoj,Yoi,Yoj,Zoi,Zoj,ppi,si,ti,ppj,sj,tj)
      implicit real*8(a-h,o-z)
      common/cpot1/qa1i,qa2i,qb1i,qb2i,qa1j,qa2j,qb1j,qb2j
      common/cpot2/AA,BB,ac,as,bc,bs,ec
      common/cpot3/cst1,cst2,cst3,cst4
c
      t1 = Xoi**2
      t4 = Xoj**2
      t5 = Yoi**2
      t8 = Yoj**2
      t9 = Zoi**2
      t12 = Zoj**2
      t13 = t1-2*Xoi*Xoj+t4+t5-2*Yoi*Yoj+t8+t9-2*Zoi*Zoj+t12
      t14 = sqrt(t13)
      t15 = t14-cst1
      t19 = 3*cst2-cst1-2*t14
      t21 = cst2-cst1
      t22 = t21**2
      t24 = 1/t22/t21
      t27 = -Yoi+Yoj
      t28 = 2*t24/t14*t27
      t30 = t15**2
      t31 = cst3*t30
      t34 = qa1i*qa1j
      t35 = cos(ppi)
      t36 = cos(si)
      t38 = cos(ti)
      t39 = sin(si)
      t40 = t38*t39
      t41 = sin(ppi)
      t43 = t35*t36-t40*t41
      t44 = t43*ac
      t46 = t38*t36
      t49 = (t35*t39+t46*t41)*as
      t50 = cos(ppj)
      t51 = cos(sj)
      t53 = cos(tj)
      t54 = sin(sj)
      t55 = t53*t54
      t56 = sin(ppj)
      t58 = t50*t51-t55*t56
      t59 = t58*ac
      t61 = t53*t51
      t64 = (t50*t54+t61*t56)*as
      t66 = (t44+t49+Xoi-t59-t64-Xoj)**2
      t69 = -t41*t36-t40*t35
      t70 = t69*ac
      t74 = (-t41*t39+t46*t35)*as
      t77 = -t56*t51-t55*t50
      t78 = t77*ac
      t82 = (-t56*t54+t61*t50)*as
      t83 = t70+t74+Yoi-t78-t82-Yoj
      t84 = t83**2
      t85 = sin(ti)
      t86 = t85*t39
      t87 = t86*ac
      t89 = t85*t36*as
      t90 = sin(tj)
      t91 = t90*t54
      t92 = t91*ac
      t94 = t90*t51*as
      t96 = (t87-t89+Zoi-t92+t94-Zoj)**2
      t98 = sqrt(t66+t84+t96)
      t102 = qa1i*qa2j
      t104 = (t44+t49+Xoi-t59+t64-Xoj)**2
      t105 = t70+t74+Yoi-t78+t82-Yoj
      t106 = t105**2
      t108 = (t87-t89+Zoi-t92-t94-Zoj)**2
      t110 = sqrt(t104+t106+t108)
      t114 = qa1i*qb1j
      t115 = t58*bc
      t117 = t56*t90*bs
      t119 = (t44+t49+Xoi-t115-t117-Xoj)**2
      t120 = t77*bc
      t122 = t50*t90*bs
      t123 = t70+t74+Yoi-t120-t122-Yoj
      t124 = t123**2
      t125 = t91*bc
      t126 = t53*bs
      t128 = (t87-t89+Zoi-t125-t126-Zoj)**2
      t130 = sqrt(t119+t124+t128)
      t134 = qa1i*qb2j
      t136 = (t44+t49+Xoi-t115+t117-Xoj)**2
      t137 = t70+t74+Yoi-t120+t122-Yoj
      t138 = t137**2
      t140 = (t87-t89+Zoi-t125+t126-Zoj)**2
      t142 = sqrt(t136+t138+t140)
      t146 = qa2i*qa1j
      t148 = (t44-t49+Xoi-t59-t64-Xoj)**2
      t149 = t70-t74+Yoi-t78-t82-Yoj
      t150 = t149**2
      t152 = (t87+t89+Zoi-t92+t94-Zoj)**2
      t154 = sqrt(t148+t150+t152)
      t158 = qa2i*qa2j
      t160 = (t44-t49+Xoi-t59+t64-Xoj)**2
      t161 = t70-t74+Yoi-t78+t82-Yoj
      t162 = t161**2
      t164 = (t87+t89+Zoi-t92-t94-Zoj)**2
      t166 = sqrt(t160+t162+t164)
      t170 = qa2i*qb1j
      t172 = (t44-t49+Xoi-t115-t117-Xoj)**2
      t173 = t70-t74+Yoi-t120-t122-Yoj
      t174 = t173**2
      t176 = (t87+t89+Zoi-t125-t126-Zoj)**2
      t178 = sqrt(t172+t174+t176)
      t182 = qa2i*qb2j
      t184 = (t44-t49+Xoi-t115+t117-Xoj)**2
      t185 = t70-t74+Yoi-t120+t122-Yoj
      t186 = t185**2
      t188 = (t87+t89+Zoi-t125+t126-Zoj)**2
      t190 = sqrt(t184+t186+t188)
      t194 = qb1i*qa1j
      t195 = t43*bc
      t197 = t41*t85*bs
      t199 = (t195+t197+Xoi-t59-t64-Xoj)**2
      t200 = t69*bc
      t202 = t35*t85*bs
      t203 = t200+t202+Yoi-t78-t82-Yoj
      t204 = t203**2
      t205 = t86*bc
      t206 = t38*bs
      t208 = (t205+t206+Zoi-t92+t94-Zoj)**2
      t210 = sqrt(t199+t204+t208)
      t214 = qb1i*qa2j
      t216 = (t195+t197+Xoi-t59+t64-Xoj)**2
      t217 = t200+t202+Yoi-t78+t82-Yoj
      t218 = t217**2
      t220 = (t205+t206+Zoi-t92-t94-Zoj)**2
      t222 = sqrt(t216+t218+t220)
      t226 = qb1i*qb1j
      t228 = (t195+t197+Xoi-t115-t117-Xoj)**2
      t229 = t200+t202+Yoi-t120-t122-Yoj
      t230 = t229**2
      t232 = (t205+t206+Zoi-t125-t126-Zoj)**2
      t234 = sqrt(t228+t230+t232)
      t238 = qb1i*qb2j
      t240 = (t195+t197+Xoi-t115+t117-Xoj)**2
      t241 = t200+t202+Yoi-t120+t122-Yoj
      t242 = t241**2
      t244 = (t205+t206+Zoi-t125+t126-Zoj)**2
      t246 = sqrt(t240+t242+t244)
      t250 = qb2i*qa1j
      t252 = (t195-t197+Xoi-t59-t64-Xoj)**2
      t253 = t200-t202+Yoi-t78-t82-Yoj
      t254 = t253**2
      t256 = (t205-t206+Zoi-t92+t94-Zoj)**2
      t258 = sqrt(t252+t254+t256)
      t262 = qb2i*qa2j
      t264 = (t195-t197+Xoi-t59+t64-Xoj)**2
      t265 = t200-t202+Yoi-t78+t82-Yoj
      t266 = t265**2
      t268 = (t205-t206+Zoi-t92-t94-Zoj)**2
      t270 = sqrt(t264+t266+t268)
      t274 = qb2i*qb1j
      t276 = (t195-t197+Xoi-t115-t117-Xoj)**2
      t277 = t200-t202+Yoi-t120-t122-Yoj
      t278 = t277**2
      t280 = (t205-t206+Zoi-t125-t126-Zoj)**2
      t282 = sqrt(t276+t278+t280)
      t286 = qb2i*qb2j
      t288 = (t195-t197+Xoi-t115+t117-Xoj)**2
      t289 = t200-t202+Yoi-t120+t122-Yoj
      t290 = t289**2
      t292 = (t205-t206+Zoi-t125+t126-Zoj)**2
      t294 = sqrt(t288+t290+t292)
      t298 = t34*ec/t98+t102*ec/t110+t114*ec/t130+t134*ec/t142+t146*ec/t
     #154+t158*ec/t166+t170*ec/t178+t182*ec/t190+t194*ec/t210+t214*ec/t2
     #22+t226*ec/t234+t238*ec/t246+t250*ec/t258+t262*ec/t270+t274*ec/t28
     #2+t286*ec/t294
      t303 = t98**2
      t309 = t110**2
      t315 = t130**2
      t321 = t142**2
      t327 = t154**2
      t333 = t166**2
      t339 = t178**2
      t345 = t190**2
      t351 = t210**2
      t357 = t222**2
      t363 = t234**2
      t369 = t246**2
      t375 = t258**2
      t381 = t270**2
      t387 = t282**2
      t393 = t294**2
      t399 = 2*t34*ec/t303/t98*t83+2*t102*ec/t309/t110*t105+2*t114*ec/t3
     #15/t130*t123+2*t134*ec/t321/t142*t137+2*t146*ec/t327/t154*t149+2*t
     #158*ec/t333/t166*t161+2*t170*ec/t339/t178*t173+2*t182*ec/t345/t190
     #*t185+2*t194*ec/t351/t210*t203+2*t214*ec/t357/t222*t217+2*t226*ec/
     #t363/t234*t229+2*t238*ec/t369/t246*t241+2*t250*ec/t375/t258*t253+2
     #*t262*ec/t381/t270*t265+2*t274*ec/t387/t282*t277+2*t286*ec/t393/t2
     #94*t289
      t401 = t13**2
      t403 = t401**2
      t413 = (cst3*t15*t19*t28-t31*t28)*t298+(t31*t19*t24+cst4)*t399/2-1
     #2*AA/t403/t401/t13*t27+6*BB/t403*t27
c
      dVdYoj=t413
c
      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      function dVdZoj(Xoi,Xoj,Yoi,Yoj,Zoi,Zoj,ppi,si,ti,ppj,sj,tj)
      implicit real*8(a-h,o-z)
      common/cpot1/qa1i,qa2i,qb1i,qb2i,qa1j,qa2j,qb1j,qb2j
      common/cpot2/AA,BB,ac,as,bc,bs,ec  
      common/cpot3/cst1,cst2,cst3,cst4
c
      t1 = Xoi**2
      t4 = Xoj**2
      t5 = Yoi**2
      t8 = Yoj**2
      t9 = Zoi**2
      t12 = Zoj**2
      t13 = t1-2*Xoi*Xoj+t4+t5-2*Yoi*Yoj+t8+t9-2*Zoi*Zoj+t12
      t14 = sqrt(t13)
      t15 = t14-cst1
      t19 = 3*cst2-cst1-2*t14
      t21 = cst2-cst1
      t22 = t21**2
      t24 = 1/t22/t21
      t27 = -Zoi+Zoj
      t28 = 2*t24/t14*t27
      t30 = t15**2
      t31 = cst3*t30
      t34 = qa1i*qa1j
      t35 = cos(ppi)
      t36 = cos(si)
      t38 = cos(ti)
      t39 = sin(si)
      t40 = t38*t39
      t41 = sin(ppi)
      t43 = t35*t36-t40*t41
      t44 = t43*ac
      t46 = t38*t36
      t49 = (t35*t39+t46*t41)*as
      t50 = cos(ppj)
      t51 = cos(sj)
      t53 = cos(tj)
      t54 = sin(sj)
      t55 = t53*t54
      t56 = sin(ppj)
      t58 = t50*t51-t55*t56
      t59 = t58*ac
      t61 = t53*t51
      t64 = (t50*t54+t61*t56)*as
      t66 = (t44+t49+Xoi-t59-t64-Xoj)**2
      t69 = -t41*t36-t40*t35
      t70 = t69*ac
      t74 = (-t41*t39+t46*t35)*as
      t77 = -t56*t51-t55*t50
      t78 = t77*ac
      t82 = (-t56*t54+t61*t50)*as
      t84 = (t70+t74+Yoi-t78-t82-Yoj)**2
      t85 = sin(ti)
      t86 = t85*t39
      t87 = t86*ac
      t89 = t85*t36*as
      t90 = sin(tj)
      t91 = t90*t54
      t92 = t91*ac
      t94 = t90*t51*as
      t95 = t87-t89+Zoi-t92+t94-Zoj
      t96 = t95**2
      t98 = sqrt(t66+t84+t96)
      t102 = qa1i*qa2j
      t104 = (t44+t49+Xoi-t59+t64-Xoj)**2
      t106 = (t70+t74+Yoi-t78+t82-Yoj)**2
      t107 = t87-t89+Zoi-t92-t94-Zoj
      t108 = t107**2
      t110 = sqrt(t104+t106+t108)
      t114 = qa1i*qb1j
      t115 = t58*bc
      t117 = t56*t90*bs
      t119 = (t44+t49+Xoi-t115-t117-Xoj)**2
      t120 = t77*bc
      t122 = t50*t90*bs
      t124 = (t70+t74+Yoi-t120-t122-Yoj)**2
      t125 = t91*bc
      t126 = t53*bs
      t127 = t87-t89+Zoi-t125-t126-Zoj
      t128 = t127**2
      t130 = sqrt(t119+t124+t128)
      t134 = qa1i*qb2j
      t136 = (t44+t49+Xoi-t115+t117-Xoj)**2
      t138 = (t70+t74+Yoi-t120+t122-Yoj)**2
      t139 = t87-t89+Zoi-t125+t126-Zoj
      t140 = t139**2
      t142 = sqrt(t136+t138+t140)
      t146 = qa2i*qa1j
      t148 = (t44-t49+Xoi-t59-t64-Xoj)**2
      t150 = (t70-t74+Yoi-t78-t82-Yoj)**2
      t151 = t87+t89+Zoi-t92+t94-Zoj
      t152 = t151**2
      t154 = sqrt(t148+t150+t152)
      t158 = qa2i*qa2j
      t160 = (t44-t49+Xoi-t59+t64-Xoj)**2
      t162 = (t70-t74+Yoi-t78+t82-Yoj)**2
      t163 = t87+t89+Zoi-t92-t94-Zoj
      t164 = t163**2
      t166 = sqrt(t160+t162+t164)
      t170 = qa2i*qb1j
      t172 = (t44-t49+Xoi-t115-t117-Xoj)**2
      t174 = (t70-t74+Yoi-t120-t122-Yoj)**2
      t175 = t87+t89+Zoi-t125-t126-Zoj
      t176 = t175**2
      t178 = sqrt(t172+t174+t176)
      t182 = qa2i*qb2j
      t184 = (t44-t49+Xoi-t115+t117-Xoj)**2
      t186 = (t70-t74+Yoi-t120+t122-Yoj)**2
      t187 = t87+t89+Zoi-t125+t126-Zoj
      t188 = t187**2
      t190 = sqrt(t184+t186+t188)
      t194 = qb1i*qa1j
      t195 = t43*bc
      t197 = t41*t85*bs
      t199 = (t195+t197+Xoi-t59-t64-Xoj)**2
      t200 = t69*bc
      t202 = t35*t85*bs
      t204 = (t200+t202+Yoi-t78-t82-Yoj)**2
      t205 = t86*bc
      t206 = t38*bs
      t207 = t205+t206+Zoi-t92+t94-Zoj
      t208 = t207**2
      t210 = sqrt(t199+t204+t208)
      t214 = qb1i*qa2j
      t216 = (t195+t197+Xoi-t59+t64-Xoj)**2
      t218 = (t200+t202+Yoi-t78+t82-Yoj)**2
      t219 = t205+t206+Zoi-t92-t94-Zoj
      t220 = t219**2
      t222 = sqrt(t216+t218+t220)
      t226 = qb1i*qb1j
      t228 = (t195+t197+Xoi-t115-t117-Xoj)**2
      t230 = (t200+t202+Yoi-t120-t122-Yoj)**2
      t231 = t205+t206+Zoi-t125-t126-Zoj
      t232 = t231**2
      t234 = sqrt(t228+t230+t232)
      t238 = qb1i*qb2j
      t240 = (t195+t197+Xoi-t115+t117-Xoj)**2
      t242 = (t200+t202+Yoi-t120+t122-Yoj)**2
      t243 = t205+t206+Zoi-t125+t126-Zoj
      t244 = t243**2
      t246 = sqrt(t240+t242+t244)
      t250 = qb2i*qa1j
      t252 = (t195-t197+Xoi-t59-t64-Xoj)**2
      t254 = (t200-t202+Yoi-t78-t82-Yoj)**2
      t255 = t205-t206+Zoi-t92+t94-Zoj
      t256 = t255**2
      t258 = sqrt(t252+t254+t256)
      t262 = qb2i*qa2j
      t264 = (t195-t197+Xoi-t59+t64-Xoj)**2
      t266 = (t200-t202+Yoi-t78+t82-Yoj)**2
      t267 = t205-t206+Zoi-t92-t94-Zoj
      t268 = t267**2
      t270 = sqrt(t264+t266+t268)
      t274 = qb2i*qb1j
      t276 = (t195-t197+Xoi-t115-t117-Xoj)**2
      t278 = (t200-t202+Yoi-t120-t122-Yoj)**2
      t279 = t205-t206+Zoi-t125-t126-Zoj
      t280 = t279**2
      t282 = sqrt(t276+t278+t280)
      t286 = qb2i*qb2j
      t288 = (t195-t197+Xoi-t115+t117-Xoj)**2
      t290 = (t200-t202+Yoi-t120+t122-Yoj)**2
      t291 = t205-t206+Zoi-t125+t126-Zoj
      t292 = t291**2
      t294 = sqrt(t288+t290+t292)
      t298 = t34*ec/t98+t102*ec/t110+t114*ec/t130+t134*ec/t142+t146*ec/t
     #154+t158*ec/t166+t170*ec/t178+t182*ec/t190+t194*ec/t210+t214*ec/t2
     #22+t226*ec/t234+t238*ec/t246+t250*ec/t258+t262*ec/t270+t274*ec/t28
     #2+t286*ec/t294
      t303 = t98**2
      t309 = t110**2
      t315 = t130**2
      t321 = t142**2
      t327 = t154**2
      t333 = t166**2
      t339 = t178**2
      t345 = t190**2
      t351 = t210**2
      t357 = t222**2
      t363 = t234**2
      t369 = t246**2
      t375 = t258**2
      t381 = t270**2
      t387 = t282**2
      t393 = t294**2
      t399 = 2*t34*ec/t303/t98*t95+2*t102*ec/t309/t110*t107+2*t114*ec/t3
     #15/t130*t127+2*t134*ec/t321/t142*t139+2*t146*ec/t327/t154*t151+2*t
     #158*ec/t333/t166*t163+2*t170*ec/t339/t178*t175+2*t182*ec/t345/t190
     #*t187+2*t194*ec/t351/t210*t207+2*t214*ec/t357/t222*t219+2*t226*ec/
     #t363/t234*t231+2*t238*ec/t369/t246*t243+2*t250*ec/t375/t258*t255+2
     #*t262*ec/t381/t270*t267+2*t274*ec/t387/t282*t279+2*t286*ec/t393/t2
     #94*t291
      t401 = t13**2
      t403 = t401**2
      t413 = (cst3*t15*t19*t28-t31*t28)*t298+(t31*t19*t24+cst4)*t399/2-1
     #2*AA/t403/t401/t13*t27+6*BB/t403*t27
c
      dVdZoj=t413
c
      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      function dVdtj(Xoi,Xoj,Yoi,Yoj,Zoi,Zoj,ppi,si,ti,ppj,sj,tj)
      implicit real*8(a-h,o-z)
      common/cpot1/qa1i,qa2i,qb1i,qb2i,qa1j,qa2j,qb1j,qb2j
      common/cpot2/AA,BB,ac,as,bc,bs,ec 
      common/cpot3/cst1,cst2,cst3,cst4     
c
      t1 = Xoi**2
      t4 = Xoj**2
      t5 = Yoi**2
      t8 = Yoj**2
      t9 = Zoi**2
      t12 = Zoj**2
      t14 = sqrt(t1-2*Xoi*Xoj+t4+t5-2*Yoi*Yoj+t8+t9-2*Zoi*Zoj+t12)
      t16 = (t14-cst1)**2
      t21 = cst2-cst1
      t22 = t21**2
      t29 = cos(ppi)
      t30 = cos(si)
      t32 = cos(ti)
      t33 = sin(si)
      t34 = t32*t33
      t35 = sin(ppi)
      t37 = t29*t30-t34*t35
      t38 = t37*ac
      t40 = t32*t30
      t43 = (t29*t33+t40*t35)*as
      t44 = cos(ppj)
      t45 = cos(sj)
      t47 = cos(tj)
      t48 = sin(sj)
      t49 = t47*t48
      t50 = sin(ppj)
      t52 = t44*t45-t49*t50
      t53 = t52*ac
      t55 = t47*t45
      t58 = (t44*t48+t55*t50)*as
      t59 = t38+t43+Xoi-t53-t58-Xoj
      t60 = t59**2
      t63 = -t35*t30-t34*t29
      t64 = t63*ac
      t68 = (-t35*t33+t40*t29)*as
      t71 = -t50*t45-t49*t44
      t72 = t71*ac
      t76 = (-t50*t48+t55*t44)*as
      t77 = t64+t68+Yoi-t72-t76-Yoj
      t78 = t77**2
      t79 = sin(ti)
      t80 = t79*t33
      t81 = t80*ac
      t83 = t79*t30*as
      t84 = sin(tj)
      t85 = t84*t48
      t86 = t85*ac
      t87 = t84*t45
      t88 = t87*as
      t89 = t81-t83+Zoi-t86+t88-Zoj
      t90 = t89**2
      t92 = sqrt(t60+t78+t90)
      t93 = t92**2
      t98 = t85*t50*ac
      t100 = t87*t50*as
      t101 = -t98+t100
      t104 = t85*t44*ac
      t106 = t87*t44*as
      t107 = -t104+t106
      t109 = t49*ac
      t110 = t55*as
      t111 = -t109+t110
      t117 = t38+t43+Xoi-t53+t58-Xoj
      t118 = t117**2
      t119 = t64+t68+Yoi-t72+t76-Yoj
      t120 = t119**2
      t121 = t81-t83+Zoi-t86-t88-Zoj
      t122 = t121**2
      t124 = sqrt(t118+t120+t122)
      t125 = t124**2
      t129 = -t98-t100
      t131 = -t104-t106
      t133 = -t109-t110
      t139 = t52*bc
      t141 = t50*t84*bs
      t142 = t38+t43+Xoi-t139-t141-Xoj
      t143 = t142**2
      t144 = t71*bc
      t146 = t44*t84*bs
      t147 = t64+t68+Yoi-t144-t146-Yoj
      t148 = t147**2
      t149 = t85*bc
      t150 = t47*bs
      t151 = t81-t83+Zoi-t149-t150-Zoj
      t152 = t151**2
      t154 = sqrt(t143+t148+t152)
      t155 = t154**2
      t160 = t85*t50*bc
      t162 = t50*t47*bs
      t163 = -t160-t162
      t166 = t85*t44*bc
      t168 = t44*t47*bs
      t169 = -t166-t168
      t171 = t49*bc
      t172 = t84*bs
      t173 = -t171+t172
      t179 = t38+t43+Xoi-t139+t141-Xoj
      t180 = t179**2
      t181 = t64+t68+Yoi-t144+t146-Yoj
      t182 = t181**2
      t183 = t81-t83+Zoi-t149+t150-Zoj
      t184 = t183**2
      t186 = sqrt(t180+t182+t184)
      t187 = t186**2
      t191 = -t160+t162
      t193 = -t166+t168
      t195 = -t171-t172
      t201 = t38-t43+Xoi-t53-t58-Xoj
      t202 = t201**2
      t203 = t64-t68+Yoi-t72-t76-Yoj
      t204 = t203**2
      t205 = t81+t83+Zoi-t86+t88-Zoj
      t206 = t205**2
      t208 = sqrt(t202+t204+t206)
      t209 = t208**2
      t220 = t38-t43+Xoi-t53+t58-Xoj
      t221 = t220**2
      t222 = t64-t68+Yoi-t72+t76-Yoj
      t223 = t222**2
      t224 = t81+t83+Zoi-t86-t88-Zoj
      t225 = t224**2
      t227 = sqrt(t221+t223+t225)
      t228 = t227**2
      t239 = t38-t43+Xoi-t139-t141-Xoj
      t240 = t239**2
      t241 = t64-t68+Yoi-t144-t146-Yoj
      t242 = t241**2
      t243 = t81+t83+Zoi-t149-t150-Zoj
      t244 = t243**2
      t246 = sqrt(t240+t242+t244)
      t247 = t246**2
      t258 = t38-t43+Xoi-t139+t141-Xoj
      t259 = t258**2
      t260 = t64-t68+Yoi-t144+t146-Yoj
      t261 = t260**2
      t262 = t81+t83+Zoi-t149+t150-Zoj
      t263 = t262**2
      t265 = sqrt(t259+t261+t263)
      t266 = t265**2
      t277 = t37*bc
      t279 = t35*t79*bs
      t280 = t277+t279+Xoi-t53-t58-Xoj
      t281 = t280**2
      t282 = t63*bc
      t284 = t29*t79*bs
      t285 = t282+t284+Yoi-t72-t76-Yoj
      t286 = t285**2
      t287 = t80*bc
      t288 = t32*bs
      t289 = t287+t288+Zoi-t86+t88-Zoj
      t290 = t289**2
      t292 = sqrt(t281+t286+t290)
      t293 = t292**2
      t304 = t277+t279+Xoi-t53+t58-Xoj
      t305 = t304**2
      t306 = t282+t284+Yoi-t72+t76-Yoj
      t307 = t306**2
      t308 = t287+t288+Zoi-t86-t88-Zoj
      t309 = t308**2
      t311 = sqrt(t305+t307+t309)
      t312 = t311**2
      t323 = t277+t279+Xoi-t139-t141-Xoj
      t324 = t323**2
      t325 = t282+t284+Yoi-t144-t146-Yoj
      t326 = t325**2
      t327 = t287+t288+Zoi-t149-t150-Zoj
      t328 = t327**2
      t330 = sqrt(t324+t326+t328)
      t331 = t330**2
      t342 = t277+t279+Xoi-t139+t141-Xoj
      t343 = t342**2
      t344 = t282+t284+Yoi-t144+t146-Yoj
      t345 = t344**2
      t346 = t287+t288+Zoi-t149+t150-Zoj
      t347 = t346**2
      t349 = sqrt(t343+t345+t347)
      t350 = t349**2
      t361 = t277-t279+Xoi-t53-t58-Xoj
      t362 = t361**2
      t363 = t282-t284+Yoi-t72-t76-Yoj
      t364 = t363**2
      t365 = t287-t288+Zoi-t86+t88-Zoj
      t366 = t365**2
      t368 = sqrt(t362+t364+t366)
      t369 = t368**2
      t380 = t277-t279+Xoi-t53+t58-Xoj
      t381 = t380**2
      t382 = t282-t284+Yoi-t72+t76-Yoj
      t383 = t382**2
      t384 = t287-t288+Zoi-t86-t88-Zoj
      t385 = t384**2
      t387 = sqrt(t381+t383+t385)
      t388 = t387**2
      t399 = t277-t279+Xoi-t139-t141-Xoj
      t400 = t399**2
      t401 = t282-t284+Yoi-t144-t146-Yoj
      t402 = t401**2
      t403 = t287-t288+Zoi-t149-t150-Zoj
      t404 = t403**2
      t406 = sqrt(t400+t402+t404)
      t407 = t406**2
      t418 = t277-t279+Xoi-t139+t141-Xoj
      t419 = t418**2
      t420 = t282-t284+Yoi-t144+t146-Yoj
      t421 = t420**2
      t422 = t287-t288+Zoi-t149+t150-Zoj
      t423 = t422**2
      t425 = sqrt(t419+t421+t423)
      t426 = t425**2
      s1 = -2*qa1i*qa1j*ec/t93/t92*(t59*t101+t77*t107+t89*t111)-2*qa1i*q
     #a2j*ec/t125/t124*(t117*t129+t119*t131+t121*t133)-2*qa1i*qb1j*ec/t1
     #55/t154*(t142*t163+t147*t169+t151*t173)-2*qa1i*qb2j*ec/t187/t186*(
     #t179*t191+t181*t193+t183*t195)-2*qa2i*qa1j*ec/t209/t208*(t201*t101
     #+t203*t107+t205*t111)-2*qa2i*qa2j*ec/t228/t227*(t220*t129+t222*t13
     #1+t224*t133)-2*qa2i*qb1j*ec/t247/t246*(t239*t163+t241*t169+t243*t1
     #73)-2*qa2i*qb2j*ec/t266/t265*(t258*t191+t260*t193+t262*t195)
      t436 = s1-2*qb1i*qa1j*ec/t293/t292*(t280*t101+t285*t107+t289*t111)
     #-2*qb1i*qa2j*ec/t312/t311*(t304*t129+t306*t131+t308*t133)-2*qb1i*q
     #b1j*ec/t331/t330*(t323*t163+t325*t169+t327*t173)-2*qb1i*qb2j*ec/t3
     #50/t349*(t342*t191+t344*t193+t346*t195)-2*qb2i*qa1j*ec/t369/t368*(
     #t361*t101+t363*t107+t365*t111)-2*qb2i*qa2j*ec/t388/t387*(t380*t129
     #+t382*t131+t384*t133)-2*qb2i*qb1j*ec/t407/t406*(t399*t163+t401*t16
     #9+t403*t173)-2*qb2i*qb2j*ec/t426/t425*(t418*t191+t420*t193+t422*t1
     #95)
      t437 = (cst3*t16*(3*cst2-cst1-2*t14)/t22/t21+cst4)*t436/2
c
      dVdtj=t437
c
      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      function dVdsj(Xoi,Xoj,Yoi,Yoj,Zoi,Zoj,ppi,si,ti,ppj,sj,tj)
      implicit real*8(a-h,o-z)
      common/cpot1/qa1i,qa2i,qb1i,qb2i,qa1j,qa2j,qb1j,qb2j
      common/cpot2/AA,BB,ac,as,bc,bs,ec   
      common/cpot3/cst1,cst2,cst3,cst4
c     
      t1 = Xoi**2
      t4 = Xoj**2
      t5 = Yoi**2
      t8 = Yoj**2
      t9 = Zoi**2
      t12 = Zoj**2
      t14 = sqrt(t1-2*Xoi*Xoj+t4+t5-2*Yoi*Yoj+t8+t9-2*Zoi*Zoj+t12)
      t16 = (t14-cst1)**2
      t21 = cst2-cst1
      t22 = t21**2
      t29 = cos(ppi)
      t30 = cos(si)
      t32 = cos(ti)
      t33 = sin(si)
      t34 = t32*t33
      t35 = sin(ppi)
      t37 = t29*t30-t34*t35
      t38 = t37*ac
      t40 = t32*t30
      t43 = (t29*t33+t40*t35)*as
      t44 = cos(ppj)
      t45 = cos(sj)
      t47 = cos(tj)
      t48 = sin(sj)
      t49 = t47*t48
      t50 = sin(ppj)
      t52 = t44*t45-t49*t50
      t53 = t52*ac
      t55 = t47*t45
      t57 = t44*t48+t55*t50
      t58 = t57*as
      t59 = t38+t43+Xoi-t53-t58-Xoj
      t60 = t59**2
      t63 = -t35*t30-t34*t29
      t64 = t63*ac
      t68 = (-t35*t33+t40*t29)*as
      t71 = -t50*t45-t49*t44
      t72 = t71*ac
      t75 = -t50*t48+t55*t44
      t76 = t75*as
      t77 = t64+t68+Yoi-t72-t76-Yoj
      t78 = t77**2
      t79 = sin(ti)
      t80 = t79*t33
      t81 = t80*ac
      t83 = t79*t30*as
      t84 = sin(tj)
      t85 = t84*t48
      t86 = t85*ac
      t87 = t84*t45
      t88 = t87*as
      t89 = t81-t83+Zoi-t86+t88-Zoj
      t90 = t89**2
      t92 = sqrt(t60+t78+t90)
      t93 = t92**2
      t97 = -t57*ac
      t98 = t52*as
      t99 = -t97-t98
      t101 = -t75*ac
      t102 = t71*as
      t103 = -t101-t102
      t105 = t87*ac
      t106 = t85*as
      t107 = -t105-t106
      t113 = t38+t43+Xoi-t53+t58-Xoj
      t114 = t113**2
      t115 = t64+t68+Yoi-t72+t76-Yoj
      t116 = t115**2
      t117 = t81-t83+Zoi-t86-t88-Zoj
      t118 = t117**2
      t120 = sqrt(t114+t116+t118)
      t121 = t120**2
      t125 = -t97+t98
      t127 = -t101+t102
      t129 = -t105+t106
      t135 = t52*bc
      t137 = t50*t84*bs
      t138 = t38+t43+Xoi-t135-t137-Xoj
      t139 = t138**2
      t140 = t71*bc
      t142 = t44*t84*bs
      t143 = t64+t68+Yoi-t140-t142-Yoj
      t144 = t143**2
      t145 = t85*bc
      t146 = t47*bs
      t147 = t81-t83+Zoi-t145-t146-Zoj
      t148 = t147**2
      t150 = sqrt(t139+t144+t148)
      t151 = t150**2
      t160 = t45*bc
      t166 = t38+t43+Xoi-t135+t137-Xoj
      t167 = t166**2
      t168 = t64+t68+Yoi-t140+t142-Yoj
      t169 = t168**2
      t170 = t81-t83+Zoi-t145+t146-Zoj
      t171 = t170**2
      t173 = sqrt(t167+t169+t171)
      t174 = t173**2
      t188 = t38-t43+Xoi-t53-t58-Xoj
      t189 = t188**2
      t190 = t64-t68+Yoi-t72-t76-Yoj
      t191 = t190**2
      t192 = t81+t83+Zoi-t86+t88-Zoj
      t193 = t192**2
      t195 = sqrt(t189+t191+t193)
      t196 = t195**2
      t207 = t38-t43+Xoi-t53+t58-Xoj
      t208 = t207**2
      t209 = t64-t68+Yoi-t72+t76-Yoj
      t210 = t209**2
      t211 = t81+t83+Zoi-t86-t88-Zoj
      t212 = t211**2
      t214 = sqrt(t208+t210+t212)
      t215 = t214**2
      t226 = t38-t43+Xoi-t135-t137-Xoj
      t227 = t226**2
      t228 = t64-t68+Yoi-t140-t142-Yoj
      t229 = t228**2
      t230 = t81+t83+Zoi-t145-t146-Zoj
      t231 = t230**2
      t233 = sqrt(t227+t229+t231)
      t234 = t233**2
      t248 = t38-t43+Xoi-t135+t137-Xoj
      t249 = t248**2
      t250 = t64-t68+Yoi-t140+t142-Yoj
      t251 = t250**2
      t252 = t81+t83+Zoi-t145+t146-Zoj
      t253 = t252**2
      t255 = sqrt(t249+t251+t253)
      t256 = t255**2
      t270 = t37*bc
      t272 = t35*t79*bs
      t273 = t270+t272+Xoi-t53-t58-Xoj
      t274 = t273**2
      t275 = t63*bc
      t277 = t29*t79*bs
      t278 = t275+t277+Yoi-t72-t76-Yoj
      t279 = t278**2
      t280 = t80*bc
      t281 = t32*bs
      t282 = t280+t281+Zoi-t86+t88-Zoj
      t283 = t282**2
      t285 = sqrt(t274+t279+t283)
      t286 = t285**2
      t297 = t270+t272+Xoi-t53+t58-Xoj
      t298 = t297**2
      t299 = t275+t277+Yoi-t72+t76-Yoj
      t300 = t299**2
      t301 = t280+t281+Zoi-t86-t88-Zoj
      t302 = t301**2
      t304 = sqrt(t298+t300+t302)
      t305 = t304**2
      t316 = t270+t272+Xoi-t135-t137-Xoj
      t317 = t316**2
      t318 = t275+t277+Yoi-t140-t142-Yoj
      t319 = t318**2
      t320 = t280+t281+Zoi-t145-t146-Zoj
      t321 = t320**2
      t323 = sqrt(t317+t319+t321)
      t324 = t323**2
      t338 = t270+t272+Xoi-t135+t137-Xoj
      t339 = t338**2
      t340 = t275+t277+Yoi-t140+t142-Yoj
      t341 = t340**2
      t342 = t280+t281+Zoi-t145+t146-Zoj
      t343 = t342**2
      t345 = sqrt(t339+t341+t343)
      t346 = t345**2
      t360 = t270-t272+Xoi-t53-t58-Xoj
      t361 = t360**2
      t362 = t275-t277+Yoi-t72-t76-Yoj
      t363 = t362**2
      t364 = t280-t281+Zoi-t86+t88-Zoj
      t365 = t364**2
      t367 = sqrt(t361+t363+t365)
      t368 = t367**2
      t379 = t270-t272+Xoi-t53+t58-Xoj
      t380 = t379**2
      t381 = t275-t277+Yoi-t72+t76-Yoj
      t382 = t381**2
      t383 = t280-t281+Zoi-t86-t88-Zoj
      t384 = t383**2
      t386 = sqrt(t380+t382+t384)
      t387 = t386**2
      t398 = t270-t272+Xoi-t135-t137-Xoj
      t399 = t398**2
      t400 = t275-t277+Yoi-t140-t142-Yoj
      t401 = t400**2
      t402 = t280-t281+Zoi-t145-t146-Zoj
      t403 = t402**2
      t405 = sqrt(t399+t401+t403)
      t406 = t405**2
      t420 = t270-t272+Xoi-t135+t137-Xoj
      t421 = t420**2
      t422 = t275-t277+Yoi-t140+t142-Yoj
      t423 = t422**2
      t424 = t280-t281+Zoi-t145+t146-Zoj
      t425 = t424**2
      t427 = sqrt(t421+t423+t425)
      t428 = t427**2
      s1 = -2*qa1i*qa1j*ec/t93/t92*(t59*t99+t77*t103+t89*t107)-2*qa1i*qa
     #2j*ec/t121/t120*(t113*t125+t115*t127+t117*t129)-2*qa1i*qb1j*ec/t15
     #1/t150*(t138*t57*bc+t143*t75*bc-t147*t84*t160)-2*qa1i*qb2j*ec/t174
     #/t173*(t166*t57*bc+t168*t75*bc-t170*t84*t160)-2*qa2i*qa1j*ec/t196/
     #t195*(t188*t99+t190*t103+t192*t107)-2*qa2i*qa2j*ec/t215/t214*(t207
     #*t125+t209*t127+t211*t129)-2*qa2i*qb1j*ec/t234/t233*(t226*t57*bc+t
     #228*t75*bc-t230*t84*t160)-2*qa2i*qb2j*ec/t256/t255*(t248*t57*bc+t2
     #50*t75*bc-t252*t84*t160)
      t441 = s1-2*qb1i*qa1j*ec/t286/t285*(t273*t99+t278*t103+t282*t107)-
     #2*qb1i*qa2j*ec/t305/t304*(t297*t125+t299*t127+t301*t129)-2*qb1i*qb
     #1j*ec/t324/t323*(t316*t57*bc+t318*t75*bc-t320*t84*t160)-2*qb1i*qb2
     #j*ec/t346/t345*(t338*t57*bc+t340*t75*bc-t342*t84*t160)-2*qb2i*qa1j
     #*ec/t368/t367*(t360*t99+t362*t103+t364*t107)-2*qb2i*qa2j*ec/t387/t
     #386*(t379*t125+t381*t127+t383*t129)-2*qb2i*qb1j*ec/t406/t405*(t398
     #*t57*bc+t400*t75*bc-t402*t84*t160)-2*qb2i*qb2j*ec/t428/t427*(t420*
     #t57*bc+t422*t75*bc-t424*t84*t160)
      t442 = (cst3*t16*(3*cst2-cst1-2*t14)/t22/t21+cst4)*t441/2
c
      DVdsj=t442
c
      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      function dVdppj(Xoi,Xoj,Yoi,Yoj,Zoi,Zoj,ppi,si,ti,ppj,sj,tj)
      implicit real*8(a-h,o-z)
      common/cpot1/qa1i,qa2i,qb1i,qb2i,qa1j,qa2j,qb1j,qb2j
      common/cpot2/AA,BB,ac,as,bc,bs,ec
      common/cpot3/cst1,cst2,cst3,cst4
c
      t1 = Xoi**2
      t4 = Xoj**2
      t5 = Yoi**2
      t8 = Yoj**2
      t9 = Zoi**2
      t12 = Zoj**2
      t14 = sqrt(t1-2*Xoi*Xoj+t4+t5-2*Yoi*Yoj+t8+t9-2*Zoi*Zoj+t12)
      t16 = (t14-cst1)**2
      t21 = cst2-cst1
      t22 = t21**2
      t29 = cos(ppi)
      t30 = cos(si)
      t32 = cos(ti)
      t33 = sin(si)
      t34 = t32*t33
      t35 = sin(ppi)
      t37 = t29*t30-t34*t35
      t38 = t37*ac
      t40 = t32*t30
      t43 = (t29*t33+t40*t35)*as
      t44 = cos(ppj)
      t45 = cos(sj)
      t47 = cos(tj)
      t48 = sin(sj)
      t49 = t47*t48
      t50 = sin(ppj)
      t52 = t44*t45-t49*t50
      t53 = t52*ac
      t55 = t47*t45
      t57 = t44*t48+t55*t50
      t58 = t57*as
      t59 = t38+t43+Xoi-t53-t58-Xoj
      t60 = t59**2
      t63 = -t35*t30-t34*t29
      t64 = t63*ac
      t68 = (-t35*t33+t40*t29)*as
      t71 = -t50*t45-t49*t44
      t72 = t71*ac
      t76 = (-t50*t48+t55*t44)*as
      t77 = t64+t68+Yoi-t72-t76-Yoj
      t78 = t77**2
      t79 = sin(ti)
      t80 = t79*t33
      t81 = t80*ac
      t83 = t79*t30*as
      t84 = sin(tj)
      t85 = t84*t48
      t86 = t85*ac
      t88 = t84*t45*as
      t90 = (t81-t83+Zoi-t86+t88-Zoj)**2
      t92 = sqrt(t60+t78+t90)
      t93 = t92**2
      t97 = -t72-t76
      t99 = -t52*ac
      t100 = -t57*as
      t101 = -t99-t100
      t107 = t38+t43+Xoi-t53+t58-Xoj
      t108 = t107**2
      t109 = t64+t68+Yoi-t72+t76-Yoj
      t110 = t109**2
      t112 = (t81-t83+Zoi-t86-t88-Zoj)**2
      t114 = sqrt(t108+t110+t112)
      t115 = t114**2
      t119 = -t72+t76
      t121 = -t99+t100
      t127 = t52*bc
      t129 = t50*t84*bs
      t130 = t38+t43+Xoi-t127-t129-Xoj
      t131 = t130**2
      t132 = t71*bc
      t134 = t44*t84*bs
      t135 = t64+t68+Yoi-t132-t134-Yoj
      t136 = t135**2
      t137 = t85*bc
      t138 = t47*bs
      t140 = (t81-t83+Zoi-t137-t138-Zoj)**2
      t142 = sqrt(t131+t136+t140)
      t143 = t142**2
      t147 = -t132-t134
      t149 = -t52*bc
      t150 = -t149+t129
      t156 = t38+t43+Xoi-t127+t129-Xoj
      t157 = t156**2
      t158 = t64+t68+Yoi-t132+t134-Yoj
      t159 = t158**2
      t161 = (t81-t83+Zoi-t137+t138-Zoj)**2
      t163 = sqrt(t157+t159+t161)
      t164 = t163**2
      t168 = -t132+t134
      t170 = -t149-t129
      t176 = t38-t43+Xoi-t53-t58-Xoj
      t177 = t176**2
      t178 = t64-t68+Yoi-t72-t76-Yoj
      t179 = t178**2
      t181 = (t81+t83+Zoi-t86+t88-Zoj)**2
      t183 = sqrt(t177+t179+t181)
      t184 = t183**2
      t194 = t38-t43+Xoi-t53+t58-Xoj
      t195 = t194**2
      t196 = t64-t68+Yoi-t72+t76-Yoj
      t197 = t196**2
      t199 = (t81+t83+Zoi-t86-t88-Zoj)**2
      t201 = sqrt(t195+t197+t199)
      t202 = t201**2
      t212 = t38-t43+Xoi-t127-t129-Xoj
      t213 = t212**2
      t214 = t64-t68+Yoi-t132-t134-Yoj
      t215 = t214**2
      t217 = (t81+t83+Zoi-t137-t138-Zoj)**2
      t219 = sqrt(t213+t215+t217)
      t220 = t219**2
      t230 = t38-t43+Xoi-t127+t129-Xoj
      t231 = t230**2
      t232 = t64-t68+Yoi-t132+t134-Yoj
      t233 = t232**2
      t235 = (t81+t83+Zoi-t137+t138-Zoj)**2
      t237 = sqrt(t231+t233+t235)
      t238 = t237**2
      t248 = t37*bc
      t250 = t35*t79*bs
      t251 = t248+t250+Xoi-t53-t58-Xoj
      t252 = t251**2
      t253 = t63*bc
      t255 = t29*t79*bs
      t256 = t253+t255+Yoi-t72-t76-Yoj
      t257 = t256**2
      t258 = t80*bc
      t259 = t32*bs
      t261 = (t258+t259+Zoi-t86+t88-Zoj)**2
      t263 = sqrt(t252+t257+t261)
      t264 = t263**2
      t274 = t248+t250+Xoi-t53+t58-Xoj
      t275 = t274**2
      t276 = t253+t255+Yoi-t72+t76-Yoj
      t277 = t276**2
      t279 = (t258+t259+Zoi-t86-t88-Zoj)**2
      t281 = sqrt(t275+t277+t279)
      t282 = t281**2
      t292 = t248+t250+Xoi-t127-t129-Xoj
      t293 = t292**2
      t294 = t253+t255+Yoi-t132-t134-Yoj
      t295 = t294**2
      t297 = (t258+t259+Zoi-t137-t138-Zoj)**2
      t299 = sqrt(t293+t295+t297)
      t300 = t299**2
      t310 = t248+t250+Xoi-t127+t129-Xoj
      t311 = t310**2
      t312 = t253+t255+Yoi-t132+t134-Yoj
      t313 = t312**2
      t315 = (t258+t259+Zoi-t137+t138-Zoj)**2
      t317 = sqrt(t311+t313+t315)
      t318 = t317**2
      t328 = t248-t250+Xoi-t53-t58-Xoj
      t329 = t328**2
      t330 = t253-t255+Yoi-t72-t76-Yoj
      t331 = t330**2
      t333 = (t258-t259+Zoi-t86+t88-Zoj)**2
      t335 = sqrt(t329+t331+t333)
      t336 = t335**2
      t346 = t248-t250+Xoi-t53+t58-Xoj
      t347 = t346**2
      t348 = t253-t255+Yoi-t72+t76-Yoj
      t349 = t348**2
      t351 = (t258-t259+Zoi-t86-t88-Zoj)**2
      t353 = sqrt(t347+t349+t351)
      t354 = t353**2
      t364 = t248-t250+Xoi-t127-t129-Xoj
      t365 = t364**2
      t366 = t253-t255+Yoi-t132-t134-Yoj
      t367 = t366**2
      t369 = (t258-t259+Zoi-t137-t138-Zoj)**2
      t371 = sqrt(t365+t367+t369)
      t372 = t371**2
      t382 = t248-t250+Xoi-t127+t129-Xoj
      t383 = t382**2
      t384 = t253-t255+Yoi-t132+t134-Yoj
      t385 = t384**2
      t387 = (t258-t259+Zoi-t137+t138-Zoj)**2
      t389 = sqrt(t383+t385+t387)
      t390 = t389**2
      s1 = -2*qa1i*qa1j*ec/t93/t92*(t59*t97+t77*t101)-2*qa1i*qa2j*ec/t11
     #5/t114*(t107*t119+t109*t121)-2*qa1i*qb1j*ec/t143/t142*(t130*t147+t
     #135*t150)-2*qa1i*qb2j*ec/t164/t163*(t156*t168+t158*t170)-2*qa2i*qa
     #1j*ec/t184/t183*(t176*t97+t178*t101)-2*qa2i*qa2j*ec/t202/t201*(t19
     #4*t119+t196*t121)-2*qa2i*qb1j*ec/t220/t219*(t212*t147+t214*t150)-2
     #*qa2i*qb2j*ec/t238/t237*(t230*t168+t232*t170)
      t399 = s1-2*qb1i*qa1j*ec/t264/t263*(t251*t97+t256*t101)-2*qb1i*qa2
     #j*ec/t282/t281*(t274*t119+t276*t121)-2*qb1i*qb1j*ec/t300/t299*(t29
     #2*t147+t294*t150)-2*qb1i*qb2j*ec/t318/t317*(t310*t168+t312*t170)-2
     #*qb2i*qa1j*ec/t336/t335*(t328*t97+t330*t101)-2*qb2i*qa2j*ec/t354/t
     #353*(t346*t119+t348*t121)-2*qb2i*qb1j*ec/t372/t371*(t364*t147+t366
     #*t150)-2*qb2i*qb2j*ec/t390/t389*(t382*t168+t384*t170)
      t400 = (cst3*t16*(3*cst2-cst1-2*t14)/t22/t21+cst4)*t399/2
c
      dVdppj=t400
c
      return
      end
