       	subroutine diffvtip3p(pot_chos,natoms,nmono,xyz,v2tot,v3tot
     &                     ,v,diffv)
	implicit real*8(a-h,o-z)
	include "param.inc"
	include "common.inc"
        common/cpot/qo,qh,AA,BB,ac,as,ec
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
	ec=ecag 
        AA=At
        BB=Ct
        pi=acos(-1d0)
        ap=alph1*pi/360d0
        as=dioa1*sin(ap)
        ac=dioa2*cos(ap)
        qo=qia1
        qh=qia2
c	write(*,*)'ecag,At,Ct,ap,as,ac,qo,qh',ecag,At,Ct,ap,as,ac,qo,qh
c	write(*,*)'ec,AA,BB,alph1,dioa1,dioa2,qia1,qia2',ec,AA,BB,alph1,dioa1,dioa2,qia1,qia2
c	read(*,*)
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
c
	      v=v+vpot3p(Xoi,Xoj,Yoi,Yoj,Zoi,Zoj,ppi,si,ti,ppj,sj,tj)
c
              dvdx=dVdXoi3p(Xoi,Xoj,Yoi,Yoj,Zoi,Zoj
     &                   ,ppi,si,ti,ppj,sj,tj)
              dvdy=dVdYoi3p(Xoi,Xoj,Yoi,Yoj,Zoi,Zoj
     &                   ,ppi,si,ti,ppj,sj,tj) 
              dvdz=dVdZoi3p(Xoi,Xoj,Yoi,Yoj,Zoi,Zoj
     &                   ,ppi,si,ti,ppj,sj,tj)
c
	      dv2dxo(i)=dv2dxo(i)+dvdx
	      dv2dyo(i)=dv2dyo(i)+dvdy
	      dv2dzo(i)=dv2dzo(i)+dvdz
	      dv2dt(i)=dv2dt(i)+dVdti3p(Xoi,Xoj,Yoi,Yoj,Zoi,Zoj
     &             ,ppi,si,ti,ppj,sj,tj)
	      dv2ds(i)=dv2ds(i)+dVdsi3p(Xoi,Xoj,Yoi,Yoj,Zoi,Zoj
     &             ,ppi,si,ti,ppj,sj,tj)
	      dv2dpp(i)=dv2dpp(i)+dVdppi3p(Xoi,Xoj,Yoi,Yoj,Zoi,Zoj
     &             ,ppi,si,ti,ppj,sj,tj)
              dv2dxo(j)=dv2dxo(j)-dvdx
	      dv2dyo(j)=dv2dyo(j)-dvdy
	      dv2dzo(j)=dv2dzo(j)-dvdz
	      dv2dt(j)=dv2dt(j)+dVdtj3p(Xoi,Xoj,Yoi,Yoj,Zoi,Zoj
     &             ,ppi,si,ti,ppj,sj,tj)
	      dv2ds(j)=dv2ds(j)+dVdsj3p(Xoi,Xoj,Yoi,Yoj,Zoi,Zoj
     &             ,ppi,si,ti,ppj,sj,tj)
	      dv2dpp(j)=dv2dpp(j)+dVdppj3p(Xoi,Xoj,Yoi,Yoj,Zoi,Zoj
     &             ,ppi,si,ti,ppj,sj,tj)
c
c
c
c              d1=dVdXoi3p(Xoi,Xoj,Yoi,Yoj,Zoi,Zoj
c     &             ,ppi,si,ti,ppj,sj,tj,as,ac)
c              d2=dVdXoj3p(Xoi,Xoj,Yoi,Yoj,Zoi,Zoj
c     &             ,ppi,si,ti,ppj,sj,tj,as,ac)
c              d3=dVdYoi3p(Xoi,Xoj,Yoi,Yoj,Zoi,Zoj
c     &             ,ppi,si,ti,ppj,sj,tj,as,ac)
c              d4=dVdYoj3p(Xoi,Xoj,Yoi,Yoj,Zoi,Zoj
c     &             ,ppi,si,ti,ppj,sj,tj,as,ac)
c              d5=dVdZoi3p(Xoi,Xoj,Yoi,Yoj,Zoi,Zoj
c     &             ,ppi,si,ti,ppj,sj,tj,as,ac)
c              d6=dVdZoj3p(Xoi,Xoj,Yoi,Yoj,Zoi,Zoj
c     &             ,ppi,si,ti,ppj,sj,tj,as,ac)
c              d7=dVdti3p(Xoi,Xoj,Yoi,Yoj,Zoi,Zoj
c     &             ,ppi,si,ti,ppj,sj,tj,as,ac)
c              d8=dVdtj3p(Xoi,Xoj,Yoi,Yoj,Zoi,Zoj
c     &             ,ppi,si,ti,ppj,sj,tj,as,ac)
c              d9=dVdsi3p(Xoi,Xoj,Yoi,Yoj,Zoi,Zoj
c     &             ,ppi,si,ti,ppj,sj,tj,as,ac)
c              d10=dVdsj3p(Xoi,Xoj,Yoi,Yoj,Zoi,Zoj
c     &             ,ppi,si,ti,ppj,sj,tj,as,ac)
c              d11=dVdppi3p(Xoi,Xoj,Yoi,Yoj,Zoi,Zoj
c     &             ,ppi,si,ti,ppj,sj,tj,as,ac)
c              d12=dVdppj3p(Xoi,Xoj,Yoi,Yoj,Zoi,Zoj
c     &             ,ppi,si,ti,ppj,sj,tj,as,ac)
c        
c              
c              write(*,*)'i,j,d1,d2',i,j,d1,d2
c              write(*,*)'d3,d4',d3,d4
c              write(*,*)'d5,d6',d5,d6
c              write(*,*)'d7,d8',d7,d8
c              write(*,*)'d9,d10',d9,d10
c              write(*,*)'d11,d12',d11,d12
c             write(*,*)'v',v
c             read(*,*)
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
      function vpot3p(Xoi,Xoj,Yoi,Yoj,Zoi,Zoj,ppi,si,ti,ppj,sj,tj)
      implicit real*8(a-h,o-z)
      common/cpot/qo,qh,AA,BB,ac,as,ec  
c      
      t1 = qo**2
      t3 = Xoi**2
      t6 = Xoj**2
      t7 = Yoi**2
      t10 = Yoj**2
      t11 = Zoi**2
      t14 = Zoj**2
      t16 = sqrt(t3-2*Xoi*Xoj+t6+t7-2*Yoi*Yoj+t10+t11-2*Zoi*Zoj+t14)
      t17 = t16+0.1E-98
      t20 = qo*qh
      t21 = cos(ppj)
      t22 = cos(sj)
      t24 = cos(tj)
      t25 = sin(sj)
      t26 = t24*t25
      t27 = sin(ppj)
      t30 = (t21*t22-t26*t27)*ac
      t32 = t24*t22
      t35 = (t21*t25+t32*t27)*as
      t37 = (Xoi-t30-t35-Xoj)**2
      t41 = (-t27*t22-t26*t21)*ac
      t45 = (-t27*t25+t32*t21)*as
      t47 = (Yoi-t41-t45-Yoj)**2
      t48 = sin(tj)
      t50 = t48*t25*ac
      t52 = t48*t22*as
      t54 = (Zoi-t50+t52-Zoj)**2
      t56 = sqrt(t37+t47+t54)
      t62 = (Xoi-t30+t35-Xoj)**2
      t64 = (Yoi-t41+t45-Yoj)**2
      t66 = (Zoi-t50-t52-Zoj)**2
      t68 = sqrt(t62+t64+t66)
      t73 = cos(ppi)
      t74 = cos(si)
      t76 = cos(ti)
      t77 = sin(si)
      t78 = t76*t77
      t79 = sin(ppi)
      t82 = (t73*t74-t78*t79)*ac
      t84 = t76*t74
      t87 = (t73*t77+t84*t79)*as
      t89 = (t82+t87+Xoi-Xoj)**2
      t93 = (-t79*t74-t78*t73)*ac
      t97 = (-t79*t77+t84*t73)*as
      t99 = (t93+t97+Yoi-Yoj)**2
      t100 = sin(ti)
      t102 = t100*t77*ac
      t104 = t100*t74*as
      t106 = (t102-t104+Zoi-Zoj)**2
      t108 = sqrt(t89+t99+t106)
      t113 = qh**2
      t114 = t113*ec
      t116 = (t82+t87+Xoi-t30-t35-Xoj)**2
      t118 = (t93+t97+Yoi-t41-t45-Yoj)**2
      t120 = (t102-t104+Zoi-t50+t52-Zoj)**2
      t122 = sqrt(t116+t118+t120)
      t127 = (t82+t87+Xoi-t30+t35-Xoj)**2
      t129 = (t93+t97+Yoi-t41+t45-Yoj)**2
      t131 = (t102-t104+Zoi-t50-t52-Zoj)**2
      t133 = sqrt(t127+t129+t131)
      t138 = (t82-t87+Xoi-Xoj)**2
      t140 = (t93-t97+Yoi-Yoj)**2
      t142 = (t102+t104+Zoi-Zoj)**2
      t144 = sqrt(t138+t140+t142)
      t150 = (t82-t87+Xoi-t30-t35-Xoj)**2
      t152 = (t93-t97+Yoi-t41-t45-Yoj)**2
      t154 = (t102+t104+Zoi-t50+t52-Zoj)**2
      t156 = sqrt(t150+t152+t154)
      t161 = (t82-t87+Xoi-t30+t35-Xoj)**2
      t163 = (t93-t97+Yoi-t41+t45-Yoj)**2
      t165 = (t102+t104+Zoi-t50-t52-Zoj)**2
      t167 = sqrt(t161+t163+t165)
      t171 = t17**2
      t172 = t171**2
      t173 = t172**2
      t180 = t1*ec/t17+t20*ec/(t56+0.1E-98)+t20*ec/(t68+0.1E-98)+t20*ec/
     #(t108+0.1E-98)+t114/(t122+0.1E-98)+t114/(t133+0.1E-98)+t20*ec/(t14
     #4+0.1E-98)+t114/(t156+0.1E-98)+t114/(t167+0.1E-98)+AA/t173/t172-BB
     #/t172/t171
c
      vpot3p=t180
c
      return
      end
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      function dVdXoi3p(Xoi,Xoj,Yoi,Yoj,Zoi,Zoj,ppi,si,ti,ppj,sj,tj)
      implicit real*8(a-h,o-z)
      common/cpot/qo,qh,AA,BB,ac,as,ec    
c
      t1 = qo**2
      t3 = Xoi**2
      t6 = Xoj**2
      t7 = Yoi**2
      t10 = Yoj**2
      t11 = Zoi**2
      t14 = Zoj**2
      t16 = sqrt(t3-2*Xoi*Xoj+t6+t7-2*Yoi*Yoj+t10+t11-2*Zoi*Zoj+t14)
      t17 = t16+0.1E-98
      t18 = t17**2
      t20 = 1/t16
      t22 = Xoi-Xoj
      t27 = qo*qh*ec
      t28 = cos(ppj)
      t29 = cos(sj)
      t31 = cos(tj)
      t32 = sin(sj)
      t33 = t31*t32
      t34 = sin(ppj)
      t37 = (t28*t29-t33*t34)*ac
      t39 = t31*t29
      t42 = (t28*t32+t39*t34)*as
      t43 = Xoi-t37-t42-Xoj
      t44 = t43**2
      t48 = (-t34*t29-t33*t28)*ac
      t52 = (-t34*t32+t39*t28)*as
      t54 = (Yoi-t48-t52-Yoj)**2
      t55 = sin(tj)
      t57 = t55*t32*ac
      t59 = t55*t29*as
      t61 = (Zoi-t57+t59-Zoj)**2
      t63 = sqrt(t44+t54+t61)
      t65 = (t63+0.1E-98)**2
      t72 = Xoi-t37+t42-Xoj
      t73 = t72**2
      t75 = (Yoi-t48+t52-Yoj)**2
      t77 = (Zoi-t57-t59-Zoj)**2
      t79 = sqrt(t73+t75+t77)
      t81 = (t79+0.1E-98)**2
      t88 = cos(ppi)
      t89 = cos(si)
      t91 = cos(ti)
      t92 = sin(si)
      t93 = t91*t92
      t94 = sin(ppi)
      t97 = (t88*t89-t93*t94)*ac
      t99 = t91*t89
      t102 = (t88*t92+t99*t94)*as
      t103 = t97+t102+Xoi-Xoj
      t104 = t103**2
      t108 = (-t94*t89-t93*t88)*ac
      t112 = (-t94*t92+t99*t88)*as
      t114 = (t108+t112+Yoi-Yoj)**2
      t115 = sin(ti)
      t117 = t115*t92*ac
      t119 = t115*t89*as
      t121 = (t117-t119+Zoi-Zoj)**2
      t123 = sqrt(t104+t114+t121)
      t125 = (t123+0.1E-98)**2
      t132 = qh**2
      t133 = t132*ec
      t134 = t97+t102+Xoi-t37-t42-Xoj
      t135 = t134**2
      t137 = (t108+t112+Yoi-t48-t52-Yoj)**2
      t139 = (t117-t119+Zoi-t57+t59-Zoj)**2
      t141 = sqrt(t135+t137+t139)
      t143 = (t141+0.1E-98)**2
      t150 = t97+t102+Xoi-t37+t42-Xoj
      t151 = t150**2
      t153 = (t108+t112+Yoi-t48+t52-Yoj)**2
      t155 = (t117-t119+Zoi-t57-t59-Zoj)**2
      t157 = sqrt(t151+t153+t155)
      t159 = (t157+0.1E-98)**2
      t166 = t97-t102+Xoi-Xoj
      t167 = t166**2
      t169 = (t108-t112+Yoi-Yoj)**2
      t171 = (t117+t119+Zoi-Zoj)**2
      t173 = sqrt(t167+t169+t171)
      t175 = (t173+0.1E-98)**2
      t182 = t97-t102+Xoi-t37-t42-Xoj
      t183 = t182**2
      t185 = (t108-t112+Yoi-t48-t52-Yoj)**2
      t187 = (t117+t119+Zoi-t57+t59-Zoj)**2
      t189 = sqrt(t183+t185+t187)
      t191 = (t189+0.1E-98)**2
      t198 = t97-t102+Xoi-t37+t42-Xoj
      t199 = t198**2
      t201 = (t108-t112+Yoi-t48+t52-Yoj)**2
      t203 = (t117+t119+Zoi-t57-t59-Zoj)**2
      t205 = sqrt(t199+t201+t203)
      t207 = (t205+0.1E-98)**2
      t214 = t18**2
      t216 = t214**2
      t220 = 2*t20*t22
      t229 = -t1*ec/t18*t20*t22-t27/t65/t63*t43-t27/t81/t79*t72-t27/t125
     #/t123*t103-t133/t143/t141*t134-t133/t159/t157*t150-t27/t175/t173*t
     #166-t133/t191/t189*t182-t133/t207/t205*t198-6*AA/t216/t214/t17*t22
     #0+3*BB/t214/t18/t17*t220
c
      dVdXoi3p=t229 
c
      return
      end
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      function dVdYoi3p(Xoi,Xoj,Yoi,Yoj,Zoi,Zoj,ppi,si,ti,ppj,sj,tj)
      implicit real*8(a-h,o-z)
      common/cpot/qo,qh,AA,BB,ac,as,ec     
c
      t1 = qo**2
      t3 = Xoi**2
      t6 = Xoj**2
      t7 = Yoi**2
      t10 = Yoj**2
      t11 = Zoi**2
      t14 = Zoj**2
      t16 = sqrt(t3-2*Xoi*Xoj+t6+t7-2*Yoi*Yoj+t10+t11-2*Zoi*Zoj+t14)
      t17 = t16+0.1E-98
      t18 = t17**2
      t20 = 1/t16
      t22 = Yoi-Yoj
      t27 = qo*qh*ec
      t28 = cos(ppj)
      t29 = cos(sj)
      t31 = cos(tj)
      t32 = sin(sj)
      t33 = t31*t32
      t34 = sin(ppj)
      t37 = (t28*t29-t33*t34)*ac
      t39 = t31*t29
      t42 = (t28*t32+t39*t34)*as
      t44 = (Xoi-t37-t42-Xoj)**2
      t48 = (-t34*t29-t33*t28)*ac
      t52 = (-t34*t32+t39*t28)*as
      t53 = Yoi-t48-t52-Yoj
      t54 = t53**2
      t55 = sin(tj)
      t57 = t55*t32*ac
      t59 = t55*t29*as
      t61 = (Zoi-t57+t59-Zoj)**2
      t63 = sqrt(t44+t54+t61)
      t65 = (t63+0.1E-98)**2
      t73 = (Xoi-t37+t42-Xoj)**2
      t74 = Yoi-t48+t52-Yoj
      t75 = t74**2
      t77 = (Zoi-t57-t59-Zoj)**2
      t79 = sqrt(t73+t75+t77)
      t81 = (t79+0.1E-98)**2
      t88 = cos(ppi)
      t89 = cos(si)
      t91 = cos(ti)
      t92 = sin(si)
      t93 = t91*t92
      t94 = sin(ppi)
      t97 = (t88*t89-t93*t94)*ac
      t99 = t91*t89
      t102 = (t88*t92+t99*t94)*as
      t104 = (t97+t102+Xoi-Xoj)**2
      t108 = (-t94*t89-t93*t88)*ac
      t112 = (-t94*t92+t99*t88)*as
      t113 = t108+t112+Yoi-Yoj
      t114 = t113**2
      t115 = sin(ti)
      t117 = t115*t92*ac
      t119 = t115*t89*as
      t121 = (t117-t119+Zoi-Zoj)**2
      t123 = sqrt(t104+t114+t121)
      t125 = (t123+0.1E-98)**2
      t132 = qh**2
      t133 = t132*ec
      t135 = (t97+t102+Xoi-t37-t42-Xoj)**2
      t136 = t108+t112+Yoi-t48-t52-Yoj
      t137 = t136**2
      t139 = (t117-t119+Zoi-t57+t59-Zoj)**2
      t141 = sqrt(t135+t137+t139)
      t143 = (t141+0.1E-98)**2
      t151 = (t97+t102+Xoi-t37+t42-Xoj)**2
      t152 = t108+t112+Yoi-t48+t52-Yoj
      t153 = t152**2
      t155 = (t117-t119+Zoi-t57-t59-Zoj)**2
      t157 = sqrt(t151+t153+t155)
      t159 = (t157+0.1E-98)**2
      t167 = (t97-t102+Xoi-Xoj)**2
      t168 = t108-t112+Yoi-Yoj
      t169 = t168**2
      t171 = (t117+t119+Zoi-Zoj)**2
      t173 = sqrt(t167+t169+t171)
      t175 = (t173+0.1E-98)**2
      t183 = (t97-t102+Xoi-t37-t42-Xoj)**2
      t184 = t108-t112+Yoi-t48-t52-Yoj
      t185 = t184**2
      t187 = (t117+t119+Zoi-t57+t59-Zoj)**2
      t189 = sqrt(t183+t185+t187)
      t191 = (t189+0.1E-98)**2
      t199 = (t97-t102+Xoi-t37+t42-Xoj)**2
      t200 = t108-t112+Yoi-t48+t52-Yoj
      t201 = t200**2
      t203 = (t117+t119+Zoi-t57-t59-Zoj)**2
      t205 = sqrt(t199+t201+t203)
      t207 = (t205+0.1E-98)**2
      t214 = t18**2
      t216 = t214**2
      t220 = 2*t20*t22
      t229 = -t1*ec/t18*t20*t22-t27/t65/t63*t53-t27/t81/t79*t74-t27/t125
     #/t123*t113-t133/t143/t141*t136-t133/t159/t157*t152-t27/t175/t173*t
     #168-t133/t191/t189*t184-t133/t207/t205*t200-6*AA/t216/t214/t17*t22
     #0+3*BB/t214/t18/t17*t220
c
      dVdYoi3p=t229
c
      return
      end
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      function dVdZoi3p(Xoi,Xoj,Yoi,Yoj,Zoi,Zoj,ppi,si,ti,ppj,sj,tj)
      implicit real*8(a-h,o-z)
      common/cpot/qo,qh,AA,BB,ac,as,ec
c
      t1 = qo**2
      t3 = Xoi**2
      t6 = Xoj**2
      t7 = Yoi**2
      t10 = Yoj**2
      t11 = Zoi**2
      t14 = Zoj**2
      t16 = sqrt(t3-2*Xoi*Xoj+t6+t7-2*Yoi*Yoj+t10+t11-2*Zoi*Zoj+t14)
      t17 = t16+0.1E-98
      t18 = t17**2
      t20 = 1/t16
      t22 = Zoi-Zoj
      t27 = qo*qh*ec
      t28 = cos(ppj)
      t29 = cos(sj)
      t31 = cos(tj)
      t32 = sin(sj)
      t33 = t31*t32
      t34 = sin(ppj)
      t37 = (t28*t29-t33*t34)*ac
      t39 = t31*t29
      t42 = (t28*t32+t39*t34)*as
      t44 = (Xoi-t37-t42-Xoj)**2
      t48 = (-t34*t29-t33*t28)*ac
      t52 = (-t34*t32+t39*t28)*as
      t54 = (Yoi-t48-t52-Yoj)**2
      t55 = sin(tj)
      t57 = t55*t32*ac
      t59 = t55*t29*as
      t60 = Zoi-t57+t59-Zoj
      t61 = t60**2
      t63 = sqrt(t44+t54+t61)
      t65 = (t63+0.1E-98)**2
      t73 = (Xoi-t37+t42-Xoj)**2
      t75 = (Yoi-t48+t52-Yoj)**2
      t76 = Zoi-t57-t59-Zoj
      t77 = t76**2
      t79 = sqrt(t73+t75+t77)
      t81 = (t79+0.1E-98)**2
      t88 = cos(ppi)
      t89 = cos(si)
      t91 = cos(ti)
      t92 = sin(si)
      t93 = t91*t92
      t94 = sin(ppi)
      t97 = (t88*t89-t93*t94)*ac
      t99 = t91*t89
      t102 = (t88*t92+t99*t94)*as
      t104 = (t97+t102+Xoi-Xoj)**2
      t108 = (-t94*t89-t93*t88)*ac
      t112 = (-t94*t92+t99*t88)*as
      t114 = (t108+t112+Yoi-Yoj)**2
      t115 = sin(ti)
      t117 = t115*t92*ac
      t119 = t115*t89*as
      t120 = t117-t119+Zoi-Zoj
      t121 = t120**2
      t123 = sqrt(t104+t114+t121)
      t125 = (t123+0.1E-98)**2
      t132 = qh**2
      t133 = t132*ec
      t135 = (t97+t102+Xoi-t37-t42-Xoj)**2
      t137 = (t108+t112+Yoi-t48-t52-Yoj)**2
      t138 = t117-t119+Zoi-t57+t59-Zoj
      t139 = t138**2
      t141 = sqrt(t135+t137+t139)
      t143 = (t141+0.1E-98)**2
      t151 = (t97+t102+Xoi-t37+t42-Xoj)**2
      t153 = (t108+t112+Yoi-t48+t52-Yoj)**2
      t154 = t117-t119+Zoi-t57-t59-Zoj
      t155 = t154**2
      t157 = sqrt(t151+t153+t155)
      t159 = (t157+0.1E-98)**2
      t167 = (t97-t102+Xoi-Xoj)**2
      t169 = (t108-t112+Yoi-Yoj)**2
      t170 = t117+t119+Zoi-Zoj
      t171 = t170**2
      t173 = sqrt(t167+t169+t171)
      t175 = (t173+0.1E-98)**2
      t183 = (t97-t102+Xoi-t37-t42-Xoj)**2
      t185 = (t108-t112+Yoi-t48-t52-Yoj)**2
      t186 = t117+t119+Zoi-t57+t59-Zoj
      t187 = t186**2
      t189 = sqrt(t183+t185+t187)
      t191 = (t189+0.1E-98)**2
      t199 = (t97-t102+Xoi-t37+t42-Xoj)**2
      t201 = (t108-t112+Yoi-t48+t52-Yoj)**2
      t202 = t117+t119+Zoi-t57-t59-Zoj
      t203 = t202**2
      t205 = sqrt(t199+t201+t203)
      t207 = (t205+0.1E-98)**2
      t214 = t18**2
      t216 = t214**2
      t220 = 2*t20*t22
      t229 = -t1*ec/t18*t20*t22-t27/t65/t63*t60-t27/t81/t79*t76-t27/t125
     #/t123*t120-t133/t143/t141*t138-t133/t159/t157*t154-t27/t175/t173*t
     #170-t133/t191/t189*t186-t133/t207/t205*t202-6*AA/t216/t214/t17*t22
     #0+3*BB/t214/t18/t17*t220
c
      dVdZoi3p=t229
c
      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      function dVdti3p(Xoi,Xoj,Yoi,Yoj,Zoi,Zoj,ppi,si,ti,ppj,sj,tj)
      implicit real*8(a-h,o-z)
      common/cpot/qo,qh,AA,BB,ac,as,ec     
c
      t2 = qh*qo*ec
      t3 = cos(ppi)
      t4 = cos(si)
      t6 = cos(ti)
      t7 = sin(si)
      t8 = t6*t7
      t9 = sin(ppi)
      t12 = (t3*t4-t8*t9)*ac
      t14 = t6*t4
      t17 = (t3*t7+t14*t9)*as
      t18 = t12+t17+Xoi-Xoj
      t19 = t18**2
      t23 = (-t9*t4-t8*t3)*ac
      t27 = (-t9*t7+t14*t3)*as
      t28 = t23+t27+Yoi-Yoj
      t29 = t28**2
      t30 = sin(ti)
      t31 = t30*t7
      t32 = t31*ac
      t33 = t30*t4
      t34 = t33*as
      t35 = t32-t34+Zoi-Zoj
      t36 = t35**2
      t38 = sqrt(t19+t29+t36)
      t40 = (t38+0.1E-98)**2
      t45 = t31*t9*ac
      t47 = t33*t9*as
      t48 = t45-t47
      t51 = t31*t3*ac
      t53 = t33*t3*as
      t54 = t51-t53
      t56 = t8*ac
      t57 = t14*as
      t58 = t56-t57
      t63 = qh**2
      t64 = t63*ec
      t65 = cos(ppj)
      t66 = cos(sj)
      t68 = cos(tj)
      t69 = sin(sj)
      t70 = t68*t69
      t71 = sin(ppj)
      t74 = (t65*t66-t70*t71)*ac
      t76 = t68*t66
      t79 = (t65*t69+t76*t71)*as
      t80 = t12+t17+Xoi-t74-t79-Xoj
      t81 = t80**2
      t85 = (-t71*t66-t70*t65)*ac
      t89 = (-t71*t69+t76*t65)*as
      t90 = t23+t27+Yoi-t85-t89-Yoj
      t91 = t90**2
      t92 = sin(tj)
      t94 = t92*t69*ac
      t96 = t92*t66*as
      t97 = t32-t34+Zoi-t94+t96-Zoj
      t98 = t97**2
      t100 = sqrt(t81+t91+t98)
      t102 = (t100+0.1E-98)**2
      t112 = t12+t17+Xoi-t74+t79-Xoj
      t113 = t112**2
      t114 = t23+t27+Yoi-t85+t89-Yoj
      t115 = t114**2
      t116 = t32-t34+Zoi-t94-t96-Zoj
      t117 = t116**2
      t119 = sqrt(t113+t115+t117)
      t121 = (t119+0.1E-98)**2
      t131 = t12-t17+Xoi-Xoj
      t132 = t131**2
      t133 = t23-t27+Yoi-Yoj
      t134 = t133**2
      t135 = t32+t34+Zoi-Zoj
      t136 = t135**2
      t138 = sqrt(t132+t134+t136)
      t140 = (t138+0.1E-98)**2
      t144 = t45+t47
      t146 = t51+t53
      t148 = t56+t57
      t153 = t12-t17+Xoi-t74-t79-Xoj
      t154 = t153**2
      t155 = t23-t27+Yoi-t85-t89-Yoj
      t156 = t155**2
      t157 = t32+t34+Zoi-t94+t96-Zoj
      t158 = t157**2
      t160 = sqrt(t154+t156+t158)
      t162 = (t160+0.1E-98)**2
      t172 = t12-t17+Xoi-t74+t79-Xoj
      t173 = t172**2
      t174 = t23-t27+Yoi-t85+t89-Yoj
      t175 = t174**2
      t176 = t32+t34+Zoi-t94-t96-Zoj
      t177 = t176**2
      t179 = sqrt(t173+t175+t177)
      t181 = (t179+0.1E-98)**2
      t192 = -t2/t40/t38*(t18*t48+t28*t54+t35*t58)-t64/t102/t100*(t80*t4
     #8+t90*t54+t97*t58)-t64/t121/t119*(t112*t48+t114*t54+t116*t58)-t2/t
     #140/t138*(t131*t144+t133*t146+t135*t148)-t64/t162/t160*(t153*t144+
     #t155*t146+t157*t148)-t64/t181/t179*(t172*t144+t174*t146+t176*t148)
c
      dVdti3p=t192
c
      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      function dVdsi3p(Xoi,Xoj,Yoi,Yoj,Zoi,Zoj,ppi,si,ti,ppj,sj,tj)
      implicit real*8(a-h,o-z)
      common/cpot/qo,qh,AA,BB,ac,as,ec     
c
      t2 = qh*qo*ec
      t3 = cos(ppi)
      t4 = cos(si)
      t6 = cos(ti)
      t7 = sin(si)
      t8 = t6*t7
      t9 = sin(ppi)
      t11 = t3*t4-t8*t9
      t12 = t11*ac
      t14 = t6*t4
      t16 = t3*t7+t14*t9
      t17 = t16*as
      t18 = t12+t17+Xoi-Xoj
      t19 = t18**2
      t22 = -t9*t4-t8*t3
      t23 = t22*ac
      t26 = -t9*t7+t14*t3
      t27 = t26*as
      t28 = t23+t27+Yoi-Yoj
      t29 = t28**2
      t30 = sin(ti)
      t31 = t30*t7
      t32 = t31*ac
      t33 = t30*t4
      t34 = t33*as
      t35 = t32-t34+Zoi-Zoj
      t36 = t35**2
      t38 = sqrt(t19+t29+t36)
      t40 = (t38+0.1E-98)**2
      t44 = -t16*ac
      t45 = t11*as
      t46 = t44+t45
      t48 = -t26*ac
      t49 = t22*as
      t50 = t48+t49
      t52 = t33*ac
      t53 = t31*as
      t54 = t52+t53
      t59 = qh**2
      t60 = t59*ec
      t61 = cos(ppj)
      t62 = cos(sj)
      t64 = cos(tj)
      t65 = sin(sj)
      t66 = t64*t65
      t67 = sin(ppj)
      t70 = (t61*t62-t66*t67)*ac
      t72 = t64*t62
      t75 = (t61*t65+t72*t67)*as
      t76 = t12+t17+Xoi-t70-t75-Xoj
      t77 = t76**2
      t81 = (-t67*t62-t66*t61)*ac
      t85 = (-t67*t65+t72*t61)*as
      t86 = t23+t27+Yoi-t81-t85-Yoj
      t87 = t86**2
      t88 = sin(tj)
      t90 = t88*t65*ac
      t92 = t88*t62*as
      t93 = t32-t34+Zoi-t90+t92-Zoj
      t94 = t93**2
      t96 = sqrt(t77+t87+t94)
      t98 = (t96+0.1E-98)**2
      t108 = t12+t17+Xoi-t70+t75-Xoj
      t109 = t108**2
      t110 = t23+t27+Yoi-t81+t85-Yoj
      t111 = t110**2
      t112 = t32-t34+Zoi-t90-t92-Zoj
      t113 = t112**2
      t115 = sqrt(t109+t111+t113)
      t117 = (t115+0.1E-98)**2
      t127 = t12-t17+Xoi-Xoj
      t128 = t127**2
      t129 = t23-t27+Yoi-Yoj
      t130 = t129**2
      t131 = t32+t34+Zoi-Zoj
      t132 = t131**2
      t134 = sqrt(t128+t130+t132)
      t136 = (t134+0.1E-98)**2
      t140 = t44-t45
      t142 = t48-t49
      t144 = t52-t53
      t149 = t12-t17+Xoi-t70-t75-Xoj
      t150 = t149**2
      t151 = t23-t27+Yoi-t81-t85-Yoj
      t152 = t151**2
      t153 = t32+t34+Zoi-t90+t92-Zoj
      t154 = t153**2
      t156 = sqrt(t150+t152+t154)
      t158 = (t156+0.1E-98)**2
      t168 = t12-t17+Xoi-t70+t75-Xoj
      t169 = t168**2
      t170 = t23-t27+Yoi-t81+t85-Yoj
      t171 = t170**2
      t172 = t32+t34+Zoi-t90-t92-Zoj
      t173 = t172**2
      t175 = sqrt(t169+t171+t173)
      t177 = (t175+0.1E-98)**2
      t188 = -t2/t40/t38*(t18*t46+t28*t50+t35*t54)-t60/t98/t96*(t76*t46+
     #t86*t50+t93*t54)-t60/t117/t115*(t108*t46+t110*t50+t112*t54)-t2/t13
     #6/t134*(t127*t140+t129*t142+t131*t144)-t60/t158/t156*(t149*t140+t1
     #51*t142+t153*t144)-t60/t177/t175*(t168*t140+t170*t142+t172*t144)
c
      dVdsi3p=t188
c
      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      function dVdppi3p(Xoi,Xoj,Yoi,Yoj,Zoi,Zoj,ppi,si,ti,ppj,sj,tj)
      implicit real*8(a-h,o-z)
      common/cpot/qo,qh,AA,BB,ac,as,ec      
c
      t2 = qh*qo*ec
      t3 = cos(ppi)
      t4 = cos(si)
      t6 = cos(ti)
      t7 = sin(si)
      t8 = t6*t7
      t9 = sin(ppi)
      t11 = t3*t4-t8*t9
      t12 = t11*ac
      t14 = t6*t4
      t16 = t3*t7+t14*t9
      t17 = t16*as
      t18 = t12+t17+Xoi-Xoj
      t19 = t18**2
      t23 = (-t9*t4-t8*t3)*ac
      t27 = (-t9*t7+t14*t3)*as
      t28 = t23+t27+Yoi-Yoj
      t29 = t28**2
      t30 = sin(ti)
      t32 = t30*t7*ac
      t34 = t30*t4*as
      t36 = (t32-t34+Zoi-Zoj)**2
      t38 = sqrt(t19+t29+t36)
      t40 = (t38+0.1E-98)**2
      t44 = t23+t27
      t46 = -t11*ac
      t47 = -t16*as
      t48 = t46+t47
      t53 = qh**2
      t54 = t53*ec
      t55 = cos(ppj)
      t56 = cos(sj)
      t58 = cos(tj)
      t59 = sin(sj)
      t60 = t58*t59
      t61 = sin(ppj)
      t64 = (t55*t56-t60*t61)*ac
      t66 = t58*t56
      t69 = (t55*t59+t66*t61)*as
      t70 = t12+t17+Xoi-t64-t69-Xoj
      t71 = t70**2
      t75 = (-t61*t56-t60*t55)*ac
      t79 = (-t61*t59+t66*t55)*as
      t80 = t23+t27+Yoi-t75-t79-Yoj
      t81 = t80**2
      t82 = sin(tj)
      t84 = t82*t59*ac
      t86 = t82*t56*as
      t88 = (t32-t34+Zoi-t84+t86-Zoj)**2
      t90 = sqrt(t71+t81+t88)
      t92 = (t90+0.1E-98)**2
      t101 = t12+t17+Xoi-t64+t69-Xoj
      t102 = t101**2
      t103 = t23+t27+Yoi-t75+t79-Yoj
      t104 = t103**2
      t106 = (t32-t34+Zoi-t84-t86-Zoj)**2
      t108 = sqrt(t102+t104+t106)
      t110 = (t108+0.1E-98)**2
      t119 = t12-t17+Xoi-Xoj
      t120 = t119**2
      t121 = t23-t27+Yoi-Yoj
      t122 = t121**2
      t124 = (t32+t34+Zoi-Zoj)**2
      t126 = sqrt(t120+t122+t124)
      t128 = (t126+0.1E-98)**2
      t132 = t23-t27
      t134 = t46-t47
      t139 = t12-t17+Xoi-t64-t69-Xoj
      t140 = t139**2
      t141 = t23-t27+Yoi-t75-t79-Yoj
      t142 = t141**2
      t144 = (t32+t34+Zoi-t84+t86-Zoj)**2
      t146 = sqrt(t140+t142+t144)
      t148 = (t146+0.1E-98)**2
      t157 = t12-t17+Xoi-t64+t69-Xoj
      t158 = t157**2
      t159 = t23-t27+Yoi-t75+t79-Yoj
      t160 = t159**2
      t162 = (t32+t34+Zoi-t84-t86-Zoj)**2
      t164 = sqrt(t158+t160+t162)
      t166 = (t164+0.1E-98)**2
      t176 = -t2/t40/t38*(t18*t44+t28*t48)-t54/t92/t90*(t70*t44+t80*t48)
     #-t54/t110/t108*(t101*t44+t103*t48)-t2/t128/t126*(t119*t132+t121*t1
     #34)-t54/t148/t146*(t139*t132+t141*t134)-t54/t166/t164*(t157*t132+t
     #159*t134)
c
      dVdppi3p=t176 
c
      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      function dVdXoj3p(Xoi,Xoj,Yoi,Yoj,Zoi,Zoj,ppi,si,ti,ppj,sj,tj)
      implicit real*8(a-h,o-z)
      common/cpot/qo,qh,AA,BB,ac,as,ec      
c
      t1 = qo**2
      t3 = Xoi**2
      t6 = Xoj**2
      t7 = Yoi**2
      t10 = Yoj**2
      t11 = Zoi**2
      t14 = Zoj**2
      t16 = sqrt(t3-2*Xoi*Xoj+t6+t7-2*Yoi*Yoj+t10+t11-2*Zoi*Zoj+t14)
      t17 = t16+0.1E-98
      t18 = t17**2
      t20 = 1/t16
      t22 = -Xoi+Xoj
      t27 = qo*qh*ec
      t28 = cos(ppj)
      t29 = cos(sj)
      t31 = cos(tj)
      t32 = sin(sj)
      t33 = t31*t32
      t34 = sin(ppj)
      t37 = (t28*t29-t33*t34)*ac
      t39 = t31*t29
      t42 = (t28*t32+t39*t34)*as
      t43 = Xoi-t37-t42-Xoj
      t44 = t43**2
      t48 = (-t34*t29-t33*t28)*ac
      t52 = (-t34*t32+t39*t28)*as
      t54 = (Yoi-t48-t52-Yoj)**2
      t55 = sin(tj)
      t57 = t55*t32*ac
      t59 = t55*t29*as
      t61 = (Zoi-t57+t59-Zoj)**2
      t63 = sqrt(t44+t54+t61)
      t65 = (t63+0.1E-98)**2
      t72 = Xoi-t37+t42-Xoj
      t73 = t72**2
      t75 = (Yoi-t48+t52-Yoj)**2
      t77 = (Zoi-t57-t59-Zoj)**2
      t79 = sqrt(t73+t75+t77)
      t81 = (t79+0.1E-98)**2
      t88 = cos(ppi)
      t89 = cos(si)
      t91 = cos(ti)
      t92 = sin(si)
      t93 = t91*t92
      t94 = sin(ppi)
      t97 = (t88*t89-t93*t94)*ac
      t99 = t91*t89
      t102 = (t88*t92+t99*t94)*as
      t103 = t97+t102+Xoi-Xoj
      t104 = t103**2
      t108 = (-t94*t89-t93*t88)*ac
      t112 = (-t94*t92+t99*t88)*as
      t114 = (t108+t112+Yoi-Yoj)**2
      t115 = sin(ti)
      t117 = t115*t92*ac
      t119 = t115*t89*as
      t121 = (t117-t119+Zoi-Zoj)**2
      t123 = sqrt(t104+t114+t121)
      t125 = (t123+0.1E-98)**2
      t132 = qh**2
      t133 = t132*ec
      t134 = t97+t102+Xoi-t37-t42-Xoj
      t135 = t134**2
      t137 = (t108+t112+Yoi-t48-t52-Yoj)**2
      t139 = (t117-t119+Zoi-t57+t59-Zoj)**2
      t141 = sqrt(t135+t137+t139)
      t143 = (t141+0.1E-98)**2
      t150 = t97+t102+Xoi-t37+t42-Xoj
      t151 = t150**2
      t153 = (t108+t112+Yoi-t48+t52-Yoj)**2
      t155 = (t117-t119+Zoi-t57-t59-Zoj)**2
      t157 = sqrt(t151+t153+t155)
      t159 = (t157+0.1E-98)**2
      t166 = t97-t102+Xoi-Xoj
      t167 = t166**2
      t169 = (t108-t112+Yoi-Yoj)**2
      t171 = (t117+t119+Zoi-Zoj)**2
      t173 = sqrt(t167+t169+t171)
      t175 = (t173+0.1E-98)**2
      t182 = t97-t102+Xoi-t37-t42-Xoj
      t183 = t182**2
      t185 = (t108-t112+Yoi-t48-t52-Yoj)**2
      t187 = (t117+t119+Zoi-t57+t59-Zoj)**2
      t189 = sqrt(t183+t185+t187)
      t191 = (t189+0.1E-98)**2
      t198 = t97-t102+Xoi-t37+t42-Xoj
      t199 = t198**2
      t201 = (t108-t112+Yoi-t48+t52-Yoj)**2
      t203 = (t117+t119+Zoi-t57-t59-Zoj)**2
      t205 = sqrt(t199+t201+t203)
      t207 = (t205+0.1E-98)**2
      t214 = t18**2
      t216 = t214**2
      t220 = 2*t20*t22
      t229 = -t1*ec/t18*t20*t22+t27/t65/t63*t43+t27/t81/t79*t72+t27/t125
     #/t123*t103+t133/t143/t141*t134+t133/t159/t157*t150+t27/t175/t173*t
     #166+t133/t191/t189*t182+t133/t207/t205*t198-6*AA/t216/t214/t17*t22
     #0+3*BB/t214/t18/t17*t220
c
      dVdXoj3p=t229
c
      return
      end
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      function dVdYoj3p(Xoi,Xoj,Yoi,Yoj,Zoi,Zoj,ppi,si,ti,ppj,sj,tj)
      implicit real*8(a-h,o-z)
      common/cpot/qo,qh,AA,BB,ac,as,ec      
c
      t1 = qo**2
      t3 = Xoi**2
      t6 = Xoj**2
      t7 = Yoi**2
      t10 = Yoj**2
      t11 = Zoi**2
      t14 = Zoj**2
      t16 = sqrt(t3-2*Xoi*Xoj+t6+t7-2*Yoi*Yoj+t10+t11-2*Zoi*Zoj+t14)
      t17 = t16+0.1E-98
      t18 = t17**2
      t20 = 1/t16
      t22 = -Yoi+Yoj
      t27 = qo*qh*ec
      t28 = cos(ppj)
      t29 = cos(sj)
      t31 = cos(tj)
      t32 = sin(sj)
      t33 = t31*t32
      t34 = sin(ppj)
      t37 = (t28*t29-t33*t34)*ac
      t39 = t31*t29
      t42 = (t28*t32+t39*t34)*as
      t44 = (Xoi-t37-t42-Xoj)**2
      t48 = (-t34*t29-t33*t28)*ac
      t52 = (-t34*t32+t39*t28)*as
      t53 = Yoi-t48-t52-Yoj
      t54 = t53**2
      t55 = sin(tj)
      t57 = t55*t32*ac
      t59 = t55*t29*as
      t61 = (Zoi-t57+t59-Zoj)**2
      t63 = sqrt(t44+t54+t61)
      t65 = (t63+0.1E-98)**2
      t73 = (Xoi-t37+t42-Xoj)**2
      t74 = Yoi-t48+t52-Yoj
      t75 = t74**2
      t77 = (Zoi-t57-t59-Zoj)**2
      t79 = sqrt(t73+t75+t77)
      t81 = (t79+0.1E-98)**2
      t88 = cos(ppi)
      t89 = cos(si)
      t91 = cos(ti)
      t92 = sin(si)
      t93 = t91*t92
      t94 = sin(ppi)
      t97 = (t88*t89-t93*t94)*ac
      t99 = t91*t89
      t102 = (t88*t92+t99*t94)*as
      t104 = (t97+t102+Xoi-Xoj)**2
      t108 = (-t94*t89-t93*t88)*ac
      t112 = (-t94*t92+t99*t88)*as
      t113 = t108+t112+Yoi-Yoj
      t114 = t113**2
      t115 = sin(ti)
      t117 = t115*t92*ac
      t119 = t115*t89*as
      t121 = (t117-t119+Zoi-Zoj)**2
      t123 = sqrt(t104+t114+t121)
      t125 = (t123+0.1E-98)**2
      t132 = qh**2
      t133 = t132*ec
      t135 = (t97+t102+Xoi-t37-t42-Xoj)**2
      t136 = t108+t112+Yoi-t48-t52-Yoj
      t137 = t136**2
      t139 = (t117-t119+Zoi-t57+t59-Zoj)**2
      t141 = sqrt(t135+t137+t139)
      t143 = (t141+0.1E-98)**2
      t151 = (t97+t102+Xoi-t37+t42-Xoj)**2
      t152 = t108+t112+Yoi-t48+t52-Yoj
      t153 = t152**2
      t155 = (t117-t119+Zoi-t57-t59-Zoj)**2
      t157 = sqrt(t151+t153+t155)
      t159 = (t157+0.1E-98)**2
      t167 = (t97-t102+Xoi-Xoj)**2
      t168 = t108-t112+Yoi-Yoj
      t169 = t168**2
      t171 = (t117+t119+Zoi-Zoj)**2
      t173 = sqrt(t167+t169+t171)
      t175 = (t173+0.1E-98)**2
      t183 = (t97-t102+Xoi-t37-t42-Xoj)**2
      t184 = t108-t112+Yoi-t48-t52-Yoj
      t185 = t184**2
      t187 = (t117+t119+Zoi-t57+t59-Zoj)**2
      t189 = sqrt(t183+t185+t187)
      t191 = (t189+0.1E-98)**2
      t199 = (t97-t102+Xoi-t37+t42-Xoj)**2
      t200 = t108-t112+Yoi-t48+t52-Yoj
      t201 = t200**2
      t203 = (t117+t119+Zoi-t57-t59-Zoj)**2
      t205 = sqrt(t199+t201+t203)
      t207 = (t205+0.1E-98)**2
      t214 = t18**2
      t216 = t214**2
      t220 = 2*t20*t22
      t229 = -t1*ec/t18*t20*t22+t27/t65/t63*t53+t27/t81/t79*t74+t27/t125
     #/t123*t113+t133/t143/t141*t136+t133/t159/t157*t152+t27/t175/t173*t
     #168+t133/t191/t189*t184+t133/t207/t205*t200-6*AA/t216/t214/t17*t22
     #0+3*BB/t214/t18/t17*t220
c
      dVdYoj3p=t229
c
      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      function dVdZoj3p(Xoi,Xoj,Yoi,Yoj,Zoi,Zoj,ppi,si,ti,ppj,sj,tj)
      implicit real*8(a-h,o-z)
      common/cpot/qo,qh,AA,BB,ac,as,ec     
c
      t1 = qo**2
      t3 = Xoi**2
      t6 = Xoj**2
      t7 = Yoi**2
      t10 = Yoj**2
      t11 = Zoi**2
      t14 = Zoj**2
      t16 = sqrt(t3-2*Xoi*Xoj+t6+t7-2*Yoi*Yoj+t10+t11-2*Zoi*Zoj+t14)
      t17 = t16+0.1E-98
      t18 = t17**2
      t20 = 1/t16
      t22 = -Zoi+Zoj
      t27 = qo*qh*ec
      t28 = cos(ppj)
      t29 = cos(sj)
      t31 = cos(tj)
      t32 = sin(sj)
      t33 = t31*t32
      t34 = sin(ppj)
      t37 = (t28*t29-t33*t34)*ac
      t39 = t31*t29
      t42 = (t28*t32+t39*t34)*as
      t44 = (Xoi-t37-t42-Xoj)**2
      t48 = (-t34*t29-t33*t28)*ac
      t52 = (-t34*t32+t39*t28)*as
      t54 = (Yoi-t48-t52-Yoj)**2
      t55 = sin(tj)
      t57 = t55*t32*ac
      t59 = t55*t29*as
      t60 = Zoi-t57+t59-Zoj
      t61 = t60**2
      t63 = sqrt(t44+t54+t61)
      t65 = (t63+0.1E-98)**2
      t73 = (Xoi-t37+t42-Xoj)**2
      t75 = (Yoi-t48+t52-Yoj)**2
      t76 = Zoi-t57-t59-Zoj
      t77 = t76**2
      t79 = sqrt(t73+t75+t77)
      t81 = (t79+0.1E-98)**2
      t88 = cos(ppi)
      t89 = cos(si)
      t91 = cos(ti)
      t92 = sin(si)
      t93 = t91*t92
      t94 = sin(ppi)
      t97 = (t88*t89-t93*t94)*ac
      t99 = t91*t89
      t102 = (t88*t92+t99*t94)*as
      t104 = (t97+t102+Xoi-Xoj)**2
      t108 = (-t94*t89-t93*t88)*ac
      t112 = (-t94*t92+t99*t88)*as
      t114 = (t108+t112+Yoi-Yoj)**2
      t115 = sin(ti)
      t117 = t115*t92*ac
      t119 = t115*t89*as
      t120 = t117-t119+Zoi-Zoj
      t121 = t120**2
      t123 = sqrt(t104+t114+t121)
      t125 = (t123+0.1E-98)**2
      t132 = qh**2
      t133 = t132*ec
      t135 = (t97+t102+Xoi-t37-t42-Xoj)**2
      t137 = (t108+t112+Yoi-t48-t52-Yoj)**2
      t138 = t117-t119+Zoi-t57+t59-Zoj
      t139 = t138**2
      t141 = sqrt(t135+t137+t139)
      t143 = (t141+0.1E-98)**2
      t151 = (t97+t102+Xoi-t37+t42-Xoj)**2
      t153 = (t108+t112+Yoi-t48+t52-Yoj)**2
      t154 = t117-t119+Zoi-t57-t59-Zoj
      t155 = t154**2
      t157 = sqrt(t151+t153+t155)
      t159 = (t157+0.1E-98)**2
      t167 = (t97-t102+Xoi-Xoj)**2
      t169 = (t108-t112+Yoi-Yoj)**2
      t170 = t117+t119+Zoi-Zoj
      t171 = t170**2
      t173 = sqrt(t167+t169+t171)
      t175 = (t173+0.1E-98)**2
      t183 = (t97-t102+Xoi-t37-t42-Xoj)**2
      t185 = (t108-t112+Yoi-t48-t52-Yoj)**2
      t186 = t117+t119+Zoi-t57+t59-Zoj
      t187 = t186**2
      t189 = sqrt(t183+t185+t187)
      t191 = (t189+0.1E-98)**2
      t199 = (t97-t102+Xoi-t37+t42-Xoj)**2
      t201 = (t108-t112+Yoi-t48+t52-Yoj)**2
      t202 = t117+t119+Zoi-t57-t59-Zoj
      t203 = t202**2
      t205 = sqrt(t199+t201+t203)
      t207 = (t205+0.1E-98)**2
      t214 = t18**2
      t216 = t214**2
      t220 = 2*t20*t22
      t229 = -t1*ec/t18*t20*t22+t27/t65/t63*t60+t27/t81/t79*t76+t27/t125
     #/t123*t120+t133/t143/t141*t138+t133/t159/t157*t154+t27/t175/t173*t
     #170+t133/t191/t189*t186+t133/t207/t205*t202-6*AA/t216/t214/t17*t22
     #0+3*BB/t214/t18/t17*t220
c
      dVdZoj3p=t229
c
      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      function dVdtj3p(Xoi,Xoj,Yoi,Yoj,Zoi,Zoj,ppi,si,ti,ppj,sj,tj)
      implicit real*8(a-h,o-z)
      common/cpot/qo,qh,AA,BB,ac,as,ec     
c
      t2 = qo*qh*ec
      t3 = cos(ppj)
      t4 = cos(sj)
      t6 = cos(tj)
      t7 = sin(sj)
      t8 = t6*t7
      t9 = sin(ppj)
      t12 = (t3*t4-t8*t9)*ac
      t14 = t6*t4
      t17 = (t3*t7+t14*t9)*as
      t18 = Xoi-t12-t17-Xoj
      t19 = t18**2
      t23 = (-t9*t4-t8*t3)*ac
      t27 = (-t9*t7+t14*t3)*as
      t28 = Yoi-t23-t27-Yoj
      t29 = t28**2
      t30 = sin(tj)
      t31 = t30*t7
      t32 = t31*ac
      t33 = t30*t4
      t34 = t33*as
      t35 = Zoi-t32+t34-Zoj
      t36 = t35**2
      t38 = sqrt(t19+t29+t36)
      t40 = (t38+0.1E-98)**2
      t45 = t31*t9*ac
      t47 = t33*t9*as
      t48 = -t45+t47
      t51 = t31*t3*ac
      t53 = t33*t3*as
      t54 = -t51+t53
      t56 = t8*ac
      t57 = t14*as
      t58 = -t56+t57
      t63 = Xoi-t12+t17-Xoj
      t64 = t63**2
      t65 = Yoi-t23+t27-Yoj
      t66 = t65**2
      t67 = Zoi-t32-t34-Zoj
      t68 = t67**2
      t70 = sqrt(t64+t66+t68)
      t72 = (t70+0.1E-98)**2
      t76 = -t45-t47
      t78 = -t51-t53
      t80 = -t56-t57
      t85 = qh**2
      t86 = t85*ec
      t87 = cos(ppi)
      t88 = cos(si)
      t90 = cos(ti)
      t91 = sin(si)
      t92 = t90*t91
      t93 = sin(ppi)
      t96 = (t87*t88-t92*t93)*ac
      t98 = t90*t88
      t101 = (t87*t91+t98*t93)*as
      t102 = t96+t101+Xoi-t12-t17-Xoj
      t103 = t102**2
      t107 = (-t93*t88-t92*t87)*ac
      t111 = (-t93*t91+t98*t87)*as
      t112 = t107+t111+Yoi-t23-t27-Yoj
      t113 = t112**2
      t114 = sin(ti)
      t116 = t114*t91*ac
      t118 = t114*t88*as
      t119 = t116-t118+Zoi-t32+t34-Zoj
      t120 = t119**2
      t122 = sqrt(t103+t113+t120)
      t124 = (t122+0.1E-98)**2
      t134 = t96+t101+Xoi-t12+t17-Xoj
      t135 = t134**2
      t136 = t107+t111+Yoi-t23+t27-Yoj
      t137 = t136**2
      t138 = t116-t118+Zoi-t32-t34-Zoj
      t139 = t138**2
      t141 = sqrt(t135+t137+t139)
      t143 = (t141+0.1E-98)**2
      t153 = t96-t101+Xoi-t12-t17-Xoj
      t154 = t153**2
      t155 = t107-t111+Yoi-t23-t27-Yoj
      t156 = t155**2
      t157 = t116+t118+Zoi-t32+t34-Zoj
      t158 = t157**2
      t160 = sqrt(t154+t156+t158)
      t162 = (t160+0.1E-98)**2
      t172 = t96-t101+Xoi-t12+t17-Xoj
      t173 = t172**2
      t174 = t107-t111+Yoi-t23+t27-Yoj
      t175 = t174**2
      t176 = t116+t118+Zoi-t32-t34-Zoj
      t177 = t176**2
      t179 = sqrt(t173+t175+t177)
      t181 = (t179+0.1E-98)**2
      t192 = -t2/t40/t38*(t18*t48+t28*t54+t35*t58)-t2/t72/t70*(t63*t76+t
     #65*t78+t67*t80)-t86/t124/t122*(t102*t48+t112*t54+t119*t58)-t86/t14
     #3/t141*(t134*t76+t136*t78+t138*t80)-t86/t162/t160*(t153*t48+t155*t
     #54+t157*t58)-t86/t181/t179*(t172*t76+t174*t78+t176*t80)
c
      dVdtj3p=t192
c
      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      function dVdsj3p(Xoi,Xoj,Yoi,Yoj,Zoi,Zoj,ppi,si,ti,ppj,sj,tj)
      implicit real*8(a-h,o-z)
      common/cpot/qo,qh,AA,BB,ac,as,ec      
c     
      t2 = qo*qh*ec
      t3 = cos(ppj)
      t4 = cos(sj)
      t6 = cos(tj)
      t7 = sin(sj)
      t8 = t6*t7
      t9 = sin(ppj)
      t11 = t3*t4-t8*t9
      t12 = t11*ac
      t14 = t6*t4
      t16 = t3*t7+t14*t9
      t17 = t16*as
      t18 = Xoi-t12-t17-Xoj
      t19 = t18**2
      t22 = -t9*t4-t8*t3
      t23 = t22*ac
      t26 = -t9*t7+t14*t3
      t27 = t26*as
      t28 = Yoi-t23-t27-Yoj
      t29 = t28**2
      t30 = sin(tj)
      t31 = t30*t7
      t32 = t31*ac
      t33 = t30*t4
      t34 = t33*as
      t35 = Zoi-t32+t34-Zoj
      t36 = t35**2
      t38 = sqrt(t19+t29+t36)
      t40 = (t38+0.1E-98)**2
      t44 = -t16*ac
      t45 = t11*as
      t46 = -t44-t45
      t48 = -t26*ac
      t49 = t22*as
      t50 = -t48-t49
      t52 = t33*ac
      t53 = t31*as
      t54 = -t52-t53
      t59 = Xoi-t12+t17-Xoj
      t60 = t59**2
      t61 = Yoi-t23+t27-Yoj
      t62 = t61**2
      t63 = Zoi-t32-t34-Zoj
      t64 = t63**2
      t66 = sqrt(t60+t62+t64)
      t68 = (t66+0.1E-98)**2
      t72 = -t44+t45
      t74 = -t48+t49
      t76 = -t52+t53
      t81 = qh**2
      t82 = t81*ec
      t83 = cos(ppi)
      t84 = cos(si)
      t86 = cos(ti)
      t87 = sin(si)
      t88 = t86*t87
      t89 = sin(ppi)
      t92 = (t83*t84-t88*t89)*ac
      t94 = t86*t84
      t97 = (t83*t87+t94*t89)*as
      t98 = t92+t97+Xoi-t12-t17-Xoj
      t99 = t98**2
      t103 = (-t89*t84-t88*t83)*ac
      t107 = (-t89*t87+t94*t83)*as
      t108 = t103+t107+Yoi-t23-t27-Yoj
      t109 = t108**2
      t110 = sin(ti)
      t112 = t110*t87*ac
      t114 = t110*t84*as
      t115 = t112-t114+Zoi-t32+t34-Zoj
      t116 = t115**2
      t118 = sqrt(t99+t109+t116)
      t120 = (t118+0.1E-98)**2
      t130 = t92+t97+Xoi-t12+t17-Xoj
      t131 = t130**2
      t132 = t103+t107+Yoi-t23+t27-Yoj
      t133 = t132**2
      t134 = t112-t114+Zoi-t32-t34-Zoj
      t135 = t134**2
      t137 = sqrt(t131+t133+t135)
      t139 = (t137+0.1E-98)**2
      t149 = t92-t97+Xoi-t12-t17-Xoj
      t150 = t149**2
      t151 = t103-t107+Yoi-t23-t27-Yoj
      t152 = t151**2
      t153 = t112+t114+Zoi-t32+t34-Zoj
      t154 = t153**2
      t156 = sqrt(t150+t152+t154)
      t158 = (t156+0.1E-98)**2
      t168 = t92-t97+Xoi-t12+t17-Xoj
      t169 = t168**2
      t170 = t103-t107+Yoi-t23+t27-Yoj
      t171 = t170**2
      t172 = t112+t114+Zoi-t32-t34-Zoj
      t173 = t172**2
      t175 = sqrt(t169+t171+t173)
      t177 = (t175+0.1E-98)**2
      t188 = -t2/t40/t38*(t18*t46+t28*t50+t35*t54)-t2/t68/t66*(t59*t72+t
     #61*t74+t63*t76)-t82/t120/t118*(t98*t46+t108*t50+t115*t54)-t82/t139
     #/t137*(t130*t72+t132*t74+t134*t76)-t82/t158/t156*(t149*t46+t151*t5
     #0+t153*t54)-t82/t177/t175*(t168*t72+t170*t74+t172*t76)
c
      dVdsj3p=t188
c
      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      function dVdppj3p(Xoi,Xoj,Yoi,Yoj,Zoi,Zoj,ppi,si,ti,ppj,sj,tj)
      implicit real*8(a-h,o-z)
      common/cpot/qo,qh,AA,BB,ac,as,ec     
c
      t2 = qo*qh*ec
      t3 = cos(ppj)
      t4 = cos(sj)
      t6 = cos(tj)
      t7 = sin(sj)
      t8 = t6*t7
      t9 = sin(ppj)
      t11 = t3*t4-t8*t9
      t12 = t11*ac
      t14 = t6*t4
      t16 = t3*t7+t14*t9
      t17 = t16*as
      t18 = Xoi-t12-t17-Xoj
      t19 = t18**2
      t23 = (-t9*t4-t8*t3)*ac
      t27 = (-t9*t7+t14*t3)*as
      t28 = Yoi-t23-t27-Yoj
      t29 = t28**2
      t30 = sin(tj)
      t32 = t30*t7*ac
      t34 = t30*t4*as
      t36 = (Zoi-t32+t34-Zoj)**2
      t38 = sqrt(t19+t29+t36)
      t40 = (t38+0.1E-98)**2
      t44 = -t23-t27
      t46 = -t11*ac
      t47 = -t16*as
      t48 = -t46-t47
      t53 = Xoi-t12+t17-Xoj
      t54 = t53**2
      t55 = Yoi-t23+t27-Yoj
      t56 = t55**2
      t58 = (Zoi-t32-t34-Zoj)**2
      t60 = sqrt(t54+t56+t58)
      t62 = (t60+0.1E-98)**2
      t66 = -t23+t27
      t68 = -t46+t47
      t73 = qh**2
      t74 = t73*ec
      t75 = cos(ppi)
      t76 = cos(si)
      t78 = cos(ti)
      t79 = sin(si)
      t80 = t78*t79
      t81 = sin(ppi)
      t84 = (t75*t76-t80*t81)*ac
      t86 = t78*t76
      t89 = (t75*t79+t86*t81)*as
      t90 = t84+t89+Xoi-t12-t17-Xoj
      t91 = t90**2
      t95 = (-t81*t76-t80*t75)*ac
      t99 = (-t81*t79+t86*t75)*as
      t100 = t95+t99+Yoi-t23-t27-Yoj
      t101 = t100**2
      t102 = sin(ti)
      t104 = t102*t79*ac
      t106 = t102*t76*as
      t108 = (t104-t106+Zoi-t32+t34-Zoj)**2
      t110 = sqrt(t91+t101+t108)
      t112 = (t110+0.1E-98)**2
      t121 = t84+t89+Xoi-t12+t17-Xoj
      t122 = t121**2
      t123 = t95+t99+Yoi-t23+t27-Yoj
      t124 = t123**2
      t126 = (t104-t106+Zoi-t32-t34-Zoj)**2
      t128 = sqrt(t122+t124+t126)
      t130 = (t128+0.1E-98)**2
      t139 = t84-t89+Xoi-t12-t17-Xoj
      t140 = t139**2
      t141 = t95-t99+Yoi-t23-t27-Yoj
      t142 = t141**2
      t144 = (t104+t106+Zoi-t32+t34-Zoj)**2
      t146 = sqrt(t140+t142+t144)
      t148 = (t146+0.1E-98)**2
      t157 = t84-t89+Xoi-t12+t17-Xoj
      t158 = t157**2
      t159 = t95-t99+Yoi-t23+t27-Yoj
      t160 = t159**2
      t162 = (t104+t106+Zoi-t32-t34-Zoj)**2
      t164 = sqrt(t158+t160+t162)
      t166 = (t164+0.1E-98)**2
      t176 = -t2/t40/t38*(t18*t44+t28*t48)-t2/t62/t60*(t53*t66+t55*t68)-
     #t74/t112/t110*(t90*t44+t100*t48)-t74/t130/t128*(t121*t66+t123*t68)
     #-t74/t148/t146*(t139*t44+t141*t48)-t74/t166/t164*(t157*t66+t159*t6
     #8)
c
      dVdppj3p=t176
c
      return
      end


















