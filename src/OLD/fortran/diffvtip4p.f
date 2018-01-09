       	subroutine diffvtip4p(pot_chos,natoms,nmono,xyz,v2tot,v3tot
     &                     ,v,diffv)
	implicit real*8(a-h,o-z)
	include "param.inc"
	include "common.inc"
        common/cpot/qm,qh,AA,BB,ac,as,bp,ec
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
        AA=At
        BB=Ct
        pi=acos(-1d0)
        ap=alph1*pi/360d0
        as=dioa1*sin(ap)
        ac=dioa2*cos(ap)
        bp=diob1
        qm=qia1
        qh=qia2
	ec=ecag
c	write(*,*)'AA,BB,ac,as,bp,qm,qh,ec',AA,BB,ac,as,bp,qm,qh,ec
c	write(*,*)'At,Ct,alph1,dioa1,dioa2,diob1,qia1,qia2,ecag'
c     & 	,At,Ct,alph1,dioa1,dioa2,diob1,qia1,qia2,ecag
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
	      v=v+vpot4p(Xoi,Xoj,Yoi,Yoj,Zoi,Zoj,ppi,si,ti,ppj,sj,tj)
c
              dvdx=dVdXoi4p(Xoi,Xoj,Yoi,Yoj,Zoi,Zoj
     &             ,ppi,si,ti,ppj,sj,tj)*4.182
              dvdy=dVdYoi4p(Xoi,Xoj,Yoi,Yoj,Zoi,Zoj
     &             ,ppi,si,ti,ppj,sj,tj)*4.182 
              dvdz=dVdZoi4p(Xoi,Xoj,Yoi,Yoj,Zoi,Zoj
     &             ,ppi,si,ti,ppj,sj,tj)*4.182
c
	      dv2dxo(i)=dv2dxo(i)+dvdx
	      dv2dyo(i)=dv2dyo(i)+dvdy
	      dv2dzo(i)=dv2dzo(i)+dvdz
	      dv2dt(i)=dv2dt(i)+dVdti4p(Xoi,Xoj,Yoi,Yoj,Zoi,Zoj
     &             ,ppi,si,ti,ppj,sj,tj)*4.182
	      dv2ds(i)=dv2ds(i)+dVdsi4p(Xoi,Xoj,Yoi,Yoj,Zoi,Zoj
     &             ,ppi,si,ti,ppj,sj,tj)*4.182
	      dv2dpp(i)=dv2dpp(i)+dVdppi4p(Xoi,Xoj,Yoi,Yoj,Zoi,Zoj
     &             ,ppi,si,ti,ppj,sj,tj)*4.182
              dv2dxo(j)=dv2dxo(j)-dvdx
	      dv2dyo(j)=dv2dyo(j)-dvdy
	      dv2dzo(j)=dv2dzo(j)-dvdz
	      dv2dt(j)=dv2dt(j)+dVdtj4p(Xoi,Xoj,Yoi,Yoj,Zoi,Zoj
     &             ,ppi,si,ti,ppj,sj,tj)*4.182
	      dv2ds(j)=dv2ds(j)+dVdsj4p(Xoi,Xoj,Yoi,Yoj,Zoi,Zoj
     &             ,ppi,si,ti,ppj,sj,tj)*4.182
	      dv2dpp(j)=dv2dpp(j)+dVdppj4p(Xoi,Xoj,Yoi,Yoj,Zoi,Zoj
     &             ,ppi,si,ti,ppj,sj,tj)*4.182
c
c
c
c              d1=dVdXoi4p(Xoi,Xoj,Yoi,Yoj,Zoi,Zoj
c     &             ,ppi,si,ti,ppj,sj,tj)
c              d2=dVdXoj4p(Xoi,Xoj,Yoi,Yoj,Zoi,Zoj
c     &             ,ppi,si,ti,ppj,sj,tj)
c              d3=dVdYoi4p(Xoi,Xoj,Yoi,Yoj,Zoi,Zoj
c     &             ,ppi,si,ti,ppj,sj,tj)
c              d4=dVdYoj4p(Xoi,Xoj,Yoi,Yoj,Zoi,Zoj
c     &             ,ppi,si,ti,ppj,sj,tj)
c              d5=dVdZoi4p(Xoi,Xoj,Yoi,Yoj,Zoi,Zoj
c     &             ,ppi,si,ti,ppj,sj,tj)
c              d6=dVdZoj4p(Xoi,Xoj,Yoi,Yoj,Zoi,Zoj
c     &             ,ppi,si,ti,ppj,sj,tj)
c              d7=dVdti4p(Xoi,Xoj,Yoi,Yoj,Zoi,Zoj
c     &             ,ppi,si,ti,ppj,sj,tj)
c              d8=dVdtj4p(Xoi,Xoj,Yoi,Yoj,Zoi,Zoj
c     &             ,ppi,si,ti,ppj,sj,tj)
c              d9=dVdsi4p(Xoi,Xoj,Yoi,Yoj,Zoi,Zoj
c     &             ,ppi,si,ti,ppj,sj,tj)
c              d10=dVdsj4p(Xoi,Xoj,Yoi,Yoj,Zoi,Zoj
c     &             ,ppi,si,ti,ppj,sj,tj)
c              d11=dVdppi4p(Xoi,Xoj,Yoi,Yoj,Zoi,Zoj
c     &             ,ppi,si,ti,ppj,sj,tj)
c              d12=dVdppj4p(Xoi,Xoj,Yoi,Yoj,Zoi,Zoj
c     &             ,ppi,si,ti,ppj,sj,tj)
c        
c              
c              write(*,*)'i,j,d1,d2',i,j,d1,d2
c              write(*,*)'d3,d4',d3,d4
c              write(*,*)'d5,d6',d5,d6
c              write(*,*)'d7,d8',d7,d8
c              write(*,*)'d9,d10',d9,d10
c              write(*,*)'d11,d12',d11,d12
c              write(*,*)'v',v
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
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      function vpot4p(Xoi,Xoj,Yoi,Yoj,Zoi,Zoj,ppi,si,ti,ppj,sj,tj)
      implicit real*8(a-h,o-z)
      common/cpot/qm,qh,AA,BB,ac,as,bp,ec      
c      
      t1 = qm**2
      t3 = sin(ppi)
      t4 = cos(si)
      t6 = cos(ti)
      t7 = sin(si)
      t8 = t6*t7
      t9 = cos(ppi)
      t11 = -t3*t4-t8*t9
      t12 = t11**2
      t13 = bp**2
      t14 = t12*t13
      t17 = t9*t4-t8*t3
      t18 = t17**2
      t19 = t18*t13
      t21 = 2*Xoi*Xoj
      t22 = Xoi**2
      t23 = Xoj**2
      t24 = Yoi**2
      t25 = sin(ti)
      t26 = t25*t7
      t27 = sin(tj)
      t29 = sin(sj)
      t33 = Zoi**2
      t34 = Yoj**2
      t35 = Zoj**2
      t36 = cos(ppj)
      t37 = cos(sj)
      t39 = cos(tj)
      t40 = t39*t29
      t41 = sin(ppj)
      t43 = t36*t37-t40*t41
      t44 = t43**2
      t45 = t44*t13
      t48 = -t41*t37-t40*t36
      t49 = t48**2
      t50 = t49*t13
      t52 = 2*Yoi*Yoj
      t54 = 2*Zoi*Zoj
      t55 = t25**2
      t56 = t7**2
      t57 = t55*t56
      t58 = t57*t13
      t59 = t14+t19-t21+t22+t23+t24-2*t26*t13*t27*t29+t33+t34+t35+t45+t5
     #0-t52-t54+t58
      t60 = t27**2
      t61 = t29**2
      t62 = t60*t61
      t63 = t62*t13
      t64 = t11*bp
      t66 = 2*t64*Yoi
      t70 = t43*bp
      t72 = 2*t70*Xoj
      t73 = Xoi*t43
      t75 = 2*t73*bp
      t77 = 2*t64*Yoj
      t78 = t48*bp
      t80 = 2*t78*Yoj
      t81 = Yoi*t48
      t83 = 2*t81*bp
      t84 = t17*bp
      t86 = 2*t84*Xoj
      t91 = 2*t84*Xoi
      t92 = Zoi*t27
      t95 = 2*t92*t29*bp
      t96 = bp*Zoj
      t97 = t26*t96
      t100 = t26*bp*Zoi
      t102 = t27*t29
      t104 = 2*t102*t96
      t105 = t63+t66-2*t11*t13*t48+t72-t75-t77+t80-t83-t86-2*t17*t13*t43
     #+t91-t95-2*t97+2*t100+t104
      t107 = sqrt(t59+t105)
      t111 = qm*qh
      t112 = t14+t19-t21+t22+t23+t24+t33+t34+t35-t52-t54
      t114 = t39*t37
      t116 = t36*t29+t114*t41
      t117 = t116*as
      t119 = 2*t117*Xoj
      t120 = ac**2
      t121 = t62*t120
      t122 = t44*t120
      t123 = t26*bp
      t124 = t102*ac
      t126 = 2*t123*t124
      t127 = t49*t120
      t130 = -t41*t29+t114*t36
      t131 = t130**2
      t132 = as**2
      t133 = t131*t132
      t134 = t58+t66-t77-t86+t119+t121+t91+t122-t126+t127+t133
      t136 = t27*t37
      t137 = t136*as
      t139 = 2*t123*t137
      t144 = 2*t60*t29*ac*t37*as
      t145 = t116**2
      t146 = t145*t132
      t147 = t37**2
      t149 = t60*t147*t132
      t151 = 2*t73*ac
      t152 = t43*ac
      t154 = 2*t152*Xoj
      t155 = t130*as
      t157 = 2*t155*Yoj
      t159 = 2*t84*t117
      t160 = t48*ac
      t162 = 2*t64*t160
      t164 = 2*t152*t117
      t167 = 2*t92*t29*ac
      t168 = t139-t144+t146+t149-t151+t154+t157-t159-t162+t164-t167
      t169 = t64*t155
      t170 = t84*t152
      t171 = t160*t155
      t173 = Yoi*t130*as
      t175 = Xoi*t116*as
      t176 = t160*Yoj
      t177 = t81*ac
      t178 = as*Zoj
      t179 = t136*t178
      t181 = t92*t37*as
      t182 = ac*Zoj
      t183 = t102*t182
      t184 = -t169-t170+t171-t173-t175+t176-t97+t100-t177-t179+t181+t183
      t187 = sqrt(t112+t134+t168+2*t184)
      t192 = t58+t66-t77-t86-t119+t121+t91+t122-t126+t127+t133
      t194 = -t139+t144+t146+t149-t151+t154-t157+t159-t162-t164-t167
      t195 = t169-t170-t171+t173+t175+t176-t97+t100-t177+t179-t181+t183
      t198 = sqrt(t112+t192+t194+2*t195)
      t203 = -t21+t22+t23+t24+t33+t34+t35+t45+t50-t52-t54
      t208 = 2*t55*t7*ac*t4*as
      t209 = t25*t4
      t210 = t209*as
      t213 = 2*t210*t102*bp
      t215 = t6*t4
      t217 = t9*t7+t215*t3
      t218 = t217**2
      t219 = t218*t132
      t222 = 2*t26*ac*Zoi
      t223 = t217*as
      t225 = 2*t223*t70
      t226 = t17*ac
      t228 = 2*t226*t223
      t229 = t63+t72-t75+t80-t83-t208+t213+t219+t222-t225+t228
      t233 = -t3*t7+t215*t9
      t234 = t233*as
      t236 = 2*t234*t78
      t238 = 2*t223*Xoj
      t240 = 2*t226*Xoj
      t241 = t233**2
      t242 = t241*t132
      t243 = t18*t120
      t244 = 2*t170
      t245 = t11*ac
      t247 = 2*t245*t234
      t250 = 2*t209*as*Zoi
      t252 = 2*t209*t178
      t253 = -t236-t126-t238-t240+t242+t243-t162-t244+t247-t250+t252
      t255 = 2*t26*t182
      t257 = 2*t245*Yoi
      t259 = 2*t226*Xoi
      t260 = t57*t120
      t262 = 2*t245*Yoj
      t263 = t12*t120
      t265 = 2*t234*Yoj
      t267 = 2*t234*Yoi
      t269 = 2*t223*Xoi
      t270 = t4**2
      t272 = t55*t270*t132
      t273 = -t255+t257+t259+t260-t262+t263-t265+t267+t269+t272-t95+t104
      t276 = sqrt(t203+t229+t253+t273)
      t281 = qh**2
      t282 = t281*ec
      t285 = 2*t11*t120*t48
      t288 = 2*t217*t132*t116
      t289 = -t21+t22+t23+t24+t33+t34+t35-t285-t52-t54-t288+t119+t121-t2
     #08+t219
      t290 = t222+t228+t122+t127+t133-t144+t146-t238+t149-t240-t151+t154
     #+t157+t242+t243+t164
      t292 = 2*t171
      t293 = 2*t173
      t294 = 2*t175
      t295 = 2*t176
      t297 = 2*t223*t152
      t299 = 2*t226*t117
      t301 = 2*t234*t160
      t302 = -t167+t292-t293-t294+t295+t247-t250+t252-t255-t297-t299-t30
     #1+t257+t259+t260-t262
      t305 = 2*t17*t120*t43
      t308 = 2*t233*t132*t130
      t311 = 2*t26*ac*t137
      t315 = 2*t209*t132*t27*t37
      t317 = 2*t210*t124
      t321 = 2*t26*t120*t27*t29
      t323 = 2*t245*t155
      t324 = 2*t177
      t325 = 2*t179
      t326 = 2*t181
      t327 = 2*t183
      t328 = -t305-t308+t263-t265+t267+t269+t272+t311-t315+t317-t321-t32
     #3-t324-t325+t326+t327
      t331 = sqrt(t289+t290+t302+t328)
      t335 = -t21+t22+t23+t24+t33+t34+t35-t285-t52-t54+t288-t119+t121-t2
     #08+t219
      t336 = t222+t228+t122+t127+t133+t144+t146-t238+t149-t240-t151+t154
     #-t157+t242+t243-t164
      t338 = -t167-t292+t293+t294+t295+t247-t250+t252-t255-t297+t299-t30
     #1+t257+t259+t260-t262
      t339 = -t305+t308+t263-t265+t267+t269+t272-t311+t315+t317-t321+t32
     #3-t324+t325-t326+t327
      t342 = sqrt(t335+t336+t338+t339)
      t346 = t63+t72-t75+t80-t83+t208-t213+t219+t222+t225-t228
      t348 = t236-t126+t238-t240+t242+t243-t162-t244-t247+t250-t252
      t349 = -t255+t257+t259+t260-t262+t263+t265-t267-t269+t272-t95+t104
      t352 = sqrt(t203+t346+t348+t349)
      t357 = -t21+t22+t23+t24+t33+t34+t35-t285-t52-t54+t288+t119+t121+t2
     #08+t219
      t358 = t222-t228+t122+t127+t133-t144+t146+t238+t149-t240-t151+t154
     #+t157+t242+t243+t164
      t360 = -t167+t292-t293-t294+t295-t247+t250-t252-t255+t297-t299+t30
     #1+t257+t259+t260-t262
      t361 = -t305+t308+t263+t265-t267-t269+t272+t311+t315-t317-t321-t32
     #3-t324-t325+t326+t327
      t364 = sqrt(t357+t358+t360+t361)
      t368 = -t21+t22+t23+t24+t33+t34+t35-t285-t52-t54-t288-t119+t121+t2
     #08+t219
      t369 = t222-t228+t122+t127+t133+t144+t146+t238+t149-t240-t151+t154
     #-t157+t242+t243-t164
      t371 = -t167-t292+t293+t294+t295-t247+t250-t252-t255+t297+t299+t30
     #1+t257+t259+t260-t262
      t372 = -t305-t308+t263+t265-t267-t269+t272-t311-t315-t317-t321+t32
     #3-t324+t325-t326+t327
      t375 = sqrt(t368+t369+t371+t372)
      t380 = sqrt(t22-t21+t23+t24-t52+t34+t33-t54+t35)
      t382 = (t380+0.1E-98)**2
      t383 = t382**2
      t384 = t383**2
      t391 = t1*ec/(t107+0.1E-98)+t111*ec/(t187+0.1E-98)+t111*ec/(t198+0
     #.1E-98)+t111*ec/(t276+0.1E-98)+t282/(t331+0.1E-98)+t282/(t342+0.1E
     #-98)+t111*ec/(t352+0.1E-98)+t282/(t364+0.1E-98)+t282/(t375+0.1E-98
     #)+AA/t384/t383-BB/t383/t382
c
      vpot4p=t391
c
      return
      end
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      function dVdXoi4p(Xoi,Xoj,Yoi,Yoj,Zoi,Zoj,ppi,si,ti,ppj,sj,tj)
      implicit real*8(a-h,o-z)
      common/cpot/qm,qh,AA,BB,ac,as,bp,ec
c
      t1 = qm**2
      t4 = 2*Yoi*Yoj
      t5 = sin(ppi)
      t6 = cos(si)
      t8 = cos(ti)
      t9 = sin(si)
      t10 = t8*t9
      t11 = cos(ppi)
      t13 = -t5*t6-t10*t11
      t14 = t13**2
      t15 = bp**2
      t16 = t14*t15
      t17 = sin(tj)
      t18 = Zoi*t17
      t19 = sin(sj)
      t22 = 2*t18*t19*bp
      t25 = t11*t6-t10*t5
      t26 = t25**2
      t27 = t26*t15
      t28 = t17*t19
      t29 = bp*Zoj
      t31 = 2*t28*t29
      t32 = sin(ppj)
      t33 = cos(sj)
      t35 = cos(tj)
      t36 = t35*t19
      t37 = cos(ppj)
      t39 = -t32*t33-t36*t37
      t40 = t39**2
      t41 = t40*t15
      t44 = t37*t33-t36*t32
      t45 = t44**2
      t46 = t45*t15
      t48 = 2*Xoi*Xoj
      t49 = Yoi*t39
      t51 = 2*t49*bp
      t52 = t17**2
      t53 = t19**2
      t54 = t52*t53
      t55 = t54*t15
      t56 = sin(ti)
      t57 = t56*t9
      t62 = t25*bp
      t64 = 2*t62*Xoi
      t65 = t56**2
      t66 = t9**2
      t67 = t65*t66
      t68 = t67*t15
      t69 = Zoj**2
      t70 = Xoi**2
      t71 = -t4+t16-t22+t27+t31+t41+t46-t48-t51+t55-2*t57*t15*t17*t19+t6
     #4+t68+t69+t70
      t74 = 2*t57*bp*Zoi
      t75 = Xoj**2
      t76 = Yoi**2
      t77 = Yoj**2
      t78 = Zoi**2
      t80 = 2*t62*Xoj
      t82 = 2*Zoi*Zoj
      t86 = Xoi*t44
      t88 = 2*t86*bp
      t89 = t13*bp
      t91 = 2*t89*Yoj
      t93 = 2*t89*Yoi
      t97 = t44*bp
      t99 = 2*t97*Xoj
      t100 = t39*bp
      t102 = 2*t100*Yoj
      t104 = 2*t57*t29
      t105 = t74+t75+t76+t77+t78-t80-t82-2*t25*t15*t44-t88-t91+t93-2*t13
     #*t15*t39+t99+t102-t104
      t107 = sqrt(t71+t105)
      t109 = (t107+0.1E-98)**2
      t118 = qm*qh*ec
      t120 = 2*t49*ac
      t121 = -t4+t16+t27-t48+t64+t68+t69+t70+t74+t75-t120
      t122 = t44*ac
      t124 = 2*t122*Xoj
      t126 = 2*t86*ac
      t128 = t35*t33
      t130 = -t32*t19+t128*t37
      t133 = 2*Yoi*t130*as
      t136 = t37*t19+t128*t32
      t137 = t136*as
      t139 = 2*t137*Xoj
      t142 = 2*Xoi*t136*as
      t143 = t39*ac
      t145 = 2*t143*Yoj
      t146 = t124+t76+t77+t78-t80-t82-t126-t133+t139-t142+t145
      t148 = t130*as
      t150 = 2*t148*Yoj
      t151 = ac**2
      t152 = t54*t151
      t153 = t33**2
      t155 = as**2
      t156 = t52*t153*t155
      t157 = t57*bp
      t158 = t17*t33
      t159 = t158*as
      t161 = 2*t157*t159
      t162 = t28*ac
      t164 = 2*t157*t162
      t165 = t130**2
      t166 = t165*t155
      t167 = t45*t151
      t172 = 2*t52*t19*ac*t33*as
      t173 = t150+t152+t156-t91+t93-t104+t161-t164+t166+t167-t172
      t175 = 2*t62*t137
      t176 = t40*t151
      t177 = t136**2
      t178 = t177*t155
      t180 = 2*t89*t143
      t182 = 2*t62*t122
      t183 = ac*Zoj
      t185 = 2*t28*t183
      t187 = 2*t89*t148
      t189 = 2*t143*t148
      t192 = 2*t18*t19*ac
      t193 = as*Zoj
      t195 = 2*t158*t193
      t197 = 2*t122*t137
      t200 = 2*t18*t33*as
      t201 = -t175+t176+t178-t180-t182+t185-t187+t189-t192-t195+t197+t20
     #0
      t204 = sqrt(t121+t146+t173+t201)
      t206 = (t204+0.1E-98)**2
      t214 = t124+t76+t77+t78-t80-t82-t126+t133-t139+t142+t145
      t216 = -t150+t152+t156-t91+t93-t104-t161-t164+t166+t167+t172
      t217 = t175+t176+t178-t180-t182+t185+t187-t189-t192+t195-t197-t200
      t220 = sqrt(t121+t214+t216+t217)
      t222 = (t220+0.1E-98)**2
      t230 = -t4-t22+t31+t41+t46-t48-t51+t55+t69+t70+t75
      t231 = t25*ac
      t233 = 2*t231*Xoj
      t235 = t8*t6
      t237 = t11*t9+t235*t5
      t238 = t237*as
      t240 = 2*t238*Xoj
      t243 = -t5*t9+t235*t11
      t244 = t243*as
      t246 = 2*t244*Yoj
      t248 = 2*t231*Xoi
      t249 = t76+t77-t233+t78-t82-t240-t88+t99+t102-t246+t248
      t251 = t67*t151
      t252 = t13*ac
      t254 = 2*t252*Yoj
      t255 = t6**2
      t257 = t65*t255*t155
      t259 = 2*t244*Yoi
      t261 = 2*t252*Yoi
      t263 = 2*t244*t100
      t264 = t56*t6
      t267 = 2*t264*as*Zoi
      t269 = 2*t252*t244
      t270 = t251-t254+t257+t259+t261-t164-t180-t182-t263-t267+t269
      t272 = 2*t57*t183
      t275 = 2*t57*ac*Zoi
      t277 = 2*t238*t97
      t279 = 2*t231*t238
      t281 = 2*t264*t193
      t282 = t264*as
      t285 = 2*t282*t28*bp
      t286 = t237**2
      t287 = t286*t155
      t292 = 2*t65*t9*ac*t6*as
      t293 = t243**2
      t294 = t293*t155
      t295 = t14*t151
      t296 = t26*t151
      t298 = 2*t238*Xoi
      t299 = -t272+t275-t277+t279+t281+t285+t287-t292+t294+t295+t296+t29
     #8
      t302 = sqrt(t230+t249+t270+t299)
      t304 = (t302+0.1E-98)**2
      t312 = qh**2
      t313 = t312*ec
      t316 = 2*t25*t151*t44
      t319 = 2*t13*t151*t39
      t322 = 2*t243*t155*t130
      t325 = 2*t237*t155*t136
      t326 = -t4-t316-t319-t48+t69+t70+t75-t322-t120+t124-t325+t76+t77-t
     #233+t78
      t327 = -t82-t126-t133+t139-t142+t145+t150+t152+t156-t240-t246+t248
     #+t251-t254+t257+t259
      t329 = t261+t166+t167-t172+t176+t178+t185+t189-t192-t195+t197-t267
     #+t200+t269-t272+t275
      t331 = 2*t252*t148
      t333 = 2*t238*t122
      t335 = 2*t231*t137
      t337 = 2*t244*t143
      t339 = 2*t282*t162
      t342 = 2*t57*ac*t159
      t346 = 2*t264*t155*t17*t33
      t350 = 2*t57*t151*t17*t19
      t351 = t279+t281-t331+t287-t292+t294-t333-t335-t337+t295+t296+t339
     #+t342-t346+t298-t350
      t354 = sqrt(t326+t327+t329+t351)
      t356 = (t354+0.1E-98)**2
      t364 = -t4-t316-t319-t48+t69+t70+t75+t322-t120+t124+t325+t76+t77-t
     #233+t78
      t365 = -t82-t126+t133-t139+t142+t145-t150+t152+t156-t240-t246+t248
     #+t251-t254+t257+t259
      t367 = t261+t166+t167+t172+t176+t178+t185-t189-t192+t195-t197-t267
     #-t200+t269-t272+t275
      t368 = t279+t281+t331+t287-t292+t294-t333+t335-t337+t295+t296+t339
     #-t342+t346+t298-t350
      t371 = sqrt(t364+t365+t367+t368)
      t373 = (t371+0.1E-98)**2
      t381 = t76+t77-t233+t78-t82+t240-t88+t99+t102+t246+t248
      t383 = t251-t254+t257-t259+t261-t164-t180-t182+t263+t267-t269
      t384 = -t272+t275+t277-t279-t281-t285+t287+t292+t294+t295+t296-t29
     #8
      t387 = sqrt(t230+t381+t383+t384)
      t389 = (t387+0.1E-98)**2
      t397 = -t82-t126-t133+t139-t142+t145+t150+t152+t156+t240+t246+t248
     #+t251-t254+t257-t259
      t399 = t261+t166+t167-t172+t176+t178+t185+t189-t192-t195+t197+t267
     #+t200-t269-t272+t275
      t400 = -t279-t281-t331+t287+t292+t294+t333-t335+t337+t295+t296-t33
     #9+t342+t346-t298-t350
      t403 = sqrt(t364+t397+t399+t400)
      t405 = (t403+0.1E-98)**2
      t413 = -t82-t126+t133-t139+t142+t145-t150+t152+t156+t240+t246+t248
     #+t251-t254+t257-t259
      t415 = t261+t166+t167+t172+t176+t178+t185-t189-t192+t195-t197+t267
     #-t200-t269-t272+t275
      t416 = -t279-t281+t331+t287+t292+t294+t333+t335+t337+t295+t296-t33
     #9-t342-t346-t298-t350
      t419 = sqrt(t326+t413+t415+t416)
      t421 = (t419+0.1E-98)**2
      t430 = sqrt(t70-t48+t75+t76-t4+t77+t78-t82+t69)
      t431 = t430+0.1E-98
      t432 = t431**2
      t433 = t432**2
      t435 = t433**2
      t441 = 2/t430*(Xoi-Xoj)
      t450 = -t1*ec/t109/t107*(-Xoj+Xoi-t97+t62)-t118/t206/t204*(-Xoj+Xo
     #i+t62-t137-t122)-t118/t222/t220*(-Xoj+Xoi+t62+t137-t122)-t118/t304
     #/t302*(-Xoj+Xoi-t97+t238+t231)-t313/t356/t354*(-Xoj+Xoi+t238+t231-
     #t137-t122)-t313/t373/t371*(-Xoj+Xoi+t238+t231+t137-t122)-t118/t389
     #/t387*(-Xoj+Xoi-t97-t238+t231)-t313/t405/t403*(-Xoj+Xoi-t238+t231-
     #t137-t122)-t313/t421/t419*(-Xoj+Xoi-t238+t231+t137-t122)-6*AA/t435
     #/t433/t431*t441+3*BB/t433/t432/t431*t441
c
      dVdXoi4p=t450
c
      return
      end
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      function dVdYoi4p(Xoi,Xoj,Yoi,Yoj,Zoi,Zoj,ppi,si,ti,ppj,sj,tj)
      implicit real*8(a-h,o-z)
      common/cpot/qm,qh,AA,BB,ac,as,bp,ec
c
      t1 = qm**2
      t4 = 2*Yoi*Yoj
      t5 = cos(ppj)
      t6 = cos(sj)
      t8 = cos(tj)
      t9 = sin(sj)
      t10 = t8*t9
      t11 = sin(ppj)
      t13 = t5*t6-t10*t11
      t14 = t13**2
      t15 = bp**2
      t16 = t14*t15
      t18 = 2*Zoi*Zoj
      t19 = cos(ppi)
      t20 = cos(si)
      t22 = cos(ti)
      t23 = sin(si)
      t24 = t22*t23
      t25 = sin(ppi)
      t27 = t19*t20-t24*t25
      t28 = t27**2
      t29 = t28*t15
      t32 = -t11*t6-t10*t5
      t33 = t32**2
      t34 = t33*t15
      t36 = 2*Xoi*Xoj
      t37 = sin(tj)
      t38 = t37*t9
      t39 = bp*Zoj
      t41 = 2*t38*t39
      t42 = t27*bp
      t44 = 2*t42*Xoj
      t46 = 2*t42*Xoi
      t47 = Zoi*t37
      t50 = 2*t47*t9*bp
      t51 = sin(ti)
      t52 = t51*t23
      t57 = t32*bp
      t59 = 2*t57*Yoj
      t60 = Yoi*t32
      t62 = 2*t60*bp
      t68 = -t25*t20-t24*t19
      t69 = t68*bp
      t71 = 2*t69*Yoj
      t72 = -t4+t16-t18+t29+t34-t36+t41-t44+t46-t50-2*t52*t15*t37*t9+t59
     #-t62-2*t27*t15*t13-t71
      t74 = 2*t69*Yoi
      t78 = t13*bp
      t80 = 2*t78*Xoj
      t81 = t37**2
      t82 = t9**2
      t83 = t81*t82
      t84 = t83*t15
      t85 = t68**2
      t86 = t85*t15
      t87 = Xoi**2
      t88 = Xoj**2
      t89 = Zoj**2
      t90 = Yoj**2
      t91 = Zoi**2
      t92 = Yoi**2
      t93 = Xoi*t13
      t95 = 2*t93*bp
      t96 = t51**2
      t97 = t23**2
      t98 = t96*t97
      t99 = t98*t15
      t101 = 2*t52*t39
      t104 = 2*t52*bp*Zoi
      t105 = t74-2*t68*t15*t32+t80+t84+t86+t87+t88+t89+t90+t91+t92-t95+t
     #99-t101+t104
      t107 = sqrt(t72+t105)
      t109 = (t107+0.1E-98)**2
      t118 = qm*qh*ec
      t120 = t8*t6
      t122 = t5*t9+t120*t11
      t123 = t122*as
      t125 = 2*t123*Xoj
      t126 = t32*ac
      t128 = 2*t126*Yoj
      t130 = 2*t60*ac
      t131 = -t4-t18+t29+t125-t36-t44+t46+t128-t71+t74-t130
      t134 = 2*t47*t6*as
      t137 = 2*t47*t9*ac
      t140 = -t11*t9+t120*t5
      t141 = t140*as
      t143 = 2*t126*t141
      t144 = t86+t87+t88+t89+t90+t91+t92+t99+t134-t137+t143
      t146 = t37*t6
      t147 = as*Zoj
      t149 = 2*t146*t147
      t151 = 2*t69*t141
      t152 = t69*t126
      t153 = 2*t152
      t154 = ac*Zoj
      t156 = 2*t38*t154
      t159 = 2*Xoi*t122*as
      t161 = 2*t93*ac
      t162 = t6**2
      t164 = as**2
      t165 = t81*t162*t164
      t166 = t13*ac
      t167 = t42*t166
      t168 = 2*t167
      t170 = 2*t166*t123
      t172 = 2*t42*t123
      t173 = ac**2
      t174 = t14*t173
      t175 = -t149-t151-t153+t156-t159-t161+t165-t168+t170-t172+t174
      t176 = t140**2
      t177 = t176*t164
      t178 = t52*bp
      t179 = t38*ac
      t180 = t178*t179
      t181 = 2*t180
      t182 = t122**2
      t183 = t182*t164
      t184 = t146*as
      t186 = 2*t178*t184
      t187 = t33*t173
      t192 = 2*t81*t9*ac*t6*as
      t193 = t83*t173
      t195 = 2*t141*Yoj
      t197 = 2*t166*Xoj
      t200 = 2*Yoi*t140*as
      t201 = t177-t181+t183+t186+t187-t192+t193+t195+t197-t200-t101+t104
      t204 = sqrt(t131+t144+t175+t201)
      t206 = (t204+0.1E-98)**2
      t214 = -t4-t18+t29-t125-t36-t44+t46+t128-t71+t74-t130
      t215 = t86+t87+t88+t89+t90+t91+t92+t99-t134-t137-t143
      t217 = t149+t151-t153+t156+t159-t161+t165-t168-t170+t172+t174
      t218 = t177-t181+t183-t186+t187+t192+t193-t195+t197+t200-t101+t104
      t221 = sqrt(t214+t215+t217+t218)
      t223 = (t221+0.1E-98)**2
      t231 = t27*ac
      t233 = 2*t231*Xoj
      t234 = -t4+t16-t18+t34-t36+t41-t50-t233+t59-t62+t80
      t235 = t68*ac
      t237 = 2*t235*Yoi
      t239 = 2*t235*Yoj
      t241 = t22*t20
      t243 = -t25*t23+t241*t19
      t244 = t243**2
      t245 = t244*t164
      t246 = t84+t237-t239+t87+t88+t89+t90+t91+t92-t95+t245
      t247 = t234+t246
      t248 = t231*Xoi
      t251 = t19*t23+t241*t25
      t252 = t251*as
      t253 = t252*t78
      t254 = t52*t154
      t255 = t51*t20
      t256 = t255*t147
      t257 = t231*t252
      t258 = t243*as
      t259 = t258*t57
      t260 = t235*t258
      t262 = t52*ac*Zoi
      t263 = -t152-t167-t180+t248-t253-t254+t256+t257-t259+t260+t262
      t268 = 2*t96*t23*ac*t20*as
      t269 = t255*as
      t272 = 2*t269*t38*bp
      t275 = 2*t255*as*Zoi
      t276 = t28*t173
      t278 = 2*t252*Xoj
      t280 = 2*t252*Xoi
      t282 = 2*t258*Yoj
      t283 = t98*t173
      t284 = t251**2
      t285 = t284*t164
      t287 = 2*t258*Yoi
      t288 = t20**2
      t290 = t96*t288*t164
      t291 = t85*t173
      t292 = -t268+t272-t275+t276-t278+t280-t282+t283+t285+t287+t290+t29
     #1
      t295 = sqrt(t247+2*t263+t292)
      t297 = (t295+0.1E-98)**2
      t305 = qh**2
      t306 = t305*ec
      t309 = 2*t68*t173*t32
      t312 = 2*t251*t164*t122
      t313 = -t4-t309-t18+t125-t36-t312-t233+t128-t130+t237-t239+t87+t88
     #+t89+t90
      t314 = t91+t92+t245+t134-t137+t143-t149+t156-t159-t161+t165+t170+t
     #174+t177+t183+t187
      t316 = 2*t248
      t317 = 2*t254
      t319 = 2*t231*t123
      t321 = 2*t252*t166
      t323 = 2*t235*t141
      t324 = 2*t256
      t325 = 2*t257
      t326 = 2*t260
      t327 = 2*t262
      t328 = -t192+t193+t195+t197+t316-t200-t317-t319-t321-t323+t324+t32
     #5+t326+t327-t268-t275
      t332 = 2*t255*t164*t37*t6
      t335 = 2*t52*ac*t184
      t337 = 2*t269*t179
      t341 = 2*t52*t173*t37*t9
      t343 = 2*t258*t126
      t346 = 2*t243*t164*t140
      t349 = 2*t27*t173*t13
      t350 = t276-t278+t280-t282+t283+t285+t287+t290+t291-t332+t335+t337
     #-t341-t343-t346-t349
      t353 = sqrt(t313+t314+t328+t350)
      t355 = (t353+0.1E-98)**2
      t363 = -t4-t309-t18-t125-t36+t312-t233+t128-t130+t237-t239+t87+t88
     #+t89+t90
      t364 = t91+t92+t245-t134-t137-t143+t149+t156+t159-t161+t165-t170+t
     #174+t177+t183+t187
      t366 = t192+t193-t195+t197+t316+t200-t317+t319-t321+t323+t324+t325
     #+t326+t327-t268-t275
      t367 = t276-t278+t280-t282+t283+t285+t287+t290+t291+t332-t335+t337
     #-t341-t343+t346-t349
      t370 = sqrt(t363+t364+t366+t367)
      t372 = (t370+0.1E-98)**2
      t380 = -t152-t167-t180+t248+t253-t254-t256-t257+t259-t260+t262
      t381 = t268-t272+t275+t276+t278-t280+t282+t283+t285-t287+t290+t291
      t384 = sqrt(t247+2*t380+t381)
      t386 = (t384+0.1E-98)**2
      t394 = -t4-t309-t18+t125-t36+t312-t233+t128-t130+t237-t239+t87+t88
     #+t89+t90
      t396 = -t192+t193+t195+t197+t316-t200-t317-t319+t321-t323-t324-t32
     #5-t326+t327+t268+t275
      t397 = t276+t278-t280+t282+t283+t285-t287+t290+t291+t332+t335-t337
     #-t341+t343+t346-t349
      t400 = sqrt(t394+t314+t396+t397)
      t402 = (t400+0.1E-98)**2
      t410 = -t4-t309-t18-t125-t36-t312-t233+t128-t130+t237-t239+t87+t88
     #+t89+t90
      t412 = t192+t193-t195+t197+t316+t200-t317+t319+t321+t323-t324-t325
     #-t326+t327+t268+t275
      t413 = t276+t278-t280+t282+t283+t285-t287+t290+t291-t332-t335-t337
     #-t341+t343-t346-t349
      t416 = sqrt(t410+t364+t412+t413)
      t418 = (t416+0.1E-98)**2
      t427 = sqrt(t87-t36+t88+t92-t4+t90+t91-t18+t89)
      t428 = t427+0.1E-98
      t429 = t428**2
      t430 = t429**2
      t432 = t430**2
      t438 = 2/t427*(Yoi-Yoj)
      t447 = -t1*ec/t109/t107*(-Yoj+Yoi-t57+t69)-t118/t206/t204*(-Yoj+Yo
     #i+t69-t141-t126)-t118/t223/t221*(-Yoj+Yoi+t69+t141-t126)-t118/t297
     #/t295*(-Yoj+Yoi-t57+t258+t235)-t306/t355/t353*(-Yoj+Yoi+t258+t235-
     #t141-t126)-t306/t372/t370*(-Yoj+Yoi+t258+t235+t141-t126)-t118/t386
     #/t384*(-Yoj+Yoi-t57-t258+t235)-t306/t402/t400*(-Yoj+Yoi-t258+t235-
     #t141-t126)-t306/t418/t416*(-Yoj+Yoi-t258+t235+t141-t126)-6*AA/t432
     #/t430/t428*t438+3*BB/t430/t429/t428*t438
c
      dVdYoi4p=t447
c
      return
      end
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      function dVdZoi4p(Xoi,Xoj,Yoi,Yoj,Zoi,Zoj,ppi,si,ti,ppj,sj,tj)
      implicit real*8(a-h,o-z)
      common/cpot/qm,qh,AA,BB,ac,as,bp,ec
c
      t1 = qm**2
      t3 = sin(ppi)
      t4 = cos(si)
      t6 = cos(ti)
      t7 = sin(si)
      t8 = t6*t7
      t9 = cos(ppi)
      t11 = -t3*t4-t8*t9
      t12 = t11*bp
      t14 = 2*t12*Yoj
      t15 = sin(tj)
      t16 = t15**2
      t17 = sin(sj)
      t18 = t17**2
      t19 = t16*t18
      t20 = bp**2
      t21 = t19*t20
      t23 = 2*t12*Yoi
      t26 = t9*t4-t8*t3
      t28 = cos(ppj)
      t29 = cos(sj)
      t31 = cos(tj)
      t32 = t31*t17
      t33 = sin(ppj)
      t35 = t28*t29-t32*t33
      t38 = t26*bp
      t40 = 2*t38*Xoi
      t42 = 2*t38*Xoj
      t46 = -t33*t29-t32*t28
      t49 = Xoi*t35
      t51 = 2*t49*bp
      t52 = Yoi*t46
      t54 = 2*t52*bp
      t55 = t46*bp
      t57 = 2*t55*Yoj
      t58 = t35*bp
      t60 = 2*t58*Xoj
      t61 = sin(ti)
      t62 = t61**2
      t63 = t7**2
      t64 = t62*t63
      t65 = t64*t20
      t67 = 2*Zoi*Zoj
      t68 = t46**2
      t69 = t68*t20
      t70 = t26**2
      t71 = t70*t20
      t72 = -t14+t21+t23-2*t26*t20*t35+t40-t42-2*t11*t20*t46-t51-t54+t57
     #+t60+t65-t67+t69+t71
      t73 = t35**2
      t74 = t73*t20
      t75 = t61*t7
      t80 = Zoj**2
      t81 = Xoi**2
      t82 = Xoj**2
      t83 = Yoi**2
      t84 = Yoj**2
      t85 = Zoi**2
      t87 = 2*Xoi*Xoj
      t89 = 2*Yoi*Yoj
      t90 = t11**2
      t91 = t90*t20
      t94 = 2*t75*bp*Zoi
      t95 = Zoi*t15
      t98 = 2*t95*t17*bp
      t99 = t15*t17
      t100 = bp*Zoj
      t102 = 2*t99*t100
      t104 = 2*t75*t100
      t105 = t74-2*t75*t20*t15*t17+t80+t81+t82+t83+t84+t85-t87-t89+t91+t
     #94-t98+t102-t104
      t107 = sqrt(t72+t105)
      t109 = (t107+0.1E-98)**2
      t113 = t99*bp
      t114 = t75*bp
      t120 = qm*qh*ec
      t121 = ac**2
      t122 = t19*t121
      t123 = t29**2
      t125 = as**2
      t126 = t16*t123*t125
      t127 = t46*ac
      t129 = 2*t127*Yoj
      t130 = t35*ac
      t132 = 2*t130*Xoj
      t134 = t31*t29
      t136 = t28*t17+t134*t33
      t139 = 2*Xoi*t136*as
      t142 = -t33*t17+t134*t28
      t143 = t142*as
      t145 = 2*t143*Yoj
      t146 = -t14+t122+t23+t40-t42+t65+t126+t129+t132-t139+t145
      t147 = t136*as
      t149 = 2*t147*Xoj
      t152 = 2*Yoi*t142*as
      t154 = 2*t49*ac
      t156 = 2*t52*ac
      t157 = t149-t152-t154-t156-t67+t71+t80+t81+t82+t83+t84
      t159 = t99*ac
      t161 = 2*t114*t159
      t163 = 2*t12*t143
      t164 = t136**2
      t165 = t164*t125
      t166 = t15*t29
      t167 = t166*as
      t169 = 2*t114*t167
      t174 = 2*t16*t17*ac*t29*as
      t175 = t85-t87-t89+t91+t94-t104-t161-t163+t165+t169-t174
      t176 = t68*t121
      t177 = t142**2
      t178 = t177*t125
      t180 = 2*t127*t143
      t181 = t73*t121
      t183 = 2*t38*t147
      t184 = as*Zoj
      t186 = 2*t166*t184
      t189 = 2*t95*t29*as
      t191 = 2*t130*t147
      t193 = 2*t12*t127
      t195 = 2*t38*t130
      t198 = 2*t95*t17*ac
      t199 = ac*Zoj
      t201 = 2*t99*t199
      t202 = t176+t178+t180+t181-t183-t186+t189+t191-t193-t195-t198+t201
      t205 = sqrt(t146+t157+t175+t202)
      t207 = (t205+0.1E-98)**2
      t215 = -t14+t122+t23+t40-t42+t65+t126+t129+t132+t139-t145
      t216 = -t149+t152-t154-t156-t67+t71+t80+t81+t82+t83+t84
      t218 = t85-t87-t89+t91+t94-t104-t161+t163+t165-t169+t174
      t219 = t176+t178-t180+t181+t183+t186-t189-t191-t193-t195-t198+t201
      t222 = sqrt(t215+t216+t218+t219)
      t224 = (t222+0.1E-98)**2
      t233 = t4*t6
      t235 = t9*t7+t233*t3
      t236 = t235**2
      t237 = t236*t125
      t240 = 2*t75*ac*Zoi
      t241 = t235*as
      t243 = 2*t241*t58
      t244 = t61*t4
      t245 = t244*as
      t247 = 2*t245*t113
      t252 = 2*t62*t7*ac*t4*as
      t253 = t90*t121
      t254 = t26*ac
      t256 = 2*t254*t241
      t259 = -t3*t7+t233*t9
      t260 = t259*as
      t262 = 2*t260*t55
      t263 = t11*ac
      t265 = 2*t263*t260
      t268 = 2*t244*as*Zoi
      t270 = 2*t244*t184
      t271 = t237+t240-t243+t247-t252+t253+t256-t262+t265-t268+t270
      t273 = 2*t75*t199
      t274 = t70*t121
      t275 = t259**2
      t276 = t275*t125
      t278 = 2*t254*Xoj
      t280 = 2*t263*Yoi
      t282 = 2*t263*Yoj
      t283 = -t273+t274+t276-t278+t280-t282+t21-t51-t54+t57+t60
      t286 = 2*t241*Xoi
      t288 = 2*t260*Yoi
      t290 = 2*t254*Xoi
      t292 = 2*t260*Yoj
      t293 = t64*t121
      t295 = 2*t241*Xoj
      t296 = t4**2
      t298 = t62*t296*t125
      t299 = t286+t288+t290-t292-t67+t69+t293-t295+t298+t74+t80
      t300 = t81+t82+t83+t84+t85-t87-t89-t98+t102-t161-t193-t195
      t303 = sqrt(t271+t283+t299+t300)
      t305 = (t303+0.1E-98)**2
      t309 = t75*ac
      t314 = qh**2
      t315 = t314*ec
      t319 = 2*t244*t125*t15*t29
      t321 = 2*t260*t127
      t323 = 2*t254*t147
      t325 = 2*t241*t130
      t327 = 2*t309*t167
      t329 = 2*t245*t159
      t333 = 2*t75*t121*t15*t17
      t335 = 2*t263*t143
      t336 = t237-t319-t321-t323-t325+t327+t329-t333-t335+t240-t252+t253
     #+t256+t265-t268
      t339 = 2*t26*t121*t35
      t342 = 2*t235*t125*t136
      t343 = t270-t273+t274+t276-t278+t280-t282+t122-t339-t342+t126+t129
     #+t132-t139+t145+t149
      t347 = 2*t11*t121*t46
      t350 = 2*t259*t125*t142
      t351 = -t152-t154-t156+t286+t288+t290-t292-t67+t293-t295-t347-t350
     #+t298+t80+t81+t82
      t352 = t83+t84+t85-t87-t89+t165-t174+t176+t178+t180+t181-t186+t189
     #+t191-t198+t201
      t355 = sqrt(t336+t343+t351+t352)
      t357 = (t355+0.1E-98)**2
      t365 = t237+t319-t321+t323-t325-t327+t329-t333+t335+t240-t252+t253
     #+t256+t265-t268
      t366 = t270-t273+t274+t276-t278+t280-t282+t122-t339+t342+t126+t129
     #+t132+t139-t145-t149
      t368 = t152-t154-t156+t286+t288+t290-t292-t67+t293-t295-t347+t350+
     #t298+t80+t81+t82
      t369 = t83+t84+t85-t87-t89+t165+t174+t176+t178-t180+t181+t186-t189
     #-t191-t198+t201
      t372 = sqrt(t365+t366+t368+t369)
      t374 = (t372+0.1E-98)**2
      t382 = t237+t240+t243-t247+t252+t253-t256+t262-t265+t268-t270
      t384 = -t286-t288+t290+t292-t67+t69+t293+t295+t298+t74+t80
      t387 = sqrt(t382+t283+t384+t300)
      t389 = (t387+0.1E-98)**2
      t397 = t237+t319+t321-t323+t325+t327-t329-t333-t335+t240+t252+t253
     #-t256-t265+t268
      t398 = -t270-t273+t274+t276-t278+t280-t282+t122-t339+t342+t126+t12
     #9+t132-t139+t145+t149
      t400 = -t152-t154-t156-t286-t288+t290+t292-t67+t293+t295-t347+t350
     #+t298+t80+t81+t82
      t403 = sqrt(t397+t398+t400+t352)
      t405 = (t403+0.1E-98)**2
      t413 = t237-t319+t321+t323+t325-t327-t329-t333+t335+t240+t252+t253
     #-t256-t265+t268
      t414 = -t270-t273+t274+t276-t278+t280-t282+t122-t339-t342+t126+t12
     #9+t132+t139-t145-t149
      t416 = t152-t154-t156-t286-t288+t290+t292-t67+t293+t295-t347-t350+
     #t298+t80+t81+t82
      t419 = sqrt(t413+t414+t416+t369)
      t421 = (t419+0.1E-98)**2
      t430 = sqrt(t81-t87+t82+t83-t89+t84+t85-t67+t80)
      t431 = t430+0.1E-98
      t432 = t431**2
      t433 = t432**2
      t435 = t433**2
      t441 = 2/t430*(Zoi-Zoj)
      t450 = -t1*ec/t109/t107*(-t113+t114-Zoj+Zoi)-t120/t207/t205*(t114-
     #Zoj+Zoi-t159+t167)-t120/t224/t222*(t114-Zoj+Zoi-t159-t167)-t120/t3
     #05/t303*(-t113-Zoj+Zoi-t245+t309)-t315/t357/t355*(-Zoj+Zoi-t245+t3
     #09-t159+t167)-t315/t374/t372*(-Zoj+Zoi-t245+t309-t159-t167)-t120/t
     #389/t387*(-t113-Zoj+Zoi+t245+t309)-t315/t405/t403*(-Zoj+Zoi+t245+t
     #309-t159+t167)-t315/t421/t419*(-Zoj+Zoi+t245+t309-t159-t167)-6*AA/
     #t435/t433/t431*t441+3*BB/t433/t432/t431*t441
c
      dVdZoi4p=t450
c
      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      function dVdti4p(Xoi,Xoj,Yoi,Yoj,Zoi,Zoj,ppi,si,ti,ppj,sj,tj)
      implicit real*8(a-h,o-z)
      common/cpot/qm,qh,AA,BB,ac,as,bp,ec
c
      t1 = qm**2
      t3 = sin(ti)
      t4 = t3**2
      t5 = sin(si)
      t6 = t5**2
      t7 = t4*t6
      t8 = bp**2
      t9 = t7*t8
      t10 = t3*t5
      t11 = bp*Zoj
      t13 = 2*t10*t11
      t14 = Xoj**2
      t15 = Yoi**2
      t17 = 2*Yoi*Yoj
      t18 = cos(ppi)
      t19 = cos(si)
      t20 = t18*t19
      t21 = cos(ti)
      t22 = t21*t5
      t23 = sin(ppi)
      t25 = t20-t22*t23
      t26 = t25*bp
      t28 = 2*t26*Xoi
      t29 = t25*t8
      t30 = cos(ppj)
      t31 = cos(sj)
      t33 = cos(tj)
      t34 = sin(sj)
      t35 = t33*t34
      t36 = sin(ppj)
      t38 = t30*t31-t35*t36
      t43 = -t36*t31-t35*t30
      t44 = Yoi*t43
      t46 = 2*t44*bp
      t47 = t43*bp
      t49 = 2*t47*Yoj
      t50 = t38*bp
      t52 = 2*t50*Xoj
      t53 = t23*t19
      t55 = -t53-t22*t18
      t56 = t55**2
      t57 = t56*t8
      t58 = t55*t8
      t61 = t55*bp
      t63 = 2*t61*Yoj
      t64 = Xoi*t38
      t66 = 2*t64*bp
      t68 = 2*t61*Yoi
      t69 = t9-t13+t14+t15-t17+t28-2*t29*t38-t46+t49+t52+t57-2*t58*t43-t
     #63-t66+t68
      t70 = t38**2
      t71 = t70*t8
      t72 = bp*Zoi
      t74 = 2*t10*t72
      t75 = sin(tj)
      t76 = Zoi*t75
      t79 = 2*t76*t34*bp
      t81 = t8*t75*t34
      t85 = 2*Zoi*Zoj
      t87 = 2*Xoi*Xoj
      t88 = t25**2
      t89 = t88*t8
      t90 = t43**2
      t91 = t90*t8
      t92 = Zoj**2
      t93 = Zoi**2
      t94 = Yoj**2
      t96 = 2*t26*Xoj
      t97 = t75**2
      t98 = t34**2
      t99 = t97*t98
      t100 = t99*t8
      t101 = t75*t34
      t103 = 2*t101*t11
      t104 = Xoi**2
      t105 = t71+t74-t79-2*t10*t81-t85-t87+t89+t91+t92+t93+t94-t96+t100+
     #t103+t104
      t107 = sqrt(t69+t105)
      t109 = (t107+0.1E-98)**2
      t113 = t22*t11
      t114 = t22*t72
      t115 = t10*t23
      t116 = t29*t115
      t117 = t10*t18
      t118 = t58*t117
      t122 = t18*bp
      t124 = t10*t122*Yoi
      t126 = t10*t122*Yoj
      t130 = t23*bp
      t132 = t10*t130*Xoi
      t134 = t10*t130*Xoj
      t135 = t3*t6
      t137 = t135*t8*t21
      t139 = -t113+t114+t116+t118-t10*t18*t8*t43+t124-t126-t10*t23*t8*t3
     #8+t132-t134+t137-t22*t81
      t143 = qm*qh*ec
      t145 = t33*t31
      t147 = -t36*t34+t145*t30
      t150 = 2*Yoi*t147*as
      t151 = t31**2
      t153 = as**2
      t154 = t97*t151*t153
      t155 = t9-t150+t154-t13+t14+t15-t17+t28+t57-t63+t68
      t157 = 2*t64*ac
      t158 = t38*ac
      t160 = 2*t158*Xoj
      t161 = ac**2
      t162 = t99*t161
      t165 = t30*t34+t145*t36
      t166 = t165*as
      t168 = 2*t166*Xoj
      t169 = t147*as
      t171 = 2*t169*Yoj
      t174 = 2*Xoi*t165*as
      t175 = t43*ac
      t177 = 2*t175*Yoj
      t178 = -t157+t160+t74+t162+t168+t171-t174-t85-t87+t89+t177
      t181 = 2*t44*ac
      t182 = t10*bp
      t183 = t75*t31
      t184 = t183*as
      t186 = 2*t182*t184
      t191 = 2*t97*t34*ac*t31*as
      t192 = t70*t161
      t193 = as*Zoj
      t195 = 2*t183*t193
      t196 = t147**2
      t197 = t196*t153
      t198 = t92+t93+t94-t96+t104-t181+t186-t191+t192-t195+t197
      t200 = 2*t158*t166
      t202 = 2*t61*t175
      t203 = t101*ac
      t205 = 2*t182*t203
      t207 = 2*t26*t158
      t209 = 2*t61*t169
      t210 = t165**2
      t211 = t210*t153
      t212 = t90*t161
      t214 = 2*t175*t169
      t217 = 2*t76*t34*ac
      t220 = 2*t76*t31*as
      t221 = ac*Zoj
      t223 = 2*t101*t221
      t225 = 2*t26*t166
      t226 = t200-t202-t205-t207-t209+t211+t212+t214-t217+t220+t223-t225
      t229 = sqrt(t155+t178+t198+t226)
      t231 = (t229+0.1E-98)**2
      t237 = t115*bp*t165*as
      t239 = t117*t47*ac
      t241 = t115*t50*ac
      t244 = t117*bp*t147*as
      t245 = t22*bp
      t246 = t245*t203
      t247 = t245*t184
      t248 = -t113+t114+t116+t118+t124-t126+t132-t134+t137-t237-t239-t24
     #1-t244-t246+t247
      t251 = t9+t150+t154-t13+t14+t15-t17+t28+t57-t63+t68
      t252 = -t157+t160+t74+t162-t168-t171+t174-t85-t87+t89+t177
      t254 = t92+t93+t94-t96+t104-t181-t186+t191+t192+t195+t197
      t255 = -t200-t202-t205-t207+t209+t211+t212-t214-t217-t220+t223+t22
     #5
      t258 = sqrt(t251+t252+t254+t255)
      t260 = (t258+0.1E-98)**2
      t264 = -t113+t114+t116+t118+t124-t126+t132-t134+t137+t237-t239-t24
     #1+t244-t246-t247
      t267 = t25*ac
      t269 = 2*t267*Xoj
      t270 = t14+t15-t269-t17-t46+t49+t52-t66+t71-t79-t85
      t272 = t21*t19
      t274 = -t23*t5+t272*t18
      t275 = t274*as
      t277 = 2*t275*Yoj
      t280 = t18*t5+t272*t23
      t281 = t280*as
      t283 = 2*t281*Xoi
      t284 = t7*t161
      t285 = -t87+t91+t92+t93+t94-t277+t100+t103+t104+t283+t284
      t288 = 2*t275*Yoi
      t289 = t19**2
      t291 = t4*t289*t153
      t293 = 2*t281*Xoj
      t294 = t55*ac
      t296 = 2*t294*Yoj
      t298 = 2*t294*Yoi
      t300 = 2*t267*Xoi
      t302 = 2*t294*t275
      t304 = 2*t275*t47
      t305 = t288+t291-t293-t296+t298+t300-t202-t205-t207+t302-t304
      t307 = 2*t267*t281
      t308 = t3*t19
      t310 = 2*t308*t193
      t311 = as*Zoi
      t313 = 2*t308*t311
      t314 = t308*as
      t315 = t101*bp
      t317 = 2*t314*t315
      t322 = 2*t4*t5*ac*t19*as
      t323 = t274**2
      t324 = t323*t153
      t325 = t88*t161
      t326 = ac*Zoi
      t328 = 2*t10*t326
      t330 = 2*t10*t221
      t331 = t280**2
      t332 = t331*t153
      t334 = 2*t281*t50
      t335 = t56*t161
      t336 = t307+t310-t313+t317-t322+t324+t325+t328-t330+t332-t334+t335
      t339 = sqrt(t270+t285+t305+t336)
      t341 = (t339+0.1E-98)**2
      t345 = t18*ac
      t348 = 2*t10*t345*Yoj
      t349 = t272*as
      t351 = 2*t349*t315
      t355 = 2*t115*ac*t280*as
      t358 = 2*t10*t345*Yoi
      t359 = t23*as
      t362 = 2*t308*t359*Xoi
      t363 = t18*as
      t366 = 2*t308*t363*Yoi
      t367 = t23*ac
      t370 = 2*t10*t367*Xoi
      t372 = 2*t272*t193
      t376 = 2*t267*t3*t53*as
      t377 = t308*t23
      t378 = as*t38
      t381 = 2*t377*t378*bp
      t382 = t308*t18
      t383 = as*t43
      t386 = 2*t382*t383*bp
      t390 = 2*t294*t3*t20*as
      t394 = 2*t117*ac*t274*as
      t395 = t10*ac
      t397 = 4*t395*t349
      t398 = -t348+t351+t355+t358-t362-t366+t370+t372-t376+t381+t386-t39
     #0+t394-t397
      t399 = t274*t153
      t400 = t399*t382
      t401 = t272*t311
      t402 = t22*t221
      t405 = t3*t289*t153*t21
      t406 = t22*t326
      t408 = t135*t161*t21
      t409 = t25*t161
      t410 = t409*t115
      t411 = t280*t153
      t412 = t411*t377
      t413 = t55*t161
      t414 = t413*t117
      t416 = t10*t367*Xoj
      t418 = t308*t359*Xoj
      t420 = t308*t363*Yoj
      t421 = -t400-t241-t239-t246-t401-t402+t405+t406+t408+t410-t412+t41
     #4-t416+t418+t420
      t425 = qh**2
      t426 = t425*ec
      t428 = 2*t409*t38
      t430 = 2*t275*t175
      t432 = 2*t395*t184
      t434 = t153*t75*t31
      t436 = 2*t308*t434
      t438 = 2*t314*t203
      t440 = t161*t75*t34
      t442 = 2*t10*t440
      t444 = 2*t411*t165
      t446 = 2*t413*t43
      t448 = 2*t399*t147
      t449 = -t428-t430+t432-t436-t150+t438-t442+t154+t14+t15-t269-t17-t
     #444-t446-t448
      t450 = -t157+t160+t162+t168+t171-t174-t85-t87+t177+t92+t93+t94-t27
     #7+t104+t283+t284
      t452 = t288+t291-t293-t296+t298+t300-t181-t191+t192-t195+t197+t200
     #+t211+t212+t214-t217
      t454 = 2*t294*t169
      t456 = 2*t267*t166
      t458 = 2*t281*t158
      t459 = t220+t223+t302+t307+t310-t313-t322+t324+t325+t328-t330+t332
     #+t335-t454-t456-t458
      t462 = sqrt(t449+t450+t452+t459)
      t464 = (t462+0.1E-98)**2
      t468 = 2*t400
      t470 = 2*t22*t440
      t474 = 2*t10*t23*t161*t38
      t478 = 2*t308*t23*t153*t165
      t482 = 2*t117*ac*t147*as
      t485 = 2*t377*t378*ac
      t486 = -t348+t355+t358-t362-t366+t370+t372-t376-t390+t394-t397-t46
     #8-t470-t474+t478-t482+t485
      t487 = t272*t434
      t490 = t10*t18*t161*t43
      t493 = t308*t18*t153*t147
      t495 = t22*ac*t184
      t496 = t349*t203
      t498 = t382*t383*ac
      t501 = t115*ac*t165*as
      t502 = -t487-t490+t493+t495+t496+t498-t501-t401-t402+t405+t406+t40
     #8+t410-t412+t414-t416+t418+t420
      t506 = -t428-t430-t432+t436+t150+t438-t442+t154+t14+t15-t269-t17+t
     #444-t446+t448
      t507 = -t157+t160+t162-t168-t171+t174-t85-t87+t177+t92+t93+t94-t27
     #7+t104+t283+t284
      t509 = t288+t291-t293-t296+t298+t300-t181+t191+t192+t195+t197-t200
     #+t211+t212-t214-t217
      t510 = -t220+t223+t302+t307+t310-t313-t322+t324+t325+t328-t330+t33
     #2+t335+t454+t456-t458
      t513 = sqrt(t506+t507+t509+t510)
      t515 = (t513+0.1E-98)**2
      t519 = -t348+t355+t358-t362-t366+t370+t372-t376-t390+t394-t397-t46
     #8-t470-t474-t478+t482+t485
      t520 = t487-t490-t493-t495+t496+t498+t501-t401-t402+t405+t406+t408
     #+t410-t412+t414-t416+t418+t420
      t524 = -t87+t91+t92+t93+t94+t277+t100+t103+t104-t283+t284
      t526 = -t288+t291+t293-t296+t298+t300-t202-t205-t207-t302+t304
      t527 = -t307-t310+t313-t317+t322+t324+t325+t328-t330+t332+t334+t33
     #5
      t530 = sqrt(t270+t524+t526+t527)
      t532 = (t530+0.1E-98)**2
      t536 = -t348-t351-t355+t358+t362+t366+t370-t372+t376-t381-t386+t39
     #0-t394+t397
      t537 = -t400-t241-t239-t246+t401-t402+t405+t406+t408+t410-t412+t41
     #4-t416-t418-t420
      t541 = -t428+t430+t432+t436-t150-t438-t442+t154+t14+t15-t269-t17+t
     #444-t446+t448
      t542 = -t157+t160+t162+t168+t171-t174-t85-t87+t177+t92+t93+t94+t27
     #7+t104-t283+t284
      t544 = -t288+t291+t293-t296+t298+t300-t181-t191+t192-t195+t197+t20
     #0+t211+t212+t214-t217
      t545 = t220+t223-t302-t307-t310+t313+t322+t324+t325+t328-t330+t332
     #+t335-t454-t456+t458
      t548 = sqrt(t541+t542+t544+t545)
      t550 = (t548+0.1E-98)**2
      t554 = -t348-t355+t358+t362+t366+t370-t372+t376+t390-t394+t397-t46
     #8-t470-t474-t478-t482-t485
      t555 = t487-t490-t493+t495-t496-t498-t501+t401-t402+t405+t406+t408
     #+t410-t412+t414-t416-t418-t420
      t559 = -t428+t430-t432-t436+t150-t438-t442+t154+t14+t15-t269-t17-t
     #444-t446-t448
      t560 = -t157+t160+t162-t168-t171+t174-t85-t87+t177+t92+t93+t94+t27
     #7+t104-t283+t284
      t562 = -t288+t291+t293-t296+t298+t300-t181+t191+t192+t195+t197-t20
     #0+t211+t212-t214-t217
      t563 = -t220+t223-t302-t307-t310+t313+t322+t324+t325+t328-t330+t33
     #2+t335+t454+t456+t458
      t566 = sqrt(t559+t560+t562+t563)
      t568 = (t566+0.1E-98)**2
      t572 = -t348-t355+t358+t362+t366+t370-t372+t376+t390-t394+t397-t46
     #8-t470-t474+t478+t482-t485
      t573 = -t487-t490+t493-t495-t496-t498+t501+t401-t402+t405+t406+t40
     #8+t410-t412+t414-t416-t418-t420
      t578 = -t1*ec/t109/t107*t139-t143/t231/t229*t248-t143/t260/t258*t2
     #64-t143/t341/t339*(t398+2*t421)/2-t426/t464/t462*(t486+2*t502)/2-t
     #426/t515/t513*(t519+2*t520)/2-t143/t532/t530*(t536+2*t537)/2-t426/
     #t550/t548*(t554+2*t555)/2-t426/t568/t566*(t572+2*t573)/2
c
      dVdti4p=t578
c
      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      function dVdsi4p(Xoi,Xoj,Yoi,Yoj,Zoi,Zoj,ppi,si,ti,ppj,sj,tj)
      implicit real*8(a-h,o-z)
      common/cpot/qm,qh,AA,BB,ac,as,bp,ec
c
      t1 = qm**2
      t3 = cos(ppi)
      t4 = cos(si)
      t6 = cos(ti)
      t7 = sin(si)
      t8 = t6*t7
      t9 = sin(ppi)
      t11 = t3*t4-t8*t9
      t12 = t11**2
      t13 = bp**2
      t14 = t12*t13
      t15 = sin(tj)
      t16 = t15**2
      t17 = sin(sj)
      t18 = t17**2
      t19 = t16*t18
      t20 = t19*t13
      t21 = sin(ti)
      t22 = t21**2
      t23 = t7**2
      t24 = t22*t23
      t25 = t24*t13
      t26 = cos(ppj)
      t27 = cos(sj)
      t29 = cos(tj)
      t30 = t29*t17
      t31 = sin(ppj)
      t33 = t26*t27-t30*t31
      t34 = t33*bp
      t36 = 2*t34*Xoj
      t39 = -t31*t27-t30*t26
      t40 = t39*bp
      t42 = 2*t40*Yoj
      t43 = Yoi*t39
      t45 = 2*t43*bp
      t48 = -t9*t4-t8*t3
      t49 = t48*bp
      t51 = 2*t49*Yoj
      t53 = 2*t49*Yoi
      t54 = Xoi*t33
      t56 = 2*t54*bp
      t57 = t11*bp
      t59 = 2*t57*Xoj
      t61 = 2*t57*Xoi
      t62 = t11*t13
      t65 = Zoj**2
      t66 = t15*t17
      t67 = bp*Zoj
      t69 = 2*t66*t67
      t70 = Zoi*t15
      t73 = 2*t70*t17*bp
      t74 = t14+t20+t25+t36+t42-t45-t51+t53-t56-t59+t61-2*t62*t33+t65+t6
     #9-t73
      t75 = t21*t7
      t76 = bp*Zoi
      t78 = 2*t75*t76
      t80 = 2*t75*t67
      t81 = Xoj**2
      t82 = Yoi**2
      t83 = Yoj**2
      t84 = Zoi**2
      t85 = Xoi**2
      t86 = t48**2
      t87 = t86*t13
      t89 = t13*t15*t17
      t92 = t33**2
      t93 = t92*t13
      t95 = 2*Zoi*Zoj
      t96 = t48*t13
      t99 = t39**2
      t100 = t99*t13
      t102 = 2*Xoi*Xoj
      t104 = 2*Yoi*Yoj
      t105 = t78-t80+t81+t82+t83+t84+t85+t87-2*t75*t89+t93-t95-2*t96*t39
     #+t100-t102-t104
      t107 = sqrt(t74+t105)
      t109 = (t107+0.1E-98)**2
      t113 = t21*t4
      t114 = t113*t67
      t115 = t113*t76
      t117 = t6*t4
      t119 = -t3*t7-t117*t9
      t120 = t62*t119
      t123 = t9*t7-t117*t3
      t124 = t96*t123
      t127 = t123*bp
      t128 = t127*Yoi
      t129 = t127*Yoj
      t132 = t119*bp
      t133 = t132*Xoi
      t134 = t132*Xoj
      t135 = t22*t7
      t137 = t135*t13*t4
      t139 = -t114+t115+t120+t124-t123*t13*t39+t128-t129-t119*t13*t33+t1
     #33-t134+t137-t113*t89
      t143 = qm*qh*ec
      t145 = t29*t27
      t147 = -t31*t17+t145*t26
      t150 = 2*Yoi*t147*as
      t151 = t27**2
      t153 = as**2
      t154 = t16*t151*t153
      t157 = t26*t17+t145*t31
      t158 = t157*as
      t160 = 2*t158*Xoj
      t161 = t33*ac
      t163 = 2*t161*Xoj
      t164 = t147*as
      t166 = 2*t164*Yoj
      t167 = t14-t150+t25+t154+t160-t51+t53-t59+t61+t163+t166
      t170 = 2*Xoi*t157*as
      t172 = 2*t54*ac
      t174 = 2*t43*ac
      t175 = -t170+t65-t172-t174+t78-t80+t81+t82+t83+t84+t85
      t177 = t75*bp
      t178 = t66*ac
      t180 = 2*t177*t178
      t185 = 2*t16*t17*ac*t27*as
      t186 = t15*t27
      t187 = t186*as
      t189 = 2*t177*t187
      t191 = 2*t57*t158
      t193 = 2*t161*t158
      t194 = t39*ac
      t196 = 2*t49*t194
      t198 = 2*t57*t161
      t200 = 2*t49*t164
      t202 = 2*t194*t164
      t205 = 2*t70*t17*ac
      t206 = t87-t180-t185+t189-t191+t193-t196-t198-t200+t202-t205
      t207 = as*Zoj
      t209 = 2*t186*t207
      t212 = 2*t70*t27*as
      t213 = ac*Zoj
      t215 = 2*t66*t213
      t216 = t157**2
      t217 = t216*t153
      t218 = ac**2
      t219 = t92*t218
      t220 = t19*t218
      t222 = 2*t194*Yoj
      t223 = t147**2
      t224 = t223*t153
      t225 = t99*t218
      t226 = -t209+t212+t215+t217+t219-t95+t220+t222+t224+t225-t102-t104
      t229 = sqrt(t167+t175+t206+t226)
      t231 = (t229+0.1E-98)**2
      t235 = t132*t158
      t236 = t127*t194
      t237 = t132*t161
      t238 = t127*t164
      t239 = t113*bp
      t240 = t239*t178
      t241 = t239*t187
      t242 = -t114+t115+t120+t124+t128-t129+t133-t134+t137-t235-t236-t23
     #7-t238-t240+t241
      t245 = t14+t150+t25+t154-t160-t51+t53-t59+t61+t163-t166
      t246 = t170+t65-t172-t174+t78-t80+t81+t82+t83+t84+t85
      t248 = t87-t180+t185-t189+t191-t193-t196-t198+t200-t202-t205
      t249 = t209-t212+t215+t217+t219-t95+t220+t222+t224+t225-t102-t104
      t252 = sqrt(t245+t246+t248+t249)
      t254 = (t252+0.1E-98)**2
      t258 = -t114+t115+t120+t124+t128-t129+t133-t134+t137+t235-t236-t23
     #7+t238-t240-t241
      t261 = t11*ac
      t263 = 2*t261*Xoj
      t264 = t48*ac
      t266 = 2*t264*Yoj
      t267 = t4**2
      t268 = t22*t267
      t269 = t268*t153
      t270 = t24*t218
      t271 = -t119*as
      t273 = 2*t271*Xoj
      t274 = t20+t36-t263+t42-t45-t56-t266+t269+t65+t270-t273
      t275 = as*Zoi
      t277 = 2*t113*t275
      t278 = t69-t73+t81+t82+t83+t84+t85-t180-t196-t198-t277
      t280 = -t123*as
      t281 = t264*t280
      t282 = t271*t34
      t283 = ac*Zoi
      t284 = t75*t283
      t285 = t75*t213
      t286 = t280*t40
      t287 = t261*t271
      t288 = t113*t207
      t289 = t271*Xoi
      t290 = t280*Yoi
      t291 = t261*Xoi
      t292 = t280*Yoj
      t293 = t281-t282+t284-t285-t286+t287+t288+t289+t290+t291-t292
      t295 = 2*t264*Yoi
      t296 = t123**2
      t297 = t296*t153
      t301 = 2*t135*ac*t4*as
      t302 = t119**2
      t303 = t302*t153
      t304 = t12*t218
      t305 = t86*t218
      t306 = t113*as
      t307 = t66*bp
      t309 = 2*t306*t307
      t310 = t295+t93-t95+t297-t301+t303+t304+t305+t100-t102+t309-t104
      t313 = sqrt(t274+t278+2*t293+t310)
      t315 = (t313+0.1E-98)**2
      t319 = -t123*t153
      t320 = t319*t48
      t321 = t11*as
      t322 = t321*Xoj
      t323 = t48*as
      t324 = t323*Yoj
      t326 = t12*ac*as
      t327 = t323*Yoi
      t328 = -t119*t153
      t329 = t328*t11
      t331 = t86*ac*as
      t332 = t123*ac
      t333 = t332*Yoj
      t334 = t11*t218
      t335 = t334*t119
      t338 = t22*t4*t153*t7
      t340 = t135*t218*t4
      t341 = t75*as
      t342 = t341*t307
      t343 = t320-t322-t324+t326-t236-t237+t327+t329+t331-t333-t240+t335
     #-t338+t340-t342
      t344 = ac*as
      t345 = t268*t344
      t346 = t24*t344
      t347 = t75*t275
      t348 = t332*t280
      t349 = t75*t207
      t350 = t113*t213
      t351 = t321*t34
      t352 = t113*t283
      t353 = t119*ac
      t354 = t353*t271
      t355 = t323*t40
      t356 = t332*Yoi
      t357 = t353*Xoj
      t358 = t48*t218
      t359 = t358*t123
      t360 = t353*Xoi
      t361 = t321*Xoi
      t362 = -t345+t346+t347+t348-t349-t350-t351+t352+t354-t355+t356-t35
     #7+t359+t360+t361
      t366 = qh**2
      t367 = t366*ec
      t369 = 2*t358*t39
      t371 = 2*t334*t33
      t373 = 2*t328*t157
      t374 = -t369-t150+t154+t160-t263-t371-t266+t269+t163+t166-t373-t17
     #0+t65+t270-t273
      t375 = -t172-t174+t81+t82+t83+t84+t85-t185+t193+t202-t205-t209+t21
     #2+t215+t217-t277
      t377 = 2*t281
      t378 = 2*t284
      t379 = 2*t285
      t380 = 2*t287
      t381 = 2*t288
      t382 = 2*t289
      t383 = 2*t290
      t384 = 2*t291
      t385 = 2*t292
      t386 = t377+t378-t379+t380+t381+t382+t383+t384-t385+t295+t219-t95+
     #t220+t222+t297-t301
      t388 = 2*t264*t164
      t390 = 2*t271*t161
      t392 = 2*t261*t158
      t394 = 2*t280*t194
      t396 = 2*t319*t147
      t398 = t218*t15*t17
      t400 = 2*t75*t398
      t402 = 2*t306*t178
      t405 = 2*t75*ac*t187
      t407 = t153*t15*t27
      t409 = 2*t113*t407
      t410 = t303+t304+t305+t224+t225-t388-t390-t392-t394-t102-t396-t400
     #-t104+t402+t405-t409
      t413 = sqrt(t374+t375+t386+t410)
      t415 = (t413+0.1E-98)**2
      t420 = t123*t218*t39
      t422 = t48*t153*t147
      t424 = t11*t153*t157
      t426 = t119*t218*t33
      t427 = t320-t322-t420-t422-t324-t424+t326-t426+t327+t329+t331-t333
     #+t335-t338+t340-t345+t346+t347
      t428 = t341*t178
      t429 = t75*t407
      t431 = t113*ac*t187
      t432 = t353*t158
      t433 = t321*t161
      t434 = t323*t194
      t435 = t332*t164
      t436 = t113*t398
      t437 = t348-t349-t350+t352+t354+t356-t357+t359+t360+t361-t428+t429
     #+t431-t432-t433-t434-t435-t436
      t441 = -t369+t150+t154-t160-t263-t371-t266+t269+t163-t166+t373+t17
     #0+t65+t270-t273
      t442 = -t172-t174+t81+t82+t83+t84+t85+t185-t193-t202-t205+t209-t21
     #2+t215+t217-t277
      t444 = t303+t304+t305+t224+t225+t388-t390+t392-t394-t102+t396-t400
     #-t104+t402-t405+t409
      t447 = sqrt(t441+t442+t386+t444)
      t449 = (t447+0.1E-98)**2
      t453 = t320-t322-t420+t422-t324+t424+t326-t426+t327+t329+t331-t333
     #+t335-t338+t340-t345+t346+t347
      t454 = t348-t349-t350+t352+t354+t356-t357+t359+t360+t361-t428-t429
     #-t431+t432-t433-t434+t435-t436
      t458 = t20+t36-t263+t42-t45-t56-t266+t269+t65+t270+t273
      t459 = t69-t73+t81+t82+t83+t84+t85-t180-t196-t198+t277
      t461 = -t281+t282+t284-t285+t286-t287-t288-t289-t290+t291+t292
      t462 = t295+t93-t95+t297+t301+t303+t304+t305+t100-t102-t309-t104
      t465 = sqrt(t458+t459+2*t461+t462)
      t467 = (t465+0.1E-98)**2
      t471 = t320+t322+t324-t326-t236-t237-t327+t329-t331-t333-t240+t335
     #-t338+t340+t342
      t472 = t345-t346-t347-t348+t349-t350+t351+t352-t354+t355+t356-t357
     #+t359+t360-t361
      t476 = -t369-t150+t154+t160-t263-t371-t266+t269+t163+t166+t373-t17
     #0+t65+t270+t273
      t477 = -t172-t174+t81+t82+t83+t84+t85-t185+t193+t202-t205-t209+t21
     #2+t215+t217+t277
      t479 = -t377+t378-t379-t380-t381-t382-t383+t384+t385+t295+t219-t95
     #+t220+t222+t297+t301
      t480 = t303+t304+t305+t224+t225-t388+t390-t392+t394-t102+t396-t400
     #-t104-t402+t405+t409
      t483 = sqrt(t476+t477+t479+t480)
      t485 = (t483+0.1E-98)**2
      t489 = t320+t322-t420+t422+t324+t424-t326-t426-t327+t329-t331-t333
     #+t335-t338+t340+t345-t346-t347
      t490 = -t348+t349-t350+t352-t354+t356-t357+t359+t360-t361+t428-t42
     #9+t431-t432+t433+t434-t435-t436
      t494 = -t369+t150+t154-t160-t263-t371-t266+t269+t163-t166-t373+t17
     #0+t65+t270+t273
      t495 = -t172-t174+t81+t82+t83+t84+t85+t185-t193-t202-t205+t209-t21
     #2+t215+t217+t277
      t497 = t303+t304+t305+t224+t225+t388+t390+t392+t394-t102-t396-t400
     #-t104-t402-t405-t409
      t500 = sqrt(t494+t495+t479+t497)
      t502 = (t500+0.1E-98)**2
      t506 = t320+t322-t420-t422+t324-t424-t326-t426-t327+t329-t331-t333
     #+t335-t338+t340+t345-t346-t347
      t507 = -t348+t349-t350+t352-t354+t356-t357+t359+t360-t361+t428+t42
     #9-t431+t432+t433+t434+t435-t436
      t512 = -t1*ec/t109/t107*t139-t143/t231/t229*t242-t143/t254/t252*t2
     #58-t143/t315/t313*(2*t343+2*t362)/2-t367/t415/t413*(2*t427+2*t437)
     #/2-t367/t449/t447*(2*t453+2*t454)/2-t143/t467/t465*(2*t471+2*t472)
     #/2-t367/t485/t483*(2*t489+2*t490)/2-t367/t502/t500*(2*t506+2*t507)
     #/2
c
      dVdsi4p=t512
c
      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      function dVdppi4p(Xoi,Xoj,Yoi,Yoj,Zoi,Zoj,ppi,si,ti,ppj,sj,tj)
      implicit real*8(a-h,o-z)
      common/cpot/qm,qh,AA,BB,ac,as,bp,ec
c
      t1 = qm**2
      t3 = sin(ppi)
      t4 = cos(si)
      t6 = cos(ti)
      t7 = sin(si)
      t8 = t6*t7
      t9 = cos(ppi)
      t11 = -t3*t4-t8*t9
      t12 = t11*bp
      t14 = 2*t12*Yoj
      t15 = cos(ppj)
      t16 = cos(sj)
      t18 = cos(tj)
      t19 = sin(sj)
      t20 = t18*t19
      t21 = sin(ppj)
      t23 = t15*t16-t20*t21
      t24 = t23*bp
      t26 = 2*t24*Xoj
      t29 = -t21*t16-t20*t15
      t30 = t29*bp
      t32 = 2*t30*Yoj
      t33 = Yoi*t29
      t35 = 2*t33*bp
      t36 = sin(tj)
      t37 = t36**2
      t38 = t19**2
      t39 = t37*t38
      t40 = bp**2
      t41 = t39*t40
      t44 = t9*t4-t8*t3
      t45 = t44*bp
      t47 = 2*t45*Xoj
      t49 = 2*t45*Xoi
      t50 = t44*t40
      t53 = sin(ti)
      t54 = t53**2
      t55 = t7**2
      t56 = t54*t55
      t57 = t56*t40
      t58 = t11*t40
      t62 = 2*t12*Yoi
      t64 = 2*Zoi*Zoj
      t65 = t23**2
      t66 = t65*t40
      t68 = 2*Yoi*Yoj
      t69 = t53*t7
      t70 = bp*Zoj
      t72 = 2*t69*t70
      t73 = -t14+t26+t32-t35+t41-t47+t49-2*t50*t23+t57-2*t58*t29+t62-t64
     #+t66-t68-t72
      t74 = t44**2
      t75 = t74*t40
      t77 = 2*Xoi*Xoj
      t78 = Xoi*t23
      t80 = 2*t78*bp
      t83 = 2*t69*bp*Zoi
      t84 = t11**2
      t85 = t84*t40
      t86 = t36*t19
      t88 = 2*t86*t70
      t89 = Zoi*t36
      t92 = 2*t89*t19*bp
      t97 = t29**2
      t98 = t97*t40
      t99 = Zoj**2
      t100 = Xoi**2
      t101 = Zoi**2
      t102 = Yoj**2
      t103 = Yoi**2
      t104 = Xoj**2
      t105 = t75-t77-t80+t83+t85+t88-t92-2*t69*t40*t36*t19+t98+t99+t100+
     #t101+t102+t103+t104
      t107 = sqrt(t73+t105)
      t109 = (t107+0.1E-98)**2
      t113 = t50*t11
      t114 = -t58*t44
      t117 = -t44*bp
      t118 = t117*Yoi
      t119 = t117*Yoj
      t121 = t12*Xoi
      t122 = t12*Xoj
      t127 = qm*qh*ec
      t129 = t18*t16
      t131 = t15*t19+t129*t21
      t134 = 2*Xoi*t131*as
      t136 = 2*t33*ac
      t137 = t23*ac
      t139 = 2*t137*Xoj
      t140 = t16**2
      t142 = as**2
      t143 = t37*t140*t142
      t144 = ac**2
      t145 = t39*t144
      t146 = -t14-t134-t47+t49+t57-t136+t139+t62+t143-t64+t145
      t149 = -t21*t19+t129*t15
      t152 = 2*Yoi*t149*as
      t153 = t131*as
      t155 = 2*t153*Xoj
      t156 = t29*ac
      t158 = 2*t156*Yoj
      t159 = t149*as
      t161 = 2*t159*Yoj
      t163 = 2*t78*ac
      t164 = -t68-t72-t152+t155+t158+t161+t75-t163-t77+t83+t85
      t167 = 2*t45*t137
      t169 = 2*t12*t156
      t171 = 2*t137*t153
      t172 = t97*t144
      t175 = 2*t89*t19*ac
      t176 = t99+t100+t101+t102+t103+t104-t167-t169+t171+t172-t175
      t178 = 2*t156*t159
      t180 = 2*t12*t159
      t181 = t131**2
      t182 = t181*t142
      t185 = 2*t89*t16*as
      t186 = ac*Zoj
      t188 = 2*t86*t186
      t189 = t36*t16
      t190 = as*Zoj
      t192 = 2*t189*t190
      t193 = t65*t144
      t195 = 2*t45*t153
      t196 = t149**2
      t197 = t196*t142
      t198 = t69*bp
      t199 = t86*ac
      t201 = 2*t198*t199
      t202 = t189*as
      t204 = 2*t198*t202
      t209 = 2*t37*t19*ac*t16*as
      t210 = t178-t180+t182+t185+t188-t192+t193-t195+t197-t201+t204-t209
      t213 = sqrt(t146+t164+t176+t210)
      t215 = (t213+0.1E-98)**2
      t219 = t12*t153
      t220 = t117*t156
      t221 = t12*t137
      t222 = t117*t159
      t226 = -t14+t134-t47+t49+t57-t136+t139+t62+t143-t64+t145
      t227 = -t68-t72+t152-t155+t158-t161+t75-t163-t77+t83+t85
      t229 = t99+t100+t101+t102+t103+t104-t167-t169-t171+t172-t175
      t230 = -t178+t180+t182-t185+t188+t192+t193+t195+t197-t201-t204+t20
     #9
      t233 = sqrt(t226+t227+t229+t230)
      t235 = (t233+0.1E-98)**2
      t243 = t6*t4
      t245 = t9*t7+t243*t3
      t246 = t245*as
      t248 = 2*t246*Xoj
      t249 = t44*ac
      t251 = 2*t249*Xoj
      t252 = t56*t144
      t253 = t4**2
      t255 = t54*t253*t142
      t258 = -t3*t7+t243*t9
      t259 = t258*as
      t261 = 2*t259*Yoi
      t262 = -t248-t251+t252+t26+t255+t32-t35+t261+t41-t64+t66
      t263 = t11*ac
      t265 = 2*t263*Yoi
      t267 = 2*t263*Yoj
      t269 = 2*t259*Yoj
      t271 = 2*t249*Xoi
      t273 = 2*t246*Xoi
      t274 = -t68+t265-t77-t80+t88-t267-t269+t271+t273-t92+t98
      t277 = 2*t263*t259
      t278 = t53*t4
      t281 = 2*t278*as*Zoi
      t282 = t99+t100+t101+t102+t103+t104-t167-t169-t201+t277-t281
      t284 = 2*t278*t190
      t286 = 2*t246*t24
      t288 = 2*t249*t246
      t290 = 2*t259*t30
      t292 = 2*t69*t186
      t295 = 2*t69*ac*Zoi
      t296 = t74*t144
      t301 = 2*t54*t7*ac*t4*as
      t302 = t258**2
      t303 = t302*t142
      t304 = t278*as
      t307 = 2*t304*t86*bp
      t308 = t84*t144
      t309 = t245**2
      t310 = t309*t142
      t311 = t284-t286+t288-t290-t292+t295+t296-t301+t303+t307+t308+t310
      t314 = sqrt(t262+t274+t282+t311)
      t316 = (t314+0.1E-98)**2
      t320 = t44*t144
      t321 = t320*t11
      t322 = t245*t142
      t323 = t322*t258
      t324 = t11*t144
      t325 = -t324*t44
      t326 = t258*t142
      t327 = -t326*t245
      t328 = t263*Xoj
      t329 = t259*Xoj
      t330 = t259*Xoi
      t331 = -t245*as
      t332 = t331*Yoi
      t333 = t263*Xoi
      t334 = t331*Yoj
      t335 = -t44*ac
      t336 = t335*Yoi
      t337 = t335*Yoj
      t338 = t335*t259
      t339 = t263*t331
      t340 = t331*t30
      t341 = t263*t246
      t342 = t249*t259
      t343 = t259*t24
      t344 = t321+t323+t325+t327-t328-t329+t330+t332+t333-t334+t336-t337
     #+t338+t339-t340+t341+t342-t343-t220-t221
      t347 = qh**2
      t348 = t347*ec
      t350 = 2*t263*t159
      t353 = 2*t69*ac*t202
      t357 = 2*t278*t142*t36*t16
      t359 = 2*t320*t23
      t363 = 2*t69*t144*t36*t19
      t365 = 2*t304*t199
      t367 = 2*t322*t131
      t368 = -t248-t251-t350+t353-t357-t359-t363+t365+t252-t134+t255-t36
     #7+t261-t136+t139
      t370 = 2*t324*t29
      t372 = 2*t326*t149
      t373 = t143-t64+t145-t68-t152+t155-t370+t158+t161+t265-t372-t163-t
     #77-t267-t269+t271
      t375 = t273+t99+t100+t101+t102+t103+t104+t171+t172-t175+t178+t182+
     #t185+t188-t192+t193
      t377 = 2*t246*t137
      t379 = 2*t249*t153
      t381 = 2*t259*t156
      t382 = t197-t209+t277-t281+t284+t288-t292+t295+t296-t301+t303-t377
     #-t379-t381+t308+t310
      t385 = sqrt(t368+t373+t375+t382)
      t387 = (t385+0.1E-98)**2
      t391 = t321+t323+t325+t327-t328-t329+t330+t332+t333-t334+t336-t337
      t392 = t335*t159
      t393 = t331*t156
      t394 = t324*t23
      t395 = t326*t131
      t397 = -t44*t144*t29
      t399 = -t245*t142*t149
      t400 = t263*t153
      t401 = t259*t137
      t402 = t338+t339+t341+t342-t392-t393-t394-t395-t397-t399-t400-t401
      t406 = -t248-t251+t350-t353+t357-t359-t363+t365+t252+t134+t255+t36
     #7+t261-t136+t139
      t407 = t143-t64+t145-t68+t152-t155-t370+t158-t161+t265+t372-t163-t
     #77-t267-t269+t271
      t409 = t273+t99+t100+t101+t102+t103+t104-t171+t172-t175-t178+t182-
     #t185+t188+t192+t193
      t410 = t197+t209+t277-t281+t284+t288-t292+t295+t296-t301+t303-t377
     #+t379-t381+t308+t310
      t413 = sqrt(t406+t407+t409+t410)
      t415 = (t413+0.1E-98)**2
      t419 = t338+t339+t341+t342+t392-t393-t394+t395-t397+t399+t400-t401
      t423 = t248-t251+t252+t26+t255+t32-t35-t261+t41-t64+t66
      t424 = -t68+t265-t77-t80+t88-t267+t269+t271-t273-t92+t98
      t426 = t99+t100+t101+t102+t103+t104-t167-t169-t201-t277+t281
      t427 = -t284+t286-t288+t290-t292+t295+t296+t301+t303-t307+t308+t31
     #0
      t430 = sqrt(t423+t424+t426+t427)
      t432 = (t430+0.1E-98)**2
      t436 = t321+t323+t325+t327-t328+t329-t330-t332+t333+t334+t336-t337
     #-t338-t339+t340-t341-t342+t343-t220-t221
      t439 = t248-t251-t350+t353+t357-t359-t363-t365+t252-t134+t255+t367
     #-t261-t136+t139
      t440 = t143-t64+t145-t68-t152+t155-t370+t158+t161+t265+t372-t163-t
     #77-t267+t269+t271
      t442 = -t273+t99+t100+t101+t102+t103+t104+t171+t172-t175+t178+t182
     #+t185+t188-t192+t193
      t443 = t197-t209-t277+t281-t284-t288-t292+t295+t296+t301+t303+t377
     #-t379+t381+t308+t310
      t446 = sqrt(t439+t440+t442+t443)
      t448 = (t446+0.1E-98)**2
      t452 = t321+t323+t325+t327-t328+t329-t330-t332+t333+t334+t336-t337
      t453 = -t338-t339-t341-t342-t392+t393-t394+t395-t397+t399-t400+t40
     #1
      t457 = t248-t251+t350-t353-t357-t359-t363-t365+t252+t134+t255-t367
     #-t261-t136+t139
      t458 = t143-t64+t145-t68+t152-t155-t370+t158-t161+t265-t372-t163-t
     #77-t267+t269+t271
      t460 = -t273+t99+t100+t101+t102+t103+t104-t171+t172-t175-t178+t182
     #-t185+t188+t192+t193
      t461 = t197+t209-t277+t281-t284-t288-t292+t295+t296+t301+t303+t377
     #+t379+t381+t308+t310
      t464 = sqrt(t457+t458+t460+t461)
      t466 = (t464+0.1E-98)**2
      t470 = -t338-t339-t341-t342+t392+t393-t394-t395-t397-t399+t400+t40
     #1
      t475 = -t1*ec/t109/t107*(t113+t114+t44*t40*t29+t118-t119-t58*t23+t
     #121-t122)-t127/t215/t213*(t113+t114+t118-t119+t121-t122-t219-t220-
     #t221-t222)-t127/t235/t233*(t113+t114+t118-t119+t121-t122+t219-t220
     #-t221+t222)-t127/t316/t314*t344-t348/t387/t385*(2*t391+2*t402)/2-t
     #348/t415/t413*(2*t391+2*t419)/2-t127/t432/t430*t436-t348/t448/t446
     #*(2*t452+2*t453)/2-t348/t466/t464*(2*t452+2*t470)/2
c
      dVdppi4p=t475 
c
      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      function dVdXoj4p(Xoi,Xoj,Yoi,Yoj,Zoi,Zoj,ppi,si,ti,ppj,sj,tj)
      implicit real*8(a-h,o-z)
      common/cpot/qm,qh,AA,BB,ac,as,bp,ec
c
      t1 = qm**2
      t3 = cos(ppi)
      t4 = cos(si)
      t6 = cos(ti)
      t7 = sin(si)
      t8 = t6*t7
      t9 = sin(ppi)
      t11 = t3*t4-t8*t9
      t12 = t11**2
      t13 = bp**2
      t14 = t12*t13
      t15 = sin(ti)
      t16 = t15**2
      t17 = t7**2
      t18 = t16*t17
      t19 = t18*t13
      t21 = 2*Yoi*Yoj
      t23 = 2*Xoi*Xoj
      t25 = 2*Zoi*Zoj
      t28 = -t9*t4-t8*t3
      t29 = t28**2
      t30 = t29*t13
      t31 = sin(ppj)
      t32 = cos(sj)
      t34 = cos(tj)
      t35 = sin(sj)
      t36 = t34*t35
      t37 = cos(ppj)
      t39 = -t31*t32-t36*t37
      t40 = t39**2
      t41 = t40*t13
      t42 = Xoj**2
      t43 = t28*bp
      t45 = 2*t43*Yoj
      t46 = sin(tj)
      t47 = t46**2
      t48 = t35**2
      t49 = t47*t48
      t50 = t49*t13
      t51 = t39*bp
      t53 = 2*t51*Yoj
      t54 = Yoi*t39
      t56 = 2*t54*bp
      t58 = 2*t43*Yoi
      t62 = t11*bp
      t64 = 2*t62*Xoj
      t65 = t14+t19-t21-t23-t25+t30+t41+t42-t45+t50+t53-t56+t58-2*t28*t1
     #3*t39-t64
      t67 = 2*t62*Xoi
      t71 = t37*t32-t36*t31
      t74 = Xoi*t71
      t76 = 2*t74*bp
      t77 = t46*t35
      t78 = bp*Zoj
      t80 = 2*t77*t78
      t81 = Zoi*t46
      t84 = 2*t81*t35*bp
      t85 = t15*t7
      t88 = 2*t85*bp*Zoi
      t93 = Zoj**2
      t94 = Zoi**2
      t95 = Yoj**2
      t96 = Xoi**2
      t97 = Yoi**2
      t98 = t71*bp
      t100 = 2*t98*Xoj
      t102 = 2*t85*t78
      t103 = t71**2
      t104 = t103*t13
      t105 = t67-2*t11*t13*t71-t76+t80-t84+t88-2*t85*t13*t46*t35+t93+t94
     #+t95+t96+t97+t100-t102+t104
      t107 = sqrt(t65+t105)
      t109 = (t107+0.1E-98)**2
      t118 = qm*qh*ec
      t119 = ac**2
      t120 = t49*t119
      t121 = t46*t32
      t122 = as*Zoj
      t124 = 2*t121*t122
      t127 = 2*t81*t35*ac
      t128 = t39*ac
      t130 = t34*t32
      t132 = -t35*t31+t130*t37
      t133 = t132*as
      t135 = 2*t128*t133
      t136 = t14+t19+t120-t21-t23-t25+t30-t124-t127+t135+t42
      t139 = 2*Yoi*t132*as
      t141 = 2*t74*ac
      t143 = 2*t54*ac
      t146 = t37*t35+t130*t31
      t149 = 2*Xoi*t146*as
      t150 = t71*ac
      t152 = 2*t150*Xoj
      t153 = t32**2
      t155 = as**2
      t156 = t47*t153*t155
      t157 = t132**2
      t158 = t157*t155
      t159 = t40*t119
      t161 = 2*t133*Yoj
      t162 = -t139-t141-t143-t149+t152+t156+t158+t159+t161-t45+t58
      t164 = t103*t119
      t165 = t146**2
      t166 = t165*t155
      t167 = -t64+t67+t164+t88+t93+t94+t95+t96+t97+t166-t102
      t168 = t146*as
      t169 = t168*Xoj
      t170 = ac*Zoj
      t171 = t77*t170
      t173 = t81*t32*as
      t174 = t43*t133
      t175 = t62*t150
      t176 = t43*t128
      t177 = t150*t168
      t178 = t62*t168
      t179 = t128*Yoj
      t180 = t85*bp
      t181 = t121*as
      t182 = t180*t181
      t183 = t77*ac
      t184 = t180*t183
      t188 = t47*t35*ac*t32*as
      t189 = t169+t171+t173-t174-t175-t176+t177-t178+t179+t182-t184-t188
      t192 = sqrt(t136+t162+t167+2*t189)
      t194 = (t192+0.1E-98)**2
      t202 = t14+t19+t120-t21-t23-t25+t30+t124-t127-t135+t42
      t203 = t139-t141-t143+t149+t152+t156+t158+t159-t161-t45+t58
      t205 = -t169+t171-t173+t174-t175-t176-t177+t178+t179-t182-t184+t18
     #8
      t208 = sqrt(t202+t203+t167+2*t205)
      t210 = (t208+0.1E-98)**2
      t218 = t18*t119
      t219 = t11*ac
      t221 = 2*t219*Xoj
      t222 = t218-t221-t21-t23-t25+t41+t42+t50+t53-t56-t76
      t224 = t6*t4
      t226 = t3*t7+t224*t9
      t227 = t226*as
      t229 = 2*t227*Xoi
      t232 = -t9*t7+t224*t3
      t233 = t232*as
      t235 = 2*t233*Yoi
      t236 = t80-t84+t93+t94+t95+t96+t97+t100+t104+t229+t235
      t239 = 2*t219*Xoi
      t241 = 2*t233*Yoj
      t242 = t12*t119
      t243 = t232**2
      t244 = t243*t155
      t245 = 2*t175
      t246 = 2*t176
      t247 = 2*t184
      t248 = t28*ac
      t250 = 2*t248*t233
      t251 = t15*t4
      t254 = 2*t251*as*Zoi
      t256 = 2*t251*t122
      t259 = 2*t85*ac*Zoi
      t260 = t239-t241+t242+t244-t245-t246-t247+t250-t254+t256+t259
      t262 = 2*t227*t98
      t264 = 2*t219*t227
      t266 = 2*t233*t51
      t268 = 2*t248*Yoi
      t270 = 2*t248*Yoj
      t271 = t4**2
      t273 = t16*t271*t155
      t275 = 2*t227*Xoj
      t276 = t251*as
      t279 = 2*t276*t77*bp
      t284 = 2*t16*t7*ac*t4*as
      t286 = 2*t85*t170
      t287 = t226**2
      t288 = t287*t155
      t289 = t29*t119
      t290 = -t262+t264-t266+t268-t270+t273-t275+t279-t284-t286+t288+t28
     #9
      t293 = sqrt(t222+t236+t260+t290)
      t295 = (t293+0.1E-98)**2
      t303 = qh**2
      t304 = t303*ec
      t307 = 2*t226*t155*t146
      t308 = t218-t307-t221+t120-t21-t23-t25-t124-t127+t135+t42-t139-t14
     #1-t143-t149
      t309 = 2*t169
      t310 = t152+t156+t158+t159+t161+t164+t93+t94+t95+t96+t97+t166+t309
     #+t229+t235+t239
      t312 = 2*t171
      t313 = 2*t173
      t314 = 2*t177
      t315 = 2*t179
      t316 = 2*t188
      t319 = 2*t85*ac*t181
      t323 = 2*t251*t155*t46*t32
      t325 = 2*t248*t133
      t329 = 2*t85*t119*t46*t35
      t331 = 2*t276*t183
      t332 = -t241+t242+t244+t312+t313+t314+t315-t316+t319-t323-t325-t32
     #9+t331+t250-t254+t256
      t334 = 2*t227*t150
      t336 = 2*t233*t128
      t338 = 2*t219*t168
      t341 = 2*t11*t119*t71
      t344 = 2*t28*t119*t39
      t347 = 2*t232*t155*t132
      t348 = t259+t264-t334-t336-t338+t268-t270+t273-t275-t341-t344-t347
     #-t284-t286+t288+t289
      t351 = sqrt(t308+t310+t332+t348)
      t353 = (t351+0.1E-98)**2
      t361 = t218+t307-t221+t120-t21-t23-t25+t124-t127-t135+t42+t139-t14
     #1-t143+t149
      t362 = t152+t156+t158+t159-t161+t164+t93+t94+t95+t96+t97+t166-t309
     #+t229+t235+t239
      t364 = -t241+t242+t244+t312-t313-t314+t315+t316-t319+t323+t325-t32
     #9+t331+t250-t254+t256
      t365 = t259+t264-t334-t336+t338+t268-t270+t273-t275-t341-t344+t347
     #-t284-t286+t288+t289
      t368 = sqrt(t361+t362+t364+t365)
      t370 = (t368+0.1E-98)**2
      t378 = t80-t84+t93+t94+t95+t96+t97+t100+t104-t229-t235
      t380 = t239+t241+t242+t244-t245-t246-t247-t250+t254-t256+t259
      t381 = t262-t264+t266+t268-t270+t273+t275-t279+t284-t286+t288+t289
      t384 = sqrt(t222+t378+t380+t381)
      t386 = (t384+0.1E-98)**2
      t394 = t218+t307-t221+t120-t21-t23-t25-t124-t127+t135+t42-t139-t14
     #1-t143-t149
      t395 = t152+t156+t158+t159+t161+t164+t93+t94+t95+t96+t97+t166+t309
     #-t229-t235+t239
      t397 = t241+t242+t244+t312+t313+t314+t315-t316+t319+t323-t325-t329
     #-t331-t250+t254-t256
      t398 = t259-t264+t334+t336-t338+t268-t270+t273+t275-t341-t344+t347
     #+t284-t286+t288+t289
      t401 = sqrt(t394+t395+t397+t398)
      t403 = (t401+0.1E-98)**2
      t411 = t218-t307-t221+t120-t21-t23-t25+t124-t127-t135+t42+t139-t14
     #1-t143+t149
      t412 = t152+t156+t158+t159-t161+t164+t93+t94+t95+t96+t97+t166-t309
     #-t229-t235+t239
      t414 = t241+t242+t244+t312-t313-t314+t315+t316-t319-t323+t325-t329
     #-t331-t250+t254-t256
      t415 = t259-t264+t334+t336+t338+t268-t270+t273+t275-t341-t344-t347
     #+t284-t286+t288+t289
      t418 = sqrt(t411+t412+t414+t415)
      t420 = (t418+0.1E-98)**2
      t429 = sqrt(t96-t23+t42+t97-t21+t95+t94-t25+t93)
      t430 = t429+0.1E-98
      t431 = t430**2
      t432 = t431**2
      t434 = t432**2
      t440 = 2/t429*(-Xoi+Xoj)
      t449 = -t1*ec/t109/t107*(-Xoi+Xoj+t98-t62)-t118/t194/t192*(-Xoi+Xo
     #j-t62+t150+t168)-t118/t210/t208*(-Xoi+Xoj-t62+t150-t168)-t118/t295
     #/t293*(-Xoi+Xoj+t98-t219-t227)-t304/t353/t351*(-Xoi+Xoj-t219-t227+
     #t150+t168)-t304/t370/t368*(-Xoi+Xoj-t219-t227+t150-t168)-t118/t386
     #/t384*(-Xoi+Xoj+t98-t219+t227)-t304/t403/t401*(-Xoi+Xoj-t219+t227+
     #t150+t168)-t304/t420/t418*(-Xoi+Xoj-t219+t227+t150-t168)-6*AA/t434
     #/t432/t430*t440+3*BB/t432/t431/t430*t440
c
      dVdXoj4p=t449
c
      return
      end
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      function dVdYoj4p(Xoi,Xoj,Yoi,Yoj,Zoi,Zoj,ppi,si,ti,ppj,sj,tj)
      implicit real*8(a-h,o-z)
      common/cpot/qm,qh,AA,BB,ac,as,bp,ec
c
      t1 = qm**2
      t3 = sin(ppi)
      t4 = cos(si)
      t6 = cos(ti)
      t7 = sin(si)
      t8 = t6*t7
      t9 = cos(ppi)
      t11 = -t3*t4-t8*t9
      t12 = t11*bp
      t14 = 2*t12*Yoi
      t15 = sin(tj)
      t16 = t15**2
      t17 = sin(sj)
      t18 = t17**2
      t19 = t16*t18
      t20 = bp**2
      t21 = t19*t20
      t24 = t9*t4-t8*t3
      t25 = t24*bp
      t27 = 2*t25*Xoj
      t29 = 2*t25*Xoi
      t31 = cos(ppj)
      t32 = cos(sj)
      t34 = cos(tj)
      t35 = t34*t17
      t36 = sin(ppj)
      t38 = t31*t32-t35*t36
      t41 = Xoi*t38
      t43 = 2*t41*bp
      t45 = 2*t12*Yoj
      t49 = -t36*t32-t35*t31
      t52 = sin(ti)
      t53 = t52**2
      t54 = t7**2
      t55 = t53*t54
      t56 = t55*t20
      t57 = t49*bp
      t59 = 2*t57*Yoj
      t60 = Yoi*t49
      t62 = 2*t60*bp
      t63 = t38*bp
      t65 = 2*t63*Xoj
      t66 = t15*t17
      t67 = bp*Zoj
      t69 = 2*t66*t67
      t70 = t52*t7
      t71 = t70*t67
      t74 = t70*bp*Zoi
      t76 = t14+t21-t27+t29-2*t24*t20*t38-t43-t45-2*t11*t20*t49+t56+t59-
     #t62+t65+t69-2*t71+2*t74
      t77 = Zoi*t15
      t80 = 2*t77*t17*bp
      t82 = 2*Yoi*Yoj
      t83 = t38**2
      t84 = t83*t20
      t85 = Yoj**2
      t86 = Yoi**2
      t87 = Xoj**2
      t88 = Xoi**2
      t89 = Zoj**2
      t94 = t24**2
      t95 = t94*t20
      t97 = 2*Zoi*Zoj
      t98 = Zoi**2
      t99 = t49**2
      t100 = t99*t20
      t101 = t11**2
      t102 = t101*t20
      t104 = 2*Xoi*Xoj
      t105 = -t80-t82+t84+t85+t86+t87+t88+t89-2*t70*t20*t15*t17+t95-t97+
     #t98+t100+t102-t104
      t107 = sqrt(t76+t105)
      t109 = (t107+0.1E-98)**2
      t118 = qm*qh*ec
      t119 = t38*ac
      t121 = 2*t25*t119
      t122 = ac**2
      t123 = t19*t122
      t125 = t34*t32
      t127 = -t36*t17+t125*t31
      t128 = t127*as
      t130 = 2*t128*Yoj
      t131 = t32**2
      t133 = as**2
      t134 = t16*t131*t133
      t136 = 2*t60*ac
      t138 = 2*t41*ac
      t139 = -t121+t14+t123-t27+t29-t45+t56+t130+t134-t136-t138
      t141 = Yoi*t127*as
      t144 = t31*t17+t125*t36
      t145 = t144*as
      t146 = t145*Xoj
      t147 = t119*Xoj
      t149 = Xoi*t144*as
      t150 = t49*ac
      t151 = t150*Yoj
      t152 = ac*Zoj
      t153 = t66*t152
      t155 = t77*t17*ac
      t156 = t15*t32
      t157 = as*Zoj
      t158 = t156*t157
      t160 = t77*t32*as
      t161 = -t141+t146+t147-t149+t151-t71+t74+t153-t155-t158+t160
      t164 = 2*t25*t145
      t166 = 2*t119*t145
      t168 = 2*t12*t150
      t170 = 2*t12*t128
      t172 = 2*t150*t128
      t173 = -t164+t166-t168-t170+t172-t82+t85+t86+t87+t88+t89
      t174 = t83*t122
      t179 = 2*t16*t17*ac*t32*as
      t180 = t70*bp
      t181 = t156*as
      t183 = 2*t180*t181
      t184 = t66*ac
      t186 = 2*t180*t184
      t187 = t144**2
      t188 = t187*t133
      t189 = t127**2
      t190 = t189*t133
      t191 = t99*t122
      t192 = t95-t97+t98+t102-t104+t174-t179+t183-t186+t188+t190+t191
      t195 = sqrt(t139+2*t161+t173+t192)
      t197 = (t195+0.1E-98)**2
      t205 = -t121+t14+t123-t27+t29-t45+t56-t130+t134-t136-t138
      t206 = t141-t146+t147+t149+t151-t71+t74+t153-t155+t158-t160
      t208 = t164-t166-t168+t170-t172-t82+t85+t86+t87+t88+t89
      t209 = t95-t97+t98+t102-t104+t174+t179-t183-t186+t188+t190+t191
      t212 = sqrt(t205+2*t206+t208+t209)
      t214 = (t212+0.1E-98)**2
      t222 = t4**2
      t224 = t53*t222*t133
      t226 = t6*t4
      t228 = t9*t7+t226*t3
      t229 = t228**2
      t230 = t229*t133
      t231 = t94*t122
      t232 = t52*t4
      t233 = t232*as
      t236 = 2*t233*t66*bp
      t238 = 2*t70*t152
      t241 = -t3*t7+t226*t9
      t242 = t241**2
      t243 = t242*t133
      t244 = t228*as
      t246 = 2*t244*t63
      t247 = t24*ac
      t249 = 2*t247*t244
      t252 = 2*t232*as*Zoi
      t257 = 2*t53*t7*ac*t4*as
      t258 = -t121+t224+t230+t231+t236-t238+t243-t246+t249-t252-t257
      t259 = t101*t122
      t260 = t241*as
      t262 = 2*t260*Yoj
      t263 = t11*ac
      t265 = 2*t263*Yoi
      t267 = 2*t244*Xoi
      t269 = 2*t244*Xoj
      t271 = 2*t247*Xoj
      t272 = t259-t262+t21-t43+t59-t62+t265+t65+t267-t269-t271
      t275 = 2*t260*Yoi
      t277 = 2*t263*Yoj
      t278 = t55*t122
      t280 = 2*t247*Xoi
      t281 = t275-t277+t278+t280+t69-t80-t168-t82+t84+t85+t86
      t283 = 2*t263*t260
      t285 = 2*t232*t157
      t288 = 2*t70*ac*Zoi
      t290 = 2*t260*t57
      t291 = t87+t88+t89+t283+t285-t97+t98+t288-t290+t100-t104-t186
      t294 = sqrt(t258+t272+t281+t291)
      t296 = (t294+0.1E-98)**2
      t304 = qh**2
      t305 = t304*ec
      t308 = 2*t228*t133*t144
      t310 = 2*t244*t119
      t313 = 2*t11*t122*t49
      t316 = 2*t24*t122*t38
      t317 = t224-t308+t230+t231-t238+t243+t249-t252-t257+t259+t123-t310
     #-t313-t316-t262
      t318 = 2*t141
      t319 = 2*t146
      t320 = 2*t147
      t321 = 2*t149
      t322 = 2*t151
      t326 = 2*t232*t133*t15*t32
      t327 = t265+t267-t269-t271+t275-t277+t130+t134-t136-t138-t318+t319
     #+t320-t321+t322-t326
      t331 = 2*t70*ac*t181
      t335 = 2*t70*t122*t15*t17
      t337 = 2*t233*t184
      t339 = 2*t247*t145
      t341 = 2*t260*t150
      t344 = 2*t241*t133*t127
      t345 = 2*t153
      t346 = 2*t155
      t347 = 2*t158
      t348 = 2*t160
      t349 = t331-t335+t337-t339-t341-t344+t278+t280+t345-t346-t347+t348
     #+t166+t172-t82+t85
      t351 = 2*t263*t128
      t352 = t86+t87+t88+t89+t283+t285-t97+t98+t288-t104+t174-t351-t179+
     #t188+t190+t191
      t355 = sqrt(t317+t327+t349+t352)
      t357 = (t355+0.1E-98)**2
      t365 = t224+t308+t230+t231-t238+t243+t249-t252-t257+t259+t123-t310
     #-t313-t316-t262
      t366 = t265+t267-t269-t271+t275-t277-t130+t134-t136-t138+t318-t319
     #+t320+t321+t322+t326
      t368 = -t331-t335+t337+t339-t341+t344+t278+t280+t345-t346+t347-t34
     #8-t166-t172-t82+t85
      t369 = t86+t87+t88+t89+t283+t285-t97+t98+t288-t104+t174+t351+t179+
     #t188+t190+t191
      t372 = sqrt(t365+t366+t368+t369)
      t374 = (t372+0.1E-98)**2
      t382 = -t121+t224+t230+t231-t236-t238+t243+t246-t249+t252+t257
      t383 = t259+t262+t21-t43+t59-t62+t265+t65-t267+t269-t271
      t385 = -t275-t277+t278+t280+t69-t80-t168-t82+t84+t85+t86
      t386 = t87+t88+t89-t283-t285-t97+t98+t288+t290+t100-t104-t186
      t389 = sqrt(t382+t383+t385+t386)
      t391 = (t389+0.1E-98)**2
      t399 = t224+t308+t230+t231-t238+t243-t249+t252+t257+t259+t123+t310
     #-t313-t316+t262
      t400 = t265-t267+t269-t271-t275-t277+t130+t134-t136-t138-t318+t319
     #+t320-t321+t322+t326
      t402 = t331-t335-t337-t339+t341+t344+t278+t280+t345-t346-t347+t348
     #+t166+t172-t82+t85
      t403 = t86+t87+t88+t89-t283-t285-t97+t98+t288-t104+t174-t351-t179+
     #t188+t190+t191
      t406 = sqrt(t399+t400+t402+t403)
      t408 = (t406+0.1E-98)**2
      t416 = t224-t308+t230+t231-t238+t243-t249+t252+t257+t259+t123+t310
     #-t313-t316+t262
      t417 = t265-t267+t269-t271-t275-t277-t130+t134-t136-t138+t318-t319
     #+t320+t321+t322-t326
      t419 = -t331-t335-t337+t339+t341-t344+t278+t280+t345-t346+t347-t34
     #8-t166-t172-t82+t85
      t420 = t86+t87+t88+t89-t283-t285-t97+t98+t288-t104+t174+t351+t179+
     #t188+t190+t191
      t423 = sqrt(t416+t417+t419+t420)
      t425 = (t423+0.1E-98)**2
      t434 = sqrt(t88-t104+t87+t86-t82+t85+t98-t97+t89)
      t435 = t434+0.1E-98
      t436 = t435**2
      t437 = t436**2
      t439 = t437**2
      t445 = 2/t434*(-Yoi+Yoj)
      t454 = -t1*ec/t109/t107*(-Yoi+Yoj+t57-t12)-t118/t197/t195*(-Yoi+Yo
     #j-t12+t128+t150)-t118/t214/t212*(-Yoi+Yoj-t12-t128+t150)-t118/t296
     #/t294*(-Yoi+Yoj+t57-t260-t263)-t305/t357/t355*(-Yoi+Yoj-t260-t263+
     #t128+t150)-t305/t374/t372*(-Yoi+Yoj-t260-t263-t128+t150)-t118/t391
     #/t389*(-Yoi+Yoj+t57+t260-t263)-t305/t408/t406*(-Yoi+Yoj+t260-t263+
     #t128+t150)-t305/t425/t423*(-Yoi+Yoj+t260-t263-t128+t150)-6*AA/t439
     #/t437/t435*t445+3*BB/t437/t436/t435*t445
c
      dVdYoj4p=t454
c
      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      function dVdZoj4p(Xoi,Xoj,Yoi,Yoj,Zoi,Zoj,ppi,si,ti,ppj,sj,tj)
      implicit real*8(a-h,o-z)
      common/cpot/qm,qh,AA,BB,ac,as,bp,ec
c
      t1 = qm**2
      t4 = 2*Xoi*Xoj
      t5 = sin(ppj)
      t6 = cos(sj)
      t8 = cos(tj)
      t9 = sin(sj)
      t10 = t8*t9
      t11 = cos(ppj)
      t13 = -t5*t6-t10*t11
      t14 = t13**2
      t15 = bp**2
      t16 = t14*t15
      t17 = sin(ppi)
      t18 = cos(si)
      t20 = cos(ti)
      t21 = sin(si)
      t22 = t20*t21
      t23 = cos(ppi)
      t25 = -t17*t18-t22*t23
      t26 = t25**2
      t27 = t26*t15
      t28 = sin(ti)
      t29 = t28*t21
      t32 = 2*t29*bp*Zoi
      t33 = Yoi**2
      t36 = t23*t18-t22*t17
      t37 = t36*bp
      t39 = 2*t37*Xoj
      t40 = t28**2
      t41 = t21**2
      t42 = t40*t41
      t43 = t42*t15
      t44 = bp*Zoj
      t46 = 2*t29*t44
      t47 = sin(tj)
      t48 = t47*t9
      t50 = 2*t48*t44
      t51 = Zoi*t47
      t54 = 2*t51*t9*bp
      t57 = t11*t6-t10*t5
      t58 = Xoi*t57
      t60 = 2*t58*bp
      t65 = 2*t37*Xoi
      t66 = t47**2
      t67 = t9**2
      t68 = t66*t67
      t69 = t68*t15
      t74 = -t4+t16+t27+t32+t33-t39+t43-t46+t50-t54-t60-2*t36*t15*t57+t6
     #5+t69-2*t29*t15*t47*t9
      t75 = Yoj**2
      t76 = Zoi**2
      t77 = Xoj**2
      t81 = t13*bp
      t83 = 2*t81*Yoj
      t84 = t57*bp
      t86 = 2*t84*Xoj
      t87 = t25*bp
      t89 = 2*t87*Yoi
      t91 = 2*t87*Yoj
      t92 = Yoi*t13
      t94 = 2*t92*bp
      t95 = Xoi**2
      t96 = Zoj**2
      t98 = 2*Zoi*Zoj
      t100 = 2*Yoi*Yoj
      t101 = t36**2
      t102 = t101*t15
      t103 = t57**2
      t104 = t103*t15
      t105 = t75+t76+t77-2*t25*t15*t13+t83+t86+t89-t91-t94+t95+t96-t98-t
     #100+t102+t104
      t107 = sqrt(t74+t105)
      t109 = (t107+0.1E-98)**2
      t113 = t48*bp
      t114 = t29*bp
      t120 = qm*qh*ec
      t122 = t8*t6
      t124 = t11*t9+t122*t5
      t127 = 2*Xoi*t124*as
      t128 = -t4+t27+t32+t33-t39+t43-t46-t127+t65+t75+t76
      t129 = ac**2
      t130 = t68*t129
      t131 = t13*ac
      t133 = 2*t131*Yoj
      t136 = 2*t51*t9*ac
      t137 = t77+t89-t91+t95+t96-t98-t100+t130+t133+t102-t136
      t141 = -t5*t9+t122*t11
      t142 = t141*as
      t144 = 2*t131*t142
      t146 = 2*t87*t142
      t147 = t57*ac
      t149 = 2*t37*t147
      t151 = 2*t87*t131
      t152 = t124*as
      t154 = 2*t147*t152
      t156 = 2*t37*t152
      t157 = ac*Zoj
      t159 = 2*t48*t157
      t162 = 2*t51*t6*as
      t163 = t47*t6
      t164 = as*Zoj
      t166 = 2*t163*t164
      t167 = t103*t129
      t169 = 2*t147*Xoj
      t170 = t144-t146-t149-t151+t154-t156+t159+t162-t166+t167+t169
      t172 = 2*t152*Xoj
      t174 = 2*t58*ac
      t176 = 2*t92*ac
      t177 = t14*t129
      t178 = t48*ac
      t180 = 2*t114*t178
      t181 = t163*as
      t183 = 2*t114*t181
      t188 = 2*t66*t9*ac*t6*as
      t189 = t141**2
      t190 = as**2
      t191 = t189*t190
      t192 = t124**2
      t193 = t192*t190
      t195 = 2*t142*Yoj
      t196 = t6**2
      t198 = t66*t196*t190
      t201 = 2*Yoi*t141*as
      t202 = t172-t174-t176+t177-t180+t183-t188+t191+t193+t195+t198-t201
      t205 = sqrt(t128+t137+t170+t202)
      t207 = (t205+0.1E-98)**2
      t215 = -t4+t27+t32+t33-t39+t43-t46+t127+t65+t75+t76
      t217 = -t144+t146-t149-t151-t154+t156+t159-t162+t166+t167+t169
      t218 = -t172-t174-t176+t177-t180-t183+t188+t191+t193-t195+t198+t20
     #1
      t221 = sqrt(t215+t137+t217+t218)
      t223 = (t221+0.1E-98)**2
      t231 = t42*t129
      t233 = t20*t18
      t235 = t23*t21+t233*t17
      t236 = t235**2
      t237 = t236*t190
      t238 = t36*ac
      t240 = 2*t238*Xoj
      t241 = t231-t4+t16+t237+t33-t240+t50-t54-t60+t69+t75
      t242 = t235*as
      t244 = 2*t242*Xoj
      t245 = t25*ac
      t247 = 2*t245*Yoj
      t248 = t76+t77+t83+t86-t244-t247-t94+t95+t96-t98-t100
      t251 = 2*t242*Xoi
      t253 = 2*t245*Yoi
      t256 = -t17*t21+t233*t23
      t257 = t256*as
      t259 = 2*t257*Yoi
      t261 = 2*t238*Xoi
      t263 = 2*t257*Yoj
      t264 = t18**2
      t266 = t40*t264*t190
      t271 = 2*t40*t21*ac*t18*as
      t272 = t251+t253+t259-t149-t151-t180+t104+t261-t263+t266-t271
      t273 = t28*t18
      t274 = t273*as
      t276 = 2*t274*t113
      t277 = t256**2
      t278 = t277*t190
      t279 = t101*t129
      t281 = 2*t29*t157
      t284 = 2*t29*ac*Zoi
      t286 = 2*t242*t84
      t288 = 2*t238*t242
      t290 = 2*t257*t81
      t292 = 2*t245*t257
      t295 = 2*t273*as*Zoi
      t297 = 2*t273*t164
      t298 = t26*t129
      t299 = t276+t278+t279-t281+t284-t286+t288-t290+t292-t295+t297+t298
      t302 = sqrt(t241+t248+t272+t299)
      t304 = (t302+0.1E-98)**2
      t308 = t29*ac
      t313 = qh**2
      t314 = t313*ec
      t317 = 2*t235*t190*t124
      t318 = t231-t4+t237+t33-t240-t127+t75+t76+t77-t317-t244-t247+t95+t
     #96-t98
      t319 = -t100+t130+t133+t251+t253+t259-t136+t144+t154+t159+t162-t16
     #6+t167+t169+t172-t174
      t323 = 2*t36*t129*t57
      t326 = 2*t25*t129*t13
      t327 = -t176+t177-t188+t191+t193+t195+t198-t201+t261-t263+t266-t32
     #3-t326-t271+t278+t279
      t329 = 2*t308*t181
      t333 = 2*t273*t190*t47*t6
      t335 = 2*t274*t178
      t339 = 2*t29*t129*t47*t9
      t342 = 2*t256*t190*t141
      t344 = 2*t245*t142
      t346 = 2*t257*t131
      t348 = 2*t242*t147
      t350 = 2*t238*t152
      t351 = -t281+t284+t329-t333+t335-t339-t342-t344-t346-t348-t350+t28
     #8+t292-t295+t297+t298
      t354 = sqrt(t318+t319+t327+t351)
      t356 = (t354+0.1E-98)**2
      t364 = t231-t4+t237+t33-t240+t127+t75+t76+t77+t317-t244-t247+t95+t
     #96-t98
      t365 = -t100+t130+t133+t251+t253+t259-t136-t144-t154+t159-t162+t16
     #6+t167+t169-t172-t174
      t367 = -t176+t177+t188+t191+t193-t195+t198+t201+t261-t263+t266-t32
     #3-t326-t271+t278+t279
      t368 = -t281+t284-t329+t333+t335-t339+t342+t344-t346-t348+t350+t28
     #8+t292-t295+t297+t298
      t371 = sqrt(t364+t365+t367+t368)
      t373 = (t371+0.1E-98)**2
      t381 = t76+t77+t83+t86+t244-t247-t94+t95+t96-t98-t100
      t383 = -t251+t253-t259-t149-t151-t180+t104+t261+t263+t266+t271
      t384 = -t276+t278+t279-t281+t284+t286-t288+t290-t292+t295-t297+t29
     #8
      t387 = sqrt(t241+t381+t383+t384)
      t389 = (t387+0.1E-98)**2
      t397 = t231-t4+t237+t33-t240-t127+t75+t76+t77+t317+t244-t247+t95+t
     #96-t98
      t398 = -t100+t130+t133-t251+t253-t259-t136+t144+t154+t159+t162-t16
     #6+t167+t169+t172-t174
      t400 = -t176+t177-t188+t191+t193+t195+t198-t201+t261+t263+t266-t32
     #3-t326+t271+t278+t279
      t401 = -t281+t284+t329+t333-t335-t339+t342-t344+t346+t348-t350-t28
     #8-t292+t295-t297+t298
      t404 = sqrt(t397+t398+t400+t401)
      t406 = (t404+0.1E-98)**2
      t414 = t231-t4+t237+t33-t240+t127+t75+t76+t77-t317+t244-t247+t95+t
     #96-t98
      t415 = -t100+t130+t133-t251+t253-t259-t136-t144-t154+t159-t162+t16
     #6+t167+t169-t172-t174
      t417 = -t176+t177+t188+t191+t193-t195+t198+t201+t261+t263+t266-t32
     #3-t326+t271+t278+t279
      t418 = -t281+t284-t329-t333-t335-t339-t342+t344+t346+t348+t350-t28
     #8-t292+t295-t297+t298
      t421 = sqrt(t414+t415+t417+t418)
      t423 = (t421+0.1E-98)**2
      t432 = sqrt(t95-t4+t77+t33-t100+t75+t76-t98+t96)
      t433 = t432+0.1E-98
      t434 = t433**2
      t435 = t434**2
      t437 = t435**2
      t443 = 2/t432*(-Zoi+Zoj)
      t452 = -t1*ec/t109/t107*(t113-t114-Zoi+Zoj)-t120/t207/t205*(-t114-
     #Zoi+Zoj-t181+t178)-t120/t223/t221*(-t114-Zoi+Zoj+t181+t178)-t120/t
     #304/t302*(t113-Zoi+Zoj+t274-t308)-t314/t356/t354*(-Zoi+Zoj+t274-t3
     #08-t181+t178)-t314/t373/t371*(-Zoi+Zoj+t274-t308+t181+t178)-t120/t
     #389/t387*(t113-Zoi+Zoj-t274-t308)-t314/t406/t404*(-Zoi+Zoj-t274-t3
     #08-t181+t178)-t314/t423/t421*(-Zoi+Zoj-t274-t308+t181+t178)-6*AA/t
     #437/t435/t433*t443+3*BB/t435/t434/t433*t443
c
      dVdZoj4p=t452
c
      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      function dVdtj4p(Xoi,Xoj,Yoi,Yoj,Zoi,Zoj,ppi,si,ti,ppj,sj,tj)
      implicit real*8(a-h,o-z)
      common/cpot/qm,qh,AA,BB,ac,as,bp,ec
c
      t1 = qm**2
      t3 = sin(ppj)
      t4 = cos(sj)
      t5 = t3*t4
      t6 = cos(tj)
      t7 = sin(sj)
      t8 = t6*t7
      t9 = cos(ppj)
      t11 = -t5-t8*t9
      t12 = t11**2
      t13 = bp**2
      t14 = t12*t13
      t15 = sin(ti)
      t16 = sin(si)
      t17 = t15*t16
      t18 = sin(tj)
      t23 = Xoj**2
      t24 = Yoj**2
      t25 = Yoi**2
      t26 = Xoi**2
      t27 = Zoj**2
      t28 = Zoi**2
      t29 = sin(ppi)
      t30 = cos(si)
      t32 = cos(ti)
      t33 = t32*t16
      t34 = cos(ppi)
      t36 = -t29*t30-t33*t34
      t37 = t36**2
      t38 = t37*t13
      t39 = t9*t4
      t41 = t39-t8*t3
      t42 = t41**2
      t43 = t42*t13
      t45 = 2*Zoi*Zoj
      t47 = 2*Yoi*Yoj
      t49 = 2*Xoi*Xoj
      t52 = t34*t30-t33*t29
      t53 = t52*bp
      t55 = 2*t53*Xoi
      t57 = 2*t53*Xoj
      t58 = t14-2*t17*t13*t18*t7+t23+t24+t25+t26+t27+t28+t38+t43-t45-t47
     #-t49+t55-t57
      t59 = t18**2
      t60 = t7**2
      t61 = t59*t60
      t62 = t61*t13
      t63 = t15**2
      t64 = t16**2
      t65 = t63*t64
      t66 = t65*t13
      t67 = t36*t13
      t70 = Xoi*t41
      t72 = 2*t70*bp
      t73 = t52*t13
      t76 = t36*bp
      t78 = 2*t76*Yoi
      t80 = 2*t76*Yoj
      t81 = Yoi*t11
      t83 = 2*t81*bp
      t84 = t11*bp
      t86 = 2*t84*Yoj
      t87 = t41*bp
      t89 = 2*t87*Xoj
      t90 = t52**2
      t91 = t90*t13
      t92 = Zoi*t18
      t93 = t7*bp
      t95 = 2*t92*t93
      t96 = bp*Zoj
      t98 = 2*t17*t96
      t101 = 2*t17*bp*Zoi
      t102 = t18*t7
      t104 = 2*t102*t96
      t105 = t62+t66-2*t67*t11-t72-2*t73*t41+t78-t80-t83+t86+t89+t91-t95
     #-t98+t101+t104
      t107 = sqrt(t58+t105)
      t109 = (t107+0.1E-98)**2
      t113 = Zoi*t6
      t114 = t113*t93
      t115 = t8*t96
      t117 = t102*t3
      t118 = t41*t13*t117
      t120 = t102*t9
      t121 = t11*t13*t120
      t122 = Yoi*t18
      t123 = t7*t9
      t124 = t123*bp
      t125 = t122*t124
      t128 = t102*t9*bp*Yoj
      t131 = t102*t3*bp*Xoj
      t133 = Xoi*t18
      t134 = t7*t3
      t135 = t134*bp
      t136 = t133*t135
      t138 = t18*t60
      t139 = t13*t6
      t140 = t138*t139
      t143 = -t114+t115+t118+t121-t125+t128+t131-t67*t120-t136-t73*t117+
     #t140-t17*t139*t7
      t147 = qm*qh*ec
      t148 = ac**2
      t149 = t12*t148
      t150 = t4*as
      t152 = 2*t92*t150
      t153 = t7*ac
      t155 = 2*t92*t153
      t156 = t18*t4
      t157 = as*Zoj
      t159 = 2*t156*t157
      t160 = t4*t6
      t162 = t123+t160*t3
      t163 = t162**2
      t164 = as**2
      t165 = t163*t164
      t166 = t42*t148
      t168 = -t134+t160*t9
      t169 = t168**2
      t170 = t169*t164
      t171 = t11*ac
      t172 = t168*as
      t174 = 2*t171*t172
      t176 = 2*t76*t172
      t177 = t41*ac
      t179 = 2*t53*t177
      t181 = 2*t76*t171
      t182 = t149+t152-t155-t159+t165+t166+t170+t174-t176-t179-t181
      t183 = t162*as
      t185 = 2*t177*t183
      t187 = 2*t53*t183
      t188 = ac*Zoj
      t190 = 2*t102*t188
      t192 = 2*t171*Yoj
      t193 = t4**2
      t195 = t59*t193*t164
      t196 = t185-t187+t190+t192+t195+t23+t24+t25+t26+t27+t28
      t198 = t17*bp
      t199 = t102*ac
      t201 = 2*t198*t199
      t202 = t156*as
      t204 = 2*t198*t202
      t209 = 2*t59*t7*ac*t4*as
      t212 = 2*Xoi*t162*as
      t213 = t38-t45-t47-t49-t201+t204-t209+t55-t57+t66-t212
      t215 = 2*t70*ac
      t217 = 2*t81*ac
      t218 = t61*t148
      t220 = 2*t172*Yoj
      t222 = 2*t177*Xoj
      t224 = 2*t183*Xoj
      t227 = 2*Yoi*t168*as
      t228 = t78-t80-t215-t217+t218+t220+t222+t224-t227+t91-t98+t101
      t231 = sqrt(t182+t196+t213+t228)
      t233 = (t231+0.1E-98)**2
      t237 = t39*as
      t238 = t122*t237
      t239 = t134*ac
      t240 = t133*t239
      t241 = t123*ac
      t242 = t122*t241
      t243 = t8*t188
      t244 = t148*t6
      t245 = t138*t244
      t247 = t164*t6
      t248 = t18*t193*t247
      t249 = t113*t150
      t250 = t160*t157
      t253 = t102*t3*ac*Xoj
      t256 = t156*t3*as*Xoj
      t259 = t156*t9*as*Yoj
      t260 = t5*as
      t261 = t133*t260
      t264 = t120*ac*t168*as
      t266 = t171*t18*t237
      t267 = t238-t240-t242+t243+t245+t248+t249-t250+t253-t256-t259+t261
     #+t264-t266
      t268 = t113*t153
      t269 = 2*t268
      t270 = t76*t18
      t272 = 2*t270*t237
      t273 = t270*t241
      t274 = 2*t273
      t275 = t53*t18
      t276 = t275*t239
      t277 = 2*t276
      t280 = t117*ac*t162*as
      t281 = 2*t280
      t283 = t156*t9
      t284 = t168*t164*t283
      t285 = 2*t284
      t288 = 2*t41*t148*t117
      t289 = t160*as
      t291 = 2*t198*t289
      t293 = 4*t199*t289
      t294 = t8*ac
      t295 = t198*t294
      t296 = 2*t295
      t298 = t156*t3
      t300 = 2*t162*t164*t298
      t303 = 2*t177*t18*t260
      t305 = 2*t275*t260
      t308 = 2*t11*t148*t120
      t312 = 2*t102*t9*ac*Yoj
      t313 = -t269+t272-t274-t277+t281-t285+t288+t291-t293-t296-t300-t30
     #3+t305+t308+t312
      t317 = t149-t152-t155+t159+t165+t166+t170-t174+t176-t179-t181
      t318 = -t185+t187+t190+t192+t195+t23+t24+t25+t26+t27+t28
      t320 = t38-t45-t47-t49-t201-t204+t209+t55-t57+t66+t212
      t321 = t78-t80-t215-t217+t218-t220+t222-t224+t227+t91-t98+t101
      t324 = sqrt(t317+t318+t320+t321)
      t326 = (t324+0.1E-98)**2
      t330 = -t238-t240-t242+t243+t245+t248-t249+t250+t253+t256+t259-t26
     #1-t264+t266
      t331 = -t269-t272-t274-t277-t281-t285+t288-t291+t293-t296-t300+t30
     #3-t305+t308+t312
      t335 = -t179-t181+t14+t23+t24+t25+t26+t27+t28+t43-t45
      t336 = t36*ac
      t338 = 2*t336*Yoi
      t340 = t32*t30
      t342 = -t29*t16+t340*t34
      t343 = t342*as
      t345 = 2*t343*Yoj
      t348 = t34*t16+t340*t29
      t349 = t348*as
      t351 = 2*t349*Xoi
      t352 = t37*t148
      t354 = 2*t336*Yoj
      t355 = -t47-t49-t201+t338+t62-t345-t72+t351+t352-t354-t83
      t357 = t15*t30
      t358 = t357*as
      t361 = 2*t358*t102*bp
      t366 = 2*t63*t16*ac*t30*as
      t367 = t90*t148
      t368 = t342**2
      t369 = t368*t164
      t371 = 2*t17*t188
      t374 = 2*t17*ac*Zoi
      t376 = 2*t349*t87
      t377 = t52*ac
      t379 = 2*t377*t349
      t381 = 2*t343*t84
      t382 = t86+t89+t361-t366+t367+t369-t371+t374-t376+t379-t381
      t384 = 2*t336*t343
      t387 = 2*t357*as*Zoi
      t389 = 2*t357*t157
      t390 = t348**2
      t391 = t390*t164
      t393 = 2*t377*Xoj
      t394 = t30**2
      t396 = t63*t394*t164
      t397 = t65*t148
      t399 = 2*t377*Xoi
      t401 = 2*t349*Xoj
      t403 = 2*t343*Yoi
      t404 = t384-t387+t389+t391-t393+t396+t397+t399-t401+t403-t95+t104
      t407 = sqrt(t335+t355+t382+t404)
      t409 = (t407+0.1E-98)**2
      t414 = t358*t8*bp
      t415 = t343*t18
      t416 = t415*t124
      t417 = t349*t18
      t418 = t417*t135
      t419 = -t114+t115+t118+t121-t125+t128+t131-t136+t140+t414-t416-t41
     #8-t273-t276-t295
      t422 = qh**2
      t423 = t422*ec
      t424 = t149+t152-t155-t159+t165+t166+t170+t174+t185+t190+t192+t195
     #+t23+t24+t25
      t425 = t348*t164
      t427 = 2*t425*t162
      t428 = t36*t148
      t430 = 2*t428*t11
      t431 = t26+t27+t28-t45-t47-t49-t209+t338-t345-t212+t351-t427-t215+
     #t352-t217-t430
      t433 = t342*t164
      t435 = 2*t433*t168
      t437 = 2*t349*t177
      t439 = 2*t377*t183
      t441 = 2*t336*t172
      t443 = 2*t343*t171
      t444 = -t354+t218-t435-t366+t367+t220+t369-t437-t439-t441-t443+t22
     #2-t371+t374+t379+t384
      t445 = t52*t148
      t447 = 2*t445*t41
      t451 = 2*t357*t164*t18*t4
      t452 = t17*ac
      t454 = 2*t452*t202
      t456 = 2*t358*t199
      t460 = 2*t17*t148*t18*t7
      t461 = -t387+t389+t224-t227+t391-t447-t451+t454+t456-t460-t393+t39
     #6+t397+t399-t401+t403
      t464 = sqrt(t424+t431+t444+t461)
      t466 = (t464+0.1E-98)**2
      t470 = t238-t240-t242+t243+t245+t248+t249-t250+t253-t256-t259+t261
     #+t264-t266-t268+t280-t284
      t472 = 2*t425*t298
      t474 = 2*t445*t117
      t477 = 2*t336*t18*t237
      t480 = 2*t17*t244*t7
      t482 = 2*t358*t294
      t484 = 2*t452*t289
      t486 = 2*t415*t241
      t489 = 2*t377*t18*t260
      t491 = 2*t417*t239
      t494 = 2*t357*t247*t4
      t496 = 2*t428*t120
      t498 = 2*t433*t283
      t499 = t288-t293-t300-t303+t308+t312+t472-t474+t477-t480+t482+t484
     #-t486+t489-t491-t494-t496+t498
      t503 = t149-t152-t155+t159+t165+t166+t170-t174-t185+t190+t192+t195
     #+t23+t24+t25
      t504 = t26+t27+t28-t45-t47-t49+t209+t338-t345+t212+t351+t427-t215+
     #t352-t217-t430
      t506 = -t354+t218+t435-t366+t367-t220+t369-t437+t439+t441-t443+t22
     #2-t371+t374+t379+t384
      t507 = -t387+t389-t224+t227+t391-t447+t451-t454+t456-t460-t393+t39
     #6+t397+t399-t401+t403
      t510 = sqrt(t503+t504+t506+t507)
      t512 = (t510+0.1E-98)**2
      t516 = -t238-t240-t242+t243+t245+t248-t249+t250+t253+t256+t259-t26
     #1-t264+t266-t268-t280-t284
      t517 = t288+t293-t300+t303+t308+t312-t472-t474-t477-t480+t482-t484
     #-t486-t489-t491+t494-t496-t498
      t521 = -t47-t49-t201+t338+t62+t345-t72-t351+t352-t354-t83
      t523 = t86+t89-t361+t366+t367+t369-t371+t374+t376-t379+t381
      t524 = -t384+t387-t389+t391-t393+t396+t397+t399+t401-t403-t95+t104
      t527 = sqrt(t335+t521+t523+t524)
      t529 = (t527+0.1E-98)**2
      t533 = -t114+t115+t118+t121-t125+t128+t131-t136+t140-t414+t416+t41
     #8-t273-t276-t295
      t536 = t26+t27+t28-t45-t47-t49-t209+t338+t345-t212-t351+t427-t215+
     #t352-t217-t430
      t538 = -t354+t218+t435+t366+t367+t220+t369+t437-t439-t441+t443+t22
     #2-t371+t374-t379-t384
      t539 = t387-t389+t224-t227+t391-t447+t451+t454-t456-t460-t393+t396
     #+t397+t399+t401-t403
      t542 = sqrt(t424+t536+t538+t539)
      t544 = (t542+0.1E-98)**2
      t548 = t288-t293-t300-t303+t308+t312-t472-t474+t477-t480-t482+t484
     #+t486+t489+t491+t494-t496-t498
      t552 = t26+t27+t28-t45-t47-t49+t209+t338+t345+t212-t351-t427-t215+
     #t352-t217-t430
      t554 = -t354+t218-t435+t366+t367-t220+t369+t437+t439+t441+t443+t22
     #2-t371+t374-t379-t384
      t555 = t387-t389-t224+t227+t391-t447-t451-t454-t456-t460-t393+t396
     #+t397+t399+t401-t403
      t558 = sqrt(t503+t552+t554+t555)
      t560 = (t558+0.1E-98)**2
      t564 = t288+t293-t300+t303+t308+t312+t472-t474-t477-t480-t482-t484
     #+t486-t489+t491-t494-t496+t498
      t569 = -t1*ec/t109/t107*t143-t147/t233/t231*(2*t267+t313)/2-t147/t
     #326/t324*(2*t330+t331)/2-t147/t409/t407*t419-t423/t466/t464*(2*t47
     #0+t499)/2-t423/t512/t510*(2*t516+t517)/2-t147/t529/t527*t533-t423/
     #t544/t542*(2*t470+t548)/2-t423/t560/t558*(2*t516+t564)/2
c
      dVdtj4p=t569
c
      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      function dVdsj4p(Xoi,Xoj,Yoi,Yoj,Zoi,Zoj,ppi,si,ti,ppj,sj,tj)
      implicit real*8(a-h,o-z)
      common/cpot/qm,qh,AA,BB,ac,as,bp,ec
c
      t1 = qm**2
      t3 = sin(tj)
      t4 = t3**2
      t5 = sin(sj)
      t6 = t5**2
      t7 = t4*t6
      t8 = bp**2
      t9 = t7*t8
      t10 = sin(ti)
      t11 = t10**2
      t12 = sin(si)
      t13 = t12**2
      t14 = t11*t13
      t15 = t14*t8
      t16 = sin(ppi)
      t17 = cos(si)
      t19 = cos(ti)
      t20 = t19*t12
      t21 = cos(ppi)
      t23 = -t16*t17-t20*t21
      t24 = t23*bp
      t26 = 2*t24*Yoi
      t29 = t21*t17-t20*t16
      t30 = t29*t8
      t31 = cos(ppj)
      t32 = cos(sj)
      t34 = cos(tj)
      t35 = t34*t5
      t36 = sin(ppj)
      t38 = t31*t32-t35*t36
      t41 = t29*bp
      t43 = 2*t41*Xoi
      t46 = -t36*t32-t35*t31
      t47 = Yoi*t46
      t49 = 2*t47*bp
      t50 = t46*bp
      t52 = 2*t50*Yoj
      t53 = t38*bp
      t55 = 2*t53*Xoj
      t56 = t23*t8
      t60 = 2*t24*Yoj
      t61 = Xoi*t38
      t63 = 2*t61*bp
      t65 = 2*t41*Xoj
      t66 = t10*t12
      t69 = 2*t66*bp*Zoi
      t70 = bp*Zoj
      t72 = 2*t66*t70
      t73 = t23**2
      t74 = t73*t8
      t75 = t9+t15+t26-2*t30*t38+t43-t49+t52+t55-2*t56*t46-t60-t63-t65+t
     #69-t72+t74
      t76 = Zoi*t3
      t79 = 2*t76*t5*bp
      t81 = 2*Xoi*Xoj
      t82 = t38**2
      t83 = t82*t8
      t85 = 2*Yoi*Yoj
      t86 = Xoi**2
      t87 = Zoj**2
      t88 = Zoi**2
      t89 = Yoj**2
      t90 = Yoi**2
      t92 = 2*Zoi*Zoj
      t93 = t46**2
      t94 = t93*t8
      t95 = t3*t5
      t97 = 2*t95*t70
      t98 = t29**2
      t99 = t98*t8
      t100 = t8*t3
      t104 = Xoj**2
      t105 = -t79-t81+t83-t85+t86+t87+t88+t89+t90-t92+t94+t97+t99-2*t66*
     #t100*t5+t104
      t107 = sqrt(t75+t105)
      t109 = (t107+0.1E-98)**2
      t114 = t76*t32*bp
      t115 = t3*t32
      t116 = t115*t70
      t119 = t34*t32
      t121 = -t31*t5-t119*t36
      t122 = t38*t8*t121
      t126 = t36*t5-t119*t31
      t127 = t46*t8*t126
      t128 = Yoi*t126
      t129 = t128*bp
      t130 = t126*bp
      t131 = t130*Yoj
      t132 = t121*bp
      t133 = t132*Xoj
      t135 = Xoi*t121
      t136 = t135*bp
      t138 = t4*t5
      t140 = t138*t8*t32
      t143 = -t114+t116+t122+t127-t129+t131+t133-t56*t126-t136-t30*t121+
     #t140-t66*t100*t32
      t147 = qm*qh*ec
      t148 = t66*bp
      t149 = t95*ac
      t151 = 2*t148*t149
      t152 = t121**2
      t153 = as**2
      t154 = t152*t153
      t155 = t38*ac
      t156 = -t121*as
      t158 = 2*t155*t156
      t160 = 2*t41*t156
      t163 = -2*Yoi*t126*as
      t164 = t115*as
      t166 = 2*t148*t164
      t167 = ac**2
      t168 = t93*t167
      t169 = ac*Zoj
      t171 = 2*t95*t169
      t172 = t126**2
      t173 = t172*t153
      t174 = ac*t32
      t177 = 2*t138*t174*as
      t178 = as*Zoj
      t180 = 2*t115*t178
      t181 = -t151+t154+t158-t160-t163+t166+t168+t171+t173-t177-t180
      t184 = 2*t76*t32*as
      t185 = t82*t167
      t187 = 2*t61*ac
      t189 = 2*t47*ac
      t190 = t46*ac
      t192 = 2*t190*Yoj
      t195 = -2*Xoi*t121*as
      t197 = 2*t155*Xoj
      t198 = -t126*as
      t200 = 2*t198*Yoj
      t201 = t7*t167
      t202 = t184+t185+t15+t26-t187-t189+t192-t195+t197+t200+t201
      t204 = t32**2
      t205 = t4*t204
      t206 = t205*t153
      t208 = 2*t156*Xoj
      t209 = t206+t208+t43-t60-t65+t69-t72+t74-t81-t85+t86
      t211 = 2*t24*t198
      t213 = 2*t41*t155
      t215 = 2*t24*t190
      t218 = 2*t76*t5*ac
      t220 = 2*t190*t198
      t221 = t87+t88+t89+t90-t92+t99+t104-t211-t213-t215-t218+t220
      t224 = sqrt(t181+t202+t209+t221)
      t226 = (t224+0.1E-98)**2
      t230 = t47*as
      t232 = t38*as*Xoj
      t233 = t29*as
      t234 = t233*t53
      t236 = t148*t115*ac
      t237 = t61*as
      t238 = t121*ac
      t239 = t238*Xoj
      t241 = t46*as*Yoj
      t243 = t82*ac*as
      t244 = t238*t156
      t246 = t76*t5*as
      t247 = t128*ac
      t248 = t23*as
      t249 = t248*t50
      t251 = t138*t167*t32
      t252 = t126*ac
      t253 = t252*t198
      t254 = ac*as
      t255 = t205*t254
      t256 = -t230+t232-t234-t236-t237+t239+t241+t243+t244-t246-t247-t24
     #9+t251+t253-t255
      t258 = t38*t167*t121
      t259 = t24*t252
      t260 = t41*t238
      t262 = -t126*t153*t46
      t264 = t46*t167*t126
      t266 = -t121*t153*t38
      t267 = t135*ac
      t268 = t252*Yoj
      t269 = t76*t174
      t270 = t95*t178
      t271 = t115*t169
      t273 = t93*ac*as
      t274 = t66*as
      t275 = t95*bp
      t276 = t274*t275
      t279 = t4*t32*t153*t5
      t280 = t7*t254
      t281 = t258-t259-t260+t262+t264+t266-t267+t268-t269+t270+t271+t273
     #-t276-t279+t280
      t285 = -t151+t154-t158+t160+t163-t166+t168+t171+t173+t177+t180
      t286 = -t184+t185+t15+t26-t187-t189+t192+t195+t197-t200+t201
      t288 = t206-t208+t43-t60-t65+t69-t72+t74-t81-t85+t86
      t289 = t87+t88+t89+t90-t92+t99+t104+t211-t213-t215-t218-t220
      t292 = sqrt(t285+t286+t288+t289)
      t294 = (t292+0.1E-98)**2
      t298 = t230-t232+t234-t236+t237+t239-t241-t243-t244+t246-t247+t249
     #+t251-t253+t255
      t299 = t258-t259-t260+t262+t264+t266-t267+t268-t269-t270+t271-t273
     #+t276-t279-t280
      t304 = t19*t17
      t306 = -t16*t12+t304*t21
      t307 = t306*as
      t309 = 2*t307*Yoi
      t312 = t21*t12+t304*t16
      t313 = t312*as
      t315 = 2*t313*Xoj
      t317 = 2*t307*Yoj
      t318 = t23*ac
      t320 = 2*t318*Yoj
      t321 = t29*ac
      t323 = 2*t321*Xoi
      t324 = -t151+t9+t309-t315-t317-t320-t49+t52+t55-t63+t323
      t325 = t17**2
      t327 = t11*t325*t153
      t329 = 2*t321*Xoj
      t331 = 2*t318*Yoi
      t332 = t14*t167
      t333 = t10*t17
      t335 = 2*t333*t178
      t336 = t98*t167
      t337 = t333*as
      t339 = 2*t337*t275
      t340 = t73*t167
      t343 = 2*t333*as*Zoi
      t345 = 2*t318*t307
      t350 = 2*t11*t12*ac*t17*as
      t351 = t327-t329+t331+t332+t335+t336+t339+t340-t343+t345-t350
      t353 = t306**2
      t354 = t353*t153
      t355 = t312**2
      t356 = t355*t153
      t358 = 2*t313*Xoi
      t360 = 2*t66*t169
      t362 = 2*t313*t53
      t365 = 2*t66*ac*Zoi
      t367 = 2*t307*t50
      t369 = 2*t321*t313
      t370 = t354+t356+t358-t360-t362+t365-t367+t369-t79-t81+t83
      t371 = -t85+t86+t87+t88+t89+t90-t92+t94+t97+t104-t213-t215
      t374 = sqrt(t324+t351+t370+t371)
      t376 = (t374+0.1E-98)**2
      t381 = t333*bp*t164
      t382 = t307*t130
      t383 = t313*t132
      t384 = -t114+t116+t122+t127-t129+t131+t133-t136+t140+t381-t382-t38
     #3-t259-t260-t236
      t387 = qh**2
      t388 = t387*ec
      t389 = t29*t167
      t391 = 2*t389*t38
      t392 = -t391+t154+t158-t163+t309+t168+t171+t173-t177-t180+t184+t18
     #5-t315-t187-t189
      t393 = t306*t153
      t395 = -2*t393*t126
      t396 = t312*t153
      t398 = -2*t396*t121
      t399 = t192-t195+t197+t200+t201+t206-t395-t317-t398+t208-t320+t323
     #+t327-t329+t331+t332
      t402 = 2*t313*t155
      t404 = 2*t321*t156
      t405 = t23*t167
      t407 = 2*t405*t46
      t409 = 2*t307*t190
      t410 = t335+t336+t340-t343+t345-t350+t354+t356+t358-t360+t365+t369
     #-t402-t404-t407-t409
      t411 = t167*t3
      t414 = 2*t66*t411*t5
      t416 = 2*t337*t149
      t417 = t153*t3
      t420 = 2*t333*t417*t32
      t423 = 2*t66*ac*t164
      t425 = 2*t318*t198
      t426 = -t414+t416-t420+t423-t425-t81-t85+t86+t87+t88+t89+t90-t92+t
     #104-t218+t220
      t429 = sqrt(t392+t399+t410+t426)
      t431 = (t429+0.1E-98)**2
      t435 = t393*t46
      t436 = t389*t121
      t437 = t396*t38
      t438 = -t230+t232-t237+t239+t241+t243-t435-t436+t244-t246-t437-t24
     #7+t251+t253-t255+t258+t262+t264
      t439 = t307*t252
      t441 = t66*t411*t32
      t443 = t333*t417*t5
      t444 = t313*t238
      t445 = t248*t190
      t447 = t333*ac*t164
      t448 = t233*t155
      t449 = t274*t149
      t450 = t405*t126
      t451 = t266-t267+t268-t269+t270+t271+t273-t279+t280-t439-t441+t443
     #-t444-t445+t447-t448-t449-t450
      t455 = -t391+t154-t158+t163+t309+t168+t171+t173+t177+t180-t184+t18
     #5-t315-t187-t189
      t456 = t192+t195+t197-t200+t201+t206+t395-t317+t398-t208-t320+t323
     #+t327-t329+t331+t332
      t458 = t335+t336+t340-t343+t345-t350+t354+t356+t358-t360+t365+t369
     #-t402+t404-t407-t409
      t459 = -t414+t416+t420-t423+t425-t81-t85+t86+t87+t88+t89+t90-t92+t
     #104-t218-t220
      t462 = sqrt(t455+t456+t458+t459)
      t464 = (t462+0.1E-98)**2
      t468 = t230-t232+t237+t239-t241-t243+t435-t436-t244+t246+t437-t247
     #+t251-t253+t255+t258+t262+t264
      t469 = t266-t267+t268-t269-t270+t271-t273-t279-t280-t439-t441-t443
     #-t444+t445+t447+t448+t449-t450
      t473 = -t151+t9-t309+t315+t317-t320-t49+t52+t55-t63+t323
      t474 = t327-t329+t331+t332-t335+t336-t339+t340+t343-t345+t350
      t476 = t354+t356-t358-t360+t362+t365+t367-t369-t79-t81+t83
      t479 = sqrt(t473+t474+t476+t371)
      t481 = (t479+0.1E-98)**2
      t485 = -t114+t116+t122+t127-t129+t131+t133-t136+t140-t381+t382+t38
     #3-t259-t260-t236
      t488 = -t391+t154+t158-t163-t309+t168+t171+t173-t177-t180+t184+t18
     #5+t315-t187-t189
      t489 = t192-t195+t197+t200+t201+t206+t395+t317+t398+t208-t320+t323
     #+t327-t329+t331+t332
      t491 = -t335+t336+t340+t343-t345+t350+t354+t356-t358-t360+t365-t36
     #9+t402-t404-t407+t409
      t492 = -t414-t416+t420+t423-t425-t81-t85+t86+t87+t88+t89+t90-t92+t
     #104-t218+t220
      t495 = sqrt(t488+t489+t491+t492)
      t497 = (t495+0.1E-98)**2
      t501 = -t230+t232-t237+t239+t241+t243+t435-t436+t244-t246+t437-t24
     #7+t251+t253-t255+t258+t262+t264
      t502 = t266-t267+t268-t269+t270+t271+t273-t279+t280+t439-t441-t443
     #+t444-t445-t447-t448-t449-t450
      t506 = -t391+t154-t158+t163-t309+t168+t171+t173+t177+t180-t184+t18
     #5+t315-t187-t189
      t507 = t192+t195+t197-t200+t201+t206-t395+t317-t398-t208-t320+t323
     #+t327-t329+t331+t332
      t509 = -t335+t336+t340+t343-t345+t350+t354+t356-t358-t360+t365-t36
     #9+t402+t404-t407+t409
      t510 = -t414-t416-t420-t423+t425-t81-t85+t86+t87+t88+t89+t90-t92+t
     #104-t218-t220
      t513 = sqrt(t506+t507+t509+t510)
      t515 = (t513+0.1E-98)**2
      t519 = t230-t232+t237+t239-t241-t243-t435-t436-t244+t246-t437-t247
     #+t251-t253+t255+t258+t262+t264
      t520 = t266-t267+t268-t269-t270+t271-t273-t279-t280+t439-t441+t443
     #+t444+t445-t447+t448+t449-t450
      t525 = -t1*ec/t109/t107*t143-t147/t226/t224*(2*t256+2*t281)/2-t147
     #/t294/t292*(2*t298+2*t299)/2-t147/t376/t374*t384-t388/t431/t429*(2
     #*t438+2*t451)/2-t388/t464/t462*(2*t468+2*t469)/2-t147/t481/t479*t4
     #85-t388/t497/t495*(2*t501+2*t502)/2-t388/t515/t513*(2*t519+2*t520)
     #/2
c
      Dvdsj4p=t525
c
      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      function dVdppj4p(Xoi,Xoj,Yoi,Yoj,Zoi,Zoj,ppi,si,ti,ppj,sj,tj)
      implicit real*8(a-h,o-z)
      common/cpot/qm,qh,AA,BB,ac,as,bp,ec
c
      t1 = qm**2
      t3 = sin(ppi)
      t4 = cos(si)
      t6 = cos(ti)
      t7 = sin(si)
      t8 = t6*t7
      t9 = cos(ppi)
      t11 = -t3*t4-t8*t9
      t12 = t11*bp
      t14 = 2*t12*Yoi
      t15 = Zoi*Zoj
      t16 = 2*t15
      t17 = Yoi*Yoj
      t18 = 2*t17
      t21 = t9*t4-t8*t3
      t22 = t21**2
      t23 = bp**2
      t24 = t22*t23
      t25 = t11**2
      t26 = t25*t23
      t27 = cos(ppj)
      t28 = cos(sj)
      t30 = cos(tj)
      t31 = sin(sj)
      t32 = t30*t31
      t33 = sin(ppj)
      t35 = t27*t28-t32*t33
      t36 = t35**2
      t37 = t36*t23
      t40 = -t33*t28-t32*t27
      t41 = t40**2
      t42 = t41*t23
      t43 = sin(ti)
      t44 = t43**2
      t45 = t7**2
      t46 = t44*t45
      t47 = t46*t23
      t49 = 2*Xoi*Xoj
      t50 = sin(tj)
      t51 = Zoi*t50
      t54 = 2*t51*t31*bp
      t55 = t50*t31
      t56 = bp*Zoj
      t58 = 2*t55*t56
      t59 = t43*t7
      t61 = 2*t59*t56
      t64 = 2*t59*bp*Zoi
      t65 = Yoi*t40
      t67 = 2*t65*bp
      t68 = t40*bp
      t70 = 2*t68*Yoj
      t71 = t14-t16-t18+t24+t26+t37+t42+t47-t49-t54+t58-t61+t64-t67+t70
      t72 = t35*bp
      t74 = 2*t72*Xoj
      t75 = t11*t23
      t79 = 2*t12*Yoj
      t80 = Xoi*t35
      t82 = 2*t80*bp
      t83 = t21*t23
      t86 = t21*bp
      t88 = 2*t86*Xoi
      t90 = 2*t86*Xoj
      t91 = t50**2
      t92 = t31**2
      t93 = t91*t92
      t94 = t93*t23
      t95 = Xoj**2
      t96 = Yoi**2
      t97 = Yoj**2
      t102 = Zoi**2
      t103 = Zoj**2
      t104 = Xoi**2
      t105 = t74-2*t75*t40-t79-t82-2*t83*t35+t88-t90+t94+t95+t96+t97-2*t
     #59*t23*t50*t31+t102+t103+t104
      t107 = sqrt(t71+t105)
      t109 = (t107+0.1E-98)**2
      t114 = t35*t23*t40
      t116 = -t35*t23*t40
      t117 = -Yoi*t35
      t118 = t117*bp
      t119 = -t35*bp
      t120 = t119*Yoj
      t121 = t68*Xoj
      t123 = Xoi*t40
      t124 = t123*bp
      t130 = qm*qh*ec
      t131 = ac*Zoj
      t133 = 2*t55*t131
      t136 = 2*t51*t28*as
      t137 = t50*t28
      t138 = as*Zoj
      t140 = 2*t137*t138
      t141 = ac**2
      t142 = t41*t141
      t144 = t30*t28
      t146 = -t33*t31+t144*t27
      t149 = 2*Yoi*t146*as
      t150 = t35*ac
      t152 = 2*t150*Xoj
      t155 = t27*t31+t144*t33
      t158 = 2*Xoi*t155*as
      t159 = t155**2
      t160 = as**2
      t161 = t159*t160
      t162 = t36*t141
      t163 = t59*bp
      t164 = t55*ac
      t166 = 2*t163*t164
      t167 = t133+t136-t140+t142-t149+t152-t158+t161+t162+t14-t166
      t172 = 2*t91*t31*ac*t28*as
      t173 = t155*as
      t175 = 2*t86*t173
      t176 = t137*as
      t178 = 2*t163*t176
      t179 = t146**2
      t180 = t179*t160
      t181 = t40*ac
      t183 = 2*t181*Yoj
      t185 = 2*t65*ac
      t186 = -t172-t175+t178+t180+t183-t185-t16-t18+t24+t26+t47
      t188 = t80*ac
      t189 = 2*t188
      t190 = t173*Xoj
      t191 = 2*t190
      t193 = 2*t150*t173
      t195 = 2*t12*t181
      t197 = 2*t86*t150
      t198 = t146*as
      t200 = 2*t198*Yoj
      t201 = t28**2
      t203 = t91*t201*t160
      t204 = -t189+t191-t49-t61+t64+t193-t195-t197+t200+t203-t79
      t205 = t93*t141
      t207 = 2*t12*t198
      t209 = 2*t181*t198
      t212 = 2*t51*t31*ac
      t213 = t88-t90+t205+t95+t96+t97+t102+t103+t104-t207+t209-t212
      t216 = sqrt(t167+t186+t204+t213)
      t218 = (t216+0.1E-98)**2
      t222 = t86*t198
      t223 = t181*t173
      t224 = t150*t198
      t225 = -t35*ac
      t226 = t12*t225
      t227 = t86*t181
      t228 = -t155*as
      t229 = t12*t228
      t230 = t225*t198
      t231 = t181*t228
      t232 = t228*Yoj
      t233 = t225*Yoj
      t235 = Xoi*t146*as
      t236 = t181*Xoj
      t237 = t198*Xoj
      t239 = -Yoi*t155*as
      t240 = t123*ac
      t241 = t117*ac
      t243 = t35*t141*t40
      t245 = -t146*t160*t155
      t247 = -t35*t141*t40
      t249 = t146*t160*t155
      t250 = -t222+t223+t224-t226-t227-t229+t230+t231+t232+t233-t235+t23
     #6+t237-t239-t240-t241+t243+t245+t247+t249
      t253 = t133-t136+t140+t142+t149+t152+t158+t161+t162+t14-t166
      t254 = t172+t175-t178+t180+t183-t185-t16-t18+t24+t26+t47
      t256 = -t189-t191-t49-t61+t64-t193-t195-t197-t200+t203-t79
      t257 = t88-t90+t205+t95+t96+t97+t102+t103+t104+t207-t209-t212
      t260 = sqrt(t253+t254+t256+t257)
      t262 = (t260+0.1E-98)**2
      t266 = t222-t223-t224-t226-t227+t229-t230-t231-t232+t233+t235+t236
     #-t237+t239-t240-t241+t243+t245+t247+t249
      t269 = t21*ac
      t271 = 2*t269*Xoj
      t272 = t4**2
      t274 = t44*t272*t160
      t276 = 2*t269*Xoi
      t277 = t11*ac
      t279 = 2*t277*Yoj
      t280 = t46*t141
      t282 = t6*t4
      t284 = -t3*t7+t282*t9
      t285 = t284*as
      t287 = 2*t285*Yoi
      t290 = t9*t7+t282*t3
      t291 = t290*as
      t293 = 2*t291*Xoj
      t294 = t25*t141
      t295 = t290**2
      t296 = t295*t160
      t297 = t43*t4
      t298 = t297*as
      t301 = 2*t298*t55*bp
      t302 = -t271+t274+t276-t166-t279+t280+t287-t293+t294+t296+t301
      t303 = t284**2
      t304 = t303*t160
      t305 = t22*t141
      t307 = 2*t291*t72
      t308 = t269*t291
      t309 = 2*t308
      t311 = 2*t285*t68
      t312 = t277*t285
      t313 = 2*t312
      t315 = t297*as*Zoi
      t316 = 2*t315
      t317 = t297*t138
      t318 = 2*t317
      t319 = t59*t131
      t320 = 2*t319
      t322 = t59*ac*Zoi
      t323 = 2*t322
      t324 = t304+t305-t307+t309-t311+t313-t316+t318-t320+t323-t16
      t329 = t44*t7*ac*t4*as
      t330 = 2*t329
      t332 = 2*t291*Xoi
      t333 = -t18+t37+t42-t330+t332-t49-t54+t58-t195-t197-t67
      t335 = 2*t277*Yoi
      t337 = 2*t285*Yoj
      t338 = t70+t74+t335-t337-t82+t94+t95+t96+t97+t102+t103+t104
      t341 = sqrt(t302+t324+t333+t338)
      t343 = (t341+0.1E-98)**2
      t347 = t285*t119
      t348 = t291*t68
      t352 = qh**2
      t353 = t352*ec
      t357 = 2*t297*t160*t50*t28
      t361 = 2*t59*t141*t50*t31
      t363 = 2*t298*t164
      t364 = t133+t136-t357-t140-t361+t363+t142-t149+t152-t158-t271+t161
     #+t162+t274+t276
      t365 = t290*t160
      t367 = 2*t365*t155
      t368 = t21*t141
      t370 = 2*t368*t35
      t371 = t11*t141
      t373 = 2*t371*t40
      t374 = t284*t160
      t376 = 2*t374*t146
      t377 = -t279-t172-t367-t370-t373-t376+t180+t280+t287-t293+t183-t18
     #5+t294+t296+t304+t305
      t379 = t277*t198
      t380 = t269*t173
      t381 = t285*t181
      t382 = t291*t150
      t384 = t59*ac*t176
      t385 = t308+t312-t315+t317-t319+t322-t15-t17-t329-t379-t188+t190-t
     #380-t381-t382+t384
      t386 = t332-t49+t193+t200+t203+t335-t337+t205+t95+t96+t97+t102+t10
     #3+t104+t209-t212
      t389 = sqrt(t364+t377+2*t385+t386)
      t391 = (t389+0.1E-98)**2
      t395 = t223+t224+t230+t231+t232+t233-t235+t236+t237-t239-t240-t241
      t396 = t277*t228
      t397 = t285*t225
      t398 = t368*t40
      t399 = t365*t146
      t400 = -t371*t35
      t401 = -t374*t155
      t402 = t269*t198
      t403 = t291*t181
      t404 = t243+t245+t247+t249-t396-t397-t398-t399-t400-t401-t402-t403
      t408 = t133-t136+t357+t140-t361+t363+t142+t149+t152+t158-t271+t161
     #+t162+t274+t276
      t409 = -t279+t172+t367-t370-t373+t376+t180+t280+t287-t293+t183-t18
     #5+t294+t296+t304+t305
      t411 = t308+t312-t315+t317-t319+t322-t15-t17-t329+t379-t188-t190+t
     #380-t381-t382-t384
      t412 = t332-t49-t193-t200+t203+t335-t337+t205+t95+t96+t97+t102+t10
     #3+t104-t209-t212
      t415 = sqrt(t408+t409+2*t411+t412)
      t417 = (t415+0.1E-98)**2
      t421 = -t223-t224-t230-t231-t232+t233+t235+t236-t237+t239-t240-t24
     #1
      t422 = t243+t245+t247+t249+t396-t397-t398+t399-t400+t401+t402-t403
      t426 = -t271+t274+t276-t166-t279+t280-t287+t293+t294+t296-t301
      t427 = t304+t305+t307-t309+t311-t313+t316-t318-t320+t323-t16
      t429 = -t18+t37+t42+t330-t332-t49-t54+t58-t195-t197-t67
      t430 = t70+t74+t335+t337-t82+t94+t95+t96+t97+t102+t103+t104
      t433 = sqrt(t426+t427+t429+t430)
      t435 = (t433+0.1E-98)**2
      t442 = t133+t136+t357-t140-t361-t363+t142-t149+t152-t158-t271+t161
     #+t162+t274+t276
      t443 = -t279-t172+t367-t370-t373+t376+t180+t280-t287+t293+t183-t18
     #5+t294+t296+t304+t305
      t445 = -t308-t312+t315-t317-t319+t322-t15-t17+t329-t379-t188+t190-
     #t380+t381+t382+t384
      t446 = -t332-t49+t193+t200+t203+t335+t337+t205+t95+t96+t97+t102+t1
     #03+t104+t209-t212
      t449 = sqrt(t442+t443+2*t445+t446)
      t451 = (t449+0.1E-98)**2
      t455 = t243+t245+t247+t249-t396+t397-t398+t399-t400+t401-t402+t403
      t459 = t133-t136-t357+t140-t361-t363+t142+t149+t152+t158-t271+t161
     #+t162+t274+t276
      t460 = -t279+t172-t367-t370-t373-t376+t180+t280-t287+t293+t183-t18
     #5+t294+t296+t304+t305
      t462 = -t308-t312+t315-t317-t319+t322-t15-t17+t329+t379-t188-t190+
     #t380+t381+t382-t384
      t463 = -t332-t49-t193-t200+t203+t335+t337+t205+t95+t96+t97+t102+t1
     #03+t104-t209-t212
      t466 = sqrt(t459+t460+2*t462+t463)
      t468 = (t466+0.1E-98)**2
      t472 = t243+t245+t247+t249+t396+t397-t398-t399-t400-t401+t402+t403
      t477 = -t1*ec/t109/t107*(t114+t116-t118+t120+t121+t75*t35-t124-t83
     #*t40)-t130/t218/t216*t250-t130/t262/t260*t266-t130/t343/t341*(t114
     #+t116-t118+t120+t121-t124-t347-t348-t226-t227)-t353/t391/t389*(2*t
     #395+2*t404)/2-t353/t417/t415*(2*t421+2*t422)/2-t130/t435/t433*(t11
     #4+t116-t118+t120+t121-t124+t347+t348-t226-t227)-t353/t451/t449*(2*
     #t395+2*t455)/2-t353/t468/t466*(2*t421+2*t472)/2
c
      dVdppj4p=t477
c
      return
      end








