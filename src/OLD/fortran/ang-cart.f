      subroutine angcart(nclust,natoms,namo,pop_clus)
      implicit real*8(a-h,o-z)
      include "param.inc" 
      include "common.inc"
c     
      real*8         pop_clus(n_max_clus,n_max_atom)
     &              ,mat_aux(n_max_clus,n_max_atom)
c
      roh1=dioa1
      roh2=dioa2
      rom1=diob1
      rom2=diob2
      alp=alph1*acos(-1d0)/180d0
      alp2=alph2*acos(-1d0)/180d0
      do i=1,nclust,1
c         write(3,*)natoms
c         write(3,*)'teste'
         do j=1,namo,1
            xo=pop_clus(i,j)
            yo=pop_clus(i,j+namo)
            zo=pop_clus(i,j+2*namo)
            ti=pop_clus(i,j+3*namo)
            si=pop_clus(i,j+4*namo)
            ppi=pop_clus(i,j+5*namo)
            sl11i=cos(ppi)*cos(si)-cos(ti)*sin(si)*sin(ppi)
            sl21i=-sin(ppi)*cos(si)-cos(ti)*sin(si)*cos(ppi)
            sl31i=sin(ti)*sin(si)
            sl12i=cos(ppi)*sin(si)+cos(ti)*cos(si)*sin(ppi)
            sl22i=-sin(ppi)*sin(si)+cos(ti)*cos(si)*cos(ppi)
            sl32i=-sin(ti)*cos(si)
            sl13i=sin(ppi)*sin(ti)
            sl23i=cos(ppi)*sin(ti)
            sl33i=cos(ti)
            x1=sl11i*roh1*cos(alp/2d0)+sl12i*roh1*sin(alp/2d0)+xo
            y1=sl21i*roh1*cos(alp/2d0)+sl22i*roh1*sin(alp/2d0)+yo 
            z1=sl31i*roh1*cos(alp/2d0)+sl32i*roh1*sin(alp/2d0)+zo
            x2=sl11i*roh2*cos(alp/2d0)-sl12i*roh2*sin(alp/2d0)+xo
            y2=sl21i*roh2*cos(alp/2d0)-sl22i*roh2*sin(alp/2d0)+yo 
            z2=sl31i*roh2*cos(alp/2d0)-sl32i*roh2*sin(alp/2d0)+zo
            x3=-sl11i*rom1*cos(alp2/2d0)+sl13i*rom1*sin(alp2/2d0)+xo
            y3=-sl21i*rom1*cos(alp2/2d0)+sl23i*rom1*sin(alp2/2d0)+yo 
            z3=-sl31i*rom1*cos(alp2/2d0)+sl33i*rom1*sin(alp2/2d0)+zo
            x4=-sl11i*rom2*cos(alp2/2d0)-sl13i*rom2*sin(alp2/2d0)+xo
            y4=-sl21i*rom2*cos(alp2/2d0)-sl23i*rom2*sin(alp2/2d0)+yo 
            z4=-sl31i*rom2*cos(alp2/2d0)-sl33i*rom2*sin(alp2/2d0)+zo
            mat_aux(i,j)=xo
            mat_aux(i,j+namo)=x1
            mat_aux(i,j+2*namo)=x2
            mat_aux(i,j+natoms)=yo
            mat_aux(i,j+natoms+namo)=y1
            mat_aux(i,j+natoms+2*namo)=y2
            mat_aux(i,j+2*natoms)=zo
            mat_aux(i,j+2*natoms+namo)=z1
            mat_aux(i,j+2*natoms+2*namo)=z2
c
c            write(3,*)'O ',xo,yo,zo
c            write(3,*)'H ',x1,y1,z1
c            write(3,*)'H ',x2,y2,z2
c            write(3,*)'cg',x3,y3,z3
c            write(3,*)'cg',x4,y4,z4
c            read(*,*)
         enddo
      enddo
c
      call matr1_matr2(nclust,3*natoms,mat_aux
     &     ,pop_clus)
cxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
c      do i=1,3*nclust
c         write(44,*)'clust',i,'      energy'
c         do j=1,natoms
c            write(44,*)j,pop_clus(i,j),pop_clus(i,j+natoms)
c     &           ,pop_clus(i,j+2*natoms)
c            if(pot_chos(1:3).eq.'tip')then
c               write(44,*)pop_clus(i,j+3*namo),pop_clus(i,j+4*namo)
c     &              ,pop_clus(i,j+5*namo)
c            endif
c         enddo
c      enddo
c      close(44)
cxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx      
c 
      return
      end
