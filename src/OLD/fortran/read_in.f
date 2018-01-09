      subroutine read_in(nclust,nat,nmono,nseed,nbet_clust,nterm_aft_rep
     &                  ,nmaxgen,per_mut,enr_nxt_m,natoms,namo,filename
     &                  ,nfilename,opt_ga,ditr_chos,opt_pred,pot_chos
     &                  ,opt_mat,opt_mut,opt_ini,inputfile,ifilelen
     &                  ,min_chos,opt_min,opt_out)
      implicit real*8(a-h,o-z)
c      
      integer        nat(10)
      character*3    el
      character*20   pot_chos,opt_mat,opt_pred,opt_ini,opt_ga
     &              ,opt_mut,ditr_chos,inputfile,filein,min_chos
     &              ,opt_min,opt_out
      character*30   filename
c     
      do i=1,10,1
         nat(i)=0
      enddo
c- read input parameters
      read(1,*)nclust
      write(*,*)nclust
      read(1,*)nat(1),nat(2),nat(3),nat(4),nat(5)
      write(*,*)nat(1),nat(2),nat(3),nat(4),nat(5)
      read(1,*)nseed
      write(*,*)nseed
      read(1,*)nbet_clust
      write(*,*)nbet_clust
      read(1,*)nterm_aft_rep
      write(*,*)nterm_aft_rep
      read(1,*)nmaxgen
      write(*,*)nmaxgen
      read(1,*)per_mut
      write(*,*)per_mut
      read(1,*)enr_nxt_m
      write(*,*)enr_nxt_m
      read(1,'(a20)')filename
      write(*,*)filename
      read(1,*)nfilename
      write(*,*)nfilename
      read(1,'(a20)')opt_ga
      write(*,*)opt_ga
      read(1,'(a20)')ditr_chos
      write(*,*)ditr_chos
      read(1,'(a20)')opt_pred
      write(*,*)opt_pred
      read(1,'(a20)')pot_chos
      write(*,*)pot_chos
      read(1,'(a20)')opt_mat
      write(*,*)opt_mat
      read(1,'(a20)')opt_mut
      write(*,*)opt_mut
      read(1,'(a20)')opt_ini
      write(*,*)opt_ini
      read(1,'(a20)')inputfile
      write(*,*)inputfile
      read(1,*)ifilelen
      write(*,*)ifilelen
      read(1,'(a20)')min_chos
      write(*,*)min_chos
      read(1,'(a20)')opt_min
      write(*,*)opt_min
      read(1,'(a20)')opt_out
      write(*,*)opt_out

c- set values 
      nmono=0
      do i=1,10,1
         if(nat(i).ne.0)nmono=nmono+1
      enddo
c
      natoms=0
      do k=1,nmono,1
         natoms=natoms+nat(k)
      enddo
      namo=natoms/nmono
c
      return
      end
