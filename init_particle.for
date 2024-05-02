!======================================================================!
!                 LAGRANGIAN PARTICLE TRACKING                         !
!----------------------------------------------------------------------!
!                           Bruño Fraga                                ! 
!                      Cardiff Uni 2013-2017                           !
!                        Stanford Uni 2018                             !
!                    Uni of Birmingham 2019-2024                       !
!======================================================================!
C#############################################################
      SUBROUTINE INIT_PARTICLE
C#############################################################
      use multidata
      use mpi
      use vars   
      use vars_pt

      implicit none   

      integer :: l,ib,f,frac1,frac_end,m
      integer :: i,j,k,nfrac,ptnr,tsnr,np_restart
      double precision :: random_number_normal,sigma_rho
      double precision :: xp,yp,zp,uop,vop,wop,Dp,sigma,rho_p
      double precision :: Wx,Wy,Wz,random_number_uniform
      double precision :: mindis,dist,distance
      logical :: random

      if (myrank.eq.0) then
            write(6,*)'................................................'
            write(6,*)'      LAGRANGIAN PARTICLE TRACKING ON'
            write(6,*)'................................................'
      endif

      open(10,file='LPT.cin')             !first read with variables applicable to all fractions
      read(10,*)                          !header
      read(10,*) PSIcell,order
      DF=.false.
         
      read(10,*) nfrac                    !how many Lag fractions you want to calculate BF2023
      read(10,*) Lcol
      read(10,*) Lcolwall

      np=0
      np_restart=0

      if (LRESTART) then                  !there is no release before the time loop
            open(20,file='final_particle.dat')
            read(20,*) np_restart
            np=np_restart
      endif

      do f=1,nfrac
            read(10,*)                    !fraction header
            read(10,*) tsnr
            read(10,*) ptnr
            if (tsnr.eq.-1) then
                  np=np+ptnr              !now we calculate the total number of particles to allocate the variables
                  if (myrank.eq.0)write(202,*) 'Initialising Lagrangian'
     &   ,' field. Releasing',ptnr,'new particles in fraction',f
            endif

            read(10,*)                    !dp
            read(10,*)                    !rhop
            read(10,*)                    !release volume
            read(10,*)random
            if (random) then 
                  read(10,*)
            else
                  do l=1,ptnr
                        read(10,*)
                  enddo
            endif
      enddo

      close(10)

      allocate (xp_pt(np),yp_pt(np),zp_pt(np))
      allocate (uop_pt(np),vop_pt(np),wop_pt(np))
      allocate (xpold(np),ypold(np),zpold(np))
      allocate (uopold(np),vopold(np),wopold(np))
      allocate (dp_pt(np),dp_old(np),rho_pt(np))
      allocate (Fu(np),Fv(np),Fw(np),rhop_old(np))
      allocate (ptsinproc(nprocs))

      open(30,file='LPT.cin')             !reopen the file to read the details of every fraction
      read(30,*)                          !header
      read(30,*)                          !PSIcell/ball
      read(30,*)                          !nfrac
      read(30,*)                          !Lcol
      read(30,*)                          !Lcolwall

      IF (myrank.eq.0) then

      IF (np_restart.gt.0) then            !there were particles before
      do l=1,np_restart
            read(20,*) xp_pt(l),yp_pt(l),zp_pt(l)
            read(20,*) uop_pt(l),vop_pt(l),wop_pt(l)
            read(20,*) dp_pt(l),rho_pt(l)
            xpold(l)=xp_pt(l)
            ypold(l)=yp_pt(l)
            zpold(l)=zp_pt(l)
            uopold(l)=uop_pt(l)
            vopold(l)=vop_pt(l)
            wopold(l)=wop_pt(l)
            dp_old(l)=dp_pt(l)
            rhop_old(l)=rho_pt(l)
            Fu(l)=0 ; Fv(l)=0 ; Fw(l)=0
      end do
      
      close (20)
      
      do ib=1,nbp
      dom(ib)%uoo = dom(ib)%u
      dom(ib)%voo = dom(ib)%v
      dom(ib)%woo = dom(ib)%w
      enddo

      ELSE  !np_restart=0

      frac1=1+np_restart                                                !starting on top of whichever particles were there before
      do f=1,nfrac
            read(30,*)                                                  !fraction header
            read(30,*) tsnr                                             !if 0 then no continuous release
            read(30,*) ptnr 
            frac_end=frac1+ptnr-1
            read(30,*) Dp,sigma
            read(30,*) rho_p,sigma_rho
            read(30,*) Wx,Wy,Wz
            read(30,*) random
            if (random) read(30,*)xp,yp,zp,uop,vop,wop
            do l=frac1,frac_end
               if (tsnr.eq.-1) then                                      !only initial release if no continuous release
                  if (random) then                                            !location
                        mindis=1.1*Dp
                        dist=0
                  do while (dist.lt.mindis)                             !avoiding overlap
                  dist=mindis
                  xp_pt(l)=random_number_uniform(xp-0.5*Wx,xp+0.5*Wx)
                  yp_pt(l)=random_number_uniform(yp-0.5*Wy,yp+0.5*Wy)
                  zp_pt(l)=random_number_uniform(zp-0.5*Wz,zp+0.5*Wz)                      
                  do m=frac1,l
                  if (l.ne.m) then 
                        distance=sqrt((xp_pt(l)-xp_pt(m))**2+(yp_pt(l)
     &                  -yp_pt(m))**2+(zp_pt(l)-zp_pt(m))**2)
                  else
                        distance=1.d9
                  endif
                        dist=min(dist,distance)
                  enddo
                  enddo
                  uop_pt(l)=uop
                  vop_pt(l)=vop
                  wop_pt(l)=wop
                  else                                                  !not random
                        read(30,*)xp_pt(l),yp_pt(l),zp_pt(l),
     &                  uop_pt(l),vop_pt(l),wop_pt(l)
                  endif

            dp_pt(l)= random_number_normal(Dp,sigma)
            rho_pt(l)= random_number_normal(rho_p,sigma_rho)


            xpold(l)=xp_pt(l)
            ypold(l)=yp_pt(l)
            zpold(l)=zp_pt(l)
            uopold(l)=uop_pt(l)
            vopold(l)=vop_pt(l)
            wopold(l)=wop_pt(l)
            dp_old(l)=dp_pt(l)
            Fu(l)=0 ; Fv(l)=0 ; Fw(l)=0

            else  !tsnr
                  if (.not.random) read(30,*)
            endif !tsnr 

!            print*,tsnr,f,np,frac1,frac_end,xp_pt(l),yp_pt(l),zp_pt(l)
!     &,wop_pt(l),dp_pt(l),rho_pt(l)       

            enddo                                                       !loop in particles within frac

            if (tsnr.eq.-1) frac1=frac1+ptnr

      enddo                                                             !loop in fracs
      end if                                                            !RESTART

      ENDIF                                                             !myrank
      
      close (30)

      RETURN
      END SUBROUTINE

C **********************************************************************
      SUBROUTINE TECPLOT(num_output)
C **********************************************************************

        use multidata
        use mpi
        use vars   
        use vars_pt

            implicit none     

            integer strlen,i,j,k,ib,ni,nj,nk,ii,idfile,num_output
            integer is,ie,js,je,ks,ke
            character(LEN=20) filename
            character(LEN=4) b_str,c_str
            double precision u_cn,v_cn,w_cn,p_cn,T_cn!,S_cn,k_cn,eps_cn,vis_cn
            

      do ib=1,nbp

!     if (dom_id(ib).eq.12.or.dom_id(ib).eq.37.or.dom_id(ib).eq.62
!     &     .or.dom_id(ib).eq.87.or.dom_id(ib).eq.112) then

        idfile=600+dom_id(ib)

        write(b_str,'(I4)') num_output
        strlen=LEN(TRIM(ADJUSTL(b_str)))
        b_str=REPEAT('0',(4-strlen))//TRIM(ADJUSTL(b_str)) ! e.g. "001"
        write(c_str,'(I4)') dom_id(ib)
        strlen=LEN(TRIM(ADJUSTL(c_str)))
        c_str=REPEAT('0',(4-strlen))//TRIM(ADJUSTL(c_str)) ! e.g. "001"

      filename='tecout_'//b_str//'_'//c_str//'.dat'

      OPEN (UNIT=idfile, FILE=filename)

      WRITE (idfile,*) 'TITLE = ', '"Eulerian field"'
      WRITE (idfile,"(A)")'VARIABLES = "X","Y","Z","U","V","W","P","T"'!,"S",
!     &"k","eps","vis"'


        is=pl+1; ie=dom(ib)%ttc_i-pl
        js=pl+1; je=dom(ib)%ttc_j-pl
        ks=pl+1; ke=dom(ib)%ttc_k-pl
        ni=ie-(is-1)+1
        nj=je-(js-1)+1
        nk=ke-(ks-1)+1

      !WRITE(idfile,*)'ZONE T="','id:',dom_id(ib),'it:',ntime,'"'
      WRITE(idfile,*)'zone ','STRANDID=', 1, 'SOLUTIONTIME=', ctime
      WRITE(idfile,*)'I=',ni,', J=',nj,', K=',nk,'F=POINT'

            do k=ks-1,ke
            do j=js-1,je
            do i=is-1,ie

                 u_cn  =0.25*(dom(ib)%u(i,j,k)+
     &dom(ib)%u(i,j+1,k)+dom(ib)%u(i,j,k+1)+
     &dom(ib)%u(i,j+1,k+1))

                 v_cn  =0.25*(dom(ib)%v(i,j,k)+
     &dom(ib)%v(i+1,j,k)+dom(ib)%v(i,j,k+1)+
     &dom(ib)%v(i+1,j,k+1))

                 w_cn  =0.25*(dom(ib)%w(i,j,k)+
     &dom(ib)%w(i+1,j,k)+dom(ib)%w(i,j+1,k)+
     &dom(ib)%w(i+1,j+1,k)) 

                 p_cn  =0.125*(dom(ib)%p(i,j,k)+
     &dom(ib)%p(i+1,j,k)    +dom(ib)%p(i,j+1,k)+
     &dom(ib)%p(i+1,j+1,k)  +dom(ib)%p(i,j,k+1)+
     &dom(ib)%p(i+1,j,k+1)  +dom(ib)%p(i,j+1,k+1)+
     &dom(ib)%p(i+1,j+1,k+1))
!                 S_cn  =0.125*(dom(ib)%S(i,j,k)+
!     &dom(ib)%S(i+1,j,k)    +dom(ib)%S(i,j+1,k)+
!     &dom(ib)%S(i+1,j+1,k)  +dom(ib)%S(i,j,k+1)+
!     &dom(ib)%S(i+1,j,k+1)  +dom(ib)%S(i,j+1,k+1)+
!     &dom(ib)%S(i+1,j+1,k+1))
!                 k_cn  =0.125*(dom(ib)%ksgs(i,j,k)+
!     &dom(ib)%ksgs(i+1,j,k)    +dom(ib)%ksgs(i,j+1,k)+
!     &dom(ib)%ksgs(i+1,j+1,k)  +dom(ib)%ksgs(i,j,k+1)+
!     &dom(ib)%ksgs(i+1,j,k+1)  +dom(ib)%ksgs(i,j+1,k+1)+
!     &dom(ib)%ksgs(i+1,j+1,k+1))
!                 eps_cn  =0.125*(dom(ib)%eps(i,j,k)+
!     &dom(ib)%eps(i+1,j,k)    +dom(ib)%eps(i,j+1,k)+
!     &dom(ib)%eps(i+1,j+1,k)  +dom(ib)%eps(i,j,k+1)+
!     &dom(ib)%eps(i+1,j,k+1)  +dom(ib)%eps(i,j+1,k+1)+
!     &dom(ib)%eps(i+1,j+1,k+1))
!                 vis_cn  =0.125*(dom(ib)%vis(i,j,k)+
!     &dom(ib)%vis(i+1,j,k)    +dom(ib)%vis(i,j+1,k)+
!     &dom(ib)%vis(i+1,j+1,k)  +dom(ib)%vis(i,j,k+1)+
!     &dom(ib)%vis(i+1,j,k+1)  +dom(ib)%vis(i,j+1,k+1)+
!     &dom(ib)%vis(i+1,j+1,k+1))
                 T_cn  =0.125*(dom(ib)%T(i,j,k)+
     &dom(ib)%T(i+1,j,k)    +dom(ib)%T(i,j+1,k)+
     &dom(ib)%T(i+1,j+1,k)  +dom(ib)%T(i,j,k+1)+
     &dom(ib)%T(i+1,j,k+1)  +dom(ib)%T(i,j+1,k+1)+
     &dom(ib)%T(i+1,j+1,k+1))

      write (idfile,'(11e14.6)') dom(ib)%x(i),dom(ib)%y(j),dom(ib)%z(k)
     & ,u_cn,v_cn,w_cn,p_cn,T_cn!S_cn,k_cn,eps_cn,vis_cn

            enddo
            enddo
            enddo
!           write (90,*) dom(ib)%isp,dom(ib)%iep,
!     & dom(ib)%jsp,dom(ib)%jep,dom(ib)%ksp,dom(ib)%kep
!     endif
      end do



      close (idfile)

!   88 FORMAT (10F15.8)

      END SUBROUTINE


C **********************************************************************
      SUBROUTINE TECPARTICLE(num_output)
C **********************************************************************
C
  
        use multidata
        use mpi
        use vars   
        use vars_pt

      implicit none     

      integer l,strlen,ib,num_output
      character(LEN=80) filename,filename2
      character(LEN=4) b_str

        write(b_str,'(I4)') num_output
        strlen=LEN(TRIM(ADJUSTL(b_str)))
        b_str=REPEAT('0',(4-strlen))//TRIM(ADJUSTL(b_str)) ! e.g. "001"

        filename='tecout_'//b_str//'_pt.dat'

        OPEN (UNIT=95, FILE=TRIM(ADJUSTL(filename)))

      WRITE (95,*) 'TITLE = ', '"Lagrangian field"'
      WRITE (95,"(A)")'VARIABLES = "X","Y","Z","U<sub>Lag<\sub>","V<sub>
     &Lag<\sub>","W<sub>Lag<\sub>","D<sub>Lag<\sub>","rho<sub>Lag<\sub>"
     &'!,"F<sub>u","F<sub>v","F<sub>w"'
      WRITE(95,*)'zone ','STRANDID=', 2, 'SOLUTIONTIME=', ctime

        do l=1,np
            WRITE (95,*) xp_pt(l),yp_pt(l),zp_pt(l)
     &                  ,uop_pt(l),vop_pt(l),wop_pt(l)
     &                  ,dp_pt(l),rho_pt(l)!,Fu(l),Fv(l),Fw(l)
        end do
        close (95)

!   88 FORMAT (10F15.8)

      END SUBROUTINE

!=======================================================================
            double precision FUNCTION random_number_uniform(a,b)
!=======================================================================

            implicit none     

            double precision :: a,b,r


            CALL RANDOM_SEED
            CALL RANDOM_NUMBER(r)

            random_number_uniform=(b-a)*r+a

            return
            end function



