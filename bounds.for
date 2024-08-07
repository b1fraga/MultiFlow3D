!##############################################################################     
        subroutine boundu
!	Note 1: prescribed inflow always comes from East
!	Note 2: wall functions apply in ALL directions/if node inside vicous layer->no slip
!	Note 3: if LSM/LSMbase are selected, inflow/outflow BCs only applied BELOW surface/lid
!##############################################################################
        use vars
        use multidata
	  		use mpi
	  		use imb
        implicit none
        integer :: i,j,k,ib,ly,im
        integer :: is,ie,js,je,ks,ke,ktop,kk,jj
	  integer :: strlen,dummy
	  double precision :: Fwallu,dyy,dzz,dxx,up
	  character (LEN=25) :: filename
	  character (LEN=29) :: fileSEM
	  character (LEN=3) :: dominio
	  character (LEN=4) :: domain
	  character (LEN=5) :: name_end

        if (PERIODIC) call exchange_bc(1,pl_ex)
        do ly=0,pl_ex
        do ib=1,nbp
           dxx=dom(ib)%dx
           dyy=dom(ib)%dy
           dzz=dom(ib)%dz
           is=dom(ib)%isu; ie=dom(ib)%ieu
           js=dom(ib)%jsu; je=dom(ib)%jeu
           ks=dom(ib)%ksu; ke=dom(ib)%keu

! Boundary Conditions for u - i.e. shear at boundaries and/or diriclet conditions
!..............................................................................
!=== West ===> ..   4=wall  ..    1=Inflow	..	7=read inflow
!..............................................................................
        if (dom(ib)%iprev.lt.0) then
           if ((dom(ib)%bc_west.eq.4).or.(dom(ib)%bc_west.ge.61)) then
              do k=ks-1,ke+1; do j=js-1,je+1
                 dom(ib)%u(is-1-ly,j,k)= 0.0	
              end do; end do

           else if (dom(ib)%bc_west.eq.1) then
              do k=ks-1,ke+1; do j=js-1,je+1
                if (L_LSM) then
			if (dom(ib)%phi(is,j,k) .ge. 0.0) then
                    dom(ib)%u(is-1-ly,j,k)= ubulk
			else
                    dom(ib)%u(is-1-ly,j,k)= 0.0
			end if
		    else									!correction
                    dom(ib)%u(is-1-ly,j,k)= ubulk
		    end if
              end do; end do

          else if (dom(ib)%bc_west.eq.7) then					!brunho2014 reading slices

        		write(name_end,'(I5)') ireadinlet
        		strlen=LEN(TRIM(ADJUSTL(name_end)))
       		name_end=REPEAT('0',(5-strlen))//
     &		TRIM(ADJUSTL(name_end)) ! e.g. "00001"				

			write(dominio,'(I3)') dom_id(ib)
		  	strlen=LEN(TRIM(ADJUSTL(dominio)))
		 	dominio=REPEAT('0',(3-strlen))//
     &		TRIM(ADJUSTL(dominio)) ! e.g. "001"

		    filename='Inlet/Inlet_'//dominio//'_'//name_end//'.dat'

		      open (unit=405, file=filename)
		      do k=ks-1,ke+1; do j=js-1,je+1
				read(405,*)dom(ib)%u(is-1-ly,j,k)
		      end do; end do	
		      close (405)	

          else if (dom(ib)%bc_west.eq.8) then					!Pablo 14/12/2015 reading SEM inlet
        		write(name_end,'(I6)') ireadinlet
        		strlen=LEN(TRIM(ADJUSTL(name_end)))
       		name_end=REPEAT('0',(6-strlen))//
     &			TRIM(ADJUSTL(name_end)) ! e.g. "000001"
			write(domain,'(I4)') dom_id(ib)
        		strlen=LEN(TRIM(ADJUSTL(domain)))
       		domain=REPEAT('0',(4-strlen))//
     &			TRIM(ADJUSTL(domain)) ! e.g. "0001"
		    fileSEM='inflow/Inlet_'//domain//'_'//name_end//'.dat'
            	open (unit=405, file=fileSEM)
			read(405,*)
			read(405,*)
             do kk=1,ke-ks+2 ; do jj=1,je-js+2 	!This is Changed!!!!
			k=kk+pl-1 ; j=jj+pl-1
	   	 	read(405,*)up,dummy,dummy
			
	  		if (UPROF_SEM.eq.12) then 		!1/7th power law inlet condition Pablo 7/12/2015 (WITH SEM)
	   		if (dom(ib)%yc(j).lt.((yen-yst)/2)) then
          		dom(ib)%u(is-1-ly,j,k) = ubulk*(1.0d0+1.0d0/7.0d0)
     &	   	*(DABS(2*dom(ib)%yc(j)/(yen-yst)))**(1./7.)
	   		else
          		dom(ib)%u(is-1-ly,j,k) =  ubulk*(1.0d0+1.0d0/7.0d0)
     &    	*(DABS(2*((yen-yst)-dom(ib)%yc(j))/(yen-yst)))**(1.d0/7.d0)
	   		endif
          dom(ib)%u(is-1-ly,j,k) = dom(ib)%u(is-1-ly,j,k)*(1.0d0+1./7.)
     &	   *(DABS(dom(ib)%zc(k)/(zen-zst)))**(1./7.) +up
	  		else 
          	dom(ib)%u(is-1-ly,j,k) = ubulk+up
	  		endif
          	end do ; end do		

	   else if (dom(ib)%bc_west.eq.12) then 		!1/7th power law inlet condition Pablo 7/12/2015 (No SEM)
           	do k=ks-1,ke+1; do j=js-1,je+1
	  	if (dom(ib)%yc(j).lt.((yen-yst)/2)) then
          	dom(ib)%u(is-1-ly,j,k) = ubulk*(1.0d0+1.0d0/7.0d0)
     &	   *(DABS(2*dom(ib)%yc(j)/(yen-yst)))**(1./7.)
	  	else
          	dom(ib)%u(is-1-ly,j,k) = ubulk*(1.0d0+1.0d0/7.0d0)
     &    	*(DABS(2*((yen-yst)-dom(ib)%yc(j))/(yen-yst)))**(1.d0/7.d0)
	  	endif
          	dom(ib)%u(is-1-ly,j,k) = dom(ib)%u(is-1-ly,j,k)*(1.0d0+1./7.)
     &	   *(DABS(dom(ib)%zc(k)/(zen-zst)))**(1./7.)
	     	enddo ; end do 

	   else if (dom(ib)%bc_west.eq.13) then 		!1/7th power law inlet condition Pablo 7/12/2015
            do k=ks-1,ke+1; do j=js-1,je+1
	  	if (dom(ib)%yc(j).lt.((yen-yst)/2)) then
          	dom(ib)%u(is-1-ly,j,k) = ubulk*(1.0d0+1.0d0/7.0d0)
     &	   *(DABS(2*dom(ib)%yc(j)/(yen-yst)))**(1./7.)
	  	else
          	dom(ib)%u(is-1-ly,j,k) = ubulk*(1.0d0+1.0d0/7.0d0)
     &   	 *(DABS(2*((yen-yst)-dom(ib)%yc(j))/(yen-yst)))**(1.d0/7.d0)
	  	endif
	     	enddo ; end do 
	   endif
        end if
!...............................................................................
!=== East ===> ..  4=wall  ..   2=Outflow
!...............................................................................
        if (dom(ib)%inext.lt.0) then
           if ((dom(ib)%bc_east.eq.4).or.(dom(ib)%bc_east.ge.61)) then
              do k=ks-1,ke+1; do j=js-1,je+1
                 dom(ib)%u(ie+1+ly,j,k)= 0.0	
              end do; end do

           else if (dom(ib)%bc_east.eq.2 .or.dom(ib)%bc_east.eq.21) then
              if (alfabc.eq.1)then
                do k=ks-1,ke+1; do j=js-1,je+1

                if (L_LSM) then

			if (dom(ib)%phi(ie,j,k) .ge. 0.0) then
                    if (dom(ib)%bc_east.eq.2) then
                       dom(ib)%u(ie+1+ly,j,k)=dom(ib)%u(ie,j,k)
                    else
                       dom(ib)%u(ie+1+ly,j,k)=dom(ib)%uoo(ie+1+ly,j,k) 
     & -dt*alfapr*(dom(ib)%uoo(ie+1+ly,j,k)-
     & dom(ib)%uoo(ie,j,k))/dom(ib)%dx
                    end if
			else
			  dom(ib)%u(ie+1+ly,j,k)=0.0
			end if

		    else		!no LSM

                    if (dom(ib)%bc_east.eq.2) then
                       dom(ib)%u(ie+1+ly,j,k)=dom(ib)%u(ie,j,k)
                    else
                       dom(ib)%u(ie+1+ly,j,k)=dom(ib)%uoo(ie+1+ly,j,k) 
     & -dt*alfapr*(dom(ib)%uoo(ie+1+ly,j,k)-
     & dom(ib)%uoo(ie,j,k))/dom(ib)%dx
                    end if

		    end if
                    end do; end do
		  endif
           end if
        end if
!...............................................................................
!=== South ===> ..   4=wall ..   3=Symmetry .. 6=Wall function
!...............................................................................
        if (dom(ib)%jprev.lt.0) then

           if (dom(ib)%bc_south.eq.4) then 
              do k=ks-1,ke+1; do i=is-1,ie+1
                 dom(ib)%u(i,js-1-ly,k)= -dom(ib)%u(i,js+ly,k)	
              end do; end do                      
           else if (dom(ib)%bc_south.eq.3) then
              do k=ks-1,ke+1; do i=is-1,ie+1
                 dom(ib)%u(i,js-1-ly,k)= dom(ib)%u(i,js+ly,k)
              end do; end do
	     else if (dom(ib)%bc_south.ge.61) then					!Wall functions Bruño2014
		if (ly.eq.0) then
	   		if (dom(ib)%bc_south.lt.63) 
     &			call log_law(3,ib)
			if (dom(ib)%bc_south.ge.63) 
     &			call wall_function(3,ib)
		      do k=ks-1,ke+1; do i=is-1,ie+1
				Fwallu = dom(ib)%tauws2(i,k)
     &					*dom(ib)%u(i,js,k)/dxx			!*Acell/Vcell
				Fwallu=SIGN(Fwallu,-dom(ib)%u(i,js,k))
				dom(ib)%u(i,js,k)=dom(ib)%u(i,js,k) 
     &						+Fwallu*alfapr*dt
		           	dom(ib)%u(i,js-1,k)= -dom(ib)%u(i,js,k)	
		      end do; end do
		else
		      do k=ks-1,ke+1; do i=is-1,ie+1
		           	dom(ib)%u(i,js-1-ly,k)= -dom(ib)%u(i,js+ly,k)	
		      end do; end do
		endif
	     endif
        end if
!.............................................................................
!=== North ===>  ..   4=wall ..   44=moving wall ..  3=Symmetry
!.............................................................................
        if (dom(ib)%jnext.lt.0) then
           if (dom(ib)%bc_north.eq.4) then
              do k=ks-1,ke+1; do i=is-1,ie+1
                 dom(ib)%u(i,je+1+ly,k)   = - dom(ib)%u(i,je-ly,k) 	
              end do; end do

           else if (dom(ib)%bc_north.eq.3) then
              do k=ks-1,ke+1; do i=is-1,ie+1
                 dom(ib)%u(i,je+1+ly,k)   =  dom(ib)%u(i,je-ly,k) 
              end do; end do
	     else if (dom(ib)%bc_north.ge.61) then					!Wall functions Bruño2014
		if (ly.eq.0) then
	   		if (dom(ib)%bc_north.lt.63) 
     &			call log_law(4,ib)
			if (dom(ib)%bc_north.ge.63) 
     &			call wall_function(4,ib)
		      do k=ks-1,ke+1; do i=is-1,ie+1
				Fwallu = dom(ib)%tauwn2(i,k)
     &					*dom(ib)%u(i,je,k)/dxx			!*Acell/Vcell
				Fwallu=sign(Fwallu,-dom(ib)%u(i,je,k))
				dom(ib)%u(i,je,k)=dom(ib)%u(i,je,k) 
     &						+Fwallu*alfapr*dt
		           	dom(ib)%u(i,je+1,k)= -dom(ib)%u(i,je,k)	
		      end do; end do
		else
		      do k=ks-1,ke+1; do i=is-1,ie+1
		           	dom(ib)%u(i,je+1+ly,k)= -dom(ib)%u(i,je-ly,k)	
		      end do; end do
		endif
	     endif
        end if
!...............................................................................
!=== Bottom ===> ..   4=wall ..   3=Symmetry
!...............................................................................
        if (dom(ib)%kprev.lt.0) then
          if (dom(ib)%bc_bottom.eq.4) then
              do j=js-1,je+1; do i=is-1,ie+1
                 dom(ib)%u(i,j,ks-1-ly)= -dom(ib)%u(i,j,ks+ly)	
              end do; end do

           else if (dom(ib)%bc_bottom.eq.1) then
              do j=js-1,je+1; do i=is-1,ie+1
                 dom(ib)%u(i,j,ks-1-ly)= 0.0
              end do; end do

           else if (dom(ib)%bc_bottom.eq.3) then
              do j=js-1,je+1; do i=is-1,ie+1
                 dom(ib)%u(i,j,ks-1-ly)= dom(ib)%u(i,j,ks+ly)
              end do; end do
	     else if (dom(ib)%bc_bottom.ge.61) then					!Wall functions Bruño2014
		if (ly.eq.0) then
	   		if (dom(ib)%bc_bottom.lt.63)
     &			call log_law(5,ib)
			if (dom(ib)%bc_bottom.ge.63) 
     &			call wall_function(5,ib)
		      do j=js-1,je+1; do i=is-1,ie+1
				Fwallu = dom(ib)%tauwb2(i,j)
     &					*dom(ib)%u(i,j,ks)/dxx			!*Acell/Vcell??
				Fwallu = sign(Fwallu,-dom(ib)%u(i,j,ks))

				dom(ib)%u(i,j,ks)=dom(ib)%u(i,j,ks) 
     &						+Fwallu*alfapr*dt

		           	dom(ib)%u(i,j,ks-1)= -dom(ib)%u(i,j,ks) 
		      end do; end do
		else
		      do j=js-1,je+1; do i=is-1,ie+1
		           	dom(ib)%u(i,j,ks-1-ly)= -dom(ib)%u(i,j,ks+ly) 
		      end do; end do
		endif
           end if
        end if
!.............................................................................
!=== Top ===>  ..   4=wall ..     3=Symmetry
!.............................................................................
	  if (L_LSMbase) then
        if (dom(ib)%z(1).le.length .and. dom(ib)%z(ke).ge.length) then		
	    k=1
          do while (dom(ib)%z(k).le.length)
	      k=k+1
	    end do
	    ktop=k-1
          do j=js-1,je+1; do i=is-1,ie+1
             dom(ib)%u(i,j,ktop) = dom(ib)%u(i,j,ktop-1)
             dom(ib)%u(i,j,ktop+1) = 0.0
          enddo; enddo
	  end if
	  end if

        if (dom(ib)%knext.lt.0) then
           if (dom(ib)%bc_top.eq.4) then
              do j=js-1,je+1; do i=is-1,ie+1
                 dom(ib)%u(i,j,ke+1+ly)   = - dom(ib)%u(i,j,ke-ly)	
              end do; end do
           else if (dom(ib)%bc_top.eq.3) then
             do j=js-1,je+1; do i=is-1,ie+1
                 dom(ib)%u(i,j,ke+1+ly)= dom(ib)%u(i,j,ke-ly)
             end do; end do
	     else if (dom(ib)%bc_top.ge.61) then					!Wall functions Bruño2014
		if (ly.eq.0) then
	   		if (dom(ib)%bc_top.lt.63) 
     &			call log_law(6,ib)
			if (dom(ib)%bc_top.ge.63) 
     &			call wall_function(6,ib)
		      do j=js-1,je+1; do i=is-1,ie+1
				Fwallu = dom(ib)%tauwt2(i,j)
     &					*dom(ib)%u(i,j,ke)/dxx			!*Acell/Vcell
				Fwallu = sign(Fwallu,-dom(ib)%u(i,j,ke))
				dom(ib)%u(i,j,ke)=dom(ib)%u(i,j,ke) 
     &						+Fwallu*alfapr*dt
                 		dom(ib)%u(i,j,ke+1) = -dom(ib)%u(i,j,ke)	
		      end do; end do
		else
		      do j=js-1,je+1; do i=is-1,ie+1
                 		dom(ib)%u(i,j,ke+1+ly) = -dom(ib)%u(i,j,ke)	
		      end do; end do
		endif
           end if
	  end if

        end do
        if(diff_sch.ne.3) call boundcoef(1)
        end do
        end subroutine
!#############################################################################
        subroutine boundv
!#############################################################################
        use vars
        use multidata
	  use mpi
	  use imb
        implicit none
        integer :: i,j,k,ib,ly
        integer :: is,ie,js,je,ks,ke,ktop,kk,jj
	  integer :: strlen
	  character (LEN=25) :: filename
	  character (LEN=3) :: dominio
	  character (LEN=29):: fileSEM
	  character (LEN=4) :: domain
	  character (LEN=5) :: name_end
	  double precision :: Fwallv,dxx,dzz,dummy,dyy

        if (PERIODIC) call exchange_bc(2,pl_ex)
 
        do ly=0,pl_ex
        do ib=1,nbp
	     dxx=dom(ib)%dx
	     dyy=dom(ib)%dy
	     dzz=dom(ib)%dz
           is=dom(ib)%isv; ie=dom(ib)%iev
           js=dom(ib)%jsv; je=dom(ib)%jev
           ks=dom(ib)%ksv; ke=dom(ib)%kev

! Boundary Conditions for V - i.e. shear at boundaries
!...............................................................................
!=== West ===>   4=wall   ..    1=Inflow
!...............................................................................
        if (dom(ib)%iprev.lt.0) then
           if (dom(ib)%bc_west.eq.4) then
              do k=ks-1,ke+1; do j=js-1,je+1
                 dom(ib)%v(is-1-ly,j,k)= -dom(ib)%v(is+ly,j,k)	
              end do; end do

           else if (dom(ib)%bc_west.eq.1.or.dom(ib)%bc_west.eq.12) then
              do k=ks-1,ke+1; do j=js-1,je+1
                 dom(ib)%v(is-1-ly,j,k)=0.0
              end do; end do
	     else if (dom(ib)%bc_west.ge.61) then					!Wall functions Bruño2014
		if (ly.eq.0) then
	   		if (dom(ib)%bc_west.lt.63) 
     &			call log_law(1,ib)
			if (dom(ib)%bc_west.ge.63) 
     &			call wall_function(1,ib)
		      do k=ks-1,ke+1; do j=js-1,je+1
				Fwallv = dom(ib)%tauww2(j,k)
     &					*dom(ib)%v(is,j,k)/dyy			!*Acell/Vcell
				Fwallv = sign(Fwallv,-dom(ib)%v(is,j,k))
				dom(ib)%v(is,j,k)=dom(ib)%v(is,j,k) 
     &						+Fwallv*alfapr*dt
                 		dom(ib)%v(is-1,j,k)= -dom(ib)%v(is,j,k)	
		      end do; end do
		else
		      do k=ks-1,ke+1; do j=js-1,je+1
                 		dom(ib)%v(is-1-ly,j,k)= -dom(ib)%v(is,j,k)	
		      end do; end do
		endif
          else if (dom(ib)%bc_west.eq.7) then					!brunho2014 reading slices

        		write(name_end,'(I5)') ireadinlet
        		strlen=LEN(TRIM(ADJUSTL(name_end)))
       		name_end=REPEAT('0',(5-strlen))//
     &		TRIM(ADJUSTL(name_end)) ! e.g. "00001"				

			write(dominio,'(I3)') dom_id(ib)
		  	strlen=LEN(TRIM(ADJUSTL(dominio)))
		 	dominio=REPEAT('0',(3-strlen))//
     &		TRIM(ADJUSTL(dominio)) ! e.g. "001"

		    filename='Inlet/Inlet_'//dominio//'_'//name_end//'.dat'
		      open (unit=405, file=filename)
		      do k=dom(ib)%ksu-1,dom(ib)%keu+1
			do j=dom(ib)%jsu-1,dom(ib)%jeu+1
				read(405,*)dummy
		      end do; end do	
		      do k=ks-1,ke+1; do j=js-1,je+1
				read(405,*)dom(ib)%v(is-1-ly,j,k)
		      end do; end do		
			close(405)	
	
          else if (dom(ib)%bc_west.eq.8) then					!Pablo 14/12/2015 reading SEM inlet
        		write(name_end,'(I5)') ireadinlet
        		strlen=LEN(TRIM(ADJUSTL(name_end)))
       		name_end=REPEAT('0',(5-strlen))//
     &			TRIM(ADJUSTL(name_end)) ! e.g. "00001"
			write(domain,'(I4)') dom_id(ib)
        		strlen=LEN(TRIM(ADJUSTL(domain)))
       		domain=REPEAT('0',(4-strlen))//
     &			TRIM(ADJUSTL(domain)) ! e.g. "0001"
		    fileSEM='inflow/Inlet_'//domain//'_'//name_end//'.dat'
            	open (unit=405, file=fileSEM)
			read(405,*)
			read(405,*)
              do kk=1,ke-ks+2 ; do jj=1,je-js+2 	!This is Changed!!!!
		   k=kk+pl-1 ; j=jj+pl-1
	   	   read(405,*)dummy,dom(ib)%v(is-1-ly,j,k),dummy
              end do; end do				
			close(405)
           end if
        end if
!...............................................................................
!=== East ===>   4=wall   ..    2=Outflow
!...............................................................................
        if (dom(ib)%inext.lt.0) then
           if (dom(ib)%bc_east.eq.4) then
              do k=ks-1,ke+1; do j=js-1,je+1
                 dom(ib)%v(ie+1+ly,j,k)= -dom(ib)%v(ie-ly,j,k)	
              end do; end do

           else if(dom(ib)%bc_east.eq.2 .or.dom(ib)%bc_east.eq.21) then
              if (alfabc.eq.1) then
                 do k=ks-1,ke+1; do j=js-1,je+1

                 if (L_LSM) then

			 if (dom(ib)%phi(ie,j,k) .ge. 0.0) then
                     if (dom(ib)%bc_east.eq.2) then
                       dom(ib)%v(ie+1+ly,j,k)= dom(ib)%v(ie,j,k)
                     else
                       dom(ib)%v(ie+1+ly,j,k)= dom(ib)%voo(ie+1+ly,j,k)
     & -dt*alfapr*(dom(ib)%voo(ie+1+ly,j,k)-
     & dom(ib)%voo(ie,j,k))/dom(ib)%dx
                     end if
			 else
			   dom(ib)%v(ie+1+ly,j,k)=0.0
			 endif

		     else									!no LSM

                     if (dom(ib)%bc_east.eq.2) then
                       dom(ib)%v(ie+1+ly,j,k)= dom(ib)%v(ie,j,k)
                     else
                       dom(ib)%v(ie+1+ly,j,k)= dom(ib)%voo(ie+1+ly,j,k)
     & -dt*alfapr*(dom(ib)%voo(ie+1+ly,j,k)-
     & dom(ib)%voo(ie,j,k))/dom(ib)%dx
                     end if

		     endif
                 end do; end do
              end if
	     else if (dom(ib)%bc_east.ge.61) then					!Wall functions Bruño2014
		if (ly.eq.0) then
	   		if (dom(ib)%bc_east.lt.63) 
     &			call log_law(2,ib)
			if (dom(ib)%bc_east.ge.63) 
     &			call wall_function(2,ib)
		      do k=ks-1,ke+1; do j=js-1,je+1
				Fwallv = dom(ib)%tauwe2(j,k)
     &					*dom(ib)%v(ie,j,k)/dyy			!*Acell/Vcell
				Fwallv = sign(Fwallv,-dom(ib)%v(ie,j,k))
				dom(ib)%v(ie,j,k)=dom(ib)%v(ie,j,k) 
     &						+Fwallv*alfapr*dt
                 		dom(ib)%v(ie+1,j,k)= -dom(ib)%v(ie,j,k)	
		      end do; end do
		else
		      do k=ks-1,ke+1; do j=js-1,je+1
                 		dom(ib)%v(ie+1+ly,j,k)= -dom(ib)%v(ie,j,k)	
		      end do; end do
		endif
           end if
        end if
!...............................................................................
!=== South ===>   4=wall   ..    3=Symmetry
!...............................................................................
        if (dom(ib)%jprev.lt.0) then 
           if (dom(ib)%bc_south.eq.3) then
              do k=ks-1,ke+1; do i=is-1,ie+1
                 dom(ib)%v(i,js-1-ly,k)=0.0
              end do; end do

	else if ((dom(ib)%bc_south.eq.4).or.(dom(ib)%bc_south.ge.61)) then
              do k=ks-1,ke+1; do i=is-1,ie+1
                 dom(ib)%v(i,js-1-ly,k)=0.0	
              end do; end do
           end if
        end if
!...............................................................................
!=== North ===>   4=wall   ..    3=Symmetry
!...............................................................................
        if (dom(ib)%jnext.lt.0) then
           if (dom(ib)%bc_north.eq.3) then
              do k=ks-1,ke+1; do i=is-1,ie+1
                 dom(ib)%v(i,je+1+ly,k)=0.0
              end do; end do

	else if ((dom(ib)%bc_north.eq.4).or.(dom(ib)%bc_north.ge.61)) then
              do k=ks-1,ke+1; do i=is-1,ie+1
                 dom(ib)%v(i,je+1+ly,k)=0.0	
              end do; end do
           end if
        end if
!...............................................................................
!=== Bottom ===> ..   4=wall ..   3=Symmetry
!...............................................................................
        if (dom(ib)%kprev.lt.0) then
          if (dom(ib)%bc_bottom.eq.4) then
              do j=js-1,je+1; do i=is-1,ie+1
                 dom(ib)%v(i,j,ks-1-ly)= -dom(ib)%v(i,j,ks+ly)	
              end do; end do

           else if (dom(ib)%bc_bottom.eq.1) then
              do j=js-1,je+1; do i=is-1,ie+1
                 dom(ib)%v(i,j,ks-1-ly)= 0.0
              end do; end do

           else if (dom(ib)%bc_bottom.eq.3) then
              do j=js-1,je+1; do i=is-1,ie+1
                 dom(ib)%v(i,j,ks-1-ly)= dom(ib)%v(i,j,ks+ly)
              end do; end do
	     else if (dom(ib)%bc_bottom.ge.61) then				!Wall functions Bruño2014
		if (ly.eq.0) then
	   		if (dom(ib)%bc_bottom.lt.63)
     &     				call log_law(5,ib)
			if (dom(ib)%bc_bottom.ge.63) 
     &			call wall_function(5,ib)
		      do j=js-1,je+1; do i=is-1,ie+1
				Fwallv = dom(ib)%tauwb2(i,j)
     &					*dom(ib)%v(i,j,ks)/dyy			!*Acell/Vcell
				Fwallv = sign(Fwallv,-dom(ib)%v(i,j,ks))
				dom(ib)%v(i,j,ks)=dom(ib)%v(i,j,ks) 
     &						+Fwallv*alfapr*dt
                 		dom(ib)%v(i,j,ks-1)= -dom(ib)%v(i,j,ks)	
		      end do; end do
		else
		      do j=js-1,je+1; do i=is-1,ie+1
                 		dom(ib)%v(i,j,ks-1-ly)= -dom(ib)%v(i,j,ks)	
		      end do; end do
		endif
           end if
        end if
!...............................................................................
!=== Top ===>   4=wall   ..    3=Symmetry
!...............................................................................
	  if (L_LSMbase) then
        if (dom(ib)%z(1).le.length .and. dom(ib)%z(ke).ge.length) then		
	    k=1
          do while (dom(ib)%z(k).le.length)
	      k=k+1
	    end do
	    ktop=k-1
          do j=js-1,je+1; do i=is-1,ie+1
            dom(ib)%v(i,j,ktop) = dom(ib)%v(i,j,ktop-1)
            dom(ib)%v(i,j,ktop+1) = 0.0
	    enddo; enddo
	  end if
	  end if

        if (dom(ib)%knext.lt.0) then
           if (dom(ib)%bc_top.eq.4) then
              do j=js-1,je+1; do i=is-1,ie+1
                 dom(ib)%v(i,j,ke+1+ly)= -dom(ib)%v(i,j,ke-ly)	
              end do; end do
           else if (dom(ib)%bc_top.eq.3) then
             do j=js-1,je+1; do i=is-1,ie+1
                 dom(ib)%v(i,j,ke+1+ly)= dom(ib)%v(i,j,ke-ly)
             end do; end do
	     else if (dom(ib)%bc_top.ge.61) then					!Wall functions Bruño2014
		if (ly.eq.0) then
	   		if (dom(ib)%bc_top.lt.63) 
     &			call log_law(6,ib)
			if (dom(ib)%bc_top.ge.63) 
     &			call wall_function(6,ib)

		      do j=js-1,je+1; do i=is-1,ie+1
				Fwallv = dom(ib)%tauwt2(i,j)
     &					*dom(ib)%v(i,j,ke)/dyy			!*Acell/Vcell
				Fwallv = sign(Fwallv,-dom(ib)%v(i,j,ke))
				dom(ib)%v(i,j,ke)=dom(ib)%v(i,j,ke) 
     &						+Fwallv*alfapr*dt
                 		dom(ib)%v(i,j,ke+1)= -dom(ib)%v(i,j,ke)	
		      end do; end do
		else
		      do j=js-1,je+1; do i=is-1,ie+1
                 		dom(ib)%v(i,j,ke+1+ly)= -dom(ib)%v(i,j,ke)	
		      end do; end do
		endif
           end if
	  end if

        end do
        if(diff_sch.ne.3) call boundcoef(2)
        end do
        end subroutine
!#############################################################################
        subroutine boundw
!#############################################################################
        use vars
        use multidata
	  use mpi
	  use imb
        implicit none
        integer :: i,j,k,ib,ly
        integer :: is,ie,js,je,ks,ke,ktop,jj,kk
	  integer :: strlen
	  character (LEN=25) :: filename
	  character (LEN=3) :: dominio
	  character (LEN=29):: fileSEM
	  character (LEN=4) :: domain
	  character (LEN=5) :: name_end
	  double precision :: Fwallw,dxx,dyy,dummy,dzz

        if (PERIODIC) call exchange_bc(3,pl_ex)

        do ly=0,pl_ex
        do ib=1,nbp
	     dxx=dom(ib)%dx
	     dyy=dom(ib)%dy
	     dzz=dom(ib)%dz
           is=dom(ib)%isw; ie=dom(ib)%iew
           js=dom(ib)%jsw; je=dom(ib)%jew
           ks=dom(ib)%ksw; ke=dom(ib)%kew

! Boundary Conditions for W - i.e. shear at boundaries
!...............................................................................
!=== West ===>   4=wall   ..    1=Inflow
!...............................................................................
        if (dom(ib)%iprev.lt.0) then
           if (dom(ib)%bc_west.eq.4) then
              do k=ks-1,ke+1; do j=js-1,je+1
                 dom(ib)%w(is-1-ly,j,k)= -dom(ib)%w(is+ly,j,k)	
              end do; end do

           else if (dom(ib)%bc_west.eq.1.or.dom(ib)%bc_west.eq.12) then
              do k=ks-1,ke+1; do j=js-1,je+1
                 dom(ib)%w(is-1-ly,j,k)= 0.0
              end do; end do
	     else if (dom(ib)%bc_west.ge.61) then					!Wall functions Bruño2014
		if (ly.eq.0) then
	   		if (dom(ib)%bc_west.lt.63) 
     &			call log_law(1,ib)
			if (dom(ib)%bc_west.ge.63) 
     &			call wall_function(1,ib)

		      do k=ks-1,ke+1; do j=js-1,je+1
				Fwallw = dom(ib)%tauww2(j,k)
     &					*dom(ib)%w(is,j,k)/dzz			!*Acell/Vcell
				Fwallw = sign(Fwallw,-dom(ib)%w(is,j,k))
				dom(ib)%w(is,j,k)=dom(ib)%w(is,j,k) 
     &						+Fwallw*alfapr*dt
                 		dom(ib)%w(is-1,j,k)= -dom(ib)%w(is,j,k)	
			enddo; enddo
		else
		      do k=ks-1,ke+1; do j=js-1,je+1
                 		dom(ib)%w(is-1-ly,j,k)= -dom(ib)%w(is,j,k)	
		      end do; end do
		endif
          else if (dom(ib)%bc_west.eq.7) then					!brunho2014 reading slices

        		write(name_end,'(I5)') ireadinlet
        		strlen=LEN(TRIM(ADJUSTL(name_end)))
       		name_end=REPEAT('0',(5-strlen))//
     &		TRIM(ADJUSTL(name_end)) ! e.g. "00001"				

			write(dominio,'(I3)') dom_id(ib)
		  	strlen=LEN(TRIM(ADJUSTL(dominio)))
		 	dominio=REPEAT('0',(3-strlen))//
     &		TRIM(ADJUSTL(dominio)) ! e.g. "001"

		    filename='Inlet/Inlet_'//dominio//'_'//name_end//'.dat'
		      open (unit=405, file=filename)
		      do k=dom(ib)%ksu-1,dom(ib)%keu+1
			do j=dom(ib)%jsu-1,dom(ib)%jeu+1
				read(405,*)dummy
		      end do; end do	
		      do k=dom(ib)%ksv-1,dom(ib)%kev+1
			do j=dom(ib)%jsv-1,dom(ib)%jev+1
				read(405,*)dummy
		      end do; end do	
		      do k=ks-1,ke+1; do j=js-1,je+1
				read(405,*)dom(ib)%w(is-1-ly,j,k)
			end do; end do	
			close(405)			
				
          else if (dom(ib)%bc_west.eq.8) then					!Pablo 14/12/2015 reading SEM inlet
        		write(name_end,'(I5)') ireadinlet
        		strlen=LEN(TRIM(ADJUSTL(name_end)))
       		name_end=REPEAT('0',(5-strlen))//
     &			TRIM(ADJUSTL(name_end)) ! e.g. "00001"
			write(domain,'(I4)') dom_id(ib)
        		strlen=LEN(TRIM(ADJUSTL(domain)))
       		domain=REPEAT('0',(4-strlen))//
     &			TRIM(ADJUSTL(domain)) ! e.g. "0001"
		    fileSEM='inflow/Inlet_'//domain//'_'//name_end//'.dat'
            	open (unit=405, file=fileSEM)
			read(405,*)
			read(405,*)
              do kk=1,ke-ks+2 ; do jj=1,je-js+2 	!This is Changed!!!!
			k=kk+pl-1 ; j=jj+pl-1
	         read(405,*)dummy,dummy,dom(ib)%w(is-1-ly,j,k)
              end do; end do			
			close(405)

           end if
        end if
!...............................................................................
!=== East ===>   4=wall   ..    2=Outflow
!...............................................................................
        if (dom(ib)%inext.lt.0) then
           if (dom(ib)%bc_east.eq.4) then
              do k=ks-1,ke+1; do j=js-1,je+1
                 dom(ib)%w(ie+1+ly,j,k)= -dom(ib)%w(ie-ly,j,k)	
              end do; end do

           elseif (dom(ib)%bc_east.eq.2 .or.dom(ib)%bc_east.eq.21) then
              if (alfabc.eq.1) then
                 do k=ks-1,ke+1; do j=js-1,je+1

                   if (L_LSM) then

			   if (dom(ib)%phi(ie,j,k) .ge. 0.0) then
                       if (dom(ib)%bc_east.eq.2) then
                       dom(ib)%w(ie+1+ly,j,k)= dom(ib)%w(ie,j,k)
                       else
                       dom(ib)%w(ie+1+ly,j,k)= dom(ib)%woo(ie+1+ly,j,k)
     & -dt*alfapr*(dom(ib)%woo(ie+1+ly,j,k)-
     & dom(ib)%woo(ie,j,k))/dom(ib)%dx
                       end if
			   else
                       dom(ib)%w(ie+1+ly,j,k)=0.0 
			   endif

			 else	!no LSM
                       if (dom(ib)%bc_east.eq.2) then
                       dom(ib)%w(ie+1+ly,j,k)= dom(ib)%w(ie,j,k)
                       else
                       dom(ib)%w(ie+1+ly,j,k)= dom(ib)%woo(ie+1+ly,j,k)
     & -dt*alfapr*(dom(ib)%woo(ie+1+ly,j,k)-
     & dom(ib)%woo(ie,j,k))/dom(ib)%dx
                       end if
			 endif

                 end do; end do
              end if
	     else if (dom(ib)%bc_east.ge.61) then					!Wall functions Bruño2014
		if (ly.eq.0) then
	   		if (dom(ib)%bc_east.lt.63) 
     &			call log_law(2,ib)
			if (dom(ib)%bc_east.ge.63) 
     &			call wall_function(2,ib)

		      do k=ks-1,ke+1; do j=js-1,je+1
				Fwallw = dom(ib)%tauwe2(j,k)
     &					*dom(ib)%w(ie,j,k)/dzz			!*Acell/Vcell
				Fwallw = sign(Fwallw,-dom(ib)%w(ie,j,k))
				dom(ib)%w(ie,j,k)=dom(ib)%w(ie,j,k) 
     &						+Fwallw*alfapr*dt
                 		dom(ib)%w(ie+1,j,k)= -dom(ib)%w(ie,j,k)	
		      end do; end do
		else
		      do k=ks-1,ke+1; do j=js-1,je+1
                 		dom(ib)%w(ie+1+ly,j,k)= -dom(ib)%w(ie,j,k)	
		      end do; end do
		endif
           end if
        end if
!...............................................................................
!=== South ===> ..   4=wall ..   3=Symmetry
!...............................................................................
        if (dom(ib)%jprev.lt.0) then
           if (dom(ib)%bc_south.eq.4) then
              do k=ks-1,ke+1; do i=is-1,ie+1
                 dom(ib)%w(i,js-1-ly,k)= -dom(ib)%w(i,js+ly,k)	
              end do; end do

           else if (dom(ib)%bc_south.eq.3) then
              do k=ks-1,ke+1; do i=is-1,ie+1
                 dom(ib)%w(i,js-1-ly,k)= dom(ib)%w(i,js+ly,k)
              end do; end do
	     else if (dom(ib)%bc_south.ge.61) then					!Wall functions Bruño2014
		if (ly.eq.0) then
	   		if (dom(ib)%bc_south.lt.63) 
     &			call log_law(3,ib)
			if (dom(ib)%bc_south.ge.63) 
     &			call wall_function(3,ib)
		      do k=ks-1,ke+1; do i=is-1,ie+1
				Fwallw = dom(ib)%tauws2(i,k)
     &					*dom(ib)%w(i,js,k)/dzz			!*Acell/Vcell
				Fwallw = sign(Fwallw,-dom(ib)%w(i,js,k))
				dom(ib)%w(i,js,k)=dom(ib)%w(i,js,k) 
     &						+Fwallw*alfapr*dt
		           	dom(ib)%w(i,js-1,k)= -dom(ib)%w(i,js,k)	
		      end do; end do
		else
		      do k=ks-1,ke+1; do i=is-1,ie+1
		           	dom(ib)%w(i,js-1-ly,k)= -dom(ib)%w(i,js,k)	
		      end do; end do
		endif
           end if
        end if
!.............................................................................
!=== North ===>  ..   4=wall ..   44=moving wall ..  3=Symmetry
!.............................................................................
        if (dom(ib)%jnext.lt.0) then
           if (dom(ib)%bc_north.eq.4) then
              do k=ks-1,ke+1; do i=is-1,ie+1
                 dom(ib)%w(i,je+1+ly,k)   = - dom(ib)%w(i,je-ly,k)	 
              end do; end do
           else if (dom(ib)%bc_north.eq.3) then
              do k=ks-1,ke+1; do i=is-1,ie+1
                 dom(ib)%w(i,je+1+ly,k)   =  dom(ib)%w(i,je-ly,k) 
              end do; end do
	     else if (dom(ib)%bc_north.ge.61) then					!Wall functions Bruño2014
		if (ly.eq.0) then
	   		if (dom(ib)%bc_north.lt.63) 
     &			call log_law(4,ib)
			if (dom(ib)%bc_north.ge.63) 
     &			call wall_function(4,ib)
		      do k=ks-1,ke+1; do i=is-1,ie+1
				Fwallw = dom(ib)%tauwn2(i,k)
     &					*dom(ib)%w(i,je,k)/dzz			!*Acell/Vcell
				Fwallw = sign(Fwallw,-dom(ib)%w(i,je,k))
				dom(ib)%w(i,je,k)=dom(ib)%w(i,je,k) 
     &						+Fwallw*alfapr*dt
		           	dom(ib)%w(i,je+1,k)= -dom(ib)%w(i,je,k)	
		      end do; end do
		else
		      do k=ks-1,ke+1; do i=is-1,ie+1
				dom(ib)%w(i,je+1+ly,k)= -dom(ib)%w(i,je,k)	
		      end do; end do
		endif
           end if
        end if
!...............................................................................
!=== Bottom ===>   4=wall   ..    3=Symmetry
!...............................................................................
        if (dom(ib)%kprev.lt.0) then 
           if (dom(ib)%bc_bottom.eq.3) then
              do j=js-1,je+1; do i=is-1,ie+1
                 dom(ib)%w(i,j,ks-1-ly)=0.0
              end do; end do

		else if ((dom(ib)%bc_bottom.eq.4).or.
     &		   (dom(ib)%bc_bottom.ge.61)) then
              do j=js-1,je+1; do i=is-1,ie+1
                 dom(ib)%w(i,j,ks-1-ly)=0.0	
              end do; end do

           else if (dom(ib)%bc_bottom.eq.1) then
              do j=js-1,je+1; do i=is-1,ie+1
                 dom(ib)%w(i,j,ks-1-ly)=1.1
              end do; end do
           end if
        end if
!...............................................................................
!=== Top ===>   4=wall   ..    3=Symmetry
!...............................................................................
        if (L_LSMbase) then
        if (dom(ib)%z(1).le.length .and. dom(ib)%z(ke).ge.length) then
          k=1
          do while (dom(ib)%z(k).le.length)
            k=k+1
	    enddo
	    ktop=k-1
          do j=js-1,je+1; do i=is-1,ie+1
             dom(ib)%w(i,j,ktop) = dom(ib)%w(i,j,ktop-1)
	       dom(ib)%w(i,j,ktop+1) = 0.0
	    enddo; enddo
	  end if
	  end if
        if (dom(ib)%knext.lt.0) then
           if (dom(ib)%bc_top.eq.3) then
             do j=js-1,je+1; do i=is-1,ie+1
                 dom(ib)%w(i,j,ke+1+ly)=0.0
             end do; end do
           else if ((dom(ib)%bc_top.eq.4).or.
     &				(dom(ib)%bc_top.ge.61)) then
              do j=js-1,je+1; do i=is-1,ie+1
                 dom(ib)%w(i,j,ke+1+ly)=0.0	
              end do; end do
           end if
	  end if

        end do
        if(diff_sch.ne.3) call boundcoef(3)
        end do
        end subroutine
!##############################################################################
        subroutine boundcoef(op)
!##############################################################################
        use vars
        use multidata
        implicit none
        integer :: i,j,k,ib,op
        integer :: is,ie,js,je,ks,ke
        double precision d,dxx,dyy,dzz
        double precision, pointer, dimension(:,:,:) :: fi


        do ib=1,nbp
           dxx=dom(ib)%dx*dom(ib)%dx
           dyy=dom(ib)%dy*dom(ib)%dy
           dzz=dom(ib)%dz*dom(ib)%dz
           
        select case (op)
           Case (1) 
           is=dom(ib)%isu; ie=dom(ib)%ieu
           js=dom(ib)%jsu; je=dom(ib)%jeu
           ks=dom(ib)%ksu; ke=dom(ib)%keu
           fi => dom(ib)%u
           Case (2) 
           is=dom(ib)%isv; ie=dom(ib)%iev
           js=dom(ib)%jsv; je=dom(ib)%jev
           ks=dom(ib)%ksv; ke=dom(ib)%kev
           fi => dom(ib)%v
           Case (3) 
           is=dom(ib)%isw; ie=dom(ib)%iew
           js=dom(ib)%jsw; je=dom(ib)%jew
           ks=dom(ib)%ksw; ke=dom(ib)%kew
           fi => dom(ib)%w
        end select
	
! Boundary Conditions for u - i.e. shear at boundaries and/or dirichlet conditions
!..............................................................................
!=== West ===> ..   4=wall  ..    1=Inflow
!..............................................................................
        if (dom(ib)%iprev.lt.0) then
           do k=ks,ke; do j=js,je
           d=fac*dom(ib)%vis(is,j,k)/dxx
           dom(ib)%ap(is,j,k)=dom(ib)%ap(is,j,k)+d
           if(diff_sch.eq.1) then
              dom(ib)%su(is,j,k)=dom(ib)%su(is,j,k)+d*fi(is-1,j,k)
           else
              dom(ib)%su(is,j,k)=dom(ib)%su(is,j,k)+2.0*d*fi(is-1,j,k)
           end if
           end do; end do
        end if
!...............................................................................
!=== East ===> ..  4=wall  ..   2=Outflow
!...............................................................................
        if (dom(ib)%inext.lt.0) then
           do k=ks,ke; do j=js,je
           d=fac*dom(ib)%vis(ie,j,k)/dxx
           dom(ib)%ap(ie,j,k)=dom(ib)%ap(ie,j,k)+d
           if(diff_sch.eq.1) then
              dom(ib)%su(ie,j,k)=dom(ib)%su(ie,j,k)+d*fi(ie+1,j,k)
           else
              dom(ib)%su(ie,j,k)=dom(ib)%su(ie,j,k)+2.0*d*fi(ie+1,j,k)
           end if
           end do; end do
        end if
!...............................................................................
!=== South ===> ..   4=wall ..   3=Symmetry ..	
!...............................................................................
        if (dom(ib)%jprev.lt.0) then         
           do k=ks,ke; do i=is,ie
           d=fac*dom(ib)%vis(i,js,k)/dyy
           dom(ib)%ap(i,js,k)=dom(ib)%ap(i,js,k)+d
           if(diff_sch.eq.1) then
              dom(ib)%su(i,js,k)=dom(ib)%su(i,js,k)+d*fi(i,js-1,k)
           else
              dom(ib)%su(i,js,k)=dom(ib)%su(i,js,k)+2.0*d*fi(i,js-1,k)
           end if
           end do; end do
        end if
!.............................................................................
!=== North ===>  ..   4=wall ..   44=moving wall ..  3=Symmetry
!.............................................................................
        if (dom(ib)%jnext.lt.0) then
           do k=ks,ke; do i=is,ie
           d=fac*dom(ib)%vis(i,je,k)/dyy
           dom(ib)%ap(i,je,k)=dom(ib)%ap(i,je,k)+d
           if(diff_sch.eq.1) then
              dom(ib)%su(i,je,k)=dom(ib)%su(i,je,k)+d*fi(i,je+1,k)
           else
              dom(ib)%su(i,je,k)=dom(ib)%su(i,je,k)+2.0*d*fi(i,je+1,k)
           end if
           end do; end do
        end if
!...............................................................................
!=== Bottom ===> ..   4=wall ..   3=Symmetry
!...............................................................................
        if (dom(ib)%kprev.lt.0) then
           do j=js,je; do i=is,ie
           d=fac*dom(ib)%vis(i,j,ks)/dzz
           dom(ib)%ap(i,j,ks)= dom(ib)%ap(i,j,ks) + d
           if(diff_sch.eq.1) then
              dom(ib)%su(i,j,ks)=dom(ib)%su(i,j,ks)+d*fi(i,j,ks-1)
           else
              dom(ib)%su(i,j,ks)=dom(ib)%su(i,j,ks)+2.0*d*fi(i,j,ks-1)
           end if 
           end do; end do
        end if
!.............................................................................
!=== Top ===>  ..   4=wall ..     3=Symmetry
!.............................................................................
        if (dom(ib)%knext.lt.0) then
           do j=js,je; do i=is,ie
           d=fac*dom(ib)%vis(i,j,ke)/dzz
           dom(ib)%ap(i,j,ke)=dom(ib)%ap(i,j,ke)+d 
           if(diff_sch.eq.1) then
              dom(ib)%su(i,j,ke)=dom(ib)%su(i,j,ke)+d*fi(i,j,ke+1)
           else
              dom(ib)%su(i,j,ke)=dom(ib)%su(i,j,ke)+2.0*d*fi(i,j,ke+1)
           end if  
           end do; end do
        end if

        end do

        end subroutine
!#############################################################################
        subroutine boundksgs_1eqn
!##############################################################################
        use vars
        use multidata
        implicit none
        integer :: i,j,k,ib,ly,ni,nj,nk
        integer :: is,ie,js,je,ks,ke

        do ly=0,2 ! due to SMART scheme
        do ib=1,nbp
           ni=dom(ib)%ttc_i; nj=dom(ib)%ttc_j; nk=dom(ib)%ttc_k
           is=dom(ib)%isp; ie=dom(ib)%iep
           js=dom(ib)%jsp; je=dom(ib)%jep
           ks=dom(ib)%ksp; ke=dom(ib)%kep
	
! Boundary Conditions for ksgs
!..............................................................................
!=== West ===> 
!..............................................................................
        if (dom(ib)%iprev.lt.0) then
           if (dom(ib)%bc_west.eq.4.or.dom(ib)%bc_west.ge.61) then
           do k=1,nk; do j=1,nj
              dom(ib)%ksgs(is-1-ly,j,k)= 0.0
           end do; end do
           else if (dom(ib)%bc_west.ne.5) then
           do k=1,nk; do j=1,nj
              dom(ib)%ksgs(is-1-ly,j,k)= dom(ib)%ksgs(is,j,k)
           end do; end do
           end if
        end if
!...............................................................................
!=== East ===> 
!...............................................................................
        if (dom(ib)%inext.lt.0) then
           if (dom(ib)%bc_east.eq.4.or.dom(ib)%bc_east.ge.61) then
           do k=1,nk; do j=1,nj
              dom(ib)%ksgs(ie+1+ly,j,k)= 0.0
           end do; end do
           else if (dom(ib)%bc_east.ne.5) then
           do k=1,nk; do j=1,nj
              dom(ib)%ksgs(ie+1+ly,j,k)= dom(ib)%ksgs(ie,j,k)
           end do; end do
           end if
        end if
!...............................................................................
!=== South ===>
!...............................................................................
        if (dom(ib)%jprev.lt.0) then
           if (dom(ib)%bc_south.eq.4.or.dom(ib)%bc_south.ge.61) then
           do k=1,nk; do i=1,ni
              dom(ib)%ksgs(i,js-1-ly,k)= 0.0	
           end do; end do
           else if (dom(ib)%bc_south.ne.5) then
           do k=1,nk; do i=1,ni
              dom(ib)%ksgs(i,js-1-ly,k)= dom(ib)%ksgs(i,js,k)
           end do; end do
           end if
        end if
!.............................................................................
!=== North ===>
!.............................................................................
        if (dom(ib)%jnext.lt.0) then
           if (dom(ib)%bc_north.eq.4.or.dom(ib)%bc_north.ge.61) then
           do k=1,nk; do i=1,ni
              dom(ib)%ksgs(i,je+1+ly,k) = 0.0	
           end do; end do
           else if (dom(ib)%bc_north.ne.5) then
           do k=1,nk; do i=1,ni
              dom(ib)%ksgs(i,je+1+ly,k) = dom(ib)%ksgs(i,je,k)	
           end do; end do
           end if
        end if
!...............................................................................
!=== Bottom ===>
!...............................................................................
        if (dom(ib)%kprev.lt.0) then
           if (dom(ib)%bc_bottom.eq.4.or.dom(ib)%bc_bottom.ge.61) then
           do j=1,nj; do i=1,ni
              dom(ib)%ksgs(i,j,ks-1-ly)= 0.0
           end do; end do
           else if (dom(ib)%bc_bottom.ne.5) then
           do j=1,nj; do i=1,ni
              dom(ib)%ksgs(i,j,ks-1-ly)= dom(ib)%ksgs(i,j,ks)
           end do; end do
           end if
        end if
!.............................................................................
!=== Top ===>
!.............................................................................
        if (dom(ib)%knext.lt.0) then
           if (dom(ib)%bc_top.eq.4.or.dom(ib)%bc_top.ge.61) then
           do j=1,nj; do i=1,ni
              dom(ib)%ksgs(i,j,ke+1+ly) = 0.0	
           end do; end do
           else if (dom(ib)%bc_top.ne.5) then
           do j=1,nj; do i=1,ni
              dom(ib)%ksgs(i,j,ke+1+ly) = dom(ib)%ksgs(i,j,ke)
           end do; end do
           end if
        end if

!==============================================================================
        end do
        end do
        end subroutine boundksgs_1eqn
!#############################################################################
