subroutine Init(natoms,lattice_pos,x,xprev,v,dt,temp)
    implicit none
    integer(8)::natoms
    real(8)::lattice_pos(natoms)
    real(8)::x(natoms),xprev(natoms) ! current and previous position
    real(8)::v(natoms)               ! velocity
    real(8)::dt,temp                 ! delta time and temperature
    !*********************************************	
    real(8)::sumv,sumv2,fs,val
    integer(8)::i
    !*********************************************
    !*** initializing all the atoms' position and velocity
    sumv =0.d0
    sumv2=0.d0
    call random_seed()
    do i=1,natoms
    	call random_number(val)
        x(i) =lattice_pos(i)
        v(i) =(val-0.5d0)
        sumv =sumv+v(i)
        sumv2=sumv2+v(i)**2
    enddo
    sumv =sumv/natoms
    sumv2=sumv2/natoms
    fs=dsqrt(3.d0*temp/sumv2)
    do i=1,natoms
    	v(i)=(v(i)-sumv)*fs
    	xprev(i)=x(i)-v(i)*dt ! from taylor expansion
    enddo
end subroutine Init