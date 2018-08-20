Program firstmd
    implicit none
    !************************************
    !*** parameters for md simulation ***
    !************************************
    real(8),allocatable::xprev(:),x(:),xnext(:) ! for last time step, current
                                                ! and next time step's position
    real(8),allocatable::v(:)                   ! velocity of atom
    real(8),allocatable::f(:)                   ! force of each atom                                            
	real(8),allocatable::lattice_pos(:)         ! position of lattice(1D)                                            
	real(8)::m             ! assume all the atoms have the same mass
	real(8)::dt            ! delta t
	real(8)::temperature   ! temperature
	real(8)::energy,totalenergy
	integer(8)::natoms     ! number of atoms

	natoms=10
	dt=1.0e-12
	temperature=100.d0
	allocate(v(natoms),lattice_pos(natoms),f(natoms))
	allocate(xprev(natoms),x(natoms),xnext(natoms))

	call CreateSample(natoms,lattice_pos)
	write(*,*) "***********************************"
	write(*,*) "*** The position of lattice is: ***"
	write(*,"(5(E14.5,1X))") lattice_pos
	write(*,*) "***********************************"
	write(*,*) "***********************************"
	write(*,"(' *** number of atoms=',I5)") natoms
	write(*,"(' *** temperature=',E14.5)") temperature
	write(*,"(' *** delt t     =',E14.5)") dt
	write(*,*) "***********************************"
	!natoms,lattice_pos,x,xprev,v,dt,temp
	call Init(natoms,lattice_pos,x,xprev,v,dt,temperature)
	call Force(natoms,x,f,energy)

	write(*,*) "***********************************"
	write(*,*) "*** The force of lattice is: ***"
	write(*,"(5(E14.5,1X))") f
	write(*,"(' *** energy is',E14.5)") energy
	write(*,*) "***********************************"

	call Integration(natoms,dt,x,xprev,f,energy,totalenergy,temperature)
	write(*,*) "***********************************"
	write(*,*) "*** The x position is: ***"
	write(*,"(5(E14.5,1X))") x
	write(*,"(' *** energy is',E14.5,1X,'total energy is:',E14.5)") energy,totalenergy
	write(*,*) "***********************************"


    deallocate(xprev,x,xnext,lattice_pos,v)
	stop
end Program firstmd