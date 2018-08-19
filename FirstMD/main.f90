Program firstmd
	implicit none
	real(8),allocatable::x(:),xm(:),v(:)
	integer(8)::nAtoms
	real(8)::temp,dt

	nAtoms=10
	allocate(x(nAtoms),xm(nAtoms),v(nAtoms))
	call CreateAtoms(nAtoms,x,0.d0,1.d0)
	write(*,"('Number of atoms:',I6)") nAtoms
	write(*,*) "Atoms information:"
	write(*,*) x
	!natoms,x,temp,dt,xm,v

	temp=273.d0;dt=1.0e-12
	call Init(nAtoms,x,temp,dt,xm,v)

	stop
end Program firstmd
	