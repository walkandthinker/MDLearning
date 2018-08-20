subroutine Force(natoms,x,f,energy)
	implicit none
	integer(8)::natoms
	real(8)::f(natoms),x(natoms),energy
	!*****************************
	integer(8)::i,j
	real(8)::xr,box,r2,r2i,r6i,ff,ecut,rc
	!****************************************
	rc=2.0e-2 ! critical distance
	ecut=4.d0*(1.d0/rc**6)*(1.d0/rc**6-1.d0)! Lennard-Jones potential
	box=1.d0
	energy=0.d0
	f=0.d0

	do i=1,natoms
		do j=i+1,natoms
			xr=x(i)-x(j)
			xr=xr-box*nint(xr/box)
			r2=xr**2
			if (r2<=rc*rc) then
				r2i=1.d0/r2
				r6i=r2i**3
				ff=48.d0*r2i*r6i*(r6i-0.5d0)
				f(i)=f(i)+ff*xr
				f(j)=f(j)-ff*xr
				energy=energy+4.d0*r6i*(r6i-1.d0)-ecut
			endif
		enddo
	enddo
end subroutine Force