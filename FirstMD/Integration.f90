subroutine Integration(natoms,dt,x,xprev,f,energy,totalenergy,temperature)
	implicit none
	integer(8)::natoms
	real(8)::dt,x(natoms),xprev(natoms),f(natoms)
	real(8)::energy,totalenergy,temperature
	!******************************************
	real(8)::sumv,sumv2,xx,vi
	integer(8)::i
	!************************
	sumv=0.d0;sumv2=0.d0
	do i=1,natoms
		xx=2.d0*x(i)-xprev(i)+dt*dt*f(i) ! taylor expansion for next step's position
		vi=(xx-xprev(i))/(2.d0*dt)
		sumv =sumv+vi
		sumv2=sumv2+vi*vi

		! update position 
		xprev(i)=x(i)
		x(i)=xx
	enddo

	temperature=sumv2/(3.d0*natoms) ! instantaneous temperature
	totalenergy=(energy+0.5*sumv2)/natoms
end subroutine Integration
