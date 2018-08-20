subroutine CreateSample(natoms,lattice_pos)
	implicit none
	integer(8)::natoms
	real(8)::lattice_pos(natoms)
	!*******************************
	!*** generate the position for 1d lattice
	!*** question: how about 2d and 3d cases,
	!***           I think it should have a special kind of 
	!***           random value generation(really random case for 10^5~10^7 atoms)
	!***           otherwise you will have trouble!!!
	!*******************************
	integer(8)::i
	real(8)::val

	call random_seed()
	do i=1,natoms
		call random_number(val)
		lattice_pos(i)=val
	enddo
end subroutine CreateSample