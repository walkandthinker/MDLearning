subroutine CreateAtoms(atomnum,x,xmin,xmax)
    implicit none
    integer(8)::atomnum
    real(8)::x(atomnum)
    real(8)::xmin,xmax
    !********************************
    integer(8)::i
    real(8)::randval
    call RANDOM_SEED()

    do i=1,atomnum
        call RANDOM_NUMBER(randval)
        x(i)=xmin+(xmax-xmin)*randval
    enddo

end subroutine CreateAtoms