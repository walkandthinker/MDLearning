subroutine init(natoms,x,temp,dt,xm,v)
    implicit none
    integer(8)::natoms
    real(8)::x(natoms),xm(natoms),v(natoms),temp,dt
    !***********************************
    integer(8)::i
    real(8)::sumv,sumv2,fs,randval

    call RANDOM_SEED()
    sumv=0.d0
    sumv2=0.d0

    do i=1,natoms
        call RANDOM_NUMBER(randval)
        v(i)=(randval-0.5d0)
        sumv=sumv+v(i)
        sumv2=sumv2+v(i)**2
    enddo

    sumv=sumv/natoms
    sumv2=sumv2/natoms
    fs=dsqrt(3.d0*temp/sumv2)
    do i=1,natoms
        v(i)=(v(i)-sumv)*fs
        xm(i)=x(i)-v(i)*dt
    enddo
end subroutine init