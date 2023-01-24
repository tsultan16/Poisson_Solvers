!Numerical Solver for 2D Poisson's Equation using 5-pt finite difference
!scheme with Conjugate-Gradient method
!Tanzid Sultan:3/17/19

program poisson2d_CG
implicit none

integer,parameter::N=64
real,parameter::tol=1.e-16

integer::i,j,counter
real*8::xk((N-1)*(N-1)),xkk((N-1)*(N-1)),f((N-1)*(N-1))
real*8::rk((N-1)*(N-1)),pk((N-1)*(N-1)),Apk((N-1)*(N-1))
real*8::A((N-1)*(N-1),(N-1)*(N-1))
real*8::alphak,betak
real*8::res,e0,ek
real*8::h
real*8::rr,pAp,rrk

character(len=6)::uniti
if(N<10)then
  write(uniti,'(I1.1)') N
else if(N>=10 .and. N<100)then
  write(uniti,'(I2.2)') N
else if(N>=100 .and. N<1000)then
  write(uniti,'(I3.3)') N
end if

open(unit=10,file=trim('CG_N=')//trim(uniti)//trim('.txt'))
open(unit=11,file=trim('exact_N=')//trim(uniti)//trim('.txt'))

!----------------------------------------------------------
!Initialization
!----------------------------------------------------------
h=1./real(N) !Cell Size
xk=0. !Initial Guess
xkk=0.
f=h*h !Source Term
rk=f !Residual
pk=rk !Search Direction

!Construct A matrix
call constructMatrix()

counter=0
res=tol+1.

!----------------------------------------------------------
!Iterations
!----------------------------------------------------------
do while(res>tol)
 rr=0.
 pAp=0.
 rrk=0.
 Apk=0.
 !Compute step size:alphak
 do i=1,(N-1)*(N-1)
  do j=1,(N-1)*(N-1)
    Apk(i)=Apk(i)+A(i,j)*pk(j) !Apk= [A]*[pk]
  end do
 end do

 do i=1,(N-1)*(N-1)
  rr=rr+rk(i)*rk(i)
  pAp=pAp+pk(i)*Apk(i) 
 end do
 alphak=rr/pAp 

 !Compute solution iterate and update residual
 do i=1,(N-1)*(N-1)
  xkk(i)=xk(i)+alphak*pk(i)
  rk(i)=rk(i)-alphak*Apk(i)
 end do

 !Compute betak
 do i=1,(N-1)*(N-1)
  rrk=rrk+rk(i)*rk(i)
 end do
 betak=rrk/rr
 
 !Update search direction
 do i=1,(N-1)*(N-1)
  pk(i)=rk(i)+betak*pk(i)
 end do
 
 xk=xkk
 
 !Compute euclidean norm of residual vector
 res=sqrt(rrk)
 
 counter=counter+1

end do

print*,'Number of iterations=',counter
print*,'Residual=',res

!save output to file
call output()
close(unit=10)


contains

subroutine constructMatrix()

integer::i,j,k,l,m,NN
real*8::D((N-1)*(N-1),(N-1)*(N-1)),Lo((N-1)*(N-1),(N-1)*(N-1))
real*8::Up((N-1)*(N-1),(N-1)*(N-1)),B(N-1,N-1),ID(N-1,N-1)

!construct matrix B and ID(identity matrix)
B=0.
ID=0.
do i=1,N-1 
 do j=1,N-1
  if(i==j)then
    B(i,j)=4.
    ID(i,j)=1.
  else if(i==j+1 .or. i==j-1)then
    B(i,j)=-1.
  end if
 end do
end do

!construct matrix D 
D=0.
!go to center of each diagonal block and fill up all the elements around it
do k=0,N-2
  l=(N/2)+k*(N-1) 
  do i=l-(N/2-1),l+(N/2-1)
   do j=l-(N/2-1),l+(N/2-1)
     D(i,j)=B(i-(l-(N/2-1))+1,j-(l-(N/2-1))+1)          
   end do
  end do
end do

!construct matrix L
Lo=0.
do k=0,N-3
  l=(N/2)+(k+1)*(N-1)
  m=(N/2)+k*(N-1) 
  do i=l-(N/2-1),l+(N/2-1)
   do j=m-(N/2-1),m+(N/2-1)
     Lo(i,j)=ID(i-(l-(N/2-1))+1,j-(m-(N/2-1))+1)          
   end do
  end do
end do

!construct matrix U
Up=0.
do k=0,N-3
  l=(N/2)+k*(N-1)
  m=(N/2)+(k+1)*(N-1) 
  do i=l-(N/2-1),l+(N/2-1)
   do j=m-(N/2-1),m+(N/2-1)
     Up(i,j)=ID(i-(l-(N/2-1))+1,j-(m-(N/2-1))+1)          
   end do
  end do
end do

!Finally, construct A
A=D-Lo-Up

end subroutine constructMatrix


subroutine output()
real*8::x,y

do i=0,N
  do j=0,N
   x=i*h
   y=j*h
   if(i==0 .or. j==0 .or. i==N .or. j==N) then
    write(10,*) x,y,0.
    write(11,*) 0._8
   else
    write(10,*) x,y,xkk(i+(j-1)*(N-1))
    write(11,*) xkk(i+(j-1)*(N-1))
   end if
  end do
end do

end subroutine output

end program poisson2d_CG
