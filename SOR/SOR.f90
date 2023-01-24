!Numerical Solver for 2D Poisson's Equation using 5-pt finite difference
!Scheme and SOR Iteration

program poisson2d_SOR
use input_mod
implicit none

integer::i,j,resflag,counter
real*8::uk(0:N,0:N),ukk(0:N,0:N),f(0:N,0:N),res
real*8::e0,ek
real*8::uexact(0:N,0:N)
real*8::h,w
character(len=6)::uniti

if(N<10)then
    write(uniti,'(I1.1)') N
  else if(N>=10 .and. N<100)then
    write(uniti,'(I2.2)') N
  else if(N>=100 .and. N<1000)then
    write(uniti,'(I3.3)') N
  else if(N>=1000 .and. N<10000)then
    write(uniti,'(I4.3)') N
  end if

!open(unit=10,file=trim('GS_N=')//trim(uniti)//trim('.txt'))

open(unit=11,file=trim('exact_N=')//trim(uniti)//trim('.txt'))


call readFile()

!Cell Size
h=1./real(N)
!Initial Guess
uk=0.
ukk=0
!Source Term
f=1.
!Optimal Weight
w=2./(1.+sin(pi*h))

if(iter_option==1)then
 print*,'Using Jacobi iteration.'
else if(iter_option==2)then
 print*,'Using Gauss-Seidel iteration.'
else if(iter_option==3)then
 print*,'Using SOR iteration.'
end if

print*,'N=',N

counter=0
resflag=0

!Initial Error
call computeError(e0,uk)
print*,'e0=',e0

ek=e0
res=tol+1.

!Compute Solution
do while((ek/e0)>1.e-6) !(resflag .eq. 0)
 !print*,'Iteration#',counter+1
 res=0.
 do i=1,N-1
  do j=1,N-1
   if(iter_option==1)then
     ukk(i,j)=0.25*(uk(i+1,j)+uk(i-1,j)+uk(i,j+1)+uk(i,j-1)+h*h*f(i,j))  
   else if(iter_option==2)then
     ukk(i,j)=0.25*(uk(i+1,j)+ukk(i-1,j)+uk(i,j+1)+ukk(i,j-1)+h*h*f(i,j))  
   else if(iter_option==3)then
    ukk(i,j)=(1.-w)*uk(i,j)+w*0.25*(uk(i+1,j)+ukk(i-1,j)+uk(i,j+1)+ukk(i,j-1)+h*h*f(i,j)) 
   end if 
   res=res+(ukk(i,j)-uk(i,j))**2
  end do
 end do
   res=sqrt(res)
 uk=ukk
 ek=0.
 call computeError(ek,ukk)
 !print*,'k=',counter+1,' ek=',ek,' res=',res
 counter=counter+1
end do

print*,'Number of iterations=',counter
print*,'u(0.5,0.5)_Numerical=',ukk(N/2,N/2)
print*,'u(0.5,0.5)_Exact=',uexact(N/2,N/2)
print*,'ek=',ek
print*,'ek/e0=',ek/e0
!save output to file
!call output()

!close(unit=10)
close(unit=11)

contains


subroutine readFile()

do i=0,N
  do j=0,N
    read(unit=11,fmt=*) uexact(i,j)
  end do
end do

end subroutine readFile

subroutine output()
real*8::x,y

do i=0,N
  do j=0,N
   x=i*h
   y=j*h
   write(10,*) x,y,ukk(i,j)
  end do
end do

end subroutine output


subroutine computeError(errorNorm,u)
!Input Variables
real*8::errorNorm,u(0:N,0:N)

do i=0,N
 do j=0,N
  errorNorm=errorNorm+(u(i,j)-uexact(i,j))**2
 end do
end do
  errorNorm=sqrt(errorNorm)
end subroutine computeError


end program poisson2d_SOR
