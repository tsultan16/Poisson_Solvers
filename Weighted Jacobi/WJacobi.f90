!Numerical Solver for 2D Poisson's Equation using 5-pt finite difference
!Scheme and Weighted Jacobi Iteration

program poisson2d_wjacobi
use input_mod
implicit none

integer::i,j,counter,counterw,wflag
real*8::uk(0:N,0:N),ukk(0:N,0:N),f(0:N,0:N),res
real*8::ukw(0:N,0:N),ukkw(0:N,0:N),resw
real*8::uexact(0:N,0:N)
real*8::h,w,w1,ek,ekw,e0
character(len=6)::uniti
real*8::alpha,alphamax
real*8::Nk


if(N<10)then
    write(uniti,'(I1.1)') N
  else if(N>=10 .and. N<100)then
    write(uniti,'(I2.2)') N
  else if(N>=100 .and. N<1000)then
    write(uniti,'(I3.3)') N
  else if(N>=1000 .and. N<10000)then
    write(uniti,'(I4.3)') N
end if
open(unit=11,file=trim('exact_N=')//trim(uniti)//trim('.txt'))

open(unit=12,file='alpha.txt')

call readFile()

!Cell Size
h=1./real(N)
!Source Term
f=1.
!Weight
w1=1.0
w=w1
alpha=0.1
alphamax=1.4

print*,'N=',N

!Initial Error
call computeError(e0,uk)
print*,'e0=',e0

do Nk=1,500

print*,'Nk=',Nk
print*,'w1,w2,=',w,alpha*w
print*,'alpha=',alpha


!Initial Guess
uk=0.
ukk=0.
ukw=0.
ukkw=0.


counter=0
counterw=0

ek=e0
ekw=e0

res=tol+1.
resw=tol+1.
wflag=0


!Compute Solution
do while((ek/e0)>etol .or. (ekw/e0)>etol)!(res>tol .or. resw>tol)
 !print*,'Iteration#',counter+1
 res=0.
 resw=0.
 do i=1,N-1
  do j=1,N-1
   ukk(i,j)=0.25*(uk(i+1,j)+uk(i-1,j)+uk(i,j+1)+uk(i,j-1)+h*h*f(i,j))   
   ukkw(i,j)=(1.-w)*ukw(i,j)+w*0.25*(ukw(i+1,j)+ukw(i-1,j)+ukw(i,j+1)+ukw(i,j-1)+h*h*f(i,j))
   !ukk(i,j)=0.25*(uk(i+1,j)+ukk(i-1,j)+uk(i,j+1)+ukk(i,j-1)+h*h*f(i,j))   
   !ukkw(i,j)=(1.-w)*ukw(i,j)+w*0.25*(ukw(i+1,j)+ukkw(i-1,j)+ukw(i,j+1)+ukkw(i,j-1)+h*h*f(i,j))


   res=res+(ukk(i,j)-uk(i,j))**2
   resw=resw+(ukkw(i,j)-ukw(i,j))**2
  end do
 end do

 res=sqrt(res)
 resw=sqrt(resw)  
 uk=ukk
 ukw=ukkw

 ek=0.
 ekw=0.
 call computeError(ek,ukk)
 call computeError(ekw,ukkw)
 !print*,'k,kw=',counter+1,counterw+1,' ek=',ek,' ekw=',ekw
 if((ek/e0)>etol)then!(res>tol)then
   counter=counter+1
 end if
 if((ekw/e0)>etol)then!(resw>tol)then
   counterw=counterw+1
 end if 

 if(wflag==0)then
   w=alpha*w1
   wflag=1
 else if(wflag==1)then
   w=w1
   wflag=0
 end if

end do

write(12,*) alpha,counterw,counter
alpha=alphamax*(Nk/500.)
!print*,'Convergence successful...'
!print*,'Total number of Jacobi iterations required=',counter
!print*,'Total number of weighted Jacobi Iterations required=',counterw

!print*,'Enter values of w1,alpha.'
!read(*,*) w,alpha

end do

close(unit=11)

contains


subroutine readFile()

do i=0,N
  do j=0,N
    read(unit=11,fmt=*) uexact(i,j)
  end do
end do

end subroutine readFile



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


end program poisson2d_wjacobi
