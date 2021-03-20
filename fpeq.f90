module fun
implicit none
contains


double precision function k(x,u,t)
double precision::x,t,u,b
b=3
!input the functional form of force as function of position x, speed u, time
k=-x-b*u+1*sin(2*t/4)
end function k


end module fun
use fun
use mathlib




implicit none
double precision,allocatable::w(:,:),x(:),u(:),k0(:,:),k1(:,:),k2(:,:),k3(:,:),w0(:,:),w1(:,:),w2(:,:),w3(:,:),e(:,:),f(:),p
double precision::t,h,hu,ht,xf,xi,uf,ui,b,tf,q,le,wmax

integer::n,i,j,l,time,timef,er,stamp

er=100
le=6.0
xi=-le
xf=le
ui=-le
uf=le
h=0.1
tf=500
ht=0.01
timef=int(tf*100)
write(*,*) timef
b=1
q=1
stamp=1
open(1,file="fpf.dat")
n=(xf-xi)/h

allocate(k0(n,n),k1(n,n),k2(n,n),k3(n,n),w(n,n),w0(n,n),w1(n,n),w2(n,n),w3(n,n),x(n),u(n),e(n,n),f(n))

do i=1,n
x(i)=xi+dble(i)*h
u(i)=ui+dble(i)*h
end do

!initial w
w=0
k0=0
k1=0
k2=0
k3=0
w0=0
w1=0
w2=0
w3=0

do i=2,n-1
	do j=2,n-1
		w(i,j)=exp(-b*(x(i))**2/(2*q)-b*(u(j))**2/(2*q))
	end do
end do

wmax=maxval(w)



!call apax2('fpf.dat',-le,le,-le,le,stamp)


time=0
do while(time.le.timef)
t=dble(time)*ht

do i=2,n-1
do j=2,n-1
if(w(i,j).gt.wmax/er .or.w(i+1,j).gt.wmax/er.or.&
w(i-1,j).gt.wmax/er .or.w(i,j+1).gt.wmax/er .or.w(i,j-1).gt.wmax/er) then	
	k0(i,j)=-u(j)*(w(i+1,j)-w(i-1,j))/(2*h)-(k(x(i),u(j+1),t)*w(i,j+1)-&
k(x(i),u(j-1),t)*w(i,j-1))/(2*h)+q*(w(i,j+1)-2*w(i,j)+w(i,j-1))/(h*h)
	w1(i,j)=w(i,j)+ht*k0(i,j)/2
	!write(*,*) time,x(i),u(j),w(i,j),(w(i+1,j)-w(i-1,j))/(2*h)

end if

end do
end do




do i=2,n-1
do j=2,n-1
if(w(i,j).gt.wmax/er .or.w(i+1,j).gt.wmax/er .or.&
w(i-1,j).gt.wmax/er .or.w(i,j+1).gt.wmax/er .or.w(i,j-1).gt.wmax/er) then	
	k1(i,j)=-u(j)*(w1(i+1,j)-w1(i-1,j))/(2*h)-(k(x(i),u(j+1),t)*w1(i,j+1)-&
k(x(i),u(j-1),t)*w1(i,j-1))/(2*h)+q*(w1(i,j+1)-2*w1(i,j)+w1(i,j-1))/(h*h)
	w2(i,j)=w1(i,j)+ht*k1(i,j)/2

end if
end do
end do


do i=2,n-1
do j=2,n-1
if(w(i,j).gt.wmax/er .or.w(i+1,j).gt.wmax/er .or.&
w(i-1,j).gt.wmax/er .or.w(i,j+1).gt.wmax/er .or.w(i,j-1).gt.wmax/er) then	
	k2(i,j)=-u(j)*(w2(i+1,j)-w2(i-1,j))/(2*h)-(k(x(i),u(j+1),t)*w2(i,j+1)&
-k(x(i),u(j-1),t)*w2(i,j-1))/(2*h)+q*(w2(i,j+1)-2*w2(i,j)+w2(i,j-1))/(h*h)
	w3(i,j)=w2(i,j)+ht*k2(i,j)/2

end if
end do
end do


do i=2,n-1
do j=2,n-1
if(w(i,j).gt.wmax/er .or.w(i+1,j).gt.wmax/er .or.&
w(i-1,j).gt.wmax/er .or.w(i,j+1).gt.wmax/er .or.w(i,j-1).gt.wmax/er) then	
	k3(i,j)=-u(j)*(w3(i+1,j)-w3(i-1,j))/(2*h)-(k(x(i),u(j+1),t)*w3(i,j+1)-&
k(x(i),u(j-1),t)*w3(i,j-1))/(2*h)+q*(w3(i,j+1)-2*w3(i,j)+w3(i,j-1))/(h*h)
end if
end do
end do


      
do i=2,n-1
do j=2,n-1
if(w(i,j).gt.wmax/er .or.w(i+1,j).gt.wmax/er .or.&
w(i-1,j).gt.wmax/er .or.w(i,j+1).gt.wmax/er .or.w(i,j-1).gt.wmax/er) then	
	w(i,j)=w(i,j)+ht*(k0(i,j)+2*k1(i,j)+2*k2(i,j)+k3(i,j))/6
	
	end if
end do
end do
if(mod(time,100).eq.0 .and. time.gt.1) then
stamp=0
do i=1,n,2
do j=1,n,2
write(1,*) time/100,x(i),u(j),w(i,j)
end do 
write(1,*)
end do
!write(*,*) w(n-1,(n-1)/2),w(0,n-1)

end if
stamp=1

!if(mod(time,200).eq.0) then

!do i=1,n,2
!open(2,file="fpe.dat")
!p=0
!do j=1,n
!p=p+w(i,j)*(0.25*x(i)**4+0.5*u(j)**2)
!end do
!write(2,*) x(i),p
!end do
close(2)
!end if
!call denplota2(x,u,w)
time=time+1
if(mod(time,100).eq.0) write(*,*) time*ht
end do
write(*,*) (time-1)/50,storage_size(w)


deallocate(x,u,w,k0,k1,k2,k3,w0,w1,w2,w3,e,f)
close(1)
end
