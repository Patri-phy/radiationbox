! 1. Parameter settings, variable type declaration

parameter (npart=10000,nx=20,ny=20,nz=20,ntend=12000)

! Geometry
integer a(nx,ny,nz), d(nx,ny,nz) 									!cavity (with window) and domain
logical condx1,condx2,condy1,condy2,condz1,condz2,cav
logical condx,condy,condz,cube,ncube
logical windx1,windx2,windy1,windy2,windz1,windz2,wind,nwind,intwind
real dx, dy, dz, prob1, prob2, al, al1, al2, th, dl

! Physical behavoir implementation
integer k,outp,cont0,cont01,cont02,cont2
integer seed
real start, finish
logical free
real mu
integer time,dt														!time dependency

! Data elaboration						 							!spacial and time evol of superf distr
! real hystx(ny,nz,ntend), hysty(nx,nz,ntend)
! real hhystx(ny,nz)
real surf(ny,nz,ntend)
integer elem
real ave(ntend),avey(ny,ntend),ave2(ntend)							!average calc
real evo(ntend,2)													!evolution of aver. (matrix for data printing)

! 1. variable initialization

pi = 3.1415926535;

!Geometry
al1 = 2; al2 = 2; 													! distances of windows edges from resp. y and z edges
al=5; th=1

! Physical behavoir implementation
dt=1; outp=0; seed = 7
escape = 0;
cont0=0;cont01=0;cont02=0;cont2=0

! Data elaboration
! hyst = 0;
surf = 0; a = 0; ave=0; mu=0

! 2. Geometry implementation with probability settings
dx = al/nx; dy = al/ny; dz = al/nz; dl = dx*.9;
prob1=.001*(1-exp(-dl)); prob2 = .004; omega0 = 1-prob2;
if (dx<.001) stop

!Window index
miniy=int(al1/dy)+1; maxiy=int((al-al1)/dy)+1;
miniz=int(al2/dz)+1; maxiz=int((al-al2)/dz)+1;

if(1.eq.2)then
	print *, "1.Initial data"
	print *,"Numb of part, nx, maxtime, elem time chang, omega0, gas scatt prob, absorb prob"
	print 1, npart,nx,ntend,dt,omega0,prob1,prob2
	print *,"dx=al/nx"
	print 2, dx
	print *,"Window indexes"
	print 3, miniy, maxiy, miniz, maxiz
	print *,"Edge length, Thickness, y-dist, z-dist, elem displ,"
	print 4, al, th, al1, al2, dl
	
	k=0
	do ix=1,nx
		condx1 = x.lt.th;
		if (condx1) then
			x = ix*dx-.5*dx;
			k=k+1;
		end if
	end do
	print *, "Num of element within a wall is "
	print 5,k
end if

!Geometric figure definition
do ix=1,nx
	do iy=1,ny
		do iz=1,nz
		x = ix*dx-.5*dx; y = iy*dy-.5*dy; z = iz*dz-.5*dz;
		
		condx1 = x.lt.th; condx2=(al-x).lt.th;
		condy1 = y.lt.th; condy2=(al-y).lt.th;
		condz1 = z.lt.th; condz2=(al-z).lt.th;
		cav = condx1.or.condx2.or.condy1.or.condy2.or.condz1.or.condz2;
		
		condx = x.lt.al; condy = y.lt.al; condz = z.lt.al;
		cube = x.lt.al.and.y.lt.al.and.z.lt.al; ncube=.not.cube

		windx1 = (al-th).lt.x ; windx2=x.lt.al;
		windy1 = al1.lt.y; windy2=y.lt.(al-al1);
		windz1 = al2.lt.z; windz2=z.lt.(al-al2);

		wind = windx1.and.windx2.and.windy1.and.windy2.and.windz1.and.windz2; nwind =.not.wind;
		
		! if(windx1.and.flag==0)then
			! flag=1; holdix=ix;
		! end if
		! intwind = wind.and.(ix==holdix)
		
		if(cube) d(ix,iy,iz)=1 								!domain
		if(wind) a(ix,iy,iz)=2 								!window
		if(cav.and.nwind) a(ix,iy,iz)=1 					!open cavity
		! if(intwind) s(ix,iy,iz)=1 						!inner layer of the window
		
		end do
	end do
end do

!Printing the 3-dim matrix
if(2.eq.3)then
	print *, "2.Cavity"
	do ix=1,nx
			print 6,((a(ix,iy,iz),iy=1,ny),iz=1,nz)
	end do
end if

! 3. Physical behavoir implementation

call srand(seed+2);
call cpu_time(start)
! ind1=0; ind2=0; ind=0
do i=1,npart   														!Particle cycle
	! ind1=ind1+1;
	x=al/2; y=al/2; z=al/2;									
	time = 0;					
	weight = 1;														!Escape factor calculation
	free = .true.;
	
	do while(free.and.time<ntend)   								!History cycle
		! ind2=ind2+1; ind = ind1+ind2-1
		time=time+dt;	if(time==ntend-1) outp=outp+1;
		
		mu = 2*rand()-1; phi = 2*pi*rand()
		st = sqrt(1-mu**2);

		deltaz = dl*mu;	deltax = dl*st*cos(phi); deltay = dl*st*sin(phi)
		xold=x;	yold=y;	zold=z;	x=x+deltax;	y=y+deltay;	z=z+deltaz	!Generating displ. and recording it
		ix = int(x/dx)+1; iy = int(y/dy)+1; iz = int(z/dz)+1;
		
		if(d(ix,iy,iz).ne.1)then									!Searching for bugs
			print *, "bug"
			stop
		end if
		
		if(a(ix,iy,iz).eq.0)then									! scattering with a gas atom
																														
			if(rand().lt.prob1)then
		
				mu = 2*rand()-1; phi = 2*pi*rand()
				st = sqrt(1-mu**2);
								
				deltaz = dl*mu;	deltax = dl*st*cos(phi); deltay = dl*st*sin(phi)
				x=x+deltax; y=y+deltay;	z=z+deltaz
				ix = int(x/dx)+1; iy = int(y/dy)+1; iz = int(z/dz)+1
				cont0=cont0+1
			end if
			
		end if
		
		if(a(ix,iy,iz).eq.1)then
			if(rand().lt.prob2)then									! absorption by an internal wall atom
				free=.false.;
				weight = weight*omega0;
				cont01=cont01+1
			else														! scattering with an internal wall atom
				flag=0;
				do while(flag.eq.0)
					x=xold;y=yold;z=zold

					mu = 2*rand()-1; phi = 2*pi*rand();
					st = sqrt(1-mu**2)
					
					deltaz = dl*mu;	deltax = dl*st*cos(phi); deltay = dl*st*sin(phi)
					x=x+deltax; y=y+deltay;	z=z+deltaz
					ix = int(x/dx)+1; iy = int(y/dy)+1; iz = int(z/dz)+1
					
				
					flag=flag+1
					if(a(ix,iy,iz).eq.1)flag=0
				end do
				cont02=cont02+1
			end if
		end if
		
		if(a(ix,iy,iz).eq.2)then									! crossing of the window
			surf(iy,iz,time)=surf(iy,iz,time)+weight	
			free=.false.
			cont2=cont2+1
		end if
		
		if(d(ix,iy,iz).ne.1)then									! searching for bugs
			print *, "bug"
			stop
		end if
		
	end do															! end of history
		
	! sel=25
	! if(i==sel*int(i/sel))then
		! call cpu_time(finish); print *, "                 ", i,finish
	! end if
	
	escape = escape + weight
end do  															! end of particle cycle

! call cpu_time(finish); print *,finish

escape=escape/npart

! hyst=hyst/npart/dx/dy/dz
surf=surf/npart/dy/dz


! hyst=hyst/npart/dx/dy/dz+.00001
! surf=surf/npart/dy/dz+.00001
if(3.eq.4)then
	print *, "3.Escape factor, number out of time and number of events"
	print *, "Escape factor: "
	print 7, escape
			
	print *, "Particles out of time: "
	print 8, outp
		    
	print *, "Num of events: gas scatt, wall absorb., wall scatt, wind cross"
	print 9, cont0,cont01,cont02,cont2
	
	if(npart==outp+cont01+cont2.and.outp.ne.0)then
		print *, "Simulation is correct because missing particles are those"
		print *, "whose time process is run out during history cycle."
	end if
	if(outp.eq.0) print *, "Simulation is correct."
end if

! 4. Data elaboration

!Average calculation
if(4.eq.5)then
	! print *, "4.Average values"
	do time=1,ntend
		do iy=miniy,maxiy
			do iz=miniz,maxiz
				avey(iy,time)=surf(iy,iz,time)+avey(iy,time)
			end do
			ave(time)=avey(iy,time)+ave(time)
		end do
	end do
	ave(time)=ave(time)/((1+maxiy-miniy)*(1+maxiz-miniz))

!finding the index of the last nonzero value of ave(time)
	flagg=1
	do time=1,ntend
		do j=ntend-time,ntend
			if(ave(j).ne.0.and.flagg.eq.1)then
				n=j; flagg=0;
			end if
		end do
	end do
	
	do j=1,n
		do i=1,2
			if(i.eq.2)then
				evo(j,i)=ave(j)
			else
				evo(j,i)=j
			end if
		end do
	end do
	
	!average print
	
	! print *, "Averages"
	! print 10, (ave(time),time=1,ntend)
	! print *, "Averages with cutoff"
	! do j=1,n
		! print 11, (evo(j,i),i=1,2)
	! end do
	
	elem=25; av=0;	ave2=0; num=ntend/elem;
	do k=0,num-1
		do j=1,n
			if(1+k*elem.le.j.and.j.le.(k+1)*elem)then
				av = av + ave(j+1);
			end if
		end do
		ave2(k+1)= av/elem;
		av=0
	end do
	! print *, "Time ev. of av. intensity distr. with ", elem,"times lower time resolution"
	! print *, "Now the num of elem is ", num
	
	
	! open(unit=1, file = 'ave21.txt' )
	! write(1,*)ave2
	
end if

if(4.eq.5)then
	print 20, (ave2(k),k=1,num)
end if

if(5.eq.6)then
	print *, "5.Surface distribution"
	! time = 3
	do time=1,100
		do iy=miniy,maxiy
			print 12,(surf(iy,iz,time),iz=miniz,maxiz)
		end do
	end do
end if

if(6.eq.7)then
	print *, "6.Data for format settings"
	print *, "Resolution, numer of window indexes"
	print 13, nx,1+maxiy-miniy
end if

! if(7.eq.7)then
	! print *, "7.Internal distribution (nx/2 plane)"
	! hhystx(iy,iz)=0
	! do time=1,ntend
		! hhystx(iy,iz)=hhystx(iy,iz)+hystx(iy,iz,time)
	! end do
	! do iy=1,ny
		! print 14,(hhystx(iy,iz),iz=1,nz)
	! end do
! end if

! if(8.eq.8)then
	! print *, "7.Internal distribution (ny/2 plane)"
	! time=15
	! do iy=1,ny
		! print 14,(hysty(ix,iz,time),iz=1,nz)
	! end do
! end if
14 format (20(F10.6)," ")

! 1.Initial data
1 format (I13,I4,I9,I17,F8.3,F14.10,F15.12)
2 format (F15.13)
3 format (I6,I6,I6,I6)
4 format (F12.8,F11.8,F8.5,F8.5,F12.8)
5 format (I32)

! 2.Cavity
6 format (20(I10)," ")

! 3.Escape factor, number out of time and number of events
7 format (F15.8)
8 format (I23)
9 format (I25,I14,I12,I12)

! 4.Average values
! 10 format (200(F10.5)," ")
 10 format (20(F11.5)," ")
 ! 11 format (2(F10.1,F10.5)," ")
 20 format (400(F11.4)," ")

! 5.Surface distribution
! check number of indexes
12 format (13(F13.3)," ")

! 6.Data for format settings
13 format (I11,I25)
end



! 123456789012345678901234567890123456789012345678901234567890123