	subroutine ckcalc (nstep,ustc,ustw,ustr,omega,lch,vk,nu,z,
     1	    dzu,dzl,capk,capks,dcapk,dcapks,zeta,KTYPE,wbl,cf)
c
	include 'maxes.inc'
c
	dimension capk(MAXZ),capks(MAXZ),dcapk(MAXZ),dcapks(MAXZ)
	dimension z(MAXZ),dzu(MAXZ),dzl(MAXZ),zeta(MAXZ)
	real nu,ustr,ustc,omega,vk,lc,lw,lch,ustw,beta,betar,gam
	character*3 KTYPE
	parameter (beta=5.4,betar=7.3,gam=1.0)
c
 	lc=ustc/(cf*6.)
 	if (lc.gt.lch) lc=lch
	if (lc.lt.250.) lc=250.
	lw=wbl*ustr/(omega*6.)
	do n=1,nstep
	    uscfz2=ustc*ustc*exp(-2.*z(n)/lc)
	    uswfz2=ustw*ustw*exp(-2.*z(n)/lw)
	    usfz=(uscfz2*uscfz2+uswfz2*uswfz2+2.*uscfz2*uswfz2)**0.25
	    if (KTYPE.eq.'NEU'.or.KTYPE.eq.'neu') then
	        capk(n)=vk*z(n)*usfz
	        capks(n)=capk(n)
	    elseif (KTYPE.eq.'STR'.or.KTYPE.eq.'str') then
	        capk(n)=vk*z(n)*usfz/(1.+beta*zeta(n))
	        capks(n)=vk*z(n)*usfz/(gam+betar*zeta(n))
	    end if
	    if (z(n).gt.1.and.capk(n).lt.nu) capk(n)=nu
	    if (z(n).gt.1.and.capks(n).lt.nu) capks(n)=nu
	enddo
	do n=2,nstep-1
	    dcapk(n)=.5*((capk(n+1)-capk(n))/dzu(n)+
     1		(capk(n)-capk(n-1))/dzl(n))
	    dcapks(n)=.5*((capks(n+1)-capks(n))/dzu(n)+
     1		(capks(n)-capks(n-1))/dzl(n))
	end do
	dcapk(1)=dcapk(2)
	dcapk(nstep)=(capk(nstep)-capk(nstep-1))/dzu(nstep-1)
	dcapks(1)=dcapks(2)
	dcapks(nstep)=(capks(nstep)-capks(nstep-1))/dzu(nstep-1)
c
	return
	end
