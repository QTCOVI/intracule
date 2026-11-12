c
c-----------------------------------------------------------------------
c
      subroutine etijcalc (m,lamx,la,lb,ce,a,b,ax,bx)
c
c-----This routine computes recursively the E_t^ij coefficients
c     resulting from the expansion of a cartesian gaussian product in
c     terms of hermite gaussians. The expressions used come from the
c     book "Molecular Electronic-Structure Theory" by T. Helgaker,
c     P.Jorgensen, and J. Olsen. All the coefficients arise from the
c     first one [E_0^00], taking also into account that E_t^ij=0 for
c     t<0 or t>i+j.
c
c-----Details
c
c     G_i (x,a,A_x) = x_A^i exp(-a x_A^2)  ;  x_A = x -  A_x
c     G_j (x,b,B_x) = x_B^j exp(-b x_B^2)  ;  x_B = x -  B_x
c
c     Omega_ij      = G_i (x,a,A_x) * G_j (x,b,B_x)
c
c     Omega_ij      = SUM [t=0,i+j] K_ab^x E_t^ij Lambda_t (x,p,P_x)
c
c     where
c
c     p                  = a+b
c     P_x                = (a A_x + b B_x)/p
c     mu                 = ab/(a+b)
c     X_AB               = A_x - B_x
c     K_ab^x             = exp(-mu X_AB^2)
c
c     The gives E_t^ij for  (0 <= i <= la; 0 <= j <= lb; 0 <= t <= i+j)
c
c     NUMERICAL NOTES (version 2, 2024):
c     --------------------------------
c     The recurrence itself is stable, but the auxiliary quantities
c
c       pa = -b*ab/p
c       pb = +a*ab/p
c       tp = 1/(2p)
c
c     can become very large when the exponents a and b are very
c     different (e.g. a=1e5, b=1e-1) AND the centers are far apart
c     (ab = ax-bx large). For the Kr atom (single center, ab=0), pa=pb=0
c     and there is no issue. For molecules, pa and pb can reach ~1e3
c     or more, and the products i*ce(...) can become large.
c
c     This version adds:
c       1. A check on the magnitude of pa, pb, tp. If any of them
c          exceeds a threshold (1e6), a warning is issued.
c       2. A check on the magnitude of ce(t,i,j,m) at the end. If any
c          value exceeds 1e200, a warning is issued (potential
c          overflow in downstream products).
c       3. In the worst case, a full log-scale version can be used,
c          but that requires changing the calling convention. Here we
c          keep the same interface and just add the warnings.
c
c-----------------------------------------------------------------------
c
      implicit none
      integer la,lb,lab,i,j,t,i1,j1,t1,m,lamx
      real(kind=8) ce(-1:2*lamx,-1:lamx,-1:lamx,3)
      real(kind=8) a,b,ax,bx,p,ab,pa,pb,tp
      real(kind=8) maxce, warn_threshold
      logical      warned_pa
c
      if (la.lt.0)    stop 'etijcalc.f: Fatal error, la < 0 '
      if (lb.lt.0)    stop 'etijcalc.f: Fatal error, lb < 0 '
      if (la.gt.lamx) stop 'etijcalc.f: Fatal error, la > lamx '
      if (lb.gt.lamx) stop 'etijcalc.f: Fatal error, lb > lamx '
      lab=la+lb
      ce(-1:lab,-1:la,-1:lb,m)=0.0D+00
      ce(0,0,0,m)=1.0D+00
      if (lab.eq.0) return
      p  = a + b
      ab = ax - bx
      pa = - b*ab/p
      pb = + a*ab/p
      tp = 1.0d+00/(2.0d+00*p)
c
c-----Check on the magnitude of pa, pb, tp
c
      warn_threshold = 1.0d6
      warned_pa = .false.
      if (abs(pa) .gt. warn_threshold .or. 
     &    abs(pb) .gt. warn_threshold .or. 
     &    abs(tp) .gt. warn_threshold) then
        if (.not. warned_pa) then
          write (*,1) '# etijcalc: WARNING large pa,pb,tp = ',pa,pb,tp
          write (*,2) '# etijcalc:   a,b = ',a,b,'  ab = ',ab
          warned_pa = .true.
        endif
      endif
c
      do i=0,la
        i1=i-1
        do j=0,lb
          j1=j-1
          do t=1,i+j
            t1=t-1
            ce(t,i,j,m)=tp*(i*ce(t1,i1,j,m)+j*ce(t1,i,j1,m))/dble(t)
          enddo
          if (i.lt.la) ce(0,i+1,j,m)=pa*ce(0,i,j,m)+ce(1,i,j,m)
          if (j.lt.lb) ce(0,i,j+1,m)=pb*ce(0,i,j,m)+ce(1,i,j,m)
        enddo
      enddo
c
c-----Check on the magnitude of ce
c
      maxce = 0.0d0
      do i = 0, la
        do j = 0, lb
          do t = 0, i+j
            if (abs(ce(t,i,j,m)) .gt. maxce) maxce = abs(ce(t,i,j,m))
          enddo
        enddo
      enddo
      if (maxce .gt. 1.0d200) then
        write (*,3) '# etijcalc: WARNING max|ce| = ',maxce, 
     &   ' for la = ',la,', lb = ',lb
        write (*,4) 
     &   '# etijcalc: potential overflow in downstream products.'
      endif
c
 1    format (1x,a,3e12.4)
 2    format (1x,a,2e12.4,a,e12.4) 
 3    format (1x,a,e12.4,a,i3,a,i3)
 4    format (1x,a)
      return
      end
