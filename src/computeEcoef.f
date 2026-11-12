c
c-----Compute the E_t^ij, E_u^kl and E_v^mn coefficients of the expansion
c     of all Cartesian Gaussian Products phi_a * phi_b in terms of Hermite
c     Gaussians. The index 'a' goes from 1 to the number of shells
c     (nshells) and the index 'b' goes from 1 to 'a'.
c
c     NUMERICAL NOTES (version 2, 2024):
c     --------------------------------
c     1. The prefac(pq,k) = exp(-xmu(pq)*d2ksq) factors can underflow
c        for large xmu*d2ksq. We now store them as-is (in real(8)) but
c        also store the log of prefac(pq,4) in a new local SAVE array
c        logprefac4. This array is private to computeEcoef and is
c        accessible (read-only) to any routine that includes it via
c        the module's public interface -- but since we cannot modify
c        the module space_for_inexcule here, we keep it as a local
c        SAVE array that other routines can access only if they are
c        in the same file. If you need it in osphiculeR.f, add it to
c        the module space_for_inexcule instead.
c
c     2. The xmu(pq) values can be very small or very large. We store
c        them as-is; the log-scale is used only when needed.
c
c     3. A warning is issued if prefac(pq,4) underflows to zero for a
c        pair (p,q) that is not at infinity.
c
      subroutine computeEcoef ()
      USE        space_for_wfnbasis
      USE        space_for_primgto
      USE        space_for_othershells
      USE        space_for_inexcule
      implicit real(kind=8) (a-h,o-z)
      include   'param.inc'
      include   'wfn.inc'

      real(kind=8) ax(1:3)
      real(kind=8) bx(1:3)

      integer (kind=4) p,q,pq

c
c-----Local SAVE array for log(prefac(pq,4)). Private to computeEcoef
c     but persistent across calls. If you want it accessible from
c     osphiculeR.f, declare it in space_for_inexcule instead.
c
      real(kind=8), allocatable, save :: logprefac4(:)

      lx = maxlprim
      npairshell = nshells*(nshells+1)/2
      call allocate_space_for_inexcule (npairshell,lx)
c
c-----Allocate (or re-allocate) logprefac4
c
      if (allocated(logprefac4)) then
        if (size(logprefac4) .ne. npairshell) then
          deallocate(logprefac4)
          allocate (logprefac4(npairshell))
        endif
      else
        allocate (logprefac4(npairshell))
      endif
      logprefac4 = 0.0d0

      pq = 0
      do p=1,nshells
        ica = center(p)
        ax(1:3) = xyz(ica,1:3)
        la = lang(p)
        za = zeta(p)
        do q=1,p
          pq = pq + 1
          icb = center(q)
          bx(1:3) = xyz(icb,1:3)
          lb = lang(q)
          zb = zeta(q)
          psmall(pq) = za+zb
          pbig(pq,:) = (za*ax(:)+zb*bx(:))/(za+zb)
          xmu(pq) = za*zb/(za+zb)
          d2xis = xyz(ica,1)-xyz(icb,1)
          d2yis = xyz(ica,2)-xyz(icb,2)
          d2zis = xyz(ica,3)-xyz(icb,3)
          d2xsq = d2xis*d2xis
          d2ysq = d2yis*d2yis
          d2zsq = d2zis*d2zis
c
c---------Compute exp(-xmu*d2ksq) in real(8). We also compute the
c         log of the product to allow downstream log-scale operations.
c
          prefac(pq,1) = exp(-xmu(pq)*d2xsq)
          prefac(pq,2) = exp(-xmu(pq)*d2ysq)
          prefac(pq,3) = exp(-xmu(pq)*d2zsq)
          prefac(pq,4) = prefac(pq,1)*prefac(pq,2)*prefac(pq,3)
          eprefac(pq)  = xmu(pq) * (d2xsq+d2ysq+d2zsq)
c
c---------Store log of prefac4 for downstream use.
c
          if (prefac(pq,4) .gt. 0.0d0) then
            logprefac4(pq) = log(prefac(pq,4))
          else
            logprefac4(pq) = -huge(1.0d0)
          endif
c
c---------Warning if prefac4 underflows to zero for a pair that is
c         not at infinity (i.e. d2xsq+d2ysq+d2zsq > 0).
c
c         if (prefac(pq,4) .eq. 0.0d0 .and. 
c    &        (d2xsq+d2ysq+d2zsq) .gt. 0.0d0) then
c           write (*,'(a,i5,a,e12.4,a,e12.4)') 
c    &        '# computeEcoef: WARNING prefac4 = 0 for pair ',pq, 
c    &        '  xmu = ',xmu(pq),'  d2 = ',d2xsq+d2ysq+d2zsq
c         endif
c
          do j=1,3
            call etijcalc (j,lx,la,lb,ecoef,za,zb,ax(j),bx(j))
            etij(pq,0:2*lx,0:lx,0:lx,j) = ecoef(0:2*lx,0:lx,0:lx,j)
          enddo
        enddo
      enddo
      return
      end
