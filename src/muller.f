!-----------------------------------------------------------------------
!     Rutina que construye las tres 2-RDM (total, Coulomb y XC) en la
!     aproximacion de Muller (equivalente al funcional BB del articulo
!     de Garcia-Revilla et al., J. Phys. Chem. A 2012, 116, 1237-1250).
!
!     Notacion:
!         ec. 3:   rho2 = rho*rho - rho2_xc
!         ec. 14:  rho2_xc = sum_{i,j} 2 sqrt(n_i n_j) chi_i(1)chi_j(1)chi_i(2)chi_j(2)
!
!     Convenio:
!         rdm2(i,k,j,l)  ->  electron 1 = (i,k),  electron 2 = (j,l)
!
!     Entrada:
!         nmo    : numero de orbitales naturales
!         occ(i) : numeros de ocupacion natural (0 <= occ(i) <= 1)
!
!     Salida:
!         rdm1    (i,k): 1-RDM total
!         rdm2    (i,k,j,l): 2-RDM total
!-----------------------------------------------------------------------
!
      subroutine muller (nmo, rdm1, rdm2)
!
      USE space_for_wfnbasis
      implicit none
!
      integer, intent(in)  :: nmo
      real(8), intent(out) :: rdm1(nmo,nmo)
      real(8), intent(out) :: rdm2(nmo,nmo,nmo,nmo)
!
      integer :: i, k, j, l
      real(8) :: coul, xc
!
      rdm1 = 0.0d+00
      rdm2 = 0.0d+00
!
      do i = 1, nmo
        do k = 1, nmo
          if (i == k) rdm1(i,k) = occ(i)
          do j = 1, nmo
            do l = 1, nmo
!
!-------------Termino de Coulomb:  delta_{ik} delta_{jl} n_i n_j
!             (ec. 3, primer termino de la derecha)
!
              coul = 0.0d+00
              if (i.eq.k .and. j.eq.l) coul = occ(i) * occ(j)
!
!-------------Termino de intercambio-correlacion (Muller / BB):
!             delta_{ik} delta_{jl} 2 sqrt(n_i n_j)
!             (ec. 14; el delta viene de integrar chi_i chi_j sobre
!             todo el espacio, ya que en la base de NOs
!             int chi_i(r) chi_j(r) dr = delta_{ij})
!
              xc = 0.0d+00
              if (i.eq.j .and. k.eq.l) xc = sqrt(occ(i)*occ(k))
!
!-------------2-RDM total: ec. 3, rho2 = rho*rho - rho2_xc
!
              rdm2(i,k,j,l) = coul - xc
            enddo
          enddo
        enddo
      enddo
!
      return
      end subroutine muller
