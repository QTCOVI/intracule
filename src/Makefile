#=======================================================================
#  Makefile para Promolden - versión corregida (v4)
#
#  NOTA IMPORTANTE:
#     El módulo "intracule_integral" está definido DENTRO del archivo
#     newintegrateIabcd.f90.  Por eso NO existe un objeto
#     "intracule_integral.o"; los archivos que hacen
#     "USE intracule_integral" deben depender de newintegrateIabcd.o.
#=======================================================================

FC      = gfortran
F77     = gfortran

FFLAGS  = -fbounds-check -g -fstack-protector -fcheck=all -fbacktrace \
          -Wall -ffree-line-length-none -fopenmp -no-pie -O2

LDFLAGS = -fopenmp -no-pie
LIBS    = -lgsl -lgslcblas -lm

PROG    = intracule

#-----------------------------------------------------------------------
#  Lista de módulos objeto (.o)
#
#  ORDEN: los módulos base (mod_prec, mod_inexcule, etc.) van primero;
#  los que dependen de ellos van después; los .o "hoja" (que nadie usa
#  pero que usan a otros) van al final.
#-----------------------------------------------------------------------
MODS = mod_prec.o \
       mod_futils.o \
       wfnbasis.o \
       primgto.o \
       wfncoef.o \
       corrarea.o \
       rdmspace.o \
       slmarea.o \
       mod_inexcule.o \
       angles_space.o \
       atdata.o \
       datatm.o \
       othershell.o \
       mod_rdm.o \
       dot4d_kahan.o \
       kummer_hypergeometric.o \
       expminus_x1F_1.o \
       newmodified_spherical_bessel.o \
       real_spherical_harmonics.o \
       slmangint.o \
       newintegrateIabcd.o \
       binary128.o \
       computeEcoef_q.o \
       radial_integ.o \
       parallel.o \
       getargs.o \
       info.o \
       initial.o \
       getdate.o \
       jacobi.o \
       lebgrid.o \
       rdwfn.o \
       setdble.o \
       setint.o \
       setword.o \
       timer.o \
       uppcase.o \
       cmcq.o \
       iqcksort.o \
       etijcalc.o \
       error.o \
       find_free_unit.o \
       writexyz.o \
       zfromSymbol.o \
       banner.o \
       computeEcoef.o \
       hermipol.o \
       clasgrid.o \
       acumula_intracule.o \
       otroiculedrv.o \
       newiculedrv.o \
       qcksort.o \
       gjrule_burkardt.o \
       newrdm2prim.o \
       iculegen1.o \
       iculegen0.o \
       iculegen2.o \
       icule0_disk.o \
       icule1_disk.o \
       icule2_disk.o \
       icule0_mem.o \
       icule1_mem.o \
       icule2_mem.o \
       hermitec.o \
       polarpow_space.o \
       gsphicule0.o \
       lebicule_claude1.o \
       lebicule.o \
       newtayloric.o \
       onecenteric.o \
       tayloric_opt2.o \
       tayloric_opt.o \
       cltayloric.o \
       osphiculeR.o \
       hermicart.o \
       ngtox.o \
       moicule.o \
       prequadpack.o \
       funbess.o \
       funbecke.o \
       rmeshgl.o \
       alltdsg.o \
       intracule.o \
       raditp.o \
       icule0grid.o \
       icule2grid.o \
       intraditp.o \
       init_yntphi_recursive.o \
       init_thetaval_recursive.o \
       thetaintegrals.o \
       phiintegrals.o \
       integrate_bessic.o \
       quadpac_bessic.o \
       quadpac_becke.o \
       safe_compute_B_l.o \
       gtogto.o \
       muller.o \
       quadpac.o \
       usage.o \
       usage_intracule.o \

#-----------------------------------------------------------------------
#  Regla principal
#-----------------------------------------------------------------------
$(PROG): $(MODS)
	$(FC) $(LDFLAGS) -o $@ $(MODS) $(LIBS)

#-----------------------------------------------------------------------
#  Reglas genéricas de compilación
#-----------------------------------------------------------------------
%.o: %.f90
	$(FC) $(FFLAGS) -c $< -o $@

%.o: %.f
	$(FC) $(FFLAGS) -c $< -o $@

# alltdsg.f holds only DATA statements (t-design grids, 63 MB): compile it
# without optimisation or run-time checks. With $(FFLAGS) gfortran needs
# hours; with -O0 about 2 minutes and ~5 GB of memory.
alltdsg.o: alltdsg.f
	$(FC) -O0 -c $< -o $@

#-----------------------------------------------------------------------
#  Reglas explícitas de compilación
#  (solo para documentar qué fuente genera cada .o)
#-----------------------------------------------------------------------
mod_prec.o:             mod_prec.f90
mod_futils.o:           mod_futils.f90
wfnbasis.o:             wfnbasis.f
wfncoef.o:              wfncoef.f
corrarea.o:             corrarea.f
parallel.o:             parallel.f
rdmspace.o:             rdmspace.f
primgto.o:              primgto.f
slmarea.o:              slmarea.f
mod_inexcule.o:         mod_inexcule.f
othershell.o:           othershell.f
angles_space.o:         angles_space.f
atdata.o:               atdata.f
datatm.o:               datatm.f
slmangint.o:            slmangint.f90
real_spherical_harmonics.o: real_spherical_harmonics.f90
kummer_hypergeometric.o:    kummer_hypergeometric.f90
expminus_x1F_1.o:       expminus_x1F_1.f90
dot4d_kahan.o:          dot4d_kahan.f90
newmodified_spherical_bessel.o: newmodified_spherical_bessel.f90
mod_rdm.o:              mod_rdm.f90
getargs.o:              getargs.f
info.o:                 info.f
initial.o:              initial.f
getdate.o:              getdate.f
jacobi.o:               jacobi.f
lebgrid.o:              lebgrid.f
rdwfn.o:                rdwfn.f
setdble.o:              setdble.f
setint.o:               setint.f
setword.o:              setword.f
timer.o:                timer.f
uppcase.o:              uppcase.f
cmcq.o:                 cmcq.f
iqcksort.o:             iqcksort.f
etijcalc.o:             etijcalc.f
error.o:                error.f
find_free_unit.o:       find_free_unit.f
writexyz.o:             writexyz.f
zfromSymbol.o:          zfromSymbol.f
banner.o:               banner.f90
computeEcoef.o:         computeEcoef.f
hermipol.o:             hermipol.f
clasgrid.o:             clasgrid.f90
acumula_intracule.o:    acumula_intracule.f90
otroiculedrv.o:         otroiculedrv.f90
newiculedrv.o:          newiculedrv.f
qcksort.o:              qcksort.f
newrdm2prim.o:          newrdm2prim.f
iculegen1.o:            iculegen1.f
iculegen0.o:            iculegen0.f
iculegen2.o:            iculegen2.f
icule0_disk.o:          icule0_disk.f
icule1_disk.o:          icule1_disk.f
icule2_disk.o:          icule2_disk.f
icule0_mem.o:           icule0_mem.f
icule1_mem.o:           icule1_mem.f
icule2_mem.o:           icule2_mem.f
hermitec.o:             hermitec.f
polarpow_space.o:       polarpow_space.f
gsphicule0.o:           gsphicule0.f
lebicule_claude1.o:     lebicule_claude1.f
lebicule.o:             lebicule.f
newtayloric.o:          newtayloric.f
onecenteric.o:          onecenteric.f
tayloric_opt2.o:        tayloric_opt2.f
tayloric_opt.o:         tayloric_opt.f
cltayloric.o:           cltayloric.f
osphiculeR.o:           osphiculeR.f
hermicart.o:            hermicart.f
ngtox.o:                ngtox.f
moicule.o:              moicule.f90
gjrule_burkardt.o:      gjrule_burkardt.f90
newintegrateIabcd.o:    newintegrateIabcd.f90
prequadpack.o:          prequadpack.f
funbess.o:              funbess.f90
funbecke.o:             funbecke.f90
rmeshgl.o:              rmeshgl.f
seltdsg.o:              seltdsg.f
intracule.o:            intracule.f
raditp.o:               raditp.f90
icule0grid.o:           icule0grid.f
icule2grid.o:           icule2grid.f
intraditp.o:            intraditp.f
init_yntphi_recursive.o: init_yntphi_recursive.f90
init_thetaval_recursive.o: init_thetaval_recursive.f
thetaintegrals.o:       thetaintegrals.f
phiintegrals.o:         phiintegrals.f
radial_integ.o:         radial_integ.f
binary128.o:            binary128.f90
computeEcoef_q.o:       computeEcoef_q.f
integrate_bessic.o:     integrate_bessic.f
quadpac_bessic.o:       quadpac_bessic.f
quadpac_becke.o:        quadpac_becke.f
safe_compute_B_l.o:     safe_compute_B_l.f90
gtogto.o:               gtogto.f
muller.o:               muller.f
quadpac.o:              quadpac.f

#-----------------------------------------------------------------------
#  Dependencias de módulos
#-----------------------------------------------------------------------
# --- Base ---
wfnbasis.o:     mod_prec.o
wfncoef.o:      mod_prec.o
corrarea.o:     mod_prec.o
parallel.o:     mod_prec.o
mod_futils.o:   mod_prec.o
rdmspace.o:     mod_prec.o
slmarea.o:      mod_prec.o
atdata.o:       mod_prec.o
datatm.o:       mod_prec.o
angles_space.o: mod_prec.o
real_spherical_harmonics.o: mod_prec.o
kummer_hypergeometric.o:    mod_prec.o
newmodified_spherical_bessel.o: mod_prec.o
rmeshgl.o:      mod_prec.o

# --- Intermedios ---
primgto.o:      wfnbasis.o
othershell.o:   primgto.o
mod_inexcule.o: mod_prec.o
slmangint.o:    angles_space.o
expminus_x1F_1.o:  mod_prec.o kummer_hypergeometric.o
mod_rdm.o:      wfnbasis.o wfncoef.o

# --- newintegrateIabcd.o define el módulo intracule_integral ---
newintegrateIabcd.o: mod_prec.o mod_inexcule.o

# --- Archivos que USAN intracule_integral ---
prequadpack.o:  mod_inexcule.o newintegrateIabcd.o
funbess.o:      mod_inexcule.o angles_space.o newintegrateIabcd.o
funbecke.o:     mod_inexcule.o angles_space.o newintegrateIabcd.o
radial_integ.o: newintegrateIabcd.o binary128.o
computeEcoef_q.o: binary128.o wfnbasis.o othershell.o
otroiculedrv.o: binary128.o

# --- integrate_bessic.o usa muchos módulos ---
integrate_bessic.o: newintegrateIabcd.o mod_inexcule.o mod_rdm.o \
                    slmangint.o real_spherical_harmonics.o \
                    newmodified_spherical_bessel.o kummer_hypergeometric.o \
                    expminus_x1F_1.o safe_compute_B_l.o dot4d_kahan.o \
                    othershell.o

# --- Cálculos con intracule ---
newiculedrv.o:  mod_rdm.o
iculegen1.o:    mod_inexcule.o
iculegen0.o:    mod_inexcule.o
iculegen2.o:    mod_inexcule.o
icule0_disk.o:  mod_inexcule.o
icule1_disk.o:  mod_inexcule.o
icule2_disk.o:  mod_inexcule.o
icule0_mem.o:   mod_inexcule.o
icule1_mem.o:   mod_inexcule.o
icule2_mem.o:   mod_inexcule.o
lebicule.o:     angles_space.o
newtayloric.o:  mod_inexcule.o mod_rdm.o angles_space.o othershell.o
onecenteric.o:  mod_inexcule.o mod_rdm.o angles_space.o othershell.o

# --- quadpac_becke/bessic usan muchos módulos ---
quadpac_becke.o:  mod_inexcule.o mod_rdm.o othershell.o angles_space.o \
                  primgto.o wfnbasis.o hermitec.o hermicart.o ngtox.o \
                  timer.o moicule.o dot4d_kahan.o computeEcoef.o
quadpac_bessic.o: mod_inexcule.o mod_rdm.o othershell.o angles_space.o \
                  primgto.o wfnbasis.o hermitec.o hermicart.o ngtox.o \
                  timer.o moicule.o dot4d_kahan.o computeEcoef.o

# --- Programa principal y utilidades ---
intracule.o:    mod_inexcule.o mod_rdm.o newiculedrv.o
moicule.o:      wfncoef.o primgto.o
gtogto.o:       primgto.o wfnbasis.o

#-----------------------------------------------------------------------
#  Limpieza
#-----------------------------------------------------------------------
.PHONY: clean cleanall
clean:
	rm -f *.o *.mod $(PROG)

cleanall: clean
	rm -f *~ *.bak
