# Ne atom, RHF/cc-pVTZ (spherical functions; the writer expands them in Cartesian primitives), written in AIMPAC .wfn format
from pyscf import gto, scf
from pyscf.tools import wfn_format
mol = gto.M(atom='Ne 0 0 0', basis='cc-pvtz', verbose=0)
mf = scf.RHF(mol).run(conv_tol=1e-12)
print('E(RHF) =', mf.e_tot)
with open('ne.wfn', 'w') as f:
    wfn_format.write_mo(f, mol, mf.mo_coeff[:, mf.mo_occ > 0], mo_occ=mf.mo_occ[mf.mo_occ > 0],
                        mo_energy=mf.mo_energy[mf.mo_occ > 0])
