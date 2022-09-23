# WT-MetaD-HREX_HJ
Well-tempered MetaDynamics with Hamiltonian Replica Exchange protocol used in Holliday Junction dynamics study.

###
Files role:
Note that these files are example used for WT-MetaD-HREX in HJ isomer I-open CV space. 
For Isomer II-open CV space, simply change the CV and the corresponding sigma value, etc.
###

run-md_gmx.sh:
  Job script to run WT-MetaD-HREX simulations in GROMACS
  It requires 6 replica tpr files, preplumed.dat and plumed.dat.
  Run WT-MetaD-HREX for 1 µs.

plumed.dat and preplumed.dat:
  PLUMED script including WT-MetaD-HREX protocol, HBfix settings and output control on-the-fly.
  preplumed.dat was used for the first 100 ns and plumed.dat was used for the rest 900 ns (in run-md_gmx.sh, we set a for-loop to run 10 100-ns simulations).
  
Check.ipynb:
  Used for averaging the biases in the final 300 ns WT-MetaD-HREX in all replicas
  Input HILLS.* files from all replicas, which are outputs from the run, defined by plumed.dat and preplumed.dat
  Output A_HILLS.* files, which would be used for the following free energies calculation.

plum-ave-reweight.dat:
  Used for calculating free-energy landscapes.
  Requires A_HILLS.* file and COLVAR.* file as input (put in the same direction).
  plumed driver --plumed plum-ave-reweight.dat --noatoms --kt 2.478957

blocks_bs.py:
  Python script used for bootstrapping analyses (standard deviation of free energy difference against open state).
  Check comment inside for details.

GHBFIX*, gHBfix*, scalingParameters.dat, typesTable.dat:
  Used for HBfix in GROMACS, executed by PLUMED

HJ_isoi.gro:
  Initial structure for WT-MetaD-HREX (isomer I-open)


