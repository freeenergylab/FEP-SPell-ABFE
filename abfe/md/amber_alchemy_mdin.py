#!/usr/bin/env python
# -*- coding:utf-8 -*-
"""
===============================================================================
CopyRight (c) By freeenergylab.
@Description:
This is the amber alchemy mdin module that can be used to construct AMBER
alchemical simulation input protocol.

@Author: Pengfei Li
@Date: May 30th, 2022
===============================================================================
"""

SOLVATED_MIN = """TI minimization stage
 &cntrl
  imin = 1,
  ntmin = 2,
  maxcyc = 10000,
  ncyc = 10000,
  drms = 0.01,
  drms = 0.0001,

  ntb = 1,
  cut = CUT,
  ntr = 0,
  restraint_wt = 0.0,

  ntp = 0,
  pres0 = 1.0,
  taup = 2.0,
  barostat = 2,

  ntc = 2,
  ntf = 1,

  ioutfm = 1,
  ntpr = 2500,
  ntwx = 0,

  icfe = 1,
  clambda = CLAMBDA,
  timask1 = TIMASK1,
  timask2 = TIMASK2,
  ifsc = 1,
  logdvdl = 0,
  scmask1 = SCMASK1,
  scmask2 = SCMASK2,
  scalpha = 0.2,
  scbeta = 50.0,
  
  gti_lam_sch = 1,
  gti_ele_sc = 1,
  gti_vdw_sc = 1,
  gti_scale_beta = 0,
  gti_cut_sc = 0,
  gti_add_sc = GTI_ADD_SC,
  tishake = 1,
  gti_syn_mass = 0,
  gti_output = 0,
  gti_cut = 1,
  gti_chg_keep = 1,
  /
"""

SOLVATED_HEAT = """TI heat MD
 &cntrl
  imin = 0,
  irest = 0,
  ntx = 1,
  nstlim = 10000,
  dt = DT,

  ntb = 1,
  cut = CUT,
  ntr = 1,
  restraint_wt = RESTRAINT_WT,
  restraintmask = "RESTRAINTMASK",

  ntp = 0,
  pres0 = 1.0,
  taup = 2.0,
  barostat = 2,

  ntc = 2,
  ntf = 1,

  ntt = 3,
  gamma_ln = 2.0,
  ig = -1,
  tempi = TEMPI,
  temp0 = TEMP0,

  ioutfm = 1,
  ntpr = 2500,
  ntwx = 0,

  icfe = 1,
  clambda = CLAMBDA,
  timask1 = TIMASK1,
  timask2 = TIMASK2,
  ifsc = 1,
  logdvdl = 0,
  scmask1 = SCMASK1,
  scmask2 = SCMASK2,
  scalpha = 0.2,
  scbeta = 50.0,
  
  gti_lam_sch = 1,
  gti_ele_sc = 1,
  gti_vdw_sc = 1,
  gti_scale_beta = 0,
  gti_cut_sc = 0,
  gti_add_sc = GTI_ADD_SC,
  tishake = 1,
  gti_syn_mass = 0,
  gti_output = 0,
  gti_cut = 1,
  gti_chg_keep = 1,

  nmropt = 1,
  /

  &wt
  type = 'TEMP_0',
  istep1 = 0, istep2 = 8000,
  value1 = TEMPI, value2 = TEMP0,
  /
  &wt
  type = 'END',
  /
"""

SOLVATED_PRESS = """TI press MD
 &cntrl
  imin = 0,
  irest = 1,
  ntx = 5,
  nstlim = 10000,
  dt = DT,

  ntb = 2,
  cut = CUT,
  ntr = 1,
  restraint_wt = RESTRAINT_WT,
  restraintmask = "RESTRAINTMASK",

  ntp = 1,
  pres0 = 1.0,
  taup = 2.0,
  barostat = 2,

  ntc = 2,
  ntf = 1,

  ntt = 3,
  gamma_ln = 2.0,
  ig = -1,
  temp0 = TEMP0,

  ioutfm = 1,
  ntpr = 2500,
  ntwx = 0,

  icfe = 1,
  clambda = CLAMBDA,
  timask1 = TIMASK1,
  timask2 = TIMASK2,
  ifsc = 1,
  logdvdl = 0,
  scmask1 = SCMASK1,
  scmask2 = SCMASK2,
  scalpha = 0.2,
  scbeta = 50.0,
  
  gti_lam_sch = 1,
  gti_ele_sc = 1,
  gti_vdw_sc = 1,
  gti_scale_beta = 0,
  gti_cut_sc = 0,
  gti_add_sc = GTI_ADD_SC,
  tishake = 1,
  gti_syn_mass = 0,
  gti_output = 0,
  gti_cut = 1,
  gti_chg_keep = 1,
 /
 """

SOLVATED_TI = """TI production MD
 &cntrl
  imin = 0,
  irest = 1,
  ntx = 5,
  nstlim = NSTLIM,
  !numexchg = NUMEXCHG,
  !gremd_acyc = 1,
  dt = DT,

  ntb = 2,
  cut = CUT,
  ntr = 0,
  restraint_wt = 0.0,

  ntp = 1,
  pres0 = 1.0,
  taup = 2.0,
  barostat = 2,

  ntc = 2,
  ntf = 1,

  ntt = 3,
  gamma_ln = 2.0,
  ig = -1,
  tempi = TEMPI,
  temp0 = TEMP0,

  ioutfm = 1,
  ntpr = 2500,
  ntwx = NTWX,

  icfe = 1,
  clambda = CLAMBDA,
  timask1 = TIMASK1,
  timask2 = TIMASK2,
  ifsc = 1,
  logdvdl = 0,
  scmask1 = SCMASK1,
  scmask2 = SCMASK2,
  scalpha = 0.2,
  scbeta = 50.0,
  
  gti_lam_sch = 1,
  gti_ele_sc = 1,
  gti_vdw_sc = 1,
  gti_scale_beta = 0,
  gti_cut_sc = 0,
  gti_add_sc = GTI_ADD_SC,
  tishake = 1,
  gti_syn_mass = 0,
  gti_output = 0,
  gti_cut = 1,
  gti_chg_keep = 1,

  ifmbar = 1,
  mbar_states = MBAR_STATES,
  mbar_lambda = MBAR_LAMBDA,
 /
 """

COMPLEX_MIN = """TI minimization
 &cntrl
  imin = 1,
  ntmin = 2,
  maxcyc = 10000,
  ncyc = 10000,
  drms = 0.01,
  drms = 0.0001,

  ntb = 1,
  cut = CUT,
  ntr = 0,
  restraint_wt = 0.0,

  ntp = 0,
  pres0 = 1.0,
  taup = 2.0,
  barostat = 2,

  ntc = 2,
  ntf = 1,

  ioutfm = 1,
  ntpr = 2500,
  ntwx = 0,

  icfe = 1,
  clambda = CLAMBDA,
  timask1 = TIMASK1,
  timask2 = TIMASK2,
  ifsc = 1,
  logdvdl = 0,
  scmask1 = SCMASK1,
  scmask2 = SCMASK2,
  scalpha = 0.2,
  scbeta = 50.0,
  
  gti_lam_sch = 1,
  gti_ele_sc = 1,
  gti_vdw_sc = 1,
  gti_scale_beta = 0,
  gti_cut_sc = 0,
  gti_add_sc = GTI_ADD_SC,
  tishake = 1,
  gti_syn_mass = 0,
  gti_output = 0,
  gti_cut = 1,
  gti_chg_keep = 1,

  nmropt = 1,
 /
 &wt
  type='END',
 &end
DISANG=restraints.inp
"""

COMPLEX_HEAT = """TI heat MD
 &cntrl
  imin = 0,
  irest = 0,
  ntx = 1,
  nstlim = 10000,
  dt = DT,
  vlimit = 10,

  ntb = 1,
  cut = CUT,
  ntr = 1,
  restraint_wt = RESTRAINT_WT,
  restraintmask = "RESTRAINTMASK",

  ntp = 0,
  pres0 = 1.0,
  taup = 2.0,
  barostat = 2,

  ntc = 2,
  ntf = 1,

  ntt = 3,
  gamma_ln = 2.0,
  ig = -1,
  tempi = TEMPI,
  temp0 = TEMP0,

  ioutfm = 1,
  ntpr = 2500,
  ntwx = 0,

  icfe = 1,
  clambda = CLAMBDA,
  timask1 = TIMASK1,
  timask2 = TIMASK2,
  ifsc = 1,
  logdvdl = 0,
  scmask1 = SCMASK1,
  scmask2 = SCMASK2,
  scalpha = 0.2,
  scbeta = 50.0,
  
  gti_lam_sch = 1,
  gti_ele_sc = 1,
  gti_vdw_sc = 1,
  gti_scale_beta = 0,
  gti_cut_sc = 0,
  gti_add_sc = GTI_ADD_SC,
  tishake = 1,
  gti_syn_mass = 0,
  gti_output = 0,
  gti_cut = 1,
  gti_chg_keep = 1,

  nmropt = 1,
 /
 &wt
 type = 'TEMP_0',
 istep1 = 0, istep2 = 8000,
 value1 = TEMPI, value2 = TEMP0,
 /
 &wt
 type = 'END',
 /
DISANG=restraints.inp
"""

COMPLEX_PRESS = """TI press MD
 &cntrl
  imin = 0,
  irest = 1,
  ntx = 5,
  nstlim = 10000,
  dt = DT,
  vlimit = 10,

  ntb = 2,
  cut = CUT,
  ntr = 1,
  restraint_wt = RESTRAINT_WT,
  restraintmask = "RESTRAINTMASK",

  ntp = 1,
  pres0 = 1.0,
  taup = 2.0,
  barostat = 2,

  ntc = 2,
  ntf = 1,

  ntt = 3,
  gamma_ln = 2.0,
  ig = -1,
  temp0 = TEMP0,

  ioutfm = 1,
  ntpr = 2500,
  ntwx = 0,

  icfe = 1,
  clambda = CLAMBDA,
  timask1 = TIMASK1,
  timask2 = TIMASK2,
  ifsc = 1,
  logdvdl = 0,
  scmask1 = SCMASK1,
  scmask2 = SCMASK2,
  scalpha = 0.2,
  scbeta = 50.0,

  gti_lam_sch = 1,
  gti_ele_sc = 1,
  gti_vdw_sc = 1,
  gti_scale_beta = 0,
  gti_cut_sc = 0,
  gti_add_sc = GTI_ADD_SC,
  tishake = 1,
  gti_syn_mass = 0,
  gti_output = 0,
  gti_cut = 1,
  gti_chg_keep = 1,

  nmropt = 1,
 /
 &wt
  type='END',
 &end
DISANG=restraints.inp
"""

COMPLEX_TI = """TI production MD
 &cntrl
  imin = 0,
  irest = 1,
  ntx = 5,
  nstlim = NSTLIM,
  !numexchg = NUMEXCHG,
  !gremd_acyc = 1,
  dt = DT,

  ntb = 2,
  cut = CUT,
  ntr = 0,
  restraint_wt = 0.0,

  ntp = 1,
  pres0 = 1.0,
  taup = 2.0,
  barostat = 2,

  ntc = 2,
  ntf = 1,

  ntt = 3,
  gamma_ln = 2.0,
  ig = -1,
  tempi = TEMPI,
  temp0 = TEMP0,

  ioutfm = 1,
  ntpr = 2500,
  ntwx = NTWX,

  icfe = 1,
  clambda = CLAMBDA,
  timask1 = TIMASK1,
  timask2 = TIMASK2,
  ifsc = 1,
  logdvdl = 0,
  scmask1 = SCMASK1,
  scmask2 = SCMASK2,
  scalpha = 0.2,
  scbeta = 50.0,
  
  gti_lam_sch = 1,
  gti_ele_sc = 1,
  gti_vdw_sc = 1,
  gti_scale_beta = 0,
  gti_cut_sc = 0,
  gti_add_sc = GTI_ADD_SC,
  tishake = 1,
  gti_syn_mass = 0,
  gti_output = 0,
  gti_cut = 1,
  gti_chg_keep = 1,

  ifmbar = 1,
  mbar_states = MBAR_STATES,
  mbar_lambda = MBAR_LAMBDA,

  nmropt = 1,
 /
 &wt
  type='DUMPFREQ', istep1=2500,
 &end
 &wt
  type='END',
 &end
DISANG=restraints.inp
DUMPAVE=restraints.out
"""

COMPLEX_RESTRAINT_MIN = """TI minimization
 &cntrl
  imin = 1,
  ntmin = 2,
  maxcyc = 10000,
  ncyc = 10000,
  drms = 0.01,
  drms = 0.0001,

  ntb = 1,
  cut = CUT,
  ntr = 0,
  restraint_wt = 0.0,

  ntp = 0,
  pres0 = 1.0,
  taup = 2.0,
  barostat = 2,

  ntc = 2,
  ntf = 1,

  ioutfm = 1,
  ntpr = 2500,
  ntwx = 0,

  icfe = 1,
  clambda = CLAMBDA,
  timask1 = ":1",
  timask2 = ":2",
  ifsc = 0,
  logdvdl = 0,

  nmropt = 1,
 /
 &wt
  type='END',
 &end
DISANG=restraints.inp
"""

COMPLEX_RESTRAINT_HEAT = """TI heat MD
 &cntrl
  imin = 0,
  irest = 0,
  ntx = 1,
  nstlim = 10000,
  dt = DT,
  vlimit = 20,

  ntb = 1,
  cut = CUT,
  ntr = 1,
  restraint_wt = RESTRAINT_WT,
  restraintmask = "RESTRAINTMASK",

  ntp = 0,
  pres0 = 1.0,
  taup = 2.0,
  barostat = 2,

  ntc = 2,
  ntf = 1,

  ntt = 3,
  gamma_ln = 2.0,
  ig = -1,
  tempi = TEMPI,
  temp0 = TEMP0,

  ioutfm = 1,
  ntpr = 2500,
  ntwx = 0,

  icfe = 1,
  clambda = CLAMBDA,
  timask1 = ":1",
  timask2 = ":2",
  ifsc = 0,
  logdvdl = 0,

  nmropt = 1,
 /
 &wt
 type = 'TEMP_0',
 istep1 = 0, istep2 = 8000,
 value1 = TEMPI, value2 = TEMP0,
 /
 &wt
 type = 'END',
 /
DISANG=restraints.inp
"""

COMPLEX_RESTRAINT_PRESS = """TI press MD
 &cntrl
  imin = 0,
  irest = 1,
  ntx = 5,
  nstlim = 10000,
  dt = DT,
  vlimit = 20,

  ntb = 2,
  cut = CUT,
  ntr = 1,
  restraint_wt = RESTRAINT_WT,
  restraintmask = "RESTRAINTMASK",

  ntp = 1,
  pres0 = 1.0,
  taup = 2.0,
  barostat = 2,

  ntc = 2,
  ntf = 1,

  ntt = 3,
  gamma_ln = 2.0,
  ig = -1,
  temp0 = TEMP0,

  ioutfm = 1,
  ntpr = 2500,
  ntwx = 0,

  icfe = 1,
  clambda = CLAMBDA,
  timask1 = ":1",
  timask2 = ":2",
  ifsc = 0,
  logdvdl = 0,

  nmropt = 1,
 /
 &wt
  type='END',
 &end
DISANG=restraints.inp
"""

COMPLEX_RESTRAINT_TI = """TI production MD
 &cntrl
  imin = 0,
  irest = 1,
  ntx = 5,
  nstlim = NSTLIM,
  !numexchg = NUMEXCHG,
  !gremd_acyc = 1,
  dt = DT,

  ntb = 2,
  cut = CUT,
  ntr = 0,
  restraint_wt = 0.0,

  ntp = 1,
  pres0 = 1.0,
  taup = 2.0,
  barostat = 2,

  ntc = 2,
  ntf = 1,

  ntt = 3,
  gamma_ln = 2.0,
  ig = -1,
  tempi = TEMPI,
  temp0 = TEMP0,

  ioutfm = 1,
  ntpr = 2500,
  ntwx = NTWX,

  icfe = 1,
  clambda = CLAMBDA,
  timask1 = ":1",
  timask2 = ":2",
  ifsc = 0,
  logdvdl = 0,

  ifmbar = 1,
  mbar_states = MBAR_STATES,
  mbar_lambda = MBAR_LAMBDA,

  nmropt = 1,
 /
 &wt
  type='DUMPFREQ', istep1=2500,
 &end
 &wt
  type='END',
 &end
DISANG=restraints.inp
DUMPAVE=restraints.out
"""