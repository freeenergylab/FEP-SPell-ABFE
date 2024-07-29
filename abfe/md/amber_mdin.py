#!/usr/bin/env python
# -*- coding:utf-8 -*-
"""
==============================================================================
CopyRight (c) By freeenergylab.
@Description:
This is the amber mdin module that can be used to construct AMBER simulation 
input protocol for equilibration stage.

@Author: Pengfei Li
@Date: May 30th, 2022
==============================================================================
"""
SOLVATED_MIN_1 = """Minimisation Stage 1
 &cntrl
  imin = 1,
  ntmin = 1,
  maxcyc = 10000,
  ncyc = 5000,
  dx0 = 0.01,
  drms = 0.0001,

  cut = CUT,
  ntr = 1,
  restraint_wt = RESTRAINT_WT,
  restraintmask = "RESTRAINTMASK",
 /
"""

SOLVATED_MIN_2 = """Minimisation Stage 2
 &cntrl
  imin = 1,
  ntmin = 1,
  maxcyc = 10000,
  ncyc = 5000,
  dx0 = 0.01,
  drms = 0.0001,

  cut = CUT,
  ntr = 0,
  restraint_wt = 0.0,
 /
"""

SOLVATED_HEAT="""Heat MD
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
  ntf = 2,

  ntt = 3,
  gamma_ln = 2.0,
  ig = -1,
  tempi = TEMPI,
  temp0 = TEMP0,

  ioutfm = 1,
  ntpr = 2500,
  ntwx = 0,

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

SOLVATED_PRESS="""Press MD
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
  ntf = 2,

  ntt = 3,
  gamma_ln = 2.0,
  ig = -1,
  temp0 = TEMP0,

  ioutfm = 1,
  ntpr = 2500,
  ntwx = 0,
 /
"""

SOLVATED_RELAX="""Relax MD
 &cntrl
  imin = 0,
  irest = 1,
  ntx = 5,
  nstlim = NSTLIM,
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
  ntf = 2,

  ntt = 3,
  gamma_ln = 2.0,
  ig = -1,
  tempi = TEMPI,
  temp0 = TEMP0,

  ioutfm = 1,
  ntpr = 2500,
  ntwx = 2500,
 /
"""

COMPLEX_MIN_1 = """Minimisation Stage 1
 &cntrl
  imin = 1,
  ntmin = 1,
  maxcyc = 10000,
  ncyc = 5000,
  dx0 = 0.01,
  drms = 0.0001,

  cut = CUT,
  ntr = 1,
  restraint_wt = RESTRAINT_WT,
  restraintmask = "RESTRAINTMASK",
 /
"""

COMPLEX_MIN_2 = """Minimisation Stage 2
 &cntrl
  imin = 1,
  ntmin = 1,
  maxcyc = 10000,
  ncyc = 5000,
  dx0 = 0.01,
  drms = 0.0001,

  cut = CUT,
  ntr = 0,
  restraint_wt = 0.0,
 /
"""

COMPLEX_HEAT = """Heat MD
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
  ntf = 2,

  ntt = 3,
  gamma_ln = 2.0,
  ig = -1,
  tempi = TEMPI,
  temp0 = TEMP0,

  ioutfm = 1,
  ntpr = 2500,
  ntwx = 0,

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

COMPLEX_PRESS = """Press MD
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
  ntf = 2,

  ntt = 3,
  gamma_ln = 2.0,
  ig = -1,
  temp0 = TEMP0,

  ioutfm = 1,
  ntpr = 2500,
  ntwx = 0,
 /
"""

COMPLEX_RELAX = """Relax MD
 &cntrl
  imin = 0,
  irest = 1,
  ntx = 5,
  nstlim = NSTLIM,
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
  ntf = 2,

  ntt = 3,
  gamma_ln = 2.0,
  ig = -1,
  tempi = TEMPI,
  temp0 = TEMP0,

  ioutfm = 1,
  ntpr = 2500,
  ntwx = 2500,
 /
"""