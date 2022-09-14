#! /bin/env python

from __future__ import division
import unittest

from solartherm import simulation
from solartherm import postproc

import numpy as np
import scipy.io as sio
import DyMat
import os

class TestFlowPathStress(unittest.TestCase):
	def setUp(self):
		fn = '../examples/GemasolarSystemOperation.mo'
		sim = simulation.Simulator(fn)
		sim.compile_model(args=['-d=nonewInst'])
		sim.compile_sim(args=['-s'])
		sim.load_init()
		T = 565 + 273.15
		h = 1396.0182*T + 0.086*pow(T,2)
		T = str(T)
		h = str(h)
		parent = os.path.join(os.path.expanduser('~'),'ownCloud/phd_update/damage/N06230/N06230_OD22.40_WT1.20_565')
		par_n = [
			'T_hot_set'
			,'T_hot_start'
			,'state_hot_set.h'
			,'controlCold.state_ref.h'
			,'opt_file'
			,'hav_file'
			,'file_dni1'
			,'file_dni2'
			,'file_dni3'
			,'file_dni4'
			,'file_mflow'
		]
		par_v = [
			T
			,T
			,h
			,h
			,os.path.join(parent,'OELTs_Solstice.motab')
			,os.path.join(parent,'HALTs_Solstice.motab')
			,os.path.join(parent,'FLUX_fp1_d0.56.motab')
			,os.path.join(parent,'FLUX_fp1_d0.87.motab')
			,os.path.join(parent,'FLUX_fp1_d1.0.motab')
			,os.path.join(parent,'FLUX_fp1_d1.39.motab')
			,os.path.join(parent,'MFLOW_Solstice_fp1.motab')
		]
		sim.update_pars(par_n,par_v)
		sim.simulate(start=0, stop='1y', step='5m', solver='dassl', nls='homotopy')
		self.res = postproc.SimResultElec(sim.res_fn)
		self.perf = self.res.calc_perf()

	def test_simple_system(self):
		print('Simulation finished')
		print('Starting post-processing')
		header = ['epy (MWh/year)','lcoe ($/MWh)','capf (%)','srev ($)']
		print(' '.join('%s'%('%s'%x + ' '*(max(len(header[i]),len('%s'%self.perf[i])) - len('%s'%x))) for i,x in enumerate(header)))
		print(' '.join('%s'%('%s'%x + ' '*(max(len(header[i]),len('%s'%self.perf[i])) - len('%s'%x))) for i,x in enumerate(self.perf)))
		os.system('cp GemasolarSystemOperation*.mat ../examples/')
		os.system('rm GemasolarSystemOperation*')


if __name__ == '__main__':
	unittest.main()

