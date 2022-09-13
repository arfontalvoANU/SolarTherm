#! /bin/env python2

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
		sim.update_pars(['T_hot_set','T_hot_start','state_hot_set.h','controlCold.state_ref.h'],
			[str(T),str(T),str(h),str(h)])
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

