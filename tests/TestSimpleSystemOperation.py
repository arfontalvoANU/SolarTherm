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
		fn = '../examples/SimpleSystemOperation.mo'
		sim = simulation.Simulator(fn)
		sim.compile_model(args=['-d=nonewInst'])
		sim.compile_sim(args=['-s'])
		sim.load_init()
		sim.simulate(start=0, stop='1y', step='1m', solver='dassl', nls='homotopy', args=['-noEventEmit'])
		self.res = postproc.SimResultElec(sim.res_fn)
		self.perf = self.res.calc_perf()

	def test_simple_system(self):
		print('Simulation finished')
		print('Starting post-processing')
		header = ['epy (MWh/year)','lcoe ($/MWh)','capf (%)','srev ($)']
		print(' '.join('%s'%('%s'%x + ' '*(max(len(header[i]),len('%s'%self.perf[i])) - len('%s'%x))) for i,x in enumerate(header)))
		print(' '.join('%s'%('%s'%x + ' '*(max(len(header[i]),len('%s'%self.perf[i])) - len('%s'%x))) for i,x in enumerate(self.perf)))
		os.system('cp SimpleSystemOperation*.mat ../examples/')
		os.system('rm SimpleSystemOperation*')


if __name__ == '__main__':
	unittest.main()

