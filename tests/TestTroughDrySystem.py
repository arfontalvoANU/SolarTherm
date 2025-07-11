# -*- coding: utf-8 -*-
from solartherm import simulation
from solartherm import postproc
from solartherm import params
import itertools
import xlsxwriter
import functools
import multiprocessing as mp
import os
import shutil
from sys import path
from time import time
import DyMat
import math
import numpy as np
import argparse

class stcompile:
	def __init__(self,wd):
		os.chdir(wd)                                                                 # Changing working directory
		sim = simulation.Simulator('%s/TroughDrySystem.mo'%(wd))                     # Compilation of the modelica model
		print('Compiling model')
		mcargs = ['-d=nonewInst']
		sim.compile_model(args=mcargs)
		sim.compile_sim(args=['-s'])
		sim.load_init()                                                              # Generating the xml file

def simulate(wd, par_n, i, par_v):
	os.chdir(wd)                                                                     # Changing working directory
	sim = simulation.Simulator('%s/TroughDrySystem.mo'%(wd), suffix=str(i))          # Create simulator
	sim.load_init()
	sim.update_pars(par_n,par_v)                                                     # Parameters to update
	sim.simulate(start=0, stop='1y', step='300s',solver='dassl', nls='homotopy')     # Simulation of the model
	res = postproc.SimResultElec(sim.res_fn)
	perf = res.calc_perf()
	res={}
	mat = DyMat.DyMatFile('%s/TroughDrySystem_res_%s.mat'%(wd,i))                    # Loading result file
	res['sm'] = mat.data('SM')[0]                                                    # Solar multiple
	res['storage'] = mat.data('t_storage')[0]                                        # Full load hours of storage
	res['T_air'] = mat.data('T_des_air_out')[0]                                       # Temperature of the air at the outlet of the collection system (inlet of rotative dryer)
	res['lcoh'] = perf[1]
	res['capf'] = perf[2]
	return res

class clean_comp_files:
	def __init__(self,wd):
		os.chdir(wd)
		os.system('rm -rf *.o *.c *.h *.json *.makefile *.xml TroughDrySystem')

def simulation_callback(i,worksheet,perf):
	k = 0
	worksheet.write(i+1,k,perf['sm']); k = k + 1
	worksheet.write(i+1,k,perf['storage']); k = k + 1
	worksheet.write(i+1,k,perf['T_air']); k = k + 1
	worksheet.write(i+1,k,perf['lcoh']); k = k + 1
	worksheet.write(i+1,k,perf['capf'])
	return None

def parameter_sweep(T_air,async):
	wd=os.path.join(os.getenv('HOME'),'PROJECTS','solartherm-arturo','examples')
	stcompile(wd)                                                                    # Compiling and initialising the simulation
	sol_multi = np.arange(1,10.1,1)                                   # Defining solar multiple to sweep
	t_storage = np.arange(5,30.1,5)                                    # Defining storage hours to sweep
	par_v = []
	for i,j in enumerate(itertools.product(sol_multi,t_storage)):
		j = map(str,list(j))
		par_v.append(j)

	par_n = ['SM','t_storage']

	workbook = xlsxwriter.Workbook('results.xlsx')                                   # Writting results in text and xlsx files
	worksheet = workbook.add_worksheet()
	col = 0
	worksheet.write(0,col,'solar_multiple'); col += 1
	worksheet.write(0,col,'t_storage'); col += 1
	worksheet.write(0,col,'T_air'); col += 1
	worksheet.write(0,col,'LCOH'); col += 1
	worksheet.write(0,col,'capacity_factor'); col += 1

	if async:
		objfunc = functools.partial(simulate, wd, par_n)
		pool = mp.Pool(processes=mp.cpu_count())
		for i,vals in enumerate(par_v):
			pool.apply_async(objfunc, args=(i, vals), callback=functools.partial(simulation_callback, i, worksheet))
		pool.close()
		pool.join()
	else:
		for i,vals in enumerate(par_v):
			perf=simulate(wd, par_n, i, vals)
			simulation_callback(i,worksheet,perf)

	workbook.close()
	os.system('cp %s/results.xlsx /home/arfontalvo/Dropbox/solarpaces_2023/decarbonisation/T%s.xlsx'%(wd,T_air))
	clean_comp_files(wd)
	return None

if __name__=="__main__":
	parser = argparse.ArgumentParser()
	parser.add_argument('T')
	parser.add_argument('--async',action='store_true')
	parser.add_argument('--sequential', dest='async', action='store_false')
	parser.set_defaults(async=True)
	args = parser.parse_args()

	tinit = time()
	parameter_sweep(args.T,args.async)
	seconds = time() - tinit
	m, s = divmod(seconds, 60)
	h, m = divmod(m, 60)
	print('Simulation time: {:d}:{:02d}:{:02d}'.format(int(h), int(m), int(s)))
