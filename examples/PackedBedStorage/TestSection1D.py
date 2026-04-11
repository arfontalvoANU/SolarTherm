#! /bin/env python

from __future__ import division
import unittest

from solartherm import simulation
import DyMat
import numpy as np
import matplotlib.pyplot as plt

from math import pi
import os
import shutil
import time

def get_segments(indices):
    segments = []
    start = indices[0]

    for i in range(1, len(indices)):
        if indices[i] != indices[i-1] + 1:
            segments.append((start, indices[i-1]))
            start = indices[i]

    segments.append((start, indices[-1]))
    return segments

class TestScheduler(unittest.TestCase):
    def setUp(self):
        self.modelname = 'Section1D'
        fn = f'{self.modelname}.mo'
        sim = simulation.Simulator(fn)
        sim.compile_model()
        sim.compile_sim(args=['-s'])
        sim.simulate(start=0, stop='10d', step='60s', solver='dassl', nls='homotopy', tolerance = '1e-06')
        self.mat = DyMat.DyMatFile(sim.res_fn)

    def test_sched(self):
        # Modelica time vector
        times = self.mat.abscissa('Tf[1]')[0]
        state = self.mat.data('state')

        chg_idx = np.where(state == 0)[0]
        dis_idx = np.where(state == 2)[0]

        chg_seg = get_segments(chg_idx)
        dis_seg = get_segments(dis_idx)

        idx = []
        for i in chg_seg:
            idx.append(i[1])

        colors = ['tab:blue','tab:orange']

        # number of axial nodes
        nz = int(self.mat.data('Nz')[0])

        # z grid from Modelica (constant over time, so read once)
        H = self.mat.data('H_tank')[0]
        z_modelica = np.zeros(nz)
        for i in range(nz):
            z_modelica[i] = self.mat.data(f'z[{i+1}]')[0]/H

        # Saving csv data
        csv = z_modelica
        headers = 'z,'

        # --- plotting ---
        fig, ax = plt.subplots(1, 1)

        charge_colors = []

        for idx, t_idx in enumerate(idx):
            # Modelica temperature profile at this time
            Tf = np.zeros(nz)
            for i in range(nz):
                Tf[i] = self.mat.data(f'Tf[{i+1}]')[t_idx]

            # --- Modelica: SOLID line ---
            line, = ax.plot(z_modelica,Tf,linestyle='-',label=f'n={idx+1}')
            charge_colors.append(line.get_color())  # store color
            csv = np.c_[csv,Tf]
            headers += f'Tchg_n{idx+1},'

        idx = []
        for i in dis_seg:
            idx.append(i[1])

        for idx, t_idx in enumerate(idx):
            # Modelica temperature profile at this time
            Tf = np.zeros(nz)
            for i in range(nz):
                Tf[i] = self.mat.data(f'Tf[{i+1}]')[t_idx]

            # --- Modelica: SOLID line ---
            ax.plot(z_modelica,Tf,linestyle='--',color=charge_colors[idx])
            csv = np.c_[csv,Tf]
            headers += f'Tdis_n{idx+1},'

        ax.set_xlabel('Non-dimensional height (-)')
        ax.set_ylabel('Fluid temperature (K)')
        ax.legend(loc='best',title='End of ith cycle')

        np.savetxt('PackedBedData.csv',csv,delimiter=',',header=headers,comments='')

        csv = times
        csv = np.c_[csv,self.mat.data('Tf[1]'),self.mat.data(f'Tf[{nz}]'),self.mat.data(f'E')]
        np.savetxt('TimeSeries.csv',csv,delimiter=',',header='t,Tb,Tt,E',comments='')

        plt.tight_layout()
        plt.show()

        os.system(f"find . -type f -name '{self.modelname}*' ! -name '*.mo' ! -name '*.mat' -delete")

if __name__ == '__main__':
    # Create working dir
    unittest.main()
