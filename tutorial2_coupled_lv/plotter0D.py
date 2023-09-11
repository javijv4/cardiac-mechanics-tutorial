#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 19 13:50:48 2022

@author: Javiera Jilberto Vallejos
"""

import matplotlib.pyplot as plt
import numpy as np

class plotter0D:
    def __init__(self, path, stepsxcycle, vol0, name =''):
        # Some names of this list might change depending on the version and the problem you're solving
        name_list = ['q_vin_r', 'p_ven_pul', 'p_at_l', 'Q_v_l', 'V_v_l', 'p_ar_sys',
                     'V_ven_pul', 'V_ar_sys', 'p_v_r', 'q_ar_pul', 'q_vin_l', 'p_at_r',
                     'q_vout_r', 'q_ven1_sys', 'V_ar_pul', 'V_at_l', 'V_ven_sys',
                     'q_vout_l', 'p_v_l', 'V_at_r', 'q_ar_sys', 'q_ven1_pul', 'V_v_r',
                     'p_ar_pul', 'p_ven_sys', 'q_arp_sys', 'p_ard_sys']

        self.compute_3Dvolume(path, vol0, name=name)
        self.vol0 = vol0

        # Load python output
        self.dict_vars = {}
        for n in name_list:
            self.dict_vars[n] = np.loadtxt(path + 'results_' + name + '_' + n + '.txt')

        # Compute number of cycles
        self.steps = len(self.dict_vars['p_v_l'])
        self.stepsxcycle = stepsxcycle
        self.ncycles = np.floor(self.steps/self.stepsxcycle).astype(int)
        if self.ncycles == 0: self.ncycles = 1

    def compute_3Dvolume(self, path, vol_n, name=''):
        ch = 'v_l'
        tmp, fluxes = np.loadtxt(''+path+'results_' + name + '_Q_'+ch+'.txt', unpack = True)
        # integrate volume (mid-point rule): Q_{n+1} = (V_{n+1} - V_{n})/dt --> V_{n+1} = Q_{n+1}*dt + V_{n}
        filename_vol = path+'/results_' + name + '_V_'+ch+'.txt'
        file_vol = open(filename_vol, 'wt')
        file_vol.write('%.16E %.16E\n' % (tmp[0], vol_n))
        for n in range(len(fluxes)-1):
            dt = tmp[n+1] - tmp[n]
            vol = -fluxes[n+1]*dt + vol_n
            file_vol.write('%.16E %.16E\n' % (tmp[n+1], vol))
            vol_n = vol
        file_vol.close()

    def get_start_stop(self, cycle):
        if cycle == 'all':
            start = 0
            stop = self.steps
        elif cycle == 'last':
            stop = self.stepsxcycle*self.ncycles
            start = stop - self.stepsxcycle
        elif type(cycle) == int:
            start = self.stepsxcycle*(cycle-1)
            stop = self.stepsxcycle*(cycle)

        return start, stop

    def plot_Wiggers(self, side, cycle, unit='KPa', offset=0., ax=None):

        if unit == 'mmHg':
            unit_mult = 7.50062
        elif unit == 'KPa':
            unit_mult = 1.

        start, stop = self.get_start_stop(cycle)

        if ax == None:
            plt.figure()
            ax = plt.gca()

        if side == 'l':
            ax.plot(self.dict_vars['p_v_l'][start:stop,0]+offset, self.dict_vars['p_v_l'][start:stop,1]*unit_mult, 'b', label = r'$p_V^\ell$')
            ax.plot(self.dict_vars['p_at_l'][start:stop,0]+offset, self.dict_vars['p_at_l'][start:stop,1]*unit_mult, 'g', ls = '-.', ms = 1,  label = r'$p_{at}^{\ell}$')
            ax.plot(self.dict_vars['p_ar_sys'][start:stop,0]+offset, self.dict_vars['p_ar_sys'][start:stop,1]*unit_mult, 'r', ls = '--', label = r'$p_{ar}^{sys}$')
            ax.plot(self.dict_vars['p_ven_sys'][start:stop,0]+offset, self.dict_vars['p_ven_sys'][start:stop,1]*unit_mult, 'y', ls = ':', ms = 1, label = r'$p_{ven}^{sys}$')

        elif side == 'r':
            ax.plot(self.dict_vars['p_v_r'][start:stop,0]+offset, self.dict_vars['p_v_r'][start:stop,1]*unit_mult, 'b', label = r'$p_V^r$')
            ax.plot(self.dict_vars['p_at_r'][start:stop,0]+offset, self.dict_vars['p_at_r'][start:stop,1]*unit_mult, 'g', ls = '-.', ms = 1,  label = r'$p_{at}^{r}$')
            ax.plot(self.dict_vars['p_ar_pul'][start:stop,0]+offset, self.dict_vars['p_ar_pul'][start:stop,1]*unit_mult, 'r', ls = '--', label = r'$p_{ar}^{pul}$')
            ax.plot(self.dict_vars['p_ven_pul'][start:stop,0]+offset, self.dict_vars['p_ven_pul'][start:stop,1]*unit_mult, 'y', ls = ':', ms = 1, label = r'$p_{ven}^{pul}$')

        ax.legend()
        ax.set_xlabel('Time [s]')
        ax.set_ylabel(r'$p$ ['+unit+']')

        return ax

    def plot_PVloop(self, side, cycle, unit='KPa', ax = None, cmap = 'viridis', **kwargs):

        if unit == 'mmHg':
            unit_mult = 7.50062
        elif unit == 'KPa':
            unit_mult = 1.

        start, stop = self.get_start_stop(cycle)
        # print(start, stop)
        if ax == None:
            plt.figure()
            ax = plt.gca()
        if cycle == 'all':
            for i in range(self.ncycles+1):
                st = start + self.stepsxcycle*i
                end = start + self.stepsxcycle*(i+1) +1
                cmap = plt.get_cmap(cmap)  # this returns a colormap
                color = cmap(1 - float(i)/(self.ncycles)) #  returns a color for each x between 0.0 and 1.0
                ax.plot(self.dict_vars['V_v_'+side][st:end,1]/1000, self.dict_vars['p_v_'+side][st:end,1]*unit_mult, color=color, label = r'Cycle ' + str(i+1))

        elif cycle == 'last':
            for i in range(self.ncycles-1,self.ncycles):
                ax.plot(self.dict_vars['V_v_'+side][start:stop,1]/1000, self.dict_vars['p_v_'+side][start:stop,1]*unit_mult,  **kwargs)

        elif type(cycle) == int:
            ax.plot(self.dict_vars['V_v_'+side][start:stop,1]/1000, self.dict_vars['p_v_'+side][start:stop,1]*unit_mult, label = r'Cycle ' + str(cycle+1))

        plt.xlabel(r'$V$ [ml]')
        plt.ylabel(r'$p$ ['+unit+']')

        return ax
