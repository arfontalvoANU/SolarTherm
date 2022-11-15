#!/usr/bin/env python3
import os
import numpy as np
import scipy.io
import matplotlib.pyplot as plt
import math
import time
import ctypes
from numpy.ctypeslib import ndpointer
from tqdm import tqdm

from srlife import receiver, solverparams, library, thermal, structural, system, damage, managers
from neml import uniaxial

strMatNormal = lambda a: [''.join(s).rstrip() for s in a]
strMatTrans  = lambda a: [''.join(s).rstrip() for s in zip(*a)]
sign = lambda x: math.copysign(1.0, x)

SODIR = os.path.join(os.getenv('HOME'),'solartherm/SolarTherm/Resources/Include/stress')

class receiver_cyl:
	def __init__(self,coolant = 'salt', Ri = 57.93/2000, Ro = 60.33/2000, T_in = 290, T_out = 565,
                      nz = 450, nt = 46, R_fouling = 0.0, ab = 0.94, em = 0.88, kp = 16.57, H_rec = 10.5, D_rec = 8.5,
                      nbins = 50, alpha = 15.6e-6, Young = 186e9, poisson = 0.31,
                      thermat = "base",defomat = "const_base",damat = "base",mat='800H',
                      debugfolder = os.path.expanduser('~'), debug = False, verification = False, Dittus=True):
		self.coolant = coolant
		self.Ri = Ri
		self.Ro = Ro
		self.thickness = Ro - Ri
		self.T_in = T_in + 273.15
		self.T_out = T_out + 273.15
		self.nz = nz
		self.nt = nt
		self.nr = 9
		self.R_fouling = R_fouling
		self.ab = ab
		self.em = em
		self.kp = kp
		self.H_rec = H_rec
		self.D_rec = D_rec
		self.dz = H_rec/nbins
		self.nbins=nbins
		self.debugfolder = debugfolder
		self.debug = debug
		self.verification = verification
		self.sigma = 5.670374419e-8
		# Discretisation parameters
		self.dt = np.pi/(nt-1)
		# Tube section diameter and area
		self.d = 2.*Ri                               # Tube inner diameter (m)
		self.area = 0.25 * np.pi * pow(self.d,2.)    # Tube flow area (m2)
		self.ln = np.log(Ro/Ri)                      # Log of Ro/Ri simplification
		#Auxiliary variables
		cosines = np.cos(np.linspace(0.0, np.pi, nt))
		self.cosines = np.maximum(cosines, np.zeros(nt))
		self.theta = np.linspace(-np.pi, np.pi,self.nt*2-2)
		self.n = 3
		self.l = alpha
		self.E = Young
		self.nu = poisson
		l = self.E*self.nu/((1+self.nu)*(1-2*self.nu));
		m = 0.5*self.E/(1+self.nu);
		props = l*np.ones((3,3)) + 2*m*np.identity(3)
		self.invprops = np.linalg.inv(props)
		# Loading dynamic library
		so_file = "%s/stress.so"%SODIR
		stress = ctypes.CDLL(so_file)
		self.fun = stress.curve_fit
		self.fun.argtypes = [ctypes.c_int,
						ctypes.c_double,
						ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"),
						ndpointer(ctypes.c_double, flags="C_CONTIGUOUS")]
		self.Dittus = Dittus
		# Choose the material models
		self.thermal_mat, self.deformation_mat, self.damage_mat = library.load_material(mat, thermat, defomat, damat)

	def density(self,T):
		"""
		   Thermal and transport properties of nitrate salt and sodium
		   Nitrate salt:  Zavoico, A. B. Solar power tower design basis document, revision 0; Topical. Sandia National Labs., 2001.
		   Liquid sodium: Fink, J. K.; Leibowitz, L. Thermodynamic and transport properties of sodium liquid and vapor. Argonne National Lab., 1995.
		"""
		if self.coolant == 'salt':
			d = 2090.0 - 0.636 * (T - 273.15)
		else:
			d = 219.0 + 275.32 * (1.0 - T / 2503.7) + 511.58 * np.sqrt(1.0 - T / 2503.7)
		return d

	def dynamicViscosity(self,T):
		"""
		   Thermal and transport properties of nitrate salt and sodium
		   Nitrate salt:  Zavoico, A. B. Solar power tower design basis document, revision 0; Topical. Sandia National Labs., 2001.
		   Liquid sodium: Fink, J. K.; Leibowitz, L. Thermodynamic and transport properties of sodium liquid and vapor. Argonne National Lab., 1995.
		"""
		if self.coolant == 'salt':
			eta = 0.001 * (22.714 - 0.120 * (T - 273.15) + 2.281e-4 * pow((T - 273.15),2) - 1.474e-7 * pow((T - 273.15),3))
			emax = 0.001 * (22.714 - 0.120 * 600. + 2.281e-4 * pow(600.,2) - 1.474e-7 * pow(600.,3))
			eta = np.maximum(eta,emax)
		else:
			eta = np.exp(-6.4406 - 0.3958 * np.log(T) + 556.835/T)
		return eta

	def thermalConductivity(self,T):
		"""
		   Thermal and transport properties of nitrate salt and sodium
		   Nitrate salt:  Zavoico, A. B. Solar power tower design basis document, revision 0; Topical. Sandia National Labs., 2001.
		   Liquid sodium: Fink, J. K.; Leibowitz, L. Thermodynamic and transport properties of sodium liquid and vapor. Argonne National Lab., 1995.
		"""
		if self.coolant == 'salt':
			k = 0.443 + 1.9e-4 * (T - 273.15)
		else:
			k = 124.67 - 0.11381 * T + 5.5226e-5 * pow(T,2) - 1.1842e-8 * pow(T,3);
		return k;

	def specificHeatCapacityCp(self,T):
		"""
		   Thermal and transport properties of nitrate salt and sodium
		   Nitrate salt:  Zavoico, A. B. Solar power tower design basis document, revision 0; Topical. Sandia National Labs., 2001.
		   Liquid sodium: Fink, J. K.; Leibowitz, L. Thermodynamic and transport properties of sodium liquid and vapor. Argonne National Lab., 1995.
		"""
		if self.coolant == 'salt':
			C = 1396.0182 + 0.172 * T
		else:
			C = 1000 * (1.6582 - 8.4790e-4 * T + 4.4541e-7 * pow(T,2) - 2992.6 * pow(T,-2))
		return C

	def phis(self,x):
		return np.exp(-1./x)/(np.exp(-1./x) + np.exp(-1./(1-x)))

	def Temperature(self, m_flow, Tf, Tamb, CG, h_ext):
		"""
		    Flow and thermal variables:
		    hf: Heat transfer coefficient due to internal forced-convection
		    mu: HTF dynamic viscosity (Pa-s)
		    kf: HTF thermal conductivity (W/m-K)
		    C:  HTF specific heat capacity (J/kg-K)
		    Re: HTF Reynolds number
		    Pr: HTF Prandtl number
		    Nu: Nusselt number due to internal forced convection
		"""

		Tf,temp = np.meshgrid(np.ones(self.nt),Tf)
		Tf = Tf*temp

		# HTF thermo-physical properties
		mu = self.dynamicViscosity(Tf)                 # HTF dynamic viscosity (Pa-s)
		kf = self.thermalConductivity(Tf)              # HTF thermal conductivity (W/m-K)
		C = self.specificHeatCapacityCp(Tf)            # HTF specific heat capacity (J/kg-K)

		m_flow,temp = np.meshgrid(np.ones(self.nt), m_flow)
		m_flow = m_flow*temp

		Tamb,temp = np.meshgrid(np.ones(self.nt), Tamb)
		Tamb = Tamb*temp

		h_ext,temp = np.meshgrid(np.ones(self.nt), h_ext)
		h_ext = h_ext*temp

		# HTF internal flow variables
		Re = m_flow * self.d / (self.area * mu)    # HTF Reynolds number
		Pr = mu * C / kf                           # HTF Prandtl number
		Re_zer = np.where(Re<=0)[0]
		Re_pos = np.where(Re>0)[0]
		Re_lam = np.where((Re>0)     & (Re<=2e3))[0]
		Re_tra = np.where((Re>2e3) & (Re<4e3))[0]
		Re_tur = np.where(Re>=4e3)[0]
		f = np.zeros_like(Re)
		Nu= 4.36*np.ones_like(Re)
		Pr_neg = np.where(Pr<=0)[0]

		if self.coolant == 'salt':
			if self.Dittus:
				Nu = 0.023 * pow(Re, 0.8) * pow(Pr, 0.4)
			else:
				f[Re_lam] = 64/Re[Re_lam]
				f[Re_tur] = pow(1.82*np.log10(Re[Re_tur]) - 1.64, -2)
				f[Re_tra] = 64/Re[Re_tra] + self.phis((Re[Re_tra]-2e3)/2e3)*(pow(1.82*np.log10(Re[Re_tra]) - 1.64, -2.) - 64/Re[Re_tra])

				Nu[Re_lam] = 4.36
				Nu[Re_tur] = (f[Re_tur]/8)*(Re[Re_tur] - 1000)*Pr[Re_tur]/(1 + 12.7*pow(f[Re_tur]/8, 0.5)*(pow(Pr[Re_tur],0.66)-1))
				Nu[Re_tra] = 4.36 + self.phis((Re[Re_tra]-2e3)/2e3)*((f[Re_tra]/8.)*(Re[Re_tra] - 1000.)*Pr[Re_tra]/(1. + 12.7*pow(f[Re_tra]/8., 0.5)*(pow(Pr[Re_tra], 0.66) -1.)) - 4.36)

		else:
			Nu = 5.6 + 0.0165 * pow(Re*Pr, 0.85) * pow(Pr, 0.01)

		# HTF internal heat transfer coefficient
		hf = Nu * kf / self.d
		if self.R_fouling>0:
			hf[Re_pos] = 1./(1./hf[Re_pos] + self.R_fouling)

		# Calculating heat flux at circumferential nodes
		cosinesm,fluxes = np.meshgrid(self.cosines,CG)
		qabs = fluxes*cosinesm 
		a = -((self.em*(self.kp + hf*self.ln*self.Ri)*self.Ro*self.sigma)/((self.kp + hf*self.ln*self.Ri)*self.Ro*(self.ab*qabs + self.em*self.sigma*pow(Tamb,4)) + hf*self.kp*self.Ri*Tf + (self.kp + hf*self.ln*self.Ri)*self.Ro*Tamb*(h_ext)))
		b = -((hf*self.kp*self.Ri + (self.kp + hf*self.ln*self.Ri)*self.Ro*(h_ext))/((self.kp + hf*self.ln*self.Ri)*self.Ro*(self.ab*qabs + self.em*self.sigma*pow(Tamb,4)) + hf*self.kp*self.Ri*Tf + (self.kp + hf*self.ln*self.Ri)*self.Ro*Tamb*(h_ext)))
		c1 = 9.*a*pow(b,2.) + np.sqrt(3.)*np.sqrt(-256.*pow(a,3.) + 27.*pow(a,2)*pow(b,4))
		c2 = (4.*pow(2./3.,1./3.))/pow(c1,1./3.) + pow(c1,1./3.)/(pow(2.,1./3.)*pow(3.,2./3.)*a)
		To = -0.5*np.sqrt(c2) + 0.5*np.sqrt((2.*b)/(a*np.sqrt(c2)) - c2)
		Ti = (To + hf*self.Ri*self.ln/self.kp*Tf)/(1 + hf*self.Ri*self.ln/self.kp)
		qnet = hf*(Ti - Tf)
		_qnet = np.concatenate((qnet[:,1:],qnet[:,::-1]),axis=1)
		Qnet = _qnet.sum(axis=1)*self.Ri*self.dt*self.dz
		net_zero = np.where(Qnet<0)[0]
		Qnet[net_zero] = 0.0
		_qnet[net_zero,:] = 0.0
		self.qnet = _qnet

		# Fourier coefficients
		indexs = np.zeros((self.nr*self.nt,2))
		stress = np.zeros((self.nr*self.nt,6))
		strain = np.zeros((self.nr*self.nt,6))
		BDp = self.Fourier(Ti)
		BDpp = self.Fourier(To)
		header='i,j,stress_rr,stress_tt,stress_zz,stress_rt,stress_rz,stress_tz,strain_rr,strain_tt,strain_zz,strain_rt,strain_rz,strain_tz'
		res = 0
		for i,rr in enumerate(np.linspace(self.Ri, self.Ro, self.nr)):
			for j,theta in enumerate(np.linspace(0.0, np.pi, self.nt)):
				tt = self.temperature(rr,theta,BDp,BDpp)
				stress[res,:],strain[res,:] = self.Thermoelastic(tt,rr,theta, BDp, BDpp)
				indexs[res,:] = i,j
				res += 1
		csv = np.c_[indexs,stress,strain]
		np.savetxt('S31609_OD33.40_WT1.651_450.csv',csv,delimiter=',',header=header)
		return Qnet

	def temperature(self,r,theta,BDp,BDpp):
		n = 1
		B0 = (BDp[0]-BDpp[0])/np.log(self.Ri/self.Ro)
		A0 = BDp[0]-B0*np.log(self.Ri)
		ss = A0 + B0*np.log(r)
		for Bpn,Bppn in zip(BDp[1:],BDpp[1:]):
			Gn=(Bpn-Bppn)/(self.Ri**n-self.Ro**n)
			Bn=(Bpn-Gn*self.Ri**n)/(self.Ri**n+self.Ro**n)*self.Ri**n*self.Ro**n
			An=Gn+Bn/(self.Ri**n*self.Ro**n)
			ss+=(An*r**n + Bn/r**n)*np.cos(n*theta)
			n+=1
		return ss

	def Fourier(self,T):
		coefs = np.empty(21)
		self.fun(self.nt, self.dt, T, coefs)
		return coefs

	def Thermoelastic(self, T, r, theta, BDp, BDpp):
		Tbar_i = BDp[0]; BP = BDp[1]; DP = BDp[2];
		Tbar_o = BDpp[0]; BPP = BDpp[1]; DPP = BDpp[2];
		a = self.Ri; b = self.Ro; a2 = a*a; b2 = b*b; r2 = r*r; r4 = pow(r,4);

		C = self.l*self.E/(2.*(1. - self.nu));
		D = 1./(2.*(1. + self.nu));
		kappa = (Tbar_i - Tbar_o)/np.log(b/a);
		kappa_theta = r*a*b/(b2 - a2)*((BP*b - BPP*a)/(b2 + a2)*np.cos(theta) + (DP*b - DPP*a)/(b2 + a2)*np.sin(theta));
		kappa_tau   = r*a*b/(b2 - a2)*((BP*b - BPP*a)/(b2 + a2)*np.sin(theta) - (DP*b - DPP*a)/(b2 + a2)*np.cos(theta));

		T_theta = T - ((Tbar_i - Tbar_o) * np.log(b/r)/np.log(b/a)) - Tbar_o;

		Qr = kappa*C*(0 -np.log(b/r) -a2/(b2 - a2)*(1 -b2/r2)*np.log(b/a) ) \
					+ kappa_theta*C*(1 - a2/r2)*(1 - b2/r2);
		Qtheta = kappa*C*(1 -np.log(b/r) -a2/(b2 - a2)*(1 +b2/r2)*np.log(b/a) ) \
					+ kappa_theta*C*(3 -(a2 +b2)/r2 -a2*b2/r4);
		Qz = kappa*self.nu*C*(1 -2*np.log(b/r) -2*a2/(b2 - a2)*np.log(b/a) ) \
					+ kappa_theta*2*self.nu*C*(2 -(a2 + b2)/r2) -self.l*self.E*T_theta;
		Qrtheta = kappa_tau*C*(1 -a2/r2)*(1 -b2/r2);

		Q_Eq = np.sqrt(0.5*(pow(Qr -Qtheta,2) + pow(Qr -Qz,2) + pow(Qz -Qtheta,2)) + 6*pow(Qrtheta,2));
		Q = np.zeros((6,))
		Q[0] = Qr; Q[1] = Qtheta; Q[2] = Qz;
		e = np.zeros((6,))
		e[0] = 1/self.E*(Qr - self.nu*(Qtheta + Qz))
		e[1] = 1/self.E*(Qtheta - self.nu*(Qr + Qz))
		e[2] = -self.l*T_theta

		if self.verification:
			print("=============== NPS Sch. 5S 1\" S31609 at 450degC ===============")
			print("Biharmonic coefficients:")
			print("Tbar_i [K]:      749.6892       %4.4f"%Tbar_i)
			print("  B'_1 [K]:      45.1191        %4.4f"%BP)
			print("  D'_1 [K]:      -0.0000        %4.4f"%DP)
			print("Tbar_o [K]:      769.7119       %4.4f"%Tbar_o)
			print(" B''_1 [K]:      79.4518        %4.4f"%BPP)
			print(" D''_1 [K]:      0.0000         %4.4f\n"%DPP)
			print("Stress at outside tube crown:")
			print("Q_r [MPa]:       0.0000         %4.4f"%(Qr/1e6))
			print("Q_rTheta [MPa]:  0.0000         %4.4f"%(Qrtheta/1e6))
			print("Q_theta [MPa]:  -101.0056       %4.4f"%(Qtheta/1e6))
			print("Q_z [MPa]:      -389.5197       %4.4f"%(Qz/1e6))
			print("Q_Eq [MPa]:      350.1201       %4.4f"%(Q_Eq/1e6))

		return Q,e

if __name__=='__main__':
	tinit = time.time()
	k = 21;alpha=20e-6;E = 165e9;nu = 0.31
	T_int = 723.15;T_amb = 293.15;CG = 0.85e6;mdot = 5
	h_ext=30           # convective loss due to wind W/(m2.K)
	a = 30.098/2e3     # inside tube radius [mm->m]
	b = 33.400/2e3     # outside tube radius [mm->m]
	model = receiver_cyl(
		Ri = a
		,Ro = b
		,alpha = alpha
		,Young = E
		,nt=91
		,ab = 0.968
		,em = 0.87
		,kp = k
		,verification=False
		,Dittus=True)
	model.Temperature(mdot,T_int,T_amb,CG,h_ext)
	seconds = time.time() - tinit
	m, s = divmod(seconds, 60)
	h, m = divmod(m, 60)
	print('Simulation time: {:d}:{:02d}:{:02d}'.format(int(h), int(m), int(s)))
