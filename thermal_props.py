##########################################
# PREAMBLE
##########################################

import numpy as np
from scipy import constants as sc
from scipy import integrate as i_
from scipy import interpolate
from scipy.signal import savgol_filter
# import skfda # (Need py 3.9)

##########################################
# THERMAL CONDUCTIVITY (NIST)
##########################################

class thermalProps:
    """
    .. versionadded:: 0.1
    
    Thermal conductivity at low temperatures, as per described at: https://trc.nist.gov/cryogenics/materials/6063_T5%20Alulminum/6063-T5Aluminum_rev.htm and https://www.matweb.com/index.aspx
    
    """
    
    # def k_mat_cal(self, z_m, rho_m, M_m):
        
    #     # number of electrons per m3
    #     n = (sc.Avogadro**2)*(z_m*rho_m/(M_m))
        
    #     # mean collision time
    #     Z = np.sqrt() # mean collision frequency
    #     tau = 1/Z
        
    #     ke = (1/3)*(np.pi**2)*(n*sc.Boltzmann*t/sc.electron_mass)
        
    #     k_mat = ke + kph
        
    #     return k_mat
    
    def cpd_cve_calc(self, T, z, ThetaD, Tm): 
        """        
        Compile the electronic contributiom to the heat capacity of a selected material as a function of temperature.
        
        :param T: Temperature range [K]
        :param z: Valence [na]
        :param ThetaD: Debye temperature [K]
        :param Tm: Melting temperature [K]
        :type T: 1d numpy array
        :type z: float
        :type ThetaD: float
        :type Tm: float
        :return cve: Specific heat capacity array, function of T [J mol-1 K-1]
        :rtype cve: array of floats
        
        .. note:: None
        
        """
        
        return (5/(24*(np.pi**3)))*z*(ThetaD**3)/((T**2)*Tm)

    def aly_cve_calc(self, cves): 
        """        
        Compile the electronic contributiom to the heat capacity of a selected alloy as a function of temperature.
        
        :param cves: Electronic contributions to heat capacity of each compiunds [K]
        :type cves: 1d numpy array
        :return cve_aly: Specific heat capacity array for a selected alloy, function of T [J mol-1 K-1]
        :rtype cve_aly: array of floats
        
        .. note:: None
        
        """
        
        cve_mult = np.zeros(np.shape(cves))
        cves_sub = cves
        
        for i in range(0, len(cves)):
            cve_mult[i] = np.multiply(cves[i], cves_sub[i+1])
            cves_sub = np.delete(cves_sub, i)
        
        cve_aly = np.sum(cves, axis = 0) + np.product(cves, axis = 0) + np.sum(cve_mult, axis = 0)
        
        return cve_aly
    
    def cpd_atmvol_calc(self, r): 
        """        
        Compile the atomic critical volume of a selected material.
        
        :param r: Atomic radius [m]
        :type r: float
        :return V: Atomic critical volume [m^3]
        :rtype V: float
        
        .. note:: None
        
        """
        
        return (4/3)*np.pi*(r**3)
    
    def aly_debye_props(self, thetaD_aly, V_aly, v_aly):
        """        
        Compile the properties related to the Debye model for a selected alloy.
        
        :param thetaD_aly: Debye temperature [K]
        :param V_aly: Critical volume [m^3]
        :param v_aly: Sound speed in alloy [m.s^-1]
        :type thetaD_aly: float
        :type V_aly: float
        :type v_aly: float
        :return w_aly: Critical frequency [rad.s^-1]
        :return D_aly: Density of state [Hz]
        :return n_aly: Atomic density in a grain of compound [na]
        :rtype w_aly: float
        :rtype D_aly: float
        :rtype n_aly: float
        
        .. note:: None
        
        """
        
        # Critical frequency of the alloy
        w_aly = sc.k*thetaD_aly/sc.hbar

        # Density of state of the compounds ('The number of different states at a particular energy level that electrons are allowed to occupy')
        D_aly = V_aly*(w_aly**2)/(2*(np.pi**2)*(v_aly**3)) # Using Al volume as it dominates, should use Gibbs Thomson relation...

        # Density of atoms in a grain of compound
        n_aly = (1/(6*(np.pi**2)))*((sc.k*thetaD_aly/(sc.hbar*w_aly))**3)
        
        return w_aly, D_aly, n_aly
        
    def aly_cv_rot_calc(self, T, x_pcts, Mmols, rs, n_aly, cpds_nr, w_aly, thetaD_aly, rho_aly, Ey_aly, Tm_aly, rho_C, Ey_C):
        
        cv_rot_cpds = np.zeros(cpds_nr)
        for j in range(1,int(n_aly)+2):
            cv_rot_cpds += x_pcts*j*(j+1)/(Mmols*(rs**2))
        cv_rot_cpds = np.sum(cv_rot_cpds)
        cv_rot_aly = (5/4)*( sc.R*(sc.hbar**3)/( (sc.k**2)*w_aly*((thetaD_aly + T)**2) ) ) * cv_rot_cpds
        cv_rot_aly = (n_aly+0.5)*( 9*cv_rot_aly + ((1-np.sqrt(Ey_aly*rho_C/(Ey_C*rho_aly)))*(sc.R*(T**3)/(thetaD_aly*(Tm_aly)**2))) )
        
        return cv_rot_aly
    
    def cv_debyemod_calc(self, T, thetaD_aly, D_aly, cve_aly, cv_rot_aly, Mmol_aly):
        
        
        i_T = np.ones(len(T))
        for i in range(0, len(T)):
            i_T[i] = i_.quad(lambda x: (x**4)*np.exp(x)/((np.exp(x)-1)**2), 0., thetaD_aly/T[i])[0]
        
        cv_debmod = ( (1+D_aly)*9*sc.R*((T/thetaD_aly)**3)*i_T*(1+cve_aly)+cv_rot_aly ) / Mmol_aly
        
        return cv_debmod
    
    def k_poly_fit_calc(self, T, k_coeffs):
    
        k = k_coeffs[0] + k_coeffs[1]*np.log10(T) + k_coeffs[2]*(np.log10(T))**2 + k_coeffs[3]*(np.log10(T))**3 + k_coeffs[4]*(np.log10(T))**4 +\
                k_coeffs[5]*(np.log10(T))**5 + k_coeffs[6]*(np.log10(T))**6 + k_coeffs[7]*(np.log10(T))**7 + k_coeffs[8]*(np.log10(T))**8
        k = 10**k
        
        return k
    
    def k_poly_fit_calc_cu(self, T, k_coeffs):
    
        k = (k_coeffs[0] + k_coeffs[2]*(np.log10(T))**0.5 + k_coeffs[4]*(np.log10(T))**1 + k_coeffs[6]*(np.log10(T))**1.5 + k_coeffs[8]*(np.log10(T))**2)/\
            (1 + k_coeffs[1]*(np.log10(T))**0.5 + k_coeffs[3]*(np.log10(T))**1 + k_coeffs[5]*(np.log10(T))**1.5 + k_coeffs[7]*(np.log10(T))**2)
        k = 10**k
        
        return k
    
    def c_poly_fit_calc(self, T, c_coeffs):
    
        c = c_coeffs[0] + c_coeffs[1]*np.log10(T) + c_coeffs[2]*(np.log10(T))**2 + c_coeffs[3]*(np.log10(T))**3 + c_coeffs[4]*(np.log10(T))**4 +\
                c_coeffs[5]*(np.log10(T))**5 + c_coeffs[6]*(np.log10(T))**6 + c_coeffs[7]*(np.log10(T))**7 + c_coeffs[8]*(np.log10(T))**8
        c = 10**c
        
        return c
    
    def my_CT(self, xy, z):
        """CT interpolator + nearest-neighbor extrapolation.

        Parameters
        ----------
        xy : ndarray, shape (npoints, ndim)
            Coordinates of data points
        z : ndarray, shape (npoints)
            Values at data points

        Returns
        -------
        func : callable
            A callable object which mirrors the CT behavior,
            with an additional neareast-neighbor extrapolation
            outside of the data range.
        """
        x = xy[:, 0]
        y = xy[:, 1]
        f = interpolate.CloughTocher2DInterpolator(xy, z)

        # this inner function will be returned to a user
        def new_f(xx, yy):
            # evaluate the CT interpolator. Out-of-bounds values are nan.
            zz = f(xx, yy)
            nans = np.isnan(zz)

            if nans.any():
                # for each nan point, find its nearest neighbor
                inds = np.argmin(
                    (x[:, None] - xx[nans])**2 +
                    (y[:, None] - yy[nans])**2
                    , axis=0)
                # ... and use its value
                zz[nans] = z[inds]
            return zz

        return new_f
    
    def ptc_heatflux_calc(self, T1, T2, PW1, PW2):
        
        T1T2 = np.c_[T1, T2]
        lut1 = interpolate.CloughTocher2DInterpolator(T1T2, PW1)
        lut11 = self.my_CT(T1T2, PW1)
        lut2 = interpolate.CloughTocher2DInterpolator(T1T2, PW2)
        lut22 = self.my_CT(T1T2, PW2)
        
        t1new = np.arange(1., 401., 1.)
        t2new = np.arange(1., 401., 1.)
        t1new, t2new = np.meshgrid(t1new, t2new)
        
        pw1_new = lut11(t1new, t2new)
        pw2_new = lut22(t1new, t2new)
        
        return t1new, t2new, pw1_new, pw2_new
    
    # def ptc_heatflux_calc_b(self, T1, T2, PW1, PW2):
    #     # t1, t2 = np.meshgrid(T1, T2)
    #     pw1 = np.meshgrid(PW1, PW1)[1]
    #     pw2 = np.meshgrid(PW2, PW2)[1]
    #     grid_points = [T1, T2]
    #     data_matrix1 = [pw1]
    #     data_matrix2 = [pw2]
    #     s1 = skfda.FDataGrid(data_matrix1, grid_points)
    #     s2 = skfda.FDataGrid(data_matrix2, grid_points)
    #     t1new = np.arange(1., 401., 1.)
    #     t2new = np.arange(1., 401., 1.)
    #     pw1_new = s1((t1new, t2new), grid=True, extrapolation="bounds")
    #     pw2_new = s2((t1new, t2new), grid=True, extrapolation="bounds")
        
    #     return t1new, t2new, pw1_new, pw2_new
    
    def ptc_heatflux_calc_c(self, T1, T2, PW1, PW2):
        
        t1new = np.arange(1., 401., 1.)
        t2new = np.arange(1., 401., 1.)
        t1new, t2new = np.meshgrid(t1new, t2new)
        rbf1 = interpolate.Rbf(T1, T2, PW1, function="multiquadric", smooth=5)
        pw1_new = rbf1(t1new, t2new)
        rbf2 = interpolate.Rbf(T1, T2, PW2, function="multiquadric", smooth=5)
        pw2_new = rbf2(t1new, t2new)
        
        return t1new, t2new, pw1_new, pw2_new
    
    def Cp_Cu_calc(sefl):
        t = np.arange(1., 26., 1.)
        a1 = 6.9434e-1
        a2 = 4.7548e-2
        a3 = 1.6314e-6
        a4 = 9.4786e-8
        a5 = -1.3639e-10
        a6 = 5.3898e-14
        Cp_Cu_25K = (a1*(t**1) + a2*(t**3) + a3*(t**5) + a4*(t**7) + a5*(t**9) + a6*(t**11))/63.54e-3/1e3
        
        t = np.arange(26., 401., 1.)
        a0 = 4.89287
        a1 = -57.51701
        a2 = 238.2039
        a3 = -345.4283
        a4 = 275.8975
        a5 = -132.5425
        a6 = 38.17399
        a7 = -6.07962
        a8 = 0.4118687
        Cp_Cu_400K = (a0*((t/100)**0) + a1*((t/100)**1) + a2*((t/100)**2) + a3*((t/100)**3) + a4*((t/100)**4) + \
        a5*((t/100)**5) + a6*((t/100)**6) + a7*((t/100)**7) + a8*((t/100)**8))/63.54e-3
        
        t = np.arange(1., 401., 1.)
        Cp_Cu = np.concatenate((Cp_Cu_25K, Cp_Cu_400K), axis=None)
        
        return t, Cp_Cu
    
    def mat_k_calc(self, T_in, k_in):
        
        t_new = np.arange(1., 401., 1.)
        # rbf1 = interpolate.CubicSpline(T_in, kappa_in, bc_type='natural')
        rbf1 = interpolate.interp1d(T_in, k_in, kind='slinear', fill_value='extrapolate')
        k_new = rbf1(t_new)
        k_new = savgol_filter(k_new, 53, 3)
        
        return t_new, k_new
    
    def mat_kappa_calc(self, T_in, kappa_in):
        
        t_new = np.arange(1., 401., 1.)
        # rbf1 = interpolate.CubicSpline(T_in, kappa_in, bc_type='natural')
        rbf1 = interpolate.interp1d(T_in, kappa_in, kind='slinear', fill_value='extrapolate')
        kappa_new = rbf1(t_new)
        kappa_new = savgol_filter(kappa_new, 53, 3)
        
        return t_new, kappa_new
    
    def mat_eps_calc(self, T_in, eps_in, pct = True):
        
        t_new = np.arange(1., 401., 1.)
        # rbf1 = interpolate.CubicSpline(T_in, kappa_in, bc_type='natural')
        rbf1 = interpolate.interp1d(T_in, eps_in, kind='slinear', fill_value='extrapolate')
        if pct:
            eps_new = rbf1(t_new)/100.
        else:
            eps_new = rbf1(t_new)
        eps_new = savgol_filter(eps_new, 53, 3)
        
        return t_new, eps_new
    
    def mat_eps_analytical_calc(self, T, k):
        
        res_mat = (((np.pi)**2)/3)*((sc.Boltzmann/sc.e)**2)*T/k
        eps_new = 0.751*(res_mat*T)**(1/2) - 0.632*(res_mat*T) + 0.67*(res_mat*T)**(3/2)-0.607*(res_mat*T)**(2)
        
        return res_mat, eps_new
    
    
    # def mat_props_calc(self, material, T):
    #     """        
    #     Compile the thermal conductivity k of a selected material as a function of temperature.
        
    #     :param material: Material of interest. Must be one of the following:
    #     "Al_6063" (T5), "Al_6061" (T6), "Al_1100"
    #     :type material: string
    #     :return k: Thermal conductivity array, function of T [W m-1 K-1]
    #     :return C: Specific heat capacity array, function of T [J kg-1 K-1]
    #     :return rho: Density, function of T [kg m-3]
    #     :rtype k: array of floats
    #     :rtype C: array of floats
    #     :rtype rho: float
        
    #     .. note::   Invalid below 4K
    #                 TO BE IMPLEMENTED: "Cu_OFHC", "G10_CR", "Brass", "Stainless_304", "Stainless_304L", "Stainless_310", "Stainless_316"
    #                 Superconducting state for specific heat capacity C and thermal conductivity k...?
        
    #     """
        
    #     rho_C = 2200 # kg.m^-3

    #     if material == r"Al_6063":
            
    #         rho = 2700 # Doesn't vary much with CTE, maybe 1% max
            
    #         a = 22.401433; b = -141.13433; c = 394.95461; d = -601.15377; e = 547.83202; f = -305.99691; g = 102.38656; h = -18.810237; i = 1.4576882; 
    #         # Debye temperatures of the different compounds
    #         ThetaD_Al = 394.; ThetaD_Mg = 318.; ThetaD_Si = 625.
            
    #         # Molar fractions and Debye temperatures
    #         x_Si = 0.004 # To be modified if needed
    #         x_Mg = 0.007 # To be modified if needed
    #         x_Al = 1. - (x_Mg + x_Si) 
    #         # min x_Si = 0.002, max = 0.006 | min x_Mg = 0.0045, max = 0.009
            
    #         # Young modulus of the compounds
    #         Ey_Al = 70e9; Ey_Mg = 42e9; Ey_Si = 168.9e9; Ey_C = 1220e9
    #         Ey_aly = x_Al*Ey_Al + x_Mg*Ey_Mg + x_Si*Ey_Si
            
    #         # Debye temperature of the alloy
    #         ThetaD = x_Al*ThetaD_Al + x_Mg*ThetaD_Mg + x_Si*ThetaD_Si
            
    #         # Melting temperatures of the different compounds
    #         Tm_Al = 660.+273.; Tm_Mg = 650.+273.; Tm_Si = 1410.+273.
    #         Tm_aly = 635.+273.
            
    #         # Valence of the different compounds
    #         z_Al = 3.; z_Mg = 2.; z_Si = 4.
            
    #         # Electronic molar contributuion to C
    #         cve_Al = (5/(24*(np.pi**3)))*z_Al*(ThetaD_Al**3)/((T**2)*Tm_Al)
    #         cve_Mg = (5/(24*(np.pi**3)))*z_Mg*(ThetaD_Mg**3)/((T**2)*Tm_Mg)
    #         cve_Si = (5/(24*(np.pi**3)))*z_Si*(ThetaD_Si**3)/((T**2)*Tm_Si)
    #         cve_aly = x_Al*cve_Al + x_Mg*cve_Mg + x_Si*cve_Si +\
    #             x_Al*cve_Al*x_Mg*cve_Mg + x_Al*cve_Al*x_Si*cve_Si + x_Mg*cve_Mg*x_Si*cve_Si +\
    #             x_Al*cve_Al*x_Mg*cve_Mg*x_Si*cve_Si
            
    #         # Atomic Radii/Volumes of the compounds
    #         r_Al = 118e-12; r_Mg = 173e-12; r_Si = 210e-12
    #         V_Al = (4/3)*np.pi*(r_Al**3)
            
    #         # Alloy Debye frequency (max phonon freq that can exist in the compound's lattice)
    #         w_aly = sc.k*ThetaD/sc.hbar
            
    #         # Sound speed in compounds (assumed constant in metal as density variation with temp is below 1%)
    #         v_Al = 3.04e3; v_Mg = 3.05e3; v_Si = 5.843e3; v_aly = x_Al*v_Al + x_Mg*v_Mg + x_Si*v_Si
            
    #         # Density of state of the compounds ('The number of different states at a particular energy level that electrons are allowed to occupy')
    #         D_aly = V_Al*(w_aly**2)/(2*(np.pi**2)*(v_aly**3)) # Using Al volume as it dominates, should use Gibbs Thomson relation...
            
    #         # Density of atoms in a grain of compound
    #         n_aly = (1/(6*(np.pi**2)))*((sc.k*ThetaD/(sc.hbar*w_aly))**3)
            
    #         # Rotational contribution to C
    #         Mmol_Al = 26.982e-3; Mmol_Mg = 24.305e-3; Mmol_Si = 28.086e-3
            
    #         cv_rot_Al = 0.; cv_rot_Mg = 0.; cv_rot_Si = 0.
    #         for j in range(1,int(n_aly)+2):
    #             cv_rot_Al += x_Al*j*(j+1)/(Mmol_Al*(r_Al**2))
    #             cv_rot_Mg += x_Mg*j*(j+1)/(Mmol_Mg*(r_Mg**2))
    #             cv_rot_Si += x_Si*j*(j+1)/(Mmol_Si*(r_Si**2))
            
    #         cv_rot = (5/4)*( sc.R*(sc.hbar**3)/( (sc.k**2)*w_aly*((ThetaD + T)**2) ) ) *\
    #             (cv_rot_Al + cv_rot_Mg + cv_rot_Si)
            
    #         # C = (12*(np.pi**4)/5)*sc.R*((T/ThetaD)**3) # Debye model, innacurate for compounds; see: https://www.scielo.br/j/mr/a/4vdpvsfLD96RJXHTQW3zsjd/?lang=en#
    #         i_T = np.ones(len(T))
    #         # C = np.ones(len(T))
    #         for i in range(0, len(i_T)):
    #             # if T[i] < ThetaD:
    #                 i_T[i] = i_.quad(lambda x: (x**4)*np.exp(x)/((np.exp(x)-1)**2), 0., ThetaD/T[i])[0]
                    
    #         c_rot = (n_aly+0.5)*( 9*cv_rot + ((1-np.sqrt(Ey_aly*rho_C/(Ey_C*rho)))*(sc.R*(T**3)/(ThetaD*(Tm_aly)**2))) )
            
            
    #         C = (1+D_aly)*9*sc.R*((T/ThetaD)**3)*i_T*(1+cve_aly) + c_rot
                    
    #                 # C_slope = (C[i] - C[i-1])/(T[i] - T[i-1])
                
    #             # else:
    #                 # C[i] = C_slope*(T[i] - T[i-1]) + C[i-1]
            
    #     if material == r"Al_6061": 
    #         a = 0.07918; b = 1.0957; c = -0.07277; d = 0.08084; e = 0.02803; f = -0.09464; g = 0.04179; h = -0.00571; i = 0; 
    #         C = 896*T # superconductivity state?!
    #         rho = 2700
        
    #     if material == r"Al_1100": 
    #         a = 23.39172; b = -148.5733; c = 422.1917; d = -653.6664; e = 607.0402; f = -346.152; g = 118.4276; h = -22.2781; i = 1.770187; 
    #         C = 904*T # superconductivity state?!
    #         rho = 2710

    #     k = a + b*np.log10(T) + c*(np.log10(T))**2 + d*(np.log10(T))**3 + e*(np.log10(T))**4 + f*(np.log10(T))**5 + g*(np.log10(T))**6 + h*(np.log10(T))**7 + i*(np.log10(T))**8
    #     k = 10**k
        
    #     return T, ThetaD, i_T, k, C, rho, cv_rot, c_rot, cve_aly