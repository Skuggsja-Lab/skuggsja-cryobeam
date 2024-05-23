##########################################
# PREAMBLE
##########################################

import numpy as np
import thermal_props
tp = thermal_props.thermalProps()

##########################################
# THERMAL FINITE DOMAIN METHOD (FDM)
##########################################

class thermalFDM:
    """
    .. versionadded:: 0.1
    
    Baseline heat transfer solver for various geometries
    
    """

    def rod_1d(self, mat, Ar, Lzr, Dz, Dt, DT, Ti, T1, T2, Qp1, Qp2):
        """        
        Compile the conductive heat transfer of a 1-dimensionnal (z-axis) rod (e.g. heat straps) over a finite period of time.
        
        :param mat: Material constituting the rod
        :param Ar: Area of the rod cross-section [m2]
        :param Lzr: Length of the rod [m]
        :param Dz: Geometrical step size [m]
        :param Dt: Time span [s]
        :param DT: Temperature step size [K]
        :param Ti: Temperature array at initial time [K]
        :param T1: Forced constant temperature at z = 0 [K]
        :param T2: Forced constant temperature at z = Lzr [K]
        :param Qp1: Constant heat flow input (from PTC) at z = 0
        :param Qp2: Constant heat flow input (from PTC) at z = Lzf
        :type mat: string
        :type Ar: float
        :type Lzr: float
        :type Dz: float
        :type Dt: float
        :type DT: float
        :type Ti: array of floats
        :type T1: float
        :type T2: float
        :type Qp1: float
        :type Qp2: float
        :return Tf: Temperature array after Dt
        :rtype Tf: float
        
        .. note:: None
        
        """

        # Dirichlet Condition (Constant T at boundary)
        if T1 != 0: Ti[0,0] = T1
        if T2 != 0: Ti[0,-1] = T2

        T, k, C, rhor = tp.mat_props_calc(mat, np.min(Ti), np.max(Ti), DT)

        zz = np.arange(0, Lzr + Dz, Dz)
        # Get value for Cr and kr at the selected temperatures from the input databases
        # Get value of Qp as a function of the temperature of the rod at connection

        l = np.zeros((1,len(zz)))

        A = np.zeros((len(zz), len(zz)))
        for i in range(0, len(zz)):
            j = np.where(T == Ti[0,i]) # make sure Ti has common values with T... to be changed later
            if i < len(zz) - 1: 
                ju = np.where(T == Ti[0,i+1])
                kr_upper = (k[ju] + k[j])/2
                kr_upper = kr_upper[0]*1e-3
            if i != 0: 
                jl = np.where(T == Ti[0,i-1])
                kr_lower = (k[jl] + k[j])/2
                kr_lower = kr_lower[0]*1e-3
            if i == 0: kr_lower = 0.
            if i == len(zz) - 1: kr_upper = 0.
            kr = k[j]
            kr = kr[0]*1e-3
            Cr = C[j]

            Aii = -(Dt)*(kr_lower + kr)/(2*Cr*rhor*(Dz)**2)
            Bii = 1 + (Dt)*(kr_lower + 2*kr + kr_upper)/(2*Cr*rhor*(Dz)**2)
            Cii = -(Dt)*(kr + kr_upper)/(2*Cr*rhor*(Dz)**2)

            if i == 0:
                A[i,i] = Bii
                A[i,i+1] = Cii
            elif i == len(zz) - 1:
                A[i,i-1] = Aii
                A[i,i] = Bii
            else:
                A[i,i-1] = Aii
                A[i,i] = Bii
                A[i,i+1] = Cii
        
        # Dirichlet Condition (Constant T at boundary)
        if T1 != 0: A[0,0] += 1; print(r"Dirichlet - T1")#Ti[0] = T1
        if T2 != 0: A[-1,-1] += 1; print(r"Dirichlet - T2") #Ti[-1] = T2

        # Neumann Condition (Constant Qp at boundary, PTC for instance)
        if Qp1 != 0:
            print(r"Neumann - Qp1")
            if len(k) == 1:
                kr = k[0]*1e-3
                A[0,0] += 1 + Dt*(2*kr)/(2*Cr*rhor*(Dz)**2)
                A[0,1] += - Dt*(2*kr)/(2*Cr*rhor*(Dz)**2)
                l[0] += -Qp1/(Dz*Ar*Cr*rhor)
            else:
                A[0,0] += 1 + Dt*(kr+k[1])/(2*Cr*rhor*(Dz)**2)
                A[0,1] += - Dt*(kr+k[1])/(2*Cr*rhor*(Dz)**2)
                l[0] += -Qp1/(Dz*Ar*Cr*rhor)
        if Qp2 != 0:
            print(r"Neumann - Qp2")
            if len(k) == 1:
                kr = k[0]*1e-3
                A[-1,-1] += 1 + Dt*(2*kr)/(2*Cr*rhor*(Dz)**2)
                A[-1,-2] += - Dt*(2*kr)/(2*Cr*rhor*(Dz)**2)
                l[-1] += -Qp2/(Dz*Ar*Cr*rhor) # again Cr should vary with temperature really
            else:
                A[-1,-1] += 1 + Dt*(k[-1]+k[-2])/(2*Cr*rhor*(Dz)**2)
                A[-1,-2] += - Dt*(k[-1]+k[-2])/(2*Cr*rhor*(Dz)**2)
                l[-1] += -Qp2/(Dz*Ar*Cr*rhor) # again Cr should vary with temperature really

        Tf = np.matmul(Ti,np.linalg.inv(A)) + np.matmul(l,np.linalg.inv(A))

        return Tf
    
    def test_kC(self, mat, Ti, DT):
        
        T, k, C, rhor = tp.mat_props_calc(mat, np.min(Ti), np.max(Ti), DT)
        
        return T, k, C, rhor
    
    def rod_1d_grid(self, mat, DT, Ti, Dz, Dt, Lz, tf):
        
        T, k, C, rhor = tp.mat_props_calc(mat, np.min(Ti), np.max(Ti), DT)
        
        Nz = int(Lz/Dz)
        Nt = int(tf/Dt)
        time=np.arange(0,(tf+Dt),Dt)
        z=np.arange(0,Lz+Dz,Dz)
        zGrid, tGrid = np.meshgrid(z, time)
        
        return zGrid, tGrid
    
    def rod_1d_matx(self, theta, mat, DT, Ti, Dz, Dt, Lz, tf, T1, T2):
        
        T, k, C, rhor = tp.mat_props_calc(mat, np.min(Ti), np.max(Ti), DT)
        
        Nz = int(Lz/Dz)
        Nt = int(tf/Dt)
        time=np.arange(0,(tf+Dt),Dt)
        z=np.arange(0,Lz+Dz,Dz)
        
        Tf=np.ones((Nz+1,Nt+1))
        b=np.zeros(Nz-1)
        
        # Boundary Condition (Dirichlet)
        for j in range (0,Nt+1):
            Tf[:,j]=Ti
            Tf[0,j]=T1
            Tf[Nz,j]=T2
        
        A=np.zeros((Nz-1,Nz-1))
        B=np.zeros((Nz-1,Nz-1))
        for i in range (0,Nz-1):
            j = np.where(T == Tf[i,0])
            kr = k[j]
            kr = kr*1e-3
            Cr = C[j]
            r = Dt*kr/(2*Cr*rhor*(Dz**2))
            print(T, Ti)
            A[i,i]=-(1+2*r*theta)
            B[i,i]=1-2*r*(1-theta)

        for i in range (0,Nz-2):
            j = np.where(T == Tf[i+1,0])
            krp1 = k[j]
            krp1 = kr[0]*1e-3
            Crp1 = C[j]
            rp1 = Dt*krp1/(2*Crp1*rhor*(Dz**2))
            if i != 0:
                j = np.where(T == Tf[i-1,0])
                krm1 = k[j]
                krm1 = kr*1e-3
                Crm1 = C[j]
                rm1 = Dt*krm1/(2*Crm1*rhor*(Dz**2))
            else:
                j = np.where(T == Tf[i,0])
                krm1 = k[j]
                krm1 = kr*1e-3
                Crm1 = C[j]
                rm1 = Dt*krm1/(2*Crm1*rhor*(Dz**2))
            A[i+1,i]=rm1*theta
            B[i+1,i]=rm1*(1-theta)
            A[i,i+1]=rp1*theta
            B[i,i+1]=rp1*(1-theta)
            
        Ainv=np.linalg.inv(A)  
        
        return Ainv, A, B, Tf
    
    def rod_1d_solve(self, theta, mat, DT, Ti, Dz, Dt, Lz, tf, T1, T2):
        
        Nz = int(Lz/Dz)
        Nt = int(tf/Dt)
        time=np.arange(0,(tf+Dt),Dt)
        z=np.arange(0,Lz+Dz,Dz)
        
        T, k, C, rhor = tp.mat_props_calc(mat, np.min(Ti), np.max(Ti), DT)
        
        Ainv, A, B, Tf = self.rod_1d_matx(theta, mat, DT, Ti, Dz, Dt, Lz, tf, T1, T2)
        
        l=np.zeros(Nz-1)
        
        for j in range (1,Nt+1):
            # At 0
            jj = np.where(T == Tf[0,j])
            kr = k[jj]
            kr = kr*1e-3
            Cr = C[jj]
            r = Dt*kr/(2*Cr*rhor*(Dz**2))
            l[0]=(1-2*r*(1-theta))*Tf[0,j-1] + (1+2*r*theta)*Tf[0,j]
            
            # At Nz
            jj = np.where(T == Tf[Nz,j])
            kr = k[jj]
            kr = kr*1e-3
            Cr = C[jj]
            r = Dt*kr/(2*Cr*rhor*(Dz**2))
            l[Nz-2]=(1-2*r*(1-theta))*Tf[Nz,j-1] + (1+2*r*theta)*Tf[Nz,j]
            
            v=np.dot(B,Tf[1:(Nz),j-1])
            Tf[1:(Nz),j]=np.dot(Ainv,v+l)
        
        return Tf, l