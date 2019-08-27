import h5py as h5
import numpy as np
import os

"""
SET PARAMETERS
"""
path = "./init/"
filename = path +"Radiation_Homogeneous.3d.hdf5"  # Name of the new file to create

EXTENDED = False
N = 32
L = 2.17080
dt_multiplier = 0.01

omega = 0.33333333
mass = 0.  # 1 or 0

def transform_PhiToChi(x): return np.exp(-2*x)

# Set components
if EXTENDED: 
    filename = path +"InitialDataPerfectFluidExt.3d.hdf5"  # Name of the new file to create
    component_names = [  # The order is important: component_0 ... component_(nth-1)
        "chi",

        "h11",    "h12",    "h13",    "h22", "h23", "h33",

        "K",

        "A11",    "A12",    "A13",    "A22", "A23", "A33",

        "Theta",

        "Gamma1", "Gamma2", "Gamma3",

        "lapse",

        "shift1", "shift2", "shift3",

        "B1",     "B2",     "B3",

        "density",  "energy", "pressure", "enthalpy",

        "u0", "u1", "u2", "u3",

        "D",  "E", "W",

        "Z0", "Z1", "Z2", "Z3",

        "V1", "V2","V3",

        "Ham",

        "Ham_ricci", "Ham_trA2", "Ham_K", "Ham_rho",  

        "Mom1",   "Mom2",   "Mom3"
    ]
else: 
    component_names = [  # The order is important: component_0 ... component_(nth-1)
        "chi",

        "h11",    "h12",    "h13",    "h22", "h23", "h33",

        "K",

        "A11",    "A12",    "A13",    "A22", "A23", "A33",

        "Theta",

        "Gamma1", "Gamma2", "Gamma3",

        "lapse",

        "shift1", "shift2", "shift3",

        "B1",     "B2",     "B3",

        "density",  "energy", "pressure", "enthalpy",

        "u0", "u1", "u2", "u3",

        "D",  "E", "W",

        "Z0", "Z1", "Z2", "Z3",

        "V1", "V2","V3",

        "Ham",

    #     "Ham_ricci", "Ham_trA2", "Ham_K", "Ham_rho",  // Extended!

        "Mom1",   "Mom2",   "Mom3"
    ]

temp_comp = np.zeros((N, N, N))   # template for components: array [Nx, Ny. Nz]
dset = dict()
# Here set the value of the components (default: to zero)

dset['chi'] = temp_comp.copy() + 1.
dset['Ham'] = temp_comp.copy()
dset['h11'] = temp_comp.copy() + 1.
dset['h22'] = temp_comp.copy() + 1.
dset['h33'] = temp_comp.copy() + 1.
dset['lapse'] = temp_comp.copy() + 1.

dset['D'] = temp_comp.copy()
dset['E'] = temp_comp.copy()
dset['density'] = temp_comp.copy()
dset['energy'] = temp_comp.copy()
dset['pressure'] = temp_comp.copy()
dset['enthalpy'] = temp_comp.copy()
dset['Z0'] = temp_comp.copy()
dset['u0'] = temp_comp.copy() + 1.
dset['W'] = temp_comp.copy() + 1.

rho_emtensor = temp_comp.copy()

# ## Constructing variables 
indices = []
for z in range(N):
    for y in range(N):
        for x in range(N):
            #wvl = 2 * np.pi * 4 / L
            ind = x + y*N + z*N**2 
            
            dset['chi'][x][y][z] = 1.
            dset['density'][x][y][z] = 0.1
            dset['energy'][x][y][z] = 3 * (dset['density'][x][y][z])**(1/3)
            dset['pressure'][x][y][z] =  omega * dset['density'][x][y][z] * dset['energy'][x][y][z]
            dset['enthalpy'][x][y][z] =  mass +  dset['energy'][x][y][z] + omega
            dset['D'][x][y][z] = dset['density'][x][y][z]
            dset['E'][x][y][z] = dset['energy'][x][y][z] * dset['density'][x][y][z]
            
            
            rho_emtensor[x][y][z] = dset['density'][x][y][z] * dset['enthalpy'][x][y][z] * dset['u0'][x][y][z] * dset['u0'][x][y][z] - dset['pressure'][x][y][z]

            indices.append(ind)
            
            
# Constructing variables  (reading from file)  # Commented out as you don't have the files
# fldata = np.loadtxt(path + 'gridofsf000.dat')
# logconfdata = np.loadtxt(path + 'gridofphi000.dat')
# indices = []
# for z in range(N):
#     for y in range(N):
#         for x in range(N):
#             ind = x + y*N + z*N**2 
#             dset['D'][x][y][z] = fldata[ind]
#             dset['density'][x][y][z] = fldata[ind]
#             dset['chi'][x][y][z] = transform_PhiToChi(logconfdata[ind]) 
#             indices.append(ind)
            
rho_mean = np.mean( np.hstack(rho_emtensor))    
K_const = - np.sqrt(rho_mean * 24 * np.pi )
dset['K'] = temp_comp.copy() + K_const

print("checks")
print("rho mean:",  rho_mean)
print("K mean:",  np.mean(dset['K']), "   ", np.sqrt(np.mean(dset['D']) *24*np.pi) )
print("E mean:",  np.mean(dset['E']))
print("D mean:",  np.mean(dset['D']))

print("Prassure mean:",  np.mean(dset['pressure']))
# print("phi mean:",  np.mean( np.hstack(dset['phi'])))





"""

CREATE HDF5

"""

if not os.path.exists(path):
    os.mkdir(path)
    print(" ! > new mkdir: ", path)

"""
Mesh and Other Params
"""
# def base attributes
base_attrb = dict()
base_attrb['time'] = 0.0
base_attrb['iteration'] = 0
base_attrb['max_level'] = 0
base_attrb['num_components'] = len(component_names)
base_attrb['num_levels'] = 1
base_attrb['regrid_interval_0'] = 1
base_attrb['steps_since_regrid_0'] = 0
for comp,  name in enumerate(component_names):
    key = 'component_' + str(comp)
    tt = 'S' + str(len(name))
    base_attrb[key] = np.array(name, dtype=tt)


# def Chombo_global attributes
chombogloba_attrb = dict()
chombogloba_attrb['testReal'] = 0.0
chombogloba_attrb['SpaceDim'] = 3

# def level0 attributes
level_attrb = dict()
level_attrb['dt'] = float(L)/N * dt_multiplier
level_attrb['dx'] = float(L)/N
level_attrb['time'] = 0.0
level_attrb['is_periodic_0'] = 1
level_attrb['is_periodic_1'] = 1
level_attrb['is_periodic_2'] = 1
level_attrb['ref_ratio']= 2
level_attrb['tag_buffer_size'] = 3
prob_dom = (0, 0, 0, N-1, N-1, N-1)
prob_dt = np.dtype([('lo_i', '<i4'), ('lo_j', '<i4'), ('lo_k', '<i4'),
                    ('hi_i', '<i4'), ('hi_j', '<i4'), ('hi_k', '<i4')])
level_attrb['prob_domain'] = np.array(prob_dom, dtype=prob_dt)
boxes = np.array([(0, 0, 0, N-1, N-1, N-1)],
      dtype=[('lo_i', '<i4'), ('lo_j', '<i4'), ('lo_k', '<i4'), ('hi_i', '<i4'), ('hi_j', '<i4'), ('hi_k', '<i4')])


""""
CREATE HDF5
"""

#TODO: if overwrite:   [...] else: raise()
if os.path.exists(filename):
    os.remove(filename)

h5file = h5.File(filename, 'w')  # New hdf5 file I want to create

# base attributes
for key in base_attrb.keys():
    h5file.attrs[key] = base_attrb[key]

# group: Chombo_global
chg = h5file.create_group('Chombo_global')
for key in chombogloba_attrb.keys():
    chg.attrs[key] = chombogloba_attrb[key]

# group: levels
l0 = h5file.create_group('level_0')
for key in level_attrb.keys():
    l0.attrs[key] = level_attrb[key]
sl0 = l0.create_group('data_attributes')
dadt = np.dtype([('intvecti', '<i4'), ('intvectj', '<i4'), ('intvectk', '<i4')])
sl0.attrs['ghost'] = np.array((3, 3, 3),  dtype=dadt)
sl0.attrs['outputGhost'] = np.array( (0, 0, 0),  dtype=dadt)
sl0.attrs['comps'] = base_attrb['num_components']
sl0.attrs['objectType'] = np.array('FArrayBox', dtype='S10')

# level datasets
dataset = np.zeros((base_attrb['num_components'], N, N, N))
for i, comp in enumerate(component_names):
    if comp in dset.keys():
        dataset[i] = dset[comp].T
fdset = []
for c in range(base_attrb['num_components']):
    fc = dataset[c].T.flatten()
    fdset.extend(fc)
fdset = np.array(fdset)

l0.create_dataset("Processors", data=np.array([0]))
l0.create_dataset("boxes",  data=boxes)
l0.create_dataset("data:offsets=0",  data=np.array([0, (base_attrb['num_components'])*N**3]))
l0.create_dataset("data:datatype=0",  data=fdset)

h5file.close()

