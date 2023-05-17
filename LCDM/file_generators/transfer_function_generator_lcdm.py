import numpy as np
import classy
from classy import Class
import scipy
import scipy.interpolate
from scipy.optimize import fsolve
from scipy.interpolate import interp1d
import math
import os
import sys

if len(sys.argv) < 3:
    print(f'Usage: python {sys.argv[0]} nmin nmax [--test]')
    print('    nmin: index of cosmology in LHS to begin generating')
    print('    nmax: index of cosmology in LHS to end generating')
    print('    --ref: Flag for generating transfer files for ref cosmology')
    exit(0)

nmin = int(sys.argv[1])
nmax = int(sys.argv[2])
ref = False
if len(sys.argv) == 4:
    ref = True

print('----------------------------------------------------------------')
print(f'Generating transfer functions for cosmologies {nmin} to {nmax}.')

# Always print which CLASS is being used!
print('Using CLASS installed in:', classy.__file__)

# Loading LHS points
lhs_path = '/home/grads/data/jonathan/cola_projects/LCDM_800/lcdm_train_800.txt'
#lhs_path = '/home/grads/data/jonathan/cola_projects/5D_emulator1/lhs.txt'
print(f'Reading LHS in {lhs_path}')
Omegam_lhs, Omegab_lhs, ns_lhs, As_lhs, h_lhs = np.loadtxt(lhs_path, unpack=True)

#############################################
# User settings controlling the figure aspect
#############################################

kvec = [float("{:g}".format(kk)) for kk in np.logspace(-5,0,5)]
z_max_pk         = 300     # highest redshift involved
k_per_decade     = 32      # number of k values, controls final resolution
k_min_tau0       = 40.     # this value controls the minimum k value in the figure (it is k_min * tau0)
P_k_max_h_Mpc    = 8.0     # this value is directly the maximum k value in the figure in Mpc
tau_num_early    = 2000    # number of conformal time values before recombination, controls final resolution
tau_num_late     = 200     # number of conformal time values after recombination, controls final resolution
tau_ini          = 10.     # first value of conformal time in Mpc
tau_label_Hubble = 20.     # value of time at which we want to place the label on Hubble crossing
tau_label_ks     = 40.     # value of time at which we want to place the label on sound horizon crossing
tau_label_kd     = 230.    # value of time at which we want to place the label on damping scale crossing
cosmo = Class()

print('Starting transfer function generation.')
for n_iter in range(nmin, nmax+1):
    if ref:
        h                = 0.67
        As               = 2.1e-9
        ns               = 0.96
        Omega_b          = 0.049
        Omega_m          = 0.319
        Omega_massivenu  = 0.058/(93.15*h**2) # JVR MODIFICATION - correct Omega_massivenu to equation (2.84) in Dodelson, Schmidt
        Omega_c          = Omega_m - Omega_b - Omega_massivenu
        #filename = f'ref_{n_iter}'
        filename = 'ref'
    else:
        # Reading cosmo parameters in LHS
        h                = h_lhs[n_iter]
        As               = As_lhs[n_iter]
        ns               = ns_lhs[n_iter]
        Omega_b          = Omegab_lhs[n_iter]
        Omega_m          = Omegam_lhs[n_iter]
        Omega_massivenu  = 0.058/(93.15*h**2) # JVR MODIFICATION - correct Omega_massivenu to equation (2.84) in Dodelson, Schmidt
        Omega_c          = Omega_m - Omega_b - Omega_massivenu
        #filename = f'lcdm_{n_iter}'
        filename = f'{n_iter}'
    
    print(f'Running cosmology {n_iter}: As = {As}, ns = {ns}, Omegab = {Omega_b}, Omegam = {Omega_m}, h = {h}.')

    # Initializing CLASS
    cosmo.empty()
    cosmo.set(
             {'output':'dTk, vTk, lTk',
              'Omega_Lambda': 0,
              'Omega_fld': 0,
              'Omega_smg': -1,
              #'a_init_nbody': 1./70.,
              'expansion_model': 'lcdm',
              'gravity_model': 'propto_omega',
              'parameters_smg':'1., 1e-9, 0., 0., 1.',
              'a_min_stability_test_smg': 1e-7,
              'radiation_streaming_approximation':3,
              'evolver': 0,
              'tol_thermo_integration':1e-6,
              'ur_fluid_approximation':2,
              'A_s':As,
              'n_s':ns,
              'z_max_pk':z_max_pk,
              'P_k_max_h/Mpc' : P_k_max_h_Mpc,
              'k_per_decade_for_pk':k_per_decade,
              'T_cmb': 2.7255,
              'background_verbose': 0,
              'perturbations_verbose': 0,
              'gauge':'Synchronous',
              #'back_integration_stepsize':1e-4,
              'perturb_sampling_stepsize':'0.25',
              'l_max_ncdm': 500,  # When creating the transfer functions for our runs uncomment these l_maxs!!!!
              'l_max_g':3500, # When creating the transfer functions for our runs uncomment these l_maxs!!!!
              'l_max_pol_g': 3500, # When creating the transfer functions for our runs uncomment these l_maxs!!!!
              'l_max_ur':3500, # When creating the transfer functions for our runs uncomment these l_maxs!!!!
              'tau_reio':0.07, 
              'h': h,
              'Omega_b':Omega_b,
              #'N_ur':3.046, 
              'N_ur':0.,#0.046,#0.00641,
              'N_ncdm':1,
              'm_ncdm':0.058/3.,
              'deg_ncdm':3.0,
              'Omega_cdm': Omega_c,
              'ncdm_fluid_approximation':3})#,
    cosmo.compute()

    background = cosmo.get_background() # load background table

    tau = np.linspace(1000,background['conf. time [Mpc]'][-1],1000)
    tau_num = len(tau)

    background_tau = background['conf. time [Mpc]'] # read conformal times in background table
    background_z = background['z'] # read redshift
    background_a = 1./(1.+background_z) # read redshift
    background_H = background['H [1/Mpc]']

    background_z_at_tau = interp1d(background_tau,background_z,fill_value="extrapolate")
    background_tau_at_z = interp1d(background_z,background_tau,fill_value="extrapolate")
    
    # Saving cosmology params
    Omega_massivenu_exact = background['(.)rho_ncdm[0]'][-1]/background['(.)rho_crit'][-1]
    Omega_c_exact = background['(.)rho_cdm'][-1]/background['(.)rho_crit'][-1]
    Omega_b_exact = background['(.)rho_b'][-1]/background['(.)rho_crit'][-1]
    Omega_m_exact = Omega_b_exact + Omega_c_exact + Omega_massivenu_exact
    #with open(f'cosmo_params/cosmo_{filename}.txt', 'w') as f:
    with open(f'/home/grads/data/jonathan/cola_projects/LCDM_800/cosmo_params/params{filename}.txt', 'w') as f:
        f.write("# ORDER: Omega_nu, Omega_c, Omega_b, Omega_m, As, ns, h\n")
        f.write(f'{Omega_massivenu_exact}\n')
        f.write(f'{Omega_c_exact}\n')
        f.write(f'{Omega_b_exact}\n')
        f.write(f'{Omega_m_exact}\n')
        f.write(f'{As}\n')
        f.write(f'{ns}\n')
        f.write(f'{h}')
    #print(f'Saved cosmological params in cosmo_params_lcdm/cosmo_{filename}.txt')

    ##### BEGIN Gui's snippet to save background for COLA to read
    H0 = cosmo.h() * 100
    c = H0 / background['H [1/Mpc]'][-1]

    H2_over_H02_COLA = background_H * background_H * c * c / H0 / H0

    H_of_a = interp1d(background_a, background_H,fill_value="extrapolate")
    dH_da_class = np.gradient(background_H, background_a)
    dlnH_dlna_clas = dH_da_class * background_a / background_H
    dlnH_dlna_clas_of_a = interp1d(background_a, dlnH_dlna_clas,fill_value="extrapolate")

    H_over_H0_COLA = np.sqrt(H2_over_H02_COLA)

    dH_da = np.gradient(H_over_H0_COLA, background_a)
    dlnH_dlna = dH_da * background_a / H_over_H0_COLA

    H_over_H0_COLA_of_a = interp1d(background_a, H_over_H0_COLA,fill_value="extrapolate")
    dlnH_dlna_of_a = interp1d(background_a, dlnH_dlna,fill_value="extrapolate")

    H_over_H0_CLASS = background_H/background_H[-1]

    as_bg = np.logspace(np.log10(10**(-4)), np.log10(10), 2000)
    as_test = np.logspace(np.log10(10**(-4)), np.log10(1), 2000)
    bg_table = np.vstack([as_test, H_over_H0_COLA_of_a(as_test), dlnH_dlna_of_a(as_test)]).T
    np.savetxt(f'/home/grads/data/jonathan/cola_projects/LCDM_800/backgrounds/bg{filename}.txt',bg_table, fmt='   %.7E   %.7E   %.7E', header = '       a       HoverH0(a)       dlogHdloga_of_a(a)')
    ##### END Gui's snippet to save background for COLA to read
    
    
    max_z_needed = background_z_at_tau(tau[0])
    if max_z_needed > z_max_pk:
        print(f'In cosmology {filename}, you must increase the value of z_max_pk to at least {max_z_needed}')
        exit(1)

    for i in range(tau_num):
        this_z = background_z_at_tau(tau[i])
        if this_z < 0:
            print(f"WARNING: interpolation error led to negative redshift {this_z} at tau = {tau[i]}. Correcting to z = 0.")
            this_z = 0
        
        one_time = cosmo.get_transfer(this_z) # transfer functions at each time tau
        if i == 0:   # if this is the first time in the loop: create the arrays (k, Theta0, phi)
            k = one_time['k (h/Mpc)']*cosmo.h()
            k_num = len(k)
            
            d_GR_CAMB = np.zeros((tau_num,k_num))
            d_m_Nb_CAMB = np.zeros((tau_num,k_num))
            d_cdm_Nb_CAMB = np.zeros((tau_num,k_num))
            d_b_Nb_CAMB = np.zeros((tau_num,k_num))
            d_g_Nb_CAMB = np.zeros((tau_num,k_num))
            d_ur_Nb_CAMB = np.zeros((tau_num,k_num))
            d_total_Nb_CAMB = np.zeros((tau_num,k_num))
            d_total_DE_Nb_CAMB = np.zeros((tau_num,k_num))


        d_GR_CAMB[i,:] = one_time['delta_GR_CAMB'][:]
        d_m_Nb_CAMB[i,:] = one_time['delta_m_Nb_CAMB'][:]
        d_cdm_Nb_CAMB[i,:] = one_time['delta_c_Nb_CAMB'][:]
        d_b_Nb_CAMB[i,:] = one_time['delta_b_Nb_CAMB'][:]
        d_g_Nb_CAMB[i,:] = one_time['delta_g_Nb_CAMB'][:]
        d_ur_Nb_CAMB[i,:] = one_time['delta_ur_Nb_CAMB'][:]
        d_total_Nb_CAMB[i,:] = one_time['delta_total_Nb_CAMB'][:]
        

    # zs_COLA = [0.000, 0.020, 0.041, 0.062, 0.085, 0.109, 0.133, 0.159, 0.186, 0.214, 0.244, 0.275, 0.308, 
    #             0.342, 0.378, 0.417, 0.457, 0.500, 0.543, 0.588, 0.636, 0.688, 0.742, 0.800, 0.862, 0.929, 
    #             1.000, 1.087, 1.182, 1.286, 1.400, 1.526, 1.667, 1.824, 2.000, 2.158, 2.333, 2.529, 2.750, 
    #             3.000, 3.289, 3.624, 4.015, 4.478, 5.036, 5.720, 6.579, 7.690, 9.182, 11.293, 14.508, 20.000, 20.5, 21]
    
    #NEW zs_COLA from z_ini=19. z>3 has different timesteps
    zs_COLA = [0.000, 0.020, 0.041, 0.062, 0.085, 0.109, 0.133, 0.159, 0.186, 0.214, 0.244, 0.275, 0.308, 
                0.342, 0.378, 0.417, 0.457, 0.500, 0.543, 0.588, 0.636, 0.688, 0.742, 0.800, 0.862, 0.929, 
                1.000, 1.087, 1.182, 1.286, 1.400, 1.526, 1.667, 1.824, 2.000, 2.158, 2.333, 2.529, 2.750, 
                3.000, 3.286, 3.615, 4.000, 4.455, 5.000, 5.667, 6.500, 7.571, 9.000, 11.000, 14.000, 19.000, 
                19.500, 20.000]
    
    zs_COLA = np.flip(zs_COLA)

    ks_COLA = one_time['k (h/Mpc)']

    T_d_GR_COLA = np.zeros((len(zs_COLA),len(ks_COLA)),'float64') 
    T_d_m_Nb_COLA = np.zeros((len(zs_COLA),len(ks_COLA)),'float64')
    T_d_cdm_Nb_COLA = np.zeros((len(zs_COLA),len(ks_COLA)),'float64')
    T_d_b_Nb_COLA = np.zeros((len(zs_COLA),len(ks_COLA)),'float64')
    T_d_g_Nb_COLA = np.zeros((len(zs_COLA),len(ks_COLA)),'float64')
    T_d_ur_Nb_COLA = np.zeros((len(zs_COLA),len(ks_COLA)),'float64')
    T_d_total_Nb_COLA = np.zeros((len(zs_COLA),len(ks_COLA)),'float64')

    T_Geff = np.zeros((len(zs_COLA),len(ks_COLA)),'float64') 


    d_GR_CAMB_at_tau = interp1d(tau,d_GR_CAMB,axis=0)
    d_m_Nb_CAMB_at_tau = interp1d(tau,d_m_Nb_CAMB,axis=0)
    d_b_Nb_CAMB_at_tau = interp1d(tau,d_b_Nb_CAMB,axis=0)
    d_cdm_Nb_CAMB_at_tau = interp1d(tau,d_cdm_Nb_CAMB,axis=0)
    d_g_Nb_CAMB_at_tau = interp1d(tau,d_g_Nb_CAMB,axis=0)
    d_ur_Nb_CAMB_at_tau = interp1d(tau,d_ur_Nb_CAMB,axis=0)
    d_total_Nb_CAMB_at_tau = interp1d(tau,d_total_Nb_CAMB,axis=0)
    # d_GR has a minus sign due to CLASS/CAMB conventions!

    for index_z in range(len(zs_COLA)):
        T_d_GR_COLA[index_z,:] = -d_GR_CAMB_at_tau(background_tau_at_z(zs_COLA[index_z]))
        T_d_m_Nb_COLA[index_z,:] = d_m_Nb_CAMB_at_tau(background_tau_at_z(zs_COLA[index_z]))
        T_d_cdm_Nb_COLA[index_z,:] = d_cdm_Nb_CAMB_at_tau(background_tau_at_z(zs_COLA[index_z]))
        T_d_b_Nb_COLA[index_z,:] = d_b_Nb_CAMB_at_tau(background_tau_at_z(zs_COLA[index_z]))
        T_d_g_Nb_COLA[index_z,:] = d_g_Nb_CAMB_at_tau(background_tau_at_z(zs_COLA[index_z]))
        T_d_ur_Nb_COLA[index_z,:] = d_ur_Nb_CAMB_at_tau(background_tau_at_z(zs_COLA[index_z]))
        T_d_total_Nb_COLA[index_z,:] = d_total_Nb_CAMB_at_tau(background_tau_at_z(zs_COLA[index_z]))
        T_Geff[index_z,:] = np.ones(len(ks_COLA))

    #T_d_no_nu is just T_d_m_Nb
    T_d_no_nu_Nb_COLA = T_d_m_Nb_COLA
    T_Weyl_de_COLA = np.zeros((len(zs_COLA),len(ks_COLA)),'float64')
    T_v_cdm_COLA = np.zeros((len(zs_COLA),len(ks_COLA)),'float64')
    T_v_b_COLA = np.zeros((len(zs_COLA),len(ks_COLA)),'float64')
    T_v_b_min_v_c_COLA = np.zeros((len(zs_COLA),len(ks_COLA)),'float64')
    T_d_total_de_COLA = np.zeros((len(zs_COLA),len(ks_COLA)),'float64')
    
    # Formatting data
    data = {}
    for index_z in range(len(zs_COLA)):
        #for index_column in range(12):
        data[zs_COLA[index_z]] = np.vstack((ks_COLA, T_d_no_nu_Nb_COLA[index_z,:], T_d_no_nu_Nb_COLA[index_z,:],
                                            T_d_g_Nb_COLA[index_z,:], T_d_ur_Nb_COLA[index_z,:],
                                            T_d_GR_COLA[index_z,:], T_d_total_Nb_COLA[index_z,:],
                                            T_d_no_nu_Nb_COLA[index_z,:], T_d_total_de_COLA[index_z,:],
                                            T_Weyl_de_COLA[index_z,:], T_v_cdm_COLA[index_z,:],
                                            T_v_b_COLA[index_z,:], T_v_b_min_v_c_COLA[index_z,:])).T

    #Routine to save a bunch of files in the COLA format
    #amypond_direc = f'/home/grads/data/Joao/COLA_projects/wCDM/transfer_functions_lcdm/{filename}/'
    amypond_direc = f'/home/grads/data/jonathan/cola_projects/LCDM_800/transfer_functions/{filename}/'
    try:
        os.mkdir(amypond_direc)
    except FileExistsError:
        pass

    #cluster_path = f'/scratch/decola/joao.reboucas2/COLA_projects/calibration/transfer_functions/{filename}' 
    #transferinfo_path = f'/home/grads/data/Joao/COLA_projects/wCDM/transferinfo_files/transferinfo{filename}.txt'

    cluster_path = f'/gpfs/projects/MirandaGroup/jonathan/cola_projects/LCDM_800/transfer_functions/{filename}/'
    transferinfo_path = f'/home/grads/data/jonathan/cola_projects/LCDM_800/transferinfo_files/transferinfo{filename}.txt'

    for i, z_COLA in enumerate(zs_COLA):
        fname = f"data_transfer_z{z_COLA:.3f}.dat"
        np.savetxt(amypond_direc+fname, data[z_COLA], fmt='   %.7E   %.7E   %.7E   %.7E   %.7E   %.7E   %.7E   %.7E   %.7E   %.7E   %.7E   %.7E   %.7E'
                , header='          k/h    CDM            baryon         photon         nu             mass_nu        total          no_nu          total_de       Weyl           v_CDM          v_b            v_b-v_c     ')

    with open(transferinfo_path, 'w') as transferinfo_file:
        transferinfo_file.write(cluster_path + f' {len(zs_COLA)}\n')
        for i in reversed(range(len(zs_COLA))):
            fname = f"data_transfer_z{zs_COLA[i]:.3f}.dat"
            transferinfo_file.write(fname + '  ' + str(np.round(zs_COLA[i],5))+'\n')
    print(f'Saved transfer files for cosmology {filename} in {amypond_direc}')
    print(f'Saved transferinfo for cosmology {filename} in {transferinfo_path}')
    cosmo.struct_cleanup()

print(f'Generation ended, check {amypond_direc}')
print('----------------------------------------------------------------')
