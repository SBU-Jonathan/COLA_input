import numpy as np
import scipy
import scipy.stats
from scipy.stats import qmc
import math
import sys, platform, os
import matplotlib
import matplotlib.pyplot as plt
import time

#Set these 4 settings:
model = 'LCDM'
num_points = 400
precision = 'default'
#need to make separate files for phase a and b to pair-fix
phase = 'b'

z_ini = 19.0

if precision == 'default':
    force_nmesh = 2048
    box_size = 1024
    prec_suffix = '1'
elif precision == 'high':
    force_nmesh = 3072
    box_size = 512
    prec_suffix = '2'

if phase == 'a':
        rev_phase = 'false'
    elif phase == 'b':
        rev_phase = 'true'

for n_iter in range(num_points):
    
    run_name = f'{n_iter}'
    
    #-----------------------------------------------------------
    #These 4 paths require attention:
    
    #the "_1" is hard-coded into params and transferinfo paths since these files are same for both precisions, so we only store them once in default-precision directory and just point 
    
    #the path on the current computer to the cosmo parameters which were output by CLASS when computing transfer functions
    cosmo_params_path = f'/home/grads/data/jonathan/cola_projects/COLA_input/{model}/{num_points}_1/cosmo_params/params{n_iter}.txt'
    
    #where on the current computer to store the lua files this script writes:
    lua_out_path = f'/home/grads/data/jonathan/cola_projects/COLA_input/{model}/{num_points}_{prec_suffix}/lua_files_{phase}/parameter_file{run_name}.lua'
    
    #where the transferinfo files will be read from on the cluster:
    #if git repo is cloned onto cluster then make it <path_to_repo>/COLA_input/{model}/{num_points}_1/transferinfo_files/transferinfor{run_name}.txt
    transferinfo_path = f'/gpfs/projects/MirandaGroup/jonathan/cola_projects/{model}/{num_points}_1/transferinfo_files/transferinfo{run_name}.txt'
    
    #where the output pk will be stored on the cluster:
    output_path = f'/gpfs/projects/MirandaGroup/jonathan/cola_projects/{model}/{num_points}_{prec_suffix}/{phase}/output/{run_name}'
    #------------------------------------------------------------
    
    params = np.loadtxt(cosmo_params_path)
    
    Om_ = params[3]
    ns_ = params[5]
    h_ = params[6]
    As_ = params[4]
    Ob_ = params[2]
    Omnu_ = params[0] 
    Oc_ = params[1]
    
    file = open(lua_out_path, 'w')
    
    file.write('all_parameters_must_be_in_file = true\n\n')
    file.write('particle_Npart_1D = 1024\n')
    file.write(f'force_nmesh = {force_nmesh}\n')
    file.write('ic_random_seed = 1234\n')
    file.write('output_redshifts = {3, 2, 1, 0.5, 0.0}\n')
    file.write('timestep_nsteps = {12, 5, 8, 9, 17}\n')
    file.write('output_particles = false\n')
    file.write('fof = false\n')
    file.write('fof_nmin_per_halo = 20\n')
    file.write('pofk = true\n')
    file.write('pofk_nmesh = 1024\n')
    file.write('ic_fix_amplitude = true\n')
    file.write(f'ic_reverse_phases = {rev_phase}\n')
    file.write('ic_type_of_input = "transferinfofile"\n')
    file.write('ic_input_filename = "' + transferinfo_path + '"\n\n')
    file.write('output_folder = "' + output_path + '"\n')
    file.write('simulation_name = "run_' + run_name + '"\n')
    file.write(f'simulation_boxsize = {box_size}\n\n')
    file.write('simulation_use_cola = true\n')
    file.write('simulation_use_scaledependent_cola = true\n\n')
    file.write('cosmology_model = "LCDM"\n')
    file.write('cosmology_OmegaCDM = ' + str(Oc_) + '\n')
    file.write('cosmology_Omegab = ' + str(Ob_) + '\n')
    file.write('cosmology_OmegaMNu = ' + str(Omnu_) + '\n')
    file.write('cosmology_OmegaLambda = 1 - cosmology_OmegaCDM - cosmology_Omegab - cosmology_OmegaMNu\n')
    file.write('cosmology_Neffective = 3.046\n')
    file.write('cosmology_TCMB_kelvin = 2.7255\n')
    file.write('cosmology_h = ' + str(h_) + '\n')
    file.write('cosmology_As = ' + str(As_) + '\n')
    file.write('cosmology_ns = ' + str(ns_) + '\n')
    file.write('cosmology_kpivot_mpc = 0.05\n')
    file.write('if cosmology_model == "w0waCDM" then \n')
    file.write('  cosmology_w0 = -1.0\n')
    file.write('  cosmology_wa = 0.0\n')
    file.write('end\n\n')
    file.write('if cosmology_model == "DGP" then \n')
    file.write('  cosmology_dgp_OmegaRC = 0.11642\n')
    file.write('end\n\n')
    file.write('gravity_model = "GR"\n\n')
    file.write('if gravity_model == "f(R)" then \n')
    file.write('  gravity_model_fofr_fofr0 = 1e-5\n')
    file.write('  gravity_model_fofr_nfofr = 1.0\n')
    file.write('  gravity_model_screening = true\n')
    file.write('  gravity_model_screening_enforce_largescale_linear = true\n')
    file.write('  gravity_model_screening_linear_scale_hmpc = 0.1\n')
    file.write('end\n\n')
    file.write('if gravity_model == "DGP" then \n')
    file.write('  gravity_model_dgp_rcH0overc = 1.0\n')
    file.write('  gravity_model_screening = true\n')
    file.write('  gravity_model_dgp_smoothing_filter = "tophat"\n')
    file.write('  gravity_model_dgp_smoothing_scale_over_boxsize = 0.0\n')
    file.write('  gravity_model_screening_enforce_largescale_linear = true\n')
    file.write('  gravity_model_screening_linear_scale_hmpc = 0.1\n')
    file.write('end\n\n')
    file.write('particle_allocation_factor = 1.25\n\n')
    file.write('output_fileformat = "GADGET"\n\n')
    file.write('timestep_method = "Quinn"\n')
    file.write('timestep_cola_nLPT = -2.5\n')
    file.write('timestep_algorithm = "KDK"\n')
    file.write('timestep_scalefactor_spacing = "linear"\n')
    file.write('if timestep_scalefactor_spacing == "powerlaw" then\n')
    file.write('  timestep_spacing_power = 1.0\n')
    file.write('end\n\n')
    file.write('ic_random_generator = "GSL"\n')
    file.write('ic_random_field_type = "gaussian"\n')
    file.write('ic_nmesh = particle_Npart_1D\n')
    file.write('ic_use_gravity_model_GR = true\n')
    file.write('ic_LPT_order = 2\n')
    file.write(f'ic_input_redshift = {z_ini}\n')
    file.write(f'ic_initial_redshift = {z_ini}\n')
    file.write('ic_sigma8_normalization = false\n')
    file.write('ic_sigma8_redshift = 0.0\n')
    file.write('ic_sigma8 = 0.83\n\n')
    file.write('if ic_random_field_type == "nongaussian" then\n')
    file.write('  ic_fnl_type = "local"\n')
    file.write('  ic_fnl = 100.0\n')
    file.write('  ic_fnl_redshift = ic_initial_redshift\n')
    file.write('end\n\n')
    file.write('if ic_random_field_type == "reconstruct_from_particles" then\n')
    file.write('  ic_reconstruct_gadgetfilepath = "output/gadget"\n')
    file.write('  ic_reconstruct_assigment_method = "CIC"\n')
    file.write('  ic_reconstruct_smoothing_filter = "sharpk"\n')
    file.write('  ic_reconstruct_dimless_smoothing_scale = 1.0 / (128 * math.pi)\n')
    file.write('  ic_reconstruct_interlacing = false\n')
    file.write('end\n\n')
    file.write('if ic_random_field_type == "read_particles" then\n')
    file.write('  ic_reconstruct_gadgetfilepath = "path/gadget"\n')
    file.write('end\n\n')
    file.write('force_density_assignment_method = "CIC"\n')
    file.write('force_kernel = "continuous_greens_function"\n')
    file.write('force_linear_massive_neutrinos = true\n\n')
    file.write('fof_linking_length = 0.2 / particle_Npart_1D\n')
    file.write('fof_nmesh_max = 0\n\n')
    file.write('pofk_interlacing = true\n')
    file.write('pofk_subtract_shotnoise = false\n')
    file.write('pofk_density_assignment_method = "PCS"\n\n')
    file.write('pofk_multipole = false\n')
    file.write('pofk_multipole_nmesh = 128\n')
    file.write('pofk_multipole_interlacing = true\n')
    file.write('pofk_multipole_subtract_shotnoise = false\n')
    file.write('pofk_multipole_ellmax = 4\n')
    file.write('pofk_multipole_density_assignment_method = "PCS"\n\n')
    file.write('bispectrum = false\n')
    file.write('bispectrum_nmesh = 128\n')
    file.write('bispectrum_nbins = 10\n')
    file.write('bispectrum_interlacing = true\n')
    file.write('bispectrum_subtract_shotnoise = false\n')
    file.write('bispectrum_density_assignment_method = "PCS"\n')

    file.close()