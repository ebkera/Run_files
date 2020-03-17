from ebk.runscripthandler import RunScriptHandler, ReadOutfiles, make_all_job_files
from ase.atoms import Atoms
from ase.io import read
import ase.io
from ebk import get_pseudopotential

# Here we have details of the pseudos we will use in this file
# All pseudopotentials are from the QE website
# ___________________________________________________
# QE_PBE_FR_NLCC_1: Sn.rel-pbe-dn-kjpaw_psl.1.0.0.UPF

# Origin: PS Library
# Author: ADC 
# Generated using "atomic" code by A. Dal Corso  v.6.3
# Pseudopotential type: PAW 
# Functional type: PBE
# Non Linear Core Correction
# Full relativistic
# ___________________________________________________
# QE_PZ_FR_NLCC_1: Sn.rel-pz-dn-kjpaw_psl.0.2.UPF

# author="ADC"
# date="11Sep2012"
# pseudo_type="PAW"
# relativistic="full"
# is_ultrasoft="T"
# is_paw="T"
# is_coulomb="F"
# has_so="T"
# core_correction="T"
# functional=" SLA  PZ   NOGX NOGC"
# total_psenergy="-4.318985078769828E+002"
# wfc_cutoff="2.554034232256697E+001"
# rho_cutoff="1.021613692902679E+002"
# ___________________________________________________
# QE_PBE_FR_NLCC_1: Sn.Sn.rel-pbe-dn-rrkjus_psl.1.0.0.UPF

# dft='PBE'
# lpaw=.false.,
# use_xsd=.FALSE.,
# pseudotype=3,
# file_pseudopw='Sn.rel-pbe-dn-rrkjus_psl.1.0.0.UPF',
# author='ADC',
# lloc=-1,
# rcloc=1.9,
# which_augfun='PSQ',
# rmatch_augfun_nc=.true.,
# nlcc=.true.,
# new_core_ps=.true.,
# rcore=1.2,
# tm=.true.

###############################################################################################################################################################
# This block creates bulk runs for convergence
IDs = ["QE_PBE_FR_NLCC_1", "QE_PBE_FR_NLCC_2", "QE_PBE_FR_NLCC_3", "QE_PBESOL_FR_NLCC_1", "QE_PBE_SR_NLCC_1", "QE_PBE_SR_NLCC_2", "QE_PBE_SR_NLCC_3", "QE_PBESOL_SR_NLCC_1"]
IDs = IDs[1:]  # since we have already done the first one we can splice the list to a smaller one.
print(IDs)
for identifier in IDs:
        pseudopotential = get_pseudopotential(identifier)
        # para = {"path":"GXWLGKL"}
        # para.update({"density": 15})
        inputs = {"calculation"   : "scf",
                "lspinorb"        : True,
                "noncolin"        : True,
                "occupations"     : 'smearing',
                "smearing"        : 'gaussian',
                "degauss"         : 0.0001,
                "mixing_beta"     : 0.7,
                "Title"           : 'Sn',
                "prefix"          : 'Sn',
                "restart_mode"    : 'from_scratch',
                "disk_io"         : 'default',
                "verbosity"       : 'high',
                "lkpoint_dir"     : False,
                "etot_conv_thr"   : 1.0e-6,
                "forc_conv_thr"   : 1.0e-4,
                "outdir"          : './',
                "structure_type"  : "bulk",
                # "partition"       : "bigmem"
                }

        Sn = RunScriptHandler(identifier = f"{identifier}", KE_cut = [75, 80, 85, 90], k = [16], a0 = [6.67], **inputs)
        Sn.job_handler = "torque"
        Sn.set_pseudopotentials(pseudopotential)
        Sn.set_pseudo_dir("carbon")

        bulk = Atoms('Sn2', [(0, 0, 0), (0.25, 0.25, 0.25)],  pbc=True)
        Sn.structure = 1

        # Printing out some run specific values
        print(f"Number of runs: {Sn.get_number_of_calculations()}")
        print(Sn.pseudopotentials)

        # times
        Sn.walltime_hours = 1
        Sn.walltime_mins = 30
        Sn.make_runs()
        Sn.create_torque_job()

make_all_job_files(IDs)






























###############################################################################################################################################################

# for identifier in IDs:
#     if identifier == "ONCV_PBE_FR1.1":
#         pseudopotential = {'Sn': f"Sn_ONCV_PBE_FR-1.1.upf"}
#         a_0 = []
#     elif identifier == "ONCV_PBE_SR1.0":
#         pseudopotential = {'Sn': f"Sn_ONCV_PBE-1.0.upf"}
#     elif identifier == "ONCV_PBE_SR1.1":
#         pseudopotential = {'Sn': f"Sn_ONCV_PBE-1.1.upf"}
#     elif identifier == "ONCV_PBE_SR1.2":
#         pseudopotential = {'Sn': f"Sn_ONCV_PBE-1.2.upf"}

# # set lattice constants here
#     if identifier == "ONCV_PBE_FR1.1":
#         a_0 = []
#     elif identifier == "ONCV_PBE_SR1.0":
#         pseudopotential = {'Sn': f"Sn_ONCV_PBE-1.0.upf"}
#     elif identifier == "ONCV_PBE_SR1.1":
#         pseudopotential = {'Sn': f"Sn_ONCV_PBE-1.1.upf"}
#     elif identifier == "ONCV_PBE_SR1.2":
#         pseudopotential = {'Sn': f"Sn_ONCV_PBE-1.2.upf"}



#     inputs = {"calculation"   : "scf",
#             "lspinorb"        : False,
#             "noncolin"        : False,
#             "occupations"     : 'smearing',
#             "smearing"        : 'gaussian',
#             "degauss"         : 0.0001,
#             "mixing_beta"     : 0.7,
#             "Title"           : 'Sn',
#             "prefix"          : 'Sn',
#             "restart_mode"    : 'from_scratch',
#             "disk_io"         : 'default',
#             "verbosity"       : 'high',
#             "lkpoint_dir"     : False,
#             "etot_conv_thr"   : 1.0e-6,
#             "forc_conv_thr"   : 1.0e-4,
#             "outdir"          : './',
#             "structure_type"  : "bulk",
#             # "partition"       : "bigmem"
#             }
#     Sn = RunScriptHandler(identifier = f"{identifier}", KE_cut = [100], k = [30], a0 = [6.4, 6.52, 6.54, 6.56, 6.58, 6.6, 6.62, 6.66, 6.68, 6.70], **inputs)
#     Sn.job_handler = "torque"
#     Sn.set_pseudopotentials(pseudopotential)
#     print(Sn.pseudopotentials)
#     Sn.set_pseudo_dir("carbon")

#     bulk = Atoms('Sn2', [(0, 0, 0), (0.25, 0.25, 0.25)],  pbc=True)
#     Sn.structure = 1

#     print(f"Number of runs: {Sn.get_number_of_calculations()}")

#     # times
#     Sn.walltime_hours = 1
#     Sn.walltime_mins = 30

#     Sn.make_runs()
#     Sn.create_torque_job()




# copy from here
# ligand.nodes = 1
# ligand.procs = ligand.nodes*48
# ligand.npool = 1
# ligand.ntasks = ligand.procs

# pseudopotentials = {'C': 'c_pbe_v1.2.uspp.F.UPF',
#                     'H': 'h_pbe_v1.4.uspp.F.UPF',
#                     'S': 's_pbe_v1.4.uspp.F.UPF'}
# ligand.set_pseudopotentials(pseudopotentials)
# print(ligand.pseudopotentials)

# ethaneDithiol = read("../../Run_files/XYZdatabase/1,2-ethaneDithiol.xyz", index=None, format="xyz")
# a = 65
# ethaneDithiol.set_cell([(a, 0, 0), (0, a, 0), (0, 0, a)], scale_atoms=False)
# ligand.structure = 0
# ligand.atoms_object = ethaneDithiol

# print(f"Number of runs: {ligand.get_number_of_calculations()}")

# times
# ligand.walltime_days = 2
# ligand.walltime_hours = 1
# ligand.walltime_mins = 30

# ligand.make_runs()
# ligand.create_bash_file()
# ligand.create_torque_job()

# reder = Read_outfiles()
# reder.read_folder_names()



###############################################################################################################################################################
# identifier = IDs[3]

# if identifier == "ONCV_PBE_FR1.1":
#         pseudopotential = {'Sn': f"Sn_ONCV_PBE_FR-1.1.upf"}
# elif identifier == "ONCV_PBE_SR1.0":
#         pseudopotential = {'Sn': f"Sn_ONCV_PBE-1.0.upf"}
# elif identifier == "ONCV_PBE_SR1.1":
#         pseudopotential = {'Sn': f"Sn_ONCV_PBE-1.1.upf"}
# elif identifier == "ONCV_PBE_SR1.2":
#         pseudopotential = {'Sn': f"Sn_ONCV_PBE-1.2.upf"}

# inputs = {"calculation"   : "bands",
#         "lspinorb"        : False,
#         "noncolin"        : False,
#         "occupations"     : 'smearing',
#         "smearing"        : 'gaussian',
#         "degauss"         : 0.0001,
#         "mixing_beta"     : 0.7,
#         "Title"           : 'Sn',
#         "prefix"          : 'Sn',
#         "restart_mode"    : 'from_scratch',
#         "disk_io"         : 'default',
#         "verbosity"       : 'high',
#         "lkpoint_dir"     : False,
#         "etot_conv_thr"   : 1.0e-6,
#         "forc_conv_thr"   : 1.0e-4,
#         "outdir"          : './',
#         "structure_type"  : "bulk",
#         # "k_path"        : {"path":"GXWLGKL", "density": 15})
#         # "partition"       : "bigmem"
#         }

# Sn = RunScriptHandler(identifier = f"{identifier}", KE_cut = [80], k = [15], a0 = [6.652], **inputs)
# Sn.job_handler = "torque"
# Sn.set_pseudopotentials(pseudopotential)
# Sn.set_pseudo_dir("carbon")

# bulk = Atoms('Sn2', [(0, 0, 0), (0.25, 0.25, 0.25)],  pbc=True)
# Sn.structure = 1

# # Printing out some run specific values
# print(f"Number of runs: {Sn.get_number_of_calculations()}")
# print(Sn.pseudopotentials)

# # times
# Sn.walltime_hours = 1
# Sn.walltime_mins = 30

# Sn.make_runs()
# Sn.create_torque_job()
