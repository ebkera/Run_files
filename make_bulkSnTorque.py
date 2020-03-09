from ebk.runscripthandler import RunScriptHandler, ReadOutfiles
from ase.atoms import Atoms
from ase.io import read
import ase.io
from ebk.kPathCreator import *
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

###############################################################################################################################################################
# This block creates bulk runs for convergence
IDs = ["QE_PBE_FR_NLCC_1"]
identifier = IDs[0]
pseudopotential = get_pseudopotential(identifier)
# para = {"path":"GXWLGKL"}
# para.update({"density": 15})

inputs = {"calculation"   : "bands",
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

Sn = RunScriptHandler(identifier = f"{identifier}", KE_cut = [76], k = [16], a0 = [6.67], **inputs)
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
