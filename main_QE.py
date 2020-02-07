import os
from ase.build import bulk
from ase.constraints import UnitCellFilter
from ase.optimize import LBFGS
from ase import Atoms
from ase.calculators.siesta import Siesta
from ase.units import Ry
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('Agg')  # no UI backend required if working in the wsl without a UIs
# from ebk import ASEInterface
import ase.io
from ebk.QE.QErunfilecreator import QERunCreator
from ebk.runscriptcreator import Runscriptcreator

PP = "Sn_ONCV_PBE_FR-1.1.upf"
pseudopotentials = {'Sn': 'Sn_ONCV_PBE_FR-1.1.upf'}
a0 = [6.6, 6.7, 6.8, 6.9]
KE_cut = [20, 40, 60, 80, 100]
E = []
KE_cut = [20, 40]
k = [3]
R = [300]
calc = f"scf"
walltime_mins = 30
nodes = 2
procs = 8
d = f"d"  # Here you can set the desired delimiter
val = f"f"

bash_file = open("run.sh", "w+")
bash_file.write(f"#!/bin/bash\n\n")

# bat_file = open("rsyn_out.bat", "w+")
# bat_file.write(f'wsl rsync -avtuz -e "ssh -p 33301" ./ rathnayake@localhost:~/Run_files')
# bat_file.close()

# bat_file_in = open("rsyn_in.bat", "w+")
# bat_file_in.write(f'wsl rsync -avtuz -e "ssh -p 33301" rathnayake@localhost:~/Run_files/ ./')
# bat_file_in.close()

for KE_cut_i in KE_cut:
    for a0_i in a0:
        for k_i in k:
            for R_i in R:
                run_name = f"QE{d}KE{val}{KE_cut_i}{d}K{val}{k_i}{d}R{val}{R_i}{d}a{val}{a0_i}{d}PP{val}{PP}{d}calc{val}{calc}"
                bulk = Atoms('Sn2', [(0, 0, 0), (0.25, 0.25, 0.25)],  pbc=True)
                b = a0_i/2.0
                bulk.set_cell([(0, b, b), (b, 0, b), (b, b, 0)], scale_atoms=True)
                ase.io.write(f"{run_name}.in", bulk, format = "espresso-in", 
                                label           = f"{run_name}",
                                pseudopotentials= pseudopotentials,
                                pseudo_dir      = "../PseudopotentialDatabase",
                                kpts            = (k_i, k_i, k_i),
                                ecutwfc         = KE_cut_i,
                                calculation     = f"{calc}",
                                lspinorb        = True,
                                noncolin        = True,
                                ecutrho         = R_i,
                                occupations     = 'smearing',
                                smearing        = 'gaussian',
                                degauss         = 0.01,
                                mixing_beta     = 0.7)

                with open (f"{run_name}.job", "w") as file:
                    file.write(f"#!/bin/bash\n")
                    file.write(f"#\n")
                    file.write(f"#  Basics: Number of nodes, processors per node (ppn), and walltime (hhh:mm:ss)\n")
                    file.write(f"#PBS -l nodes={nodes}:ppn={procs}\n")
                    file.write(f"#PBS -l walltime=0:{walltime_mins}:00\n")
                    file.write(f"#PBS -N {run_name}\n")
                    file.write(f"#PBS -A cnm66441\n")
                    file.write(f"#\n")
                    file.write(f"#  File names for stdout and stderr.  If not set here, the defaults\n")
                    file.write(f"#  are <JOBNAME>.o<JOBNUM> and <JOBNAME>.e<JOBNUM>\n")
                    file.write(f"#PBS -o job.out\n")
                    file.write(f"#PBS -e job.err\n")
                    file.write("\n")
                    file.write(f"# Send mail at begin, end, abort, or never (b, e, a, n). Default is 'a'.\n")
                    file.write(f"#PBS -m bea erathnayake@sivananthanlabs.us\n")
                    file.write("\n")
                    file.write(f"# change into the directory where qsub will be executed\n")
                    file.write(f"cd $PBS_O_WORKDIR\n")
                    file.write("\n")
                    file.write(f"# start MPI job over default interconnect; count allocated cores on the fly.\n")
                    file.write(f"mpirun -machinefile  $PBS_NODEFILE -np $PBS_NP pw.x -in {run_name}.in > {run_name}.out\n")
                os.mkdir(f"{run_name}")
                os.rename(f"{run_name}.in", f"./{run_name}/{run_name}.in")
                os.rename(f"{run_name}.job", f"./{run_name}/{run_name}.job")
                bash_file.write(f"qsub -A cnm66441 ./{run_name}/{run_name}.job\n\n")
bash_file.close()