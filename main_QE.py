import os
import shutil
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
pseudopotentials = {'Sn': PP}
a0 = [6.5, 6.6, 6.7, 6.8, 6.9, 7.0, 7.1]
a0 = [6.7]
KE_cut = [60, 80, 100]
E = []
k = [10, 15, 20]
R = [1]
print(f"r:{R}")
calc = f"scf"
walltime_mins = 30
nodes = 2
procs = 8
d = f"^"  # Here you can set the desired delimiter
equals = f"+"

bash_file = open("run.sh", "w+")
bash_file.write(f"#!/bin/bash\n\n")
bash_file.write(f"dir_list=(")

bat_file = open("rsyn_out.bat", "w+")
bat_file.write(f'wsl rsync -avtuz -e "ssh -p 33301" ./ rathnayake@localhost:~/Run_files')
bat_file.close()

bat_file_in = open("rsyn_in.bat", "w+")
bat_file_in.write(f'wsl rsync -avtuz -e --min-size=5m "ssh -p 33301" rathnayake@localhost:~/Run_files/ ./')
bat_file_in.close()

#Set naming parameters here
dict = {"Calc":"QE",  # (Calculator) QE-Quantum Espresso, SI-siesta
        "Specie":"Sn",
        "XC":"pbe",
        "Struct":"bulk"
        }
run_name2 = ""
for key, val in dict.items():
    run_name2 = f"{run_name2}{key}{equals}{val}{d}"

for KE_cut_i in KE_cut:
    for a0_i in a0:
        for k_i in k:
            R_i = 4*KE_cut_i
            run_name = f"{run_name2}KE{equals}{KE_cut_i}{d}K{equals}{k_i}{d}R{equals}{R_i}{d}a{equals}{a0_i}{d}PP{equals}{PP}{d}type{equals}{calc}"
            bulk = Atoms('Sn2', [(0, 0, 0), (0.25, 0.25, 0.25)],  pbc=True)
            b = a0_i/2.0
            bulk.set_cell([(0, b, b), (b, 0, b), (b, b, 0)], scale_atoms=True)
            ase.io.write(f"run.in", bulk, format = "espresso-in", 
                            label           = f"{run_name}",
                            pseudopotentials= pseudopotentials,
                            pseudo_dir      = "../PseudopotentialDatabase",
                            kpts            = (k_i, k_i, k_i),
                            ecutwfc         = KE_cut_i,
                            calculation     = f"{calc}",
                            lspinorb        = True,
                            noncolin        = True,
                            # ecutrho         = R_i,
                            occupations     = 'smearing',
                            smearing        = 'gaussian',
                            degauss         = 0.01,
                            mixing_beta     = 0.7)

            with open (f"run.job", "w") as file:
                file.write(f"#!/bin/bash\n")
                file.write(f"#\n")
                file.write(f"#  Basics: Number of nodes, processors per node (ppn), and walltime (hhh:mm:ss)\n")
                file.write(f"#PBS -l nodes={nodes}:ppn={procs}\n")
                file.write(f"#PBS -l walltime=1:{walltime_mins}:00\n")
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
                file.write(f"mpirun -machinefile  $PBS_NODEFILE -np $PBS_NP pw.x -in run.in > run.out\n")
            if os.path.exists(run_name):
                shutil.rmtree(run_name)
                print("Path exists!! Overwriting")
            os.mkdir(f"{run_name}")
            os.rename(f"run.in", f"./{run_name}/run.in")
            os.rename(f"run.job", f"./{run_name}/run.job")
            bash_file.write(f" '{run_name}'")
            # bash_file.write(f"dos2uix ./{run_name}/run.job\n")
            # bash_file.write(f"qsub -A cnm66441 ./{run_name}/run.job\n\n")
bash_file.write(")\n")
bash_file.write('for dir in "${dir_list[@]}"\n')
bash_file.write(f"do\n")
bash_file.write(f'  dos2unix "$dir"/run.job | qsub -w "$dir" -A cnm66441\ndone\n')
bash_file.close()