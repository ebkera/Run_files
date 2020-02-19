#!/bin/bash
# Submit jobs from explicitly specified directories;
# stern, 2020-02-18

shopt -s extglob	# handle "+()" patterns

# Let's use a shell loop to read a list of tasks from a "Here-Document" (at the
# end of the loop).  See also: http://www.tldp.org/LDP/abs/html/here-docs.html
while read dir
do
    # Skip empty lines and comments
    [[ $dir == +(""|"#"*) ]] && continue

    # Let's use some basic shell variable string operations to remove the
    # leading and trailing parts of the dir names that are the same, isolating
    # the changing middle section as a useful short name.
    job_name=${dir#*KE+}
    job_name=${job_name%^a*}

    # Use another here-doc to read the job script, to obviate individual files
    # in each data directory.
    #
    # A here-doc with leading "-" will get leading TABs removed.  (unwise to
    # use, however, if your editor munges TABs.)
    #
    # Double-quote $dir to avoid parameter substitution.  (This is not strictly
    # necessary here, though, because the chars "^+" are not special -- in the
    # circumstances used here.)
    qsub -w "$PWD/$dir" -N "$job_name" <<-END_JOB_SCRIPT
	#!/bin/bash
	#PBS -l nodes=2:ppn=8
	#PBS -l walltime=1:30:00
	#PBS -A cnm66441
	#PBS -e $PWD/$dir/job.err
	#PBS -j eo
	#PBS -m bea

	cd \$PBS_O_WORKDIR

	# use a per-job lineup of modules; stern
	module purge
	module load intel
	module load openmpi/1.10/intel-17
	module load quantum-espresso/5.4/openmpi-1.10
	module list

	mpirun pw.x -in run.in -out run.out
END_JOB_SCRIPT


done <<'END_TASKLIST'
    # Single quoting the limit string 'EOT' will pass strings without shell variable and execution expansion.
    # Comments and empty line are fine because we explicitly skip them.

	 Calc+QE^Specie+Sn^XC+pbe^Struct+bulk^KE+60^K+10^R+240^a+6.7^PP+Sn_ONCV_PBE_FR-1.1.upf^type+scf
	 Calc+QE^Specie+Sn^XC+pbe^Struct+bulk^KE+60^K+15^R+240^a+6.7^PP+Sn_ONCV_PBE_FR-1.1.upf^type+scf
	 Calc+QE^Specie+Sn^XC+pbe^Struct+bulk^KE+60^K+20^R+240^a+6.7^PP+Sn_ONCV_PBE_FR-1.1.upf^type+scf
	 Calc+QE^Specie+Sn^XC+pbe^Struct+bulk^KE+80^K+10^R+320^a+6.7^PP+Sn_ONCV_PBE_FR-1.1.upf^type+scf
	 Calc+QE^Specie+Sn^XC+pbe^Struct+bulk^KE+80^K+15^R+320^a+6.7^PP+Sn_ONCV_PBE_FR-1.1.upf^type+scf
	 Calc+QE^Specie+Sn^XC+pbe^Struct+bulk^KE+80^K+20^R+320^a+6.7^PP+Sn_ONCV_PBE_FR-1.1.upf^type+scf
	Calc+QE^Specie+Sn^XC+pbe^Struct+bulk^KE+100^K+10^R+400^a+6.7^PP+Sn_ONCV_PBE_FR-1.1.upf^type+scf
	Calc+QE^Specie+Sn^XC+pbe^Struct+bulk^KE+100^K+15^R+400^a+6.7^PP+Sn_ONCV_PBE_FR-1.1.upf^type+scf
	Calc+QE^Specie+Sn^XC+pbe^Struct+bulk^KE+100^K+20^R+400^a+6.7^PP+Sn_ONCV_PBE_FR-1.1.upf^type+scf
	Calc+QE^Specie+Sn^XC+pbe^Struct+bulk^KE+100^K+20^R+400^a+6.7^PP+Sn_ONCV_PBE_FR-1.1.upf^type+scf

END_TASKLIST
