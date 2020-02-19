#!/bin/bash
# Submit jobs from explicitly specified directories;
# stern, 2020-02-18

shopt -s extglob	# handle "+()" patterns

# Let's use a "Here-Document", http://www.tldp.org/LDP/abs/html/here-docs.html
while read dir
do
  # skip empty lines and comments
  [[ $dir == +(""|"#"*) ]] && continue

  # Run in a sub-shell "( ... )" to avoid "cd .." getting lost in hierarchy,
  # which would happen if a $dir does not exist.
  #
  # Double-quote $dir to avoid parameter substitution (before being passed to
  # "cd") due to the unusual file names.  (This is, however, not strictly
  # necessary here because the chars "^+" are not special -- in the
  # circumstances used here.)
  #
  ( cd "$dir"; qsub *job )

done <<'EOT'
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

EOT
