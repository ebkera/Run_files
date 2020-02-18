#!/bin/bash

dir_list=("Calc=QE^Specie=Sn^XC=pbe^Struct=bulk^KE=60^K=10^R=240^a=6.7^PP=Sn_ONCV_PBE_FR-1.1.upf^type=scf" "Calc=QE^Specie=Sn^XC=pbe^Struct=bulk^KE=60^K=15^R=240^a=6.7^PP=Sn_ONCV_PBE_FR-1.1.upf^type=scf" "Calc=QE^Specie=Sn^XC=pbe^Struct=bulk^KE=60^K=20^R=240^a=6.7^PP=Sn_ONCV_PBE_FR-1.1.upf^type=scf" "Calc=QE^Specie=Sn^XC=pbe^Struct=bulk^KE=80^K=10^R=320^a=6.7^PP=Sn_ONCV_PBE_FR-1.1.upf^type=scf" "Calc=QE^Specie=Sn^XC=pbe^Struct=bulk^KE=80^K=15^R=320^a=6.7^PP=Sn_ONCV_PBE_FR-1.1.upf^type=scf" "Calc=QE^Specie=Sn^XC=pbe^Struct=bulk^KE=80^K=20^R=320^a=6.7^PP=Sn_ONCV_PBE_FR-1.1.upf^type=scf" "Calc=QE^Specie=Sn^XC=pbe^Struct=bulk^KE=100^K=10^R=400^a=6.7^PP=Sn_ONCV_PBE_FR-1.1.upf^type=scf" "Calc=QE^Specie=Sn^XC=pbe^Struct=bulk^KE=100^K=15^R=400^a=6.7^PP=Sn_ONCV_PBE_FR-1.1.upf^type=scf" "Calc=QE^Specie=Sn^XC=pbe^Struct=bulk^KE=100^K=20^R=400^a=6.7^PP=Sn_ONCV_PBE_FR-1.1.upf^type=scf" )
for dir in $dir_list
do
  dos2unix "$dir"/run.job | qsub -w "$dir" -A cnm66441
done
