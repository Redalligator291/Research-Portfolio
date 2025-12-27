#!/bin/bash
export OMPI_ALLOW_RUN_AS_ROOT=1
export OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1


echo 8 | gmx pdb2gmx -f  -water tip3p
gmx editconf -f conf.gro -o newbox.gro -c -d 1.0 -bt cubic
gmx solvate -cp newbox.gro -cs spc216.gro -o solv.gro -p topol.top
gmx grompp -f ions.mdp -c solv.gro -p topol.top -o ions.tpr

echo 13 | gmx genion -s ions.tpr -o solv_ions.gro -p topol.top -conc 0.02 -pname NA -nname CL -neutral

echo q! | gmx make_ndx -f solv_ions.gro -o index.ndx