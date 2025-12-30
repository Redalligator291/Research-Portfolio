echo 0 | gmx_mpi trjconv -f md_2.xtc -s md_2.tpr -n index.ndx -o md_solu.xtc -pbc whole
echo 0 | gmx_mpi trjconv -f npt.gro -s npt.tpr -n index.ndx -o md_solu.pdb -pbc whole
#echo 2 | gmx_mpi trjconv -f md_2.xtc -s md_2.tpr -n index.ndx -o md.gro -dump 1500000 
echo 0 | gmx_mpi trjconv -f md_2.gro -s md_2.tpr -n index.ndx -o md_final.pdb -pbc whole
