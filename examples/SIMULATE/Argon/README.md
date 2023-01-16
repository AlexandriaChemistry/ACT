Simple test of the alexandria simulate utility
----------------------------------------------
First minimize the input coordinates using
```
alexandria simulate -minimize -f tetramer.xyz -ff noblegas.xml -c minimized.pdb -maxiter 10000 -minalg Steep -toler 1e-4 -g em.log
```
This will give you an output struture minimized.pdb.
Then you can run a short simulation with starting temperature 50 K.
```
alexandria simulate -f minimized.pdb -ff noblegas.xml -temp 50 -nsteps 10000 -nstxout 100 -g md.log
```
Check the energies using the ACT script viewxvg
```
viewxvg -f energy.xvg
```
and please check the trajectory using your favorite pdb viewer.

Please do also inspect the log files.
