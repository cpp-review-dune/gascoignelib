../../bin/Examples/oscill_alt/para_box3d 1 0.51 0.01 0.501 0.001 1> /dev/null 2>> timelog.txt
mkdir Plots_and_data/box3d/C01tc51I1
cp Results/p.*.bup Plots_and_data/box3d/C01tc51I1
cp Results/p.*.vtk Plots_and_data/box3d/C01tc51I1
cp para_functionalC01tc51I1.txt Plots_and_data/box3d/C01tc51I1

sleep 30
../../bin/Examples/oscill_alt/para_box3d 2 0.51 0.01 0.501 0.001 1> /dev/null 2>> timelog.txt
mkdir Plots_and_data/box3d/C01tc51I2
cp Results/p.*.bup Plots_and_data/box3d/C01tc51I2
cp Results/p.*.vtk Plots_and_data/box3d/C01tc51I2
cp para_functionalC51tc51I2.txt Plots_and_data/box3d/C01tc51I2

sleep 30
../../bin/Examples/oscill_alt/para_box3d 3 0.51 0.01 0.501 0.001 1> /dev/null 2>> timelog.txt
mkdir Plots_and_data/box3d/C01tc51I3
cp Results/p.*.bup Plots_and_data/box3d/C01tc51I3
cp Results/p.*.vtk Plots_and_data/box3d/C01tc51I3
cp para_functionalC51tc51I3.txt Plots_and_data/box3d/C01tc51I3
