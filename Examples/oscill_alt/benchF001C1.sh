../../bin/Examples/oscill_alt/para_box3d 1 0.6 0.1 0.501 0.001 1> /dev/null 2>> timelog.txt
mkdir Plots_and_data/box3d/C1tc6I1
cp Results/p.*.bup Plots_and_data/box3d/C1tc6I1
cp Results/p.*.vtk Plots_and_data/box3d/C1tc6I1
cp para_functionalC1tc6I1.txt Plots_and_data/box3d/C1tc6I1

sleep 30
../../bin/Examples/oscill_alt/para_box3d 2 0.6 0.1 0.501 0.001 1> /dev/null 2>> timelog.txt
mkdir Plots_and_data/box3d/C1tc6I2
cp Results/p.*.bup Plots_and_data/box3d/C1tc6I2
cp Results/p.*.vtk Plots_and_data/box3d/C1tc6I2
cp para_functionalC1tc6I2.txt Plots_and_data/box3d/C1tc6I2

sleep 30
../../bin/Examples/oscill_alt/para_box3d 3 0.6 0.1 0.501 0.001 1> /dev/null 2>> timelog.txt
mkdir Plots_and_data/box3d/C1tc6I3
cp Results/p.*.bup Plots_and_data/box3d/C1tc6I3
cp Results/p.*.vtk Plots_and_data/box3d/C1tc6I3
cp para_functionalC1tc6I3.txt Plots_and_data/box3d/C1tc6I3
