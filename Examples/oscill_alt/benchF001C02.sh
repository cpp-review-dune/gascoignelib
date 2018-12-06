../../bin/Examples/oscill_alt/para_box3d 1 0.52 0.02 0.501 0.001 1> /dev/null 2>> timelog.txt
mkdir Plots_and_data/box3d/C02tc52I1
cp Results/p.*.bup Plots_and_data/box3d/C02tc52I1
cp Results/p.*.vtk Plots_and_data/box3d/C02tc52I1
cp para_functionalC02tc52I1.txt Plots_and_data/box3d/C02tc52I1

sleep 30
../../bin/Examples/oscill_alt/para_box3d 2 0.52 0.02 0.501 0.001 1> /dev/null 2>> timelog.txt
mkdir Plots_and_data/box3d/C02tc52I2
cp Results/p.*.bup Plots_and_data/box3d/C02tc52I2
cp Results/p.*.vtk Plots_and_data/box3d/C02tc52I2
cp para_functionalC02tc52I2.txt Plots_and_data/box3d/C02tc52I2

sleep 30
../../bin/Examples/oscill_alt/para_box3d 3 0.52 0.02 0.501 0.001 1> /dev/null 2>> timelog.txt
mkdir Plots_and_data/box3d/C02tc52I3
cp Results/p.*.bup Plots_and_data/box3d/C02tc52I3
cp Results/p.*.vtk Plots_and_data/box3d/C02tc52I3
cp para_functionalC02tc52I3.txt Plots_and_data/box3d/C02tc52I3
