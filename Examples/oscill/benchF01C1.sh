sleep 30
../../bin/paraoscill 1 0.6 0.1 0.51 0.01 1> /dev/null 2>> timelog.txt
mkdir Plots_and_data/wall_mount/F01C1tc6I1
cp Results/p.?.* Plots_and_data/wall_mount/F01C1tc6I1
cp Results/p.0000{0,1,2,3,4,5,6,7}.vtk Plots_and_data/wall_mount/F01C1tc6I1
cp para_functionalC1tc6I1.txt Plots_and_data/wall_mount/F01C1tc6I1

sleep 30
../../bin/paraoscill 2 0.6 0.1 0.51 0.01 1> /dev/null 2>> timelog.txt
mkdir Plots_and_data/wall_mount/F01C1tc6I2
cp Results/p.?.* Plots_and_data/wall_mount/F01C1tc6I2
cp Results/p.0000{0,1,2,3,4,5,6,7}.vtk Plots_and_data/wall_mount/F01C1tc6I2
cp para_functionalC1tc6I2.txt Plots_and_data/wall_mount/F01C1tc6I2

sleep 30
../../bin/paraoscill 3 0.6 0.1 0.51 0.01 1> /dev/null 2>> timelog.txt
mkdir Plots_and_data/wall_mount/F01C1tc6I3
cp Results/p.?.* Plots_and_data/wall_mount/F01C1tc6I3
cp Results/p.0000{0,1,2,3,4,5,6,7}.vtk Plots_and_data/wall_mount/F01C1tc6I3
cp para_functionalC1tc6I3.txt Plots_and_data/wall_mount/F01C1tc6I3
