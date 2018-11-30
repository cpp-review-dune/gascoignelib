sleep 30
../../bin/paraoscill 1 0.55 0.05 0.505 0.005 1> /dev/null 2>> timelog.txt
mkdir Plots_and_data/wall_mount/F005C05tc55I1
cp Results/p.?.* Plots_and_data/wall_mount/F005C05tc55I1
cp Results/p.0000{0,1,2,3,4,5,6,7}.vtk Plots_and_data/wall_mount/F005C05tc55I1
cp para_functionalC05tc55I1.txt Plots_and_data/wall_mount/F005C05tc55I1

sleep 30
../../bin/paraoscill 2 0.55 0.05 0.505 0.005 1> /dev/null 2>> timelog.txt
mkdir Plots_and_data/wall_mount/F005C05tc55I2
cp Results/p.?.* Plots_and_data/wall_mount/F005C05tc55I2
cp Results/p.0000{0,1,2,3,4,5,6,7}.vtk Plots_and_data/wall_mount/F005C05tc55I2
cp para_functionalC05tc55I2.txt Plots_and_data/wall_mount/F005C05tc55I2

sleep 30
../../bin/paraoscill 3 0.55 0.05 0.505 0.005 1> /dev/null 2>> timelog.txt
mkdir Plots_and_data/wall_mount/F005C05tc55I3
cp Results/p.?.* Plots_and_data/wall_mount/F005C05tc55I3
cp Results/p.0000{0,1,2,3,4,5,6,7}.vtk Plots_and_data/wall_mount/F005C05tc55I3
cp para_functionalC05tc55I3.txt Plots_and_data/wall_mount/F005C05tc55I3
