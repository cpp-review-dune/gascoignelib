#!/bin/zsh
../../bin/paraoscill 1 0.52 0.02 0.502 0.002 1> /dev/null 2>> timelog.txt
mkdir functionals_etc/wall_mount/F002C02tc52I1
cp Results/p.?.* functionals_etc/wall_mount/F002C02tc52I1
cp Results/p.0000{0,1,2,3,4,5,6,7}.vtk functionals_etc/wall_mount/F002C02tc52I1
cp para_functionalC02tc52I1.txt functionals_etc/wall_mount/F002C02tc52I1

sleep 30
../../bin/paraoscill 1 0.525 0.025 0.502 0.002 1> /dev/null 2>> timelog.txt
mkdir functionals_etc/wall_mount/F002C025tc525I1
cp Results/p.?.* functionals_etc/wall_mount/F002C025tc525I1
cp Results/p.0000{0,1,2,3,4,5,6,7}.vtk functionals_etc/wall_mount/F002C025tc525I1
cp para_functionalC025tc525I1.txt functionals_etc/wall_mount/F002C025tc525I1

sleep 30
../../bin/paraoscill 1 0.54 0.04 0.502 0.002 1> /dev/null 2>> timelog.txt
mkdir functionals_etc/wall_mount/F002C04tc54I1
cp Results/p.?.* functionals_etc/wall_mount/F002C04tc54I1
cp Results/p.0000{0,1,2,3,4,5,6,7}.vtk functionals_etc/wall_mount/F002C04tc54I1
cp para_functionalC04tc54I1.txt functionals_etc/wall_mount/F002C04tc54I1

sleep 30
../../bin/paraoscill 2 0.52 0.02 0.502 0.002 1> /dev/null 2>> timelog.txt
mkdir functionals_etc/wall_mount/F002C02tc52I2
cp Results/p.?.* functionals_etc/wall_mount/F002C02tc52I2
cp Results/p.0000{0,1,2,3,4,5,6,7}.vtk functionals_etc/wall_mount/F002C02tc52I2
cp para_functionalC02tc52I2.txt functionals_etc/wall_mount/F002C02tc52I2

sleep 30
../../bin/paraoscill 2 0.525 0.025 0.502 0.002 1> /dev/null 2>> timelog.txt
mkdir functionals_etc/wall_mount/F002C025tc525I2
cp Results/p.?.* functionals_etc/wall_mount/F002C025tc525I2
cp Results/p.0000{0,1,2,3,4,5,6,7}.vtk functionals_etc/wall_mount/F002C025tc525I2
cp para_functionalC025tc525I2.txt functionals_etc/wall_mount/F002C025tc525I2

sleep 30
../../bin/paraoscill 2 0.54 0.04 0.502 0.002 1> /dev/null 2>> timelog.txt
mkdir functionals_etc/wall_mount/F002C04tc54I2
cp Results/p.?.* functionals_etc/wall_mount/F002C04tc54I2
cp Results/p.0000{0,1,2,3,4,5,6,7}.vtk functionals_etc/wall_mount/F002C04tc54I2
cp para_functionalC04tc54I2.txt functionals_etc/wall_mount/F002C04tc54I2

sleep 30
../../bin/paraoscill 3 0.52 0.02 0.502 0.002 1> /dev/null 2>> timelog.txt
mkdir functionals_etc/wall_mount/F002C02tc52I3
cp Results/p.?.* functionals_etc/wall_mount/F002C02tc52I3
cp Results/p.0000{0,1,2,3,4,5,6,7}.vtk functionals_etc/wall_mount/F002C02tc52I3
cp para_functionalC02tc52I3.txt functionals_etc/wall_mount/F002C02tc52I3

sleep 30
../../bin/paraoscill 3 0.525 0.025 0.502 0.002 1> /dev/null 2>> timelog.txt
mkdir functionals_etc/wall_mount/F002C025tc525I3
cp Results/p.?.* functionals_etc/wall_mount/F002C025tc525I3
cp Results/p.0000{0,1,2,3,4,5,6,7}.vtk functionals_etc/wall_mount/F002C025tc525I3
cp para_functionalC025tc525I3.txt functionals_etc/wall_mount/F002C025tc525I3

sleep 30
../../bin/paraoscill 3 0.54 0.04 0.502 0.002 1> /dev/null 2>> timelog.txt
mkdir functionals_etc/wall_mount/F002C04tc54I3
cp Results/p.?.* functionals_etc/wall_mount/F002C04tc54I3
cp Results/p.0000{0,1,2,3,4,5,6,7}.vtk functionals_etc/wall_mount/F002C04tc54I3
cp para_functionalC04tc54I3.txt functionals_etc/wall_mount/F002C04tc54I3

#fine step 0.004
sleep 30
../../bin/paraoscill 1 0.525 0.025 0.504 0.004 1> /dev/null 2>> timelog.txt
mkdir functionals_etc/wall_mount/F004C02tc52I1
cp Results/p.?.* functionals_etc/wall_mount/F004C02tc52I1
cp Results/p.0000{0,1,2,3,4,5,6,7}.vtk functionals_etc/wall_mount/F004C02tc52I1
cp para_functionalC02tc52I1.txt functionals_etc/wall_mount/F004C02tc52I1

sleep 30
../../bin/paraoscill 1 0.54 0.04 0.504 0.004 1> /dev/null 2>> timelog.txt
mkdir functionals_etc/wall_mount/F004C04tc54I1
cp Results/p.?.* functionals_etc/wall_mount/F004C04tc54I1
cp Results/p.0000{0,1,2,3,4,5,6,7}.vtk functionals_etc/wall_mount/F004C04tc54I1
cp para_functionalC04tc54I1.txt functionals_etc/wall_mount/F004C04tc54I1

sleep 30
../../bin/paraoscill 2 0.525 0.025 0.504 0.004 1> /dev/null 2>> timelog.txt
mkdir functionals_etc/wall_mount/F004C025tc525I2
cp Results/p.?.* functionals_etc/wall_mount/F004C025tc525I2
cp Results/p.0000{0,1,2,3,4,5,6,7}.vtk functionals_etc/wall_mount/F004C025tc525I2
cp para_functionalC025tc525I2.txt functionals_etc/wall_mount/F004C025tc525I2

sleep 30
../../bin/paraoscill 2 0.54 0.04 0.504 0.004 1> /dev/null 2>> timelog.txt
mkdir functionals_etc/wall_mount/F004C04tc54I2
cp Results/p.?.* functionals_etc/wall_mount/F004C04tc54I2
cp Results/p.0000{0,1,2,3,4,5,6,7}.vtk functionals_etc/wall_mount/F004C04tc54I2
cp para_functionalC04tc54I2.txt functionals_etc/wall_mount/F004C04tc54I2

sleep 30
../../bin/paraoscill 3 0.525 0.025 0.504 0.004 1> /dev/null 2>> timelog.txt
mkdir functionals_etc/wall_mount/F004C025tc525I3
cp Results/p.?.* functionals_etc/wall_mount/F004C025tc525I3
cp Results/p.0000{0,1,2,3,4,5,6,7}.vtk functionals_etc/wall_mount/F004C025tc525I3
cp para_functionalC025tc525I3.txt functionals_etc/wall_mount/F004C025tc525I3

sleep 30
../../bin/paraoscill 3 0.54 0.04 0.504 0.004 1> /dev/null 2>> timelog.txt
mkdir functionals_etc/wall_mount/F004C04tc54I3
cp Results/p.?.* functionals_etc/wall_mount/F004C04tc54I3
cp Results/p.0000{0,1,2,3,4,5,6,7}.vtk functionals_etc/wall_mount/F004C04tc54I3
cp para_functionalC04tc54I3.txt functionals_etc/wall_mount/F004C04tc54I3
