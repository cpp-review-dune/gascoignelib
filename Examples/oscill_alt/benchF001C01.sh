../../bin/Examples/oscill_alt/paraoscill 1 0.55 0.05 0.501 0.001 1> /dev/null 2>> timelog.txt
mkdir Plots_and_data/paraoscill/C05tc55I1
cp Results/p.*.bup Plots_and_data/paraoscill/C05tc55I1
cp Results/p.*.vtk Plots_and_data/paraoscill/C05tc55I1
cp para_functionalC05tc55I1.txt Plots_and_data/paraoscill/C05tc55I1

sleep 10
../../bin/Examples/oscill_alt/para_paraoscill 2 0.51 0.01 0.501 0.001 1> /dev/null 2>> timelog.txt
mkdir Plots_and_data/paraoscill/C01tc51I2
cp Results/p.*.bup Plots_and_data/paraoscill/C01tc51I2
cp Results/p.*.vtk Plots_and_data/paraoscill/C01tc51I2
cp para_functionalC51tc51I2.txt Plots_and_data/paraoscill/C01tc51I2

sleep 10
../../bin/Examples/oscill_alt/para_paraoscill 3 0.51 0.01 0.501 0.001 1> /dev/null 2>> timelog.txt
mkdir Plots_and_data/paraoscill/C01tc51I3
cp Results/p.*.bup Plots_and_data/paraoscill/C01tc51I3
cp Results/p.*.vtk Plots_and_data/paraoscill/C01tc51I3
cp para_functionalC51tc51I3.txt Plots_and_data/paraoscill/C01tc51I3
