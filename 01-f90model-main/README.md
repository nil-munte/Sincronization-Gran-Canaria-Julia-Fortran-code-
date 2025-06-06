# power-grid


This repository contains the program used for the simulations presented in:

María Martínez-Barbeito, Damià Gomila, & Pere Colet (2022). Dynamical model for power grid frequency fluctuations: application to islands with high penetration of wind generation. Submitted for publication.

In what follows, we give an explanation for the files. For the details on the equations and the specific parameter values considered, we refer the reader to the manuscript.



Main program: grid_implicit_sparse.f90 

Program written by María Martínez-Barbeito and Pere Colet.
- Integrates power grid dynamics using a first order semi-implicit Euler method: x(t+dt)=x(t)+dt*F(x(t+dt)) approximated as x(t+dt)=x(t)+dx with dx solution of the linear set of equations (I-dt*J)*dx=F(x)*dt where J is the Jacobian.
- Takes advantage of the sparse structure of the grid (<2% non-zero terms in I-dtJ). The program uses PARDISO routines included in Intel Math Kernel Libraries (MKL) to solve the sparse linear set of equations. There is also an open source version as part of the PARDISO project: https://www.pardiso-project.org/
- Uses RANDOM_NUMBER and Box-Muller algorithm to obtain Gaussian random numbers.



Input file: parameters.dat 

- Contains all parameter values to be used in the main program, as well as file names.


Input file: dispatch/dispatch_dd-mm-yyyy.dat

- dispatch_07-02-2018.dat and dispatch_08-02-2018.dat given as examples. These files can be used to reproduce Figure 3 from the paper.
- Dispatch files were created using 10-minute power data from Red Eléctrica de España (https://demanda.ree.es/visiona/home) for Gran Canaria.
- Each file contains 156 dispatches. There is 1 dispatch/10 minutes. Dispatch 1 is equivalent to 23:00 (11 pm) from the previous day, and dispatch 156 is equivalent to 01:00 (1 am) from the next day (UTC time).
- Each dispatch contains:
    - The load per node in MW.
    - The operation set point (Pref) of each generating group in MW.
    - The timescale of Pref (lambda) for each plant in minutes.
    - The power provided by renewables (assets) in MW, which is assumed to be 0 if it is not in the file.


Output file: results/frequency_dd-mm-yyyy.dat

- frequency_07-02-2018.dat and frequency_08-02-2018.dat given as examples. These files can be used to reproduce Figure 3 from the paper.
- Frequency files contain the frequency values from simulations written every 1 second. Timestamp is not given explicitly, but it can be easily deduced.

