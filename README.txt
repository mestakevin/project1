The main focus of this project is to display different techniques for solving a second order ordinary differential equation and numerically integrating an analytical function. For this project the ODE was a damped spring mass system and the analytical function the probability density for a particle in a 1D box. In order to run the code to display the correct results the following Python libraries are required:

matplotlib
scipy
math

Also ensure that the following files are in the same directory:

main.py
ODE_methods.py
Integral_methods.py

In order to run the code, run the main script by calling "python3 main.py". This will then call both compare() methods within the ODE_methods module and Integral_methods module. Following this call you will be asked to enter the mass of the spring mass system, the spring constant, the damping constant, the initial position and the initial velocity. Finally you will be asked to enter the time step to use to plot the displacement of the spring as a function of time. 

A sample input looks like this:
#########################################################################################################
Please enter a mass: 1
Please enter a spring constant: 1
Please enter a damping constant: 1
Please enter the initial position: 1
Please enter the initial velocity: 1
Please enter the time step: 0.01
#########################################################################################################
The allowed values for mass are positive floats
The allowed values for spring constants are positive floats
The allowed values for damping constants are positive floats
The allowed values for initial position are positive or negative floats
The allowed values for initial velocity are positive or negative floats
The allowed values for time steps are positive floats, preferably less than 0.1

Afterwards a series of 5 plots will be displayed showing how the displacement vs time plots appear for the two manual methods, explicit and RK4, as well as an analytically determined solution and an RK4 plot using the scipy library.

After closing the 5th plot the terminal will begin asking the user for inputs for the probability density problem. The user will be asked to enter the length of the one-dimensional box, the energy level, the number of steps(sub-intervals) for integration, the lower bound of the integration and the upper bound of the integration. 

A sample input looks like this:
#########################################################################################################
Please enter the length of the box: 1
Please enter the energy level: 1
Please enter the number of steps for integration: 1
Please enter from where you'd like to start the integral: 0
Please enter where you'd like to end the integral: 1
#########################################################################################################

The allowed values for length are positive floats
The allowed values for energy level are positive integers
The allowed values for number of steps are positive integers
The allowed values for starting point of the integral is any floating number between 0 and the length, inclusive
The allowed values for ending point of the integral is any floating point greater than or equal to the starting point of the integral but less than or equal to the length

Afterwards a plot of the probability density will appear to help visual what the function that you are integrating looks like. Other than that it does not represent anything.

Once you close the plot, the terminal will display calculated integrals using the three manual rules, midpoint, trapezoidal, and simpson's, as well as the trapezoidal and simpson's rule from the scipy library.

At this point the program will be terminated and you can call it again changing the initial conditions for the ODE or the parameters of the box for the analytical function to obtain new results
