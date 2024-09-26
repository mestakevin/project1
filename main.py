# import modules here

import matplotlib.pyplot as plt
import ODE_methods as ode
import Integral_methods as integ
import math

#purpose: runs the main compare functions of ODE_methods and Integral_methods and produces desired comparions for different methods
#assumptions: the modules math, scipy.integrate and matplotlib.plot have been installed, ODE_methods.py and Integral_methods.py are in the same directory as main.py
#inputs: none
#outputs: plots the displacement vs time plots for the ODE problem and prints the area of the integral of the analytical function using different integration techniques

def main():
    ode.compare()
    integ.compare()


main()
