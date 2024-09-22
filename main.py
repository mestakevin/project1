#import modules here

import matplotlib.pyplot as plt
import ODE_methods as ode
import Integral_methods as integ
import math

def numInput(prompt):
        try:
            num = float(input(prompt))
        except ValueError:
            num = numInput("Invalid input, please enter a number: ")
        return num

def main():
        ode.compare()
        

main()
