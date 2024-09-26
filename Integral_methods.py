# Integral Methods

import math
import matplotlib.pyplot as plt
from scipy import integrate as scint


#purpose: obtains a float or integer number from the user, if it is not a number it will ask again until a number is entered
#assumptions: number entered is acceptable for the given prompt
#inputs: prompt - message displayed to the user instructing them what type of number to enter
#outputs: returns the number entered as a float
def num_input(prompt):
    try:
        num = float(input(prompt))
    except ValueError:
        num = num_input("Invalid input, please enter a number: ")
    return num

#purpose: determines the value of P(x), the probability density, of finding a particle at a certain position in a box of specific length and energy level 
#assumptions: the user has inputted constants for setting up parameters for particle in a box
#inputs: length - length of box region
#       eng_lvl - quantized energy level of the particle in the box
#       pos - position of particle in box within the box
#outputs:
def particle_box_prob(length, eng_lvl, pos):
    return (2 / length) * ((math.sin((eng_lvl * math.pi * pos) / length)) ** 2)


#purpose: generates two lists to represent the probability density as a function of position within a box of specifc length and energy level
#assumptions: the user has inputted constants for setting up parameters for particle in a box
#inputs: length - length of box region
#       eng_lvl - quantized energy level of the particle in the box
#       start - position within box to start list from
#       stop - position within box to stop list at
#outputs: returns list of probability density and corresponding position within box
def plot_particle_prob(length, eng_lvl, start, stop):
    prob_list = []
    pos_list = []
    pos = start
    while pos <= stop:
        pos_list.append(pos)
        prob_list.append(particle_box_prob(length, eng_lvl, pos))
        pos += length / 10000

    return pos_list, prob_list


#purpose: calculates the riemann sum using the midpoint rule over a given range for a box of specifc length and energy level
#assumptions: the user has inputted constants for setting up parameters for particle in a box
#inputs: num_steps - number of subintervals to break up the region of the box into 
#       length - length of the box region
#       eng_lvl - quantized energy level of the particle in the box
#       start - position within box to start sum from
#       stop - position within box to stop sum at
#outputs: returns the calculated area according to the midpoint rule
def simpleRiemann(num_steps, length, eng_lvl, start, stop):
    area = 0
    step_size = (stop - start) / num_steps
    pos = start + step_size
    while pos < stop:
        area += particle_box_prob(length, eng_lvl, pos) * step_size
        pos += step_size
    return area


#purpose: calculates the riemann sum using the trapezoid rule over a given range for a box of specifc length and energy level
#assumptions: the user has inputted constants for setting up parameters for particle in a box
#inputs: num_steps - number of subintervals to break up the region of the box into 
#       length - length of the box region
#       eng_lvl - quantized energy level of the particle in the box
#       start - position within box to start sum from
#       stop - position within box to stop sum at
#outputs: returns the calculated area according to the trapezoid rule
def trapezoidRule(num_steps, length, eng_lvl, start, stop):
    area = 0
    step_size = (stop - start) / num_steps
    pos = start
    i = 0
    while i <= num_steps:
        if i == 0 or i == num_steps:
            area += particle_box_prob(length, eng_lvl, pos)
            i += 1
            pos += step_size
        else:
            area += 2 * particle_box_prob(length, eng_lvl, pos)
            i += 1
            pos += step_size

    area = (step_size / 2) * area
    return area


#purpose: calculates the riemann sum using the trapezoid rule from scipy.integrate over a given range for a box of specifc length and energy level
#assumptions: the user has inputted constants for setting up parameters for particle in a box
#inputs: num_steps - number of subintervals to break up the region of the box into 
#       length - length of the box region
#       eng_lvl - quantized energy level of the particle in the box
#       start - position within box to start sum from
#       stop - position within box to stop sum at
#outputs: returns the calculated area according to the trapezoid rule from scipy
def scipyTrapezoid(num_steps, length, eng_lvl, start, stop):
    step_size = (stop - start) / num_steps
    pos_list, prob_list = plot_particle_prob(length, eng_lvl, start, stop)
    return scint.trapezoid(prob_list, pos_list, step_size)


#purpose: calculates the riemann sum using simpsons rule over a given range for a box of specifc length and energy level
#assumptions: the user has inputted constants for setting up parameters for particle in a box
#inputs: num_steps - number of subintervals to break up the region of the box into 
#       length - length of the box region
#       eng_lvl - quantized energy level of the particle in the box
#       start - position within box to start sum from
#       stop - position within box to stop sum at
#outputs: returns the calculated area according to simpsons rule
def simpsonsRule(num_steps, length, eng_lvl, start, stop):
    area = 0
    step_size = (stop - start) / num_steps
    pos = start
    i = 0
    while i <= num_steps:
        if i == 0 or i == num_steps:
            area += particle_box_prob(length, eng_lvl, pos)
            i += 1
            pos += step_size
        elif i % 2 == 1:
            area += 4 * particle_box_prob(length, eng_lvl, pos)
            i += 1
            pos += step_size
        else:
            area += 2 * particle_box_prob(length, eng_lvl, pos)
            i += 1
            pos += step_size
    area = (1 / 3) * step_size * area
    return area


#purpose: calculates the riemann sum using simpsons rule from scipy.integrate over a given range for a box of specifc length and energy level
#assumptions: the user has inputted constants for setting up parameters for particle in a box
#inputs: num_steps - number of subintervals to break up the region of the box into 
#       length - length of the box region
#       eng_lvl - quantized energy level of the particle in the box
#       start - position within box to start sum from
#       stop - position within box to stop sum at
#outputs: returns the calculated area according to simpsons rule from scipy
def scipySimpson(num_steps, length, eng_lvl, start, stop):
    step_size = (stop - start) / num_steps
    pos_list, prob_list = plot_particle_prob(length, eng_lvl, start, stop)
    return scint.simpson(prob_list, pos_list, step_size)


#purpose: displays the graph of probability density vs position within the box as well as integrating across a region specified by the user using the midpoint rule, trapezoid rule, simpsons rule and the scipy equivalents
#assumptions: user inputs appropiate values for specified inputs such as positive integers for energy levels or start and stop values within the length of the box
#inputs: none
#outputs: displays graph of probability density and prints the integral calculated according to different integration methods
def compare():
    length = num_input("Please enter the length of the box: ")
    eng_lvl = num_input("Please enter the energy level: ")
    num_steps = num_input("Please enter the number of steps for integration: ")
    start = num_input("Please enter from where you'd like to start the integral: ")
    stop = num_input("Please enter where you'd like to end the integral: ")

    plt.figure()
    pos_list, prob_list = plot_particle_prob(length, eng_lvl, 0, length)
    plt.plot(pos_list, prob_list)
    plt.xlabel("Particle Position")
    plt.ylabel("Probability Density")
    plt.show()

    print(
        "The calculated area using Midpoint Riemann sum is:",
        simpleRiemann(num_steps, length, eng_lvl, start, stop),
    )

    print(
        "The calculated area using Trapezoid rule is:",
        trapezoidRule(num_steps, length, eng_lvl, start, stop),
    )

    print(
        "The calculated area using Trapezoid rule from Scipy is:",
        scipyTrapezoid(num_steps, length, eng_lvl, start, stop),
    )

    print(
        "The calculated area using Simpsons rule is:",
        simpsonsRule(num_steps, length, eng_lvl, start, stop),
    )

    print(
        "The calculated area using Simpsons rule from Scipy is:",
        scipySimpson(num_steps, length, eng_lvl, start, stop),
    )
