#Integral Methods

import math
import matplotlib.pyplot as plt
from scipy import integrate as scint

def numInput(prompt):
        try:
            num = float(input(prompt))
        except ValueError:
            num = numInput("Invalid input, please enter a number: ")
        return num


def particleBoxProb(length, eng_lvl, pos):
        return (2 / length) * ( (math.sin( (eng_lvl * math.pi * pos) / length)) ** 2)


def plotParticleProb(length, eng_lvl, start, stop):
        Prob_list = []
        pos_list = []
        pos = start 
        while pos <= stop:
            pos_list.append(pos)
            Prob_list.append(particleBoxProb(length, eng_lvl, pos))
            pos += length / 10000


        return pos_list, Prob_list


def simpleRiemann(num_steps, length, eng_lvl, start, stop):
        area = 0 
        step_size = (stop - start) / num_steps
        pos = start + step_size
        while pos < stop: 
            area += particleBoxProb(length, eng_lvl, pos) * step_size
            pos += step_size
        return area

def trapezoidRule(num_steps, length, eng_lvl, start, stop):
        area = 0 
        step_size = (stop - start) / num_steps
        pos = start
        i = 0 
        while i <= num_steps:
            if i == 0 or i == num_steps:
                area += particleBoxProb(length, eng_lvl, pos)
                i += 1
                pos += step_size
            else:
                area += 2 * particleBoxProb(length, eng_lvl, pos)
                i += 1
                pos += step_size

        area = (step_size / 2) * area
        return area

def scipyTrapezoid(num_steps, length, eng_lvl, start, stop):
        step_size = (stop - start) / num_steps
        pos_list, Prob_list = plotParticleProb(length, eng_lvl,start,stop)
        return scint.trapezoid(Prob_list, pos_list, step_size)


def simpsonsRule(num_steps, length, eng_lvl, start, stop):
        area = 0 
        step_size = (stop - start) / num_steps
        pos = start
        i = 0 
        while i <= num_steps:
            if i == 0 or i == num_steps:
                area += particleBoxProb(length, eng_lvl, pos)
                i += 1
                pos += step_size
            elif i%2 == 1:
                area += 4 * particleBoxProb(length, eng_lvl, pos)
                i += 1
                pos += step_size
            else:
                area += 2* particleBoxProb(length, eng_lvl, pos)
                i += 1
                pos += step_size
        area = (1/3) * step_size * area
        return area

def scipySimpson(num_steps, length, eng_lvl, start, stop):
        step_size = (stop - start) / num_steps
        pos_list, Prob_list = plotParticleProb(length, eng_lvl, start, stop)
        return scint.simpson(Prob_list, pos_list, step_size)

def compare():
        length = numInput("Please enter the length of the box: ")
        eng_lvl = numInput("Please enter the energy level: ")
        num_steps = numInput("Please enter the number of steps for integration: ")
        start = numInput("Please enter from where you'd like to start the integral: ")
        stop = numInput("Please enter where you;d like to end the integral: ")

        plt.figure()
        
        pos_list, Prob_list = plotParticleProb(length, eng_lvl, 0,length)

        plt.plot(pos_list, Prob_list)
        plt.xlabel('Particle Position')
        plt.ylabel('Probability Density')
        plt.show()


        print("The calculated area using Midpoint Riemann sum is:", simpleRiemann(num_steps, length, eng_lvl, start, stop) ) 

        print("The calculated area using Trapezoid rule is:", trapezoidRule(num_steps, length, eng_lvl, start, stop) ) 
            
        print("The calculated area using Trapezoid rule from Scipy is:", scipyTrapezoid(num_steps, length, eng_lvl, start, stop) ) 

        print("The calculated area using Simpsons rule is:", simpsonsRule(num_steps, length, eng_lvl, start, stop) ) 
        
        print("The calculated area using Simpsons rule from Scipy is:", scipySimpson(num_steps, length, eng_lvl, start, stop) ) 
