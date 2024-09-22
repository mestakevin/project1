#Integral Methods

import math
import matplotlib.pyplot as plt

def numInput(prompt):
        try:
            num = float(input(prompt))
        except ValueError:
            num = numInput("Invalid input, please enter a number: ")
        return num


def particleBoxProb(length, eng_lvl, pos):
        return (2/length) * (math.sin( eng_lvl * math.pi * pos / length) ** 2)

#Simple Riemann sum

def plotParticleProb(num_steps, length, eng_lvl):
        Prob_list = []
        pos_list = []
        pos = 0 
        while pos <= length:
            pos_list.append(pos)
            Prob_list.append(particleBoxProb(length, eng_lvl, pos))
            pos += length / num_steps

        return pos_list, Prob_list


def simpleRiemann(num_steps, domain, eng_lvl):
        area = 0 
        step_size = domain / num_steps
        pos = step_size / 2
        while pos < domain: 
            area += particleBoxProb(domain, eng_lvl, pos) * step_size
            pos += step_size
        return area

def integrateProb():
        length = numInput("Please enter the length of the box: ")
        eng_lvl = numInput("Please enter the energy level: ")
        num_steps = numInput("Please enter the number of steps for integration: ")

        plt.figure()
        
        pos_list, Prob_list = plotParticleProb(num_steps, length, eng_lvl)
        plt.plot(pos_list, Prob_list)
        plt.xlabel('Particle Position')
        plt.ylabel('Probability Density')
        plt.show()


        print("The calculated area using these inputs is:", simpleRiemann(num_steps, length, eng_lvl) ) 


integrateProb()
