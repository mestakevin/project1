#import modules here
import matplotlib.pyplot as plt
import ODE_methods as ode
import math

#calculates force
def getForce(cur_x, cur_vel, k_cons, damp_cons):
        force = (-k_cons * cur_x) - (damp_cons * cur_vel)
        return force


def oscillateSpring(mass, x0, v0, k_cons, damp_cons, dt):
        pos_list = []
        time_list = []
        iteration = 0
        time = 0
        cur_x = x0
        cur_vel = v0

        while iteration < 800:
            
            pos_list.append(cur_x)
            time_list.append(time)
            
            force = getForce(cur_x, cur_vel, k_cons, damp_cons)

            cur_x, cur_vel = ode.explictMethod(mass, cur_x, cur_vel, force, dt)
        
            time += dt
            iteration += 1

        return pos_list, time_list

def analyticalSpring(mass, x0, v0, k_cons, damp_cons, dt): 
        pos_list = []
        time_list = []
        iteration = 0
        time = 0

        cur_x = x0         
        alp =  (-damp_cons / (2* mass))
        beta = ( ( (4 * mass * k_cons - (damp_cons ** 2) ) / (4* ( mass ** 2) ) ) ** 0.5 )
    
      
        while iteration < 800:
            pos_list.append(cur_x)
            time_list.append(time)
            
            time += dt
            cur_x = (math.exp(alp) )*(x0 * math.cos(beta * time ) + ((v0 - alp * x0) / beta) * math.sin(beta * time) )
            iteration += 1        
      

        return pos_list, time_list


def main():
        mass = ode.floatInput("Please enter a mass: ")         
        k_cons = ode.floatInput("Please enter a spring constant: ")
        damp_cons  = ode.floatInput("Please enter a damping constant: ")

        x_init = ode.floatInput("Please enter the initial position: ")
        v_init = ode.floatInput("Please enter the initial velocity: ")
        
        dt = ode.floatInput("Please enter the time step: ")
        
        #pos_list, time_list = oscillateSpring(mass, x_init, v_init, k_cons, damp_cons, dt)
        #pos_list, time_list = analyticalSpring(mass, x_init, v_init, k_cons, damp_cons, dt)

        plt.figure()
        
        pos_list, time_list = oscillateSpring(mass, x_init, v_init, k_cons, damp_cons, dt)
        plt.plot(time_list, pos_list, label = 'Explict')
       
        pos_list, time_list = analyticalSpring(mass, x_init, v_init, k_cons, damp_cons, dt)
        plt.plot(time_list, pos_list, label = 'Analytical')
        
        plt.xlabel('Time')
        plt.ylabel('X displacement')
        plt.show()

#run main
main()
