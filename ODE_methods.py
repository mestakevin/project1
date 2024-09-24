#ODE Integrators
import math
import matplotlib.pyplot as plt
import scipy.integrate as scpy


def numInput(prompt):
        try:
            num = float(input(prompt))
        except ValueError:
            num = numInput("Invalid input, please enter a number: ")
        return num


def getForce(cur_x, cur_vel, k_cons, damp_cons):
        force = (-k_cons * cur_x) - (damp_cons * cur_vel)
        return force

#Explict Method
def explictMethod(mass, cur_x, cur_vel, force, dt):
        new_vel = cur_vel + (force/mass) * dt
        new_x   = cur_x   + cur_vel * dt
        return new_x, new_vel

#Fourth Order Runge-Kutta

def RK4(mass, cur_x, cur_vel, k_cons, damp_cons, dt):
        k1x = func1RK4(cur_x, cur_vel)
        k1v = func2RK4(cur_x, cur_vel, k_cons, damp_cons, mass)
        
        k2x = func1RK4(cur_x + (dt * k1x / 2), cur_vel + (dt * k1v / 2 ) )
        k2v = func2RK4(cur_x + (dt * k1x / 2), cur_vel + (dt * k1v / 2 ), k_cons, damp_cons, mass )

        k3x = func1RK4(cur_x + (dt * k2x / 2), cur_vel + (dt * k2v / 2 ) )
        k3v = func2RK4(cur_x + (dt * k2x / 2), cur_vel + (dt * k2v / 2 ), k_cons, damp_cons, mass )

        k4x = func1RK4(cur_x + (dt * k3x / 2), cur_vel + (dt * k3v / 2 ) )
        k4v = func2RK4(cur_x + (dt * k3x / 2), cur_vel + (dt * k3v / 2 ), k_cons, damp_cons, mass )

        new_x = cur_x + (dt/6) * (k1x + 2*k2x + 2*k3x + k4x) 
        new_v = cur_vel + (dt/6) * (k1v + 2*k2v + 2*k3v + k4v)
        return new_x, new_v

def func1RK4(x, v):
        return v

def func2RK4(x, v, k_cons, damp_cons, mass):
        return ( ( (-k_cons/mass) * x) - ( (damp_cons/mass) * v) )



def RK4Spring(mass, x0, v0, k_cons, damp_cons, dt):
        pos_list = []
        time_list = []
        iteration = 0
        time = 0
        cur_x = x0
        cur_vel = v0

        while iteration < 5000:
            
            pos_list.append(cur_x)
            time_list.append(time)
            

            cur_x, cur_vel = RK4(mass, cur_x, cur_vel, k_cons, damp_cons, dt)
        
            time += dt
            iteration += 1
        
        return pos_list, time_list

def explicitSpring(mass, x0, v0, k_cons, damp_cons, dt):
        pos_list = []
        time_list = []
        iteration = 0
        time = 0
        cur_x = x0
        cur_vel = v0

        while iteration < 5000:
            
            pos_list.append(cur_x)
            time_list.append(time)
            
            force = getForce(cur_x, cur_vel, k_cons, damp_cons)

            cur_x, cur_vel = explictMethod(mass, cur_x, cur_vel, force, dt)
        
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
        #beta = ( ( (4 * mass * k_cons - (damp_cons ** 2) ) / (4* ( mass ** 2) ) ) ** 0.5 )
        beta = ( (k_cons/mass) - (alp**2) ) ** 2 
      
        while iteration < 5000:
            pos_list.append(cur_x)
            time_list.append(time)
            
            time += dt
            cur_x = (math.exp(alp * time) )*(x0 * math.cos(beta * time ) + ((v0 - alp * x0) / beta) * math.sin(beta * time) )
            iteration += 1        
      
        return pos_list, time_list

def compare():
        mass = numInput("Please enter a mass: ")         
        k_cons = numInput("Please enter a spring constant: ")
        damp_cons  = numInput("Please enter a damping constant: ")

        x_init = numInput("Please enter the initial position: ")
        v_init = numInput("Please enter the initial velocity: ")
        
        dt = numInput("Please enter the time step: ")
        

        plt.figure()
        
        pos_list, time_list = explicitSpring(mass, x_init, v_init, k_cons, damp_cons, dt)
        plt.plot(time_list, pos_list, label = 'Explict')
       
        pos_list, time_list = RK4Spring(mass, x_init, v_init, k_cons, damp_cons, dt)
        plt.plot(time_list, pos_list, label = 'RK4')

        pos_list, time_list = analyticalSpring(mass, x_init, v_init, k_cons, damp_cons, dt)
        plt.plot(time_list, pos_list, label = 'Analytical')


        plt.legend()
        plt.xlabel('Time')
        plt.ylabel('X displacement')
        plt.show()

