# ODE Integrators
import math
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

#purpose: obtains a float or integer number from the user, if it is not a number it will ask again until a number is entered
#assumptions: number entered is acceptable for the given prompt
#inputs: prompt - message displayed to the user instructing them what type of number to enter
#outputs: returns the number entered as a float
object on spring has begun moving and constants have been initalized
def num_input(prompt):
    try:
        num = float(input(prompt))
    except ValueError:
        num = num_input("Invalid input, please enter a number: ")
    return num

#purpose: calculates the force experienced by the object on the spring
#assumptions: object on spring has begun moving and constants have been initalized  
#inputs: cur_x - current displacement of object from equilibrium position
#       cur_vel - current velocity of the object relative to equilibium position
#       K_CONS - spring constant
#       DAMP_CONS - damping constant
#outputs: returns the force experienced by the object

def get_force(cur_x, cur_vel, K_CONS, DAMP_CONS):
    force = (-K_CONS * cur_x) - (DAMP_CONS * cur_vel)
    return force


# Explict Method
#purpose: updates the velocity and position of the object on the spring using the explicit method
#assumption: object on spring has begun moving and constants have been initalized
#inputs: MASS - mass of object on the spring
#       cur_x - current displacement of object from equilibrium position
#       cur_vel - current velocity of the object relative to equilibrium position
#       force - current force experienced by object
#       dt - time step
#outputs: returns the new displacement and new velocity
def explict_method(MASS, cur_x, cur_vel, force, dt):
    new_vel = cur_vel + (force / MASS) * dt
    new_x = cur_x + cur_vel * dt
    return new_x, new_vel


# Fourth Order Runge-Kutta
#purpose: updates the velocity and the position of the object on the spring using the RK4 method
#assumption: object on spring has begun moving and constants have been initalized, functions for RK4 method have been created
#inputs: MASS - mass of object on the spring
#       cur_x - current displacement of object from equilibrium position
#       cur_vel - current velocity of the object relative to equilibirum position
#       force - current force experienced by object
#       dt - time step
#outputs: returns the new displacement and new velocity
def kevin_RK4(MASS, cur_x, cur_vel, K_CONS, DAMP_CONS, dt):
    k1x = func1_RK4(cur_x, cur_vel)
    k1v = func2_RK4(cur_x, cur_vel, K_CONS, DAMP_CONS, MASS)
    k2x = func1_RK4(cur_x + (dt * k1x / 2), cur_vel + (dt * k1v / 2))
    k2v = func2_RK4(
        cur_x + (dt * k1x / 2), cur_vel + (dt * k1v / 2), K_CONS, DAMP_CONS, MASS
    )
    k3x = func1_RK4(cur_x + (dt * k2x / 2), cur_vel + (dt * k2v / 2))
    k3v = func2_RK4(
        cur_x + (dt * k2x / 2), cur_vel + (dt * k2v / 2), K_CONS, DAMP_CONS, MASS
    )
    k4x = func1_RK4(cur_x + (dt * k3x / 2), cur_vel + (dt * k3v / 2))
    k4v = func2_RK4(
        cur_x + (dt * k3x / 2), cur_vel + (dt * k3v / 2), K_CONS, DAMP_CONS, MASS
    )

    new_x = cur_x + (dt / 6) * (k1x + 2 * k2x + 2 * k3x + k4x)
    new_v = cur_vel + (dt / 6) * (k1v + 2 * k2v + 2 * k3v + k4v)
    return new_x, new_v


#purpose: first function for the kevin_RK4 function that returns value of dx/dt
#assumption: object on spring has begun moving and constants have been initalized and kevin_RK4 method is being run
#inputs: x - current displacement of object from equilibrium position
#       v - current velocity of object relative to equilibrium position
#outputs: returns the value of dx/dt which is velocity
def func1_RK4(x, v):
    return v


#purpose: second function for the kevin_RK4 function that returns value of dv/dt
#assumption: object on spring has begun moving and constants have been initalized and kevin_RK4 method is being run
#inputs: x - current displacement of object from equilibrium position
#       v - current velocity of object relative to equilibrium position       
#       K_CONS - spring constant
#       DAMP_CONS - damping constant
#       MASS - mass of object on the spring
#outputs: returns the value of dv/dt which is acceleration
def func2_RK4(x, v, K_CONS, DAMP_CONS, MASS):
    return ((-K_CONS / MASS) * x) - ((DAMP_CONS / MASS) * v)


#purpose: calculates dx/dt and dv/dt for the solve_ivp function from scipy.integrate written in format for RK4 method 
#assumption: scipy.integrate has been imported, constants have also been intialized
#inputs: t - current time
#       y - vector with displacement of mass from equilibrium and velocity of mass
#       K_CONS - spring constant
#       DAMP_CONS - damping constant
#       MASS - mass of object on the spring
#outputs: returns dx/dt and dv/dt in a list
def scipy_RK45(t, y, K_CONS, DAMP_CONS, MASS):
    x, v = y
    dxdt = v
    dvdt = ((-K_CONS / MASS) * x) - ((DAMP_CONS / MASS) * v)
    return [dxdt, dvdt]


#purpose: oscillates spring using own RK4 method
#assumption: user has inputted constants for spring simulation
#inputs: MASS - mass of object on the spring
#       x0 - initial displacement of object from equilibrium position
#       v0 - initial velocity of object relative to equilibrium position
#       K_CONS - spring constant
#       DAMP_CONS - damping constant
#       dt - time step
#outputs: returns list of displacements of object from equilibirium position and corresponding time list
def RK4Spring(MASS, x0, v0, K_CONS, DAMP_CONS, dt):
    pos_list = []
    time_list = []
    iteration = 0
    time = 0
    cur_x = x0
    cur_vel = v0
    while iteration < 5000:
        pos_list.append(cur_x)
        time_list.append(time)
        cur_x, cur_vel = kevin_RK4(MASS, cur_x, cur_vel, K_CONS, DAMP_CONS, dt)
        time += dt
        iteration += 1

    return pos_list, time_list


#purpose: oscillates spring using own explicit method
#assumption: user has inputted constants for spring simulation
#inputs: MASS - mass of object on the spring
#       x0 - initial displacement of object from equilibrium position
#       v0 - initial velocity of object relative to equilibrium position
#       K_CONS - spring constant
#       DAMP_CONS - damping constant
#       dt - time step
#outputs: returns list of displacements of object from equilibirium position and corresponding time list
def explicitSpring(MASS, x0, v0, K_CONS, DAMP_CONS, dt):
    pos_list = []
    time_list = []
    iteration = 0
    time = 0
    cur_x = x0
    cur_vel = v0

    while iteration < 5000:
        pos_list.append(cur_x)
        time_list.append(time)
        force = get_force(cur_x, cur_vel, K_CONS, DAMP_CONS)
        cur_x, cur_vel = explict_method(MASS, cur_x, cur_vel, force, dt)
        time += dt
        iteration += 1

    return pos_list, time_list


#purpose: oscillates spring using analytically determined solution for four cases, harmonic oscillator, underdamped, overdamped or critically damped
#assumption: user has inputted constants for spring simulation
#inputs: MASS - mass of object on the spring
#       x0 - initial displacement of object from equilibrium position
#       v0 - initial velocity of object relative to equilibrium position
#       K_CONS - spring constant
#       DAMP_CONS - damping constant
#outputs: returns list of displacements of object from equilibirium position and corresponding time list
def analyticalSpring(MASS, x0, v0, K_CONS, DAMP_CONS, dt):
    pos_list = []
    time_list = []
    iteration = 0
    time = 0
    cur_x = x0
    alp = -DAMP_CONS / (2 * MASS)
    beta = ((K_CONS / MASS) - (alp**2)) ** 0.5

    while iteration < 5000:
        pos_list.append(cur_x)
        time_list.append(time)
        time += dt
        if beta == 0:  # Critically damped
            cur_x = (math.exp(alp * time)) * (x0)
        elif (K_CONS / MASS - alp**2) < 0:  # Over damped
            root1 = (-DAMP_CONS - (DAMP_CONS**2 - 4 * MASS * K_CONS) ** 0.5) / (
                2 * MASS
            )
            root2 = (-DAMP_CONS + (DAMP_CONS**2 - 4 * MASS * K_CONS) ** 0.5) / (
                2 * MASS
            )
            b_cons = (v0 - x0 * root1) / (root2 - root1)
            a_cons = x0 - b_cons
            cur_x = a_cons * math.exp(root1 * time) + b_cons * math.exp(root2 * time)
        else:  # Under damped
            cur_x = (math.exp(alp * time)) * (
                x0 * math.cos(beta * time)
                + ((v0 - alp * x0) / beta) * math.sin(beta * time)
            )
        iteration += 1

    return pos_list, time_list

#purpose: plots two different methods on the same figure
#assumption: at least two different spring functions have been run and the outputs stored
#inputs: pos_list1 - displacement of spring from equilibrium position using method 1
#       time_list1 - associated time list for pos_list1
#       label1 - name of method 1
#       pos_list1 - displacement of spring from equilibrium position using method 2
#       time_list - associated time list for pos_list2
#       label2 - name of method 2
#outputs: plot with both methods calculated displacements labeled accordingly

def plotDataSets(pos_list1, time_list1, label1, pos_list2, time_list2, label2):
    plt.figure()
    plt.plot(time_list1, pos_list1, label=label1)
    plt.plot(time_list2, pos_list2, label=label2)
    plt.legend()
    plt.xlabel("Time")
    plt.ylabel("X displacement")
    plt.show()


#purpose: performs all four methods to determine object on a damped spring's position over time and displays comparison of explicit method with analytical solution and own RK4 method with scipy's RK4 method
#assumption: user will only input values that are appropiate, such as a positive mass/spring constant
#inputs: none
#outputs: returns two plots, one with the explicit and analytical methods and another with own RK4 and scipy RK4 methods
def compare():
    MASS = num_input("Please enter a mass: ")
    K_CONS = num_input("Please enter a spring constant: ")
    DAMP_CONS = num_input("Please enter a damping constant: ")
    x_init = num_input("Please enter the initial position: ")
    v_init = num_input("Please enter the initial velocity: ")
    dt = num_input("Please enter the time step: ")

    ex_pos, ex_time = explicitSpring(MASS, x_init, v_init, K_CONS, DAMP_CONS, dt)
    an_pos, an_time = analyticalSpring(MASS, x_init, v_init, K_CONS, DAMP_CONS, dt)
    rk_pos, rk_time = RK4Spring(MASS, x_init, v_init, K_CONS, DAMP_CONS, dt)

    plotDataSets(ex_pos, ex_time, "Explicit", an_pos, an_time, "Analytical")

    t_span = [0.0, dt * 5000]
    y0 = [x_init, v_init]
    scipy_rk4 = solve_ivp(
        scipy_RK45,
        t_span,
        y0,
        method="RK45",
        args=(K_CONS, DAMP_CONS, MASS),
        max_step=dt,
    )
    scipy_time = scipy_rk4.t
    scipy_pos = scipy_rk4.y[0]

    plotDataSets(rk_pos, rk_time, "Kevin's RK4", scipy_pos, scipy_time, "Scipy's RK4")
