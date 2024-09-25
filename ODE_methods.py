# ODE Integrators
import math
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp


def num_input(prompt):
    try:
        num = float(input(prompt))
    except ValueError:
        num = num_input("Invalid input, please enter a number: ")
    return num


def get_force(cur_x, cur_vel, K_CONS, DAMP_CONS):
    force = (-K_CONS * cur_x) - (DAMP_CONS * cur_vel)
    return force


# Explict Method
def explict_method(MASS, cur_x, cur_vel, force, dt):
    new_vel = cur_vel + (force / MASS) * dt
    new_x = cur_x + cur_vel * dt
    return new_x, new_vel


# Fourth Order Runge-Kutta


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


def func1_RK4(x, v):
    return v


def func2_RK4(x, v, K_CONS, DAMP_CONS, MASS):
    return ((-K_CONS / MASS) * x) - ((DAMP_CONS / MASS) * v)


def scipy_RK45(t, y, K_CONS, DAMP_CONS, MASS):
    x, v = y
    dxdt = v
    dvdt = ((-K_CONS / MASS) * x) - ((DAMP_CONS / MASS) * v)
    return [dxdt, dvdt]


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


def plotDataSets(pos_list1, time_list1, label1, pos_list2, time_list2, label2):
    plt.figure()
    plt.plot(time_list1, pos_list1, label=label1)
    plt.plot(time_list2, pos_list2, label=label2)
    plt.legend()
    plt.xlabel("Time")
    plt.ylabel("X displacement")
    plt.show()


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
