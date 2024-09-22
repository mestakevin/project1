#ODE Integrators


def floatInput(prompt):
        num = float(input(prompt))
        return num




#Explict Method
def explictMethod(mass, cur_x, cur_vel, force, dt):
        new_vel = cur_vel + (force/mass) * dt
        new_x   = cur_x   + cur_vel * dt
        return new_x, new_vel

#Fourth Order Runge-Kutta


