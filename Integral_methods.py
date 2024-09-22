#Integral Methods

#Simple Riemann sum

def simpleRiemann(num_steps, domain, function):
        area = 0 
        step_size = domain / num_steps
        pos = step_size / 2
        while pos < domain: 
            area += function(pos) * step_size
            pos += step_size
        return area



#Trapezoidal Rule

#Simpsons Rule


