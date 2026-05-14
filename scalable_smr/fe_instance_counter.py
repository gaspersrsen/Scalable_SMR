import openmc
import numpy as np

def filled(x,y,r):
    return x**2 + y**2 - r**2 <= r-1

def fatfilled(x,y,radius):
	return filled(x, y, radius) and not (
		filled(x + 1, y, radius) and
		filled(x - 1, y, radius) and
		filled(x, y + 1, radius) and
		filled(x, y - 1, radius) and
		filled(x + 1, y + 1, radius) and
		filled(x + 1, y - 1, radius) and
		filled(x - 1, y - 1, radius) and
		filled(x - 1, y + 1, radius))

def get_n_fe(d):
    n_instances = 0
    list_i=[]
    r = (d-1) / 2
    j = -int(np.floor(d/2))
    for y in np.linspace(-r,r,d):
        i = -int(np.floor(d/2))
        found_edge = False
        for x in np.linspace(-r,r,d):
            if x >= 0 and y >= 0 and y <= x:
                if filled(x,y,r):
                    n_fei = 8
                    if x==0: n_fei/=2
                    if y==0: n_fei/=2
                    if x==y: n_fei/=2
                    n_instances += int(n_fei)
    return n_instances