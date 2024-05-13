#Brownian Motion
#Forbes Li and Thomas Nguyen
import math
import os
import matplotlib
import matplotlib.pyplot as plt
matplotlib.use('TkAgg',force=True)
import numpy as np
from scipy.stats import rayleigh
from scipy.optimize import curve_fit

f1 = open(os.path.dirname(os.path.abspath(__file__))+"/Bead1.txt", "r")
f2 = open(os.path.dirname(os.path.abspath(__file__))+"/Bead2.txt", "r")
f3 = open(os.path.dirname(os.path.abspath(__file__))+"/Bead3.txt", "r")
f4 = open(os.path.dirname(os.path.abspath(__file__))+"/Bead4.txt", "r")
f5 = open(os.path.dirname(os.path.abspath(__file__))+"/Bead5.txt", "r")
f6 = open(os.path.dirname(os.path.abspath(__file__))+"/Bead6.txt", "r")
f7 = open(os.path.dirname(os.path.abspath(__file__))+"/Bead7.txt", "r")
f8 = open(os.path.dirname(os.path.abspath(__file__))+"/Bead8.txt", "r")
f9 = open(os.path.dirname(os.path.abspath(__file__))+"/Bead9.txt", "r")
f10 = open(os.path.dirname(os.path.abspath(__file__))+"/Bead10.txt", "r")
total_data = [f1,f2,f3,f4,f5,f6,f7,f8,f9,f10]
#print(f1)

def generate_ri(total_data):
    ri_2 = []
    d_r = []
    
    dr = 0
    dx = 1.9 * 10 **(-7)
    m1 = 0
    m2 = 0
    print(len(total_data))
    r_2_ns = []
    for i in range(len(total_data)):
        a1 = 0
        x = 0
        y = 0
        r = 0
        #D = 0
        r_2 = 0
        #r_2_ns = []
        x_a = []
        y_a = []
        
        
        a = str(total_data[i].read())
        a = a.replace('X (pixels)  	Y (pixels)','')
        a = a[a.find('\n') + 1:]
        a = a[a.find('\n') + 1:]
        while(a.count('\n') > 0):
            a1 = a[0:a.find('\n')]
            a = a[a.find('\n') + 1:]
            x = float(a1[:a1.find('\t')]) * 0.1155 * 10**(-6)
            y = float(a1[a1.find('\t') + 1:]) * 0.1155 * 10**(-6)
            x_a.append(x)
            y_a.append(y)
        a1 = a[0:a.find('\n')]
        x = float(a1[:a1.find('\t')]) * 0.1155 * 10**(-6)
        y = float(a1[a1.find('\t') + 1:]) * 0.1155 * 10**(-6)
        x_a.append(x)
        y_a.append(y)
        for n in range(len(x_a) - 1):
            x = x_a[n + 1] - x_a[n]
            y = y_a[n + 1] - y_a[n]
            r = x**2 + y**2
            r_2 = r**(1/2)
            m1 = 2*(2*x*dx)**(2)
            m2 = 2*(2*y*dx)**(2)
            r_2_ns.append(r_2)
            ri_2.append(r)
            
    
    
    
    #print(ri_2)
    print(len(r_2_ns))
    print(len(ri_2))
    
    maximum_likelihood(ri_2)
    mean_squared_displacement_graph(ri_2)
    histogram(r_2_ns)
    error_maximmum_likelihood(ri_2)

def calculate_D(ri):
    total = 0
    est = 0
    for  i in range(len(ri)):
        total += ri[i] 
    est = total/(2*len(ri))
    return est/2


def maximum_likelihood(ri):
    k = 0
    vis = 10 ** (-3)
    d =  1.9 * 10 ** (-6)
    D = calculate_D(ri)
    # D = est # D = est/2t where t = 60s
    g = 6* math.pi * vis * d
    k = (g * D)/296.5
    print('the 1st value of k is: ', k)

# def error_maximmum_likelihood(ri):
#     delta_k = 0
#     vis = 10 ** (-3)
#     d = 

   

def mean_squared_displacement_graph(ri_2):
    x = len(ri_2) * [0]
    y = len(ri_2) * [0]
    y_a = len(ri_2) * [0]
    y_a[0] = ri_2[0]
    D = 0
    for i in range(len(ri_2)):
        x[i] = i/2
        if i >= 1:
            y_a[i] = ri_2[i] + y_a[i-1]
    m, b = np.polyfit(x, y_a, 1)
    for i in range(len(y_a)):
        y[i] = m*x[i] + b
    plt.scatter(x, y_a, color= "blue", marker= ".", s=10)
    D = abs(m/4)
    g = 6* math.pi * 0.93*10 ** (-3) * 0.95 * 10 ** (-6)
    k = (g * D)/296.5
    print('the 2nd value of k is: ',k)
    plt.plot(x, y)
    plt.xlabel('time (s)')
    plt.ylabel('mean squared displacement (m^2)')
    plt.title('Mean Squared Displacement as a Function of Time')
    plt.show()


def histogram(ri_2):
    
    D = 0
    bins = np.linspace(0,1.2e-6,24)
    y,x,_ = plt.hist(ri_2, bins,density = True, color = 'green',histtype = 'bar', rwidth = 1)
    x = (x[1:] + x[:-1])/2
    plt.xlabel('Step Size (m)')
    plt.ylabel('Number of occurences')
    plt.title('Histogram of step sizes')
    D,_ = curve_fit(func,x,y, p0 = 1.1662771482352923e-15)
    #print(D)
    plt.plot(x,func(x,D),'y-')
    plt.plot(x,func(x,2.29*10**-13),'r-')
    plt.plot(x,func(x,2.46*10**-13),'b-')
    g = 6* math.pi * 0.93*10 ** (-3) * 0.95 * 10 ** (-6)
    k = (g * D)/296.5
    print('the 3rd value of k is: ',k)
    plt.show()

'''
    x_values = np.linspace(0, np.max(ri_2), 100)
    plt.plot(x_values, rayleigh.pdf(x_values, scale=maximum_likelihood(ri_2)), 'r-', label='Rayleigh Fit')

    plt.xlabel('Step Size')
    plt.ylabel('Probability Density')
    plt.title('Step Size Distribution')
    plt.legend()
    plt.show()
'''

def func(x,D):
    return (x/D)* np.exp(-(x**2)/(2*D))

generate_ri(total_data)
