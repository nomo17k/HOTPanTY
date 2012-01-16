#!/usr/bin/env python2.6
"""
Plot convolution kernel functions.
"""
import numpy as np
import matplotlib.pyplot as plt


x=np.arange(-10., +10., 0.01)


def g(x, s, i):
    y = np.exp(-x**2/2./s**2) * x**i
    return y / y.max()

factor = 1.8
factor = 2.0
#factor = 3.
print(3., 3./factor, 3./factor**2)

for i in range(7):
    plt.plot(x, g(x, 3./factor**2, i), 'b:')
for i in range(5):
    plt.plot(x, g(x, 3./factor, i), 'r--')
for i in range(3):
    plt.plot(x, g(x, 3., i), 'g-')
for i in range(1):
    plt.plot(x, g(x, 0.167, i), 'y-')

plt.show()
