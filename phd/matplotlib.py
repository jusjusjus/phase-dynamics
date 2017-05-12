
import numpy as np
from .core import period as p
import matplotlib.pyplot as plt

def xticks(*args, **kwargs):
    plt.xticks(p*np.arange(5)/4., (r'$0$', r'$\pi/2$', r'$\pi$', r'$3\pi/2$', r'$2\pi$'), *args, **kwargs)
    plt.xlim(0., p)

def yticks(*args, **kwargs):
    plt.yticks(p*np.arange(5)/4., (r'$0$', r'$\pi/2$', r'$\pi$', r'$3\pi/2$', r'$2\pi$'), *args, **kwargs)
    plt.ylim(0., p)

def xticks_centered(*args, **kwargs):
	plt.xticks(p*np.arange(5)/4.-0.5*p, (r'$-\pi$', r'$-\pi/2$', r'$0$', r'$\pi/2$', r'$\pi$'), *args, **kwargs)
	plt.xlim(-0.5*p, 0.5*p)

def yticks_centered(*args, **kwargs):
	plt.yticks(p*np.arange(5)/4.-0.5*p, (r'$-\pi$', r'$-\pi/2$', r'$0$', r'$\pi/2$', r'$\pi$'), *args, **kwargs)
	plt.ylim(-0.5*p, 0.5*p)
