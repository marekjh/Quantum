import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as ani

def squared_modulus(c):
    return np.real(c)**2 + np.imag(c)**2

a = np.sqrt(1/8)
V0 = 10
L = 5
j = np.arange(-30, 31, 1)
dv = 6 * a / (len(j) - 1)

v0 = 2.584
v = v0 + j * dv
k = np.sqrt(v**2 + V0)
R = np.exp(-2j*L*v) * ((1j*v*np.tan(k*L) + k) / (1j*v*np.tan(k*L) - k))
A = -2j*v*np.exp(-1j*L*v) / (k*np.cos(k*L) - 1j*v*np.sin(k*L))
G = np.exp(-(v - v0)**2 / (2*a**2)) # G for Gaussian

x_in = np.linspace(0, L, L*10)
x_out = np.linspace(L, 35, (35 - L)*10)
t = np.linspace(-2.8, 2, 24)

def psi_out_squared(t):
    def psi_out(x, t):
        return np.sum(dv*G*(np.exp(-1j*v*x) + R*np.exp(1j*v*x))*np.exp(-1j*v**2*t))

    return np.array([squared_modulus(psi_out(x, t)) for x in x_out])

def psi_in_squared(t):
    def psi_in(x, t):
        return np.sum(dv*G*A*np.sin(k*x)*np.exp(-1j*v**2*t))

    return np.array([squared_modulus(psi_in(x, t)) for x in x_in])

fig = plt.figure()
def draw_plot(tt):
    plt.clf()
    plt.plot(x_in, psi_in_squared(tt), color="blue")
    plt.plot(x_out, psi_out_squared(tt), color="blue")
    plt.title(fr"$\nu_0={v0}$, $t={round(tt, 3)}$") 
    plt.xlabel(r"$x$")
    plt.ylabel(r"$|\Psi(x,t)|^2$")
    plt.ylim(0, 2)
animator = ani.FuncAnimation(fig, draw_plot, t, interval=40)
plt.show()
