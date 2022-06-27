import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import root, newton
from matplotlib.widgets import Slider, Button

g = 9.81
def x(t, v, beta):
    return v*t*np.cos(beta)

def y(t, v, beta):
    return v*t*np.sin(beta) - g/2*t**2

def NST(v, beta):
    return 2*v*np.sin(beta) / g

init_beta = np.pi / 4
init_v = 100

t = np.linspace(0, NST(init_v, init_beta), 1000)
# t = np.linspace(0, 10, 1000)

fig, ax = plt.subplots()
line, = plt.plot(x(t, init_v, init_beta), y(t, init_v, init_beta), lw=2)

plt.subplots_adjust(left=0.2, bottom=0.25)

axangle = plt.axes([0.1, 0.25, 0.0225, 0.63])
angle_slider = Slider(
    ax=axangle,
    label='Angle',
    valmin=10**-30,
    valmax=np.pi / 2,
    valinit=init_beta,
    orientation="vertical"
)

axspeed = plt.axes([0.25, 0.12, 0.65, 0.03])
speed_slider = Slider(
    ax=axspeed,
    label='Speed',
    valmin=10**-16,
    valmax=110,
    valinit=init_v,
)

def update(val):
    t = np.linspace(0, NST(speed_slider.val, angle_slider.val), 1000)
    line.set_xdata(x(t, speed_slider.val, angle_slider.val))
    line.set_ydata(y(t, speed_slider.val, angle_slider.val))
    fig.canvas.draw_idle()

speed_slider.on_changed(update)
angle_slider.on_changed(update)

# Create a `matplotlib.widgets.Button` to reset the sliders to initial values.
resetax = plt.axes([0.05, 0.125, 0.1, 0.04])
button = Button(resetax, 'Reset', hovercolor='0.975')


def reset(event):
    speed_slider.reset()
    angle_slider.reset()
button.on_clicked(reset)

plt.show()
