import numpy as np

data = np.loadtxt("out.dat", delimiter=",").T

import matplotlib.pyplot as plt


# Create plotting routine
def diff_plot(t, y, fig_name=None):
    
    # Plot
    fig, ax = plt.subplots()
    ax.plot(t, y[0], 'r', label='$y_{0}$')
    ax.plot(t, y[1], 'b', label='$y_{1}$')
    ax.set(xlabel='$t$', ylabel='$y$')
    ax.legend(loc='best')
    
    # Show figure
    plt.show()
    
    # Save figure
    if fig_name is not None:
        fig.savefig(f'{fig_name}.pdf')
      
      
diff_plot(data[0, :], data[1:, :])

print('C Shape:', data.shape)

def func(t, y):
    return [(1. - 0.01 * y[1]) * y[0], (0.02 * y[0] - 1.) * y[1]]
t_span = (0.0, 500.0)
y0 = (20.0, 20.0)
rtol = 1.0e-7
atol = 1.0e-8

from scipy.integrate import solve_ivp

solution = solve_ivp(func, t_span, y0, method='RK45', rtol=rtol, atol=atol)

print('Solve IVP Shape', solution.y.shape)
diff_plot(solution.t, solution.y)

# breakpoint()