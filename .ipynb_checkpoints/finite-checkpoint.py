import numpy as np
import matplotlib.pyplot as plt
#plt.style.use('classic')

V_0 = 40 # in [eV]
hc = 197.3 # in [eV nm]
mc2= 0.511e6 # in [eV]
a = 0.05 # in [nm]

    
def F(E: float) -> float:
    """
    Defining the function whose roots are to be found
    """
    l = np.sqrt((2 * mc2 / hc**2) * (E + V_0))
    kappa = np.sqrt((-2 * mc2 / hc**2) * E)
    
    return l * np.tan(l * a) - kappa


def G(v_0: float) -> float:
    E = -13.6
    l = np.sqrt((2 * mc2 / hc**2) * (E + v_0))
    kappa = np.sqrt((-2 * mc2 / hc**2) * E)
    
    return l * np.tan(l * a) - kappa


def plot_full(*energies:list):
    """
    Optionally accepts list of allowed energies for plotting wavefunctions
    """
    # Define the x values
    x = np.linspace(-2*a, 2*a, 1000)
    
    # Create an array for the potential energy
    V = np.zeros_like(x)

    # Set the potential energy inside the well region
    inside_well = (abs(x) <= a)
    V[inside_well] = -V_0

    plt.plot(x, V)
    plt.xlabel('Position (x)')
    plt.ylabel('Potential Energy (V)')
    plt.title('Finite Square Well Potential')
    plt.grid(True)
    plt.show()
    

def zero_crossings2(func, x_l, x_u):
    threshold = 1e-3
    if abs(x_u - x_l) < threshold:
        print('Root Constrained!')
        return(x_l, x_u)
    else:
        x_mid = (x_l + x_u) / 2
        if func(x_l) * func(x_mid) < 0:
            print('going left', x_l, x_u)
            return zero_crossings2(x_l, x_mid)
            
        elif func(x_mid) * func(x_u) < 0:
            print('going right', x_l, x_u)
            return zero_crossings2(x_mid, x_u)
        else:
            print('Root Not Found')
            return None
    
    

def zero_crossings1(x:np.numarray):
    crossings = []
    for i in range(len(x) - 1):
        prod = x[i] * x[i + 1]
        if prod < 0:
            crossings.append((i * V_0 / len(x), (x[i] + x[i+1]) / 2))
    
    print('Done!')
    return crossings


def wavefunction(x:float):
    """
    Define the wavefunction in the three potentials
    """
    if x > a:
        pass
    elif x < a:
        pass
    else:
        pass


def find_energy_root():
    e = np.linspace(-V_0*0.9, -10, 1000) # Eilon values go from 0 to V_0
    f = F(e)
    try:
        xLower, xUpper = zero_crossings2(F, -V_0+1, -10)
        print(xLower, xUpper)
    except TypeError:
        print('Rootfinding Aborted')
        pass
 
    fig, ax = plt.subplots(1, figsize=(8, 5))
    ax.axhline(y=0, color='k', linestyle='--')
    ax.plot(e, f, 'r-', label='stuff')
    #ax.set(xlim=(-40, 0), ylim=(-10, 10))
    fig.show()


def find_potential_root():
    v = np.linspace(13.6, 40, 1000)
    g = G(v)
    try:
        xLower, xUpper = zero_crossings2(G, 35, 20)
        print(xLower, xUpper)
    except TypeError:
        print('Rootfinding Aborted')
        pass

    fig, ax = plt.subplots(1, figsize=(8, 5))
    ax.axhline(y=0, color='k', linestyle='--')
    ax.plot(v, g, 'r-', label='stuff')
    #ax.set(xlim=(-40, 0), ylim=(-10, 10))
    fig.show()
    

def main():
    #plot_explicit()
    #plot_full()
    find_potential_root()
    
if __name__ == "__main__":
    main()
    
