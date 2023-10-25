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


def plot_full(energies:dict, scale=1.0):
    """
    Optionally accepts list of allowed energies for plotting wavefunctions
    
    parameters
    ----------
    energies : dict
    dictionary with keys 'even',and 'odd' corresponding
    corresponding to the parityy of the wave function solutions.

    scale : float
    scale factor for the wave function for better visualisation.
    Cannot,as yet, confirm the appropirateness
    given the normalisation constraints.
    """

    # Define extent of potential
    l = 13
    # Define the x values
    x = np.linspace(-l*a, l*a, 1000)
    
    # Create an array for the potential energy
    V = np.zeros_like(x)

    # Set the potential energy inside the well region
    inside_well = (abs(x) <= a)
    V[inside_well] = -V_0

    # Wave Functions
    # psi_0 = wave_function(x, energies[0])

    fig, ax = plt.subplots(1, figsize=(8, 5))

    # Wave functions
    for parity, Energy in energies.items():
        ax.plot(x, scale * wave_function(x, Energy[0], parity) + Energy[0], label=f'Energy: {Energy[0]:.2f} eV')
        
    ax.plot(x, V, 'k')
    [ax.axhline(y=i, alpha=0.5, color='r', linestyle='--') for i in energies.values()]
    ax.set(xlabel='x [nm]', 
           ylabel='Potential Energy [eV]', 
           title='Finite Square Well Potential')
    ax.grid(True)
    ax.legend()
    #plt.savefig('wave_function.png', dpi=600)
    plt.show()
    

def wave_function(x, E:float, parity:str) -> float:
    """
    Takes an array of positions and energy for a bound state.
    Returns an array of wave function values for those positions.
    checks parity of bound state and adjusts output accordingly
    """
    l = np.sqrt((2 * mc2 / hc**2) * (E + V_0))
    kappa = np.sqrt((-2 * mc2 / hc**2) * E)
    A = 1 / np.sqrt(a + 1 / kappa) 

    # Region IIII
    mask_IIII = x > a
    # Region III
    mask_III = np.logical_and(x >= 0, x <= a)
    # Region II (By Symmetry)
    mask_II = np.logical_and(x >= -a, x <= 0)
    # Region I (By Symmetry)
    mask_I = x < -a
    
    if parity == 'even':
        return np.piecewise(x=x, 
            condlist=[mask_IIII, mask_III, mask_II, mask_I], 
            funclist=[lambda x: A * np.cos(l * a) * np.exp(-kappa * (x - a)),
                        lambda x: A * np.cos(l * x), 
                        lambda x: A * np.cos(l * -x), 
                        lambda x: A * np.cos(l * a) * np.exp(-kappa * (-x - a))])
    elif parity == 'odd':
        return np.piecewise(x=x, 
            condlist=[mask_IIII, mask_III, mask_II, mask_I], 
            funclist=[lambda x: A * np.sin(l * a) * np.exp(-kappa * (x - a)),
                        lambda x: A * np.sin(l * x), 
                        lambda x: -A * np.sin(l * -x), 
                        lambda x: -A * np.sin(l * a) * np.exp(-kappa * (-x - a))])
    else:
        print("Neither???")
        return None


def zero_crossings2(func, x_l, x_u) -> tuple:
    threshold = 1e-3
    if abs(x_u - x_l) < threshold:
        print('Root Constrained!')
        return(x_l, x_u)
    else:
        x_mid = (x_l + x_u) / 2
        if func(x_l) * func(x_mid) < 0:
            print('going left', x_l, x_u)
            return zero_crossings2(func, x_l, x_mid)
            
        elif func(x_mid) * func(x_u) < 0:
            print('going right', x_l, x_u)
            return zero_crossings2(func, x_mid, x_u)
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


def find_energy_root():
    e = np.linspace(-39, -5, 1000) # Energy values go from 0 to V_0
    f = F(e)
    try:
        xLower, xUpper = zero_crossings2(F, -V_0+1, -10)
        midpoint = (xLower + xUpper) / 2
        dist = abs(xLower - xUpper)
        print(f"Estimated root: {midpoint}")
        
    except TypeError:
        print('Rootfinding Aborted')
        pass
 
    fig, ax = plt.subplots(1, figsize=(8, 5))
    ax.axhline(y=0, color='k', linestyle='--')
    ax.plot(e, f, 'r-', label='stuff')
    #ax.set(xlim=(-40, 0), ylim=(-10, 10))
    plt.show()


def find_potential_root():
    v = np.linspace(13.6, 40, 1000)
    g = G(v)
    try:
        xLower, xUpper = zero_crossings2(G, 35, 20)
        midpoint = (xLower + xUpper) / 2
        print(f"Estimated root: {midpoint}")
    except TypeError:
        print('Rootfinding Aborted')
        pass

    fig, ax = plt.subplots(1, figsize=(8, 5))
    ax.axhline(y=0, color='k', linestyle='--')
    ax.plot(v, g, 'r-', label='stuff')
    #ax.set(xlim=(-40, 0), ylim=(-10, 10))
    plt.show()
    

def main():
    #plot_explicit()
    plot_full(energies={'even': [-26.34322554814030], 
                        'odd': [-0.090842234843586]
                        },
              scale=5)
    # find_energy_root()
    # find_potential_root()
    # print(f"{hc**2 / (2 * mc2) * (np.pi / (2 * a))**2}")
if __name__ == "__main__":
    main()
    
