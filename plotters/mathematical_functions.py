import matplotlib.patheffects as pe
import numpy as np
from scipy.special import jv 

from amuse.lab import constants, units


def GW_interferometer_semi_ecc(ax, m1, m2):
    """
    Over plot the LISA and muAres frequency range in a vs. (1-e) plots
    
    Args:
        ax (object):  Axes object to plot on
        m1 (float):  Mass of the primary object
        m2 (float):  Mass of the secondary object
    Returns:
        ax (object):  Axes object with the plot
    """
    ecc_range = np.linspace(0.0001, (1.-1e-8), 50)

    max_Ares = max_LISA = 1 | units.Hz
    min_Ares = 1e-7 | units.Hz
    min_LISA = 1e-5 | units.Hz
    peak_Ares = 1e-3 | units.Hz
    peak_LISA = 1e-2 | units.Hz
    
    Ares_semimaj_max = gw_cfreq_semi(ecc_range[1:], max_Ares, m1, m2)
    Ares_semimaj_min = gw_cfreq_semi(ecc_range[1:], min_Ares, m1, m2)
    Ares_semimaj_peak = gw_cfreq_semi(ecc_range[1:], peak_Ares, m1, m2)
    Ares_data = [Ares_semimaj_min, Ares_semimaj_peak, Ares_semimaj_max]

    LISA_semimaj_max = gw_cfreq_semi(ecc_range[1:], max_LISA, m1, m2)
    LISA_semimaj_min = gw_cfreq_semi(ecc_range[1:], min_LISA, m1, m2)
    LISA_semimaj_peak = gw_cfreq_semi(ecc_range[1:], peak_LISA, m1, m2)
    LISA_data = [LISA_semimaj_min, LISA_semimaj_peak, LISA_semimaj_max]

    ecc_range = [np.log10(1.-i) for i in ecc_range[1:]]
    text_angle = -45
    
    ax.text(1.6, -5.6, r"LISA ($f_{\rm peak} = 10^{-2}$ Hz)", 
            path_effects=[pe.withStroke(linewidth=1, foreground="white")],
            color='red', fontsize=10, rotation=text_angle)
    ax.text(3.5, -6.1, r"$\mu$Ares ($f_{\rm peak} = 10^{-3}$ Hz)", 
            path_effects=[pe.withStroke(linewidth=1, foreground="white")],
            color='blue', fontsize=10, rotation=text_angle)

    for data, c in zip([Ares_data, LISA_data], ["blue", "red"]):
        ax.plot(data[0], ecc_range, linestyle=':', color=c, zorder=2)
        ax.plot(data[1], ecc_range, linestyle='-.', color=c, zorder=3)
        ax.plot(data[2], ecc_range, linestyle=':', color=c, zorder=4)
        ax.fill_between(np.append(data[0], data[2]), 
                        np.append(ecc_range[:], ecc_range[::-1]), 
                        alpha=0.25, color=c, zorder=1
                        )

    return ax

def gfunc(ecc):
    """
    Compute the g(e) function for the gravitational wave strain.
    Uses eq. A1 of Peters & Matthews 1963.
    
    Args:
        ecc (float):  Eccentricity of the binary system
    Returns:
        g (float):  The g(e) function"""
    nharm = gw_harmonic_mode(ecc)
    return nharm**4/32.*((jv(nharm-2., nharm*ecc)-2.*ecc*jv(nharm-1., nharm*ecc) + 2./nharm*jv(nharm, nharm*ecc) \
            + 2.*ecc*jv(nharm+1., nharm*ecc) - jv(nharm+2., nharm*ecc))**2 + (1.-ecc**2)*(jv(nharm-2., nharm*ecc) \
            - 2.*jv(nharm, nharm*ecc) + jv(nharm+2., nharm*ecc))**2 + 4./(3.*nharm**2)*(jv(nharm, nharm*ecc)**2))

def gw_cfreq_semi(ecc_arr, freq_val, m1, m2):
    """
    Function to get constant frequency curves (Eqn. 43 of Samsing et al. 2014).
    LISA (1e-2 Hz) peak sensitivity with range 1e-5<f<1 [https://lisa.nasa.gov/]
    muAres (1e-3 Hz) peak sensitivity with range 1e-7<f<1 [arXiv:1908.11391]
    
    Args:
        ecc_arr (array):  The eccentricity array
        freq_val (float):  The constant frequency wishing to plot
        m1/m2 (float):  The binary mass
    Returns:
        semi_maj (array):  The semi-major axes of the system
    """
    term1 = np.sqrt(constants.G*(m1+m2))/np.pi
    semi_maj = [np.log10(((term1*(1.+i)**1.1954/(1.-i**2)**1.5*freq_val**-1)**(2./3.)).value_in(units.au)) for i in ecc_arr]
    return semi_maj

def gw_freq(semi, ecc, m1, m2):
    """
    Frequency equation is based on Samsing et al. 2014 eqn (43). 
    
    Args:
        semi (float):  The semi-major axes of the system
        ecc (float):  The eccentricity of the binary system
        m1/m2 (float):  The binary mass
    Returns:
        freq (float):  The gravitational wave
    """
    nharm = gw_harmonic_mode(ecc)
    freq =  (2.*np.pi)**-1*np.sqrt(constants.G*(m1+m2)/abs(semi)**3)*nharm
    return freq

def gw_harmonic_mode(ecc):
    """
    Finding the peak harmonic of gravitational frequency for a given eccentric orbit.
    Equation 36 of Wen (2003)
    
    Args:
        ecc (float):  The eccentricity of the binary system
    Returns:
        nharm (float):  The peak harmonic mode
    """ 
    nharm = 2.*(1.+ecc)**1.1954/(1.-ecc**2)**1.5
    return nharm

def gw_dfreq(ecc, m1, m2, semi, chirp_mass, ecc_func):
    """
    Function to take into account the limited LISA observation time, T ~ 5yrs
    Based on equation (6) of Kremer et al. 2019.
    The redshift is calculated from Ned Wright's calculator with data from arXiv:1807.06209
    (omegaM = 0.315, omegaL = 0.685, H0 = 67.4 km/s/Mpc)
    
    Args:
        ecc (float):  The eccentricity of the binary system
        m1/m2 (float):  The binary mass
        semi (float):  The semi-major axes of the system
        chirp_mass (float):  The chirp mass of the binary system
        ecc_func (float):  The eccentricity function
    Returns:
        dfreq (float):  The change in frequency over some obs. time
    """
    redshift = 0.116 
    nharm = gw_harmonic_mode(ecc)
    forb = np.sqrt(constants.G*(m1+m2))/(2*np.pi)*abs(semi)**-1.5*(redshift+1.)**-1
    dfreq = (96.*nharm)/(10.*np.pi)*(constants.G*chirp_mass)**(5./3.)/(constants.c**5)*(2.*np.pi*forb)**(11./3.)*abs(ecc_func)
    return dfreq

def gw_strain(semi, ecc, m1, m2):
    """
    Calculate GW strain using eqn (7) Kremer et al. 2018 and eqn (20) of Peters and Matthews (1963).
    
    Args:
        semi (float):  The semi-major axes of the system
        ecc (float):  The eccentricity of the binary system
        m1/m2 (float):  The individual mass components of the binary system
    Returns:
        strain (float):  The gravitational wave strain
    """

    dist = 0.5 | units.Gpc
    redshift = 0.116  # arXiv:1807.06209 --> OmegaM=0.315, OmegaL=0.685, H0=67.4 km/s/Mpc
    ecc = abs(ecc)

    chirp_mass = (m1*m2)**0.6/(m1+m2)**0.2*(1.+redshift)**-1
    cfactor = 2./(3.*np.pi**(4./3.))*(constants.G**(5./3.))/(constants.c**3)*(dist*(1.+redshift))**-2.
    ecc_func = (1.+(73./24.)*ecc**2+(37./96.)*ecc**4.)*(1.-ecc**2.)**-3.5

    nharm = gw_harmonic_mode(ecc)
    freq = gw_freq(semi, ecc, m1, m2)
    dfreq = gw_dfreq(ecc, m1, m2, semi, chirp_mass, ecc_func)
    factor = dfreq*(5. | units.yr)/freq

    strain = min(1., factor)*cfactor*chirp_mass**(5./3.)*freq**(-1./3.)*(2./nharm)**(2./3.)*(gfunc(ecc)/ecc_func)
    return (strain.value_in(units.s**-1.6653345369377348e-16))**0.5
    
def tGW_peters(m1, m2, semi, ecc):
    """
    Calculate the gravitational wave timescale using the Peters formula
    
    Args:
        m1 (float):  Mass of the primary object
        m2 (float):  Mass of the secondary object
        semi (float):  Semi-major axis of the binary
        ecc (float):  Eccentricity of the binary
    Returns:
        tgw (float):  Gravitational wave inspiral time
    """
    reduced_mass = (m1*m2)/(m1+m2)
    total_mass = m1 + m2
    tgw = (5./256.)*(constants.c)**5./(constants.G**3.)*(semi**4.*(1.-ecc**2.)**3.5)/(reduced_mass*total_mass**2.)
    return tgw

def peters_orb_param_evolution(semi_major, eccentricity, mass1, mass2, dt, sim_len):
    """
    Track the predicted orbital parameter evolution of a binary system
    Uses equation 5.6 and 5.7 of Peters (1964) to calculate the semi-major axis 
    and eccentricity of the system.
    
    Args:  
        semi_major (float):  The initial semi-major axis of the binary system
        eccentricity (float):  The initial eccentricity of the binary system
        mass1 (float):  The mass of the primary object
        mass2 (float):  The mass of the secondary object
        dt (float):  The simulation timestep
        sim_len (int): The number of iterations in simulation
    Returns:
        predicted_sem_evol (array):  The predicted semi-major axis evolution
        predicted_ecc_evol (array):  The predicted eccentricity evolution
    """
    system_mass = mass1 + mass2
    predicted_sem_evol = [ ]
    predicted_ecc_evol = [ ]
    
    rperi = semi_major*(1-eccentricity)
    iteration = 0
    while rperi > (6.*constants.G*system_mass)/(constants.c**2.) \
        and iteration < sim_len:
            if not (predicted_sem_evol and predicted_ecc_evol):
                predicted_sem_evol.append(semi_major.value_in(units.au))
                predicted_ecc_evol.append(eccentricity)
                
            else:
                coefficient = (constants.G**3.*mass1*mass2*system_mass)/(constants.c**5.)
                ecc_func = (1.+(73./24.)*eccentricity**2.+(37./96.)*eccentricity**4.)
                
                da = -(64./5.)*coefficient/(semi_major**3.*(1.-eccentricity**2.)**(7./2.))*ecc_func
                de = -(304./15.)*eccentricity*coefficient/(semi_major**4.*(1.-eccentricity**2.)**(5./2.))*(1.+(121./304.)*eccentricity**2.)
                
                semi_major += da*dt
                eccentricity += de*dt
                rperi = semi_major*(1.-eccentricity)
                
                predicted_sem_evol.append(abs(semi_major).value_in(units.au))
                predicted_ecc_evol.append(eccentricity)
        
            iteration += 1
            
    return predicted_sem_evol, predicted_ecc_evol

def power_law_fit(x, slope, alpha, yint):
    """
    Function to fit a power law to the data
    
    Args:
        x (array):  The x-axis data
        slope (float):  The slope of the power law
        alpha (float):  The power of the power law
        yint (float):  The y-intercept of the power law
    Returns:
        y (array):  The y-axis line-of-best-fit data
    """
    return slope*x**alpha+yint