import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
import numpy as np
import statsmodels.api as sm


class PlotterSetup(object):
    def __init__(self):
        plt.rcParams["font.family"] = "Times New Roman"
        plt.rcParams["mathtext.fontset"] = "cm"
        self.font_size = 14
        self.tick_size = 14
        
        left = .45
        width = .5
        bottom = .58
        height = .5
        self.right = left + width
        self.top = bottom + height
        
        self.pop_colours = ["red", "darkorange", "gold",
                            "mediumturquoise", "dodgerblue",
                            "blue", "purple"
                            ]

    def tickers(self, ax):
        """
        Function to setup axis
        
        Args:
            ax (object):  Axes object to format
        Returns:
            ax (object):  Re-formatted axes
        """
        ax.yaxis.set_ticks_position('both')
        ax.xaxis.set_ticks_position('both')
        ax.xaxis.set_minor_locator(mtick.AutoMinorLocator())
        ax.yaxis.set_minor_locator(mtick.AutoMinorLocator())

        ax.tick_params(
            axis="y", which='both', direction="in", 
            labelsize=self.tick_size
        )
        ax.tick_params(
            axis="x", which='both', direction="in", 
            labelsize=self.tick_size
        )
        
        return ax

    def cdf_maker(self, data):
        """
        Function to make a CDF of the data
        
        Args:
            data (array):  The data to make the CDF of
        Returns:
            x (array):  The x-axis CDF data
            y (array):  The y-axis CDF data
        """
        x = np.sort(data)
        y = np.linspace(0, 1, len(x))
        return x, y
    
    def kde_maker(self, data):
        """
        Function to make a KDE of the data
        
        Args:
            data (array):  The data to make the KDE of
        Returns:
            x (array):  The x-axis KDE data
            y (array):  The y-axis KDE data
        """
        kde = sm.nonparametric.KDEUnivariate(data).fit()
        kde.density /= kde.density.max()
        return kde.support, kde.density