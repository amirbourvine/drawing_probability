import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import norm
from time import time
from tqdm import tqdm
from scipy import stats

from sampling_script import create_graph,sample_strands

# assumption: sampling with replacement!
"""
This script is is for plotting a graph of the slope of the linear growth of the std to the M/n ratio.
"""


def plot_xy(x,y):
    """
    Plots the given array where the x-axis represents the index and the y-axis represents the array values.
    
    :param arr: List or NumPy array of numerical values to plot.
    """

    plt.figure(figsize=(8, 5))
    plt.plot(x, y, marker='o', linestyle='-', color='b')
    plt.xlabel('M/n ratio')
    plt.ylabel('slope')
    plt.legend()
    plt.grid(True)
    plt.show()



def simulate_sampling(ns, Ms, iterations):
    """Simulate the iterative sampling process."""
    
    ratios = []
    slopes = []

    for n,M in tqdm(zip(ns, Ms), total=len(ns)):
        probabs = np.ones(n) / n

        strand_counts = np.zeros(n)

        strand_counts = sample_strands(strand_counts, probabs, M)
    
        stds = []

        for i in range(iterations):
            x,y = create_graph(strand_counts)
            _, std = stats.norm.fit(np.repeat(x, y))

            stds.append(std)

            probabs = strand_counts / ((i+1)*M)

            if( ((i+1)*M) != strand_counts.sum()):
                print("Error")
                exit(1)
            
            strand_counts = sample_strands(strand_counts, probabs, M)
        
        ratios.append(float(M)/n)
        slopes.append(np.mean(np.abs(np.diff(stds))))

    
    sorted_indices = np.argsort(ratios)  # Get sorting indices
    sorted_ratios = np.array(ratios)[sorted_indices]  # Sort ratios
    sorted_slopes = np.array(slopes)[sorted_indices]  # Reorder slopes
    plot_xy(sorted_ratios,sorted_slopes)

    np.savez("slope_to_ratio", x=sorted_ratios, y=sorted_slopes)

    
 
    
    
if __name__ == "__main__":
    n = int(1e3)
    iterations = 10

    Ms = list(range(50 * n, 1000 * n + 1, 100 * n))
    ns = [n]*len(Ms)

    simulate_sampling(ns, Ms, iterations)
    