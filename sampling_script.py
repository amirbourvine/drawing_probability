import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import norm
from time import time
from tqdm import tqdm
from scipy import stats

# assumption: sampling with replacement!
"""
This script is the base experiment setup. 
You can use:
    plot_distribution to plot the histogram
    plot_multiple_distributions to plot the progression of the normal curve
    plot_array_for_index to see the progression of stds/mus
"""

def create_graph(new_counts):
    """
    gets new_counts- new_counts[i] is the amount of time strand i was picked
    return x,y where:
    x is number of copies
    y is number of strands sampled x times
    """
    x = np.unique(new_counts)
    x = np.sort(x)

    # st = time()
    # y = np.zeros(x.shape[0])
    # for i in range(y.shape[0]):
    #     y[i] = np.sum(new_counts == x[i])
    # print(f"time regular {time()-st}")
    # st = time()

    unique_vals, counts = np.unique(new_counts, return_counts=True)
    y = np.zeros_like(x, dtype=int)  # Initialize result array

    # Match counts to corresponding values in x
    indices = np.searchsorted(unique_vals, x)
    mask = unique_vals[indices] == x  # Ensure only exact matches are counted
    y[mask] = counts[indices[mask]]

    # print(f"time otimized {time()-st}")
    # print(f"create_graph: {np.all(y_try==y)}")
    return x,y
    

def sample_strands(strand_counts, probabilities, M):
    """Take a sample of size M from the strands, based on their relative abundance."""

    sampled_strands = np.random.choice(len(strand_counts), size=M, p=probabilities)
    new_counts = np.bincount(sampled_strands, minlength=len(strand_counts))

    summed = strand_counts + new_counts

    return summed

def plot_distribution(x,y, iteration):
    """Plot the the requested graph"""
    plt.figure(figsize=(16, 10))

    plt.bar(x, y, color='blue', edgecolor='blue', alpha=0.7)
    plt.xlabel('Number of Copies')
    plt.ylabel('Number of Samples')
    plt.title(f'Samples per number of Copies - Iteration {iteration}')
    

    # Fit normal distribution to the data
    mu, std = stats.norm.fit(np.repeat(x, y))
    print(f"for iteration #{iteration}: mu = {mu} and std = {std}")
    
    # Generate points for the fitted normal distribution
    xmin, xmax = plt.xlim()
    x_fit = np.linspace(xmin, xmax, 100)
    y_fit = stats.norm.pdf(x_fit, mu, std) * sum(y)
    
    # Plot the fitted normal distribution
    plt.plot(x_fit, y_fit, 'r-', linewidth=2, label='Fitted Normal Distribution')


    plt.legend()
    plt.show()


def plot_multiple_distributions(mus, stds):
    """
    inputs are mus and stds and plots in the same plot all of the normal distributions
    """
    plt.figure(figsize=(16, 10))
    x = np.linspace(min(mus) - 3*max(stds), max(mus) + 3*max(stds), 1000)
    
    for i, (mu, std) in enumerate(zip(mus, stds)):
        y = norm.pdf(x, mu, std)
        plt.plot(x, y, label=f'Distribution {i}')
        
        # Add index number next to the peak of each curve
        peak_x = mu
        peak_y = norm.pdf(peak_x, mu, std)
        plt.text(peak_x, peak_y, str(i), fontsize=12, verticalalignment='bottom')
    
    # plt.xlabel('X-axis')
    # plt.ylabel('Probability Density')
    # plt.title('Multiple Normal Distributions')
    plt.legend()
    plt.grid(True)
    plt.show()


def plot_array_for_index(arr):
    """
    Plots the given array where the x-axis represents the index and the y-axis represents the array values.

    :param arr: List or NumPy array of numerical values to plot.
    """

    plt.figure(figsize=(8, 5))
    plt.plot(range(len(arr)), arr, marker='o', linestyle='-', color='b')
    plt.xlabel('iteration')
    plt.ylabel('std')
    # plt.title('Array Plot')
    plt.legend()
    plt.grid(True)
    plt.show()



def simulate_sampling(n, M, iterations):
    """Simulate the iterative sampling process."""
    
    probabs = np.ones(n) / n

    strand_counts = np.zeros(n)

    strand_counts = sample_strands(strand_counts, probabs, M)

    # print(new_counts.sum())
    
    mus = []
    stds = []

    for i in range(iterations):
        x,y = create_graph(strand_counts)
        plot_distribution(x,y, i)

        mu, std = stats.norm.fit(np.repeat(x, y))
        print(f"for iteration #{i}: mu = {mu} and std = {std}")

        mus.append(mu)
        stds.append(std)

        probabs = strand_counts / ((i+1)*M)

        if( ((i+1)*M) != strand_counts.sum()):
            print("Error")
            exit(1)
        
            
        strand_counts = sample_strands(strand_counts, probabs, M)

    # plot_multiple_distributions(mus, stds)

    # plot_array_for_index(stds)
    
    
if __name__ == "__main__":
    # Parameters
    n = int(1e4)      # Number of unique strands
    M = int(1e7)        # Sample size
    # print(M)
    iterations = 10    # Number of iterations

    simulate_sampling(n, M, iterations)
    