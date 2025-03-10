# Analysis of the Drawing Probability

Assume there are n strands, each with billions of copies. 
A sample of size M is taken from these n strands so each draw is chosen uniformly at random. 
The plot that describes the number of strands (y axis) with a given number of copies (x axis) should have a normal distribution.
Then, another sample is taken according to the new distribution and this process repeats several time. 
The goal is to analyze how this plot changes with the number of iterations.


### Installing

Follow the following steps in order to run our project:

Clone the git to your local machine:
```
git clone https://github.com/amirbourvine/drawing_probability.git
```

Install the requirements:
```
pip install -r requirements.txt
```

Run any of the scripts. 
Documentation could be found in the beginning of each scrtip.



### Details about the Scripts

* sampling_script.py: This script is the base experiment setup. 
You can use:
    plot_distribution to plot the histogram
    plot_multiple_distributions to plot the progression of the normal curve
    plot_array_for_index to see the progression of stds/mus

* exp_vs_analytical.py: This script is analytical vs experiments.
Plots the analytically calculated std values with the stds observed in the experiment.

* slope_for_ration.py: This script is is for plotting a graph of the slope of the linear growth of the std to the M/n ratio.

* curve_fitting.py: This script is fitting various functions to the slope for ratio curve.
first run slope_for_ration.py to create slope_to_ratio.npz, that this script uses.


