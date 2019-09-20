# Moment-based analysis of interacting cell populations
This repository contains the code used for the paper 'Moment-based analysis of biochemical networks in a heterogeneous population of communicating cells' (https://arxiv.org/abs/1905.02053). All code and plots are provided in the form of Python notebooks (https://jupyter.org/). The provided code can be used to reproduce the case studies shown in the paper and is not meant to be a general tool. For further information, contact David Gonzales (gonzales@mpi-cbg.de) & Christoph Zechner (zechner@mpi-cbg.de).

# Installation
The provided code was written using Python v3.6.1 and uses the following libraries:
- Numpy v1.13.3 (https://www.numpy.org/)
- Scipy v1.1.0 (https://www.scipy.org/)
- Tellurium v2.1.0 (https://tellurium.readthedocs.io/en/latest/)
- Pandas v0.23.4 (https://pandas.pydata.org/)
- Matplotlib v2.0.2 (https://matplotlib.org/)

In case you do not have anything installed, the easiest way to start is to install using Anaconda. Instructions to install Python and Jupyter can be found in https://jupyter.readthedocs.io/en/latest/install.html. Numpy, Scipy, and Sympy are already included in the Anaconda installation. To install Tellurium:
```
pip install tellurium
```

# Usage
The python scripts (.py) contain functions used in the Python notebooks. To recreate the plots in the paper, run each of the notebooks for each case study. The current number of SSA samples are set to 200, but can be changed. We used 1000 stochastic simulations in the paper.