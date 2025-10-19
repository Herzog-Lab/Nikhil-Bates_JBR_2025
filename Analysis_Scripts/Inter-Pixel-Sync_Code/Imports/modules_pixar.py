# -*- coding: utf-8 -*-
"""
Created on Fri Oct 29 18:16:02 2021

@author: Nikhil
"""

import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
import matplotlib as mpl
from matplotlib import ticker
import matplotlib.gridspec as gridspec
import matplotlib.patches as mpatches
from mpl_toolkits import mplot3d
from matplotlib.ticker import PercentFormatter, MultipleLocator, AutoMinorLocator
import numpy as np
import pandas as pd
from pyboat import WAnalyzer
from skimage import io
from skimage.transform import rescale, resize
from skimage.filters import gaussian
from skimage import measure
import statistics as stat
from scipy import ndimage
from scipy.stats import pearsonr
from astropy.timeseries import LombScargle
import seaborn as sns
import random
import os
import sys
import shutil
import warnings
from datetime import datetime
import platform

if platform.system() == 'Windows':
    from rpy2 import robjects as ro
elif platform.system() == 'Darwin':
    import subprocess
else:
    sys.exit("The code supports only windows and osx")
