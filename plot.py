#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jul 21 04:03:13 2018

@author: katerina

PLOTS
"""

import matplotlib.pyplot as plt
import numpy as np

x, y = np.loadtxt('QW_NORM_qplate_n10_T1_steps2_niter-success5_only_NORM.txt', delimiter=',', unpack=True)

plt.plot(x,y,'g^')

plt.xlabel('x')
plt.ylabel('y')
plt.title('Schmidt norm after n=2 steps')
plt.legend()
plt.savefig('QW_NORM_qplate_n10_T1_steps2_niter-success5_only_NORM.png')