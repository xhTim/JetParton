#!/usr/bin/env python3

import argparse
import collections
import re
import math
import numpy as np
import os
import sys

playground = os.listdir("../../../Playground")
QAQBall = np.zeros((11))
Neventall = np.zeros((11))
user_input = int(sys.argv[1])
for ifolder in range(user_input, user_input + 1):
    print(ifolder)
    folderpath = os.listdir("../../../Playground/job-{}/results/Nch".format(ifolder))
    for ifile in range(len(folderpath)):
        QAQBtemp = np.loadtxt("../../../Playground/job-{}/results/Nch/{}".format(ifolder, folderpath[ifile]))
        for ii in range(len(QAQBtemp)):
            nchb = int(QAQBtemp[ii]/10)
            if (nchb >10):
                nchb = 10
            QAQBall[nchb] = QAQBtemp[ii] + QAQBall[nchb]
            Neventall[nchb] = Neventall[nchb] +1
aaa = np.array([Neventall, QAQBall])
np.savetxt("Nch_{}".format(user_input), aaa.T)


