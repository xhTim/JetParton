#!/usr/bin/env python3

import argparse
import collections
import re
import math
import numpy as np
import os

playground = os.listdir("../Playground")
QAQBall = 0
for ifolder in range(50):
    print(ifolder)
    folderpath = os.listdir("../Playground/job-{}/results/QAQB".format(ifolder))
    for ifile in range(len(folderpath)):
        QAQBtemp = np.loadtxt("../Playground/job-{}/results/QAQB/{}".format(ifolder, folderpath[ifile]))
        QAQBall = QAQBtemp + QAQBall
np.savetxt("QAQBall", QAQBall)




