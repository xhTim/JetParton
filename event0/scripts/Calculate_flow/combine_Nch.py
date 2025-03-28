#!/usr/bin/env python3

import argparse
import collections
import re
import math
import numpy as np
import os

playground = os.listdir("./Playground")
QAQBall = 0
for ifolder in range(50):
    print(ifolder)
    #folderpath = os.listdir("../Playground/job-{}/results/QAQiB".format(ifolder))
    #for ifile in range(len(folderpath)):
    QAQBtemp = np.loadtxt("./Playground/job-{}/Nch_{}".format(ifolder, ifolder))
    QAQBall = QAQBtemp + QAQBall
np.savetxt("Nchall", QAQBall)

