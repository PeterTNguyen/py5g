#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt

class phy:
    def __init__(self):
        print('test')

    def ofdm_modulate(self):
        pass

    def gen_pss(self, n_cell_id):
        # TS 38.211, 7.4.2.2
        x = np.array([0,1,1,0,1,1,1]);
        for i in range(127 - len(x)):
            x = np.append(x, (x[i] + x[i+4])%2)
        m = 43 * (n_cell_id % 3);
        pss = 1 - 2*np.roll(x, -m);
        return pss;

    def gen_sss(self, n_cell_id):
        # TS 38.211, 7.4.2.3
        x0 = np.array([1, 0, 0, 0, 0, 0, 0]);
        x1 = np.array([1, 0, 0, 0, 0, 0, 0]);
        for i in range(127 - len(x1)):
            x0 = np.append(x0, (x0[i] + x0[i + 4]) % 2)
            x1 = np.append(x1, (x1[i] + x1[i + 1]) % 2)
        n1 = int(n_cell_id)/3;
        n2 = int(n_cell_id) % 3;
        m0 = 15 * (n1/112) + 5*n2;
        m1 = n1 % 112;
        s0 = 1 - 2*np.roll(x0, -m0);
        s1 = 1 - 2*np.roll(x1, -m1);
        sss = np.multiply(s0, s1);
        return sss

class phy_consts:
    def __init__(self):
        self.pss_ind = np.arange(127)+56;
        self.sss_ind = 2*240 + np.arange(127)+56;


class gnb:
    def __init__(self):
        self.DLNRB = 25;


phy_test = phy();
pss = phy_test.gen_pss(200);
print(pss);
sss = phy_test.gen_sss(200);
print(sss);
phy_const = phy_consts();
print(phy_const.pss_ind);
