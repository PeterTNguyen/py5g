#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt

class phy:
    def __init__(self):
        self.dpss = [1, -1, -1, 1, -1, -1, -1, -1, 1, 1, -1, -1, -1, 1, 1, -1, 1, -1, 1, -1, -1, 1, 1, -1, -1,
                     1, 1, 1, 1, 1, -1, -1, 1, -1, -1, 1, -1, 1, -1, -1, -1, 1, -1, 1, 1, 1, -1, -1, 1, 1,
                     -1, 1, 1, 1, -1, 1, 1, 1, 1, 1, 1, -1, 1, 1, -1, 1, 1, -1, -1, 1, -1, 1, 1, -1, -1,
                     -1, -1, 1, -1, -1, -1, 1, 1, 1, 1, -1, -1, -1, -1, -1, -1, -1, 1, 1, 1, -1, -1, -1, 1, -1,
                     -1, 1, 1, 1, -1, 1, -1, 1, 1, -1, 1, -1, -1, -1, -1, -1, 1, -1, 1, -1, 1, -1, 1, 1, 1, 1, -1];
        print('test')

    def ofdm_modulate(self):
        pass

    def gen_pss(self, n_cell_id):
        n2_off = 43 * (n_cell_id % 3)
        pss = self.dpss[n2_off:] + self.dpss[:n2_off];
        return pss;

class phy_consts:
    def __init__(self):
        self.pss_ind = 1;


class gnb:
    def __init__(self):
        self.DLNRB = 25;


phy_test = phy();
pss = phy_test.gen_pss(200);
print(pss);
phy_const = phy_consts();
