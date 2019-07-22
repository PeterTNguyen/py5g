#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import math
import utils

class phy:
    def __init__(self):
        pass

    def ofdm_modulate(self, sym_grid, nr_info):
        [n_sc,n_syms] = np.shape(sym_grid);
        n_sc = nr_info.n_sc;
        n_fft = nr_info.n_fft;

        include_zsc = 1; # TODO: dc subcarrier parameters
        waveform = np.empty(0);

        # Modulate Symbols
        for sym in range(n_syms):

            # prepare IFFT, zeropad if necessary
            if(n_sc == n_fft):
                ifft_in = grid[:][sym];
            else:
                first_sc = (n_fft / 2) - (n_dl_rb * 6);
                ifft_in = np.zeros(n_fft);
                ifft_in[first_sc:first_sc+n_sc/2] = grid[0:n_sc/2,sym];
                mid_sc = first_sc + n_sc / 2 + (include_zsc == 0);
                ifft_in[mid_sc:mid_sc+n_sc/2] = grid[n_sc / 2:, sym];

            ifft_out = np.fft.ifft(np.fft.ifftshift(ifft_in));

            # Add Cyclic Prefix
            cp_len = nr_info.cp_lengths[sym%nr_info.syms_per_slot];
            cp_syms = ifft_out[(n_fft - cp_len):];
            ifft_out = np.append(cp_syms,ifft_out);
            waveform = np.append(waveform, ifft_out);

        return waveform

    def ofdm_demodulate(self, syms_in, nr_info):
        [n_samps, n_ants] = np.shape(syms_in);
        n_sc = nr_info.n_sc;
        n_fft = nr_info.n_fft;
        n_subframes = math.floor(n_samps/nr_info.samps_per_subframe);
        n_syms = n_subframes * nr_info.syms_per_subframe;
        sym_grid = np.zeros((n_sc, n_syms ));
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

class nr_numerology:
    def __init__(self, sc_spacing=15, n_dl_rb=30, e_cp=0):
        self.n_dl_rb = n_dl_rb;
        self.n_rb = n_dl_rb;
        self.sc_spacing = sc_spacing;
        self.e_cp = e_cp;
        self.n_sc = n_dl_rb * 12;
        self.n_fft = 2 ** (math.ceil(math.log(n_dl_rb * 12 / 0.85, 2.0)));
        self.n_fft = max(128, self.n_fft);
        self.samp_rate = self.n_fft * 15e3;
        if (e_cp):
            cp_lengths = 512 * np.ones(12);
        else:
            cp_lengths = 144 * np.ones(14);
            cp_lengths[0] = 160;
            cp_lengths[7] = 160;
        cp_lengths = cp_lengths * self.n_fft / 2048.0;
        self.cp_lengths = cp_lengths;
        self.syms_per_slot = len(cp_lengths);
        self.sym_lengths = cp_lengths + self.n_fft;
        scs_config = math.floor(math.log(sc_spacing/15.0, 2.0));
        scs_scale = 2**scs_config;

        self.syms_per_slot = len(self.cp_lengths);
        self.slots_per_subframe = math.floor(14/self.syms_per_slot)*scs_scale;
        self.syms_per_subframe = self.syms_per_slot * self.slots_per_subframe;

        self.samps_per_subframe = np.sum(cp_lengths) + self.n_fft * self.syms_per_subframe;
        self.subframe_period = self.samps_per_subframe / self.samp_rate;


phy_test = phy();

sss = phy_test.gen_sss(200);
phy_const = phy_consts();

n_dl_rb = 36;
n_sc = n_dl_rb * 12;
grid = np.zeros((n_dl_rb*12, 14));
filename = "data/nr_sync.bin"
wave = utils.read_data(filename);
pss_ind = phy_const.pss_ind;
sc_offset = n_sc/2 - 120;

pss_syms_ind = np.arange(4)+4;

plt.figure();

nr_num = nr_numerology(n_dl_rb = n_dl_rb, sc_spacing=30);

print(nr_num);

pss_corr = np.zeros(3);
pss_corr_ind = np.zeros(3);


for cell_id in range(3):
    pss = phy_test.gen_pss(cell_id);
    pss_grid = np.zeros((240, 4));
    pss_grid[pss_ind,0] = pss;

    grid[sc_offset:sc_offset+240, pss_syms_ind] = pss_grid;
    pss_t = phy_test.ofdm_modulate(grid, nr_num);
    pss_t = pss_t[np.nonzero(pss_t)];
    corr  = np.correlate(wave[:,0],pss_t)
    corr += np.correlate(wave[:, 1], pss_t)
    pss_corr_ind[cell_id] = np.argmax(np.abs(corr));
    pss_corr = np.max(np.abs(corr));
    print('ind: ' + str(np.argmax(np.abs(corr))))
    print('max: ' + str(np.max(np.abs(corr))));



    plt.plot(np.abs(corr));
pss_max_corr = np.argmax(pss_corr);
pss_offset = pss_corr_ind[pss_max_corr];

print(pss_offset);
plt.show()



