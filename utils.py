#!/usr/bin/env python
import numpy as np

filename = "data/nr5g_dl_float.bin"

def read_data(filename):
	with open(filename, 'rb') as f:
		data = np.fromfile(f, dtype=np.float32)
		real_syms = data[::2];
		imag_syms = data[1::2];

		num_syms = np.size(real_syms);
		syms = np.zeros(shape=(num_syms,1), dtype=complex);
		#for real, imag in zip(real_syms, imag_syms):
		for i in range(num_syms):
			syms[i] = real_syms[i] + imag_syms[i]*1j;

		syms_all = syms.reshape(2,num_syms/2).transpose();
	return syms_all;

read_data(filename);