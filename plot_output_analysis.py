import numpy as np
import matplotlib.pyplot as plt
from scipy import signal
import csv
from scipy.signal import find_peaks
fs = 19.23 # Sampling frequency
# Generate the time vector properly

input_signal = []
output_signal= []
freq_scale   = []
mag_spec_input=[]
phase_spec_input = []
impulse_resp = []
freq_resp_mag = []
freq_resp_phase = []
mag_spec_out = []
mag_phase_out = []

with open('outputsAnalysis.csv', 'rU') as csvfile:
	spamreader = csv.reader(csvfile, delimiter=',', quotechar='|')
	for row in spamreader:
		try:
			#print row
			input_signal.append(float(row[0]))
			output_signal.append(float(row[1]))
			freq_scale.append(float(row[2]))
			mag_spec_input.append(float(row[3]))
			phase_spec_input.append(float(row[4]))
			impulse_resp.append(float(row[5]))
			freq_resp_mag.append(float(row[6]))
			freq_resp_phase.append(float(row[7]))
			mag_spec_out.append(float(row[8]))
			mag_phase_out.append(float(row[9]))
		except:
			pass
			
			
			
			
			


plt.plot(freq_scale, mag_spec_out, label='output frequency spectrum')
plt.plot(freq_scale, mag_spec_input, label='input frequency spectrum')



#plt.plot(freq_scale, input_signal, label='input signal')
#plt.plot(freq_scale, output_signal, label='output signal')

#plt.plot(freq_scale, impulse_resp, label='impulse response')
plt.ylabel('Magnitude')
plt.xlabel('Frequency (Hz)')
plt.legend()
plt.show()
