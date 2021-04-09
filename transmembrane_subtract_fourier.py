
# plt.figure()
# plt.plot(t, host_cell.step_response / np.linalg.norm(host_cell.step_response))
# plt.plot(t, virus.step_response / np.linalg.norm(virus.step_response))
# plt.show()
#
# n = 1000000
# host_cell_fft = np.fft.fft(host_cell.step_response / np.linalg.norm(host_cell.step_response),n)
# virus_fft = np.fft.fft(virus.step_response / np.linalg.norm(virus.step_response),n)
# freq = np.fft.fftfreq(n,dt)
#
# plt.figure()
# plt.plot(freq, virus_fft-host_cell_fft)
# plt.show()
# plt.figure()
#
# shaped_pulse = np.fft.ifft(virus_fft-host_cell_fft)
# shaped_dt = dt*(len(t) / len(shaped_pulse)) # to get high freq res, need a large fft n; but this messes up the timestep
# # shaped_pulse = shaped_pulse[0:len(shaped_pulse) //3]
# shaped_pulse -= np.min(np.real(shaped_pulse)) #add a DC offset
# shaped_pulse /= np.max(shaped_pulse) #normalize
# # shaped_pulse = np.pad(shaped_pulse,[0,500000])
# # shaped_pulse = np.tile(shaped_pulse, 20)
#
# shaped_times = np.arange(len(shaped_pulse)) * shaped_dt
#
# plt.plot(shaped_times,np.real(shaped_pulse))
# plt.show()
#
# plt.figure()
# plt.plot(shaped_times,convolve_output(shaped_pulse, host_cell, shaped_dt))
# plt.plot(shaped_times,convolve_output(shaped_pulse, virus, shaped_dt))
# plt.show()
#
