import numpy as np
import matplotlib.pyplot as plt

SAMPLE_RATE = 1000
def damped_harmonic_oscillator(a, v, t, w, theta):
    '''
    a = amplitude > 0 (constant)
    v = strength of damping (constant)
    t = time
    w = angular frequency > 0 (constant)
    theta = phase angle (constant)
        Starting condition, where does the signal start in degrees relative to the peak?
        0 is trough, 180 is peak
    '''

    return a * np.exp(-v * t / 2) * np.cos(w * t - theta)

def add_waves(t1, w1, t2, w2):
    t1_end_ind = int(t1[-1] * SAMPLE_RATE)
    t2_start_ind = int(t2[0] * SAMPLE_RATE)
    return np.concatenate([w1[:t2_start_ind], w1[t2_start_ind:] + w2[:t1_end_ind-t2_start_ind]])

def phase_to_time_given_period(phase, period):
    ratio = phase / np.pi
    return ratio * period

amplitude = 1
time = np.linspace(0, 1, SAMPLE_RATE)
# convert from Hz to angular freq
period = .032
hertz = 1/period
w = 2 * np.pi * hertz
theta = np.deg2rad(-90)
# critical damping * scalar
v = 2 * w * 0.10
# v = 0

plt.subplot(2, 1, 1)
offset = phase_to_time_given_period(np.pi, period)
w1 = damped_harmonic_oscillator(amplitude, v, time, w, theta)
w2 = -1 * damped_harmonic_oscillator(amplitude, v, time, w, theta)
w1_p_w2 = add_waves(time, w1, time+offset, w2)
plt.plot(time, w1, label='w1')
plt.plot(time+offset, w2, label='w2')
plt.plot(time, w1_p_w2, label='w1+w2_{}'.format(np.pi))
plt.xlim([0, 0.3])
plt.legend()
# plt.show()
print('amplitude to cancel at t: a{} t{}'.format(w1_p_w2[int(SAMPLE_RATE*period*1.5)], period*1.5))


# create two more waves at some time period
t_period = int(period*SAMPLE_RATE)
print('orig error', np.sum(np.abs(w1_p_w2[t_period:])), 'total', np.sum(np.abs(w1_p_w2[t_period:])))
# for a_ind, amplitude in enumerate(np.linspace(0.01, 0.06, 10)):
plt.subplot(2, 1, 2)
best = [0, 0, 0, np.inf]
# for a_ind, amplitude in enumerate(np.linspace(0.01, 0.06, 20)):
    # # plt.subplot(1, 20, a_ind+1)
    # errors = []
    # for offset_scale in np.linspace(0, np.pi*0.75, 100):
        # pi_remain = np.pi-offset_scale
        # for len_scale in np.linspace(0, pi_remain, 100):
            # ratio = offset_scale / np.pi
            # t_offset = ratio * period
            # ratio = len_scale / np.pi
            # t_len = ratio * period
            # w1 = damped_harmonic_oscillator(amplitude, v, time, w, theta)
            # w2 = -1 * damped_harmonic_oscillator(amplitude, v, time, w, theta)
            # t_w1_p_w2 = add_waves(time, w1, time+t_len, w2)
            # new = add_waves(time, w1_p_w2, time+t_offset, t_w1_p_w2)
            # error = np.sum(np.abs(new[t_period:]))
            # errors.append([offset_scale, len_scale, error])

            # # plt.plot(time+t_offset, w1, label='w1')
            # # plt.plot(time+t_offset+t_len, w2, label='w2')
            # # plt.plot(time, w1_p_w2, label='orig')
            # # plt.plot(time+t_offset, t_w1_p_w2, label='w1+w2')
            # # plt.plot(time, new, label='w1+w2+dip')
            # # plt.legend()
            # # plt.xlim([0, 0.3])
            # # plt.show()
            # if error < best[-1]:
                # best = [amplitude, offset_scale, len_scale, error]

    # errors = np.asarray(errors)
    # print('amp', amplitude, 'min', errors[np.argmin(errors[:, 2])])
    # plt.scatter(errors[:, 0], errors[:, 1], c=errors[:, 2])
    # plt.colorbar()

# plot best
amplitude, offset_scale, len_scale = [0.020526315789473684, 0.30939928, 1.08710453]
# amplitude, len_scale, offset_scale = [amplitude*1.3756, 0.271, 0.0191]
# amplitude, offset_scale, len_scale, e = best
# ratio = offset_scale / np.pi
# t_offset = ratio * period
# ratio = len_scale / np.pi
# t_len = ratio * period
t_len = period * len_scale
t_offset = period * offset_scale
print(t_len, t_offset)
w1 = damped_harmonic_oscillator(amplitude, v, time, w, theta)
w2 = -1 * damped_harmonic_oscillator(amplitude, v, time, w, theta)
t_w1_p_w2 = add_waves(time, w1, time+t_len, w2)
new = add_waves(time, w1_p_w2, time+t_offset, t_w1_p_w2)
error = np.sum(np.abs(new[t_period:]))
print('orig error', np.sum(np.abs(w1_p_w2[t_period:])), 'total', np.sum(np.abs(w1_p_w2)))
print('best new', error, 'total', np.sum(np.abs(new)))

plt.plot(time+t_offset, w1, label='w1')
plt.plot(time+t_offset+t_len, w2, label='w2')
plt.plot(time+t_offset, t_w1_p_w2, label='w1+w2')
# plt.plot(time, new, label='w1+w2+dip')
plt.legend()
plt.xlim([0, 0.3])
plt.show()
