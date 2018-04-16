from tables import open_file
import matplotlib.pyplot as plt

def get_hist(data, x_lab, title):
    plt.figure()
    plt.title(title)
    ax = plt.gca()
    ref_val = data[0]
    ax.hist(data)
    plt.axvline(ref_val, color='red')
    ax.set_xlabel(x_lab)
    ax.set_ylabel('Counts')
    ax.get_xaxis().get_major_formatter().set_useOffset(False)


f = open_file('./data/decoherence_test.h5', 'r')

bunches = f.root.bunch

Wx_dK = bunches.dK.Wx_hist.read()
Wy_dK = bunches.dK.Wy_hist.read()
tilt_dK = bunches.dK.tilt_hist.read()

Wx_dK_means = Wx_dK.mean(axis=1)
Wx_dK_std = Wx_dK.std
tilt_dK_means = tilt_dK.mean(axis=1)


plt.rc('text', usetex=True)
font = {'family' : 'normal',
        'weight' : 'bold',
        'size'   : 16}

plt.rc('font', **font)

plt.plot(tilt_dK_means, Wy_dK[:,0], '.')
plt.xlabel('Mean tilt angle, $\Theta$ [rad]')
plt.ylabel('Reference horizontal spin precession frequency, $\Omega_y$ [rad/sec]')

i = tilt_dK_means.argmax()
tilt = tilt_dK_means[i]
get_hist(Wx_dK[i], r'$\Omega_x$ [rad]', r'$\Omega_x$ histogram for mean tilt angle {:4.2e} rad'.format(tilt))
get_hist(Wy_dK[i], r'$\Omega_y$ [rad]', r'$\Omega_y$ histogram for mean tilt angle {:4.2e} rad'.format(tilt))
i = tilt_dK_means.argmin()
tilt = tilt_dK_means[i]
get_hist(Wx_dK[i], r'$\Omega_x$ [rad]', r'$\Omega_x$ histogram for mean tilt angle {:4.2e} rad'.format(tilt))
## didn't include Wy because it's the same distribution for all tilts

plt.plot(tilt_dK_means, Wy_dK[:, 0], '.')
plt.xlabel('Mean tilt angle, $\Theta$ [rad]')
plt.ylabel('Reference horizontal spin precession frequency, $\Omega_y$ [rad/sec]')

