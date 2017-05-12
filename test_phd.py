
# Can only be set before all other imports.
from phd.core import set_period
# set_period(1.0)

import numpy as np
import phd.matplotlib as pht
import matplotlib.pyplot as plt
from phd.core import unmod, threshold_data, period

# plt.plot(np.arange(0., 10.))
# pht.xticks_centered(fontsize=20)
# pht.yticks_centered(fontsize=20)
# plt.show()


def test_unmod():
    x = period * np.arange(0.0, 10., 0.1)
    x_mod = np.mod(x, period)
    y = unmod(x_mod)
    assert np.all(x-y == 0.0)


def test_threshold_data():
    dt = 0.1
    t = np.arange(0., 20., dt)
    x = np.sin(t)
    threshold = 0.9
    segments = threshold_data(x, threshold)
    for seg in segments:
        assert all(x[seg] > threshold)
    # plt.plot(t, x, 'ko-')
    # for seg in segments:
    #     plt.plot(t[seg], x[seg], 'ro-', lw=1.5)
    #     plt.axvline(x=dt*seg.start)
    #     plt.axvline(x=dt*seg.stop)
    # plt.axhline(y=threshold)
    # plt.show()


test_unmod()
test_threshold_data()
