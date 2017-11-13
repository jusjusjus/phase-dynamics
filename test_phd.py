
# Can only be set before all other imports.
from phd.core import set_period
# set_period(1.0)

import phd
import numpy as np
from phd.core import period
import phd.matplotlib as pht
import matplotlib.pyplot as plt


def test_matplotlib():
    plt.plot(np.arange(0., 10.))
    pht.xticks_centered(fontsize=20)
    pht.yticks_centered(fontsize=20)
    plt.show()


def test_unmod():
    x = period * np.arange(0.0, 10., 0.1)
    x_mod = np.mod(x, period)
    y = phd.unmod(x_mod)
    assert np.all(x-y == 0.0)


def test_threshold_data():
    dt = 0.1
    t = np.arange(0., 20., dt)
    x = np.sin(t)
    threshold = 0.9
    segments = phd.threshold_data(x, threshold)
    for seg in segments:
        assert all(x[seg] > threshold)
    # plt.plot(t, x, 'ko-')
    # for seg in segments:
    #     plt.plot(t[seg], x[seg], 'ro-', lw=1.5)
    #     plt.axvline(x=dt*seg.start)
    #     plt.axvline(x=dt*seg.stop)
    # plt.axhline(y=threshold)
    # plt.show()


def test_poincare_times():
    assert period == 2.0*np.pi
    dt = 0.1
    T  = 1.0
    x0 = 0.#3.0*np.pi/2

    w = 2.0*np.pi/T
    t = T * np.arange(0., 10., dt)
    Wt = 0.2 * np.sqrt(dt) * np.random.randn(t.size).cumsum()
    x = w * t + Wt
    x_mod = phd.mod(x)

    idx, ti = phd.poincare_times(x, x0=x0, interp=True)
    Ti = dt * ti
    
    assert all(Ti > 0.9 * T)

    # plt.plot(t, x_mod, 'ko-')
    # for i in idx:
    #     plt.axvline(x=dt*i, lw=2.0)
    # plt.axhline(y=x0)
    # idx = idx.astype(int)
    # plt.plot(t[idx], x[idx], 'ro')
    # pht.yticks()
    # plt.show()


def test_gradient():
    dt = 0.1
    T  = 1.0
    x0 = 3.0*np.pi/2

    w = 2.0*np.pi/T
    t = T * np.arange(0., 10., dt).astype(np.float64)
    x = w * t
    x_mod = phd.mod(x)
    empirical = {
        md: phd.core.gradient(x_mod, mode=md)/dt
        for md in phd.core._gradient
    }
    for md, emp in empirical.items():
        assert all(abs(w-emp)<10.0**-12), 'error in mode {}'.format(md)


test_unmod()
test_threshold_data()
test_poincare_times()
test_gradient()
