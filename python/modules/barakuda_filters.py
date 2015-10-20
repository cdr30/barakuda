
# A collection of filters for low-pass filtering ans spectral-analysis stuff...

# Laurent Brodeau, Jan. 2015

import numpy as nmp

from scipy.interpolate import UnivariateSpline
from scipy.signal      import wiener, filtfilt, butter, gaussian, freqz
from scipy.ndimage     import filters

#import scipy.optimize as op

# For spectral analysis:
#from scipy.fftpack import fft



def ssqe(sm, s, npts):
    return nmp.sqrt(nmp.sum(nmp.power(s-sm,2)))/npts


 
def testGauss(x, y):
    b = gaussian(39, 10)
    #ga = filtfilt(b/b.sum(), [1.0], y)
    ga = filters.convolve1d(y, b/b.sum())
    return ga
 
def testButterworth(nyf, x, y):
    b, a = butter(4, 1.5/nyf)
    fl = filtfilt(b, a, y)
    return fl
 
def testWiener(x, y):
    #wi = wiener(y, mysize=29, noise=0.5)
    wi = wiener(y, mysize=29)
    return wi
 
def testSpline(x, y, rS):
    sp = UnivariateSpline(x, y, k=4, s=rS)
    return sp(x)

#def plotPowerSpectrum(y, w):
#    ft = nmp.fft.rfft(y)
#    ps = nmp.real(ft*nmp.conj(ft))*nmp.square(dt)








def Amp_Spctrm(vy, rdt=1., lwin=False):

    # L. Brodeau, 2015
    # Computes a Single-Sided Amplitude Spectrum of vy(t)
    #
    # In general, you should apply a window function to your time domain data
    # before calculating the FFT, to reduce spectral leakage. The Hann window is
    # almost never a bad choice. You can also use the rfft function to skip the
    # -FSample/2, 0 part of the output.

    # Depending on the implementation, your FFT routine might output yout =
    # yin/N, or yout = yin/sqrt(N), or if you applied a window, you have to
    # replace N with a coefficient that depends on the chosen window.
    
    N = len(vy) # length of the signal

    if lwin:
        sY = 2.*nmp.fft.rfft(vy*nmp.hanning(N))/(float(N)/2.) # fft computing and normalization (Fast Fourier Transform)
    else:
        sY = nmp.fft.rfft(vy)/(float(N)/2.)      # fft computing and normalization (Fast Fourier Transform)

    sY = nmp.abs(sY[:N/2])          # Taking real part, and symetric around 0, only keep half!

    vfreq = nmp.fft.fftfreq(N, d=rdt) # in hours, or d=1.0/24 in days
    vfreq = vfreq[:N/2]

    # The inverse Fourier Transform would be: numpy.fft.ifft

    return vfreq, sY


def Pow_Spctrm(vy, rdt=1., lwin=False):

    [ vfreq, spY ] = Amp_Spctrm(vy, lwin=lwin, rdt=rdt)
    spY = spY*spY
    return vfreq, spY


