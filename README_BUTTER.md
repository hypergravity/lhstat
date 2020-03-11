## The Butter Worth filter

### The usage of the function:
``` python
from scipy import signal
signal.butter(N, Wn, btype='low', analog=False, output='ba')
```
- N: the order of the filter
- Wn: normalized cutoff frequency
- btype: type of filter {'lowpass', 'highpass', 'bandpass', 'bandstop'},

### The SQM filter is a low pass BW filter:
```python
import numpy as np
from scipy import signal
def lfilter(x, y, N=8, Wn=0.05, btype='low', analog=False, output='ba'):
    b, a = signal.butter(N, Wn, btype=btype, analog=analog, output=output)
    padlen = np.int(0.3 * len(x))
    yf = signal.filtfilt(b, a, y, padlen=padlen)
    return yf
```

### The Parameters of the filter
The `Wn` is determined using `Wn=2*F_cutoff/F_sampling`.
The `F_sampling`, the sampling frequency is typically 1/1min.
`Wn=0.05` means `F_cutoff=Wn*F_sampling/2=0.05*(1/1min)/2=1/40min`.
Therefore, this function filters the signals whose timescale are below 40 minutes.

### The sunny/cloudy criterion
A moving standard deviation method calculates the std at a given time within (plus/minus) 10 min.
This moving std is compared to a threshold, which is set to 0.03 mag in our code. 