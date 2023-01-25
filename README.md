# Udacity Sensor Fusion Nanodegree - Project 3: Radar Target Detection

In this project a single matlab code is completed to simulate transmission and receive of the reflected radar waves from an object at a distance of 110 m approaching at a velocity of 20 m/s.

The parameters for the Frequency Modulated Continuous Wave radar are to be found such that we can have a range resolution of 1 m and a velocity resolution of 3 m/s. 

## Project rubric

### Using the given system requirements, design a FMCW waveform. Find its Bandwidth (B), chirp time (Tchirp) and slope of the chirp.
The sweep time, `Tchirp` is found from the relation `Tchirp = 5.5*(2*Rmax)/c` and is found to be 7.33 microseconds. The bandwidth is related to the required resolution and is found from the relation `Bandwidth = c / (2*Rresolution)` to be 150 MHz. These two values give the slope found from `Slope = Bandwidth / Tchirp` to be 2.0455e+13 Hz/second.

### Simulate Target movement and calculate the beat or mixed signal for every timestamp.

The given code is completed as below to first calculate the current range and then the delay in the reflected signal. Thereafter, the two signals are mixed as shown below. Note that these operations can be written as element-wise operations that will be executed faster in Matlab. 

```Matlab
%% Signal generation and Moving Target simulation
% Running the radar scenario over the time.
for i=1:length(t)
    r_t(i) = R + V * t(i);  %Range of the Target for constant velocity.
    td(i) = 2 * r_t(i) / c; %The delay in signal

    Tx(i) = cos(2*pi*(fc*(t(i)      ) + Slope* t(i)^2       /2)); %transmitted
    Rx(i) = cos(2*pi*(fc*(t(i)-td(i)) + Slope*(t(i)-td(i))^2/2)); %received signal.

    Mix(i) = Tx(i) * Rx(i); %beat signal by mixing the Transmit and Receive
end
```


### Implement the Range FFT on the Beat or Mixed Signal and plot the result.
After mixing, the vector is reshaped and a fourier transform in the range direction is performed.

```Matlab
Mix_reshaped = reshape(Mix, [Nr, Nd]); %reshape the vector into Nr*Nd array
fftMix = fft(Mix_reshaped)/Nr; %FFT on the beat signal and normalize.

absfftMix2 = abs(fftMix);  % The absolute value of FFT output
absfftMix1 = absfftMix2(1:Nr/2+1,:); % Throw out half of the samples.
```

The figure below shows the one vector output of 1D FFT in the range dimension. The total simulation time is 0.94 ms. The target has moved only about 2 cm during this time.

<center>
  <img src="media/fft1.bmp" width="779" height="414" />
</center>


### Implement the 2D CFAR process on the output of 2D FFT operation, i.e the Range Doppler Map.
The 2D FFT in both range and doppler dimensions was implemented in the code. After that the FFT is performed, first, the positive half of the ranges is taken and then the matrix is shifted in the doppler dimension so that zero will be in the middle of the axis. Finally, the result is expressed in the Db scale and plotted.  

```Matlab
sig_fft2 = fft2(Mix,Nr,Nd);

% Taking just one side of signal from Range dimension.
sig_fft2 = sig_fft2(1:Nr/2,1:Nd);
sig_fft2 = fftshift (sig_fft2);
RDM = abs(sig_fft2);
RDM = 10*log10(RDM) ;
```
Figure below presents `RDM` in a range-doppler map.

<center>
  <img src="media/fft2.bmp" width="779" height="414" />
</center>

The final step is to apply CFAR in order to subpress the noise and find the target. This is done by first defining a training band around the cell under test with a distance to the cell in order to avoid influence of the signal on estimate of the noise. In the below code, the width of the training band in the range dimension is denoted by `Tr` and in the doppler dimensino by `Td`. The distance from the training band to the cell under test is denoted by `Gr` and `Gd` standing for Guard cells in range and doppler dimensions.

```Matlab
Tr = 5; % Range training
Td = 5; % Doppler training

Gr = 5; % Guard cells in range
Gd = 5; % Guard cells in doppler
```

The method is to average the signal in the training band for each cell and then use the average value as a measure for noise level. In order to average the noise level, one has to convert the decibel values to power values first and then once the average value is calculated, convert the level back to decibel values.

The algorithm aims to find the cells that have a value larger than the local noise level, but before making the comparison, the decibel values must be offset but a value related to the signal noise ratio.


```Matlab
result = zeros(size(RDM));
pRDM = db2pow(RDM);
NumberOfCells = (2*(1+2*Tr+2*Gr)*Td+2*Td*Tr);
NLim = size(RDM);

for r = 1 + Tr + Gr : NLim(1) - Tr - Gr
  for d = 1 + Td + Gd : NLim(2) - Td - Gd
     noise_level = (sum(sum(pRDM(r-Tr-Gr:r+Tr+Gr,d+Gd+1 :d+Gd+Td)))+...
                        sum(sum(pRDM(r-Tr-Gr:r+Tr+Gr,d-Gd-Td:d-Gd-1 )))+...
                        sum(sum(pRDM(r-Tr-Gr:r-Gr-1 ,d-Gd-Td:d-Gd-1 )))+...
                        sum(sum(pRDM(r+Gr+1 :r+Gr+Tr,d+Gd+1 :d+Gd+Td))))/NumberOfCells;
     if (RDM(r,d) > pow2db(noise_level)+offset)
       result(r,d) = 1;
     end
  end
end
```

Once this is done, the cells with a value higher than the local noise value multiplied by the signal-noise ratio are selected and the result is set to 1 for those cells. Note that the mentioned product becomes an additive contribution in the decibel scale, i.e. the comparison becomes 

```Matlab
if (RDM(r,d) > pow2db(noise_level)+offset)
   result(r,d) = 1;
end
```
The figure below shows the final result in the range-doppler map with `Tr = 5, Td = 5, Gr = 5, Gd = 5` .
<center>
  <img src="media/cfar.bmp" width="779" height="414" />
</center>

With `Tr = 10, Td = 8, Gr = 4, Gd = 4` we get 
<center>
  <img src="media/cfar2.bmp" width="779" height="414" />
</center>

and finally we will have 
<center>
  <img src="media/cfar3.bmp" width="779" height="414" />
</center>

with  `Tr = 10, Td = 8, Gr = 1, Gd = 1` .
