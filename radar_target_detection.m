clear all
clc;
close all

%% Radar Specifications
%%%%%%%%%%%%%%%%%%%%%%%%%%%
fc= 77e9;             %Operating carrier frequency of Radar = 77GHz
Rmax = 200;           %Max Range Rmax = 200m
Rresolution = 1;      %Range resolution = 1 m
Vmax = 70;            % Max Velocity = 70 m/s
Vresolution = 3;      % Velocity resolution = 3 m/s
c= 3e8;               %speed of light = 3e8
%%%%%%%%%%%%%%%%%%%%%%%%%%%

R = 110;              %% User Defined Range of target
V = -20;              %% User Defined Velocity of target

%% FMCW Waveform Generation

%Design the FMCW waveform by giving the specs of each of its parameters.
% Calculate the Bandwidth (B), Chirp Time (Tchirp) and Slope (slope) of the FMCW
Tchirp = 5.5*(2*Rmax)/c; %Sweep time according to the rule = 5.5*(2*Rmax)/c;
Bandwidth = c / (2*Rresolution); % Bandwidth for given resolution - Range resolution = c/(2*Band_width);

Slope = Bandwidth / Tchirp;

%The number of chirps in one sequence  FFT
%for Doppler Estimation.
Nd=128; % #of doppler cells OR #of sent periods % number of chirps Power of 2 for faster FFT

%The number of samples on each chirp.
Nr=1024;                  %for length of time OR # of range cells

% Timestamp for running the displacement scenario for every sample on each chirp
t=linspace(0,Nd*Tchirp,Nr*Nd); %total time for samples

%Creating the vectors for Tx, Rx and Mix based on the total samples input.
Tx=zeros(1,length(t)); %transmitted signal
Rx=zeros(1,length(t)); %received signal
Mix = zeros(1,length(t)); %beat signal

%Similar vectors for range_covered and time delay.
r_t=zeros(1,length(t));
td=zeros(1,length(t));

%% Signal generation and Moving Target simulation
% Running the radar scenario over the time.
for i=1:length(t)
    r_t(i) = R + V * t(i);  %Range of the Target for constant velocity.
    td(i) = 2 * r_t(i) / c; %The delay in signal

    Tx(i) = cos(2*pi*(fc*(t(i)      ) + Slope* t(i)^2       /2)); %transmitted
    Rx(i) = cos(2*pi*(fc*(t(i)-td(i)) + Slope*(t(i)-td(i))^2/2)); %received signal.

    Mix(i) = Tx(i) * Rx(i); %beat signal by mixing the Transmit and Receive
end

%% RANGE MEASUREMENT

Mix_reshaped = reshape(Mix, [Nr, Nd]); %reshape the vector into Nr*Nd array
fftMix = fft(Mix_reshaped)/Nr; %FFT on the beat signal and normalize.

absfftMix2 = abs(fftMix);  % The absolute value of FFT output
absfftMix1 = absfftMix2(1:Nr/2+1,:); % Throw out half of the samples.

%plotting the range
figure ('Name','Range from First FFT')
plot(absfftMix1(:,10))
axis ([0 200 0 1]);

% Range Doppler Map Generation.

% The output of the 2D FFT is an image that has reponse in the range and
% doppler FFT bins. So, it is important to convert the axis from bin sizes
% to range and doppler based on their Max values.

Mix=reshape(Mix,[Nr,Nd]);

% 2D FFT using the FFT size for both dimensions.
sig_fft2 = fft2(Mix,Nr,Nd);

% Taking just one side of signal from Range dimension.
sig_fft2 = sig_fft2(1:Nr/2,1:Nd);
sig_fft2 = fftshift (sig_fft2);
RDM = abs(sig_fft2);
RDM = 10*log10(RDM) ;

%use the surf function to plot the output of 2DFFT and to show axis in both
%dimensions
doppler_axis = linspace(-100,100,Nd);
range_axis = linspace(-200,200,Nr/2)*((Nr/2)/400);
figure,surf(doppler_axis,range_axis,RDM);

%% CFAR implementation

Tr = 5; % Range training
Td = 5;  % Doppler training

Gr = 5; % Guard cells in range
Gd = 5; % Guard cells in doppler

offset = 12; % offset the threshold by SNR value in dB

%Create a vector to store noise_level for each iteration on training cells
noise_level = zeros(1,1);

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

figure, surf(doppler_axis,range_axis,result); %display the CFAR
colorbar;

