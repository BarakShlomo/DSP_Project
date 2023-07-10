clc;
clear;
[y,Fs] = audioread('File_1.ogg'); % Read the data from the audio file.
Dt = 1/Fs; % Sampling period
N = length(y); % Number of data samples.
t = (0:N-1)*Dt ; % Time vector.
figure(1);
plot (t,y); % Audio signal in plane of time.
grid on ;
xlabel('Time [sec]');% x label name.
ylabel('|S(t)|');% y label name.
title('Signal In The Time Plane'); % graph title.
hold on ;
Delta_f = Fs/N;
Y = fft(y); % FFT on the audio file data.
freq = (0:Delta_f:Fs-Delta_f); % Frequency vector.
figure(2);
subplot(1,2,1);
plot(freq ,abs(Y)); % Audio signal in plane of frequency.
grid on ;
xlabel('Frequency [Hz]');% x label name.
ylabel('|S(f)|');% y label name.
title('Signal In The Frequency Plane'); % graph title.
Max_y = max(y);%maximum signal's value

%Centering code:
center1 = fftshift(Y); % fft shift
f_center = (-N/2:N/2-1)*(Fs/N); % zero-centered frequency range
powercenter1 = abs(center1);     % zero-centered power
subplot(1,2,2);
plot(f_center,powercenter1) % Audio signal in plane of frequency centering.
grid on ;
xlabel('Frequency [Hz]');% x label name.
ylabel('|S(f)|');% y label name.
title('Signal In The Frequency Plane - Centering'); % graph title.

%Pre-warping code:
Freq_c_l = 85; % given in the assignment
Freq_c_h = 155 ; % given in the assignment
F_c_l =(2/Dt)*tan(Freq_c_l*Dt/2);% New low cuttoff frequency
F_c_h = (2/Dt)*tan(Freq_c_h*Dt/2); % New high cuttoff frequency

[a,b] = butter(3,[F_c_l/(Fs/2) , F_c_h/(Fs/2)],'bandpass');% Butterworth filter.
y_har_1 = filter(a,b,y);%Signal affter filtering.

%Normalize code :
Max_y_1=max(y_har_1);%Maximum harmony's value
Norm_y_1 = Max_y/Max_y_1;%Normalization constant
y_har_1 = Norm_y_1*y_har_1;
Y_har_1 = fft(y_har_1);% New signal affter FFT

figure(3);
subplot(1,2,1);
plot (freq,abs(Y_har_1));
grid on; 
xlabel('Frequency [Hz]');% x label name.
ylabel('|S(f)|');% y label name.
title('Signal In The Frequency Plane Affter BPF'); % graph title.

%Centering code:
center_har_1 = fftshift(Y_har_1);%fft shift
powercenter_har_1 = abs(center_har_1);% zero-centered power
subplot(1,2,2);
plot(f_center,powercenter_har_1)
grid on ;
xlabel('Frequency [Hz]');% x label name.
ylabel('|S(f)|');% y label name.
title('Signal In The Frequency Plane Affter BPF - Centering'); % graph title.

%sound(y_har_1,Fs); %Playing the new signal.
%filename = 'File_2.ogg';% Name for the new signal.
%audiowrite(filename,y_har_1,Fs) % Saving the new signal.

%Mean frequency code:
abs_Y_new= abs(Y_har_1);
L_new = length(abs_Y_new);
F_new_Squared = abs_Y_new.*abs_Y_new;
sum1 = 0;
for i = 1:L_new/2
       sum1 = sum1 + ((i*F_new_Squared(i))/(L_new*2))*2*Fs;
end
 sum2 = sum(F_new_Squared);
mean_freq= (sum1/sum2);
center_freq = mean_freq+F_c_l; %zero point in BPF c.o.l frequency.

SD = std(Y_har_1);% Standard deviation code

%Pre-warping code :
Freq_c_l_2 = 2*center_freq - SD; % given in the assignment
Freq_c_h_2 = 2*center_freq + SD ; % given in the assignment
F_c_l_2 =(2/Dt)*tan(Freq_c_l_2*Dt/2);% New low cuttoff frequency
F_c_h_2 = (2/Dt)*tan(Freq_c_h_2*Dt/2); % New high cuttoff frequency

[c,d] = butter(3,[F_c_l_2/(Fs/2) , F_c_h_2/(Fs/2)],'bandpass');% Butterworth filter.
y_har_2 = filter(c,d,y);%Signal affter filtering.
%Normalize code :
Max_y_2=max(y_har_2);%Maximum harmony's value
Norm_y_2 = Max_y/Max_y_2;%Normalization constant
y_har_2 = Norm_y_2*y_har_2;
Y_har_2 = fft(y_har_2);% New signal affter FFT
figure(4);
subplot(1,2,1);
plot (freq,abs(Y_har_2));
grid on; 
xlabel('Frequency [Hz]');% x label name.
ylabel('|S(f)|');% y label name.
title('Signal Second Harmony'); % graph title.

%Centering code:
center_har_2 = fftshift(Y_har_2); %fft shift
powercenter_har_2 = abs(center_har_2);% zero-centered power
subplot(1,2,2);
plot(f_center,powercenter_har_2)
grid on ;
xlabel('Frequency [Hz]');% x label name.
ylabel('|S(f)|');% y label name.
title('Signal Second Harmony - Centering'); % graph title.

%sound(y_har_2,Fs); %Playing the new signal.
%filename = 'File_3.ogg';% Name for the new signal.
%audiowrite(filename,y_har_2,Fs) % Saving the new signal.

%Pre-warping code :
Freq_c_l_3 = 3*center_freq - SD; % given in the assignment
Freq_c_h_3 = 3*center_freq + SD ; % given in the assignment
F_c_l_3 =(2/Dt)*tan(Freq_c_l_3*Dt/2);% New low cuttoff frequency
F_c_h_3 = (2/Dt)*tan(Freq_c_h_3*Dt/2); % New high cuttoff frequency

[e,f] = butter(3,[F_c_l_3/(Fs/2) , F_c_h_3/(Fs/2)],'bandpass');% Butterworth filter.
y_har_3 = filter(e,f,y_har_2);%Signal affter filtering.

%Normailze code:
Max_y_3=max(y_har_3);% Maximum harmony's value 
Norm_y_3 = Max_y/Max_y_3;%Normalization constant
y_har_3 = Norm_y_3*y_har_3;
Y_har_3 = fft(y_har_3);% New signal affter FFT

figure(5);
subplot(1,2,1);
plot (freq,abs(Y_har_3));%Harmony's signal in frequency plane.
grid on; 
xlabel('Frequency [Hz]');% x label name.
ylabel('|S(f)|');% y label name.
title('Signal Third Harmony'); % graph title.

%Centering code:
center_har_3 = fftshift(Y_har_3);%fft shift
powercenter_har_3 = abs(center_har_3);% zero-centered power
subplot(1,2,2);
plot(f_center,powercenter_har_3);
grid on ;
xlabel('Frequency [Hz]');% x label name.
ylabel('|S(f)|');% y label name.
title('Signal Third Harmony - Centering'); % graph title.

%sound(y_har_3,Fs); %Playing the new signal.
%filename = 'File_4.ogg';% Name for the new signal.
%audiowrite(filename,y_har_3,Fs) % Saving the new signal.

cosine = cos(2*pi*center_freq*t); % cosine wave
y_1 = y_har_1.*cosine'; % Harmony 1 After multiplying by cosine

%Normailize code :
Max_y1=max(y_1); % Maximum harmony 1 affter shift value
Norm_y1 = Max_y/Max_y1;%Normalization constant
y_1 = Norm_y1*y_1;
Y_1 = fft(y_1);%New signal affter fft
figure(6);
subplot(1,2,1);
plot (freq,abs(Y_1));
grid on;
xlabel('Frequency [Hz]');% x label name.
ylabel('|S(f)|');% y label name.
title('First Harmony Shift'); % graph title.

%Centering code:
center_Y_1 = fftshift(Y_1);%fft shift
powercenter_Y_1 = abs(center_Y_1);% zero-centered power
subplot(1,2,2);
plot(f_center,powercenter_Y_1)
grid on ;
xlabel('Frequency [Hz]');% x label name.
ylabel('|S(f)|');% y label name.
title('First Harmony Shift - Centering'); % graph title.

%Pre-warping code:
Freq_c_l_shift = 85 + center_freq; % given in the assignment
Freq_c_h_shift = 155 + center_freq ; % given in the assignment
F_c_l_shift =(2/Dt)*tan(Freq_c_l_shift*Dt/2);% New low cuttoff frequency
F_c_h_shift = (2/Dt)*tan(Freq_c_h_shift*Dt/2); % New high cuttoff frequency

[A,B] = butter(3,[F_c_l_shift/(Fs/2) , F_c_h_shift/(Fs/2)],'bandpass');% Butterworth filter.
y_har1 = filter(A,B,y_1);%Signal affter filtering.

%Normailize code :
Max_y1=max(y_har1);%Maximum harmony's value
Norm_y1 = Max_y/Max_y1;%Normalization constant
y_1 = Norm_y1*y_1;
Y_har1 = fft(y_har1);% New signal affter fft

figure(7);
subplot(1,2,1);
plot (freq,abs(Y_har1));
grid on;
xlabel('Frequency [Hz]');% x label name.
ylabel('|S(f)|');% y label name.
title('First Harmony Shift Affter BPF'); % graph title.

%Centering code:
center_har1 = fftshift(Y_har1);%fft shift
powercenter_har1 = abs(center_har1);% zero-centered power
subplot(1,2,2);
plot(f_center,powercenter_har1)
grid on ;
xlabel('Frequency [Hz]');% x label name.
ylabel('|S(f)|');% y label name.
title('First Harmony Shift Affter BPF - Centering'); % graph title.

%sound(y_har1,Fs); %Playing the new signal.
%filename = 'File_5.ogg';% Name for the new signal.
%audiowrite(filename,y_har1,Fs) % Saving the new signal.

y_2 = y_har_2.*cosine';% Harmony 2 After multiplying by cosine

%Normailize code:
Max_y2=max(y_2);% Maximum harmony 2 affter shift value
Norm_y2 = Max_y/Max_y2;%Normalization constant
y_2 = Norm_y2*y_2;
Y_2 = fft(y_2);%New harmony shift affter fft
figure(8);
subplot(1,2,1);
plot (freq,abs(Y_2));
grid on;
xlabel('Frequency [Hz]');% x label name.
ylabel('|S(f)|');% y label name.
title('Second Harmony Shift'); % graph title.

%Centering code:
center_Y_2 = fftshift(Y_2);%fft shift
powercenter_Y_2 = abs(center_Y_2);% zero-centered power
subplot(1,2,2);
plot(f_center,powercenter_Y_2)
grid on ;
xlabel('Frequency [Hz]');% x label name.
ylabel('|S(f)|');% y label name.
title('Second Harmony Shift - Centering'); % graph title.

%Pre-warping code:
Freq_c_l_shift2 = 2*center_freq-SD; % given in the assignment
Freq_c_h_shift2 = 2*center_freq+SD ; % given in the assignment
F_c_l_shift2 =(2/Dt)*tan(Freq_c_l_shift2*Dt/2);% New low cuttoff frequency
F_c_h_shift2 = (2/Dt)*tan(Freq_c_h_shift2*Dt/2); % New high cuttoff frequency

[C,D] = butter(3,[F_c_l_shift2/(Fs/2) , F_c_h_shift2/(Fs/2)],'bandpass');% Butterworth filter.
y_har2 = filter(C,D,y_1);%Signal affter filtering.

%Normailze code:
Max_y2=max(y_har2);%Maximum harmony 2 affter shift and bpf value
Norm_y2 = Max_y/Max_y2;%Normalization constant
y_har2 = Norm_y2*y_har2;
Y_har2 = fft(y_har2);% New harmony 2 affter shift and bpf

figure(9);
subplot(1,2,1);
plot (freq,abs(Y_har2));
grid on;
xlabel('Frequency [Hz]');% x label name.
ylabel('|S(f)|');% y label name.
title('Second Harmony Shift Affter BPF'); % graph title.

%Centering code:
center_har2 = fftshift(Y_har2);%fft shift
powercenter_har2 = abs(center_har2);% zero-centered power
subplot(1,2,2);
plot(f_center,powercenter_har2)
grid on ;
xlabel('Frequency [Hz]');% x label name.
ylabel('|S(f)|');% y label name.
title('Second Harmony Shift Affter BPF - Centering'); % graph title.

%sound(y_har2,Fs); %Playing the new signal.
%filename = 'File_6.ogg';% Name for the new signal.
%audiowrite(filename,y_har2,Fs) % Saving the new signal.

y_3 = y_har_3.*cosine';% Harmony 2 After multiplying by cosine

%Normalize code : 
Max_y3=max(y_3);%Maximum harmony 3 affter shift value
Norm_y3 = Max_y/Max_y3;%Normalization constant
y_3 = Norm_y3*y_3;
Y_3 = fft(y_3);% New harmony 3 affter shift and bpf

figure(10);
subplot(1,2,1);
plot (freq,abs(Y_3));
grid on;
xlabel('Frequency [Hz]');% x label name.
ylabel('|S(f)|');% y label name.
title('Third Harmony Shift'); % graph title.
%centering code:
center_Y_3 = fftshift(Y_3);%fft shift
powercenter_Y_3 = abs(center_Y_3);% zero-centered power
subplot(1,2,2);
plot(f_center,powercenter_Y_3)
grid on ;
xlabel('Frequency [Hz]');% x label name.
ylabel('|S(f)|');% y label name.
title('Third Harmony Shift - Centering'); % graph title.

y_har3 = filter(A,B,y_3);%Signal affter filtering.

%Normailze code:
Max_y3=max(y_har3);%Maximum harmony 3 affter shift and bpf value
Norm_y3 = Max_y/Max_y3;%Normailzation constant
y_3 = Norm_y3*y_har3;
Y_har3 = fft(y_har3);% New signal affter FFT

figure(11);
subplot(1,2,1);
plot (freq,abs(Y_har3));
grid on;
xlabel('Frequency [Hz]');% x label name.
ylabel('|S(f)|');% y label name.
title('Third Harmony Shift Affter BPF'); % graph title.

%Centering code:
center_har3 = fftshift(Y_har3);%fft shift
powercenter_har3 = abs(center_har3);% zero-centered power
subplot(1,2,2);
plot(f_center,powercenter_har3)
grid on ;
xlabel('Frequency [Hz]');% x label name.
ylabel('|S(f)|');% y label name.
title('Third Harmony Shift Affter BPF - Centering'); % graph title.

%sound(y_har3,Fs); %Playing the new signal.
%filename = 'File_7.ogg';% Name for the new signal.
%audiowrite(filename,y_har3,Fs) % Saving the new signal.

y_new = y_har1 + y_har2 + y_har3;% The sum of the harmonics after shifting

%Normalize code :
Max_y_new=max(y_new);%Maximum sum of harmonics value
Norm_y_new = Max_y/Max_y_new;%Normailzation constant
y_1 = Norm_y_new*y_new;
Y_new = fft(y_new);% New signal affter FFT

figure(12);
subplot(1,2,1);
plot (freq,abs(Y_new));
grid on;
xlabel('Frequency [Hz]');% x label name.
ylabel('|S(f)|');% y label name.
title('New Signal'); % graph title.

%Centering code:
center_y_new = fftshift(Y_new);%fft shift
powercenter_Y_new = abs(center_y_new);% zero-centered power
subplot(1,2,2);
plot(f_center,powercenter_Y_new)
grid on ;
xlabel('Frequency [Hz]');% x label name.
ylabel('|S(f)|');% y label name.
title('New Signal - Centering'); % graph title.

%sound(y_new,Fs); %Playing the new signal.
%filename = 'File_8.ogg';% Name for the new signal.
%audiowrite(filename,y_new,Fs) % Saving the new signal.

%RLS Filter:
N_RLS=2;%Filter's order
p0 = 2 * eye(N_RLS); % Initial Inverse Covariance
lambda = 0.91; % Forgetting Factor
rls = dsp.RLSFilter(N_RLS,'ForgettingFactor',lambda,...
   'InitialInverseCovariance',p0);
[signal,err] = rls(y_new,y);

%filename = 'File_9.ogg';% Name for the new signal.
%audiowrite(filename,signal,Fs) % Saving the new signal.

figure(13); freqz(err,signal,[],Fs);

SIGNAL = fft(signal);
figure(14);
subplot(1,2,1);
plot (freq,abs(SIGNAL));
grid on;
xlabel('Frequency [Hz]');% x label name.
ylabel('|S(f)|');% y label name.
title('Signal Affter RLS Filter'); % graph title.

%centering code:
center_SIGNAL = fftshift(SIGNAL);%fft shift
powercenter_SIGNAL = abs(center_SIGNAL);% zero-centered power
subplot(1,2,2);
plot(f_center,powercenter_SIGNAL)
grid on ;
xlabel('Frequency [Hz]');% x label name.
ylabel('|S(f)|');% y label name.
title('Signal Affter RLS Filter - Centering'); % graph title.
