%-- test cases : fe_curve --%
comgui('close all clear');

%-- --%

fs=1024;  %-- sampling rate
aqc_length=4;  %-- time signals length (in seconds)

%-- spectrum definition for noise generation --%

%-- Define a spectrum model --%
Spec=[10*ones(1,80) 10:50 50*ones(1,50) ...
        30*ones(1,50) 10*ones(1,50) 10:-.5:0 0:2:50 50:-1:0];

Noise=fe_curve('Noise',fs,fs,Spec);
Time=Noise.X{1};
Noise_t=Noise.Y;

Noise_w=fft(Noise_t);
Noise_w=Noise_w(1:length(Time)/2+1);
Freq=fs*(0:length(Time)-1)/length(Time)*2*pi;
Freq=Freq(1:length(Time)/2+1);

%-- Compute noise power spectrum density --%

Pyy = Noise_w .* conj(Noise_w) / length(Noise_w);

figure(1); subplot(2,1,1);plot(Time,Noise_t);
title('Noise');xlabel('s');
subplot(2,1,2);
plot(Freq/2/pi,Pyy);title('Power spectrum density of generated noise');
xlabel('Hz');

%-- 3 DOF oscillator definition --%

Puls = [20 150 270]*2*pi;  %-- natural frequencies
Damp = [.005 .003 .008]*5;  %-- damping
Amp = [1 2 -1;2 -1 1;-1 1 2];  %-- "mode shapes"
Amp=Amp./det(Amp);

C=[1 0 0];  %-- Observation matrix
B=[0 1 0]';  %-- Command matrix

%-- Compute transfert function of oscillator 

Trans_w=[];

for i1=1:length(Freq)
  Trans_w(i1)=C*Amp*...
    (diag(1./(Puls.^2-Freq(i1).^2*[1 1 1]+2*j*Freq(i1)*Damp.*Puls))...
    *Amp'*B);
end;

%-- Compute time response of oscillator to noise

Rep_t=fe_curve('TimeFreq',Noise,[Trans_w zeros(1,length(Trans_w)-2)].');

%-- Compute FRF estimations from time signals

frame{1}.X=Time';
frame{1}.Y=[Noise.X{1}(:) Rep_t.Y];
out=fe_curve('H1H2',frame,'hanning');
figure(2);
semilogy(out.X(2:end),abs(out.H1(2:end)),':b','linewidth',1.5);
hold on
semilogy(out.X(2:end),abs(out.H2(2:end)),'-.m','linewidth',1);
legend('H1 estimator','H2 estimator');

%-- Compute filtered response using true bandpass filter --%

out5=fe_curve('BandPass Hz 10 30',frame);
figure(3);
plot(frame{1}.X,frame{1}.Y(:,2));
hold on;
plot(out5.X,out5.Y(:,2),'r');
legend('Initial response','Filtered response');
title('Bandpass filtered signal between 10 & 30 Hz')

%-- Curve handling utilities --%

clear model

%-- definition of a curve --%

model.Stack{1,1}='curve';
model.Stack{1,2}='Time response';
model.Stack{1,3}.ID=8;
model.Stack{1,3}.X=Time;
model.Stack{1,3}.Y=Rep_t;
model.Stack{1,3}.xunit.label='time';
model.Stack{1,3}.xunit.unit='(s.)';
model.Stack{1,3}.yunit.label='displ.';
model.Stack{1,3}.yunit.unit='(m.)';
model.Stack{1,3}.name='Response of DOF 1 to a noisy input at DOF 2';

%-- plots time response using fe_curve --%

%fe_curve('plot',model,'Time response');

%-- compute time reponse using fe_curve --%
%
% /!\ /!\ /!\ /!\ /!\ /!\ /!\ /!\
%
% Beware of the sampling frequency of the time signal 
% to be coherent with frequency definition of given transfert
%
% /!\ /!\ /!\ /!\ /!\ /!\ /!\ /!\ 

val=frame{1};
val.Y=val.Y(:,2);
val.yunit=fe_curve('datatype','Displacement');

freq=linspace(0,1000,100);
damp=.01:.01:.05;

out1=fe_curve('ResSpectrum T R A',val,freq,damp);
view(0,0);
out2=fe_curve('ResSpectrum P A A',val,freq,damp);
view(0,0);

wd=fileparts(which('fe_curve'))

st=sprintf('read %s',fullfile(wd,'test','bagnol_ns.cyt'));

bagnol_ns=fe_curve(st);
bagnol_ns.yunit=fe_curve('datatype','Acceleration');
st=sprintf('read %s',fullfile(wd,'test','bagnol_ns_rspec_pa.cyt'));

bagnol_ns_rspec_pa=fe_curve(st);

cur=fe_curve('ResSpectrum True Rel. Acc.',bagnol_ns,bagnol_ns_rspec_pa.X/2/pi,.05);

fe_curve('plot',cur);
hold on;
plot(cur.X,bagnol_ns_rspec_pa.Y,'r');
legend('fe\_curve','cyberquake');

figure(3);view(2)

