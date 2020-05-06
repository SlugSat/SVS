clc;
clear;
close all;

%Retrieves Solar Spectrum and creates a matrix with the visible spectrum
filename='ASTM_SolarIrradiance_AM0.csv';
Solar_reference = csvread(filename,1,0,[1 0 1697 1]);
total=trapz(Solar_reference(:,1),Solar_reference(:,2));
slot = Solar_reference(Solar_reference(:,1)>0.4 & Solar_reference(:,1)<0.5,:);
first = trapz(slot(:,1),slot(:,2));
slot = Solar_reference(Solar_reference(:,1)>0.5 & Solar_reference(:,1)<0.6,:);
second = trapz(slot(:,1),slot(:,2));
slot = Solar_reference(Solar_reference(:,1)>0.6 & Solar_reference(:,1)<0.7,:);
third = trapz(slot(:,1),slot(:,2));
slot = Solar_reference(Solar_reference(:,1)>0.7 & Solar_reference(:,1)<0.8,:);
fourth = trapz(slot(:,1),slot(:,2));
slot = Solar_reference(Solar_reference(:,1)>0.8 & Solar_reference(:,1)<0.9,:);
fifth = trapz(slot(:,1),slot(:,2));
slot = Solar_reference(Solar_reference(:,1)>0.9 & Solar_reference(:,1)<1.1,:);
six = trapz(slot(:,1),slot(:,2));

LENGTH = 0.01;
AREA = LENGTH.^2;
%figure (3);
x=0:LENGTH/100:LENGTH;
y=0:LENGTH/100:LENGTH;
DISTANCE =0.115;
[X Y LED_W] = LED_irradiance([0 0 DISTANCE],0.2835,LENGTH, 120);
W= trapz(y,trapz(x,LED_W,2));
W=W./AREA;

%figure (2);
[X Y LED_B] = LED_irradiance([0 0 DISTANCE],0.072188,LENGTH, 120);
B= trapz(y,trapz(x,LED_B,2));
B=B./AREA;

%figure(3);
[X Y LED_DP] = LED_irradiance([0 0 DISTANCE],0.21867,LENGTH, 120);
DP= trapz(y,trapz(x,LED_DP,2));
DP=DP./AREA;

%figure(4);
[X Y LED_G] = LED_irradiance([0 0 DISTANCE],0.089,LENGTH, 120);
G= trapz(y,trapz(x,LED_G,2));
G=G./AREA;

%figure(5);
[X Y LED_R] = LED_irradiance([0 0 DISTANCE],0.08976,LENGTH, 120);
R= trapz(y,trapz(x,LED_R,2));
R=R./AREA;

%figure(6);
[X Y LED_UV] = LED_irradiance([0 0 DISTANCE],0.1273,LENGTH, 120);
UV= trapz(y,trapz(x,LED_UV,2));
UV=UV./AREA;

[X Y LED_V] = LED_irradiance([0 0 DISTANCE],0.1089,LENGTH, 165);
V= trapz(y,trapz(x,LED_V,2));
V=V./AREA;

[X Y LED_C] = LED_irradiance([0 0 DISTANCE],0.1344,LENGTH, 120);
C= trapz(y,trapz(x,LED_C,2));
C=C./AREA;


%Retrieves Blue LED portion of curve and converts it to scaled intensity
filename='white_blue_v2_spectrum.csv';
LED_W_B =csvread(filename,1,0,[1 0 66 1]);

blue_intensity_multiplier = (W.*0.2368) ./ (trapz(LED_W_B(:,1), LED_W_B(:,2)));
LED_W_B (:,2) = LED_W_B(:,2) .*blue_intensity_multiplier; %should be in w/nm

%Retrieves yellow LED portion of curve and converts it to scaled intensity
filename='white_yellow_v2_spectrum.csv';
LED_W_Y =csvread(filename,1,0,[1 0 136 1]);

yellow_intensity_multiplier = (W.*0.7632) ./ (trapz(LED_W_Y(:,1), LED_W_Y(:,2)));
LED_W_Y (:,2) = LED_W_Y(:,2) .*yellow_intensity_multiplier; %should be in w/nm

LED_W= [LED_W_B(1:end-1,:);LED_W_Y];


%Reads Deep Blue from data
filename='GD_CSXPM1_14_20160712_spectrum.csv';
LED_B = csvread(filename,1,0,[1 0 66 1]);

B_intensity_multiplier = DP ./ trapz(LED_B(:,1), LED_B(:,2));
LED_B (:,2) = LED_B(:,2) .*B_intensity_multiplier; %should be in w/nm

%Getting Far red spectrum from data
filename='far_red_spectrum.csv';
LED_R = csvread(filename,1,0,[1 0 136 1]);

red_intensity_multiplier = R ./ trapz(LED_R(:,1), LED_R(:,2));
LED_R (:,2) = LED_R(:,2) .*red_intensity_multiplier; %should be in w/nm

%creates an estimate for the UV LED using a guasing
LED_UV= guass_estimate(405,20);
UV_intensity_multiplier = UV ./ trapz(LED_UV(:,1), LED_UV(:,2));
LED_UV (:,2) = LED_UV(:,2) .*UV_intensity_multiplier; %should be in w/nm

%creates an estimate for the Violet LED using a guasing
LED_V= guass_estimate(425,20);
V_intensity_multiplier = V ./ trapz(LED_V(:,1), LED_V(:,2));
LED_V (:,2) = LED_V(:,2) .*V_intensity_multiplier; %should be in w/nm

%creates an estimate for the Cyan LED using a guasing
LED_C= guass_estimate(490,25);
C_intensity_multiplier = C ./ trapz(LED_C(:,1), LED_C(:,2));
LED_C (:,2) = LED_C(:,2) .*C_intensity_multiplier; %should be in w/nm

%creates an estimate for the real blue using a guasing
LED_BB= guass_estimate(475,20);
Bblue_intensity_multiplier = DP ./ (trapz(LED_BB(:,1),LED_BB(:,2)));
LED_BB (:,2) = LED_BB(:,2) .*Bblue_intensity_multiplier; %should be in w/nm

%Retrieves Green LED portion of curve and converts it to scaled intensity
% filename='GREEN_LED.csv';
% LED_G = csvread(filename,1,0,[1 0 100 1]);
filename='true_green_spectrum.csv';
LED_G = csvread(filename,1,0,[1 0 150 1]);
green_intensity_multiplier= G ./ (trapz(LED_G(:,1),LED_G(:,2)));
LED_G (:,2) = (LED_G(:,2) .*green_intensity_multiplier); %should be in w/nm

%Retrieves IR data
filename='QTH_V1.csv';
QTH = csvread(filename,1,0,[1 0 211 1]);
QTH(:,2)= QTH(:,2)/4; %Converting mW values to W values and multiply by 25 due to emitting closer

N_W_LEDS=9; %number of white leds
N_G_LEDS=4; %number of green leds
N_UV_LEDS=4; 
N_QTH=1;
N_B_LEDS=1;
N_BB_LEDS=1;
N_R_LEDS=3;
N_V_LEDS=2;
N_C_LEDS=3;

totalIrradiance = combineSpectrum(LED_W,LED_UV,N_W_LEDS,N_UV_LEDS);
totalIrradiance = combineSpectrum(totalIrradiance,LED_G,1,N_G_LEDS);
totalIrradiance = combineSpectrum(totalIrradiance,LED_B,1,N_B_LEDS);
totalIrradiance = combineSpectrum(totalIrradiance,LED_BB,1,N_BB_LEDS);
totalIrradiance = combineSpectrum(totalIrradiance,LED_R,1,N_R_LEDS);
totalIrradiance = combineSpectrum(totalIrradiance,LED_V,1,N_V_LEDS);
totalIrradiance = combineSpectrum(totalIrradiance,LED_C,1,N_C_LEDS);
totalIrradiance=combineSpectrum(totalIrradiance,QTH,1,N_QTH);
figure (1);
plot(Solar_reference(:,1).*1000,(Solar_reference(:,2)./1000));
hold on
plot(LED_W(:,1),(LED_W(:,2).*N_W_LEDS));
hold on
plot(LED_G(:,1),(LED_G(:,2).*N_G_LEDS));
hold on
plot(LED_UV(:,1),(LED_UV(:,2).*N_UV_LEDS));
hold on
plot(LED_BB(:,1),(LED_BB(:,2).*N_BB_LEDS));
hold on
plot(LED_B(:,1),(LED_B(:,2).*N_B_LEDS));
hold on
plot(LED_R(:,1),(LED_R(:,2).*N_R_LEDS));
hold on
plot(LED_V(:,1),(LED_V(:,2).*N_V_LEDS));
hold on
plot(LED_C(:,1),(LED_C(:,2).*N_C_LEDS));
hold on
plot(QTH(:,1),QTH(:,2).*N_QTH);
hold off

legend({'Solar Reference','White','Green','UV','BB','B','R'})
grid on;
axis([400 1100 0 2.5]);
title('Plot of Spectral Irradiance @ AM0 with Light Sources');
xlabel('Wavelength (nm)'); 
ylabel('Spectral Irradiance (W/m^2 /nm )');

figure (2);
plot(Solar_reference(:,1).*1000,Solar_reference(:,2)./1000);

hold on
axis([400 1100 0 2.5]);
section=totalIrradiance(totalIrradiance(:,1)>=300 & totalIrradiance(:,1)<=400,:);
LED_sum=trapz(section(:,1),section(:,2));
area(section(:,1),section(:,2));
%fprintf('300-400 %.1f\n',((LED_sum)/f).*100);

hold on
section=totalIrradiance(totalIrradiance(:,1)>=400 & totalIrradiance(:,1)<=500,:);
LED_sum=trapz(section(:,1),section(:,2));
area(section(:,1),section(:,2));
fprintf('400-500 %.1f\n',((LED_sum)/first).*100);

hold on
section=totalIrradiance(totalIrradiance(:,1)>=500 & totalIrradiance(:,1)<=600,:);
LED_sum=trapz(section(:,1),section(:,2));
area(section(:,1),section(:,2));
fprintf('500-600 %.1f\n',((LED_sum)/second).*100);

hold on
section=totalIrradiance(totalIrradiance(:,1)>=600 & totalIrradiance(:,1)<=700,:);
LED_sum=trapz(section(:,1),section(:,2));
area(section(:,1),section(:,2));
fprintf('600-700 %.1f\n',((LED_sum)/third).*100);

hold on
section=totalIrradiance(totalIrradiance(:,1)>=700 & totalIrradiance(:,1)<=800,:);
LED_sum=trapz(section(:,1),section(:,2));
area(section(:,1),section(:,2));
fprintf('700-800 %.1f\n',((LED_sum)/fourth).*100);

IR_axis = 800:1:1500;
QTH_interp=interp1(totalIrradiance(:,1),totalIrradiance(:,2),IR_axis, 'spine');
IR = [IR_axis.' QTH_interp.'];
hold on
section=IR(IR(:,1)>=800 & IR(:,1)<=900,:);
LED_sum=trapz(section(:,1),section(:,2));
area(section(:,1),section(:,2));
fprintf('800-900 %.1f\n',((LED_sum)/fifth).*100);

hold on
section=IR(IR(:,1)>=900 & IR(:,1)<=1100,:);
LED_sum=trapz(section(:,1),section(:,2));
area(section(:,1),section(:,2));
fprintf('900-1100 %.1f\n',((LED_sum)/six).*100);

hold on
section=IR(IR(:,1)>=1100 & IR(:,1)<=1200,:);
LED_sum=trapz(section(:,1),section(:,2));
area(section(:,1),section(:,2));
%fprintf('1100-1300 %.1f\n',((LED_sum)/total).*100);

hold on
section=IR(IR(:,1)>=1200 & IR(:,1)<=1300,:);
LED_sum=trapz(section(:,1),section(:,2));
area(section(:,1),section(:,2));
%fprintf('1200-1300 %.1f\n',((LED_sum)/total).*100);

hold on
section=IR(IR(:,1)>=1300 & IR(:,1)<=1400,:);
LED_sum=trapz(section(:,1),section(:,2));
area(section(:,1),section(:,2));
%fprintf('1200-1300 %.1f\n',((LED_sum)/total).*100);

hold on
section=IR(IR(:,1)>=1400 & IR(:,1)<=1500,:);
LED_sum=trapz(section(:,1),section(:,2));
area(section(:,1),section(:,2));
%fprintf('1200-1300 %.1f\n',((LED_sum)/total).*100);

hold on
plot(Solar_reference(:,1).*1000,Solar_reference(:,2)./1000);

hold off
grid on;
title('Plot of AM0 Spectral Match Groups');
xlabel('Wavelength (nm)'); 
ylabel('Spectral Irradiance (W/m^2 /nm )');