%%
% 
%  PREFORMATTED
%  TEXT
% 
%Spectral_Match_v1
%Calculates the Spectrum Match of LEDs and QTH based on datasheets relative
%intensity spectrum and radiant flux\luminous flux
%Simulates with least amount of components using White, UV, Green, red,
%and different waveband blues
clc;
clear;
close all;

%Retrieves ASTM Solar Reference Irradiance Spectrum
filename='ASTM_SolarIrradiance_AM0.csv';
Solar_reference = csvread(filename,1,0,[1 0 1697 1]);
total=trapz(Solar_reference(:,1),Solar_reference(:,2));
fprintf('Solar Reference Total Irradiance %.2f\n',total);

%Retrieves Photopic curve
filename='PhotopicCurve.csv';
photopic = csvread(filename,1,0,[1 0 471 1]);

%Retrieves QTH Irradiance Spectrum
filename='QTH_v1.csv';
QTH = csvread(filename,1,0,[1 0 212 1]);

%Retrieving OSRAM Deep Blue LED Relative Intensity Spectrum from CSV
watt=.660; %at 350mA
filename='GD_CSXPM1_14_20160712_spectrum.csv';
LED_B = csvread(filename,1,0,[1 0 66 1]);

%Calculates absolute Radiant flux Spectrum
B_intensity_multiplier = (watt ./ trapz(LED_B(:,1), LED_B(:,2)));
LED_B (:,2) = LED_B(:,2) .*B_intensity_multiplier;

%Retrieving OSRAM Far Red LED Relative Intensity Spectrum from CSV
watt = 0.383; % at 350mA
filename='far_red_spectrum.csv';
LED_R = csvread(filename,1,0,[1 0 136 1]);

%Calculates absolute Radiant flux Spectrum
red_intensity_multiplier = watt ./ trapz(LED_R(:,1), LED_R(:,2));
LED_R (:,2) = LED_R(:,2) .*red_intensity_multiplier; %should be in w/nm

%Creates relative intensity spectrum for UV LED from Gaussian curve
watt=.930*3/4; %with current being near 0.35A
LED_UV= guass_estimate(405,20);

%Calculates absolute Radiant flux Spectrum
UV_intensity_multiplier = watt ./ trapz(LED_UV(:,1), LED_UV(:,2));
LED_UV (:,2) = LED_UV(:,2) .*UV_intensity_multiplier; %should be in w/nm

%Creates relative intensity spectrum for Blue LED from Gaussian curve
lumen=45.7; %with current near 0.35A
LED_BB= guass_estimate(475,20);

%Retrieves photopic band that alignes with Blue LED band
photopic_band = photopic(photopic(:,1)>LED_BB(1,1) & photopic(:,1)<LED_BB(end,1),:);
LED_BB_interp=interp1(LED_BB(:,1),LED_BB(:,2),photopic_band(:,1)); % aligns data sets

%Calculates absolute Radiant flux Spectrum
Bblue_intensity_multiplier = (lumen) ./ (683*trapz(photopic_band(:,1), LED_BB_interp.*photopic_band(:,2)));
LED_BB (:,2) = LED_BB(:,2) .*Bblue_intensity_multiplier; %should be in w/nm

%Retrieves Blue spectrum portion of White LED and converts it to scaled intensity
lumen = 140;
filename='white_blue_v2_spectrum.csv';
LED_W_B =csvread(filename,1,0,[1 0 66 1]);

photopic_band = photopic(photopic(:,1)>LED_W_B(1,1) & photopic(:,1)<LED_W_B(end,1),:);

LED_W_interp=interp1(LED_W_B(:,1),LED_W_B(:,2),photopic_band(:,1));
blue_intensity_multiplier = (lumen.*0.0426) ./ (683*trapz(photopic_band(:,1), LED_W_interp.*photopic_band(:,2)));
LED_W_B (:,2) = LED_W_B(:,2) .*blue_intensity_multiplier; %should be in w/nm

%Retrieves yellow LED portion of white LED and converts it to scaled intensity
filename='white_yellow_v2_spectrum.csv';
LED_W_Y =csvread(filename,1,0,[1 0 136 1]);
photopic_band = photopic (photopic(:,1)>LED_W_Y(1,1) & photopic(:,1)<LED_W_Y(end,1),:);
LED_W_interp=interp1(LED_W_Y(:,1),LED_W_Y(:,2),photopic_band(:,1));
yellow_intensity_multiplier = (lumen.*0.9574) ./ (683*trapz(photopic_band(:,1), LED_W_interp.*photopic_band(:,2)));
LED_W_Y (:,2) = LED_W_Y(:,2) .*yellow_intensity_multiplier; %should be in w/nm
LED_W= [LED_W_B(1:end-1,:);LED_W_Y];

%Retrieves Green LED portion of curve and converts it to scaled intensity
filename='true_green_spectrum.csv';
LED_G = csvread(filename,1,0,[1 0 150 1]);
lumen=143;
photopic_band = photopic (photopic(:,1)>LED_G(1,1) & photopic(:,1)<LED_G(end,1),:);
LED_G_interp=interp1(LED_G(:,1),LED_G(:,2),photopic_band(:,1));
green_intensity_multiplier= (lumen) ./ (683*trapz(photopic_band(:,1), LED_G_interp.*photopic_band(:,2)));
LED_G (:,2) = (LED_G(:,2) .*green_intensity_multiplier); %should be in w/nm

%Number of Components
N_W_LEDS=9; %number of white LEDS
N_G_LEDS=4; %number of green LEDS
N_UV_LEDS=2; %number of UV LEDS
N_QTH=0.25; %Quartz tungsten lamp
N_B_LEDS=1; %number of Royal Blue LED
N_BB_LEDS=1;%number of Blue Leds
N_R_LEDS=3;%number of Red LEDs

%Add Absolute Spectrum Radiant Flux of all components
totalLED = combineSpectrum(LED_W,LED_UV,N_W_LEDS,N_UV_LEDS);
totalLED = combineSpectrum(totalLED,LED_G,1,N_G_LEDS);
totalLED = combineSpectrum(totalLED,LED_B,1,N_B_LEDS);
totalLED = combineSpectrum(totalLED,LED_BB,1,N_BB_LEDS);
totalLED = combineSpectrum(totalLED,LED_R,1,N_R_LEDS);

%Calculate the totalIrradiance assuming all radiant flux is emitted on test
%area
AREA=0.02; %Test area in m^2
totalIrradiance=totalLED;
totalIrradiance(:,2) = totalLED(:,2)./AREA;
totalIrradiance=combineSpectrum(totalIrradiance,QTH,1,N_QTH); %Add the irradiance of QTH

%Plots the ASTM Solar Reference Spectrum
figure (1);
tiledlayout(2,1)
nexttile
plot(Solar_reference(:,1),Solar_reference(:,2));
grid on;
axis([0 4 0 2300]);
title('Linear Plot of Spectral Irradiance @ AM0');
xlabel('Wavelength (\mum)'); 
ylabel('Spectral Irradiance (W/m^2 -\mum )');
nexttile
loglog(Solar_reference(:,1),Solar_reference(:,2));
grid on;
title('log-log Plot of Spectral Irradiance @ AM0');
xlabel('Wavelength (\mum)'); 
ylabel('Spectral Irradiance (W/m^2 -\mum )');

%Plots Irradiance Spectrum of LEDs
figure (2);
plot(Solar_reference(:,1).*1000,Solar_reference(:,2)./1000);
hold on
plot(LED_W(:,1),(LED_W(:,2).*N_W_LEDS)./AREA);
hold on
plot(LED_G(:,1),(LED_G(:,2).*N_G_LEDS)./AREA);
hold on
plot(LED_UV(:,1),(LED_UV(:,2).*N_UV_LEDS)./AREA);
hold on
plot(LED_BB(:,1),(LED_BB(:,2).*N_B_LEDS)./AREA);
hold on
plot(LED_B(:,1),(LED_B(:,2).*N_B_LEDS)./AREA);
hold on
plot(LED_R(:,1),(LED_R(:,2).*N_R_LEDS)./AREA);
hold on
plot(QTH(:,1),QTH(:,2).*N_QTH);
hold off

grid on;
axis([400 1100 0 2.5]);
title('Plot of Spectral Irradiance @ AM0 with Light Sources');
xlabel('Wavelength (nm)'); 
ylabel('Spectral Irradiance (W/m^2 /nm )');

%Creates graph with ASTM waveband sections color coded
figure(3);
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
fprintf('400-500 %.1f\n',((LED_sum)/total).*100);

hold on
section=totalIrradiance(totalIrradiance(:,1)>=500 & totalIrradiance(:,1)<=600,:);
LED_sum=trapz(section(:,1),section(:,2));
area(section(:,1),section(:,2));
fprintf('500-600 %.1f\n',((LED_sum)/total).*100);

hold on
section=totalIrradiance(totalIrradiance(:,1)>=600 & totalIrradiance(:,1)<=700,:);
LED_sum=trapz(section(:,1),section(:,2));
area(section(:,1),section(:,2));
fprintf('600-700 %.1f\n',((LED_sum)/total).*100);

hold on
section=totalIrradiance(totalIrradiance(:,1)>=700 & totalIrradiance(:,1)<=800,:);
LED_sum=trapz(section(:,1),section(:,2));
area(section(:,1),section(:,2));
fprintf('700-800 %.1f\n',((LED_sum)/total).*100);

IR_axis = 800:1:1500;
QTH_interp=interp1(totalIrradiance(:,1),totalIrradiance(:,2),IR_axis, 'spine');
IR = [IR_axis.' QTH_interp.'];
hold on
section=IR(IR(:,1)>=800 & IR(:,1)<=900,:);
LED_sum=trapz(section(:,1),section(:,2));
area(section(:,1),section(:,2));
fprintf('800-900 %.1f\n',((LED_sum)/total).*100);

hold on
section=IR(IR(:,1)>=900 & IR(:,1)<=1100,:);
LED_sum=trapz(section(:,1),section(:,2));
area(section(:,1),section(:,2));
fprintf('900-1100 %.1f\n',((LED_sum)/total).*100);

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