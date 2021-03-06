clc;
clear;
close all;

%Retrieves Solar Spectrum and creates a matrix with the visible spectrum
filename='ASTM_SolarIrradiance_AM0.csv';
Solar_reference = csvread(filename,1,0,[1 0 1697 1]);

colorSpectrum = Solar_reference(Solar_reference(:,1)>0.3 & Solar_reference(:,1)<0.9,:);
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
colorIrradiance=trapz(colorSpectrum(:,1),colorSpectrum(:,2));
percentageOfVisibleLight= (colorIrradiance/total)*100;

fprintf('Total Irradiance %.2f\n',total);
fprintf('Visible Light Irradiance %.2f\n', colorIrradiance);
fprintf('Visible Light Irradiance percentage of total %.2f%%\n', percentageOfVisibleLight);

%Retrieves Photopic curve
filename='PhotopicCurve.csv';
photopic = csvread(filename,1,0,[1 0 471 1]);

%Retrieves IR data
% filename='QTH.csv';
% QTH = csvread(filename,1,0,[1 0 625 1]);
filename='QTH_v1.csv';
QTH = csvread(filename,1,0,[1 0 212 1]);

%creates an estimate for the Deep Blue LED using a guasing
watt=.660/2; %with current being near 0.25A
filename='GD_CSXPM1_14_20160712_spectrum.csv';
LED_B = csvread(filename,1,0,[1 0 66 1]);

B_intensity_multiplier = watt ./ trapz(LED_B(:,1), LED_B(:,2));
LED_B (:,2) = LED_B(:,2) .*B_intensity_multiplier; %should be in w/nm

%Getting Far red spectrum from data
watt = 0.383/2; % at 200mA
filename='far_red_spectrum.csv';
LED_R = csvread(filename,1,0,[1 0 136 1]);

red_intensity_multiplier = watt ./ trapz(LED_R(:,1), LED_R(:,2));
LED_R (:,2) = LED_R(:,2) .*red_intensity_multiplier; %should be in w/nm

%creates an estimate for the UV LED using a guasing
watt=.930/2; %with current being near 0.25A
LED_UV= guass_estimate(405,20);

UV_intensity_multiplier = watt ./ trapz(LED_UV(:,1), LED_UV(:,2));
LED_UV (:,2) = LED_UV(:,2) .*UV_intensity_multiplier; %should be in w/nm

%creates an estimate for the real blue using a guasing
lumen=45.7/2; %with current being near 0.20A
LED_BB= guass_estimate(475,20);

photopic_band = photopic(photopic(:,1)>LED_BB(1,1) & photopic(:,1)<LED_BB(end,1),:);

LED_BB_interp=interp1(LED_BB(:,1),LED_BB(:,2),photopic_band(:,1));
Bblue_intensity_multiplier = (lumen) ./ (683*trapz(photopic_band(:,1), LED_BB_interp.*photopic_band(:,2)));
LED_BB (:,2) = LED_BB(:,2) .*Bblue_intensity_multiplier; %should be in w/nm

%Retrieves Blue LED portion of curve and converts it to scaled intensity
cd = 95.23;
% filename='GW_CS8PM1_EM__blue__spectrum.csv';
% LED_W_B = csvread(filename,1,0,[1 0 66 1]);
filename='white_blue_v2_spectrum.csv';
LED_W_B =csvread(filename,1,0,[1 0 66 1]);

photopic_band = photopic(photopic(:,1)>LED_W_B(1,1) & photopic(:,1)<LED_W_B(end,1),:);

LED_W_interp=interp1(LED_W_B(:,1),LED_W_B(:,2),photopic_band(:,1));
blue_intensity_multiplier = (cd.*0.0426) ./ (683*trapz(photopic_band(:,1), LED_W_interp.*photopic_band(:,2)));
LED_W_B (:,2) = LED_W_B(:,2) .*blue_intensity_multiplier; %should be in w/nm

%Retrieves yellow LED portion of curve and converts it to scaled intensity
% filename='GW_CS8PM1_EM_yellow_spectrum.csv';
% LED_W_Y =csvread(filename,1,0,[1 0 135 1]);
filename='white_yellow_v2_spectrum.csv';
LED_W_Y =csvread(filename,1,0,[1 0 136 1]);

photopic_band = photopic (photopic(:,1)>LED_W_Y(1,1) & photopic(:,1)<LED_W_Y(end,1),:);

LED_W_interp=interp1(LED_W_Y(:,1),LED_W_Y(:,2),photopic_band(:,1));
yellow_intensity_multiplier = (cd.*0.9574) ./ (683*trapz(photopic_band(:,1), LED_W_interp.*photopic_band(:,2)));
LED_W_Y (:,2) = LED_W_Y(:,2) .*yellow_intensity_multiplier; %should be in w/nm

LED_W= [LED_W_B(1:end-1,:);LED_W_Y];
Estimation=trapz(LED_W(:,1),LED_W(:,2));
disp(Estimation);


%Retrieves Green LED portion of curve and converts it to scaled intensity
% filename='GREEN_LED.csv';
% LED_G = csvread(filename,1,0,[1 0 100 1]);
filename='true_green_spectrum.csv';
LED_G = csvread(filename,1,0,[1 0 150 1]);

lumen=143; %
photopic_band = photopic (photopic(:,1)>LED_G(1,1) & photopic(:,1)<LED_G(end,1),:);
LED_G_interp=interp1(LED_G(:,1),LED_G(:,2),photopic_band(:,1));
green_intensity_multiplier= (lumen) ./ (683*trapz(photopic_band(:,1), LED_G_interp.*photopic_band(:,2)));
LED_G (:,2) = (LED_G(:,2) .*green_intensity_multiplier); %should be in w/nm

N_W_LEDS=6; %number of white leds
N_G_LEDS=1; %number of green leds
N_UV_LEDS=1; %1 led will be taking up two spots
N_QTH=0.25;
N_B_LEDS=1;
N_BB_LEDS=1;
N_R_LEDS=1;

totalLED = combineSpectrum(LED_W,LED_UV,N_W_LEDS,N_UV_LEDS);
totalLED = combineSpectrum(totalLED,LED_G,1,N_G_LEDS);
totalLED = combineSpectrum(totalLED,LED_B,1,N_B_LEDS);
totalLED = combineSpectrum(totalLED,LED_BB,1,N_BB_LEDS);
totalLED = combineSpectrum(totalLED,LED_R,1,N_R_LEDS);
AREA=0.01; %Area the light is immited in m^2
totalIrradiance=totalLED;
totalIrradiance(:,2) = totalLED(:,2)./AREA;

totalIrradiance=combineSpectrum(totalIrradiance,QTH,1,N_QTH);

figure (1);
tiledlayout(2,1)

nexttile
plot(Solar_reference(:,1),Solar_reference(:,2));
hold on
area(colorSpectrum(:,1),colorSpectrum(:,2));
hold off
grid on;
axis([0 4 0 2300]);
title('Linear Plot of Spectral Irradiance @ AM0');
xlabel('Wavelength (\mum)'); 
ylabel('Spectral Irradiance (W/m^2 -\mum )');

nexttile
loglog(Solar_reference(:,1),Solar_reference(:,2));
hold on
area(colorSpectrum(:,1),colorSpectrum(:,2));
hold off
grid on;
title('log-log Plot of Spectral Irradiance @ AM0');
xlabel('Wavelength (\mum)'); 
ylabel('Spectral Irradiance (W/m^2 -\mum )');

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