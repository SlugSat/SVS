%Using Values from Spectral_Matchv1, the intensity and flux can be calculated
%interchangeably. This Script places LED lambertian light sources 118mm
%above test area of 2500mm^2.
clc;
clear;
close all;

%Assigning test area and distance from
LENGTH = 0.05;
AREA = LENGTH.^2;
x=0:LENGTH/10:LENGTH;
y=0:LENGTH/10:LENGTH;
DISTANCE =0.118;

LED_W_Positions = [ 
    0 0 DISTANCE; 
    0.035 0 DISTANCE;
    0.035 0.035 DISTANCE; 
    0.035 -0.035 DISTANCE;
    -0.035 0.035 DISTANCE;
    -0.035 -0.035 DISTANCE;
    -0.035 0 DISTANCE;
    0 -0.035 DISTANCE;
    0 0.035 DISTANCE
    ];
sum = zeros(length(x),length(x));
for i = 1:9
    [X Y LED] = LED_irradiance(LED_W_Positions(i,:),0.2835,LENGTH, 120);
    sum =sum+LED;
end
W= trapz(y,trapz(x,sum,2));
W=W./AREA;
% figure (1);
% nexttile;
% surf(X,Y,sum);
% view(2);
% colorbar
% xlabel('Cubesat X Coordinate m');
% ylabel('Cubesat Y Coordinate m');
% zlabel('Irradiance W/m^2');
% zlim([180 200]);

% figure(2);
LED_DP_Positions = [ 
     -0.00551 0 DISTANCE
];
sumBuffer = zeros(length(x),length(x));
for i = 1:1
    [X Y LED] = LED_irradiance(LED_DP_Positions(i,:),0.21867,LENGTH, 120);
    sumBuffer =sumBuffer+LED;
end
DP= trapz(y,trapz(x,sumBuffer,2));
DP=DP./AREA;
% surf(X,Y,sumBuffer);
% colorbar
% view(2);
% xlabel('Cubesat X Coordinate m');
% ylabel('Cubesat Y Coordinate m');
% zlabel('Irradiance W/m^2');

sum = sum+sumBuffer;

% figure(3);
LED_G_Positions = [ 
    0.02183 0.01863 DISTANCE;
    0.02183 -0.01863 DISTANCE;
    -0.02183 0.01863 DISTANCE;
    -0.02183 -0.01863 DISTANCE;
];
sumBuffer = zeros(length(x),length(x));
for i = 1:4
    [X Y LED] = LED_irradiance(LED_G_Positions(i,:),0.089,LENGTH, 120);
    sumBuffer =sumBuffer+LED;
end
G= trapz(y,trapz(x,sumBuffer,2));
G=G./AREA;
% surf(X,Y,sumBuffer);
% colorbar
% view(2);
% xlabel('Cubesat X Coordinate m');
% ylabel('Cubesat Y Coordinate m');
% zlabel('Irradiance W/m^2');

sum = sum+sumBuffer;

% figure(4);
LED_R_Positions = [ 
    0.035 -0.01398 DISTANCE;
    -0.01207 -0.035 DISTANCE;
    -0.035 0.02356 DISTANCE;
];
sumBuffer = zeros(length(x),length(x));
for i = 1:3
    [X Y LED] = LED_irradiance(LED_R_Positions(i,:),0.08976,LENGTH, 120);
    sumBuffer =sumBuffer+LED;
end
R= trapz(y,trapz(x,sumBuffer,2));
R=R./AREA;
% surf(X,Y,sumBuffer);
% colorbar
% view(2);
% xlabel('Cubesat X Coordinate m');
% ylabel('Cubesat Y Coordinate m');
% zlabel('Irradiance W/m^2');

sum = sum+sumBuffer;

% figure(5);
LED_UV_Positions = [
    0.01313 0.01313 DISTANCE;
    0.01313 -0.01313 DISTANCE;
    -0.01313 0.01313 DISTANCE;
    -0.01313 -0.01313 DISTANCE;
];
sumBuffer = zeros(length(x),length(x));
for i = 1:3
    [X Y LED] = LED_irradiance(LED_UV_Positions(i,:),0.1273,LENGTH, 120);
    sumBuffer =sumBuffer+LED;
end
UV= trapz(y,trapz(x,sumBuffer,2));
UV=UV./AREA;
% surf(X,Y,sumBuffer);
% colorbar
% view(2);
% xlabel('Cubesat X Coordinate m');
% ylabel('Cubesat Y Coordinate m');
% zlabel('Irradiance W/m^2');

sum=sum+sumBuffer;

% figure(6);
LED_V_Positions = [ 
    -0.01313 0.02123 DISTANCE;
    0.01313 -0.02123 DISTANCE;
];
sumBuffer = zeros(length(x),length(x));
for i = 1:2
    [X Y LED] = LED_irradiance(LED_V_Positions(i,:),0.1089,LENGTH, 165);
    sumBuffer =sumBuffer+LED;
end
V= trapz(y,trapz(x,sumBuffer,2));
V=V./AREA;
% surf(X,Y,sumBuffer);
% colorbar
% view(2);
% xlabel('Cubesat X Coordinate m');
% ylabel('Cubesat Y Coordinate m');
% zlabel('Irradiance W/m^2');

sum =sum+sumBuffer;

% figure(7);
LED_C_Positions = [ 
    -0.035 -0.00551 DISTANCE;
    0.00551 0.035 DISTANCE;
    0.035 -0.029 DISTANCE;
];
sumBuffer = zeros(length(x),length(x));
for i = 1:3
    [X Y LED] = LED_irradiance(LED_C_Positions(i,:),0.1344,LENGTH, 120);
    sumBuffer =sumBuffer+LED;
end
C= trapz(y,trapz(x,sumBuffer,2));
C=C./AREA;
% surf(X,Y,sumBuffer);
% colorbar
% view(2);
% xlabel('Cubesat X Coordinate m');
% ylabel('Cubesat Y Coordinate m');
% zlabel('Irradiance W/m^2');

sum=sum+sumBuffer;

% figure(8);
LED_B_Positions = [ 
    %%
    % 
    %  PREFORMATTED
    %  TEXT
    % 
    0.00551 0 DISTANCE;

];
sumBuffer = zeros(length(x),length(x));
for i = 1:1
    [X Y LED] = LED_irradiance(LED_B_Positions(i,:),0.072188,LENGTH, 120);
    sumBuffer =sumBuffer+LED;
end
B= trapz(y,trapz(x,sumBuffer,2));
B=B./AREA;
% surf(X,Y,sumBuffer);
% colorbar
% view(2);
% xlabel('Cubesat X Coordinate m');
% ylabel('Cubesat Y Coordinate m');
% zlabel('Irradiance W/m^2');

sum=sum+sumBuffer;

figure(9);
surf(X,Y,sum);
cbar =colorbar;
cbar.Label.String = 'Irradiance W/m^2';

view(2);
xlabel('Cubesat X Coordinate (m)');
ylabel('Cubesat Y Coordinate (m)');
zlabel('Irradiance W/m^2');
title('Irradiance on CubeSat Surface with 0.118m Distance');

max = max(max(sum));
min=min(min(sum));
fprintf('The Non-Uniformity is %.1f%%\n',((max-min)./(max+min)).*100);

figure (10);
plot(LED_C_Positions(:,1),LED_C_Positions(:,2),'cX','Markersize',20);
hold on
plot(LED_DP_Positions(:,1),LED_DP_Positions(:,2),'X','Markersize',20);
hold on
plot(LED_G_Positions(:,1),LED_G_Positions(:,2),'gX','Markersize',20);
hold on
plot(LED_R_Positions(:,1),LED_R_Positions(:,2),'rX','Markersize',20);
hold on
plot(LED_UV_Positions(:,1),LED_UV_Positions(:,2),'X','Markersize',20);
hold on
plot(LED_V_Positions(:,1),LED_V_Positions(:,2),'mX','Markersize',20);
hold on
plot(LED_W_Positions(:,1),LED_W_Positions(:,2),'kX','Markersize',20);
hold on
plot(LED_B_Positions(:,1),LED_B_Positions(:,2),'bX','Markersize',20);
%hold on
%circle([LED_C_Positions(:,1);LED_DP_Positions(:,1);LED_G_Positions(:,1);LED_R_Positions(:,1);LED_UV_Positions(:,1);LED_V_Positions(:,1);LED_W_Positions(:,1);LED_B_Positions(:,1)],[LED_C_Positions(:,2);LED_DP_Positions(:,2);LED_G_Positions(:,2);LED_R_Positions(:,2);LED_UV_Positions(:,2);LED_V_Positions(:,2);LED_W_Positions(:,2);LED_B_Positions(:,2)], 0.0198);
hold off
axis equal;
axis([-0.07/2 0.07/2 -0.07/2 0.07/2])
legend('Cyan','Deep Blue', 'Green', 'Red', 'UV','Violet','White','Blue');
title('LED Bulb Placement on Panel');
xlabel('X position (m)'); 
ylabel('Y position (m)'); 

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

%creates an estimate for the royal blue using a guasing
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

totalIrradiance = combineSpectrum(LED_W,LED_UV,1,1);
totalIrradiance = combineSpectrum(totalIrradiance,LED_G,1,1);
totalIrradiance = combineSpectrum(totalIrradiance,LED_B,1,1);
totalIrradiance = combineSpectrum(totalIrradiance,LED_BB,1,1);
totalIrradiance = combineSpectrum(totalIrradiance,LED_R,1,1);
totalIrradiance = combineSpectrum(totalIrradiance,LED_V,1,1);
totalIrradiance = combineSpectrum(totalIrradiance,LED_C,1,1);
totalIrradiance(:,2)=totalIrradiance(:,2).*1.5;
totalIrradiance=combineSpectrum(totalIrradiance,QTH,1,1);
figure (1);
plot(Solar_reference(:,1).*1000,(Solar_reference(:,2)./1000));
hold on
plot(LED_W(:,1),(LED_W(:,2)));
hold on
plot(LED_G(:,1),(LED_G(:,2)));
hold on
plot(LED_UV(:,1),(LED_UV(:,2)));
hold on
plot(LED_BB(:,1),(LED_BB(:,2)));
hold on
plot(LED_B(:,1),(LED_B(:,2)));
hold on
plot(LED_R(:,1),(LED_R(:,2)));
hold on
plot(LED_V(:,1),(LED_V(:,2)));
hold on
plot(LED_C(:,1),(LED_C(:,2)));
hold on
plot(QTH(:,1),QTH(:,2));
hold off

legend({'ASTM Solar Reference','White','Green','UV','Deep Blue','Blue','Red','Violet','Cyan','QTH'})
grid on;
axis([400 1100 0 2.5]);
title('Spectral Irradiance of Each component');
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
legend({'ASTM Solar Reference'});
hold off
grid on;
title('Plot of AM0 Spectral Match Groups');
xlabel('Wavelength (nm)'); 
ylabel('Spectral Irradiance (W/m^2 /nm )');
