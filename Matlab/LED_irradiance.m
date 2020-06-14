%Function LED_irradiance has 4 inputs and 3 outputs. The 4 inputs are
%position, intensity, length of test area, and apex angle. With these input
%the Irradiance model calcuates the irradiance on the designated surface.
%The ouput is an X,Y meshgrid and surface being the irradiance
function [x_surf y_surf surface] = LED_irradiance( position,intensity,length, apex)
    [x_surf y_surf] = meshgrid(0:length/10:length);
    M=-1.*(log(2)./log(cosd(apex/2)));
    position = position + [length./2 length./2 0];
    R=sqrt((x_surf-position(1)).^2 + (y_surf-position(2)).^2 + position(3).^2);
    surface = ((position(3).^(M+1)).*intensity)./(R.^(M+3));

end

