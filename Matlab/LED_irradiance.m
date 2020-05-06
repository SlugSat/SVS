function [x_surf y_surf surface] = LED_irradiance( position,intensity,length, apex)

    [x_surf y_surf] = meshgrid(0:length/1000:length);
    M=-log(2)./log(cosd(apex/2));
    R=( (x_surf-position(1)).^2 + (y_surf-position(2)).^2 + position(3).^2);
    position = position + [length./2 length./2 0];
    surface = (position(3).^(M+1).*intensity)./(R).^(M+5);

end
