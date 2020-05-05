function [surface] = irradiance( position,intensity, resolution)

    [x_surf y_surf] = meshgrid(0:resolution:50);
    
    surface = (position(3).^2.*intensity)./(((x_surf-position(1)).^2 + (y_surf-position(2)).^2 + position(3).^2)).^2;

end
