function [x_surf y_surf surface] = LED_irradiance( position,intensity, resolution)

    [x_surf y_surf] = meshgrid(0:resolution:0.05);
    position = position + [0.025 0.025 0];
    surface = (position(3).^2.*intensity)./(((x_surf-position(1)).^2 + (y_surf-position(2)).^2 + position(3).^2)).^2;

end
