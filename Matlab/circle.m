function h = circle(x,y,r)
    th = 0:pi/50:2*pi;
    for l = 1:length(x)
        xunit = r * cos(th) + x(l);
        yunit = r * sin(th) + y(l);
        h = plot(xunit, yunit);
        hold on;
    end
end