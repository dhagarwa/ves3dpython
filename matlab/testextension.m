function [] = testextension()

theta = -2*pi:0.02*pi:2*pi;
x_mapfrom_theta = cos(theta);
plot(theta, fun(x_mapfrom_theta));


end


function val = fun(x)

    val = x;

end