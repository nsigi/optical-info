clear;
n = 1000; # точки на входе
a = 5; #border
hx = 2 * a/n; # step
x = -a:hx:(a - hx/2) #
m = 200; #точки на выходе
b = 5;
hxi = 2 * b/m;
xi = (-b:hxi:(b - hxi/2)).'; # транспонирование
f = sin(x);
# x - ветор  строка (суммирование по ней)
[X, XI] = meshgrid(xi, x); # произведение
Kernel = exp(-2*pi*1i*X.*XI);
F = Kernel * f * hx;
#plot(xi, arg(F));
plot(xi, abs(F));

