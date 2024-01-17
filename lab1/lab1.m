clear;
n = 1000;
m = 1000;
alpha = 1;
beta = 0.1;
a = 5;
p = 10;
hx = 2 * a/n;
x = -a:hx:(a - hx/2)
hxi = 2 * p/m;
xi = (-p:hxi:(p - hxi/2));
f = exp(1i * x * beta);
[X, XI] = meshgrid(x, xi);
Kernel = 1i * exp(-alpha * abs(X.-XI));
F = Kernel * f.' * hx;
figure(1);
plot(x, abs(f)); # Амплитуда входного сигнала
figure(2);
plot(x, arg(f)); # Фаза входного сигнала
figure(3);
plot(xi, abs(F));# Амплитуда выходного сигнала
figure(4);
plot(xi, arg(F));# Фаза выходного сигнала
