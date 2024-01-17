n = 1000;
beta = 0.1;
m = 1000;
a = -5;
b = 5;
p = -10;
q = 10;
alpha = 1;
h = (b-a)/m;
x = linspace(a, b, n);
xi = linspace(p, q, m);
f1 = 0i;
f = exp(1i*x*beta);
for i = 1:m
  f1 = f1 + 1i * exp(-alpha*abs(x(i) - x))*exp(1i*x(i)*beta)*h;
 endfor
[X, XI] = meshgrid(x, xi);
figure(3);
plot(xi, abs(f1)); %%Амплитуда
figure(4);
plot(xi, arg(f1)); %%Фаза
figure(1);
plot(x, abs(f)); %%Амплитуда входного сигнала
figure(2);
plot(x, arg(f)); %%Фаза входного сигнала

