clear;
# БПФ
function F = DFT(func, M, hx)
  s = (M - length(func)) / 2;

  func = [zeros(s, 1).', func, zeros(s, 1).'];
  func = [func(1 + M/2:end), func(1:M / 2)];

  F = fft(func) * hx;
  F = [F(1 + M/2:end), F(1:M / 2)];
  F = F(s + 1:M - s);
end

# Аналитическое решение
function solve = AnaliticSolve(u, a)
  #e = exp(2i*a*pi.*u);
  #solve = (1i/pi).*(4*u - 2)./(4*u.^2 - 4*u - 3).*(e.^2-1)./e;
  solve = (1i/pi).*(4*u - 2)./(4*u.^2 - 4*u - 3).*(exp(2i*a*pi.*u)-exp(-2i*a*pi.*u));
end

# Аналитическое решение для двумерного случая
function solve = AnaliticSolve2D(u, v, a)
  #e = exp(2i*a*pi.*u);
  #solve = (1i/pi).*(4*u - 2)./(4*u.^2 - 4*u - 3).*(e.^2-1)./e;
  solve = (-1/(pi^2)).*(4*u - 2)./(4*u.^2 - 4*u - 3).*(exp(2i*a*pi.*u)-exp(-2i*a*pi.*u)).*(4*v - 2)./(4*v.^2 - 4*v - 3).*(exp(2i*a*pi.*v)-exp(-2i*a*pi.*v));
end

N = 64;
a = 5;
hx = 2 * a / N;
x = -a:hx:(a - hx / 2);

# Гауссов пучок
f = exp(-(x.^2));
#figure(1);
#plot(x, abs(f)); # Амплитуда
#figure(2);
#plot(x, arg(f)); # Фаза

# БПФ гауссова пучка
M = 256;
F = DFT(f, M, hx);
b = N ^ 2 / (4 * a * M);
hxi = 2 * b / N;
xi = -b:hxi:(b - hxi / 2);
#figure(1);
#plot(xi, abs(F), "b"); # Амплитуда
#figure(2);
#plot(xi, arg(F), "b"); # Фаза

# Стандартный метод численного интегрирования
[X, XI] = meshgrid(x, xi);
Kernel = exp(-2 * pi * 1i * X.*XI);
F2 = Kernel * f.' * hx;
#figure(1);
#plot(xi, abs(F2), "k"); # Амплитуда
#figure(2);
#plot(xi, arg(F2), "k"); # Фаза

# Сравнение БПФ и стандартного метода численного интегрирования
#figure(3);
#plot(xi, abs(F)); # Амплитуда
#hold on;
#plot(xi, abs(F2), "k"); # Амплитуда
#figure(4);
#plot(xi, arg(F)); # Фаза
#hold on;
#plot(xi, arg(F2), "k"); # Фаза

# Cветовое поле
LF = exp(-pi * 1i * x) + exp(3 * pi * 1i * x)
#figure(1);
#plot(x, abs(LF));  # Амплитуда
#figure(2);
#plot(x, arg(LF)); # Фаза

# Преобразование светового поля
FLF = DFT(LF, M, hx);
#figure(1);
#plot(xi, abs(FLF)); # Амплитуда
#figure(2);
#plot(xi, arg(FLF)); # Фаза

# Аналитическое решение
ASolve = AnaliticSolve(xi, a);
#figure(1);
#plot(xi, abs(FLF)); # Амплитуда
#hold on;
#plot(xi, abs(ASolve), "r"); # Фаза
#figure(2);
#plot(xi, arg(FLF)); # Амплитуда
#hold on;
#plot(xi, arg(ASolve), "r"); # Фаза

# Двумерный гауссов пучок
y = x.';
f = exp(-(x.^2)-(y.^2));
intervalA = [-a, a];
#figure(1);
#imagesc(intervalA, intervalA, abs(f)); # Амплитуда
#colorbar;
#figure(2);
#imagesc(intervalA, intervalA, arg(f)); # Фаза
#colorbar;

# Двумерное преобразование гауссова пучка
[X, Y] = meshgrid(x, y);
F2D = exp(-(X.^2)-(Y.^2));
intervalB = [-b, b];
for i = 1:rows(F2D)
  F2D(i, :) = DFT(F2D(i, :), M, hx);
endfor
for j = 1:columns(F2D)
  F2D(:, j) = DFT(F2D(:, j).', M, hx).';
endfor
#figure(1);
#imagesc(intervalB, intervalB, abs(F2D)); # Амплитуда
#colorbar;
#figure(2);
#imagesc(intervalB, intervalB, arg(F2D)); # Фаза
#colorbar;

# Двумерное световое поле
y = x.';
LF2D = (exp(-pi * 1i * x) + exp(3 * pi * 1i * x)).*(exp(-pi * 1i * y) + exp(3 * pi * 1i * y));
#figure(1);
#imagesc(intervalA, intervalA, abs(LF2D)); # Амплитуда
#colorbar;
#figure(2);
#imagesc(intervalA, intervalA, arg(LF2D)); # Фаза
#colorbar;

# Двумерное преобразование светового поля
[X, Y] = meshgrid(x, y);
FLF2D = (exp(-pi * 1i * x) + exp(3 * pi * 1i * x)).*(exp(-pi * 1i * y) + exp(3 * pi * 1i * y));
intervalB = [-b, b];
for i = 1:rows(FLF2D)
  FLF2D(i, :) = DFT(FLF2D(i, :), M, hx);
endfor
for j = 1:columns(FLF2D)
  FLF2D(:, j) = DFT(FLF2D(:, j).', M, hx).';
endfor
#figure(1);
#imagesc(intervalB, intervalB, abs(FLF2D)); # Амплитуда
#colorbar;
#figure(2);
#imagesc(intervalB, intervalB, arg(FLF2D)); # Фаза
#colorbar;

# Аналитическое решение двумерного случая
[X, Y] = meshgrid(x, y);
ASolve2D = AnaliticSolve2D(X, Y, a);
figure(1);
imagesc(intervalB, intervalB, abs(ASolve2D)); # Амплитуда
colorbar;
figure(2);
imagesc(intervalB, intervalB, arg(ASolve2D)); # Фаза
colorbar;


