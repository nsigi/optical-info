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
# Обобщенный полином Лагерра
function poly = L(n, p, r)
    poly = zeros(size(r));
    for j = 0:n
        poly += power(-1, j) * nchoosek(n + p, n - j) * r.^ j / factorial(j);
    end
end
# Мода Гаусса-Лагерра
function moda = GL(n, p, r)
  p = abs(p);
  moda = exp(r.^ 2 / -2) .* r.^ p .* L(n, p, r.^ 2);
end

m = 2;
n = 2;
p = -3;

N = 256;
R = 5;
interval = [-R, R];
hx = R / N;
x = 0:hx:(R - hx / 2);

F = GL(n, p, x);

# Графики входной функции
figure(1);
plot(x, abs(F)); # Амплитуда
figure(2);
plot(x, arg(F)); # Фаза

# Восстанавление изображения в двумерный массив
rec_img = zeros(2 * N, 2 * N);

for i = 1:rows(rec_img)
  for j = 1:columns(rec_img)
    alpha = round(sqrt((i - N)^2 + (j - N)^2)) + 1;
    if (alpha <= N)
      rec_img(i, j) = F(alpha).* exp(1i * m * atan2(j - N, i - N));
    endif
  endfor
endfor

figure(3);
imagesc(interval, interval, abs(rec_img)); # Амплитуда
colorbar;
figure(4);
imagesc(interval, interval, arg(rec_img)); # Фаза
colorbar;


# Преобразования Ханкеля
tic();

[X, XI] = meshgrid(x, x);
A = (2 * pi / 1i^ m)  * besselj(m, 2 * pi * X.* XI) .* X;
H = A * F.' * hx;

printf('Time H: %f sec\n', toc());

# График преобразования Ханкеля
figure(5);
plot(x, abs(H)); # Амплитуда
figure(6);
plot(x, arg(H)); # Фаза

# Восстанавливаем изображение для Ханкеля
rec_img_H = zeros(2 * N - 1, 2 * N - 1);

for i = 1:rows(rec_img_H)
  for j = 1:columns(rec_img_H)
    alpha = round(sqrt((i - N)^2 + (j - N)^2)) + 1;
    if (alpha <= N)
      rec_img_H(i, j) = H(alpha).* exp(1i * m * atan2(j - N, i - N));
    endif
  endfor
endfor

figure(7);
imagesc(interval, interval, abs(rec_img_H)); # Амплитуда
colorbar;
figure(8);
imagesc(interval, interval, arg(rec_img_H)); # Фаза
colorbar;


# Двумерное преобразованье Фурье через БПФ
M = 2048;
b = N ^ 2 / (4 * R * M);

tic();
for i = 1:rows(rec_img)
  rec_img(i, :) = DFT(rec_img(i, :), M, hx);
endfor

for j = 1:columns(rec_img)
  rec_img(:, j) = DFT(rec_img(:, j).', M, hx).';
endfor

printf('Time F: %f sec\n', toc());

figure(9);
imagesc(interval, interval, abs(rec_img)); # Амплитуда
colorbar;
figure(10);
imagesc(interval, interval, arg(rec_img)); # Фаза
colorbar;
