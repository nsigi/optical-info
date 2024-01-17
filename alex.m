clear;

function F = dft(f, M, hx)
  N = length(f);
  pad_size = (M - N) / 2;
  f = [zeros(pad_size, 1).', f, zeros(pad_size, 1).'];
  f = [f(M/2+1:end), f(1:M/2)];
  F = fft(f) * hx;
  F = [F(M/2+1:end), F(1:M/2)];
  F = F(pad_size + 1:M-pad_size);

end
N = 256;
R = 5;
hx = R / N;
x = 0:hx:(R - hx / 2);

function result = L(n, p, r)
    result = zeros(size(r));
    for j = 0:n
        result += power(-1, j) * r.^ j * nchoosek(n + p, n - j) / factorial(j);
    end
endfunction
function result = GL(n, p, r)
  p = abs(p);
  result = exp(r.^ 2 / -2) .* r.^ p .* L(n, p, r.^ 2);
endfunction

m = -3;
n = 3;
p = -2;

F = GL(n, p, x);


# График входной функции
figure(1);
subplot(2, 1, 1);
plot(x, abs(F));
title("Амплитуда");
subplot(2, 1, 2);
plot(x, arg(F));
title("Фаза");



# Восстанавливаем изображение в двумерный массив
image_input = zeros(2 * N, 2 * N);

for row = 1:rows(image_input)
  for col = 1:columns(image_input)
    alpha = round(sqrt((row - N)^2 + (col - N)^2)) + 1;
    if (alpha <= N)
      image_input(row, col) = F(alpha).* exp(1i * m * atan2(col - N, row - N));
    endif
  endfor
endfor

figure(2);
subplot(2, 1, 1);
imagesc([-R, R], [-R, R], abs(image_input));
colorbar;
title("2D Амплитуда");
subplot(2, 1, 2);
imagesc([-R, R], [-R, R], arg(image_input));
colorbar;
title("2D Фаза");


# Численная реализация преобразования Ханкеля
tic();

[X, XI] = meshgrid(x, x);
A = (2 * pi / 1i^ m)  * besselj(m, 2 * pi * X.* XI) .* X;
H = A * F.' * hx;

printf('Total cpu time: %f seconds\n', toc());

# График преобразования Ханкеля
figure(3);
subplot(2, 1, 1);
plot(x, abs(H));
title("Амплитуда Хенкеля");
subplot(2, 1, 2);
plot(x, arg(H));
title("Фаза Хенкеля");

# Восстанавливаем изображение в двумерный массив
image_output_hankel = zeros(2 * N - 1, 2 * N - 1);

for row = 1:rows(image_output_hankel)
  for col = 1:columns(image_output_hankel)
   alpha = round(sqrt((row - N)^2 + (col - N)^2)) + 1;
   if (alpha <= N)
    image_output_hankel(row, col) = H(alpha).* exp(1i * m * atan2(col - N, row - N));
   endif
  endfor
endfor

figure(4);
subplot(2, 1, 1);
imagesc([-R, R], [-R, R], abs(image_output_hankel));
colorbar;
title("2D Амплитуда Хенкеля");
subplot(2, 1, 2);
imagesc([-R, R], [-R, R], arg(image_output_hankel));
colorbar;
title("2D Фаза Хенкеля");

# Двумерное преобразованье Фурье через БПФ
M = 4096;

b = N ^ 2 / (4 * R * M);

tic();

for row = 1:rows(image_input)
    image_input(row, :) = dft(image_input(row, :), M, hx);
endfor

for col = 1:columns(image_input)
    image_input(:, col) = dft(image_input(:, col).', M, hx).';
endfor

printf('Total cpu time: %f seconds\n', toc());

figure(5);
subplot(2, 1, 1);
imagesc([-b, b], [-b, b], abs(image_input));
colorbar;
title("2D Амплитуда Фурье");
subplot(2, 1, 2);
imagesc([-b, b], [-b, b], arg(image_input));
colorbar;
title("2D Фаза Фурье");

