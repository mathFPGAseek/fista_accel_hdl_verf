% blur experiment
% zero out the bottom-right quadrant of x0
x_test = x0;
x_test(17:end,17:end) = 0;

blur_test = ifft2(H .* fft2(x_test));
figure(5); colormap gray;
imagesc(abs(blur_test(1:20,1:20)));
axis image off;
title('Experiment to show cropped area contains mixture')
