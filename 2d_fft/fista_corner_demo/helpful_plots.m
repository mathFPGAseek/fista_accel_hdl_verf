% Helpful plots to understand corner layout coordinates to desing
% FPGA

figure(2); imagesc(x0); axis image off; title('True x moved to corner layout'); colormap gray;
figure(3); imagesc(blur_full); axis image off; title('blur full in corner layout'); colormap gray;
figure(4); imagesc(b); axis image off; title('Measurement is corner cropped'); colormap gray;