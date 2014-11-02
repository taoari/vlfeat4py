I = single(load('lena.txt'));
Is = vl_imsmooth(I, 5.0);
dlmwrite('lena_smooth_matlab.txt', Is, '\t');
disp(size(Is));

