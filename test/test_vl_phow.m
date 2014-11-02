load lena.mat

I = single(lena_gray);
[f,d] = vl_phow(I, 'step', 50, 'floatDescriptors', 1, 'verbose', 1);
dlmwrite('phow_f_matlab.txt', f', '\t');
dlmwrite('phow_d_matlab.txt', d', '\t');
