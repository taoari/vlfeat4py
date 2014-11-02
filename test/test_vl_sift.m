I = single(load('lena.txt'));
[f,d] = vl_sift(I, 'floatDescriptors', 'verbose');
dlmwrite('sift_f_matlab.txt', f, '\t');
dlmwrite('sift_d_matlab.txt', d, '\t');
