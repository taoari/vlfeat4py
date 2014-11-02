I = single(load('lena.txt'));
[f,d] = vl_dsift(I, 'step', 100, 'fast', 'norm', 'floatDescriptors', 'verbose');
dlmwrite('dsift_f_matlab.txt', f, '\t');
dlmwrite('dsift_d_matlab.txt', d, '\t');

tic
for i=1:10
	[f,d] = vl_sift(I, 'floatDescriptors');
end
toc
