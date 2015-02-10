I = single(load('lena.txt'));
[f,d] = vl_dsift(I, 'step', 5, 'fast', 'norm', 'floatDescriptors');
numClusters = 30 ;
[means, covariances, priors] = vl_gmm(d, numClusters);
enc = vl_fisher(d, means, covariances, priors);
