function r = randgputest(n)
rand1 = randn(n,3*n);
rand2 = randn(3*n, 3*n);
rand3 = randn(3*n,1,2*n-1);
disp('ugh.')
tic
invtest = rand1*inv(rand2);
ult = zeros(n,1,2*n-1);
for i = 1:(2*n-1)
ult(:,1,i) = invtest*rand3(:,:,i);
end
disp(toc)

gpu1 = gpuArray(rand1);
gpu2 = gpuArray(rand2);
gpu3 = gpuArray(rand3);
tic
invg = gpu1*inv(gpu2);
ultg = pagefun(@mtimes, invg, gpu3);
wait(gpuDevice)
disp(toc)
end 

