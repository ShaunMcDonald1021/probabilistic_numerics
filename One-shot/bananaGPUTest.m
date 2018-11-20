clear all; close all; clc;
gpuDevice(2);
num_iter = 15;
dim_list = [25 50 100 250 500];
num_dim = length(dim_list);
logCPUtime = zeros(num_iter, num_dim);
logGPUtime = zeros(num_iter, num_dim);

syms x z real
b = 0.5;
sigsq = 3;


for q = num_dim:-1:1
n = dim_list(q);

y = sym('y', [1 n]);
assumeAlso(y, 'real');

f(y) = mvnpdf([y(1) y(2)+b*(y(1)^2-sigsq) y(3:n)], zeros(1,n), diag([sigsq ones(1,n-1)]));

%Fortunately we already know what the gradient and Hessian are
%so in order to get this done we'll "cheat" a little and simply plug them in
%(not sure how much this will mirror the "real world" but oh well)
mle = [0 b*sigsq zeros(1,n-2)];
laplacemat = blkdiag(double(subs(hessian(log(f),y(1:2)), y(1:2), mle(1:2))), -eye(n-2));



epsilon = 1.5; %rule of thumb
sd = 1./sqrt(-diag(laplacemat));

m_t(y) = subs(f,y,mle)*exp(0.5*(y-mle)*laplacemat*(y-mle)');

f_func = matlabFunction(f, 'Vars', {[y]});
m_t_func = matlabFunction(m_t, 'Vars', {[y]});
disp('got everything set up except for the integral functions')
int_func_cell = cell(n,1);

%Again, knowing a bit about the structure of the banana allows us
%to cheat. Since all of the components of the L.A. are independent,
%there's really no point sitting around and waiting for integrals
%when we already know we just have to integrate one component

w = sym('y', [1 n-1]);
assumeAlso(w, 'real');
for i = 1:2
m0sym(y([1:(i-1) (i+1):n]),x,z) = int(m_t, y(i), x, z);
int_func_cell{i} = matlabFunction(m0sym, 'Vars', {[y([1:(i-1) (i+1):n]) x z]});
disp(strcat(num2str(n), ' dim, integral function ', num2str(i), ' done'))
end


syms rho real %running out of letters here
mle_red = mle([1:2, 4:end]);
laplacemat_red = laplacemat([1:2, 4:end], [1:2, 4:end]);
template_func(w,x,z) = subs(f,y,mle)*exp(0.5*(w-mle_red)*laplacemat_red*...
(w-mle_red)')*int(normpdf(rho,0,1),rho,x,z);

filename = strcat('intfunc',num2str(n),'.mat');

for i = 3:n
argcell = sym2cell([y([1:(i-1) (i+1):n]) x z]);

int_func_cell{i} = matlabFunction(template_func(argcell{:}), 'Vars', {[y([1:(i-1) (i+1):n]) x z]});
disp(strcat(num2str(n), ' dim, integral function ', num2str(i), ' done'))
save(filename,'int_func_cell')

end




QRcell = cell(n,1);
RRcell = cell(n,1);
QQcell = cell(n,1);
diffcell = cell(n,1);

disp('got EVERYTHING set up now')
for t = 1:num_iter
 

m0_arr = zeros(n, 1, 2*n-1);
diff_point = zeros(n,1);
diff_mat = zeros(3,1,2*n-1);

%CPU
tic
for i = 1:n
    lambda = 1.5*epsilon*sd(i);
    alpha = 1/(sd(i)*epsilon);
    s1 = mle(i) - epsilon*sd(i);
    s3 = mle(i) + epsilon*sd(i);
    s = [s1, mle(i), s3]';
     
    m0func = int_func_cell{i};
    
    cond_points = vertcat(repmat(mle([1:(i-1) (i+1):n]),n-1,1) - ...
        epsilon*diag(sd([1:(i-1) (i+1):n])), mle([1:(i-1) (i+1):n]), ...
        repmat(mle([1:(i-1) (i+1):n]),n-1,1) + ...
        epsilon*diag(sd([1:(i-1) (i+1):n])));
    
    for j = 1:3
        diff_point(i) = s(j);
        for k = 1:(2*n-1);
            diff_point([1:(i-1) (i+1):n]) = cond_points(k,:);
           % diff_pointcell = num2cell(diff_point);
            diff_mat(j,1,k) = f_func(diff_point') - m_t_func(diff_point');
        end
    end
    diffcell{i} = diff_mat;
    
    for k = 1:(2*n-1)
        %val_cell = num2cell([cond_points(k,:), s1, s3]);
        m0_arr(i,1,k) = m0func([cond_points(k,:), s1, s3]);
    end
    QRcell{i} = QR1d_se(s3,s,lambda,s1,s3);
    RRcell{i} = RR1d_se(s, s, lambda, s1, s3);
    QQcell{i} = QQ1d_se(s3, s3, lambda, s1, s3);
end

bias_coefs = blkdiag(QRcell{:})/blkdiag(RRcell{:});
variance = diag(blkdiag(QQcell{:}) - bias_coefs*blkdiag(QRcell{:})');
bias = zeros(n,1,2*n-1);
diffs = vertcat(diffcell{:});
for i = 1:(2*n-1)
    bias(:,1,i) = bias_coefs*diffs(:,:,i);
end

logCPUtime(t,q) = log(toc);

disp(strcat('CPU, ', num2str(n), ' dim, ', num2str(exp(logCPUtime(t,q))), ' s'))



%GPU
%Function evaluations are MUCH slower on the GPU, so we'll just push things there at the end
%for computations
tic
sdGPU = gpuArray(sd);
mleGPU = gpuArray(mle);
disp('GPU arrays set up')
for i = 1:n
    lambda = 1.5*1.5*sdGPU(i);
    alpha = 1/(sdGPU(i)*1.5);
    s1 = mleGPU(i) - 1.5*sdGPU(i);
    s3 = mleGPU(i) + 1.5*sdGPU(i);
    s = [s1, mleGPU(i), s3]';

    sCPU = [mle(i) - 1.5*sd(i) mle(i) mle(i) + 1.5*sd(i)]; 

    m0func = int_func_cell{i};
    
    cond_points = vertcat(repmat(mle([1:(i-1) (i+1):n]),n-1,1) - ...
        1.5*diag(sd([1:(i-1) (i+1):n])), mle([1:(i-1) (i+1):n]), ...
        repmat(mle([1:(i-1) (i+1):n]),n-1,1) + ...
        1.5*diag(sd([1:(i-1) (i+1):n])));
    
    for j = 1:3
        diff_point(i) = sCPU(j);
        for k = 1:(2*n-1);
   diff_point([1:(i-1) (i+1):n]) = cond_points(k,:);
           % diff_pointcell = num2cell(diff_point);
            diff_mat(j,1,k) = f_func(diff_point') - m_t_func(diff_point');

end
    end
    diffcell{i} = diff_mat;
    
    for k = 1:(2*n-1)
m0_arr(i,1,k) = m0func([cond_points(k,:), sCPU(1), sCPU(3)]);
    end
    QRcell{i} = QR1d_se(s3,s,lambda,s1,s3);
    RRcell{i} = RR1d_se(s, s, lambda, s1, s3);
    QQcell{i} = QQ1d_se(s3, s3, lambda, s1, s3);
%disp(strcat(num2str(i), num2str(n), ' GPU'))
end

bias_coefs = blkdiag(QRcell{:})/blkdiag(RRcell{:});
disp('bias_coefs done GPU')
variance = diag(blkdiag(QQcell{:}) - bias_coefs*blkdiag(QRcell{:})');
disp('variance done GPU')
diffs = gpuArray(vertcat(diffcell{:}));

%Looks like a for-loop may actually be faster than pagefun here??
%probably because pagefun wastes time making copies
bias = zeros(n,1,2*n-i, 'gpuArray');
for i = 1:(2*n-1)
    bias(:,1,i) = bias_coefs*diffs(:,:,i);
end
%bias = pagefun(@mtimes, bias_coefs, diffs);

disp('bias done GPU')

biasgat = gather(bias);
disp('gathered bias')
variancegat = gather(variance);
disp('gathered variance')
%m0gat = gather(m0_arrGPU);
%disp('gathered posterior mean')
%wait(gpuDevice)
logGPUtime(t,q) = log(toc);
disp(strcat('GPU, ', num2str(n), ' dim, ', num2str(exp(logGPUtime(t,q))), ' s'))
reset(gpuDevice)
save('times_everythingelse.mat', 'logCPUtime', 'logGPUtime')

end
end

quit
