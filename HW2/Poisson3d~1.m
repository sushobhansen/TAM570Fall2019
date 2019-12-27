clear all; close all; clc;

N = [2:2:10]';
times = zeros(size(N,1),1);

for i=1:size(N,1)
    for j=1:5
    tstart = tic;
    U = poisson3d(N(i));
    times(i) = times(i) + toc(tstart);
    end
end

semilogy(N,times/5,'b',N,1e-6*N.^3,'r');
xlabel('N'); ylabel('Time (sec)'); 
legend('Time','O(N^3)','Location','northwest');
axis square;
