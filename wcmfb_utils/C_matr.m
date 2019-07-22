function C = C_matr(w,N)
%% Function generate cosine matrix:
% C(w) = [2cos(w/2)  2cos(3w/2)  ...  2cos(|N-1|w/2)] for N even
% C(w) = [1  2cos(w)  ...  2cos(|N-1|w/2)] for N odd
N2=N/2;
C = zeros(1);
if (mod(N,2)==1),
    for n=1:N2-1,
        C(n+1,:) = [1 2*cos(w*n)];
    end
else
    n=(1:N2)'- 1/2;
    nw = n*w;
    C = 2*cos(nw);
end
