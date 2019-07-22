function [Nsf Lsf eq_coeff] = allpass_char(alpha,N)

ap_a=[1 -alpha];        % Coefficients of  allpass filter
ap_b=[-alpha 1];

x = [1 zeros(1,11000)];
tmp=x;
for i=2:N,
    tmp=filter(ap_b,ap_a,tmp);
end

[Offset Nsf] = IndFinder(tmp,0.1e-4);
eq_coeff = flipud(tmp(Offset:Nsf)');
Lsf = Nsf - Offset;
if nargout == 1
    Nsf = eq_coeff;
end
return

%% Function IndFinder()
function [ResInd1 ResInd2]=IndFinder(Arr,Threshold)
n=length(Arr);

for N=1:n
    if(abs(Arr(N))>Threshold)
        break
    end
end
ResInd1=N;

for N=n:-1:1
    if(abs(Arr(N))>Threshold)
        break
    end
end
ResInd2=N+1;
return
