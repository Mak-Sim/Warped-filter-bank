function G = my_zerophase(h,f,N)
% Zero-phase amplitude response G(f)
if (mod(N,2)==1)
    %% WORKING EXAMPLE
    N2=1+(N-1)/2;
    G=h(N2)*ones(length(f),1);      
    for m=1:length(f),
        w=2*pi*f(m);
        for k = 1:N2-1,
            G(m)=G(m)+2*h(N2-k)*cos(w*k);
        end
    end
else
    %% WORKING EXAMPLE
    N2=N/2;
    G=zeros(length(f),1); 
    for m=1:length(f),
        w=2*pi*f(m);
        for k = 1:N2,
            G(m)=G(m)+2*h(N2-k+1)*cos(w*(k-0.5));
        end
    end
end;