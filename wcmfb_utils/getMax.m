function [Mx,n] = getMax (E,w)

Mx(1) = E(1);
n(1) = 1;
i_mx = 2;
for i=2:length(E)-1,
    if (E(i-1)< E(i)) && (E(i+1)<E(i)),
        Mx(i_mx) = E(i);
        n(i_mx) = i;
        i_mx = i_mx+1;
    end
end
Mx(i_mx) = E(end);
n(i_mx) = length(E);

if (Mx(1)<Mx(2) && (Mx(1)<0)),
    Mx(1) = (Mx(2) - Mx(3))/( w(n(2)) - w(n(3)) ) * (w(n(1)) - w(n(2))) + Mx(2);
end

if Mx(i_mx)<0.1*Mx(i_mx-1),
    Mx(i_mx) = (Mx(i_mx-1) - Mx(i_mx-2))/( w(n(i_mx-1)) - w(n(i_mx-2)) ) * (w(n(i_mx)) - w(n(i_mx-1))) + Mx(i_mx-1);
end
Mx(i_mx)= abs(Mx(i_mx));

% check the condition V(wl) <= 0.01*min[V(wl-1), V(wl+1)]
for i=2:length(Mx)-1
    m = min(Mx(i-1), Mx(i+1))*0.01;
    if (m>Mx(i)) Mx(i)= m; end;
end