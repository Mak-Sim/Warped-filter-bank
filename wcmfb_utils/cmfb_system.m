function xr=cmfb_system(h,M,x,a)
%% This function pass input signal x through analysis and than synthesis
% filters of cmfb. 
% h  - Coefficients of filter prototype
% M  - Number of channel
% x  - Input signal
% xr - Reconstructed signal
% a  - Coefficient of allpas transform

ap_a=[1 -a];        % Coefficients of  allpass filter
ap_b=[-a 1];        % Coefficients of  allpass filter
N = length (h);     % Order of FIR filter prototype
XMIN = -95;         % Parament for output plot lower bound -100 dB
poly_ord = N/(2*M);
Npt = length(x);
h = reshape(h,2*M,length(h)/(2*M));
h(:,2:2:size(h,2))=-h(:,2:2:size(h,2));     %G(z)=G(-z)
h = h(:);

X=zeros(length(h),length(x));     %Channel signals
X(1,:)=x;
tmp=x;
% 1. Allpass chain
if (a~=0),
    for i=2:N,
        tmp=filter(ap_b,ap_a,tmp);
        X(i,1:Npt) = tmp;       
    end
elseif (a==0)
    for i=2:N,
        X(i,1:Npt) = tmp;
        tmp = [0 tmp(1:end-1)];   
    end
end

% 2.Polyphase filtering
I_2M = eye(2*M);
polyphase_sum = sparse(kron(ones(1,poly_ord),I_2M));
for i=1:N,
    X(i,:)=h(i)*X(i,:);
end
Xk = polyphase_sum*X;

% 3. Cosime-modulation matrix
c = zeros(M,2*M);
for k =0:M-1,
    for l=0:2*M-1,
        c(k+1,l+1) = 2*cos( (2*k+1)*(pi/(2*M))*(l - (N-1)/2) + (-1)^(k)*(pi/4) );
    end;
end;

% 4. Cosine modulation 
X = c*Xk;

figure;
axes('Position',[0.12, 0.1, 0.85, 0.85]);
colors=[0.2 0.2 0.3; 0.2 0.5 0.2; 0 0.4 1; 1 0 0; 1 0 1; 0.1 0.6 0.1; 0 0 0; 0 0 0.5; 0 0.5 1; 0 1 0.5; 1 0.5 0; 0.5 0.25 1; 0.1 0.2 0.3; 0.3 0.1 0.2; 0.8 0.1 0.1; 0.1 0.5 0.2];
for k=1:M,
    Hk=X(k,:);
    plot((0:length(Hk)-1)/length(Hk),20*log10(abs(fft(Hk))),'LineWidth',2,'Color',colors(mod(k,16)+1,:) );
    ylabel('$|H_{k}(e^{j\omega})|,\mathrm{dB}$', 'Interpreter', 'Latex','FontSize',14);
    xlabel('$\omega/2\pi$', 'Interpreter', 'Latex','FontSize',14);
    hold on;
end; 
grid on;
xlim([0 0.5]);ylim([XMIN 5]);

%Decimation/Interpolation

% Synthesis
%6. Inverse modulation
X=c'*X;

%7. Polyphase filtering
Xk = repmat(flipud(X),poly_ord,1);
h = -h;
for i=1:N,
    Xk(i,:)=Xk(i,:)*h(i);
end;

% 8. Phase compensation
if (a~=0),
    for i=N:-1:2,
        Xk(i,:)=filter(ap_b,ap_a,Xk(i,:));
        Xk(i-1,:)=Xk(i,:)+Xk(i-1,:);
    end
    x=Xk(1,:);
else
    for i=2:N,
        for j=1:(i-1),
            Xk(i,:) =[0 Xk(i,1:size(Xk,2)-1)];
        end
    end
    x = sum(Xk);
end

xr = x;

X=fft(x);

figure;
plot((0:length(x)-1)./length(x),20*log10(abs(X)),'LineWidth',2,'Color',[0.6 0.1 0.6]);
grid on;xlim([0 0.5]);
ylabel('$|T_0(e^{j\omega})|,dB$', 'Interpreter', 'Latex','FontSize',14);
xlabel('$\omega/2\pi$', 'Interpreter', 'Latex','FontSize',14);

figure;
plot((0:length(x)-1)./length(x),unwrap(angle(X)),'LineWidth',2,'Color',[0.1 0.6 0.1]);
grid on;hold on;xlim([0 0.5]);
ylabel('$\varphi(\omega), \mathrm{deg}$', 'Interpreter', 'Latex','FontSize',14);
xlabel('$\omega/2\pi$', 'Interpreter', 'Latex','FontSize',14);

figure;
plot(0:length(xr)-1,xr,'LineWIdth',2,'color',[0 0 0]); xlim([0 length(xr)]);
xlabel('Samples');
ylabel('Normilized amplitude');
title('Filter bank response');
grid on;
set(gca,'LineWidth',2.5);
figure;
[gd,w]=grpdelay(x,1);ylabel('samples');xlim([0 0.5]);
plot(w/pi,gd,'LineWidth',2.5,'Color',[0.1 0.6 0.6]); 
grid on; xlim([0 0.5]);
ylabel('Samples','FontSize',14);
xlabel('$\omega / 2\pi$', 'Interpreter', 'Latex','FontSize',14);
end