function [ X ] = warped_cmfb_analysis(wcmfb_parms, x)
%warped_cmfb_analysis - реализуют неравнополосный банк фильтров анализа.
%   wcmfb_parms  --  структура, содержащая параметры банка фильтров
%   x            --  анализируемый сигнал.

% 1. Параметры
M = wcmfb_parms.NumberOfChannel;    % Число каналов банка фильтров
N = wcmfb_parms.OrderOfFP;          % Порядок фильтра-прототипа (ФП)
m = wcmfb_parms.OrderOfPPComponents;% Порядок полифазных компонент ФП
alpha = wcmfb_parms.Alpha;
h = wcmfb_parms.CoeffOfFP;          % Коэффициенты фильтра-прототипа
Npt = length(x);
X = zeros(M,Npt);                   % матрица канальный сигналов
allpass_delay_chain = zeros(N,1);

% 2. Коэффициенты обобщенной полифазной структуры
h = reshape(h,2*M,length(h)/(2*M));
h(:,2:2:size(h,2))=-h(:,2:2:size(h,2));     %G(z)=G(-z)
h = h(:);

% 3. Подготовка модуляционной матрицы
c = zeros(M,2*M);
for k =0:M-1,
    for l=0:2*M-1,
        c(k+1,l+1) = 2*cos( (2*k+1)*(pi/(2*M))*(l - (N-1)/2) + (-1)^(k)*(pi/4) );
    end;
end;
    
% 4. Реализация банка фильтров анализа
I_2M = eye(2*M);
polyphase_sum = sparse(kron(ones(1,m),I_2M));
for n=1:Npt
    % 4.1 Цепочка фазовых звеньев
    R1 = x(n);
    R2 = allpass_delay_chain(1);
    allpass_delay_chain(1) = R1;
    for i=2:N
        R3 = R2;
        R2 = allpass_delay_chain(i);
        R1 = (R2-R1)*alpha + R3;
        allpass_delay_chain(i) = R1; 
    end
    
    % 4.2 Умножение на коэффициенты фильтра-прототипа
    PP_out = allpass_delay_chain.*h;
    
    Xk = polyphase_sum*PP_out;

    % 4.3 Косинусная модуляция 
    X(:,n) = (c*Xk)';
end
end

