function [Sk] = subsampling_selection(M, alpha, rho)
%% Функция осуществляет выбор коэффициентов децимации/интерполяции для
%  неравнополосного косинусно-модулированного банка фильтров
% M     -- число каналов БФ.
% alpha -- параметр деформации частотной оси

f_arr=2*pi*(0:1/(2*M):0.5);     % Находим границы полос равнополосного БФ
ws=(1+rho)*pi/(2*M);

% Находим границы полос на деформированной оси
for i=1:M+1,
    w=f_arr(i);
    f_def(i)= w - 2*atan((alpha*sin(w))/(alpha*cos(w)-1) );
end
ws_def = (f_arr(2)+rho*pi/(2*M)) - 2*atan((alpha*sin(f_arr(2)+rho*pi/(2*M)))/(alpha*cos(f_arr(2)+rho*pi/(2*M))-1) );
f_def=f_def/(2*pi);            % перходим к нормированным частотам
fs_def = ws_def/(2*pi);
% Рассчитаем коэффициенты передискретизации так, чтобы уровень наложения
% спектров в банке фильтров определялся только ослаблением в полосе
% непропускания фильтра-прототипа.
f_arr=f_arr/(2*pi);             % перходим к нормированным частотам
Sk = ones(1,M);
for k=1:M
    Sk(k)=1;
    if (k==1)
        fU(k)=fs_def;
%         fU(k)=f_def(k+2);
        fL(k)=f_def(k);
    elseif (k==M)
            fU(k)=f_def(M+1);
            fL(k)=f_def(M-1);
    else
        fU(k)=f_def(k+2);
        fL(k)=f_def(k-1);
    end;
    B = fU(k)-fL(k);
%     nk_max = floor(1/(2*B));    % максимальное количество копий спектра
    nk_max = floor(fU(k)/(fU(k)-fL(k)));    % максимальное количество копий спектра
    for nk=1:nk_max,
        ub = floor(nk/(2*fU(k)));
        lb = max(ceil((nk-1)/(2*fL(k))),1);
        if (ub>=lb)
            for s=lb:ub,
                if (s*B<=0.5)
                    Sk(k)= max(Sk(k),s);
                end
            end
        end
    end
end
disp(Sk);
fprintf(1,'Total oversampling ratio (TOR) = %1.3f\n',sum(1./Sk));
end

