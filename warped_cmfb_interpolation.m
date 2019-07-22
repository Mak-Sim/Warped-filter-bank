function [Y] = warped_cmfb_interpolation(wcmfb_params, X)
%warped_cmfb_decimation ��������� ������������ ��������� �������� �����
%��������

sb_pt = wcmfb_params.SubbandSamples;
Sk = wcmfb_params.SubsamplingFactors;   % ������������ ���������
M = wcmfb_params.NumberOfChannel;       % ����� ������� ����� ��������
Y = zeros(size(X));
for i=1:M
    tmp = Sk(i)*upsample(X(i,1:sb_pt(i)),Sk(i));
    len_tmp = length(tmp);
    len_tmp = min(len_tmp, size(Y,2));
    Y(i,1:len_tmp) = tmp(1:len_tmp);
end

end