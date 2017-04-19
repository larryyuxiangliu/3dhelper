function cc = coeffBORto3D(c, Q)
% coeffBORto3D changes the coefficients c_jm ~ back to the 3D space, this
% allows us to evalute the fields more efficiently using double sum or FMM;
%
% cc = coeffBORto3D(c, Q)
%
% Inputs:   c: M-by-P matrix: stores the coefficients for each Fourier
%              mode;
%           Q: number of trap quadrature pts;
%
% Outputs: cc: MQ-by-1, stores the coefficients of each MFS source in the
%              3D space;
%
% Larry Liu, 07/30/2014

[N, P] = size(c);  % P is the number of Fourier modes; M is the number of rings
if Q>=P
    cQ = zeros(N, Q);
    cQ(:,1:P/2) = c(:, 1:P/2);
    cQ(:,end-P/2+1:end) = c(:,end-P/2+1:end);
else
    cQ = [c(:,1:Q/2), c(:,end-Q/2+1:end)];
end

cc = zeros(N*Q, 1);
for i=1:N
    cc((i-1)*Q+1:i*Q) = ifft(cQ(i,:)) * 2*pi;
end
