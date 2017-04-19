function y = imagshiftcurve(x,d)
% IMAGSHIFTCURVE - displace a set of points on curve in imag param direction
%
% y = imagshiftcurve(x,d)
% Inputs:
% x : list of N points in complex plane which uniformly sample a smooth
%     closed curve
% d : distance to move in imaginary direction (positive goes inside)
% Outputs:
% y : list of N displaced points in complex plane
%
% imagshiftcurve('test') runs tests
%
% To do: make maxgro an input option?
%
% Barnett 10/8/14, for MFS project with Larry Liu.
if strcmp(x,'test'), testimagshiftcurve, return; end
N = numel(x);
f = fft(x);
m = [0:N/2-1, 0, -N/2+1:-1];  % Fourier mode indices
m = m'; 
maxgro = 1e8;         % lose log10 of this many digits for large N limit
y = ifft(f.*min(maxgro,exp(-m*d)));  % applying a spectral filter

function testimagshiftcurve
d = 0.2;
for N=[30 100 1e3]   % testing large N crucial since maxgro cuts in
    s = (1:N)/N*2*pi;
    Z = @(s) exp(1i*s);  % circle...
    y = imagshiftcurve(Z(s),d); ye = Z(s+1i*d); norm(y-ye) % ans & exact pts
    Z = @(s) (1+0.3*cos(5*s)).*exp(1i*s); % starfish...
    y = imagshiftcurve(Z(s),d); ye = Z(s+1i*d); norm(y-ye)
end
figure;plot(Z(s),'k-');hold on;plot(ye,'b+');plot(y,'r.');axis equal;%keyboard
