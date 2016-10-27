% Use single-layer potential to evolve a front

dt = 1e-3;
% time step size

N = 128;
% number of points on front
theta = (0:N-1)'*2*pi/N;
modes = (-N/2:N/2-1)';

x = cos(theta);
y = sin(theta)/10;
z0 = x + 1i*y;
% initial shape

bgFlow_x = @(x,y) 4*ones(size(x));
bgFlow_y = @(x,y) 10*ones(size(x));

eps2 = 1e-1;
% regularization parameter

for k = 1:200
  z = x + 1i*y;
  z0 = z;
  zh = fftshift(fft(z));
  dzh = 1i*modes.*zh;
  dz = ifft(ifftshift(dzh));
  % derivative of shape
  sa = abs(dz);
  % arclength of shape

  nx = imag(dz)./sa;
  ny = -real(dz)./sa;
  % outward normal of shape


  vel_normal = zeros(N,1);

  for k = 1:N
    dist2 = (x(k) - x).^2 + (y(k) - y).^2;
    vel_normal(k) = sum(1./sqrt(dist2 + eps2).*sa)*2*pi/N;
  end
%  for k = 1:2:N
%    ind = (2:2:N);
%    dist2 = (x(k) - x(ind)).^2 + (y(k) - y(ind)).^2;
%    vel_normal(k) = sum(1./sqrt(dist2).*sa(ind))*2*pi/N;
%  end
%  for k = 2:2:N
%    ind = (1:2:N);
%    dist2 = (x(k) - x(ind)).^2 + (y(k) - y(ind)).^2;
%    vel_normal(k) = sum(1./sqrt(dist2).*sa(ind))*2*pi/N;
%  end
  vel_normal = vel_normal*2;
  figure(2)
  clf
  plot(vel_normal)
%  semilogy(modes,abs(fftshift(fft(vel_normal/N))))
%  ylim([1e-16 1e2])
  ylim([10 20])
%  pause(1e-1)
  pause

  x = x + dt*(vel_normal.*nx + bgFlow_x(x,y));
  y = y + dt*(vel_normal.*ny + bgFlow_y(x,y));

  figure(1);
  clf; 
%  subplot(2,1,1)
  plot(z0,'b-')
  hold on
  plot(x,y,'r-')
  axis equal;
  axis([-3 3 -1 3])
%  subplot(2,1,2)
%  semilogy(abs(sa-sa(1)))
%  pause
end
