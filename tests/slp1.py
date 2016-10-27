import numpy as np
from matplotlib import pyplot as plt

########################################################################
def bgFlow(x,y):
  # background flow
  u = 0*np.ones(len(x))
  v = 0*np.ones(len(x))
  return u,v

########################################################################
def main():

  T = 1 # time horizon
  m = 1000 # number of time steps
  dt = T/m # time step size
  lag = 4
  N = 30 # number of points on front
  print( "Number of points on front is %d" % N)

  dtheta = 2*np.pi/N
  # parameter spacing

  theta = np.linspace(0,2*np.pi,N,endpoint=False)
  # linspace for the praameter variable

  modes1 = np.arange(0,N/2+1,1)
  # positive Fourier modes
  modes2 = np.arange(-N/2+1,0,1)
  # negative Fourier modes
  modes = np.concatenate((modes1,modes2))
  # Fourier modes consistent with FFT output


  x = np.cos(theta)
  # x-coordinate of the front
  y = np.sin(theta)/5
  # y-coordinate of the front
  xs=np.zeros((lag+1,len(x)))
  ys = np.zeros((lag+1,len(y)))
  #Set initial history as starting values
  for i in range(lag+1):
      xs[i,:]=x.copy()*.9
      ys[i,:]=y.copy()*.9
  eps2 = 1e-1

  time = 0. # current time
  plt.ion()
  plt.plot(xs[0,:], ys[0,:], 'r')
  plt.show()
  for i in range(0,m):
    #Lag all data

    time = time + dt
    z = xs[0,:] + 1j*ys[0,:]
    # coordinates in complex[0,:] plane
    zh = np.fft.fft(z)
    # Fourier coefficients of the shape
    Dzh = zh*modes*1j
    # Fourier coefficients of the derivative of the shape
    Dz = np.fft.ifft(Dzh)
    # spectral derivative of the shape
    sa = abs(Dz)
    Dx = np.real(Dz)
    # x-coordinate of the derivative of the shape
    Dy = np.imag(Dz)
    # y-coordinate of the derivative of the shape
    speed = np.sqrt(Dx*Dx + Dy*Dy)
    # Jacobian or arclength of the curve

    nx = Dy/speed
    # x-coordinate of outward normal
    ny = -Dx/speed
    # y-coordinate of outward normal

    vel_normal = np.zeros(N)
    # normal velocity

    for k in range(0, N):
      dist2 = (xs[0,k] - xs[0,:]) ** 2 + (ys[0,k] - ys[0,:]) ** 2
      thickness = np.sqrt((xs[0,k]-xs[lag,k])**2+(ys[0,k]-ys[lag,k])**2)
      vel_normal[k] = 10*thickness*sum(1. / np.sqrt(dist2 + eps2) * sa) * 2 * np.pi / N
      print(xs[0,k],xs[lag,k])

    """
    for k in range (0,N,2):
      # loop over odd indexed points
      ind = np.arange(1,N,2)

      # only count contribution from even indexed points
      dist2 = np.square(x[k]-x[ind]) + np.square(y[k]-y[ind])

      vel_normal[k] = np.sum(np.divide(speed[ind],np.sqrt(dist2)))*2*np.pi/N
      print(np.sum(speed[ind])*2*np.pi/N)
      # normal velocity is the single-layer potential

    for k in range(1,N,2):
      # loop over even indexed points
      ind = np.arange(0,N,2)

      # only count contribution from even indexed points
      dist2 = np.square(x[k]-x[ind]) + np.square(y[k]-y[ind])
      print(np.sum(speed[ind]) * 2 * np.pi / N)
      vel_normal[k] = np.sum(np.divide(speed[ind],np.sqrt(dist2)))*2*np.pi/N
      # normal velocity is the single-layer potential
    """
    vel_normal = 2*vel_normal
    # multiply by 2 since we only looked at half of the indices

    velx,vely = bgFlow(xs[0,:],ys[0,:])
    # background velocity
    xs[0,:] = xs[0,:] + dt*(vel_normal*nx + velx)
    ys[0,:] = ys[0,:] + dt*(vel_normal*ny + vely)
    for j in range(0,lag):
        xs[lag-j,:]=xs[lag-(j+1),:]
        ys[lag-j,:]=ys[lag-(j+1),:]


    plt.clf()

    #plt.autoscale(False)
    plt.plot(xs[0,:],ys[0,:],'r')
    plt.plot(xs[lag, :], ys[lag, :], 'b')
    plt.draw()



    #plt.pause(0.0001)
    plt.quiver(xs[0,:],ys[0,:],(nx + velx),(ny + vely))
    plt.xlim([-1.5,1.5])
    plt.ylim([-.4, 1.])
    plt.waitforbuttonpress()




main()



