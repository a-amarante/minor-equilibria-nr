c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c      EQUILIBRIUM.INC    (22 May 2020)
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c Author:
c
c Andre Amarante Luiz (A. Amarante) - a.amarante.br@gmail.com
c
c Parameters that you may want to alter at some point:
c
c novmax  = maximum number of vertices
c nopmax  = maximum number of faces
c noed = maximum number of vertices per face
c nocen= maximum number of centroids
c NMESS = maximum number of messages in message.in
c NEQ = maximum number of equilibrium points
c NBO = maximum number of boxes
c TINY  = a small number
c HUGE  = an implausibly large number
c
      integer novmax, nopmax, noed, nocen, NMESS, NEQ, NBO
      real*8 TINY,HUGE
c
      parameter (novmax = 500000)
      parameter (nopmax = 500000)
      parameter (noed = 10)
      parameter (nocen = 500000)
      parameter (NMESS = 60)
      parameter (NEQ = 10000)
      parameter (NBO = 10000)
      parameter (TINY = 4.D-15)
      parameter (HUGE = 9.9d29)
c
c------------------------------------------------------------------------------
c
c Constants:
c
c DR = conversion factor from degrees to radians
c gc = gravitational constant in km^3/s^2/kg
c
      real*8 pi,DR,TWOPI
c
      parameter (pi = 3.141592653589793d0)
      parameter (TWOPI = PI * 2.d0)
      parameter (DR = pi / 180.d0)
c      parameter (gc = 6.6719842229637d-20)
