import math
import random
import sys

#################################################################################

def read_data(filename):
  data = [] 				# declare data container  

  try:
    infile = open(filename, "r")	# open file
  except:
    return data

  while True:
    line = infile.readline()
    if not line:
      break
    line = line.split(',')

    i=0
    while i<len(line):
      if line[i] == '' or line[i] == ' ':
        line[i] = '0'
      i += 1

    if line[0] != 'GM':
      data.append(line)  		# put data in container

  infile.close()			# close file

  data = zip(*data)			# inverse rows and columns

  data = [list(x) for x in data]        # make list type    
  return data				# return data

def read_data_float(filename):
  data = read_data(filename)

  N = len(data)
  i=0
  while i<N:
    M = len(data[i])	  
    j=0
    while j<M:
      try:
        data[i][j] = float(data[i][j])	    
      except:
        data[i][j] = 0.0
      j += 1
    i += 1

  return data

#################################################################################

def get_random_power_law_mass(m_max, m_index):
  while True:
    m_trial = random.random()*-10.
    fm = -m_trial*m_index
    p = random.random()*10.*m_index
    if p < fm:
      return 10.**m_trial*m_max

def get_power_law_mass(data, m_disk, m_index, m_max):
  G = 1.3274935e11 # km^3 / (Msun*s^2)

  # convert known mass values
  m = [x/G for x in data]

  # generate random masses for unknown mass values
  i=0
  while i<len(m):
    if m[i] == 0:
      m[i] = get_random_power_law_mass(m_max, m_index)
    i += 1

  return m

#################################################################################

def convert_mean_to_eccentric_binary(a, e, i, RA, Arp, A):
    delta = 1e-12
    EA = A
    diff = 1.0
    if(e > 0.4):
      while diff > delta:
        EA0 = EA
        f = 1.0 - e*math.cos(EA0)
        if f == 0.:
          EA = 1.
        else:
          EA = EA0 + (A + e*math.sin(EA0) - EA0) / (1.0 - e*math.cos(EA0))
          diff = abs(EA-EA0)
    else:
      while diff > delta:
        EA0 = EA
        EA = A + e*math.sin(EA0)
        diff = abs(EA-EA0)
    return EA

def convert_mean_to_eccentric_hyperbola(a, e, i, RA, Arp, A):
    delta = 1e-12
    if e == 1.:
      e += delta
    EA = A
    diff = 1.0
    while diff > delta:
      EA0 = EA
      f = e*math.cosh(EA0) - 1.
      if f == 0.:
        EA = 1.
      else:
        EA = EA0 + (A - e*math.sinh(EA0) + EA0) / (e*math.cosh(EA0) - 1.)
        diff = abs(EA-EA0)
    return EA

def convert_mean_to_eccentric_anomaly(a, e, i, RA, Arp, A):
    if a > 0. and e >= 1.:
      return -1
    elif a > 0. and e > 0:
      return convert_mean_to_eccentric_binary(a, e, i, RA, Arp, A)
    elif a > 0. and e <= 0.:
      e = 1e-12
      return convert_mean_to_eccentric_binary(a, e, i, RA, Arp, A)
    elif a <= 0. and e >= 1.:
      return convert_mean_to_eccentric_hyperbola(a, e, i, RA, Arp, A)
    elif a <= 0. and e > 0.:
      return -1
    else:
      return -1

def make_cartesian_binary(mu, a, e, i, RA, Arp, E):
    px = math.cos(Arp) * math.cos(RA) - math.sin(Arp) * math.cos(i) * math.sin(RA)
    py = math.cos(Arp) * math.sin(RA) + math.sin(Arp) * math.cos(i) * math.cos(RA)
    pz = math.sin(i) * math.sin(Arp)

    qx = -math.sin(Arp) * math.cos(RA) - math.cos(Arp) * math.cos(i) * math.sin(RA)
    qy = -math.sin(Arp) * math.sin(RA) + math.cos(Arp) * math.cos(i) * math.cos(RA)
    qz = math.sin(i) * math.cos(Arp)

    C1 = a*(math.cos(E)-e)
    C2 = a*math.sqrt(1.0-e*e)*math.sin(E)
    x = C1*px + C2*qx
    y = C1*py + C2*qy
    z = C1*pz + C2*qz  

    Edot = math.sqrt(mu/(a*a*a)) / (1.0-e*math.cos(E))
    C3 = -a*math.sin(E)*Edot
    C4 = a*math.sqrt(1.0-e*e)*math.cos(E)*Edot
    vx = C3*px + C4*qx
    vy = C3*py + C4*qy
    vz = C3*pz + C4*qz
    return x, y, z, vx, vy, vz

def make_cartesian_hyperbola(mu, a, e, i, RA, Arp, E):
    px = math.cos(Arp) * math.cos(RA) - math.sin(Arp) * math.cos(i) * math.sin(RA)
    py = math.cos(Arp) * math.sin(RA) + math.sin(Arp) * math.cos(i) * math.cos(RA)
    pz = math.sin(i) * math.sin(Arp)

    qx = -math.sin(Arp) * math.cos(RA) - math.cos(Arp) * math.cos(i) * math.sin(RA)
    qy = -math.sin(Arp) * math.sin(RA) + math.cos(Arp) * math.cos(i) * math.cos(RA)
    qz = math.sin(i) * math.cos(Arp)

    C1 = a*(math.cos(E)-e)
    C2 = a*math.sqrt(1.0-e*e)*math.sin(E)
    x = C1*px + C2*qx
    y = C1*py + C2*qy
    z = C1*pz + C2*qz  

    Edot = math.sqrt(mu/(a*a*a)) / (1.0-e*math.cos(E))
    C3 = -a*math.sin(E)*Edot
    C4 = a*math.sqrt(1.0-e*e)*math.cos(E)*Edot
    vx = C3*px + C4*qx
    vy = C3*py + C4*qy
    vz = C3*pz + C4*qz
    return x, y, z, vx, vy, vz

def get_cartesian_coordinates(mu, a, e, i, RA, Arp, E):
    if a > 0. and e >= 1.:
      return 0, 0, 0, 0, 0, 0
    elif a > 0. and e > 0:
      return make_cartesian_binary(mu, a, e, i, RA, Arp, E)
    elif a > 0. and e <= 0.:
      return make_cartesian_binary(mu, a, e, i, RA, Arp, E)
    elif a <= 0. and e >= 1.:
      return make_cartesian_hyperbola(mu, a, e, i, RA, Arp, E)
    elif a <= 0. and e > 0.:
      return 0, 0, 0, 0, 0, 0
    else:
      return 0, 0, 0, 0, 0, 0

#################################################################################

if __name__=="__main__":

  #fileName  = 'jpl_asteroids.csv'
  fileName  = 'jpl_comets.csv'

  m_disk  = 1e-9 # Msun
  m_index = 2.   # powerlaw index
  m_max = 1e-13  # maximum mass

  #################################################################################

  data = read_data_float(fileName)

  #################################################################################

  m = get_power_law_mass(data[0], m_disk, m_index, m_max)

  #################################################################################

  x = []
  y = []  
  z = []
  vx = []
  vy = []
  vz = []

  Nhyp = 0
  Nbin = 0
  Ne = 0
  Nflyby = 0
  Nflybin = 0
  Nunknown = 0

  q=0
  while q<len(m):
    print >> sys.stderr, q+1, '/', len(m)

    mu = 1.0 + m[q]
    a = data[2][q]
    e = data[3][q]
    i = data[4][q]
    RA = data[5][q]
    Arp = data[6][q]
    A = data[7][q]

    if a > 0. and e >= 1.:
      Nhyp += 1
    elif a > 0. and e > 0:
      Nbin += 1
    elif a > 0. and e <= 0.:
      Ne += 1
    elif a <= 0. and e >= 1.:
      Nflyby += 1
    elif a <= 0. and e > 0.:
      Nflybin += 1
    else:
      Nunknown += 1

    E = convert_mean_to_eccentric_anomaly(a, e, i, RA, Arp, A)
    if E != -1:
      myx, myy, myz, myvx, myvy, myvz = get_cartesian_coordinates(mu, a, e, i, RA, Arp, E)

      if myx != 0.:
        x.append(myx)
        y.append(myy)
        z.append(myz)
        vx.append(myvx)
        vy.append(myvy)
        vz.append(myvz)
    q += 1

  print >> sys.stderr, 'disk mass =', sum(m)
  print >> sys.stderr, 'category numbers', Nhyp, Nbin, Ne, Nflyby, Nflybin, Nunknown

  #################################################################################

  t = 0.0
  N = len(x)
  t_cpu = 0.0

  print '#', t, N, t_cpu
  
  i=0
  while i<N:

    print m[i], x[i], y[i], z[i], vx[i], vy[i], vz[i]

    i += 1

  #################################################################################


