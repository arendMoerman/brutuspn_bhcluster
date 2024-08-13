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

if __name__=="__main__":

  #fileName  = 'jpl_asteroids.csv'
  fileName  = 'jpl_comets.csv'

  #################################################################################

  data = read_data_float(fileName)

  #################################################################################


