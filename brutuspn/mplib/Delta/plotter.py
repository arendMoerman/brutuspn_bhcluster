################################################################################

from exp1 import *
from tb_stats import *
from pylab import *

################################################################################

dirroot = "/host/data/exps/brutus/"
numroot = "brutus_N3/"

index_first = 0
index_last = 10000

eta = "0.015625" # 0.5, 0.25, 0.125, 0.0625, 0.03125, 0.015625

fileName_tcpu = "N3_tcpu"

################################################################################

root = dirroot + numroot

rootlf = root + "leapfrog2/"
roothe = root + "hermite4/"
rootbr = root + "brutus/"

################################################################################

header = False
toZip = False

tcpu_lf = []
i=index_first
while i<index_last:
  fileName = rootlf + "file_eta" + eta + "_" + str(i) + ".log"
  data = read_data(fileName, header, toZip)
  M = len(data)
  j=0
  while j<M:
    if len(data[j]) > 0:

      if data[j][0] == "t_cpu":
        tcpu_lf.append(float(data[j][2]))

    j += 1
  i += 1

tcpu_he = []
i=index_first
while i<index_last:
  fileName = roothe + "file_eta" + eta + "_" + str(i) + ".log"
  data = read_data(fileName, header, toZip)
  M = len(data)
  j=0
  while j<M:
    if len(data[j]) > 0:

      if data[j][0] == "t_cpu":
        tcpu_he.append(float(data[j][2]))

    j += 1
  i += 1

tcpu_br = []
i=index_first
while i<index_last:
  fileName = rootbr + "file_" + str(i) + ".log"
  data = read_data(fileName, header, toZip)
  M = len(data)
  j=0
  while j<M:
    if len(data[j]) > 0:

      if data[j][0] == "t_cpu":
        tcpu_br.append(float(data[j][2]))

    j += 1
  i += 1

print len(tcpu_lf), len(tcpu_he), len(tcpu_br)

################################################################################

u_lf, v_lf = get_cdf(tcpu_lf)
u_he, v_he = get_cdf(tcpu_he)
u_br, v_br = get_cdf(tcpu_br)

################################################################################

## plot tcpu
axis([1e-3, 1e4, 0.0, 1.0])
xlabel("t$_{cpu}$ [s]")
ylabel("f")
semilogx()
plot(u_lf, v_lf, color='red', linestyle=':', linewidth=3)
plot(u_he, v_he, color='blue', linestyle='-.', linewidth=3)
plot(u_br, v_br, color='black', linestyle='-', linewidth=3)

savefig(fileName_tcpu+'.png')
savefig(fileName_tcpu+'.pdf')




