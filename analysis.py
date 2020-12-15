import numpy as np
import matplotlib.pyplot as plt
import sys

fileName = sys.argv[1]

f = open(fileName)

temperature     = 0.0
cutoff          = 0.0
sigma           = 0.0
temperature     = []
time            = []
epot            = []
ekin            = []
etot            = []
density         = 0.0
gr              = [[],[]]

line = f.readline()

while line.strip().split()[0] != "PRINTING":
    line = f.readline()+'.'

while line.strip().split()[0] != "Cutoff":
    line = f.readline()+'.'

cutoff = float(line.strip().split()[-3])

line = f.readline().strip().split()
line = f.readline().strip().split()
density = float(line[-2])

line = f.readline()
line = f.readline()
line = f.readline()
line = f.readline()
line = f.readline()
line = f.readline()
line = f.readline()
line = f.readline()
line = f.readline()
while line.strip().split()[0] != 'Computation':
    line = line.strip().split()
    if line[0] == "Equilibration":
        line = f.readline()
        continue
    time.append(float(line[1]))
    epot.append(float(line[2]))
    ekin.append(float(line[3]))
    etot.append(float(line[4]))
    temperature.append(float(line[5]))
    line = f.readline()

line = f.readline()
line = f.readline()
line = f.readline()
line = f.readline()
line = f.readline()
line = f.readline()
while line.strip().split()[0] != '='*150:
    line = line.strip().split()
    gr[0].append(float(line[0]))
    gr[1].append(float(line[1]))
    line = f.readline()


fig, (ax1,ax2) = plt.subplots(2,1)

ax1.plot(time,epot,label='Potential energy')
ax1.plot(time,ekin,label='Kinetic energy')
ax1.plot(time,etot,label='Total energy')
ax2.plot(time,temperature)
ax1.set_xlabel('Time (fs)')
ax1.legend()
ax2.set_xlabel('Time (fs)')
ax1.set_ylabel('Energy (Dal.bohr$^2$.fs$^{-2}$)')
ax2.set_ylabel('Temperature (K)')
ax1.set_title("Evolution of the energy over time")
ax2.set_title("Evolution of the temperature over time")

fig2, ax3 = plt.subplots(1,1)
ax3.plot(gr[0],gr[1])
ax3.set_title("Radial pair distribution function")
ax3.set_ylabel("g(r)")
ax3.set_xlabel("Distance (bohr)")


plt.show()
