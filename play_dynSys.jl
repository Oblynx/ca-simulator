using DynamicalSystems
using Plots
pyplot()

sys= Systems.lorenz96(4, rand(4); F=0.01);
trj= trajectory(sys,10);

plot3d(trj[:,1],trj[:,2],trj[:,3])
