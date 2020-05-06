#unsteady diffusion equation using fvm
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
np.set_printoptions(linewidth=np.inf,precision=2)
from matplotlib import pyplot as plt
import time

m=100
n=100              #grid points 100*100
nps=m*n
dx=1.0/(m); dy=1.0/(n)
dt =1.0e-5                  #time_step
gamma=1.0                   #diffusion coefficient
t=0.0;  error=1.0;    errortime=0.0
volume=dx*dy
areae=areaw=distancey=dy
arean=areas=distancex=dx
ae=(gamma*areae)/distancex
aw=(gamma*areaw)/distancex
an=(gamma*arean)/distancey
asth =(gamma*areas)/distancey
ap0=volume/dt
anb=ae+aw+an+asth
ap=ap0
x = np.linspace(0, m+1,n+1)
y = np.linspace(0, m+1, n+1)
x_1, y_1 = np.meshgrid(x, y)

tem=np.zeros((n+1,m+1),dtype='double')
tem_old=np.zeros((n+1,m+1),dtype='double')
#print(tem_new)

#Initialization and Boundary conditions
tem[0,:]=1.0
#print(tem_new)

iterations=0; timestep=0
starttime=time.time()
#EXPLICIT METHOD
while True: #OUTER TIME LOOP
    tem_old[:,:]=tem[:,:]       #COPYING PRESENT TO PREVIOUS
    iterations=0
    while error>0.00001:    #INNER LOOP
        error=0.0
        for j in range(1,n):
            for i in range(1,m):
                temp=tem[j,i]
                tem[j,i]=((ap0-anb)*tem_old[j,i]+ae*tem[j,i+1]+aw*tem[j,i-1]+an*tem[j+1,i]+asth*tem[j-1,i]) /ap
                error=error+pow((tem[j,i]-temp),2.0)

        error=pow((error/nps),0.5)
        print('Iterations=',iterations,'  Error=',error)
        iterations+=1

    error=1.0
    errortime=0.0
    for j in range(1,n):
        for i in range(1,m):
            errortime=errortime+pow((tem[j,i]-tem_old[j,i]),2.0)
    errortime=pow((errortime/nps),0.5)
    t=t+dt
    timestep+=1
    print('Time=',t,'Timestep=',timestep,'ErrorTime=',errortime)

    # SAVING TEMP CONTOUR IMAGES at every 10 timestep FOR ANIMATION
    if timestep %20==0:
        #plt.contourf(x_1, y_1, tem, cmap='hsv') #CONTOUR PLOT
        plt.contour(x_1, y_1, tem, cmap='hsv')  #CONTOUR LINES
        # plt.colorbar()
        plt.title("UDE")
        # plt.xlabel("TimeStep= %d Iteration=%d"%(timestep,iterations))
        plt.xlabel("TimeStep= %d" % (timestep))
        name = str(timestep) + '.png'
        plt.savefig(name, format="png")

    if errortime<1e-5:break


#POTTING TEMPERATURE CONTURE OF FINAL TIMESTEP
#reversed_tem = tem[::-1]
#print(reversed_tem)
#print('maxtemp= ',tem_new.max())
#print('mintemp= ',tem_new.min())

'''fig = plt.figure(figsize = (11,7), dpi=100)
plt.contourf(x_1, y_1, reversed_tem, cmap='hsv')
plt.title("UDE")
plt.xlabel("TimeStep=",timestep,"Iteration=",iterations)
plt.colorbar()
plt.show()'''
endtime=time.time()-starttime
fig,ax=plt.subplots()
#ax=pyplot.subplot(111)
cs=ax.contourf(x_1, y_1,tem, cmap='hsv')
fig.colorbar(cs)
plt.title("UnsteadyDiffusionEquation")
plt.xlabel("EndTimeStep= %d total time= %f" % (timestep,endtime))
plt.savefig('UDE.png')
plt.show()
