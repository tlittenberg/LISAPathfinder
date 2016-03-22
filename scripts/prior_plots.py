#!/usr/bin/env python
import numpy as np
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
import math
from matplotlib.backends.backend_pdf import PdfPages

C1_X = -9.260000e-01 # x coordinate of spacecraft bottom deck corner 1 [m]
C1_Y = -2.168000e-01 # y coordinate of spacecraft bottom deck corner 1 [m]
C2_X = -9.260000e-01 # x coordinate of spacecraft bottom deck corner 2 [m]
C2_Y = 2.048000e-01  # y coordinate of spacecraft bottom deck corner 2 [m]
C3_X = -5.263000e-01 # x coordinate of spacecraft bottom deck corner 3 [m]
C3_Y = 8.970000e-01  # y coordinate of spacecraft bottom deck corner 3 [m]
C4_X = 5.163000e-01  # x coordinate of spacecraft bottom deck corner 4 [m]
C4_Y = 8.970000e-01  # y coordinate of spacecraft bottom deck corner 4 [m]
C5_X = 9.160000e-01  # x coordinate of spacecraft bottom deck corner 5 [m]
C5_Y = 2.048000e-01  # y coordinate of spacecraft bottom deck corner 5 [m]
C6_X = 9.160000e-01  # x coordinate of spacecraft bottom deck corner 6 [m]
C6_Y = -2.168000e-01 # y coordinate of spacecraft bottom deck corner 6 [m]
C7_X = 5.163000e-01  # x coordinate of spacecraft bottom deck corner 7 [m]
C7_Y = -9.090000e-01 # y coordinate of spacecraft bottom deck corner 7 [m]
C8_X = -5.263000e-01 # x coordinate of spacecraft bottom deck corner 8 [m]
C8_Y = -9.090000e-01 # y coordinate of spacecraft bottom deck corner 8 [m]
H = 8.315000e-01

def area(phi_off):
    phi_off=phi_off
    x=np.zeros(8)
    x[0]=math.sin(math.atan2(C1_Y,C1_X)+phi_off);
    x[1]=math.sin(math.atan2(C2_Y,C2_X)+phi_off);
    x[2]=math.sin(math.atan2(C3_Y,C3_X)+phi_off);
    x[3]=math.sin(math.atan2(C4_Y,C4_X)+phi_off);
    x[4]=math.sin(math.atan2(C5_Y,C5_X)+phi_off);
    x[5]=math.sin(math.atan2(C6_Y,C6_X)+phi_off);
    x[6]=math.sin(math.atan2(C7_Y,C7_X)+phi_off);
    x[7]=math.sin(math.atan2(C8_Y,C8_X)+phi_off);
    w=max(x)-min(x)
    return w*H

def area2(phi_off):
    sum=0
    for i in range(8):
        cth=math.cos(phi_norm[i]-phi_off)
        if(cth>0):
            sum+=cth*fa[i]
    return sum

def area_3d(theta): 
    return mean_area*math.sin(theta)+top_area*abs(math.cos(theta))

def triangle_area(x0,y0,x1,y1,x2,y2):
    return abs((y2-y0)*(x1-x0)-(x2-x0)*(y1-y0))/2.0

def side_area(x0,y0,x1,y1):
    return math.sqrt((y1-y0)**2+(x1-x0)**2)*H

def side_normal(x0,y0,x1,y1):
    return np.array([y1-y0,-(x1-x0)])/math.sqrt((y1-y0)**2+(x1-x0)**2)

top_area=0
top_area+=triangle_area(C1_X,C1_Y,C2_X,C2_Y,C3_X,C3_Y)
top_area+=triangle_area(C1_X,C1_Y,C4_X,C4_Y,C3_X,C3_Y)
top_area+=triangle_area(C1_X,C1_Y,C4_X,C4_Y,C5_X,C5_Y)
top_area+=triangle_area(C1_X,C1_Y,C6_X,C6_Y,C5_X,C5_Y)
top_area+=triangle_area(C1_X,C1_Y,C6_X,C6_Y,C7_X,C7_Y)
top_area+=triangle_area(C1_X,C1_Y,C8_X,C8_Y,C7_X,C7_Y)
fa=np.zeros(10)
fa[0]=side_area(C1_X,C1_Y,C2_X,C2_Y)
fa[1]=side_area(C3_X,C3_Y,C2_X,C2_Y)
fa[2]=side_area(C3_X,C3_Y,C4_X,C4_Y)
fa[3]=side_area(C5_X,C5_Y,C4_X,C4_Y)
fa[4]=side_area(C5_X,C5_Y,C6_X,C6_Y)
fa[5]=side_area(C7_X,C7_Y,C6_X,C6_Y)
fa[6]=side_area(C7_X,C7_Y,C8_X,C8_Y)
fa[7]=side_area(C1_X,C1_Y,C8_X,C8_Y)
fa[8]=fa[9]=top_area;
for iface in range(10):
    print "Area of face ",iface," is ",fa[iface]
fn=np.zeros((10,2))
fn[0]=side_normal(C1_X,C1_Y,C2_X,C2_Y)
fn[1]=-side_normal(C3_X,C3_Y,C2_X,C2_Y)
fn[2]=side_normal(C3_X,C3_Y,C4_X,C4_Y)
fn[3]=-side_normal(C5_X,C5_Y,C4_X,C4_Y)
fn[4]=side_normal(C5_X,C5_Y,C6_X,C6_Y)
fn[5]=-side_normal(C7_X,C7_Y,C6_X,C6_Y)
fn[6]=side_normal(C7_X,C7_Y,C8_X,C8_Y)
fn[7]=-side_normal(C1_X,C1_Y,C8_X,C8_Y)

phi_norm=np.zeros(8)
for iface in range(8):
    phi_norm[iface]=math.atan2(fn[iface][1],fn[iface][0])
    print "Normal to face ",iface," is ",fn[iface],",  phi_norm=",phi_norm[iface]

phi=math.atan2(C2_Y-C1_Y,C2_X-C1_X);print("x-sect at phi=",phi,": ",area(phi)," or ",area2(phi))
phi=math.atan2(C3_Y-C2_Y,C3_X-C2_X);print("x-sect at phi=",phi,": ",area(phi)," or ",area2(phi))
phi=math.atan2(C4_Y-C3_Y,C4_X-C3_X);print("x-sect at phi=",phi,": ",area(phi)," or ",area2(phi))
phi=math.atan2(C5_Y-C4_Y,C5_X-C4_X);print("x-sect at phi=",phi,": ",area(phi)," or ",area2(phi))
phi=math.atan2(C6_Y-C5_Y,C6_X-C5_X);print("x-sect at phi=",phi,": ",area(phi)," or ",area2(phi))
phi=math.atan2(C7_Y-C6_Y,C7_X-C6_X);print("x-sect at phi=",phi,": ",area(phi)," or ",area2(phi))
phi=math.atan2(C8_Y-C7_Y,C8_X-C7_X);print("x-sect at phi=",phi,": ",area(phi)," or ",area2(phi))
phi=math.atan2(C1_Y-C8_Y,C1_X-C8_X);print("x-sect at phi=",phi,": ",area(phi)," or ",area2(phi))
    
data=np.loadtxt("impactchain.dat")
#data=np.loadtxt("mini_impactchain.dat")
skydata=data[:,6:8]
print np.shape(data)
print data
print np.shape(skydata)
print skydata

pp = PdfPages('histograms.pdf')
print "histo err ~ ",50/math.sqrt(np.shape(data)[0]/50),"%"
n, bins, patches = plt.hist(skydata[:,1], 150, normed=1, facecolor='green', alpha=0.75)
plt.xlabel('phi')
plt.ylabel('Probability')
plt.grid(True)
pp.savefig()
plt.clf()

phineareq=[t[1] for t in skydata if abs(t[0])<0.1]
print "histo err ~ ",50/math.sqrt(np.shape(phineareq)[0]/50),"%"
n, bins, patches = plt.hist(phineareq, 50, normed=1, facecolor='green', alpha=0.75)
y=np.array([ area(b) for b in bins ])
mean_area=np.mean(y)
print "mean_area=",mean_area
scale=np.mean(n)/mean_area
y=y*scale
l = plt.plot(bins, y, 'r--', linewidth=1)
y=np.array([ area2(b) for b in bins ])
mean_area=np.mean(y)
print "mean_area=",mean_area
scale=np.mean(n)/mean_area
y=y*scale
l = plt.plot(bins, y, 'b--', linewidth=1)
plt.xlabel('phi-near-equator')
plt.ylabel('Probability')
plt.grid(True)
#plt.show()
pp.savefig()
plt.clf()

n, bins, patches = plt.hist(skydata[:,0], 50, normed=1, facecolor='green', alpha=0.75)
y=np.array([ area_3d(math.acos(b)) for b in bins ])
mean_area3d=np.mean(y)
scale=np.mean(n)/mean_area3d
y=y*scale
l = plt.plot(bins, y, 'r--', linewidth=1)
plt.xlabel('cos(theta)')
plt.ylabel('Probability')
plt.grid(True)
#plt.show()
pp.savefig()
plt.clf()

bins=np.arange(11)-0.5
n, bins, patches = plt.hist(data[:,8], bins, normed=1, facecolor='green', alpha=0.75)
l = plt.plot(np.arange(10)+0.0, fa/sum(fa), '--ro', linewidth=1)
plt.xlabel('face')
plt.ylabel('Probability')
plt.grid(True)
#plt.show()
pp.savefig()
plt.clf()

facecth=[t[6] for t in data if abs(t[8])==8]
nbins=50
n, bins, patches = plt.hist(facecth, nbins, normed=1, facecolor='green', alpha=0.75)
l = plt.plot(bins, -2*bins, '--r', linewidth=1)
plt.xlabel('cth on face 8')
plt.ylabel('Probability')
plt.grid(True)
#plt.show()
pp.savefig()
plt.clf()

facecth=[t[6] for t in data if abs(t[8])==9]
nbins=50
n, bins, patches = plt.hist(facecth, nbins, normed=1, facecolor='green', alpha=0.75)
l = plt.plot(bins, 2*bins, '--r', linewidth=1)
plt.xlabel('cth on face 9')
plt.ylabel('Probability')
plt.grid(True)
#plt.show()
pp.savefig()
plt.clf()

print "checking incidence angles:"

for t in data:
    iface=t[8]
    if(iface<8):
        if(iface<0 or math.cos(t[7]-phi_norm[iface])<0):
            print "err:"+str(t.tolist())

for iface in range(8):
    facephi=[t[7] for t in data if t[8]==iface]
    nbins=50
    n, bins, patches = plt.hist(facephi, nbins, normed=1, facecolor='green', alpha=0.75)
    l = plt.plot(bins, np.cos(bins-phi_norm[iface])/2, '--r', linewidth=1)
    plt.xlabel('phi on face '+str(iface))
    plt.ylabel('Probability')
    plt.grid(True)
    #plt.show()
    pp.savefig()
    plt.clf()

for iface in range(8):
    facecth=[t[6] for t in data if abs(t[8])==iface]
    nbins=50
    n, bins, patches = plt.hist(facecth, nbins, normed=1, facecolor='green', alpha=0.75)
    l = plt.plot(bins, np.sqrt(1-bins*bins)/(math.pi/2), '--r', linewidth=1)
    plt.xlabel('cth on face '+str(iface))
    plt.ylabel('Probability')
    plt.grid(True)
    #plt.show()
    pp.savefig()
    plt.clf()

pp.close()

