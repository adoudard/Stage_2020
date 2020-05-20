# -*- coding: utf-8 -*-

import os
import matplotlib.pyplot as plt
import numpy as np
from scipy import interpolate, optimize
import centers


##### Méthodes calculatoires #####

def centre_of_mass(list_x,list_y):
    
    assert len(list_x)==len(list_y), "Les deux listes de coordonnées doivent être de la même longueur."
    A=0
    Cx=0
    Cy=0
    
    for i in range (len(list_x)):
        B=(list_x[i-1] * list_y[i] - list_x[i] * list_y[i-1])
        
        A+= B/2
        
        Cx+= (1/6) * B * (list_x[i-1] + list_x[i])
        
        Cy+= (1/6) * B * (list_y[i-1] + list_y[i])
    
    return(Cx/A, Cy/A)
        




#Lecture des acquis C=270 et C=90 du CC4 :
    
x4_270,y4_270,r4_270,sx,sy = centers.Hough_center('E:\\Stage 2020\\Acquisitions 30-04', 'WL CC4-10-012.dcm', None, True, sig=0.5, pui=2, thr=0.7)
x4_90,y4_90,r4_90,sx,sy = centers.Hough_center('E:\\Stage 2020\\Acquisitions 30-04', 'WL CC4-10-013.dcm', None, True, sig=0.5, pui=2, thr=0.7)

ex = (x4_270-x4_90)*sx
ey = (y4_270-y4_90)*sy

cx = 0.5*(x4_270+x4_90)
cy = 0.5*(y4_270+y4_90)

# Ces valeurs nous serviront à générer des listes écart bille-champ virtuels avec le CC4, avec les suffixes v.

# Ouverture en lecture seule
analyse=open('E:\\Stage 2020\\Scripts\\30_04.txt',"r")
lines = analyse.readlines()[1:]  #Skip du header
Gs=[]  # List des Gantry Angles
Cs=[]  # List des Table Angles
Ts=[]  # List des Collimator Angles
Fs=[]  # List des Field Coordinates
Bs=[]  # List des Ball Coordinates

# Les informations sont séparées par des tabulations, codées par \t

for line in lines:
    
    Gs.append(eval(line.split("\t")[1]))
    Ts.append(eval(line.split("\t")[2]))
    Cs.append(eval(line.split("\t")[3]))
    Fs.append(eval(line.split("\t")[4]))
    Bs.append(eval(line.split("\t")[5]))
    
analyse.close()

# G0, T0, C variable
XFG0=np.array([Fs[i][0] for i in range (2)])
YFG0=np.array([Fs[i][1] for i in range (2)])
XBG0=np.array([Bs[i][0] for i in range (2)])
YBG0=np.array([Bs[i][1] for i in range (2)])

# G90, T0, C variable
XFG90=np.array([Fs[i+2][0] for i in range (2)])
YFG90=np.array([Fs[i+2][1] for i in range (2)])
XBG90=np.array([Bs[i+2][0] for i in range (2)])
YBG90=np.array([Bs[i+2][1] for i in range (2)])

# G180, T0, C variable
XFG180=np.array([Fs[i+4][0] for i in range (2)])
YFG180=np.array([Fs[i+4][1] for i in range (2)])
XBG180=np.array([Bs[i+4][0] for i in range (2)])
YBG180=np.array([Bs[i+4][1] for i in range (2)])

# G270, T0, C variable
XFG270=np.array([Fs[i+6][0] for i in range (2)])
YFG270=np.array([Fs[i+6][1] for i in range (2)])
XBG270=np.array([Bs[i+6][0] for i in range (2)])
YBG270=np.array([Bs[i+6][1] for i in range (2)])


#G0, C270, T variable
XFG0Ts=np.array([Fs[i+8][0] for i in range (4)])
YFG0Ts=np.array([Fs[i+8][1] for i in range (4)])
XBG0Ts=np.array([Bs[i+8][0] for i in range (4)])
YBG0Ts=np.array([Bs[i+8][1] for i in range (4)])

XFG0Ts=np.append(XFG0Ts,Fs[0][0])
YFG0Ts=np.append(YFG0Ts,Fs[0][1])
XBG0Ts=np.append(XBG0Ts,Bs[0][0])
YBG0Ts=np.append(YBG0Ts,Bs[0][1])

(x0, y0) = (np.mean(XFG0),np.mean(YFG0))
(x90, y90) = (np.mean(XFG90),np.mean(YFG90))
(x180, y180) = (np.mean(XFG180),np.mean(YFG180))
(x270, y270) = (np.mean(XFG270),np.mean(YFG270))
(x0t, y0t) = (np.mean(XFG0Ts),np.mean(YFG0Ts))

dx270 = x4_270 - XFG0[0]
dy270 = y4_270 - YFG0[0]
dx90 = x4_90 - XFG0[1]
dy90 = y4_90 - YFG0[1]

XFG0v = XFG0 + np.array([dx270, dx90])
YFG0v = YFG0 + np.array([dy270, dy90])

XFG90v = XFG90 + np.array([dx90, dx270])
YFG90v = YFG90 + np.array([dy90, dy270])

XFG180v = XFG180 + np.array([dx270, dx90])
YFG180v = YFG180 + np.array([dy270, dy90])

XFG270v = XFG270 + np.array([dx90, dx270])
YFG270v = YFG270 + np.array([dy90, dy270])

XFG0Tsv = XFG0Ts + np.array([dx270, dx270, dx270, dx270, dx270])
YFG0Tsv = YFG0Ts + np.array([dy270, dy270, dy270, dy270, dy270])

(x0v, y0v) = (np.mean(XFG0v),np.mean(YFG0v))
(x90v, y90v) = (np.mean(XFG90v),np.mean(YFG90v))
(x180v, y180v) = (np.mean(XFG180v),np.mean(YFG180v))
(x270v, y270v) = (np.mean(XFG270v),np.mean(YFG270v))
(x0tv, y0tv) = (np.mean(XFG0Tsv),np.mean(YFG0Tsv))


s=(2/3)*0.336 #scale


# Système de coordonnées : 
    # z = axe de la table
    # y = axe antéro postérieur
    # x = axe pour que x,y,z soit othogonal et orienté positivement.
    

AX0=s*(XBG0-XFG0)  # [270,90]
AY0=s*(YBG0-YFG0)

AX0v=s*(XBG0-XFG0v)  # [270,90]
AY0v=s*(YBG0-YFG0v)

AX90=s*(XBG90-XFG90)  # [90,270]
AY90=s*(YBG90-YFG90)

AX90v=s*(XBG90-XFG90v)  # [90,270]
AY90v=s*(YBG90-YFG90v)

AX180=s*(XBG180-XFG180)  # [270,90]
AY180=s*(YBG180-YFG180)

AX180v=s*(XBG180-XFG180v)  # [270,90]
AY180v=s*(YBG180-YFG180v)

AX270=s*(XBG270-XFG270)  # [90,270]
AY270=s*(YBG270-YFG270)

AX270v=s*(XBG270-XFG270v)  # [90,270]
AY270v=s*(YBG270-YFG270v)

AX0T=s*(XBG0Ts-XFG0Ts)  # [270]
AY0T=s*(YBG0Ts-YFG0Ts)

AX0Tv=s*(XBG0Ts-XFG0Tsv)  # [270]
AY0Tv=s*(YBG0Ts-YFG0Tsv)

#Valeurs mesurées par RIT

AX0r=np.array([-0.46,-0.07])
AY0r=np.array([0.27,0.18])

AX90r=np.array([0.08,-0.28])
AY90r=np.array([0.01,0.03])

AX180r=np.array([0.02,0.36])
AY180r=np.array([-0.2,-0.23])

AX270r=np.array([0.22,-0.14])
AY270r=np.array([-0.04,-0.01])

AX0Tr=np.array([-0.05,-0.13,-0.5,-0.43,-0.46])
AY0Tr=np.array([0.02,0.17,0.48,0.64,0.27])

AY180T = AY0T - AY0[0] + AY180[0]
AY180Tr = AY0Tr - AY0r[0] + AY180r[0]

AY180Tp = AY0T - AY0[0] + AY180[1]
AY180Tpr = AY0Tr - AY0r[0] + AY180r[1]

AX180T = AX0T - AX0[0] - AX180[0]
AX180Tr = AX0Tr - AX0r[0] - AX180r[0]

AX180Tp = AX0T - AX0[0] - AX180[1]
AX180Tpr = AX0Tr - AX0r[0] - AX180r[1]



angles=[(0,0),(+90,0),(+180,0),(+270,0),(0,-90),(0,-45),(0,45),(0,90)]



def Low(AX0,AY0,AX90,AY90,AX180,AY180,AX270,AY270,AX0T,AY0T):
    
    A=np.array([[1,0,0],[0,1,0],[1,0,0],[0,0,-1],[1,0,0],[0,-1,0],[1,0,0],[0,0,1],[0,1,0],[-1,0,0],[0.5*np.sqrt(2),0.5*np.sqrt(2),0],[-0.5*np.sqrt(2),0.5*np.sqrt(2),0],[0.5*np.sqrt(2),-0.5*np.sqrt(2),0],[0.5*np.sqrt(2),0.5*np.sqrt(2),0],[0,-1,0],[1,0,0]])
    
    B=np.dot(np.linalg.inv(np.dot(A.transpose(),A)),A.transpose())

    ksi=np.array([AY0[0],AX0[0],AY90[1],AX90[1],AY180[0],AX180[0],AY270[1],AX270[1],AY0T[0],AX0T[0],AY0T[1],AX0T[1],AY0T[2],AX0T[2],AY0T[3],AX0T[3]])
    
    delta=np.dot(B,ksi.transpose())
    return delta, A, ksi, B

delta, A, ksi, B=Low(AX0,AY0,AX90,AY90,AX180,AY180,AX270,AY270,AX0T,AY0T)
# deltav=Low(AX0v,AY0v,AX90v,AY90v,AX180v,AY180v,AX270v,AY270v,AX0Tv,AY0Tv)
deltar, Ar, ksir, Br=Low(AX0r,AY0r,AX90r,AY90r,AX180r,AY180r,AX270r,AY270r,AX0Tr,AY0Tr)

print("Optimisation de Low : La bille est excentrée de dX = ", round(delta[0],3)," mm, dY = ", round(delta[1],3), " mm et de dZ = ", round(delta[2],3)," mm dans le système standard, avec les seules mesures en CC270, soit un dR = ",round(np.sqrt(delta[0]**2+delta[1]**2+0**2),3), "mm.")    
      
plt.figure()

plt.subplot(234)
plt.text(0,0,"Avant correction de la position de la bille")
plt.xlim(-1,5)
plt.ylim(-0.5,0.5)
plt.xticks([])
plt.yticks([])

plt.subplot(232)
shift0=(np.mean(s*XBG0),np.mean(s*YBG0))
ray=s*np.sqrt(np.max((x0-XFG0)**2+(y0-YFG0)**2))
circle_f=plt.Circle((s*x0-shift0[0],s*y0-shift0[1]),radius=ray,color='red',fill=False)
plt.scatter(s*XFG0-shift0[0],s*YFG0-shift0[1],color='red',label='Champ')
plt.scatter(s*XBG0-shift0[0],s*YBG0-shift0[1],color='orange',label='Bille')
plt.scatter(s*x0-shift0[0],s*y0-shift0[1],color='brown',label='Centroïde Champ')
plt.xlim(-1,1)
plt.ylim(-1,1)
fig2 = plt.gcf()
ax2 = fig2.gca()
ax2.add_artist(circle_f)
plt.text(s*(x0-1)-shift0[0],s*(y0-1)-shift0[1], "Excentrage : "+str(round(s*ray,3))+ "mm")
plt.text(s*np.mean(XBG0)-shift0[0],s*np.mean(YBG0)-shift0[1],"Incertitude : "+str(round(s*np.sqrt(np.std(XBG0)**2+np.std(YBG0)**2),3))+ "mm")
plt.title('G=0, C variable')
plt.legend(loc='upper left')
plt.xticks(np.linspace(-1,1,9))
plt.grid(linestyle='--')
plt.xlabel('Y,-v, mm')
plt.ylabel('X, u, mm')
plt.gca().set_aspect('equal', adjustable='box')

plt.subplot(233)
shift90=(np.mean(s*XBG90),np.mean(s*YBG90))
ray=s*np.sqrt(np.max((x90-XFG90)**2+(y90-YFG90)**2))
circle_f=plt.Circle((s*x90-shift90[0],s*y90-shift90[1]),radius=ray,color='red',fill=False)
plt.scatter(s*XFG90-shift90[0],s*YFG90-shift90[1],color='red')
plt.scatter(s*XBG90-shift90[0],s*YBG90-shift90[1],color='orange')
plt.scatter(s*x90-shift90[0],s*y90-shift90[1],color='brown')
plt.xlim(-1,1)
plt.ylim(-1,1)
fig2 = plt.gcf()
ax2 = fig2.gca()
ax2.add_artist(circle_f)
plt.text(s*(x90-1)-shift90[0],s*(y90-1)-shift90[1], "Excentrage : "+str(round(s*ray,3))+ "mm")
plt.text(s*np.mean(XBG90)-shift90[0],s*np.mean(YBG90)-shift90[1],"Incertitude : "+str(round(s*np.sqrt(np.std(XBG90)**2+np.std(YBG90)**2),3))+ "mm")
plt.title('G=90, C variable')
plt.xlabel('-Z, -v, mm')
plt.ylabel('X, u, mm')
plt.xticks(np.linspace(-1,1,9))
plt.grid(linestyle='--')
plt.gca().set_aspect('equal', adjustable='box')

plt.subplot(235)
shift180=(np.mean(s*XBG180),np.mean(s*YBG180))
ray=s*np.sqrt(np.max((x180-XFG180)**2+(y180-YFG180)**2))
circle_f=plt.Circle((s*x180-shift180[0],s*y180-shift180[1]),radius=ray,color='red',fill=False)
plt.scatter(s*XFG180-shift180[0],s*YFG180-shift180[1],color='red')
plt.scatter(s*XBG180-shift180[0],s*YBG180-shift180[1],color='orange')
plt.scatter(s*x180-shift180[0],s*y180-shift180[1],color='brown')
plt.xlim(-1,1)
plt.ylim(-1,1)
fig2 = plt.gcf()
ax2 = fig2.gca()
ax2.add_artist(circle_f)
plt.text(s*(x180-1)-shift180[0],s*(y180-1)-shift180[1], "Excentrage : "+str(round(s*ray,3))+ "mm")
plt.text(s*np.mean(XBG180)-shift180[0],s*np.mean(YBG180)-shift180[1],"Incertitude : "+str(round(0.224*np.sqrt(np.std(XBG180)**2+np.std(YBG180)**2),3))+ "mm")
plt.title('G=180, C variable')
plt.xlabel('-Y, -v, mm')
plt.ylabel('X, u, mm')
plt.xticks(np.linspace(-1,1,9))
plt.grid(linestyle='--')
plt.gca().set_aspect('equal', adjustable='box')

plt.subplot(231)
shift=(np.mean(s*XFG0Ts),np.mean(s*YFG0Ts))
plt.scatter(s*XFG0Ts-shift[0],s*YFG0Ts-shift[1],color='red')
plt.scatter(s*XBG0Ts-shift[0],s*YBG0Ts-shift[1],color='orange')
plt.scatter(s*x0t-shift[0],s*y0t-shift[1],color='brown')
plt.xlim(-1,1)
plt.ylim(-1,1)
plt.title('G=0, T variable')
plt.xlabel('mm')
plt.ylabel('mm')
plt.xticks(np.linspace(-1,1,9))
plt.grid(linestyle='--')
plt.gca().set_aspect('equal', adjustable='box')

plt.subplot(236)
shift270=(np.mean(s*XBG270),np.mean(s*YBG270))
ray=s*np.sqrt(np.max((x270-XFG270)**2+(y270-YFG270)**2))
circle_f=plt.Circle((s*x270-shift270[0],s*y270-shift270[1]),radius=ray,color='red',fill=False)
plt.scatter(s*XFG270-shift270[0],s*YFG270-shift270[1],color='red')
plt.scatter(s*XBG270-shift270[0],s*YBG270-shift270[1],color='orange')
plt.scatter(s*x270-shift270[0],s*y270-shift270[1],color='brown')
plt.xlim(-1,1)
plt.ylim(-1,1)
fig3 = plt.gcf()
ax3 = fig3.gca()
ax3.add_artist(circle_f)
plt.text(s*(x270-1)-shift270[0],s*(y270-1)-shift270[1], "Excentrage : "+str(round(s*ray,3))+ "mm")
plt.text(s*np.mean(XBG270)-shift270[0],s*np.mean(YBG270)-shift270[1],"Incertitude : "+str(round(0.224*np.sqrt(np.std(XBG270)**2+np.std(YBG270)**2),3))+ "mm")
plt.title('G=270, C variable')
plt.xlabel('Z, -v, mm')
plt.ylabel('X, u, mm')
plt.xticks(np.linspace(-1,1,9))
plt.grid(linestyle='--')
plt.gca().set_aspect('equal', adjustable='box')


# "Avant correction de la position de la bille\nWL virtuels avec CC4"

plt.figure()

plt.subplot(234)
plt.text(0,0,"Avant correction de la position de la bille\nWL virtuels avec CC4")
plt.xlim(-1,5)
plt.ylim(-0.5,0.5)
plt.xticks([])
plt.yticks([])

plt.subplot(232)
shift0=(np.mean(s*XBG0),np.mean(s*YBG0))
ray=s*np.sqrt(np.max((x0v-XFG0v)**2+(y0v-YFG0v)**2))
circle_f=plt.Circle((s*x0v-shift0[0],s*y0v-shift0[1]),radius=ray,color='red',fill=False)
plt.scatter(s*XFG0v-shift0[0],s*YFG0v-shift0[1],color='red',label='Champ')
plt.scatter(s*XBG0-shift0[0],s*YBG0-shift0[1],color='orange',label='Bille')
plt.scatter(s*x0v-shift0[0],s*y0v-shift0[1],color='brown',label='Centroïde Champ')
plt.xlim(-1,1)
plt.ylim(-1,1)
fig2 = plt.gcf()
ax2 = fig2.gca()
ax2.add_artist(circle_f)
plt.text(s*(x0v-1)-shift0[0],s*(y0v-1)-shift0[1], "Excentrage : "+str(round(s*ray,3))+ "mm")
plt.text(s*np.mean(XBG0)-shift0[0],s*np.mean(YBG0)-shift0[1],"Incertitude : "+str(round(s*np.sqrt(np.std(XBG0)**2+np.std(YBG0)**2),3))+ "mm")
plt.title('G=0, C variable')
plt.legend(loc='upper left')
plt.xticks(np.linspace(-1,1,9))
plt.grid(linestyle='--')
plt.xlabel('Y,-v, mm')
plt.ylabel('X, u, mm')
plt.gca().set_aspect('equal', adjustable='box')

plt.subplot(233)
shift90=(np.mean(s*XBG90),np.mean(s*YBG90))
ray=s*np.sqrt(np.max((x90v-XFG90v)**2+(y90v-YFG90v)**2))
circle_f=plt.Circle((s*x90v-shift90[0],s*y90v-shift90[1]),radius=ray,color='red',fill=False)
plt.scatter(s*XFG90v-shift90[0],s*YFG90v-shift90[1],color='red')
plt.scatter(s*XBG90-shift90[0],s*YBG90-shift90[1],color='orange')
plt.scatter(s*x90v-shift90[0],s*y90v-shift90[1],color='brown')
plt.xlim(-1,1)
plt.ylim(-1,1)
fig2 = plt.gcf()
ax2 = fig2.gca()
ax2.add_artist(circle_f)
plt.text(s*(x90v-1)-shift90[0],s*(y90v-1)-shift90[1], "Excentrage : "+str(round(s*ray,3))+ "mm")
plt.text(s*np.mean(XBG90)-shift90[0],s*np.mean(YBG90)-shift90[1],"Incertitude : "+str(round(s*np.sqrt(np.std(XBG90)**2+np.std(YBG90)**2),3))+ "mm")
plt.title('G=90, C variable')
plt.xlabel('-Z, -v, mm')
plt.ylabel('X, u, mm')
plt.xticks(np.linspace(-1,1,9))
plt.grid(linestyle='--')
plt.gca().set_aspect('equal', adjustable='box')

plt.subplot(235)
shift180=(np.mean(s*XBG180),np.mean(s*YBG180))
ray=s*np.sqrt(np.max((x180v-XFG180v)**2+(y180v-YFG180v)**2))
circle_f=plt.Circle((s*x180v-shift180[0],s*y180v-shift180[1]),radius=ray,color='red',fill=False)
plt.scatter(s*XFG180v-shift180[0],s*YFG180v-shift180[1],color='red')
plt.scatter(s*XBG180-shift180[0],s*YBG180-shift180[1],color='orange')
plt.scatter(s*x180v-shift180[0],s*y180v-shift180[1],color='brown')
plt.xlim(-1,1)
plt.ylim(-1,1)
fig2 = plt.gcf()
ax2 = fig2.gca()
ax2.add_artist(circle_f)
plt.text(s*(x180v-1)-shift180[0],s*(y180v-1)-shift180[1], "Excentrage : "+str(round(s*ray,3))+ "mm")
plt.text(s*np.mean(XBG180)-shift180[0],s*np.mean(YBG180)-shift180[1],"Incertitude : "+str(round(0.224*np.sqrt(np.std(XBG180)**2+np.std(YBG180)**2),3))+ "mm")
plt.title('G=180, C variable')
plt.xlabel('Y, -v, mm')
plt.ylabel('X, u, mm')
plt.xticks(np.linspace(-1,1,9))
plt.grid(linestyle='--')
plt.gca().set_aspect('equal', adjustable='box')

plt.subplot(231)
shift=(np.mean(s*XFG0Tsv),np.mean(s*YFG0Tsv))
plt.scatter(s*XFG0Tsv-shift[0],s*YFG0Tsv-shift[1],color='red')
plt.scatter(s*XBG0Ts-shift[0],s*YBG0Ts-shift[1],color='orange')
plt.scatter(s*x0tv-shift[0],s*y0tv-shift[1],color='brown')
plt.xlim(-1,1)
plt.ylim(-1,1)
plt.title('G=0, T variable')
plt.xlabel('mm')
plt.ylabel('mm')
plt.xticks(np.linspace(-1,1,9))
plt.grid(linestyle='--')
plt.gca().set_aspect('equal', adjustable='box')

plt.subplot(236)
shift270=(np.mean(s*XBG270),np.mean(s*YBG270))
ray=s*np.sqrt(np.max((x270v-XFG270v)**2+(y270v-YFG270v)**2))
circle_f=plt.Circle((s*x270v-shift270[0],s*y270v-shift270[1]),radius=ray,color='red',fill=False)
plt.scatter(s*XFG270v-shift270[0],s*YFG270v-shift270[1],color='red')
plt.scatter(s*XBG270-shift270[0],s*YBG270-shift270[1],color='orange')
plt.scatter(s*x270v-shift270[0],s*y270v-shift270[1],color='brown')
plt.xlim(-1,1)
plt.ylim(-1,1)
fig3 = plt.gcf()
ax3 = fig3.gca()
ax3.add_artist(circle_f)
plt.text(s*(x270v-1)-shift270[0],s*(y270v-1)-shift270[1], "Excentrage : "+str(round(s*ray,3))+ "mm")
plt.text(s*np.mean(XBG270)-shift270[0],s*np.mean(YBG270)-shift270[1],"Incertitude : "+str(round(0.224*np.sqrt(np.std(XBG270)**2+np.std(YBG270)**2),3))+ "mm")
plt.title('G=270, C variable')
plt.xlabel('Z, -v, mm')
plt.ylabel('X, u, mm')
plt.xticks(np.linspace(-1,1,9))
plt.grid(linestyle='--')
plt.gca().set_aspect('equal', adjustable='box')

##Methode d'optimisation de Hancock##

theta=np.array([np.radians(0),np.radians(90),np.radians(180),np.radians(270)])

def calcul_PM(AY0, AY90, AY180, AY270, AX0, AX90, AX180, AX270):

    #1.1 Isocentre bras avec C=90,T=0
    
    ug=np.array([AY0[1],AY90[0],AY180[1],AY270[0]])
    vg=np.array([AX0[1],AX90[0],AX180[1],AX270[0]])
    
    xg=ug
    yg=vg*np.cos(theta)
    zg=-1*vg*np.sin(theta)
    
    PG=0.5*np.array([np.max(xg)+np.min(xg), np.max(yg)+np.min(yg), np.max(zg)+np.min(zg)])
    
    #1.2 Isocentre bras avec C=270,T=0
    
    uc=np.array([AY0[0],AY90[1],AY180[0],AY270[1]])
    vc=np.array([AX0[0],AX90[1],AX180[0],AX270[1]])
    
    xc=uc
    yc=vc*np.cos(theta)
    zc=-1*vc*np.sin(theta)

    PC=0.5*np.array([np.max(xc)+np.min(xc), np.max(yc)+np.min(yc), np.max(zc)+np.min(zc)])
    
    #1.3 Isocentre moyen
    
    PM=0.5*(PG+PC)
    
    return PM,PC,PG

PM,PC,PG = calcul_PM(AY0, AY90, AY180, AY270, AX0, AX90, AX180, AX270)

PMv,PCv,PGv = calcul_PM(AY0v, AY90v, AY180v, AY270v, AX0v, AX90v, AX180v, AX270v)

PMr,PCr,PGr = calcul_PM(AY0r, AY90r, AY180r, AY270r, AX0r, AX90r, AX180r, AX270r)

# On va soustraire PM (ses projections) à la bille

def calculs_GT_AB(AX0,AY0,AX90,AY90,AX180,AY180,AX270,AY270,PM):

    GT90=np.array([AY0[1]-PM[0], AY90[0]-PM[0], AY180[1]-PM[0], AY270[0]-PM[0], AY0[1]-PM[0]])
    GT270=np.array([AY0[0]-PM[0], AY90[1]-PM[0], AY180[0]-PM[0], AY270[1]-PM[0], AY0[0]-PM[0]])
    AB90=np.array([AX0[1]-PM[1], AX90[0]+PM[2], AX180[1]+PM[1], AX270[0]-PM[2], AX0[1]-PM[1]])   #  *-1 car le poly RIT n'est pas cohérent
    AB270=np.array([AX0[0]-PM[1], AX90[1]+PM[2], AX180[0]+PM[1], AX270[1]-PM[2], AX0[0]-PM[1]])  #  avec lui-même ... Même si ça ne change fondamentalement rien
    
    return (GT90,GT270,AB90,AB270)

(GT90,GT270,AB90,AB270) = calculs_GT_AB(AX0,AY0,AX90,AY90,AX180,AY180,AX270,AY270,PM)

(GT90v,GT270v,AB90v,AB270v) = calculs_GT_AB(AX0v,AY0v,AX90v,AY90v,AX180v,AY180v,AX270v,AY270v,PMv)

(GT90r,GT270r,AB90r,AB270r) = calculs_GT_AB(AX0r,AY0r,AX90r,AY90r,AX180r,AY180r,AX270r,AY270r,PMr)

plt.figure()
x=np.append(np.degrees(theta),360.)
plt.scatter(x, AB90, color='blue')
plt.scatter(x, AB270,color='blue')
plt.scatter(x, GT90, color='lime')
plt.scatter(x, GT270, color='lime')
plt.scatter(x, np.sqrt(GT90**2+AB90**2), color='red')
plt.scatter(x, np.sqrt(GT270**2+AB270**2), color='red')
tAB90 = interpolate.splrep(x, AB90, s=0)
tAB270 = interpolate.splrep(x, AB270, s=0)
tGT90 = interpolate.splrep(x, GT90, s=0)
tGT270 = interpolate.splrep(x, GT270, s=0)
xnew = np.arange(0, 360.5, 0.5)
yAB90 = interpolate.splev(xnew, tAB90, der=0)
yAB270 = interpolate.splev(xnew, tAB270, der=0)
yGT90 = interpolate.splev(xnew, tGT90, der=0)
yGT270 = interpolate.splev(xnew, tGT270, der=0)
plt.plot(xnew,yAB270,color='blue', linestyle='-',label='AB270')
plt.plot(xnew,yAB90,color='blue', linestyle=':',label='AB90')
plt.plot(xnew,yGT270,color='lime', linestyle='-',label='GT270')
plt.plot(xnew,yGT90,color='lime', linestyle=':',label='GT90')
plt.plot(xnew,np.sqrt(yGT270**2+yAB270**2),color='red', linestyle='-',label='Tot270')
plt.plot(xnew,np.sqrt(yGT90**2+yAB90**2),color='red', linestyle=':',label='Tot90')
plt.plot(xnew,np.ones(len(xnew)),color='black')
plt.plot(xnew,-1*np.ones(len(xnew)),color='black')
plt.ylim(-0.3,0.4)
plt.xlim(0,360)
plt.legend(loc='upper left')
plt.ylabel('mm')
plt.xlabel('Angle Bras')
plt.title("Déplacements latéraux(AB) et avant-arrière(GT) \n en fonction de l'angle bras")
plt.grid(linestyle='--')

plt.figure()

plt.subplot(234)
plt.text(0,0,"Après correction de la position de la bille\nsans prise en compte de la table")
plt.xlim(-1,5)
plt.ylim(-0.5,0.5)
plt.xticks([])
plt.yticks([])

plt.subplot(232)
shift0=(np.mean(s*XBG0-PM[1]),np.mean(s*YBG0-PM[0]))
ray=s*np.sqrt(np.max((x0-XFG0)**2+(y0-YFG0)**2))
circle_f=plt.Circle((s*x0-shift0[0],s*y0-shift0[1]),radius=ray,color='red',fill=False)
plt.scatter(s*XFG0-shift0[0],s*YFG0-shift0[1],color='red',label='Champ')
plt.scatter(s*XBG0-PM[1]-shift0[0],s*YBG0-PM[0]-shift0[1],color='orange',label='Bille')
plt.scatter(s*x0-shift0[0],s*y0-shift0[1],color='brown',label='Centroïde Champ')
plt.xlim(-1,1)
plt.ylim(-1,1)
fig2 = plt.gcf()
ax2 = fig2.gca()
ax2.add_artist(circle_f)
plt.text(s*(x0-1)-shift0[0],s*(y0-1)-shift0[1], "Excentrage : "+str(round(s*ray,3))+ "mm")
plt.text(s*np.mean(XBG0)-shift0[0],s*np.mean(YBG0)-shift0[1],"Incertitude : "+str(round(s*np.sqrt(np.std(XBG0)**2+np.std(YBG0)**2),3))+ "mm")
plt.title('G=0, C variable')
plt.legend(loc='upper left')
plt.xticks(np.linspace(-1,1,9))
plt.grid(linestyle='--')
plt.xlabel('Y,-v, mm')
plt.ylabel('X, u, mm')
plt.gca().set_aspect('equal', adjustable='box')

plt.subplot(233)
shift90=(np.mean(s*XBG90+PM[2]),np.mean(s*YBG90-PM[0]))
ray=s*np.sqrt(np.max((x90-XFG90)**2+(y90-YFG90)**2))
circle_f=plt.Circle((s*x90-shift90[0],s*y90-shift90[1]),radius=ray,color='red',fill=False)
plt.scatter(s*XFG90-shift90[0],s*YFG90-shift90[1],color='red')
plt.scatter(s*XBG90+PM[2]-shift90[0],s*YBG90-PM[0]-shift90[1],color='orange')
plt.scatter(s*x90-shift90[0],s*y90-shift90[1],color='brown')
plt.xlim(-1,1)
plt.ylim(-1,1)
fig2 = plt.gcf()
ax2 = fig2.gca()
ax2.add_artist(circle_f)
plt.text(s*(x90-1)-shift90[0],s*(y90-1)-shift90[1], "Excentrage : "+str(round(s*ray,3))+ "mm")
plt.text(s*np.mean(XBG90)-shift90[0],s*np.mean(YBG90)-shift90[1],"Incertitude : "+str(round(s*np.sqrt(np.std(XBG90)**2+np.std(YBG90)**2),3))+ "mm")
plt.title('G=90, C variable')
plt.xlabel('-Z, -v, mm')
plt.ylabel('X, u, mm')
plt.xticks(np.linspace(-1,1,9))
plt.grid(linestyle='--')
plt.gca().set_aspect('equal', adjustable='box')

plt.subplot(235)
shift180=(np.mean(s*XBG180+PM[1]),np.mean(s*YBG180-PM[0]))
ray=s*np.sqrt(np.max((x180-XFG180)**2+(y180-YFG180)**2))
circle_f=plt.Circle((s*x180-shift180[0],s*y180-shift180[1]),radius=ray,color='red',fill=False)
plt.scatter(s*XFG180-shift180[0],s*YFG180-shift180[1],color='red')
plt.scatter(s*XBG180+PM[1]-shift180[0],s*YBG180-PM[0]-shift180[1],color='orange')
plt.scatter(s*x180-shift180[0],s*y180-shift180[1],color='brown')
plt.xlim(-1,1)
plt.ylim(-1,1)
fig2 = plt.gcf()
ax2 = fig2.gca()
ax2.add_artist(circle_f)
plt.text(s*(x180-1)-shift180[0],s*(y180-1)-shift180[1], "Excentrage : "+str(round(s*ray,3))+ "mm")
plt.text(s*np.mean(XBG180)-shift180[0],s*np.mean(YBG180)-shift180[1],"Incertitude : "+str(round(0.224*np.sqrt(np.std(XBG180)**2+np.std(YBG180)**2),3))+ "mm")
plt.title('G=180, C variable')
plt.xlabel('-Y, -v, mm')
plt.ylabel('X, u, mm')
plt.xticks(np.linspace(-1,1,9))
plt.grid(linestyle='--')
plt.gca().set_aspect('equal', adjustable='box')

plt.subplot(231)
shift=(np.mean(s*XFG0Ts-PM[1]),np.mean(s*YFG0Ts+PM[0]))
plt.scatter(s*XFG0Ts-PM[1]-shift[0],s*YFG0Ts+PM[0]-shift[1],color='red')
plt.scatter(s*XBG0Ts-shift[0],s*YBG0Ts-shift[1],color='orange')
plt.scatter(s*x0t-PM[1]-shift[0],s*y0t+PM[0]-shift[1],color='brown')
plt.xlim(-1,1)
plt.ylim(-1,1)
plt.title('G=0, T variable')
plt.xlabel('mm')
plt.ylabel('mm')
plt.xticks(np.linspace(-1,1,9))
plt.grid(linestyle='--')
plt.gca().set_aspect('equal', adjustable='box')

plt.subplot(236)
shift270=(np.mean(s*XBG270-PM[2]),np.mean(s*YBG270-PM[0]))
ray=s*np.sqrt(np.max((x270-XFG270)**2+(y270-YFG270)**2))
circle_f=plt.Circle((s*x270-shift270[0],s*y270-shift270[1]),radius=ray,color='red',fill=False)
plt.scatter(s*XFG270-shift270[0],s*YFG270-shift270[1],color='red')
plt.scatter(s*XBG270-PM[2]-shift270[0],s*YBG270-PM[0]-shift270[1],color='orange')
plt.scatter(s*x270-shift270[0],s*y270-shift270[1],color='brown')
plt.xlim(-1,1)
plt.ylim(-1,1)
fig3 = plt.gcf()
ax3 = fig3.gca()
ax3.add_artist(circle_f)
plt.text(s*(x270-1)-shift270[0],s*(y270-1)-shift270[1], "Excentrage : "+str(round(s*ray,3))+ "mm")
plt.text(s*np.mean(XBG270)-shift270[0],s*np.mean(YBG270)-shift270[1],"Incertitude : "+str(round(0.224*np.sqrt(np.std(XBG270)**2+np.std(YBG270)**2),3))+ "mm")
plt.title('G=270, C variable')
plt.xlabel('Z, -v, mm')
plt.ylabel('X, u, mm')
plt.xticks(np.linspace(-1,1,9))
plt.grid(linestyle='--')
plt.gca().set_aspect('equal', adjustable='box')

print("Optimisation de Hancock sans table : La bille est excentrée de dX = ", round(PM[0],3)," mm, dY = ", round(PM[1],3), " mm et de dZ = ", round(PM[2],3)," mm dans le système de coorodonnées de base.")

plt.figure()
x=np.append(np.degrees(theta),360.)
plt.scatter(x, AB90v, color='blue')
plt.scatter(x, AB270v,color='blue')
plt.scatter(x, GT90v, color='lime')
plt.scatter(x, GT270v, color='lime')
plt.scatter(x, np.sqrt(GT90v**2+AB90v**2), color='red')
plt.scatter(x, np.sqrt(GT270v**2+AB270v**2), color='red')
tAB90v = interpolate.splrep(x, AB90v, s=0)
tAB270v = interpolate.splrep(x, AB270v, s=0)
tGT90v = interpolate.splrep(x, GT90v, s=0)
tGT270v = interpolate.splrep(x, GT270v, s=0)
xnew = np.arange(0, 360, 0.5)
yAB90v = interpolate.splev(xnew, tAB90v, der=0)
yAB270v = interpolate.splev(xnew, tAB270v, der=0)
yGT90v = interpolate.splev(xnew, tGT90v, der=0)
yGT270v = interpolate.splev(xnew, tGT270v, der=0)
plt.plot(xnew,yAB270v,color='blue', linestyle='-',label='AB270 virtuel')
plt.plot(xnew,yAB90v,color='blue', linestyle=':',label='AB90 virtuel')
plt.plot(xnew,yGT270v,color='lime', linestyle='-',label='GT270 virtuel')
plt.plot(xnew,yGT90v,color='lime', linestyle=':',label='GT90 virtuel')
plt.plot(xnew,np.sqrt(yGT270v**2+yAB270v**2),color='red', linestyle='-',label='Tot270 virtuel')
plt.plot(xnew,np.sqrt(yGT90v**2+yAB90v**2),color='red', linestyle=':',label='Tot90 virtuel')
plt.plot(xnew,np.ones(len(xnew)),color='black')
plt.plot(xnew,-1*np.ones(len(xnew)),color='black')
plt.ylim(-0.3,0.4)
plt.xlim(0,360)
plt.legend(loc='upper left')
plt.ylabel('mm')
plt.xlabel('Angle Bras')
plt.title("Déplacements latéraux(AB) et avant-arrière(GT) \n en fonction de l'angle bras\n avec mesures virtuelles CC4")
plt.grid(linestyle='--')

plt.figure()

plt.subplot(234)
plt.text(0,0,"Après correction de la position de la bille\nsans prise en compte de la table\nen virtuel CC4")
plt.xlim(-1,5)
plt.ylim(-0.5,0.5)
plt.xticks([])
plt.yticks([])

plt.subplot(232)
shift0=(np.mean(s*XBG0-PMv[1]),np.mean(s*YBG0-PMv[0]))
ray=s*np.sqrt(np.max((x0v-XFG0v)**2+(y0v-YFG0v)**2))
circle_f=plt.Circle((s*x0v-shift0[0],s*y0v-shift0[1]),radius=ray,color='red',fill=False)
plt.scatter(s*XFG0v-shift0[0],s*YFG0v-shift0[1],color='red',label='Champ')
plt.scatter(s*XBG0-PMv[1]-shift0[0],s*YBG0-PMv[0]-shift0[1],color='orange',label='Bille')
plt.scatter(s*x0v-shift0[0],s*y0v-shift0[1],color='brown',label='Centroïde Champ')
plt.xlim(-1,1)
plt.ylim(-1,1)
fig2 = plt.gcf()
ax2 = fig2.gca()
ax2.add_artist(circle_f)
plt.text(s*(x0v-1)-shift0[0],s*(y0v-1)-shift0[1], "Excentrage : "+str(round(s*ray,3))+ "mm")
plt.text(s*np.mean(XBG0)-shift0[0],s*np.mean(YBG0)-shift0[1],"Incertitude : "+str(round(s*np.sqrt(np.std(XBG0)**2+np.std(YBG0)**2),3))+ "mm")
plt.title('G=0, C variable')
plt.legend(loc='upper left')
plt.xticks(np.linspace(-1,1,9))
plt.grid(linestyle='--')
plt.xlabel('Y,-v, mm')
plt.ylabel('X, u, mm')
plt.gca().set_aspect('equal', adjustable='box')

plt.subplot(233)
shift90=(np.mean(s*XBG90+PMv[2]),np.mean(s*YBG90-PMv[0]))
ray=s*np.sqrt(np.max((x90v-XFG90v)**2+(y90v-YFG90v)**2))
circle_f=plt.Circle((s*x90v-shift90[0],s*y90v-shift90[1]),radius=ray,color='red',fill=False)
plt.scatter(s*XFG90v-shift90[0],s*YFG90v-shift90[1],color='red')
plt.scatter(s*XBG90+PMv[2]-shift90[0],s*YBG90-PMv[0]-shift90[1],color='orange')
plt.scatter(s*x90v-shift90[0],s*y90v-shift90[1],color='brown')
plt.xlim(-1,1)
plt.ylim(-1,1)
fig2 = plt.gcf()
ax2 = fig2.gca()
ax2.add_artist(circle_f)
plt.text(s*(x90v-1)-shift90[0],s*(y90v-1)-shift90[1], "Excentrage : "+str(round(s*ray,3))+ "mm")
plt.text(s*np.mean(XBG90)-shift90[0],s*np.mean(YBG90)-shift90[1],"Incertitude : "+str(round(s*np.sqrt(np.std(XBG90)**2+np.std(YBG90)**2),3))+ "mm")
plt.title('G=90, C variable')
plt.xlabel('-Z, -v, mm')
plt.ylabel('X, u, mm')
plt.xticks(np.linspace(-1,1,9))
plt.grid(linestyle='--')
plt.gca().set_aspect('equal', adjustable='box')

plt.subplot(235)
shift180=(np.mean(s*XBG180+PMv[1]),np.mean(s*YBG180-PMv[0]))
ray=s*np.sqrt(np.max((x180v-XFG180v)**2+(y180v-YFG180v)**2))
circle_f=plt.Circle((s*x180v-shift180[0],s*y180v-shift180[1]),radius=ray,color='red',fill=False)
plt.scatter(s*XFG180v-shift180[0],s*YFG180v-shift180[1],color='red')
plt.scatter(s*XBG180+PMv[1]-shift180[0],s*YBG180-PMv[0]-shift180[1],color='orange')
plt.scatter(s*x180v-shift180[0],s*y180v-shift180[1],color='brown')
plt.xlim(-1,1)
plt.ylim(-1,1)
fig2 = plt.gcf()
ax2 = fig2.gca()
ax2.add_artist(circle_f)
plt.text(s*(x180v-1)-shift180[0],s*(y180v-1)-shift180[1], "Excentrage : "+str(round(s*ray,3))+ "mm")
plt.text(s*np.mean(XBG180)-shift180[0],s*np.mean(YBG180)-shift180[1],"Incertitude : "+str(round(0.224*np.sqrt(np.std(XBG180)**2+np.std(YBG180)**2),3))+ "mm")
plt.title('G=180, C variable')
plt.xlabel('-Y, -v, mm')
plt.ylabel('X, u, mm')
plt.xticks(np.linspace(-1,1,9))
plt.grid(linestyle='--')
plt.gca().set_aspect('equal', adjustable='box')

plt.subplot(231)
shift=(np.mean(s*XFG0Tsv-PMv[1]),np.mean(s*YFG0Tsv+PMv[0]))
plt.scatter(s*XFG0Tsv-PMv[1]-shift[0],s*YFG0Tsv+PMv[0]-shift[1],color='red')
plt.scatter(s*XBG0Ts-shift[0],s*YBG0Ts-shift[1],color='orange')
plt.scatter(s*x0tv-PMv[1]-shift[0],s*y0tv+PMv[0]-shift[1],color='brown')
plt.xlim(-1,1)
plt.ylim(-1,1)
plt.title('G=0, T variable')
plt.xlabel('mm')
plt.ylabel('mm')
plt.xticks(np.linspace(-1,1,9))
plt.grid(linestyle='--')
plt.gca().set_aspect('equal', adjustable='box')

plt.subplot(236)
shift270=(np.mean(s*XBG270-PMv[2]),np.mean(s*YBG270-PMv[0]))
ray=s*np.sqrt(np.max((x270v-XFG270v)**2+(y270v-YFG270v)**2))
circle_f=plt.Circle((s*x270v-shift270[0],s*y270v-shift270[1]),radius=ray,color='red',fill=False)
plt.scatter(s*XFG270v-shift270[0],s*YFG270v-shift270[1],color='red')
plt.scatter(s*XBG270-PMv[2]-shift270[0],s*YBG270-PMv[0]-shift270[1],color='orange')
plt.scatter(s*x270v-shift270[0],s*y270v-shift270[1],color='brown')
plt.xlim(-1,1)
plt.ylim(-1,1)
fig3 = plt.gcf()
ax3 = fig3.gca()
ax3.add_artist(circle_f)
plt.text(s*(x270v-1)-shift270[0],s*(y270v-1)-shift270[1], "Excentrage : "+str(round(s*ray,3))+ "mm")
plt.text(s*np.mean(XBG270)-shift270[0],s*np.mean(YBG270)-shift270[1],"Incertitude : "+str(round(0.224*np.sqrt(np.std(XBG270)**2+np.std(YBG270)**2),3))+ "mm")
plt.title('G=270, C variable')
plt.xlabel('Z, -v, mm')
plt.ylabel('X, u, mm')
plt.xticks(np.linspace(-1,1,9))
plt.grid(linestyle='--')
plt.gca().set_aspect('equal', adjustable='box')


# Graphes AB GT en fonction des excentrages relevés par RIT

plt.figure()
x=np.append(np.degrees(theta),360.)
plt.scatter(x, AB90r, color='blue')
plt.scatter(x, AB270r,color='blue')
plt.scatter(x, GT90r, color='lime')
plt.scatter(x, GT270r, color='lime')
plt.scatter(x, np.sqrt(GT90r**2+AB90r**2), color='red')
plt.scatter(x, np.sqrt(GT270r**2+AB270r**2), color='red')
tAB90r = interpolate.splrep(x, AB90r, s=0)
tAB270r = interpolate.splrep(x, AB270r, s=0)
tGT90r = interpolate.splrep(x, GT90r, s=0)
tGT270r = interpolate.splrep(x, GT270r, s=0)
xnew = np.arange(0, 360, 0.5)
yAB90r = interpolate.splev(xnew, tAB90r, der=0)
yAB270r = interpolate.splev(xnew, tAB270r, der=0)
yGT90r = interpolate.splev(xnew, tGT90r, der=0)
yGT270r = interpolate.splev(xnew, tGT270r, der=0)
plt.plot(xnew,yAB270r,color='blue', linestyle='-',label='AB270 RIT')
plt.plot(xnew,yAB90r,color='blue', linestyle=':',label='AB90 RIT')
plt.plot(xnew,yGT270r,color='lime', linestyle='-',label='GT270 RIT')
plt.plot(xnew,yGT90r,color='lime', linestyle=':',label='GT90 RIT')
plt.plot(xnew,np.sqrt(yGT270r**2+yAB270r**2),color='red', linestyle='-',label='Tot270 RIT')
plt.plot(xnew,np.sqrt(yGT90r**2+yAB90r**2),color='red', linestyle=':',label='Tot90 RIT')
plt.plot(xnew,np.ones(len(xnew)),color='black')
plt.plot(xnew,-1*np.ones(len(xnew)),color='black')
plt.ylim(-0.3,0.4)
plt.xlim(0,360)
plt.legend(loc='upper left')
plt.ylabel('mm')
plt.xlabel('Angle Bras')
plt.title("Déplacements latéraux(AB) et avant-arrière(GT) \n en fonction de l'angle bras\n avec mesures RIT")
plt.grid(linestyle='--')

#2. Projections en fonction de la table

# Il va falloir chercher à optimiser la position de l'axe de rotation de la table. Puisque PM est un point autour duquel on veut se stabiliser, on va devoir rapprocher l'axe de rotation de PM.

# Dans la série initiale, on peut déjà calculer 12 distances : les XB-XF + YB-YF.

# On déduit 8 autres distances par linéarité : celle pour tous les angles de table à bras 180 deg,
# et pour des angles colli de 90 et 270 deg. Calculons ce tableau:
    
# Explicitons les 20 distances, au format D_G_C_T:

def calcul_Dmax(AY0, AY90, AY180, AY270, AX0, AX90, AX180, AX270, AX0T, AY0T):
    
    AY180T = AY0T - AY0[0] + AY180[0]

    AY180Tp = AY0T - AY0[0] + AY180[1]

    AX180T = AX0T - AX0[0] - AX180[0]

    AX180Tp = AX0T - AX0[0] - AX180[1]
    
    D_0_270_0 = np.sqrt((AX0)**2+(AY0)**2)[0]
    D_0_90_0 = np.sqrt((AX0)**2+(AY0)**2)[1]
    D_90_90_0 = np.sqrt((AX90)**2+(AY90)**2)[0]
    D_90_270_0 = np.sqrt((AX90)**2+(AY90)**2)[1]
    D_180_270_0 = np.sqrt((AX180)**2+(AY180)**2)[0]
    D_180_90_0 = np.sqrt((AX180)**2+(AY180)**2)[1]
    D_270_90_0 = np.sqrt((AX270)**2+(AY270)**2)[0]
    D_270_270_0 = np.sqrt((AX270)**2+(AY270)**2)[1]
    D_0_270_270 = np.sqrt((AX0T)**2+(AY0T)**2)[0]
    D_0_270_315 = np.sqrt((AX0T)**2+(AY0T)**2)[1]
    D_0_270_45 = np.sqrt((AX0T)**2+(AY0T)**2)[2]
    D_0_270_90 = np.sqrt((AX0T)**2+(AY0T)**2)[3]
    
    D_180_90_270 = np.sqrt(AX180Tp[0]**2 + AY180Tp[0]**2)
    D_180_90_315 = np.sqrt(AX180Tp[1]**2 + AY180Tp[1]**2)
    D_180_90_45 = np.sqrt(AX180Tp[2]**2 + AY180Tp[2]**2)
    D_180_90_90 = np.sqrt(AX180Tp[3]**2 + AY180Tp[3]**2)
    
    D_180_270_270 = np.sqrt(AX180T[0]**2 + AY180T[0]**2)
    D_180_270_315 = np.sqrt(AX180T[1]**2 + AY180T[1]**2)
    D_180_270_45 = np.sqrt(AX180T[2]**2 + AY180T[2]**2)
    D_180_270_90 = np.sqrt(AX180T[3]**2 + AY180T[3]**2)
    
    Dmax=max(D_0_270_0,D_0_90_0,D_90_90_0,D_90_270_0,D_180_270_0,D_180_90_0,D_270_90_0,D_270_270_0,D_0_270_270,D_0_270_315,D_0_270_45,D_0_270_90,D_180_90_270,D_180_90_315,D_180_90_45,D_180_90_90,D_180_270_270,D_180_270_315,D_180_270_45,D_180_270_90)
    
    return Dmax

Dmax = calcul_Dmax(AY0, AY90, AY180, AY270, AX0, AX90, AX180, AX270, AX0T, AY0T)
Dmaxr = calcul_Dmax(AY0r, AY90r, AY180r, AY270r, AX0r, AX90r, AX180r, AX270r, AX0Tr, AY0Tr)
Dmaxv = calcul_Dmax(AY0v, AY90v, AY180v, AY270v, AX0v, AX90v, AX180v, AX270v, AX0Tv, AY0Tv)
# Vérification sur schéma et test avec mes valeurs : ok

def centre_rot_table(AX0T, AY0T):  # Least-Squares Circle Fit de Randy Bullock
    
    x=AX0T
    y=AY0T
    # coordinates of the barycenter
    x_m = np.mean(x)
    y_m = np.mean(y)
    
    # calculation of the reduced coordinates
    u = x - x_m
    v = y - y_m
    
    # linear system defining the center (uc, vc) in reduced coordinates:
    #    Suu * uc +  Suv * vc = (Suuu + Suvv)/2
    #    Suv * uc +  Svv * vc = (Suuv + Svvv)/2
    Suv  = sum(u*v)
    Suu  = sum(u**2)
    Svv  = sum(v**2)
    Suuv = sum(u**2 * v)
    Suvv = sum(u * v**2)
    Suuu = sum(u**3)
    Svvv = sum(v**3)
    
    # Solving the linear system
    A = np.array([ [ Suu, Suv ], [Suv, Svv]])
    B = np.array([ Suuu + Suvv, Svvv + Suuv ])/2.0
    uc, vc = np.linalg.solve(A, B)
    
    xc_1 = x_m + uc
    yc_1 = y_m + vc
    
    # Calcul des distances au centre (xc_1, yc_1)
    Ri_1     = np.sqrt((x-xc_1)**2 + (y-yc_1)**2)
    R_1      = np.mean(Ri_1)
    
    return (R_1,Ri_1,xc_1,yc_1)

def plot_controle(Y0, Y90, Y180, Y270, X0, X90, X180, X270, X0T, Y0T):
    
    plt.figure()
    
    plt.subplot(234)
    plt.text(0,0,"Après corr bille et table de la position de la bille")
    plt.xlim(-1,5)
    plt.ylim(-0.5,0.5)
    plt.xticks([])
    plt.yticks([])
    
    plt.subplot(232)
    ray=s*np.sqrt(np.max((x0-XFG0)**2+(y0-YFG0)**2))
    circle_f=plt.Circle((-np.mean(X0),-np.mean(Y0)),radius=ray,color='red',fill=False)
    plt.scatter(-X0,-Y0,color='red',label='Champ')
    plt.scatter(0,0,color='orange',label='Bille')
    plt.scatter(-np.mean(X0),-np.mean(Y0),color='brown',label='Centroïde Champ')
    plt.xlim(-1,1)
    plt.ylim(-1,1)
    fig2 = plt.gcf()
    ax2 = fig2.gca()
    ax2.add_artist(circle_f)
    plt.text(-0.7,-0.7, "Excentrage : "+str(round(s*ray,3))+ "mm")
    plt.title('G=0, C variable')
    plt.legend(loc='upper left')
    plt.xticks(np.linspace(-1,1,9))
    plt.grid(linestyle='--')
    plt.xlabel('Y,-v, mm')
    plt.ylabel('X, u, mm')
    plt.gca().set_aspect('equal', adjustable='box')
    
    plt.subplot(233)
    ray=s*np.sqrt(np.max((x90-XFG90)**2+(y90-YFG90)**2))
    circle_f=plt.Circle((-np.mean(X90),-np.mean(Y90)),radius=ray,color='red',fill=False)
    plt.scatter(-X90,-Y90,color='red',label='Champ')
    plt.scatter(0,0,color='orange',label='Bille')
    plt.scatter(-np.mean(X90),-np.mean(Y90),color='brown',label='Centroïde Champ')
    plt.xlim(-1,1)
    plt.ylim(-1,1)
    fig2 = plt.gcf()
    ax2 = fig2.gca()
    ax2.add_artist(circle_f)
    plt.text(-0.7,-0.7, "Excentrage : "+str(round(s*ray,3))+ "mm")
    plt.title('G=90, C variable')
    plt.legend(loc='upper left')
    plt.xticks(np.linspace(-1,1,9))
    plt.grid(linestyle='--')
    plt.xlabel('-Z,-v, mm')
    plt.ylabel('X, u, mm')
    plt.gca().set_aspect('equal', adjustable='box')
    
    plt.subplot(235)
    ray=s*np.sqrt(np.max((x180-XFG180)**2+(y180-YFG180)**2))
    circle_f=plt.Circle((-np.mean(X180),-np.mean(Y180)),radius=ray,color='red',fill=False)
    plt.scatter(-X180,-Y180,color='red',label='Champ')
    plt.scatter(0,0,color='orange',label='Bille')
    plt.scatter(-np.mean(X180),-np.mean(Y180),color='brown',label='Centroïde Champ')
    plt.xlim(-1,1)
    plt.ylim(-1,1)
    fig2 = plt.gcf()
    ax2 = fig2.gca()
    ax2.add_artist(circle_f)
    plt.text(-0.7,-0.7, "Excentrage : "+str(round(s*ray,3))+ "mm")
    plt.title('G=180, C variable')
    plt.legend(loc='upper left')
    plt.xticks(np.linspace(-1,1,9))
    plt.grid(linestyle='--')
    plt.xlabel('-Y,-v, mm')
    plt.ylabel('X, u, mm')
    plt.gca().set_aspect('equal', adjustable='box')
    
    plt.subplot(231)
    shift=(np.mean(s*XFG0Ts),np.mean(s*YFG0Ts))
    plt.scatter(s*XFG0Ts-shift[0],s*YFG0Ts-shift[1],color='red')
    plt.scatter(X0T[0],Y0T[0],color='orangered',label='T270')
    plt.scatter(X0T[1],Y0T[1],color='orange',label='T315')
    plt.scatter(X0T[2],Y0T[2],color='chartreuse',label='T45')
    plt.scatter(X0T[3],Y0T[3],color='turquoise',label='T90')
    plt.scatter(X0T[4],Y0T[4],color='yellow',label='T0')
    plt.scatter(0,0,color='brown')
    plt.xlim(-1,1)
    plt.ylim(-1,1)
    plt.title('G=0, T variable')
    plt.xlabel('mm')
    plt.ylabel('mm')
    plt.legend()
    plt.xticks(np.linspace(-1,1,9))
    plt.grid(linestyle='--')
    plt.gca().set_aspect('equal', adjustable='box')
    
    plt.subplot(236)
    ray=s*np.sqrt(np.max((x270-XFG270)**2+(y270-YFG270)**2))
    circle_f=plt.Circle((-np.mean(X270),-np.mean(Y270)),radius=ray,color='red',fill=False)
    plt.scatter(-X270,-Y270,color='red',label='Champ')
    plt.scatter(0,0,color='orange',label='Bille')
    plt.scatter(-np.mean(X270),-np.mean(Y270),color='brown',label='Centroïde Champ')
    plt.xlim(-1,1)
    plt.ylim(-1,1)
    fig2 = plt.gcf()
    ax2 = fig2.gca()
    ax2.add_artist(circle_f)
    plt.text(-0.7,-0.7, "Excentrage : "+str(round(s*ray,3))+ "mm")
    plt.title('G=270, C variable')
    plt.legend(loc='upper left')
    plt.xticks(np.linspace(-1,1,9))
    plt.grid(linestyle='--')
    plt.xlabel('Z,-v, mm')
    plt.ylabel('X, u, mm')
    plt.gca().set_aspect('equal', adjustable='box')
    
def optimisation_2(AY0, AY90, AY180, AY270, AX0, AX90, AX180, AX270, AX0T, AY0T):
    
    Lis_dis=np.zeros((101,101))
    Lis_PM=np.zeros((101,101),dtype=list)
    Lis_ST=np.zeros((101,101),dtype=list)
    
    
    
    for i in range (101):
        for j in range (101):
            
            Y0=np.copy(AY0)
            Y90=np.copy(AY90)
            Y180=np.copy(AY180)
            Y270=np.copy(AY270)
            Y0T=np.copy(AY0T)
            X0=np.copy(AX0)
            X90=np.copy(AX90)
            X180=np.copy(AX180)
            X270=np.copy(AX270)
            X0T=np.copy(AX0T)
            
            ST = np.array([(i-50)*0.01, (j-50)*0.01])
            
            theta0 = np.arctan2(ST[1],ST[0])
            
            R = np.sqrt(ST[1]**2+ST[0]**2)
            
            X0=X0-R*np.cos(theta0)
            Y0=Y0-R*np.sin(theta0)
            X90=X90
            Y90=Y90-R*np.sin(theta0)
            X180=X180+R*np.cos(theta0)
            Y180=Y180-R*np.sin(theta0)
            X270=X270
            Y270=Y270-R*np.sin(theta0)
            X0T=X0T-R*np.array([np.cos(theta0+np.radians(+90)),np.cos(theta0+np.radians(+45)),np.cos(theta0+np.radians(-45)),np.cos(theta0+np.radians(-90)),np.cos(theta0)])
            Y0T=Y0T-R*np.array([np.sin(theta0+np.radians(+90)),np.sin(theta0+np.radians(+45)),np.sin(theta0+np.radians(-45)),np.sin(theta0+np.radians(-90)),np.sin(theta0)])
            
            M,C,G = calcul_PM(Y0, Y90, Y180, Y270, X0, X90, X180, X270)
            
            X0=X0-M[1]
            Y0=Y0-M[0]
            X90=X90+M[2]
            Y90=Y90-M[0]
            X180=X180+M[1]
            Y180=Y180-M[0]
            X270=X270-M[2]
            Y270=Y270-M[0]
            X0T=X0T-M[1]
            Y0T=Y0T-M[0]
            
            D = calcul_Dmax(Y0, Y90, Y180, Y270, X0, X90, X180, X270, X0T, Y0T)
            
            Lis_dis[i,j]=D
            Lis_PM[i,j]=M
            Lis_ST[i,j]=ST
            
            
    plt.figure()
    plt.imshow(Lis_dis,origin='lower')
    
    return
    
# optimisation_2(AY0, AY90, AY180, AY270, AX0, AX90, AX180, AX270, AX0T, AY0T)

def optimisation_3(AY0, AY90, AY180, AY270, AX0, AX90, AX180, AX270, AX0T, AY0T, PM, Dmax, epsilon):
    
    
    
    Y0=np.copy(AY0)
    Y90=np.copy(AY90)
    Y180=np.copy(AY180)
    Y270=np.copy(AY270)
    Y0T=np.copy(AY0T)
    X0=np.copy(AX0)
    X90=np.copy(AX90)
    X180=np.copy(AX180)
    X270=np.copy(AX270)
    X0T=np.copy(AX0T)
    
    
    X0=X0-PM[1]
    Y0=Y0-PM[0]
    X90=X90+PM[2]
    Y90=Y90-PM[0]
    X180=X180+PM[1]
    Y180=Y180-PM[0]
    X270=X270-PM[2]
    Y270=Y270-PM[0]
    X0T=X0T-PM[1]
    Y0T=Y0T-PM[0]
    
    SB=np.copy(PM)
    ST=np.array([0.,0.])
    
    D = calcul_Dmax(AY0, AY90, AY180, AY270, AX0, AX90, AX180, AX270, AX0T, AY0T)

    x=[]
    y=[]
    
    STX=[50]
    STY=[50]
    i=0
    
    loop=True
    def gradT(Y0, Y90, Y180, Y270, X0, X90, X180, X270, X0T, Y0T, D, xTrueoryFalse):
        
        if xTrueoryFalse:
            
            X0=X0+0.005
            X180=X180-0.005
            X0T=X0T+np.array([0,0.0025*np.sqrt(2),-0.0025*np.sqrt(2),0,0.005])
            Y0T=Y0T+np.array([0.005,0.0025*np.sqrt(2),-0.0025*np.sqrt(2),-0.005,0])
            
            D1 = calcul_Dmax(Y0, Y90, Y180, Y270, X0, X90, X180, X270, X0T, Y0T)
        
        else:
        
            Y0=Y0+0.005
            Y90=Y90+0.005
            Y180=Y180+0.005
            Y270=Y270+0.005
            X0T=X0T+np.array([-0.005,-0.0025*np.sqrt(2),0.0025*np.sqrt(2),0.005,0])
            Y0T=Y0T+np.array([0,0.0025*np.sqrt(2),0.0025*np.sqrt(2),0,0.005])
    
            D1 = calcul_Dmax(Y0, Y90, Y180, Y270, X0, X90, X180, X270, X0T, Y0T)
        return ((D1-D)/0.005)
    
    while(loop and i<100):
        
        x.append(i)
        y.append(D)
       
        G=(gradT(Y0, Y90, Y180, Y270, X0, X90, X180, X270, X0T, Y0T, D, True),gradT(Y0, Y90, Y180, Y270, X0, X90, X180, X270, X0T, Y0T, D, False))
        
        step = np.sqrt(G[0]**2+G[1]**2)
        if np.sqrt(G[0]**2+G[1]**2)<epsilon:
            loop=False
            
        else:
            
            theta0=np.pi+np.arctan2(G[1],G[0])
            
       
            ST+=step*0.005*np.array([np.cos(theta0),np.sin(theta0)])
    
            X0=X0-step*0.005*np.cos(theta0)
            Y0=Y0-step*0.005*np.sin(theta0)
            X90=X90
            Y90=Y90-step*0.005*np.sin(theta0)
            X180=X180+step*0.005*np.cos(theta0)
            Y180=Y180-step*0.005*np.sin(theta0)
            X270=X270
            Y270=Y270-step*0.005*np.sin(theta0)
            X0T=X0T+step*0.005*np.array([np.cos(theta0+np.radians(+90)),np.cos(theta0+np.radians(+45)),np.cos(theta0+np.radians(-45)),np.cos(theta0+np.radians(-90)),np.cos(theta0)])
            Y0T=Y0T+step*0.005*np.array([np.sin(theta0+np.radians(+90)),np.sin(theta0+np.radians(+45)),np.sin(theta0+np.radians(-45)),np.sin(theta0+np.radians(-90)),np.sin(theta0)])
            
            M,C,G = calcul_PM(Y0, Y90, Y180, Y270, X0, X90, X180, X270)
            
            X0=X0-M[1]
            Y0=Y0-M[0]
            X90=X90+M[2]
            Y90=Y90-M[0]
            X180=X180+M[1]
            Y180=Y180-M[0]
            X270=X270-M[2]
            Y270=Y270-M[0]
            X0T=X0T-M[1]
            Y0T=Y0T-M[0]
        
            SB = SB + M
            
            D = calcul_Dmax(Y0, Y90, Y180, Y270, X0, X90, X180, X270, X0T, Y0T)
            i+=1
            STX.append(ST[0]*100 + 50)
            STY.append(ST[1]*100 + 50)
    
    plt.figure()
    plt.plot(x,y)
    plot_controle(Y0, Y90, Y180, Y270, X0, X90, X180, X270, X0T, Y0T)
    plt.figure()
    plt.plot(STX, STY)
    plt.xlim(0,100)
    plt.ylim(0,100)   
    plt.gca().set_aspect('equal', adjustable='box')     
    return(SB, ST, Y0, Y90, Y180, Y270, X0, X90, X180, X270, X0T, Y0T,D)

# (SB, ST, Y0, Y90, Y180, Y270, X0, X90, X180, X270, X0T, Y0T,D)=optimisation_3(AY0, AY90, AY180, AY270, AX0, AX90, AX180, AX270, AX0T, AY0T, PM, Dmax, 0.001)

def optimisation_table(AY0, AY90, AY180, AY270, AX0, AX90, AX180, AX270, AX0T, AY0T, PM, epsilon=0,step=0.02):
    
    D0 = calcul_Dmax(AY0, AY90, AY180, AY270, AX0, AX90, AX180, AX270, AX0T, AY0T)
    
    Y0=np.copy(AY0)
    Y90=np.copy(AY90)
    Y180=np.copy(AY180)
    Y270=np.copy(AY270)
    Y0T=np.copy(AY0T)
    X0=np.copy(AX0)
    X90=np.copy(AX90)
    X180=np.copy(AX180)
    X270=np.copy(AX270)
    X0T=np.copy(AX0T)
    
    
    X0=X0-PM[1]
    Y0=Y0-PM[0]
    X90=X90+PM[2]
    Y90=Y90-PM[0]
    X180=X180+PM[1]
    Y180=Y180-PM[0]
    X270=X270-PM[2]
    Y270=Y270-PM[0]
    X0T=X0T-PM[1]
    Y0T=Y0T-PM[0]
    
    R,Ri,xc,yc=centre_rot_table(X0T, Y0T)
    
    #On détermine les erreurs systématiques liées à l'imperfection de la rotation de table.
    # On compare donc, en fonction des rayons mesurés;, la position des billes.
    # L'angle 0 théorique est alors confondu avec l'angle 0 de la table, car c'est
    # l'angle de référence en traitement.
    
    # d0=[(X0T[4]-xc)*(1-R/Ri[4]), (Y0T[4]-yc)*(1-R/Ri[4])]
    
    # l'angle de référence, l'angle 0, est donné par :
        
    theta0=np.arctan2((Y0T[4]-yc),X0T[4]-xc)
    
    # d90=[X0T[3]-(xc+R*np.cos(theta0+np.radians(-90))), Y0T[3]-(yc+R*np.sin(theta0+np.radians(-90)))]
    # d45=[X0T[2]-(xc+R*np.cos(theta0+np.radians(-45))), Y0T[2]-(yc+R*np.sin(theta0+np.radians(-45)))]
    # d315=[X0T[1]-(xc+R*np.cos(theta0+np.radians(+45))), Y0T[1]-(yc+R*np.sin(theta0+np.radians(45)))]
    # d270=[X0T[0]-(xc+R*np.cos(theta0+np.radians(+90))), Y0T[0]-(yc+R*np.sin(theta0+np.radians(+90)))]

    # Ces décalages systématiques seront toujours appliqués. Ils sont en fait inutiles
    # apriori, mais je les garde au cas ou j'en aurais besoin pour remanier le code.
    
    
    SB=np.copy(PM) # Shift bille
    
    ST = step*np.array([-Ri[4]*np.cos(theta0), -Ri[4]*np.sin(theta0)])
    
    X0=X0-step*Ri[4]*np.cos(theta0)
    Y0=Y0-step*Ri[4]*np.sin(theta0)
    X90=X90
    Y90=Y90-step*Ri[4]*np.sin(theta0)
    X180=X180+step*Ri[4]*np.cos(theta0)
    Y180=Y180-step*Ri[4]*np.sin(theta0)
    X270=X270
    Y270=Y270-step*Ri[4]*np.sin(theta0)
    X0T=X0T+step*np.array([-Ri[4]*np.cos(theta0+np.radians(+90)),-Ri[4]*np.cos(theta0+np.radians(+45)),-Ri[4]*np.cos(theta0+np.radians(-45)),-Ri[4]*np.cos(theta0+np.radians(-90)),-Ri[4]*np.cos(theta0)])
    Y0T=Y0T+step*np.array([-Ri[4]*np.sin(theta0+np.radians(+90)),-Ri[4]*np.sin(theta0+np.radians(+45)),-Ri[4]*np.sin(theta0+np.radians(-45)),-Ri[4]*np.sin(theta0+np.radians(-90)),-Ri[4]*np.sin(theta0)])
    
    D1 = calcul_Dmax(Y0, Y90, Y180, Y270, X0, X90, X180, X270, X0T, Y0T)
    
    x = np.arange(0, 202, 1)
    y=np.zeros(202)
    
    y[0]=D0
    y[1]=D1
    STX=[50]
    STY=[50]
    n=0
    while (n<200 ): #and abs((D1-D0)/D1)>2*epsilon*step
        
        D0=D1
        n+=1
        
        M,C,G = calcul_PM(Y0, Y90, Y180, Y270, X0, X90, X180, X270)
        
        X0=X0-M[1]
        Y0=Y0-M[0]
        X90=X90+M[2]
        Y90=Y90-M[0]
        X180=X180+M[1]
        Y180=Y180-M[0]
        X270=X270-M[2]
        Y270=Y270-M[0]
        X0T=X0T-M[1]
        Y0T=Y0T-M[0]
        
        SB = SB + M
        D0 = D1
        D1 = calcul_Dmax(Y0, Y90, Y180, Y270, X0, X90, X180, X270, X0T, Y0T)
        
        R,Ri,xc,yc=centre_rot_table(X0T, Y0T)
        
        theta0=np.arctan2((Y0T[4]-yc),X0T[4]-xc)
        
        X0=X0-step*Ri[4]*np.cos(theta0)
        Y0=Y0-step*Ri[4]*np.sin(theta0)
        X90=X90
        Y90=Y90-step*Ri[4]*np.sin(theta0)
        X180=X180+step*Ri[4]*np.cos(theta0)
        Y180=Y180-step*Ri[4]*np.sin(theta0)
        X270=X270
        Y270=Y270-step*Ri[4]*np.sin(theta0)
        X0T=X0T+step*np.array([-Ri[4]*np.cos(theta0+np.radians(+90)),-Ri[4]*np.cos(theta0+np.radians(+45)),-Ri[4]*np.cos(theta0+np.radians(-45)),-Ri[4]*np.cos(theta0+np.radians(-90)),-Ri[4]*np.cos(theta0)])
        Y0T=Y0T+step*np.array([-Ri[4]*np.sin(theta0+np.radians(+90)),-Ri[4]*np.sin(theta0+np.radians(+45)),-Ri[4]*np.sin(theta0+np.radians(-45)),-Ri[4]*np.sin(theta0+np.radians(-90)),-Ri[4]*np.sin(theta0)])
        
        SB+=M
        ST+=step*np.array([-Ri[4]*np.cos(theta0), -Ri[4]*np.sin(theta0)])
        
        D1 = calcul_Dmax(Y0, Y90, Y180, Y270, X0, X90, X180, X270, X0T, Y0T)
        y[n+1]=D1
        
        STX.append(ST[0]*100 + 50)
        STY.append(ST[1]*100 + 50)
        
    plt.figure()
    plt.plot(x,y)
    plot_controle(Y0, Y90, Y180, Y270, X0, X90, X180, X270, X0T, Y0T)
    plt.figure()
    plt.plot(STX,STY)
    plt.xlim(0,100)
    plt.ylim(0,100)
    plt.gca().set_aspect('equal', adjustable='box')
        
    return(SB, ST, Y0, Y90, Y180, Y270, X0, X90, X180, X270, X0T, Y0T, D1)
        
# (SB, ST, Y0, Y90, Y180, Y270, X0, X90, X180, X270, X0T, Y0T) = optimisation_table(AY0, AY90, AY180, AY270, AX0, AX90, AX180, AX270, AX0T, AY0T, PM, epsilon=0.01)
    # Ensuite, on décale virtuellement l'isocentre de table sur le point mesuré à table 0.
# (SBv, STv, Y0v, Y90v, Y180v, Y270v, X0v, X90v, X180v, X270v, X0Tv, Y0Tv) = optimisation_table(AY0v, AY90v, AY180v, AY270v, AX0v, AX90v, AX180v, AX270v, AX0Tv, AY0Tv, PMv, epsilon=0.01)
    # Cela revient à pousser la table pour que le point 0 


#################
# YXr=np.concatenate((AX0r,AX90r,AX180r,AX270r,AX0Tr))

# YX=np.concatenate((AX0,AX90,AX180,AX270,AX0T))

# YY=np.concatenate((AY0,AY90,AY180,AY270,AY0T))

# YYr=np.concatenate((AY0r,AY90r,AY180r,AY270r,AY0Tr))

# x = np.arange(1, 14, 1)

# plt.figure()

# plt.plot(x,YX,label='Horizontal mesuré')
# plt.plot(x,YY,label='Vertical mesuré')
# plt.plot(x,YXr, label='Horizontal RIT')
# plt.plot(x,YYr, label='Vertical RIT')
# plt.legend()