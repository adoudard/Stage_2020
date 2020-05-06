# -*- coding: utf-8 -*-

import os
import matplotlib.pyplot as plt
import numpy as np
from scipy import interpolate

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
        

# dossier ou se trouve le fichier txt
os.chdir('E:\\Stage 2020\\Scripts')

# Ouverture en lecture seule
analyse=open("30_04.txt","r")
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


s=(2/3)*0.336 #scale


# Système de coordonnées : 
    # z = axe de la table
    # y = axe antéro postérieur
    # x = axe pour que x,y,z soit othogonal et orienté positivement.
    

AX0=s*(XBG0-XFG0)
AY0=s*(YBG0-YFG0)

AX90=s*(XBG90-XFG90)
AY90=s*(YBG90-YFG90)

AX180=s*(XBG180-XFG180)
AY180=s*(YBG180-YFG180)

AX270=s*(XBG270-XFG270)
AY270=s*(YBG270-YFG270)

AX0T=s*(XBG0Ts-XFG0Ts)
AY0T=s*(YBG0Ts-YFG0Ts)

angles=[(0,0),(90,0),(180,0),(270,0),(0,-90),(0,-45),(0,45),(0,90)]

A=np.zeros((2*len(angles),3))

for i in range (len(angles)):
    
    A[2*i]=[round(np.cos(np.radians(angles[i][1])),5), -1*round(np.sin(np.radians(angles[i][1])),5), 0]
    A[2*i+1]=[round(np.cos(np.radians(angles[i][0]))*np.sin(np.radians(angles[i][1])),5), round(np.cos(np.radians(angles[i][0]))*np.cos(np.radians(angles[i][1])),5), round(np.sin(np.radians(angles[i][0])),5)]

B=np.dot(np.linalg.inv(np.dot(A.transpose(),A)),A.transpose())

B0=(np.mean(AX0),np.mean(AY0))
B90=(np.mean(AX90),np.mean(AY90))
B180=(np.mean(AX180),np.mean(AY180))
B270=(np.mean(AX270),np.mean(AY270))

ksi=np.array([B0[0],B0[1],B90[0],B90[1],B180[0],B180[1],B270[0],B270[1],AX0T[0],AY0T[0],AX0T[1],AY0T[1],AX0T[2],AY0T[2],AX0T[3],AY0T[3]])

delta=np.dot(B,ksi.transpose())

print("Optimisation de Low : La bille est excentrée de dX = ",-1*round(delta[1],3)," mm, dY = ",-1*round(delta[0],3), " mm et de dZ = ", -1* round(delta[2],3)," mm dans le système de coorodonnées Varian IEC 1217.")
              
plt.figure()

plt.subplot(232)
ray=np.sqrt(np.max((x0-XFG0)**2+(y0-YFG0)**2))
circle_f=plt.Circle((x0,y0),radius=ray,color='red',fill=False)
plt.scatter(XFG0,YFG0,color='red',label='Champ')
plt.scatter(XBG0,YBG0,color='orange',label='Bille')
plt.scatter(x0,y0,color='brown',label='Centroïde Champ')
plt.xlim(636,644)
plt.ylim(636,644)
fig1 = plt.gcf()
ax1 = fig1.gca()
ax1.add_artist(circle_f)
plt.text(x0-1,y0-1, "Excentrage : "+str(round(0.224*ray,3))+ "mm")
plt.text(np.mean(XBG0),np.mean(YBG0),"Incertitude : "+str(round(0.224*np.sqrt(np.std(XBG0)**2+np.std(YBG0)**2),3))+ "mm")
plt.title('G=0, C variable')
plt.legend(loc='upper left')
plt.gca().set_aspect('equal', adjustable='box')

plt.subplot(233)
ray=np.sqrt(np.max((x90-XFG90)**2+(y90-YFG90)**2))
circle_f=plt.Circle((x90,y90),radius=ray,color='red',fill=False)
plt.scatter(XFG90,YFG90,color='red')
plt.scatter(XBG90,YBG90,color='orange')
plt.scatter(x90,y90,color='brown')
plt.xlim(636,644)
plt.ylim(636,644)
fig4 = plt.gcf()
ax4 = fig4.gca()
ax4.add_artist(circle_f)
plt.text(x90-1,y90-1, "Excentrage : "+str(round(0.224*ray,3))+ "mm")
plt.text(np.mean(XBG90),np.mean(YBG90),"Incertitude : "+str(round(0.224*np.sqrt(np.std(XBG90)**2+np.std(YBG90)**2),3))+ "mm")
plt.title('G=90, C variable')
plt.gca().set_aspect('equal', adjustable='box')

plt.subplot(235)
ray=np.sqrt(np.max((x180-XFG180)**2+(y180-YFG180)**2))
circle_f=plt.Circle((x180,y180),radius=ray,color='red',fill=False)
plt.scatter(XFG180,YFG180,color='red')
plt.scatter(XBG180,YBG180,color='orange')
plt.scatter(x180,y180,color='brown')
plt.xlim(636,644)
plt.ylim(636,644)
fig2 = plt.gcf()
ax2 = fig2.gca()
ax2.add_artist(circle_f)
plt.text(x180-1,y180-1, "Excentrage : "+str(round(0.224*ray,3))+ "mm")
plt.text(np.mean(XBG180),np.mean(YBG180),"Incertitude : "+str(round(0.224*np.sqrt(np.std(XBG180)**2+np.std(YBG180)**2),3))+ "mm")
plt.title('G=180, C variable')
plt.gca().set_aspect('equal', adjustable='box')

plt.subplot(231)
plt.scatter(XFG0Ts,YFG0Ts,color='red')
plt.scatter(XBG0Ts,YBG0Ts,color='orange')
plt.scatter(x0t,y0t,color='brown')
plt.xlim(636,644)
plt.ylim(636,644)
plt.title('G=0, T variable')
plt.gca().set_aspect('equal', adjustable='box')

plt.subplot(236)
ray=np.sqrt(np.max((x270-XFG270)**2+(y270-YFG270)**2))
circle_f=plt.Circle((x270,y270),radius=ray,color='red',fill=False)
plt.scatter(XFG270,YFG270,color='red')
plt.scatter(XBG270,YBG270,color='orange')
plt.scatter(x270,y270,color='brown')
plt.xlim(636,644)
plt.ylim(636,644)
fig3 = plt.gcf()
ax3 = fig3.gca()
ax3.add_artist(circle_f)
plt.text(x270-1,y270-1, "Excentrage : "+str(round(0.224*ray,3))+ "mm")
plt.text(np.mean(XBG270),np.mean(YBG270),"Incertitude : "+str(round(0.224*np.sqrt(np.std(XBG270)**2+np.std(YBG270)**2),3))+ "mm")
plt.title('G=270, C variable')
plt.gca().set_aspect('equal', adjustable='box')


plt.show()

plt.figure()

xs=[0,90,180,270]
ys=[s*(np.mean(XBG0)-x0),s*(np.mean(XBG90)-x90),s*(np.mean(XBG180)-x180),s*(np.mean(XBG270)-x270)]
plt.scatter(xs,ys)


##Methode d'optimisation de Hancock##

#1.1 Isocentre bras avec C=90,T=0

ug=np.array([AY0[1],AY90[0],AY180[1],AY270[0]])
vg=np.array([AX0[1],AX90[0],AX180[1],AX270[0]])

theta=np.array([np.radians(0),np.radians(90),np.radians(180),np.radians(270)])

xg=ug
yg=-vg*np.cos(theta)
zg=vg*np.sin(theta)

PG=0.5*np.array([np.max(xg)+np.min(xg), np.max(yg)+np.min(yg), np.max(zg)+np.min(zg)])

#1.2 Isocentre bras avec C=270,T=0

uc=np.array([AY0[0],AY90[1],AY180[0],AY270[1]])
vc=np.array([AX0[0],AX90[1],AX180[0],AX270[1]])

xc=uc
yc=-vc*np.cos(theta)
zc=vc*np.sin(theta)

PC=0.5*np.array([np.max(xc)+np.min(xc), np.max(yc)+np.min(yc), np.max(zc)+np.min(zc)])

#1.3 Isocentre moyen

PM=0.5*(PG+PC)

# On va soustraire PM (ses projections) à la bille

GT90=np.array([AY0[1]-PM[0], AY90[0]-PM[0], AY180[1]-PM[0], AY270[0]-PM[0], AY0[1]-PM[0]])
GT270=np.array([AY0[0]-PM[0], AY90[1]-PM[0], AY180[0]-PM[0], AY270[1]-PM[0], AY0[0]-PM[0]])
AB90=np.array([AX0[1]+PM[1], AX90[0]-PM[2], AX180[1]-PM[1], AX270[0]+PM[2], AX0[1]+PM[1]]) #*-1 ???
AB270=np.array([AX0[0]+PM[1], AX90[1]-PM[2], AX180[0]-PM[1], AX270[1]+PM[2], AX0[0]+PM[1]]) #*-1 ???

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
xnew = np.arange(0, 360, 0.5)
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

plt.subplot(232)
shift0=(np.mean(s*XBG0+PM[1]),np.mean(s*YBG0-PM[0]))
ray=s*np.sqrt(np.max((x0-XFG0)**2+(y0-YFG0)**2))
circle_f=plt.Circle((s*x0-shift0[0],s*y0-shift0[1]),radius=ray,color='red',fill=False)
plt.scatter(s*XFG0-shift0[0],s*YFG0-shift0[1],color='red',label='Champ')
plt.scatter(s*XBG0+PM[1]-shift0[0],s*YBG0-PM[0]-shift0[1],color='orange',label='Bille')
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
plt.gca().set_aspect('equal', adjustable='box')

plt.subplot(233)
shift90=(np.mean(s*XBG90-PM[2]),np.mean(s*YBG90-PM[0]))
ray=s*np.sqrt(np.max((x90-XFG90)**2+(y90-YFG90)**2))
circle_f=plt.Circle((s*x90-shift90[0],s*y90-shift90[1]),radius=ray,color='red',fill=False)
plt.scatter(s*XFG90-shift90[0],s*YFG90-shift90[1],color='red')
plt.scatter(s*XBG90-PM[2]-shift90[0],s*YBG90-PM[0]-shift90[1],color='orange')
plt.scatter(s*x90-shift90[0],s*y90-shift90[1],color='brown')
plt.xlim(-1,1)
plt.ylim(-1,1)
fig2 = plt.gcf()
ax2 = fig2.gca()
ax2.add_artist(circle_f)
plt.text(s*(x90-1)-shift90[0],s*(y90-1)-shift90[1], "Excentrage : "+str(round(s*ray,3))+ "mm")
plt.text(s*np.mean(XBG90)-shift90[0],s*np.mean(YBG90)-shift90[1],"Incertitude : "+str(round(s*np.sqrt(np.std(XBG90)**2+np.std(YBG90)**2),3))+ "mm")
plt.title('G=90, C variable')
plt.gca().set_aspect('equal', adjustable='box')

plt.subplot(235)
shift180=(np.mean(s*XBG180-PM[1]),np.mean(s*YBG180-PM[0]))
ray=s*np.sqrt(np.max((x180-XFG180)**2+(y180-YFG180)**2))
circle_f=plt.Circle((s*x180-shift180[0],s*y180-shift180[1]),radius=ray,color='red',fill=False)
plt.scatter(s*XFG180-shift180[0],s*YFG180-shift180[1],color='red')
plt.scatter(s*XBG180-PM[1]-shift180[0],s*YBG180-PM[0]-shift180[1],color='orange')
plt.scatter(s*x180-shift180[0],s*y180-shift180[1],color='brown')
plt.xlim(-1,1)
plt.ylim(-1,1)
fig2 = plt.gcf()
ax2 = fig2.gca()
ax2.add_artist(circle_f)
plt.text(s*(x180-1)-shift180[0],s*(y180-1)-shift180[1], "Excentrage : "+str(round(s*ray,3))+ "mm")
plt.text(s*np.mean(XBG180)-shift180[0],s*np.mean(YBG180)-shift180[1],"Incertitude : "+str(round(0.224*np.sqrt(np.std(XBG180)**2+np.std(YBG180)**2),3))+ "mm")
plt.title('G=180, C variable')
plt.gca().set_aspect('equal', adjustable='box')

plt.subplot(231)
shift=(np.mean(s*XFG0Ts-PM[1]),np.mean(s*YFG0Ts+PM[0]))
plt.scatter(s*XFG0Ts-shift[0],s*YFG0Ts-shift[1],color='red')
plt.scatter(s*XBG0Ts+PM[1]-shift[0],s*YBG0Ts-PM[0]-shift[1],color='orange')
plt.scatter(s*x0t-shift[0],s*y0t-shift[1],color='brown')
plt.xlim(-1,1)
plt.ylim(-1,1)
plt.title('G=0, T variable')
plt.gca().set_aspect('equal', adjustable='box')

plt.subplot(236)
shift270=(np.mean(s*XBG270+PM[2]),np.mean(s*YBG270-PM[0]))
ray=s*np.sqrt(np.max((x270-XFG270)**2+(y270-YFG270)**2))
circle_f=plt.Circle((s*x270-shift270[0],s*y270-shift270[1]),radius=ray,color='red',fill=False)
plt.scatter(s*XFG270-shift270[0],s*YFG270-shift270[1],color='red')
plt.scatter(s*XBG270+PM[2]-shift270[0],s*YBG270-PM[0]-shift270[1],color='orange')
plt.scatter(s*x270-shift270[0],s*y270-shift270[1],color='brown')
plt.xlim(-1,1)
plt.ylim(-1,1)
fig3 = plt.gcf()
ax3 = fig3.gca()
ax3.add_artist(circle_f)
plt.text(s*(x270-1)-shift270[0],s*(y270-1)-shift270[1], "Excentrage : "+str(round(s*ray,3))+ "mm")
plt.text(s*np.mean(XBG270)-shift270[0],s*np.mean(YBG270)-shift270[1],"Incertitude : "+str(round(0.224*np.sqrt(np.std(XBG270)**2+np.std(YBG270)**2),3))+ "mm")
plt.title('G=270, C variable')
plt.gca().set_aspect('equal', adjustable='box')

#2. Projections en fonction de la table

