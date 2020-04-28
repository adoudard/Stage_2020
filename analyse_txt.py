# -*- coding: utf-8 -*-
"""
Created on Mon Apr 13 10:24:10 2020

@author: alex2
"""


import os
import matplotlib.pyplot as plt
import numpy as np

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
analyse=open("28043b2.txt","r")
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
XFG0=[Fs[i][0] for i in range (5)]
YFG0=[Fs[i][1] for i in range (5)]
XBG0=[Bs[i][0] for i in range (5)]
YBG0=[Bs[i][1] for i in range (5)]

# G90, T0, C variable
XFG90=[Fs[i+5][0] for i in range (5)]
YFG90=[Fs[i+5][1] for i in range (5)]
XBG90=[Bs[i+5][0] for i in range (5)]
YBG90=[Bs[i+5][1] for i in range (5)]

# G180, T0, C variable
XFG180=[Fs[i+10][0] for i in range (5)]
YFG180=[Fs[i+10][1] for i in range (5)]
XBG180=[Bs[i+10][0] for i in range (5)]
YBG180=[Bs[i+10][1] for i in range (5)]

# G270, T0, C variable
XFG270=[Fs[i+15][0] for i in range (5)]
YFG270=[Fs[i+15][1] for i in range (5)]
XBG270=[Bs[i+15][0] for i in range (5)]
YBG270=[Bs[i+15][1] for i in range (5)]

# G180E, T0, C variable
XFG180E=[Fs[i+20][0] for i in range (5)]
YFG180E=[Fs[i+20][1] for i in range (5)]
XBG180E=[Bs[i+20][0] for i in range (5)]
YBG180E=[Bs[i+20][1] for i in range (5)]

#G180, C288, T variable
XFG180Ts=[Fs[i+25][0] for i in range (4)]
YFG180Ts=[Fs[i+25][1] for i in range (4)]
XBG180Ts=[Bs[i+25][0] for i in range (4)]
YBG180Ts=[Bs[i+25][1] for i in range (4)]

# On rajoute également aux listes au dessus l'acquisition pour T=0
XFG180Ts.append(Fs[14][0])
YFG180Ts.append(Fs[14][1])
XBG180Ts.append(Bs[14][0])
YBG180Ts.append(Bs[14][1])

(x0, y0) = centre_of_mass(XFG0,YFG0)
(x90, y90) = centre_of_mass(XFG90,YFG90)
(x180, y180) = centre_of_mass(XFG180,YFG180)
(x180e, y180e) = centre_of_mass(XFG180E,YFG180E)
(x180t, y180t) = centre_of_mass(XFG180Ts,YFG180Ts)
(x270, y270) = centre_of_mass(XFG270,YFG270)

s=(2/3)*0.336 #scale

XFG0 = np.asarray(XFG0)
YFG0 = np.asarray(YFG0)
XBG0 = np.asarray(XBG0)
YBG0 = np.asarray(YBG0)

# G90, T0, C variable
XFG90 = np.asarray(XFG90)
YFG90 = np.asarray(YFG90)
XBG90 = np.asarray(XBG90)
YBG90 = np.asarray(YBG90)

# G180, T0, C variable
XFG180 = np.asarray(XFG180)
YFG180 = np.asarray(YFG180)
XBG180 = np.asarray(XBG180)
YBG180 = np.asarray(YBG180)

# G270, T0, C variable
XFG270 = np.asarray(XFG270)
YFG270 = np.asarray(YFG270)
XBG270 = np.asarray(XBG270)
YBG270 = np.asarray(YBG270)

# G180E, T0, C variable
XFG180E = np.asarray(XFG180E)
YFG180E = np.asarray(YFG180E)
XBG180E = np.asarray(XBG180E)
YBG180E = np.asarray(YBG180E)

#G180, C288, T variable
XFG180Ts = np.asarray(XFG180Ts)
YFG180Ts = np.asarray(YFG180Ts)
XBG180Ts = np.asarray(XBG180Ts)
YBG180Ts = np.asarray(YBG180Ts)

Dis0 = np.sqrt((XFG0-np.float64(x0))**2+(YFG0-np.float64(y0))**2)
Dis90 = np.sqrt((XFG90-np.float64(x90))**2+(YFG90-np.float64(y90))**2)
Dis180 = np.sqrt((XFG180-np.float64(x180))**2+(YFG180-np.float64(y180))**2)
Dis180E = np.sqrt((XFG180E-np.float64(x180e))**2+(YFG180E-np.float64(y180e))**2)
Dis270 = np.sqrt((XFG270-np.float64(x270))**2+(YFG270-np.float64(y270))**2)

Dist_champ_rot = np.mean([Dis0, Dis90, Dis180, Dis180E, Dis270])

# Dis0 = np.sqrt((np.mean(XBG0)-np.mean(XFG0))**2+(np.mean(YBG0)-np.mean(YFG0))**2)
# Dis90 = np.sqrt((np.mean(XBG90)-np.mean(XFG90))**2+(np.mean(YBG90)-np.mean(YFG90))**2)
# Dis180 = np.sqrt((np.mean(XBG180)-np.mean(XFG180))**2+(np.mean(YBG180)-np.mean(YFG180))**2)
# Dis180E = np.sqrt((np.mean(XBG180E)-np.mean(XFG180E))**2+(np.mean(YBG180E)-np.mean(YFG180E))**2)
# Dis270 = np.sqrt((np.mean(XBG270)-np.mean(XFG270))**2+(np.mean(YBG270)-np.mean(YFG270))**2)

# Dist_isocentres = np.mean([Dis0, Dis90, Dis180, Dis180E, Dis270])

# Contrib_table = np.sqrt((np.mean(XBG180Ts)-np.mean(XFG180Ts))**2+(np.mean(YBG180Ts)-np.mean(YFG180Ts))**2)

print("Ecart entre l'axe de rotation du collimateur et l'axe de la source de : ", (2/3)*0.336*Dist_champ_rot, " mm.")
# print("Ecart entre l'isocentre machine et celui d'irradiation de : ", (2/3)*0.336*Dist_isocentres, " mm.")
# print("Ecart supplémentaire lié à l'isocentre table de : ", (2/3)*0.336*Contrib_table, " mm.")


# Système de coordonnées : 
    # z = axe de la table
    # y = axe antéro postérieur
    # x = axe pour que x,y,z soit othogonal et orienté positivement.
    

A1=s*(XBG0-XFG0)
A2=s*(YBG0-YFG0)
A4=s*(YBG90-YFG90)
A3=s*(XBG90-XFG90)
A5=s*(XBG180-XFG180)
A6=s*(YBG180-YFG180)
A8=s*(YBG270-YFG270)
A7=s*(XBG270-XFG270)
A9=s*(XBG180E-XFG180E)
A10=s*(YBG180E-YFG180E)
A11=s*(YBG180Ts-YFG180Ts)
A12=s*(XBG180Ts-XFG180Ts)


angles=[(0,0),(90,0),(180,0),(270,0),(179.99,0),(180,-90),(180,-45),(180,90),(180,45)]
anglesp=[]
for i in range(len(Gs)):
    anglesp.append((Gs[i],Ts[i],Cs[i]))
    pass


A=np.zeros((2*len(angles),3))
Ap=np.zeros((2*len(anglesp),3))

for i in range (len(angles)):
    
    A[2*i]=[np.cos(np.radians(angles[i][1])), np.sin(np.radians(angles[i][1])), 0]
    A[2*i+1]=[-1*np.cos(np.radians(angles[i][0]))*np.sin(np.radians(angles[i][1])), np.cos(np.radians(angles[i][0]))*np.cos(np.radians(angles[i][1])), np.sin(np.radians(angles[i][0]))]

for i in range (len(anglesp)):
    Ap[2*i]=[np.cos(np.radians(anglesp[i][2]))*np.cos(np.radians(anglesp[i][1]))-np.sin(np.radians(anglesp[i][2]))*np.cos(np.radians(anglesp[i][0]))*np.sin(np.radians(anglesp[i][1])), np.cos(np.radians(anglesp[i][2]))*np.sin(np.radians(anglesp[i][1]))+np.sin(np.radians(anglesp[i][2]))*np.cos(np.radians(anglesp[i][1]))*np.cos(np.radians(anglesp[i][0])), np.sin(np.radians(anglesp[i][2]))*np.sin(np.radians(anglesp[i][0]))]
    Ap[2*i+1]=[-1*np.sin(np.radians(anglesp[i][2]))*np.cos(np.radians(anglesp[i][1]))-np.cos(np.radians(anglesp[i][2]))*np.cos(np.radians(anglesp[i][0]))*np.sin(np.radians(anglesp[i][1])), -1*np.sin(np.radians(anglesp[i][2]))*np.sin(np.radians(anglesp[i][1]))+np.cos(np.radians(anglesp[i][2]))*np.cos(np.radians(anglesp[i][1]))*np.cos(np.radians(anglesp[i][0])), np.cos(np.radians(anglesp[i][2]))*np.sin(np.radians(anglesp[i][0]))]

B=np.dot(np.linalg.inv(np.dot(A.transpose(),A)),A.transpose())
Bp=np.dot(np.linalg.inv(np.dot(Ap.transpose(),Ap)),Ap.transpose())

B12=centre_of_mass(A1,A2)
B34=centre_of_mass(A4,A3)
B56=centre_of_mass(A5,A6)
B78=centre_of_mass(A8,A7)
B910=centre_of_mass(A9,A10)
ksi=np.array([B12[0],B12[1],B34[0],B34[1],B56[0],B56[1],B78[0],B78[1],B910[0],B910[1],A12[0],A11[0],A12[1],A11[1],A12[2],A11[2],A12[3],A11[3]])
ksip=np.array([item for pair in zip(A1, A2) for item in pair]+[item for pair in zip(A4, A3) for item in pair]+[item for pair in zip(A5, A6) for item in pair]+[item for pair in zip(A8, A7) for item in pair]+[item for pair in zip(A9, A10) for item in pair]+[item for pair in zip(A12, A11) for item in pair][:-2])

delta=np.dot(B,ksi.transpose())
deltap=np.dot(Bp,ksip.transpose())

print("La bille doit être optimalement décalée de ",-1*delta[0]," mm vers le bras, de ",delta[1], " mm vers la droite du bras et de ",delta[2]," mm vers le haut.")
              
plt.figure()

plt.subplot(231)
ray=np.sqrt(np.mean((x0-XFG0)**2+(y0-YFG0)**2))
circle_f=plt.Circle((x0,y0),radius=ray,color='red',fill=False)
plt.scatter(XFG0,YFG0,color='red',label='Champ')
plt.scatter(XBG0,YBG0,color='orange',label='Bille')
plt.scatter(x0,y0,color='brown',label='Centroïde Champ')
plt.xlim(637,641)
plt.ylim(637,641)
fig1 = plt.gcf()
ax1 = fig1.gca()
ax1.add_artist(circle_f)
plt.text(x0-1,y0-1, "Excentrage : "+str(round(0.224*ray,3))+ "mm")
plt.text(np.mean(XBG0),np.mean(YBG0),"Incertitude : "+str(round(0.224*np.sqrt(np.std(XBG0)**2+np.std(YBG0)**2),3))+ "mm")
plt.title('G=0, C variable')
plt.legend(loc='upper left')
plt.gca().set_aspect('equal', adjustable='box')

plt.subplot(234)
ray=np.sqrt(np.mean((x90-XFG90)**2+(y90-YFG90)**2))
circle_f=plt.Circle((x90,y90),radius=ray,color='red',fill=False)
plt.scatter(XFG90,YFG90,color='red')
plt.scatter(XBG90,YBG90,color='orange')
plt.scatter(x90,y90,color='brown')
plt.xlim(637,641)
plt.ylim(637,641)
fig4 = plt.gcf()
ax4 = fig4.gca()
ax4.add_artist(circle_f)
plt.text(x90-1,y90-1, "Excentrage : "+str(round(0.224*ray,3))+ "mm")
plt.text(np.mean(XBG90),np.mean(YBG90),"Incertitude : "+str(round(0.224*np.sqrt(np.std(XBG90)**2+np.std(YBG90)**2),3))+ "mm")
plt.title('G=90, C variable')
plt.gca().set_aspect('equal', adjustable='box')

plt.subplot(232)
ray=np.sqrt(np.mean((x180-XFG180)**2+(y180-YFG180)**2))
circle_f=plt.Circle((x180,y180),radius=ray,color='red',fill=False)
plt.scatter(XFG180,YFG180,color='red')
plt.scatter(XBG180,YBG180,color='orange')
plt.scatter(x180,y180,color='brown')
plt.xlim(637,641)
plt.ylim(637,641)
fig2 = plt.gcf()
ax2 = fig2.gca()
ax2.add_artist(circle_f)
plt.text(x180-1,y180-1, "Excentrage : "+str(round(0.224*ray,3))+ "mm")
plt.text(np.mean(XBG180),np.mean(YBG180),"Incertitude : "+str(round(0.224*np.sqrt(np.std(XBG180)**2+np.std(YBG180)**2),3))+ "mm")
plt.title('G=180, C variable')
plt.gca().set_aspect('equal', adjustable='box')

plt.subplot(235)
plt.scatter(XFG180Ts,YFG180Ts,color='red')
plt.scatter(XBG180Ts,YBG180Ts,color='orange')
plt.xlim(637,641)
plt.ylim(637,641)
plt.title('G=180, T variable')
plt.gca().set_aspect('equal', adjustable='box')

plt.subplot(233)
ray=np.sqrt(np.mean((x270-XFG270)**2+(y270-YFG270)**2))
circle_f=plt.Circle((x270,y270),radius=ray,color='red',fill=False)
plt.scatter(XFG270,YFG270,color='red')
plt.scatter(XBG270,YBG270,color='orange')
plt.scatter(x270,y270,color='brown')
plt.xlim(637,641)
plt.ylim(637,641)
fig3 = plt.gcf()
ax3 = fig3.gca()
ax3.add_artist(circle_f)
plt.text(x270-1,y270-1, "Excentrage : "+str(round(0.224*ray,3))+ "mm")
plt.text(np.mean(XBG270),np.mean(YBG270),"Incertitude : "+str(round(0.224*np.sqrt(np.std(XBG270)**2+np.std(YBG270)**2),3))+ "mm")
plt.title('G=270, C variable')
plt.gca().set_aspect('equal', adjustable='box')

plt.subplot(236)
ray=np.sqrt(np.mean((x180e-XFG180E)**2+(y180e-YFG180E)**2))
circle_f=plt.Circle((x180e,y180e),radius=ray,color='red',fill=False)
plt.scatter(XFG180E,YFG180E,color='red',label='Champ')
plt.scatter(XBG180E,YBG180E,color='orange',label='Bille')
plt.scatter(x180e,y180e,color='brown', label='Centroïde Champ')
plt.xlim(637,641)
plt.ylim(637,641)
fig6 = plt.gcf()
ax6 = fig6.gca()
ax6.add_artist(circle_f)
plt.text(x180e-1,y180e-1, "Excentrage : "+str(round(0.224*ray,3))+ "mm")
plt.text(np.mean(XBG180E),np.mean(YBG180E),"Incertitude : "+str(round(0.224*np.sqrt(np.std(XBG180E)**2+np.std(YBG180E)**2),3))+ "mm")
plt.title('G=180E, C variable')
plt.legend(loc='lower right')
plt.gca().set_aspect('equal', adjustable='box')

plt.show()

plt.figure()

xs=[0,90,179.9,180,270]
ys=[s*(np.mean(XBG0)-x0),s*(np.mean(XBG90)-x90),s*(np.mean(XBG180E)-x180e),s*(np.mean(XBG180)-x180),s*(np.mean(XBG270)-x270)]
plt.scatter(xs,ys)