import pydicom as dcm
import SimpleITK as sitk
import numpy as np
from math import *
import os
import matplotlib.pyplot as plt
from scipy import ndimage
import cv2
from skimage import feature, filters, morphology
import time

#
# Mettre en place le remplissage pondéré de l'EH
# 

###----------------NOTE POUR L'UTILISATION-----------------###

#Pour utiliser le fichier, bien se placer dans le bon répertoire de travail. Par exemple :

# import os
# os.chdir('E:\\Stage 2020\\Scripts') : répertoire où se trouve le fichier centers.py
# ddd='E:\\Stage 2020\\Acquisitions 10-03' : dossier des acquisitions
# import centers

# Puis appeler les fonctions simplement de la façon suivante :

# centers.analyse_1(ddd)
# centers.analyse_2(ddd)
# centers.analyse_3a(ddd)
# centers.analyse_3b(ddd)   #1min 30 si save_fig=False

# Pour travailler avec les simulations, soit rentrer les paramètres selon ses envies, soit utiliser les paramètres stables: 

# centers.simul(1)
# centers.simul_WL(1,1)  ~30s)

###--------------------------------------------------------###


def extract_ROI(data_dcm_directory,Image_file,Image,Is_a_file,diam):
    
    if Is_a_file:
        #Lecture des données pixel du fichier DICOM
        Lec = sitk.ReadImage(os.path.join(data_dcm_directory, Image_file))
                
        #Conversion en array
        Image = sitk.GetArrayFromImage(Lec[:,:,0])
        space_y=Lec.GetSpacing()[0]
        space_x=Lec.GetSpacing()[1]
    else:
        space_x=0.336
        space_y=0.336
        
    ind_max_y=np.argmax(Image,axis=0)
    ind_max_x=np.argmax(Image,axis=1)
    
    margin_y = (3/2)*2*diam/space_y
    margin_x = (3/2)*2*diam/space_x
    
    ROI=Image[ind_max_y-margin_y:ind_max_y+margin_y,ind_max_x-margin_x:ind_max_x+margin_x]
    
    return (ROI,ind_max_x,ind_max_y,margin_x,margin_y)


def Hough_center(data_dcm_directory, Image_file, Image, Is_a_file, sig, pui, thr, resol_r=0.5, resol_x=1, resol_y=1):
    
    if Is_a_file:
        #Lecture des données pixel du fichier DICOM
        Lec = sitk.ReadImage(os.path.join(data_dcm_directory, Image_file))
        
        #Conversion en array
        Image = sitk.GetArrayFromImage(Lec[:,:,0])
        space_y=Lec.GetSpacing()[0]
        space_x=Lec.GetSpacing()[1]
        
    #Conversion de l'image de non signé 16bits vers le type float
    Image = Image.astype(float)
    #Lissage de l'image par un filtre gaussien
    IG8 = filters.gaussian(Image,sigma=sig)
    # IG8 = filters.median(Image,selem=morphology.square(5))
    #Augmentation du contraste de l'image filtrée
    IG8 = IG8/(0.7*np.max(IG8))
    IG8 = IG8**pui
    
    
    # Obtention de l'image de gradient par filtre de Scharr
    Scharr = filters.scharr(IG8)
   
    #Définition du seuil pour l'image de gradient
    B = thr*np.max(Scharr)
    Thresh = Scharr*(Scharr>B)
    
    #Par soucis d'economie en temps de calcul et gestion de mémoire, on va rogner notre espace de travail. En effet, si on travaille avec une forme circulaire, travailler sur l'ensemble de l'image n'a pas d'intérêt : l'information est localisée.
    #On cherche donc les indices minimums et maximums selon les deux axes de l'image, et on travaille donc sur un espace de Hough réduit, et un remplissage depuis une image tronquée.
    
    y_min = np.min(np.nonzero(Thresh)[0])
    y_max = np.max(np.nonzero(Thresh)[0])
    x_min = np.min(np.nonzero(Thresh)[1])
    x_max = np.max(np.nonzero(Thresh)[1])
    
    Dx = x_max-x_min
    Dy = y_max-y_min
    
    Dxt = int(round(0.3*(Dx)))   # Marges en x et y au delà de xmin et xmax
    Dyt = int(round(0.3*(Dy)))

    #Figure de contrôle : Image de gradient seuillé
    plt.figure()
    plt.title("Image de seuil de " + str(Image_file))
    plt.imshow(Thresh[y_min-Dyt:y_max+Dyt+1,x_min-Dxt:x_max+Dxt+1])
    
    # #Figure de contrôle : Image de base après pré-traitement
    # plt.figure()
    # plt.title("Image filtrée de " + str(Image_file))
    # plt.imshow(IG8[y_min-Dyt:y_max+Dyt+1,x_min-Dxt:x_max+Dxt+1])


    D = max(Dx,Dy)
    Dt = max(Dxt,Dyt)
    
    #Valeurs de contrôles affichées.
    # print(y_min,y_max,x_min,x_max,Dx,Dy,Dxt,Dyt,Nrt)
    
    #Création de l'espace de Hough
    
    Hough_space = np.zeros((int(Dt/resol_r),int((2*D)/resol_y)+1,int((2*D)/resol_x)+1),dtype=float)                                                               # Création de l'espace de Hough
    
    beta_x = (1/2)*(1-resol_x)
    beta_y = (1/2)*(1-resol_y)
    
    # n=0
    for yi in range (Dy+1):               # Les xi et yi (image) ne sont non nuls
        for xi in range (Dx+1):           # que sur les True de l'image de seuil.
            T=Thresh[yi+y_min,xi+x_min]
            if (T!=0):
                
                for rh in range(0,int((Dt/resol_r))):
                    
                    xh = np.arange(int((2*D)/resol_y)+1)   #On crée un tableau vide de même
                    yh = np.arange(int((2*D)/resol_x)+1)   #taille que les dim x et y de
                                                              #l'espace de Hough

                    #Masque circulaire de centre (xi, yi), de rayon (D/2-Dt) + rh
                    
                    # Calcul du coefficient epsilon de l'équation (x-xc)**2 + ( y-yc)**2 - (r-r0)**2 < epsilon**2
                    
                    alpha = np.sqrt((resol_x/2)**2+(resol_y/2)**2)
                    
                    epsilon2 = abs(2*alpha*(rh*resol_r + (D/2) - Dt) + alpha**2)
                    
                    # Matrice conditionnelle sur epsilon : on calcule tous les membres de gauche de l'eq au dessus
                    
                    cond = abs((xh[np.newaxis,:]*resol_x-beta_x-(xi+D/2))**2+(yh[:,np.newaxis]*resol_y-beta_y-(yi+D/2))**2-(rh*resol_r + D/2 - Dt)**2)
                    
                    # Test : Soit on est inférieur à epsilon2 et on garde la valeur, sachant que plus on est proche de epsilon moins on a de poids
                    # Soit on est supérieur et on prend la valeur 0. np clip assure que toutes les valeurs sont comprises entre 0 et 1.
                    # La condition <1 est inutile ici, mais un argument doit être passé à np clip. On met 1.1 pour ne pas exclure les cas
                    # parfaits ou cond = 0.
                    
                    w = np.clip((epsilon2-cond)/epsilon2,0,1.1) # Matrice des poids
                    
                    # if n==30:
                    #     plt.figure()
                    #     plt.subplot(121)
                    #     plt.imshow(w)
                    #     plt.subplot(122)
                    #     plt.imshow(cond)               # Boucle à activer pour démonstration. Penser à décommenter au dessus n=0 et en
                    #     print(epsilon2)                # dessous n+=1.
                    
                    Hough_space[int(rh)]+=w*T
                    # n+=1
 

    
    #Augmentation du contraste de l'espace de Hough
    Hough_space = Hough_space/np.percentile(Hough_space, 95)                         
    Hough_space_10 = Hough_space**5
    
    #Valeur max de l'espace de Hough : le cdm n'a de sens qu'a proximité
    Max_Hough = np.max(Hough_space_10)
    #On retire donc toutes les valeurs trop éloignées du max
    Hough_trunc = np.where(Hough_space_10>0.01*Max_Hough,Hough_space_10,0)
    #Centre de masse de l'espace de Hough
    cdm_10 = ndimage.measurements.center_of_mass(Hough_trunc)
    #Figures de contrôle
    # print(cdm_10[0],D,Dt)
    # plt.figure()
    # plt.title("Espace de Hough proche du max de "+str(Image_file))
    # plt.subplot(121)
    # plt.imshow(Hough_space_10[int(np.ceil(cdm_10[0])),:,:])
    # plt.subplot(122)
    # plt.imshow(Hough_space_10[int(np.floor(cdm_10[0])),:,:])
    
    
    r_est_10=cdm_10[0]*resol_r + D/2 - Dt
    y_est_10 = (cdm_10[1])*resol_y + y_min-(D/2) - beta_y
    x_est_10 = (cdm_10[2])*resol_x + x_min-(D/2) - beta_x
    
    if not Is_a_file:
        space_y = 0.336
        space_x = 0.336
    
    print("x_est = ",round(x_est_10,3),", y_est = ",round(y_est_10,3),", r_est = ",round(r_est_10,3))   # Coordonnées en pixel dans l'image
    
    return (x_est_10, y_est_10, r_est_10, space_x, space_y)



def Hough_edge_center(data_dcm_directory, Image_file, Image, Is_a_file, sig, pui, thr, resol_r=0.25, resol_x=0.25, resol_y=0.25):
    
    if Is_a_file:
        #Lecture des données pixel du fichier DICOM
        Lec = sitk.ReadImage(os.path.join(data_dcm_directory, Image_file))
        
        #Conversion en array
        Image = sitk.GetArrayFromImage(Lec[:,:,0])
        space_y = Lec.GetSpacing()[0]
        space_x = Lec.GetSpacing()[1]
        
    #Conversion de l'image de non signé 16bits vers le type float
    Image = Image.astype(float)
    #Lissage de l'image par un filtre gaussien
    IG8 = filters.gaussian(Image,sigma=sig)
    
    #Augmentation du contraste de l'image filtrée
    IG8 = IG8/(0.7*np.max(IG8))
    IG8 = IG8**pui
    
    
    # Obtention de l'image de gradient par filtre de Scharr
    Scharr = filters.scharr(IG8)
    # On recupère également l'information sur les gradients horizontaux et verticaux. scharr_v donne les contours verticaux, donc les gradients
    # horizontaux. Vice-versa pour scharr_h
    Gx = filters.scharr_v(IG8)
    Gy = filters.scharr_h(IG8)
    
    # Figures de contrôle
    # plt.figure()
    # plt.subplot(131)
    # plt.imshow(Scharr)
    # plt.subplot(132)
    # plt.imshow(Schx)
    # plt.subplot(133)
    # plt.imshow(Schy)
   
    #Définition du seuil pour l'image de gradient
    B = thr*np.max(Scharr)
    Thresh = Scharr*(Scharr>B)
    
    #Par soucis d'economie en temps de calcul et gestion de mémoire, on va rogner notre espace de travail. En effet, si on travaille avec une forme circulaire, travailler sur l'ensemble de l'image n'a pas d'intérêt : l'information est localisée.
    #On cherche donc les indices minimums et maximums selon les deux axes de l'image, et on travaille donc sur un espace de Hough réduit, et un remplissage depuis une image tronquée.
    
    y_min = np.min(np.nonzero(Thresh)[0])
    y_max = np.max(np.nonzero(Thresh)[0])
    x_min = np.min(np.nonzero(Thresh)[1])
    x_max = np.max(np.nonzero(Thresh)[1])
    
    Dx = x_max-x_min
    Dy = y_max-y_min
    
    Dxt = int(round(0.3*(Dx)))   # Marges en x et y au delà de xmin et xmax
    Dyt = int(round(0.3*(Dy)))

    # #Figure de contrôle : Image de gradient seuillé
    # plt.figure()
    # plt.title("Image de seuil de " + str(Image_file))
    # plt.imshow(Thresh[y_min-Dyt:y_max+Dyt+1,x_min-Dxt:x_max+Dxt+1])
    
    # #Figure de contrôle : Image de base après pré-traitement
    # plt.figure()
    # plt.title("Image filtrée de " + str(Image_file))
    # plt.imshow(IG8[y_min-Dyt:y_max+Dyt+1,x_min-Dxt:x_max+Dxt+1])


    D = max(Dx,Dy)
    Dt = max(Dxt,Dyt)
    
    #Valeurs de contrôles affichées.
    # print(y_min,y_max,x_min,x_max,Dx,Dy,Dxt,Dyt,Nrt)
    
    #Création de l'espace de Hough
    
    Hough_space = np.zeros((int(Dt/resol_r),int((2*D)/resol_y)+1,int((2*D)/resol_x)+1),dtype=float)                                                               # Création de l'espace de Hough
    
    beta_x = (1/2)*(1-resol_x)
    beta_y = (1/2)*(1-resol_y)

    for yi in range (Dy+1):               # Les xi et yi (image) ne sont non nuls
        for xi in range (Dx+1):           # que sur les True de l'image de seuil.
            T = Thresh[yi+y_min,xi+x_min]
            if (T!=0):
                
                for rh in range(0,int((Dt/resol_r))):
                    
                    #theta est l'angle entre l'orientation du gradient et l'horizontale. Géométriquement, on a tan(theta)=Gy/Gx
                    theta = np.arctan2(Gy[yi+y_min,xi+x_min],Gx[yi+y_min,xi+x_min])
                    # theta est en plus signé, grâce aux méthodes scharr_v et scharr_h. On a donc seulement un seul point à considérer,
                    # à l'intérieur de l'image de grad, et sur l'axe défini par l'orientation du gradient
                    xh = ((xi+D/2))/resol_x + ((rh*resol_r + (D/2) - Dt)*np.cos(theta))/resol_x + beta_x
                    yh = ((yi+D/2))/resol_y + ((rh*resol_r + (D/2) - Dt)*np.sin(theta))/resol_y + beta_y
                    #Enfin, puisque cette méthode est sensible au bruit, on n'ajoute pas un seul point mais un kernel 3*3 centré sur ce point.
                    Hough_space[int(rh),int(round(yh))-1:int(round(yh))+2,int(round(xh))-1:int(round(xh))+2]+=0.001

    
    # Augmentation du contraste de l'espace de Hough
    # Hough_space=Hough_space/np.percentile(Hough_space, 97)                         
    Hough_space_10=Hough_space**5
    
    #Valeur max de l'espace de Hough : le cdm n'a de sens qu'a proximité
    Max_Hough = np.max(Hough_space_10)
    #On retire donc toutes les valeurs trop éloignées du max
    Hough_trunc = np.where(Hough_space_10>0.01*Max_Hough,Hough_space_10,0)
    #Centre de masse de l'espace de Hough
    cdm_10 = ndimage.measurements.center_of_mass(Hough_trunc)
    
    #Figures de contrôle
    # plt.figure()
    # plt.title("Espace de Hough proche du max de "+str(Image_file))
    # plt.subplot(121)
    # plt.imshow(Hough_space_10[int(np.ceil(cdm_10[0])),:,:])
    # plt.subplot(122)
    # plt.imshow(Hough_space_10[int(np.floor(cdm_10[0])),:,:])
    
    
    r_est_10=cdm_10[0]*resol_r + D/2 - 2*Dt
    y_est_10 = (cdm_10[1])*resol_y + y_min-(D/2) - beta_y
    x_est_10 = (cdm_10[2])*resol_x + x_min-(D/2) - beta_x
    
    if not Is_a_file:
        space_y = 0.336
        space_x = 0.336
    
    print("x_est = ",round(x_est_10,3),", y_est = ",round(y_est_10,3),", r_est = ",round(r_est_10,3))   # Coordonnées en pixel dans l'image
    
    return (x_est_10, y_est_10, space_x, space_y)

def cdm_sans_seuil(data_dcm_directory, Image_file, sig):
   
   
    #Lecture des données pixel du fichier DICOM
    Lec = sitk.ReadImage(os.path.join(data_dcm_directory, Image_file))
        
    #Conversion en array
    Image = sitk.GetArrayFromImage(Lec[:,:,0])    

    #Conversion de l'image de non signé 16bits vers le type float
    Image = Image.astype(float)
    #Lissage de l'image par un filtre gaussien
    IG8 = filters.gaussian(Image,sigma=sig)
    #Augmentation du contraste de l'image filtrée
    IG8 = IG8/(0.7*np.max(IG8))

    cdm=ndimage.measurements.center_of_mass(IG8)
    
    space_y = Lec.GetSpacing()[0]
    space_x = Lec.GetSpacing()[1]
    y_est = cdm[0]
    x_est = cdm[1]
    
    print ("x_est = ",round(x_est,3)," y_est = ",round(y_est,3))   # Coordonnées en pixel dans l'image
    return(x_est, y_est, space_x, space_y)

def cdm_avec_seuil(data_dcm_directory, Image_file, sig, pui, thr):
    

    #Lecture des données pixel du fichier DICOM
    Lec = sitk.ReadImage(os.path.join(data_dcm_directory, Image_file))
        
    #Conversion en array
    Image = sitk.GetArrayFromImage(Lec[:,:,0])    

    #Conversion de l'image de non signé 16bits vers le type float
    Image = Image.astype(float)
    #Lissage de l'image par un filtre gaussien
    IG8 = filters.gaussian(Image,sigma=sig)
    #Augmentation du contraste de l'image filtrée
    IG8=IG8/(0.7*np.max(IG8))
    IG8 = IG8
    
    
    # Obtention de l'image de gradient par filtre de Scharr
    Scharr = (filters.scharr(IG8))**pui
   
    #Définition du seuil pour l'image de gradient
    B = thr*np.max(Scharr)
    Thresh = Scharr>B
    
    cdm = ndimage.measurements.center_of_mass(Thresh)
    
    space_y = Lec.GetSpacing()[0]
    space_x = Lec.GetSpacing()[1]
    y_est = cdm[0]
    x_est = cdm[1]
    
    print ("x_est = ",round(x_est,3)," y_est = ",round(y_est,3))  # Coordonnées en pixel dans l'image
    return(x_est, y_est, space_x, space_y)

def morpho_split(data_dcm_directory, Image_file, Image, Is_a_file, sig, diam):
    
    if Is_a_file:
        
        #Lecture des données pixel du fichier DICOM
        Lec = sitk.ReadImage(os.path.join(data_dcm_directory, Image_file))
        
        #Conversion en array
        Image = sitk.GetArrayFromImage(Lec[:,:,0])    
    
    #Conversion de l'image de non signé 16bits vers le type float
    Image = Image.astype(float)
    #Lissage de l'image par un filtre gaussien
    IG8 = Image
    #IG8 = filters.gaussian(Image,sigma=sig)
    #Augmentation du contraste de l'image filtrée
    IG8 = IG8/(0.7*np.max(IG8))
    
    indices = np.unravel_index(np.argmax(IG8,axis=None),IG8.shape)
    max_y=int(np.median(indices[0],overwrite_input=True))
    max_x=int(np.median(indices[1],overwrite_input=True))
    margin_y = int(2*(3/2)*diam/0.336)  #Marge de sécurité de deux fois le diamètre
    margin_x = int(2*(3/2)*diam/0.336)  # en +/-x et en +/-y.
    
    print(max_y,margin_y,max_x,margin_x)
        
    #Définition du Kernel pour le filtrage morphologique
    Kerb = morphology.disk(40)
    Field=IG8.copy()
    Filled_Field=IG8.copy()
    #Kerf = morphology.disk(60)
    #Field = morphology.white_tophat(IG8,Kerf)    # A priori inutile de faire du top-hat pour avoir le champ. Mathématiquement, pour
    #                                             # enlever la bille, il suffit d'ajouter le bottom-hat au top hat, qui est ici l'image de base.
    #                                             # Cela peut être faux, attention!
    Ball = morphology.black_tophat(IG8[max_y-margin_y:max_y+margin_y,max_x-margin_x:max_x+margin_x],Kerb)
    
    B = np.zeros((1280,1280))
    B[max_y-margin_y:max_y+margin_y,max_x-margin_x:max_x+margin_x] = Ball
    Filled_Field[max_y-margin_y:max_y+margin_y,max_x-margin_x:max_x+margin_x]+=Ball
    
    Ball=B
    #Figures de contrôle
    
    # plt.figure()
    # plt.subplot(131)
    # plt.imshow(Filled_Field)
    # plt.subplot(132)
    # plt.imshow(Field)
    # plt.subplot(133)
    # plt.imshow(Ball)
    
    return(Filled_Field, Ball)


##### Méthodes calculatoires #####

def centre_of_mass(list_x,list_y):
    
    # Application de la définition du centroïde d'un polygone.
    
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

##########################################################
###############        SIMULATIONS        ################
##########################################################
    
def simul(pnb=3, I1=600, I0=0, rc=20, sig=0.5, pui=2, thr=0.7, xc=None, yc=None):
    
    # pnb : taille de la pénombre, en pixels
    # I1 : Intensité du champ
    # I0 : Intensité du fond
    # rc : rayon du cercle
    Sim = np.zeros((200,200)) # Création de l'image simulée
    if ((xc==None and yc==None)):
        xc = 100.0 + np.random.normal(scale=2)  # Génération des points centraux du champ, non entiers
        yc = 100.0 + np.random.normal(scale=2)
    print(" -------- Coordonnées simulées des centres --------")
    print("xc = ",round(xc,3),", yc = ",round(yc,3))
    print ("--------- Génération de l'image simulée ----------")
    for xi in range (200):
        for yi in range (200):
            dij = np.sqrt((xi-xc)**2+(yi-yc)**2)
            Sim[yi,xi] = int(max((I0-I1)/(1+np.exp(-(dij-rc)/pnb)) + I1 + np.random.normal(scale=0.3),0))
    print ("---------- Estimation du centre champ -----------")
    x_est,y_est,r_est,sx,sy = Hough_center(None, None, Sim, False, sig, pui, thr, resol_r=0.5, resol_x=0.5, resol_y=0.5)
    # plt.figure()
    # plt.imshow(Sim)
    # circle=plt.Circle((x_est,y_est),radius=r_est,color='red',fill=False)
    # plt.scatter(x_est,y_est,color='red')
    # fig = plt.gcf()
    # ax = fig.gca()
    # ax.add_artist(circle)
    # plt.xlim(60,140)
    # plt.ylim(60,140)
    # plt.show()
    return(xc,yc,x_est,y_est)

def simul_shift(pnb=3, I1=600, I0=0, rc=20, sig=0.5, pui=2, thr=0.7,steps=5):
    
    x0,y0,x0_est,y0_est=simul(pnb,I1,I0,rc,sig,pui,thr)
    Ex=np.zeros((steps,steps))
    Ey=np.zeros((steps,steps))
    for i in range(steps):
        for j in range(steps):
            S1,S2,S3,S4=simul(pnb,I1,I0,rc,sig,pui,thr,xc=(x0-0.5+i/steps),yc=(y0-0.5+j/steps))
            Ex[j,i]=S3-(x0-0.5+i/steps)
            Ey[i,j]=S4-(y0-0.5+j/steps)
    return (Ex,Ey)

def simul_WL(pnbc=1.2, pnbb=0.5, I1=600, I0=0, rext=18, I2=170, rint=7.5, sig=0.5, pui=2, thr=0.7, xc=None, yc=None, xb=None, yb=None):
    
    # pnb : taille de la pénombre, pour le champ (c) et la bille (b)
    # I1 : Intensité du champ
    # I0 : Intensité du fond
    # rext : rayon du champ
    # rint : rayon de la bille
    # I2 : intensité max perdue du champ sous la bille
    Sim = np.zeros((1280,1280)) # Création de l'image simulée
    if ((xc==None or yc==None or xb==None or yb==None)):
        xc = 640 + np.random.normal(scale=2)  # Génération des points centraux du champ, non entiers
        yc = 640+ np.random.normal(scale=2)
        xb = 640 + np.random.normal(scale=1)  # Génération des points centraux du champ, non entiers
        yb = 640 + np.random.normal(scale=1)
    print(" -------- Coordonnées simulées des centres --------")
    print("xc = ",xc,", yc = ",yc,", xb = ",xb,", yb = ",yb)
    print ("--------- Génération de l'image simulée ----------")
    for xi in range (1280):
        for yi in range (1280):
            dijc = np.sqrt((xi-xc)**2+(yi-yc)**2)
            dijb = np.sqrt((xi-xb)**2+(yi-yb)**2)
            Sim[yi,xi] = int(max((I0-I1)/(1+np.exp(-(dijc-rext)/pnbc)) + I1 + np.random.normal(scale=0.3) - ((I0-I2)/(1+np.exp(-(dijb-rint)/pnbb)) + I2 + np.random.normal(scale=0.3)),0))
     #Conversion de l'image de non signé 16bits vers le type float
    Sim = Sim.astype(float)
    #Lissage de l'image par un filtre gaussien
    # plt.figure()
    # plt.imshow(Sim)
    print ("---------- Méthodes morphologiques ... -----------")
    F,B = morpho_split(None, None, Sim, False, 0.1,10)
    print ("---------- Estimation du centre champ -----------")
    x_est_c, y_est_c, r_est_c, sx, sy = Hough_center(None, None, F, False, sig, pui, thr)
    print ("---------- Estimation du centre bille -----------")
    x_est_b, y_est_b, r_est_b, sx, sy = Hough_center(None, None, B, False, sig, pui, 0.5*thr)
    
    plt.figure()
    plt.imshow(Sim)
    circle_f=plt.Circle((x_est_c,y_est_c),radius=r_est_c,color='red',fill=False)
    circle_b=plt.Circle((x_est_b,y_est_b),radius=r_est_b,color='blue',fill=False)
    plt.scatter(x_est_c,y_est_c,color='red',label='Centre champ')
    plt.scatter(x_est_b,y_est_b,color='blue',label='Centre bille')
    fig = plt.gcf()
    ax = fig.gca()
    ax.add_artist(circle_f)
    ax.add_artist(circle_b)
    plt.xlim(600,680)
    plt.ylim(600,680)
    plt.title("Centres bille et champ estimés à partir de l'image simulée.")
    plt.legend()
    plt.show() 

    return(xc,yc,xb,yb,x_est_c,y_est_c,x_est_b,y_est_b)

def simul_WL_shift(pnbc=2, pnbb=1, I1=600, I0=0, rext=15, I2=200, rint=7.5, sig=0.5, pui=2, thr=0.7,steps=3):
    
    xc0,yc0,xb0,yb0,x0_est_c,y0_est_c,x0_est_b,y0_est_b=simul_WL(pnbc,pnbb,I1,I0,rext,I2,rint,sig,pui,thr)
    Exc=np.zeros((steps,steps))
    Eyc=np.zeros((steps,steps))
    Exb=np.zeros((steps,steps))
    Eyb=np.zeros((steps,steps))
    for i in range(steps):
        for j in range(steps):
            S1,S2,S3,S4,S5,S6,S7,S8=simul_WL(pnbc,pnbb,I1,I0,rext,I2,rint,sig,pui,thr,xc=(xc0-0.5+i/steps),yc=(yc0-0.5+j/steps),xb=(xb0-0.5+i/steps),yb=(yb0-0.5+j/steps))
            Exc[j,i]=S5-(xc0-0.5+i/steps)
            Eyc[i,j]=S6-(yc0-0.5+j/steps)
            Exb[i,j]=S7-(xb0-0.5+i/steps)
            Eyb[i,j]=S8-(yb0-0.5+j/steps)
    return (Exc,Eyc,Exb,Eyb)
##########################################################
################        ANALYSES        ##################
##########################################################
    
def analyse_1(data_dcm_directory,sig=0.5, pui=2, thr=0.7):
    file_list = ['G0T0C0cc10-000.dcm','G0T0C0cc10-001.dcm','G0T0C0cc10-002.dcm','G0T0C0cc10-003.dcm','G0T0C0cc10-004.dcm']
    #,'G0T0C0cc10-001.dcm','G0T0C0cc10-002.dcm','G0T0C0cc10-003.dcm','G0T0C0cc10-004.dcm'
    print (" 1 - Etude de la reproductibilité du positionnement des CC")
    print ("--------Etude via cdm sans seuil--------")
    x_cdss = []
    y_cdss = []
    x_cdm = []
    y_cdm = []
    x_h10 = []
    y_h10 = []
    for i in range(len(file_list)):
        print (" -------------------------------------- ")
        print (" Acquisition ",i )
        x,y,sx,sy = cdm_sans_seuil(data_dcm_directory, file_list[i], sig)
        x_cdss.append(x*sx)
        y_cdss.append(y*sy)
    print(" -------------------------------------- ")
    print ("l'écart type de la mesure est de ",round((2/3)*np.std(x_cdss),3), " mm en x et de ",round((2/3)*np.std(y_cdss),3), " mm en y." )
    print ("l'écart max de la mesure est de ",round((2/3)*np.max(x_cdss)-(2/3)*np.min(x_cdss),3), " mm en x et de ",round((2/3)*np.max(y_cdss)-(2/3)*np.min(y_cdss),3), " mm en y." )
    print ("--------Etude via cdm avec seuil--------")
    for i in range(len(file_list)):
        print (" -------------------------------------- ")
        print (" Acquisition ",i )
        x,y,sx,sy = cdm_avec_seuil(data_dcm_directory, file_list[i], sig, pui, thr)
        x_cdm.append(x*sx)
        y_cdm.append(y*sy)
    print(" -------------------------------------- ")
    print ("l'écart type de la mesure est de ",round((2/3)*np.std(x_cdm),3), " mm en x et de ",round((2/3)*np.std(y_cdm),3), " mm en y." )
    print ("l'écart max de la mesure est de ",round((2/3)*np.max(x_cdm)-(2/3)*np.min(x_cdm),3), " mm en x et de ",round((2/3)*np.max(y_cdm)-(2/3)*np.min(y_cdm),3), " mm en y." )
    print ("--------Etude via TH sur filtre--------")
    for i in range(len(file_list)):
        print (" -------------------------------------- ")
        print (" Acquisition ",i )
        x10,y10,r10,sx,sy = Hough_center(data_dcm_directory, file_list[i], None, True, sig, pui, thr)
        x_h10.append(x10*sx)
        y_h10.append(y10*sy)
        Lec = sitk.ReadImage(os.path.join(data_dcm_directory, file_list[i]))
        Image = sitk.GetArrayFromImage(Lec[:,:,0])
        plt.figure()
        circle=plt.Circle((x10,y10),radius=r10,color='red',fill=False)
        plt.imshow(Image)
        plt.xlim(600,680)
        plt.ylim(600,680)
        fig = plt.gcf()
        ax = fig.gca()
        ax.add_artist(circle)
        plt.scatter(x10,y10,color='red',label='Centre champ')
        
        plt.title("Centres bille et champ estimés à partir de l'acquisition " + str(i))
        plt.legend()
    print(" -------------------------------------- ")
    print ("L'écart type de la mesure est de ",round((2/3)*np.std(x_h10),3), " mm en x et de ",round((2/3)*np.std(y_h10),3), " mm en y." )
    print ("L'écart max de la mesure est de ",round((2/3)*np.max(x_h10)-(2/3)*np.min(x_h10),3), " mm en x et de ",round((2/3)*np.max(y_h10)-(2/3)*np.min(y_h10),3), " mm en y." )
    
def analyse_2(data_dcm_directory,sig=0.5, pui=2, thr=0.7):
    file_list = ['G0T0C0cc-000.dcm','G0T0C0cc-001.dcm','G0T0C0cc-002.dcm','G0T0C0cc-003.dcm','G0T0C0cc-004.dcm','G0T0C0cc-005.dcm']
    print (" 1 - Etude de la reproductibilité du positionnement des CC")
    print ("--------Etude via cdm sans seuil--------")
    x_cdss = []
    y_cdss = []
    x_cdm = []
    y_cdm = []
    x_h10 = []
    y_h10 = []
    for i in range(len(file_list)):
        print (" -------------------------------------- ")
        print (" Acquisition ",i )
        x,y,sx,sy = cdm_sans_seuil(data_dcm_directory, file_list[i], sig)
        x_cdss.append(x*sx)
        y_cdss.append(y*sy)
    print(" -------------------------------------- ")
    print ("l'écart type de la mesure est de ",round((2/3)*np.std(x_cdss),3), " mm en x et de ",round((2/3)*np.std(y_cdss),3), " mm en y." )
    print ("l'écart max de la mesure est de ",round((2/3)*np.max(x_cdss)-(2/3)*np.min(x_cdss),3), " mm en x et de ",round((2/3)*np.max(y_cdss)-(2/3)*np.min(y_cdss),3), " mm en y." )
    print ("---------Etude via cdm avec seuil---------")
    for i in range(len(file_list)):
        print (" -------------------------------------- ")
        print (" Acquisition ",i )
        x,y,sx,sy = cdm_avec_seuil(data_dcm_directory, file_list[i], sig, pui, thr)
        x_cdm.append(x*sx)
        y_cdm.append(y*sy)
    print(" -------------------------------------- ")
    print ("l'écart type de la mesure est de ",round((2/3)*np.std(x_cdm),3), " mm en x et de ",round((2/3)*np.std(y_cdm),3), " mm en y." )
    print ("l'écart max de la mesure est de ",round((2/3)*np.max(x_cdm)-(2/3)*np.min(x_cdm),3), " mm en x et de ",round((2/3)*np.max(y_cdm)-(2/3)*np.min(y_cdm),3), " mm en y." )
    print ("---------Etude via TH---------")
    for i in range(len(file_list)):
        print (" -------------------------------------- ")
        print (" Acquisition ",i )
        x10,y10,r10,sx,sy = Hough_center(data_dcm_directory, file_list[i], None, True, sig, pui, thr)
        x_h10.append(x10*sx)
        y_h10.append(y10*sy)
        Lec = sitk.ReadImage(os.path.join(data_dcm_directory, file_list[i]))
        Image = sitk.GetArrayFromImage(Lec[:,:,0])
        plt.figure()
        circle=plt.Circle((x10,y10),radius=r10,color='red',fill=False)
        plt.imshow(Image)
        plt.xlim(600,680)
        plt.ylim(600,680)
        fig = plt.gcf()
        ax = fig.gca()
        ax.add_artist(circle)
        plt.scatter(x10,y10,color='red',label='Centre champ')
        
        plt.title("Centres bille et champ estimés à partir de l'acquisition " + str(i))
        plt.legend()
    print(" -------------------------------------- ")
    print ("L'écart type de la mesure est de ",round((2/3)*np.std(x_h10),3), " mm en x et de ",round((2/3)*np.std(y_h10),3), " mm en y." )
    print ("L'écart max de la mesure est de ",round((2/3)*np.max(x_h10)-(2/3)*np.min(x_h10),3), " mm en x et de ",round((2/3)*np.max(y_h10)-(2/3)*np.min(y_h10),3), " mm en y." )
    
def analyse_3a(data_dcm_directory,sig=1.1, pui=2, thr=0.5):
    file_list_1 = ['G0T0Ccc10-000.dcm','G0T0Ccc10-001.dcm']
    file_list_2 = ['G0T0Ccc10-002.dcm','G0T0Ccc10-003.dcm','G0T0Ccc10-004.dcm','G0T0Ccc10-005.dcm','G0T0Ccc10-006.dcm']
    print("----- Détermination du centre de rotation collimateur avec deux points -----")
    x1 = []
    y1 = []
    xs1 = []
    ys1 = []
    for i in range(len(file_list_1)):
        print (" -------------------------------------- ")
        print (" Acquisition ",i )
        x,y,r,sx,sy = Hough_center(data_dcm_directory, file_list_1[i], None, True, sig, pui, thr)
        xs1.append(x*sx)
        ys1.append(y*sy)
        x1.append(x)
        y1.append(y)
    print(" -------------------------------------- ")
    print ("Le centre se situe en position (pixels) ",round(np.mean(x1),3), " en x et de ",round(np.mean(y1),3), " en y." )
    print("----- Détermination du centre de rotation collimateur avec cinq points -----")
    x51 = []
    y51 = []
    xs51 = []
    ys51 = []
    for i in range(len(file_list_2)):
        print (" -------------------------------------- ")
        print (" Acquisition ",i )
        x,y,r,sx,sy = Hough_center(data_dcm_directory, file_list_2[i], None, True, sig, pui, thr)
        xs51.append(x*sx)
        ys51.append(y*sy)
        x51.append(x)
        y51.append(y)
    (xc,yc)=centre_of_mass(x51,y51)
    (px,py)=(xc*sx,yc*sy)
    print(" -------------------------------------- ")
    print ("Le centre se situe en position (pixels) ",round(xc,3), " en x et de ",round(yc,3), " en y." )
    print(" Soit un décalage entre les deux positions de : ",round((2/3)*abs(px-np.mean(xs1)),3), " mm en x et de ",round((2/3)*abs(py-np.mean(ys1)),3), " mm en y.")
    print ("L'erreur de centrage du cône sur l'axe de rotation du collimateur est estimée à : ",round((2/3)*np.sqrt((np.mean(xs1)-xs1[0])**2 + (np.mean(ys1)-ys1[0])**2),3)," avec la méthode à 2 points et à : ",round((2/3)*np.sqrt(np.mean((xs1-px)**2+(ys1-py)**2)),3), " mm selon la méthode à 5 points")

def analyse_3b(data_dcm_directory,txtfilename,save_fig=True,sig=0.5,pui=2,thr=0.7):
    
    file_list = []
    for i in range(29):
        file_list.append('GTCcc10-'+str(i).zfill(3)+'.dcm')
        # Liste des acquisitions WL
        
    res_file = open(txtfilename,"w")
    Data = []
    Header = ['File Name', 'Gantry Angle', 'Collimator Angle', 'Table Angle', 'Field Coordinates', 'Ball Coordinates']
    Data.append(Header)
    t = time.time()  # Contrôle de la durée de l'analyse
    t1 = t
    i = 1 
    for file in file_list:
        
        print('Acq '+str(i)+' sur '+str(len(file_list)))
        # Filtrage morphologique de l'image
        F,B = morpho_split(data_dcm_directory, file, None, True, 0.1, 10)
        # Détermination des centres bille et champ
        xf, yf, rf, sx, sy = Hough_center(data_dcm_directory, None, F, False, sig, pui, thr)
        xb, yb, rb, sx, sy = Hough_center(data_dcm_directory, None, B, False, sig, pui, 0.7*thr)
        # Récupération des tags DICOM
        tags = dcm.dcmread(os.path.join(data_dcm_directory,file))
        # Récupération du contenu des tags avec le mot clé correspondant.
        Data.append([file,str(tags['GantryAngle'].value),str(tags['BeamLimitingDeviceAngle'].value),str(tags['PatientSupportAngle'].value),str((xf,yf)),str((xb,yb))])
        t2=time.time()
        print("traitement de l'acquistion en : ",round(t2-t1,1)," s.")
        print("Temps total : ", round(t2-t,1), " s.")
        t1 = t2
        i+=1
        #Lecture des données pixel du fichier DICOM
        if save_fig:
            
            Lec = sitk.ReadImage(os.path.join(data_dcm_directory, file))
            Image = sitk.GetArrayFromImage(Lec[:,:,0])
            plt.figure()
            plt.imshow(Image)
            circle_f=plt.Circle((xf,yf),radius=rf,color='red',fill=False)
            circle_b=plt.Circle((xb,yb),radius=rb,color='blue',fill=False)
            plt.xlim(600,680)
            plt.ylim(600,680)
            plt.scatter(xf,yf,color='red',label='Centre champ')
            plt.scatter(xb,yb,color='blue',label='Centre bille')
            fig = plt.gcf()
            ax = fig.gca()
            ax.add_artist(circle_f)
            ax.add_artist(circle_b)
            plt.title("G = " + str(round(tags['GantryAngle'].value)) + ", T = " + str(round(tags['PatientSupportAngle'].value)) + " C = " + str(round(tags['BeamLimitingDeviceAngle'].value)) )
            plt.legend()
            plt.savefig("Acq "+str(i-1)+".png")
            plt.close()
            
    for line in Data :
        res_file.write("\t".join(line) + "\n")
    
    res_file.close()
    return None

##### ANALYSE DU FICHIER TXT CREE ####
##analyse_txt.py##
    
