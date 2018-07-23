# -*- coding: utf-8 -*-
"""
Created on Fri Jun 29 12:15:59 2018

@author: kgratsea
"""

# -*- coding: utf-8 -*-
"""
Quantum walk with GDT. 3 free parameters for each coin operator

function to be minimzed norm

@author: kgratsea
"""

import numpy as np
import cmath
import math
from scipy import optimize
global final
import random

metrhths=0

def tensor(vectorA,vectorB) :
    m = np.size(vectorA,0)
    n = np.size(vectorB,0)
    tens=np.zeros((m,n))
    for i in range(m) :
        for j in range(n) :
            tens[i][j] = vectorA[i]*vectorB[j]
    return (tens);

def func(z) :  
    global metrhths 
    metrhths += 1
    n=2 #number of steps
    k=2*n+1 #number of sites at the final state
    
    initial = np.zeros((2*k,1),dtype=complex)
    #localised on one site
    initial[2*n][0]= 1.
    initial[2*n+1][0]= 1.5
    initial/= np.linalg.norm(initial)
    
    #Initial = initial
    #print (Initial)
   
    #definition of matrixS  
    #the one I use
    qplate = np.zeros((2*k,2*k),dtype=complex)
    #(m,up)--> (m+1,up) (m,down)--> (m-1,down)
    
    i=1
    while (i+2) < 2*k :
        qplate[i][i+1] = 1.0
        i += 2
       
    j=1
    while (j+2) < 2*k :
        qplate[j+1][j] = 1.0
        j += 2
        
    matrixS = qplate    
        
    listSt = []
    listc = []
    listC = []

    listSt.append (initial)
    
    #Define coin operators with gdt
    
    l = 0 # for corresponding the correct coin parameters at each step n
    for j in range (0,n,+1) : 
        #print ("n",j)
        c=np.zeros((2,2),dtype=complex)
        theta=z[0+l]
        ksi=z[1+l]
        zeta=z[2+l]
        
        c[0][0]= (math.cos(ksi*math.pi) + math.sin(ksi*math.pi)*1j)*math.cos(theta*math.pi/2)
        c[0][1]= (math.cos(zeta*math.pi) + math.sin(zeta*math.pi)*1j)*math.sin(theta*math.pi/2) 
        c[1][0]= (math.cos(zeta*math.pi) - math.sin(zeta*math.pi)*1j)*math.sin(theta*math.pi/2)         
        c[1][1]= - (math.cos(ksi*math.pi) - math.sin(ksi*math.pi)*1j)*math.cos(theta*math.pi/2)  
        
        listc.append(c)
        matrixC = np.zeros((2*k,2*k),dtype=complex)
        #print (c)
        
        for i in range (0,2*k,2):
            matrixC[0+i][0+i] = c[0][0]
            matrixC[1+i][1+i] = c[1][1]
            matrixC[0+i][1+i] = c[0][1]          
            matrixC[1+i][0+i] = c[1][0]   
         
        listC.append (matrixC)    
        
        #print (initial)
        m1 = np.dot(matrixC,initial)
        m2 = np.dot(matrixS,m1)   #next state
        #print ("steps")
        #print (m1)
        #print (m2)
        listSt.append (m2)
        initial = m2/np.linalg.norm(m2)
        l += 3 # moving to the next coin parameters
        
    Phi=initial    
    Phi_reshaped =np.zeros((2,k),dtype=complex)
    for i in range(0,2,1):
        q=0
        for j in range(0,k,1): 
            Phi_reshaped[i][j] = Phi[i+q][0]
            q +=2
    #Phi_internal = np.delete(Phi,[0,1,2*k*2*k-2,2*k*2*k-1],None) #delete outer sites
    #print ("Phi_reshaped",Phi_reshaped)
    
    psiA, l, psiB = np.linalg.svd(Phi_reshaped,full_matrices=1) #decomposition of initial matrix
    #print ("l",l)
    
    NORM=0.0
    p=1.0 # p has to be larger or equal than 1 for the algorithm to work
    
    m = np.size(Phi_reshaped,0)  #number of rows of initial matrix
    n = np.size(Phi_reshaped,1)  #number of columns of initial matrix
    sum=np.zeros((m,n))
    
    for i in range(2) :
         tens= tensor(psiA[i],psiB[i])
         sum += math.pow(l[i],p-1)*tens
         NORM = NORM + math.pow(l[i],p) 

    NORM = math.pow(NORM,1./p)

    #print (NORM)
    
    with open('QW_NORM_qplate_n10_T1_steps2_niter-success5_only_NORM.txt', 'a+') as f:
        print (metrhths,",",NORM,file=f)
    f.close()
    
    if (-NORM+math.sqrt(2)<0.0000001) :
        f = open("QW_NORM_qplate_n10_T1_steps2_niter-success5.txt","a+")
        f.write("initial")
        f.close()
        with open('QW_NORM_qplate_n10_T1_steps2_niter-success5.txt', 'a+') as f:
            print (initial,file=f)
        f.close()
        
        f = open("QW_NORM_qplate_n10_T1_steps2_niter-success5.txt","a+")
        f.write("l,NORM")
        f.close()
        with open('QW_NORM_qplate_n10_T1_steps2_niter-success5.txt', 'a+') as f:
            print (l,NORM,file=f)
        f.close()
    
        f = open("QW_NORM_qplate_n10_T1_steps2_niter-success5.txt","a+")
        f.write("z")
        f.close()
        with open('QW_NORM_qplate_n10_T1_steps2_niter-success5.txt', 'a+') as f:
            print (z,file=f)
        f.close()
        
        
        f = open("QW_NORM_qplate_n10_T1_steps2_niter-success5.txt","a+")
        f.write("listc")
        f.close()
        with open('QW_NORM_qplate_n10_T1_steps2_niter-success5.txt', 'a+') as f:
            print (listc,file=f)
        f.close()
        
        
        f = open("QW_NORM_qplate_n10_T1_steps2_niter-success5.txt","a+")
        f.write("Phi")
        f.close()
        with open('QW_NORM_qplate_n10_T1_steps2_niter-success5.txt', 'a+') as f:
            print (Phi,file=f)
        f.close()
        
        f = open("QW_NORM_qplate_n10_T1_steps2_niter-success5.txt","a+")
        f.write("Phi_reshaped")
        f.close()
        with open('QW_NORM_qplate_n10_T1_steps2_niter-success5.txt', 'a+') as f:
            print (Phi_reshaped,file=f)
        f.close()    
    
    return (-NORM+ math.sqrt(2))

    #max.ent.over sites i ={2,3}
    #Phi_target= np.array([[ 0.0066874],       [ 0.       ],       [ 0.6148   ],       [ 0.3492   ],       [ 0.3493   ],       [-0.6148   ],       [ 0.       ],       [ 0.0067874]])
    
    '''
    #max.ent.over sites i ={2,5}    
    #d = [0.48107031, 0.63099229]  #solution of orthogonality condition
    #Phi_target = np.array([ [0.4811], [0.52619533], [0.27592482], [ 0.26093188],  [0.2510644],[0 ], [0], [0.2290105] , [ 0.48391584], [ 0.3214155], [ 0.35487592],[0.6310]])
    Phi_target= np.array([  [0.4811], [0] ,[0.52619533],  [0.2290105] ,  [0.27592482],  [ 0.48391584], [ 0.26093188],  [ 0.3214155], [0.2510644], [ 0.35487592] ,[0 ], [0.6310] ])
    Phi_target /= np.linalg.norm(Phi_target)
    Phi_target_internal = np.delete(Phi_target,[0,1,2*k-2,2*k-1],None)
    print ("Phi_target",Phi_target) '''

   
   
#initial_coin_parameters=[1/2,0,math.pi,1/4,math.pi/2.,math.pi/2.,1/8,2*math.pi,0]  #n=3
#initial_coin_parameters=[1/2,0,math.pi,1/4,math.pi/2.,math.pi/2.,1/8,2*math.pi,0,1/2,3*math.pi,math.pi/4,1/3,2*math.pi/3,math.pi/3] #n=5
#Initial_coin_par=initial_coin_parameters  
    
f=2
    
my_randoms=[]
for i in range (3*f):
    my_randoms.append(random.randrange(1,10,1))

print (my_randoms)    
    
initial_coin_parameters=my_randoms

        
minimizer_kwargs = {"method": "BFGS"}


ret = optimize.basinhopping(func,initial_coin_parameters, minimizer_kwargs=minimizer_kwargs,niter=10, T=1.0, disp = True,niter_success=5 )
  
l=0
listc=[]
for j in range (0,f,+1) : 
        print ("j",j)
        c=np.zeros((2,2),dtype=complex)
        theta=ret.x[0+l]
        ksi=ret.x[1+l]
        zeta=ret.x[2+l]
        
        c[0][0]= (math.cos(ksi*math.pi) + math.sin(ksi*math.pi)*1j)*math.cos(theta*math.pi/2)
        c[0][1]= (math.cos(zeta*math.pi) + math.sin(zeta*math.pi)*1j)*math.sin(theta*math.pi/2) 
        c[1][0]= (math.cos(zeta*math.pi) - math.sin(zeta*math.pi)*1j)*math.sin(theta*math.pi/2)         
        c[1][1]= - (math.cos(ksi*math.pi) - math.sin(ksi*math.pi)*1j)*math.cos(theta*math.pi/2)  
        
        listc.append(c)
        l+=3