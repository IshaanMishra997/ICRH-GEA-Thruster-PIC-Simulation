#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov  9 10:20:40 2020
author: Ishaan Mishra
Hybrid Cyclotron Plasma & Gridded Ion Thruster
"""
import pylab as pl
import numpy as np 
import random as rd
import math

eps0 = 8.854187e-12
e = 1.602176e-19
amu=1.67377e-27
me = 9.109e-31
Kb = 1.3806e-23 # Boltzmann constant
eV = 11604.52 # eV in Kelvin; QE/K
J_to_eV=6.242e+18
mol=6.02214076e23

#Info of Helium, Argon, Xenon
He_m= 4.002602 * amu
Ar_m= 39.948 * amu
Xe_m= 131.293 * amu

He_BP=5193.1*He_m*160*J_to_eV #eV
Ar_BP=520.33*Ar_m*160*J_to_eV #eV
Xe_BP=158.32*Xe_m*160*J_to_eV #eV

He_ie = 2372.3e3*J_to_eV/mol
Ar_ie = 1520.6e3*J_to_eV/mol
Xe_ie = 1170.4e3*J_to_eV/mol

mi=Xe_m
starttemp=Xe_BP
temp=Xe_ie + Xe_BP

'''
Algorithm: 
    1) Setup - where particles enter, walls, exit
    2) Magnetic fields - magnetic mirror, nozzle
    
    LOOP:
        3) Helicon Plasma Generator-refer to papers for ionisation %
        4) Ion, Electron density solver          (Fully kinetic/Hybrid PIC???)
        5) Potential solver
        6) Electric field solver
        7) Move particles
        8) Record exit velocity
'''

#---------------------------------DEFINITION---------------------------------#

dr = dz = 1e-2
dt = 5e-9
mdot = 10e-3 #mass flow rate in kg 
plasma_den = n0 = 10e18#m^-3 

spwt = 1e9
r1=1e-2#radis of particles entering chamber
#nump=plasma_den*a1*dz*nz1/2/spwt
nump = int(mdot/mi*dt/spwt)

phi0=0

#particle definition
class particle:
    def __init__(self,x,y,vz,vr):
        self.x=x
        self.y=y
        self.vz=vz
        self.vr=vr
        self.KE=1/2*(vz**2 + vr**2)
    def is_at_wall(self):
        global nz,nr, nodes, i 
        if i<3: 
            return False
        elif self.x<0: 
            return True 
        elif self.y<0: 
            return True
        else: 
            node_valx=int(self.x/dz)
            node_valy=int(self.y/dr)
            if node_valy>nodes[node_valx][node_valy]:
                return True
            else: 
                return False
            
        
class species:      
    def __init__(self, mass, charge, spwt, starttemp):
        self.mass=mass
        self.charge=charge*e
        self.spwt=spwt
        self.starttemp=starttemp
        self.arr=[]
        self.nump=len(self.arr)


#-----------------------------------SETUP-----------------------------------#


#defining thruster nodes, region 
l1 = 3e-2
l2 = l1 + 45e-2
l3 = l2 + 18e-2
l4 = l3+ 2e-2
nz = np.zeros(int(l4/dz))
nr= np.zeros(int(l4/dz))
i=0
c=0
numcells=int(l4/dz)-1 

while c<int(l4/dz):
    nz[c]=i
    if i<l1: 
        nr[c]=(5/3*nz[c])
        k=nr[c]
    elif i<l2: 
        nr[c]=(k)
        #k=x[c]-y[c]
    elif i<l3: 
       nr[c]=((nz[c]-l2)/np.sqrt(6)+k)
       k2=nr[c]
    elif i<=l4: 
        nr[c]=(k2)
    i+=dz
    c+=1

    
pl.figure()
pl.plot (nz,nr, 'r')
pl.xlabel('z(m)')
pl.ylabel('r(m)')
   

nodes=[]
le=0
c=0
for i in nr: 
    nodes.append([j*dr for j in range(int(i/dr)+1)])
    #nodes.append(np.zeros(int(i/dr)))
    le+=int(i/dr)
    c+=1

for i in range(len(nodes)): 
    k=[]
    for j in range(len(nodes[i])):
        k.append(i*dz)
    pl.scatter(k,nodes[i], s=1)
print(le)
#print(nodes)

#setting up node volume 
node_volume=[]
for i in range(0,len(nodes)):
    node_volume.append([])
    for j in range(0,len(nodes[i])):
        j_min = j-0.5
        j_max = j+0.5
        if j_min<0: 
            j_min=0
        if j_max>len(nodes[i])-1: 
            j_max=len(nodes[i])-1
        
        if (i==0 or i==(numcells)-1):
            a=0.5
        else:
            a=1.0
        #note, this is r*dr for non-boundary nodes
        node_volume[i].append(a*dz*((j_max*dr)**2-(j_min*dr)**2))
        #print(a*dz*((j_max*dr)**2-(j_min*dr)**2))
node_volume[0][0]=1/2*node_volume[1][0]


#initialising ions, electrons 
ions=species(mi,1,spwt,starttemp)
electrons=species(me, -1, spwt, starttemp)

empty_node_set=[]

for i in range(len(nodes)): 
    tp=[]
    for i in range(len(nodes[i])):
        tp.append(0)
    empty_node_set.append(tp)
#print(empty_node_set)

EF=phi=rho=empty_node_set


#-----------------------------REQUIRED_FUNCTIONS-----------------------------#


def generate(species, nump):
    for i in range(nump):
        v_th= (1/2*Kb*species.starttemp/species.mass)**(1/2)
        r=rd.random()
        #print(r)
        if r<0.5:
            r=1-r
        r1= rd.random()
        r2= rd.random()
        r3= rd.random()
        v= 2**(1/2) * v_th * (r1+r2+r3- 1.5)
        vz=np.sqrt((v**2)*r)
        vr=np.sqrt((v**2)*(1-r))
        r=rd.random()
        y=r*len(nodes[0])*dr

        species.arr += [particle(0.025,y+0.001,vz,vr)]
        species.nump=len(species.arr)        

#def solve_density_electrons():
 #   global phi
#print(node_volume)
def solve_density(species):
    global node_volume
    #identifying 4 nodes
    #print('doin')
    for i in range(len(species.arr)):
        p=species.arr[i]
        if p.x>0:
            znode=int(p.x/dz)
        #print(p.y, dr)
        if p.y>0:
            rnode=int(p.y/dr)
        
        #print(p.y)
        #print('ay yo')
        n1=[znode,rnode,0] # znode, rnode, charge
        n2=[0,rnode,0]
        n3=[znode,0,0]
        n4=[0,0,0]    
        def rnodecheck():
            #global n3,n4,rnode
            if rnode==0 and p.y>0:
                n3[1]=n4[1]=rnode+1
            elif rnode==len(rho[len(rho)-1])-1 and p.y<rnode*dr:
                n3[1]=n4[1]=rnode-1
            elif rnode*dr>p.y and rnode>0 and rnode<len(nodes[znode]):
                n3[1]=n4[1]=rnode-1
            elif rnode*dr<p.y and rnode>0 and rnode<len(nodes[znode]):
                n3[1]=n4[1]=rnode+1
        
        if znode==0 and p.x>0: 
            n2[0]=n4[0]=znode+1
            rnodecheck()
        elif znode==len(nodes)-1 and p.x<znode*dz:
            n2[0]=n4[0]=znode-1
            rnodecheck()           
        elif znode*dz<p.x:
            n2[0]=n4[0]=znode+1
            rnodecheck()
        elif znode*dz>p.x:
            n2[0]=n4[0]=znode-1
            rnodecheck()
            
        #Assigning Density weightage to local nodes 
        
        lx=(1-np.abs(p.x/dz-znode))#*species.charge*species.spwt

        ly=(1-np.abs(p.y/dr-rnode*dr))#*species.charge*species.spwt
        #print(lx, ly)
        
        n1[2]=(lx*ly*species.charge*species.spwt)
        n2[2]=((1-lx)*ly*species.charge*species.spwt)
        n3[2]=(lx*(1-ly)*species.charge*species.spwt)
        n4[2]=((1-lx)*(1-ly)*species.charge*species.spwt)
        #print(n1)
        #print(node_volume[n1[0]][n1[1]])
        #print(node_volume[0][0])
        rho[n1[0]][n1[1]]+=n1[2]#/node_volume[n1[0]][n1[1]]
        rho[n2[0]][n2[1]]+=n2[2]#/node_volume[n2[0]][n2[1]]
        #print(n3)
        rho[n3[0]][n3[1]]+=n3[2]#/node_volume[n3[0]][n3[1]]
        rho[n4[0]][n4[1]]+=n4[2]#/node_volume[n4[0]][n4[1]]

    for i in range(len(node_volume)): 
        for j in range(len(node_volume[i])): 
            
            rho[i][j]/=node_volume[i][j]
            print(rho[i][j])
    
    #print('YEYEYEYE')
  

#simple Jacobian solver, does not do any convergence checking
def solvePotential():
    global phi
    max_it=100
    #make copy of dirichlet nodes
    P = phi 
    
    g = empty_node_set
    dz2 = dz*dz
    dr2 = dr*dr

    #rho_e = empty_node_set
    
    #set radia
    r = empty_node_set
    for i in range(len(r)):
        for j in range(len(r[i])):
            r[i][j] = j*dr
    #b=empty_node_set
    for i in range (max_it):
        for j in range(len(rho)):
            for k in range(len(rho[j])):
                rho[j][k]=rho[j][k]/eps0
                
        for j in range(1,len(g)):

            for k in range(1,len(g[j])):
                if k>=2:
                    p1=phi[j][k]
                else:
                    p1=0
                if k<len(g[j])-1:
                    p2=p1=phi[j][k]
                else: 
                    p2=0
                if j>=2:
                    p3=phi[j][k]
                else: 
                    p3=0
                if j<=len(g)-1:
                    p4=phi[j][k]
                else: 
                    p4=0
                    
                g[j][k]=(rho[j][k]+p1+p2)/dr2 + (p2-p1)/(2*dr*r[j][k]) + (p3+p4/dz2)/(2/dr2 + 2/dz2)

        #neumann boundaries
        g[0] = g[1]       #left
        g[-1] = g[-2]     #right
        g[:][-1] = g[:][-2] #top
        g[:][0] = g[:][1]
        
        #dirichlet nodes
        phi = np.where(True,P,g)
    return phi

#computes electric field  
efz=efr=empty_node_set                  
def computeEF():
    global phi, efz, efr
    last = len(efz[len(nodes)-1])-1
    
    
    for k in range(len(efz[0])):
        efz[0][k]=(phi[0][k] - phi[1][k])/dz
        print('0 efz is ',efz[0][k])
    for k in range(len(phi[last-1])):
        efz[last][k]=(phi[last-1][k] - phi[last][k])/dz
        
        
    for k in range(len(efz[0])):
        efr[k][0]=(phi[k][0] - phi[k][1])/dr
    for k in range(len(phi[last-1])):
        efr[k][len(efr[k])-1]=(phi[k][len(efr[k])-2] - phi[k][len(efr[k])-1])/dr
        
    
    
    
    for i in range(1,len(efz)-1): 
        for j in range(1,len(efz[i])-1):

            efz[i][j]= (efz[i][j-1]-efz[i][j+1])/(2*dz)
            
        for j in range(1,len(efr[i])-1):
            #print(j)
            if j>len(efr[i-1])-1:
                efr[i][j]=(phi[i][j] - phi[i+1][j])/dr
            else:
                efr[i][j]= (efr[i-1][j]-efr[i+1][j])/(2*dr)
            

    
def rebound(particle):
    particle.vr=-1*particle.vr
    particle.x+=particle.vz*dt
    particle.y+=particle.vr*dt




vbar=0
vbarz=0
vbarr=0
pcount=0    
def log_velocity(p):
    global vbar, vbarz, vbarr, pcount
    rat=species.mass/ions.mass
    vbar+=(pcount*vbar + (((p.vr**2 + p.vz**2**(1/2)))*rat))/(pcount+rat)
    vbarz+=(pcount*vbarz + p.vz*rat)/(pcount+rat)
    vbarr+=(pcount*vbarr + p.vr*rat)/(pcount+rat)
    pcount+=rat

def mover(species): 
    i=1
    while i < len(species.arr):
        p=species.arr[i]
        
        if p.x>len(efz)*dz:
            log_velocity(p)
            species.arr[i]=species.arr[len(species.arr)-1]
            del species.arr[len(species.arr)-1]
            i-=1
        
        if p.is_at_wall(): 
            rebound(p)
        else:
            if p.x>0 and p.y>0:
                znode=int(p.x/dz)
                rnode=int(p.y/dr)
            else: 
                break
            
            n1=[znode,rnode] # znode, rnode, charge
            n2=[0,rnode]
            n3=[znode,0]
            n4=[0,0]    
            def rnodecheck():
                if rnode==0 and p.y>0:
                    n3[1]=n4[1]=rnode+1
                elif rnode==len(nodes[len(rho[len(rho)-1])])-1 and p.y<rnode*dr:
                    n3[1]=n4[1]=rnode-1
                elif rnode*dr>p.y and rnode>0 and rnode<len(nodes[znode]):
                    n3[1]=n4[1]=rnode-1
                elif rnode*dr<p.y and rnode>0 and rnode<len(nodes[znode]):
                    n3[1]=n4[1]=rnode+1
            if znode==0 and p.x>0: 
                n2[0]=n4[0]=znode+1
                rnodecheck()
            elif znode==len(nodes)-1 and p.x<znode*dz:
                n2[0]=n4[0]=znode-1
                rnodecheck()           
            elif znode*dz<p.x:
                n2[0]=n4[0]=znode+1
                rnodecheck()
            elif znode*dz>p.x:
                n2[0]=n4[0]=znode-1
                rnodecheck()
            
            
            lx=(1-np.abs(p.x/dz-znode))
            ly=(1-np.abs(p.y/dr-rnode*dr))
            p.vz+=species.charge/species.mass*dt*(lx*efz[n1[0]][n1[1]] - (1-lx)*efz[n2[0]][n2[1]])
            p.vr+=species.charge/species.mass*dt*(ly*efr[n1[0]][n1[1]] - (1-ly)*efr[n3[0]][n3[1]]) 
            p.x+=p.vz*dt
            p.y+=p.vr*dt
            if math.isnan(p.x):
                pass
            species.arr[i]=p                                     

        i+=1


def ICRH(): 
    pass
#------------------------------------LOOP------------------------------------#
it=100
generate(ions,nump)
#generate(electrons,nump)
for i in range(it):
    print('iteration:', i+1)
    EF=phi=rho=empty_node_set
    solve_density(ions)
    #solve_density(electrons)
    solvePotential()
    computeEF()
    mover(ions)
    #mover(electrons)
    #print(len(ions.arr), len(electrons.arr))
    print(len(ions.arr))
    print(ions.arr[100].vz,ions.arr[100].vr)
    #print(electrons.arr[100].vz,electrons.arr[100].vr)
print("DONE HA")
print(vbar)
