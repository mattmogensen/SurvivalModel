#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 18 12:20:56 2020

@author: mogensen
"""

import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
import numpy as np
import numpy.random
from scipy.stats import norm
from scipy.stats import skewnorm
from scipy.stats import beta
from scipy.optimize import curve_fit
import pandas as pd
import seaborn as sns

def propagate(chemical,sp,dc):
    newchemical = chemical
    for i in range(size):
        for j in range(size):
            #If chemical present
            if chemical[i][j]>0:
                if i > 0:
                    #propagate up
                    if newchemical[i-1][j] == 0:
                        newchemical[i-1][j]=sp*chemical[i][j]
                if i < size-1:
                    #propagate down
                    if newchemical[i+1][j] == 0:
                        newchemical[i+1][j]=sp*chemical[i][j]
                if j > 0:
                    #propagate left
                    if newchemical[i][j-1] == 0:
                        newchemical[i][j-1]=sp*chemical[i][j]
                if j < size - 1:
                    #propagate right
                    if newchemical[i][j+1] == 0:
                        newchemical[i][j+1]=sp*chemical[i][j]
                        
    for i in range(size):
        for j in range(size):
            newchemical[i][j]=dc*newchemical[i][j]
            
    return newchemical
            
def moveup(model,vec):
    x,y = vec
    z = model[x][y]
    w = max(0,min(z,np.random.poisson(lam3), capacity - model[x-1][y]))
    model[x][y] = model[x][y] - w
    model[x-1][y] = model[x-1][y] + w - p*w
    return model, 0, p*w
    
def movedown(model,vec):
    x,y = vec
    z = model[x][y]
    w = max(0,min(z,np.random.poisson(lam3), capacity - model[x+1][y]))
    model[x][y] = model[x][y] - w
    model[x+1][y] = model[x+1][y] + w - p*w
    return model, 0, p*w
    
def moveleft(model,vec):
    x,y = vec
    z = model[x][y]
    w = max(0,min(z,np.random.poisson(lam2),capacity - model[x][y-1]))
    model[x][y] = model[x][y] - w
    model[x][y-1] = model[x][y-1] + w - p*w
    return model, 0, p*w
    
def moveright(model,vec):
    x,y = vec
    z = model[x][y]
    w = max(0,min(z,np.random.poisson(lam2),capacity - model[x][y+1]))
    model[x][y] = model[x][y] - w
    model[x][y+1] = model[x][y+1] + w - p*w
    return model, 0, p*w

def moveleftandup(model,vec):
    x,y = vec
    z = model[x][y]
    if np.random.uniform(1) < 0.5:
        w1 = max(0,min(z,np.random.poisson(lam2),capacity - model[x][y-1]))
        w2 = max(0,min(z-w1,np.random.poisson(lam3),capacity - model[x-1][y]))
        model[x][y] = model[x][y] - w1 - w2
        model[x][y-1] = model[x][y-1] + w1 - p*w1
        model[x-1][y] = model[x-1][y] + w2 - p*w2
        return model, 0, p*w1+p*w2
    else:
        w2 = max(0,min(z,np.random.poisson(lam3),capacity - model[x-1][y]))
        w1 = max(0,min(z-w2,np.random.poisson(lam2),capacity - model[x][y-1]))
        model[x][y] = model[x][y] - w1 - w2
        model[x][y-1] = model[x][y-1] + w1 - p*w1
        model[x-1][y] = model[x-1][y] + w2 - p*w2
        return model, 0, p*w1+p*w2
 
def moveleftanddown(model,vec):
    x,y = vec
    z = model[x][y]
    if np.random.uniform(1) < 0.5:
        w1 = max(0,min(z,np.random.poisson(lam2),capacity - model[x][y-1]))
        w2 = max(0,min(z-w1,np.random.poisson(lam3),capacity - model[x+1][y]))
        model[x][y] = model[x][y] - w1 - w2
        model[x][y-1] = model[x][y-1] + w1 - p*w1
        model[x+1][y] = model[x+1][y] + w2 - p*w2
        return model, 0, p*w1+p*w2
    else:
        w2 = max(0,min(z,np.random.poisson(lam3),capacity - model[x+1][y]))
        w1 = max(0,min(z-w2,np.random.poisson(lam2),capacity - model[x][y-1]))
        model[x][y] = model[x][y] - w1 - w2
        model[x][y-1] = model[x][y-1] + w1 - p*w1
        model[x+1][y] = model[x+1][y] + w2 - p*w2
        return model, 0, p*w1+p*w2
     
def evacuate(model,vec):
    x,y = vec
    z = model[x][y]
    
    #if no people present, do nothing
    if z==0:
        return model,0,0
    
    #If at evacuation point, evacuate
    if x == 0 and y == 0:
        w = max(0,min(z,np.random.poisson(lam1)))
        model[x][y] = model[x][y] - w
        return model, w, 0
    
    #If on top row, move left
    if (x == 0 and y > 0):
        return moveleft(model,vec)
    
    #If on left column, move up
    if x > 0 and y == 0:
        return moveup(model,vec)
    
    #If on bottom row, move left or up
    if x == size-1 and y <= chem[1]:
        return moveleftandup(model,vec)
    
    #If on bottom row, move left
    if x == size-1:
        return moveleft(model,vec)
    
    #Immediately right of chemical, move right
    if x == chem[0] and y == chem[0]+1:
        return moveright(model,vec)
    
    #if right of chemical, move up
    if x == chem[0] and y > chem[1]:
        return moveup(model,vec)
    
    #If immediately below chemical, move left
    if x == chem[0]+1 and y == chem[1]:
        return moveleft(model,vec)
    
    #In upper right quadrant: move left or up
    if x <= chem[0] and y > chem[1]:
        return moveleftandup(model,vec)
    
    #In lower right quadrant: move left or down
    if x > chem[0] and y > chem[1]:
        return moveleftanddown(model,vec)
    
    #If left of chemical:
    if y <= chem[1]:
        return moveleftandup(model,vec)
    
    return model,0,0

def clean(vec):
    vec2 = []
    for i in range(len(vec)):
        if vec[i] >= 0 and vec[i] < 0.001:
            vec2.append(0.001)
        if vec[i] <= 1 and vec[i] > 0.999:
            vec2.append(0.999)
        if vec[i] >= 0.001 and vec[i] <= 0.999:
            vec2.append(vec[i])
    return vec2

spreadfactor = [0.5, 0.7, 0.9,0.0] #lethality
decayfactor = [0.7, 0.9,0.0] #persistence
city = [100,100]

splabel=['Low','Medium','High','Zero (FP)']
dclabel=['low','high','zero (FP)']
pplabel=['rural','urban']

#evacuation parameters
lam1 = 800
lam2 = 300
lam3 = 100
capacity = 500
evac = 1
evaclbl ='Yes'
loops = 10

finaldata = {'Area':[],\
            'Lethality':[],\
            'Persistence':[],\
            'Delay':[],\
            'Evacuate':[],\
            'SafeA':[],\
            'SafeB':[],\
            'StampedeA':[],\
            'StampedeB':[],\
            'ChemicalA':[],\
            'ChemicalB':[]}

df = pd.DataFrame(data=finaldata)

print(df)

for z in range(2):
    q = 3
    w = 2
    #rural or urban
    #delay
    delaywarning = 0

#fixed parameters
    size = 10
    num = city[z]
    chem = [int(np.random.uniform(2,8)),int(np.random.uniform(2,8))]
    pop = size*size*num
    spl = splabel[q]
    dcl = dclabel[w]
    pp = pplabel[z]
    DatasetSafe = []
    DatasetStampede = []
    DatasetChemical = []
    
    print(str(q)+' '+str(w)+' '+str(z)+' '+str(delaywarning)+' '+str(evac))
    for n in range(loops):
        
        #display current progress
        Model = [[num for col in range(size)] for row in range(size)]
        Chemical = [[0 for col in range(size)] for row in range(size)]
        Chemical[chem[0]][chem[1]] = 1
        safe = 0
        stampededead = 0
        chemicaldead = 0
        p = max(0.001,np.random.normal(0.005,0.0001))
        sp = 0
        dc = 0
        spl = splabel[q]
    
        count = 0
        while count < delaywarning:
            Chemical = propagate(Chemical,sp,dc)
            count = count + 1
            
        count = 0
        while chemicaldead + stampededead + safe < pop-10**(-10) and count < 1000:
            
            count = count + 1
            #sns.heatmap(Model)
            #sns.heatmap(Chemical)
            #plt.show()
            
            #print(chemicaldead+stampededead+safe)
        
            #update people movement and people crowd deaths
            for i in range(size):
                for j in range(size):
                    if evac == 1:
                        Model, s, d = evacuate(Model,[i,j])
                        safe = safe + s
                        stampededead = stampededead + d
        
            #Update chemical deaths      
            Chemical = propagate(Chemical,sp,dc)
            for i in range(size):
                for j in range(size):
                    chemicaldead = chemicaldead + Chemical[i][j]*Model[i][j]
                    Model[i][j] = Model[i][j] - Chemical[i][j]*Model[i][j]
    
        DatasetSafe.append(safe/pop)
        DatasetStampede.append(stampededead/pop)

    #clean data
    DatasetSafe = clean(DatasetSafe)
    DatasetStampede = clean(DatasetStampede)

    # Fit a Beta distribution to the data:
    aSa, bSa , locSa, scaleSa = beta.fit(DatasetSafe, floc=0,fscale=1)
    aSt, bSt , locSt, scaleSt = beta.fit(DatasetStampede, floc=0,fscale=1)

    fig1,ax1 = plt.subplots()
    fig1.suptitle(str(spl) +' lethality, '+str(dcl)+' persistence, '+str(pp) + ' area (pop ' + str(pop) + ')')
    fig2,ax2 = plt.subplots()
    fig2.suptitle(str(spl) +' lethality, '+str(dcl)+' persistence, '+str(pp) + ' area (pop ' + str(pop) + ')')
        
    ax1.hist(DatasetSafe, bins=50, density=True, alpha=0.3, color='darkblue')
    xmin1, xmax1 = min(DatasetSafe), max(DatasetSafe)
    x1 = np.linspace(xmin1, xmax1, 100)
    p1 = beta.pdf(x1, aSa,bSa,locSa,scaleSa)
    ax1.plot(x1, p1, 'k', linewidth=2)
    title1 = " Delay warning: " + str(delaywarning) + ". Safely evacuated ~ Beta(%.2f,%.2f)" % (aSa, bSa)
    ax1.set_title(title1)
    fig1.savefig('Chemical_deaths_'+ 'None' +'-lethality_' + 'None' + '-persistence_' + 'None' + '-area_' + str(evaclbl) + '-evacuation.png')
    plt.close(fig1)   
        
    ax2.hist(DatasetStampede, bins=50, density=True, alpha=0.3, color='darkorange')
    xmin2, xmax2 = min(DatasetStampede), max(DatasetStampede)
    x2 = np.linspace(xmin2, xmax2, 100)
    p2 = beta.pdf(x2, aSt,bSt,locSt,scaleSt)
    ax2.plot(x2, p2, 'k', linewidth=2)
    title2 = " Delay warning: " + str(delaywarning) + ". Stampede deaths ~ Beta(%.2f,%.2f)" % (aSt, bSt)
    ax2.set_title(title2)
    fig2.savefig('Chemical_deaths_'+ 'None' +'-lethality_' + 'None' + '-persistence_' + 'None' + '-area_' + str(evaclbl) + '-evacuation.png')
    plt.close(fig2)
    
    df2 = pd.DataFrame(data = {'Area':[str(pp)],'Lethality':[str(spl)],'Persistence':[str(dcl)],'Delay':[str(delaywarning)],'Evacuate': [str(evaclbl)],'SafeA':[str(aSa)],\
          'SafeB':[str(bSa)],'StampedeA':[str(aSt)],'StampedeB':[str(bSt)],'ChemicalA':'NN/A','ChemicalB':'N/A','MeanSafe':[str(np.average(DatasetSafe))],\
          'MeanStampede':[str(np.average(DatasetStampede))],'MeanChemical':[str(0)]})
    
    df = df.append(df2)
        
    
for q in range(3):
    for w in range(2):
        for z in range(2):
            for delaywarning in range(7):
                evac = 1
                evaclbl = 'Yes'
                
                if delaywarning == 6:
                    evac = 0
                    evaclbl = 'No'
                    
                #display current progress
                print(str(q)+' '+str(w)+' '+str(z)+' '+str(delaywarning)+' '+str(evac))
                
                #fixed parameters
                size = 10
                num = city[z]
                chem = [int(np.random.uniform(2,8)),int(np.random.uniform(2,8))]
                pop = size*size*num
                spl = splabel[q]
                dcl = dclabel[w]
                pp = pplabel[z]
                DatasetSafe = []
                DatasetStampede = []
                DatasetChemical = []
                
                for n in range(loops):
                
                    Model = [[num for col in range(size)] for row in range(size)]
                    Chemical = [[0 for col in range(size)] for row in range(size)]
                    Chemical[chem[0]][chem[1]] = 1
                    safe = 0
                    stampededead = 0
                    chemicaldead = 0
                    p = max(np.random.normal(0.005,0.001),0.001)
                    #print(p)
                    sp = min(0.99,np.random.normal(spreadfactor[q],0.1))
                    #print(sp)
                    dc = min(0.99,np.random.normal(decayfactor[w],0.1))
                    #print(dc)
                    spl = splabel[q]
                    
                
                    count = 0
                    while count < delaywarning:
                        Chemical = propagate(Chemical,sp,dc)
                        count = count + 1
                        
                    count = 0
                    while (chemicaldead + stampededead + safe < pop-10**(-10)) and (count < 500):
                        
                        count = count + 1
                        sns.heatmap(Model)
                        #sns.heatmap(Chemical)
                        plt.show()
                    
                        #update people movement and people crowd deaths
                        for i in range(size):
                            for j in range(size):
                                if evac == 1:
                                    Model, s, d = evacuate(Model,[i,j])
                                    safe = safe + s
                                    stampededead = stampededead + d
                    
                        #Update chemical deaths      
                        Chemical = propagate(Chemical,sp,dc)
                        for i in range(size):
                            for j in range(size):
                                if Chemical[i][j]<0:
                                    print('error!!!')
                                if Chemical[i][j]>=1:
                                    print('erroororororor')
                                chemicaldead = chemicaldead + Chemical[i][j]*Model[i][j]
                                Model[i][j] = Model[i][j] - Chemical[i][j]*Model[i][j]
                        if chemicaldead> pop:
                            print('chemicaldead is ' + str(chemicaldead))
                            print('safe is ' + str(safe))
                            print('stampede is '+str(stampededead))
                
                    if evac == 1:
                        DatasetSafe.append(safe/pop)
                        DatasetStampede.append(stampededead/pop)
                    else:
                        DatasetSafe.append((pop-chemicaldead)/pop)
                        
                    DatasetChemical.append(chemicaldead/pop)               
                    
    

                if evac == 1:
                    DatasetSafe = clean(DatasetSafe)
                    DatasetStampede = clean(DatasetStampede)
                    DatasetChemical = clean(DatasetChemical)
                    aSa, bSa , locSa, scaleSa = beta.fit(DatasetSafe, floc=0,fscale=1)
                    aSt, bSt , locSt, scaleSt = beta.fit(DatasetStampede, floc=0,fscale=1)
                    aCh, bCh , locCh, scaleCh = beta.fit(DatasetChemical, floc=0,fscale=1)
                    
                    # evacuation, true positive
                    fig1,ax1 = plt.subplots()
                    fig1.suptitle(str(spl) +' lethality, '+str(dcl)+' persistence, '+str(pp) + ' area (pop ' + str(pop) + ')')
                    fig2,ax2 = plt.subplots()
                    fig2.suptitle(str(spl) +' lethality, '+str(dcl)+' persistence, '+str(pp) + ' area (pop ' + str(pop) + ')')
                    fig3,ax3 = plt.subplots()
                    fig3.suptitle(str(spl) +' lethality, '+str(dcl)+' persistence, '+str(pp) + ' area (pop ' + str(pop) + ')')
                    
                    ax1.hist(DatasetSafe, bins=50, density=True, alpha=0.3, color='darkblue')
                    xmin1, xmax1 = min(DatasetSafe), max(DatasetSafe)
                    x1 = np.linspace(xmin1, xmax1, 100)
                    p1 = beta.pdf(x1, aSa,bSa,locSa,scaleSa)
                    ax1.plot(x1, p1, 'k', linewidth=2)
                    title1 = " Delay warning: " + str(delaywarning) + ". Safely evacuated ~ Beta(%.2f,%.2f)" % (aSa, bSa)
                    ax1.set_title(title1)
                    fig1.savefig('Safely_Evacuated_'+ spl +'-lethality_' + dcl + '-persistence_' + pp + '-area_' + str(delaywarning) + '-delay.png')
                    plt.close(fig1)
                    
                    ax2.hist(DatasetStampede, bins=50, density=True, alpha=0.3, color='darkorange')
                    xmin2, xmax2 = min(DatasetStampede), max(DatasetStampede)
                    x2 = np.linspace(xmin2, xmax2, 100)
                    p2 = beta.pdf(x2, aSt,bSt,locSt,scaleSt)
                    ax2.plot(x2, p2, 'k', linewidth=2)
                    title2 = " Delay warning: " + str(delaywarning) + ". Stampede deaths ~ Beta(%.2f,%.2f)" % (aSt, bSt)
                    ax2.set_title(title2)
                    fig2.savefig('Stampede_deaths_'+ spl +'-lethality_' + dcl + '-persistence_' + pp + '-area_' + str(delaywarning) + '-delay.png')
                    plt.close(fig2)
                    
                    ax3.hist(DatasetChemical, bins=50, density=True, alpha=0.3, color='darkred')
                    xmin3, xmax3 = min(DatasetChemical), max(DatasetChemical)
                    x3 = np.linspace(xmin3, xmax3, 100)
                    p3 = beta.pdf(x3, aCh,bCh,locCh,scaleCh)
                    ax3.plot(x3, p3, 'k', linewidth=2)
                    title3 = " Delay warning: " + str(delaywarning) + ". Chemical deaths ~ Beta(%.2f,%.2f)" % (aCh, bCh)
                    ax3.set_title(title3)
                    fig3.savefig('Chemical_deaths_'+ spl +'-lethality_' + dcl + '-persistence_' + pp + '-area_' + str(delaywarning) + '-delay.png')
                    plt.close(fig3)
                    
                    df2 = pd.DataFrame(data = {'Area':[str(pp)],'Lethality':[str(spl)],'Persistence':[str(dcl)],'Delay':[str(delaywarning)],'Evacuate': [str(evaclbl)],'SafeA':[str(aSa)],\
                                                       'SafeB':[str(bSa)],'StampedeA':[str(aSt)],'StampedeB':[str(bSt)],'ChemicalA':[str(aCh)],'ChemicalB':[str(bCh)],'MeanSafe':[str(np.average(DatasetSafe))],\
                                                       'MeanStampede':[str(np.average(DatasetStampede))],'MeanChemical':[str(np.average(DatasetChemical))]})
                    df = df.append(df2)
                    
                else:
                    # no evacuation,false negative
                    DatasetSafe = clean(DatasetSafe)
                    DatasetChemical = clean(DatasetChemical)
                    aSa, bSa , locSa, scaleSa = beta.fit(DatasetSafe, floc=0,fscale=1)
                    aCh, bCh , locCh, scaleCh = beta.fit(DatasetChemical, floc=0,fscale=1)
                    
                    fig1,ax1 = plt.subplots()
                    fig1.suptitle(str(spl) +' lethality, '+str(dcl)+' persistence, '+str(pp) + ' area (pop ' + str(pop) + ')')
                    fig3,ax3 = plt.subplots()
                    fig3.suptitle(str(spl) +' lethality, '+str(dcl)+' persistence, '+str(pp) + ' area (pop ' + str(pop) + ')')
  
                    ax1.hist(DatasetSafe, bins=50, density=True, alpha=0.3, color='darkblue')
                    xmin1, xmax1 = min(DatasetSafe), max(DatasetSafe)
                    x1 = np.linspace(xmin1, xmax1, 100)
                    p1 = beta.pdf(x1, aSa,bSa,locSa,scaleSa)
                    ax1.plot(x1, p1, 'k', linewidth=2)
                    title1 = " Evacuate: " + str(evaclbl) + ". Safely evacuated ~ Beta(%.2f,%.2f)" % (aSa, bSa)
                    ax1.set_title(title1)
                    fig1.savefig('Safe_'+ spl +'-lethality_' + dcl + '-persistence_' + pp + '-area_' + str(evaclbl) + '-evacuatiion.png')
                    plt.close(fig1)
                    
                    ax3.hist(DatasetChemical, bins=50, density=True, alpha=0.3, color='darkred')
                    xmin3, xmax3 = min(DatasetChemical), max(DatasetChemical)
                    x3 = np.linspace(xmin3, xmax3, 100)
                    p3 = beta.pdf(x3, aCh,bCh,locCh,scaleCh)
                    ax3.plot(x3, p3, 'k', linewidth=2)
                    title3 = " Evacuate: " + str(evaclbl)  + ". Chemical deaths ~ Beta(%.2f,%.2f)" % (aCh, bCh)
                    ax3.set_title(title3)
                    fig3.savefig('Chemical_deaths_'+ spl +'-lethality_' + dcl + '-persistence_' + pp + '-area_' + str(evaclbl) + '-evacuation.png')
                    plt.close(fig3)

                    df2 = pd.DataFrame(data = {'Area':[str(pp)],'Lethality':[str(spl)],'Persistence':[str(dcl)],'Delay':[str(delaywarning)],'Evacuate': [str(evaclbl)],'SafeA':[str(aSa)],\
                                                       'SafeB':[str(bSa)],'StampedeA':'N/A','StampedeB':'N/A','ChemicalA':[str(aCh)],'ChemicalB':[str(bCh)],'MeanSafe':[str(np.average(DatasetSafe))],\
                                                       'MeanStampede':[str(0)],'MeanChemical':[str(np.average(DatasetChemical))]})
   
                    df = df.append(df2)

print(df)
df.to_excel('Evacuation_Data.xlsx')

