from __future__ import print_function
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import webbrowser
from numpy import trapz
from scipy.integrate import simps
from scipy.optimize import curve_fit
from mpl_toolkits.axes_grid1 import host_subplot
import mpl_toolkits.axisartist as AA


Speed_Crystal_T80=[]
Speed_Crystal_T90=[]
Speed_Crystal_T100=[]
Speed_Crystal_T110=[]
Speed_Crystal_T120=[]
Speed_Crystal_T130=[]
Speed_Crystal_T140=[]
Speed_Time_T80=[]
Speed_Time_T90=[]
Speed_Time_T100=[]
Speed_Time_T110=[]
Speed_Time_T120=[]
Speed_Time_T130=[]
Speed_Time_T140=[]
Rel_Cryst_T80=[]
Rel_Cryst_T90=[]
Rel_Cryst_T100=[]
Rel_Cryst_T110=[]
Rel_Cryst_T120=[]
Rel_Cryst_T130=[]
Rel_Cryst_T140=[]
Rel_TCryst_T80=[]
Rel_TCryst_T90=[]
Rel_TCryst_T100=[]
Rel_TCryst_T110=[]
Rel_TCryst_T120=[]
Rel_TCryst_T130=[]
Rel_TCryst_T140=[]



###Data thermodynamic simulation
Data = pd.read_csv('Temperature_LH02.csv.')
Temperature_LH02=Data['Temperature']

Data = pd.read_csv('Temperature_LH01.csv.')
Temperature_LH01=Data['Temperature']

Data = pd.read_csv('Temperature_LH005.csv.')
Temperature_LH005=Data['Temperature']

Data = pd.read_csv('Time_LH02.csv.')
Time_LH02=Data['Time']

Data = pd.read_csv('Time_LH01.csv.')
Time_LH01=Data['Time']

Data = pd.read_csv('Time_LH005.csv.')
Time_LH005=Data['Time']

###crystallinity speed calculation
for i in range(0,5):
    if i == 0: 
        Degree='80'
    if i == 1: 
        Degree='90'
    if i == 2: 
        Degree='100'
    if i == 3: 
        Degree='110'
    if i == 4: 
        Degree='120'

    Data = pd.read_csv('crystallization '+Degree+' deg.csv.')
    Time=Data['Time']*60
    HF=-Data['Heat Flow Endo Down']
    
    index=len(HF)
    offset=HF[index-1]
    
    for j in range(0,len(HF)):
        if offset>0:
            HF[j]=HF[j]-offset
        else:
            HF[j]=HF[j]-offset

    Area=[]
    TimeArea=[]
    AreaDx=[]
    Compteur=0
    Crystal=[]
    for j in range(0,len(HF)-1):
        AreaDx.append(HF[j])
        AreaDx.append(HF[j+1])
        dx=Time[j+1]-Time[j]
        ProvArea=trapz(AreaDx, dx=dx)
        if ProvArea>0:
            Area.append(Compteur+ProvArea)
            Compteur=Compteur+ProvArea
        else:
            Area.append(Compteur-ProvArea)
            Compteur=Compteur-ProvArea
        TimeArea.append(Time[j])
        AreaDx=[]
    
    ##relative crystalinity
    for j in range(0,len(Area)):
        Crystal.append(Area[j]/Compteur) 


    ##Cacul of speed of crystallization
    for j in range(0,len(Crystal)-1):
        dt=TimeArea[j+1]-TimeArea[j]
        dXc=Crystal[j+1]-Crystal[j]

        if i ==0:
            Speed_Crystal_T80.append((0.56*dXc)/(dt))
            Speed_Time_T80.append(TimeArea[j])
            Rel_Cryst_T80.append(Crystal[j])
            Rel_TCryst_T80.append(TimeArea[j])
           
                
        elif i ==1: 
            Speed_Crystal_T90.append((0.56*dXc)/(dt))
            Speed_Time_T90.append(TimeArea[j])
            Rel_Cryst_T90.append(Crystal[j])
            Rel_TCryst_T90.append(TimeArea[j])
       
                
        elif i ==2: 
            Speed_Crystal_T100.append((0.56*dXc)/(dt))
            Speed_Time_T100.append(TimeArea[j])
            Rel_Cryst_T100.append(Crystal[j])
            Rel_TCryst_T100.append(TimeArea[j])

                
        elif i ==3: 
            Speed_Crystal_T110.append((0.56*dXc)/(dt))
            Speed_Time_T110.append(TimeArea[j])
            Rel_Cryst_T110.append(Crystal[j])
            Rel_TCryst_T110.append(TimeArea[j])

                
        elif i ==4: 
           # if (j+1)%3==0 or j==0:
            Speed_Crystal_T120.append((0.56*dXc)/(dt))
            Speed_Time_T120.append(TimeArea[j])
            Rel_Cryst_T120.append(Crystal[j])
            Rel_TCryst_T120.append(TimeArea[j])
  
    

for T in range(0,3):

    Crystallization=[]
    Time_Crystallization=[]
    AVC=[]
    Crystal_EF=0
    interval=1
    Time_Tot=0
    Compteur=0
    N=0.01
    
    if T ==0:
        G=len(Temperature_LH02)
        X=Time_LH02
        Y=Temperature_LH02
        K='0.2 mm'
    elif T==1:
        G=len(Temperature_LH01)
        X=Time_LH01
        Y=Temperature_LH01
        K='0.1 mm'
    else:
        G=len(Temperature_LH005)
        X=Time_LH01
        Y=Temperature_LH005
        K='0.05 mm'
     
    for i in range(0,G-interval):
        Contact=0
        AVc=0
        
        if T ==0:
            G=len(Temperature_LH02)
            label='LH0.2 mm'
            ##Definition of the input 
            Input_Temperature=(Temperature_LH02[i+interval]+Temperature_LH02[i])/2
            Input_Time=Time_LH02[i+interval]-Time_LH02[i]
            color='b'
            
        elif T==1:
            G=len(Temperature_LH01)
            label='LH 0.1mm'
            ##Definition of the input 
            Input_Temperature=(Temperature_LH01[i+interval]+Temperature_LH01[i])/2
            Input_Time=Time_LH01[i+interval]-Time_LH01[i]
            color='r'
            
        else:
            G=len(Temperature_LH005)
            label='LH 0.05 mm'
            ##Definition of the input 
            Input_Temperature=(Temperature_LH005[i+interval]+Temperature_LH005[i])/2
            Input_Time=Time_LH005[i+interval]-Time_LH005[i]
            color='k'
       
        ##Definition of the function
        #Selection of the curve
        if Input_Temperature>83 and Input_Temperature<85:
            Curve_SC=Speed_Crystal_T80
            Curve_Time=Speed_Time_T80
            Curve_Crystal=Rel_Cryst_T80
            index1=min(range(len(Rel_Cryst_T80)), key=lambda j: abs(Rel_Cryst_T80[j]-Crystal_EF))
            index2=index1+1
            Compteur=Compteur+1
            Contact=1
        
        if Input_Temperature>=85 and Input_Temperature<95:
            Curve_SC=Speed_Crystal_T90
            Curve_Time=Speed_Time_T90
            index1=min(range(len(Rel_Cryst_T90)), key=lambda j: abs(Rel_Cryst_T90[j]-Crystal_EF))
            index2=index1+1
            Compteur=Compteur+1
            Contact=1
        elif Input_Temperature>=95 and Input_Temperature<105:
            Curve_SC=Speed_Crystal_T100
            Curve_Time=Speed_Time_T100
            index1=min(range(len(Rel_Cryst_T100)), key=lambda j: abs(Rel_Cryst_T100[j]-Crystal_EF))
            index2=index1+1
            Compteur=Compteur+1
            Contact=1
        elif Input_Temperature>=105 and Input_Temperature<115:
            Curve_SC=Speed_Crystal_T110
            Curve_Time=Speed_Time_T110
            index1=min(range(len(Rel_Cryst_T110)), key=lambda j: abs(Rel_Cryst_T110[j]-Crystal_EF))
            index2=index1+1
            Compteur=Compteur+1
            Contact=1
        elif Input_Temperature>=115 and Input_Temperature<125:
            Curve_SC=Speed_Crystal_T120
            Curve_Time=Speed_Time_T120
            index1=min(range(len(Rel_Cryst_T120)), key=lambda j: abs(Rel_Cryst_T120[j]-Crystal_EF))
            index2=index1+1
            Compteur=Compteur+1
            Contact=1
        
        
        #Find the data for dt
        if Compteur==1 and Contact==1 and Input_Temperature<125 and Input_Temperature>83:
            
            Vc1=Curve_SC[index1]
            Vc2=Curve_SC[index2]
            AVc=(Vc2+Vc1)/2
            AVC.append(AVc)
            Crystal_EF=Crystal_EF+(Input_Time*AVc)
            Crystallization.append(Crystal_EF*100)
            
            Time_Tot=Time_Tot + Input_Time
            Time_Crystallization.append(Time_Tot)
            
        elif Contact==1 and Compteur>1 and Input_Temperature<125 and Input_Temperature>83:
                 
            Vc1=Curve_SC[index1]
            Vc2=Curve_SC[index2]
            AVc=(Vc2+Vc1)/2
            AVC.append(AVc)
            
            Crystal_EF=Crystal_EF+(Input_Time*AVc)
            Crystallization.append(Crystal_EF*100)
            
            Time_Tot=Time_Tot + Input_Time
            Time_Crystallization.append(Time_Tot)
        
            
        if Input_Temperature<83 or Input_Temperature>125:
            Crystallization.append(Crystal_EF*100)
            Time_Tot=Time_Tot + Input_Time
            Time_Crystallization.append(Time_Tot)
            AVC.append(AVc)
 
    fig, ax1 = plt.subplots()
    TX1=[-1,60]
    TY1=[83,83]
    TY2=[125,125]
    TX3=[4,4]
    TX3=[8,8]
    TX3=[12,12]
    TX3=[16,16]
    
    TY3=[60,230]
    color = 'tab:blue'
    ax1.set_xlabel('time (s)')
    ax1.set_ylabel('Temperature (TL) (°C)', color=color)
    ax1.set_title(K+' Layer height')
    ax1.set_ylim(50,225)
    ax1.set_xlim(-1,20)
    ax1.set_xticks(np.arange(0, 21, 4))
    ax1.plot(X,Y,color,linewidth=0.8, label='Temperature')
    ax1.plot(TX1,TY1,'k',linewidth=0.8)
    ax1.plot(TX1,TY2,'k',linewidth=0.8)
    ax1.tick_params(axis='y', labelcolor=color)
    
    ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis
    
    color = 'tab:red'
    ax2.set_ylabel('Crystallization (%)', color=color)  # we already handled the x-label with ax1
    ax2.set_ylim(-1,20)
    ax2.plot(Time_Crystallization,Crystallization,color,linewidth=0.8,label=label)
    ax2.tick_params(axis='y', labelcolor=color)
    
    fig.tight_layout()  # otherwise the right y-label is slightly clipped
    plt.savefig('Crystallization'+str(T)+'.pdf',bbox_inches='tight')
    plt.show()
    if i ==G-2:
        print(K+'='+str(Crystal_EF))



plt.plot( Speed_Time_T110, Speed_Crystal_T110)
