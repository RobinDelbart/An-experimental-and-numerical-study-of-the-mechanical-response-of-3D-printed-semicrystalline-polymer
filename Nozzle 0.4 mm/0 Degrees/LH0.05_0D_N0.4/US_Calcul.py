import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import statistics
UStress=[]
for i in range(1,4):
    Data = pd.read_csv(f'Test {i}.csv.')
    W0=13
    T0=3.2
    S0=W0*T0
    Stress=((Data['Load ']*1000)/S0)
    
    UStress.append(max(Stress))

T=(UStress[0]+UStress[1]+UStress[2])/3
print(T)

