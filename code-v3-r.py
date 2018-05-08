# -*- coding: utf-8 -*-
"""
Created on Sat Mar 17 18:17:57 2018
模拟不同浓度ATP的Kinesin的运动（MC Simulation）

@author: gg
"""

#导入需要的包
import numpy as np
import matplotlib.pyplot as plt

#设置随机数种子myseed
myseed = 91

#kinesin的状态用K表示
K = 113 
#K是一个三位数字，百位代表后头的状态（1-ATP;2-ADP;3-空）
#十位代表Kinesin的状态（1-two head bound;2-中间态）
#个位代表前头的状态（1-ATP;2-ADP;3-空）

#设置参数值
kplus = 104 #后头ATP->ADP速率；s-1
kminus = 8 #前头ATP->ADP速率；s-1
kb = 1.2 #ATP结合二阶速率；s-1uM-1
ATP_conc = 2000 #ATP浓度；uM
F = 7.5 #加载力；pN
F_all = [3.5]
r0 = 900 #0pN step ratio
Fs = 8 #Stop Force;pN
r = r0 ** (1-F/Fs) #step ratio
Pe = r / (kplus/kminus + r) #向前概率
d = 8 #step length;nm

#kinesin运动的时间由T表示
T = 0.0
Time = []

#kinesin的后头位置由X表示
X = 0
P = 0.0#两头平均位置
Position = []

#走的步数
Nstep = 0

#消耗的ATP个数
NofATP = 0
#Number = []

#模拟次数N
N = 1000

#记录每次模拟的runlength--R
#R = []

#记录每次模拟的速度V
V = []
V_F = []

#记录每次模拟步数
NATP = []
NATP_F = []

#每次走的步数
Number = []
Number_F = []

#模拟开始

for F_index in range(len(F_all)):
    F = F_all[F_index]
    r = r0 ** (1-F/Fs) #step ratio
    Pe = r / (kplus/kminus + r) #向前概率
    print(F)
    V = []
    NATP = []
    Number = []
    np.random.seed(myseed)
    
    for index in range(N):
        K = 113 #初始状态：T--0
        T = 0.0
        X = 0
        Time = []
        Position = []
        NofATP = 0
        Nstep = 0
        Time.append(T)
        P = X+d/2
        Position.append(P)
        while 1:        
            K_now = K #记录当前状态
            if K_now == 113: #T--0
                kT = kplus + kb*ATP_conc
                T = T + (np.random.exponential(1.0/kT))
                decision_ch = np.random.random() #选择
                if decision_ch <= kplus / kT: #后头水解
                    K = 3221 #0D;尾数1代表向前
                    X = X+d
                    P = X
                    NofATP = NofATP+1
                else: #前头结合ATP
                    K = 111 #T--T
                    X = X
                    P = X+d/2
                Time.append(T)
                Position.append(P)
                    
            if K_now == 3221: #0D;尾数1代表向前
                kT = kb*ATP_conc
                T = T + (np.random.exponential(1.0/kT))
                decision_ch = np.random.random() #选择
                if decision_ch <= Pe: #Forward
                    K = 113 #T--0
                    X = X
                    P = X+d/2
                    Nstep = Nstep+1
                else: #Back
                    K = 311 #0--T
                    X = X-d
                    P = X+d/2
                Time.append(T)
                Position.append(P)
                    
            if K_now == 3222: #0D;尾数2代表向后
                kT = kb*ATP_conc
                T = T + (np.random.exponential(1.0/kT))
                decision_ch = np.random.random() #选择
                if decision_ch <= Pe: #Forward
                    K = 113 #T--0
                    X = X
                    P = X+d/2
                else: #Back
                    K = 311 #0--T
                    X = X-d
                    P = X+d/2
                    Nstep = Nstep+1
                Time.append(T)
                Position.append(P)
                    
            if K_now == 111: #T--T
                kT = kplus + kminus
                T = T + (np.random.exponential(1.0/kT))
                decision_ch = np.random.random() #选择
                if decision_ch <= kplus / kT: #后头水解
                    K = 1221 #TD;尾数1代表向前
                    X = X+d
                    P = X
                    NofATP = NofATP+1
                else: #前头水解
                    K = 1222 #TD;尾数2代表向后
                    X = X
                    P = X
                    NofATP = NofATP+1
                Time.append(T)
                Position.append(P)
                    
            if K_now == 311: #0--T
                kT = kminus + kb*ATP_conc
                T = T + (np.random.exponential(1.0/kT))
                decision_ch = np.random.random() #选择
                if decision_ch <= kminus / kT: #前头水解
                    K = 3222 #0D;尾数2代表向后
                    X = X
                    P = X
                    NofATP = NofATP+1
                else: #后头结合ATP
                    K = 111 #T--T
                    X = X
                    P = X+d/2
                Time.append(T)
                Position.append(P)
                    
            if K_now == 1221: #TD
                decision_ch = np.random.random() #选择
                if decision_ch <= Pe: #Forward
                    K = 113 #T--0
                    X = X
                    P = X+d/2
                    Nstep = Nstep+1
                else: #Back
                    K = 311 #0--T
                    X = X-d
                    P = X+d/2
                Time.append(T)
                Position.append(P)
                    
            if K_now == 1222: #TD
                decision_ch = np.random.random() #选择
                if decision_ch <= Pe: #Forward
                    K = 113 #T--0
                    X = X
                    P = X+d/2
                else: #Back
                    K = 311 #0--T
                    X = X-d
                    P = X+d/2
                    Nstep = Nstep+1
                Time.append(T)
                Position.append(P)
                    
            if T > 100: #每条轨迹1000s
     #       if NofATP > 10000: #每条轨迹10000ATP
                break
        np.savetxt('Time%d.dat' % (index),Time)
        np.savetxt('Position%d.dat' % (index),Position)
                
        V.append(P/T)
        NATP.append(P)
            
        if Nstep != 0:
            Number.append(NofATP/Nstep)
    #    R.append(np.sum(X))
    V_F.append(np.mean(V))
    print(np.mean(V))
    NATP_F.append(np.var(NATP)/d/np.mean(NATP))
    print(np.var(NATP)/d/np.mean(NATP))
    Number_F.append(np.mean(Number))
    print(np.mean(Number))
    TT =[]
    PP =[]
    Var=[]
    
    for index in range(50):
        t=index*2
        TT.append(t)
        PP = []
        for N_index in range(1000):
            Time = np.loadtxt('Time%d.dat' % (N_index))
            Position = np.loadtxt('Position%d.dat' % (N_index))
            for i in range(len(Time)):
                if Time[i]>=t:
                    PP.append(Position[i])
                    break
        Var.append(np.var(PP))
        
    plt.plot(TT,Var)
    np.savetxt('TT.dat',TT)
    np.savetxt('Var.dat',Var)
    
    z1 = np.polyfit(TT, Var, 1)  #一次多项式拟合，相当于线性拟合
    p1 = np.poly1d(z1)
    print(p1)
    print(p1[1]/8/np.mean(V))
#print(np.mean(V))
#print(np.std(V))
#print(np.mean(NATP))
#print(np.std(NATP))
