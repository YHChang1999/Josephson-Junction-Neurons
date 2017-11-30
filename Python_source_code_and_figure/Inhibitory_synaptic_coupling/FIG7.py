#引入標頭檔
import matplotlib.pyplot as plt     #畫圖用
import math     #三角函數用
#

#參數
start=0     #開始時間
finish=1000     #停止時間
I_b1=-1.9
I_b2=1.76
Lambda=0.1
Lambda_s1=0.5
Lambda_p1=0.5
Lambda_s2=0.65
Lambda_p2=0.35
Lambda_c=0
Gamma=2.0
eta=1
Omega=1.0
Q=0.05
Lambda_syn=0.3
r12=0.6
startval=[0,0,0,0,0,0,0,0,0,0,0]    #initial conditions
step=0.005  #每一個數值解的間距
def I_in(t):
    return 0.33*((t>=175)and(t<=675))
#論文上面是用0.3，可能因為迭代之後會有些許誤差，因此使用0.3時圖會長的不一樣，此處使用0.33來代替。
#

#微分方程
def fy1(t,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11):
    return y2
def fy2(t,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11):
    return -Gamma*y2-math.sin(y1)-Lambda*(y1+y3)+Lambda_s1*I_in(t)+(1-Lambda_p1)*I_b1-(y11+Lambda*y9/(Lambda_syn*Omega*Omega))  #這裡有加上第五頁所提及的back-action
def fy3(t,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11):
    return y4
def fy4(t,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11):
    return -Gamma*y4-math.sin(y3)+(-Lambda*(y1+y3)+Lambda_s1*I_in(t)-Lambda_p1*I_b1)/eta
def fy5(t,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11):
    return y6
def fy6(t,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11):
    return -Gamma*y6-math.sin(y5)-Lambda*(y5+y7)+Lambda_s2*y11+(1-Lambda_p2)*I_b2
def fy7(t,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11):
    return y8
def fy8(t,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11):
    return -Gamma*y8-math.sin(y7)+(-Lambda*(y5+y7)+Lambda_s2*y11-Lambda_p2*I_b2)/eta
def fy9(t,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11):
    return y10
def fy10(t,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11):
    return (Omega*Omega)*((-Q*y10/Omega)-y9+y2-(Q*Omega*Lambda_syn/Lambda)*y11-(1/(1-Lambda_syn))*(-r12*y11/Gamma+y9-Lambda_syn*(y6+y8)))
def fy11(t,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11):
    return (Lambda/(Lambda_syn*(1-Lambda_syn)))*((-r12*y11/Gamma)+y9-Lambda_syn*(y6+y8))
#

#Runge-Kutta 係數初始化
#參照書上的方法
b=[1/6,1/3,1/3,1/6]
c=[[0,0,0,0],[0.5,0,0,0],[0,0.5,0,0],[0,0,1,0]]
d=[0,0.5,0.5,1]
k=[0,0,0,0]
order=4
#

#數值初始化
y1=startval[0]
y2=startval[1]
y3=startval[2]
y4=startval[3]
y5=startval[4]
y6=startval[5]
y7=startval[6]
y8=startval[7]
y9=startval[8]
y10=startval[9]
y11=startval[10]
t=start
ck=0    #迭代用的參數
#

#最後呈現的曲線陣列
Flux1=[(y1+y3)*Lambda-0.7]      #這裡有做平移讓所有曲線可以在同一張圖裡呈現
Flux2=[(y5+y7)*Lambda]
tvals=[start]
I_of_t=[I_in(start)/3-0.9]      #這裡有做縮放及平移讓所有曲線可以在同一張圖裡呈現

#

#Runge-Kutta Classical
#之前書上給的方法中使用最經典的方法
j=start+step
while j<=finish:
    #對1
    k[0]=step*fy1(t,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11)
    for i in range (1,order):
        ck=0
        for l in range(0,i):
            ck+=c[i][l]*k[l]
        k[i]=step*fy1(t+step*d[i],y1+ck,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11)
    y1tmp=y1
    for i in range (0,order):
        y1tmp+=b[i]*k[i]

    #對2
    k[0]=step*fy2(t,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11)
    for i in range (1,order):
        ck=0
        for l in range(0,i):
            ck+=c[i][l]*k[l]
        k[i]=step*fy2(t+step*d[i],y1,y2+ck,y3,y4,y5,y6,y7,y8,y9,y10,y11)
    y2tmp=y2
    for i in range (0,order):
        y2tmp+=b[i]*k[i]

    #對3
    k[0]=step*fy3(t,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11)
    for i in range (1,order):
        ck=0
        for l in range(0,i):
            ck+=c[i][l]*k[l]
        k[i]=step*fy3(t+step*d[i],y1,y2,y3+ck,y4,y5,y6,y7,y8,y9,y10,y11)
    y3tmp=y3
    for i in range (0,order):
        y3tmp+=b[i]*k[i]

    #對4
    k[0]=step*fy4(t,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11)
    for i in range (1,order):
        ck=0
        for l in range(0,i):
            ck+=c[i][l]*k[l]
        k[i]=step*fy4(t+step*d[i],y1,y2,y3,y4+ck,y5,y6,y7,y8,y9,y10,y11)
    y4tmp=y4
    for i in range (0,order):
        y4tmp+=b[i]*k[i]

    #對5
    k[0]=step*fy5(t,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11)
    for i in range (1,order):
        ck=0
        for l in range(0,i):
            ck+=c[i][l]*k[l]
        k[i]=step*fy5(t+step*d[i],y1,y2,y3,y4,y5+ck,y6,y7,y8,y9,y10,y11)
    y5tmp=y5
    for i in range (0,order):
        y5tmp+=b[i]*k[i]

    #對6
    k[0]=step*fy6(t,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11)
    for i in range (1,order):
        ck=0
        for l in range(0,i):
            ck+=c[i][l]*k[l]
        k[i]=step*fy6(t+step*d[i],y1,y2,y3,y4,y5,y6+ck,y7,y8,y9,y10,y11)
    y6tmp=y6
    for i in range (0,order):
        y6tmp+=b[i]*k[i]

    #對7
    k[0]=step*fy7(t,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11)
    for i in range (1,order):
        ck=0
        for l in range(0,i):
            ck+=c[i][l]*k[l]
        k[i]=step*fy7(t+step*d[i],y1,y2,y3,y4,y5,y6,y7+ck,y8,y9,y10,y11)
    y7tmp=y7
    for i in range (0,order):
        y7tmp+=b[i]*k[i]

    #對8
    k[0]=step*fy8(t,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11)
    for i in range (1,order):
        ck=0
        for l in range(0,i):
            ck+=c[i][l]*k[l]
        k[i]=step*fy8(t+step*d[i],y1,y2,y3,y4,y5,y6,y7,y8+ck,y9,y10,y11)
    y8tmp=y8
    for i in range (0,order):
        y8tmp+=b[i]*k[i]

    #對9
    k[0]=step*fy9(t,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11)
    for i in range (1,order):
        ck=0
        for l in range(0,i):
            ck+=c[i][l]*k[l]
        k[i]=step*fy9(t+step*d[i],y1,y2,y3,y4,y5,y6,y7,y8,y9+ck,y10,y11)
    y9tmp=y9
    for i in range (0,order):
        y9tmp+=b[i]*k[i]

    #對10
    k[0]=step*fy10(t,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11)
    for i in range (1,order):
        ck=0
        for l in range(0,i):
            ck+=c[i][l]*k[l]
        k[i]=step*fy10(t+step*d[i],y1,y2,y3,y4,y5,y6,y7,y8,y9,y10+ck,y11)
    y10tmp=y10
    for i in range (0,order):
        y10tmp+=b[i]*k[i]

    #對11
    k[0]=step*fy11(t,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11)
    for i in range (1,order):
        ck=0
        for l in range(0,i):
            ck+=c[i][l]*k[l]
        k[i]=step*fy11(t+step*d[i],y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11+ck)
    y11tmp=y11
    for i in range (0,order):
        y11tmp+=b[i]*k[i]
    #上面都在做迭代
    t1=t+step
    #將迭代後的數值加入最後的陣列
    tvals.append(t1)
    I_of_t.append(I_in(t1)/3-0.9)
    t=t1
    y1=y1tmp
    y2=y2tmp
    y3=y3tmp
    y4=y4tmp
    y5=y5tmp
    y6=y6tmp
    y7=y7tmp
    y8=y8tmp
    y9=y9tmp
    y10=y10tmp
    y11=y11tmp
    j+=step
    Flux1.append((y1+y3)*Lambda-0.7)
    Flux2.append((y5+y7)*Lambda)
#

#畫圖

plt.plot(tvals,Flux1)
plt.plot(tvals,Flux2)
plt.plot(tvals,I_of_t)
plt.xlabel("Time(arb.units)")
plt.ylabel("Flux(arb.units)")
plt.ylim(-1,0.6)
plt.title("FIG7")

plt.show()
#plt.savefig("FIG7.png",dpi=300,format="png")
#如果要存檔就把上一行的井字號去掉
#

