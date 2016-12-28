using ODE
using PyCall
# PyCall package to call Scipy integrate function
@pyimport scipy.integrate as si


#Baseline Parameters
C_m  =   1.0
#membrane capacitance, in uF/cm^2"""
g_Na = 120.0
#Sodium (Na) maximum conductances, in mS/cm^2"""
g_K  =  36.0
#Postassium (K) maximum conductances, in mS/cm^2"""
g_L  =   0.3
#Leak maximum conductances, in mS/cm^2"""
E_Na =  50.0
#Sodium (Na) Nernst reversal potentials, in mV"""
E_K  = -77.0
#Postassium (K) Nernst reversal potentials, in mV"""
E_L  = -54.387
#Leak channels Nernst reversal potentials, in mV"""


# Time steps
step = collect(range(0,0.01,45000))

#Channel Dynamics functions
alpha_m(V) = 0.1*(V+40.0)/(1.0 - exp(-(V+40.0) / 10.0))

beta_m(V) =  4.0 * exp(-(V+65.0) / 18.0)

alpha_h(V) = 0.07 * exp(-(V+65.0) / 20.0)

beta_h(V) = 1.0/(1.0 + exp(-(V+35.0) / 10.0))

alpha_n(V) = 0.01*(V+55.0)/(1.0 - exp(-(V+55.0) / 10.0))

beta_n(V) = 0.125 * exp(-(V+65) / 80.0)

# To Calculate Inhibitory Sodium current 
I_Na(V, m, h) = g_Na * m^3 * h * (V - E_Na)

# To Calculate Inhibitory Potassium current 
I_K(V, n) = g_K  * n^4 * (V - E_K)

# To Calculate Inhibitory Leak current 
I_L(V) = g_L * (V - E_L)

function fcheck(x,t)
    if t > x
        return 1
    else
        return 0
    end
end
    

I_inj(t) = 10*(fcheck(100,t)) - 10*(fcheck(200,t)) + 35*(fcheck(300,t)) - 35*(fcheck(400,t))

function dALLdt(X,t)
    V,m,h,n = X 
    dVdt = (I_inj(t) - I_Na(V, m, h) - I_K(V, n) - I_L(V)) / C_m
    dmdt = alpha_m(V)*(1.0-m) - beta_m(V)*m
    dhdt = alpha_h(V)*(1.0-h) - beta_h(V)*h
    dndt = alpha_n(V)*(1.0-n) - beta_n(V)*n
#    println(I_inj(t))
#    print(I_Na(V,m,h))
#    print(I_K(V,n))
#    print(I_L(V))
#    print(t+"\t"+I_inj(t)+"\t"+I_Na(V, m, h)+"\t"+I_K(V, n)+"\t"+I_L(V)+"\t"+C_m+"\t"+V+"\t"+m+"\t"+h+"\t"+n+"\t"+alpha_m(V)+"\t"+alpha_h(V)+"\t"+alpha_n(V)+"\t"+beta_m(V)+"\t"+beta_h(V)+"\t"+beta_n(V)+"\t"+dVdt+"\t"+dmdt+"\t"+dhdt+"\t"+dndt)  
#    print(dVdt)
#    print(dmdt)
#    print(dhdt)
#    print(dndt)
#    print("\n")
    return dVdt, dmdt, dhdt, dndt
end

#t,y = ODE.ode23s(dALLdt,[-65, 0.05, 0.6, 0.32],step)
y = [[]]
y = si.odeint(dALLdt, [-65, 0.05, 0.6, 0.32],step)

using Plots

Vh = y[:,1]
mh = y[:,2]
hh = y[:,3]
nh = y[:,4]

# To calculate inhibitory currents for Sodium channel, Potassium channel and Leak currents
xlen = length(Vh)
iNa = Array(Float64,45000) 
iK = Array(Float64,45000) 
iNa = fill(0.0,45000)
iK = fill(0.0,45000)
iL = I_L(Vh)                            # Leak Current
for i = 1:45000
    iNa[i] = I_Na(Vh[i], mh[i], hh[i])  # for Sodium channel
    iK[i] = I_K(Vh[i],nh[i])            # for Potassium channel
end


plot(Vh,label = "V")

plot(mh,label = "m")
plot!(hh,label = "h")
plot!(nh,label = "n")

plot(iNa,label = "INa")
plot!(iK,label = "IK")
plot!(iL,label = "IL")


