#4/6/23 EJZ
#simple 0D NPBD model for ocean surface
#assume no light limitation for time being for pp (since in surface)

#State variables (uM N):
#N -- nutrient (inorganic)
#P -- phytoplankton (photosynthetic primary producer)
#B -- bacteria (heterotrophic remineralizer)
#D -- detritus

#Model

#parameters as struct:
struct Prms
    T::Int64
    dt::Float64
    umax_p::Float64 
    k_n::Float64 
    vmax_b::Float64 
    y_b::Float64 
    k_d::Float64 
    m_q::Float64 
    S_n::Float64
    rsink::Float64
    nIC::Float64
    pIC::Float64
    bIC::Float64
    dIC::Float64
end

#ecosystem model:
function ecof(n, p, b, d, prms)

    u_p = prms.umax_p.*n/(n + prms.k_n) 
    
    u_b = prms.y_b*prms.vmax_b.*d/(d + prms.k_d) 
    
    dndt = prms.S_n - u_p.*p + (1. ./ prms.y_b .- 1.)*u_b.*b
    dpdt = u_p.*p - prms.m_q.*p.^2 
    dbdt = u_b.*b - prms.m_q.*b.^2 
    dddt = prms.m_q.*(p.^2 + b.^2) - (1. ./ prms.y_b)*u_b.*b - prms.rsink.*d

    return dndt, dpdt, dbdt, dddt

end #end ecof function

#running everything:
function runNPBD(prms)

    nt = Int(T/dt) #number of time points

    #set up empty arrays:
    nrec = Int(nt + 1) # # of points to record (includes 0)
    timet = Array{Float64,1}(undef, nrec)
    n = Array{Float64,1}(undef, nrec)
    p = Array{Float64,1}(undef, nrec)
    b = Array{Float64,1}(undef, nrec)
    d = Array{Float64,1}(undef, nrec)

    #Initial conditions (IC):
    #record 0th time step 
    timet[1] = 0. 
    n[1] = prms.nIC
    p[1] = prms.pIC
    b[1] = prms.bIC
    d[1] = prms.dIC
    #Set tracers at initial conditions
    ntemp = prms.nIC
    ptemp = prms.pIC
    btemp = prms.bIC
    dtemp = prms.dIC

    for t = 1:nt #start time loop

        dndt, dpdt, dbdt, dddt = ecof(ntemp, ptemp, btemp, dtemp, prms)

        #simplest possible forward integration
        ntemp2 = ntemp .+ dt.*dndt 
        ptemp2 = ptemp .+ dt.*dpdt 
        btemp2 = btemp .+ dt.*dbdt 
        dtemp2 = dtemp .+ dt.*dddt 

        #save new time points:
        ntemp = ntemp2
        ptemp = ptemp2
        btemp = btemp2
        dtemp = dtemp2

        #record:
        timet[t+1] = t.*dt 
        n[t+1] = ntemp 
        p[t+1] = ptemp 
        b[t+1] = btemp 
        d[t+1] = dtemp 

        if mod(t, 100) == 0
            println(t.*dt)
        end

    end #end time loop

    return timet, n, p, b, d

end #end runNPBD function

###########################################
#Run the model:

#Time stepping:
T = 100
dt = 0.01

#Microbial parameters:
umax_p = 1
k_n = 0.1
vmax_b = 1
y_b = 0.3
k_d = 0.1
m_q = 1

#Nutrient supply rate and sinking rate:
S_n = 0.1
rsink = 1

#Initial conditions:
nIC = 1
pIC = 0.1
bIC = 0.1
dIC = 0.1

params = Prms(T, dt, umax_p, k_n, vmax_b, y_b, k_d, m_q, S_n, rsink, nIC, pIC, bIC, dIC)

timet, n, p, b, d = runNPBD(params)

###########################################
#Plot:

using Plots

p1 = plot(timet, n, title = "N")
p2 = plot(timet, p, title = "P")
p3 = plot(timet, b, title = "B")
p4 = plot(timet, d, title = "D")

plot(p1, p2, p3, p4,
     linewidth = 2,
     layout = [1,1,1,1],
)


