using Plots
gr()
using DataFrames,CSV

using LsqFit
r2=0

#%% Load the data
include("/home/kevin/Documents/1 - PhD Projects/2 - Confidence Attractor Network/ConfianceAttractorNetwork/src/decision_making.jl")
#include("/home/kevin/Documents/1 - PhD Projects/2 - Confidence Attractor Network/src/reaction_time_analysis_functions.jl")
include("/home/kevin/Documents/1 - PhD Projects/2 - Confidence Attractor Network/ConfianceAttractorNetwork/src/distributions.jl")

function sim_mf_decisionmaking_pattern_fastv2(coh_temp::Array{Float64},Ir::Float64,If::Float64,
  param::Dict{String,Float64},fixed_param::Dict{String,Float64},simu_param::Dict{String,Float64})



  # FI curve parameters - adapted from Wanf fit
  a=fixed_param["a"]
  b=fixed_param["b"]
  d=fixed_param["d"]
  falset= param["false_time"]
  thr = param["Threshold"]


  # time variables of the  simulation

  dt=simu_param["dt"]
  time_wind=simu_param["time_wind"] # mean window of the variables
  slide_wind = simu_param["slide_wind"]
  false_time=param["false_time"]
  falsetIf=param["false_timeIf"]


  #---- Intialise and vectorise variables to be used in loops below ------

  # Final Variables
  res=zeros(length(coh_temp))
  reaction_time=zeros(length(coh_temp))
  diffS=zeros(length(coh_temp))
  S1=zeros(length(coh_temp))
  S2=zeros(length(coh_temp))
  R1=zeros(length(coh_temp))
  R2=zeros(length(coh_temp))

  ###### Simu variables aprameters
  noise_amp=param["noise_amp"]

  mu0=  param[  "mu0"]

  I0E1p=  param[  "I0E1"]
  I0E2p =  param[  "I0E2"]
  JN11 =  param[ "JN11"]
  JN22  = param[ "JN22"]
  JN12 =  param[ "JN12"]
  JN21 =  param[ "JN21"]
  JAext  = param[ "JAext"]
  Tnmda = param[  "Tnmda"]
  gamma = param[  "gamma"]
  Tstim = param[  "Tstim"]
  Tampa = param[  "Tampa"]
  coh=  param[ "coh"]

  false_timeIf=  param[  "false_timeIf"]

  #####

  # Temporary variables

    Isyn1=0.0
  Isyn2=0.0
  Trest=param["T_rest"]
  s1=0.1
  s2=0.1
  phi1=0.0
  phi2=0.0
  I_eta1=noise_amp * randn()
  I_eta2=noise_amp * randn()

  nu1=Float64[2.0]
  nu2=Float64[2.0]
  k=0 # variable for mod 4 in the mean
  th=false
  result=true

  for (j,coh) in enumerate(coh_temp) # loop on trials number

    th=false
    result=true
    t=0

    nu1=Float64[nu1[end]]
    nu2=Float64[nu2[end]]
    k=0
    while th==false # while we do not cross the threshold
      t=t+dt
      k=k+1

      if t<Trest && result==true
        I0E1=I0E1p -Ir*exp(-t/falset)
        I0E2=I0E2p -Ir*exp(-t/falset)

        I_stim1=0
        I_stim2=0

      elseif t<Trest
        I0E1=I0E1p -Ir*exp(-t/falset) -If*exp(-t/falsetIf)
        I0E2=I0E2p -Ir*exp(-t/falset) -If*exp(-t/falsetIf)
        I_stim1=0
        I_stim2=0

      elseif t<Trest+dt/2 && t>Trest-dt/2
        I0E1=I0E1p
        I0E2=I0E2p
        I_stim1 = JAext*mu0*(1+coh/100)
        I_stim2 = JAext*mu0*(1-coh/100)


      else
        I0E1=I0E1p
        I0E2=I0E2p
        I_stim1 = JAext*mu0*(1+coh/100)
        I_stim2 = JAext*mu0*(1-coh/100)
      end


      Isyn1=JN11*s1-JN12*s2+I_stim1 + I_eta1
      Isyn2=JN22*s2-JN21*s1+I_stim2 + I_eta2

      phi1=(a*Isyn1-b)/(1-exp(-d*(a*Isyn1-b)))
      phi2=(a*Isyn2-b)/(1-exp(-d*(a*Isyn2-b)))

      #---- Dynamical equations -------------------------------------------
      # Mean NMDA-mediated synaptic dynamics updating
      s1 = s1 + dt*(-s1/Tnmda + (1-s1)*gamma*nu1[end]/1000);


      s2 = s2 + dt*(-s2/Tnmda + (1-s2)*gamma*nu2[end]/1000);


      # Ornstein-Uhlenbeck generation of noise in pop1 and 2
      I_eta1 = I_eta1 + (dt/Tampa)*(I0E1-I_eta1) + sqrt(dt/Tampa)*noise_amp*randn() ;
      I_eta2 = I_eta2 + (dt/Tampa)*(I0E2-I_eta2) + sqrt(dt/Tampa)*noise_amp*randn() ;

      # To ensure firing rates are always positive. Large noise amplitude
      # may result in unwanted negative values
      if phi1 < 0
        push!(nu1,0)
        phi1 = 0;
      else
        push!(nu1,phi1)
      end;
      if phi2 < 0
        push!(nu2,0)
        phi2 = 0;
      else
        push!(nu2,phi2)
      end;

      if mod(k,4)==0

        # Test the crossing of a threshold
        @views m1= mean(nu1[k-3:k])
        @views m2 = mean(nu2[k-3:k])
        if  t>Trest&& m1>thr
          th=true
          if coh>=0
          @inbounds  res[j]=1
            result=true
          else
            @inbounds  res[j]=-1
            result=false
          end

        elseif t>Trest&& m2>thr
          th = true
          if coh>=0
            @inbounds  res[j]=-1
            result=false
          else
            @inbounds  res[j]=1
            result=true
          end

        end
      end
    end
  @inbounds  reaction_time[j]=t-Trest
  @inbounds    diffS[j]=s2-s1
  @inbounds    R1[j]=nu1[end]
  @inbounds  R2[j]=nu2[end]
  @inbounds    S1[j]=s1
  @inbounds  S2[j]=s2
    th=false
  end

  result=DataFrame()
  result[:RTs]=reaction_time
  result[:NotError]=res
  result[:Distri]=coh_temp
  result[:DiffS]=diffS
  result[:S1]=S1
  result[:S2]=S2
  result[:R1]=R1
  result[:R2]=R2
  result[:LastTrial]=vcat([0],res[1:end-1])
  result[:LastCoh]=vcat([0],coh_temp[1:end-1])
  return result
end

function list_elements(l)
      a = Float64[]
        i=1
        while (length(a) <20)
           if any(a.==l[i])

            else
                push!(a,l[i])
            end
            i=i+1
        end
        return sort(a)
    end


    using Distributions, StatPlots,Colors

colorsT=[:red,:blue,:darkorange,:violet,:green]
colorsT=colormap("Blues",5)
o=4

param0=zeros(40)
paramR0=zeros(40)
paramR1=zeros(40)
param1=zeros(40)
xAcc=zeros(160)
xTr=zeros(160)
f=0
u=1
fig=plot(legend=false)

for Tr in [5000]#[200,500,1500,2500,5000]
    o=o
    c_min=-51
     c_max = 51
     number = 20

       distri = make_distributions(c_min,c_max,number,10000)

     param["T_rest"] = 1000
     Ir=0.035
     If=0.0

      r2= sim_mf_decisionmaking_pattern_fastv2(distri,Ir,If,
        param,fixed_param,simu_param)
        j=0
Accrep=Float64[]
Accalt=Float64[]
for a in list_elements(r2[:Distri])
    j=j+1

for i=2:10000
    if r2[:Distri][i]==a
        if r2[:NotError][i-1]*sign(r2[:Distri][i-1])==r2[:NotError][i]*sign(r2[:Distri][i])

            push!(Accrep,r2[:RTs][i])
        else
            push!(Accalt,r2[:RTs][i])
        end
    end
end


end
histogram!(Accalt,color=:forestgreen,alpha=0.4,bins=100,normed=true)
histogram!(Accrep,color=:darkorange,alpha=0.4,bins=100,normed=true)
print(mean(Accalt)-mean(Accrep))

#savefig(fig,"./LeftRightIndecAcc$Tr.svg")

end
fig
xlims!((0,1000))
savefig(fig,"./Hist0355000.svg")

fig



#### sauvegarde RSI energy package


# /home/kevin/Documents/1 - PhD Projects/1 - Sequential Decision making/3-SequentialEffects/2-Data/testsAltRep
param = Dict(
    "mu0"=>30.0,
    "noise_amp"=>0.02,
    "I0E1"=>0.3255,
    "I0E2"=>0.3255,
    "JN11"=>0.2609,
    "JN22"=>0.2609,
    "JN12"=>0.0497,
    "JN21"=>0.0497,
    "JAext" => 0.00052,
    "Tnmda" => 100.0,
    "Tampa" => 2.0,
    "gamma" => 0.641,
    "Tstim"=>5000.0,
    "T_rest"=>700.0,
    "coh"=>5.0,
    "false_time"=>200.0,
    "false_timeIf"=>500.0,
    "Threshold"=>20.0
    )
    fixed_param=Dict("a"=>270.0,"b"=>108.0,"d"=>0.1540)
  simu_param=Dict("dt"=>0.5,"time_wind"=>4.0,"slide_wind"=>4.0) #4=2/dt
Icdlist=[0.035,0.04,0.05,0.06,0.07,0.08,0.09,0.1,0.2,0.3]
#Icdlist=[0.2,0.3]#
Trlist=[500,1000,1500,2000,2500,5000,10000]

for Icd in Icdlist
  for Tr in Trlist

    c_min=-51
     c_max = 51
     number = 20

       distri = make_distributions(c_min,c_max,number,10000)

     param["T_rest"] = Tr
     Ir=Icd
     If=0.0

      r2= sim_mf_decisionmaking_pattern_fastv2(distri,Ir,If,
        param,fixed_param,simu_param)
        CSV.write("/home/kevin/Documents/1 - PhD Projects/1 - Sequential Decision making/3-SequentialEffects/2-Data/testsAltRep/EnergyTest/10000RepAlt-T$Tr-Ir$Icd.csv",r2)

  end
end
