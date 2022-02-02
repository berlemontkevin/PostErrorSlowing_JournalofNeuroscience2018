#### Simulations for "Post-error with attractor network" Berlemont & Nadal###
# Author: Kevin Berlemont
# 10 July 2018
# Shared privately with the reviewers.
########################################################
using Distributions,CSV,Query
using DataFrames
using Plots
gr()
srand(1234)

param = Dict("T_rest"=>500.0,
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
"coh"=>5.0,
"false_time"=>200.0,
"Threshold"=>20.0
)

fixed_param=Dict("a"=>270.0,"b"=>108.0,"d"=>0.1540)
#Parameters of the simulation for EUler methods
simu_param=Dict("dt"=>0.5,"time_wind"=>4.0,"slide_wind"=>4.0) #4=2/dt

using FileIO
using DataFrames
using CSV

Tr=500
ICD=0.047
println(Tr)


for i=1:50
    print(i)
    print("\n")
    for Tr in [200,300,400,500,600,700,800,900,1000,1100,1200,1300,1400,1500,1600,
        1700,1800,1900,2000,3000,4000,5000]
        for c in 0.0:1.0:20.0

        c_min=-c
        c_max = c
        number = 2

        distri = make_distributions(c_min,c_max,number,5000)

        param["T_rest"] = Tr
        If=0.0

        result= sim_mf_decisionmaking_pattern_fast(distri,ICD,0.0,
        param,fixed_param,simu_param)

        CSV.write("./I47/Coh$c-Final200Exp$i-T$Tr-Ir$Ir.csv",result)


    end

end
