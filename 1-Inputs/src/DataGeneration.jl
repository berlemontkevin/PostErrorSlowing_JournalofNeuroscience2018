

# Simulations parameters
include("./decision_making.jl")
include("./reaction_time_analysis_functions.jl")
include("./distributions.jl")


param = Dict("T_rest"=>10000,
"mu0"=>30,
"noise_amp"=>0.02,
"I0E1"=>0.3255,
"I0E2"=>0.3255,
"JN11"=>0.2609,
"JN22"=>0.2609,
"JN12"=>0.0497,
"JN21"=>0.0497,
"JAext" => 0.00052,
"Tnmda" => 100,
"Tampa" => 2,
"gamma" => 0.641,
"Tstim"=>5000,
"coh"=>5,
"false_time"=>200,
"Threshold"=>20
)

fixed_param=Dict("a"=>270,"b"=>108,"d"=>0.1540)
#Parameters of the simulation for EUler methods
simu_param=Dict("dt"=>0.5,"time_wind"=>4,"slide_wind"=>4) #4=2/dt

using FileIO
using DataFrames
using CSV

#a=2
Tr=500
println(Tr)
# de 0.02 à 0.025 déjàa sauvegardé (voir 0.02,0.3 50)
df=DataFrame()
df[:Ir]=linspace(0.03,0.1,30)
for j=2:10
#CSV.write("Irsmall.csv",df)
 for i=1:25

for Ir in [0.035,0.047]
println(Ir)
Cohlist=[2,10,20]


 #Need to do the quantile analysis after with r s etc ...




for c in Cohlist[1:end]
c_min=-c
 c_max = c
 number = 2

   distri = make_distributions(c_min,c_max,number,3000)

 param["T_rest"] = Tr
 If=0

  result,r2= sim_mf_decisionmaking_pattern(distri,Ir,If,
    param,fixed_param,simu_param)

 CSV.write("./datatempAll/V$j-Coh$c-Final200Exp$i-T$Tr-Ir$Ir.csv",result)
 CSV.write("./datatempAll/V$j-NondynCoh$c-Final200Exp$i-T$Tr-Ir$Ir.csv",r2)

 end
end

end
end
