# Simulations parameters
include("/home/kevin/Documents/1 - PhD Projects/1 - Sequential Decision making/3-SequentialEffects/1-Inputs/src/decision_making.jl")
include("/home/kevin/Documents/1 - PhD Projects/1 - Sequential Decision making/3-SequentialEffects/1-Inputs/src/reaction_time_analysis_functions.jl")
include("/home/kevin/Documents/1 - PhD Projects/1 - Sequential Decision making/3-SequentialEffects/1-Inputs/src/distributions.jl")


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
"false_time"=>300,
"Threshold"=>20
)

fixed_param=Dict("a"=>270,"b"=>108,"d"=>0.1540)
#Parameters of the simulation for EUler methods
simu_param=Dict("dt"=>0.5,"time_wind"=>4,"slide_wind"=>4) #4=2/dt


using DataFrames
using CSV

#a=2
Tr=500
println(Tr)
# de 0.02 à 0.025 déjàa sauvegardé (voir 0.02,0.3 50)


Irlist=[0.035,0.047]
Cohlist=[2,10,20]

using Plots
gr()



############################
## Fig neuronale PES
############################
a=300
PEsameL=zeros(251+a)#
PEdiffL=zeros(251+a)
PEsameW=zeros(251+a)
PEdiffW=zeros(251+a)
PCsameL=zeros(251+a)#
PCdiffL=zeros(251+a)
PCsameW=zeros(251+a)
PCdiffW=zeros(251+a)
numPEsame=0
numPEdiff=0
numPEsame=0
numPEdiff=0
numPCsame=0
numPCdiff=0
numPCsame=0
numPCdiff=0
Ir=Irlist[2]
c=Cohlist[3]

numr=0
numr2=1

for i=1:25
    r=CSV.read("./datatempAll/Coh$c-Final200Exp$i-T$Tr-Ir$Ir.csv")
    r2=CSV.read("./datatempAll/NondynCoh$c-Final200Exp$i-T$Tr-Ir$Ir.csv")

temp = 250
for i=1:2500
temp= temp + Int(r2[:RTs][i]/2)
if r2[:NotError][i] == 1 && r2[:Distri][i] >0 && r2[:Distri][i+1] >0
PCsameL=PCsameL + r[:R2][temp:(temp+250+a)]
PCsameW=PCsameW + r[:R1][temp:(temp+250+a)]

numPCsame=numPCsame+1
#rlist2=rlist2 + r[:R2][temp:(temp+250)]

elseif r2[:NotError][i] == -1 && r2[:Distri][i] <0 && r2[:Distri][i+1] >0
    #rlist=rlist + r[:R2][temp:(temp+250)]
    PEsameL=PEsameL + r[:R2][temp:(temp+250+a)]
    PEsameW=PEsameW + r[:R1][temp:(temp+250+a)]

    numPEsame=numPEsame+1

elseif r2[:NotError][i] == -1 && r2[:Distri][i] >0 && r2[:Distri][i+1] <0
    #rlist=rlist + r[:R2][temp:(temp+250)]
    PEsameL=PEsameL + r[:R1][temp:(temp+250+a)]
    PEsameW=PEsameW + r[:R2][temp:(temp+250+a)]

    numPEsame=numPEsame+1
elseif r2[:NotError][i] == 1 && r2[:Distri][i] <0 && r2[:Distri][i+1] <0
    #rlist=rlist + r[:R2][temp:(temp+250)]
    PCsameL=PCsameL + r[:R1][temp:(temp+250+a)]
    PCsameW=PCsameW + r[:R2][temp:(temp+250+a)]

    numPCsame=numPCsame+1









elseif r2[:NotError][i] == 1 && r2[:Distri][i] >0 && r2[:Distri][i+1] <0 && r2[:NotError][i+1] == 1
    PCdiffL=PCdiffL + r[:R2][temp:(temp+250+a)]
    PCdiffW=PCdiffW + r[:R1][temp:(temp+250+a)]

    numPCdiff=numPCdiff+1
    #rlist2=rlist2 + r[:R2][temp:(temp+250)]

elseif r2[:NotError][i] == -1 && r2[:Distri][i] <0 && r2[:Distri][i+1] <0 && r2[:NotError][i+1] == 1
        #rlist=rlist + r[:R2][temp:(temp+250)]
        PEdiffL=PEdiffL + r[:R2][temp:(temp+250+a)]
        PEdiffW=PEdiffW + r[:R1][temp:(temp+250+a)]

        numPEdiff=numPEdiff+1

    elseif r2[:NotError][i] == -1 && r2[:Distri][i] >0 && r2[:Distri][i+1] >0 && r2[:NotError][i+1] == 1
        #rlist=rlist + r[:R2][temp:(temp+250)]
        PEdiffL=PEdiffL + r[:R1][temp:(temp+250+a)]
        PEdiffW=PEdiffW + r[:R2][temp:(temp+250+a)]

        numPEdiff=numPEdiff+1
    elseif r2[:NotError][i] == 1 && r2[:Distri][i] <0 && r2[:Distri][i+1] >0 && r2[:NotError][i+1] == 1
        #rlist=rlist + r[:R2][temp:(temp+250)]
        PCdiffL=PCdiffL + r[:R1][temp:(temp+250+a)]
        PCdiffW=PCdiffW + r[:R2][temp:(temp+250+a)]

        numPCdiff=numPCdiff+1

end
temp = temp+250
end

end

PEsameL=PEsameL./numPEsame
PEdiffL=PEdiffL./numPEdiff
PEsameW=PEsameW./numPEsame
PEdiffW=PEdiffW./numPEdiff
PCsameL=PCsameL./numPCsame
PCdiffL=PCdiffL./numPCdiff
PCsameW=PCsameW./numPCsame
PCdiffW=PCdiffW./numPCdiff




plot(PEsameL,linewidth=3)
plot!(PCsameL,linewidth=3)

plot(PEdiffL[5:250],linewidth=4,color=:darkorange,line=:dash)
plot!(PCdiffL[5:250],linewidth=4,color=:blue,line=:dash)
plot!(PEdiffW[5:250],linewidth=4,color=:darkorange)
plot!(PCdiffW[5:250],linewidth=4,color=:blue)
plot!(legend=:bottomleft)


savefig("./DiffPEQTrue.svg")


mean(r2[:RTs])

plot!(PEsameW,linewidth=3)
plot!(PCsameW,linewidth=3)


plot(PEsameL[5:250],linewidth=4,color=:darkorange,line=:dash)
plot!(PCsameL[5:250],linewidth=4,color=:blue,line=:dash)
plot!(PEsameW[5:250],linewidth=4,color=:darkorange)
plot!(PCsameW[5:250],linewidth=4,color=:blue)
plot!(legend=:topleft)

savefig("./SamePEQ.svg")


plot(PEdiffL[30:320]-PCdiffL[30:320],linewidth=3,color=:lightgray)
plot!(-PEsameW[30:320]+PCsameW[30:320],linewidth=3,color=:black)
vline!([220],linewidth=4,color=:black,line=:dash)

savefig("./DifferencePES.svg")






############################
## Fig neuronale PES: time from sacccade
############################
a=80
PEsameL=zeros(a+1)#
PEdiffL=zeros(a+1)
PEsameW=zeros(a+1)
PEdiffW=zeros(a+1)
PCsameL=zeros(a+1)#
PCdiffL=zeros(a+1)
PCsameW=zeros(a+1)
PCdiffW=zeros(a+1)
numPEsame=0
numPEdiff=0
numPEsame=0
numPEdiff=0
numPCsame=0
numPCdiff=0
numPCsame=0
numPCdiff=0
Ir=Irlist[2]
c=Cohlist[1]

numr=0
numr2=1

for i=1:25
    r=CSV.read("./datatempAll/Coh$c-Final200Exp$i-T$Tr-Ir$Ir.csv")
    r2=CSV.read("./datatempAll/NondynCoh$c-Final200Exp$i-T$Tr-Ir$Ir.csv")

temp = 250
for i=1:2500
temp= temp + Int(r2[:RTs][i]/2) + 250
if r2[:NotError][i] == 1 && r2[:Distri][i] >0 && r2[:Distri][i+1] >0
PCsameL=PCsameL + r[:R2][temp:temp+a]
PCsameW=PCsameW + r[:R1][temp:temp+a]

numPCsame=numPCsame+1
#rlist2=rlist2 + r[:R2][temp:(temp+250)]

elseif r2[:NotError][i] == -1 && r2[:Distri][i] <0 && r2[:Distri][i+1] >0
    #rlist=rlist + r[:R2][temp:(temp+250)]
    PEsameL=PEsameL + r[:R2][temp:temp+a]
    PEsameW=PEsameW + r[:R1][temp:temp+a]

    numPEsame=numPEsame+1

elseif r2[:NotError][i] == -1 && r2[:Distri][i] >0 && r2[:Distri][i+1] <0
    #rlist=rlist + r[:R2][temp:(temp+250)]
    PEsameL=PEsameL + r[:R1][temp:temp+a]
    PEsameW=PEsameW + r[:R2][temp:temp+a]

    numPEsame=numPEsame+1
elseif r2[:NotError][i] == 1 && r2[:Distri][i] <0 && r2[:Distri][i+1] <0
    #rlist=rlist + r[:R2][temp:(temp+250)]
    PCsameL=PCsameL + r[:R1][temp:temp+a]
    PCsameW=PCsameW + r[:R2][temp:temp+a]

    numPCsame=numPCsame+1









elseif r2[:NotError][i] == 1 && r2[:Distri][i] >0 && r2[:Distri][i+1] <0 && r2[:NotError][i+1] == 1
    PCdiffL=PCdiffL + r[:R2][temp:temp+a]
    PCdiffW=PCdiffW + r[:R1][temp:temp+a]

    numPCdiff=numPCdiff+1
    #rlist2=rlist2 + r[:R2][temp:(temp+250)]

elseif r2[:NotError][i] == -1 && r2[:Distri][i] <0 && r2[:Distri][i+1] <0 && r2[:NotError][i+1] == 1
        #rlist=rlist + r[:R2][temp:(temp+250)]
        PEdiffL=PEdiffL + r[:R2][temp:temp+a]
        PEdiffW=PEdiffW + r[:R1][temp:temp+a]

        numPEdiff=numPEdiff+1

    elseif r2[:NotError][i] == -1 && r2[:Distri][i] >0 && r2[:Distri][i+1] >0 && r2[:NotError][i+1] == 1
        #rlist=rlist + r[:R2][temp:(temp+250)]
        PEdiffL=PEdiffL + r[:R1][temp:temp+a]
        PEdiffW=PEdiffW + r[:R2][temp:temp+a]

        numPEdiff=numPEdiff+1
    elseif r2[:NotError][i] == 1 && r2[:Distri][i] <0 && r2[:Distri][i+1] >0 && r2[:NotError][i+1] == 1
        #rlist=rlist + r[:R2][temp:(temp+250)]
        PCdiffL=PCdiffL + r[:R1][temp:temp+a]
        PCdiffW=PCdiffW + r[:R2][temp:temp+a]

        numPCdiff=numPCdiff+1

end

end

end

PEsameL=PEsameL./numPEsame
PEdiffL=PEdiffL./numPEdiff
PEsameW=PEsameW./numPEsame
PEdiffW=PEdiffW./numPEdiff
PCsameL=PCsameL./numPCsame
PCdiffL=PCdiffL./numPCdiff
PCsameW=PCsameW./numPCsame
PCdiffW=PCdiffW./numPCdiff


using Plots
gr()
plot(PEsameL,linewidth=3,color=:darkorange,line=:dash)
plot!(PCsameL,linewidth=3,color=:blue,line=:dash)


plot!(PEsameW,color=:darkorange,linewidth=3)
plot!(PCsameW,color=:blue,linewidth=3)
plot!(legend=:topleft)

#savefig("./PEQSameMotiononset.svg")
savefig("./NoSameMotiononset.svg")

plot(PEdiffL,linewidth=3,color=:darkorange,line=:dash)
plot!(PCdiffL,linewidth=3,color=:blue,line=:dash)
plot!(PEdiffW,color=:darkorange,linewidth=3)
plot!(PCdiffW,color=:blue,linewidth=3)
plot!(legend=:topleft)
savefig("./NoDiffTrueMotiononset.svg")




c_min=-30
 c_max = 30
 number = 40

   distri = make_distributions(c_min,c_max,number,20000)

 param["T_rest"] = 1500
 Ir=0.035
 If=0

  result,r2= sim_mf_decisionmaking_pattern(distri,Ir,If,
    param,fixed_param,simu_param)



Accrept=[]
Accaltt=[]

for a in linspace(-30,30,40)
    Accrep=Float64[]
    Accalt=Float64[]
for i=2:1000
    if r2[:Distri][i]==a
        if r2[:NotError][i-1]*sign(r2[:Distri][i-1])==1

            push!(Accrep,r2[:NotError][i]*sign(r2[:Distri][i]))
        else
            push!(Accalt,r2[:NotError][i]*sign(r2[:Distri][i]))
        end
    end
end

push!(Accrept,mean(Accrep))
push!(Accaltt,mean(Accalt))

end
r2

scatter(linspace(-30,30,40),0.5+0.5.*Accrept,color=:blue)
scatter!(linspace(-30,30,40),0.5+0.5.*Accaltt,color=:darkorange)

using LsqFit

model(x, p) = 1./(1+exp.(-(p[1]+p[2].*x)))

xdata = linspace(-30,30,40)
ydata = 0.5+0.5.*Accrept
p0 = [0.5, 0.5]
fit = curve_fit(model, xdata, ydata, p0)

fit.param


plot!(linspace(-30,30,40),1./(1+exp.(-(fit.param[1]+fit.param[2].*linspace(-30,30,40)))),color=:blue,label="Repetition",linewidth=3)

fit = curve_fit(model, xdata, 0.5+0.5.*Accaltt, p0)

fit.param


plot!(linspace(-30,30,40),1./(1+exp.(-(fit.param[1]+fit.param[2].*linspace(-30,30,40)))),color=:darkorange,label="Alternation",linewidth=3)
ylabel!("Percentage Rightward choice")
plot!(legend=:topleft)
savefig("./accuracy.pdf")






c_min=-30
 c_max = 30
 number = 20

   distri = make_distributions(c_min,c_max,number,10000)

 param["T_rest"] = 500
 Ir=0.04
 If=0

  result= sim_mf_decisionmaking_pattern_fast(distri,Ir,If,
    param,fixed_param,simu_param)

using LsqFit
r2=0

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

colorsT=[:red,:blue,:darkorange,:violet,:green]
o=1

param0=zeros(10)
paramR0=zeros(10)

paramR1=zeros(10)
param1=zeros(10)
fig=plot(legend=false)
for Tr in [200,500,1500,2500,5000]
    fig=plot()
    o=2
    Accrept=zeros(20,10)
    Accaltt=zeros(20,10)
for u=1:10
    r2 = CSV.read("/home/kevin/Documents/1 - PhD Projects/1 - Sequential Decision making/3-SequentialEffects/2-Data/AccuracySequential/Accuracyc51Tau200Exp$u-T$Tr-Ir0.04.csv")
j=0
for a in list_elements(r2[:Distri])
    j=j+1
    Accrep=Float64[]
    Accalt=Float64[]
for i=2:10000
    if r2[:Distri][i]==a
        if r2[:NotError][i-1]*sign(r2[:Distri][i-1])==1

            push!(Accrep,r2[:NotError][i]*sign(r2[:Distri][i]))
        else
            push!(Accalt,r2[:NotError][i]*sign(r2[:Distri][i]))
        end
    end
end

Accrept[j,u]=mean(Accrep)
Accaltt[j,u]=mean(Accalt)

end
end

scatter!(list_elements(r2[:Distri]),0.5+0.5.*mean(Accrept,2),yerr=sqrt.(var(Accrept,2)),color=colorsT[o])
scatter!(list_elements(r2[:Distri]),0.5+0.5.*mean(Accaltt,2),yerr=sqrt.(var(Accaltt,2)),color=colorsT[o])


model(x, p) = 1./(1+exp.(-(p[1]+p[2].*x)))
for u=1:10
xdata = linspace(-51,51,20)
ydata = 0.5+0.5.*Accrept[:,u]
p0 = [0.5, 0.5]
fit = curve_fit(model, xdata, ydata, p0)

paramR0[u]=fit.param[1]
paramR1[u]=fit.param[2]
plot!(linspace(-51,51,40),1./(1+exp.(-(fit.param[1]+fit.param[2].*linspace(-51,51,40)))),alpha=0.3,label="Repetition",linewidth=3,color=colorsT[o])
ydata = 0.5+0.5.*Accaltt[:,u]
p0 = [0.5, 0.5]
fit = curve_fit(model, xdata, ydata, p0)
plot!(linspace(-51,51,40),1./(1+exp.(-(fit.param[1]+fit.param[2].*linspace(-51,51,40)))),alpha=0.3,label="Repetition",linewidth=3,line=:dash,color=colorsT[o])
param0[u]=fit.param[1]
param1[u]=fit.param[2]
end
plot!(legend=false)
fig
savefig(fig,"./Acc$Tr.svg")
end

#### Left right indecision stuff
# compute param0/param1

histogram(+paramR0./paramR1-param0./param1,bins=10)
-param0./param1



scatter(linspace(-30,30,20),0.5+0.5.*Accrept,color=:blue)
scatter!(linspace(-30,30,20),0.5+0.5.*Accaltt,color=:darkorange)

using LsqFit

model(x, p) = 1./(1+exp.(-(p[1]+p[2].*x)))
for u=1:10
xdata = linspace(-20,20,12)
ydata = 0.5+0.5.*Accrept[:,u]
p0 = [0.5, 0.5]
fit = curve_fit(model, xdata, ydata, p0)

fit.param
plot!(linspace(-20,20,40),1./(1+exp.(-(fit.param[1]+fit.param[2].*linspace(-20,20,40)))),alpha=0.1,color=:blue,label="Repetition",linewidth=3)
end

fig

fit = curve_fit(model, xdata, 0.5+0.5.*Accaltt, p0)

fit.param


plot!(linspace(-30,30,40),1./(1+exp.(-(fit.param[1]+fit.param[2].*linspace(-30,30,40)))),alpha=0.5,color=:darkorange,label="Alternation",linewidth=3)
ylabel!("Percentage Rightward choice")
plot!(legend=:topleft)
savefig("./accuracy.pdf")



#######################################################
##########################################
#Mean stuff for neuronal activity
################################################
#################################################


############################
## Fig neuronale PES
############################
a=300
PEsameL=zeros(251+a,25)#
PEdiffL=zeros(251+a,25)
PEsameW=zeros(251+a,25)
PEdiffW=zeros(251+a,25)
PCsameL=zeros(251+a,25)#
PCdiffL=zeros(251+a,25)
PCsameW=zeros(251+a,25)
PCdiffW=zeros(251+a,25)

j=2
for u=2:25
numPEsame=0
numPEdiff=0
numPEsame=0
numPEdiff=0
numPCsame=0
numPCdiff=0
numPCsame=0
numPCdiff=0
Ir=Irlist[1]
c=Cohlist[2]

numr=0
numr2=1


    r=CSV.read("/home/kevin/Documents/1 - PhD Projects/1 - Sequential Decision making/3-SequentialEffects/2-Data/FiringRates/datatempAll/V$j-Coh$c-Final200Exp$u-T$Tr-Ir$Ir.csv")
    r2=CSV.read("/home/kevin/Documents/1 - PhD Projects/1 - Sequential Decision making/3-SequentialEffects/2-Data/FiringRates/datatempAll/V$j-NondynCoh$c-Final200Exp$u-T$Tr-Ir$Ir.csv")

temp = 250
for i=1:2990
temp= temp + Int(r2[:RTs][i]/2)
if r2[:NotError][i] == 1 && r2[:Distri][i] >0 && r2[:Distri][i+1] >0
PCsameL[:,u]=PCsameL[:,u] + r[:R2][temp:(temp+250+a)]
PCsameW[:,u]=PCsameW[:,u] + r[:R1][temp:(temp+250+a)]

numPCsame=numPCsame+1
#rlist2=rlist2 + r[:R2][temp:(temp+250)]

elseif r2[:NotError][i] == -1 && r2[:Distri][i] <0 && r2[:Distri][i+1] >0
    #rlist=rlist + r[:R2][temp:(temp+250)]
    PEsameL[:,u]=PEsameL[:,u] + r[:R2][temp:(temp+250+a)]
    PEsameW[:,u]=PEsameW[:,u] + r[:R1][temp:(temp+250+a)]

    numPEsame=numPEsame+1

elseif r2[:NotError][i] == -1 && r2[:Distri][i] >0 && r2[:Distri][i+1] <0
    #rlist=rlist + r[:R2][temp:(temp+250)]
    PEsameL[:,u]=PEsameL[:,u] + r[:R1][temp:(temp+250+a)]
    PEsameW[:,u]=PEsameW[:,u] + r[:R2][temp:(temp+250+a)]

    numPEsame=numPEsame+1
elseif r2[:NotError][i] == 1 && r2[:Distri][i] <0 && r2[:Distri][i+1] <0
    #rlist=rlist + r[:R2][temp:(temp+250)]
    PCsameL[:,u]=PCsameL[:,u] + r[:R1][temp:(temp+250+a)]
    PCsameW[:,u]=PCsameW[:,u] + r[:R2][temp:(temp+250+a)]

    numPCsame=numPCsame+1









elseif r2[:NotError][i] == 1 && r2[:Distri][i] >0 && r2[:Distri][i+1] <0 && r2[:NotError][i+1] == 1
    PCdiffL[:,u]=PCdiffL[:,u] + r[:R2][temp:(temp+250+a)]
    PCdiffW[:,u]=PCdiffW[:,u] + r[:R1][temp:(temp+250+a)]

    numPCdiff=numPCdiff+1
    #rlist2=rlist2 + r[:R2][temp:(temp+250)]

elseif r2[:NotError][i] == -1 && r2[:Distri][i] <0 && r2[:Distri][i+1] <0 && r2[:NotError][i+1] == 1
        #rlist=rlist + r[:R2][temp:(temp+250)]
        PEdiffL[:,u]=PEdiffL[:,u] + r[:R2][temp:(temp+250+a)]
        PEdiffW[:,u]=PEdiffW[:,u] + r[:R1][temp:(temp+250+a)]

        numPEdiff=numPEdiff+1

    elseif r2[:NotError][i] == -1 && r2[:Distri][i] >0 && r2[:Distri][i+1] >0 && r2[:NotError][i+1] == 1
        #rlist=rlist + r[:R2][temp:(temp+250)]
        PEdiffL[:,u]=PEdiffL[:,u] + r[:R1][temp:(temp+250+a)]
        PEdiffW[:,u]=PEdiffW[:,u] + r[:R2][temp:(temp+250+a)]

        numPEdiff=numPEdiff+1
    elseif r2[:NotError][i] == 1 && r2[:Distri][i] <0 && r2[:Distri][i+1] >0 && r2[:NotError][i+1] == 1
        #rlist=rlist + r[:R2][temp:(temp+250)]
        PCdiffL[:,u]=PCdiffL[:,u] + r[:R1][temp:(temp+250+a)]
        PCdiffW[:,u]=PCdiffW[:,u] + r[:R2][temp:(temp+250+a)]

        numPCdiff=numPCdiff+1

end
temp = temp+250
end


PEsameL[:,u]=PEsameL[:,u]./numPEsame
PEdiffL[:,u]=PEdiffL[:,u]./numPEdiff
PEsameW[:,u]=PEsameW[:,u]./numPEsame
PEdiffW[:,u]=PEdiffW[:,u]./numPEdiff
PCsameL[:,u]=PCsameL[:,u]./numPCsame
PCdiffL[:,u]=PCdiffL[:,u]./numPCdiff
PCsameW[:,u]=PCsameW[:,u]./numPCsame
PCdiffW[:,u]=PCdiffW[:,u]./numPCdiff

end

#############################"
##########" Same with S
#################


a=300
SPEsameL=zeros(251+a,25)#
SPEdiffL=zeros(251+a,25)
SPEsameW=zeros(251+a,25)
SPEdiffW=zeros(251+a,25)
SPCsameL=zeros(251+a,25)#
SPCdiffL=zeros(251+a,25)
SPCsameW=zeros(251+a,25)
SPCdiffW=zeros(251+a,25)

j=2
for u=2:25
numPEsame=0
numPEdiff=0
numPEsame=0
numPEdiff=0
numPCsame=0
numPCdiff=0
numPCsame=0
numPCdiff=0
Ir=Irlist[1]
c=Cohlist[2]

numr=0
numr2=1


    r=CSV.read("/home/kevin/Documents/1 - PhD Projects/1 - Sequential Decision making/3-SequentialEffects/2-Data/FiringRates/datatempAll/V$j-Coh$c-Final200Exp$u-T$Tr-Ir$Ir.csv")
    r2=CSV.read("/home/kevin/Documents/1 - PhD Projects/1 - Sequential Decision making/3-SequentialEffects/2-Data/FiringRates/datatempAll/V$j-NondynCoh$c-Final200Exp$u-T$Tr-Ir$Ir.csv")

temp = 250
for i=1:2990
temp= temp + Int(r2[:RTs][i]/2)
if r2[:NotError][i] == 1 && r2[:Distri][i] >0 && r2[:Distri][i+1] >0
SPCsameL[:,u]=SPCsameL[:,u] + r[:S2][temp:(temp+250+a)]
SPCsameW[:,u]=SPCsameW[:,u] + r[:S1][temp:(temp+250+a)]

numPCsame=numPCsame+1
#rlist2=rlist2 + r[:R2][temp:(temp+250)]

elseif r2[:NotError][i] == -1 && r2[:Distri][i] <0 && r2[:Distri][i+1] >0
    #rlist=rlist + r[:R2][temp:(temp+250)]
    SPEsameL[:,u]=SPEsameL[:,u] + r[:S2][temp:(temp+250+a)]
    SPEsameW[:,u]=SPEsameW[:,u] + r[:S1][temp:(temp+250+a)]

    numPEsame=numPEsame+1

elseif r2[:NotError][i] == -1 && r2[:Distri][i] >0 && r2[:Distri][i+1] <0
    #rlist=rlist + r[:R2][temp:(temp+250)]
    SPEsameL[:,u]=SPEsameL[:,u] + r[:S1][temp:(temp+250+a)]
    SPEsameW[:,u]=SPEsameW[:,u] + r[:S2][temp:(temp+250+a)]

    numPEsame=numPEsame+1
elseif r2[:NotError][i] == 1 && r2[:Distri][i] <0 && r2[:Distri][i+1] <0
    #rlist=rlist + r[:R2][temp:(temp+250)]
    SPCsameL[:,u]=SPCsameL[:,u] + r[:S1][temp:(temp+250+a)]
    SPCsameW[:,u]=SPCsameW[:,u] + r[:S2][temp:(temp+250+a)]

    numPCsame=numPCsame+1









elseif r2[:NotError][i] == 1 && r2[:Distri][i] >0 && r2[:Distri][i+1] <0 && r2[:NotError][i+1] == 1
    SPCdiffL[:,u]=SPCdiffL[:,u] + r[:S2][temp:(temp+250+a)]
    SPCdiffW[:,u]=SPCdiffW[:,u] + r[:S1][temp:(temp+250+a)]

    numPCdiff=numPCdiff+1
    #rlist2=rlist2 + r[:R2][temp:(temp+250)]

elseif r2[:NotError][i] == -1 && r2[:Distri][i] <0 && r2[:Distri][i+1] <0 && r2[:NotError][i+1] == 1
        #rlist=rlist + r[:R2][temp:(temp+250)]
        SPEdiffL[:,u]=SPEdiffL[:,u] + r[:S2][temp:(temp+250+a)]
        SPEdiffW[:,u]=SPEdiffW[:,u] + r[:S1][temp:(temp+250+a)]

        numPEdiff=numPEdiff+1

    elseif r2[:NotError][i] == -1 && r2[:Distri][i] >0 && r2[:Distri][i+1] >0 && r2[:NotError][i+1] == 1
        #rlist=rlist + r[:R2][temp:(temp+250)]
        SPEdiffL[:,u]=SPEdiffL[:,u] + r[:S1][temp:(temp+250+a)]
        SPEdiffW[:,u]=SPEdiffW[:,u] + r[:S2][temp:(temp+250+a)]

        numPEdiff=numPEdiff+1
    elseif r2[:NotError][i] == 1 && r2[:Distri][i] <0 && r2[:Distri][i+1] >0 && r2[:NotError][i+1] == 1
        #rlist=rlist + r[:R2][temp:(temp+250)]
        SPCdiffL[:,u]=SPCdiffL[:,u] + r[:S1][temp:(temp+250+a)]
        SPCdiffW[:,u]=SPCdiffW[:,u] + r[:S2][temp:(temp+250+a)]

        numPCdiff=numPCdiff+1

end
temp = temp+250
end


SPEsameL[:,u]=SPEsameL[:,u]./numPEsame
SPEdiffL[:,u]=SPEdiffL[:,u]./numPEdiff
SPEsameW[:,u]=SPEsameW[:,u]./numPEsame
SPEdiffW[:,u]=SPEdiffW[:,u]./numPEdiff
SPCsameL[:,u]=SPCsameL[:,u]./numPCsame
SPCdiffL[:,u]=SPCdiffL[:,u]./numPCdiff
SPCsameW[:,u]=SPCsameW[:,u]./numPCsame
SPCdiffW[:,u]=SPCdiffW[:,u]./numPCdiff

end

using Bootstrap
PEL = Float64[]
PELB=Float64[]
PELU = Float64[]

PCL = Float64[]
PCLB=Float64[]
PCLU = Float64[]

PEW = Float64[]
PEWB=Float64[]
PEWU = Float64[]

PCW = Float64[]
PCWB=Float64[]
PCWU = Float64[]


for i=5:300
    n_boot = 2000

      ## basic bootstrap
      #bs1 = bootstrap(PEsameL[i,2:end], mean, BasicSampling(n_boot))
      bs1 = bootstrap(PEdiffL[i,2:end], mean, BasicSampling(n_boot))

      cil = 0.95;

      ## basic CI
      bci1 = ci(bs1, BasicConfInt(cil));
      push!(PEL,bci1[1][1])
      push!(PELB,-bci1[1][3]+bci1[1][1])
      push!(PELU,bci1[1][2]-bci1[1][1])

     # bs1 = bootstrap(PCsameL[i,2:end], mean, BasicSampling(n_boot))
      bs1 = bootstrap(PCdiffL[i,2:end], mean, BasicSampling(n_boot))
      cil = 0.95;

      ## basic CI
      bci1 = ci(bs1, BasicConfInt(cil));
      push!(PCL,bci1[1][1])
      push!(PCLB,-bci1[1][3]+bci1[1][1])
      push!(PCLU,bci1[1][2]-bci1[1][1])


     # bs1 = bootstrap(PEsameW[i,2:end], mean, BasicSampling(n_boot))
     bs1 = bootstrap(PEdiffW[i,2:end], mean, BasicSampling(n_boot))

      cil = 0.95;

      ## basic CI
      bci1 = ci(bs1, BasicConfInt(cil));
      push!(PEW,bci1[1][1])
      push!(PEWB,-bci1[1][3]+bci1[1][1])
      push!(PEWU,bci1[1][2]-bci1[1][1])

    #  bs1 = bootstrap(PCsameW[i,2:end], mean, BasicSampling(n_boot))

      bs1 = bootstrap(PCdiffW[i,2:end], mean, BasicSampling(n_boot))
      cil = 0.95;

      ## basic CI
      bci1 = ci(bs1, BasicConfInt(cil));
      push!(PCW,bci1[1][1])
      push!(PCWB,-bci1[1][3]+bci1[1][1])
      push!(PCWU,bci1[1][2]-bci1[1][1])
end

fig=plot()
plot!(PCW,ribbon=(PCWB,PCWU),linewidth=1,color=:blue)
plot!(PEW,ribbon=(PEWB,PEWU),linewidth=1,color=:darkorange)

plot!(PEL,ribbon=(PELB,PELU),linewidth=1,line=:dash,color=:darkorange)
plot!(PCL,ribbon=(PCLB,PCLU),linewidth=1,line=:dash,color=:blue)

plot!(legend=false)


plot!(mean(PEsameL[5:300,2:end],2),line=:dash,ribbon=sqrt.(var(PEsameL[5:300,2:end],2)),color=:darkorange)
plot!(mean(PCsameL[5:300,2:end],2),line=:dash,ribbon=sqrt.(var(PCsameL[5:300,2:end],2)),color=:blue)

plot(linspace(5,300,296),mean(PEsameW[5:300,2:end],2),ribbon=sqrt.(var(PEsameW[5:300,2:end],2)),color=:darkorange,linewidth=3)
plot!(linspace(5,300,296),mean(PCsameW[5:300,2:end],2),ribbon=sqrt.(var(PCsameW[5:300,2:end],2)),color=:blue,linewidth=3)

plot!(legend=false)
savefig("./PESsame.svg")


plot(mean(log.(PEsameW[200:300,2:end]),2),ribbon=sqrt.(var(log.(PEsameW[200:300,2:end]),2)),color=:darkorange,linewidth=3)
plot!(mean(log.(PCsameW[200:300,2:end]),2),ribbon=sqrt.(var(log.(PCsameW[200:300,2:end]),2)),color=:blue,linewidth=3)

plot!(legend=false)

n_boot=1000
MSPE=Float64[]
MSPC=Float64[]
VUSPE=Float64[]
VDSPE=Float64[]
VUSPC=Float64[]
VDSPC=Float64[]
DMSPE=Float64[]
DVUSPE=Float64[]
DVDSPE=Float64[]

using Bootstrap
for i=10:350
    bs1 = bootstrap(SPEsameW[i,2:end], mean, BasicSampling(n_boot))

    cil = 0.95;

    ## basic CI
    bci1 = ci(bs1, BasicConfInt(cil));
    push!(MSPE,bci1[1][1])
    push!(VUSPE,-bci1[1][3]+bci1[1][1])
    push!(VDSPE,bci1[1][2]-bci1[1][1])


    bs1 = bootstrap(SPCsameW[i,2:end], mean, BasicSampling(n_boot))

    cil = 0.95;

    ## basic CI
    bci1 = ci(bs1, BasicConfInt(cil));
    push!(MSPC,bci1[1][1])
    push!(VUSPC,-bci1[1][3]+bci1[1][1])
    push!(VDSPC,bci1[1][2]-bci1[1][1])

    bs1 = bootstrap(SPEsameW[i,2:end]-SPCsameW[i,2:end], mean, BasicSampling(n_boot))

    cil = 0.95;

    ## basic CI
    bci1 = ci(bs1, BasicConfInt(cil));
    push!(DMSPE,bci1[1][1])
    push!(DVUSPE,-bci1[1][3]+bci1[1][1])
    push!(DVDSPE,bci1[1][2]-bci1[1][1])



end
MSPE
gr()


p1=plot(linspace(10,300,141),MSPE[1:141],ribbon=(VDSPE[1:141],VUSPE[1:141]),color=:darkorange,linewidth=2)
plot!(p1,linspace(10,300,141),MSPC[1:141],ribbon = (VDSPC[1:141],VUSPC[1:141]),color=:blue,linewidth=2)
ps=twinx(p1)
plot!(ps,linspace(600,700,51),MSPE[291:end],ribbon=(VDSPE[291:end],VUSPE[291:end]),color=:darkorange,linewidth=2)
plot!(ps,linspace(600,700,51),MSPC[291:end],ribbon=(VDSPC[291:end],VUSPC[291:end]),color=:blue,linewidth=2)
plot!(legend=false)
vline!([500],line=:dash,color=:black)
plot!(legend=false)



savefig("./PESsameS.svg")


plot(linspace(10,700,341),DMSPE./MSPE.*100,ribbon=(DVDSPE./MSPE.*100,DVUSPE./MSPE.*100),linewidth=3)

plot!(twinx(),linspace(10,700,341),DMSPE./MSPE+sqrt.(linspace(10,700,341)),ribbon=(DVDSPE./MSPE,DVUSPE./MSPE),linewidth=3)

plot!(legend=false)
savefig("./PESsameDiff.svg")


plot(mean(log.(SPEsameW[200:300,2:end]),2),ribbon=sqrt.(var(log.(SPEsameW[200:300,2:end]),2)),color=:darkorange,linewidth=3)
plot!(mean(log.(SPCsameW[200:300,2:end]),2),ribbon=sqrt.(var(log.(SPCsameW[200:300,2:end]),2)),color=:blue,linewidth=3)

plot!(legend=false)

plot(MSPE,ribbon=(VDSPE,VUSPC))
plot!(MSPC,ribbon = (VDSPC,VUSPC))


n_boot=1000
MSPE=Float64[]
MSPC=Float64[]
DMSPE=Float64[]
DVUSPE=Float64[]
DVDSPE=Float64[]

VUSPE=Float64[]
VDSPE=Float64[]
VUSPC=Float64[]
VDSPC=Float64[]
using Bootstrap
for i=10:350
    bs1 = bootstrap(PEsameW[i,2:end], mean, BasicSampling(n_boot))

    cil = 0.95;

    ## basic CI
    bci1 = ci(bs1, BasicConfInt(cil));
    push!(MSPE,bci1[1][1])
    push!(VUSPE,-bci1[1][3]+bci1[1][1])
    push!(VDSPE,bci1[1][2]-bci1[1][1])


    bs1 = bootstrap(PCsameW[i,2:end], mean, BasicSampling(n_boot))

    cil = 0.95;

    ## basic CI
    bci1 = ci(bs1, BasicConfInt(cil));
    push!(MSPC,bci1[1][1])
    push!(VUSPC,-bci1[1][3]+bci1[1][1])
    push!(VDSPC,bci1[1][2]-bci1[1][1])


    bs1 = bootstrap(PEsameW[i,2:end]-PCsameW[i,2:end], mean, BasicSampling(n_boot))

    cil = 0.95;

    ## basic CI
    bci1 = ci(bs1, BasicConfInt(cil));
    push!(DMSPE,bci1[1][1])
    push!(DVUSPE,-bci1[1][3]+bci1[1][1])
    push!(DVDSPE,bci1[1][2]-bci1[1][1])


end


p1=plot(linspace(10,300,141),MSPE[1:141],ribbon=(VDSPE[1:141],VUSPE[1:141]),color=:darkorange,linewidth=2)
plot!(p1,linspace(10,300,141),MSPC[1:141],ribbon = (VDSPC[1:141],VUSPC[1:141]),color=:blue,linewidth=2)
ps=twinx(p1)
plot!(ps,linspace(420,700,141),MSPE[201:end],ribbon=(VDSPE[201:end],VUSPE[201:end]),color=:darkorange,linewidth=2)
plot!(ps,linspace(420,700,141),MSPC[201:end],ribbon=(VDSPC[201:end],VUSPC[201:end]),color=:blue,linewidth=2)
vline!([500],line=:dash,color=:black)
plot!(legend=false)
savefig("./PESBootstrapsame.svg")


plot(linspace(10,700,341),MSPE,ribbon=(VDSPE,VUSPC),linewidth=3,color=:darkorange)
plot!(linspace(10,700,341),MSPC,ribbon = (VDSPC,VUSPC),linewidth=3,color=:blue)
plot!(legend=false)
savefig("./PESBootstrapsame.svg")


plot(linspace(10,700,341),DMSPE./MSPE.*100,ribbon=(DVDSPE./MSPE.*100,DVUSPE./MSPE.*100),linewidth=3)

plot(linspace(10,700,341),DMSPE,ribbon=(DVDSPE,DVUSPE),linewidth=3)
plot!(legend=false)
savefig("./DiffPESsame.svg")




savefig(fig,"./PEQdiffdec.png")



plot(mean(PEdiffL[5:300,2:end],2),ribbon=sqrt.(var(PEdiffL[5:300,2:end],2)))
plot!(mean(PCdiffL[5:300,2:end],2),ribbon=sqrt.(var(PCdiffL[5:300,2:end],2)))

plot!(legend=false)

plot!(mean(PEdiffW[5:300,2:end],2),ribbon=sqrt.(var(PEdiffW[5:300,2:end],2)))
plot!(mean(PCdiffW[5:300,2:end],2),ribbon=sqrt.(var(PCdiffW[5:300,2:end],2)))

plot!(legend=false)


mean(PEsameL[:,2:3],2)
var(PEsameL[:,2:3],2)


#### Indiv stuffs
i=1
r=CSV.read("/home/kevin/Documents/1 - PhD Projects/1 - Sequential Decision making/3-SequentialEffects/2-Data/FiringRates/datatempAll/Coh$c-Final200Exp$i-T$Tr-Ir$Ir.csv")
r2=CSV.read("/home/kevin/Documents/1 - PhD Projects/1 - Sequential Decision making/3-SequentialEffects/2-Data/FiringRates/datatempAll/NondynCoh$c-Final200Exp$i-T$Tr-Ir$Ir.csv")


for i=1:2300
if r2[:NotError][i]==-1
    print(i)
end
end

r2[:NotError][35]
r2[:Distri][35]
r2[:Distri][36]

r2[:NotError][270]
r2[:Distri][144]
r2[:Distri][145]
tinit = sum(r2[:RTs][1:357]./2)+250*357
fig=plot()
#plot(r[:R1][Int(tinit)-100:Int(tinit)+450],linewidth=2,color=:darkorange)
plot!(linspace(0,1000,501),r[:R2][Int(tinit)-50:Int(tinit)+450],linewidth=2,color=:darkorange)
r2[:NotError][103]
r2[:Distri][108]
r2[:Distri][119]
tinit = sum(r2[:RTs][1:118]./2)+250*118
plot!(linspace(0,1000,501),r[:R1][Int(tinit)-50:Int(tinit)+450],linewidth=2,color=:blue)
#plot!(r[:R2][Int(tinit)-100:Int(tinit)+450],linewidth=2,color=:blue)
plot!(legend=false)
vline!([100],line=:dash,linewidth=2,color=:black)
vline!([600],line=:dash,linewidth=2,color=:black)

fig
savefig(fig,"./PEQdiffFR.svg")



r2[:NotError][521]
r2[:Distri][521]
r2[:Distri][522]
tinit = sum(r2[:RTs][1:521]./2)+250*521
fig=plot()
#plot(r[:R1][Int(tinit)-100:Int(tinit)+450],linewidth=2,color=:darkorange)
plot!(linspace(0,900,451),r[:R2][Int(tinit)-50:Int(tinit)+400],linewidth=2,color=:darkorange)
r2[:NotError][103]
r2[:Distri][108]
r2[:Distri][111]
tinit = sum(r2[:RTs][1:111]./2)+250*111
plot!(linspace(0,900,451),r[:R2][Int(tinit)-50:Int(tinit)+400],linewidth=2,color=:blue)
#plot!(r[:R2][Int(tinit)-100:Int(tinit)+450],linewidth=2,color=:blue)
plot!(legend=false)
vline!([100],line=:dash,linewidth=2,color=:black)
vline!([600],line=:dash,linewidth=2,color=:black)

fig
savefig(fig,"./neeffFR.svg")




#### Indiv stuffs
i=1
c=Cohlist[2]
Ir=Irlist[1]
r=CSV.read("/home/kevin/Documents/1 - PhD Projects/1 - Sequential Decision making/3-SequentialEffects/2-Data/FiringRates/datatempAll/Coh$c-Final200Exp$i-T$Tr-Ir$Ir.csv")
r2=CSV.read("/home/kevin/Documents/1 - PhD Projects/1 - Sequential Decision making/3-SequentialEffects/2-Data/FiringRates/datatempAll/NondynCoh$c-Final200Exp$i-T$Tr-Ir$Ir.csv")


for i=1:2300
if r2[:NotError][i]==-1
    print(i)
end
end

r2[:NotError][53]
r2[:Distri][52]
r2[:Distri][53]

r2[:NotError][270]
r2[:Distri][144]
r2[:Distri][145]
fig=plot()
tinit = sum(r2[:RTs][1:499]./2)+250*499
#plot(r[:R1][Int(tinit)-100:Int(tinit)+450],linewidth=2,color=:darkorange)
#linspace(0,1000,501),
plot!(linspace(0,1000,501),r[:R2][Int(tinit)-50:Int(tinit)+450],linewidth=2,color=:darkorange)
r2[:NotError][161]
r2[:Distri][160]
r2[:Distri][161]
tinit = sum(r2[:RTs][1:365]./2)+250*365
plot!(linspace(0,1000,501),r[:R2][Int(tinit)-50:Int(tinit)+450],linewidth=2,color=:blue)
#plot!(r[:R2][Int(tinit)-100:Int(tinit)+450],linewidth=2,color=:blue)
plot!(legend=false)
vline!([100],line=:dash,linewidth=2,color=:black)
vline!([600],line=:dash,linewidth=2,color=:black)
savefig(fig,"./PESsameFR.svg")


for i=1:length(r2[:NotError])-1
    if r2[:NotError][i]==1 && r2[:NotError][i+1]==1 &&
         sign(r2[:Distri][i])==sign(r2[:Distri][i+1])
print(i)
print("\n")
end
end


for i=1:length(r2[:NotError])-1
    if r2[:NotError][i]==-1 && r2[:NotError][i+1]==1 &&
         sign(r2[:Distri][i])!=sign(r2[:Distri][i+1])
print(i)
print("\n")
end
end


fig



savefig(fig,"./PESdiffFR.svg")
