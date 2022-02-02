# Implementation of a decision making model

# Ref :Wang 2006 : A recurrent network mechanisms of time integrtion in perceptual decision
#

# Some part of this implementation are extended from the original model.
# The idea is to model decision making with two competitve population
# in Mean field framework
using Distributions
using DataFrames

function sim_mf_decisionmaking(N_trials,param,S1,S2,s_save=false)

  # Args :
  #  N_trials : TOtal number of trials
  # Tstim : total duration stimulus
  # Tstart : start of the stimulus
  # mu0 : external stimulus strength
  # noise_amp : Noise amplitude into selective populations
  # coh : coherence of the stiimulus

  # Return a list of s1,s2,r1,r2 for the dynamical system
  # Return a DataFrame object

  print("Start of the simulation")

  # some constants of the article
  Tnmda = 100
  Tampa = 2
  gamma = 0.641
	th = param["Threshold"]
  # FI curve parameters
  a = 270
  b=108
  d=0.1540

  s_list1=[]
  s_list2=[]

  r1_traj=[]
  r2_traj=[]
  s1_traj=[]
  s2_traj=[]

  nums1=[]
  rt_list=Float64[]

  for n=1:N_trials
    trials_no = n

    nu1_wind = []
    nu2_wind=[]
    s1_wind=[]
    s2_wind=[]

    s1_in = S1
    s2_in=S2

    nu1_in = 2
    nu2_in=2

    I_eta1_in = param["noise_amp"] * randn()
    I_eta2_in = param["noise_amp"] * randn()

    dt=0.5
    time_wind=2/dt # mean window of the variables
    T_total = (time_wind*dt+param["Tstim"])/dt
    slide_wind = time_wind

    #---- Intialise and vectorise variables to be used in loops below ------

    s1 = s1_in*ones(T_total); s2 = s2_in*ones(T_total);
    nu1 = nu1_in*ones(T_total); nu2 = nu2_in*ones(T_total);
    phi1 = nu1_in*ones(T_total); phi2 = nu2_in*ones(T_total);
    I_eta1 = I_eta1_in*ones(T_total); I_eta2 = I_eta2_in*ones(T_total);
    Isyn1=zeros(T_total);
    Isyn2=zeros(T_total);

    for t=1:Int64(T_total)-1

      if 0<t && t<(param["Tstim"]+param["T_rest"])/dt
        I_stim1 = param["JAext"]*param["mu0"]*(1+param["coh"]/100)

      else
        I_stim1 = 0

      end

      if 0<t && t<(param["Tstim"]+param["T_rest"])/dt
        I_stim2 = param["JAext"]*param["mu0"]*(1-param["coh"]/100)

      else
        I_stim2 = 0

      end

      Isyn1[t]=param["JN11"]*s1[t]-param["JN12"]*s2[t]+I_stim1 + I_eta1[t]
      Isyn2[t]=param["JN22"]*s2[t]-param["JN21"]*s1[t]+I_stim2 + I_eta2[t]

      phi1[t]= (a*Isyn1[t]-b)/(1-exp(-d*(a*Isyn1[t]-b)))
      phi2[t]= (a*Isyn2[t]-b)/(1-exp(-d*(a*Isyn2[t]-b)))

      #---- Dynamical equations -------------------------------------------

      # Mean NMDA-mediated synaptic dynamics updating
      s1[t+1] = s1[t] + dt*(-s1[t]/Tnmda + (1-s1[t])*gamma*nu1[t]/1000); # Because we are in mS and Hz simultaneously
      s2[t+1] = s2[t] + dt*(-s2[t]/Tnmda + (1-s2[t])*gamma*nu2[t]/1000); # Because we are in mS and Hz simultaneously

      # Ornstein-Uhlenbeck generation of noise in pop1 and 2
      I_eta1[t+1] = I_eta1[t] + (dt/Tampa)*(param["I0E1"]-I_eta1[t]) + sqrt(dt/Tampa)*param["noise_amp"]*randn() ;
      I_eta2[t+1] = I_eta2[t] + (dt/Tampa)*(param["I0E2"]-I_eta2[t]) + sqrt(dt/Tampa)*param["noise_amp"]*randn() ;

      # To ensure firing rates are always positive. Large noise amplitude
      # may result in unwanted negative values
      if phi1[t] < 0
        nu1[t+1] = 0;
        phi1[t] = 0;
      else
        nu1[t+1] = phi1[t];
      end;
      if phi2[t] < 0
        nu2[t+1] = 0;
        phi2[t] = 0;
      else
        nu2[t+1] = phi2[t];
      end;

    end

    #---- Calculating the mean rates and gating variables with sliding window -----
    mean(nu1[1:Int(time_wind)])
    push!(nu1_wind,mean(nu1[1:Int(time_wind)]))
    push!(nu2_wind,mean(nu2[1:Int(time_wind)]))
    push!(s1_wind,mean(s1[1:Int(time_wind)]))
    push!(s2_wind,mean(s2[1:Int(time_wind)]))

    for t = 1:Int64((T_total-time_wind)/slide_wind)-1
      push!(nu1_wind,mean(nu1[Int(slide_wind*t):Int(slide_wind*t+time_wind)]))
      push!(nu2_wind,mean(nu2[Int(slide_wind*t):Int(slide_wind*t+time_wind)]))
      push!(s1_wind,mean(s1[Int(slide_wind*t):Int(slide_wind*t+time_wind)]))
      push!(s2_wind,mean(s2[Int(slide_wind*t):Int(slide_wind*t+time_wind)]))
    end

    push!(r1_traj,nu1_wind)
    push!(r2_traj,nu2_wind)

    push!(s1_traj,s1_wind)
    push!(s2_traj,s2_wind)

    r1=0
    r2=0
    t=0
    while r1 <th && r2<th && t<2750
      t=t+1

      r1=nu1_wind[t]
      r2=nu2_wind[t]
    end
    if r1>th && param["coh"] >=0
      push!(nums1,1)
    elseif r1>th && param["coh"]<0

      push!(nums1,-1)
    elseif r2>th && param["coh"]<0
      push!(nums1,1)
    else
      push!(nums1,-1)
    end
    push!(rt_list,(t*slide_wind+time_wind/2)*dt-param["T_rest"])

  end


  if s_save==true
    result = Dict("s1"=>s1_traj,
    "s2"=>s2_traj,
    "r1"=>r1_traj,
    "r2"=> r2_traj,
    "NotError"=>nums1,
    "RTs"=>rt_list)
  else
    result = Dict("NotError"=>nums1,
    "RTs"=>rt_list)

  end
  return result
end


function sim_mf_decisionmaking_Huk(N_trials,List_coh,List_rest,List_time1,
  dic_param=("I0E1"=>0.3255,"I0E2"=>0.3255),mu0=30,noise_amp=0.009,
  I0E1=0.3297,I0E2=0.3297,JN11=0.3725,JN22=0.3725,JN12=0.1137,
  JN21=0.1137)
  # Description : Return the trajectory of the dynamical system in the case of
  # decision making with two competitive populations
  # and for a serie of stimulus. Between stimuli the input I0 is cut.
  # I0 correpsonds to selective and non selective background in our model
  # New synaptics couplings : stronger for better interactions between attractor state ???

  # Args:
  # List_rest : list of the time of reset between stimuli
  # List_coh : coherence of the several stimuli
  # List_time : duration of the stimuli
  # Dic_param : general param of the model
  # Dic_spec : specific param to the simulation

  # Output :
  # Reaction time with half direction selectivity
  # Wrong or Right list
  # S1 S2 R1 R2
  Tnoise = 60 # new fit
  Tampa = 2
  gamma = 0.641
  # FI curve parameters : never change
  a = 270
  b=108
  d=0.1540

  JAext = 0.0011

  r1_traj=[]
  r2_traj=[]
  s1_traj=[]
  s2_traj=[]
  t_tot=[]
  I_list=[]

  dt=0.5
  time_wind=2/dt # mean window of the variables
  T_total = (time_wind*dt+sum(List_time1)+sum(List_rest))/dt
  slide_wind = time_wind

  List_time = ones(length(List_time1)+length(List_rest))
  for i =1:length(List_time1)
    List_time[2*i]=List_time1[i]
  end
  for i =1:length(List_rest)
    List_time[2*i-1]=List_rest[i]
  end

  for n=1:N_trials
    println(n)
    trials_no = n

    nu1_wind = []
    nu2_wind=[]
    s1_wind=[]
    s2_wind=[]
    t_list=[]
    I_stim=[]
    s1_in = 0.1
    s2_in=0.1

    nu1_in = 2
    nu2_in=2

    I_eta1_in = noise_amp * randn()
    I_eta2_in = noise_amp * randn()


    #---- Intialise and vectorise variables to be used in loops below ------
    s1 = s1_in*ones(T_total); s2 = s2_in*ones(T_total);
    nu1 = nu1_in*ones(T_total); nu2 = nu2_in*ones(T_total);
    phi1 = nu1_in*ones(T_total); phi2 = nu2_in*ones(T_total);
    I_eta1 = I_eta1_in*ones(T_total); I_eta2 = I_eta2_in*ones(T_total);
    Isyn1=zeros(T_total);
    Isyn2=zeros(T_total);
    I_stimtemp=zeros(T_total);

    for t=1:Int64(T_total)-1

      indice = 1
      if  t<List_time[1]/dt
        I_stim1 = 0
        I_stim2=0
        I0E1 = 0.3255 # 0 current external
        I0E2 = 0.3255
      else
        for i=2:length(List_time)
          indice = indice + 1
          if  t<sum(List_time[1:i])/dt&& t>=sum(List_time[1:(i-1)])/dt
            break
          end

        end

      end

      # search the coh value
      if mod(indice,2)==1
        # no stimuli
        I_stim1=0
        I_stim2=0
        I0E1=0.30
        I0E2 = 0.30
        if  t<List_time[1]/dt
          I_stim1 = 0
          I_stim2=0
          I0E1 = 0.3255 # 0 current external
          I0E2 = 0.3255
        end
        I_stimtemp[t]=0
      else
        coh = List_coh[Int64(indice/2)]
        I_stim1 = JAext*mu0*(1+0.45*coh/100)
        I_stim2 = JAext*mu0*(1-0.45*coh/100)
        I0E1=0.3255#dic_param["I0E1"]
        I0E2=0.3255#dic_param["I0E2"]
        I_stimtemp[t]=0.5*sign(coh)
      end


      Isyn1[t]=JN11*s1[t]-JN12*s2[t]+I_stim1 + I_eta1[t]
      Isyn2[t]=JN22*s2[t]-JN21*s1[t]+I_stim2 + I_eta2[t]

      phi1[t]= (a*Isyn1[t]-b)/(1-exp(-d*(a*Isyn1[t]-b)))
      phi2[t]= (a*Isyn2[t]-b)/(1-exp(-d*(a*Isyn2[t]-b)))

      #---- Dynamical equations -------------------------------------------

      # Mean NMDA-mediated synaptic dynamics updating
      s1[t+1] = s1[t] + dt*(-s1[t]/Tnmda + (1-s1[t])*gamma*nu1[t]/1000);
      s2[t+1] = s2[t] + dt*(-s2[t]/Tnmda + (1-s2[t])*gamma*nu2[t]/1000);

      # Ornstein-Uhlenbeck generation of noise in pop1 and 2
      I_eta1[t+1] = I_eta1[t] + (dt/Tampa)*(I0E1-I_eta1[t]) + sqrt(dt/Tampa)*noise_amp*randn() ;
      I_eta2[t+1] = I_eta2[t] + (dt/Tampa)*(I0E2-I_eta2[t]) + sqrt(dt/Tampa)*noise_amp*randn() ;

      # To ensure firing rates are always positive. Large noise amplitude
      # may result in unwanted negative values
      if phi1[t] < 0
        nu1[t+1] = 0;
        phi1[t] = 0;
      else
        nu1[t+1] = phi1[t];
      end;
      if phi2[t] < 0
        nu2[t+1] = 0;
        phi2[t] = 0;
      else
        nu2[t+1] = phi2[t];
      end;

    end

    #---- Calculating the mean rates and gating variables with sliding window -----
    mean(nu1[1:Int(time_wind)])
    push!(nu1_wind,mean(nu1[1:Int(time_wind)]))
    push!(nu2_wind,mean(nu2[1:Int(time_wind)]))
    push!(s1_wind,mean(s1[1:Int(time_wind)]))
    push!(s2_wind,mean(s2[1:Int(time_wind)]))
    push!(t_list,time_wind/2*dt)
    push!(I_stim,mean(I_stimtemp[1:Int(time_wind)]))

    for t = 1:Int64((T_total-time_wind)/slide_wind)-1
      push!(nu1_wind,mean(nu1[Int(slide_wind*t):Int(slide_wind*t+time_wind)]))
      push!(nu2_wind,mean(nu2[Int(slide_wind*t):Int(slide_wind*t+time_wind)]))
      push!(s1_wind,mean(s1[Int(slide_wind*t):Int(slide_wind*t+time_wind)]))
      push!(s2_wind,mean(s2[Int(slide_wind*t):Int(slide_wind*t+time_wind)]))
      push!(t_list,(slide_wind*t+time_wind/2)*dt)
      push!(I_stim,mean(I_stimtemp[1:Int(slide_wind*t):Int(slide_wind*t+time_wind)]))
    end

    push!(r1_traj,nu1_wind)
    push!(r2_traj,nu2_wind)

    push!(s1_traj,s1_wind)
    push!(s2_traj,s2_wind)
    if n==1
      push!(t_tot,t_list)
    end
    if n==1
      push!(I_list,I_stim)
    end
  end

  return s1_traj,s2_traj,r1_traj,r2_traj,t_tot,I_list
end


function sim_mf_decisionmaking_pattern_fast(coh_temp,Ir,IF,
  param,fixed_param,simu_param)
  # Description : Return the trajectory of the dynamical system in the case of
  # decision making with two competitive populations
  # and for a serie of stimulus.


  # Args:
  # List_rest : list of the time of reset between stimuli
  # List_coh : coherence of the several stimuli
  # List_time : duration of the stimuli
  #time : dictionnary with all the time specificity of our model
  #simu_param : Dictionary of the specifity of simulation
  # param : general param of the model
  # Dic_spec : specific param to the simulation

  # Output :
  # Reaction time with half direction selectivity
  # Return 1 or 2 r1 r2
  # Save param et time
  # return timing spike (for example fixed at the value in Wang)
  # return the true list of decision, and maybe rest distribution ????
  # Wrong or Right list
  # S1 S2 R1 R2


  # FI curve parameters - adapted from Wanf fit
  a=fixed_param["a"]
  b=fixed_param["b"]
  d=fixed_param["d"]
  th = 20


  # time variables of the  simulation

  dt=simu_param["dt"]
  time_wind=simu_param["time_wind"] # mean window of the variables

  slide_wind = simu_param["slide_wind"]

  false_time=param["false_time"]

  deltaSstim=[]

  #---- Intialise and vectorise variables to be used in loops below ------

  # Final Variables

  res=[]

  reaction_time=[]
  Istim1 = []
  # Temporary variables

  Isyn1=[0.0]
  Isyn2=[0.0]
  Trest=param["T_rest"]
  s1=[0.1]
  s2=[0.1]
  phi1=[0.0]
  phi2=[0.0]
  I_eta1=[param["noise_amp"] * randn()]
  I_eta2=[param["noise_amp"] * randn()]
  nu1=[2.0]
  nu2=[2.0]
  I=[0.0]
  diffS=[]
S1=[]
S2=[]
R1=[]
R2=[]
  k=0 # variable for mod 4 in the mean
  th=false
  result=true
  for j=1:length(coh_temp) # loop on trials number
    coh=coh_temp[j]
    t=0

    # reinit
    Isyn1=[Isyn1[end]]
    Isyn2=[Isyn2[end]]
    s1=[s1[end]]
    s2=[s2[end]]
    phi1=[phi1[end]]
    phi2=[phi2[end]]
    I_eta1=[I_eta1[end]]
    I_eta2=[I_eta2[end]]
    nu1=[nu1[end]]
    nu2=[nu2[end]]
    k=0
    while th==false # while we do not cross the threshold
      t=t+dt
      k=k+1
      if t<Trest && result==true
        I0E1=0.3255 -Ir*exp(-t/false_time) #100ms time duration
        I0E2=0.3255 -Ir*exp(-t/false_time) #100ms time duration
        #push!(I,I0E1)
        I_stim1=0
        I_stim2=0
        cht=0
      elseif t<Trest
        I0E1=0.3255 -Ir*exp(-t/false_time) -IF #100ms time duration
        I0E2=0.3255 -Ir*exp(-t/false_time) -IF#100ms time duration
        #push!(I,I0E1)
        I_stim1=0
        I_stim2=0
        cht=0
      elseif t<Trest+dt/2 && t>Trest-dt/2
        I0E1=0.3255#dic_param["I0E1"]
        #push!(I,I0E1) #assuming input on 1 and 2 are the same
        I0E2=0.3255#dic_param["I0E2"]
        I_stim1 = param["JAext"]*param["mu0"]*(1+coh/100)
        I_stim2 = param["JAext"]*param["mu0"]*(1-coh/100)
        cht=coh
        push!(deltaSstim,s1[end] - s2[end])

      else
        I0E1=0.3255#dic_param["I0E1"]
        #push!(I,I0E1) #assuming input on 1 and 2 are the same
        I0E2=0.3255#dic_param["I0E2"]
        I_stim1 = param["JAext"]*param["mu0"]*(1+coh/100)
        I_stim2 = param["JAext"]*param["mu0"]*(1-coh/100)
        cht=coh

      end


      push!(Isyn1,param["JN11"]*s1[end]-param["JN12"]*s2[end]+I_stim1 + I_eta1[end])
      push!(Isyn2,param["JN22"]*s2[end]-param["JN21"]*s1[end]+I_stim2 + I_eta2[end])

      push!(phi1,(a*Isyn1[end]-b)/(1-exp(-d*(a*Isyn1[end]-b))))
      push!(phi2,(a*Isyn2[end]-b)/(1-exp(-d*(a*Isyn2[end]-b))))
      #---- Dynamical equations -------------------------------------------

      # Mean NMDA-mediated synaptic dynamics updating
      s1t = s1[end] + dt*(-s1[end]/param["Tnmda"] + (1-s1[end])*param["gamma"]*nu1[end]/1000);
      push!(s1,s1t)
      s2t = s2[end] + dt*(-s2[end]/param["Tnmda"] + (1-s2[end])*param["gamma"]*nu2[end]/1000);
      push!(s2,s2t)
      # Ornstein-Uhlenbeck generation of noise in pop1 and 2
      I_eta1t = I_eta1[end] + (dt/param["Tampa"])*(I0E1-I_eta1[end]) + sqrt(dt/param["Tampa"])*param["noise_amp"]*randn() ;
      I_eta2t = I_eta2[end] + (dt/param["Tampa"])*(I0E2-I_eta2[end]) + sqrt(dt/param["Tampa"])*param["noise_amp"]*randn() ;
      push!(I_eta1,I_eta1t)
      push!(I_eta2,I_eta2t)
      # To ensure firing rates are always positive. Large noise amplitude
      # may result in unwanted negative values
      if phi1[end] < 0
        push!(nu1,0)
        phi1[end] = 0;
      else
        push!(nu1,phi1[end])
      end;
      if phi2[end] < 0
        push!(nu2,0)
        phi2[end] = 0;
      else
        push!(nu2,phi2[end])
      end;

      if mod(k,4)==0

        # Test the crossing of a threshold
        if mean(nu1[k-3:k])>20 && t>Trest
          th=true
          if coh>=0
            push!(res,1) # push the result
            result=true
          else
            push!(res,-1)
            result=false
          end

        elseif mean(nu2[k-3:k])[end]>20 && t>Trest
          th = true
          if coh>=0
            push!(res,-1)
            result=false
          else
            push!(res,1)
            result=true
          end

        end
      end
    end

    push!(reaction_time,t-Trest)
    push!(diffS,s2[end]-s1[end])
   # push!(S1,nu1[end])
   #push!(S2,nu2[end])
     push!(S1,s1[end])
    push!(S2,s2[end])
     push!(R1,nu1[end])
    push!(R2,nu2[end])


    th=false
  end
  result=DataFrame()
result[:RTs]=reaction_time
result[:NotError]=res
result[:Distri]=coh_temp
result[:DiffS]=diffS
result[:S1]=S1
result[:S2]=S2

result[:DiffSstim]=deltaSstim
result[:LastTrial]=vcat([0],res[1:end-1])
result[:LastCoh]=vcat([0],coh_temp[1:end-1])

  return result
end

function sim_mf_decisionmaking_pattern(coh_temp,Ir,IF,
  param,fixed_param,simu_param)
  # Description : Return the trajectory of the dynamical system in the case of
  # decision making with two competitive populations
  # and for a serie of stimulus.


  # Args:
  # List_rest : list of the time of reset between stimuli
  # List_coh : coherence of the several stimuli
  # List_time : duration of the stimuli
  #time : dictionnary with all the time specificity of our model
  #simu_param : Dictionary of the specifity of simulation
  # param : general param of the model
  # Dic_spec : specific param to the simulation

  # Output :
  # Reaction time with half direction selectivity
  # Return 1 or 2 r1 r2
  # Save param et time
  # return timing spike (for example fixed at the value in Wang)
  # return the true list of decision, and maybe rest distribution ????
  # Wrong or Right list
  # S1 S2 R1 R2


  # FI curve parameters - adapted from Wanf fit
  a=fixed_param["a"]
  b=fixed_param["b"]
  d=fixed_param["d"]
  th = 20

  false_time=param["false_time"]
  # time variables of the  simulation

  dt=simu_param["dt"]
  time_wind=simu_param["time_wind"] # mean window of the variables

  slide_wind = simu_param["slide_wind"]





  #---- Intialise and vectorise variables to be used in loops below ------

  # Final Variables
  s1_mean=[]
  s2_mean=[]
  r1_mean=[]
  r2_mean=[]
  res=[]

  t_mean=[]
  reaction_time=[]

  # Temporary variables

  Isyn1=[]
  Isyn2=[]
  Trest=param["T_rest"]
  s1=[0.1]
  s2=[0.1]
  phi1=[]
  phi2=[]
  I_eta1=[param["noise_amp"] * randn()]
  I_eta2=[param["noise_amp"] * randn()]
  nu1=[2.0]
  nu2=[2.0]
  I=[]


  k=0 # variable for mod 4 in the mean
  th=false
  result=true
  for j=1:length(coh_temp) # loop on trials number
    coh=coh_temp[j]
    t=0
    while th==false # while we do not cross the threshold
      t=t+dt
      k=k+1
      if t<Trest && result==true
        I0E1=0.3255 -Ir*exp(-t/false_time) #100ms time duration
        I0E2=0.3255 -Ir*exp(-t/false_time) #100ms time duration

        I_stim1=0
        I_stim2=0
        cht=0
      elseif t<Trest
        I0E1=0.3255 -Ir*exp(-t/false_time) -IF #100ms time duration
        I0E2=0.3255 -Ir*exp(-t/false_time) -IF#100ms time duration

        I_stim1=0
        I_stim2=0
        cht=0
      else
        I0E1=0.3255#dic_param["I0E1"]
        push!(I,I0E1) #assuming input on 1 and 2 are the same
        I0E2=0.3255#dic_param["I0E2"]
        I_stim1 = param["JAext"]*param["mu0"]*(1+coh/100)
        I_stim2 = param["JAext"]*param["mu0"]*(1-coh/100)
        cht=coh

      end


      push!(Isyn1,param["JN11"]*s1[end]-param["JN12"]*s2[end]+I_stim1 + I_eta1[end])
      push!(Isyn2,param["JN22"]*s2[end]-param["JN21"]*s1[end]+I_stim2 + I_eta2[end])

      push!(phi1,(a*Isyn1[end]-b)/(1-exp(-d*(a*Isyn1[end]-b))))
      push!(phi2,(a*Isyn2[end]-b)/(1-exp(-d*(a*Isyn2[end]-b))))
      #---- Dynamical equations -------------------------------------------

      # Mean NMDA-mediated synaptic dynamics updating
      s1t = s1[end] + dt*(-s1[end]/param["Tnmda"] + (1-s1[end])*param["gamma"]*nu1[end]/1000);
      push!(s1,s1t)
      s2t = s2[end] + dt*(-s2[end]/param["Tnmda"] + (1-s2[end])*param["gamma"]*nu2[end]/1000);
      push!(s2,s2t)
      # Ornstein-Uhlenbeck generation of noise in pop1 and 2
      I_eta1t = I_eta1[end] + (dt/param["Tampa"])*(I0E1-I_eta1[end]) + sqrt(dt/param["Tampa"])*param["noise_amp"]*randn() ;
      I_eta2t = I_eta2[end] + (dt/param["Tampa"])*(I0E2-I_eta2[end]) + sqrt(dt/param["Tampa"])*param["noise_amp"]*randn() ;
      push!(I_eta1,I_eta1t)
      push!(I_eta2,I_eta2t)
      # To ensure firing rates are always positive. Large noise amplitude
      # may result in unwanted negative values
      if phi1[end] < 0
        push!(nu1,0)
        phi1[end] = 0;
      else
        push!(nu1,phi1[end])
      end;
      if phi2[end] < 0
        push!(nu2,0)
        phi2[end] = 0;
      else
        push!(nu2,phi2[end])
      end;

      if mod(k,4)==0
        push!(s1_mean,mean(s1[k-3:k]))
        push!(s2_mean,mean(s2[k-3:k]))
        push!(r1_mean,mean(nu1[k-3:k]))
        push!(r2_mean,mean(nu2[k-3:k]))

        push!(t_mean,(k/4*slide_wind+time_wind/2)*dt) # new time

        # Test the crossing of a threshold
        if r1_mean[end]>20 && t>Trest#+100
          th=true
          if coh>=0
            push!(res,1) # push the result
            result=true
          else
            push!(res,-1)
            result=false
          end

        elseif r2_mean[end]>20 && t>Trest#+100
          th = true
          if coh>=0
            push!(res,-1)
            result=false
          else
            push!(res,1)
            result=true
          end

        end
      end
    end

    push!(reaction_time,t-Trest)
    th=false
  end
  result=DataFrame()
    result2=DataFrame()
	result2[:RTs]=reaction_time
	result2[:NotError]=res
	result2[:Distri]=coh_temp
	result[:R1]=r1_mean
	result[:S1]=s1_mean
	result[:S2]=s2_mean
	result[:R2]=r2_mean
	result[:Temps]=t_mean
	result2[:LastTrial]=vcat([0],res[1:end-1])


  return result,result2
end


function sim_mf_decisionmaking_withoutnoise(N_trials,param,S0,S1,R1,R2,Ir)

  # Args :
  #  N_trials : TOtal number of trials
  # Tstim : total duration stimulus
  # Tstart : start of the stimulus
  # mu0 : external stimulus strength
  # noise_amp : Noise amplitude into selective populations
  # coh : coherence of the stiimulus

  # Return a list of s1,s2,r1,r2 for the dynamical system
  # Return a DataFrame object

  print("Start of the simulation")

  # some constants of the article
  Tnmda = 100
  Tampa = 2
  gamma = 0.641
	th = param["Threshold"]
  # FI curve parameters
  a = 270
  b=108
  d=0.1540

  s_list1=[]
  s_list2=[]

  r1_traj=[]
  r2_traj=[]
  s1_traj=[]
  s2_traj=[]

  nums1=[]
  rt_list=Float64[]

  for n=1:N_trials
    trials_no = n

    nu1_wind = []
    nu2_wind=[]
    s1_wind=[]
    s2_wind=[]

    s1_in = S0
    s2_in=S1

    nu1_in = R1
    nu2_in=R2

    I_eta1_in = 0
    I_eta2_in = 0

    dt=0.5
    time_wind=2/dt # mean window of the variables
    T_total = (time_wind*dt+param["T_rest"])/dt
    slide_wind = time_wind

    #---- Intialise and vectorise variables to be used in loops below ------

    s1 = s1_in*ones(T_total); s2 = s2_in*ones(T_total);
    nu1 = nu1_in*ones(T_total); nu2 = nu2_in*ones(T_total);
    phi1 = nu1_in*ones(T_total); phi2 = nu2_in*ones(T_total);

    Isyn1=zeros(T_total);
    Isyn2=zeros(T_total);

    for t=1:Int64(T_total)-1


        I_stim1 = 0


        I_stim2 = 0


      Isyn1[t]=param["JN11"]*s1[t]-param["JN12"]*s2[t]+I_stim1 +0.3255-Ir*exp(-t/200)
      Isyn2[t]=param["JN22"]*s2[t]-param["JN21"]*s1[t]+I_stim2 +0.3255-Ir*exp(-t/200)
      phi1[t]= (a*Isyn1[t]-b)/(1-exp(-d*(a*Isyn1[t]-b)))
      phi2[t]= (a*Isyn2[t]-b)/(1-exp(-d*(a*Isyn2[t]-b)))

      #---- Dynamical equations -------------------------------------------

      # Mean NMDA-mediated synaptic dynamics updating
      s1[t+1] = s1[t] + dt*(-s1[t]/Tnmda + (1-s1[t])*gamma*nu1[t]/1000); # Because we are in mS and Hz simultaneously
      s2[t+1] = s2[t] + dt*(-s2[t]/Tnmda + (1-s2[t])*gamma*nu2[t]/1000); # Because we are in mS and Hz simultaneously


      # To ensure firing rates are always positive. Large noise amplitude
      # may result in unwanted negative values
      if phi1[t] < 0
        nu1[t+1] = 0;
        phi1[t] = 0;
      else
        nu1[t+1] = phi1[t];
      end;
      if phi2[t] < 0
        nu2[t+1] = 0;
        phi2[t] = 0;
      else
        nu2[t+1] = phi2[t];
      end;

    end

    #---- Calculating the mean rates and gating variables with sliding window -----
    mean(nu1[1:Int(time_wind)])
    push!(nu1_wind,mean(nu1[1:Int(time_wind)]))
    push!(nu2_wind,mean(nu2[1:Int(time_wind)]))
    push!(s1_wind,mean(s1[1:Int(time_wind)]))
    push!(s2_wind,mean(s2[1:Int(time_wind)]))

    for t = 1:Int64((T_total-time_wind)/slide_wind)-1
      push!(nu1_wind,mean(nu1[Int(slide_wind*t):Int(slide_wind*t+time_wind)]))
      push!(nu2_wind,mean(nu2[Int(slide_wind*t):Int(slide_wind*t+time_wind)]))
      push!(s1_wind,mean(s1[Int(slide_wind*t):Int(slide_wind*t+time_wind)]))
      push!(s2_wind,mean(s2[Int(slide_wind*t):Int(slide_wind*t+time_wind)]))
    end

    push!(r1_traj,nu1_wind)
    push!(r2_traj,nu2_wind)

    push!(s1_traj,s1_wind)
    push!(s2_traj,s2_wind)

    r1=0
    r2=0
    t=0

   # push!(rt_list,(t*slide_wind+time_wind/2)*dt-param["T_rest"])

  end



    result = Dict("s1"=>s1_traj,
    "s2"=>s2_traj,
    "r1"=>r1_traj,
    "r2"=> r2_traj
    )

  return result
end


function sim_mf_decisionmaking_relaxation(N_trials,param,S0,S1,R1,R2,Ir)

  # Args :
  #  N_trials : TOtal number of trials
  # Tstim : total duration stimulus
  # Tstart : start of the stimulus
  # mu0 : external stimulus strength
  # noise_amp : Noise amplitude into selective populations
  # coh : coherence of the stiimulus

  # Return a list of s1,s2,r1,r2 for the dynamical system
  # Return a DataFrame object

  print("Start of the simulation")

  # some constants of the article
  Tnmda = 100
  Tampa = 2
  gamma = 0.641
	th = param["Threshold"]
  # FI curve parameters
  a = 270
  b=108
  d=0.1540

  s_list1=[]
  s_list2=[]

  r1_traj=[]
  r2_traj=[]
  s1_traj=[]
  s2_traj=[]

  nums1=[]
  rt_list=Float64[]

  for n=1:N_trials
    trials_no = n

    nu1_wind = []
    nu2_wind=[]
    s1_wind=[]
    s2_wind=[]

    s1_in = S0
    s2_in=S1

    nu1_in = R1
    nu2_in=R2

    I_eta1_in = 0
    I_eta2_in = 0

    dt=0.5
    time_wind=2/dt # mean window of the variables
    T_total = (time_wind*dt+param["T_rest"])/dt
    slide_wind = time_wind

    #---- Intialise and vectorise variables to be used in loops below ------

    s1 = s1_in*ones(T_total); s2 = s2_in*ones(T_total);
    nu1 = nu1_in*ones(T_total); nu2 = nu2_in*ones(T_total);
    phi1 = nu1_in*ones(T_total); phi2 = nu2_in*ones(T_total);

    Isyn1=zeros(T_total);
    Isyn2=zeros(T_total);

    for t=1:Int64(T_total)-1


        I_stim1 = 0


        I_stim2 = 0


      Isyn1[t]=param["JN11"]*s1[t]-param["JN12"]*s2[t]+I_stim1 +0.3255-Ir*exp(-t/200)
      Isyn2[t]=param["JN22"]*s2[t]-param["JN21"]*s1[t]+I_stim2 +0.3255-Ir*exp(-t/200)
      phi1[t]= (a*Isyn1[t]-b)/(1-exp(-d*(a*Isyn1[t]-b)))
      phi2[t]= (a*Isyn2[t]-b)/(1-exp(-d*(a*Isyn2[t]-b)))

      #---- Dynamical equations -------------------------------------------

      # Mean NMDA-mediated synaptic dynamics updating
      s1[t+1] = s1[t] + dt*(-s1[t]/Tnmda + (1-s1[t])*gamma*nu1[t]/1000); # Because we are in mS and Hz simultaneously
      s2[t+1] = s2[t] + dt*(-s2[t]/Tnmda + (1-s2[t])*gamma*nu2[t]/1000); # Because we are in mS and Hz simultaneously


      # To ensure firing rates are always positive. Large noise amplitude
      # may result in unwanted negative values
      if phi1[t] < 0
        nu1[t+1] = 0;
        phi1[t] = 0;
      else
        nu1[t+1] = phi1[t];
      end;
      if phi2[t] < 0
        nu2[t+1] = 0;
        phi2[t] = 0;
      else
        nu2[t+1] = phi2[t];
      end;

    end

    #---- Calculating the mean rates and gating variables with sliding window -----
    mean(nu1[1:Int(time_wind)])
    push!(nu1_wind,mean(nu1[1:Int(time_wind)]))
    push!(nu2_wind,mean(nu2[1:Int(time_wind)]))
    push!(s1_wind,mean(s1[1:Int(time_wind)]))
    push!(s2_wind,mean(s2[1:Int(time_wind)]))

    for t = 1:Int64((T_total-time_wind)/slide_wind)-1
      push!(nu1_wind,mean(nu1[Int(slide_wind*t):Int(slide_wind*t+time_wind)]))
      push!(nu2_wind,mean(nu2[Int(slide_wind*t):Int(slide_wind*t+time_wind)]))
      push!(s1_wind,mean(s1[Int(slide_wind*t):Int(slide_wind*t+time_wind)]))
      push!(s2_wind,mean(s2[Int(slide_wind*t):Int(slide_wind*t+time_wind)]))
    end

    push!(r1_traj,nu1_wind)
    push!(r2_traj,nu2_wind)

    push!(s1_traj,s1_wind)
    push!(s2_traj,s2_wind)

    r1=0
    r2=0
    t=0

   # push!(rt_list,(t*slide_wind+time_wind/2)*dt-param["T_rest"])

  end



    result = Dict("s1"=>s1_traj,
    "s2"=>s2_traj,
    "r1"=>r1_traj,
    "r2"=> r2_traj
    )

  return result
end
