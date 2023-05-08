"""
################################################################################
# Code for "Grid-optimal energy community planning from a systems perspective"                                         #
Arpan, Selina                                   #
################################################################################
# Based on StochasticPowerModels.jl                                            #
# An extention package of PowerModels.jl for Stochastic (Optimal) Power Flow   #
# See http://github.com/timmyfaraday/StochasticPowerModels.jl                  #
################################################################################
# This example is for Numerical Illustration A where powerflow evaluation is done for rooftop PV capacity
#In this case study power flow is done when using all possible rooftop area.
################################################################################
"""

using Pkg
Pkg.activate(".")
using JuMP
using Ipopt
using PowerModels
using StochasticPowerModels
using PowerModelsDistribution
using JSON
using DataFrames
using CSV
using Statistics
using Plots
# constants 
const PM = PowerModels
const SPM = StochasticPowerModels

# solvers
ipopt_solver = Ipopt.Optimizer

# input
deg  = 2
aux  = true
red  = false


datadir= "C:/Users/karpan/Documents/PhD/Collaboration_selina/Inputfiles_Arpan/Inputfiles_Arpan/" #This need to be changed according tho the data folder, you can make it more easier once we have final data used for paper
#if those data can be open then a new folder can be made in the data section
load_data=CSV.read(datadir*"CaseSpain/Load_data.tab", DataFrame;header=1)
f=[ "All_feeder/65019_73796_configuration.json", #Feeders for Selina
"All_feeder/65019_74430_configuration.json",
"All_feeder/65019_74469_configuration.json",
"All_feeder/65019_74478_configuration.json",
"All_feeder/65019_74559_configuration.json",
"All_feeder/65019_74572_configuration.json",
]

feeder = "All_feeder/65019_74478_configuration.json"  #The curent implementation is only for one feeder


#these are module to rea all the input data
pv_data_S30=CSV.read(datadir*"CaseSpain/PVprofiles/PVprofile_S_30.csv", DataFrame;header=1)
pv_data_SE=CSV.read(datadir*"CaseSpain/PVprofiles/PVprofile_SE.csv", DataFrame;header=1)
pv_data_SW=CSV.read(datadir*"CaseSpain/PVprofiles/PVprofile_SW.csv", DataFrame;header=1)
pv_data_E=CSV.read(datadir*"CaseSpain/PVprofiles/PVprofile_E.csv", DataFrame;header=1)
pv_data_W=CSV.read(datadir*"CaseSpain/PVprofiles/PVprofile_W.csv", DataFrame;header=1)
pv_data=DataFrame(S=pv_data_S30[!,"0"], SE= pv_data_SE[!,"0"], SW=pv_data_SW[!,"0"],E=pv_data_E[!,"0"],W=pv_data_W[!,"0"])
inst_data=CSV.read(datadir*"ConsumerBuildingData.csv", DataFrame;header=1)



load_dist= "beta_lm_2016_8_6.csv"
pv_dist = "beta_pm_2016_8_6.csv"
pov_feeder=[]
# for feeder in f
#feeder = "All_feeder/86315_785381_configuration.json" #50)% error feeder
    file  = joinpath(BASE_DIR, "test/data/Spanish/")   

    data  = SPM.build_mathematical_model_single_phase(file, feeder,load_dist,pv_dist, t_s= 0)

    mn_network= SPM.mn_data(data,load_data,pv_data, inst_data, time_steps=8759)

    r,v=SPM.naive_time_series_analysis(mn_network,8759)
    #6This is to check percentage of the time where voltage in network above 1.05
    #this counts all violations if there are multiple in the same time-period
    ov= sum(v.>1.05)/sum(v.>0)*100
    #however, other interesting way is just to compute the maximum voltage per timestamp as below; this gives actual percentage wher violations occur
    ov_t=sum([maximum(v[j,:]) for j=1:length(v[:,1])].>1.02)/length(v[:,1])
    # the following module is to get branch current using (V2-V1)/Z however it is not working now 
    # Don't really know the reason 
    d=[abs2(im*(r["$t"]["solution"]["bus"]["$(data["branch"]["$b"]["t_bus"])"]["vi"]-r["$t"]["solution"]["bus"]["$(data["branch"]["$b"]["f_bus"])"]["vi"])+
    (r["$t"]["solution"]["bus"]["$(data["branch"]["$b"]["t_bus"])"]["vr"]-r["$t"]["solution"]["bus"]["$(data["branch"]["$b"]["f_bus"])"]["vr"]))/abs2(data["branch"]["$b"]["br_r"]+im*data["branch"]["$b"]["br_x"])
    for t in 12:length(r), b in 1:length(data["branch"])] 
    #For this also I propose an alternative where only the main feeder current is checked.
    #Thats also how it is done in actual practice.
    # Usually a limit of 200 Amp is fixed and 
    current_slack=[r["$j"]["solution"]["gen"]["1"]["pg"]*500/0.23 for j=1:length(r)]
    oc=(sum(current_slack.>200)+sum(current_slack.<-200))/length(r)
    
#this one is one time run to get 30 highest irradiance points
    # max_30=sort(current_slack)[1:30]
    # time_opf=[findall(x->x==c,current_slack)[1] for c in max_30]
    # writedlm("time_opf.csv",time_opf)
#Once it is done is saved in input_files
    
  ## I left the part to repeat to other feeers for you just that you will also get an idea what is happening