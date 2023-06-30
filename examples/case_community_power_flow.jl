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
Pkg.instantiate()
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


datadir= "C:/Users/arpan/OneDrive/Documents/PhD/Collaboration_selina/Inputfiles_Arpan/Inputfiles_Arpan/"#if those data can be open then a new folder can be made in the data section
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
pv_data_S30=CSV.read(datadir*"CaseSpain/PVprofiles/PVprofile_S.csv", DataFrame;header=1)
pv_data_SE=CSV.read(datadir*"CaseSpain/PVprofiles/PVprofile_SE.csv", DataFrame;header=1)
pv_data_SW=CSV.read(datadir*"CaseSpain/PVprofiles/PVprofile_SW.csv", DataFrame;header=1)
pv_data_E=CSV.read(datadir*"CaseSpain/PVprofiles/PVprofile_E.csv", DataFrame;header=1)
pv_data_W=CSV.read(datadir*"CaseSpain/PVprofiles/PVprofile_W.csv", DataFrame;header=1)
pv_data=DataFrame(S=pv_data_S30[!,"0"], SE= pv_data_SE[!,"0"], SW=pv_data_SW[!,"0"],E=pv_data_E[!,"0"],W=pv_data_W[!,"0"])
inst_data=CSV.read(datadir*"ConsumerBuildingData.csv", DataFrame;header=1)

case_data=CSV.read(datadir*"CaseCompare_panel405W_case2b_only.csv", DataFrame; header=1)
cases= ["case1","case2a","case2b","case4"]
load_dist= "beta_lm_2016_8_6.csv"
pv_dist = "beta_pm_2016_8_6.csv"
f_name=[]
OV=[]
OC=[]
EOV=[]
case_name=[]
max_mv=[]
mean_mv=[]
std_mv=[]
pv=[]
global it=0
all_feeder=DataFrame()
for case in cases
    for feeder in f
    # feeder="All_feeder/65019_73796_configuration.json"
    #feeder = "All_feeder/86315_785381_configuration.json" #50)% error feeder
            file  = joinpath(BASE_DIR, "test/data/Spanish/")


            data  = SPM.build_mathematical_model_single_phase(file, feeder,load_dist,pv_dist, t_s= 0)
            [data["bus"]["$i"]["vmax"]=1.15 for i in 1:length(data["bus"])]
            [data["bus"]["$i"]["vmin"]=0.0 for i in 1:length(data["bus"])]

            data_PV=SPM.add_PV_from_case(data,case_data, case)

            mn_network= SPM.mn_data_pf(data_PV,load_data,pv_data, inst_data, time_steps=8760)

            r,v=SPM.naive_time_series_analysis(mn_network,8760)
            #6This is to check percentage of the time where voltage in network above 1.05
            #this counts all violations if there are multiple in the same time-period
            ov= sum(v.>1.05)/sum(v.>0)*100
            #however, other interesting way is just to compute the maximum voltage per timestamp as below; this gives actual percentage wher violations occur
            ov_t=sum([maximum(v[j,:]) for j=1:length(v[:,1])].>1.05)/length(v[:,1])
            ov_th2=sum([maximum(v[j,:]) for j=1:length(v[:,1])].>1.1)/length(v[:,1])
            max_max_ov=maximum([maximum(v[j,:]) for j=1:length(v[:,1])])
            mean_max_ov=mean([maximum(v[j,:]) for j=1:length(v[:,1])])
            std_max_ov=std([maximum(v[j,:]) for j=1:length(v[:,1])])
            # the following module is to get branch current using (V2-V1)/Z however it is not working now 
            # Don't really know the reason 
            # d=[abs2(im*(r["$t"]["solution"]["bus"]["$(data["branch"]["$b"]["t_bus"])"]["vi"]-r["$t"]["solution"]["bus"]["$(data["branch"]["$b"]["f_bus"])"]["vi"])+
            # (r["$t"]["solution"]["bus"]["$(data["branch"]["$b"]["t_bus"])"]["vr"]-r["$t"]["solution"]["bus"]["$(data["branch"]["$b"]["f_bus"])"]["vr"]))/abs2(data["branch"]["$b"]["br_r"]+im*data["branch"]["$b"]["br_x"])
            # for t in 12:length(r), b in 1:length(data["branch"])] 
            #For this also I propose an alternative where only the main feeder current is checked.
            #Thats also how it is done in actual practice.
            # Usually a limit of 200 Amp is fixed and 
            current_slack=[r["$j"]["solution"]["gen"]["1"]["pg"]*500/0.23 for j=1:length(r)]
            oc=(sum(current_slack.>200)+sum(current_slack.<-200))/length(r)
            
            # violations=Dict()
            push!(f_name,feeder[13:end-19])
            push!(case_name,case)
            push!(OV, ov_t)
            push!(EOV, ov_th2)
            push!(OC,oc) 
            push!(max_mv, max_max_ov)
            push!(mean_mv, mean_max_ov)
            push!(std_mv, std_max_ov)
            push!(pv,sum([data_PV["PV"]["$j"]["p_size"][1] for j=1:length(data["PV"])]))
            # df=DataFrame(violations)
            # append!(all_feeder, df)
            
    end
    violations=Dict()
    violations["feederId"]=f_name
    violations["case"]=case_name
    violations["OV"]=OV
    violations["EOV"]=EOV
    violations["OC"]=OC
    violations["max_max_OV"]=max_mv
    violations["mean_max_OV"]=mean_mv
    violations["std_max_OV"]=std_mv
    df=DataFrame(violations)
    CSV.write("HC_violation_v02.csv",df)

end

            violations=Dict()
            violations["feederId"]=f_name
            violations["case"]=case_name
            violations["OV"]=OV
            violations["EOV"]=EOV
            violations["OC"]=OC
            violations["max_max_OV"]=max_mv
            violations["mean_max_OV"]=mean_mv
            violations["std_max_OV"]=std_mv
            df=DataFrame(violations)
            append!(all_feeder, df)
            CSV.write("HC_violation_new_v02.csv",all_feeder)
#this one is one time run to get 30 highest irradiance points
    # max_30=sort(current_slack)[1:30]
    # time_opf=[findall(x->x==c,current_slack)[1] for c in max_30]
    # writedlm("time_opf.csv",time_opf)
#Once it is done is saved in input_files
    
  ## I left the part to repeat to other feeers for you just that you will also get an idea what is happening