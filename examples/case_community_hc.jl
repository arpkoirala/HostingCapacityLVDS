"""
################################################################################
# Code for "Grid-optimal energy community planning from a systems perspective"                                         #
Arpan, Selina                                   #
################################################################################
# Based on StochasticPowerModels.jl                                            #
# An extention package of PowerModels.jl for Stochastic (Optimal) Power Flow   #
# See http://github.com/timmyfaraday/StochasticPowerModels.jl                  #
################################################################################
# This example is for Numerical Illustration for OPF based HC considering the grid limitation
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

voltage_base = 0.230  # (kV)
power_base = 0.5  # (MW)
Z_base = voltage_base^2/power_base # (Ohm)
current_base = power_base/(voltage_base*1e-3) # (A)

# solvers
ipopt_solver = Ipopt.Optimizer

# input
deg  = 2
aux  = true
red  = false

#To be in sync with previous study I propose doing deterministic OPF for 30 timeperiod where possibility of the violation is very high
## the 30 time schedule with highest current is obtained from the  

datadir= "C:/Users/karpan/Documents/PhD/Collaboration_selina/Inputfiles_Arpan/Inputfiles_Arpan/"

feeder = "All_feeder/65019_74478_configuration.json" 

load_data=CSV.read(datadir*"CaseSpain/Load_data.tab", DataFrame;header=1)
f=[ "All_feeder/65019_73796_configuration.json", #Feeders for Selina
"All_feeder/65019_74430_configuration.json",
"All_feeder/65019_74469_configuration.json",
"All_feeder/65019_74478_configuration.json",
"All_feeder/65019_74559_configuration.json",
"All_feeder/65019_74572_configuration.json",
]

pv_data_S30=CSV.read(datadir*"CaseSpain/PVprofiles/PVprofile_S_30.csv", DataFrame;header=1)
pv_data_SE=CSV.read(datadir*"CaseSpain/PVprofiles/PVprofile_SE.csv", DataFrame;header=1)
pv_data_SW=CSV.read(datadir*"CaseSpain/PVprofiles/PVprofile_SW.csv", DataFrame;header=1)
pv_data_E=CSV.read(datadir*"CaseSpain/PVprofiles/PVprofile_E.csv", DataFrame;header=1)
pv_data_W=CSV.read(datadir*"CaseSpain/PVprofiles/PVprofile_W.csv", DataFrame;header=1)
pv_data=DataFrame(S=pv_data_S30[!,"0"], SE= pv_data_SE[!,"0"], SW=pv_data_SW[!,"0"],E=pv_data_E[!,"0"],W=pv_data_W[!,"0"])
inst_data=CSV.read(datadir*"ConsumerBuildingData.csv", DataFrame;header=1)


time_opf=CSV.read(datadir*"time_opf.csv", DataFrame;header=0)
time_opf = time_opf[!,"Column1"]




load_dist= "beta_lm_2016_8_6.csv"
pv_dist = "beta_pm_2016_8_6.csv"
pov_feeder=[]


all_feeder=DataFrame()

for feeder in f[1]
#feeder = "All_feeder/86315_785381_configuration.json" #50)% error feeder
    file  = joinpath(BASE_DIR, "test/data/Spanish/")   

    data  = SPM.build_mathematical_model_single_phase(file, feeder,load_dist,pv_dist, t_s= 0)

    mn_network= SPM.mn_data_opf(data,load_data,pv_data, inst_data, time_opf)
    [mn_network["nw"]["$j"]["branch"]["1"]["rate_a"]=200/current_base for j=1:length(mn_network["nw"])] #limiting the first branch current to 200

    result_hc,p_size= SPM.time_series_hc(mn_network, PM.IVRPowerModel, ipopt_solver, aux=aux, deg=deg, red=red, stochastic=false, time=[1,29])

    #TO have a dictionary with no of panels and device id (meter name).

    no_panels=Dict()
    no_panels["meterId"]=[data["load"]["$i"]["source_id"] for i =1:length(data["PV"])]
    no_panels["p_max"]=[mn_network["nw"]["1"]["PV"]["$i"]["p_max"] for i =1:length(data["PV"])]
    no_panels["p_size_worst"]=[result_hc["1"]["solution"]["PV"]["$i"]["p_size"] for i =1:length(data["PV"])]
    no_panels["p_size_99_5"]=[result_hc["29"]["solution"]["PV"]["$i"]["p_size"] for i =1:length(data["PV"])]
    no_panels["no_of_installation"]=[result_hc["1"]["solution"]["PV"]["$i"]["p_size"]/0.220 for i =1:length(data["PV"])]

    df=DataFrame(no_panels)
    append!(all_feeder, df)
end


