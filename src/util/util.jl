################################################################################
#  Copyright 2021, Tom Van Acker, Arpan Koirala                                #
################################################################################
# StochasticPowerModels.jl                                                     #
# An extention package of PowerModels.jl for Stochastic (Optimal) Power Flow   #
# See http://github.com/timmyfaraday/StochasticPowerModels.jl                  #
################################################################################

# input data
""
function parse_dst(dst, pa, pb, deg)
    dst == "Beta"    && return _PCE.Beta01OrthoPoly(deg, pa, pb; Nrec=5*deg)
    dst == "Normal"  && return _PCE.GaussOrthoPoly(deg; Nrec=5*deg)
    dst == "Uniform" && return _PCE.Uniform01OrthoPoly(deg; Nrec=5*deg)
end


""
function parse_dst_beta(dst, pa, pb, deg)
    dst == "Beta"    && return _PCE.Beta01OrthoPoly(deg, pa, pb; Nrec=5*deg)
    #dst == "Normal"  && return _PCE.GaussOrthoPoly(deg; Nrec=5*deg)
    #dst == "Uniform" && return _PCE.Uniform01OrthoPoly(deg; Nrec=5*deg)
end
"""
    StochasticPowerModels.build_stochastic_data(data::Dict{String,Any}, deg::Int)

Function to build the multi-network data representative of the polynomial chaos
expansion of a single-network data dictionary.
"""
function build_stochastic_data(data::Dict{String,Any}, deg::Int)
    # add maximum current
    for (nb, branch) in data["branch"]
        f_bus = branch["f_bus"]
        branch["cmax"] = branch["rate_a"] / data["bus"]["$f_bus"]["vmin"]
    end

    # build mop
    opq = [parse_dst(ns[2]["dst"], ns[2]["pa"], ns[2]["pb"], deg) for ns in data["sdata"]]
    mop = _PCE.MultiOrthoPoly(opq, deg)

    # build load matrix
    Nd, Npce = length(data["load"]), mop.dim
    pd, qd = zeros(Nd, Npce), zeros(Nd, Npce)
    for nd in 1:Nd 
        # reactive power
        qd[nd,1] = data["load"]["$nd"]["qd"]
        # active power
        nb = data["load"]["$nd"]["load_bus"]
        ni = data["bus"]["$nb"]["dst_id"]
        if ni == 0
            pd[nd,1] = data["load"]["$nd"]["pd"]
        else
            base = data["baseMVA"]
            μ, σ = data["bus"]["$nb"]["μ"] / base, data["bus"]["$nb"]["σ"] / base
            if mop.uni[ni] isa _PCE.GaussOrthoPoly
                pd[nd,[1,ni+1]] = _PCE.convert2affinePCE(μ, σ, mop.uni[ni])
            else
                pd[nd,[1,ni+1]] = _PCE.convert2affinePCE(μ, σ, mop.uni[ni], kind="μσ")
            end
        end
    end

    # replicate the data
    data = _PM.replicate(data, Npce)

    # add the stochastic data 
    data["T2"] = _PCE.Tensor(2,mop)
    data["T3"] = _PCE.Tensor(3,mop)
    data["T4"] = _PCE.Tensor(4,mop)
    data["mop"] = mop
    for nw in 1:Npce, nd in 1:Nd
        data["nw"]["$nw"]["load"]["$nd"]["pd"] = pd[nd,nw]
        data["nw"]["$nw"]["load"]["$nd"]["qd"] = qd[nd,nw]
    end

    return data
end

# output data
"""
    StochasticPowerModels.pce_coeff(result, element::String, id::Int, var::String)

Returns all polynomial chaos coefficients associated with the variable `var` of 
the `id`th element `element`.
"""
pce_coeff(result, element::String, id::Int, var::String) =
    [nw[2][element]["$id"][var] for nw in sort(collect(result["solution"]["nw"]), by=x->parse(Int,x[1]))]

"""
    StochasticPowerModels.sample(sdata, result, element::String, id::Int, var::String; sample_size::Int=1000)

Return an `sample_size` sample of the variable `var` of the `id`th element 
`element`.
"""
sample(result, element::String, id::Int, var::String; sample_size::Int=1000) =
    _PCE.samplePCE(sample_size, pce_coeff(result, element, id, var), result["mop"])

"""
    StochasticPowerModels.density(sdata, result, element::String, id::Int, var::String; sample_size::Int=1000)

Return an kernel density estimate of the variable `var` of the `id`th element 
`element`.
"""
density(result, element::String, id::Int, var::String; sample_size::Int=1000) =
    _KDE.kde(sample(result, element, id, var; sample_size=sample_size))

function print_summary(obj::Dict{String,<:Any}; kwargs...)
    if _IM.ismultinetwork(obj)
        for (n,nw) in obj["nw"]
            println("----------------")
            println("PCE index $n")
            _PM.summary(stdout, nw; kwargs...)
        end
    end
end



"Converts JSON file of three phase DN to single phase equivalent"
function build_mathematical_model_single_phase(dir, config_file_name, load_dist_csv, pv_dist_csv; t_s=52, pd = 0.0, qd = 0.0, scale_factor = 1.0, curt=0.0, cross_area_fact=1.0)
#configuration = "star"

"""
Specify voltage and power base, voltage base should be the phase-to-ground voltage
of the feeder studied (kV), the power base can be arbitrairly chosen (MW)
"""
voltage_base = 0.230  # (kV)
power_base = 0.5  # (MW)
Z_base = voltage_base^2/power_base # (Ohm)
current_base = power_base/(voltage_base*1e-3) # (A)

mwpu = 1/power_base
kwpu = (1e-3)/power_base

network_model = Dict{String,Any}()
configuration_json_dict = Dict{Any,Any}()
device_df=CSV.read(dir*config_file_name[1:length(config_file_name)-19]*".csv", DataFrame)

dist_lv=CSV.read(dir*load_dist_csv, DataFrame)
dist_pv=CSV.read(dir*pv_dist_csv, DataFrame)
# dist_pv=CSV.read(dir*"beta_pm_2022_181"*".csv", DataFrame)
dist_pv_ts= dist_pv[in([t_s]).(dist_pv.timeslot),:]
dist_lv_ts=dist_lv[in([t_s]).(dist_lv.timeslot),:]

dist_lv_ts_feeder = dist_lv_ts[in(unique(device_df.category)).(dist_lv_ts.cluster),:]
s_dict=Dict()
i=1
for dist in eachrow(dist_lv_ts_feeder)
    s=Dict()
    s["dst"]= "Beta"
    s["dst_id"] = dist["cluster"]
    s["pa"]= dist["alpha"]
    s["pb"]= dist["beta"]
    s["pc"]= dist["lower"]
    s["pd"]= dist["lower"]+dist["upper"]
    s_dict[string(i)] = s
    i=i+1
end


##add Irradiance if day time or there is some Irradiance
if dist_pv_ts.upper[1]>0
    s=Dict()
    s["dst"]= "Beta"
    s["dst_id"] = 55
    s["pa"]= dist_pv_ts[!,"alpha"][1]
    s["pb"]= dist_pv_ts[!,"beta"][1]
    s["pc"]= dist_pv_ts[!,"lower"][1]
    s["pd"]= dist_pv_ts[!,"lower"][1]+dist_pv_ts[!,"upper"][1]
    s_dict[string(i)] = s
end


#network_model["is_kron_reduced"] = true
network_model["p_factor"] = 0.95
network_model["q_factor"] = sqrt(1-0.95^2)
network_model["dcline"] = Dict{String,Any}()
network_model["switch"] = Dict{String,Any}()
#network_model["is_projected"] = true
network_model["per_unit"] = true
#network_model["data_model"] = MATHEMATICAL
network_model["shunt"] = Dict{String,Any}()
network_model["transformer"] = Dict{String,Any}()
network_model["bus"] = Dict{String,Any}()
network_model["map"] = Dict{String,Any}()
#network_model["conductors"] = 1
network_model["baseMVA"] =  power_base
network_model["basekv"] =  voltage_base
network_model["bus_lookup"] = Dict{Any,Int64}()
network_model["run_type"] = 1
network_model["load"] = Dict{String,Any}()
network_model["gen"] = Dict{String,Any}("1" => Dict{String,Any}(
"pg"            => 0.2,
"model"         => 2,
#"connections"   => [1, 2, 3],
"shutdown"      => 0.0,
"startup"       => 0.0,
#"configuration" => WYE,
"name"          => "virtual_generator",
"qg"            => 0.0,
"gen_bus"       => 1,
"vbase"         =>  voltage_base,
"source_id"     => Any["gen",1],
"index"         => 1,
"cost"          => [20000.0, 1400.0, 0.0],
"gen_status"    => 1,
"qmax"          => 1.275,
"qmin"          => -1.275,
"pmax"          => 1.5,
"pmin"          => -1.5,
"ncost"         => 3,
"λpmin"         => 1.65, #1.03643 ,
"λpmax"         => 1.65, #1.03643 ,
"λqmin"         => 1.65, #1.03643 ,
"λqmax"         => 1.65 #1.03643
))
network_model["settings"] = Dict{String,Any}(
"sbase_default"        => power_base,
"vbases_default"       => Dict{String,Any}(), #No default is specified for now, since default is never used
"voltage_scale_factor" => 1E3, #Voltages are thus expressed in kV
"sbase"                => power_base,
"power_scale_factor"   => 1E6, #Power is expressed in MW
"base_frequency"       => 50.0 #Hertz
)
network_model["branch"] = Dict{String,Any}()
network_model["storage"] = Dict{String,Any}()
open(dir * config_file_name,"r") do io
configuration_json_dict = JSON.parse(io)
end;
#voltage_base = configuration_json_dict["gridConfig"]["basekV"]
#power_base = configuration_json_dict["gridConfig"]["baseMVA"]
sub_dir=splitpath(config_file_name)[1]
configuration = configuration_json_dict["gridConfig"]["connection_configuration"]
branches_file_name = sub_dir*"/"*splitpath(configuration_json_dict["gridConfig"]["branches_file"])[2]
buses_file_name = sub_dir*"/"*splitpath(configuration_json_dict["gridConfig"]["buses_file"])[2]
devices_file_name = sub_dir*"/"*splitpath(configuration_json_dict["gridConfig"]["devices_file"])[2]


open(dir * buses_file_name,"r") do io
buses_json_dict = JSON.parse(io)
    for bus in buses_json_dict
        id = bus["busId"] + 1 #Indexing starts at one in Julia
        id_s = string(id)
        network_model["bus_lookup"][id_s] = id
        network_model["settings"]["vbases_default"][id_s] =  voltage_base

        if id == 1 #Settings for slack bus
            network_model["bus"][id_s] = Dict{String,Any}(
                "name"      => "slack",
                "bus_type"  => 3,
                ##"grounded"  => Bool[0, 0, 0],
                #"terminals" => [1, 2, 3],
                "vbase"     =>  voltage_base,
                "index"     => id,
                "bus_i"     => id,
                "λvmin"     => 1.65, #1.03643,
                "λvmax"     => 1.65, #1.03643,
                "vmin"      => 1.0,
                "vmax"      => 1,
                "va"        => 0.0,
                "vm"        => 1, 
                #"LPp"       => Dict(c => Dict(b => 0.0 for b in 1:length(keys(buses_json_dict))) for c in 1:3),
                #"LPq"       => Dict(c => Dict(b => 0.0 for b in 1:length(keys(buses_json_dict))) for c in 1:3),
                #"LQp"       => Dict(c => Dict(b => 0.0 for b in 1:length(keys(buses_json_dict))) for c in 1:3),
                #"LQq"       => Dict(c => Dict(b => 0.0 for b in 1:length(keys(buses_json_dict))) for c in 1:3),
                #"run_type"  => 1
                )
        else
            network_model["bus"][id_s] = Dict{String,Any}(
                "name"      => id_s,
                "bus_type"  => 1,
                #"grounded"  => Bool[0, 0, 0],
                #"terminals" => [1, 2, 3],
                "vbase"     =>  voltage_base,
                "index"     => id,
                "bus_i"     => id,
                "λvmin"     => 1.65,#1.03643,
                "λvmax"     => 1.65, #1.03643,
                "vmin"      =>  0.95,
                "vmax"      => 1.05, 
                #"LPp"       => Dict(c => Dict(b => 0.0 for b in 1:length(keys(buses_json_dict))) for c in 1:3),
                #"LPq"       => Dict(c => Dict(b => 0.0 for b in 1:length(keys(buses_json_dict))) for c in 1:3),
                #"LQp"       => Dict(c => Dict(b => 0.0 for b in 1:length(keys(buses_json_dict))) for c in 1:3),
                #"LQq"       => Dict(c => Dict(b => 0.0 for b in 1:length(keys(buses_json_dict))) for c in 1:3),
                "run_type"  => 1
                )
        end;
    end;
end;

#print(device_df)
open(dir * devices_file_name,"r") do io
devices_json_dict = JSON.parse(io)
  for device in devices_json_dict["LVcustomers"]
    id = device["deviceId"] + 1 #Indexing starts at one in Julia
    d=device_df[in(id-1).(device_df.dev_id),:]
    id_s = string(id)
    μ = dist_lv_ts_feeder[in(d[!,"category"][1]).(dist_lv_ts_feeder.cluster),:][!,"lower"][1]
    σ  = dist_lv_ts_feeder[in(d[!,"category"][1]).(dist_lv_ts_feeder.cluster),:][!,"upper"][1] 
    cons = convert(Float64,device["yearlyNetConsumption"])
    network_model["load"][id_s] = Dict{String,Any}(
        #"connections"   => vec(Int.(device["phases"])),
        "name"          => id_s*"-"*device["coded_ean"],
        "status"        => 1,
        "vbase"         =>  voltage_base,
        "vnom_kv"       => 1.0,
        "source_id"     => device["coded_ean"],
        "load_bus"      => device["busId"] + 1,
        "dispatchable"  => 0,
        "index"         => id,
        "yearlyNetConsumption" => cons,
        #"phases"        => device["phases"],
        "pd"            => max(0.1, μ)/1e3/power_base/ 3,
        "qd"            =>  max(0.01,μ)/1e3/ power_base/ 3/10,
        "p_inj"         => 0.0,
        "q_inj"         => 0.0,
        "conn_cap_kW"   => device["connectionCapacity"],
        "dst_id" => d[!,"category"][1],
        "cluster_id"  => findall(x->x==1,[s_dict["$i"]["dst_id"]==d[!,"category"][1] for i=1:length(s_dict)])[1],
        "μ"  => μ,
        "σ"  => σ 

    )
    end;
end;

open(dir * branches_file_name,"r") do io
branches_json_dict = JSON.parse(io)
impedance_dict = Dict{String,Any}(
"BT - Desconocido BT" => [0.21, 0.075],
"BT - MANGUERA" => [0.3586, 0.089], #[1.23, 0.08],
"BT - RV 0,6/1 KV 2*16 KAL" => [2.14, 0.09], #16 = 20 A
"BT - RV 0,6/1 KV 2*25 KAL" => [1.34, 0.097],
"BT - RV 0,6/1 KV 3(1*150 KAL) + 1*95 KAL" => [0.2309, 0.085],
"BT - RV 0,6/1 KV 3(1*240 KAL) + 1*150 KAL" => [0.1602, 0.079],
"BT - RV 0,6/1 KV 3(1*240 KAL) + 1*95 KAL" => [0.1602, 0.079],
"BT - RV 0,6/1 KV 4*25 KAL" => [1.34, 0.097],
"BT - RV 0,6/1 KV 4*50 KAL" => [0.71849, 0.093],
"BT - RV 0,6/1 KV 4*95 KAL" => [0.3586, 0.089],
"BT - RX 0,6/1 KV 2*16 Cu" => [1.23, 0.08],
"BT - RX 0,6/1 KV 2*2 Cu" => [9.9, 0.075],
"BT - RX 0,6/1 KV 2*4 Cu" => [4.95, 0.075],
"BT - RX 0,6/1 KV 2*6 Cu" => [3.3, 0.075],
"BT - RZ 0,6/1 KV 2*16 AL" => [2.14, 0.09],
"BT - RZ 0,6/1 KV 3*150 AL/80 ALM" => [0.2309, 0.85],
"BT - RZ 0,6/1 KV 3*150 AL/95 ALM" => [0.2309, 0.85],
"BT - RZ 0,6/1 KV 3*25 AL/54,6 ALM" => [1.34, 0.097],
"BT - RZ 0,6/1 KV 3*35 AL/54,6 ALM" => [0.9073, 0.095],
"BT - RZ 0,6/1 KV 3*50 AL/54,6 ALM" => [0.718497, 0.093],
"BT - RZ 0,6/1 KV 3*70 ALM/54,6 AL" => [0.4539, 0.091],
"BT - RZ 0,6/1 KV 3*95 AL/54,6 ALM" => [0.3586, 0.089],
"BT - RZ 0,6/1 KV 4*16 AL" => [2.14, 0.09],
"aansluitkabel" => [1.15, 0.150]
)

open(dir * branches_file_name,"r") do io
branches_json_dict = JSON.parse(io)
currentmax_dict = Dict{String,Any}(
"BT - Desconocido BT" => 200,
"BT - MANGUERA" => 150, #200 certain  40.18#150
"BT - RV 0,6/1 KV 2*16 KAL" => 75,
"BT - RV 0,6/1 KV 2*25 KAL" => 100,
"BT - RV 0,6/1 KV 3(1*150 KAL) + 1*95 KAL" => 305,
"BT - RV 0,6/1 KV 3(1*240 KAL) + 1*150 KAL" => 344,
"BT - RV 0,6/1 KV 3(1*240 KAL) + 1*95 KAL" => 344,
"BT - RV 0,6/1 KV 4*25 KAL" => 100,
"BT - RV 0,6/1 KV 4*50 KAL" => 150,
"BT - RV 0,6/1 KV 4*95 KAL" => 230,
"BT - RX 0,6/1 KV 2*16 Cu" => 75,#95,
"BT - RX 0,6/1 KV 2*2 Cu" => 40, #30,
"BT - RX 0,6/1 KV 2*4 Cu" => 60,#40,
"BT - RX 0,6/1 KV 2*6 Cu" => 80, #50,
"BT - RZ 0,6/1 KV 2*16 AL" => 20, #75,
"BT - RZ 0,6/1 KV 3*150 AL/80 ALM" => 305, #264,
"BT - RZ 0,6/1 KV 3*150 AL/95 ALM" => 305,#264,
"BT - RZ 0,6/1 KV 3*25 AL/54,6 ALM" => 100, #78.98,
"BT - RZ 0,6/1 KV 3*35 AL/54,6 ALM" => 120,
"BT - RZ 0,6/1 KV 3*50 AL/54,6 ALM" => 150, #118.47,
"BT - RZ 0,6/1 KV 3*70 ALM/54,6 AL" => 160,
"BT - RZ 0,6/1 KV 3*95 AL/54,6 ALM" => 230, # 182.21,
"BT - RZ 0,6/1 KV 4*16 AL" => 75,
"aansluitkabel" => 120 #200
)


    for branch in branches_json_dict
        id = branch["branchId"] +1
        id_s = string(id)
        network_model["branch"][id_s] = Dict{String,Any}(
            "shift"         => 0.0,
            #"f_connections" => [1, 2, 3],
            "name"          => id_s,
            "switch"        => false,
            "g_to"          => 0.0,
            "c_rating_a"    => currentmax_dict[branch["cableType"]],
            "vbase"         =>  voltage_base,
            "g_fr"          => 0.0,
            #"t_connections" => [1, 2, 3],
            "f_bus"         => branch["upBusId"]+1,
            "b_fr"          => 0.0,
            "c_rating_b"    => currentmax_dict[branch["cableType"]],
            "br_status"     => 1,
            "t_bus"         => branch["downBusId"]+1,
            "b_to"          => 0.0,
            "index"         => id,
            "angmin"        => -1.0472,
            "angmax"        => 1.0472,
            "transformer"   => false,
            "tap"           => 1.0,
            "c_rating_c"    => currentmax_dict[branch["cableType"]], 
            "λcmax"         => 1.65 #1.65 #2.5 #1.03643    
            )

        if haskey(impedance_dict,branch["cableType"])
            # network_model["branch"][id_s]["br_x"] = (cross_area_fact^(1/4)).* impedance_dict[branch["cableType"]][2] .* (branch["cableLength"]/1000+1E-6)./  Z_base 
            # network_model["branch"][id_s]["br_r"] = (1/cross_area_fact).* impedance_dict[branch["cableType"]][1] .* (branch["cableLength"]/1000+1E-6)./  Z_base
            network_model["branch"][id_s]["br_x"] = (1/cross_area_fact).* impedance_dict[branch["cableType"]][2] .* (branch["cableLength"]/1000+1E-6)./  Z_base 
            network_model["branch"][id_s]["br_r"] = (1/cross_area_fact).* impedance_dict[branch["cableType"]][1] .* (branch["cableLength"]/1000+1E-6)./  Z_base  
        end;
        

        if haskey(currentmax_dict,branch["cableType"])
            network_model["branch"][id_s]["rate_a"] = cross_area_fact.*((currentmax_dict[branch["cableType"]]*voltage_base)/1e3)/power_base
            
            #network_model["branch"][id_s]["I_rating"] = currentmax_dict[branch["cableType"]]/current_base

        end;
    end;
end;
end;


network_model["sdata"]= s_dict 
network_model["curt"]= curt

network_model["PV"]=deepcopy(network_model["load"]);
[network_model["PV"][d]["μ"]=s_dict[string(length(s_dict))]["pc"] for d in   keys(network_model["PV"])]
[network_model["PV"][d]["σ"]=s_dict[string(length(s_dict))]["pd"] for d in   keys(network_model["PV"])]
[network_model["PV"][d]["pd"]=s_dict[string(length(s_dict))]["pd"]/1e6/ power_base / 3 for d in   keys(network_model["PV"])]
return network_model
end;


"""
    StochasticPowerModels.build_stochastic_data_hc(data::Dict{String,Any}, deg::Int)

Function to build the multi-network data representative of the polynomial chaos
expansion of a single-network data dictionary for HC problem in DN.
"""
function build_stochastic_data_hc(data::Dict{String,Any}, deg::Int, t_s=50)
    # add maximum current
    curt=data["curt"]
    for (nb, branch) in data["branch"]
        f_bus = branch["f_bus"]
        branch["cmax"] = branch["rate_a"] / data["bus"]["$f_bus"]["vmin"]
    end


    # build mop
    opq = [parse_dst_beta(ns[2]["dst"], ns[2]["pa"], ns[2]["pb"], deg) for ns in sort(data["sdata"])]
    mop = _PCE.MultiOrthoPoly(opq, deg)

    # build load matrix
    Nd, Npce = length(data["load"]), mop.dim
    pd, qd = zeros(Nd, Npce), zeros(Nd, Npce)
    pd_g, qd_g = zeros(Nd, Npce), zeros(Nd, Npce)
    for nd in 1:Nd 
        # reactive power
        qd[nd,1] = data["load"]["$nd"]["qd"]
        # active power
        nb = data["load"]["$nd"]["load_bus"]
        ni = data["load"]["$nd"]["cluster_id"]
        if ni == 55
            pd[nd,1] = data["load"]["$nd"]["pd"]
        else
            base = data["baseMVA"]
            μ, σ = data["load"]["$nd"]["μ"] /1e3/ base/ 3, data["load"]["$nd"]["σ"] /1e3/ base/3
            if mop.uni[ni] isa _PCE.GaussOrthoPoly
                pd[nd,[1,ni+1]] = _PCE.convert2affinePCE(μ, σ, mop.uni[ni])
            else
                pd[nd,[1,ni+1]] = _PCE.convert2affinePCE(μ, σ, mop.uni[ni])
            end
        end
        np = length(opq)
        base = data["baseMVA"]
        μ, σ = data["PV"]["1"]["μ"]/1e6 / base / 3, data["PV"]["1"]["σ"] /1e6/ base / 3
        
            if mop.uni[np] isa _PCE.GaussOrthoPoly
                pd_g[nd,[1,np+1]] = _PCE.convert2affinePCE(μ, σ, mop.uni[np])
            else
                pd_g[nd,[1,np+1]] = _PCE.convert2affinePCE(μ, σ, mop.uni[np])
            end
        
    end

    # replicate the data
    data = _PM.replicate(data, Npce)

    # add the stochastic data 
    data["T2"] = _PCE.Tensor(2,mop)
    data["T3"] = _PCE.Tensor(3,mop)
    data["T4"] = _PCE.Tensor(4,mop)
    data["mop"] = mop
    data["curt"]= curt
    for nw in 1:Npce, nd in 1:Nd
       
        data["nw"]["$nw"]["load"]["$nd"]["pd"] = pd[nd,nw]
        data["nw"]["$nw"]["load"]["$nd"]["qd"] = qd[nd,nw]
    end
    
    for nw in 1:Npce, nd in 1:Nd

        data["nw"]["$nw"]["PV"]["$nd"]["pd"] = pd_g[nd,nw]
        data["nw"]["$nw"]["PV"]["$nd"]["qd"] = qd_g[nd,nw]
    end

    return data
end




"Converts JSON file of three phase DN to single phase equivalent"
function build_mathematical_model_single_phase_market(dir, config_file_name,load_dist_csv, pv_dist_csv, load_dir, load_folder; grid_only=true, t_s=52, pd = 0.0, qd = 0.0, scale_factor = 1.0, curt=0.0, cross_area_fact=1.0)
#configuration = "star"

"""
Specify voltage and power base, voltage base should be the phase-to-ground voltage
of the feeder studied (kV), the power base can be arbitrairly chosen (MW)
"""
voltage_base = 0.230  # (kV)
power_base = 0.5  # (MW)
Z_base = voltage_base^2/power_base # (Ohm)
current_base = power_base/(voltage_base*1e-3) # (A)

mwpu = 1/power_base
kwpu = (1e-3)/power_base

PVins= CSV.read(load_dir*"/PVinstallations_$load_folder.csv", DataFrame)

demand_p=CSV.read(load_dir*"/demand_P_kW_$load_folder.csv", DataFrame;header=1)
rename!(demand_p, string.(PVins.consumer_index))
demand_q=CSV.read(load_dir*"/demand_Q_kVAR_$load_folder.csv", DataFrame;header=1)
rename!(demand_q, string.(PVins.consumer_index))
generation_p=CSV.read(load_dir*"/generation_P_kW_$load_folder.csv", DataFrame;header=1)
rename!(generation_p, string.(PVins.consumer_index))
generation_q=CSV.read(load_dir*"/generation_Q_kVAR_$load_folder.csv", DataFrame;header=1)
rename!(generation_q, string.(PVins.consumer_index))

solar_AF=CSV.read(load_dir*"/solar_AF.csv", DataFrame;header=1)

PVins=PVins[PVins.EAN.!="[\"none\"]",:]
# PVsize=Dict(["$(PVins.EAN[i])" => parse.(Float64, split(chop(PVins.installed_PVcap_kW[i]; head=1, tail=1), ',')) for i=1:size(PVins)[1]])
# PVorien= Dict([if length(PVins.PVorientations[i])>5 "$(PVins.EAN[i])" =>(PVins.PVorientations[i][3],PVins.PVorientations[i][8]) else "$(PVins.EAN[i])" =>(PVins.PVorientations[i][3]) end for i=1:size(PVins)[1]])
dev_ean=[[split(PVins.EAN[i][1:length(PVins.EAN[i])-1],",")[j][3:15] for j in 1:length(split(PVins.EAN[i][2:length(PVins.EAN[i])-1],","))] for i in 1:length(PVins.EAN)]
network_model = Dict{String,Any}()
configuration_json_dict = Dict{Any,Any}()
device_df=CSV.read(dir*config_file_name[1:length(config_file_name)-19]*".csv", DataFrame)
PVsize=Dict() 
solar_cost = Dict()
for i=1:length(dev_ean)
    for j=1:length(dev_ean[i])
        if dev_ean[i][j][1] =="?"
            dev_ean[i][j]= "??"*dev_ean[i][j]
        end
       merge!(PVsize, Dict("$(dev_ean[i][j])" => PVins.installed_PVcap_kW[i]/length(dev_ean[i]))) 
       merge!(solar_cost, Dict("$(dev_ean[i][j])" => PVins.solar_cost_eurkW[i]))   
    end
end
dist_lv=CSV.read(dir*load_dist_csv, DataFrame)
dist_pv=CSV.read(dir*pv_dist_csv, DataFrame)
# dist_pv=CSV.read(dir*"beta_pm_2022_181"*".csv", DataFrame)
dist_pv_ts= dist_pv[in([t_s]).(dist_pv.timeslot),:]
dist_lv_ts=dist_lv[in([t_s]).(dist_lv.timeslot),:]

dist_lv_ts_feeder = dist_lv_ts[in(unique(device_df.category)).(dist_lv_ts.cluster),:]

#uncertainties based on day 10
# for i=1:nrow(dist_lv_ts_feeder)
#     dev=device_df[device_df.category.==dist_lv_ts_feeder.cluster[i],:].consumer_id.+1
#     dist_lv_ts_feeder.upper[i]=maximum(maximum.([demand_p[!,"$(i)"] for i in dev]))
#     dist_lv_ts_feeder.lower[i]= minimum(minimum.([demand_p[!,"$(i)"] for i in dev]))
# end

s_dict=Dict()
i=1
for dist in eachrow(dist_lv_ts_feeder)
    s=Dict()
    s["dst"]= "Beta"
    s["dst_id"] = dist["cluster"]
    s["pa"]= dist["alpha"]
    s["pb"]= dist["beta"]
    s["pc"]= dist["lower"]
    s["pd"]= dist["upper"]
    s_dict[string(i)] = s
    i=i+1
end


##add Irradiance if day time or there is some Irradiance
if dist_pv_ts.upper[1]>0
    s=Dict()
    s["dst"]= "Beta"
    s["dst_id"] = 55
    s["pa"]= dist_pv_ts[!,"alpha"][1]
    s["pb"]= dist_pv_ts[!,"beta"][1]
    s["pc"]= dist_pv_ts[!,"lower"][1]#minimum(solar_AF[!,"AF_E"][222:230])
    s["pd"]= dist_pv_ts[!,"lower"][1]+dist_pv_ts[!,"upper"][1] #maximum(solar_AF[!,"AF_E"][221:231])
    s_dict[string(i)] = s
end


#network_model["is_kron_reduced"] = true
network_model["p_factor"] = 0.95
network_model["q_factor"] = sqrt(1-0.95^2)
network_model["dcline"] = Dict{String,Any}()
network_model["switch"] = Dict{String,Any}()
#network_model["is_projected"] = true
network_model["per_unit"] = true
#network_model["data_model"] = MATHEMATICAL
network_model["shunt"] = Dict{String,Any}()
network_model["transformer"] = Dict{String,Any}()
network_model["bus"] = Dict{String,Any}()
network_model["map"] = Dict{String,Any}()
#network_model["conductors"] = 1
network_model["baseMVA"] =  power_base
network_model["basekv"] =  voltage_base
network_model["bus_lookup"] = Dict{Any,Int64}()
network_model["run_type"] = 1
network_model["load"] = Dict{String,Any}()
network_model["gen"] = Dict{String,Any}("1" => Dict{String,Any}(
"pg"            => 0.2,
"model"         => 2,
#"connections"   => [1, 2, 3],
"shutdown"      => 0.0,
"startup"       => 0.0,
#"configuration" => WYE,
"name"          => "virtual_generator",
"qg"            => 0.0,
"gen_bus"       => 1,
"vbase"         =>  voltage_base,
"source_id"     => Any["gen",1],
"index"         => 1,
"cost"          => [20000.0, 1400.0, 0.0],
"gen_status"    => 1,
"qmax"          => 1.275,
"qmin"          => -1.275,
"pmax"          => 1.5,
"pmin"          => -1.5,
"ncost"         => 3,
"λpmin"         => 1.65, #1.03643 ,
"λpmax"         => 1.65, #1.03643 ,
"λqmin"         => 1.65, #1.03643 ,
"λqmax"         => 1.65 #1.03643
))
network_model["settings"] = Dict{String,Any}(
"sbase_default"        => power_base,
"vbases_default"       => Dict{String,Any}(), #No default is specified for now, since default is never used
"voltage_scale_factor" => 1E3, #Voltages are thus expressed in kV
"sbase"                => power_base,
"power_scale_factor"   => 1E6, #Power is expressed in MW
"base_frequency"       => 50.0 #Hertz
)
network_model["branch"] = Dict{String,Any}()
network_model["storage"] = Dict{String,Any}()
open(dir * config_file_name,"r") do io
configuration_json_dict = JSON.parse(io)
end;
#voltage_base = configuration_json_dict["gridConfig"]["basekV"]
#power_base = configuration_json_dict["gridConfig"]["baseMVA"]
sub_dir=splitpath(config_file_name)[1]
configuration = configuration_json_dict["gridConfig"]["connection_configuration"]
branches_file_name = sub_dir*"/"*splitpath(configuration_json_dict["gridConfig"]["branches_file"])[2]
buses_file_name = sub_dir*"/"*splitpath(configuration_json_dict["gridConfig"]["buses_file"])[2]
devices_file_name = sub_dir*"/"*splitpath(configuration_json_dict["gridConfig"]["devices_file"])[2]


open(dir * buses_file_name,"r") do io
buses_json_dict = JSON.parse(io)
    for bus in buses_json_dict
        id = bus["busId"] + 1 #Indexing starts at one in Julia
        id_s = string(id)
        network_model["bus_lookup"][id_s] = id
        network_model["settings"]["vbases_default"][id_s] =  voltage_base

        if id == 1 #Settings for slack bus
            network_model["bus"][id_s] = Dict{String,Any}(
                "name"      => "slack",
                "bus_type"  => 3,
                ##"grounded"  => Bool[0, 0, 0],
                #"terminals" => [1, 2, 3],
                "vbase"     =>  voltage_base,
                "index"     => id,
                "bus_i"     => id,
                "λvmin"     => 1.65, #1.03643,
                "λvmax"     => 1.65, #1.03643,
                "vmin"      => 1.0,
                "vmax"      => 1,
                "va"        => 0.0,
                "vm"        => 1, 
                #"LPp"       => Dict(c => Dict(b => 0.0 for b in 1:length(keys(buses_json_dict))) for c in 1:3),
                #"LPq"       => Dict(c => Dict(b => 0.0 for b in 1:length(keys(buses_json_dict))) for c in 1:3),
                #"LQp"       => Dict(c => Dict(b => 0.0 for b in 1:length(keys(buses_json_dict))) for c in 1:3),
                #"LQq"       => Dict(c => Dict(b => 0.0 for b in 1:length(keys(buses_json_dict))) for c in 1:3),
                #"run_type"  => 1
                )
        else
            network_model["bus"][id_s] = Dict{String,Any}(
                "name"      => id_s,
                "bus_type"  => 1,
                #"grounded"  => Bool[0, 0, 0],
                #"terminals" => [1, 2, 3],
                "vbase"     =>  voltage_base,
                "index"     => id,
                "bus_i"     => id,
                "λvmin"     => 1.65,#1.03643,
                "λvmax"     => 1.65, #1.03643,
                "vmin"      =>  0.90,
                "vmax"      => 1.05, 
                #"LPp"       => Dict(c => Dict(b => 0.0 for b in 1:length(keys(buses_json_dict))) for c in 1:3),
                #"LPq"       => Dict(c => Dict(b => 0.0 for b in 1:length(keys(buses_json_dict))) for c in 1:3),
                #"LQp"       => Dict(c => Dict(b => 0.0 for b in 1:length(keys(buses_json_dict))) for c in 1:3),
                #"LQq"       => Dict(c => Dict(b => 0.0 for b in 1:length(keys(buses_json_dict))) for c in 1:3),
                "run_type"  => 1
                )
        end;
    end;
end;

#print(device_df)
open(dir * devices_file_name,"r") do io
devices_json_dict = JSON.parse(io)
  for device in devices_json_dict["LVcustomers"]
    id = device["deviceId"] + 1 #Indexing starts at one in Julia
    d=device_df[in(id-1).(device_df.dev_id),:]
    id_s = string(id)
    μ = dist_lv_ts_feeder[in(d[!,"category"][1]).(dist_lv_ts_feeder.cluster),:][!,"lower"][1]
    σ  = dist_lv_ts_feeder[in(d[!,"category"][1]).(dist_lv_ts_feeder.cluster),:][!,"upper"][1] 
    cons = convert(Float64,device["yearlyNetConsumption"])
    
    if startswith(device["coded_ean"], '?')
        device["coded_ean"]= "??"*device["coded_ean"]
    end
    network_model["load"][id_s] = Dict{String,Any}(
        #"connections"   => vec(Int.(device["phases"])),
        "name"          => id_s*"-"*device["coded_ean"],
        "status"        => 1,
        "vbase"         =>  voltage_base,
        "vnom_kv"       => 1.0,
        "source_id"     => device["coded_ean"],
        "load_bus"      => device["busId"] + 1,
        "dispatchable"  => 0,
        "index"         => id,
        "yearlyNetConsumption" => cons,
        #"phases"        => device["phases"],
        "pd"            => max(0.1, μ)/1e3/power_base/ 3,
        "qd"            =>  max(0.01,μ)/1e3/ power_base/ 3/10,
        "p_inj"         => 0.0,
        "q_inj"         => 0.0,
        "conn_cap_kW"   => device["connectionCapacity"],
        "dst_id" => d[!,"category"][1],
        "cluster_id"  => findall(x->x==1,[s_dict["$i"]["dst_id"]==d[!,"category"][1] for i=1:length(s_dict)])[1],
        "μ"  => μ,
        "σ"  => σ 

    )
    end;
end;

open(dir * branches_file_name,"r") do io
branches_json_dict = JSON.parse(io)
impedance_dict = Dict{String,Any}(
"BT - Desconocido BT" => [0.21, 0.075],
"BT - MANGUERA" => [0.3586, 0.089], #[1.23, 0.08],
"BT - RV 0,6/1 KV 2*16 KAL" => [2.14, 0.09], #16 = 20 A
"BT - RV 0,6/1 KV 2*25 KAL" => [1.34, 0.097],
"BT - RV 0,6/1 KV 3(1*150 KAL) + 1*95 KAL" => [0.2309, 0.085],
"BT - RV 0,6/1 KV 3(1*240 KAL) + 1*150 KAL" => [0.1602, 0.079],
"BT - RV 0,6/1 KV 3(1*240 KAL) + 1*95 KAL" => [0.1602, 0.079],
"BT - RV 0,6/1 KV 4*25 KAL" => [1.34, 0.097],
"BT - RV 0,6/1 KV 4*50 KAL" => [0.71849, 0.093],
"BT - RV 0,6/1 KV 4*95 KAL" => [0.3586, 0.089],
"BT - RX 0,6/1 KV 2*16 Cu" => [1.23, 0.08],
"BT - RX 0,6/1 KV 2*2 Cu" => [9.9, 0.075],
"BT - RX 0,6/1 KV 2*4 Cu" => [4.95, 0.075],
"BT - RX 0,6/1 KV 2*6 Cu" => [3.3, 0.075],
"BT - RZ 0,6/1 KV 2*16 AL" => [2.14, 0.09],
"BT - RZ 0,6/1 KV 3*150 AL/80 ALM" => [0.2309, 0.85],
"BT - RZ 0,6/1 KV 3*150 AL/95 ALM" => [0.2309, 0.85],
"BT - RZ 0,6/1 KV 3*25 AL/54,6 ALM" => [1.34, 0.097],
"BT - RZ 0,6/1 KV 3*35 AL/54,6 ALM" => [0.9073, 0.095],
"BT - RZ 0,6/1 KV 3*50 AL/54,6 ALM" => [0.718497, 0.093],
"BT - RZ 0,6/1 KV 3*70 ALM/54,6 AL" => [0.4539, 0.091],
"BT - RZ 0,6/1 KV 3*95 AL/54,6 ALM" => [0.3586, 0.089],
"BT - RZ 0,6/1 KV 4*16 AL" => [2.14, 0.09],
"aansluitkabel" => [1.15, 0.150]
)

open(dir * branches_file_name,"r") do io
branches_json_dict = JSON.parse(io)
currentmax_dict = Dict{String,Any}(
"BT - Desconocido BT" => 200,
"BT - MANGUERA" => 150, #200 certain  40.18#150
"BT - RV 0,6/1 KV 2*16 KAL" => 75,
"BT - RV 0,6/1 KV 2*25 KAL" => 100,
"BT - RV 0,6/1 KV 3(1*150 KAL) + 1*95 KAL" => 305,
"BT - RV 0,6/1 KV 3(1*240 KAL) + 1*150 KAL" => 344,
"BT - RV 0,6/1 KV 3(1*240 KAL) + 1*95 KAL" => 344,
"BT - RV 0,6/1 KV 4*25 KAL" => 100,
"BT - RV 0,6/1 KV 4*50 KAL" => 150,
"BT - RV 0,6/1 KV 4*95 KAL" => 230,
"BT - RX 0,6/1 KV 2*16 Cu" => 75,#95,
"BT - RX 0,6/1 KV 2*2 Cu" => 40, #30,
"BT - RX 0,6/1 KV 2*4 Cu" => 60,#40,
"BT - RX 0,6/1 KV 2*6 Cu" => 80, #50,
"BT - RZ 0,6/1 KV 2*16 AL" => 20, #75,
"BT - RZ 0,6/1 KV 3*150 AL/80 ALM" => 305, #264,
"BT - RZ 0,6/1 KV 3*150 AL/95 ALM" => 305,#264,
"BT - RZ 0,6/1 KV 3*25 AL/54,6 ALM" => 100, #78.98,
"BT - RZ 0,6/1 KV 3*35 AL/54,6 ALM" => 120,
"BT - RZ 0,6/1 KV 3*50 AL/54,6 ALM" => 150, #118.47,
"BT - RZ 0,6/1 KV 3*70 ALM/54,6 AL" => 160,
"BT - RZ 0,6/1 KV 3*95 AL/54,6 ALM" => 230, # 182.21,
"BT - RZ 0,6/1 KV 4*16 AL" => 75,
"aansluitkabel" => 120 #200
)


    for branch in branches_json_dict
        id = branch["branchId"] +1
        id_s = string(id)
        network_model["branch"][id_s] = Dict{String,Any}(
            "shift"         => 0.0,
            #"f_connections" => [1, 2, 3],
            "name"          => id_s,
            "switch"        => false,
            "g_to"          => 0.0,
            "c_rating_a"    => currentmax_dict[branch["cableType"]],
            "vbase"         =>  voltage_base,
            "g_fr"          => 0.0,
            #"t_connections" => [1, 2, 3],
            "f_bus"         => branch["upBusId"]+1,
            "b_fr"          => 0.0,
            "c_rating_b"    => currentmax_dict[branch["cableType"]],
            "br_status"     => 1,
            "t_bus"         => branch["downBusId"]+1,
            "b_to"          => 0.0,
            "index"         => id,
            "angmin"        => -1.0472,
            "angmax"        => 1.0472,
            "transformer"   => false,
            "tap"           => 1.0,
            "c_rating_c"    => currentmax_dict[branch["cableType"]], 
            "λcmax"         => 1.65 #1.65 #2.5 #1.03643    
            )

        if haskey(impedance_dict,branch["cableType"])
            # network_model["branch"][id_s]["br_x"] = (cross_area_fact^(1/4)).* impedance_dict[branch["cableType"]][2] .* (branch["cableLength"]/1000+1E-6)./  Z_base 
            # network_model["branch"][id_s]["br_r"] = (1/cross_area_fact).* impedance_dict[branch["cableType"]][1] .* (branch["cableLength"]/1000+1E-6)./  Z_base
            network_model["branch"][id_s]["br_x"] = (1/cross_area_fact).* impedance_dict[branch["cableType"]][2] .* (branch["cableLength"]/1000+1E-6)./  Z_base 
            network_model["branch"][id_s]["br_r"] = (1/cross_area_fact).* impedance_dict[branch["cableType"]][1] .* (branch["cableLength"]/1000+1E-6)./  Z_base  
        end;
        

        if haskey(currentmax_dict,branch["cableType"])
            network_model["branch"][id_s]["rate_a"] = cross_area_fact.*((currentmax_dict[branch["cableType"]]*voltage_base)/1e3)/power_base
            
            #network_model["branch"][id_s]["I_rating"] = currentmax_dict[branch["cableType"]]/current_base

        end;
    end;
end;
end;


network_model["sdata"]= s_dict 
network_model["curt"]= curt

network_model["PV"]=deepcopy(network_model["load"]);
[network_model["PV"][d]["μ"]=s_dict[string(length(s_dict))]["pc"] for d in   keys(network_model["PV"])]
[network_model["PV"][d]["σ"]=s_dict[string(length(s_dict))]["pd"] for d in   keys(network_model["PV"])]
[network_model["PV"][d]["pd"]=s_dict[string(length(s_dict))]["pd"]/1e6/ power_base / 3 for d in   keys(network_model["PV"])]
[network_model["PV"][d]["p_min"]=0 for d in   keys(network_model["PV"])] 
if grid_only == true
    for d in keys(network_model["PV"])
        print(solar_cost[network_model["PV"][d]["source_id"]])
        if solar_cost[network_model["PV"][d]["source_id"]] == 0
            network_model["PV"][d]["p_max"] =0
            print("0")
        else
            network_model["PV"][d]["p_max"] = 15
            print("15")
        end
    end
else 
    # [network_model["PV"][d]["p_max"]=PVins.installed_PVcap_kW_total[parse.(Int,d)] for d in   keys(network_model["PV"])] PVsize[data["PV"]["1"]["source_id"]]
    [network_model["PV"][d]["p_max"] = PVsize[network_model["PV"][d]["source_id"]] for d in   keys(network_model["PV"])]
end
return network_model
end;


function naive_time_series_analysis(mn_model, time_steps)
    solver = Ipopt.Optimizer
    len_nodes=length(mn_model["nw"]["1"]["bus"])
    result_pf = Dict{String,Any}()
    power_1=zeros(time_steps,len_nodes-1)
    voltage_5= zeros(time_steps,len_nodes)
    current_1=zeros(time_steps,len_nodes-1)

    s1 = Dict("output" => Dict("branch_flows" => true))
    s2 = Dict("output" => Dict("duals" => true)) 
    #branc_details = PowerModels.component_table(network_data, "branch", ["f_bus", "t_bus", "rate_a","br_status","br_r", "br_x"])
    branc_details = PowerModels.component_table(mn_model["nw"]["1"], "branch", ["f_bus", "t_bus", "rate_a","br_status","br_r", "br_x"])
    for (n,network) in mn_model["nw"]
        network["per_unit"] = true
        a=parse(Int,n)
        #print(n)
        # result[n] = run_ac_opf(network, solver; setting = s2);
        # d0=PowerModels.component_table(result[n]["solution"], "bus", ["lam_kcl_i", "lam_kcl_r"]) ;
       
        result_pf[n] = PowerModels.run_pf(network, _PM.IVRPowerModel, solver; setting = s1);
        m0 = PowerModels.component_table(result_pf[n]["solution"], "bus", ["vi", "vr"])
        voltage_5[a,:]= sqrt.(m0[:,2].^2+m0[:,3].^2)
        # current_1[a,:]=v[1,data["branch"]["2"]["t_bus"]]-v[1,data["branch"]["2"]["f_bus"]]/sqrt.(data["branch"]["2"]["br_r"]^2+data["branch"]["2"]["br_x"]^2)
        # #power_1[n]=result[n]["solution"]["gen"]["1"]["pg"];
        
        # Injection_ac_orig = PowerModels.component_table(result_pf[n]["solution"], "branch", ["qf", "pf", "qt","pt"])
        # thermal_ac_original = max(sqrt.(Injection_ac_orig[:,2].^2 + Injection_ac_orig[:,3].^2), sqrt.(Injection_ac_orig[:,4].^2 + Injection_ac_orig[:,5].^2))
        # loss_1[a,:] = 100 * abs.(Injection_ac_orig[:,3] + Injection_ac_orig[:,5]) ./ branc_details[:,4]
        # power_1[a,:] = 100 * thermal_ac_original ./ branc_details[:,4] 
    end;
    return result_pf, voltage_5
end;


function time_series_hc(mn_model, model, solver; aux=true, deg=2, red=false, stochastic=false, time=time)
    result = Dict{String,Any}()
    len_pv=length(mn_model["nw"]["1"]["PV"])
    pv_size= zeros(length(mn_model["nw"]),len_pv)
    for (n,network) in mn_model["nw"]
        a=parse(Int,n) 
        if a in time
            network["per_unit"] = true
            display(n)
            result[n] = _PM.run_model(network, model, solver, _SPM.build_sopf_hc_deterministic; multinetwork=false)
            pv_size[a,:]=[result[n]["solution"]["PV"]["$j"]["p_size"] for j=1:len_pv]
        end
    end
    return result,pv_size

end