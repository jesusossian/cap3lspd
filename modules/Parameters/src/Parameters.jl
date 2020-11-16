module Parameters

struct ParameterData
	instName::String
	form::String
	balanced::Int
	solver::String
	maxtime::Int
	tolgap::Float64
	printsol::Int
	capacity::Float64
	disablesolver::Int
	dpseed::Int
	dpalpha::Float64
	dptype::Int
	dpmaxiter::Int
	tolsep::Float64
	cpmaxrounds::Int
	maxnodes::Int
	dpredperiods::Float64
	dpcostreductionret::Float64
	dpcostreductionware::Float64
end

export ParameterData, readInputParameters

function readInputParameters(ARGS)

    println("Running Parameters.readInputParameters")

    ### Set standard values for the parameters ###
	instName="instances/N50T15/N50T15P1W10SD_DF1.dat"
	form="std"
	balanced = 1
	solver = "Gurobi"
	maxtime = 3600
	tolgap = 0.000001
	printsol = 0
	capacity = 1.5
	disablesolver = 0
	dpseed = 1234
	dpalpha = 0.20
	dptype = 1
	dpmaxiter = 500
	tolsep = 10.0
	cpmaxrounds = 1000
	maxnodes = 10000.0

	dpredperiods = 0.7
	dpcostreductionret = 0.7
	dpcostreductionware = 0.85

    ### Read the parameters and set correct values whenever provided ###
    for param in 1:length(ARGS)
	    if ARGS[param] == "--inst"
            instName = ARGS[param+1]
            param += 1
        elseif ARGS[param] == "--solver"
            solver = ARGS[param+1]
            param += 1
        elseif ARGS[param] == "--maxtime"
            maxtime = parse(Int,ARGS[param+1])
            param += 1
        elseif ARGS[param] == "--tolgap"
            tolgap = parse(Float64,ARGS[param+1])
            param += 1
		elseif ARGS[param] == "--capacity"
			capacity = parse(Float64,ARGS[param+1])
			param += 1
        elseif ARGS[param] == "--printsol"
            printsol = parse(Int,ARGS[param+1])
            param += 1
		elseif ARGS[param] == "--balanced"
            balanced = parse(Int,ARGS[param+1])
            param += 1
		elseif ARGS[param] == "--disablesolver"
            disablesolver = parse(Int,ARGS[param+1])
            param += 1
        elseif ARGS[param] == "--form"
            form = ARGS[param+1]
            param += 1
		elseif ARGS[param] == "--dpseed"
            dpseed = parse(Int,ARGS[param+1])
            param += 1
		elseif ARGS[param] == "--dpalpha"
            dpalpha = parse(Float64,ARGS[param+1])
            param += 1
		elseif ARGS[param] == "--dptype"
            dptype = parse(Int,ARGS[param+1])
            param += 1
		elseif ARGS[param] == "--dpmaxiter"
			dpmaxiter = parse(Int,ARGS[param+1])
			param += 1
		elseif ARGS[param] == "--tolsep"
            tolsep = parse(Float64,ARGS[param+1])
            param += 1
		elseif ARGS[param] == "--cpmaxrounds"
			cpmaxrounds = parse(Int,ARGS[param+1])
			param += 1
		elseif ARGS[param] == "--maxnodes"
			maxnodes = parse(Float64,ARGS[param+1])
			param += 1
		elseif ARGS[param] == "--dpredperiods"
			dpredperiods = parse(Float64,ARGS[param+1])
			param += 1
		elseif ARGS[param] == "--dpcostreductionret"
			dpcostreductionret = parse(Float64,ARGS[param+1])
			param += 1
		elseif ARGS[param] == "--dpcostreductionware"
			dpcostreductionware = parse(Float64,ARGS[param+1])
			param += 1
        end
    end

    params = ParameterData(instName,form,balanced,solver,maxtime,tolgap,printsol,capacity,disablesolver,dpseed,dpalpha,dptype,dpmaxiter,tolsep,cpmaxrounds,maxnodes,dpredperiods,dpcostreductionret,dpcostreductionware)

    return params

end ### end readInputParameters

end ### end module
