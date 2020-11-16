module DPHeuristics

using JuMP
using Gurobi
using CPLEX
using Data
using Parameters
using Random
#using FixAndOptimize

mutable struct stdFormVars
	xp
	xr
	xw
	yp
	yr
	yw
	sp
	sr
	sw
end

export threelevelFormulation, stdFormVars

function RandomizedDPHeuristicBottomUp(inst,params)
	#G= zeros(Float64,inst.NR,inst.NT)
	#GK= zeros(Int,inst.NR,inst.NT)

	bestSETR = Vector{Vector{Int}}(undef,inst.NR)
	bestSETW = Vector{Vector{Int}}(undef,inst.NW)
	bestSETP = Vector{Int}

	t1 = time_ns()

	bestcost = 999999999.0
	alpha = params.dpalpha
	rng = MersenneTwister(params.dpseed);

	greedycost = 0.0

	randSCR =zeros(Int,inst.NR,inst.NT)
	randSCW = zeros(Int,inst.NW,inst.NT)

	sumcosts = 0.0

	for iter in 1:params.dpmaxiter

		if iter != 1

			if params.dptype == 1
				randValuesSCR = rand!(rng,zeros(inst.NR,inst.NT))*params.dpalpha
				#println("randValuesSCR = ",randValuesSCR)
				randValuesSCW = rand!(rng,zeros(inst.NW,inst.NT))*params.dpalpha
				#println("randValuesSCW = ",randValuesSCW)
			elseif params.dptype == 2
				randValuesSCR = rand!(rng,zeros(inst.NR,inst.NT))*2*params.dpalpha
				for r in 1:inst.NR, t in 1:inst.NT
					randValuesSCR[r,t] -= params.dpalpha
				end
				#println("randValuesSCR = ",randValuesSCR)
				randValuesSCW = rand!(rng,zeros(inst.NW,inst.NT))*2*params.dpalpha
				for w in 1:inst.NW, t in 1:inst.NT
					randValuesSCW[w,t] -= params.dpalpha
				end
				#println("randValuesSCW = ",randValuesSCW)
			end

			for r in 1:inst.NR,t in 1:inst.NT
				randSCR[r,t] = inst.SCR[r,t] + round(inst.SCR[r,t]*(randValuesSCR[r,t]))
			end
			for w in 1:inst.NW,t in 1:inst.NT
				randSCW[w,t] = inst.SCW[w,t] + round(inst.SCW[w,t]*(randValuesSCW[w,t]))
			end


		else
			for r in 1:inst.NR,t in 1:inst.NT
				randSCR[r,t] = inst.SCR[r,t]
			end
			for w in 1:inst.NW,t in 1:inst.NT
				randSCW[w,t] = inst.SCW[w,t]
			end
		end

		G,GK = DPRetailer(inst.HCR,randSCR,inst.D,inst.NR,inst.NT)
		newDW = calculateDW(inst,GK)
		GW, GWK = DPWarehouse(inst.HCW,randSCW,newDW,inst.NW,inst.NT)

		newDP = calculateDP(inst,GWK,newDW)
		GP,GPK = DPPlant(inst.HCP,inst.SCP,newDP,inst.NT)

		solcost = sum(G[:,inst.NT]) + sum(GW[:,inst.NT]) + GP[inst.NT]
		#println("solcost = ",solcost)



		SETR,SETW,SETP,calcsolcost = RandomizedDPrecoversolution(inst,GK,GWK,GPK,newDW,newDP)

		sumcosts += calcsolcost

		if calcsolcost < bestcost
			bestcost = calcsolcost
			bestSETR = SETR
			bestSETW =SETW
			bestSETP = SETP
		end
		if iter == 1
			greedycost = calcsolcost
		end

		#println("solcostcalculated = ",calcsolcost)
		#println("\n\n\n")

	end

	t2 = time_ns()
	elapsedtime = (t2-t1)/1.0e9
	#println("Elapsed = ", elapsedtime)
	#println("Greedycost = ",greedycost)
	#println("Bestcost = ",bestcost)

	averagecost = sumcosts/params.dpmaxiter

	if params.form == "randdpheur"
		open("saida.txt","a") do f
			write(f,";$bestcost;$averagecost;$elapsedtime;$(params.dptype);$(params.dpalpha) \n")
		end
	elseif params.form == "randdpheurmc"

		#println("SETR = ",SETR)
		#println("SETW = ",SETW)
		#println("SETP = ",SETP)

		optvalue = partialmulticommodityRounding(inst,params,bestSETR,bestSETW,bestSETP,bestcost)

	end

	#println("SETR = ",SETR)
	#println("SETW = ",SETW)
	#println("SETP = ",SETP)

	#optvalue = partialmulticommodityRounding(inst,params,bestSETR,bestSETW,bestSETP)

	#FixAndOptimize.FixAndOptimizeStandardFormulation(inst,params,bestSETR,bestSETW,bestSETP,bestcost)

	#gap = 100*(bestcost-optvalue)/bestcost
	#println("Heuristic gap = ",gap)

	return bestSETR,bestSETW,bestSETP

end


function DPHeuristicBottomUp(inst,params)
	#G= zeros(Float64,inst.NR,inst.NT) #
	#GK= zeros(Int,inst.NR,inst.NT)

	t1 = time_ns()

	G,GK = DPRetailer(inst.HCR,inst.SCR,inst.D,inst.NR,inst.NT)
	newDW = calculateDW(inst,GK)
	GW, GWK = DPWarehouse(inst.HCW,inst.SCW,newDW,inst.NW,inst.NT)

	newDP = calculateDP(inst,GWK,newDW)
	GP,GPK = DPPlant(inst.HCP,inst.SCP,newDP,inst.NT)

	solcost = sum(G[:,inst.NT]) + sum(GW[:,inst.NT]) + GP[inst.NT]
	#println("solcost = ",solcost)

	t2 = time_ns()
	elapsedtime = (t2-t1)/1.0e9
	#println("Elapsed = ", elapsedtime)

	SETR,SETW,SETP = DPrecoversolution(inst,GK,GWK,GPK)

	#println("SETR = ",SETR)
	#println("SETW = ",SETW)
	#println("SETP = ",SETP)

	#optvalue = partialmulticommodityRounding(inst,params,SETR,SETW,SETP)

	#gap = 100*(solcost-optvalue)/solcost
	#println("Heuristic gap = ",gap)

end


function DPrecoversolution(inst,GK,GWK,GPK)

	SETR = Vector{Vector{Int}}(undef,inst.NR)
	for r in 1:inst.NR
		SETR[r] = []
	end

	SETW = Vector{Vector{Int}}(undef,inst.NW)
	for w in 1:inst.NW
		SETW[w] = []
	end

	SETP = Vector{Int}
	SETP = []

	for r in 1:inst.NR
		k = GK[r,inst.NT]
		push!(SETR[r],k)
		while k > 1
			oldk = k
			k = GK[r,k-1]
			push!(SETR[r],k)
		end
	end

	for w in 1:inst.NW
		k = GWK[w,inst.NT]
		push!(SETW[w],k)
		while k > 1
			oldk = k
			k = GWK[w,k-1]
			push!(SETW[w],k)
		end
	end

	k = GPK[inst.NT]
	push!(SETP,k)
	while k > 1
		oldk = k
		k = GPK[k-1]
		push!(SETP,k)
	end

	return SETR,SETW,SETP

end #

function RandomizedDPrecoversolution(inst,GK,GWK,GPK,newDW,newDP)

	DPHCP = preparePlantAggregatedHoldingCosts(inst.HCP,newDP,inst.NT)
	DPHCW = prepareWarehouseAggregatedHoldingCosts(inst.HCW,newDW,inst.NW,inst.NT)
	DPHCR = prepareRetailerAggregatedHoldingCosts(inst.HCR,inst.D,inst.NR,inst.NT)

	SETR = Vector{Vector{Int}}(undef,inst.NR)
	for r in 1:inst.NR
		SETR[r] = []
	end

	SETW = Vector{Vector{Int}}(undef,inst.NW)
	for w in 1:inst.NW
		SETW[w] = []
	end

	SETP = Vector{Int}
	SETP = []

	solcost = 0

	for r in 1:inst.NR
		k = GK[r,inst.NT]
		solcost += inst.SCR[r,k] + DPHCR[r,k,inst.NT]
		push!(SETR[r],k)
		while k > 1
			oldk = k
			k = GK[r,k-1]
			solcost += inst.SCR[r,k] + DPHCR[r,k,oldk-1]
			push!(SETR[r],k)
		end
	end

	for w in 1:inst.NW
		k = GWK[w,inst.NT]
		solcost += inst.SCW[w,k] + DPHCW[w,k,inst.NT]
		push!(SETW[w],k)
		while k > 1
			oldk = k
			k = GWK[w,k-1]
			solcost += inst.SCW[w,k] + DPHCW[w,k,oldk-1]
			push!(SETW[w],k)
		end
	end

	k = GPK[inst.NT]
	solcost += inst.SCP[1,k] + DPHCP[k,inst.NT]
	push!(SETP,k)
	while k > 1
		oldk = k
		k = GPK[k-1]
		solcost += inst.SCP[1,k] + DPHCP[k,oldk-1]
		push!(SETP,k)
	end

	return SETR,SETW,SETP,solcost

end #


function calculateDP(inst,GWK,newDW)

	newDP = zeros(Int,inst.NT)

	#println("Checking plant")
	for w in 1:inst.NW
	#	println("\t Checking warehouse = ",w)
	#	println("\t\t GK = ",GWK[w,:])
	#	println("\t\t D = ",newDW[w,:])
		k = GWK[w,inst.NT]
		demk = sum(newDW[w,k:inst.NT])
	#	println("\t\t k = ",k,",     demk = ",demk)
		newDP[k] += demk
		while k > 1
			oldk = k
			k = GWK[w,k-1]
			#if k>= 1
				demk = sum(newDW[w,k:oldk-1])
	#			println("\t\t k = ",k,",     demk = ",demk)
				newDP[k] += demk
			#end
		end
	#	println("finished")
	end

	return newDP

end #calculateDP


function calculateDW(inst,GK)

	newDW = zeros(Int,inst.NW,inst.NT)

	for w in 1:inst.NW
		#println("Checking warehouse = ",w)
		for rw in 1:length(inst.DeltaW[w])
			r = inst.DeltaW[w][rw]
			#println("\t Checking retailer = ",r)
			#println("\t\t GK = ",GK[r,:])
			#println("\t\t D = ",inst.D[r,:])
			k = GK[r,inst.NT]
			demk = sum(inst.D[r,k:inst.NT])
			#println("\t\t k = ",k,",     demk = ",demk)
			newDW[w,k] += demk
			while k > 1
				oldk = k
				k = GK[r,k-1]
				#if k>= 1
					demk = sum(inst.D[r,k:oldk-1])
					#println("\t\t k = ",k,",     demk = ",demk)
					newDW[w,k] += demk
				#end
			end
			#println("finished")
		end

	end

	return newDW

end #calculateDW


function DPPlant(HCP,SCP,newDP,NT)

	DPHCP = preparePlantAggregatedHoldingCosts(HCP,newDP,NT)

	GP= zeros(Float64,NT)
	GPK= zeros(Int,NT)

	GP[1] = SCP[1,1] + DPHCP[1,1]
	GPK[1] = 1
	for t in 2:NT
		minvalue =  SCP[1,1] + DPHCP[1,t]
		minindex = 1
		for k in 2:t
			kvalue = GP[k-1] + SCP[1,k] + DPHCP[k,t]
			if kvalue < minvalue
				minvalue = kvalue
				minindex = k
			end
		end
		GP[t] = minvalue
		GPK[t] = minindex
	end

	totalcost = sum(GP[NT])
	#println("Totalcost = ",totalcost)

	return GP, GPK

end # DPPlant


function DPWarehouse(HCW,SCW,newDW,NW,NT)

	DPHCW = prepareWarehouseAggregatedHoldingCosts(HCW,newDW,NW,NT)

	GW= zeros(Float64,NW,NT) #
	GWK= zeros(Int,NW,NT)
	for w in 1:NW
		GW[w,1] = SCW[w,1] + DPHCW[w,1,1]
		GWK[w,1] = 1
		for t in 2:NT
			minvalue =  SCW[w,1] + DPHCW[w,1,t]
			minindex = 1
			for k in 2:t
				kvalue = GW[w,k-1] + SCW[w,k] + DPHCW[w,k,t]
				if kvalue < minvalue
					minvalue = kvalue
					minindex = k
				end
			end
			GW[w,t] = minvalue
			GWK[w,t] = minindex
		end
	end

	totalcost = sum(GW[:,NT])
	#println("Totalcost = ",totalcost)

	return GW, GWK

end # DPRetailer


function DPRetailer(HCR,SCR,D,NR,NT)

	DPHCR = prepareRetailerAggregatedHoldingCosts(HCR,D,NR,NT)

	G= zeros(Float64,NR,NT) #
	GK= zeros(Int,NR,NT)
	for r in 1:NR
		G[r,1] = SCR[r,1] + DPHCR[r,1,1]
		GK[r,1] = 1
		for t in 2:NT
			minvalue =  SCR[r,1] + DPHCR[r,1,t]
			minindex = 1
			for k in 2:t
				kvalue = G[r,k-1] + SCR[r,k] + DPHCR[r,k,t]
				if kvalue < minvalue
					minvalue = kvalue
					minindex = k
				end
			end
			G[r,t] = minvalue
			GK[r,t] = minindex
		end
	end

	totalcost = sum(G[:,NT])
	#println("Totalcost = ",totalcost)

	return G, GK

end # DPRetailer


function prepareRetailerAggregatedHoldingCosts(HCR,D,NR,NT)

	DPHCR = zeros(Float64,NR,NT,NT)
	for r in 1:NR
		for t in 1:NT
			DPHCR[r,t,t] = 0
			for k in t+1:NT
				DPHCR[r,t,k] = DPHCR[r,t,k-1] + (k-t)*HCR[r]*D[r,k]
			end
		end
	end
	return DPHCR

end #prepareRetailerAggregatedHoldingCosts


function prepareWarehouseAggregatedHoldingCosts(HCW,newDW,NW,NT)

	DPHCW = zeros(Float64,NW,NT,NT)
	for w in 1:NW
		for t in 1:NT
			DPHCW[w,t,t] = 0
			for k in t+1:NT
				DPHCW[w,t,k] = DPHCW[w,t,k-1] + (k-t)*HCW[w]*newDW[w,k]
			end
		end
	end
	return DPHCW

end #prepareRetailerAggregatedHoldingCosts


function preparePlantAggregatedHoldingCosts(HCP,newDP,NT)

	DPHCP = zeros(Float64,NT,NT)

	for t in 1:NT
		DPHCP[t,t] = 0
		for k in t+1:NT
			DPHCP[t,k] = DPHCP[t,k-1] + (k-t)*HCP[1]*newDP[k]
		end
	end

	return DPHCP

end #preparePlantAggregatedHoldingCosts


function DPHeuristicBottomUpCC(inst,params)
	#G= zeros(Float64,inst.NR,inst.NT) #
	#GK= zeros(Int,inst.NR,inst.NT)

	t1 = time_ns()

	G,GK = DPRetailer(inst.HCR,inst.SCR,inst.D,inst.NR,inst.NT)
	newDW = calculateDW(inst,GK)
	GW, GWK = DPWarehouse(inst.HCW,inst.SCW,newDW,inst.NW,inst.NT)

	newDP = calculateDP(inst,GWK,newDW)
	GP,GPK = DPPlant(inst.HCP,inst.SCP,newDP,inst.NT)

	GP,GPK = DPPlantCC(inst.HCP,inst.SCP,newDP,inst.NT,inst.C[1])

	#solcost = sum(G[:,inst.NT]) + sum(GW[:,inst.NT]) + GP[inst.NT]
	#println("solcost = ",solcost)

	t2 = time_ns()
	elapsedtime = (t2-t1)/1.0e9
	println("Elapsed = ", elapsedtime)


	#SETR,SETW,SETP = DPrecoversolution(inst,GK,GWK,GPK)

	#println("SETR = ",SETR)
	#println("SETW = ",SETW)
	#println("SETP = ",SETP)

	#optvalue = partialmulticommodityRounding(inst,params,SETR,SETW,SETP)

	#gap = 100*(solcost-optvalue)/solcost
	#println("Heuristic gap = ",gap)

end


function DPPlantCC(HCP,SCP,newDP,NT,C)
 	#solve constant capacity lot-sizing for the plant
	#println("Running DPPlantCC")

	GP= zeros(Float64,NT) #
	GPK= zeros(Int,NT)

	regenerationIntervalCost(HCP,SCP,newDP,1,NT,C)

	#for a in 1:inst.NT, b in a:inst.NT
	#	solveproblemQ(HCP,SCP,newDP,NT,C,j)

	#end

	return GP,GPK

end # DPPlant

function regenerationIntervalCost(HCP,SCP,newDP,a,b,C)
	#println("Running regenerationIntervalCost [$(a),$(b)]")

	deltat = zeros(Int,b-a+1)
	for t in a:b
		deltat[t-a+1] = t - ceil( sum(newDP[a:a+t-1])/C )
	end
	println("delta = ",deltat)

end #regenerationIntervalCost


function solveproblemQ(HCP,SCP,newDP,NT,C,j)

end #solveproblemQ

function partialmulticommodityRounding(inst::InstanceData, params::ParameterData,SETR,SETW,SETP,dpheurvalue)
	#println("Running NewFormulations.partialmulticommodityFormulation")
	if params.solver == "Gurobi"
		if params.disablesolver == 1 #Disable gurobi cuts and presolve
		   if params.maxnodes < 999.0
			   model = Model(solver = GurobiSolver(TimeLimit=params.maxtime,MIPGap=params.tolgap,CliqueCuts=0, CoverCuts=0, FlowCoverCuts=0,FlowPathCuts=0,MIRCuts=0,NetworkCuts=0,GomoryPasses=0, PreCrush=1,NodeLimit=params.maxnodes,Method=2))
		   else
			   model = Model(solver = GurobiSolver(TimeLimit=params.maxtime,MIPGap=params.tolgap,CliqueCuts=0, CoverCuts=0, FlowCoverCuts=0,FlowPathCuts=0,MIRCuts=0,NetworkCuts=0,GomoryPasses=0, PreCrush=1,Method=2))
		   end
	   else
		   if params.maxnodes < 999.0
			   model = Model(solver = GurobiSolver(TimeLimit=params.maxtime,MIPGap=params.tolgap,PreCrush=1,NodeLimit=params.maxnodes,Method=2))
		   else
			   model = Model(solver = GurobiSolver(TimeLimit=params.maxtime,MIPGap=params.tolgap,PreCrush=1,Method=2))
		   end
	   end
	elseif params.solver == "Cplex"
        env = Cplex.Env()
		 model = Model(solver = CplexSolver(CPX_PARAM_TILIM=params.maxtime,CPX_PARAM_EPGAP=params.tolgap))
	else
		println("No solver selected")
		return 0
	end

	### Defining variables ###
	@variable(model,0 <= xxp[r=1:inst.NR,t=1:inst.NT,k=t:inst.NT] <= Inf) #sum(inst.DP[i,t] for t=1:inst.NT))
	@variable(model,0 <= xxr[r=1:inst.NR,t=1:inst.NT,k=t:inst.NT] <= Inf) #sum(inst.DR[i,t] for t=1:inst.NT))
	@variable(model,0 <= xxw[r=1:inst.NR,t=1:inst.NT,k=t:inst.NT] <= Inf) #sum(inst.D[i,t] for t=1:inst.NT))
	@variable(model, yp[p=1:inst.NP,t=1:inst.NT], Bin)
	@variable(model, yr[r=1:inst.NR,t=1:inst.NT], Bin)
	@variable(model, yw[w=1:inst.NW,t=1:inst.NT], Bin)
	@variable(model,0 <= ssp[r=1:inst.NR,t=0:inst.NT-1,k=t+1:inst.NT] <= Inf) #sum(inst.DP[i,t] for t=1:inst.NT))
	@variable(model,0 <= ssr[r=1:inst.NR,t=0:inst.NT-1,k=t+1:inst.NT] <= Inf) #sum(inst.DR[i,t] for t=1:inst.NT))
	@variable(model,0 <= ssw[r=1:inst.NR,t=0:inst.NT-1,k=t+1:inst.NT] <= Inf)

	ypsol = zeros(inst.NP,inst.NT)
	ywsol = zeros(inst.NW,inst.NT)
	yrsol = zeros(inst.NR,inst.NT)

	for t in 1:length(SETP)
		ypsol[1,SETP[t]] = 1
	end
	for w in 1:inst.NW, t in 1:length(SETW[w])
		ywsol[w,SETW[w][t]] = 1
	end
	for r in 1:inst.NR, t in 1:length(SETR[r])
		yrsol[r,SETR[r][t]] = 1
	end

	for t = 1:inst.NT
		JuMP.setvalue(yp[1,t],ypsol[1,t])
	end
	for w=1:inst.NW, t=1:inst.NT
		JuMP.setvalue(yw[w,t], ywsol[w,t])
	end
	for r=1:inst.NR, t=1:inst.NT
		JuMP.setvalue(yr[r,t],yrsol[r,t])
	end

	#if params.reducefl != 0
		useX = ones(Int,inst.NR, inst.NT, inst.NT)
		for r in 1:inst.NR, t in 1:inst.NT-1, k in t+1:inst.NT
		#for r in 6:6, t in 10:10, k in t+1:t+2 #inst.NT
			if (useX[r,t,k-1] == 0) || ( inst.D[r,k]*(k-t)*inst.HCR[r] >= inst.D[r,k]*(k-t)*inst.HCW[inst.DeltamR[r]] + inst.SCR[r,k])
			#if  ( inst.D[r,k]*(k-t)*inst.HCR[r] > inst.D[r,k]*(k-t)*inst.HCW[inst.DeltamR[r]] + inst.SCR[r])
				#println("r = $(r), [$(t),$(k)]")
				#println(inst.D[r,k]*(k-t)*inst.HCR[r],"  >=  ",inst.D[r,k]*(k-t)*inst.HCW[inst.DeltamR[r]] + inst.SCR[r,k])
				#println("dem $(inst.D[r,k]), hcr $(inst.HCR[r]), hcw $(inst.HCW[inst.DeltamR[r]]), scr $(inst.SCR[r,k])\n")
				useX[r,t,k] = 0
			end
		end
	#end
	#println("useX = ",useX)
	totalVariables = 0
	totalFixed = 0
	for r in 1:inst.NR, t in 1:inst.NT-1, k in t+1:inst.NT
		totalVariables += 1
		if useX[r,t,k] == 0
			totalFixed +=1
			#JuMP.fix(xxr[r,t,k],0)
			@constraint(model,xxr[r,t,k] == 0)
		end
	end
	percentagefixed = 100*(totalFixed/totalVariables)
	#println("Percentage of fixed = ",percentagefixed)

	#println("Defined variables")
	### Objective function ###
	@objective(model, Min,
		sum(inst.SCP[p,t]*yp[p,t] for p=1:inst.NP, t=1:inst.NT) + sum(inst.HCP[1]*ssp[r,k,t] for r=1:inst.NR,k=1:inst.NT-1,t= k+1:inst.NT)
		+ sum(inst.SCW[w,t]*yw[w,t]  for w=1:inst.NW, t=1:inst.NT) + sum(inst.HCW[inst.DeltamR[r]]*ssw[r,k,t] for r=1:inst.NR,k=1:inst.NT-1,t= k+1:inst.NT)
		+ sum(inst.SCR[r,t]*yr[r,t]  for r=1:inst.NR, t=1:inst.NT) + sum(inst.HCR[r]*ssr[r,k,t] for r=1:inst.NR,k=1:inst.NT-1,t= k+1:inst.NT)
	)

	#println("Finished objective function")

	### Setup constraints ###
	@constraint(model,
				setupP[r=1:inst.NR,t=1:inst.NT,k=t:inst.NT], xxp[r,t,k] <= inst.D[r,k]*yp[1,t]
				)
	@constraint(model,
				setupW[r=1:inst.NR, t=1:inst.NT,k=t:inst.NT], xxw[r,t,k] <= inst.D[r,k]*yw[inst.DeltamR[r],t]
				)
	@constraint(model,
				setupR[r=1:inst.NR,t=1:inst.NT,k=t:inst.NT], xxr[r,t,k] <= inst.D[r,k]*yr[r,t]
				)

	if params.capacity != 0
		@constraint(model,
					capP[t=1:inst.NT], sum(xxp[r,t,k] for r=1:inst.NR,k=t:inst.NT) <= inst.C[t]*yp[1,t]
					)
	end

	#println("Finished setup constraints")

	### no initial inventory
	@constraint(model,
				zeroinvP[r=1:inst.NR,t=1:inst.NT],
				ssp[r,0,t] == 0
				)
	@constraint(model,
				zeroinvW[r=1:inst.NR,t=1:inst.NT],
				ssw[r,0,t] == 0
				)
	@constraint(model,
				zeroinvR[r=1:inst.NR,t=1:inst.NT],
				ssr[r,0,t] == 0
				)

	#println("Finished zero inventory constraints")

	### Balance constraints ###
	@constraint(model,
				balanceP[r=1:inst.NR, t=1:inst.NT, k=t:inst.NT],
				ssp[r,t-1,k] + xxp[r,t,k] == xxw[r,t,k] + sum(ssp[r,t,z] for z in k:k if t<z)
				)
				#println("h")
	@constraint(model,
				balanceW[r=1:inst.NR,t=1:inst.NT,k=t:inst.NT],
				ssw[r,t-1,k] + xxw[r,t,k] == xxr[r,t,k] + sum(ssw[r,t,z] for z in k:k if t<z)
				)
				#println("h1")
	@constraint(model,
				balanceR[r=1:inst.NR,t=1:inst.NT,k=t:inst.NT],
				ssr[r,t-1,k] + xxr[r,t,k] == sum(inst.D[r,z] for z in k:k if t==k) + sum(ssr[r,t,z] for z in k:k if t<k)
				)

	#writeLP(model,"modelo.lp",genericnames=false)

	t1 = time_ns()
	solve(model)

	#bestsol = getobjectivevalue(model)
	#bestbound = getobjbound(model)

	t2 = time_ns()
	elapsedtime = (t2-t1)/1.0e9

	bestsol = getobjectivevalue(model)
	bestbound = getobjbound(model)
	numnodes = getnodecount(model)
	time = getsolvetime(model)
	gap = 100*(bestsol-bestbound)/bestsol
	#println("bestsol = ", bestsol)
	#println("bestbound = ", bestbound)
    #println("gap = ", gap)
    #println("time = ", time)
    #println("nodes = ", numnodes)

	open("saida.txt","a") do f
		write(f,"$(params.form);$dpheurvalue;$bestbound;$bestsol;$gap;$time;$numnodes;$percentagefixed \n")
	end

	return bestsol

	#if params.printsol == 1
	#	printStandardFormulationSolution(inst,xp,xr,xw,yp,yr,yw,sp,sr,sw)
	#end

end #function multicommodityFormulation()

end
