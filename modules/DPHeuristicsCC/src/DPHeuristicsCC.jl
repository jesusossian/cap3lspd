module DPHeuristicsCC

using JuMP
using Gurobi
using CPLEX
using Data
using Parameters
using Random
using FixAndOptimize

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
	println("Running DPHeuristicsCC.RandomizedDPHeuristicBottomUp")
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





	FixAndOptimize.FixAndOptimizeStandardFormulation(inst,params,bestSETR,bestSETW,bestSETP,bestcost)

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


	plantcost = 0
	k = GPK[inst.NT]
	solcost += inst.SCP[1,k] + DPHCP[k,inst.NT]
	plantcost += inst.SCP[1,k] + DPHCP[k,inst.NT]
	push!(SETP,k)
	while k > 1
		oldk = k
		k = GPK[k-1]
		solcost += inst.SCP[1,k] + DPHCP[k,oldk-1]
		plantcost += inst.SCP[1,k] + DPHCP[k,oldk-1]
		push!(SETP,k)
	end

	println("DP plant cost = $(plantcost)")

	return SETR,SETW,SETP,solcost

end #



function RandomizedDPCCrecoversolution(inst,GK,GWK,GPK,newDW,newDP)

	#println("C: ",inst.C)

	#println("newDP: ",newDP)

	#println("hcp: ",inst.HCP)

	#println("scp: ",inst.SCP)


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

	retcost = 0
	for r in 1:inst.NR
		k = GK[r,inst.NT]
		solcost += inst.SCR[r,k] + DPHCR[r,k,inst.NT]
		retcost += inst.SCR[r,k] + DPHCR[r,k,inst.NT]
		push!(SETR[r],k)
		while k > 1
			oldk = k
			k = GK[r,k-1]
			solcost += inst.SCR[r,k] + DPHCR[r,k,oldk-1]
			retcost += inst.SCR[r,k] + DPHCR[r,k,oldk-1]
			push!(SETR[r],k)
		end
	end
	#println("DP retailers cost = $(retcost)")

	warecost = 0
	for w in 1:inst.NW
		k = GWK[w,inst.NT]
		solcost += inst.SCW[w,k] + DPHCW[w,k,inst.NT]
		warecost += inst.SCW[w,k] + DPHCW[w,k,inst.NT]
		push!(SETW[w],k)
		while k > 1
			oldk = k
			k = GWK[w,k-1]
			solcost += inst.SCW[w,k] + DPHCW[w,k,oldk-1]
			warecost += inst.SCW[w,k] + DPHCW[w,k,oldk-1]
			push!(SETW[w],k)
		end
	end

	#println("DP warehouses cost = $(warecost)")

	plantcost = 0
	k = GPK[inst.NT]
	#println("interval [$(k),$(inst.NT)]")
	intSETP,intcost = recovergenerationintervalsolutionDLSCC(inst.HCP,inst.SCP,newDP,k,inst.NT,inst.C[1])
	#println("SETP = $(intSETP)")
	#println("cost = $(intcost)")
	solcost += intcost
	plantcost += intcost
	for k in 1:length(intSETP)
		push!(SETP,intSETP[k])
	end
	while k > 1
		oldk = k
		k = GPK[k-1]

		#println("interval [$(k),$(oldk-1)]")
		intSETP,intcost = recovergenerationintervalsolutionDLSCC(inst.HCP,inst.SCP,newDP,k,oldk-1,inst.C[1])
		#println("SETP = $(intSETP)")
		#println("cost = $(intcost)")
		solcost += intcost
		plantcost += intcost
		for k in 1:length(intSETP)
			push!(SETP,intSETP[k])
		end

	end

	#println("DP plant cost dd = $(plantcost)")


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






function RandomizedDPHeuristicBottomUpCC(inst,params)
	println("RandomizedDPHeuristicBottomUpCC")
	#G= zeros(Float64,inst.NR,inst.NT)
	#GK= zeros(Int,inst.NR,inst.NT)

	bestSETR = Vector{Vector{Int}}(undef,inst.NR)
	bestSETW = Vector{Vector{Int}}(undef,inst.NW)
	bestSETP = Vector{Int}

	t1 = time_ns()

	bestcost = 999999999.0
	alpha = params.dpalpha
	dpredperiods = params.dpredperiods
	dpcostreductionret = params.dpcostreductionret
	dpcostreductionware = params.dpcostreductionware
	rng = MersenneTwister(params.dpseed);

	greedycost = 0.0

	randSCR =zeros(Int,inst.NR,inst.NT)
	randSCW = zeros(Int,inst.NW,inst.NT)

	sumcosts = 0.0

	feasiblerounds = 0

	for iter in 1:params.dpmaxiter

		if params.dptype == 1
			randValuesSCR = rand!(rng,zeros(inst.NR,inst.NT))*params.dpalpha*(-1)
			for r in 1:inst.NR, t in 1:ceil(Int,dpredperiods*inst.NT)
				randValuesSCR[r,t] -= dpcostreductionret
				#randValuesSCR[r,t] += 0.03*(ceil(inst.NT/2)-t)
			end
			#println("randValuesSCR = ",randValuesSCR)
			randValuesSCW = rand!(rng,zeros(inst.NW,inst.NT))*params.dpalpha*(-1)
			#println("randValuesSCW = ",randValuesSCW)
			for w in 1:inst.NW, t in 1:ceil(Int,dpredperiods*inst.NT)
				#randValuesSCW[w,t] -= 0.01*(ceil(inst.NT/2)-t)
				randValuesSCW[w,t] -= dpcostreductionware
			end
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



		G,GK = DPRetailer(inst.HCR,randSCR,inst.D,inst.NR,inst.NT)
		newDW = calculateDW(inst,GK)
		GW, GWK = DPWarehouse(inst.HCW,randSCW,newDW,inst.NW,inst.NT)

		newDP = calculateDP(inst,GWK,newDW)


		if feasibleCapacitatedProblem(newDP,inst.NT,inst.C[1])

			feasiblerounds += 1

			GP,GPK = DPPlantCC(inst.HCP,inst.SCP,newDP,inst.NT,inst.C[1])


			solcost = sum(G[:,inst.NT]) + sum(GW[:,inst.NT]) + GP[inst.NT]
			#println("solcost = ",solcost)

			#return

			SETR,SETW,SETP,calcsolcost = RandomizedDPCCrecoversolution(inst,GK,GWK,GPK,newDW,newDP)

			sumcosts += calcsolcost

			if calcsolcost < bestcost
				bestcost = calcsolcost
				bestSETR = SETR
				bestSETW =SETW
				bestSETP = SETP
			end

 			#println("##################### solcostcalculated = ",calcsolcost)
			#println("\n\n\n")



		end

	end

	t2 = time_ns()
	elapsedtime = (t2-t1)/1.0e9
	println("Elapsed = ", elapsedtime)
	#println("Greedycost = ",greedycost)
	println("Bestcost = ",bestcost)

	println("Feasible rounds = ",100*feasiblerounds/params.dpmaxiter)

	averagecost = sumcosts/params.dpmaxiter

	bestcostdpheur = bestcost
	elapsedtimefo = 0.0

	if params.form == "randdpheurCCfo" || params.form == "lr"

		#FixAndOptimize.lpbasedimprovement(inst,params,bestSETR,bestSETW,bestSETP,bestcost)
		t1 = time_ns()

		bestcost = FixAndOptimize.FixAndOptimizeEchelonStockFormulation2(inst,params,bestSETR,bestSETW,bestSETP,bestcost)
		#FixAndOptimize.FixAndOptimizeStandardFormulation(inst,params,bestSETR,bestSETW,bestSETP,bestcost)
		t2 = time_ns()
		elapsedtimefo = (t2-t1)/1.0e9
	end


	if params.form == "randdpheurCC"
		open("saida.txt","a") do f
			write(f,";$bestcost;$averagecost;$elapsedtime;$(params.dptype);$(params.dpalpha);$(params.dpredperiods);$(params.dpcostreductionware);$(params.dpcostreductionret);$(100*feasiblerounds/params.dpmaxiter) \n")
		end
	elseif params.form == "randdpheurCCfo"
			open("saida.txt","a") do f
				write(f,";$bestcostdpheur;$averagecost;$elapsedtime;$(params.dptype);$(params.dpalpha);$(params.dpredperiods);$(params.dpcostreductionware);$(params.dpcostreductionret);$(100*feasiblerounds/params.dpmaxiter);$bestcost;$elapsedtimefo \n")
			end
	elseif params.form == "randdpheurCCmc"

		#println("SETR = ",SETR)
		#println("SETW = ",SETW)
		#println("SETP = ",SETP)

		#optvalue = partialmulticommodityRounding(inst,params,bestSETR,bestSETW,bestSETP,bestcost)

		optvalue = standardFormulation(inst,params,bestSETR,bestSETW,bestSETP,bestcost)

	end

	#println("SETR = ",SETR)
	#println("SETW = ",SETW)
	#println("SETP = ",SETP)

	#optvalue = partialmulticommodityRounding(inst,params,bestSETR,bestSETW,bestSETP)

	#FixAndOptimize.FixAndOptimizeStandardFormulation(inst,params,bestSETR,bestSETW,bestSETP,bestcost)

	#gap = 100*(bestcost-optvalue)/bestcost
	#println("Heuristic gap = ",gap)

	return bestSETR,bestSETW,bestSETP,bestcost

end



function RandomizedDPHeuristicBottomUpCCAdaptative(inst,params)
	println("RandomizedDPHeuristicBottomUpCCAdaptative")
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

	feasiblerounds = 0




	for iter in 1:params.dpmaxiter

		if params.dptype == 1
			randValuesSCR = rand!(rng,zeros(inst.NR,inst.NT))*params.dpalpha*(-1)
			#for r in 1:inst.NR, t in 1:ceil(Int,0.7*inst.NT)
			#	randValuesSCR[r,t] *= 0.25
				#randValuesSCR[r,t] += 0.03*(ceil(inst.NT/2)-t)
			#end
			#println("randValuesSCR = ",randValuesSCR)
			randValuesSCW = rand!(rng,zeros(inst.NW,inst.NT))*params.dpalpha*(-1)
			#println("randValuesSCW = ",randValuesSCW)
			#for w in 1:inst.NW, t in 1:ceil(Int,0.7*inst.NT)
			#	#randValuesSCW[w,t] -= 0.01*(ceil(inst.NT/2)-t)
			#	randValuesSCW[w,t] *= 0.36
			#end
		end

		for r in 1:inst.NR,t in 1:inst.NT
			randSCR[r,t] = inst.SCR[r,t] + round(inst.SCR[r,t]*(randValuesSCR[r,t]))
			if t <= 0.7*inst.NT
				randSCR[r,t] = round(randSCR[r,t]*0.9)
			end
		end
		for w in 1:inst.NW,t in 1:inst.NT
			randSCW[w,t] = inst.SCW[w,t] + round(inst.SCW[w,t]*(randValuesSCW[w,t]))
			if t <= 0.7*inst.NT
				randSCW[w,t] = round(randSCW[w,t]*0.9)
			end
		end



		G,GK = DPRetailer(inst.HCR,randSCR,inst.D,inst.NR,inst.NT)
		newDW = calculateDW(inst,GK)
		GW, GWK = DPWarehouse(inst.HCW,randSCW,newDW,inst.NW,inst.NT)

		newDP = calculateDP(inst,GWK,newDW)


		if feasibleCapacitatedProblem(newDP,inst.NT,inst.C[1])

			feasiblerounds += 1

			GP,GPK = DPPlantCC(inst.HCP,inst.SCP,newDP,inst.NT,inst.C[1])


			solcost = sum(G[:,inst.NT]) + sum(GW[:,inst.NT]) + GP[inst.NT]
			#println("solcost = ",solcost)

			#return

			SETR,SETW,SETP,calcsolcost = RandomizedDPCCrecoversolution(inst,GK,GWK,GPK,newDW,newDP)

			sumcosts += calcsolcost

			if calcsolcost < bestcost
				bestcost = calcsolcost
				bestSETR = SETR
				bestSETW =SETW
				bestSETP = SETP
			end

 			#println("##################### solcostcalculated = ",calcsolcost)
			#println("\n\n\n")



		end

	end

	t2 = time_ns()
	elapsedtime = (t2-t1)/1.0e9
	println("Elapsed = ", elapsedtime)
	#println("Greedycost = ",greedycost)
	println("Bestcost = ",bestcost)

	println("Feasible rounds = ",100*feasiblerounds/params.dpmaxiter)

	averagecost = sumcosts/params.dpmaxiter

	if params.form == "randdpheurCCfo"

		FixAndOptimize.lpbasedimprovement(inst,params,bestSETR,bestSETW,bestSETP,bestcost)

		FixAndOptimize.FixAndOptimizeEchelonStockFormulation2(inst,params,bestSETR,bestSETW,bestSETP,bestcost)
		#FixAndOptimize.FixAndOptimizeStandardFormulation(inst,params,bestSETR,bestSETW,bestSETP,bestcost)
	end


	if params.form == "randdpheurCC"
		open("saida.txt","a") do f
			write(f,";$bestcost;$averagecost;$elapsedtime;$(params.dptype);$(params.dpalpha) \n")
		end
	elseif params.form == "randdpheurCCmc"

		#println("SETR = ",SETR)
		#println("SETW = ",SETW)
		#println("SETP = ",SETP)

		#optvalue = partialmulticommodityRounding(inst,params,bestSETR,bestSETW,bestSETP,bestcost)

		optvalue = standardFormulation(inst,params,bestSETR,bestSETW,bestSETP,bestcost)

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



function DPPlantCC(HCP,SCP,newDP,NT,C)
 	#solve constant capacity lot-sizing for the plant
	#println("Running DPPlantCC")

	GP= zeros(Float64,NT) #
	GPK= zeros(Int,NT)




	REGCOST = regenerationIntervalCost(HCP,SCP,newDP,NT,C)



	GP[1] = REGCOST[1,1]
	GPK[1] = 1
	for t in 2:NT
		minvalue =  REGCOST[1,t]
		minindex = 1
		for k in 2:t
			kvalue = GP[k-1] + REGCOST[k,t]
			if kvalue < minvalue
				minvalue = kvalue
				minindex = k
			end
		end
		GP[t] = minvalue
		GPK[t] = minindex
	end

	#println("GPK = ",GPK)
	totalcost = sum(GP[NT])
	#println("Totalcost = ",totalcost)

	return GP, GPK



end # DPPlant

function regenerationIntervalCost(HCP,SCP,newDP,NT,C)
	#println("Running regenerationIntervalCost [$(a),$(b)]")

	REGCOST = zeros(Float64,NT,NT)

	for a in 1:NT, b in a:NT
		REGCOST[a,b] = solveDLSCC(HCP,SCP,newDP,a,b,C)
	end

	return REGCOST

end #regenerationIntervalCost


function solveDLSCC(HCP,SCP,newDP,a,b,C)
	#println("Running solveDLSCC [$(a),$(b)]")

	origicost = 0.0

	#println("Capacity = $(C)")
	#println("newDP = $(newDP)")

	#println("SCP = ",SCP)

	newSCP = zeros(Float64,b)

	for t in a:b
		newSCP[t] = SCP[t] + C*HCP[1]*(b-t)
	end
	#println("newSCP = ",newSCP)
	#println("HCP = ",HCP)

	solvector = zeros(Int,b-a+1)

	deltat = zeros(Int,b-a+1)
	for t in a:b
		deltat[t-a+1] = t-a+1 - ceil( sum(newDP[a:t])/C )
		if deltat[t-a+1] < 0
			cost = 999999.0
			#println("cost = ",cost)
			return  cost
		end
	end
	#println("delta = ",deltat)


	sortedperiods = Vector{Tuple{Int,Float64}}()
	for t in a:b
		push!(sortedperiods,(t,newSCP[t]))
	end
	sort!(sortedperiods,by = last,rev = true)

	#println("here")

	#println(sortedperiods)

	cost = 0.0
	for k in 1:b-a+1
		jk = sortedperiods[k][1]
		#println("1")
		#println("jk = ",jk)
		mind = minimum(deltat[jk-a+1:b-a+1])
		if mind > 0
			#println("2")
			for t in a+jk-1:b
				deltat[t-a+1] = max(0,deltat[t-a+1]-1)
			end
		else
			#println("3")
			solvector[jk-a+1] = 1
			cost+= newSCP[jk]
			origicost += SCP[jk]
		end
		#println("4")
	end

	#for k in b:-1:a

		#if solvector[jk-a+1] == 1

	newProdCost = zeros(Float64,b)
	toremove = 0.0
	for t in a:b
		newProdCost[t] = HCP[1]*(b-t)
		if t<b
			toremove += HCP[1]*sum(newDP[a:t])
		end
	end

	#println("newProdCost = ",newProdCost)

 	toserve = 0.0
	for t in b:-1:a
		toserve += newDP[t]

		if solvector[t-a+1] == 1
			#println("period: $(t)")
			origicost += newProdCost[t]*min(C,toserve)
			toserve -= min(C,toserve)
		end

	end

	origicost -= toremove

	#end
	#println("const now 1: ",cost)
	for t in a:b-1
		cost += HCP[1]*sum(newDP[t+1:b])
	end
	#println("const now 2: ",cost)
	cost -= HCP[1]*max(0,b-a)*C*ceil(sum(newDP[a:b])/C)
	#println("const now 3: ",cost)

	#println(solvector)
	#println("cost = ",cost)


	#println("origicost = ",origicost)

	return origicost

end #regenerationIntervalCost


function recovergenerationintervalsolutionDLSCC(HCP,SCP,newDP,a,b,C)
	#println("Running recovergenerationintervalsolutionDLSCC [$(a),$(b)]")
	#println("Capacity = ",C)
	#println("newDP = ",newDP)

	setupperiods = Vector{Int}()

	deltat = zeros(Int,b-a+1)
	for t in a:b
		deltat[t-a+1] = t-a+1 - ceil( sum(newDP[a:t])/C )
		if deltat[t-a+1] < 0
			cost = 999999.0
	#		println("cost = ",cost)
			return  cost
		end
	end
	#println("delta = ",deltat)


	newSCP = zeros(Float64,b)

	for t in a:b
		newSCP[t] = SCP[t] + C*HCP[1]*max(0,b-t)
	end
	#println("newSCP = ",newSCP)

	sortedperiods = Vector{Tuple{Int,Float64}}()
	for t in a:b
		push!(sortedperiods,(t,newSCP[t]))
	end
	sort!(sortedperiods,by = last,rev = true)

	#println("here")

	#println(sortedperiods)

	origicost = 0.0

	cost = 0.0
	for k in 1:b-a+1
		jk = sortedperiods[k][1]
		#println("1")
		#println("jk = ",jk)
		#println("delta from $(jk-a+1) to $(b-a+1) = $(deltat)")
		mind = minimum(deltat[jk-a+1:b-a+1])
		if mind > 0
			#println("2")
			for t in jk-a+1:b-a+1
				deltat[t] = max(0,deltat[t]-1)
			end
		else
			#println("3")
			push!(setupperiods,jk)
			cost+= newSCP[jk]
			origicost += SCP[jk]
		end
		#println("4")
	end


	newProdCost = zeros(Float64,b)
	toremove = 0.0
	for t in a:b
		newProdCost[t] = HCP[1]*(b-t)
		if t<b
			toremove += HCP[1]*sum(newDP[a:t])
		end
	end

	#println("newProdCost = ",newProdCost)

 	toserve = 0.0
	for t in b:-1:a
		toserve += newDP[t]

		if t in setupperiods
	#		println("period: $(t)")
			origicost += newProdCost[t]*min(C,toserve)
			toserve -= min(C,toserve)
		end

	end

	origicost -= toremove


	for t in a:b-1
		cost += HCP[1]*sum(newDP[t+1:b])
	end
	#println("const now 2: ",cost)
	cost -= HCP[1]*max(0,b-1-a+1)*C*ceil(sum(newDP[a:b])/C)

	#println("setupperiods = ",setupperiods)
	#println("cost = ",cost)

	return setupperiods,origicost

end #regenerationIntervalCost

function feasibleCapacitatedProblem(newDP,NT,C)
	b = NT
	a = 1
	deltat = zeros(Int,b-a+1)
	for t in a:b
		deltat[t-a+1] = t-a+1 - ceil( sum(newDP[a:t])/C )
		if deltat[t-a+1] < 0
			infeas = sum(newDP[a:t]) - C*(t-a+1)
			#println("Infeasibility in period $(t-a+1): $(deltat[t-a+1])  -->  $(infeas)")
			return false
		end
	end


	return true
end


function solveproblemQ(HCP,SCP,newDP,NT,C,j)

end #solveproblemQ




function standardFormulation(inst::InstanceData, params::ParameterData,SETR,SETW,SETP,dpheurvalue)
    #println("Running Formulations.standardFormulation")

	if params.solver == "Gurobi"
		if params.disablesolver == 1 #Disable gurobi cuts and presolve
			if params.maxnodes < 999.0
				model = Model(solver = GurobiSolver(TimeLimit=params.maxtime,MIPGap=params.tolgap,CliqueCuts=0, CoverCuts=0, FlowCoverCuts=0,FlowPathCuts=0,MIRCuts=0,NetworkCuts=0,GomoryPasses=0, PreCrush=1,NodeLimit=params.maxnodes))
			else
				model = Model(solver = GurobiSolver(TimeLimit=params.maxtime,MIPGap=params.tolgap,CliqueCuts=0, CoverCuts=0, FlowCoverCuts=0,FlowPathCuts=0,MIRCuts=0,NetworkCuts=0,GomoryPasses=0, PreCrush=1))
			end
		else
			if params.maxnodes < 999.0
        		model = Model(solver = GurobiSolver(TimeLimit=params.maxtime,MIPGap=params.tolgap,PreCrush=1,NodeLimit=params.maxnodes))
			else
				model = Model(solver = GurobiSolver(TimeLimit=params.maxtime,MIPGap=params.tolgap,PreCrush=1))
			end
		end
	elseif params.solver == "Cplex"
        model = Model(solver = CplexSolver(CPX_PARAM_TILIM=params.maxtime,CPX_PARAM_EPGAP=params.tolgap))
	else
		println("No solver selected")
		return 0
	end

	### Defining variables ###
    @variable(model,0 <= xp[i=1:inst.NP,t=1:inst.NT] <= Inf)
    @variable(model,0 <= xr[i=1:inst.NR,t=1:inst.NT] <= Inf)
    @variable(model,0 <= xw[i=1:inst.NW,t=1:inst.NT] <= Inf)
    @variable(model, yp[i=1:inst.NP,t=1:inst.NT], Bin)
    @variable(model, yr[i=1:inst.NR,t=1:inst.NT], Bin)
    @variable(model, yw[i=1:inst.NW,t=1:inst.NT], Bin)
    @variable(model,0 <= sp[i=1:inst.NP,t=0:inst.NT] <= Inf)
    @variable(model,0 <= sr[i=1:inst.NR,t=0:inst.NT] <= Inf)
    @variable(model,0 <= sw[i=1:inst.NW,t=0:inst.NT] <= Inf)


	println("SETP = ",SETP)
	println("SETW = ",SETW)
	println("SETR = ",SETR)

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










	#println("Defined variables")
	### Objective function ###
	@objective(model, Min,
        sum(inst.SCP[i,t]*yp[i,t] + inst.HCP[i]*sp[i,t] for i=1:inst.NP, t=1:inst.NT)
        + sum(inst.SCR[i,t]*yr[i,t] + inst.HCR[i]*sr[i,t] for i=1:inst.NR, t=1:inst.NT)
        + sum(inst.SCW[i,t]*yw[i,t] + inst.HCW[i]*sw[i,t] for i=1:inst.NW, t=1:inst.NT)
	)

	#println("Finished objective function")

	### Setup constraints ###
	@constraint(model,
				setupP[t=1:inst.NT], xp[1,t] <= sum(inst.DP[1,k] for k in t:inst.NT)*yp[1,t]
				)
	@constraint(model,
				setupW[i=1:inst.NW, t=1:inst.NT], xw[i,t] <= sum(inst.DW[i,k] for k in t:inst.NT)*yw[i,t]
				)
	@constraint(model,
				setupR[i=1:inst.NR, t=1:inst.NT], xr[i,t] <= sum(inst.D[i,k] for k in t:inst.NT)*yr[i,t]
				)

	#println("Finished setup constraints")

	### no initial inventory
	@constraint(model,
				zeroinvP[p=1:inst.NP],
				sp[p,0] == 0
				)
	@constraint(model,
				zeroinvW[w=1:inst.NW],
				sw[w,0] == 0
				)
	@constraint(model,
				zeroinvR[r=1:inst.NR],
				sr[r,0] == 0
				)

	#println("Finished zero inventory constraints")

	### Balance constraints ###
	@constraint(model,
				balanceP[i=1:inst.NP, t=1:inst.NT],
				sp[i,t-1] + xp[i,t] == sum(xw[w,t] for w in 1:inst.NW) + sp[i,t]
				)
	@constraint(model,
				balanceW[w=1:inst.NW,t=1:inst.NT],
				sw[w,t-1] + xw[w,t] == sum(xr[inst.DeltaW[w][k],t]   for k in 1:length(inst.DeltaW[w])) + sw[w,t]
				)
	@constraint(model,
				balanceR[i=1:inst.NR,t=1:inst.NT],
				sr[i,t-1] + xr[i,t] == inst.D[i,t] + sr[i,t]
				)

	if params.capacity != 0.0
		@constraint(model,
					capP[t=1:inst.NT], xp[1,t] <= min(inst.C[t],sum(inst.DP[1,k] for k in t:inst.NT))*yp[1,t]
					)
	end

	### Capacity constraints ###
	#@constraint(model, capacity[i=1:inst.NP, t=1:inst.NT], xp[i,t] <= min(inst.DP[i,t],inst.C[t])*yp[i,t])

	#writeLP(model,"modelo.lp",genericnames=false)

	t1 = time_ns()
	status = solve(model)
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
		write(f,";$(params.form);$bestbound;$bestsol;$gap;$time;$numnodes;$(params.disablesolver) \n")
	end

#	if params.printsol == 1
#		printStandardFormulationSolution(inst,xp,xr,xw,yp,yr,yw,sp,sr,sw)
#	end

end #function standardFormulation()








function partialmulticommodityRounding(inst::InstanceData, params::ParameterData,SETR,SETW,SETP,dpheurvalue)
	#println("Running NewFormulations.partialmulticommodityFormulation")
	if params.solver == "Gurobi"
		if params.disablesolver == 1 #Disable gurobi cuts and presolve
		   if params.maxnodes < 999.0
			   model = Model(solver = GurobiSolver(TimeLimit=params.maxtime,MIPGap=params.tolgap,CliqueCuts=0, CoverCuts=0, FlowCoverCuts=0,FlowPathCuts=0,MIRCuts=0,NetworkCuts=0,GomoryPasses=0, PreCrush=1,NodeLimit=params.maxnodes))
		   else
			   model = Model(solver = GurobiSolver(TimeLimit=params.maxtime,MIPGap=params.tolgap,CliqueCuts=0, CoverCuts=0, FlowCoverCuts=0,FlowPathCuts=0,MIRCuts=0,NetworkCuts=0,GomoryPasses=0, PreCrush=1))
		   end
	   else
		   if params.maxnodes < 999.0
			   model = Model(solver = GurobiSolver(TimeLimit=params.maxtime,MIPGap=params.tolgap,PreCrush=1,NodeLimit=params.maxnodes))
		   else
			   model = Model(solver = GurobiSolver(TimeLimit=params.maxtime,MIPGap=params.tolgap,PreCrush=1))
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


	println("SETP = ",SETP)
	println("SETW = ",SETW)
	println("SETR = ",SETR)

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
