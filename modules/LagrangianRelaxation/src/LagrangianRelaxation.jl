module LagrangianRelaxation

using JuMP
using Gurobi
using CPLEX
using Data
using Parameters
using DPHeuristicsCC

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

function lagrangianRelaxation(inst::InstanceData, params::ParameterData)

	lrmaxiter = 120

	lowerbounds = Vector{Float64}()
	upperbounds = Vector{Float64}()

	alpha = 1.6
	lambda = zeros(Float64,inst.NT)
	G = zeros(Float64,inst.NT)

	bestupperbound = Inf
	bestlowerbound = -Inf

	bestSETR,bestSETW,bestSETP,bestupperbound = DPHeuristicsCC.RandomizedDPHeuristicBottomUpCC(inst,params)

	noimprovements = 0

	for iter in 1:lrmaxiter
		println("iteration $(iter)")

		value_current,xp_vals,yp_vals = standardFormulation(inst,params,lambda)
		#value_current, xp_vals, yp_vals = multicommodityFormulation(inst,params,lambda)
		println(value_current,"   X   ",bestlowerbound)

		if value_current > bestlowerbound + 0.000001
			bestlowerbound = value_current
			noimprovements = 0
		elseif iter == 1
			bestlowerbound = value_current
			noimprovements = 0
		else
			noimprovements += 1
		end

		#bestSETR,bestSETW,bestSETP,upperbound_current = DPHeuristicsCC.RandomizedDPHeuristicBottomUpCC(inst,params)

		upperbound_current = bestupperbound

		#push!(lowerbounds,value_current)
		push!(lowerbounds,bestlowerbound)
		push!(upperbounds,upperbound_current)

		if upperbound_current < bestupperbound
			bestupperbound = upperbound_current
		end

		for t in 1:inst.NT
			G[t] = xp_vals[1,t] - min(inst.C[t],sum(inst.DP[1,k] for k in t:inst.NT))*yp_vals[1,t]
		end
		println("Gradient vector = ",G)

		stepsize = (bestupperbound*1.05-bestlowerbound)*alpha/sum(G.^2)
		println("stepsize = ",stepsize)

		for t in 1:inst.NT
			lambda[t] = max(0.0,lambda[t]+stepsize*G[t])
		end
		println("lambda = ",lambda)

		if noimprovements > 3
			#println("BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB")
			alpha = alpha/2
			noimprovements = 0
			if alpha <= 0.005
				break
			end
		end

		open("lagrange.txt","a") do f
			write(f,"$(iter);$bestlowerbound;$bestupperbound;$stepsize; \n")
		end

	end

	print("upperbounds = ",upperbounds)
	println("lowerbounds = ",lowerbounds)

end


function standardFormulation(inst::InstanceData, params::ParameterData,lambda)

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
    #@variable(model, yp[i=1:inst.NP,t=1:inst.NT], Bin)
    #@variable(model, yr[i=1:inst.NR,t=1:inst.NT], Bin)
    #@variable(model, yw[i=1:inst.NW,t=1:inst.NT], Bin)
	@variable(model, 0 <= yp[p=1:inst.NP,t=1:inst.NT] <= 1)
	@variable(model, 0 <= yr[r=1:inst.NR,t=1:inst.NT] <= 1)
	@variable(model, 0 <= yw[w=1:inst.NW,t=1:inst.NT] <= 1)
	@variable(model,0 <= sp[i=1:inst.NP,t=0:inst.NT] <= Inf)
    @variable(model,0 <= sr[i=1:inst.NR,t=0:inst.NT] <= Inf)
    @variable(model,0 <= sw[i=1:inst.NW,t=0:inst.NT] <= Inf)

	### Objective function ###
	@objective(model, Min,
        sum(inst.SCP[i,t]*yp[i,t] + inst.HCP[i]*sp[i,t] for i=1:inst.NP, t=1:inst.NT)
        + sum(inst.SCR[i,t]*yr[i,t] + inst.HCR[i]*sr[i,t] for i=1:inst.NR, t=1:inst.NT)
        + sum(inst.SCW[i,t]*yw[i,t] + inst.HCW[i]*sw[i,t] for i=1:inst.NW, t=1:inst.NT)
		+ sum( lambda[t]*(xp[1,t] - min(inst.C[t],sum(inst.DP[1,k] for k in t:inst.NT))*yp[1,t]) for t=1:inst.NT)
	)

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

	#writeLP(model,"modelo.lp",genericnames=false)

	#t1 = time_ns()
	status = solve(model)
	#t2 = time_ns()
	#elapsedtime = (t2-t1)/1.0e9

	#bestsol = getobjectivevalue(model)
	bestbound = getobjbound(model)
	#numnodes = getnodecount(model)
	#time = getsolvetime(model)
	#gap = 100*(bestsol-bestbound)/bestsol
	#println("bestsol = ", bestsol)
	#println("bestbound = ", bestbound)
    #println("gap = ", gap)
    #println("time = ", time)
    #println("nodes = ", numnodes)

	#open("saida.txt","a") do f
	#	write(f,";$(params.form);$bestbound;$bestsol;$gap;$time;$numnodes;$(params.disablesolver) \n")
	#end

	#if params.printsol == 1
	#	printStandardFormulationSolution(inst,xp,xr,xw,yp,yr,yw,sp,sr,sw)
	#end

	xp_vals = getvalue(xp)
	yp_vals = getvalue(yp)

	#println("xp_vals = ",xp_vals)
	#println("xp_vals = ",yp_vals)

	return bestbound,xp_vals,yp_vals

end #function standardFormulation()


function multicommodityFormulation(inst::InstanceData, params::ParameterData,lambda)

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
        model = Model(solver = CplexSolver(CPX_PARAM_TILIM=params.maxtime,CPX_PARAM_EPGAP=params.tolgap))
	else
		println("No solver selected")
		return 0
	end

	### Defining variables ###
	@variable(model,0 <= xxp[r=1:inst.NR,t=1:inst.NT,k=t:inst.NT] <= Inf)
	@variable(model,0 <= xxr[r=1:inst.NR,t=1:inst.NT,k=t:inst.NT] <= Inf)
	@variable(model,0 <= xxw[r=1:inst.NR,t=1:inst.NT,k=t:inst.NT] <= Inf)
	#@variable(model, yp[p=1:inst.NP,t=1:inst.NT], Bin)
	#@variable(model, yr[r=1:inst.NR,t=1:inst.NT], Bin)
	#@variable(model, yw[w=1:inst.NW,t=1:inst.NT], Bin)
	@variable(model, 0 <= yp[p=1:inst.NP,t=1:inst.NT] <= 1)
	@variable(model, 0 <= yr[r=1:inst.NR,t=1:inst.NT] <= 1)
	@variable(model, 0 <= yw[w=1:inst.NW,t=1:inst.NT] <= 1)
	@variable(model,0 <= ssp[r=1:inst.NR,t=0:inst.NT-1,k=t+1:inst.NT] <= Inf)
	@variable(model,0 <= ssr[r=1:inst.NR,t=0:inst.NT-1,k=t+1:inst.NT] <= Inf)
	@variable(model,0 <= ssw[r=1:inst.NR,t=0:inst.NT-1,k=t+1:inst.NT] <= Inf)
	@variable(model,0 <= xp[i=1:inst.NP,t=1:inst.NT] <= Inf)

	### Objective function ###
	@objective(model, Min,
		sum(inst.SCP[p,t]*yp[p,t] for p=1:inst.NP, t=1:inst.NT) + sum(inst.HCP[1]*ssp[r,k,t] for r=1:inst.NR,k=1:inst.NT-1,t= k+1:inst.NT)
		+ sum(inst.SCW[w,t]*yw[w,t]  for w=1:inst.NW, t=1:inst.NT) + sum(inst.HCW[inst.DeltamR[r]]*ssw[r,k,t] for r=1:inst.NR,k=1:inst.NT-1,t= k+1:inst.NT)
		+ sum(inst.SCR[r,t]*yr[r,t]  for r=1:inst.NR, t=1:inst.NT) + sum(inst.HCR[r]*ssr[r,k,t] for r=1:inst.NR,k=1:inst.NT-1,t= k+1:inst.NT)
		+ sum( lambda[t]*(xp[1,t] - min(inst.C[t],sum(inst.DP[1,k] for k in t:inst.NT))*yp[1,t]) for t=1:inst.NT)
	)

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

	@constraint(model,
				prodVariables[t=1:inst.NT], xp[1,t] == sum(xxp[r,t,k] for r=1:inst.NR,k=t:inst.NT)
	)

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

	### Balance constraints ###
	@constraint(model,
				balanceP[r=1:inst.NR, t=1:inst.NT, k=t:inst.NT],
				ssp[r,t-1,k] + xxp[r,t,k] == xxw[r,t,k] + sum(ssp[r,t,z] for z in k:k if t<z)
				)
	@constraint(model,
				balanceW[r=1:inst.NR,t=1:inst.NT,k=t:inst.NT],
				ssw[r,t-1,k] + xxw[r,t,k] == xxr[r,t,k] + sum(ssw[r,t,z] for z in k:k if t<z)
				)
	@constraint(model,
				balanceR[r=1:inst.NR,t=1:inst.NT,k=t:inst.NT],
				ssr[r,t-1,k] + xxr[r,t,k] == sum(inst.D[r,z] for z in k:k if t==k) + sum(ssr[r,t,z] for z in k:k if t<k)
				)

	#writeLP(model,"modelo.lp",genericnames=false)

	#t1 = time_ns()
	status = solve(model)
	#t2 = time_ns()
	#elapsedtime = (t2-t1)/1.0e9

	#bestsol = getobjectivevalue(model)
	bestbound = getobjbound(model)
	#numnodes = getnodecount(model)
	#time = getsolvetime(model)
	#gap = 100*(bestsol-bestbound)/bestsol

	#println("bestsol = ", bestsol)
	#println("bestbound = ", bestbound)
    #println("gap = ", gap)
    #println("time = ", time)
    #println("nodes = ", numnodes)

#	open("saida.txt","a") do f
#		write(f,";$(params.form);$bestbound;$bestsol;$gap;$time;$numnodes;$(params.disablesolver) \n")
#	end

#	if params.printsol == 1
#		printStandardFormulationSolution(inst,xp,xr,xw,yp,yr,yw,sp,sr,sw)
#	end

	xp_vals = getvalue(xp)
	yp_vals = getvalue(yp)
	#println("xp_vals = ",xp_vals)
	#println("xp_vals = ",yp_vals)

	return bestbound,xp_vals,yp_vals

end #function multicommodityFormulation()

end
