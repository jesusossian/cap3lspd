module Formulations

using JuMP
using Gurobi
using CPLEX
using Data
using Parameters

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

export standardFormulation, stdFormVars, schelonStockFormulation

function standardFormulation(inst::InstanceData, params::ParameterData)
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


function multicommodityFormulation(inst::InstanceData, params::ParameterData)
	#println("Running Formulations.multicommodityFormulation")

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
	@variable(model, yp[p=1:inst.NP,t=1:inst.NT], Bin)
	@variable(model, yr[r=1:inst.NR,t=1:inst.NT], Bin)
	@variable(model, yw[w=1:inst.NW,t=1:inst.NT], Bin)
	@variable(model,0 <= ssp[r=1:inst.NR,t=0:inst.NT-1,k=t+1:inst.NT] <= Inf)
	@variable(model,0 <= ssr[r=1:inst.NR,t=0:inst.NT-1,k=t+1:inst.NT] <= Inf)
	@variable(model,0 <= ssw[r=1:inst.NR,t=0:inst.NT-1,k=t+1:inst.NT] <= Inf)

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

    if params.capacity != 0.0
        @constraint(model,
            capP[t=1:inst.NT], sum(xxp[r,t,k] for r=1:inst.NR,k=t:inst.NT) <= min(inst.C[t],sum(inst.DP[1,k] for k in t:inst.NT))*yp[1,t]
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

end #function multicommodityFormulation()


function echelonStockFormulation(inst::InstanceData, params::ParameterData)
    #println("Running Formulations.schelonStockFormulation")

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
    @variable(model,0 <= isp[i=1:inst.NP,t=0:inst.NT] <= Inf)
    @variable(model,0 <= isr[i=1:inst.NR,t=0:inst.NT] <= Inf)
    @variable(model,0 <= isw[i=1:inst.NW,t=0:inst.NT] <= Inf)

	#println("Defined variables")

	### Objective function ###
	@objective(model, Min,
        sum(inst.SCP[i,t]*yp[i,t] for i=1:inst.NP, t=1:inst.NT)
        + sum(inst.SCR[i,t]*yr[i,t] for i=1:inst.NR, t=1:inst.NT)
        + sum(inst.SCW[i,t]*yw[i,t] for i=1:inst.NW, t=1:inst.NT)
        + sum(inst.HCP[i]*isp[i,t] for i=1:inst.NP, t=1:inst.NT)
        + sum((inst.HCW[i]- inst.HCP[1])*isw[i,t] for i=1:inst.NW, t=1:inst.NT)
        + sum((inst.HCR[i]- inst.HCW[inst.DeltamR[i]])*isr[i,t] for i=1:inst.NR, t=1:inst.NT)
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

    if params.capacity != 0.0
        @constraint(model,
                    capP[t=1:inst.NT], xp[1,t] <= min(inst.C[t],sum(inst.DP[1,k] for k in t:inst.NT))*yp[1,t]
                    )
    end

	#println("Finished setup constraints")

	### no initial inventory
	@constraint(model, zeroinvP[p=1:inst.NP], isp[p,0] == 0)
	@constraint(model, zeroinvW[w=1:inst.NW], isw[w,0] == 0)
	@constraint(model, zeroinvR[r=1:inst.NR], isr[r,0] == 0)

	#println("Finished zero inventory constraints")

	### Balance constraints ###
	@constraint(model,
				balanceP[i=1:inst.NP, t=1:inst.NT],
				isp[i,t-1] + xp[i,t] == inst.DP[i,t] + isp[i,t]
				)
	@constraint(model,
				balanceW[w=1:inst.NW,t=1:inst.NT],
				isw[w,t-1] + xw[w,t] == inst.DW[w,t] + isw[w,t]
				)
	@constraint(model,
				balanceR[i=1:inst.NR,t=1:inst.NT],
				isr[i,t-1] + xr[i,t] == inst.D[i,t] + isr[i,t]
				)

    # echelon stock
	@constraint(model,
				echelonP[i=1:inst.NP, t=1:inst.NT],
				isp[i,t] >= sum(isw[w,t] for w in 1:inst.NW)
				)

	@constraint(model,
				echeloW[w=1:inst.NW,t=1:inst.NT],
				isw[w,t] >= sum(isr[inst.DeltaW[w][k],t]   for k in 1:length(inst.DeltaW[w]))
				)

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

end #function echelonStockFormulation()


function echelonStockFormulationESN(inst::InstanceData, params::ParameterData)
    #println("Running Formulations.schelonStockFormulationESN")

	if params.solver == "Gurobi"
        model = Model(solver = GurobiSolver(TimeLimit=params.maxtime,MIPGap=params.tolgap))
	elseif params.solver == "Cplex"
        model = Model(solver = CplexSolver(CPX_PARAM_TILIM=params.maxtime,CPX_PARAM_EPGAP=params.tolgap))
	else
		println("No solver selected")
		return 0
	end

	### Defining variables ###
    @variable(model, yp[i=1:inst.NP,t=1:inst.NT], Bin)
    @variable(model, yr[i=1:inst.NR,t=1:inst.NT], Bin)
    @variable(model, yw[i=1:inst.NW,t=1:inst.NT], Bin)

    @variable(model,0 <= isp[i=1:inst.NP,t=0:inst.NT] <= Inf)
    @variable(model,0 <= isr[i=1:inst.NR,t=0:inst.NT] <= Inf)
    @variable(model,0 <= isw[i=1:inst.NW,t=0:inst.NT] <= Inf)

    @variable(model,0 <= zp[i=1:inst.NP,t=1:inst.NT,k=1:inst.NT] <= Inf)
    @variable(model,0 <= zr[i=1:inst.NR,t=1:inst.NT,k=1:inst.NT] <= Inf)
    @variable(model,0 <= zw[i=1:inst.NW,t=1:inst.NT,k=1:inst.NT] <= Inf)

    sdp = zeros(Int,inst.NP,inst.NT,inst.NT)
    sdw = zeros(Int,inst.NW,inst.NT,inst.NT)
    sdr = zeros(Int,inst.NR,inst.NT,inst.NT)

    for i in 1:inst.NP
        for t in 1:inst.NT
            sdp[i,t,t] = inst.DP[i,t]
            for k in t+1:inst.NT
                sdp[i,t,k] = sdp[i,t,k-1] + inst.DP[i,k]
            end
        end
    end

    for i in 1:inst.NW
        for t in 1:inst.NT
            sdw[i,t,t] = inst.DW[i,t]
            for k in t+1:inst.NT
                sdw[i,t,k] = sdw[i,t,k-1] + inst.DW[i,k]
            end
        end
    end

    for i in 1:inst.NR
        for t in 1:inst.NT
            sdr[i,t,t] = inst.D[i,t]
            for k in t+1:inst.NT
                sdr[i,t,k] = sdr[i,t,k-1] + inst.D[i,k]
            end
        end
    end

	#println("Defined variables")

	### Objective function ###
	@objective(model, Min,
        sum(inst.SCP[i,t]*yp[i,t] for i=1:inst.NP, t=1:inst.NT)
        + sum(inst.SCR[i,t]*yr[i,t] for i=1:inst.NR, t=1:inst.NT)
        + sum(inst.SCW[i,t]*yw[i,t] for i=1:inst.NW, t=1:inst.NT)
        + sum(inst.HCP[i]*isp[i,t] for i=1:inst.NP, t=1:inst.NT)
        + sum((inst.HCW[i]- inst.HCP[1])*isw[i,t] for i=1:inst.NW, t=1:inst.NT)
        + sum((inst.HCR[i]- inst.HCW[inst.DeltamR[i]])*isr[i,t] for i=1:inst.NR, t=1:inst.NT)
	)

	#println("Finished objective function")

	### Setup constraints ###
    @constraint(model,
                setupP[i=1:inst.NP, t=1:inst.NT], sum(zp[i,t,k] for k in t:inst.NT if sdp[i,t,k] > 0) <= yp[i,t]
                )
    @constraint(model,
                setupW[i=1:inst.NW, t=1:inst.NT], sum(zw[i,t,k] for k in t:inst.NT if sdw[i,t,k] > 0) <= yw[i,t]
                )
    @constraint(model,
                setupR[i=1:inst.NR, t=1:inst.NT], sum(zr[i,t,k] for k in t:inst.NT if sdr[i,t,k] > 0) <= yr[i,t]
                )

    #### capacity ####
	if params.capacity != 0.0
		@constraint(model,
                    capP[t=1:inst.NT],
                    sum(zp[1,t,k]*sdp[1,t,k] for k in t:inst.NT) <= min(inst.C[t],sum(inst.DP[1,k] for k in t:inst.NT))*yp[1,t]
                    )
	end

	#println("Finished setup constraints")

	### Balance constraints ###
    @constraint(model,
                balanceP[i=1:inst.NP],
                sum(zp[i,1,k] for k in 1:inst.NT) == 1
                )
    @constraint(model,
                balanceW[i=1:inst.NW],
                sum(zw[i,1,k] for k in 1:inst.NT) == 1
                )
    @constraint(model,
                balanceR[i=1:inst.NR],
                sum(zr[i,1,k] for k in 1:inst.NT) == 1
                )

#########

    @constraint(model,
                balanceP1[i=1:inst.NP,t=2:inst.NT],
                sum(zp[i,l,t-1] for l in 1:(inst.NT-1)) ==  sum(zp[i,l,t] for l in t:inst.NT)
                )
    @constraint(model,
                balanceW1[i=1:inst.NW,t=2:inst.NT],
                sum(zw[i,l,t-1] for l in 1:(inst.NT-1)) == sum(zw[i,l,t] for l in t:inst.NT)
                )
    @constraint(model,
                balanceR1[i=1:inst.NP,t=2:inst.NT],
                sum(zr[i,l,t-1]  for l in 1:(inst.NT-1)) ==  sum(zr[i,l,t] for l in t:inst.NT)
                )

	@constraint(model,
				balanceP2[i=1:inst.NP, t=1:inst.NT],
				isp[i,t] == sum(sdp[i,l,k] for l=1:t,k=l:inst.NT)*zp[i,l,k] - sdp[i,1,t]
				)
	@constraint(model,
				balanceW2[i=1:inst.NW, t=1:inst.NT],
				isw[i,t] == sum(sdw[i,l,k] for l=1:t,k=l:inst.NT)*zw[i,l,k] - sdw[i,1,t]
				)
	@constraint(model,
				balanceR2[i=1:inst.NR, t=1:inst.NT],
				isr[i,t] == sum(sdr[i,l,k] for l=1:t,k=l:inst.NT)*zr[i,l,k] - sdr[i,1,t]
				)

    # echelon stock
	@constraint(model,
				echelonP[i=1:inst.NP, t=1:inst.NT],
				isp[i,t] >= sum(isw[w,t] for w in 1:inst.NW)
				)

	@constraint(model,
				echeloW[w=1:inst.NW,t=1:inst.NT],
				isw[w,t] >= sum(isr[inst.DeltaW[w][k],t]   for k in 1:length(inst.DeltaW[w]))
				)

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

end #function echelonStockFormulationESN()


function echelonStockFormulationESTP(inst::InstanceData, params::ParameterData)
    #println("Running Formulations.schelonStockFormulationESTP")

	if params.solver == "Gurobi"
        model = Model(solver = GurobiSolver(TimeLimit=params.maxtime,MIPGap=params.tolgap))
	elseif params.solver == "Cplex"
        model = Model(solver = CplexSolver(CPX_PARAM_TILIM=params.maxtime,CPX_PARAM_EPGAP=params.tolgap))
	else
		println("No solver selected")
		return 0
	end

	### Defining variables ###
    @variable(model, yp[i=1:inst.NP,t=1:inst.NT], Bin)
    @variable(model, yr[i=1:inst.NR,t=1:inst.NT], Bin)
    @variable(model, yw[i=1:inst.NW,t=1:inst.NT], Bin)

    @variable(model,0 <= isp[i=1:inst.NP,t=0:inst.NT] <= Inf)
    @variable(model,0 <= isr[i=1:inst.NR,t=0:inst.NT] <= Inf)
    @variable(model,0 <= isw[i=1:inst.NW,t=0:inst.NT] <= Inf)

    @variable(model,0 <= sxp[i=1:inst.NP,t=1:inst.NT,k=1:inst.NT] <= Inf)
    @variable(model,0 <= sxr[i=1:inst.NR,t=1:inst.NT,k=1:inst.NT] <= Inf)
    @variable(model,0 <= sxw[i=1:inst.NW,t=1:inst.NT,k=1:inst.NT] <= Inf)

    sdp = zeros(Int,inst.NP,inst.NT,inst.NT)
    sdw = zeros(Int,inst.NW,inst.NT,inst.NT)
    sdr = zeros(Int,inst.NR,inst.NT,inst.NT)

    for i in 1:inst.NP
        for t in 1:inst.NT
            sdp[i,t,t] = inst.DP[i,t]
            for k in t+1:inst.NT
                sdp[i,t,k] = sdp[i,t,k-1] + inst.DP[i,k]
            end
        end
    end

    for i in 1:inst.NW
        for t in 1:inst.NT
            sdw[i,t,t] = inst.DW[i,t]
            for k in t+1:inst.NT
                sdw[i,t,k] = sdw[i,t,k-1] + inst.DW[i,k]
            end
        end
    end

    for i in 1:inst.NR
        for t in 1:inst.NT
            sdr[i,t,t] = inst.D[i,t]
            for k in t+1:inst.NT
                sdr[i,t,k] = sdr[i,t,k-1] + inst.D[i,k]
            end
        end
    end

	#println("Defined variables")

	### Objective function ###
	@objective(model, Min,
        sum(inst.SCP[i,t]*yp[i,t] for i=1:inst.NP, t=1:inst.NT)
        + sum(inst.SCR[i,t]*yr[i,t] for i=1:inst.NR, t=1:inst.NT)
        + sum(inst.SCW[i,t]*yw[i,t] for i=1:inst.NW, t=1:inst.NT)
        + sum(inst.HCP[i]*isp[i,t] for i=1:inst.NP, t=1:inst.NT)
        + sum((inst.HCW[i]- inst.HCP[1])*isw[i,t] for i=1:inst.NW, t=1:inst.NT)
        + sum((inst.HCR[i]- inst.HCW[inst.DeltamR[i]])*isr[i,t] for i=1:inst.NR, t=1:inst.NT)
	)

 	#println("Finished objective function")

	### Setup constraints ###
    @constraint(model,
                setupP[i=1:inst.NP, k=1:inst.NT, t=1:k], sxp[i,t,k] <= inst.DP[i,k]*yp[i,t]
                )
    @constraint(model,
                setupW[i=1:inst.NW, k=1:inst.NT, t=1:k], sxw[i,t,k] <= inst.DW[i,k]*yw[i,t]
                )
    @constraint(model,
                setupR[i=1:inst.NR, k=1:inst.NT, t=1:k], sxr[i,t,k] <= inst.D[i,k]*yr[i,t]
                )

    #### capacity ####
    if params.capacity != 0.0
        @constraint(model,
                    capP[t=1:inst.NT], sum(sxp[1,t,k] for k in t:inst.NT) <= min(inst.C[t],sum(inst.DP[1,k] for k in t:inst.NT))*yp[1,t]
                    )
    end

#	println("Finished setup constraints")

	### no initial inventory
	@constraint(model, zeroinvP[p=1:inst.NP], isp[p,0] == 0)
	@constraint(model, zeroinvW[w=1:inst.NW], isw[w,0] == 0)
	@constraint(model, zeroinvR[r=1:inst.NR], isr[r,0] == 0)

	### Balance constraints ###
	@constraint(model,
				balanceP[i=1:inst.NP, t=1:inst.NT],
				isp[i,t-1] + sum(sxp[i,t,k] for k in t:inst.NT) == inst.DP[i,t] + isp[i,t]
				)
	@constraint(model,
				balanceW[w=1:inst.NW,t=1:inst.NT],
				isw[w,t-1] + sum(sxw[w,t,k] for k in t:inst.NT) == inst.DW[w,t] + isw[w,t]
				)
	@constraint(model,
				balanceR[i=1:inst.NR,t=1:inst.NT],
				isr[i,t-1] + sum(sxr[i,t,k] for k in t:inst.NT) == inst.D[i,t] + isr[i,t]
				)

	### Balance constraints ###
    @constraint(model,
                balanceP1[i=1:inst.NP,t=1:inst.NT],
                sum(sxp[i,t,k] for k in 1:t) == inst.DP[i,t]
                )
    @constraint(model,
                balanceW1[i=1:inst.NW,t=1:inst.NT],
                sum(sxw[i,t,k] for k in 1:t) == inst.DW[i,t]
                )
    @constraint(model,
                balanceR1[i=1:inst.NR,t=1:inst.NT],
                sum(sxr[i,t,k] for k in 1:t) == inst.D[i,t]
                )

    # echelon stock
	@constraint(model,
				echelonP[i=1:inst.NP, t=1:inst.NT],
				isp[i,t] >= sum(isw[w,t] for w in 1:inst.NW)
				)

	@constraint(model,
				echeloW[w=1:inst.NW,t=1:inst.NT],
				isw[w,t] >= sum(isr[inst.DeltaW[w][k],t]   for k in 1:length(inst.DeltaW[w]))
				)

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

end #function echelonStockFormulationESTP()

function echelonStockFormulationESLS(inst::InstanceData, params::ParameterData)
    #println("Running Formulations.schelonStockFormulationESLS")

	if params.solver == "Gurobi"
        model = Model(solver = GurobiSolver(TimeLimit=params.maxtime,MIPGap=params.tolgap))
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
    @variable(model,0 <= isp[i=1:inst.NP,t=0:inst.NT] <= Inf)
    @variable(model,0 <= isr[i=1:inst.NR,t=0:inst.NT] <= Inf)
    @variable(model,0 <= isw[i=1:inst.NW,t=0:inst.NT] <= Inf)

	#println("Defined variables")

	### Objective function ###
	@objective(model, Min,
        sum(inst.SCP[i,t]*yp[i,t] for i=1:inst.NP, t=1:inst.NT)
        + sum(inst.SCR[i,t]*yr[i,t] for i=1:inst.NR, t=1:inst.NT)
        + sum(inst.SCW[i,t]*yw[i,t] for i=1:inst.NW, t=1:inst.NT)
        + sum(inst.HCP[i]*isp[i,t] for i=1:inst.NP, t=1:inst.NT)
        + sum((inst.HCW[i]- inst.HCP[1])*isw[i,t] for i=1:inst.NW, t=1:inst.NT)
        + sum((inst.HCR[i]- inst.HCW[inst.DeltamR[i]])*isr[i,t] for i=1:inst.NR, t=1:inst.NT)
	)

 #   println("Finished objective function")

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

    if params.capacity != 0.0
        @constraint(model,
                    capP[t=1:inst.NT], xp[1,t] <= min(inst.C[t],sum(inst.DP[1,k] for k in t:inst.NT))*yp[1,t]
                    )
    end

	#println("Finished setup constraints")

	### no initial inventory
	@constraint(model, zeroinvP[p=1:inst.NP], isp[p,0] == 0)
	@constraint(model, zeroinvW[w=1:inst.NW], isw[w,0] == 0)
	@constraint(model, zeroinvR[r=1:inst.NR], isr[r,0] == 0)

	#println("Finished zero inventory constraints")

	### Balance constraints ###
	@constraint(model,
				balanceP[i=1:inst.NP, t=1:inst.NT],
				isp[i,t-1] + xp[i,t] == inst.DP[i,t] + isp[i,t]
				)
	@constraint(model,
				balanceW[w=1:inst.NW,t=1:inst.NT],
				isw[w,t-1] + xw[w,t] == inst.DW[w,t] + isw[w,t]
				)
	@constraint(model,
				balanceR[i=1:inst.NR,t=1:inst.NT],
				isr[i,t-1] + xr[i,t] == inst.D[i,t] + isr[i,t]
				)

    #echelon stock
	@constraint(model,
				echelonP[i=1:inst.NP, t=1:inst.NT],
				isp[i,t] >= sum(isw[w,t] for w in 1:inst.NW)
				)

	@constraint(model,
				echeloW[w=1:inst.NW,t=1:inst.NT],
				isw[w,t] >= sum(isr[inst.DeltaW[w][k],t]   for k in 1:length(inst.DeltaW[w]))
				)

	#echelon stock ls inequalities
	@constraint(model,
	            echelon1P[i=1:inst.NP, l=1:inst.NT, k=1:l],
	            isp[i,k-1] >= sum(inst.DP[i,j] for j in k:l) - sum(inst.DP[i,j]*yp[i,u] for j=k:l, u=k:j)
	            )

	@constraint(model,
	            echelon1W[i=1:inst.NW, l=1:inst.NT, k=1:l],
	            isw[i,k-1] >= sum(inst.DW[i,j] for j in k:l) - sum(inst.DW[i,j]*yw[i,u] for j=k:l, u=k:j)
	            )

	@constraint(model,
	            echelon1R[i=1:inst.NR, l=1:inst.NT, k=1:l],
	            isr[i,k-1] >= sum(inst.D[i,j] for j in k:l) - sum(inst.D[i,j]*yr[i,u] for j=k:l, u=k:j)
	            )

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

end #function echelonStockFormulationESLS()


function multicommodityFormulationMCE(inst::InstanceData, params::ParameterData)
	#println("Running Formulations.multicommodityFormulation")

	if params.solver == "Gurobi"
		 model = Model(solver = GurobiSolver(TimeLimit=params.maxtime,MIPGap=params.tolgap))
	elseif params.solver == "Cplex"
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

    if params.capacity != 0.0
        @constraint(model,
            capP[t=1:inst.NT], sum(xxp[r,t,k] for r=1:inst.NR,k=t:inst.NT) <= min(inst.C[t],sum(inst.DP[1,k] for k in t:inst.NT))*yp[1,t]
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

end #function multicommodityFormulationMCE()

end
