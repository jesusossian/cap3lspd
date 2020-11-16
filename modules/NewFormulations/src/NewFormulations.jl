module NewFormulations

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

export threelevelFormulation, stdFormVars

function threelevelFormulation(inst::InstanceData, params::ParameterData)
	#println("Running NewFormulations.threelevelFormulation")

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
	@variable(model,0 <= x0[i=1:inst.NR,t=1:inst.NT] <= Inf)
	@variable(model,0 <= x2[i=1:inst.NR,t=1:inst.NT] <= Inf)
	@variable(model,0 <= x1[i=1:inst.NR,t=1:inst.NT] <= Inf)
	@variable(model, yp[i=1:inst.NP,t=1:inst.NT], Bin)
	@variable(model, yr[i=1:inst.NR,t=1:inst.NT], Bin)
	@variable(model, yw[i=1:inst.NW,t=1:inst.NT], Bin)
	@variable(model,0 <= s0[i=1:inst.NR,t=0:inst.NT] <= Inf)
	@variable(model,0 <= s2[i=1:inst.NR,t=0:inst.NT] <= Inf)
	@variable(model,0 <= s1[i=1:inst.NR,t=0:inst.NT] <= Inf)

	#println("Defined variables")
	### Objective function ###
	@objective(model, Min,
		sum(inst.SCP[i,t]*yp[i,t]  for i=1:inst.NP, t=1:inst.NT)
		+ sum(inst.SCR[i,t]*yr[i,t] for i=1:inst.NR, t=1:inst.NT)
		+ sum(inst.SCW[i,t]*yw[i,t] for i=1:inst.NW, t=1:inst.NT)
		+ sum( inst.HCP[1]*s0[r,t] for r=1:inst.NR, t=1:inst.NT)
		+ sum( inst.HCR[r]*s2[r,t] for r=1:inst.NR, t=1:inst.NT)
		+ sum( inst.HCW[inst.DeltamR[r]]*s1[r,t] for r=1:inst.NR, t=1:inst.NT)
	)

	#println("Finished objective function")

	### Setup constraints ###
	@constraint(model,
				setupP[r=1:inst.NR,t=1:inst.NT], x0[r,t] <= sum(inst.D[r,k] for k in t:inst.NT)*yp[1,t]
				)
	@constraint(model,
				setupW[r=1:inst.NR, t=1:inst.NT], x1[r,t] <= sum(inst.D[r,k] for k in t:inst.NT)*yw[inst.DeltamR[r],t]
				)
	@constraint(model,
				setupR[r=1:inst.NR, t=1:inst.NT], x2[r,t] <= sum(inst.D[r,k] for k in t:inst.NT)*yr[r,t]
				)

	if params.capacity != 0
		@constraint(model,
					capP[t=1:inst.NT], sum(x0[r,t] for r in 1:inst.NR) <= inst.C[t]*yp[1,t]
					)
	end

	#println("Finished setup constraints")

	### no initial inventory
	@constraint(model,
				zeroinvP[r=1:inst.NR],
				s0[r,0] == 0
				)
	@constraint(model,
				zeroinvW[r=1:inst.NR],
				s1[r,0] == 0
				)
	@constraint(model,
				zeroinvR[r=1:inst.NR],
				s2[r,0] == 0
				)

	#println("Finished zero inventory constraints")

	### Balance constraints ###
	@constraint(model,
				balanceP[r=1:inst.NR, t=1:inst.NT],
				s0[r,t-1] + x0[r,t] == x1[r,t] + s0[r,t]
				)
	@constraint(model,
				balanceW[r=1:inst.NR,t=1:inst.NT],
				s1[r,t-1] + x1[r,t] == x2[r,t] + s1[r,t]
				)
	@constraint(model,
				balanceR[r=1:inst.NR,t=1:inst.NT],
				s2[r,t-1] + x2[r,t] == inst.D[r,t] + s2[r,t]
				)

######################################3
# WW inequalities
######################################

	#@constraint(model,
	#			wwineq0[r=1:inst.NR,t=1:inst.NT,l=t:inst.NT],
	#			s0[r,t-1] + s1[r,t-1] + s2[r,t-1] + sum(yp[1,k]*inst.cumdem[r,k,l] for k in t:l) >= inst.cumdem[r,t,l]
	#			)

	#@constraint(model,
	#			wwineq1[r=1:inst.NR,t=1:inst.NT,l=t:inst.NT],
	#			 s1[r,t-1] + s2[r,t-1] + sum(yw[inst.DeltamR[r],k]*inst.cumdem[r,k,l] for k in t:l) >= inst.cumdem[r,t,l]
	#			)

	#@constraint(model,
	#			wwineq2[r=1:inst.NR,t=1:inst.NT,l=t:inst.NT],
	#			 s2[r,t-1] + sum(yr[r,k]*inst.cumdem[r,k,l] for k in t:l) >= inst.cumdem[r,t,l]
	#			)

	#@constraint(model,
	#			wwineq01[r=1:inst.NR,t=1:inst.NT-1,j=t:inst.NT-1,l=j+1:inst.NT],
	#			s0[r,t-1] + s1[r,t-1] + s2[r,t-1] + sum(yp[1,k]*inst.cumdem[r,k,l] for k in t:j) + sum(yw[inst.DeltamR[r],k]*inst.cumdem[r,k,l] for k in j+1:l) >= inst.cumdem[r,t,l]
	#			)

	#@constraint(model,
	#			wwineq12[r=1:inst.NR,t=1:inst.NT-1,j=t:inst.NT-1,l=j+1:inst.NT],
	#			 s1[r,t-1] + s2[r,t-1] + sum(yw[inst.DeltamR[r],k]*inst.cumdem[r,k,l] for k in t:j) + sum(yr[r,k]*inst.cumdem[r,k,l] for k in j+1:l) >= inst.cumdem[r,t,l]
	#			)

	### Capacity constraints ###
	#@constraint(model, capacity[i=1:inst.NP, t=1:inst.NT], x0[i,t] <= min(inst.DP[i,t],inst.C[t])*yp[i,t])

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
		write(f,";$(params.form);$bestbound;$bestsol;$gap;$time;$numnodes \n")
	end

#	if params.printsol == 1
#		printthreelevelFormulationSolution(inst,xp,xr,xw,yp,yr,yw,sp,sr,sw)
#	end

end #function threelevelFormulation()


function partialmulticommodityFormulation(inst::InstanceData, params::ParameterData)
	#println("Running NewFormulations.partialmulticommodityFormulation")

	if params.capacity == 0
		TKP=ceil(inst.NT/1)
		TKW=ceil(inst.NT/1)
		TKR=ceil(inst.NT/1)
	else
		TKP=ceil(inst.NT/2)
		TKW=ceil(inst.NT/2)
		TKR=ceil(inst.NT/2)
	end

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
	@variable(model,0 <= xxp[r=1:inst.NR,t=1:inst.NT,k=t:inst.NT] <= Inf)
	@variable(model,0 <= xxr[r=1:inst.NR,t=1:inst.NT,k=t:inst.NT] <= Inf)
	@variable(model,0 <= xxw[r=1:inst.NR,t=1:inst.NT,k=t:inst.NT] <= Inf)
	@variable(model, yp[p=1:inst.NP,t=1:inst.NT], Bin)
	@variable(model, yr[r=1:inst.NR,t=1:inst.NT], Bin)
	@variable(model, yw[w=1:inst.NW,t=1:inst.NT], Bin)
	@variable(model,0 <= ssp[r=1:inst.NR,t=0:inst.NT-1,k=t+1:inst.NT] <= Inf)
	@variable(model,0 <= ssr[r=1:inst.NR,t=0:inst.NT-1,k=t+1:inst.NT] <= Inf)
	@variable(model,0 <= ssw[r=1:inst.NR,t=0:inst.NT-1,k=t+1:inst.NT] <= Inf)

	if params.capacity == 0
		useX = ones(Int,inst.NR, inst.NT, inst.NT)
		for r in 1:inst.NR, t in 1:inst.NT-1, k in t+1:inst.NT
		#for r in 6:6, t in 10:10, k in t+1:t+2 #inst.NT
			if (useX[r,t,k-1] == 0) || ( inst.D[r,k]*(k-t)*inst.HCR[r] >= inst.D[r,k]*(k-t)*inst.HCW[inst.DeltamR[r]] + inst.SCR[r,k])
			#if  ( inst.D[r,k]*(k-t)*inst.HCR[r] > inst.D[r,k]*(k-t)*inst.HCW[inst.DeltamR[r]] + inst.SCR[r])
#				println("r = $(r), [$(t),$(k)]")
#				println(inst.D[r,k]*(k-t)*inst.HCR[r],"  >=  ",inst.D[r,k]*(k-t)*inst.HCW[inst.DeltamR[r]] + inst.SCR[r,k])
#				println("dem $(inst.D[r,k]), hcr $(inst.HCR[r]), hcw $(inst.HCW[inst.DeltamR[r]]), scr $(inst.SCR[r,k])\n")
				useX[r,t,k] = 0
			end
		end
		#println("useX = ",useX)
		totalVariables = 0
		totalFixed = 0
		for r in 1:inst.NR, t in 1:inst.NT-1, k in t+1:inst.NT
			totalVariables += 1
			if useX[r,t,k] == 0
				totalFixed +=1
				@constraint(model,
							xxr[r,t,k] == 0
							)
			end
		end
		#println("Percentage of fixed = ",100*(totalFixed/totalVariables))
	end
for r=1:inst.NR,t=1:inst.NT,k=t+TKP:inst.NT
	JuMP.fix(xxp[r,t,k], 0)
end
for r=1:inst.NR,t=1:inst.NT,k=t+TKW:inst.NT
	JuMP.fix(xxw[r,t,k], 0)
end
for r=1:inst.NR,t=1:inst.NT,k=t+TKR:inst.NT
	JuMP.fix(xxr[r,t,k], 0)
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

	if params.capacity != 0.0
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
		write(f,";$(params.form);$bestbound;$bestsol;$gap;$time;$numnodes \n")
	end

#	for r=1:inst.NR,t=1:inst.NT,k=t:inst.NT
#		if(getvalue(xxw[r,t,k])>0.001)
#			println("xxw[$(r),$(t),$(k)] = ",trunc(Int,round(getvalue(xxw[r,t,k]))))
#		end
#	end
#	for r=1:inst.NR,t=1:inst.NT,k=t:inst.NT
#		if(getvalue(xxr[r,t,k])>0.00000001)
#			println("xxr[$(r),$(t),$(k)] = ",trunc(Int,round(getvalue(xxr[r,t,k]))))
#		end
#	end


#	for w=1:inst.NW, t=1:inst.NT
#		if(getvalue(yw[w,t])>0.000000001)
#			println("yw[$(w),$(t)] = ",getvalue(yw[w,t]))
#		end
#	end
#	for r=1:inst.NR, t=1:inst.NT
#		if(getvalue(yr[r,t])>0.000000001)
#			println("yr[$(r),$(t)] = ",getvalue(yr[r,t]))
#		end
#	end

#	if params.printsol == 1
#		printStandardFormulationSolution(inst,xp,xr,xw,yp,yr,yw,sp,sr,sw)
#	end


end #function multicommodityFormulation()


function enhancedthreelevelFormulation(inst::InstanceData, params::ParameterData)
	#println("Running NewFormulations.threelevelFormulation")

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
	@variable(model,0 <= x0[i=1:inst.NR,t=1:inst.NT] <= Inf)
	@variable(model,0 <= x1[i=1:inst.NR,t=1:inst.NT] <= Inf)
	@variable(model,0 <= xxr[r=1:inst.NR,t=1:inst.NT,k=t:inst.NT] <= Inf)
	@variable(model, yp[i=1:inst.NP,t=1:inst.NT], Bin)
	@variable(model, yr[i=1:inst.NR,t=1:inst.NT], Bin)
	@variable(model, yw[i=1:inst.NW,t=1:inst.NT], Bin)
	@variable(model,0 <= s0[i=1:inst.NR,t=0:inst.NT] <= Inf)
	@variable(model,0 <= ssr[r=1:inst.NR,t=0:inst.NT-1,k=t+1:inst.NT] <= Inf)
	@variable(model,0 <= s1[i=1:inst.NR,t=0:inst.NT] <= Inf)

	cumdem = zeros(Int,inst.NR,inst.NT,inst.NT)

	for r in 1:inst.NT, t in 1:inst.NT
		cumdem[r,t,t] = inst.D[r,t]
		for k in t+1:inst.NT
			cumdem[r,t,k] = cumdem[r,t,k-1] + inst.D[r,k]
		end
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
			JuMP.fix(xxr[r,t,k],0)
		end
	end
	#println("Percentage of fixed = ",100*(totalFixed/totalVariables))

	#println("Defined variables")
	### Objective function ###
	@objective(model, Min,
		sum(inst.SCP[i,t]*yp[i,t]  for i=1:inst.NP, t=1:inst.NT)
		+ sum(inst.SCR[i,t]*yr[i,t] for i=1:inst.NR, t=1:inst.NT)
		+ sum(inst.SCW[i,t]*yw[i,t] for i=1:inst.NW, t=1:inst.NT)
		+ sum( inst.HCP[1]*s0[r,t] for r=1:inst.NR, t=1:inst.NT)
		+ sum(inst.HCR[r]*ssr[r,k,t] for r=1:inst.NR,k=1:inst.NT-1,t= k+1:inst.NT)
		+ sum( inst.HCW[inst.DeltamR[r]]*s1[r,t] for r=1:inst.NR, t=1:inst.NT)
	)

	#println("Finished objective function")

	### Setup constraints ###
	@constraint(model,
				setupP[r=1:inst.NR,t=1:inst.NT], x0[r,t] <= sum(inst.D[r,k] for k in t:inst.NT)*yp[1,t]
				)
	@constraint(model,
				setupW[r=1:inst.NR, t=1:inst.NT], x1[r,t] <= sum(inst.D[r,k] for k in t:inst.NT)*yw[inst.DeltamR[r],t]
				)
		@constraint(model,
					setupR[r=1:inst.NR,t=1:inst.NT,k=t:inst.NT], xxr[r,t,k] <= inst.D[r,k]*yr[r,t]
					)

	if params.capacity != 0.0
		@constraint(model,
					capP[t=1:inst.NT], sum(x0[r,t] for r in 1:inst.NR) <= inst.C[t]*yp[1,t]
					)
	end

	#println("Finished setup constraints")

	### no initial inventory
	@constraint(model,
				zeroinvP[r=1:inst.NR],
				s0[r,0] == 0
				)
	@constraint(model,
				zeroinvW[r=1:inst.NR],
				s1[r,0] == 0
				)
	@constraint(model,
				zeroinvR[r=1:inst.NR,t=1:inst.NT],
				ssr[r,0,t] == 0
				)

	#println("Finished zero inventory constraints")

	### Balance constraints ###
	@constraint(model,
				balanceP[r=1:inst.NR, t=1:inst.NT],
				s0[r,t-1] + x0[r,t] == x1[r,t] + s0[r,t]
				)
	@constraint(model,
				balanceW[r=1:inst.NR,t=1:inst.NT],
				s1[r,t-1] + x1[r,t] == sum(xxr[r,t,k] for k in t:inst.NT) + s1[r,t]
				)

	@constraint(model,
				balanceR[r=1:inst.NR,t=1:inst.NT,k=t:inst.NT],
				ssr[r,t-1,k] + xxr[r,t,k] == sum(inst.D[r,z] for z in k:k if t==k) + sum(ssr[r,t,z] for z in k:k if t<k)
				)

######################################
# WW inequalities
######################################

	#@constraint(model,
	#			wwineq0[r=1:inst.NR,t=1:inst.NT,l=t:inst.NT],
	#			s0[r,t-1] + s1[r,t-1] + s2[r,t-1] + sum(yp[1,k]*cumdem[r,k,l] for k in t:l) >= cumdem[r,t,l]
	#			)



	#@constraint(model,
	#			wwineq01[r=1:inst.NR,t=1:inst.NT-1,j=t:inst.NT-1,l=j+1:inst.NT],
	#			s0[r,t-1] + s1[r,t-1] + s2[r,t-1] + sum(yp[1,k]*cumdem[r,k,l] for k in t:j) + sum(yw[inst.DeltamR[r],k]*cumdem[r,k,l] for k in j+1:l) >= cumdem[r,t,l]
	#			)


	### Capacity constraints ###
#	@constraint(model, capacity[i=1:inst.NP, t=1:inst.NT], x0[i,t] <= min(inst.DP[i,t],inst.C[t])*yp[i,t])

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
		write(f,";$(params.form);$bestbound;$bestsol;$gap;$time;$numnodes \n")
	end

	#if params.printsol == 1
		#printthreelevelFormulationSolution(inst,xp,xr,xw,yp,yr,yw,sp,sr,sw)
	#end

end #function enhancedthreelevelFormulation()

end
