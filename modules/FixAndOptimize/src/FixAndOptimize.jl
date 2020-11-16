module FixAndOptimize

using JuMP
using Gurobi
using CPLEX
using Data
using Parameters

export FixAndOptimizeStandardFormulation

function FixAndOptimizeStandardFormulation(inst::InstanceData, params,SETR,SETW,SETP,bestsol)

	fixsizefo = 5
	horsizefo = 7
	tolgapfo = 0.00001
	maxfixsizefo = fixsizefo
	maxhorsizefo = horsizefo
	maxtimefo = 60

	nbsubprob = ceil(inst.NT/fixsizefo)
	maxtimesubprob = maxtimefo/nbsubprob


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

	println("Running FixAndOptimize.FixAndOptimizeStandardFormulation")
	if params.solver == "Gurobi"
		env = Gurobi.Env()
		model = Model(solver=GurobiSolver(TimeLimit=maxtimesubprob,MIPGap=tolgapfo))
	elseif params.solver == "Cplex"
		#env = Cplex.Env()
		model = Model(solver=CplexSolver(CPX_PARAM_TILIM=maxtimesubprob,CPX_PARAM_EPGAP=tolgapfo))
	else
		println("No solver selected")
		return 0
	end

	### Defining variables ###
    @variable(model,0 <= xp[i=1:inst.NP,t=1:inst.NT] <= Inf) #sum(inst.DP[i,t] for t=1:inst.NT))
    @variable(model,0 <= xr[i=1:inst.NR,t=1:inst.NT] <= Inf) #sum(inst.DR[i,t] for t=1:inst.NT))
    @variable(model,0 <= xw[i=1:inst.NW,t=1:inst.NT] <= Inf) #sum(inst.D[i,t] for t=1:inst.NT))
    @variable(model, yp[i=1:inst.NP,t=1:inst.NT], Bin)
    @variable(model, yr[i=1:inst.NR,t=1:inst.NT], Bin)
    @variable(model, yw[i=1:inst.NW,t=1:inst.NT], Bin)
    @variable(model,0 <= sp[i=1:inst.NP,t=0:inst.NT] <= Inf) #sum(inst.DP[i,t] for t=1:inst.NT))
    @variable(model,0 <= sr[i=1:inst.NR,t=0:inst.NT] <= Inf) #sum(inst.DR[i,t] for t=1:inst.NT))
    @variable(model,0 <= sw[i=1:inst.NW,t=0:inst.NT] <= Inf)

println("Defined variables")
	### Objective function ###
	@objective(model, Min,
        sum(inst.SCP[i,t]*yp[i,t] + inst.HCP[i]*sp[i,t] for i=1:inst.NP, t=1:inst.NT)
        + sum(inst.SCR[i,t]*yr[i,t] + inst.HCR[i]*sr[i,t] for i=1:inst.NR, t=1:inst.NT)
        + sum(inst.SCW[i,t]*yw[i,t] + inst.HCW[i]*sw[i,t] for i=1:inst.NW, t=1:inst.NT)
	)

println("Finished objective function")

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

	println("Finished setup constraints")

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

println("Finished zero inventory constraints")

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





	### Fix and optimize ###
	println("\n\n\n FIX AND OPTIMIZE")
	#setparam!(env, MIPGap, 0.000001)
	k = horsizefo
	kprime = fixsizefo
	maxk = maxhorsizefo
	maxkprime = maxfixsizefo

	elapsedtime = 0

	alpha = 1
	beta = min(alpha + k - 1, inst.NT)


	for t in 1:inst.NT
		setvalue(yp[1,t],ypsol[1,t])
	end
	for w in 1:inst.NW, t in 1:inst.NT
		setvalue(yw[w,t],ywsol[w,t])
	end
	for r in 1:inst.NR, t in 1:inst.NT
		setvalue(yr[r,t],yrsol[r,t])
	end


	for it2 in kprime:maxkprime
		println(kprime,"  ",maxkprime)
		println(k," ",kprime+1," <? ",maxk)
		for it in max(k,kprime+1):maxk
			println("\n\n ########### RESTARTING FIX AND OPTIMIZE with k=$(it) and kprime=$(it2)")
			proceed = true
			while (proceed && elapsedtime + maxtimesubprob <= maxtimefo )
				proceed = false



				alpha = 1
				beta = min(alpha + k - 1, inst.NT)
				while(beta <= inst.NT && elapsedtime + maxtimesubprob <= maxtimefo)

					println("FIX AND OPTIMIZE [$(alpha), $(beta)]")

					t1 = time_ns()

					for t in 1:inst.NT
						if ypsol[1,t] == 1
							setlowerbound(yp[1,t],1)
						else
							setlowerbound(yp[1,t],0)
						end

						for w in 1:inst.NW
							if ywsol[w,t] == 1
								setlowerbound(yw[w,t],1)
							else
								setlowerbound(yw[w,t],0)
							end
						end

						for r in 1:inst.NR
							if yrsol[r,t] == 1
								setlowerbound(yr[r,t],1)
							else
								setlowerbound(yr[r,t],0)
							end
						end

					end


					for t in alpha:beta
						setlowerbound(yp[1,t],0)
						for w in 1:inst.NW
							setlowerbound(yw[w,t],0)
						end
						for r in 1:inst.NR
							setlowerbound(yr[r,t],0)
						end
					end


					solve(model)
					if getobjectivevalue(model) < bestsol - 0.01
						if getobjectivevalue(model) < bestsol - 0.5
							proceed = true
						end
						bestsol = getobjectivevalue(model)
						for t in 1:inst.NT
							if getvalue(yp[1,t]) >= 0.99
								ypsol[1,t] = 1
							else
								ypsol[1,t] = 0
							end
							for w in 1:inst.NW
								if getvalue(yw[w,t]) >= 0.99
									ywsol[w,t] = 1
								else
									ywsol[w,t] = 0
								end
							end
							for r in 1:inst.NR
								if getvalue(yr[r,t]) >= 0.99
									yrsol[r,t] = 1
								else
									yrsol[r,t] = 0
								end
							end
						end

					end




					alpha = alpha + kprime
					if beta == inst.NT
						beta = inst.NT+1
					else
						beta = min(alpha + k -1,inst.NT)
					end

					t2 = time_ns()
					elapsedtime += (t2-t1)/1.0e9
					println("Elapsed ",elapsedtime)

				end

			end
		end
	end


	return ypsol,ywsol,yrsol

end




function lpbasedimprovement(inst::InstanceData, params,SETR,SETW,SETP,bestsol)


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

	println("Running FixAndOptimize.lpbasedimprovement")
	if params.solver == "Gurobi"
		env = Gurobi.Env()
		model = Model(solver=GurobiSolver(TimeLimit=params.maxtime))
	elseif params.solver == "Cplex"
		#env = Cplex.Env()
		model = Model(solver=CplexSolver(CPX_PARAM_TILIM=params.maxtime))
	else
		println("No solver selected")
		return 0
	end

	### Defining variables ###
    @variable(model,0 <= xp[i=1:inst.NP,t=1:inst.NT] <= Inf) #sum(inst.DP[i,t] for t=1:inst.NT))
    @variable(model,0 <= xr[i=1:inst.NR,t=1:inst.NT] <= Inf) #sum(inst.DR[i,t] for t=1:inst.NT))
    @variable(model,0 <= xw[i=1:inst.NW,t=1:inst.NT] <= Inf) #sum(inst.D[i,t] for t=1:inst.NT))
    @variable(model, 0 <=  yp[i=1:inst.NP,t=1:inst.NT] <= 1)
    @variable(model, 0 <=yr[i=1:inst.NR,t=1:inst.NT] <=1)
    @variable(model, 0 <=  yw[i=1:inst.NW,t=1:inst.NT]<=1)
    @variable(model,0 <= sp[i=1:inst.NP,t=0:inst.NT] <= Inf) #sum(inst.DP[i,t] for t=1:inst.NT))
    @variable(model,0 <= sr[i=1:inst.NR,t=0:inst.NT] <= Inf) #sum(inst.DR[i,t] for t=1:inst.NT))
    @variable(model,0 <= sw[i=1:inst.NW,t=0:inst.NT] <= Inf)

println("Defined variables")
	### Objective function ###
	@objective(model, Min,
        sum(inst.SCP[i,t]*yp[i,t] + inst.HCP[i]*sp[i,t] for i=1:inst.NP, t=1:inst.NT)
        + sum(inst.SCR[i,t]*yr[i,t] + inst.HCR[i]*sr[i,t] for i=1:inst.NR, t=1:inst.NT)
        + sum(inst.SCW[i,t]*yw[i,t] + inst.HCW[i]*sw[i,t] for i=1:inst.NW, t=1:inst.NT)
	)

println("Finished objective function")

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

	println("Finished setup constraints")

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

println("Finished zero inventory constraints")

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





	for t in 1:inst.NT
		@constraint(model,yp[1,t]==ypsol[1,t])
	end
	for w in 1:inst.NW, t in 1:inst.NT
		@constraint(model,yw[w,t]==ywsol[w,t])
	end
	for r in 1:inst.NR, t in 1:inst.NT
		@constraint(model,yr[r,t]==yrsol[r,t])
	end


	solve(model)


	plantcost = 0
	for t in 1:inst.NT
		plantcost += inst.SCP[1,t]*getvalue(yp[1,t]) + inst.HCP[1]*getvalue(sp[1,t])
 	end
	println("Plant cost = $(plantcost)")

	warecosts = 0
	for w in 1:inst.NW, t in 1:inst.NT
		warecosts += inst.SCW[w,t]*getvalue(yw[w,t]) + inst.HCW[w]*getvalue(sw[w,t])
	end
	println("Warehouses cost = $(warecosts)")

	retcosts = 0
	for r in 1:inst.NR, t in 1:inst.NT
		retcosts += inst.SCR[r,t]*getvalue(yr[r,t]) + inst.HCR[r]*getvalue(sr[r,t])
	end
	println("Retailers cost = $(retcosts)")

	total = plantcost + warecosts + retcosts

	println("Total cost = $(total)")

	println("Test unused setups setups")
	for t in 1:inst.NT
		if getvalue(yp[1,t]) > 0.9999 && getvalue(xp[1,t]) < 0.0001
			println("plant, period $(t)")
		end
 	end

	for w in 1:inst.NW, t in 1:inst.NT
		if getvalue(yw[w,t]) > 0.9999 && getvalue(xw[w,t]) < 0.0001
			println("warehouse $(w), period $(t)")
		end
	end

	for r in 1:inst.NR, t in 1:inst.NT
		if getvalue(yr[r,t]) > 0.9999 && getvalue(xr[r,t]) < 0.0001
			println("retailer $(r), period $(t)")
		end
	end

	return ypsol,ywsol,yrsol

end




function FixAndOptimizeEchelonStockFormulation(inst::InstanceData, params,SETR,SETW,SETP,bestsol)

	fixsizefo = 4
	horsizefo = 6
	tolgapfo = 0.0001
	maxfixsizefo = fixsizefo + 1
	maxhorsizefo = horsizefo + 1

	maxtimefo = 120

	nbsubprob = ceil(inst.NT/fixsizefo)
	maxtimesubprob = floor((maxtimefo-5)/nbsubprob)


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

	println("Running FixAndOptimize.FixAndOptimizeStandardFormulation")
	if params.solver == "Gurobi"
		env = Gurobi.Env()
		model = Model(solver=GurobiSolver(TimeLimit=maxtimesubprob,MIPGap=tolgapfo))
	elseif params.solver == "Cplex"
		#env = Cplex.Env()
		model = Model(solver=CplexSolver(CPX_PARAM_TILIM=maxtimesubprob,CPX_PARAM_EPGAP=tolgapfo))
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





	### Fix and optimize ###
	println("\n\n\n FIX AND OPTIMIZE")
	#setparam!(env, MIPGap, 0.000001)
	k = horsizefo
	kprime = fixsizefo
	maxk = maxhorsizefo
	maxkprime = maxfixsizefo

	elapsedtime = 0

	alpha = 1
	beta = min(alpha + k - 1, inst.NT)


	for t in 1:inst.NT
		setvalue(yp[1,t],ypsol[1,t])
	end
	for w in 1:inst.NW, t in 1:inst.NT
		setvalue(yw[w,t],ywsol[w,t])
	end
	for r in 1:inst.NR, t in 1:inst.NT
		setvalue(yr[r,t],yrsol[r,t])
	end


	for it2 in kprime:maxkprime
		println(kprime,"  ",maxkprime)
		println(k," ",kprime+1," <? ",maxk)
		for it in max(k,kprime+1):maxk
			println("\n\n ########### RESTARTING FIX AND OPTIMIZE with k=$(it) and kprime=$(it2)")
			proceed = true
			while (proceed && elapsedtime + maxtimesubprob <= maxtimefo )
				proceed = false



				alpha = 1
				beta = min(alpha + k - 1, inst.NT)
				println("sum: $(elapsedtime + maxtimesubprob)")
				println("max: $(maxtimefo)")
				while(beta <= inst.NT && elapsedtime + maxtimesubprob <= maxtimefo)

					println("FIX AND OPTIMIZE [$(alpha), $(beta)]")

					t1 = time_ns()

					for t in 1:inst.NT
						if ypsol[1,t] == 1
							setlowerbound(yp[1,t],1)
						else
							setlowerbound(yp[1,t],0)
						end

						for w in 1:inst.NW
							if ywsol[w,t] == 1
								setlowerbound(yw[w,t],1)
							else
								setlowerbound(yw[w,t],0)
							end
						end

						for r in 1:inst.NR
							if yrsol[r,t] == 1
								setlowerbound(yr[r,t],1)
							else
								setlowerbound(yr[r,t],0)
							end
						end

					end


					for t in alpha:beta
						setlowerbound(yp[1,t],0)
						for w in 1:inst.NW
							setlowerbound(yw[w,t],0)
						end
						for r in 1:inst.NR
							setlowerbound(yr[r,t],0)
						end
					end


					solve(model)
					if getobjectivevalue(model) < bestsol - 0.01
						if getobjectivevalue(model) < bestsol - 0.5
							proceed = true
						end
						bestsol = getobjectivevalue(model)
						for t in 1:inst.NT
							if getvalue(yp[1,t]) >= 0.99
								ypsol[1,t] = 1
							else
								ypsol[1,t] = 0
							end
							for w in 1:inst.NW
								if getvalue(yw[w,t]) >= 0.99
									ywsol[w,t] = 1
								else
									ywsol[w,t] = 0
								end
							end
							for r in 1:inst.NR
								if getvalue(yr[r,t]) >= 0.99
									yrsol[r,t] = 1
								else
									yrsol[r,t] = 0
								end
							end
						end

					end




					alpha = alpha + kprime
					if beta == inst.NT
						beta = inst.NT+1
					else
						beta = min(alpha + k -1,inst.NT)
					end

					t2 = time_ns()
					elapsedtime += (t2-t1)/1.0e9
					println("Elapsed ",elapsedtime)

					println("sum: $(elapsedtime + maxtimesubprob)")
					println("max: $(maxtimefo)")

				end

			end
		end
	end


	return ypsol,ywsol,yrsol

end





function FixAndOptimizeEchelonStockFormulation2(inst::InstanceData, params,SETR,SETW,SETP,bestsol)


	toltime = 0

	fixsizefo = 4
	horsizefo = fixsizefo+2
	tolgapfo = 0.0001
	#maxfixsizefo = fixsizefo + 2
	#maxhorsizefo = horsizefo + 2

	maxfixsizefo = inst.NT - 2
	maxhorsizefo = inst.NT

	if params.maxtime < 3600
		maxtimefo = params.maxtime
	else
		maxtimefo = 300
	end

	nbsubprob = ceil(inst.NT/fixsizefo)





	maxtimesubprob = floor(maxtimefo/(2*nbsubprob)) #try to guarantee at least two passes


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

	println("Running FixAndOptimize.FixAndOptimizeStandardFormulation")
	if params.solver == "Gurobi"
		env = Gurobi.Env()
		model = Model(solver=GurobiSolver(TimeLimit=maxtimesubprob,MIPGap=tolgapfo))
	elseif params.solver == "Cplex"
		#env = Cplex.Env()
		model = Model(solver=CplexSolver(CPX_PARAM_TILIM=maxtimesubprob,CPX_PARAM_EPGAP=tolgapfo))
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







	### Fix and optimize ###
	println("\n\n\n FIX AND OPTIMIZE")
	#setparam!(env, MIPGap, 0.000001)
	k = horsizefo
	kprime = fixsizefo
	maxk = maxhorsizefo
	maxkprime = maxfixsizefo

	elapsedtime = 0

	alpha = 1
	beta = min(alpha + k - 1, inst.NT)


	for t in 1:inst.NT
		setvalue(yp[1,t],ypsol[1,t])
	end
	for w in 1:inst.NW, t in 1:inst.NT
		setvalue(yw[w,t],ywsol[w,t])
	end
	for r in 1:inst.NR, t in 1:inst.NT
		setvalue(yr[r,t],yrsol[r,t])
	end

	tolimprov = 5

	it2 = kprime
	rounds = 0

	keepfo = true
	previousbestsol = Inf
	while keepfo == true
		keepfo = false

		rounds += 1

		#break
		println(kprime,"  ",maxkprime)
		println(k," ",kprime+1," <? ",maxk)
		it = it2+2
		#for it in max(k,kprime+1):maxk
			println("\n\n ########### RESTARTING FIX AND OPTIMIZE with k=$(it) and kprime=$(it2)")
			proceed = true
			while (proceed && elapsedtime + maxtimesubprob <= maxtimefo+toltime )
				proceed = false



				alpha = 1
				beta = min(alpha + k - 1, inst.NT)
				println("sum: $(elapsedtime + maxtimesubprob)")
				println("max: $(maxtimefo)")
				while(beta <= inst.NT && elapsedtime + maxtimesubprob <= maxtimefo+toltime)

					println("FIX AND OPTIMIZE [$(alpha), $(beta)]")

					t1 = time_ns()

					println("bestsolbefore = ",bestsol)
					bestsol = FixAndOptimizeInterval(inst,model,yp,yw,yr,ypsol,ywsol,yrsol,bestsol,alpha,beta)
					println("bestsolafter = ",bestsol)



					alpha = alpha + kprime
					if beta == inst.NT
						beta = inst.NT+1
					else
						beta = min(alpha + k -1,inst.NT)
					end

					t2 = time_ns()
					elapsedtime += (t2-t1)/1.0e9
					println("Elapsed ",elapsedtime)

					println("sum: $(elapsedtime + maxtimesubprob)")
					println("max: $(maxtimefo)")
					println("maxtimesubprob: $(maxtimesubprob)")

				end

			end

		if bestsol < previousbestsol - tolimprov
			previousbestsol = bestsol
			if elapsedtime+maxtimesubprob < maxtimefo + toltime
				keepfo = true
			end
		elseif elapsedtime+maxtimesubprob < maxtimefo + toltime && it2 < maxkprime
			keepfo = true
			it2 += 1
		end

		#end
	end


	#println("********************************")
	#println("Fix and optimize for warehouses")
	#for w in 1:inst.NW

	#	FixAndOptimizeWarehouse(inst,model,yp,yw,yr,ypsol,ywsol,yrsol,bestsol,w)

	#end

	#return ypsol,ywsol,yrsol
	return bestsol

end



function FixAndOptimizeInterval(inst,model,yp,yw,yr,ypsol,ywsol,yrsol,bestsol,alpha,beta)


	for t in 1:inst.NT
		if ypsol[1,t] == 1
			setlowerbound(yp[1,t],1)
		else
			setlowerbound(yp[1,t],0)
		end

		for w in 1:inst.NW
			if ywsol[w,t] == 1
				setlowerbound(yw[w,t],1)
			else
				setlowerbound(yw[w,t],0)
			end
		end

		for r in 1:inst.NR
			if yrsol[r,t] == 1
				setlowerbound(yr[r,t],1)
			else
				setlowerbound(yr[r,t],0)
			end
		end

	end


	for t in alpha:beta
		setlowerbound(yp[1,t],0)
		for w in 1:inst.NW
			setlowerbound(yw[w,t],0)
		end
		for r in 1:inst.NR
			setlowerbound(yr[r,t],0)
		end
	end


	solve(model)
	if getobjectivevalue(model) < bestsol - 0.01
		if getobjectivevalue(model) < bestsol - 0.5
			proceed = true
		end
		bestsol = getobjectivevalue(model)
		for t in 1:inst.NT
			if getvalue(yp[1,t]) >= 0.99
				ypsol[1,t] = 1
			else
				ypsol[1,t] = 0
			end
			for w in 1:inst.NW
				if getvalue(yw[w,t]) >= 0.99
					ywsol[w,t] = 1
				else
					ywsol[w,t] = 0
				end
			end
			for r in 1:inst.NR
				if getvalue(yr[r,t]) >= 0.99
					yrsol[r,t] = 1
				else
					yrsol[r,t] = 0
				end
			end
		end

	end

	return bestsol

end




function FixAndOptimizeWarehouse(inst,model,yp,yw,yr,ypsol,ywsol,yrsol,bestsol,ware)


	for t in 1:inst.NT
		if ypsol[1,t] == 1
			setlowerbound(yp[1,t],1)
		else
			setlowerbound(yp[1,t],0)
		end

		for w in 1:inst.NW
			if ywsol[w,t] == 1
				setlowerbound(yw[w,t],1)
			else
				setlowerbound(yw[w,t],0)
			end
		end

		for r in 1:inst.NR
			if yrsol[r,t] == 1
				setlowerbound(yr[r,t],1)
			else
				setlowerbound(yr[r,t],0)
			end
		end

	end


	for t in 1:inst.NT
		setlowerbound(yp[1,t],0)
		setlowerbound(yw[ware,t],0)
		for k in 1:length(inst.DeltaW[ware])
			r = inst.DeltaW[ware][k]
			setlowerbound(yr[r,t],0)
		end
	end


	solve(model)
	if getobjectivevalue(model) < bestsol - 0.01
		if getobjectivevalue(model) < bestsol - 0.5
			proceed = true
		end
		bestsol = getobjectivevalue(model)
		for t in 1:inst.NT
			if getvalue(yp[1,t]) >= 0.99
				ypsol[1,t] = 1
			else
				ypsol[1,t] = 0
			end
			for w in 1:inst.NW
				if getvalue(yw[w,t]) >= 0.99
					ywsol[w,t] = 1
				else
					ywsol[w,t] = 0
				end
			end
			for r in 1:inst.NR
				if getvalue(yr[r,t]) >= 0.99
					yrsol[r,t] = 1
				else
					yrsol[r,t] = 0
				end
			end
		end

	end


end





end
