module Heuristics

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


function partialmulticommodityRounding(inst::InstanceData, params::ParameterData)
	println("Running NewFormulations.partialmulticommodityFormulation")
TK=ceil(inst.NT/2)

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
	println("Percentage of fixed = ",100*(totalFixed/totalVariables))

for r=1:inst.NR,t=1:inst.NT,k=t+TK:inst.NT
	JuMP.fix(xxp[r,t,k],0)
	JuMP.fix(xxw[r,t,k], 0)
	JuMP.fix(xxr[r,t,k], 0)
end

println("Defined variables")
	### Objective function ###
	@objective(model, Min,
		sum(inst.SCP[p,t]*yp[p,t] for p=1:inst.NP, t=1:inst.NT) + sum(inst.HCP[1]*ssp[r,k,t] for r=1:inst.NR,k=1:inst.NT-1,t= k+1:inst.NT)
		+ sum(inst.SCW[w,t]*yw[w,t]  for w=1:inst.NW, t=1:inst.NT) + sum(inst.HCW[inst.DeltamR[r]]*ssw[r,k,t] for r=1:inst.NR,k=1:inst.NT-1,t= k+1:inst.NT)
		+ sum(inst.SCR[r,t]*yr[r,t]  for r=1:inst.NR, t=1:inst.NT) + sum(inst.HCR[r]*ssr[r,k,t] for r=1:inst.NR,k=1:inst.NT-1,t= k+1:inst.NT)
	)

println("Finished objective function")

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

	println("Finished setup constraints")

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

println("Finished zero inventory constraints")

	### Balance constraints ###
	@constraint(model,
				balanceP[r=1:inst.NR, t=1:inst.NT, k=t:inst.NT],
				ssp[r,t-1,k] + xxp[r,t,k] == xxw[r,t,k] + sum(ssp[r,t,z] for z in k:k if t<z)
				)
				println("h")
	@constraint(model,
				balanceW[r=1:inst.NR,t=1:inst.NT,k=t:inst.NT],
				ssw[r,t-1,k] + xxw[r,t,k] == xxr[r,t,k] + sum(ssw[r,t,z] for z in k:k if t<z)
				)
				println("h1")
	@constraint(model,
				balanceR[r=1:inst.NR,t=1:inst.NT,k=t:inst.NT],
				ssr[r,t-1,k] + xxr[r,t,k] == sum(inst.D[r,z] for z in k:k if t==k) + sum(ssr[r,t,z] for z in k:k if t<k)
				)

	#writeLP(model,"modelo.lp",genericnames=false)

	for t in 1:inst.NT
		setcategory(yp[1,t],:Cont)
		for r in 1:inst.NR
			setcategory(yr[r,t],:Cont)
		end
		for w in 1:inst.NW
			setcategory(yw[w,t],:Cont)
		end
	end

	t1 = time_ns()
	solve(model)

	bestsol = getobjectivevalue(model)
	#bestbound = getobjbound(model)

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

	ypsol = zeros(inst.NP,inst.NT)
	ywsol = zeros(inst.NW,inst.NT)
	yrsol = zeros(inst.NR,inst.NT)

	epsilon = 0.00000001
	for t = 1:inst.NT
		if(getvalue(yp[1,t])>epsilon)
			println("yp[1,$(t)] = ",getvalue(yp[1,t]))
		end
		if getvalue(yp[1,t]) > 0.5
			ypsol[1,t] = 1
		else
			ypsol[1,t] = 0
		end
	end
	for w=1:inst.NW, t=1:inst.NT
		if getvalue(yw[w,t])>epsilon
			println("yw[$(w),$(t)] = ",getvalue(yw[w,t]))
		end
		if getvalue(yw[w,t]) > 0.3
			ywsol[w,t] = 1
		else
			ywsol[w,t] = 0
		end
	end
	for r=1:inst.NR, t=1:inst.NT
		if getvalue(yr[r,t])>epsilon
			println("yr[$(r),$(t)] = ",getvalue(yr[r,t]))
		end
		if getvalue(yr[r,t]) > epsilon
			yrsol[r,t] = 1
		else
			yrsol[r,t] = 0
		end
	end

	for t in 1:inst.NT
		setcategory(yp[1,t],:Bin)
		for r in 1:inst.NR
			setcategory(yr[r,t],:Bin)
		end
		for w in 1:inst.NW
			setcategory(yw[w,t],:Bin)
		end
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

	solve(model)

	t2 = time_ns()
	elapsedtime = (t2-t1)/1.0e9

	bestsol = getobjectivevalue(model)
	println("bestsol = ", bestsol)
	println("Elapsed = ", elapsedtime)

	open("saida.txt","a") do f
		write(f,";$bestsol;$elapsedtime \n")
	end

#	if params.printsol == 1
#		printStandardFormulationSolution(inst,xp,xr,xw,yp,yr,yw,sp,sr,sw)
#	end

end #function multicommodityFormulation()

end
