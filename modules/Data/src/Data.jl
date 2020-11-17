module Data

struct InstanceData
	NR::Int #number of retailers
	NT::Int #number of periods
	NI::Int #number of the instance
	NP::Int #number of plants
	NW::Int #number of warehouses
	#NS::Int #number of sucessor

	# For the production plant:
	HCP 	#holding_cost_for_plant
	SCP 	#setup costs for each period
	DP		# demand for each period

	# For each warehouse:
	HCW 	# holding_cost_of_warehouse
	SCW  	# setup costs for each period
	DW		#demand for each period
	# For each retailer:
	HCR 	# holding_cost_of_retailer
	SCR 	# setup costs for each period
	D  		# demand for each period

	# Capacity
	C

	#Set delta
	DeltaW # a vector for the production and for the warehouses
	DeltamR # an integer for each retailer
	cumdem
end

export InstanceData, readData

function readData(instanceFile,params)

	println("Running Data.readData with file $(instanceFile)")
	file = open(instanceFile)
	fileText = read(file, String)
	tokens = split(fileText) #tokens will have all the tokens of the input file in a single vector. We will get the input token by token

	# read the problem's dimensions NR, NT, NI, NP and NW
	aux = 1
	NR = parse(Int,tokens[aux])
	aux = aux+1
	NT = parse(Int,tokens[aux])
	aux = aux+1
	NI = parse(Int,tokens[aux])
	aux = aux+1
    NP = parse(Int,tokens[aux])
	aux = aux+1
    NW = parse(Int,tokens[aux])

	#print("NR = $(NR) NT = $(NT) NI = $(NI) NP = $(NP) NW = $(NW) \n")

	#resize data structures according to NI, NT and NP
	HCP = zeros(Float64,NP)	# holding_cost_of_plant
	SCP = zeros(Int,NP,NT)	# setup costs for each period
	SomaDP = zeros(Int,NP,NT)	# setup costs for each period
	HCW = zeros(Float64,NW)	# holding_cost_of_warehouse
	SCW = zeros(Int,NW,NT)	# setup costs for each period
	HCR = zeros(Float64,NR) # holding_cost_of_retailer
  	SCR = zeros(Int,NR,NT) 	# setup costs for each period
  	D = zeros(Int,NR,NT)  	# retailer demand for each period
  	DP = zeros(Int,NP,NT)  	# plant demand for each period
	DW = zeros(Int,NR,NT) #retailer demand for each period
	C = zeros(Float64, NT)	# capacity for each period

	for p in 1:NP
		aux = aux+1
		#plant number in instance ignored
		aux = aux+1
		HCP[p] = parse(Float64,tokens[aux])

		for q in 1:NT
			aux = aux + 1
			SCP[p,q] = parse(Int,tokens[aux])
		end
	end

	for p in 1:NW
		aux = aux + 1
		#string w? in instance ignored
		aux = aux + 1
		HCW[p] = parse(Float64,tokens[aux])
		for q in 1:NT
			aux = aux + 1
			SCW[p,q] = parse(Int,tokens[aux])
		end
	end

	for p in 1:NR
		aux = aux + 1
		#string r? in instance ignored
		aux = aux + 1
		HCR[p] = parse(Float64,tokens[aux])
		for q in 1:NT
			aux = aux + 1
			SCR[p,q] = parse(Int,tokens[aux])
		end
		for q in 1:NT
			aux = aux + 1
			D[p,q] = parse(Int,tokens[aux])
		end
	end

	close(file)

	for t in 1:NT
		for r=1:NR
			DP[1,t] += D[r,t]
		end
	end

	aux = sum(D[:,:])
	for t in 1:NT
		C[t] = params.capacity*aux/NT
	end

	DeltaW,DeltamR = calculateDeltas(NW,NR,params.balanced)

	for w in 1:NW
		for t in 1:NT
			for k in 1:length(DeltaW[w])
				DW[w,t] += D[DeltaW[w][k],t]
			end
		end
	end

	open("saida.txt","a") do f
		write(f,"$instanceFile")
	end

	cumdem = zeros(Int,NR,NT,NT)

	for r in 1:NR, t in 1:NT
		cumdem[r,t,t] = D[r,t]
		for k in t+1:NT
			cumdem[r,t,k] = cumdem[r,t,k-1] + D[r,k]
		end
	end

	inst = InstanceData(NR,NT,NI,NP,NW,HCP,SCP,DP,HCW,SCW,DW,HCR,SCR,D,C,DeltaW,DeltamR,cumdem)

	return inst

end


function calculateDeltas(NW,NR,balanced)

	DeltaW = Vector{Vector{Int}}(undef,NW)
	DeltamR = zeros(Int,NR)

	for w in 1:NW
		DeltaW[w] = []
	end

	retperwar = calculateretperwar(NW,NR,balanced)

	r = 0
	for w in 1:NW
		for k in 1:retperwar[w]
			r+=1
			push!(DeltaW[w],r)
			DeltamR[r] = w
		end
	end

	#println("DeltaW: ")
	#for w in 1:NW
	#	println(DeltaW[w])
	#end

	#println("DeltamR ")
	#println(DeltamR)

	return DeltaW,DeltamR

end #function calculateDeltas()


function calculateretperwar(NW,NR,balanced)

	retperwar = zeros(Int,NW) #stores the number of retailers in each warehouse

	if balanced == 1
		if mod(NR,NW) == 0
			for w in 1:NW
				retperwar[w] = NR/NW
			end
		elseif NW == 15 && NR == 50
			for w in 1:10
				retperwar[w] = 3
			end
			for w in 11:15
				retperwar[w] = 4
			end
		elseif NW == 15 && NR == 100
			for w in 1:5
				retperwar[w] = 6
			end
			for w in 6:15
				retperwar[w] = 7
			end
		elseif NW==15 && NR == 200
			for w in 1:10
				retperwar[w] = 14
			end
			for w in 11:15
				retperwar[w] = 12
			end
		elseif NW == 20 && NR == 50
			for w in 1:10
				retperwar[w] = 3
			end
			for w in 11:20
				retperwar[w] = 2
			end
		end

	elseif balanced == 0
		if NW == 5 && NR == 50
			retperwar[1] = 40
			for w in 2:3
				retperwar[w] = 3
			end
			for w in 4:5
				retperwar[w] = 2
			end
		elseif NW == 5 && NR == 100
			retperwar[1] = 80
			for w in 2:5
				retperwar[w] = 5
			end
		elseif NW == 5 && NR == 200
			retperwar[1] = 160
			for w in 2:5
				retperwar[w] = 10
			end
		elseif NW == 10 && NR == 50
			for w in 1:2
				retperwar[w] = 17
			end
			for w in 3:10
				retperwar[w] = 2
			end
		elseif NW == 10 && NR == 100
			for w in 1:2
				retperwar[w] = 38
			end
			for w in 3:10
				retperwar[w] = 3
			end
		elseif NW == 10 && NR == 200
			for w in 1:2
				retperwar[w] = 80
			end
			for w in 3:10
				retperwar[w] = 5
			end
		elseif NW == 15 && NR == 50
			for w in 1:2
				retperwar[w] = 9
			end
			retperwar[3] = 8
			for w in 4:15
				retperwar[w] = 2
			end
		elseif NW == 15 && NR == 100
			for w in 1:2
				retperwar[w] = 25
			end
			retperwar[3] = 26
			for w in 4:15
				retperwar[w] = 2
			end
		elseif NW == 15 && NR == 200
			for w in 1:2
				retperwar[w] = 54
			end
			retperwar[3] = 56
			for w in 4:15
				retperwar[w] = 3
			end
		elseif NW == 20 && NR == 50
			for w in 1:2
				retperwar[w] = 5
			end
			for w in 3:4
				retperwar[w] = 4
			end
			for w in 5:20
				retperwar[w] = 2
			end
		elseif NW == 20 && NR == 100
			for w in 1:4
				retperwar[w] = 17
			end
			for w in 5:20
				retperwar[w] = 2
			end
		elseif NW == 20 && NR == 200
			for w in 1:4
				retperwar[w] = 38
			end
			for w in 5:20
				retperwar[w] = 3
			end

		end

	end

#	println("retperwar = ",retperwar)

	return retperwar

end #function calculateretperwar(NW,NR,balanced)

end
