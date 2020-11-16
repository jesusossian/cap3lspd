#!/bin/bash
version=cap2 # cap2 cap175 cap150   
style=bal # unbal bal(balanced) 
for model in esn #std mc 3level partialmc echelon esn estp esls
do
    for N in 50 #50 100 200
    do
	    for T in 15 #15 30
	    do
		    for P in 1
	    	do
		    	for W in 5 #5 10 15 20
		    	do
		    		for V in DD_SF #SD_SF SD_DF DD_SF DD_DF
		    		do
		    			echo "instance;bestboud;opt;gap;t(s);nodes" >> saida.txt
		    			for id in $(seq 5)
		    			do
		    		   		julia threeplsp.jl --inst instances/N${N}T${T}/N${N}T${T}P${P}W${W}${V}${id}.dat --form ${model} >> report/${version}_${style}_${model}_out_N${N}T${T}P${P}W${W}${V}${id}.txt
		    			done
		    		    mv saida.txt result/${version}_${style}_${model}_N${N}_T${T}_result.txt		
		    		done
		    	done
		    done
        done
    done
done
