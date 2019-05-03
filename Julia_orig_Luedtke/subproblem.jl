function subprob(a,scenind)
# a in B*T
# used in second stage to check fw
T=3



       loadval=[0 ;21; 94; 47; 7 ;11; 0; 0; 29; 9; 3; 6; 13; 14]*ones(1,T) #T*busnum

        for Tind=1:T


            fwtemp=f2nd(a[:,Tind],loadval[:,Tind]) #simpler version of 2nd stage
				#println(" ===scenind: ",scenind," ===T: ",Tind)

				if fwtemp==0
                    println(" =====scenind: ",scenind," fail the current design at t=",Tind,"=========")
                    println()
				return 0,Tind
				end

		end
        println(" =====scenind: ",scenind," pass the current design =========")
        println()
		return 1,T
end
