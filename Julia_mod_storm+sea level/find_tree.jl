#called in load redistribution to find the connecting tree for a particular bus
#change B
function findtree(busind)
        B=14

        treenet=zeros(Int,B,B) # indicator * level
        treenet[busind,1]=1

        treenetind=zeros(Int,B,B)# busnum * indicator
        treenetind[1,1]=busind

        nodenum=zeros(Int,1,B)# track node number in each level
        nodenum[1]=1


        #identify the load first level connection
        mpc = loadcase("case14")
        branchnum=size(mpc["branch"],1)
        fromto=round.(Int,mpc["branch"][:,1:2]) #integer


          for level=1:B-1
                  for nodei=1:nodenum[level]
                          if treenetind[nodei,level]>0
                          for j=1:branchnum
                                if fromto[j,1]==treenetind[nodei,level]
                                   nodenum[level+1]=nodenum[level+1]+1
                                   treenet[fromto[j,2],level+1]=1
                                   treenetind[nodenum[level+1],level+1]=fromto[j,2]
                                   fromto[j,:]=[0 0]
                                end

                                if fromto[j,2]==treenetind[nodei,level]
                                    nodenum[level+1]=nodenum[level+1]+1
                                    treenet[fromto[j,1],level+1]=1
                                    treenetind[nodenum[level+1],level+1]=fromto[j,1]
                                    fromto[j,:]=[0 0]
                                end
                         end
                         end
                 end
          end

        # there are duplication need to fix
        treenetfix=zeros(Int,B,B)
        for i=1:B
                for j=1:B
                        if treenet[i,j]==1
                                treenet[i,:]=0
                                treenetfix[i,j]=1
                        end
                end
        end

        if sum(treenetfix)!=B
                println()
                println()
                println("ERROR on identifying Treenetwork for bus:",i)
        end

       return treenetfix



end
