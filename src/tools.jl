"""
This file contains helpful tools for spectral analysis such as a combination difference checker
"""

function TriangleFinder(linds)
   #finds all the combination differences triangles
      
   initial = linds[:,1:3] #initial quantum states
   final = linds[:,4:6] #final quantum states
   
   initialtog = fill("0",length(initial[:,1]))
   
   for i in 1:length(initial[:,1]) #smashes initial quantum states together
     initialtog[i] = string(initial[i,1], initial[i,2], initial[i,3])
   end
   
   finaltog = fill("0",length(final[:,1]))
   
   for i in 1:length(final[:,1]) #smashes final quantum states together
     finaltog[i] = string(final[i,1], final[i,2], final[i,3])
   end
   
   array = zeros(Int,length(initial)^2,3) #zeroes array that should have more lines than possible triangles
   j = 1 #index initializer
   
   for c in 1:length(initialtog) #for a given initial state
     initmatch = findall(isequal(initialtog[c]),initialtog) #find all the lines with the initial state (pt A)
     finalmatch = finaltog[initmatch] #find all the final states that match (the B in AB lines)
   
     if length(initmatch) == 1 #if there's only one match
     else #if not there's a triangle
         for b in 1:length(initmatch) #for each of those lines with that initial state
             finalfinal = findall(isequal(finalmatch[b]),initialtog) #find the lines with the B as their initial state
             finalfinalmatch = finaltog[finalfinal] #the lower states of *those* lines (C in BC lines)
             for i in 1:length(finalmatch) #for all those AB lines
                 tempval = finalmatch[i]
                 for a in 1:length(finalfinalmatch) #and for all those C's
                     if finalfinalmatch[a] == tempval #is there an AC line?
                     array[j,1] = initmatch[b]
                     array[j,2] = initmatch[i]
                     array[j,3] = finalfinal[a]
                     j += 1
                     else
                     end
                 end
             end
         end
     end
   end
   
   triangles0 = unique(array, dims = 1) #sort out unique triangles only
   triangles = triangles0[1:end-1,:] #get rid of the stray zero at the end
   
   return triangles

end

function TriangleTester(triangles, ofreqs)
      #makes sure all the triangles that were previously found sum to zero
   
   
   for i in 1:length(triangles[:,1])
       A = ofreqs[triangles[i,1]]
       B = ofreqs[triangles[i,2]]
       C = ofreqs[triangles[i,3]]
       line = [A,B,C] #writes the three energies to an array
       max = findmax(line)
       line[max[2]] *= -1 #negates the smallest one
       if sum(line) >= 0.1 #checks that they sum to zero
           println("The transitions on lines ", triangles[i,1], ", ", triangles[i,2], ", and ", triangles[i,3], " don't form a triangle.")
       else
       end
   end
   
end
