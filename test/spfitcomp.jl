#prepare var file
#cp blank.var test.var
##copy template file
##set spin degeneracy
#sed -i 's/SD/$SD/g' test.var
##set K max
##replace parameter values
#prepare int file
#cp blank.int test.inp
##set temperature
##set Fmax
##set νmax
#run simulation
#grab partion function value
#replace partition function value
#rerun simulation
#convert cat file to csv
##should be similar to my pred2csv script
#convert csv to lne
#build westerfit parameter input
#fit spcat output in westerfit
#print results summary
#rm test.*
