using MAT

file = matopen("abm_mat/data/labmatfile.mat")
read(file, "varname") # note that this does NOT introduce a variable ``varname`` into scope
close(file)P