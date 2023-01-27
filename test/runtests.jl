using socialpde: testabm, testpde, testensemble, testcontrol, testoptimization

cd("..")  # in order to find the output directories

testabm("4inf")
testpde("4inf")
testpde("uniform")

testensemble()
testcontrol()
testoptimization()
