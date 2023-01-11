################# influencer moving to right corner
#followersum, masscorner,speedpenalty, sols, Ps = solutionfixedtargets([1.5 1.5];  p=PDEconstruct(),q=parameters_control(), r=parameters_optcont(ntarg = 1 ,start="zero"), scenario="optimalcontrol",countercontrol = "no",stubborntarget=[1.5 1.5])
#plotsnapshots(sols, Ps, [5. 7. 8.5 10.]; scenario="infmax",followercount=true)


############## influencer counteraction
targets =  [ 0.5604084588036462
0.5455196695338128
0.642358069258012
0.6527252366140318
0.7454806015269124
0.734952127347543
0.8060867535176255
0.7718781212141405]


#followersum, masscorner,speedpenalty, sols, Ps = solutionfixedtargets(targets;  p=PDEconstruct(),q=parameters_control(), r=parameters_optcont(ntarg = 4 ,start="zero"), scenario="optimalcontrol",countercontrol = "inf",stubborntarget=[1.5 1.5])
#plotsnapshots(sols, Ps, [5. 6. 7.5 10.]; scenario="counterinf",followercount=true)

########### media counteraction
targets = [  1.8106679226546316
-0.9805088659711532
 1.7642225907551847
-1.7710551914852335
 1.6338751208236737
-1.8358368944346843
 1.5683845474974514
-1.7962804090565507]

#followersum, masscorner,speedpenalty, sols, Ps = solutionfixedtargets(targets;  p=PDEconstruct(),q=parameters_control(), r=parameters_optcont(ntarg = 4 ,start="zero"), scenario="optimalcontrol",countercontrol = "med",stubborntarget=[1.5 1.5])
#plotsnapshots(sols, Ps, [5. 6. 7.5 10.];  scenario="countermed",followercount=true)