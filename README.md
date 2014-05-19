TTree File Task Code Description:

This is a repository for my work on Christina's TTree File Task where we were asked to produce a series of histograms from a simulated data file of pion production. The macro reads in the simulated data file and finds all the variables held in the tree file. Each event is then looped over in a for loop where TLorentzVector is used to calculate the values used in the various histograms. The histograms are then printed, and TSPectrum is used on the pion and proton mass histograms to print out the peak's x-position in the terminal.

How to Run Pi0TreeTask.C:

Before loading the macro into ROOT you will have to go into the macro and change the location of the simulated file. After this is done simply compile the macro in ROOT and run Pi0TreeTask().
