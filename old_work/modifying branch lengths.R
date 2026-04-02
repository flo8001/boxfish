setwd("/Users/grahamslater/Dropbox/boxfish/boxfish trees");

t <- read.tree("/Users/grahamslater/Dropbox/boxfish/boxfish trees/newboxcomplete.phy")[[2]];


plot(t);

t$node.labels <- NULL;

t$edge.length

nodelabels()
rm(tt)
tt <- chronopl(t, lambda = 0, age.min = c(36, 8), age.max = c(50, 40), node = c("root", 36));

#plot(tt);

#axisPhylo()

tt <- chronopl(t, lambda = 100, age.min = 50, age.max = 55, node = "root"); # root only

tt <- chronopl(t, lambda = 10, age.min = 8, age.max = 40, node = 36); # ostraciidae only

plot(tt);
axisPhylo()

nodelabels(round(branching.times(tt),2));

write.tree(tt, "USE_THIS_TINA_timeCalibratedFull.phy");
?nodelabels


d <- read.csv("/Users/grahamslater/Dropbox/boxfish/contour plot practice/allData.csv");

colnames(d);

species <- as.character(d$Species);
rm(names)
tt$tip.label