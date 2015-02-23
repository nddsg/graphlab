library(ggplot2)
library(Rmisc)
library(Hmisc)
library(scales)
library(reshape2)
library(data.table)

dataraw <- read.table("stats.csv", header=TRUE, sep=",")
dataraw[, "percent"] <- as.numeric(gsub("%", "", as.character(dataraw[, "percent"])))
datafiltered <- dataraw[dataraw$exit == 0, -c(which(colnames(dataraw) == "exit"))]
datafiltered[, 'percent'] = datafiltered[, 'percent'] / 58 / 100

breaks <- c("gl", "nx", "glp", "gtp", "snp")
labels <- c("VCED", "APSP", "HADI", "Reverse Cuthill-McKee", "Random-BFS")
xbreaks <- 10^(1:5)
ybreaks <- 10^(-3:5)

gtype = c('er', 'er', 'sf')
gdirected = c(0, 1, 0)
gsfmodifier = 2 + ((0:10) * .1)
germodifier = 2^(1:5)
vars = c()
for (i in c(1:length(gtype))) {
  if (gtype[i] == 'er') {
    modifier = germodifier
  }
  else {
    modifier = gsfmodifier
  }
  for (mod in modifier) {
    vars = rbind(vars, c(gtype[i], gdirected[i], mod))
  }
}

outgraph <- function(data, filterby, type, directed, modifier, title, yaxis, ybreaks, ylabel, log = T) {
  summary <- summarySE(data, measurevar=filterby, groupvars=c("program", "type", "directed", "nodes", "modifier"))
  colnames(summary)[7] <- 'val'
  g <- ggplot(summary[summary$type==type & summary$directed==directed & summary$modifier==modifier, ], aes(x=nodes, y=val, group=program, linetype=program)) +
    geom_errorbar(aes(ymin=val-se, ymax=val+se), color='grey40', width=.1, linetype=1) +
    geom_line(size=.75) +
    scale_x_log10(breaks=xbreaks, labels=comma) +
    xlab("Number of Vertices") +
    ylab(yaxis) +
    ggtitle(paste0(title, '(', toupper(type), ', ', ifelse(directed == 0, 'undirected', 'directed'), ')')) +
    scale_linetype_manual(name="Program", breaks=breaks, labels=labels, values=c("solid", "52", "1252", "dotted", "125212"))
  if (log) {
    g <- g + scale_y_log10(breaks=ybreaks, labels=ylabel)
  }
  else {
    g <- g + scale_y_continuous(breaks = ybreaks, labels=ylabel)
  }
  g <- g + theme_classic() +
    theme(legend.direction ="vertical", legend.position="bottom", legend.key.width=unit(3,"line")) +
    guides(linetype = guide_legend(nrow = 2, byrow = T))
  show (g)
  
  ggsave(paste0("/home/cpennycu/diameter/graphs/eps/", paste(filterby, type, directed, modifier, sep="-"), ".eps"), g, width=6.5, height=5, units="in")
}

for (v in 1:length(vars)) {
  outgraph(datafiltered, 'time',        vars[v, 1], vars[v,2], vars[v,3], 'Wall Time Usage',              'Seconds of Execution', ybreaks,      comma)
  outgraph(datafiltered, 'clocktime',   vars[v, 1], vars[v,2], vars[v,3], 'User Time',                    'Seconds of Execution', ybreaks,      comma)
  outgraph(datafiltered, 'percent',     vars[v, 1], vars[v,2], vars[v,3], 'CPU Usage',                    'Percentage of Use',    ybreaks = (0:10) * .1, percent, log=F)
  ticks <- 2^(0:9)
  kbytes <- c()
  for (i in 0:2) {
    kbytes <- c(kbytes, ticks * 1024^i)
  }
  klabel <- paste0(ticks, c(rep('Kb', 10), rep('Mb', 10), rep('Gb', 10)))
  outgraph(datafiltered, 'memory',      vars[v, 1], vars[v,2], vars[v,3], 'Maximum RAM Usage',            'RAM Usage',               kbytes,  klabel)
  outgraph(datafiltered, 'voluntary',   vars[v, 1], vars[v,2], vars[v,3], 'Voluntary Context Switches',   'Number of Switches',   c(10^(0:20)),  comma)
  outgraph(datafiltered, 'involuntary', vars[v, 1], vars[v,2], vars[v,3], 'Involuntary Context Switches', 'Number of Switches',   c(10^(0:20)),  comma)
}

# Show Average Error of Pseudo Diameter methods
d <- dcast(datafiltered, type + directed + batch + nodes + modifier ~ program, value.var='diameter')
d <- d[d['gl']!= 0,]
dtemp <- c()
for (col in c('glp', 'gtp', 'snp')) {
  dfiltered <- d[d[col]!=0 & complete.cases(d[col]),]
  dfiltered[col] = (dfiltered[col] - dfiltered['gl']) / dfiltered[col]
  summary <- summarySE(dfiltered[dfiltered$nodes>10,], measurevar=col, groupvars=c("nodes"))
  setnames(summary, col, 'val')
  dtemp <- rbind(dtemp, cbind(summary, 'program'=col))
}

g <- ggplot(dtemp, aes(x=nodes, y=val, group=program, linetype=program)) +
  geom_errorbar(aes(ymin=val-se, ymax=val+se), color='grey40', width=.1, linetype=1) +
  geom_line(size=.75) +
  scale_x_log10(breaks=xbreaks, labels=comma) +
  xlab("Number of Vertices") +
  scale_y_continuous(breaks=c(-10:10)/10, labels=percent) +
  ylab("Percent Error") +
  ggtitle('Percent Error of Pseudo Diameter Algorithms') +
  scale_linetype_manual(name="Program", breaks=breaks, labels=labels, values=c("52", "1252", "125212")) +
  theme_classic() +
  theme(legend.key.width=unit(3,"line"))
show (g)

ggsave(paste0("/home/cpennycu/diameter/graphs/pseudo-error.eps"), g, width=8, height=5, units="in")

# Percent error by graph type
title <- cbind('glp'='GraphLab Pseudo', 'gtp'='Graph-Tool Pseudo', 'snp'='Snap Pseudo')
graphtype <- cbind('er0'='ER Undirected', 'er1'='ER Directed', 'sf0'='SF')
graphs <- c()
for (col in c('glp', 'gtp', 'snp')) {
  d <- dcast(datafiltered, type + directed + batch + nodes + modifier ~ program, value.var='diameter')
  dfiltered <- d[d['gl']!= 0 & d[col]!=0 & complete.cases(d[col]),]
  dfiltered[col] = (dfiltered[col] - dfiltered['gl']) / dfiltered[col]
  summary <- summarySE(dfiltered[dfiltered$nodes>10,], measurevar=col, groupvars=c('nodes', 'type', 'directed'))
  setnames(summary, col, 'val')
  summary <- cbind(summary, 'munged'=paste0(levels(summary$type)[summary$type], summary$directed))
  
  g <- ggplot(summary, aes(x=nodes, y=val, group=munged, linetype=munged)) +
    geom_errorbar(aes(ymin=val-se, ymax=val+se), color='grey40', width=.1, linetype=1) +
    geom_line(size=.75) +
    scale_x_log10(breaks=xbreaks, labels=comma) +
    xlab("Number of Vertices") +
    scale_y_continuous(breaks=c(-10:10)/10, labels=percent, limits = c(-.4, .7)) +
    ylab("Percent Error") +
    ggtitle(paste0(title[1,col], ' Error by Graph Type')) +
    scale_linetype_manual(name="Graph Type", breaks=colnames(graphtype), labels=graphtype[1,], values=c("solid", "5252", "12")) +
    theme_classic() +
    theme(legend.direction ="vertical", legend.position="bottom", legend.key.width=unit(3,"line")) +
    guides(linetype = guide_legend(nrow = 1, byrow = T))
  graphs <- rbind(graphs, g)
  show(g)
  ggsave(paste0("/home/cpennycu/diameter/graphs/algoerror",col,".eps"), g, width=5, height=5, units="in")
}
png("/home/cpennycu/diameter/graphs/algoerror.png")
multiplot(plotlist=graphs, cols=3)
dev.off()

# Percent error by graph type (all on one graph)
d <- dcast(datafiltered, type + directed + batch + nodes + modifier ~ program, value.var='diameter')
dtemp <- c()
for (col in c('glp', 'gtp', 'snp')) {
  dfiltered <- d[d$nodes>10 & d['gl']!=0 & d[col]!=0 & complete.cases(d[col]),]
  dfiltered[col] = (dfiltered[col] - dfiltered['gl']) / dfiltered[col]
  summary <- summarySE(dfiltered, measurevar=col, groupvars=c('nodes', 'type', 'directed'))
  setnames(summary, col, 'val')
  dtemp <- rbind(dtemp, cbind(summary, 'program'=col))
}
dtemp <- cbind(dtemp, 'munged'=paste0(dtemp$program, dtemp$type, dtemp$directed))

g <- ggplot(dtemp, aes(x=nodes, y=val, group=munged, linetype=munged)) +
  geom_errorbar(aes(ymin=val-se, ymax=val+se), color='grey40', width=.1, linetype=1) +
  geom_line(size=.75) +
  scale_x_log10(breaks=xbreaks, labels=comma) +
  xlab("Number of Vertices") +
  scale_y_continuous(breaks=c(-10:10)/10, labels=percent) +
  ylab("Percent Error") +
  ggtitle('Percent Error by Graph Type') +
  theme_classic() +
  theme(legend.key.width=unit(3,"line"))
show (g)

