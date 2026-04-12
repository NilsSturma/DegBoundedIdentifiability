
# Read data
#str = readLines("results/4-nodes.m2")
str = readLines("results/random-n10-0.2-10s.m2")
#str = readLines("results/5-HTC-inconclusive.m2")

# Clean 
str = gsub("matrix ", "", str)
str = gsub("[{}]", "", str)
str = as.numeric(strsplit(str, ", ")[[1]])
df <- as.data.frame(matrix(str, ncol = 6, byrow = TRUE))
colnames(df) <- c("GP", "DegBd", "GPdeg", "DegBddeg", "timeGP", "timeDegBd")
head(df)

# Check whether both method give the same result on all graphhs
all(df[df$GP==1, "DegBd"] == 1)
all(df[df$DegBd==1, "GP"] == 1) 

nrow(df[df$GP==1, ])
nrow(df[df$DegBd==1, ])

# Select only graphs which are identifiable according to both algorithms
dfID = df[df$GP==1, ] 
nrow(dfID)
tail(dfID)

##############
# Plot times #
##############
summary(dfID[, c("timeGP", "timeDegBd")])

pdf("results/boxplot.pdf", width = 5, height = 6)
boxplot(dfID[, c("timeGP", "timeDegBd")],
        main = "",
        ylab = "Seconds",
        col = c("lightblue", "lightgreen"),
        names = c("GP", "DegBd"), 
        log = "y"
       )
dev.off()

# Histograms
all_vals <- c(dfID$timeGP, dfID$timeDegBd)
nbins = 50
breaks <- seq(min(log(all_vals)), max(log(all_vals)), length.out = nbins)
#breaks <- exp(log_breaks)

hist(log(dfID$timeGP), breaks = breaks, col = rgb(0, 0, 1, 0.4), main = "", xlab="", xaxt = "n")
hist(log(dfID$timeDegBd),
     breaks = breaks,
     col = rgb(1, 0, 0, 0.4),   # transparent red 
     add = TRUE)

# Add a nice log x-axis
axis_ticks <- round(exp(floor(log(min(all_vals))) : ceiling(log(max(all_vals)))),2)
axis(1, at = (floor(log(min(all_vals))) : ceiling(log(max(all_vals)))), labels = axis_ticks)



########################
# Plot maximal degrees #
########################

# Histograms
all_vals <- c(dfID$GPdeg, dfID$DegBddeg)
#nbins = 50
breaks <- seq(min(all_vals), max(all_vals))+0.5
breaks <- c(min(breaks) - 1, breaks)
#breaks <- exp(log_breaks)

hist(dfID$GPdeg, breaks = breaks, col = rgb(0, 0, 1, 0.4), main = "", xlab="", ylim = c(0,300))
hist(dfID$DegBddeg,
     breaks = breaks,
     col = rgb(1, 0, 0, 0.4),   # transparent red 
     add = TRUE)

# Add a nice log x-axis
axis_ticks <- round(exp(floor(log(min(all_vals))) : ceiling(log(max(all_vals)))),2)
axis(1, at = (floor(log(min(all_vals))) : ceiling(log(max(all_vals)))), labels = axis_ticks)

###################################
# Create table by number of edges #
###################################
library(rjson)
graphs <- fromJSON(file = "graphs/random-n10-0.2.json")
#graphs <- fromJSON(file = "graphs/acyclic-mixed-graphs-4-nodes.json")

df$nedges = sapply(seq_len(nrow(df)), function(i) {graphs[[i]]$ndir+graphs[[i]]$nbidir})
#df$nedges = sapply(seq_len(nrow(df)), function(i) {sum(graphs[[i]]$D)+(sum(graphs[[i]]$B)/2)})
unique(df$nedges)
length(unique(df$nedges))
max(unique(df$nedges))
min(unique(df$nedges))
head(df)

table <- matrix(0, nrow=length(unique(df$nedges)), ncol=6)
rownames(table) = unique(df$nedges)
table = table[order(as.integer(rownames(table))), ]

for (i in 1:nrow(df)){
  table[as.character(df[i, "nedges"]), 1] = table[as.character(df[i, "nedges"]), 1] + 1
  if (df[i,"GP"]==1){
    table[as.character(df[i, "nedges"]), 2] = table[as.character(df[i, "nedges"]), 2] + 1
    table[as.character(df[i, "nedges"]), 4] = table[as.character(df[i, "nedges"]), 4] + df[i,"timeGP"]
    #table[as.character(df[i, "nedges"]), 7] = table[as.character(df[i, "nedges"]), 7] + df[i,"GPdeg"]
  }
  if (df[i,"DegBd"]==1){
    table[as.character(df[i, "nedges"]), 3] = table[as.character(df[i, "nedges"]), 3] + 1
    table[as.character(df[i, "nedges"]), 5] = table[as.character(df[i, "nedges"]), 5] + df[i,"timeDegBd"]
    #table[as.character(df[i, "nedges"]), 8] = table[as.character(df[i, "nedges"]), 8] + df[i,"DegBddeg"]
  }
  if (graphs[[i]]$htcID){
    table[as.character(df[i, "nedges"]), 6] = table[as.character(df[i, "nedges"]), 6] + 1
  }
}
table[,4] = ifelse(table[, 2]==0, NA, table[,4] / table[, 2])
#table[,7] = ifelse(table[, 2]==0, NA, table[,7] / table[, 2])
table[,5] = ifelse(table[, 3]==0, NA, table[,5] / table[, 3])
#table[,8] = ifelse(table[, 3]==0, NA, table[,8] / table[, 3])

colnames(table) <- c("Total", "GP-ID", "DegBd-ID", "Avg. Time GP", "Avg. Time DegBd", "HTC-ID") #, "Avg. Deg. GP", "Avg. Deg. DegBD")
round(table, 2)

colSums(table)
colMeans(table, na.rm=TRUE)

write.table(table, file = "results/table-random-n10-0.2-10s.txt", sep = "\t", row.names = TRUE, col.names = TRUE)
#write.table(table, file = "results/table-4-nodes.txt", sep = "\t", row.names = TRUE, col.names = TRUE)
