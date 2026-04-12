library(rjson)

# Read data
str = readLines("experiments/results/random-n10-0.2-10s.m2")

# Clean 
str = gsub("matrix ", "", str)
str = gsub("[{}]", "", str)
str = as.numeric(strsplit(str, ", ")[[1]])
df <- as.data.frame(matrix(str, ncol = 6, byrow = TRUE))
colnames(df) <- c("GP", "DegBd", "GPdeg", "DegBddeg", "timeGP", "timeDegBd")

head(df)
nrow(df[df$GP==1, ])
nrow(df[df$DegBd==1, ])

# Create table by number of edges
graphs <- fromJSON(file = "graphs/random-n10-0.2.json")

df$nedges = sapply(seq_len(nrow(df)), function(i) {graphs[[i]]$ndir+graphs[[i]]$nbidir})
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
  }
  if (df[i,"DegBd"]==1){
    table[as.character(df[i, "nedges"]), 3] = table[as.character(df[i, "nedges"]), 3] + 1
    table[as.character(df[i, "nedges"]), 5] = table[as.character(df[i, "nedges"]), 5] + df[i,"timeDegBd"]
  }
  if (graphs[[i]]$htcID){
    table[as.character(df[i, "nedges"]), 6] = table[as.character(df[i, "nedges"]), 6] + 1
  }
}
table[,4] = ifelse(table[, 2]==0, NA, table[,4] / table[, 2])
table[,5] = ifelse(table[, 3]==0, NA, table[,5] / table[, 3])

colnames(table) <- c("Total", "GP-ID", "DegBd-ID", "Avg. Time GP", "Avg. Time DegBd", "HTC-ID") 
round(table, 2)

colSums(table)

write.table(table, file = "experiments/results/table-random-n10-0.2-10s.txt", sep = "\t", row.names = TRUE, col.names = TRUE)
