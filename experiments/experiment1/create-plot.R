# Read data
str = readLines("experiments/results/4-nodes.m2")

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

pdf("experiments/results/boxplot.pdf", width = 5, height = 6)
boxplot(dfID[, c("timeGP", "timeDegBd")],
        main = "",
        ylab = "Seconds",
        col = c("lightblue", "lightgreen"),
        names = c("GP", "DegBd"), 
        log = "y"
       )
dev.off()
