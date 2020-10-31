library('VennDiagram')

GFLasso <-  read.csv("C:/Users/vrrani/Documents/SNF/AllSimilarities/Ovarian/GFLasso.csv",header = TRUE)
Lasso <-  read.csv("C:/Users/vrrani/Documents/SNF/AllSimilarities/Ovarian/Lasso.csv",header = TRUE)
SGL <-  read.csv("C:/Users/vrrani/Documents/SNF/AllSimilarities/Ovarian/SGL.csv",header = TRUE)
SIOL <-  read.csv("C:/Users/vrrani/Documents/SNF/AllSimilarities/Ovarian/SIOL.csv",header = TRUE)
Fused <-  read.csv("C:/Users/vrrani/Documents/SNF/AllSimilarities/Ovarian/Fused.csv",header = TRUE)



head(GFLasso)
head(Lasso)
head(SGL)
head(SIOL)
head(Fused)

GFLassogenes <- GFLasso[,3]
Lassogenes <- Lasso[,3]
SGLgenes <- SGL[,3]
SIOLgenes <- SIOL[,3]
Fusedgenes <- Fused[,3]

data = list(GFLasso = GFLassogenes, Lasso = Lassogenes, SGL = SGLgenes, SIOL = SIOLgenes)

commonPairs <- Reduce(intersect,  data)

venn.plot <- venn.diagram(
x = data,
filename = "CommonMetylationPairs-Ovarian.tiff",
col = "transparent",
fill = c("cornflowerblue", "green", "yellow", "darkorchid1"),
alpha = 0.50,
label.col = c("orange", "white", "darkorchid4", "white",
"white", "white", "white", "white", "darkblue", "white",
"white", "white", "white", "darkgreen", "white"),
cex = 1.5,
fontfamily = "serif",
fontface = "bold",
cat.col = c("darkblue", "darkgreen", "orange", "darkorchid4"),
cat.cex = 1.5,
cat.pos = 0,
cat.dist = 0.07,
cat.fontfamily = "serif",
rotation.degree = 270,
margin = 0.2
);


GFLasso.Fused <- Reduce(intersect,  list(GFLassogenes, Fusedgenes))
Lasso.Fused <- Reduce(intersect,  list(Lassogenes, Fusedgenes))
SGL.Fused <- Reduce(intersect,  list(SGLgenes, Fusedgenes))
SIOL.Fused <- Reduce(intersect,  list(SIOLgenes, Fusedgenes))

length(GFLasso.Fused)
length(Lasso.Fused)
length(SGL.Fused)
length(SIOL.Fused)