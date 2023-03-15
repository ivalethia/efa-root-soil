# ALGUNAS REFERENCIAS USADAS: https://statsandr.com/blog/correlogram-in-r-how-to-highlight-the-most-correlated-variables-in-a-dataset/ https://cran.r-project.org/web/packages/corrplot/vignettes/corrplot-intro.html https://r-graphics.org/recipe-output-vector-svg https://www.r-bloggers.com/2013/02/exporting-nice-plots-in-r/

library(readr)

T1_CON <- read_csv("T1_CON.csv")
T1_ORG <- read_csv("T1_ORG.csv")
T2_CON <- read_csv("T2_CON.csv")
T2_ORG <- read_csv("T2_ORG.csv")

corrplot2 <- function(data,
                      method = "pearson",
                      sig.level = 0.05,
                      order = "original",
                      diag = FALSE,
                      type = "upper",
                      tl.srt = 90,
                      tl.cex = 1,
                      cl.cex = 1,
                      number.font = 1,
                      number.cex = 0.0,
                      pch.cex = 0.5,
                      insig = "label_sig",
                      mar = c(0, 0, 0, 0)) {
  library(corrplot)
  library(svglite)
  svglite(filename = paste0("Heatmap.svg"), width = 5, height = 5)
  data_incomplete <- data
  data <- data[complete.cases(data), ]
  mat <- cor(data, method = method)
  cor.mtest <- function(mat, method) {
    mat <- as.matrix(mat)
    n <- ncol(mat)
    p.mat <- matrix(NA, n, n)
    diag(p.mat) <- 0
    for (i in 1:(n - 1)) {
      for (j in (i + 1):n) {
        tmp <- cor.test(mat[, i], mat[, j], method = method)
        p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
      }
    }
    colnames(p.mat) <- rownames(p.mat) <- colnames(mat)
    p.mat
  }
  p.mat <- cor.mtest(data, method = method)
  col <- colorRampPalette(c("#474747", "#FFFFFF", "#4472C4"))
  corrplot(mat,
           method = "color", col = col(10), number.font = number.font,
           addgrid.col = "darkgray",
           mar = mar, number.cex = number.cex,
           type = type, order = order, tl.cex = tl.cex, cl.cex = cl.cex, pch.cex = pch.cex,
           addCoef.col = "black", # add correlation coefficient
           tl.col = "black", tl.srt = tl.srt, # rotation of text labels
           p.mat = p.mat, insig = insig, sig.level = c(0.001, 0.01, 0.05), # combine with significance level
           diag = diag # hide correlation coefficients on the diagonal
  )
  dev.off()}

# edit from here
corrplot2(
  data = T1_CON,
  tl.cex = 0.25,
  cl.cex = 0.25, 
  pch.cex = 0.5,
  insig = "label_sig",
  number.cex = 0.001,
  method = "pearson",
  sig.level = 0.05,
  order = "original",
  diag = FALSE,
  type = "upper",
  tl.srt = 90,
)
corrplot2(
  data = T1_ORG,
  tl.cex = 0.25,
  cl.cex = 0.25, 
  pch.cex = 0.5,
  insig = "label_sig",
  number.cex = 0.001,
  method = "pearson",
  sig.level = 0.05,
  order = "original",
  diag = FALSE,
  type = "lower",
  tl.srt = 90,
)
corrplot2(
  data = T2_CON,
  tl.cex = 0.25,
  cl.cex = 0.25, 
  pch.cex = 0.5,
  insig = "label_sig",
  number.cex = 0.001,
  method = "pearson",
  sig.level = 0.05,
  order = "original",
  diag = FALSE,
  type = "upper",
  tl.srt = 90,
)
corrplot2(
  data = T2_ORG,
  tl.cex = 0.25,
  cl.cex = 0.25, 
  pch.cex = 0.5,
  insig = "label_sig",
  number.cex = 0.001,
  method = "pearson",
  sig.level = 0.05,
  order = "original",
  diag = FALSE,
  type = "lower",
  tl.srt = 90,
)
