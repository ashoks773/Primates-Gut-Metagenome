Wilcoxon_pairwise <- function(x)
{
  WILCOXON_TEST_TABLE_Gorilla_Human <- data.frame()
  
  CLASS <- colnames(x)
  for (i in seq_along(CLASS[-(length(CLASS))]))
  {
    Gorilla <- x[which(x[,ncol(x)] == "Gorilla"),i]
    Human <- x[which(x[,ncol(x)] == "Human"), i]
    
    Gorilla_mean <- mean(Gorilla)
    Human_mean <- mean(Human)
    Condition <- x[,ncol(x)]
    Values <- c(Gorilla, Human)
    TEST <- pairwise.wilcox.test(x = Values, g = Condition, paired = FALSE, p.adjust.method="fdr", alternative = c("two.sided"), conf.level = 0.95)
    P_Value <- TEST[[3]]
    newline <- data.frame(t(c(paste(CLASS[[i]]), P_Value[1,1], Gorilla_mean, Human_mean)))
    
    WILCOXON_TEST_TABLE_Gorilla_Human <- rbind(WILCOXON_TEST_TABLE_Gorilla_Human, newline)
    
  }
  colnames(WILCOXON_TEST_TABLE_Gorilla_Human)  <- c("Genes", "P-Value", "Gorilla_Mean", "Human_Mean")
  
  write.table(WILCOXON_TEST_TABLE_Gorilla_Human, file = "Wilcoxon.txt", sep = "\t", col.names = TRUE, row.names = FALSE)
  
}