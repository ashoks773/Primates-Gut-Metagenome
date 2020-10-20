Kruskall_Test <- function (x)
{
  TEST_TABLE <- data.frame()
  CLASS <- colnames(x)
  for (i in seq_along(CLASS[-(length(CLASS))]))
  {
    BaAka <- x[which(x[,ncol(x)] == "BaAka"),i]
    Bantu <- x[which(x[,ncol(x)] == "Bantu"), i]
    Dry <- x[which(x[,ncol(x)] == "Dry"),i]
    Wet <- x[which(x[,ncol(x)] == "Wet"), i]
    BaAka_Mean <- mean(BaAka)
    Bantu_Mean <- mean(Bantu)
    Dry_Mean <- mean(Dry)
    Wet_Mean <- mean(Wet)
    Values = x[,i]
    Group = x[,ncol(x)]
    
    Kruskall_Test <- kruskal.test(Values ~ Group, data = x)
    KRUS_P_Value <- Kruskall_Test[[3]]
    Wilcoxon_Test <- pairwise.wilcox.test(x = Values, g = Group, paired = FALSE, p.adjust.method="fdr", alternative = c("two.sided"), conf.level = 0.95, correct = TRUE, exact = FALSE)
    P_VALUE <- Wilcoxon_Test[[3]]
    BaAka_Bantu_P_Value <- P_VALUE[1,1]
    BaAka_Wet_P_Value <- P_VALUE[2,1]
    BaAka_Dry_P_Value <- P_VALUE[3,1]
    newline <- data.frame(t(c(paste(CLASS[[i]]), KRUS_P_Value, BaAka_Bantu_P_Value, BaAka_Wet_P_Value, BaAka_Dry_P_Value, BaAka_Mean, Bantu_Mean, Dry_Mean, Wet_Mean)))
    TEST_TABLE <- rbind(TEST_TABLE, newline)
    
  }
  colnames(TEST_TABLE)  <- c("Species", "Kruskall_Wallis_P_Value", "BaAka_Bantu_P_Value", "BaAka_Wet_P_Value", "BaAka_Dry_P_Value", "BaAka_Mean", "Bantu_Mean", "Dry_Mean", "Wet_Mean")
  write.table(TEST_TABLE, file = "COMPARATIVE_TEST_TABLE", sep = "\t", col.names = TRUE, row.names = FALSE, append = FALSE)
}