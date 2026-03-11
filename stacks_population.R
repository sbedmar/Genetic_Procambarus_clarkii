library(tidyverse)

# read rows 2 to 38 from a tab-separated values file
sumstats <- read.table(
  "data/stacks_populations/pclarkii.qc_ac_bial.lowmiss_maf.p.sumstats_summary.tsv",
  header = FALSE, sep = "\t", stringsAsFactors = FALSE, skip = 2, nrows = 36
  )
# assign column names
# Pop_ID	Private	Num_Indv	Num_Indv_Var	Num_Indv_StdErr	P	P_Var	P_StdErr	Obs_Het	Obs_Het_Var	Obs_Het_StdErr	Obs_Hom	Obs_Hom_Var	Obs_Hom_StdErr	Exp_Het	Exp_Het_Var	Exp_Het_StdErr	Exp_Hom	Exp_Hom_Var	Exp_Hom_StdErr	Pi	Pi_Var	Pi_StdErr	Fis	Fis_Var	Fis_StdErr
colnames(sumstats) <- c(
  "Pop_ID", "Private", "Num_Indv", 
  "Num_Indv_Var", "Num_Indv_StdErr", 
  "P", "P_Var", "P_StdErr", 
  "Obs_Het", "Obs_Het_Var", "StdErr_Obs_Het", 
  "Obs_Hom", "Obs_Hom_Var", "Obs_Hom_StdErr", 
  "Exp_Het", "Exp_Het_Var", "Exp_Het_StdErr", 
  "Exp_Hom", "Exp_Hom_Var", "Exp_Hom_StdErr", 
  "Pi", "Pi_Var", "Pi_StdErr", 
  "Fis", "Fis_Var", "Fis_StdErr"
)

# ggplot of Obs_Het by Pop_ID sorted by Obs_Het
sumstats <- sumstats |> arrange(desc(Obs_Het))

ggplot(sumstats, aes(x = factor(Pop_ID, levels = Pop_ID), y = Obs_Het)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  geom_errorbar(aes(ymin = Obs_Het - StdErr_Obs_Het, ymax = Obs_Het + StdErr_Obs_Het), width = 0.2) +
  theme_minimal() +
  labs(title = "Observed Heterozygosity by Population",
       x = "Population ID",
       y = "Observed Heterozygosity") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# sumstats 2
sumstats2 <- read.table(
  "data/stacks_populations/pclarkii.qc_ac_bial.lowmiss_maf.p.sumstats_summary.tsv",
  header = FALSE, sep = "\t", stringsAsFactors = FALSE, skip = 38, nrows = 36
)
# assign column names
# Pop_ID	Private_Sites	Variant_Sites	Polymorphic_Sites	%Polymorphic_Loci	Num_Indv	Num_Indv_Var	Num_Indv_StdErr	P	P_Var	P_StdErr	Obs_Het	Obs_Het_Var	Obs_Het_StdErr	Obs_Hom	Obs_Hom_Var	Obs_Hom_StdErr	Exp_Het	Exp_Het_Var	Exp_Het_StdErr	Exp_Hom	Exp_Hom_Var	Exp_Hom_StdErr	Pi	Pi_Var	Pi_StdErr	Fis	Fis_Var	Fis_StdErr
colnames(sumstats2) <- c(
  "Pop_ID", "Private", "Sites", "Variant_Sites", 
  "Polymorphic_Sites", "Percent_Polymorphic_Loci", 
  "Num_Indv", "Num_Indv_Var", "Num_Indv_StdErr", 
  "P", "P_Var", "P_StdErr", 
  "Obs_Het", "Obs_Het_Var", "Obs_Het_StdErr", 
  "Obs_Hom", "Obs_Hom_Var", "Obs_Hom_StdErr", 
  "Exp_Het", "Exp_Het_Var", "Exp_Het_StdErr", 
  "Exp_Hom", "Exp_Hom_Var", "Exp_Hom_StdErr", 
  "Pi", "Pi_Var", "Pi_StdErr", 
  "Fis", "Fis_Var", "Fis_StdErr"
)

# ggplot of Pi by Pop_ID
sumstats2 <- sumstats2 |> arrange(desc(Pi))

ggplot(sumstats, aes(x = factor(Pop_ID, levels = Pop_ID), y = Pi)) +
  geom_bar(stat = "identity", fill = "darkgreen") +
  geom_errorbar(aes(ymin = Pi - Pi_StdErr, ymax = Pi + Pi_StdErr), width = 0.2) +
  theme_minimal() +
  labs(title = "Nucleotide Diversity (Pi) by Population",
       x = "Population ID",
       y = "Nucleotide Diversity (Pi)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggplot(sumstats2, aes(x = Pop_ID, y = Pi)) +
  geom_bar(stat = "identity", fill = "darkgreen") +
  geom_errorbar(aes(ymin = Pi - Pi_StdErr, ymax = Pi + Pi_StdErr), width = 0.2) +
  theme_minimal() +
  labs(title = "Nucleotide Diversity (Pi) by Population",
       x = "Population ID",
       y = "Nucleotide Diversity (Pi)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# ggplot of Obs_Het by Pop_ID
sumstats2 <- sumstats2 |> arrange(desc(Obs_Het))
ggplot(sumstats2, aes(x =  factor(Pop_ID, levels = Pop_ID), y = Obs_Het)) +
  geom_bar(stat = "identity", fill = "orange") +
  geom_errorbar(aes(ymin = Obs_Het - Obs_Het_StdErr, ymax = Obs_Het + Obs_Het_StdErr), width = 0.2) +
  theme_minimal() +
  labs(title = "Observed Heterozygosity by Population",
       x = "Population ID",
       y = "Observed Heterozygosity") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

######

# sumstats 2
sumstats2 <- read.table(
  "data/stacks_populations/pclarkii.qc_ac_bial.lowmiss_maf.p.sumstats_summary.tsv",
  header = FALSE, sep = "\t", stringsAsFactors = FALSE, skip = 38, nrows = 36
)

# assign column names
# Pop_ID	Private_Sites	Variant_Sites	Polymorphic_Sites	%Polymorphic_Loci	Num_Indv	Num_Indv_Var	Num_Indv_StdErr	P	P_Var	P_StdErr	Obs_Het	Obs_Het_Var	Obs_Het_StdErr	Obs_Hom	Obs_Hom_Var	Obs_Hom_StdErr	Exp_Het	Exp_Het_Var	Exp_Het_StdErr	Exp_Hom	Exp_Hom_Var	Exp_Hom_StdErr	Pi	Pi_Var	Pi_StdErr	Fis	Fis_Var	Fis_StdErr
colnames(sumstats2) <- c(
  "Pop_ID", "Private", "Sites", "Variant_Sites", 
  "Polymorphic_Sites", "Percent_Polymorphic_Loci", 
  "Num_Indv", "Num_Indv_Var", "Num_Indv_StdErr", 
  "P", "P_Var", "P_StdErr", 
  "Obs_Het", "Obs_Het_Var", "Obs_Het_StdErr", 
  "Obs_Hom", "Obs_Hom_Var", "Obs_Hom_StdErr", 
  "Exp_Het", "Exp_Het_Var", "Exp_Het_StdErr", 
  "Exp_Hom", "Exp_Hom_Var", "Exp_Hom_StdErr", 
  "Pi", "Pi_Var", "Pi_StdErr", 
  "Fis", "Fis_Var", "Fis_StdErr"
)

# 
pop_ordered_list <- c("IB", "JB", "HB", "CS", "AM", "DM", "AB", "EA", "AA", "BB", "CQ", "AR", "EC", "FB", "BA", "BM", "CB", "EB", "EM", "DB", "IM", "JM", "IA", "JA", "GB", "GM", "FM", "GA", "CA", "CM", "HA", "HM", "DA", "FA", "EX", "LA")
sumstats2$Pop_ID <- factor(sumstats2$Pop_ID, levels = pop_ordered_list)

# plot Obs_Het and Pi side by side in horizontal barplot
library(ggplot2)
library(tidyr)
library(cowplot)
sumstats2_long <- sumstats2 |> select(Pop_ID, Obs_Het, Pi) |>
  pivot_longer(cols = c(Obs_Het, Pi), names_to = "Metric", values_to = "Value")
ggplot(sumstats2_long, aes(x = Value, y = Pop_ID, fill = Metric)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_minimal() +
  labs(title = "Observed Heterozygosity and Nucleotide Diversity by Population",
       x = "Value",
       y = "Population ID") +
  scale_fill_manual(values = c("Obs_Het" = "orange", "Pi" = "darkgreen")) +
  theme(legend.title = element_blank())
# plot Pi only with scale starting at 0.15
ggplot(sumstats2, aes(x = Pi, y = Pop_ID)) +
  geom_bar(stat = "identity", fill = "darkgreen") +
  theme_minimal() +
  labs(title = "Nucleotide Diversity (Pi) by Population",
       x = "Nucleotide Diversity (Pi)",
       y = "Population ID") +
  xlim(0.15, 0.25)

bars <- ggplot(sumstats2, aes(y = Pop_ID)) +
  geom_rect(aes(xmin = 0.15,
                xmax = Pi,
                ymin = as.numeric(Pop_ID) - 0.45,
                ymax = as.numeric(Pop_ID) + 0.45),
            fill = "darkgreen") +
  theme_minimal() +
  labs(x = "Nucleotide Diversity (Pi)",
       y = "Population ID") + theme_void()
bars
k=6
print(k)
bars2 <- bars
for (pop in rev(pop_ordered_list)){
  filename <- paste0("plots/admixture/", k, "/", pop, ".png")
  ggdraw(bars) +
    draw_image(filename,
               x = 0.18,        # center horizontally
               y = 0.5,        # center vertically
               hjust = 0.5,
               vjust = 0.5,
               scale = 0.05)
  bars + draw_image(filename,
                    x = 0.25, y = 10)
    }

####
