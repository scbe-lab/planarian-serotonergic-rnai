elena_motifs <- read.delim2(
  file = "~/colabos/elena_sero/rnaseq_reanalysis/outputs/homer/out.txt",
  dec = ".",
  header = TRUE
)

pitx = "h1SMcG0012776"
lhx = "h1SMcG0003793"

elena_motifs_pitx <- elena_motifs[elena_motifs$PositionID == pitx,]
elena_motifs_lhx <- elena_motifs[elena_motifs$PositionID == lhx,]

# plot(density(elena_motifs$MotifScore), col = "gray", lwd = 2)

hist(
  main = "Histogram of Lhx/Pitx motif likelihood in promoters",
  elena_motifs$MotifScore,
  breaks = 100,
  col = c("#232b2b"), border = "#232b2b",
  xlab = "logOddsPWM"
  )
text(6.44,15000,"*",cex = 2)
text(6.99,8000,"*",cex = 2)
abline(
  v=quantile(elena_motifs$MotifScore),
  col = colorspace::sequential_hcl(5,"GnBu", rev = TRUE),
  lwd = 2
  )