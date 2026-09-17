d <- read.table("parent_depths.tsv", header = TRUE, sep = "\t")

# histograms, capped at 60x so the long tail doesn't flatten everything
par(mfrow = c(1, 2))
hist(pmin(d$p1_11AD, 60), breaks = 60, main = "11-AD parental depth",
     xlab = "AD-sum depth (capped at 60)", col = "steelblue")
abline(v = 10, col = "red", lwd = 2)   # your confidence threshold

hist(pmin(d$p2_7AD, 60), breaks = 60, main = "7-AD parental depth",
     xlab = "AD-sum depth (capped at 60)", col = "darkorange")
abline(v = 10, col = "red", lwd = 2)

# quick summary
summary(d)
mean(d$p1_11AD >= 10 & d$p2_7AD >= 10)   # fraction of sites confident in BOTH


thr <- 0:40
n_both <- sapply(thr, function(t) sum(d$p1_11AD >= t & d$p2_7AD >= t))
plot(thr, n_both, type = "l", lwd = 2,
     xlab = "min depth in BOTH parents", ylab = "# confident sites",
     main = "Usable sites vs parental depth threshold")
abline(v = 10, col = "red", lty = 2)