library(ggplot2)
library(gridExtra)
stouffer = function(b1, b2, epsilon, n) {
    zvals = replicate(n = n, {
        obs = rbeta(n = 1, shape1 = b1 + epsilon, shape2 = b2)
        pval = pbeta(q = obs, shape1 = b1, shape2 = b2, lower.tail = FALSE)
        zval = qnorm(p = 1 - pval)
    })
    stouffer = sum(zvals) / sqrt(n)
    pnorm(stouffer, lower.tail = FALSE)
    return(list(pval = pnorm(stouffer, lower.tail = FALSE),
                stouffer_z = stouffer,
                mean_null = b1 / (b1 + b2),
                mean_true = (b1 + epsilon) / (b1 + b2 + epsilon)))
}

b1 = 2; b2 = 2; epsilon = .15
zs = replicate(10000, stouffer(b1,b2,epsilon,2000)$stouffer_z)
ps = replicate(10000, stouffer(b1,b2,epsilon,2000)$pval)
theme_set(theme_bw())
z_plot = qplot(zs, geom = "histogram") + xlab("z-value") +
    ggtitle("Distribution of Stouffer's Z\nwith n = 2000")
p_plot = qplot(ps, geom = "histogram") +
    scale_x_log10("p-value") + ggtitle("Distribution of p-values from\nStouffer's Z with n = 2000")

pdf("stouffer_distributions.pdf", width=6,height=2.5)
grid.arrange(z_plot, p_plot, ncol = 2)
dev.off()
xgrid = seq(0,1,length.out=1000)
beta_df = data.frame(density = c(dbeta(xgrid, b1, b2), dbeta(xgrid, b1 + epsilon, b2)),
                     x = c(xgrid, xgrid),
                     distribution = rep(c("null", "true"), each = length(xgrid)))
pdf("two_betas.pdf",width=4,height=2.25)
ggplot(beta_df) +
    geom_line(aes(x = x, y = density, linetype = distribution)) +
    xlab("") +
    ggtitle("Distribution of betas\nused for simulation")
dev.off()
