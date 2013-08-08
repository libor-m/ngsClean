
### Given a tabulated file this script will try to fit a mixture of two distributions (Gamma and NegBionm) to the data (currently only to the first column)


library(optparse)
library(ggplot2)
library(reshape2)
library(BB)

option_list <- list(make_option(c('-i', '--in_file'), action='store', type='character', default=NULL, help='Input distrib file'),
                    make_option(c('-o', '--out_file'), action='store', type='character', default=NULL, help='Output PDF file'),
                    make_option(c('-r', '--rnd_sample'), action='store', type='numeric', default=NULL, help='Number of observations to sample, to speed up computation'),
                    make_option(c('-d', '--debug'), action='store', type='numeric', default=0, help='Run on debug mode. 0: disabled; 1: plots sepparate distributions; 2: ignores input file and simulates dataset')
                    )
opt <- parse_args(OptionParser(option_list = option_list))

############################################################################################

#opt <- list(in_file = "../nivara-LC_chrALL.sdepth", rnd_sample=10000, debug=FALSE)
#opt <- list(in_file = "/home/fgvieira/nivara-LC_chrALL.sdepth", debug=FALSE, rnd_sample=10000)
plot_in_data <- c()
plot_sim_data <- c()
plot_out_data <- c()
options(warn=-1)


#################
### Functions ###
#################

# Maximizing function
f <- function(shape, scale, size, mu, mix_prop) {
  g <- dgamma(data, shape=shape, scale=scale)
  nb <- dnbinom(data, size=size, mu=mu)

  # Parameter check
  if( is.nan(sum(g)) || is.infinite(sum(g)) || is.nan(sum(nb)) || is.infinite(sum(nb)) ) {
    print(c(shape, scale, size, mu, mix_prop, sum(g), sum(nb)))
  }
  return( -sum(log(mix_prop*g + (1-mix_prop)*nb)) )
}

# Adapter function
fcn <- function(x) f(x[1], x[2], x[3], x[4], x[5])

# Simulate mixed distribution
gamma.nb.mix <- function(n, shape=1, scale=2, size=10, mu=30, mix_prop=0.5) {
  u <- runif(n)
  apply( as.matrix(u), 1, function(x) ifelse(x<=mix_prop, rgamma(1,shape=shape,scale=scale), rnbinom(1,size=size,mu=mu)) )
}



############
### Main ###
############

if ( opt$debug != 2) {
  #################
  ### Read file ###
  #################
  cat("# Reading data from file", opt$in_file, fill=TRUE)
  in_data <- read.table(opt$in_file, sep="\t", header=TRUE, row.names = 1)#, nrows=100000)
} else {
  #############
  ### Debug ###
  #############
  cat("# Simulating data...", fill=TRUE)
  # generate test data to check model efficiency
  shape <- 1; scale <- 100; size <- 5; mu <- 500; mix_prop <- 0.5

  # Simulate distributions separately (round and remove zeros from gamma)
  sim_dist <- round(rgamma(100000,shape=shape,scale=scale))
  sim_dist <- sim_dist[sim_dist != 0]
  plot_sim_data <- rbind(plot_sim_data, melt(data.frame("sim_gamma"=sim_dist)))
  sim_dist <- rnbinom(100000,size=size,mu=mu)
  plot_sim_data <- rbind(plot_sim_data, melt(data.frame("sim_nbinom"=sim_dist)))
  
  # Simulate mixed distribution
  sim_data <- round(gamma.nb.mix(1000000, shape, scale, size, mu, mix_prop))
  sim_data <- sim_data[sim_data != 0]
  in_data <- data.frame("sim_mix"=sim_data)
}

# Sample data for analysis speed-up
if( is.null(opt$rnd_sample) ) {
  data <- in_data[,1]
} else {
  seed <- sample(1:1e6, 1)
  cat("# Using", seed, "as seed for random sampling", fill=TRUE)
  set.seed(seed)
  data <- sample(in_data[,1], opt$rnd_sample)
  plot_in_data <- rbind(plot_in_data, melt(data.frame("rnd_sample"=data)))
}

# Plot data (either real or simulated)
plot_in_data <- rbind(plot_in_data, melt(data.frame("in_data"=in_data)))
rm(in_data)



###################
### Fit distrib ###
###################

M <- mean(data)
V <- var(data)

# Starting points (probability parameter is important to remain as close as possible)
start <- list("shape"=(M^2)/V, "scale"=V/M, "size"=ceiling((M^2)/(V-M)), "mu"=M, "mix_prop"=0.5)
cat("# Setting starting points to:", as.numeric(start), fill=TRUE)

# Fit
cat("# Fiting data...", fill=TRUE)
params <- spg( par=as.numeric(start), fn=fcn, lower=c(1e-6,1e-6,1,0,0), upper=c(Inf,Inf,Inf,Inf,1), control=list(maxit=100000, maxfeval = 1e6, ftol=1e-6, triter=100) )
cat(params$message, "after", params$iter, "iters:", params$par, fill=TRUE)


##############
### Output ###
##############

cat("# Preparing plots...", fill=TRUE)
# Get estimated parameters
shape <- params$par[1]
scale <- params$par[2]
size <- params$par[3]
mu <- params$par[4]
prob <- mu/(size+mu)
mix_prop <- params$par[5]

# Plot estimated distributions
cat("# Estimated mixed distribution...", fill=TRUE)
mix <- round(gamma.nb.mix(100000, shape, scale, size, mu, mix_prop))
plot_mixout_data <- melt(data.frame("estimated_mix"=mix))
cat("# Estimated Gamma distribution...", fill=TRUE)
gamma <- round(rgamma(100000, shape=shape, scale=scale))
plot_out_data <- rbind(plot_out_data, melt(data.frame("gamma"=gamma)))
cat("# Estimated NBinom distribution...", fill=TRUE)
nbinom <- rnbinom(100000, size=size, mu=mu)
plot_out_data <- rbind(plot_out_data, melt(data.frame("nbinom"=nbinom)))

cat("# Determining upper and lower cutting points...", fill=TRUE)
l_lim <- qnbinom(0.00001, size=size, mu=mu, lower.tail = TRUE)
u_lim <- qnbinom(0.99999, size=size, mu=mu, lower.tail = TRUE)

dist_overlap <- vector(mode="numeric", length=max(gamma, nbinom))
for ( i in 1:length(dist_overlap) ) {
  dist_overlap[i] <- dgamma(i, shape=shape, scale=scale) - dnbinom(i, size=size, mu=mu)
  if(dist_overlap[i] < 0) {
    break
  }
}
gamma_lim <- ifelse( abs(dist_overlap[i-1]) < abs(dist_overlap[i]) || i==1, i-1, i)
l_lim <- max(gamma_lim, l_lim)
u_lim <- max(gamma_lim, u_lim)

if(l_lim == u_lim) {
  cat("ERROR: Lower and upper limit are the same!\n", fill=TRUE)
  quit()
}

cat("# Outputing results...", fill=TRUE)
# Output
pdf( ifelse(is.null(opt$out_file), paste(opt$in_file, "pdf", sep="."), opt$out_file), width=20 )

mean <- mean(c(plot_sim_data$value, plot_in_data$value, plot_out_data$value, plot_mixout_data$value))
max <- max(plot_sim_data$value, plot_in_data$value, plot_out_data$value, plot_mixout_data$value)
min <- min(plot_sim_data$value, plot_in_data$value, plot_out_data$value, plot_mixout_data$value)

plot <- ggplot() +
  geom_density(data=rbind(plot_in_data, plot_mixout_data), aes(x=value, y=..density.., color=variable, fill=variable), adjust=2, alpha=0.2, size=1) +
  scale_colour_brewer(palette="Spectral") +
  scale_x_continuous(limits=c(0, 2*u_lim)) +
  geom_vline(xintercept=l_lim, colour="blue") +
  geom_vline(xintercept=u_lim, colour="blue") +
#  facet_wrap(~ variable) +
  ggtitle( paste("Fit of Gamma and NBinomial to data (n = ", opt$rnd_sample, " / seed = ", seed, ")\n",
          "shape = ", shape, " / scale = ", scale,
          " / size = ", size, " / mu = ", mu, " / p = ", prob,
          " / mix_prop = ", mix_prop,
          "\nmin = ", min, " / max = ", max, " / mean = ", mean,
          " / l_lim = ", l_lim, " / u_lim = ", u_lim, sep="") )

if( opt$debug > 0 ){
  plot <- plot + geom_histogram(data=rbind(plot_sim_data, plot_out_data), aes(x=value, y=..density.., fill=variable), binwidth=max/1000, position="dodge")
}

print(plot)
dev.off()

cat(l_lim, u_lim, fill=TRUE)

### Memory profilling
#sort( sapply(ls(),function(x){object.size(get(x))}))
