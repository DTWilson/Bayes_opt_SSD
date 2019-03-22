library(DiceKriging)
library(ggplot2)
library(plyr)

###########################################################################################
### Plot of a GP constraint, annotating some expected improvemets at a couple of points ###
###########################################################################################

get_power <- function(n)
{
  return(power.t.test(n=n, delta=0.3)$power)
}


df <- data.frame(n = c(10, 110, 260, 290))
df$f<- apply(df, 1, get_power)


model <- km(design = df[,1, drop=FALSE], response = df[,2])

x <- data.frame(n=seq(100, 300))

p <- predict.km(model, newdata=x, type="SK")

x$f <- p$mean
x$sd <- p$sd
x$min <- x$f - 1.96*x$sd
x$max <- x$f + 1.96*x$sd

get_p <- function(input)
{
  return(1-pnorm(0.8, input$f, input$sd))
}

probs <- adply(x, 1, get_p)
probs$EFI <- (260-probs$n)*probs$V1/450 + 0.5
probs$EFI <- ifelse(probs$EFI < 0.5 , 0.5, probs$EFI)

ggplot(x, aes(n, f, ymin=min, ymax=max)) + geom_ribbon(alpha=0.3) + geom_line() +
  geom_line(data=probs, aes(n, V1)) +
  geom_line(data=probs, aes(n, EFI))

# smallest probably feasible point is n=203 with m=0.8605029, s=0.0306694814
m <- 0.8605029
s <- 0.0306694814
# expected improvemet here is (1-pnorm(0.8, m, s))*(260-203) = 55.62

x2 <- data.frame(f=seq(0,m,0.001))
x2$n <- apply(x2, 1, getn, m=m, s=s, n=203)
x2 <- subset(x2, f>0.5)

x3 <- data.frame(f=seq(m,1,0.001))
x3$n <- apply(x3, 1, getn, m=m, s=s, n=203)

# Another example point, n=176  with m=0.8112502 and s=0.0357850962, 
# (1-pnorm(0.8, m, s))*(260-176) = 52.36
m <- 0.8112502
s <- 0.0357850962

x4 <- data.frame(f=seq(0,m,0.001))
x4$n <- apply(x4, 1, getn, m=m, s=s, n=176)
x4 <- subset(x4, f>0.5)

x5 <- data.frame(f=seq(m,1,0.001))
x5$n <- apply(x5, 1, getn, m=m, s=s, n=176)

# Optimal expected feasible improvement is at 190, with
# m = 0.8388252, sd = 3.464424e-02

m <- 0.8388252
s <- 3.464424e-02

x6 <- data.frame(f=c(seq(0,m,0.001), m))
x6$n <- apply(x6, 1, getn, m=m, s=s, n=190)
x6 <- subset(x6, f>0.5)

x7 <- data.frame(f=seq(m,1,0.001))
x7$n <- apply(x7, 1, getn, m=m, s=s, n=190)

ggplot(x6, aes(n, f)) + geom_line(colour="#009E73") + geom_line(data=x7, colour="#009E73") + 
  geom_line(data=x, colour="#D55E00") + geom_ribbon(data=x, aes(ymin=min, ymax=max), alpha=0.1) + 
  #geom_line(data=x4, colour="red") + geom_line(data=x5, colour="red") +
  geom_hline(yintercept = 0.8, linetype=2) + geom_point(data=subset(df, n>100)) +
  scale_y_continuous(limits = c(0.5, 1)) +
  geom_line(data=probs, aes(n, EFI), colour="#CC79A7") +
  theme_minimal() + ylab("Power") + xlab("Sample size")

#ggsave("EGO.png", width=7, height=4.5, dpi=2000)

# For paper
ggplot(x6, aes(n, f)) + geom_line(colour="#009E73", linetype=2) + geom_line(data=x7, colour="#009E73", linetype=2) + 
  geom_line(data=x, colour="#D55E00") + geom_ribbon(data=x, aes(ymin=min, ymax=max), alpha=0.1) + 
  #geom_line(data=x4, colour="red") + geom_line(data=x5, colour="red") +
  geom_point(data=subset(df, n>100)) +
  scale_y_continuous(limits = c(0.5, 1)) +
  geom_line(data=probs, aes(n, EFI), colour="#CC79A7", linetype=3) +
  theme_minimal() + ylab("Power") + xlab("Sample size")

#ggsave("GP_example.pdf", width=5, height=3.3)

####################################################################################
### Chart showing an MC estimte improving as N increases, with 95 and 99% bounds ###
####################################################################################

df <- data.frame(N = seq(1, 22026))
df$x <- rbinom(22026, 1, 0.8)

estimate <- function(input, df)
{
  z <- input$N
  return(data.frame(y = mean(subset(df, N <= z)$x)))
}

df <- adply(df, 1, estimate, df=df)
df$up <- df$y + qnorm(0.99)*sqrt(0.2*0.8/df$N)
df$lo <- df$y - qnorm(0.99)*sqrt(0.2*0.8/df$N)
df$logN <- log(df$N)
df$time <- df$N*112/1000

df2 <- data.frame(logN = c(5.52, 5.52), y= c(0.756, 0.7))
df3 <- data.frame(logN = c(6.91, 6.91), y= c(0.777, 0.9))
df4 <- data.frame(logN = c(8.52, 8.52), y= c(0.788, 0.7))

ggplot(df, aes(logN, y)) + geom_ribbon(aes(ymin=lo, ymax=up), alpha=0.1, fill="#D55E00") + geom_line() +
  scale_y_continuous(limits = c(0.65, 1)) + scale_x_continuous(limits = c(4, 10)) +
  theme_minimal() + 
  geom_line(data=df2, linetype=2) + annotate("text", x = 5.52, y = 0.69, label = "N = 250, 0.5 mins") +
  geom_line(data=df3, linetype=2) + annotate("text", x = 6.91, y = 0.915, label = "N = 1000, 1.9 mins") +
  geom_line(data=df4, linetype=2) + annotate("text", x = 8.52, y = 0.69, label = "N = 5000, 9.3 mins") +
  ylab("Estimated power") + xlab("log(N)")

#ggsave("MC_run.png", width=5, height=3.3, dpi=2000)

# N = 250, logN = 5.52, time = 0.5 mins
# N = 1000, logN = 6.91, time = 1.9
# N = 5000, logN = 8.52, time = 9.3


############################################################
### Plot a 3d scatter of an initial space filling design ###
############################################################

library(plotly)

plot_ly(DoE[1:30,1:3], x=~m, y=~k, z=~j, marker=list(size=3))

