library('rjags')
library('coda')

mchain <- function(data, prior = 'dbeta(1,1)', niter = 15000 , nburnin = 500) {
    # Load the data:
    y = data            # The y values are in the column named y.
    Ntotal = length(y)  # Compute the total number of "flips".
    dataList = list(    # Put the information into a list.
    y = y ,
    Ntotal = Ntotal 
    )

    # Define the model:
    modelString = paste(
    "model {
    for ( i in 1:Ntotal ) {
        y[i] ~ dbern( theta )
    }
    theta ~ ",prior,"}")
    # close quote for modelString
    writeLines( modelString , con = "TEMPmodel.txt" )

    # Initialize the chains based on MLE of data.
    # Option: Use single initial value for all chains:
    # thetaInit = sum(y)/length(y)
    # initsList = list( theta=thetaInit )
    # Option: Use function that generates random values for each chain:
    initsList <- function() {
        resampledY = sample( y , replace=TRUE )
        thetaInit = sum(resampledY)/length(resampledY)
        #   thetaInit = 0.001+0.998*thetaInit # keep away from 0,1
        return(list(theta = thetaInit))
    }

    # Run the chains:
    jagsModel <- jags.model(file="TEMPmodel.txt" , data = dataList , inits = initsList , 
                            n.chains = 4 , n.adapt = 4000, quiet = TRUE)
    update(jagsModel, n.iter = nburnin)
    codaSamples <- coda.samples(jagsModel, variable.names = c("theta"),
                               n.iter=niter)
    return(codaSamples)
}

# Pfizer old (7 days after 2nd dose)
rownames <- c('tot_old', 'tot_new', '12-15', '16-64', '65-74', 'more75')
vpos <- c(8, 77, 0, 7, 1, 0)
ppos <- c(162, 850, 18, 143, 14, 5)
vtot <- c(18198, 18198, 1131, 14216, 3176, 804)
ptot <- c(18325, 18325, 1129, 14299, 3226, 812)
pfizer.df <- data.frame(row.names = rownames, vpos = vpos, ppos = ppos, vtot = vtot, ptot = ptot)

pfizer.old <- pfizer.df['tot_old',]
pfizer.data.pv.old <- c(rep(1, len = pfizer.old[[1]]), 
                        rep(0, len = pfizer.old[[3]] - pfizer.old[[1]]))

pfizer.data.pp.old <- c(rep(1, len = pfizer.old[[2]]),
                        rep(0, len = pfizer.old[[4]] - pfizer.old[[2]]))

pfizer.new <- pfizer.df['tot_new',]
pfizer.data.pv.new <- c(rep(1, len = pfizer.new[[1]]), 
                        rep(0, len = pfizer.new[[3]] - pfizer.new[[1]]))

pfizer.data.pp.new <- c(rep(1, len = pfizer.new[[2]]),
                        rep(0, len = pfizer.new[[4]] - pfizer.new[[2]]))

print('inferring pv for Pfizer old')
pv.pfizer.old <- mchain(pfizer.data.pv.old, prior = 'dbeta(2.88,93.12)')
print('inferring pp for Pfizer old')
pp.pfizer.old <- mchain(pfizer.data.pp.old, prior = 'dbeta(1.41,66.11)')

eff.pfizer.old <- 1 - as.matrix(pv.pfizer.old)/as.matrix(pp.pfizer.old)

# Pfizer new (6 months after 2nd dose)
print('inferring pv for Pfizer new')
pv.pfizer.new <- mchain(pfizer.data.pv.new, prior = 'dbeta(2.88,93.12)')
print('inferring pp for Pfizer new')
pp.pfizer.new <- mchain(pfizer.data.pp.new, prior = 'dbeta(1.41,66.11)')

eff.pfizer.new <- 1 - as.matrix(pv.pfizer.new)/as.matrix(pp.pfizer.new)

# Moderna
young.pop= 3732
vratio = 2/3
moderna.df <- data.frame(
              row.names=c('tot', '12-17', '18-64', 'more65'),
              vpos=c(11, 0, 7, 4),
              ppos=c(185, 4, 156, 29),
              vtot=c(14134, young.pop*vratio, 10551, 3583),
              ptot=c(14073, young.pop*(1-vratio), 10521, 3552))

moderna <- moderna.df['tot',]
moderna.data.pv <- c(rep(1, len = moderna[[1]]), 
                        rep(0, len = moderna[[3]] - moderna[[1]]))

moderna.data.pp <- c(rep(1, len = moderna[[2]]),
                        rep(0, len = moderna[[4]] - moderna[[2]]))

print('inferring pv for Moderna')
pv.moderna <- mchain(moderna.data.pv, prior = 'dbeta(2.88,93.12)')
print('inferring pp for Moderna')
pp.moderna <- mchain(moderna.data.pp, prior = 'dbeta(1.41,66.11)')

eff.moderna <- 1 - as.matrix(pv.moderna)/as.matrix(pp.moderna)

# Astrazeneca
rownames = c('tot')
vpos = c(84)
ppos = c(248)
vtot = c(8597)
ptot = c(8581)
astrazeneca.df =  data.frame(row.names = rownames, vpos=vpos, ppos=ppos, vtot=vtot, ptot=ptot)

astrazeneca <- astrazeneca.df['tot',]
astrazeneca.data.pv <- c(rep(1, len = astrazeneca[[1]]), 
                        rep(0, len = astrazeneca[[3]] - astrazeneca[[1]]))

astrazeneca.data.pp <- c(rep(1, len = astrazeneca[[2]]),
                        rep(0, len = astrazeneca[[4]] - astrazeneca[[2]]))

print('inferring pv for Astrazeneca')
pv.astrazeneca <- mchain(astrazeneca.data.pv, prior = 'dbeta(2.88,93.12)')
print('inferring pp for Astrazeneca')
pp.astrazeneca <- mchain(astrazeneca.data.pp, prior = 'dbeta(1.41,66.11)')

eff.astrazeneca <- 1 - as.matrix(pv.astrazeneca)/as.matrix(pp.astrazeneca)

# Janssen
janssen.df <- data.frame(
              row.names=c('tot', '18-59', 'more60'),
              vpos=c(116, 95, 21),
              ppos=c(348, 260, 88),
              vtot=c(19630, 12830, 6800),
              ptot=c(19691, 12881, 6810))

janssen <- janssen.df['tot',]
janssen.data.pv <- c(rep(1, len = janssen[[1]]), 
                        rep(0, len = janssen[[3]] - janssen[[1]]))

janssen.data.pp <- c(rep(1, len = janssen[[2]]),
                        rep(0, len = janssen[[4]] - janssen[[2]]))

print('inferring pv for Janssen')
pv.janssen <- mchain(janssen.data.pv, prior = 'dbeta(2.88,93.12)')
print('inferring pp for Janssen')
pp.janssen <- mchain(janssen.data.pp, prior = 'dbeta(1.41,66.11)')

eff.janssen <- 1 - as.matrix(pv.janssen)/as.matrix(pp.janssen)

save(
    pv.pfizer.old, pp.pfizer.old, eff.pfizer.old, pv.pfizer.new, pp.pfizer.new, eff.pfizer.new, 
    pv.moderna, pp.moderna, eff.moderna, pv.astrazeneca, pp.astrazeneca, eff.astrazeneca, 
    pv.janssen, pp.janssen, eff.janssen, file = './mcmc_data/mcmc_data.RData'
)