pfizer.old <- pfizer.df[1,][-1]#['tot_7d',] # rimuovo il nome con [-1]
pfizer.data.pv.old <- c(rep(1, len = pfizer.old[[1]]), 
                        rep(0, len = pfizer.old[[3]] - pfizer.old[[1]]))

pfizer.data.pp.old <- c(rep(1, len = pfizer.old[[2]]),
                        rep(0, len = pfizer.old[[4]] - pfizer.old[[2]]))

pfizer.new <- pfizer.df[2,][-1] #['tot_6m',]
pfizer.data.pv.new <- c(rep(1, len = pfizer.new[[1]]), 
                        rep(0, len = pfizer.new[[3]] - pfizer.new[[1]]))

pfizer.data.pp.new <- c(rep(1, len = pfizer.new[[2]]),
                        rep(0, len = pfizer.new[[4]] - pfizer.new[[2]]))

moderna <- moderna.df[1,][-1] #['tot',]
moderna.data.pv <- c(rep(1, len = moderna[[1]]), 
                        rep(0, len = moderna[[3]] - moderna[[1]]))

moderna.data.pp <- c(rep(1, len = moderna[[2]]),
                        rep(0, len = moderna[[4]] - moderna[[2]]))

astrazeneca <- astrazeneca.df[1,][-1] #['tot',]
astrazeneca.data.pv <- c(rep(1, len = astrazeneca[[1]]), 
                        rep(0, len = astrazeneca[[3]] - astrazeneca[[1]]))

astrazeneca.data.pp <- c(rep(1, len = astrazeneca[[2]]),
                        rep(0, len = astrazeneca[[4]] - astrazeneca[[2]]))

janssen <- janssen.df[1,][-1] #['tot',]
janssen.data.pv <- c(rep(1, len = janssen[[1]]), 
                        rep(0, len = janssen[[3]] - janssen[[1]]))

janssen.data.pp <- c(rep(1, len = janssen[[2]]),
                        rep(0, len = janssen[[4]] - janssen[[2]]))