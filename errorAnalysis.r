library(dplyr)

getProportionFromLastLine <- function(fileName){
    if (!file.exists(fileName)){
        prop = NULL
    } else {
        prop = system(paste("tail -n1", fileName), intern=T) %>%
                    strsplit(., "\t") %>% unlist %>% as.numeric
    }
    return(prop)
}


computeEffectiveK <- function (prop){
    sumsq = sum(prop^2)
    if (sumsq == 0){
        return(0)
    } else {    
        return(1/sumsq)
    }
}


computeInferred.k <- function(prop){
    return(sum(prop > .01))
}


chooseSeed <- function(seedAndInferredK){ # seedAndInferredK is a dataframe
    inferred.k = seedAndInferredK$inferred.k
    seed = seedAndInferredK$seed
    idx = which(inferred.k == as.numeric(names(which.max(table(inferred.k)))))[1]
    if (length(idx) == 0){
        return(NULL)
    }
    return(seed[idx])
}


dEploidOutError_3<-function(h.pair, h.pair.true, rel.cost.switch=2,  do.plot=FALSE) {
    l <- ncol(h.pair);
    n.hap <- nrow(h.pair)
    possible.permn = combinat::permn(1:n.hap)
    current.len = length(possible.permn)

    possible.permn[[7]] = c(1,2,1)
    possible.permn[[8]] = c(1,1,2)
    possible.permn[[9]] = c(2,1,1)

    possible.permn[[10]] = c(1,3,1)
    possible.permn[[11]] = c(1,1,3)
    possible.permn[[12]] = c(3,1,1)

    possible.permn[[13]] = c(2,3,2)
    possible.permn[[14]] = c(2,2,3)
    possible.permn[[15]] = c(3,2,2)

    possible.permn[[16]] = c(2,1,2)
    possible.permn[[17]] = c(2,2,1)
    possible.permn[[18]] = c(1,2,2)

    possible.permn[[19]] = c(3,1,3)
    possible.permn[[20]] = c(3,3,1)
    possible.permn[[21]] = c(1,3,3)

    possible.permn[[22]] = c(3,2,3)
    possible.permn[[23]] = c(3,3,2)
    possible.permn[[24]] = c(2,3,3)

    n.permn = length(possible.permn)
    print(n.permn)
    v<-rep(0, n.permn);
    vn<-v;

    tb<-array(0, c(n.permn, l));

    for ( j in 1:n.permn){
        v[j] = sum(h.pair[,1]!=h.pair.true[possible.permn[[j]],1]);
    }

    ee <- rep(0, n.permn)
    same.path <- rep(0, n.permn)

    for (i in 2:l) {
        for ( j in 1:n.permn){
            ee[j] = sum(h.pair[,i]!=h.pair.true[possible.permn[[j]],i]);
            ones = rep(1, n.permn)
            ones[j] = 0
            drop.ones = rep(1, n.permn)
            drop.ones[j] = 0
            if (j<7){
                drop.ones[1:6] = 0
            } else {
                drop.ones[7:24] = 0
            }


            tmp <- v + rel.cost.switch * ones + same.path *drop.ones

            vn[j] <- min(tmp) + ee[j]
            transit.to = which.min(tmp)

            if ( j == transit.to ){
                same.path = same.path + 1
            } else {
                same.path = 0
            }

            tb[j, i] <- transit.to
        }
        v<-vn;
    }

    #decode
    wm<-which.min(v);
    op<-array(0, c(1,l));
    n.s<-0;

    if (wm!=0){
        n.gt<-sum(h.pair[,l]!=h.pair.true[possible.permn[[wm]],l]);
    }


    op[l]<-wm;
    for (i in l:2) {
        wmp<-tb[wm,i];
        if (wmp!=wm) n.s<-n.s+1;

        if (wmp!=0){
            n.gt <- n.gt + sum(h.pair[,i-1] != h.pair.true[possible.permn[[wmp]],i-1]);
        }

        wm<-wmp;
        op[i-1]<-wm;
    }

    if (do.plot) {
        par(mfrow=c(2,1))
        plot(0,0,type="n", xlab="Position", ylab="", yaxt="n", xlim=c(0,ncol(h.pair)), bty="n", ylim=c(-0.5,1.5));
        image(x=1:ncol(h.pair), z=t(rbind(h.pair, rep(0, ncol(h.pair)), h.pair.true)), col=c("white", "black"),
            add=T);
        del<-which(diff(op[1,])!=0);
        if (length(del)>0) points(x=del, y=rep(0.5, length(del)), pch="|", col="red");

        for ( j in 1:n.hap ){
            wd1<-which(h.pair[j,] != h.pair.true[j,]);
    #
            if (length(wd1)>0) {
                points(x=wd1, y=rep(1.2, length(wd1)), pch=25, col=j+1, cex=0.5);
            }

        }
        image(t(op))
    }

    k.eff.permn = rep(0, n.permn)
    for (j in 1:n.permn){
        k.eff.permn[j] = length(unique(possible.permn[[j]]))
    }
    k.eff = k.eff.permn[op]

    drop.strain.permn = rep(0, n.permn)
    for (j in 1:n.permn){
        dropped = which(!c(1,2,3) %in% possible.permn[[j]])
        if ( length(dropped) > 0 ){
            drop.strain.permn[j] = dropped
        }
    }
    drop.strain =  drop.strain.permn[op]

    dropTimes = sum(diff(drop.strain) != 0)
    dropError = sum(drop.strain != 0)
    cat("\nDecoding gives:\nNo. switches:\t", n.s-dropTimes, "\nNo. GT errs:\t", n.gt, "\nNo. Drop errs:\t", dropError,"\n");
    return (list(switchError = n.s, # - dropTimes,
             mutError = n.gt,
             dropError = dropError,
             op = op, drop.strain = drop.strain,
             possible.permn = possible.permn) )
}
