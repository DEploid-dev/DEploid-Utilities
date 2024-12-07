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
    return (list(switchError = n.s - dropTimes,
             mutError = n.gt,
             dropError = dropError,
             op = op, drop.strain = drop.strain,
             possible.permn = possible.permn) )
}


dEploidOutError_2 <-function(h.pair, h.pair.true, rel.cost.switch=2, do.plot=FALSE) {
    l <- ncol(h.pair);
    n.hap <- nrow(h.pair)
    possible.permn = combinat::permn(1:n.hap)

    possible.permn[[3]] = c(1,1)
    possible.permn[[4]] = c(2,2)

    #possible.permn = list()
    #possible.permn[[1]] = c(1,1)
    #possible.permn[[2]] = c(2,2)
    #possible.permn[[3]] = c(1,2)
    #possible.permn[[4]] = c(2,1)

    n.permn = length(possible.permn)
    v<-rep(0, n.permn);
    vn<-v;

    tb<-array(0, c(n.permn, l));

    for ( j in 1:n.permn){
        v[j] = sum(h.pair[,1]!=h.pair.true[possible.permn[[j]],1]);
    }

#rel.cost.drop.current = 0
    ee <- rep(0, n.permn)
    same.path <- rep(0, n.permn)

    for (i in 2:l) {
#        rel.cost.drop.current = rel.cost.drop.current + 1
        for ( j in 1:n.permn){
            ee[j] = sum(h.pair[,i]!=h.pair.true[possible.permn[[j]],i]);

            ones = rep(1, n.permn)
            ones[j] = 0
            drop.ones = rep(1, n.permn)
            drop.ones[j] = 0
            if (j<3){
                drop.ones[1:2] = 0
            } else {
                drop.ones[3:4] = 0
            }
#            tmp <- v + rel.cost.switch * ones
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


#        if ( rel.cost.drop.current > rel.cost.drop){
#            rel.cost.drop.current = 0
#        }
        v<-vn;
#        cat ("site ", i, " cost: ", v,"\n")
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

    k.eff.permn = rep(0, n.permn)
    for (j in 1:n.permn){
        k.eff.permn[j] = length(unique(possible.permn[[j]]))
    }
    k.eff = k.eff.permn[op]

    changed.permn.at = which(diff(op[1,]) != 0)
    changed.k.eff.at = which(diff(k.eff) != 0)

#print("changed.permn.at")
#print(changed.permn.at)
#print("changed.k.eff.at")
#print(changed.k.eff.at)


    drop.strain.permn = rep(0, n.permn)
    for (j in 1:n.permn){
        dropped = which(!c(1,2) %in% possible.permn[[j]])
        if ( length(dropped) > 0 ){
            drop.strain.permn[j] = dropped
        }
    }
    drop.strain =  drop.strain.permn[op]

#    definite.switch = sum(!changed.permn.at %in% changed.k.eff.at)
    dropTimes = sum(diff(drop.strain) != 0) #length(changed.k.eff.at)
    dropError = sum(drop.strain != 0)
    cat("\nDecoding gives:\nNo. switches:\t", n.s - dropTimes, "\nNo. GT errs:\t", n.gt, "\nNo. Drop switches", dropTimes, "\nNo. Drop errs:\t", dropError,"\n");


    hap.pair.true.idx = array(unlist(possible.permn[op]), c(2, length(op)))
    genotype.error = rep(0, n.hap)
    switch.error = rep(0, n.hap)
    drop.error = rep(0, n.hap)
    for ( j in 1:n.hap ){
        wd1 = c()
        for ( i in 1:l ){
            if ( h.pair[j,i] != h.pair.true[hap.pair.true.idx[j,i],i] ){
                wd1 = c(wd1, i)
            }
        }
        genotype.error[j] = length(wd1)
    }

    switch.error[1] = (sum(diff(hap.pair.true.idx[1,])!=0))
    switch.error[2] = (sum(diff(hap.pair.true.idx[2,])!=0))
    drop.error[1] = (sum(drop.strain == 1))
    drop.error[2] = (sum(drop.strain == 2))

	if (do.plot) {
#        par(mfrow=c(2,1))
        layout(matrix(c(1,1,1,1,2,2), 3, 2, byrow = T))
		plot(0,0,type="n", xlab="Position", ylab="", yaxt="n", xlim=c(0,ncol(h.pair)), bty="n", ylim=c(-0.5,1.5));
		image(x=1:ncol(h.pair), z=t(rbind(h.pair, rep(0, ncol(h.pair)), h.pair.true)), col=c("white", "black"),
			add=T);
		del<-which(diff(op[1,])!=0);
		if (length(del)>0) {points(x=del, y=rep(0.5, length(del)), pch="|", col="red");}


        for ( j in 1:n.hap ){
            wd1 = c()
            for ( i in 1:l ){
                if ( h.pair[j,i] != h.pair.true[hap.pair.true.idx[j,i],i] ){
                    wd1 = c(wd1, i)
                }
            }

            if (length(wd1)>0) {
                points(x=wd1, y=rep(1.2, length(wd1)), pch=25, col=j+1, cex=0.5);
#                print(length(wd1))
            }

        }
		plot(0,0,type="n", xlab="Position", ylab="", yaxt="n", xlim=c(0,ncol(h.pair)), bty="n", ylim=c(-0.5,1.5));
#        image(x=1:ncol(h.pair), z = t(rbind(op, array(drop.strain,c(1,l)))), col = heat.colors(4), add=T)
        image(x=1:ncol(h.pair), z = t(op), col = heat.colors(4), add=T)
    }


        return (list(switchError = n.s - dropTimes,
                 mutError = n.gt,
                 dropError = dropError,
                 op = op, drop.strain = drop.strain,
                 possible.permn = possible.permn, tb = tb,
                 switch.error = switch.error,
                 drop.error = drop.error,
                 genotype.error = genotype.error )
                 )
}
