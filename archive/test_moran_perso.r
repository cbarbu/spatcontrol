## get clean data
data.base <- read.csv("DB_mm_blocks_05May2011.csv",header=T);
# data<-data.base[data.base$northing>8185000&data.base$northing<8185500&data.base$easting>233250&data.base$easting<233750,]
data<-data.base[data.base$northing&data.base$northing&data.base$easting&data.base$easting,]
library("spdep");
library("spam")
threshold<-100;
dnb <- dnearneigh(as.matrix(data[,c("easting","northing")]), 0, threshold)
isolated<-which(card(dnb)<1)

if(length(isolated)>0){
	data<-data[-isolated,];
}

## functions
moran_perso<-function(dmt,values){
	# dmt a spam matrix with the weights between houses
	# values the values for each location
	# easily 20 and probably up to 100 times faster than moran.test

	# correction of dmt for correct calculus of I
	dim_mp<-length(data_non_NA[,1])
	one_vect<-rep(1,dim_mp)
	#need 0 on the diagonal of weight
	diag(dmt)<-0
	#need weight normalized by row
	rsum<-dmt%*%one_vect
	rsum[rsum==0]<-1;
	rsum_mat<-dmt*0;
	diag(rsum_mat)<-1/rsum
	dmt<-rsum_mat%*%dmt

	# matrix of residuals
	res_z<- values-mean(values)
	res_z_mat<-diag.spam(1,length(values));
	diag.spam(res_z_mat)<-res_z

	### Moran's I
	# covariance term
	int<-dmt%*%res_z_mat;
	int<-res_z_mat%*%int

	MI<-length(res_z_mat)*one_vect%*%int%*%one_vect/(sum(dmt@entries)*res_z%*%res_z)
	return(drop(MI))
}
moran.perso <- function (x, listw, n, S0, zero.policy = NULL, NAOK = FALSE) 
{
    if (is.null(zero.policy)) 
        zero.policy <- get("zeroPolicy", env = .spdepOptions)
    stopifnot(is.logical(zero.policy))
    n1 <- length(listw$neighbours)
    x <- c(x)
    if (n1 != length(x)) 
        stop("objects of different length")
    xx <- mean(x, na.rm = NAOK)
    z <- x - xx
    zz <- sum(z^2, na.rm = NAOK)
    K <- (length(x) * sum(z^4, na.rm = NAOK))/(zz^2)
    lz <- lag.listw(listw, z, zero.policy = zero.policy, NAOK = NAOK)
    I <- (n/S0) * ((sum(z * lz, na.rm = NAOK))/zz)
    print(I)
    res <- list(I = I, K = K)
    res
}


## moran official
moran.test.perso <- function (x, listw, randomisation = TRUE, zero.policy = NULL, 
    alternative = "greater", rank = FALSE, na.action = na.fail, 
    spChk = NULL, adjust.n = TRUE) 
{
    alternative <- match.arg(alternative, c("greater", "less", 
        "two.sided"))
    if (!inherits(listw, "listw")) 
        stop(paste(deparse(substitute(listw)), "is not a listw object"))
    if (!is.numeric(x)) 
        stop(paste(deparse(substitute(x)), "is not a numeric vector"))
    if (is.null(zero.policy)) 
        zero.policy <- get("zeroPolicy", env = .spdepOptions)
    stopifnot(is.logical(zero.policy))
    if (is.null(spChk)) 
        spChk <- get.spChkOption()
    if (spChk && !chkIDs(x, listw)) 
        stop("Check of data and weights ID integrity failed")
    xname <- deparse(substitute(x))
    wname <- deparse(substitute(listw))
    NAOK <- deparse(substitute(na.action)) == "na.pass"
    x <- na.action(x)
    na.act <- attr(x, "na.action")
    if (!is.null(na.act)) {
        subset <- !(1:length(listw$neighbours) %in% na.act)
        listw <- subset(listw, subset, zero.policy = zero.policy)
    }
    n <- length(listw$neighbours)
    if (n != length(x)) 
        stop("objects of different length")
    wc <- spweights.constants(listw, zero.policy = zero.policy, 
        adjust.n = adjust.n)
    S02 <- wc$S0 * wc$S0
    res <- moran.perso(x, listw, wc$n, wc$S0, zero.policy = zero.policy, 
        NAOK = NAOK)
    I <- res$I
    K <- res$K
    if (rank) 
        K <- (3 * (3 * wc$n^2 - 7))/(5 * (wc$n^2 - 1))
    EI <- (-1)/wc$n1
    if (randomisation) {
        VI <- wc$n * (wc$S1 * (wc$nn - 3 * wc$n + 3) - wc$n * 
            wc$S2 + 3 * S02)
        tmp <- K * (wc$S1 * (wc$nn - wc$n) - 2 * wc$n * wc$S2 + 
            6 * S02)
        VI <- (VI - tmp)/(wc$n1 * wc$n2 * wc$n3 * S02)
        VI <- VI - EI^2
    }
    else {
        VI <- (wc$nn * wc$S1 - wc$n * wc$S2 + 3 * S02)/(S02 * 
            (wc$nn - 1))
    VI <- VI - EI^2
    }
    ZI <- (I - EI)/sqrt(VI)
    statistic <- ZI
    names(statistic) <- "Moran I statistic standard deviate"
    if (alternative == "two.sided") 
	    PrI <- 2 * pnorm(abs(ZI), lower.tail = FALSE)
    else if (alternative == "greater") 
    PrI <- pnorm(ZI, lower.tail = FALSE)
    else PrI <- pnorm(ZI)
    if (!is.finite(PrI) || PrI < 0 || PrI > 1) 
	    warning("Out-of-range p-value: reconsider test arguments")
    vec <- c(I, EI, VI)
    names(vec) <- c("Moran I statistic", "Expectation", "Variance")
    method <- paste("Moran's I test under", ifelse(randomisation, 
		    "randomisation", "normality"))
    data.name <- paste(xname, ifelse(rank, "using rank correction", 
		    ""), "\nweights:", wname, ifelse(is.null(na.act), "", 
		    paste("\nomitted:", paste(na.act, collapse = ", "))), 
	    "\n")
    res <- list(statistic = statistic, p.value = PrI, estimate = vec, 
	    alternative = alternative, method = method, data.name = data.name)
    if (!is.null(na.act)) 
	    attr(res, "na.action") <- na.act
    class(res) <- "htest"
    res
}

## carefull this one doesn't normalize weight as it "should"
Moran.I.perso <- function (x, weight, scaled = FALSE, na.rm = FALSE, alternative = "two.sided") 
{
    if (dim(weight)[1] != dim(weight)[2]) 
        stop("'weight' must be a square matrix")
    n <- length(x)
    if (dim(weight)[1] != n) 
        stop("'weight' must have as many rows as observations in 'x'")
    ei <- -1/(n - 1)
    nas <- is.na(x)
    if (any(nas)) {
        if (na.rm) {
            x <- x[!nas]
            n <- length(x)
            weight <- weight[!nas, !nas]
        }
        else {
            warning("'x' has missing values: maybe you wanted to set na.rm=TRUE?")
            return(list(observed = NA, expected = ei, sd = NA, 
                p.value = NA))
        }
    }
    ROWSUM <- rowSums(weight)
    ROWSUM[ROWSUM == 0] <- 1
    weight <- weight/ROWSUM
    s <- sum(weight)
    m <- mean(x)
    y <- x - m
    cv <- sum(weight * y %o% y)
    v <- sum(y^2)
    obs <- (n/s) * (cv/v)
    if (scaled) {
        i.max <- (n/s) * (sd(rowSums(weight) * y)/sqrt(v/(n - 
            1)))
        obs <- obs/i.max
    }
    S1 <- 0.5 * sum((weight + t(weight))^2)
    S2 <- sum((apply(weight, 1, sum) + apply(weight, 2, sum))^2)
    s.sq <- s^2
    k <- (sum(y^4)/n)/(v/n)^2
    sdi <- sqrt((n * ((n^2 - 3 * n + 3) * S1 - n * S2 + 3 * s.sq) - 
        k * (n * (n - 1) * S1 - 2 * n * S2 + 6 * s.sq))/((n - 
        1) * (n - 2) * (n - 3) * s.sq) - 1/((n - 1)^2))
    alternative <- match.arg(alternative, c("two.sided", "less", 
        "greater"))
    pv <- pnorm(obs, mean = ei, sd = sdi)
    if (alternative == "two.sided") 
        pv <- if (obs <= ei) 
            2 * pv
        else 2 * (1 - pv)
    if (alternative == "greater") 
        pv <- 1 - pv
    list(observed = obs, expected = ei, sd = sdi, p.value = pv)
}

# prep data for autocorrelation
data_non_NA<-data[data$status!=9,]
dim_mp<-length(data_non_NA[,1])

dmtmp <-nearest.dist(x=data_non_NA[,c("easting","northing")], y=NULL, method="euclidian", delta=50, upper=NULL);          
dmtmp@entries<-rep(1,length(dmtmp@entries))# [dmt@entries!=0]<-1 # 1 only when dist_mat not 0
SBmp<- nearest.dist(x=cbind(data_non_NA$block_num,rep(0,dim_mp)), method="euclidian", upper=NULL,delta=0.1)
ASmp<-dmtmp-SBmp; # get 1 whereever the distances matrix is defined(under threshold) and not same block
ASmp<-as.spam(ASmp)

# autocorrelation with classic method
houses_XY2<- data_non_NA[c("easting","northing")]; # not unicode and number of bugs
values <- data_non_NA$status # get presence absence data 
temp=system.time(dnb <- dnearneigh(as.matrix(houses_XY2), 0, 50));
temp=system.time(lw <- nb2listw(dnb, zero.policy=TRUE));
temp=system.time(mt <- moran.test.perso(values, lw, zero.policy=TRUE,adjust.n=FALSE));
print(mt)

# # matrix of weight
# dmt1<-dmtmp
# dmt1@entries[dmt1@entries<100&dmt1@entries>120]<-0
# dmt1<-as.spam(dmt1)
# dmt2<-dmt1-SBmp;


temp2<-system.time(mI<-moran_perso(dmtmp,values));
# library(ape)
# Moran.I.perso(values,as.matrix(dmtmp)) # carefull, don't normalize weights

