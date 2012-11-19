determinant.test <-
function (x, logarithm = TRUE, pivot = "MMD", method = "NgPeyton", 
    memory = list(), eps = 0.0000001, ...) 
{
    # if (eps < .Machine$double.eps) 
    #     stop("'eps' should not be smaller than machine precision", 
    #         call. = FALSE)
    logdet <- list()
    nrow <- x@dimension[1]
    nnzA <- as.integer(x@rowpointers[nrow + 1] - 1)
    if (nrow != x@dimension[2]) 
        stop("non-square matrix in 'chol'", call. = FALSE)
    # if (.Spam$cholsymmetrycheck) {
    #     test <- all.equal.spam(x, t.spam(x), tolerance = eps * 
    #         100)
    #     if (!identical(TRUE, test)) 
    #         stop("Input matrix to 'chol' not symmetric (up to 100*eps)", 
    #             call. = FALSE)
    # }
    if (method != "NgPeyton") 
        warning(gettextf("method = '%s' is not supported. Using 'NgPeyton'", 
            method), domain = NA)
    if (length(pivot) == nrow) {
        doperm <- as.integer(0)
        pivot <- as.vector(pivot, "integer")
        # if (.Spam$cholpivotcheck) {
        #     checkpivot(pivot, nrow)
        # }
    }
    else if (length(pivot) == 1) {
        if (pivot == FALSE) {
            doperm <- as.integer(0)
            pivot <- seq_len(nrow)
        }
        else if (pivot == TRUE) {
            doperm <- as.integer(1)
            pivot <- vector("integer", nrow)
        }
        else {
            doperm <- as.integer(switch(match.arg(pivot, c("MMD", 
                "RCM")), MMD = 1, RCM = 2))
            pivot <- vector("integer", nrow)
        }
    }
    else stop("'pivot' should be 'MMD', 'RCM' or a permutation")
    nnzcfact <- c(5, 1, 5)
    nnzRfact <- c(5, 1, 2)
    if (is.null(memory$nnzcolindices)) {
        nnzcolindices <- ifelse((nnzA/nrow < 5), max(1000, nnzA * 
            (1.05 * nnzA/nrow - 3.8)), nnzA) * nnzcfact[doperm + 
            1]
        nnzcolindices <- max(nnzcolindices, nnzA)
    }
    else {
        nnzcolindices <- max(memory$nnzcolindices, nnzA)
        memory$nnzcolindices <- NULL
    }
    if (is.null(memory$nnzR)) 
        nnzR <- min(max(4 * nnzA, floor(0.2 * nnzA^1.3)) * nnzRfact[doperm + 
            1], nrow * (nrow + 1)/2)
    else {
        nnzR <- memory$nnzR
        memory$nnzR <- NULL
    }
    if (is.null(memory$cache)) 
        cache <- 64
    else {
        cache <- memory$cache
        memory$cache <- NULL
    }
    if (length(memory) > 0) 
        warning("The component(s) ", paste("'", names(memory), 
            "'", sep = "", collapse = ","), " of the argument 'memory'\npassed to function 'chol' not meaningful and hence ignored.", 
            call. = FALSE)
    z <- .Fortran("cholstepwise", nrow = nrow, nnzA = as.integer(x@rowpointers[nrow + 
        1] - 1), d = dcheck(x@entries), jd = x@colindices, id = x@rowpointers, 
        doperm = doperm, invp = vector("integer", nrow), perm = pivot, 
        nnzlindx = vector("integer", 1), nnzcolindices = as.integer(nnzcolindices), 
        lindx = vector("integer", nnzcolindices), xlindx = vector("integer", 
            nrow + 1), nsuper = vector("integer", 1), nnzR = as.integer(nnzR), 
        lnz = vector("double", nnzR), xlnz = vector("integer", 
            nrow + 1), snode = vector("integer", nrow), xsuper = vector("integer", 
            nrow + 1), cachesize = as.integer(cache), ierr = as.integer(0), 
        NAOK = TRUE, DUP = FALSE)
    if (z$ierr == 1) 
        stop("Singularity problem when calculating the Cholesky factor.")
    if (z$ierr == 6) 
        stop("Inconsitency in the input", call. = FALSE)
    while (z$ierr > 1) {
        if (z$ierr == 4) {
            warning("Increased 'nnzR' with 'NgPeyton' method\n", 
                "(currently set to ", nnzR, " from ", ceiling(nnzR * 
                  .Spam$cholpar[1]), ")", call. = FALSE)
            nnzR <- ceiling(nnzR * .Spam$nnzRinc)
        }
        if (z$ierr == 5) {
            warning("Increased 'nnzcolindices' with 'NgPeyton' method\n", 
                "(currently set to ", nnzcolindices, " from ", 
                ceiling(nnzcolindices * .Spam$cholpar[2]), ")", 
                call. = FALSE)
            nnzcolindices <- ceiling(nnzcolindices * .Spam$cholpar[2])
        }
        z <- .Fortran("cholstepwise", nrow = nrow, nnzA = as.integer(x@rowpointers[nrow + 
            1] - 1), d = dcheck(x@entries), jd = x@colindices, 
            id = x@rowpointers, doperm = doperm, invp = vector("integer", 
                nrow), perm = pivot, nnzlindx = vector("integer", 
                1), nnzcolindices = as.integer(nnzcolindices), 
            lindx = vector("integer", nnzcolindices), xlindx = vector("integer", 
                nrow + 1), nsuper = vector("integer", 1), nnzR = as.integer(nnzR), 
            lnz = vector("double", nnzR), xlnz = vector("integer", 
                nrow + 1), snode = vector("integer", nrow), xsuper = vector("integer", 
                nrow + 1), cachesize = as.integer(cache), ierr = as.integer(0), 
            NAOK = !.Spam$safemode[3], DUP = FALSE)
    }
    if (z$ierr == 1) {
        warning("singularity problem or matrix not positive definite", 
            call. = FALSE)
        logdet$modulus <- NA
    }
    else {
        tmp <- 2 * sum(log(z$lnz[z$xlnz[-(z$nrow + 1)]]))
        if (logarithm) 
            logdet$modulus <- tmp
        else logdet$modulus <- exp(tmp)
    }
    attr(logdet$modulus, "logarithm") <- logarithm
    logdet$sign <- ifelse(z$ierr == 1, NA, 1)
    attr(logdet, "class") <- "det"
    return(logdet)
}
