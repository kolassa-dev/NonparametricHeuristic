#' Construct a analysis of variance table for multiple regression, highlighting
#' the null hypothesis of no effect from any of the covariates.  This table
#' is common in statistical text books.
#'
#' @param object A fit as created by lm.
#' @return A matrix with three rows, corresponing to model, error and total,
#'   and five columns, representing sum of squares, degrees of freedom, mean
#'   squares, the F statistic, and the p-value.
#'
#' @importFrom stats pf coef var
#' @export
traditionalanova <- function (object) {
    z <- object
    p <- z$rank
    rdf <- z$df.residual
    if (p == 0) {
        r <- z$residuals
        n <- length(r)
        w <- z$weights
        if (is.null(w)) {
            rss <- sum(r^2)
        } else {
            rss <- sum(w * r^2)
            r <- sqrt(w) * r
        }
        resvar <- rss/rdf
        ans <- z[c("call", "terms", if(!is.null(z$weights)) "weights")]
        class(ans) <- "summary.lm"
        ans$aliased <- is.na(coef(object))  # used in print method
        ans$residuals <- r
        ans$df <- c(0L, n, length(ans$aliased))
        ans$coefficients <- matrix(NA_real_, 0L, 4L, dimnames =
   list(NULL, c("Estimate", "Std. Error", "t value", "Pr(>|t|)")))
        ans$sigma <- sqrt(resvar)
        ans$r.squared <- ans$adj.r.squared <- 0
        ans$cov.unscaled <- matrix(NA_real_, 0L, 0L)
        return(ans)
    }
    if (is.null(z$terms)) stop("invalid 'lm' object:  no 'terms' component")
    if(!inherits(object, "lm")) warning("calling summary.lm(<fake-lm-object>) ...")
    Qr <- object$qr
    n <- NROW(Qr$qr)
    if(is.na(z$df.residual) || n - p != z$df.residual)
        warning("residual degrees of freedom in object suggest this is not an \"lm\" fit")
    ## do not want missing values substituted here
    r <- z$residuals
    f <- z$fitted.values
    w <- z$weights
    if (is.null(w)) {
        mss <- if (attr(z$terms, "intercept"))
            sum((f - mean(f))^2) else sum(f^2)
        rss <- sum(r^2)
    } else {
        mss <- if (attr(z$terms, "intercept")) {
            m <- sum(w * f /sum(w))
            sum(w * (f - m)^2)
        } else sum(w * f^2)
        rss <- sum(w * r^2)
        r <- sqrt(w) * r
    }
    resvar <- rss/rdf
    ## see thread at https://stat.ethz.ch/pipermail/r-help/2014-March/367585.html
    if (is.finite(resvar) &&
        resvar < (mean(f)^2 + var(c(f))) * 1e-30)  # a few times .Machine$double.eps^2
        warning("essentially perfect fit: summary may be unreliable")
    p1 <- 1L:p
    R <- chol2inv(Qr$qr[p1, p1, drop = FALSE])
    se <- sqrt(diag(R) * resvar)
    est <- z$coefficients[Qr$pivot[p1]]
    tval <- est/se
    if (p != attr(z$terms, "intercept")) {
 df.int <- if (attr(z$terms, "intercept")) 1L else 0L
    } else ans$r.squared <- ans$adj.r.squared <- 0
    out<-array(NA,c(3,5))
    dimnames(out)<-list(c("Model","Residual","Total"),
       c("Sum of Squares","Degrees of Freedom","Mean Squares","F","P"))
    out[,1]<-c(mss,rss,mss+rss)
    out[,2]<-c((p - df.int),rdf,rdf+p-df.int)
    out[1:2,3]<- out[1:2,1]/out[1:2,2]
    out[1,4]<-  (mss/(p - df.int))/resvar
    out[1,5]<-  pf(out[1,4],out[1,2],out[2,2],lower.tail=FALSE)
    return(out)
}
