\name{pmh}
\alias{pmh}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
A measure of departure from partial marginal homogeneity for square contingency tables
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
When we want to measure departure from partial marginal homogeneity for square contingency tables, we can use this function.
}
\usage{
pmh(data = j_sight, lambda = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
%%     ~~Describe \code{data} here~~
A square contingency table
}
  \item{lambda}{
%%     ~~Describe \code{lambda} here~~
lambda
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
Again, when we want to measure departure from partial marginal homogeneity for square contingency tables, we can use this function.
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
Ytr
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (data = j_sight, lambda = 0)
{
    if (1 == 2) {
        lambda = 0
    }
    log0 = function(x) {
        if (is.infinite(log(x))) {
            res = 0
        }
        else {
            res = log(x)
        }
        return(res)
    }
    p_ij = dat/sum(dat)
    n.row = n.col = nrow(p_ij)
    dp_ij = matrix(NA, n.row, n.col)
    A = (2^lambda)/((2^lambda) - 1)
    A2 = 1/log0(2)
    p_id = c()
    p_di = c()
    p1rat = p2rat = c()
    him = enm = lenm = pistar = c()
    alln = sum(dat)
    toriihim = toriienm = 1
    for (i in 1:n.row) {
        p_id[i] = sum(p_ij[i, ])
        p_di[i] = sum(p_ij[, i])
        p1rat[i] = p_id[i]/(p_id[i] + p_di[i])
        p2rat[i] = p_di[i]/(p_id[i] + p_di[i])
        him[i] = 1 - A * (1 - ((p1rat[i])^(lambda + 1)) - ((p2rat[i])^(lambda +
            1)))
        enm[i] = 1 + A2 * (p1rat[i] * log0(p1rat[i]) + p2rat[i] *
            log0(p2rat[i]))
        pistar[i] = (p_id[i] + p_di[i])/2
        toriihim = toriihim * (him[i]^pistar[i])
        toriienm = toriienm * (enm[i]^pistar[i])
    }
    if (lambda != 0) {
        for (i in 1:n.row) {
            for (j in 1:n.col) {
                if (i == j) {
                  dp_ij[i, j] = toriihim * (log0(him[i]) + A *
                    (lambda + 1) * (p_id[i]^lambda - p_di[i]^lambda) *
                    (p_di[i] - p_id[i])/(2 * him[i] * (p_id[i] +
                    p_di[i])^(lambda + 1)))
                }
                else {
                  dp_ij[i, j] = (1/2) * toriihim * (log0(him[i]) +
                    A * (lambda + 1) * (p_id[i]^lambda - p_di[i]^lambda) *
                      p_di[i]/(him[i] * (p_id[i] + p_di[i])^(lambda +
                      1)) + log0(him[j]) + A * (lambda + 1) *
                    (p_di[j]^lambda - p_id[j]^lambda) * p_id[j]/(him[j] *
                    (p_id[j] + p_di[j])^(lambda + 1)))
                }
            }
        }
        Sigma1 = Sigma22 = 0
        for (i in 1:n.row) {
            for (j in 1:n.col) {
                Sigma1 = Sigma1 + p_ij[i, j] * (dp_ij[i, j]^2)
                Sigma22 = Sigma22 + p_ij[i, j] * dp_ij[i, j]
            }
        }
        Sigma2 = Sigma22^2
        Sigma = Sigma1 - Sigma2
        CI.low = toriihim - qnorm(1 - 0.05/2) * sqrt(Sigma/alln)
        CI.upper = toriihim + qnorm(1 - 0.05/2) * sqrt(Sigma/alln)
    }
    if (lambda == 0) {
        for (i in 1:n.row) {
            for (j in 1:n.col) {
                if (i == j) {
                  dp_ij[i, j] = toriienm * (log0(enm[i]) + A2 *
                    (p_di[i] - p_id[i]) * log0(p_id[i]/p_di[i])/(2 *
                    (p_id[i] + p_di[i]) * enm[i]))
                }
                else {
                  dp_ij[i, j] = toriienm * (1/2) * (log0(enm[i]) +
                    A2 * p_di[i] * log0(p_id[i]/p_di[i])/((p_id[i] +
                      p_di[i]) * enm[i]) + log0(enm[j]) - A2 *
                    p_id[j] * log0(p_id[j]/p_di[j])/((p_id[j] +
                    p_di[j]) * enm[j]))
                }
            }
        }
        Sigma1 = Sigma22 = 0
        for (i in 1:n.row) {
            for (j in 1:n.col) {
                Sigma1 = Sigma1 + p_ij[i, j] * (dp_ij[i, j]^2)
                Sigma22 = Sigma22 + p_ij[i, j] * dp_ij[i, j]
            }
        }
        Sigma2 = Sigma22^2
        Sigma = Sigma1 - Sigma2
        CI.low = toriienm - qnorm(1 - 0.05/2) * sqrt(Sigma/alln)
        CI.upper = toriienm + qnorm(1 - 0.05/2) * sqrt(Sigma/alln)
        toriihim = toriienm
    }
    return(list(dep = toriihim, vari = Sigma, CI = c(CI.low,
        CI.upper)))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
