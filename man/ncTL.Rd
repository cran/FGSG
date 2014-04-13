\name{ncTL}
\alias{ncTL}
\title{Non Convex Truncated L1 Feature Grouping and Selection}
\description{
The ncTL function implements a variant of non-convex truncated feature grouping and selection algorithm in [2], which only uses the truncated \eqn{l_1} regularization.
Given \eqn{A = \{a_1,\dots,a_n\}}, the response \eqn{y}, and a set of edges \eqn{E}, this function aims to solves
\deqn{min\frac{1}{2}||Ax-y||^2 + \lambda_1 \sum\limits_i(J_\tau(|x_i|)) + \lambda_2 \sum\limits_{(i,j)\in E}w_{(i,j)}|x_i-r_{(i,j)}x_j| }
where \eqn{J_\tau(x)=min(x/\tau,1)} is a surrogate of the L0 norm, \eqn{w_{(i,j)}} is the weight of the edge \eqn{(i,j)}, and \eqn{r_{(i,j)}} is the sign of the correlation between features \eqn{a_i} and \eqn{a_j}. The weight and sign can be specified in \eqn{Rwt}: \eqn{w=abs(Rwt)}, and \eqn{r=sign(Rwt)}.
}

\usage{
ncTL(A, y, tp, s1, s2, RmaxIter = 100, RaMaxIter = 1000, Rrho = 5, Rtau = 0.15,
Rwt = rep(1, length(tp)), Rtol = 0.001, RaTol = 0.001, x0 = rep(0, ncol(A)))
}
\arguments{
  \item{A}{
\code{A} The data matrix of size \eqn{n \times p}, each row corresponds to one sample.
}
  \item{y}{
\code{y} The response vector of length n.
}
  \item{tp}{
\code{tp} The edges vector of length 2g (eg. (1,2,3,4) means an edge between 1 and 2, and an edge between 3 and 4, g=2 is the number of edges).
}
  \item{s1}{
\code{s1} The \eqn{l_1} regularization parameter, \eqn{s1 >=0}.
}
  \item{s2}{
\code{s2} Tge grouping penatly parameter, \eqn{s2 >=0}.
}
  \item{RmaxIter}{
\code{RmaxIter} The maximum number of iterations in DC programming (default 100).
}
  \item{RaMaxIter}{
\code{RaMaxIter} The maximum number of iterations in ADMM (default 1000).
}
  \item{Rrho}{
\code{Rrho} The dual update length ofor ADMM (default 5).
}
  \item{Rtau}{
\code{Rtau} The tuning parameter for non-convex penalty (default 0.15).
}
  \item{Rwt}{
\code{Rwt} The weight and signs of edges (default rep(1,g)).
}
  \item{Rtol}{
\code{Rtol} The tolerance for convergence in DC programming (default 1e-3).
}
  \item{RaTol}{
\code{RaTol} The tolerance for convergence in ADMM (default 1e-3).
}
  \item{x0}{
\code{x0} The returned weight vector (default rep(0,p)).
}
}
\value{
  Returned value x0 is the solution to the optimizaiton problem.
}
\examples{
A<-matrix(rnorm(25),5,5)
y<-rnorm(5)
tp<-c(1,2,2,3,3,4,4,5)
ncTL(A,y,tp,0,0)
}
\references{
 S.Yang, L.Yuan, Y.Lai, X.Shen, P.Wonka, and J.Ye. Feature grouping and selection over an undirected graph. KDD, 2012.
}
\author{
  Yiwen Sun
}