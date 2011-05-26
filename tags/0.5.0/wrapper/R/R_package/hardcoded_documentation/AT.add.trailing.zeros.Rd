\name{AT.add.trailing.zeros}
\alias{AT.add.trailing.zeros}
\title{AT.add.trailing.zeros}
\description{Adds trailing zeros to a number and represent it in a string of fixed length.}
\usage{AT.add.trailing.zeros( x, digits = 5)
}
\arguments{
  \item{x}{ number}
  \item{digits}{length of resulting string.}
}
\value{
  Character string with \code{digits} characters.
}

\examples{
# Represent 99.1 as '99.100'
AT.add.trailing.zeros(x = 99.1, digits = 5)
}
