\name{AT.SGR.leading.zeros}
\alias{AT.SGR.leading.zeros}
\title{AT.SGR.leading.zeros}
\description{Adds leading zeros to a number and represent it in a string of fixed length.}
\usage{AT.SGR.leading.zeros( x, digits = 5)
}
\arguments{
  \item{x}{ number}
  \item{digits}{length of resulting string. Should be larger that log10(x).}
}
\value{
  Character string with \code{digits} characters.
}

\examples{
# Represent 99 as '00099'
AT.SGR.leading.zeros(x = 9, digits = 5)
}
