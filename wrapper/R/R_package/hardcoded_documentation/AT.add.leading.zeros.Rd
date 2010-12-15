\name{AT.add.leading.zeros}
\alias{AT.add.leading.zeros}
\title{AT.add.leading.zeros}
\description{Adds leading zeros to a number and represent it in a string of fixed length.}
\usage{AT.add.leading.zeros( x, digits = 5)
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
AT.add.leading.zeros(x = 9, digits = 5)
}
