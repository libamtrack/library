\name{AT.inverse.gamma.response}
\alias{AT.inverse.gamma.response}
\title{AT.inverse.gamma.response}
\description{Computes the inverse gamma response, i.e. the dose for a given response.

!!THIS WORKS UPTONOW FOR THE LIN-QUAD-MODEL ONLY, BUT SHOULD BE DEVELOPED TO A GENERAL ROUTINE!!}
\usage{AT.inverse.gamma.response( surv, alpha, beta)
}
\arguments{
  \item{surv}{ Response (survival)}
  \item{alpha}{ Alpha parameter in (1/Gy)}
  \item{beta}{ Beta parameter in (1/Gy2)}
}
\value{Dose in Gy}

\examples{
# Compute dose for 10% survival for HSG cells
AT.inverse.gamma.response( surv = 0.1, alpha = 0.2, beta = 0.05)
}
