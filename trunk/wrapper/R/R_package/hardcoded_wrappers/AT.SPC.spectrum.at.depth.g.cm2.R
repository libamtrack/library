AT.SPC.spectrum.at.depth.g.cm2 <- function(spc, depth.g.cm2)
{
    depth.g.cm2.in.spc    <- unique(spc$depth.g.cm2)
    closest.depth.idx     <- which(abs(depth.g.cm2.in.spc - depth.g.cm2) == min(abs(depth.g.cm2.in.spc - depth.g.cm2)))
    closest.depth.g.cm2   <- depth.g.cm2.in.spc[closest.depth.idx]
    depth.step            <- unique(spc$depth.step)[closest.depth.idx]

    cat( paste( "Closest depth [g/cm2] is ", 
                sprintf("%4.3f", closest.depth.g.cm2), 
                " (depth step ", 
                depth.step, 
                ") with offset ", 
                sprintf("%4.3f", depth.g.cm2 - closest.depth.g.cm2), 
                ".\n", 
                sep = ""))

   return(AT.SPC.spectrum.at.depth.step(spc, depth.step))
}
