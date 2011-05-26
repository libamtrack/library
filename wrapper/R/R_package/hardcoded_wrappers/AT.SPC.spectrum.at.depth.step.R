AT.SPC.spectrum.at.depth.step <- function(spc, depth.step)
{
    if(depth.step %in% unique(spc$data$depth.step)){
        ii     <- spc$data$depth.step == depth.step
        df     <- data.frame( E.MeV.u     = spc$data$E.mid.MeV.u[ii],
                              particle.no = AT.particle.no.from.Z.and.A(Z = spc$data$Z[ii], A = spc$data$A[ii])$particle.no,
                              fluence.cm2 = spc$data$H[ii])
        #cat(paste("Returning", sum(ii), "entries in data frame.\n"))
        return(df)
    }else{
        cat("Depth step not found in data. Skipping.")
        return(NULL)
    }
}
