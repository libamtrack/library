require(libamtrack)
rm(list = ls())

path.to.files <- "/home/greilich/workspace/eclipse/libamtrack/data/SPC/auswahl/*.spc"

spc.list <- AT.SPC.get.list(path.to.files = path.to.files)

spc <- AT.SPC.get(spc.list = spc.list, energy.MeV.u = 255)

spc.lower <- AT.SPC.read("/home/greilich/workspace/eclipse/libamtrack/data/SPC/auswahl/libamtrack.12C.H2O.active3.MeV24000.spc", endian = "little")
spc.upper <- AT.SPC.read(spc.list$file.name[3], endian = "little")

ddd.lower <- AT.SPC.tapply(spc = spc.lower$spc, INDEX = c("depth.g.cm2"), 
                           FUN = AT.total.D.Gy, additional.arguments = list(c("material.no", 
            "AT.material.no.from.material.name('Water, Liquid')", 
            FALSE), c("stopping.power.source.no", "2", FALSE)), 
        names.results = "D.Gy")

ddd.interp <- AT.SPC.tapply(spc = spc$spc, INDEX = c("depth.g.cm2"), 
                           FUN = AT.total.D.Gy, additional.arguments = list(c("material.no", 
            "AT.material.no.from.material.name('Water, Liquid')", 
            FALSE), c("stopping.power.source.no", "2", FALSE)), 
        names.results = "D.Gy")

ddd.upper <- AT.SPC.tapply(spc = spc.upper$spc, INDEX = c("depth.g.cm2"), 
                           FUN = AT.total.D.Gy, additional.arguments = list(c("material.no", 
            "AT.material.no.from.material.name('Water, Liquid')", 
            FALSE), c("stopping.power.source.no", "2", FALSE)), 
        names.results = "D.Gy")

ddd.lower$which <- "lower"
ddd.interp$which <- "interp"
ddd.upper$which <- "upper"
ddd <- rbind.data.frame(ddd.lower, ddd.interp, ddd.upper)

require(lattice)

xyplot(D.Gy ~ depth.g.cm2,
       ddd,
       groups = which,
       auto.key = TRUE,
       type = 'o')

##################
# Check intermediate spectrum of interpolated spc -> no more double peaks 
# if energy / z grid is dense enough
spectrum <- AT.SPC.spectrum.at.depth.g.cm2(spc = spc$spc, depth.g.cm2 = 20.0) 

xyplot(log10(fluence.cm2) ~ E.MeV.u,
       spectrum,
       groups = particle.no,
       type = 's')

##################
# Check intermediate spectrum of original spc

first.depth <- unique(spc.upper$spc$depth.g.cm2[spc.upper$spc$depth.step == 2])
sec.depth <- unique(spc.upper$spc$depth.g.cm2[spc.upper$spc$depth.step == 3])
spc.depth.g.cm2 <- c(seq(0,first.depth,length.out = 6),
                     seq(first.depth*1.16, sec.depth, length.out = 6))
df <- NULL
for (i in 1:length(spc.depth.g.cm2)){
  spectrum <- AT.SPC.spectrum.at.depth.g.cm2(spc = spc.upper$spc, 
                                             depth.g.cm2 = spc.depth.g.cm2[i])
  df <- rbind(df, 
              data.frame(depth.g.cm2 = spc.depth.g.cm2[i],
                         spectrum))
}

custom.panel <- function(x,y,subscripts,...){
  panel.grid(h = -1, v = -1)
  panel.xyplot(x=x,y=y,subscripts=subscripts,...)
  sum.fluence <- sum(y, na.rm = TRUE)
  panel.text(sprintf("%3.2f",sum.fluence), x = 800, y = 0.3)
}

xyplot( fluence.cm2 ~ E.MeV.u|sprintf("depth %04.3f g/cm2", depth.g.cm2),
        df,
        groups = particle.no,
        auto.key = list(space = 'right'),
        layout = c(4,3),
        as.table = TRUE,
        panel = custom.panel,
        type = 's')

xyplot( log10(fluence.cm2) ~ E.MeV.u|sprintf("depth %04.3f g/cm2", depth.g.cm2),
        df,
        groups = particle.no,
        auto.key = list(space = 'right'),
        layout = c(4,3),
        as.table = TRUE,
        panel = custom.panel,
        type = 's')


##################

# Try spc mmap

require(libamtrack)
file <- "/home/greilich/workspace/eclipse/libamtrack/data/SPC/libamtrack.12C.H2O.active3.MeV27000.spc"
spc <- AT.SPC.read(file, 
                   endian = "little")
spc <- AT.SPC.read(file, 
                   flavour = "C",
                   endian = "little")