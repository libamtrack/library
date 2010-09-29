types <- c("bool", "long", "int", "long", "double")
conversion.c <- c("int*\t", "int*\t", "int*\t", "long*\t", "float*")
conversion.R <- c("as.integer", "as.integer", "as.integer", "as.integer", "as.single")
extension <- c("_long", "_long", "_long", "_long", "_double")

conversion.table <- data.frame(	type = c(types, paste("const", types)),
						conversion.c = c(paste(conversion.c, "\t", sep = ""),
									 paste("const", conversion.c)),
						conversion.R = c(conversion.R, conversion.R),
						extension = c(extension, extension),
						stringsAsFactors = F
					)

# add additional "\t\t" to "float*"
conversion.table$conversion.c[conversion.table$conversion.c == "float*\t"] <- "float*\t\t" 


convert.c <- function(type){
				ii <- type == conversion.table$type
#				cat("type=_",type,"_ ",conversion.table$type,ii,"\n",sep="")
				if( sum(ii) > 0 )
					return(conversion.table$conversion.c[ii])
				else
					return(type)				
			}

convert.R <- function(type){
				ii <- type == conversion.table$type
				return(conversion.table$conversion.R[ii])
			}