library(dplyr)


file <-  "scalability/scalability_gene_cor_celltype.R"

for ( i in c(1:3)){
  
  # "200_200", "200_500", "200_1000", "500_200", "500_500", "500_1000", "1000_200", "1000_500", "1000_1000", 
  # "3000_200", "3000_500", "3000_1000", "5000_200", "5000_500", "5000_1000"
  for ( size in c("3000_200", "3000_500", "3000_1000", "5000_200", "5000_500", "5000_1000")){
    model <- "scDesign3_gau_ri"
    message("start ", model, ":", size, " repeat: ", i)
    setwd("/home/xiaoqi/spatial_simulationV2")
    res <- try(processx::run("scalability/memusg.py", 
                             args = c("Rscript", file,  size ))
               , silent = TRUE
    )
    
    if(assertthat::is.error(res)){
      stop("R code execution resulted in an error")
    } else {
      res <- stringr::str_split(res[["stderr"]], "\\n") %>% unlist()
      pos <- which(stringr::str_detect(res, "memusg: vmpeak: "))
      if(length(pos) == 1){
        val <- res[pos] %>%
          stringr::str_split("\\s") %>%
          unlist() %>%
          `[`(3) %>%
          as.numeric() 
        #  %>%  magrittr::divide_by(2^20)
      } else {
        stop("Something went wrong with memusg. Try debug mode.")
      }
    }
    
    setwd("scalability/result")
    write(c(val,size, model), ncolumns = 3,  file="memory.txt",append=TRUE, sep = ",")
    message("finish ", model, ":", size, " repeat: ", i)
  }
  
}
