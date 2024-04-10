

table_BIC <- function(BIC_vector, C_values){
  border_string <- "border-bottom: 2px solid;border-top: 2px solid;border-left: 2px solid;border-right: 2px solid;"
  ff = as.matrix(cbind(Clusters = paste("C =", C_values), BIC = round(BIC_vector,2)))
  cn <- colnames(ff)
  idx_best <- which.min(BIC_vector)
  colnames(ff) <- NULL
  ff %>% kbl(align = c("c", "c")) %>%
    kable_styling(font_size = 20)  %>%
    row_spec(idx_best, color = "black",
           background = "rgba(0, 255, 0, 0.6)", bold = F) %>%
    row_spec(seq(1,length(BIC_vector)), extra_css = "border-bottom: 2px solid;") %>%
    column_spec(1:2, extra_css = "border-left: 2px solid;") %>%
    add_header_above(cn, extra_css = border_string, 
                   background = alpha("lightblue",0.2))
  
}
