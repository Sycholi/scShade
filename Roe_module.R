meta = seu@meta.data
library(Startrac)
Roe = calTissueDist(meta,
                    byPatient = F,
                    colname.cluster = "type", 
                    colname.patient = "Sample", 
                    colname.tissue = "Group", 
                    method = "chisq", # "chisq", "fisher", and "freq" 
                    min.rowSum = 0)  

roe_res = data.frame(T = Roe[,2], P = Roe[,1])
rownames(roe_res) = rownames(Roe)

p1 = pheatmap(roe_res,
              color = colorRampPalette(c("#ADC6A9", "#7DB391",'#469F80',
                                         '#16887E','#136C72','#1B5465',
                                         '#1B3955','#151D44'))(100) %>% rev(),
              display_numbers = T,
              border_color = NA,
              treeheight_row = 0,
              treeheight_col = 0)
ggsave(p1,
       filename = 'Roe_of_Cell_type.pdf',
       height = 1000,
       width = 1000,
       units = 'px',
       dpi = 300)
