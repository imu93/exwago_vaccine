setwd("~/storage/Data/vaccine_exwago/analyses/ss3/fn_cls/")
pacman::p_load(ggplot2, ggpubr, pals, reshape2, scales, dplyr, purrr, edgeR, tibble)

inFile="heligmosomoides_bakeri.fn_cls_ustrand.1827.Rds"
outFile = paste0(sub("\\..*", "", inFile),".exWAGO_IP_raw_cpm.Rds")
mtx = readRDS(inFile)

fn  = mtx$fn
cls = mtx$cls 

groups = names(fn) %>% strsplit("_") %>%
  map(function(x){paste(x[c(1,2,3)], collapse = "_")}) %>%  unlist() %>% 
  unique()

names(groups) = groups

prop_l_cpm = list()
mean_l_cpm = list()
mean_tab_cpm = list()
group_lst_cpm = list()
mean_fn_cpm_lis = list()
simple_list = list()
for (i in names(groups)) {
  print(i)
  # Get libraries from the same group
  libs = fn[grepl(groups[i], names(fn))]
  counts_lts = list()
  # I need a mtx with counts by fn and length 
  for (fn_lib in names(libs)) {
    tmp.lib = libs[[fn_lib]] %>% as.data.frame()
    tmp.lib$len = rownames(tmp.lib)
    tmp.mlt = melt(tmp.lib) # Melt and collapse
    rownames(tmp.mlt) = paste0(tmp.mlt$variable, "_",tmp.mlt$len)
    tmp.mlt = select(tmp.mlt, c("value"))
    colnames(tmp.mlt) = fn_lib
    counts_lts[[fn_lib]] = tmp.mlt
  }
  # Now I need to merge and transform into CPMs
  # LogCPMs does not make sense in this context
  names(counts_lts) = NULL
  tmp.counts = do.call(cbind, counts_lts)
  tmp.counts[is.na(tmp.counts)] = 0
  mtx_cpm = cpm(tmp.counts)
  
  # Now I need to tranform into a mtx format
  libs = list()
  for (exp in colnames(mtx_cpm)) {
    x = mtx_cpm[,exp]
    cpm_mtx = as.data.frame(x) %>% 
      rownames_to_column(var = "rownames") %>%
      mutate(first_char = substr(rownames, 1, 1)) %>%
      group_split(first_char) %>%
      setNames(c("df_A", "df_C", "df_G", "df_T")) %>% 
      lapply(., function(x){select(x,c(1,2))}) %>% 
      do.call(cbind, .) 
    
    rownames(cpm_mtx) = sub(".*_", "", cpm_mtx$df_A.rownames)
    cpm_mtx = cpm_mtx[,!grepl("rownames", colnames(cpm_mtx))]
    colnames(cpm_mtx) =  gsub("^df_", "", substr( colnames(cpm_mtx), 4, 4)) 
    libs[[exp]] = cpm_mtx
  }
  
  # libs = lapply(libs, function(x){x[c(1:10),]}) # Edit this line to select specific rows
  # For very long matrices some rows could be emth, thus I need a way to have
  # all the same avlues in each table. To do so, first I need all possible values 
  vals = lapply(libs, rownames) %>%  melt() %>%  select(value) %>% pull() %>% factor() %>% 
    levels() %>%  as.numeric() %>%  sort() %>% as.character()
  n_libs = list()
  
  for (x in 1:length(libs)) {
    l = libs[[x]]
    ns = setdiff(vals, rownames(l))
    vs =  as.data.frame(matrix(0, length(ns), ncol(l)))
    colnames(vs) =  colnames(l)
    rownames(vs) = ns
    l = rbind(l, vs)
    l =  l[vals,]
    n_libs[[x]] = l
  }
  
  libs = n_libs
  # I'll use Reduce to obtain the mean of the matrices within the group list
  mean_fn = Reduce("+", libs) / length(libs)
  # Change T for U (these are sRNAs no DNA)
  colnames(mean_fn) = sub("T", "U", colnames(mean_fn))
  index_fn = 18:27
  n_size_list = list()
  for (read_size in index_fn) {
    if (read_size %in% rownames(mean_fn) == F) {
      n_size = c("C"=0,"U"=0, "A"=0, "G"=0)
      n_size_list[[as.character(read_size)]] = n_size
    }else{
      next
    }
    
  }
  ids = names(n_size_list)
  n_size = do.call(rbind, n_size_list) %>% data.frame()
  
  mean_fn = rbind(n_size, mean_fn)
  mean_fn = mean_fn[as.character(18:27),]
  rownames(mean_fn) = as.numeric(rownames(mean_fn))
  mean_fn[is.na(mean_fn)] = 0
  mean_fn_cpm_lis[[i]] = mean_fn
  
  mean_fn_long = mean_fn
  mean_fn_long$length = as.factor(rownames(mean_fn_long))
  mean_fn_long =  melt(mean_fn_long, id.vars = "length")
  # I need row means to estimate the average count number per base
  prop_fn = mean_fn/rowSums(mean_fn)*100
  prop_fn_long = prop_fn
  prop_fn_long$length = as.factor(rownames(prop_fn_long))
  # I'll use melt to format the mtx by length 
  prop_fn_long=  melt(prop_fn_long, id.vars = "length")
  colnames(mean_fn_long) = c("Length", "First_nt", "Reads")
  colnames(prop_fn_long) = c("Length", "First_nt", "Reads")
  # Now I need to order the nt using factors 
  # Note, edit in case of Ns
  mean_fn_long$First_nt = factor(mean_fn_long$First_nt, levels = rev(c("C", "A", "U", "G"))) 
  prop_fn_long$First_nt = factor(prop_fn_long$First_nt, levels = rev(c("C", "A", "U", "G")))
  
  mean_fn_long$Length = factor(mean_fn_long$Length, levels = 18:27)
  prop_fn_long$Length = factor(prop_fn_long$Length, levels = 18:27)
  
  # Add colors
  colors = rev(c("#00008AFF","#09891A", "#E6110E","#FFE400"))
  #colors = rev(plasma(4))
  
  mean_fn_long$IP = groups[[i]]
  mean_p = ggplot(mean_fn_long, aes(x=Length, y=Reads, fill=First_nt)) + 
    geom_bar(stat = "identity", show.legend = F) + theme_light() +  scale_fill_manual(values = colors) + 
    ylab("CPM") + 
    theme(axis.title.x = element_text(face = "bold", size = 16), axis.text.x = element_text( size = 16)) +
    theme(axis.title.y = element_text(face = "bold", size = 16),
          axis.text.y = element_text(size = 16, angle = 90, hjust = .5)) +
    scale_y_continuous(labels= comma) + theme( strip.text = element_text(size = 18, face = "bold", colour = "black")) +
    theme(plot.title = element_text(hjust = 0.5, size = 24)) + theme(legend.position="bottom") + facet_wrap(~groups[[i]]) +
    scale_y_continuous(labels = unit_format(unit = "", scale = 1e-6)) +  ylim(c(0, .4e6)) +
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank())
  
  mean_p = mean_p + labs(fill = "5' nt")
  mean_l_cpm[[groups[i]]] = mean_p
  mean_tab_cpm[[groups[i]]] = mean_fn_long
  
  prop_p = ggplot(prop_fn_long, aes(x=Length, y=Reads, fill=First_nt)) + 
    geom_bar(stat = "identity") + theme_light() +  scale_fill_manual(values = colors) + 
    ylab("% of CPM") + xlab("Length distribution (nt)") + 
    theme(axis.title.x = element_text(face = "bold", size = 16), axis.text.x = element_text(face = "bold", size = 16, angle = 0, hjust = .5, vjust = 1)) +
    theme(axis.title.y = element_text(face = "bold", size = 16),
          axis.text.y = element_text(face = "bold", size = 16, angle = 90, hjust = .5)) +
    theme(plot.title = element_text(hjust = 0.5, size = 24)) + scale_y_continuous(breaks = c(0, 50, 100)) 
  
  prop_p = prop_p + labs(fill = "5' nt")
  prop_l_cpm[[groups[i]]] = prop_p
  # Save the merged plot in a list 
  group_lst_cpm[[groups[i]]] = ggarrange(mean_p, ggplot() + theme_minimal(), prop_p, ncol = 1, nrow = 3, heights = c(2, .05 , 1), common.legend = T)  
  #group_lst[[groups[i]]] = ggarrange(mean_p,  ncol = 1, nrow = 1, common.legend = T) 
  mean_p = ggplot(mean_fn_long, aes(x=Length, y=Reads, fill=First_nt)) + 
    geom_bar(stat = "identity") + theme_light() +  scale_fill_manual(values = colors) + 
    ylab("CPM") + xlab("Length distribution (nt)") + 
    theme(axis.title.x = element_text(face = "bold", size = 16), axis.text.x = element_text(face = "bold", size = 16)) +
    theme(axis.title.y = element_text(face = "bold", size = 16),
          axis.text.y = element_text(face = "bold", size = 16, angle = 90, hjust = .5)) +
    scale_y_continuous(labels= comma) + theme( strip.text = element_text(size = 18, face = "bold", colour = "black")) +
    theme(plot.title = element_text(hjust = 0.5, size = 24)) + theme(legend.position="bottom") + facet_wrap(~IP) + 
    scale_y_continuous(labels = unit_format(unit = "", scale = 1e-6)) + ylim(c(0, .5e6))
  simple_list[[i]] = mean_p
}

mean_pl = list()
prop_pl = list()
mean_dfl = list()
prop_dfl = list()
for (i in names(groups)) {
  print(i)
  libs = cls[grepl(groups[i], names(cls))]
  # I'll use Reduce to obtain the mean of the matrices within the group list
  mean_fn = Reduce("+", libs) / length(libs) 
  mean_fn = mean_fn %>% as.data.frame()
  
  mean_fn_long = as.data.frame(mean_fn)
  mean_fn_long$length= as.factor(rownames(mean_fn_long))
  mean_fn_long =  melt(mean_fn_long, id.vars = "length")
  prop_fn = mean_fn/rowSums(mean_fn)*100
  prop_fn_long = as.data.frame(prop_fn)
  prop_fn_long$length = as.factor(rownames(prop_fn_long))
  # I'll use melt to format the mtx by length 
  prop_fn_long=  melt(prop_fn_long, id.vars = "length")
  colnames(mean_fn_long) = c("Length", "Class", "Reads")
  colnames(prop_fn_long) = c("Length", "Class", "Reads")
  
  #mean_fn_long$Class = factor(mean_fn_long$Class, levels = c("Intergenic", "Other_repeats", "Other_ncRNA",
  #                                                          "miRNA", "lincRNA", "rRNA","yRNA","tRNA","Unknown", "Retrotransposon", "Transposon",
  #                                                         "Intron","Exon"))
  
  #prop_fn_long$Class = factor(prop_fn_long$Class, levels = c("Intergenic", "Other_repeats", "Other_ncRNA",
  #                                                          "miRNA", "lincRNA", "rRNA","yRNA","tRNA", "Unknown", "Retrotransposon", "Transposon",
  #                                                         "Intron","Exon"))
  
  #colors = rev(c("#0000C9","#0051FFFF","#0092EFFF", "#13F240FF", "#9AFF16FF", "#FFBA00FF", "#FF7F00FF" ,"#FF3E00", "#BF0000", "#601573", "#4C4C4F","#76767A", "#DBDBDB"))
  mean_fn_long$Class = factor(mean_fn_long$Class, levels =rev(c("Protein-coding", "Transposon","ncRNA", "Other_repeat", "Unannotated")))
  
  prop_fn_long$Class = factor(prop_fn_long$Class, levels =rev(c("Protein-coding", "Transposon","ncRNA", "Other_repeat", "Unannotated")))
  
  colors = rev(c("#0B3BFA","#90DB3F", "#601573", "#A8A8A8", "#DBDBDB"))
  exp = groups[i]
  mean_fn_long$exp =  exp
  mean_p = ggplot(mean_fn_long, aes(x=Length, y=Reads, fill=Class)) + 
    geom_bar(stat = "identity") + theme_test() +  scale_fill_manual(values = colors, guide = guide_legend(reverse = TRUE)) + 
    ylab("Number of aligned reads") + xlab("Length distribution (nt)") + 
    theme(axis.title.x = element_text(face = "bold", size = 14), axis.text.x = element_text(size = 14)) +
    theme(axis.title.y = element_text(face = "bold", size = 14),
          axis.text.y = element_text(size = 14, angle = 90, hjust = .5)) + facet_wrap(~exp) +
    theme( strip.text = element_text(size = 14, colour = "black")) +
    theme(plot.title = element_text(hjust = 0.5, size = 24)) + theme(legend.position="bottom", legend.text = element_text(size = 14)) +
    scale_y_continuous(labels = unit_format(unit = "M",  scale = 1e-6)) 
  
  
  mean_p = mean_p +  xlab("") + theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())  +
    theme(legend.position = "none") 
  
  mean_pl[[groups[i]]] = mean_p
  
  
  prop_p = ggplot(prop_fn_long, aes(x=Length, y=Reads, fill=Class)) + 
    geom_bar(stat = "identity", show.legend = F) + theme_light() +  scale_fill_manual(values = colors, guide = guide_legend(reverse = TRUE)) + 
    theme(legend.position="bottom") + ylab("% of reads") + xlab("Length distribution (nt)")+ 
    theme(axis.title.x = element_text(face = "bold", size = 18), axis.text.x = element_text(size = 16)) +
    theme(axis.title.y = element_text(face = "bold", size = 18), axis.text.y = element_text(size = 16, angle = 90, hjust = .5)) +
    scale_y_continuous(breaks = c(0, 50, 100)) +  theme(plot.title = element_text(hjust = 0.5, size = 24))  +
    theme(legend.key.size = unit(.9, 'cm')) 
  
  
  prop_pl[[groups[i]]] = prop_p

}

plt_lst = list() 
for (i in names(groups)) {
  tmp.fn = mean_l_cpm[[i]]
  tmp.cls = prop_pl[[i]]
  plt_lst[[i]] = ggarrange(tmp.fn, ggplot() + theme_minimal(),tmp.cls, ncol = 1, heights = c(1.2, .1, .7), align = c("hv"))
}


#saveRDS(plt_lst$Adult_IP_polyP, "hbakeri_adult_IP.Rds")


#ms_plot = ggarrange(plt_lst$Sup_IP_polyP, plt_lst$EV_IP_polyP, common.legend = T)
#ggsave(filename = "hbakeri_exwago_secreted_ips_fn.pdf", device = "pdf", plot = ms_plot, 
#       width = 10, height = 7, dpi = 300,
#       path = "~/storage/Data/vaccine_exwago/analyses/figures/")
