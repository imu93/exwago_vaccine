# This script aims to produce fn plots normalized by cpm
# depens on 01.bam2Rds and 02.get_fn_mtx.R

# For Neophytou et al., 2024 the whole script is not 100% necesary but 
# I'll save some lines to be able to accses to the code for non CPM norm plots
pacman::p_load(ggplot2, ggpubr, pals, reshape2, scales, dplyr, purrr, edgeR, tibble)
mtx = list.files(pattern = "acey.exwago.fn.all_lengths.Rds")
lst = lapply(mtx, readRDS)
lst = lst[[1]] 

groups = names(lst) %>% strsplit("_") %>%
  map(function(x){paste(x[c(1,2,3,4)], collapse = "_")}) %>%  unlist() %>% 
  unique()
names(groups) = groups

prop_l = list()
mean_l = list()
mean_tab = list()
group_lst = list()
simple_list = list()
for (i in names(groups)) {
  print(i)
  libs = lst[grepl(groups[i], names(lst))]
  #libs = lapply(libs, function(x){x[c(1:10),]}) # Edit this line to select specific rows
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
  index_fn = 15:40
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
  mean_fn = mean_fn[as.character(15:40),]
  rownames(mean_fn) = as.numeric(rownames(mean_fn))
  mean_fn[is.na(mean_fn)] = 0
  
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
  
  mean_fn_long$Length = factor(mean_fn_long$Length, levels = 15:40)
  prop_fn_long$Length = factor(prop_fn_long$Length, levels = 15:40)
  
  # Add colors
  colors = rev(c("#00008AFF","#09891A", "#E6110E","#FFE400"))
  #colors = rev(plasma(4))
  
  mean_fn_long$IP = groups[[i]]
  mean_p = ggplot(mean_fn_long, aes(x=Length, y=Reads, fill=First_nt)) + 
    geom_bar(stat = "identity") + theme_light() +  scale_fill_manual(values = colors) + 
    ylab("Number of reads") + xlab("") + 
    theme(axis.title.x = element_text(face = "bold", size = 16), axis.text.x = element_text(face = "bold", size = 16)) +
    theme(axis.title.y = element_text(face = "bold", size = 16),
          axis.text.y = element_text(face = "bold", size = 16, angle = 90, hjust = .5)) +
    scale_y_continuous(labels= comma) + theme( strip.text = element_text(size = 18, face = "bold", colour = "black")) +
    theme(plot.title = element_text(hjust = 0.5, size = 24)) + theme(legend.position="bottom") + facet_wrap(~groups[[i]]) + 
    theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
    scale_y_continuous(labels = unit_format(unit = "M", scale = 1e-6)) 
  
  mean_p = mean_p + labs(fill = "5' nt")
  mean_l[[groups[i]]] = mean_p
  mean_tab[[groups[i]]] = mean_fn_long
  
  prop_p = ggplot(prop_fn_long, aes(x=Length, y=Reads, fill=First_nt)) + 
    geom_bar(stat = "identity") + theme_light() +  scale_fill_manual(values = colors) + 
    ylab("% of reads") + xlab("Length distribution (nt)") + 
    theme(axis.title.x = element_text(face = "bold", size = 16), axis.text.x = element_text(face = "bold", size = 16, angle = 0, hjust = .5, vjust = 1)) +
    theme(axis.title.y = element_text(face = "bold", size = 16),
          axis.text.y = element_text(face = "bold", size = 16, angle = 90, hjust = .5)) +
    theme(plot.title = element_text(hjust = 0.5, size = 24)) + scale_y_continuous(breaks = c(0, 50, 100)) 
  
  prop_p = prop_p + labs(fill = "5' nt")
  prop_l[[groups[i]]] = prop_p
  # Save the merged plot in a list 
  group_lst[[groups[i]]] = ggarrange(mean_p, ggplot() + theme_minimal(), prop_p, ncol = 1, nrow = 3, heights = c(2, .05 , 1), common.legend = T)  
  #group_lst[[groups[i]]] = ggarrange(mean_p,  ncol = 1, nrow = 1, common.legend = T)  
}

ubs = list(group_lst$Adult_IP_polyP)
ggarrange(plotlist = ubs, ncol = 1, nrow = 1, common.legend = T, legend = "bottom")

prop_l_cpm = list()
mean_l_cpm = list()
mean_tab_cpm = list()
group_lst_cpm = list()
mean_fn_cpm_lis = list()
for (i in names(groups)) {
  print(i)
  # Get libraries from the same group
  libs = lst[grepl(groups[i], names(lst))]
  counts_lts = list()
  # I need a mtx with counts by fn and length 
  for (fn_lib in names(libs)) {
    tmp.lib = libs[[fn_lib]]
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
  index_fn = 15:40
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
  mean_fn = mean_fn[as.character(15:40),]
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
  
  mean_fn_long$Length = factor(mean_fn_long$Length, levels = 15:40)
  prop_fn_long$Length = factor(prop_fn_long$Length, levels = 15:40)
  
  # Add colors
  colors = rev(c("#00008AFF","#09891A", "#E6110E","#FFE400"))
  #colors = rev(plasma(4))
  
  mean_fn_long$IP = groups[[i]]
  mean_p = ggplot(mean_fn_long, aes(x=Length, y=Reads, fill=First_nt)) + 
    geom_bar(stat = "identity") + theme_light() +  scale_fill_manual(values = colors) + 
    ylab("CPM") + xlab("") + 
    theme(axis.title.x = element_text(face = "bold", size = 16), axis.text.x = element_text(face = "bold", size = 16)) +
    theme(axis.title.y = element_text(face = "bold", size = 16),
          axis.text.y = element_text(face = "bold", size = 16, angle = 90, hjust = .5)) +
    scale_y_continuous(labels= comma) + theme( strip.text = element_text(size = 18, face = "bold", colour = "black")) +
    theme(plot.title = element_text(hjust = 0.5, size = 24)) + theme(legend.position="bottom") + facet_wrap(~groups[[i]]) +
    scale_y_continuous(labels = unit_format(unit = "", scale = 1e-6)) +  ylim(c(0, .45e6))
  
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

# For manuscript figures I will only save the cpm table to then plot all libs
saveRDS(mean_fn_cpm_lis, "hBake_exWAGO_IP_raw_cpm.Rds")

