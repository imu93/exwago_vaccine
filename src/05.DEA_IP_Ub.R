setwd("~/storage/Data/vaccine_paper/analyses/ss3/dea/")
pacman::p_load(edgeR, ggplot2, ggpubr, ggplotify, dplyr, purrr, gplots, rtracklayer, plyranges, reshape, kableExtra, ggExtra,
               ggdist, ggbreak, Rmisc, stringr, pals, ggdist)
list.files()
#####################################################################################################################################################################
conts = list(c("EV_IP", "EV_Ub"), c("Sup_IP", "Sup_Ub"), c("Adult_IP_polyP", "Adult_total_polyP"))
names(conts) = lapply(conts, function(x){paste(x, collapse = "-")})
l_lfc = c(1, 1, 2)
l_fdr = c(.05,.05,.01)
names(l_lfc) = names(conts)
names(l_fdr) = names(conts)
dt_lst_nonorm = list()
dt_lst_normrrna = list()
control_regs_lst = list()
control_norm_regs_lst = list()
nonorm_rrna_ma = list()
vplots = list()
tmm_ma = list()
rrna_ma = list()
boxlist = list()
norm_rrna_simp_list = list()
norm_rrna_linc_list = list()
contro_norm_ma_with_dens = list()
rrna_trna_mas = list()
cameraRes_lis = list()
cameraRes_fam_lis = list()
for ( cont in names(conts)) {
  #cont = "EV_IP-EV_Ub" 
  print(cont)
  contrast_id = cont
  libs = conts[[cont]] %>% paste(collapse = "|")
  # Read table and create groups 
  counts = read.delim("nxHelBake1.1.counts.1827.txt.gz")
  counts = counts[,grepl(libs, colnames(counts))]
  colnames(counts) = sub("\\.trim.*", "", colnames(counts))
  groups = colnames(counts) %>% strsplit("_") %>%
    map(function(x){paste(x[c(1,2)], collapse = "_")}) %>%  unlist() %>% 
    factor()
  ################################################################################################################################################  
  lib.size = colSums(counts)
  # How many reads do we have?
  round((lib.size)/1e6, 2)
  # Create DGE and filter by expression
  dge = DGEList(counts=counts, group=groups, lib.size = lib.size)
  keep = filterByExpr(dge, min.count= 5)
  dge = dge[keep, , keep.lib.sizes=FALSE]
  print(table(keep))
  f.exp.counts = colSums(dge$counts)
  dge = calcNormFactors(dge, method = "TMM")
  ##############################################################################################################################################
  if (grepl("EV|Sup", libs)) {
    # As Kyriaki mentioned we should use a batch effect since this are paired worms
    batch = factor(paste0("batch",sub(".*(r.*)","\\1",rownames(dge$samples))))
    design = model.matrix(~0+dge$samples$group+batch)
    colnames(design) =  c(levels(dge$samples$group),levels(batch)[-1])
  }else{
    design = model.matrix(~0+dge$samples$group)
    colnames(design) =  c(levels(dge$samples$group))
  }
  ##############################################################################################################################################
  dge = estimateDisp(dge, design=design, robust=TRUE)
  plotBCV(dge)
  # As expected, we have variance in low expressed regions and low in highly expressed 
  # Let's fit glm
  ncont = colnames(design)[c(1,2)] %>% paste(collapse = "-")
  fit = glmFit(dge, design = dge$design, dispersion = dge$common.dispersion)
  contrast = makeContrasts(ncont,levels = dge$design)
  cont_fc = l_lfc[cont]
  cont_fdr = l_fdr[cont]
  gt = glmTreat(fit,  contrast=contrast, lfc = cont_fc )
  topTags(gt)
  dt = decideTests(gt, adjust.method="BH",  p.value=cont_fdr, lfc = cont_fc)
  tb_nonorm = table(dt)
  dt_lst_nonorm[[ncont]] = tb_nonorm
  
  topTable = topTags(gt, n=Inf)$table
  topTable = topTable[rownames(dt),]
  de_tbl = data.frame('Ave_CPM'=topTable$logCPM, 'log-FC'= topTable$logFC, 'cl'= dt)
  colnames(de_tbl) = c('Ave_CPM', 'logFC', 'cl')
  de_tbl$col = ifelse(de_tbl$cl == 0, "#C0C0C0", ifelse(de_tbl$cl == 1, "#5E2129", "#333333"))
  de_tbl$ord = ifelse(de_tbl$col == "#C0C0C0", 1, ifelse(de_tbl$col == "#333333", 2, 3))
  de_tbl = de_tbl[order(de_tbl$ord),]
  outDEtbl = "nxHelBake1.1.other_TMM_de_tbl.txt"
  outDEtbl = sub("other", ncont, outDEtbl)
  write.table(de_tbl,outDEtbl, quote = F, sep = "\t")
  ##############################################################################################################################################
  de_tbl$cond = cont 
  de_tbl$cond = ifelse(grepl("Adult", cont), "Adult exWAGO IP vs Adult total", ifelse(grepl("EV", cont), "Vesicular exWAGO IP vs Unbound", "Non-vesicular exWAGO IP vs Unbound"))
  de_tbl$col = factor(de_tbl$col, levels = c("#5E2129", "#C0C0C0", "#333333"))
  
  p1 = ggplot(de_tbl, aes(x=Ave_CPM, y=logFC, colour = col)) +
    geom_point() +  xlab(bquote(~log[2]~' CPM')) + theme_light() +
    ylab(bquote(~log[2]~ 'FC IP/Ub')) + theme(axis.title.x = element_text(size = 20),
                                              axis.text.x = element_text( size = 24)) +
    theme(axis.title.y = element_text(size = 20),axis.text.y = element_text(size = 24)) +
    geom_hline(aes(yintercept = 0), colour = "blue", linewidth = 1, linetype="dashed", alpha=0.4) +
    theme(plot.title = element_text(hjust = 0.5, size = 24, face = "bold")) + 
    scale_color_manual(values = adjustcolor(c("#002A86", "#C0C0C0", "#333333"), alpha.f = .6), 
                       labels=c('IP enriched', 'not-DE', 'Unbound'), name=NULL) +
    theme(legend.position = "bottom", legend.text = element_text(size = 13)) +  
    guides(color = guide_legend(override.aes = list(size=10))) + ylim(c(-18,12)) +
    theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), 
                             "inches"))
  p1 
  
  ##############################################################################################################################################
  cts = c("miRNA_S", "rRNA_S", "lincRNA_S", "tRNA_S")
  df_lst = list()
  for (i in cts) {
    df = de_tbl
    df$col = ifelse(grepl(i, rownames(df)), "#144675", "#C0C0C0")
    df$ord = ifelse(df$col == "#C0C0C0", 1,2)
    df$class = i
    df = df[order(df$ord),]
    df$ord = factor(df$ord)
    df_lst[[i]] = df
  }
  names(df_lst) = NULL
  df = do.call(rbind, df_lst)
  control = ggplot(df, aes(x=Ave_CPM, y=logFC, colour = col)) +
    geom_point(show.legend = F) +  xlab(bquote('Average'~log[2]~' CPM')) + theme_minimal() + facet_wrap(~class, ncol = 2) + theme_light() +
    theme( strip.text = element_text(size = 11, face = "bold", colour = "black")) +
    ylab(bquote(~log[2]~ ' Fold-change')) + theme(axis.title.x = element_text(face = "bold", size = 24), axis.text.x = element_text(face = "bold", size = 24)) +
    theme(axis.title.y = element_text(face = "bold", size = 24),axis.text.y = element_text(face = "bold", size = 24)) +
    geom_hline(aes(yintercept = 0), colour = "blue", linewidth = 1, linetype="dashed", alpha=0.4) +
    theme(plot.title = element_text(hjust = 0.5, size = 20, face = "bold")) + 
    scale_color_manual(values = adjustcolor(c("#144675", "#BCBBBD"), alpha.f = .6)) + ylim(c(-18,18))
  
  control_regs_lst[[ncont]] = control
  
  df_rrna = df[df$class == "rRNA_S",]
  df_rrna$col = ifelse(grepl("28s_rRNA_S",rownames(df_rrna)), "#22B6E3", "#BCBBBD")
  df_rrna$col = ifelse(grepl("18s_rRNA_S",rownames(df_rrna)), "#223CE3", df_rrna$col)
  df_rrna$col = ifelse(grepl("^8s_rRNA_S",rownames(df_rrna)), "#7879E3", df_rrna$col)
  df_rrna$col = factor(df_rrna$col, levels = c("#BCBBBD", "#22B6E3", "#223CE3", "#7879E3"))
  
  rrna_nonorm = ggplot(df_rrna, aes(x=Ave_CPM, y=logFC, colour = col)) +
    geom_point(show.legend = T) +  xlab(bquote('Average'~log[2]~' CPM')) + theme_minimal() + facet_wrap(~class, ncol = 2) + theme_light() +
    theme( strip.text = element_text(size = 11, face = "bold", colour = "black")) +
    ylab(bquote(~log[2]~ ' Fold-change IP/Ub')) + theme(axis.title.x = element_text(face = "bold", size = 24), axis.text.x = element_text(face = "bold", size = 24)) +
    theme(axis.title.y = element_text(face = "bold", size = 24),axis.text.y = element_text(face = "bold", size = 24)) +
    geom_hline(aes(yintercept = 0), colour = "blue", linewidth = 1, linetype="dashed", alpha=0.4) +
    theme(plot.title = element_text(hjust = 0.5, size = 20, face = "bold")) + 
    scale_color_manual(values = adjustcolor(c("#BCBBBD", "#22B6E3", "#223CE3", "#7879E3"), alpha.f = .6), labels=c("other","28s", "18s", "8s")) + 
    ylim(c(-15,18)) +
    theme(legend.position = "bottom", legend.text = element_text(size = 14)) +
    guides(color = guide_legend(override.aes = list(size=10))) + ylim(c(-18, 18))
  
  nonorm_rrna_ma[[ncont]] = rrna_nonorm
  #####################################################################################################################################################################
  # Since rRNA_S are at the top of the second pic I'll use this category to estimate new normfactors
  normCounts = dge$counts[grepl("rRNA_S", rownames(dge$counts)),]
  sumNormRNAs = colSums(normCounts)
  totalLibsize = colSums(dge$counts)
  normFactors = sumNormRNAs / totalLibsize
  normFactors = normFactors / (prod(normFactors)^(1/length(normFactors)))
  dge$samples$norm.factors = normFactors
  dge = estimateDisp(dge, design=design, robust=TRUE)
  plotBCV(dge)
  #####################################################################################################################################################################
  fit = glmFit(dge, design=design, robust=T, dispersion = dge$common.dispersion)
  contrast = makeContrasts(ncont,levels = dge$design)
  cont_fc = l_lfc[cont]
  cont_fdr = l_fdr[cont]
  gt = glmTreat(fit,  contrast=contrast, lfc = cont_fc )
  topTags(gt)
  dt = decideTests(gt, adjust.method="BH", p.value=cont_fdr, lfc = cont_fc)
  normrrna_dt = table(dt)
  dt_lst_normrrna[[ncont]] = normrrna_dt
  topTable = topTags(gt, n=Inf)$table
  topTable = topTable[rownames(dt),]
  outTab = "f_tables/nxHelBake1.1.exwago_other_1827nt_lfc1_topTable.txt"
  outTab = sub("other", ncont, outTab)
  write.table(topTable, outTab , 
              sep = "\t", quote = F, col.names = T, row.names = T)
  #####################################################################################################################################################################
  de_tbl_norm = data.frame('Ave_CPM'=topTable$logCPM, 'log-FC'= topTable$logFC, 'cl'= dt)
  colnames(de_tbl_norm) = c('Ave_CPM', 'logFC', 'cl')
  de_tbl_norm$col = ifelse(de_tbl_norm$cl == 0, "#C0C0C0", ifelse(de_tbl_norm$cl == 1, "#15C1FB", "#C0C0C0"))
  de_tbl_norm$col = ifelse(grepl("^8s_rRNA_S", rownames(de_tbl_norm)), "#DB1A27", de_tbl_norm$col)
  de_tbl_norm$col = ifelse(grepl("18s_rRNA_S", rownames(de_tbl_norm)), "#DB971A", de_tbl_norm$col)
  de_tbl_norm$col = ifelse(grepl("28s_rRNA_S", rownames(de_tbl_norm)), "#0B5727", de_tbl_norm$col)
  de_tbl_norm$col = factor(de_tbl_norm$col, levels = c("#0B5727","#DB971A", "#DB1A27", "#15C1FB", "#C0C0C0"))
  de_tbl_norm$ord = ifelse(de_tbl_norm$col == "#C0C0C0", 1, ifelse(de_tbl_norm$col == "#333333", 2, 3))
  de_tbl_norm$ord = ifelse(grepl("rRNA_S", rownames(de_tbl_norm)), 4, de_tbl_norm$ord)
  de_tbl_norm = de_tbl_norm[order(de_tbl_norm$ord),]
  recol = c("#22B6E3", "#0064FA", "#7879E3", "#002A86","#C0C0C0")
  outTab2 = "f_tables/nxHelBake1.1.exwago_other_1827nt_lfc1_deTblnorm.txt"
  
  outTab2 = sub("other", ncont, outTab2)
  write.table(de_tbl_norm, outTab2, 
              sep = "\t", quote = F, col.names = T, row.names = T)
  ######################################################################################################################################################################
  de_tbl_norm_simp = de_tbl_norm
  de_tbl_norm_simp$cond = ifelse(grepl("Adult", cont), "Adult exWAGO IP vs Adult total", ifelse(grepl("EV", cont), "Vesicular exWAGO IP vs Unbound", "Non-vesicular exWAGO IP vs Unbound"))
  de_tbl_norm_simp$col = ifelse(grepl("rRNA_S", rownames(de_tbl_norm_simp)), "#DB971A", de_tbl_norm_simp$col)
  de_tbl_norm_simp$col = factor(de_tbl_norm_simp$col, levels=c(4,5,"#DB971A") )
  
  recol_2 = c("#002A86", "#C0C0C0", "#FA0C57")
  ma_smp = ggplot(de_tbl_norm_simp, aes(x=Ave_CPM, y=logFC, colour = col)) +
    geom_point() +  xlab(bquote(~log[2]~' CPM') ) +
    scale_color_manual(values = adjustcolor(recol_2, alpha.f = .9) , name=NULL, labels =c("IP enriched","Unbound", "rRNA_S")) + theme_light() +
    ylab(bquote(~log[2]~ 'FC IP/Ub')) + theme(axis.title.x = element_text(size = 20), axis.text.x = element_text(size = 24)) +
    theme(axis.title.y = element_text(size = 20),axis.text.y = element_text( size = 24)) +
    geom_hline(aes(yintercept = 0), colour = "blue", linewidth = 1, linetype="dashed", alpha=0.4) +
    theme(plot.title = element_text(hjust = 0.5, size = 24, face = "bold")) +  
    theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold")) + theme(legend.position = "bottom", legend.text = element_text(size = 14)) +
    guides(color = guide_legend(override.aes = list(size=10))) +  ylim(c(-18,18)) +
    theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), 
                             "inches"))
    #facet_wrap(~cond) + theme( strip.text = element_text(size = 16, face = "bold", colour = "black"))
  
  
  norm_rrna_simp_list[[ncont]] = ma_smp
  #gm1 = ggMarginal(ma_smp, groupColour = T, groupFill = T)
  
  #####################################################################################################################################################################
  de_tbl_norm_tRNA =  data.frame('Ave_CPM'=topTable$logCPM, 'log-FC'= topTable$logFC, 'cl'= dt)
  de_tbl_norm_tRNA$cond = ifelse(grepl("Adult", cont), "Adult exWAGO IP vs Adult total", ifelse(grepl("EV", cont), "Vesicular exWAGO IP vs Unbound", "Non-vesicular exWAGO IP vs Unbound"))
  colnames(de_tbl_norm_tRNA) = c('Ave_CPM', 'logFC', 'cl')
  de_tbl_norm_tRNA$col = ifelse(de_tbl_norm_tRNA$cl == 0, "#C0C0C0", ifelse(de_tbl_norm_tRNA$cl == 1, "#15C1FB", "#C0C0C0"))
  de_tbl_norm_tRNA$col = ifelse(grepl("tRNA_S",rownames(de_tbl_norm_tRNA)), "#DB1A27", de_tbl_norm_tRNA$col)
  de_tbl_norm_tRNA$col = factor(de_tbl_norm_tRNA$col, levels = c("#15C1FB", "#C0C0C0", "#DB1A27"))
  de_tbl_norm_tRNA$ord = ifelse(de_tbl_norm_tRNA$col == "#C0C0C0", 1, ifelse(de_tbl_norm_tRNA$col == "#15C1FB", 2,
                                                                             ifelse(de_tbl_norm_tRNA$col == "#DB1A27", 3, 0)))
  
  de_tbl_norm_tRNA = de_tbl_norm_tRNA[order(de_tbl_norm_tRNA$ord),]
  recol_2 = c("#002A86","#C0C0C0", "#E08707")
  ma_tRNA = ggplot(de_tbl_norm_tRNA, aes(x=Ave_CPM, y=logFC, colour = col)) +
    geom_point() +  xlab(bquote(~log[2]~' CPM') ) +
    scale_color_manual(values = adjustcolor(recol_2, alpha.f = .9) , name=NULL, labels =c("IP enriched", "Unbound", "tRNA_S")) + theme_light() +
    ylab(bquote(~log[2]~ 'FC IP/Ub')) + theme(axis.title.x = element_text(face = "bold", size = 20), axis.text.x = element_text(face = "bold", size = 24)) +
    theme(axis.title.y = element_text(face = "bold", size = 20),axis.text.y = element_text(face = "bold", size = 24)) +
    geom_hline(aes(yintercept = 0), colour = "blue", linewidth = 1, linetype="dashed", alpha=0.4) +
    theme(plot.title = element_text(hjust = 0.5, size = 24, face = "bold")) +  
    theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold")) + theme(legend.position = "bottom", legend.text = element_text(size = 14)) +
    guides(color = guide_legend(override.aes = list(size=10))) + ylim(c(-18, 18)) 
  gm2 = ggMarginal(ma_tRNA, groupColour = T, groupFill = T)
  
  #rrna_trna_mas[[ncont]] =  ggarrange(gm1, gm2, ncol = 2, align = "hv")
  
  ###################################################################################################################################################################### 
  de_tbl_norm_linc =  data.frame('Ave_CPM'=topTable$logCPM, 'log-FC'= topTable$logFC, 'cl'= dt)
  colnames(de_tbl_norm_linc) = c('Ave_CPM', 'logFC', 'cl')
  de_tbl_norm_linc$col = ifelse(de_tbl_norm_linc$cl == 0, "#C0C0C0", ifelse(de_tbl_norm_linc$cl == 1, "#15C1FB", "#C0C0C0"))
  de_tbl_norm_linc$col = ifelse(grepl("lincRNA_As",rownames(de_tbl_norm_linc)), "#DB1A27", de_tbl_norm_linc$col)
  de_tbl_norm_linc$col = factor(de_tbl_norm_linc$col, levels = c("#DB1A27","#15C1FB", "#C0C0C0"))
  de_tbl_norm_linc$ord = ifelse(de_tbl_norm_linc$col == "#C0C0C0", 1, ifelse(de_tbl_norm_linc$col == "#15C1FB", 2,
                                                                             ifelse(de_tbl_norm_linc$col == "#DB1A27", 3, 0)))
  
  de_tbl_norm_linc = de_tbl_norm_linc[order(de_tbl_norm_linc$ord),]
  recol_2 = c("#0FFA00","#002A86","#C0C0C0")
  ma_linc = ggplot(de_tbl_norm_linc, aes(x=Ave_CPM, y=logFC, colour = col)) +
    geom_point() +  xlab(bquote('Average'~log[2]~' CPM') ) +
    scale_color_manual(values = adjustcolor(recol_2, alpha.f = .9) , name=NULL, labels =c("lincRNA_As", "IP enriched", "Unbound")) + theme_light() +
    ylab(bquote(~log[2]~ ' Fold-change IP/Ub')) + theme(axis.title.x = element_text(face = "bold", size = 20), axis.text.x = element_text(face = "bold", size = 24)) +
    theme(axis.title.y = element_text(face = "bold", size = 20),axis.text.y = element_text(face = "bold", size = 24)) +
    geom_hline(aes(yintercept = 0), colour = "blue", linewidth = 1, linetype="dashed", alpha=0.4) +
    theme(plot.title = element_text(hjust = 0.5, size = 24, face = "bold")) +  
    theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold")) + theme(legend.position = "bottom", legend.text = element_text(size = 14)) +
    guides(color = guide_legend(override.aes = list(size=10))) + ylim(c(-18, 18))
  norm_rrna_linc_list[[ncont]]= ma_linc
  ###################################################################################################################################################################### 
  volcanoData = data.frame(topTable$logFC, -log10(topTable$FDR))
  volcanoData$dt = as.factor(dt)
  rownames(volcanoData) = rownames(topTable)
  colnames(volcanoData) = c("logFC", "neglogPval","dt")
  volcanoData$cat =  ifelse(volcanoData$dt == 0, "#C0C0C0", ifelse(volcanoData$dt == 1, "#15C1FB", "#C0C0C0"))
  volcanoData$cat = ifelse(grepl("rRNA_S", rownames(volcanoData)), "#0B5727", volcanoData$cat)
  volcanoData$ord = ifelse(volcanoData$cat == "#C0C0C0", 1, ifelse(volcanoData$cat == "#333333", 2, 3))
  volcanoData$ord = ifelse(grepl("rRNA_S", rownames(volcanoData)), 4, volcanoData$ord)
  volcanoData = volcanoData[order(volcanoData$ord),]
  volcanoData$cond = ifelse(grepl("Adult", cont), "Adult exWAGO IP vs Adult total", ifelse(grepl("EV", cont), "Vesicular exWAGO IP vs Unbound", "Non-vesicular exWAGO IP vs Unbound"))
  volcanoData$cat = factor(volcanoData$cat, levels = c("#15C1FB","#C0C0C0","#0B5727"))
  
  vp = ggplot(volcanoData, aes(logFC,neglogPval, color = cat)) + theme_light() +
    geom_point(size = 2) +
    xlab(expression("log"[2]*"FC IP/Ub")) + 
    ylab(expression("-log"[10]*"FDR")) +
    scale_color_manual(values = adjustcolor(c("#002A86", "#C0C0C0", "#FA0C57"), alpha.f = .6),
                       labels=factor(c("IP enriched", "Unbound", "rRNA_S"),
                                     levels = c("IP enriched","Unbound","rRNA_S")), name=NULL) +
    theme(axis.title.x = element_text(face = "bold", size = 20), 
          axis.text.x = element_text(face = "bold", size = 24)) +
    theme(axis.title.y = element_text(face = "bold", size = 20),
          axis.text.y = element_text(face = "bold", size = 24)) +
    theme(legend.position = "bottom", legend.text = element_text(size = 14)) +
    guides(color = guide_legend(override.aes = list(size=10))) +
    geom_vline(aes(xintercept = 2), linewidth=1.5, linetype = "dashed", color=adjustcolor("blue", alpha.f = .3)) +
    geom_vline(aes(xintercept = -2), linewidth=1.5, linetype = "dashed", color=adjustcolor("blue", alpha.f = .3)) +
    facet_wrap(~cond) + theme( strip.text = element_text(size = 16, face = "bold", colour = "black"))
  
  
  vplots[[ncont]] = vp
  ###################################################################################################################################################################### 
  de_tbl_norm$cond = ifelse(grepl("Adult", cont), "Adult exWAGO IP vs Adult total", ifelse(grepl("EV", cont), "Vesicular exWAGO IP vs Unbound", "Non-vesicular exWAGO IP vs Unbound"))
  ma = ggplot(de_tbl_norm, aes(x=Ave_CPM, y=logFC, colour = col)) +
    geom_point() +  xlab(bquote('Average'~log[2]~' CPM') ) +
    scale_color_manual(values = adjustcolor(recol, alpha.f = .7) , name=NULL, labels =c("28s", "18s", "8s", "IP enriched", "Unbound")) + theme_light() +
    ylab(bquote(~log[2]~ ' Fold-change')) + theme(axis.title.x = element_text(face = "bold", size = 20), axis.text.x = element_text(face = "bold", size = 24)) +
    theme(axis.title.y = element_text(face = "bold", size = 20),axis.text.y = element_text(face = "bold", size = 24)) +
    geom_hline(aes(yintercept = 0), colour = "blue", linewidth = 1, linetype="dashed", alpha=0.4) +
    theme(plot.title = element_text(hjust = 0.5, size = 24, face = "bold")) +  
    theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold")) + theme(legend.position = "bottom", legend.text = element_text(size = 14)) +
    guides(color = guide_legend(override.aes = list(size=10))) + ylim(c(-14, 14)) + 
    facet_wrap(~cond) + theme( strip.text = element_text(size = 16, face = "bold", colour = "black"))
  
  tmm_ma[[ncont]] = p1
  rrna_ma[[ncont]] = ma
  ########################################################################################################################################################################
  cts = c("LTR.*_As", "intergenic", "LINE.*_As", "exons_As")
  df_lst = list()
  for (i in cts) {
    df = de_tbl_norm
    df$col = ifelse(grepl(i, rownames(df)) & df$cl == 1, "#144675", "#C0C0C0")
    df$ord = ifelse(df$col == "#C0C0C0", 1,2)
    df$class = i
    df = df[order(df$ord),]
    df_lst[[i]] = df
  }
  names(df_lst) = NULL
  df = do.call(rbind, df_lst)
  control_norm = ggplot(df, aes(x=Ave_CPM, y=logFC, colour = col)) +
    geom_point(show.legend = F) +  xlab(bquote('Average'~log[2]~' CPM')) + theme_minimal() + facet_wrap(~class, ncol = 2) + theme_light() +
    theme( strip.text = element_text(size = 11, face = "bold", colour = "black")) +
    ylab(bquote(~log[2]~ ' Fold-change')) + theme(axis.title.x = element_text(face = "bold", size = 24), axis.text.x = element_text(face = "bold", size = 24)) +
    theme(axis.title.y = element_text(face = "bold", size = 24),axis.text.y = element_text(face = "bold", size = 24)) +
    geom_hline(aes(yintercept = 0), colour = "blue", linewidth = 1, linetype="dashed", alpha=0.4) +
    theme(plot.title = element_text(hjust = 0.5, size = 20, face = "bold")) + 
    scale_color_manual(values = adjustcolor(c("#144675", "#BCBBBD"), alpha.f = .6)) + ylim(c(-15,18))
  control_norm_regs_lst[[ncont]] = control_norm
  
  #######################################################################################################################################################################
  cts = c("LTR.*_As", "intergenic", "LINE.*_As", "exons_As")
  plots = list()
  for (i in cts) {
    df = de_tbl_norm
    df$col = ifelse(grepl(i, rownames(df)), "#144675", "#C0C0C0")
    df$ord = ifelse(df$col == "#C0C0C0", 1,2)
    df = df[order(df$ord),]
    df$col = factor(df$col, levels= c("#C0C0C0", "#144675"))
    p = ggplot(df, aes(x=Ave_CPM, y=logFC, colour = df$col)) +
      geom_point() +  xlab(bquote('Average'~log[2]~' CPM')) + theme_light() +
      ylab(bquote(~log[2]~ ' Fold-change')) + theme(axis.title.x = element_text(face = "bold", size = 20), axis.text.x = element_text(face = "bold", size = 18)) +
      theme(axis.title.y = element_text(face = "bold", size = 20),axis.text.y = element_text(face = "bold", size = 18)) +
      geom_hline(aes(yintercept = 0), colour = "blue", linewidth = 1, linetype="dashed", alpha=0.4) + ggtitle(i) +
      theme(plot.title = element_text(hjust = 0.5, size = 18, face = "bold")) +
      scale_color_manual(values = c("#C0C0C0", "#144675"), name=NULL, labels =c("Other", i)) + theme(legend.position = "bottom", legend.text = element_text(size = 14)) +
      guides(color = guide_legend(override.aes = list(size=10)))
    p1 = ggMarginal(p, groupColour = T)
    plots[[i]] = p1
  }
  ar_p = ggarrange(plotlist = plots, ncol = 2, nrow = 2)
  contro_norm_ma_with_dens[[ncont]] =  ar_p
  ########################################################################################################################################################################
  
  rRNA_S = topTable[grepl("^rRNA_S", rownames(topTable)),]$logFC
  yRNA_S = topTable[grepl("yRNA_S", rownames(topTable)),]$logFC
  miRNA_S = topTable[grepl("miRNA_S", rownames(topTable)),]$logFC
  tRNA_S = topTable[grepl("^tRNA_S", rownames(topTable)),]$logFC
  lincRNA_As = topTable[grepl("^lincRNA_As", rownames(topTable)),]$logFC
  lincRNA_S = topTable[grepl("^lincRNA_S", rownames(topTable)),]$logFC
  genes_As = topTable[grepl("exons_As|introns_As", rownames(topTable)),]$logFC
  retro_As = topTable[grepl("LINE.*_As|SINE.*_As|LTR.*_As", rownames(topTable)),]$logFC
  trans_As = topTable[grepl("DNA.*_As|RC.*_As", rownames(topTable)),]$logFC
  rRNA_As = topTable[grepl("^rRNA_As", rownames(topTable)),]$logFC
  yRNA_As = topTable[grepl("yRNA_As", rownames(topTable)),]$logFC
  miRNA_As = topTable[grepl("miRNA_As", rownames(topTable)),]$logFC
  tRNA_As = topTable[grepl("^tRNA_As", rownames(topTable)),]$logFC
  genes_S = topTable[grepl("exons_S|introns_S", rownames(topTable)),]$logFC
  retro_S = topTable[grepl("LINE.*_S|SINE.*_S|LTR.*_S", rownames(topTable)),]$logFC
  trans_S = topTable[grepl("DNA.*_S|RC.*_S", rownames(topTable)),]$logFC
  #unk = topTable[grepl("Unk", rownames(topTable)),]$logFC
  
  fc = list("miRNA_S"=miRNA_S, "miRNA_As"= miRNA_As, "lincRNA_S"=lincRNA_S, "lincRNA_As"=lincRNA_As ,"yRNA_S"=yRNA_S, "yRNA_As"=yRNA_As,
            "tRNA_S"=tRNA_S, "tRNA_As"=tRNA_As, "Genes_S"=genes_S, "Genes_As"=genes_As,
            "Retrotransposon_S"=retro_S, "Retrotransposon_As"=retro_As,
            "Transposon_S"=trans_S, "Transposon_As"=trans_As)
  
  df= reshape::melt(fc)
  colnames(df) = c("value", "class")
  df$class = factor(df$class, levels = c("Retrotransposon_S", "Retrotransposon_As",
                                         "Transposon_S",  "Transposon_As", "Genes_S","Genes_As",
                                         "miRNA_S", "miRNA_As", "lincRNA_S", "lincRNA_As", "yRNA_S", "yRNA_As",
                                         "tRNA_S", "tRNA_As"))
  fc = list("miRNA_S"=miRNA_S, "lincRNA_As"=lincRNA_As, "yRNA_S"=yRNA_S, 
            "tRNA_S"=tRNA_S,  "Genes_As"=genes_As,
            "Retrotransposon_As"=retro_As,
            "Transposon_As"=trans_As)
  #######################################################################################################################################################################
  df= reshape::melt(fc)
  colnames(df) = c("value", "class")
  df$class = factor(df$class, levels = c( "Retrotransposon_As",
                                          "Transposon_As", "Genes_As",
                                          "miRNA_S", "lincRNA_As","yRNA_S", 
                                          "tRNA_S"))
  df$cond = ifelse(grepl("Adult", cont), "Adult exWAGO IP vs Adult total", ifelse(grepl("EV", cont), "Vesicular exWAGO IP vs Unbound", "Non-vesicular exWAGO IP vs Unbound"))
  bp_mexp_nout = ggplot(df, aes(x =class, y=value, fill = class)) + 
    geom_boxplot(outlier.shape = NA) + theme_light() +
    scale_fill_manual(values = c(rep("gray", 12)))+
    ylab(bquote(~log[2]~ 'FC IP/Ub')) + xlab(NULL) +
    theme(axis.title.x = element_text( size = 20), axis.text.x = element_text(size = 14, angle = 90, hjust = 1)) +
    theme(axis.title.y = element_text(size = 20),axis.text.y = element_text(size = 18)) +
    theme(legend.position = "bottom") + guides(fill=guide_legend(nrow=2, byrow=TRUE)) +
    geom_hline(aes(yintercept = 0), colour = "blue", linewidth = 2, linetype="dashed", alpha=0.4) +
    theme(legend.position = "none") + ylim(c(-15, 15)) +
    facet_wrap(~cond) + theme( strip.text = element_text(size = 16, face = "bold", colour = "black"))
  
  boxlist[[ncont]] = bp_mexp_nout
  dgeFile = "dge_rds/nxHelBake1.1.exwago_other_lfc1_rrna_dge.Rds"
  dgeFile = sub("other", ncont, dgeFile)
  saveRDS(dge, dgeFile)
  ################################################################################################################################################################################
  categories = c("DNA.*As|RC.*As", "DNA.*S|RC.*S", "LINE.*As|PLE.*As|Retrop*As", "LINE.*S|PLE.*S|Retrop*S", "LTR.*As", "LTR.*S",
                 "Unk", "Low|Simple|Sat|SINE|MITE", "lincRNA_As", "miRNA_S", 
                 "^tRNA_S", "_rRNA|yRNA|snRNA|snoRNA|other_ncRNA|piRNA|lncRNA|miRNA_As|lincRNA_S|^tRNA_As", 
                 "exons_As", "exons_S", "introns_As", "introns_S", "intergenic")
  
  names(categories) = c("DNA_As", "DNA_S", "LINE_As", "LINE_S", "LTR_As", "LTR_S", "Unknown", 
                        "Other_repeat", "lincRNA_As", "miRNA_S", "tRNA_S", 
                        "other_ncRNA", "exons_As", "exons_S", "introns_As", "introns_S", "intergenic")
  allFams = list()
  for (i in names(categories)) {
    tmp.cat = categories[i]
    tmp.sqs = dge$counts[grepl(tmp.cat, rownames(dge$counts)),] %>% rownames()
    allFams[[i]] =tmp.sqs 
  }
  design = dge$design
  allFamsIdx = lapply(allFams, function(x) rownames(dge$counts) %in% x)
  allFamsIdx = allFamsIdx[sapply(allFamsIdx, sum) > 5]
  cameraRes = camera(dge, allFamsIdx, design, contrast)
  cameraRes_lis[[ncont]] = cameraRes
  
  ##########################################################################################################################################################################
  aveCPM = 2 ** aveLogCPM(dge, normalized.lib.sizes = T, prior.count = .01)
  names(aveCPM) = rownames(dge$counts)
  
  aveCPMByfam = split(aveCPM, sub(":.*", "", names(aveCPM)))
  
  categories =  names(aveCPMByfam)[grepl("DNA|LINE|LTR|Unk",names(aveCPMByfam))]
  names(categories) = categories
  
  allFams = list()
  for (i in names(categories)) {
    tmp.cat = categories[i]
    tmp.sqs = dge$counts[grepl(tmp.cat, rownames(dge$counts)),] %>% rownames()
    allFams[[i]] =tmp.sqs 
  }
  
  allFamsIdx = lapply(allFams, function(x) rownames(dge$counts) %in% x)
  allFamsIdx = allFamsIdx[sapply(allFamsIdx, sum) > 5]
  cameraRes = camera(dge, allFamsIdx, design, contrast)
  cameraRes = cameraRes[match(names(categories), rownames(cameraRes)),]
  rownames(cameraRes) = names(categories)
  cameraRes$Feature = ifelse(grepl("DNA|RC", rownames(cameraRes)), "DNA", 
                             ifelse(grepl("LINE|PLE", rownames(cameraRes)), "LINE",
                                    ifelse(grepl("LTR", rownames(cameraRes)), "LTR", "Unknown")))
  cameraRes_fam_lis[[ncont]] = cameraRes
}
#######################################################################################################################################################################
# plot and save
# This could be a mapply but I use this way just if i need to change a particular plot
ma_tmm = ggarrange(tmm_ma$`Adult_IP-Adult_total`, tmm_ma$`EV_IP-EV_Ub`,
                   tmm_ma$`Sup_IP-Sup_Ub`,
                   ncol = 3, common.legend = T, legend = "bottom")



#ma_thesis_fig2 = ggarrange(tmm_ma$`Adult_IP-Adult_total`, norm_rrna_simp_list$`Adult_IP-Adult_total`, labels = c('A)', 'B)'), font.label = list(size=24))
#ggsave(filename = "nxHelBake1.1.ma_tmm_vs_rRNA_adult.png", plot = ma_thesis_fig2, device = "png", 
#       path = "figures/", 
#       dpi = 300, width = 12, height = 5)



ggsave(filename = "nxHelBake1.1.ma_tmm.png", plot = ma_tmm, device = "png", 
       path = "figures/", 
       dpi = 300, width = 18, height = 5)

ma_norm = ggarrange(rrna_ma$`Adult_IP-Adult_total`, rrna_ma$`EV_IP-EV_Ub`,
                    rrna_ma$`Sup_IP-Sup_Ub`, 
                    ncol = 3, common.legend = T, legend = "bottom")

ggsave(filename = "nxHelBake1.1.ma_rrna.png", plot = ma_norm, device = "png", 
       path = "figures/", 
       dpi = 300, width = 16, height = 5)


ma_norm_simp = ggarrange(norm_rrna_simp_list$`Adult_IP-Adult_total`, 
                         norm_rrna_simp_list$`EV_IP-EV_Ub`,
                         norm_rrna_simp_list$`Sup_IP-Sup_Ub`,
                         ncol = 3, common.legend = T, legend = "bottom")

ggsave(filename = "nxHelBake1.1.ma_rrna_simp.png", plot = ma_norm_simp, device = "png", 
       path = "figures/", 
       dpi = 300, width = 16, height = 5)


ma_norm_linc = ggarrange(norm_rrna_linc_list$`Adult_IP-Adult_total`, 
                         norm_rrna_linc_list$`EV_IP-EV_Ub`,
                         norm_rrna_linc_list$`Sup_IP-Sup_Ub`,
                         ncol = 3, common.legend = T, legend = "bottom")

ggsave(filename = "nxHelBake1.1.ma_rrna_linc.png", plot = ma_norm_linc, device = "png", 
       path = "figures/", 
       dpi = 300, width = 16, height = 5)


v_plot = ggarrange(vplots$`Adult_IP-Adult_total`, vplots$`EV_IP-EV_Ub`,
                   vplots$`Sup_IP-Sup_Ub`,
                   ncol = 3, common.legend = T, legend = "bottom")

ggsave(filename = "nxHelBake1.1.vp_rrna.png", plot = v_plot, device = "png", 
       path = "figures/", 
       dpi = 300, width = 16, height = 5)

b_plots = ggarrange(boxlist$`Adult_IP-Adult_total`, boxlist$`EV_IP-EV_Ub`,
                    boxlist$`Sup_IP-Sup_Ub`,
                    ncol = 3)
ggsave(filename = "nxHelBake1.1.bp_rrna.png", plot = b_plots, device = "png", 
       path = "figures/", 
       dpi = 300, width = 16, height = 5)


rr_nonorm = ggarrange(nonorm_rrna_ma$`Adult_IP-Adult_total`, 
                      nonorm_rrna_ma$`EV_IP-EV_Ub`, 
                      nonorm_rrna_ma$`Sup_IP-Sup_Ub`,
                      ncol = 3, common.legend = T, legend = "bottom")

ggsave(filename = "nxHelBake1.1.ma_rrna_nonorm.png", plot = rr_nonorm, device = "png", 
       path = "figures/", 
       dpi = 300, width = 16, height = 5)

cont_regs = ggarrange(control_regs_lst$`Adult_IP-Adult_total`, 
                      control_regs_lst$`EV_IP-EV_Ub`, 
                      control_regs_lst$`Sup_IP-Sup_Ub`,
                      nrow = 3)

ggsave(filename = "nxHelBake1.1.mas_cont-regs_nonorm.png", plot = cont_regs, device = "png", 
       path = "figures/", 
       dpi = 300, width = 14, height = 26)


control_adult_enr =  contro_norm_ma_with_dens$`Adult_IP-Adult_total`
ggsave(filename = "nxHelBake1.1.adult_cont-regs_nonorm_dens.png", plot = control_adult_enr, device = "png", 
       path = "figures/", 
       dpi = 300, width = 7, height = 8)

ma_cont_rrna_trna =  ggarrange(rrna_trna_mas$`Adult_IP-Adult_total`, rrna_trna_mas$`EV_IP-EV_Ub`, rrna_trna_mas$`Sup_IP-Sup_Ub`, nrow=3)
ggsave(filename = "nxHelBake1.1.mas_rrna_trna.png", plot = ma_cont_rrna_trna, device = "png", 
       path = "figures/", 
       dpi = 300, width = 14, height = 16)

########################################################################################################
#path_vp2 = "/Users/isaacmartinezugalde/storage/Data/evo_exwago/analyses/hBake_exwago_guides/vaccine_paper_v2"

#ma_tmm = ggarrange(tmm_ma$`EV_IP-EV_Ub`,
#                   tmm_ma$`Sup_IP-Sup_Ub`,
#                   ncol = 2, common.legend = T, legend = "bottom")

#ggsave(filename = "nxHelBake1.1.ma_tmm.png", plot = ma_tmm, device = "png", 
#       path = path_vp2, 
#       dpi = 300, width = 11, height = 5)


#ma_norm_simp = ggarrange(norm_rrna_simp_list$`EV_IP-EV_Ub`,
#                         norm_rrna_simp_list$`Sup_IP-Sup_Ub`,
#                         ncol = 2, common.legend = T, legend = "bottom")

#ggsave(filename = "nxHelBake1.1.ma_rrna_simp.png", plot = ma_norm_simp, device = "png", 
#       path = path_vp2, 
#       dpi = 300, width = 11, height = 5)

#b_plots = ggarrange(boxlist$`EV_IP-EV_Ub`,
#                    boxlist$`Sup_IP-Sup_Ub`,
#                    ncol = 2)
#ggsave(filename = "nxHelBake1.1.bp_rrna.png", plot = b_plots, device = "png", 
#       path = path_vp2, 
#       dpi = 300, width = 11, height = 5)

#v_plot = ggarrange(vplots$`EV_IP-EV_Ub`,
#                   vplots$`Sup_IP-Sup_Ub`,
#                   ncol = 2, common.legend = T, legend = "bottom")

#ggsave(filename = "nxHelBake1.1.vp_rrna.png", plot = v_plot, device = "png", 
#       path = path_vp2, 
#       dpi = 300, width = 11, height = 5)

saveRDS(cameraRes_lis,"nxHelBake1.1_exWAGO_all_IPs_CAMERA.Rds")
saveRDS(cameraRes_fam_lis, "nxHelBake1.1_exWAGO_all_IPs_CAMERA_fams.Rds")
