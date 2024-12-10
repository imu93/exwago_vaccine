setwd("~/storage/Data/vaccine_exwago/analyses/ss3/dea/")
pacman::p_load(edgeR, rtracklayer, ggplot2, ggpubr, kableExtra, 
               gplots, RColorBrewer, grDevices, dplyr, purrr, ggrepel, ggExtra)
counts = read.table("nxHelBake1.1.counts.1827.txt.gz")
colnames(counts) = sub("\\.trim.*", "", colnames(counts))
colnames(counts)
counts = counts[,grepl("EV_IP_polyP|Sup_IP_polyP", colnames(counts)),]
###################################################################################################################################################################################################################

# Now let's check the number of counts to have a grasp of the analysis
lib.size = colSums(counts)
round((lib.size)/1e6, 2)
# I'll use a MDS plot to have an idea of the influence of the batch effect
# But first I need to build my dge object
group = colnames(counts) %>% strsplit("_") %>%
  map(function(x){paste(x[c(1,2)], collapse = "_")}) %>%  unlist() %>% 
  factor()

dge = DGEList(counts=counts, group=group, lib.size = lib.size)
keep = filterByExpr(dge, min.count= 5)
dge = dge[keep, , keep.lib.sizes=FALSE]
print(table(keep))
f.exp.counts = colSums(dge$counts)


batch = factor(paste0("batch",sub(".*(r.*)","\\1",rownames(dge$samples))))
design = model.matrix(~0+dge$samples$group+batch)
colnames(design) =  c(levels(dge$samples$group),levels(batch)[-1])
###################################################################################################################################################################################################################
dge = calcNormFactors(dge, method = "TMM")
# As a techincal note: use alpha in rich colors  to have transparency
colors = rich.colors(length(levels(group)), alpha = .5)
names(colors) = levels(group)
par(mfrow=c(1,2))
plotMDS(dge, col=colors[dge$samples$group], pch = 16, cex=3,main="MDS plot of H. bakeri exWAGO EV vs Sup")
plot.new()
legend("topright", legend = names(colors), pch = 16, col = unique(colors))
###################################################################################################################################################################################################################

# Let's estimate de dispersion
dge = estimateDisp(dge, design=design, robust=TRUE)
plotBCV(dge)
fit = glmFit(dge, design, dispersion = dge$common.dispersion)
contrast = makeContrasts("EV_IP-Sup_IP",levels = dge$design)
gt = glmTreat(fit, contrast=contrast, lfc = 1)
topTags(gt)
dt = decideTestsDGE(gt, adjust.method="BH", p.value=0.05, lfc = 1)
table(dt)
topTable = topTags(gt, n=Inf)$table
topTable = topTable[rownames(dt),]
de_tbl = data.frame('Ave_CPM'=topTable$logCPM, 'log-FC'= topTable$logFC, 'cl'= dt)
colnames(de_tbl) = c('Ave_CPM', 'logFC', 'cl')
de_tbl$col = ifelse(de_tbl$cl == 0, "#C0C0C0", ifelse(de_tbl$cl == 1, "#5E2129", "#333333"))
de_tbl$ord = ifelse(de_tbl$col == "#C0C0C0", 1, ifelse(de_tbl$col == "#333333", 2, 3))
de_tbl = de_tbl[order(de_tbl$ord),]
#write.table(de_tbl, "~/storage/Data/evo_exwago/analyses/hBake_exwago_guides/ss3/thesis_prop_tables/Hb_exWAGO_EV-Sup_IP_deTable.txt", quote = F, sep = "\t")
#######################################################################################################################################################################
volcanoData = data.frame(topTable$logFC, -log10(topTable$FDR))
volcanoData$dt = as.factor(dt)
rownames(volcanoData) = rownames(topTable)
colnames(volcanoData) = c("logFC", "negLogPval","dt")
volcanoData$exp = "Vesicular vs Non-vesicular exWAGO IPs"
p1= ggplot(volcanoData, aes(logFC,negLogPval, color = dt)) + theme_light() +
  geom_point(size = 2) +
  xlab(expression("log"[2]*"FC IP/IP")) + 
  ylab(expression("-log"[10]*"FDR")) +
  scale_color_manual(values = adjustcolor(c("#E8252C","#BCBBBD", "#25A0E8"), alpha.f = .6),
                     labels=factor(c("EV-depleted HES", "non-DE","EVs"),
                                   levels = c("EV-depleted HES", "non-DE","EVs")), name=NULL) +
  theme(axis.title.x = element_text(face = "bold", size = 20), 
        axis.text.x = element_text( size = 24)) +
  theme(axis.title.y = element_text(face = "bold", size = 20),
        axis.text.y = element_text(size = 24)) +
  theme(legend.position = "bottom", legend.text = element_text(size = 14)) +
  guides(color = guide_legend(override.aes = list(size=14))) +
  geom_vline(aes(xintercept = 2), linewidth=1.5, linetype = "dashed", color=adjustcolor("#1B2CC2", alpha.f = .3)) +
  geom_vline(aes(xintercept = -2), linewidth=1.5, linetype = "dashed", color=adjustcolor("#1B2CC2", alpha.f = .3)) +
  theme(strip.text = element_text(size = 14, face = "bold", colour = "black"))  +
  theme(plot.title = element_text(hjust = 0.5, size = 24, face = "bold")) 



ggsave(filename = "nxHelBake1.ewago_ip_ev_vs_sup_vplot.png", plot = p1, device = "png", path = "~/storage/Data/vaccine_exwago/analyses/figures/", 
       width = 6, height = 5, dpi = 300)
#######################################################################################################################################################################
de_tbl$col = factor(de_tbl$col, levels = c("#5E2129","#C0C0C0", "#333333"))

p = ggplot(de_tbl, aes(x=Ave_CPM, y=logFC, colour = col)) +
  geom_point() +  xlab(bquote(~log[2]~' CPM')) + theme_light() +
  ylab(bquote(~log[2]~ 'FC IP/IP')) + theme(axis.title.x = element_text(face = "bold", size = 24), axis.text.x = element_text(face = "bold", size = 24)) +
  theme(axis.title.y = element_text(face = "bold", size = 24),axis.text.y = element_text(face = "bold", size = 24)) +
  geom_hline(aes(yintercept = 0), colour = "#1B2CC2", linewidth = 2, linetype="dashed", alpha=0.4) + theme( strip.text = element_text(size = 14, face = "bold", colour = "black"))  +
  theme(plot.title = element_text(hjust = 0.5, size = 24, face = "bold")) + scale_color_manual(values = adjustcolor(c("#E8252C","#BCBBBD", "#25A0E8"), alpha.f = .6),labels=c("EVs", "non-DE", "EV-depleted HES"), name=NULL)

p = p + theme(legend.position = "bottom", legend.text = element_text(size = 12)) +
  guides(color = guide_legend(override.aes = list(size=12)))

#######################################################################################################################################################################
evs = dge$counts[,grepl("EV_IP", colnames(dge$counts))]
# But only those DE
des = evs[match(rownames(de_tbl[de_tbl$cl == 1,]), rownames(evs)),]
# I'll use the mean
des_mean = rowMeans(des)
names(des_mean) = rownames(des)
sum_fam = split(des_mean, gsub(":.*", "", names(des_mean)))
options(scipen = 999)
fams = unlist(lapply(sum_fam, function(x){sum(x)*100/sum(des_mean)}))
EV_fams = fams 
EV_fams = EV_fams[!grepl("exons|introns|miRNA|other|inter|^tRNA|^rRNA|miRNA|piRNA|yRNA|snRNA|snoRNA|other_ncRNA|Low|Sat|Simple|inter|linc", names(EV_fams))]
EV_top = EV_fams %>% sort(decreasing = T) %>% head()


sup = dge$counts[,grepl("Sup_IP", colnames(dge$counts))]
# But only those DE
des = evs[match(rownames(de_tbl[de_tbl$cl == -1,]), rownames(sup)),]
# I'll use the mean
des_mean = rowMeans(des)
names(des_mean) = rownames(des)
sum_fam = split(des_mean, gsub(":.*", "", names(des_mean)))
options(scipen = 999)
fams = unlist(lapply(sum_fam, function(x){sum(x)*100/sum(des_mean)}))
Sup_fams = fams 
Sup_fams = Sup_fams[!grepl("exons|introns|miRNA|other|inter|^tRNA|^rRNA|miRNA|piRNA|yRNA|snRNA|snoRNA|other_ncRNA|Low|Sat|Simple|inter|linc", names(Sup_fams))]
Sup_top = Sup_fams %>% sort(decreasing = T) %>% head()


s_fams = c(EV_top,Sup_top) %>% names() %>% unique()

sub_ev = EV_fams[s_fams]
names(sub_ev) = s_fams
sub_ev = ifelse(is.na(sub_ev) == T, 0, sub_ev)
sub_sup = Sup_fams[s_fams]


lst2plt = list(sub_ev, sub_sup)
names(lst2plt) = c("EVs", "EV-depleted HES")
color = c(rich.colors(11))
f_lst = list()
for (i in names(lst2plt)) {
  x = as.data.frame(lst2plt[[i]])
  rownames(x) = s_fams
  x$class = rownames(x)
  x$Treatment = i
  x$col = color
  colnames(x) = c("Library_prop", "Class","Treatment", "col")
  f_lst[[i]] = x
}

names(f_lst) = NULL

df  = do.call(rbind, f_lst)
###################################################################################################################################################################################################################
df$Treatment = factor(df$Treatment, levels = c("EV-depleted HES", "EVs"))
df$col = factor(df$col, levels = color)
df$exp = "exWAGO TE family enrichment comparison"
p2 = ggplot(df, aes(x=Treatment , y=Library_prop, group=Class, fill=Class)) +
  geom_line(linetype = "dashed") + geom_point(size=8, shape=21) + 
  scale_fill_manual(values = adjustcolor(rev(color), alpha.f = .7), name="Family") + 
  theme_light() + 
  theme(strip.text = element_text(size = 14, face = "bold", colour = "black"))  +
  theme(plot.title = element_text(hjust = 0.5, size = 24, face = "bold")) +
  ylab("% of sRNAs") + xlab("") +
  theme(axis.title.x = element_text(face = "bold", size = 20), axis.text.x = element_text(size = 18)) +
  theme(axis.title.y = element_text(face = "bold", size = 20),axis.text.y = element_text(size = 18)) +
  theme(plot.title = element_text(hjust = 0, size = 24, face = "bold")) + theme(legend.position = "bottom") + 
  guides(color = guide_legend(override.aes = list(size=8)), fill=guide_legend(ncol=1)) + theme( legend.text = element_text(size = 14), legend.title = element_text(size = 14))

ggsave(filename = "nxHelBake1.ewago_ip_ev_vs_sup_topFams.pdf", plot = p2, device = "pdf", path = "~/storage/Data/vaccine_exwago/analyses/figures/", 
       width = 6, height = 10, dpi = 300)
