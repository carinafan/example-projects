# this script uses the same workspace that is created by the main analysis script

# function to create shared legend across multiple ggplots
# from http://rpubs.com/sjackman/grid_arrange_shared_legend

grid_arrange_shared_legend <- function(...) {
  plots <- list(...)
  g <- ggplotGrob(plots[[1]] + theme(legend.position="bottom"))$grobs
  legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
  lheight <- sum(legend$height)
  grid.arrange(
    do.call(arrangeGrob, lapply(plots, function(x)
      x + theme(legend.position="none"))),
    legend,
    ncol = 1,
    heights = grid::unit.c(unit(1, "npc") - lheight, lheight))
}

#---- figure 1: scree plot and dim 1 vs 2 for PCA on SAM ----

# scree plot

tiff("figures/sam_scree.tiff", width = 6, height = 4, units = "in", res = 600)

PTCA4CATA::PlotScree(
  ev = pca.sam$ExPosition.Data$pdq$tau,
  max.ev = max(pca.sam$ExPosition.Data$eig)
) +
  coord_fixed(ratio = .5)

dev.off()

# set up for PCA plots

Fj.sam = pca.sam$ExPosition.Data$fj

map_labels = createxyLabels.gen(
  x_axis = 1, y_axis = 2,
  lambda = pca.sam$ExPosition.Data$eigs,
  tau = pca.sam$ExPosition.Data$t)

jmap.sam = PTCA4CATA::createFactorMap(
  Fj.sam,
  col.points = colour.sam,
  col.labels = colour.sam,
  alpha.points. = .3,
  text.cex = 3)

jaggmap.sam =
  jmap.sam$zeMap_background +
  jmap.sam$zeMap_dots +
  jmap.sam$zeMap_text +
  coord_fixed(ratio = 1.2) +
  map_labels

Fi.sam = pca.sam$ExPosition.Data$fi

imap.sam = PTCA4CATA::createFactorMap(
  Fi.sam,
  display.labels = FALSE,
  alpha.points = .3,
  cex = .3,
  text.cex = 3)

iaggmap.sam =
  imap.sam$zeMap_background +
  imap.sam$zeMap_dots +
  coord_equal(ratio = .8)

# PCA plot of variables

tiff("figures/sam_pca_dim1v2.tiff", width = 3, height = 4, units = "in", res = 600)

print(jaggmap.sam)

dev.off()

# PCA plot of subjects 

tiff("figures/sam_pca_dim1v2_subjects.tiff", width = 3, height = 2, units = "in", res = 600)

print(iaggmap.sam)

dev.off()

#---- figure 2: univariate relationships between SAM and OSIQ ----

scatter_samE_osiq =   
  df.sam.osiq %>%
  select(sam_epi_raw, osiq_mini_o, osiq_mini_s) %>%
  melt(id = "sam_epi_raw", measured = c("osiq_mini_o", "osiq_mini_s")) %>% 
  ggplot(aes(x = value, y = sam_epi_raw, colour = variable)) +
  geom_point(alpha = .1, position = "jitter") +
  geom_smooth(method = "lm", size = 2) +
  theme_bw() +
  labs(x = "OSIQ score", y = "SAM-episodic score") +
  scale_colour_manual(name = "OSIQ domain", 
                      values = colourCode.osiq,
                      labels = c("OSIQ-object", "OSIQ-spatial"))

scatter_samS_osiq =   
  df.sam.osiq %>%
  select(sam_spa_raw, osiq_mini_o, osiq_mini_s) %>%
  melt(id = "sam_spa_raw", measured = c("osiq_mini_o", "osiq_mini_s")) %>% 
  ggplot(aes(x = value, y = sam_spa_raw, colour = variable)) +
  geom_point(alpha = .1, position = "jitter") +
  geom_smooth(method = "lm", size = 2) +
  theme_bw() +
  labs(x = "OSIQ score", y = "SAM-spatial score") +
  scale_colour_manual(name = "OSIQ domain", 
                      values = colourCode.osiq,
                      labels = c("OSIQ-object", "OSIQ-spatial"))

scatter_samF_osiq =   
  df.sam.osiq %>%
  select(sam_fut_raw, osiq_mini_o, osiq_mini_s) %>%
  melt(id = "sam_fut_raw", measured = c("osiq_mini_o", "osiq_mini_s")) %>% 
  ggplot(aes(x = value, y = sam_fut_raw, colour = variable)) +
  geom_point(alpha = .1, position = "jitter") +
  geom_smooth(method = "lm", size = 2) +
  theme_bw() +
  labs(x = "OSIQ score", y = "SAM-future score") +
  scale_colour_manual(name = "OSIQ domain",
                      values = colourCode.osiq,
                      labels = c("OSIQ-object", "OSIQ-spatial"))

scatter_samM_osiq =
  df.sam.osiq %>%
  select(sam_sem_raw, osiq_mini_o, osiq_mini_s) %>%
  melt(id = "sam_sem_raw", measured = c("osiq_mini_o", "osiq_mini_s")) %>% 
  ggplot(aes(x = value, y = sam_sem_raw, colour = variable)) +
  geom_point(alpha = .1, position = "jitter") +
  geom_smooth(method = "lm", size = 2) +
  theme_bw() +
  labs(x = "OSIQ score", y = "SAM-semantic score") +
  scale_colour_manual(name = "OSIQ domain",
                      values = colourCode.osiq,
                      labels = c("OSIQ-object", "OSIQ-spatial"))

png("figures/sam_osiq_scatter.png", width = 10, height = 10, units = "in", res = 300)

grid_arrange_shared_legend(scatter_samE_osiq, scatter_samS_osiq, 
                           scatter_samF_osiq, scatter_samM_osiq)

dev.off()

#---- figure 3: PLSC on SAM and OSIQ ----

# scree plot

png("figures/sam_osiq_scree.png", width = 5, height = 5, units = "in", res = 600)

PTCA4CATA::PlotScree(
  ev = plsc.sam.osiq$TExPosition.Data$eigs,
  max.ev = max(plsc.sam.osiq$TExPosition.Data$eigs)
)

dev.off()

# variable plot

loadings = rbind(plsc.sam.osiq$TExPosition.Data$fi,
                 plsc.sam.osiq$TExPosition.Data$fj)

colnames(loadings) = paste0("Dimension ", 1:ncol(loadings))

colour.plsc = c(colour.osiq, colour.sam)

map.loadings = PTCA4CATA::createFactorMap(
  loadings,
  col.points = colour.plsc,
  alpha.points = .2,
  cex = 2,
  col.labels = colour.plsc)

map_labels = createxyLabels.gen(
  x_axis = 1, y_axis = 2,
  lambda = plsc.sam.osiq$TExPosition.Data$eigs,
  tau = plsc.sam.osiq$TExPosition.Data$t)

aggmap.loadings =
  map.loadings$zeMap_background +
  map.loadings$zeMap_text +
  map.loadings$zeMap_dots +
  map_labels

png("figures/sam_osiq_plsc.png", width = 6, height = 4, units = "in", res = 300)

print(aggmap.loadings)

dev.off()

#---- figure 4: contributions to PLSC LVs ----

# LV1 OSIQ

Fi = plsc.sam.osiq$TExPosition.Data$fi

contrib.i = plsc.sam.osiq$TExPosition.Data$ci

contrib.i.signed = contrib.i* sign(Fi)

png("figures/plsc_osiq_lv1.png", width = 5, height = 3, units = "in", res = 600)

PrettyBarPlot2(
  bootratio = round(100 * contrib.i.signed[,1]),
  threshold = 100/ nrow(contrib.i.signed),
  ylim = NULL,
  color4bar = gplots::col2hex(colour.osiq),
  color4ns = "gray75",
  plotnames = TRUE,
  main = "Important Contributions to LV1 from OSIQ items",
  ylab = "Signed Contributions")

dev.off()

# LV1 SAM

Fj = plsc.sam.osiq$TExPosition.Data$fj

contrib.j = plsc.sam.osiq$TExPosition.Data$cj

contrib.j.signed = contrib.j * sign(Fj)

png("figures/plsc_sam_lv1.png", width = 5, height = 3, units = "in", res = 600)

PrettyBarPlot2(
  bootratio = round(100 * contrib.j.signed[,1]),
  threshold = 100/ nrow(contrib.j.signed),
  ylim = NULL,
  color4bar = gplots::col2hex(colour.sam),
  color4ns = "gray75",
  plotnames = TRUE,
  main = "Important Contributions to LV1 from SAM items",
  ylab = "Signed Contributions")

dev.off()

# LV2 OSIQ

png("figures/plsc_osiq_lv2.png", width = 5, height = 3, units = "in", res = 600)

PrettyBarPlot2(
  bootratio = round(100 * contrib.i.signed[,2]),
  threshold = 100/ nrow(contrib.i.signed),
  ylim = NULL,
  color4bar = gplots::col2hex(colour.osiq),
  color4ns = "gray75",
  plotnames = TRUE,
  main = "Important Contributions to LV2 from OSIQ items",
  ylab = "Signed Contributions")

dev.off()

# LV2 SAM

png("figures/plsc_sam_lv2.png", width = 5, height = 3, units = "in", res = 600)

PrettyBarPlot2(
  bootratio = round(100 * contrib.j.signed[,2]),
  threshold = 100/ nrow(contrib.j.signed),
  ylim = NULL,
  color4bar = gplots::col2hex(colour.sam),
  color4ns = "gray75",
  plotnames = TRUE,
  main = "Important Contributions to LV2 from SAM items",
  ylab = "Signed Contributions")

dev.off()

#---- supplemental figure 1: SAM PCA dim 3 through 5 ----

# set up PCA plots

jmap.sam.3v4 = PTCA4CATA::createFactorMap(
  Fj.sam,
  axis1 = 3, axis2 = 4,
  col.points = colour.sam,
  col.labels = colour.sam,
  alpha.points. = .3)

jaggmap.sam.3v4 =
  jmap.sam.3v4$zeMap_background +
  jmap.sam.3v4$zeMap_dots +
  jmap.sam.3v4$zeMap_text

imap.sam.3v4 = PTCA4CATA::createFactorMap(
  Fi.sam,
  axis1 = 3, axis2 = 4,
  cex = 1,
  display.labels = FALSE,
  alpha.points = .1)

iaggmap.sam.3v4 =
  imap.sam.3v4$zeMap_background +
  imap.sam.3v4$zeMap_dots

jmap.sam.4v5 = PTCA4CATA::createFactorMap(
  Fj.sam,
  axis1 = 4, axis2 = 5,
  col.points = colour.sam,
  col.labels = colour.sam,
  alpha.points. = .1)

jaggmap.sam.4v5 =
  jmap.sam.4v5$zeMap_background +
  jmap.sam.4v5$zeMap_dots +
  jmap.sam.4v5$zeMap_text

imap.sam.4v5 = PTCA4CATA::createFactorMap(
  Fi.sam,
  axis1 = 4, axis2 = 5,
  cex = 1,
  display.labels = FALSE,
  alpha.points = .1)

iaggmap.sam.4v5 =
  imap.sam.4v5$zeMap_background +
  imap.sam.4v5$zeMap_dots

# PCA plot of variables (3 vs 4)

png("figures/sam_pca_dim3v4.png", width = 4, height = 4, units = "in", res = 300)

print(jaggmap.sam.3v4)

dev.off()

# PCA plot of subjects (3 vs 4)

png("figures/sam_pca_dim3v4_subjects.png", width = 4, height = 4, units = "in", res = 300)

print(iaggmap.sam.3v4)

dev.off()

# PCA plot of variables (4 vs 5)

png("supplemental/figures/sam_pca_dim4v5.png", 
    width = 4, height = 4, units = "in", res = 300)

print(jaggmap.sam.4v5)

dev.off()

# PCA plot of subjects (4 vs 5)

png("figures/sam_pca_dim4v5_subjects.png", width = 4, height = 4, units = "in", res = 300)

print(iaggmap.sam.4v5)

dev.off()

#---- supplemental figure 2: SAM PCA on independent sample (2013 dataset) ----

# scree plot

png("supplemental/figures/scree_samOG.png", 
    width = 6, height = 4, units = "in", res = 600)

PTCA4CATA::PlotScree(
  ev = pca.sam.og$ExPosition.Data$pdq$tau,
  max.ev = max(pca.sam.og$ExPosition.Data$eig)
) +
  coord_fixed(ratio = .5)

dev.off()

# set up PCA plots

Fj.sam.og = pca.sam.og$ExPosition.Data$fj

jmap.sam.og.1v2 = PTCA4CATA::createFactorMap(
  Fj.sam.og,
  axis1 = 1, axis2 = 2,
  col.points = colour.sam,
  col.labels = colour.sam,
  alpha.points. = .3)

jaggmap.sam.og.1v2 =
  jmap.sam.og.1v2$zeMap_background +
  jmap.sam.og.1v2$zeMap_dots +
  jmap.sam.og.1v2$zeMap_text

Fi.sam.og = pca.sam.og$ExPosition.Data$fi

imap.sam.og.1v2 = PTCA4CATA::createFactorMap(
  Fi.sam.og,
  axis1 = 1, axis2 = 2,
  cex = 1,
  display.labels = FALSE,
  alpha.points = .1)

iaggmap.sam.og.1v2 =
  imap.sam.og.1v2$zeMap_background +
  imap.sam.og.1v2$zeMap_dots

jmap.sam.og.3v4 = PTCA4CATA::createFactorMap(
  Fj.sam.og,
  axis1 = 3, axis2 = 4,
  col.points = colour.sam,
  col.labels = colour.sam,
  alpha.points. = .1)

jaggmap.sam.og.3v4 =
  jmap.sam.og.3v4$zeMap_background +
  jmap.sam.og.3v4$zeMap_dots +
  jmap.sam.og.3v4$zeMap_text

imap.sam.og.3v4 = PTCA4CATA::createFactorMap(
  Fi.sam.og,
  axis1 = 3, axis2 = 4,
  cex = 1,
  display.labels = FALSE,
  alpha.points = .1)

iaggmap.sam.og.3v4 =
  imap.sam.og.3v4$zeMap_background +
  imap.sam.og.3v4$zeMap_dots

# PCA plot of variables (1 vs 2)

png("supplemental/figures/pca_samOG_1v2.png", 
    width = 4, height = 4, units = "in", res = 300)

print(jaggmap.sam.og.1v2)

dev.off()

# PCA plot of subjects (1 vs 2)

png("supplemental/figures/pca_samOG_1v2_subjects.png", 
    width = 4, height = 4, units = "in", res = 300)

print(iaggmap.sam.og.1v2)

dev.off()

# PCA plot of variables (3 vs 4)

png("supplemental/figures/pca_samOG_3v4.png", 
    width = 4, height = 4, units = "in", res = 300)

print(jaggmap.sam.og.3v4)

dev.off()

# PCA plot of subjects (3 vs 4)

png("supplemental/figures/pca_samOG_3v4_subjects.png", 
    width = 4, height = 4, units = "in", res = 300)

print(iaggmap.sam.og.3v4)

dev.off()

#---- supplemental figure 3: SAM PCA on independent sample (2019 dataset) ----

# scree plot

png("supplemental/figures/scree_samET.png", 
    width = 6, height = 4, units = "in", res = 600)

PTCA4CATA::PlotScree(
  ev = pca.sam.et$ExPosition.Data$pdq$tau,
  max.ev = max(pca.sam.et$ExPosition.Data$eig)
) +
  coord_fixed(ratio = .5)

dev.off()

# set up PCA plots

Fj.sam.et = pca.sam.et$ExPosition.Data$fj

jmap.sam.et.1v2 = PTCA4CATA::createFactorMap(
  Fj.sam.et,
  axis1 = 1, axis2 = 2,
  col.points = colour.sam,
  col.labels = colour.sam,
  alpha.points. = .3)

jaggmap.sam.et.1v2 =
  jmap.sam.et.1v2$zeMap_background +
  jmap.sam.et.1v2$zeMap_dots +
  jmap.sam.et.1v2$zeMap_text

Fi.sam.et = pca.sam.et$ExPosition.Data$fi

imap.sam.et.1v2 = PTCA4CATA::createFactorMap(
  Fi.sam.et,
  axis1 = 1, axis2 = 2,
  cex = 1,
  display.labels = FALSE,
  alpha.points = .1)

iaggmap.sam.et.1v2 =
  imap.sam.et.1v2$zeMap_background +
  imap.sam.et.1v2$zeMap_dots

jmap.sam.et.3v4 = PTCA4CATA::createFactorMap(
  Fj.sam.et,
  axis1 = 3, axis2 = 4,
  col.points = colour.sam,
  col.labels = colour.sam,
  alpha.points. = .1)

jaggmap.sam.et.3v4 =
  jmap.sam.et.3v4$zeMap_background +
  jmap.sam.et.3v4$zeMap_dots +
  jmap.sam.et.3v4$zeMap_text

imap.sam.et.3v4 = PTCA4CATA::createFactorMap(
  Fi.sam.et,
  axis1 = 3, axis2 = 4,
  cex = 1,
  display.labels = FALSE,
  alpha.points = .1)

iaggmap.sam.et.3v4 =
  imap.sam.et.3v4$zeMap_background +
  imap.sam.et.3v4$zeMap_dots

# PCA plot of variables (1 vs 2)

png("supplemental/figures/pca_samET_1v2.png", 
    width = 4, height = 4, units = "in", res = 300)

print(jaggmap.sam.et.1v2)

dev.off()

# PCA plot of subjects (1 vs 2)

png("supplemental/figures/pca_samET_1v2_subjects.png", 
    width = 4, height = 4, units = "in", res = 300)

print(iaggmap.sam.et.1v2)

dev.off()

# PCA plot of variables (3 vs 4)

png("supplemental/figures/pca_samET_3v4.png", 
    width = 4, height = 4, units = "in", res = 300)

print(jaggmap.sam.et.3v4)

dev.off()

# PCA plot of subjects (3 vs 4)

png("supplemental/figures/pca_samET_3v4_subjects.png", 
    width = 4, height = 4, units = "in", res = 300)

print(iaggmap.sam.et.3v4)

dev.off()

#---- supplemental figure 4: projection of OSIQ onto SAM ----

osiq.sup = supplementaryCols(
  SUP.DATA = expo.scale(df.osiq, scale = FALSE),
  res = pca.sam,
  center = FALSE,
  scale = FALSE)

Fj.sam.osiq = rbind(Fj.sam, osiq.sup$fjj)

colour.sam.osiq = c(colour.sam, colour.osiq)

constraints.sup = lapply(minmaxHelper(Fj.sam, Fj.osiq.sup), "*", 1.1)

jmap.sup = PTCA4CATA::createFactorMap(
  Fj.sam.osiq,
  col.points = colour.sam.osiq,
  col.labels = colour.sam.osiq,
  alpha.points. = .3,
  cex = 1,
  constraints = constraints.sup)

jaggmap.sup = 
  jmap.sup$zeMap_background +
  jmap.sup$zeMap_dots +
  jmap.sup$zeMap_text

png("supplemental/figures/osiq_projection.png", width = 5, height = 6, units = "in", res = 300)

print(jaggmap.sup)

dev.off()


#---- supplemental figure 5: correlation heat map of SAM and OSIQ ----


cor.sam.osiq = cor(df.sam, df.osiq)

png("supplemental/figures/cor_sam_osiq.png", width = 8, height = 8, units = "in", res = 300)

heatmap.2(cor.sam.osiq,
          cellnote = round(cor.sam.osiq, 2),
          notecol = "black",
          trace = "none",
          scale = "none",
          labRow = rownames(cor.sam.osiq),
          labCol = colnames(cor.sam.osiq),
          Rowv = FALSE,
          Colv = FALSE, # no tree
          dendrogram = "none",
          ColSideColors = as.vector(colour.osiq),
          RowSideColors = as.vector(colour.sam),
          colCol = as.vector(colour.osiq),
          colRow = as.vector(colour.sam),
          col = colourCode.cor,
          srtCol = 60,
          xlab = "OSIQ items",
          ylab = "SAM items",
          key = TRUE)

dev.off()


#---- supplemental figure 6: univariate relationships between SAM and OSIQ in an independent sample (2019 dataset) ----

scatter_samE_osiq_et =   
  df.sam.osiq.et %>%
  select(sam_epi_raw, osiq_mini_o, osiq_mini_s) %>%
  melt(id = "sam_epi_raw", measured = c("osiq_mini_o", "osiq_mini_s")) %>% 
  ggplot(aes(x = value, y = sam_epi_raw, colour = variable)) +
  geom_point(alpha = .4, position = "jitter") +
  geom_smooth(method = "lm", size = 2) +
  theme_bw() +
  labs(x = "OSIQ score", y = "SAM-episodic score") +
  scale_colour_manual(name = "OSIQ domain", 
                      values = colourCode.osiq,
                      labels = c("OSIQ-object", "OSIQ-spatial"))

scatter_samS_osiq_et =   
  df.sam.osiq.et %>%
  select(sam_spa_raw, osiq_mini_o, osiq_mini_s) %>%
  melt(id = "sam_spa_raw", measured = c("osiq_mini_o", "osiq_mini_s")) %>% 
  ggplot(aes(x = value, y = sam_spa_raw, colour = variable)) +
  geom_point(alpha = .4, position = "jitter") +
  geom_smooth(method = "lm", size = 2) +
  theme_bw() +
  labs(x = "OSIQ score", y = "SAM-spatial score") +
  scale_colour_manual(name = "OSIQ domain", 
                      values = colourCode.osiq,
                      labels = c("OSIQ-object", "OSIQ-spatial"))

scatter_samF_osiq_et =   
  df.sam.osiq.et %>%
  select(sam_fut_raw, osiq_mini_o, osiq_mini_s) %>%
  melt(id = "sam_fut_raw", measured = c("osiq_mini_o", "osiq_mini_s")) %>% 
  ggplot(aes(x = value, y = sam_fut_raw, colour = variable)) +
  geom_point(alpha = .4, position = "jitter") +
  geom_smooth(method = "lm", size = 2) +
  theme_bw() +
  labs(x = "OSIQ score", y = "SAM-future score") +
  scale_colour_manual(name = "OSIQ domain",
                      values = colourCode.osiq,
                      labels = c("OSIQ-object", "OSIQ-spatial"))

scatter_samM_osiq_et =
  df.sam.osiq.et %>%
  select(sam_sem_raw, osiq_mini_o, osiq_mini_s) %>%
  melt(id = "sam_sem_raw", measured = c("osiq_mini_o", "osiq_mini_s")) %>% 
  ggplot(aes(x = value, y = sam_sem_raw, colour = variable)) +
  geom_point(alpha = .4, position = "jitter") +
  geom_smooth(method = "lm", size = 2) +
  theme_bw() +
  labs(x = "OSIQ score", y = "SAM-semantic score") +
  scale_colour_manual(name = "OSIQ domain",
                      values = colourCode.osiq,
                      labels = c("OSIQ-object", "OSIQ-spatial"))

png("supplemental/figures/sam_osiq_scatter_et.png", width = 10, height = 10, units = "in", res = 600)

grid_arrange_shared_legend(scatter_samE_osiq_et, scatter_samS_osiq_et, 
                           scatter_samF_osiq_et, scatter_samM_osiq_et)

dev.off()


#---- supplemental figure 7: PLSC on SAM and OSIQ in an independent sample (2019 dataset) ----


# scree plot

png("supplemental/figures/ET_plsc_scree.png", 
    width = 5, height = 5, units = "in", res = 600)

PTCA4CATA::PlotScree(
  ev = plsc.sam.osiq.et$TExPosition.Data$eigs,
  max.ev = max(plsc.sam.osiq$TExPosition.Data$eigs)
)

dev.off()

# variable plot

loadings.et = rbind(plsc.sam.osiq.et$TExPosition.Data$fi,
                    plsc.sam.osiq.et$TExPosition.Data$fj)

colnames(loadings.et) = paste0("Dimension ", 1:ncol(loadings.et))

map.loadings.et = PTCA4CATA::createFactorMap(
  title = "All item loadings on dimensions 1 and 2",
  loadings.et,
  col.points = colour.plsc,
  alpha.points = .4,
  cex = 2,
  text.cex = 3,
  col.labels = colour.plsc)

aggmap.loadings.et =
  map.loadings.et$zeMap_background +
  map.loadings.et$zeMap_text +
  map.loadings.et$zeMap_dots +
  coord_fixed(ratio = 2.3)

png("supplemental/figures/ET_plsc_loadings.png", 
    width = 5, height = 5, units = "in", res = 600)

print(aggmap.loadings.et)

dev.off()

# LV1 OSIQ

Fi.et = plsc.sam.osiq.et$TExPosition.Data$fi

contrib.i.et = plsc.sam.osiq.et$TExPosition.Data$ci

contrib.i.signed.et = contrib.i.et * sign(Fi.et)

png("supplemental/figures/ET_plsc_lv1_osiq.png", 
    width = 5, height = 3, units = "in", res = 600)

PrettyBarPlot2(
  bootratio = round(100 * contrib.i.signed.et[,1]),
  threshold = 100 / nrow(contrib.i.signed.et),
  ylim = NULL,
  color4bar = gplots::col2hex(colour.osiq),
  color4ns = "gray75",
  plotnames = TRUE,
  main = "Important Contributions to LV1 from OSIQ items",
  ylab = "Signed Contributions")

dev.off()

# LV1 SAM

Fj.et = plsc.sam.osiq.et$TExPosition.Data$fj

contrib.j.et = plsc.sam.osiq.et$TExPosition.Data$cj

contrib.j.signed.et = contrib.j.et * sign(Fj.et)

png("supplemental/figures/ET_plsc_lv1_sam.png", 
    width = 5, height = 3, units = "in", res = 600)

PrettyBarPlot2(
  bootratio = round(100 * contrib.j.signed.et[,1]),
  threshold = 100/ nrow(contrib.j.signed.et),
  ylim = NULL,
  color4bar = gplots::col2hex(colour.sam),
  color4ns = "gray75",
  plotnames = TRUE,
  main = "Important Contributions to LV1 from SAM items",
  ylab = "Signed Contributions")

dev.off()

# LV2 OSIQ

png("supplemental/figures/ET_plsc_lv2_osiq.png", 
    width = 5, height = 3, units = "in", res = 600)

PrettyBarPlot2(
  bootratio = round(100 * contrib.i.signed.et[,2]),
  threshold = 100 / nrow(contrib.i.signed.et),
  ylim = NULL,
  color4bar = gplots::col2hex(colour.osiq),
  color4ns = "gray75",
  plotnames = TRUE,
  main = "Important Contributions to LV2 from OSIQ items",
  ylab = "Signed Contributions")

dev.off()

# LV2 SAM

png("supplemental/figures/ET_plsc_lv2_sam.png", 
    width = 5, height = 3, units = "in", res = 600)

PrettyBarPlot2(
  bootratio = round(100 * contrib.j.signed.et[,2]),
  threshold = 100/ nrow(contrib.j.signed.et),
  ylim = NULL,
  color4bar = gplots::col2hex(colour.sam),
  color4ns = "gray75",
  plotnames = TRUE,
  main = "Important Contributions to LV2 from SAM items",
  ylab = "Signed Contributions")

dev.off()


#---- supplemental figure 8: PCA on SAM split by gender ----

# female 

Fj.sam.female = pca.sam.female$ExPosition.Data$fj

jmap.sam.female = PTCA4CATA::createFactorMap(
  Fj.sam.female,
  col.points = colour.sam,
  col.labels = colour.sam,
  alpha.points. = .3)

jaggmap.sam.female =
  jmap.sam.female$zeMap_background +
  jmap.sam.female$zeMap_dots +
  jmap.sam.female$zeMap_text

png("supplemental/figures/sam_pca_female.png", width = 5, height = 6, units = "in", res = 300)

print(jaggmap.sam.female)

dev.off()

# male  

Fj.sam.male = pca.sam.male$ExPosition.Data$fj

jmap.sam.male = PTCA4CATA::createFactorMap(
  Fj.sam.male,
  col.points = colour.sam,
  col.labels = colour.sam,
  alpha.points. = .3)

jaggmap.sam.male =
  jmap.sam.male$zeMap_background +
  jmap.sam.male$zeMap_dots +
  jmap.sam.male$zeMap_text

png("supplemental/figures/sam_pca_male.png", width = 5, height = 6, units = "in", res = 300)

print(jaggmap.sam.male)

dev.off()


#---- supplemental figure 9: projection of BFI onto SAM ----


sup.bfi = supplementaryCols(
  SUP.DATA = expo.scale(df.bfi, scale = FALSE),
  res = pca.sam,
  center = FALSE,
  scale = FALSE)

Fj.sam.bfi = rbind(Fj.sam, sup.bfi$fjj)

colour.sam.bfi = c(colour.sam, colour.bfi)

constraints.bfi.sup = 
  lapply(minmaxHelper(Fj.sam, Fj.sup.bfi, axis1 = 1, axis2 = 2), "*", 1.1) 

jmap.bfi.sup = PTCA4CATA::createFactorMap(
  Fj.sam.bfi,
  col.points = colour.sam.bfi,
  col.labels = colour.sam.bfi,
  alpha.points. = .3,
  cex = 1,
  constraints = constraints.bfi.sup)

jaggmap.bfi.sup = 
  jmap.bfi.sup$zeMap_background +
  jmap.bfi.sup$zeMap_dots +
  jmap.bfi.sup$zeMap_text

png("supplemental/figures/bfi_projection.png", width = 5, height = 6, units = "in", res = 300)

print(jaggmap.bfi.sup)

dev.off()

#---- supplemental figure 10: PLSC on SAM and BFI ----

loadings.sam.bfi = 
  rbind(plsc.sam.bfi$TExPosition.Data$fi,
        plsc.sam.bfi$TExPosition.Data$fj)

colnames(loadings.sam.bfi) = paste0("Dimension ", 1:ncol(loadings.sam.bfi))

colour.sam.bfi = c(colour.sam, colour.bfi)

map.sam.bfi.loadings = createFactorMap(
  loadings.sam.bfi,
  col.points = colour.sam.bfi,
  col.labels = colour.sam.bfi,
  alpha.points = .2,
  cex = 1)

aggmap.sam.bfi.loadings = 
  map.sam.bfi.loadings$zeMap_background +
  map.sam.bfi.loadings$zeMap_text +
  map.sam.bfi.loadings$zeMap_dots

png("supplemental/figures/bfi_plsc.png", width = 6, height = 4, units = "in", res = 300)

print(aggmap.sam.bfi.loadings)

dev.off()

#---- supplemental figure 11: projection of PHQ onto SAM ----

sup.phq = supplementaryCols(
  SUP.DATA = expo.scale(df.phq, scale = FALSE),
  res = pca.sam,
  center = FALSE,
  scale = FALSE)

Fj.sup.phq = sup.phq$fjj

sup.constraints.phq= lapply(
  minmaxHelper(Fj.sam, Fj.sup.phq, axis1 = 1, axis2 = 2), "*", 1.1) 

jmap.sup.sam4phq = createFactorMap(
  Fj.sam,
  axis1 = 1, axis2 = 2,
  col.points = colour.sam,
  col.labels = colour.sam,
  alpha.points. = .3,
  cex = 1,
  constraints = sup.constraints.phq)

jaggmap.sup.sam4phq = 
  jmap.sup.sam4phq$zeMap_background +
  jmap.sup.sam4phq$zeMap_dots +
  jmap.sup.sam4phq$zeMap_text

jmap.sup.phq = createFactorMap(
  Fj.sup.phq,
  axis1 = 1, axis2 = 2,
  constraints = sup.constraints.phq,
  cex = 1,
  col.points = "lavenderblush4",
  col.labels = "lavenderblush4")

jaggmap.sup.phq= 
  jmap.sup.phq$zeMap_background +
  jmap.sup.phq$zeMap_text +
  jmap.sup.phq$zeMap_dots

jaggmap.sam.phq = 
  jaggmap.sup.sam4phq+
  jmap.sup.phq$zeMap_text +
  jmap.sup.phq$zeMap_dots

sup.phq = supplementaryCols(
  SUP.DATA = expo.scale(df.phq, scale = FALSE),
  res = pca.sam,
  center = FALSE,
  scale = FALSE)

Fj.sam.phq = rbind(Fj.sam, sup.phq$fjj)

colour.sam.phq = c(colour.sam, rep("lavenderblush4", 9))

constraints.phq.sup = 
  lapply(minmaxHelper(Fj.sam, Fj.sup.phq, axis1 = 1, axis2 = 2), "*", 1.1)

jmap.phq.sup = PTCA4CATA::createFactorMap(
  Fj.sam.phq,
  col.points = colour.sam.phq,
  col.labels = colour.sam.phq,
  alpha.points. = .3,
  cex = 1,
  constraints = constraints.phq.sup)

jaggmap.phq.sup = 
  jmap.phq.sup$zeMap_background +
  jmap.phq.sup$zeMap_dots +
  jmap.phq.sup$zeMap_text

png("supplemental/figures/phq_projection.png", width = 5, height = 6, units = "in", res = 300)

print(jaggmap.phq.sup)

dev.off()

#---- supplemental figure 12: PLSC on SAM and PHQ ----

loadings.sam.phq = 
  rbind(plsc.sam.phq$TExPosition.Data$fi,
        plsc.sam.phq$TExPosition.Data$fj)

colnames(loadings.sam.phq) = paste0("Dimension ", 1:ncol(loadings.sam.phq))

colour.sam.phq = c(colour.sam, rep("lavenderblush4", ncol(df.phq)))

map.sam.phq.loadings = createFactorMap(
  loadings.sam.phq,
  col.points = colour.sam.phq,
  col.labels = colour.sam.phq,
  alpha.points = .2,
  cex = 1)

aggmap.sam.phq.loadings = 
  map.sam.phq.loadings$zeMap_background +
  map.sam.phq.loadings$zeMap_text +
  map.sam.phq.loadings$zeMap_dots

png("supplemental/figures/phq_plsc.png", width = 6, height = 3, units = "in", res = 300)

print(aggmap.sam.phq.loadings)

dev.off()
