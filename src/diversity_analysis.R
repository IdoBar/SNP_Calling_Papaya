# load custom functions from github
devtools::source_gist("7f63547158ecdbacf31b54a58af0d1cc", 
                      filename = "util.R")
# install.packages("pak")
library(pacman)
# pak::pkg_install(c("green-striped-gecko/dartR", 'YuLab-SMU/ggtree', "hemstrow/snpR"))
# p_load_gh("Mikata-Project/ggthemr")
# p_load_gh("green-striped-gecko/dartR", update=FALSE)
# remotes::install_github('YuLab-SMU/ggtree')
# remotes::install_github("green-striped-gecko/dartR")

pacs <- c("tidyverse", 'RColorBrewer', 'paletteer', 'ggrepel', 'glue','janitor', 
          'hierfstat' , "qvalue",  'poppr', 'tidytree', 'StAMPP', 'mmod','phangorn', 'phytools',
          "here", 'showtext', 'ape', 'readxl', 'ggtree', "ggstar", "dendextend", 
          'pheatmap', 'igraph', 'adegraphics', "ggplotify", 'SNPfiltR',
          'hemstrow/snpR', "vcfR", "Mikata-Project/ggthemr")
# git_pacs <- c('hemstrow/snpR', "green-striped-gecko/dartR", "Mikata-Project/ggthemr")
pak::pak(c(pacs))
pacman::p_load(char = basename(pacs), install = FALSE)

# p_load(char = pacs, install = FALSE) 
# p_load_gh(git_pacs, install = FALSE)
# quadprog, raster, mvtnorm,pheatmap, ComplexHeatmap, strataG, diveRsity, DECIPHER, cowplot, ggtree
# library(dartR)
# BiocManager::install("qvalue")
# pacman::p_load_gh("Mikata-Project/ggthemr")
# remotes::install_github("green-striped-gecko/dartR")
# pacman::p_load_gh("green-striped-gecko/dartR")


pale_theme <- ggthemr::ggthemr("pale", text_size = 18, set_theme = FALSE)

#### Read DArT Data ####
populations_file <- here("data/Complete_coffee_genotype_list.xlsx")
samples_table <- read_excel(populations_file) %>% 
 #  select(id=Genotype, pop=Identification, Type, Origin) %>% 
  mutate(vcf_samples=paste0(Sequencing_ID, "A"),
        # across(where(is.character),.fns = ~replace_na(.x, "Unknown")),
         across(where(is.character),.fns = str_trim))
# read in vcf data
coffee_vcf <- read.vcfR("../BGI_WGRS/FB_SNP_calling_03_11_2025/Coffea_BGI_WGRS_ET39.fb.diploid.snps.Q30.mac5.nomiss0.8.poly.recode.vcf.gz")
gl <- adegenet::import2genind()
Mam00s#&J1nj1e&M@Y@l00k
samples_metadata <- tibble(id=indNames(gl)) %>% left_join(samples_table) %>%
  write_csv("data/Coffee_wgrs_metadata.csv")
# 

# metadata <- read_csv(here("data/DArT_metadata.csv"))
# double check matching names and length
length(indNames(gl)) == nrow(samples_metadata)
all(indNames(gl) == samples_metadata$id)

# check which samples names don't appear in the genotyping results
samples_table %>% filter(!id %in% samples_metadata$id) 

# set consistent strata
strata_factors <- gl@other$ind.metrics %>%
  # select(indNames, Location, State, Region, Year, Host, Pathotype, Path_rating, Haplotype, Taxa) %>%
  mutate(across(c("pop", "Type", "Origin"), factor))
summary(strata_factors) 

#### Exploratory Data Analysis ####
# ploidy(Arab_gl) <- 1
gl.report.monomorphs(gl)  # all loci are polymorphic
# check missing rates
gl.report.callrate(gl, method = "loc", v=3)
gl.report.callrate(gl, method = "ind", v = 3) # the 2 individuals from Spain are the most remote
# check reproducibility (RepAVG)
gl.report.RepAvg(gl)
# check multiple loci per tag and reproducibility
gl.report.secondaries(gl)

# add strata information 
gl@strata <- strata_factors 
setPop(gl) <- ~pop
# apply filters for export #####
# find loci with high number of SNPs
snps_in_tag_thres <- 4
drop_loci <- gl@other$loc.metrics %>% count(CloneID) %>% 
  filter(n>=snps_in_tag_thres)
# keep_loci <- Arab_gl@other$loc.metrics %>% count(CloneID) %>% filter(n<snps_in_tag_thres)
# idebtifying which SNPs came from loci with high number of SNPs
drop_snps <- str_subset(locNames(gl), paste(drop_loci$CloneID, collapse = "|"))
# drop_snps <- locNames(gl)[sub("^(\\d+)-.+", "\\1", locNames(gl)) %in% drop_loci$CloneID]
# drop pooled individuals
drop_ind <- str_subset(indNames(gl), "pool")

# Entire population ####
missing_thres <- 0.2

glsub <- gl %>% gl.compliance.check() %>%  gl.drop.loc(drop_snps, verbose=3) %>% 
  gl.drop.ind(drop_ind, verbose = 3) %>% 
  gl.filter.reproducibility(., threshold = 0.95, verbose = 3) %>% 
  gl.filter.callrate(., method = "loc", threshold = 0.95, verbose = 3, mono.rm = TRUE, recalc=TRUE) %>% 
  # gl.recalc.metrics() %>% 
  gl.filter.callrate(., method = "ind", threshold = 0.9, verbose = 3, mono.rm = TRUE, recalc=TRUE)  %>%  
  # gl.recalc.metrics() %>% 
  gl.filter.callrate(., method = "loc", threshold = 0.95, verbose = 3, mono.rm = TRUE, recalc=TRUE)

# Subset strata table according to our filtered dataset, remove empty factor levels and assign back to filtered dataset
glsub_strata <- strata_factors %>% 
  # mutate(Origin=fct_recode(State, ICARDA = "SPAIN")) %>% 
  filter(id %in% indNames(glsub))  %>% mutate(across(where(is.factor), fct_drop))
glsub@strata <- glsub_strata


load("data/Papaya_DArT_analysis_03_02_2023.RData")

# save tags as fasta for homology search
dart_for_fasta <- glsub@other$loc.metrics %>% as_tibble() %>% 
  mutate(tag_id=paste0("tag_", as.character(CloneID)),
         tag_seq= as.character(AlleleSequenceRef)) %>% 
  select(tag_seq, SnpPosition, SNP, tag_id) 

tag_seqs <- map(dart_for_fasta$tag_seq, .f = seqinr::s2c) %>% 
  set_names(dart_for_fasta$tag_id)
# write fasta using my function (better and easier)!
dart_for_fasta %>% select(tag_id, tag_seq) %>% 
  write_fasta(., filename = "output/papaya_dart_tags.fasta")

# save as vcf
# pak::pkg_install("thierrygosselin/radiator")
# radiator::genomic_converter(glsub, output = "vcf", filename = "data/Papaya_DArT_filtered_SNPs.vcf")
gl2vcf(glsub, "Papaya_DArT_filtered_SNPs", outpath = "data", 
       plink_path = "C:\\ProgramData\\plink_win64_20230116",
       snp_pos = "ChromPosSnp_Papaya_v1.0", snp_chr = "Chrom_Papaya_v1.0")


# Popgen analysis #####
papaya_genind <- gl2gi(glsub, probar = TRUE)
papaya_genind@strata <- glsub@strata

# calculate within-pop Statistics with dartR (need to filter each pop to remove missing loci)
gl_by_pop <- seppop(glsub)
popgen_basic_stats <- gl_by_pop %>%
  imap_dfr(~ {
    LogMsg(glue("Processing pop {.y}"))
    popgl_sub <- gl.filter.callrate(.x, method = "loc", threshold = 0.65, verbose = 0, recalc=TRUE) 
    as_tibble(as.list(gl.basic.stats(popgl_sub)$overall)) %>% mutate(pop=.y)
  })
# test <- gl.basic.stats(gl_by_pop$`1B`)
# as_tibble(as.list(gl.basic.stats(gl_by_pop[["Holland_5"]] %>% 
#                  gl.filter.callrate(., method = "loc", threshold = 0.65, verbose = 0, recalc=TRUE))$overall)) %>% mutate(pop="Holland_5")
# 
popgen_basic_stats %>% pivot_wider(values_from = values, names_from = pop)
popgen_basic_stats %>% relocate(pop) %>% 
  write_xlsx(., "output/papaya_popgen.xlsx",
           sheet = "basic_genotype_stats",
           overwritesheet = TRUE)

# with poppr (takes time!)
papaya_poppr_stats <- poppr(papaya_genind)
write_xlsx(papaya_poppr_stats, "output/papaya_popgen.xlsx",
           sheet = "poppr_genotype_stats",
           overwritesheet = TRUE)
save.image(filedate("Papaya_DArT_analysis", ".RData", outdir = "data"))
# Pairwise Fst
papaya_pairFst <- stamppFst(glsub, nboots = 1000, nclusters = 4)
papaya_pairFst_results <- expand.grid(rownames(papaya_pairFst$Fsts), 
                                      colnames(papaya_pairFst$Fsts),
                                      stringsAsFactors = FALSE) %>% 
  mutate(Fst = as.vector(papaya_pairFst$Fsts)) %>% 
  filter(!is.na(Fst)) %>% rename(Pop1=Var1, Pop2=Var2) %>% 
  arrange(desc(Fst))
write_xlsx(papaya_pairFst_results, "output/papaya_popgen.xlsx", 
           sheet = "pop_Fst",
           overwritesheet = TRUE)
papaya_pairFst_results %>% pivot_wider(names_from = "Pop2", values_from = "Fst") %>% 
  .[c(3,1,6,4,2,5),c("Pop1", "Fairbairn", "Waruma", "Cania", "Somerset", "Callide", "ACH")] %>%
  write_xlsx(., "output/papaya_popgen.xlsx", 
             sheet = "pop_Fst_mat",
             overwritesheet = TRUE)
# population level with mmod
# test <- seppop(papaya_genind)
papaya_mmod <- diff_stats(papaya_genind)
as_tibble(as.list(papaya_mmod$global)) %>% mutate(pop="All") %>% 
  write_xlsx(., "output/papaya_popgen.xlsx", 
             sheet = "pop_Gst", overwritesheet = TRUE)

# calculate inbreeding 
pak::pkg_install("mastoffel/inbreedR")
# library(inbreedR)
# overall
papaya_gt <- genind2df(papaya_genind, oneColPerAll = TRUE)
papaya_snps <- popNames(glsub) %>% 
  map(~inbreedR::convert_raw(papaya_gt %>% filter(pop==.x) %>% .[-1]) )
names(papaya_snps) <- popNames(glsub)

# calculate sMLH
papaya_mlh <- inbreedR::sMLH(papaya_snps$`papaya park`)
# r2 between inbreeding and heterozygosity
papaya_r2_hf <- papaya_snps %>% 
  imap_dfr(~ inbreedR::r2_hf(.x, 
                             type = "snps",nboot = 100, 
                             CI = 0.95, parallel = TRUE, 
                             ncores = 4)[c("r2_hf_full", "CI_boot")] %>% as_tibble(.) %>% mutate(CI=c("CI_min", "CI_max"), pop=.y) %>% 
             pivot_wider(values_from = "CI_boot", names_from = "CI"))
papaya_r2_hf %>% relocate(pop) %>% write_xlsx(., "output/papaya_popgen.xlsx", sheet = "pop_inbreeding_MLH_r2",
                                              overwritesheet = TRUE)
# calculate g2
papaya_g2 <- papaya_snps %>% 
  imap_dfr(~ as_tibble(inbreedR::g2_snps(.x, 
                                         nperm = 100, nboot = 100, 
                                         CI = 0.95, parallel = TRUE, 
                                         ncores = 4)[c("nloc", "g2", "g2_se", "p_val")]) %>% 
             mutate(pop=.y))

write_xlsx(papaya_g2, "output/papaya_popgen.xlsx", sheet = "pop_inbreeding",
           overwritesheet = TRUE)

# define set colour palettes for each factor
# Population
# Sex

# specify colours
# paletteer_d("RColorBrewer::RdYlGn", direction = -1)
paletteer_d("RColorBrewer::Set1")
my_colours = list(
  pop = levels(glsub@strata$pop) %>% 
    setNames(as.character(paletteer_d("ggsci::default_igv", length(.))), .),
  Type = levels(glsub@strata$Type) %>% 
    setNames(c("#E05723", "#FFCC22"), .), #YEllow/Red
  Origin = levels(glsub@strata$Origin) %>% 
    setNames(as.character(paletteer_d("ggsci::category10_d3", length(.))), .)
)

# snpR Analysis ######
pak::pkg_install("hemstrow/snpR@fb2904f")
save(glsub, file = "data/papaya_glsub.RDS")

glsub_loc_info <- tibble(snpID=locNames(glsub), 
                         Chrom=glsub@other$loc.metrics$Chrom_Papaya_v1.0,
                         Pos = glsub@other$loc.metrics$ChromPosSnp_Papaya_v1.0) %>% 
  mutate(across(.fns = ~str_replace_all(string = .x, replacement = "_", pattern = "\\.")))

glsub_sample_metadata <- glsub$strata %>% as_tibble() %>% 
  mutate(across(.fns = ~str_replace_all(string = .x, replacement = "_", pattern = "\\.")))

papaya_snpr <- snpR::convert_genlight(glsub)
papaya_snpr <- snpR::read_vcf("data/Papaya_DArT_filtered_SNPs.vcf", 
                              snp.meta = glsub_loc_info,
                              sample.meta = glsub_sample_metadata)

summarize_facets(papaya_snpr, facets = "pop")
parentage_analysis <- run_sequoia(papaya_snpr, facets = "pop")


# PCA ####
# Convert to allele frequencies (manage missing data)
X <- adegenet::tab(glsub, freq = TRUE, NA.method = "mean")
pca1 <- ade4::dudi.pca(X, scale = FALSE, scannf = FALSE, nf=20)

pca_data <- pca1$li %>% rownames_to_column("id")  %>%
  left_join(glsub_strata) # %>% mutate(alpha=ifelse(Origin=="ICARDA", 0.65, 0.35))
# %>% mutate(Host=fct_relevel(fct_infreq(fct_lump_min(Host, min=6)), "Other",
#                                                   after = Inf),
#                                  State=factor(State, levels = c("QLD", "NSW", "VIC", "SA", "WA", "SPAIN")))
# Calculate variance for each component (columns starting with 'C')
pca_var <- pca1$li %>% summarise_at(vars(starts_with("Axis")), var)
# Calculate the percentage of each component of the total variance
percentVar <- pca_var/sum(pca_var)
# Define colour palette, but get rid of the awful yellow - number 6
# pal <- brewer.pal(9, "Set1")
# pal2 <- brewer.pal(8, "Dark2")
# text_cols <- adjustcolor( c("grey15", "dodgerblue3"), alpha.f = 0.8)
# seq_face <- setNames(c("bold", "plain"),unique(pca_data$Sequenced))
shapes <- c(16,15)
# Isolate to annotate
# outliers <- pca_data %>% filter(Origin=="ICARDA")
# 
#### PCA plot 1-2 geographic ####
# Create the plot for C1 and C2
plot_comps <- 1:2
plot_axes <- paste0("Axis", plot_comps)

pca_12 <- ggplot(pca_data, aes(x=!!sym(plot_axes[1]), y=!!sym(plot_axes[2]), 
                               color=Origin)) + # shape=State,    alpha=alpha)
  aes(label=pop, shape=Type) + 
  #  geom_jitter(alpha = 0.65, size = 4, width = 0.95, height = 0.95) +
  geom_point(alpha = 0.65, size=4) +
  pale_theme$theme + 
  # scale_color_paletteer_d("ggsci::category10_d3") + 
  scale_colour_manual(values = my_colours$Origin) +
  guides(colour = guide_legend(override.aes = list(alpha = 1)),
         shape="none") + 
  scale_shape_manual(values = shapes) +
  # scale_color_manual(values=seq_cols) +
  # # scale_size_continuous(range = c(3,6)) +
  # guides(shape = guide_legend(override.aes = list(size = 5 , fill=pal[2], color="grey15"), order = 1),  #, #
  #        fill = guide_legend(override.aes = list(shape = 21))) +   #  guide_legend(override.aes = list(size = 5, pch=shapes[1]), order = 2)
  labs( x=glue::glue("C{plot_comps[1]}: {round(percentVar[plot_comps[1]]*100,2)}% variance"),
        y=glue::glue("C{plot_comps[2]}: {round(percentVar[plot_comps[2]]*100,2)}% variance"), colour="Genotype")
# Save plot to a pdf file
ggsave(filedate(glue::glue("All_samples_PCA_Origin_PC{paste(plot_comps, collapse='-')}"),
                ext = ".pdf", outdir = here("output/plots")), width = 8, height=6)

ggplotly(pca_12) %>% hide_legend()

plot_comps <- 3:4
plot_axes <- paste0("Axis", plot_comps)
# ggthemr("pale", text_size = 20)
pca_34 <- ggplot(pca_data, aes(x=!!sym(plot_axes[1]), y=!!sym(plot_axes[2]), 
                               colour=pop, shape = Type)) + # shape=State,
  # geom_point(alpha = 0.65, size=4) +
  geom_point(alpha = 0.65, size=4) +
  pale_theme$theme + 
  # scale_color_paletteer_d("ggsci::category10_d3") + 
  # scale_fill_brewer(palette = "RdYlGn", direction = -1) +
  # scale_fill_paletteer_d("RColorBrewer", "RdYlGn", direction = -1 ) + # ggsci, category10_d3
  scale_shape_manual(values = shapes) +
  scale_color_manual(values=my_colours$pop) +
  # # scale_size_continuous(range = c(3,6)) +
  # guides(shape = guide_legend(override.aes = list(size = 5 , fill=pal[2], color="grey15"), order = 1),  #, #
  #        fill = guide_legend(override.aes = list(shape = 21))) +   #  guide_legend(override.aes = list(size = 5, pch=shapes[1]), order = 2)
  # plot_theme(baseSize = 20) + #  size="Year",
  labs( x=glue::glue("C{plot_comps[1]}: {round(percentVar[plot_comps[1]]*100,2)}% variance"),
        y=glue::glue("C{plot_comps[2]}: {round(percentVar[plot_comps[2]]*100,2)}% variance"),
        colour="Genotype")

cowplot::plot_grid(pca_12, pca_34, labels = c('A', 'B'), label_size = 22)
# Save plot to a pdf file
ggsave(filedate(glue::glue("All_samples_PCA_Origin_composite"),
                ext = ".pdf", outdir = here("output/plots")), width = 14, height=6)

#### DAPC analysis ####
# adegenet::adegenetTutorial("dapc")
# https://github.com/thibautjombart/adegenet/blob/master/tutorials/tutorial-dapc.pdf
# https://grunwaldlab.github.io/Population_Genetics_in_R/DAPC.html
# papaya_genind <- gl2gi(glsub)
# ploidy(papaya_genind) <- 1
# papaya_genind@strata <- strata_factors %>% 
#   mutate(Origin=fct_recode(State, ICARDA = "SPAIN")) %>% 
#   filter(id %in% indNames(Arab_genind)) %>% mutate_if(is.factor, ~fct_drop(.))


papaya_dapc <- dapc(papaya_genind, var.contrib = TRUE, scale = FALSE, 
                    n.pca = 20, n.da = nPop(papaya_genind) - 1)
scatter(papaya_dapc, cell = 0, cstar = 0, mstree = TRUE, lwd = 2, lty = 2)
# cross-validating the number of PCs to retain
myInset <- function(dapc_obj, text_x=5, text_y=90, text_cex=0.85){
  # dapc_obj=Year_dapc$DAPC
  temp <- dapc_obj$pca.eig
  temp <- 100* cumsum(temp)/sum(temp)
  plot(temp, col=rep(c("black","lightgrey"),
                     c(dapc_obj$n.pca,1000)), ylim=c(0,100),
       xlab="PCA axis", ylab="Cumulated variance (%)",
       cex=1, pch=20, type="h", lwd=2)
  text(text_x, text_y,  sprintf("PCs=%s", dapc_obj$n.pca),
       cex=text_cex, pos=4)
}

# now repeat 100 times
# plot_dapc <- function(gen_obj, strata_var, plot_filename, plot_widt=8, plot_heig=7,
#                       pal_pack="rcartocolor", pal="Bold",
#                       pca_range=6:18, reps=500, shapes=15:19, label_groups=FALSE, leg_pos="topright",
#                       pca_pos="bottomleft", ncores=parallel::detectCores()-1){
#   # calculate optimal number of PCs to retain
#   Arabx <- xvalDapc(tab(gen_obj, NA.method = "mean"), strata(gen_obj)[[strata_var]],
#                         n.pca = pca_range, n.rep = reps,
#                         parallel="snow", ncpus=ncores)
# 
#   # save plot
#   pdf(file = plot_filename, width = plot_widt, height=plot_heig)
#   scatter(Arabx$DAPC,  cex = 2, legend = TRUE, pch=shapes, 
#           col = paletteer_d(!!pal_pack, !!pal), scree.da=FALSE, 
#           clabel = label_groups, posi.leg = leg_pos, scree.pca = FALSE, 
#           cleg = 0.75, xax = 1, yax = 2, inset.solid = 1)
#   add.scatter(myInset(Arabx$DAPC), posi=pca_pos,
#               inset=c(-0.03,-0.01), ratio=.25,
#               bg=transp("white"))
#   dev.off()
#   return(invisible(Arabx))
# }

# set constants
pca_range=5:10
reps=500
shapes=16 # c(15:18, 9,25)
ncores=parallel::detectCores()-1
# cross-validating the number of PCs to retain (rep=1000)
# cluster by population
strata_var <- "pop"
# pal <- "Bold"
pop_dapc <- xvalDapc(tab(papaya_genind, NA.method = "mean"), papaya_genind@strata[[strata_var]],
                     n.pca = pca_range, n.rep = reps,
                     parallel="snow", ncpus=ncores)

# save plot
library(patchwork)
# scatter_12 <- ggplotify::as.ggplot(~{scatter(pop_dapc$DAPC,  cex = 2, legend = TRUE, pch=shapes, 
#                                              col = my_colours[[strata_var]], scree.da=FALSE, 
#                                              clabel = FALSE, posi.leg = "bottomright", 
#                                              scree.pca = FALSE, ylab = "PC2", xlab = "PC1", 
#                                              cleg = 0.85, xax = 1, yax = 2, inset.solid = 1); 
#   add.scatter(myInset(pop_dapc$DAPC), posi="bottomleft",
#               inset=c(0.01,-0.01), ratio=.2,
#               bg=transp("white")) })
# scatter_12 + labs(x="PC1")
wrap_elements(panel = ~{scatter(pop_dapc$DAPC,  cex = 2, legend = TRUE, pch=shapes, 
                                col = my_colours[[strata_var]], scree.da=FALSE, 
                                clabel = FALSE, posi.leg = "bottomright", 
                                scree.pca = FALSE, ylab = "PC2", xlab = "PC1", 
                                cleg = 0.85, xax = 1, yax = 2, inset.solid = 1); 
  add.scatter(myInset(pop_dapc$DAPC), posi="bottomleft",
              inset=c(0.01,-0.01), ratio=.2,
              bg=transp("white")) }) + ggtitle('A') +
  wrap_elements(panel = ~scatter(pop_dapc$DAPC,  cex = 2, legend = TRUE, pch=shapes, 
                                 col = my_colours[[strata_var]], scree.da=FALSE, 
                                 clabel = FALSE, posi.leg = "bottomright", 
                                 scree.pca = FALSE, 
                                 cleg = 0.85, xax = 3, yax = 4, inset.solid = 1)) + ggtitle('B')
pdf(file = filedate(sprintf("papaya_DArT_DAPC_%s_analysis", strata_var), 
                    ".pdf", outdir = "output/plots"), width = 7, height=6)
scatter(pop_dapc$DAPC,  cex = 2, legend = TRUE, pch=shapes, 
        col = my_colours[[strata_var]], scree.da=FALSE, 
        clabel = FALSE, posi.leg = "bottomright", scree.pca = FALSE, 
        cleg = 0.85, xax = 1, yax = 2, inset.solid = 1)
add.scatter(myInset(pop_dapc$DAPC), posi="bottomleft",
            inset=c(0.01,-0.01), ratio=.2,
            bg=transp("white"))

dev.off()
# plot 3rd and 4th PCs
pdf(file = filedate(sprintf("papaya_DArT_DAPC_%s_PC3-4", strata_var), 
                    ".pdf", outdir = "output/plots"), width = 7, height=6)
scatter(pop_dapc$DAPC,  cex = 2, legend = TRUE, pch=shapes, 
        col = my_colours[[strata_var]], scree.da=FALSE, 
        clabel = FALSE, posi.leg = "bottomright", scree.pca = FALSE, 
        cleg = 0.85, xax = 3, yax = 4, inset.solid = 1)
dev.off()

# produce structure-like plots
# compoplot(pop_dapc$DAPC,col = my_colours[[strata_var]], posi = 'top')
dapc.results <- as_tibble(pop_dapc$DAPC$posterior) %>% 
  mutate(Original_Pop=fct_relevel(pop(glsub), "Solo"),
         Sample = indNames(glsub)) %>% 
  pivot_longer(where(is.numeric), names_to = "Assigned_Pop", 
               values_to = "Posterior_membership_probability") %>% 
  mutate(Assigned_Pop=fct_relevel(factor(Assigned_Pop), "Solo"))
# plot with ggplot2
ggplot(dapc.results, aes(x=Sample, y=Posterior_membership_probability, fill=Assigned_Pop)) +
  geom_bar(stat='identity') +
  scale_fill_manual(values = my_colours$pop) + 
  facet_wrap(~Original_Pop, scales = "free", nrow = 2) +
  theme_minimal(base_size = 18) + 
  labs(y="Membership Prob.", fill = "Assigned Pop.") +
  # pale_theme$theme + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8),
        panel.grid = element_blank())
ggsave("output/plots/papaya_pops_compoplot.pdf", width = 12, height = 5)

# Calculate Tree, heatmap and clusters ####

# papaya_dist <- gl.dist.pop(glsub)
tree <- aboot(glsub, tree = "upgma", distance = "bitwise.dist", sample = 1000, 
              showtree = T, rroot = T, 
              cutoff = 50, quiet = T,
              missing="ignore")
D <- bitwise.dist(glsub)
x <- as.vector(D)
y <- as.vector(as.dist(cophenetic(tree)))
plot(x, y, xlab="original pairwise distances", ylab="pairwise distances on the tree",
     main="Is UPGMA  appropriate?", pch=20, col=transp("black",.1), cex=3)
abline(lm(y~x), col="red")
# Test Nj tree #
nj_tre <- nj(D)

y_nj <- as.vector(as.dist(cophenetic(nj_tre)))
plot(x, y_nj, xlab="original pairwise distances", ylab="pairwise distances on the tree", 
     main="Is NJ appropriate?", pch=20, col=transp("black",.1), cex=3)
abline(lm(y~x), col="red")

# Try ML tree #####
# read DArT data and export as concatenated fasta
gl.read.dart("data/Report_DCari22-7639_SNP_2.csv",
             ind.metafile = "data/Papaya_DArT_metadata.csv" ) %>%
# indNames(gl)  
  gl.keep.ind(., ind.list = indNames(glsub), recalc = T, mono.rm = T, verbose = 3) %>%
  gl.filter.callrate(., threshold = 1, mono.rm=TRUE, recalc=TRUE, verbose = 3 ) %>% 
  # gl.keep.loc(., loc.list = locNames(glsub), verbose = 3) %>%
gl2fasta(., outfile = "Papaya_DArT_no_missing_selected_inds_concat.fasta", outpath = "data",
         verbose = 3)
# import data as fasta
papaya_dart_DNA <- fasta2DNAbin(file="data/Papaya_DArT_tags_concat.fasta")
dna2 <- as.phyDat(papaya_dart_DNA) #assign the original dna sequences data as a phyDat object...
# test which distance model to use
mt <- modelTest(dna2)
print(mt)

# create initial NJ tree
tre.ini <- nj(dist.dna(papaya_dart_DNA,model="TN93"))

str(dna2)
pml(tre.ini, dna2)
# export as Phylip
fasta <- ape::read.FASTA("data/Papaya_DArT_tags_concat.fasta")
colWidth <- max(sapply(fasta, length))
ape::write.dna(fasta, file="data/Papaya_DArT_tags_concat.phylip", 
               colsep="", format="sequential",
               colw=colWidth)


strata_factors$Type[strata_factors$pop=="F5-67"] <- "Yellow"
tree_data <- as_tibble(tree) %>% 
  left_join(strata_factors, by = c("label"="id")) %>% 
  mutate(id=label, label=pop) %>% 
  as.treedata()
tree_data %>% as_tibble() %>% arrange(pop)
# visualise the tree
ggtree(tree_data, layout="daylight") + #geom_treescale() #+ 
  geom_tiplab(aes(angle=angle, colour=Type), size=4) + 
  # geom_tippoint(aes(colour=pop, shape=Type)) +  # , shape=sex)
  # scale_shape_manual(values = shapes) +
  scale_color_manual(values=my_colours$Type) +
  # plot_theme("bw", 18) + 
  theme_tree2(text = element_text(size=8)) + 
  theme(legend.position="right") +
  # guides(colour = guide_legend(override.aes = list(size = 6))) +
  guides(colour = "none") +
  labs(colour = "Genotype", x = "Genetic distance (proportion of loci that are different)") #  shape="Sex"

# save the plot
ggsave("output/plots/papaya_ind_tree.pdf", width= 12, height = 14)

save.image(filedate("Papaya_DArT_analysis", ext = ".RData", outdir = "data"))
# population-level Tree ####

# Analysis
set.seed(999)
pop_tree <- papaya_genind   %>%
  genind2genpop(pop = ~pop) %>%
  aboot(cutoff = 50, quiet = TRUE, sample = 1000, distance = provesti.dist)
# pop_tree_data <- as_tibble(tree) %>% left_join(strata_factors, by = c("label"="id")) %>% 
#   as.treedata()
pop_tree_data <- as_tibble(pop_tree) %>% 
  left_join(strata_factors %>% group_by(pop) %>% 
              slice(1), by = c("label"="pop")) %>% 
  # mutate(id=label, label=pop) %>% 
  as.treedata()
# visualise the tree
ggtree(pop_tree_data, layout="daylight") + #geom_treescale() #+ 
  geom_tippoint(aes(colour=Type, shape=Type), size=0, alpha=0.8) +  # , shape=sex)
  scale_shape_manual(values = shapes) +
  geom_highlight(node=MRCA(pop_tree_data, c("Hybrid_7", "C3-15")), 
                           fill="coral", alpha=0.3) +
  geom_highlight(node=MRCA(pop_tree_data, c("First_Lady", "Solo")), 
                 fill="coral2", alpha=0.3) +
  geom_highlight(node=MRCA(pop_tree_data, c("1B", "ML2-8-17")), 
                 fill=my_colours$Type["Yellow"], alpha=0.3) +
  
  # scale_color_manual(values=my_colours$pop) +
  # plot_theme("bw", 18) + 
  geom_tiplab(aes(colour=Type), size=3, fontface =2) +
  scale_color_manual(values=my_colours$Type) +
  # theme_tree2(text = element_text(size=8)) + 
  theme(legend.position="right") +
  coord_cartesian(clip="off") + 
  # ggplot2::xlim(0, 0.06) + 
  guides(color = guide_legend(override.aes = list(size = 4, shape=16, label=""))) +
  labs(colour = "Type") #  shape="Sex"

# save the plot
ggsave("output/plots/papaya_pop_tree_shaded_clades.pdf", width= 8, height = 7)

# visualise the tree (circular)
ggtree(pop_tree_data, layout = 'circular') + #geom_treescale() #+  , branch.length='none'
  geom_tippoint(aes(colour=Type, shape=Type), size=0, alpha=0.8) +  # , shape=sex)
  scale_shape_manual(values = shapes) +
  geom_highlight(node=MRCA(pop_tree_data, c("Hybrid_7", "C3-15")), 
                 fill="coral", alpha=0.3) +
  geom_highlight(node=MRCA(pop_tree_data, c("First_Lady", "Solo")), 
                 fill="coral2", alpha=0.3) +
  geom_highlight(node=MRCA(pop_tree_data, c("1B", "ML2-8-17")), 
                 fill=my_colours$Type["Yellow"], alpha=0.3) +
  
  # scale_color_manual(values=my_colours$pop) +
  # plot_theme("bw", 18) + 
  geom_tiplab(aes(colour=Type), size=3, fontface =2) +
  scale_color_manual(values=my_colours$Type) +
  # theme_tree2(text = element_text(size=8)) + 
  theme(legend.position="right") +
  # coord_cartesian(clip="off") + 
  # ggplot2::xlim(0, 0.06) + 
  guides(color = guide_legend(override.aes = list(size = 4, shape=16, label=""))) +
  labs(colour = "Type") #  shape="Sex"

# save the plot
ggsave("output/plots/papaya_pop_tree_circ.pdf", width= 8, height = 7)

# Network analysis ####
dist_matrix <- bitwise.dist(glsub, euclidean = TRUE)
msn <- poppr.msn(glsub, dist_matrix, showplot = FALSE, include.ties = T)
