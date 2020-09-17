#********************************************************************#
#                                                                    #
#                              R v3.5.1                              #
#                                                                    #
#         2016 16S rRNA Multivariate Analysis of 78 Samples          #
#                                                                    #
#                       Author: Shawn Kroetsch                       #
#                                                                    #
#********************************************************************#

## 1) Set Working Directory
setwd("/data/shawn/2016_Multivariate_Analysis")


## 2) Install Packages
install.packages("vegan")
install.packages("devtools")
install_github("vqv/ggbiplot")

library(ggplot2)
library(plyr)
library(car)
library(MASS)
library(vegan)
library(devtools)
library(ggbiplot)


## 3) Presence & Absence - Host Invertebrate Taxonomy (Genus)

# Generating Sorensen Dissimilarity Matrix:
presence_absence_2016_genus <- read.csv("/data/shawn/2016_Multivariate_Analysis/presence_absence_2016_genus.csv")
presence_absence_2016_genus <- as.data.frame(presence_absence_2016_genus)
row.names(presence_absence_2016_genus) <- presence_absence_2016_genus$Genus
presence_absence_2016_genus <- presence_absence_2016_genus[,-1]
presence_absence_2016_genus_DM <- vegdist(presence_absence_2016_genus, method = "bray", binary = TRUE)

# Unconstrained Principal Coordinates Analysis (PCoA)
presence_absence_2016_genus.ord <- dbrda(presence_absence_2016_genus_DM ~ 1, add = TRUE)
summary(presence_absence_2016_genus.ord)
PCoA_scores <- scores(presence_absence_2016_genus.ord, display="si")
SampleID_scores <- data.frame(PCoA_scores)

# Assign Taxonomy
variables_2016 <- read.csv("/data/shawn/2016_Multivariate_Analysis/variables_2016.csv")
variables_2016 <- as.data.frame(variables_2016)
row.names(variables_2016) <- variables_2016$SampleID
variables_2016 <- variables_2016[,-1]
Genus <- variables_2016$Genus
SampleID_scores$Genus <- Genus

# Calculate Multivariate Dispersions Using Betadisper
mod <- betadisper(presence_absence_2016_genus_DM, SampleID_scores$Genus)

# Permutation Test
permutest(mod, pairwise = TRUE, permutations = 999)

# Tukey's Honest Significant Differences
(mod.HSD <- TukeyHSD(mod))
plot(mod.HSD)


## 4) Presence & Absence - Host Invertebrate Taxonomy (Family)

# Generating Sorensen Dissimilarity Matrix:
presence_absence_2016_family <- read.csv("/data/shawn/2016_Multivariate_Analysis/presence_absence_2016_family.csv")
presence_absence_2016_family <- as.data.frame(presence_absence_2016_family)
row.names(presence_absence_2016_family) <- presence_absence_2016_family$Family
presence_absence_2016_family <- presence_absence_2016_family[,-1]
presence_absence_2016_family_DM <- vegdist(presence_absence_2016_family, method = "bray", binary = TRUE)

# Unconstrained Principal Coordinates Analysis (PCoA)
presence_absence_2016_family.ord <- dbrda(presence_absence_2016_family_DM ~ 1, add = TRUE)
summary(presence_absence_2016_family.ord)
PCoA_scores <- scores(presence_absence_2016_family.ord, display="si")
SampleID_scores <- data.frame(PCoA_scores)

# Assign Taxonomy
variables_2016 <- read.csv("/data/shawn/2016_Multivariate_Analysis/variables_2016.csv")
variables_2016 <- as.data.frame(variables_2016)
row.names(variables_2016) <- variables_2016$SampleID
variables_2016 <- variables_2016[,-1]
Family <- variables_2016$Family
SampleID_scores$Family <- Family

# Calculate Multivariate Dispersions Using Betadisper
mod <- betadisper(presence_absence_2016_genus_DM, SampleID_scores$Family)

# Permutation Test
permutest(mod, pairwise = TRUE, permutations = 999)

# Tukey's Honest Significant Differences
(mod.HSD <- TukeyHSD(mod))
plot(mod.HSD)


## 5) Presence & Absence - Host Invertebrate Taxonomy (Order)

# Generating Sorensen Dissimilarity Matrix:
presence_absence_2016_order <- read.csv("/data/shawn/2016_Multivariate_Analysis/presence_absence_2016_order.csv")
presence_absence_2016_order <- as.data.frame(presence_absence_2016_order)
row.names(presence_absence_2016_order) <- presence_absence_2016_order$Order
presence_absence_2016_order <- presence_absence_2016_order[,-1]
presence_absence_2016_order_DM <- vegdist(presence_absence_2016_order, method = "bray", binary = TRUE)

# Unconstrained Principal Coordinates Analysis (PCoA)
presence_absence_2016_order.ord <- dbrda(presence_absence_2016_family_DM ~ 1, add = TRUE)
summary(presence_absence_2016_order.ord)
PCoA_scores <- scores(presence_absence_2016_order.ord, display="si")
SampleID_scores <- data.frame(PCoA_scores)

# Assign Taxonomy
variables_2016 <- read.csv("/data/shawn/2016_Multivariate_Analysis/variables_2016.csv")
variables_2016 <- as.data.frame(variables_2016)
row.names(variables_2016) <- variables_2016$SampleID
variables_2016 <- variables_2016[,-1]
Order <- variables_2016$Order
SampleID_scores$Order <- Order

# Calculate Multivariate Dispersions Using Betadisper
mod <- betadisper(presence_absence_2016_order_DM, SampleID_scores$Order)

# Permutation Test
permutest(mod, pairwise = TRUE, permutations = 999)

# Tukey's Honest Significant Differences
(mod.HSD <- TukeyHSD(mod))
plot(mod.HSD)


## 6) Relative Abundance - Host Invertebrate Taxonomy (Genus)

# Generating Bray-Curtis Dissimilarity Matrix:
relative_abundance_2016_genus <- read.csv("/data/shawn/2016_Multivariate_Analysis/relative_abundance_2016_genus.csv")
relative_abundance_2016_genus <- as.data.frame(relative_abundance_2016_genus)
row.names(relative_abundance_2016_genus) <- relative_abundance_2016_genus$Genus
relative_abundance_2016_genus <- relative_abundance_2016_genus[,-1]
relative_abundance_2016_genus_DM <- vegdist(relative_abundance_2016_genus, method = "bray")

# Unconstrained Principal Coordinates Analysis (PCoA)
relative_abundance_2016_genus.ord <- dbrda(relative_abundance_2016_genus_DM ~ 1, add = TRUE)
summary(relative_abundance_2016_genus.ord)
PCoA_scores <- scores(relative_abundance_2016_genus.ord, display="si")
SampleID_scores <- data.frame(PCoA_scores)

# Assign Taxonomy
variables_2016 <- read.csv("/data/shawn/2016_Multivariate_Analysis/variables_2016.csv")
variables_2016 <- as.data.frame(variables_2016)
row.names(variables_2016) <- variables_2016$SampleID
variables_2016 <- variables_2016[,-1]
Genus <- variables_2016$Genus
SampleID_scores$Genus <- Genus

# Calculate Multivariate Dispersions Using Betadisper
mod <- betadisper(relative_abundance_2016_genus_DM, SampleID_scores$Genus)

# Permutation Test
permutest(mod, pairwise = TRUE, permutations = 999)

# Tukey's Honest Significant Differences
(mod.HSD <- TukeyHSD(mod))
plot(mod.HSD)


## 7) Relative Abundance - Host Invertebrate Taxonomy (Family)

# Generating Bray-Curtis Dissimilarity Matrix:
relative_abundance_2016_family <- read.csv("/data/shawn/2016_Multivariate_Analysis/relative_abundance_2016_family.csv")
relative_abundance_2016_family <- as.data.frame(relative_abundance_2016_family)
row.names(relative_abundance_2016_family) <- relative_abundance_2016_family$Family
relative_abundance_2016_family <- relative_abundance_2016_family[,-1]
relative_abundance_2016_family_DM <- vegdist(relative_abundance_2016_family, method = "bray")

# Unconstrained Principal Coordinates Analysis (PCoA)
relative_abundance_2016_family.ord <- dbrda(relative_abundance_2016_family_DM ~ 1, add = TRUE)
summary(relative_abundance_2016_family.ord)
PCoA_scores <- scores(relative_abundance_2016_family.ord, display="si")
SampleID_scores <- data.frame(PCoA_scores)

# Assign Taxonomy
variables_2016 <- read.csv("/data/shawn/2016_Multivariate_Analysis/variables_2016.csv")
variables_2016 <- as.data.frame(variables_2016)
row.names(variables_2016) <- variables_2016$SampleID
variables_2016 <- variables_2016[,-1]
Family <- variables_2016$Family
SampleID_scores$Family <- Family

# Calculate Multivariate Dispersions Using Betadisper
mod <- betadisper(relative_abundance_2016_family_DM, SampleID_scores$Family)

# Permutation Test
permutest(mod, pairwise = TRUE, permutations = 999)

# Tukey's Honest Significant Differences
(mod.HSD <- TukeyHSD(mod))
plot(mod.HSD)


## 8) Relative Abundance - Host Invertebrate Taxonomy (Order)

# Generating Bray-Curtis Dissimilarity Matrix:
relative_abundance_2016_order <- read.csv("/data/shawn/2016_Multivariate_Analysis/relative_abundance_2016_order.csv")
relative_abundance_2016_order <- as.data.frame(relative_abundance_2016_order)
row.names(relative_abundance_2016_order) <- relative_abundance_2016_order$Order
relative_abundance_2016_order <- relative_abundance_2016_order[,-1]
relative_abundance_2016_order_DM <- vegdist(relative_abundance_2016_order, method = "bray")

# Unconstrained Principal Coordinates Analysis (PCoA)
relative_abundance_2016_order.ord <- dbrda(relative_abundance_2016_order_DM ~ 1, add = TRUE)
summary(relative_abundance_2016_order.ord)
PCoA_scores <- scores(relative_abundance_2016_order.ord, display="si")
SampleID_scores <- data.frame(PCoA_scores)

# Assign Taxonomy
variables_2016 <- read.csv("/data/shawn/2016_Multivariate_Analysis/variables_2016.csv")
variables_2016 <- as.data.frame(variables_2016)
row.names(variables_2016) <- variables_2016$SampleID
variables_2016 <- variables_2016[,-1]
Order <- variables_2016$Order
SampleID_scores$Order <- Order

# Calculate Multivariate Dispersions Using Betadisper
mod <- betadisper(relative_abundance_2016_order_DM, SampleID_scores$Order)

# Permutation Test
permutest(mod, pairwise = TRUE, permutations = 999)

# Tukey's Honest Significant Differences
(mod.HSD <- TukeyHSD(mod))
plot(mod.HSD)


## 9) Presence & Absence - Discriminant Function Analysis for Functional Feeding Group

# Discriminant Function Analysis Using Jacknifed Prediction
presence_absence_2016_FFG <- read.csv("/data/shawn/2016_Multivariate_Analysis/presence_absence_2016_FFG.csv")
presence_absence_2016_FFG <- as.data.frame(presence_absence_2016_FFG)
row.names(presence_absence_2016_FFG) <- presence_absence_2016_FFG$SampleID
presence_absence_2016_FFG <- presence_absence_2016_FFG[,-1]

DFA_fit_JK <- lda(FFG ~ denovo36 + denovo122 + denovo133 + denovo140 + denovo162 + denovo216 + denovo255 + 
						denovo332 + denovo367 + denovo391 + denovo402 + denovo435 + denovo448 + denovo465 + denovo478 +
						denovo511 + denovo535 + denovo544 + denovo567 + denovo584 + denovo624 + denovo680 +
						denovo707 + denovo728 + denovo730 + denovo747 + denovo837 + denovo922 + denovo924 + denovo946 +
						denovo964 + denovo991 + denovo994 + denovo1013 + denovo1027 + denovo1041 + denovo1047 + denovo1048 +
						denovo1054 + denovo1077 + denovo1106 + denovo1118 + denovo1142 + denovo1182 + denovo1211 +
						denovo1223 + denovo1229 + denovo1292 + denovo1332 + denovo1335 + denovo1339 + denovo1365 +
						denovo1407 + denovo1408 + denovo1487 + denovo1524 + denovo1533 + denovo1540 + denovo1545 + denovo1553 +
						denovo1621 + denovo1625 + denovo1632 + denovo1642 + denovo1664 + denovo1669 + denovo1707 + denovo1767 +
						denovo1773 + denovo1836 + denovo1839 + denovo1887 + denovo1941 + denovo1944 + denovo1989 + denovo1992 +
						denovo2011 + denovo2029 + denovo2120 + denovo2147 + denovo2179 + denovo2189 + denovo2245 + denovo2269 +
						denovo2290 + denovo2297 + denovo2351 + denovo2359 + denovo2367 + denovo2418 + denovo2465 + denovo2479 +
						denovo2491 + denovo2505 + denovo2549 + denovo2554 + denovo2556 + denovo2567 + denovo2581 + denovo2606 +
						denovo2619 + denovo2677 + denovo2713 + denovo2800 + denovo2808 + denovo2824 + denovo2835 + denovo2839 +
						denovo2848 + denovo2867 + denovo2952 + denovo3030 + denovo3065 + denovo3079 + denovo3151 +
						denovo3198 + denovo3246 + denovo3319 + denovo3335 + denovo3346 + denovo3363 + denovo3368 +
						denovo3406 + denovo3426 + denovo3432 + denovo3449 + denovo3461 + denovo3480 + denovo3489 + denovo3536 +
						denovo3599 + denovo3643 + denovo3652 + denovo3654 + denovo3673 + denovo3679 + denovo3700 + denovo3745 +
						denovo3768 + denovo3777 + denovo3838 + denovo3853 + denovo3858 + denovo3882 + denovo3929 + denovo3959 +
						denovo4025 + denovo4040 + denovo4074 + denovo4079 + denovo4082 + denovo4114 + denovo4154 +
						denovo4164 + denovo4205 + denovo4228 + denovo4236 + denovo4247 + denovo4283 + denovo4289 + denovo4312 +
						denovo4319 + denovo4325 + denovo4340 + denovo4352 + denovo4378 + denovo4408 + denovo4419 + denovo4424 +
						denovo4437 + denovo4440 + denovo4499 + denovo4501 + denovo4533 + denovo4575 + denovo4609 + denovo4631 +
						denovo4636 + denovo4711 + denovo4746 + denovo4759 + denovo4762 + denovo4770 + denovo4805 + denovo4834 +
						denovo4836 + denovo4859 + denovo4875 + denovo4900 + denovo4904 + denovo4984 + denovo4988 + denovo5067 +
						denovo5133 + denovo5168 + denovo5198 + denovo5204 + denovo5226 + denovo5247 + denovo5263 +
						denovo5269 + denovo5288 + denovo5293 + denovo5324 + denovo5414 + denovo5448 + denovo5509 + denovo5550 +
						denovo5581 + denovo5647 + denovo5690 + denovo5741 + denovo5776 + denovo5801 + denovo5859 +
						denovo5971 + denovo5975 + denovo5980 + denovo5986 + denovo6019 + denovo6138 + denovo6148 +
						denovo6150 + denovo6209 + denovo6214 + denovo6226 + denovo6255 + denovo6316 + denovo6402 + denovo6411 +
						denovo6438 + denovo6439 + denovo6451 + denovo6453 + denovo6472 + denovo6473 + denovo6493 + denovo6496 +
						denovo6499 + denovo6504 + denovo6520 + denovo6521 + denovo6527 + denovo6550 + denovo6573 + denovo6577 +
						denovo6583 + denovo6608 + denovo6660 + denovo6689 + denovo6709 + denovo6735 + denovo6749 + denovo6811 +
						denovo6816 + denovo6833 + denovo6835 + denovo6856 + denovo6872 + denovo6878 + denovo6909 + denovo6918 +
						denovo7003 + denovo7005 + denovo7007 + denovo7053 + denovo7060 + denovo7074 + denovo7086 + denovo7124 +
						denovo7136 + denovo7141 + denovo7162 + denovo7171 + denovo7186 + denovo7187 + denovo7195 + denovo7197 +
						denovo7204 + denovo7206 + denovo7238 + denovo7239 + denovo7251 + denovo7257 + denovo7288 + denovo7303 +
						denovo7334 + denovo7336 + denovo7342 + denovo7347 + denovo7359 + denovo7361 + denovo7369 + denovo7414 +
						denovo7417 + denovo7432 + denovo7456 + denovo7459 + denovo7465 + denovo7470 + denovo7473 + denovo7495 +
						denovo7508 + denovo7522 + denovo7528 + denovo7541 + denovo7551 + denovo7559 + denovo7672 + denovo7682 +
						denovo7729 + denovo7741 + denovo7753 + denovo7771 + denovo7794 + denovo7800 + denovo7841 + denovo7866 +
						denovo7880 + denovo7895 + denovo7908 + denovo7939 + denovo7953 + denovo7958 + denovo7966 + denovo7995 +
						denovo8002 + denovo8017 + denovo8034 + denovo8057 + denovo8087 + denovo8121 + denovo8126 + denovo8155 +
						denovo8181 + denovo8186 + denovo8188 + denovo8230 + denovo8263 + denovo8288 + denovo8303 + denovo8350 +
						denovo8352 + denovo8357 + denovo8385 + denovo8386 + denovo8405 + denovo8428 + denovo8463 + denovo8485 +
						denovo8502 + denovo8563 + denovo8582 + denovo8599 + denovo8602 + denovo8606 + denovo8696 +
						denovo8704 + denovo8714 + denovo8735 + denovo8741 + denovo8790 + denovo8817 + denovo8861 + denovo8879 +
						denovo8889 + denovo8906 + denovo8908 + denovo8929 + denovo8970 + denovo8978 + denovo8986 + denovo8991 +
						denovo9005 + denovo9017 + denovo9030 + denovo9075 + denovo9097 + denovo9115 + denovo9118 + denovo9181 +
						denovo9204 + denovo9206 + denovo9216 + denovo9227 + denovo9248 + denovo9265 + denovo9271 + denovo9276 +
						denovo9341 + denovo9350 + denovo9366 + denovo9409 + denovo9431 + denovo9457 + denovo9464 + denovo9476 +
						denovo9501 + denovo9502 + denovo9606 + denovo9615 + denovo9636 + denovo9678 + denovo9684 +
						denovo9687 + denovo9711 + denovo9726 + denovo9758 + denovo9761 + denovo9771 + denovo9802 + denovo9809 +
						denovo9812 + denovo9900 + denovo9914 + denovo9920 + denovo9932 + denovo9953 + denovo9957 +
						denovo9967 + denovo10061 + denovo10127 + denovo10155 + denovo10157 + denovo10176 + denovo10181 +
						denovo10213 + denovo10226 + denovo10227 + denovo10241 + denovo10264 + denovo10282 + denovo10315 + denovo10329 +
						denovo10333 + denovo10352 + denovo10366 + denovo10376 + denovo10383 + denovo10395 + denovo10447 + denovo10452 +
						denovo10464 + denovo10487 + denovo10516 + denovo10559 + denovo10580 + denovo10581 + denovo10602 + denovo10618 +
						denovo10621 + denovo10634 + denovo10639 + denovo10646 + denovo10647 + denovo10648 + denovo10650 + denovo10675 +
						denovo10689 + denovo10718 + denovo10732 + denovo10752 + denovo10831 + denovo10842 + denovo10866 + denovo10889 +
						denovo10898 + denovo10992 + denovo11010 + denovo11038 + denovo11050 + denovo11051 + denovo11073 + denovo11095 +
						denovo11152 + denovo11191 + denovo11209 + denovo11229 + denovo11234 + denovo11235 + denovo11245 + denovo11247 +
						denovo11249 + denovo11286 + denovo11290 + denovo11329 + denovo11367 + denovo11403 + denovo11410 +
						denovo11419 + denovo11430 + denovo11505 + denovo11556 + denovo11574 + denovo11588 + denovo11597 + denovo11623 +
						denovo11629 + denovo11646 + denovo11689 + denovo11730 + denovo11772 + denovo11778 + denovo11784 + denovo11789 +
						denovo11814 + denovo11817 + denovo11820 + denovo11844 + denovo11866 + denovo11873 + denovo11890 + denovo11958 +
						denovo11971 + denovo11980 + denovo12022 + denovo12027 + denovo12030 + denovo12052 + denovo12062 + denovo12087 +
						denovo12137 + denovo12185 + denovo12192 + denovo12221 + denovo12242 + denovo12248 + denovo12255 + denovo12265 +
						denovo12274 + denovo12275 + denovo12407 + denovo12420 + denovo12474 + denovo12497 + denovo12548 + denovo12609 +
						denovo12621 + denovo12641 + denovo12649 + denovo12658 + denovo12690 + denovo12693 + denovo12728 + denovo12736 +
						denovo12765 + denovo12777 + denovo12783 + denovo12808 + denovo12842 + denovo12870 + denovo12873 +
						denovo12881 + denovo12890 + denovo12913 + denovo12943 + denovo12952 + denovo12960 + denovo12963 + denovo13006 +
						denovo13017 + denovo13068 + denovo13075 + denovo13081 + denovo13096 + denovo13120 + denovo13121 + denovo13148 +
						denovo13183 + denovo13197 + denovo13265 + denovo13269 + denovo13312 + denovo13334 + denovo13336 + denovo13348 +
						denovo13417 + denovo13466 + denovo13515 + denovo13517 + denovo13524 + denovo13529 + denovo13536 +
						denovo13547 + denovo13569 + denovo13641 + denovo13669 + denovo13673 + denovo13708 + denovo13710 + denovo13767 +
						denovo13781 + denovo13800 + denovo13838 + denovo13877 + denovo13900 + denovo13928 + denovo13935 + denovo13948 +
						denovo13952 + denovo13996 + denovo14012 + denovo14019 + denovo14029 + denovo14044 + denovo14089 +
						denovo14125 + denovo14153 + denovo14202 + denovo14214 + denovo14249 + denovo14250 + denovo14251 + denovo14255 +
						denovo14273 + denovo14278 + denovo14294 + denovo14308 + denovo14329 + denovo14332 + denovo14346 +
						denovo14358 + denovo14400 + denovo14411 + denovo14544 + denovo14565 + denovo14640 + denovo14642 + denovo14666 +
						denovo14679 + denovo14703 + denovo14705 + denovo14713 + denovo14781 + denovo14802 + denovo14803 + denovo14809 +
						denovo14818 + denovo14820 + denovo14821 + denovo14971 + denovo14978 + denovo14982 + denovo14991 +
						denovo15000 + denovo15015 + denovo15066 + denovo15076 + denovo15101 + denovo15150 + denovo15151 +
						denovo15185 + denovo15212 + denovo15284 + denovo15307 + denovo15356 + denovo15359 + denovo15368 + denovo15371 +
						denovo15376 + denovo15430 + denovo15438 + denovo15449 + denovo15450 + denovo15509 + denovo15558 + denovo15581 +
						denovo15601 + denovo15711 + denovo15713 + denovo15714 + denovo15729 + denovo15764 + denovo15801 + denovo15811 +
						denovo15835 + denovo15849 + denovo15855 + denovo15895 + denovo15899 + denovo15968 + denovo15998 + denovo16021 +
						denovo16058 + denovo16064 + denovo16075 + denovo16081 + denovo16177 + denovo16182 + denovo16194 + denovo16353 +
						denovo16407 + denovo16437 + denovo16466 + denovo16470 + denovo16504 + denovo16516 + denovo16524 + denovo16528 +
						denovo16552 + denovo16577 + denovo16584 + denovo16603 + denovo16678 + denovo16722 + denovo16737 + denovo16746 +
						denovo16763 + denovo16782 + denovo16812 + denovo16819 + denovo16828 + denovo16857 + denovo16859 + denovo16861 +
						denovo16914 + denovo16915 + denovo16918 + denovo16923 + denovo16979 + denovo17022 + denovo17052 + denovo17078 +
						denovo17092 + denovo17190 + denovo17194 + denovo17232 + denovo17277 + denovo17292 + denovo17303 + denovo17304 +
						denovo17324 + denovo17336 + denovo17365 + denovo17373 + denovo17402 + denovo17436 + denovo17515 + denovo17562 +
						denovo17628 + denovo17633 + denovo17634 + denovo17645 + denovo17646 + denovo17683 + denovo17691 + denovo17713 +
						denovo17735 + denovo17762 + denovo17769 + denovo17836 + denovo17846 + denovo17856 + denovo17859 + denovo17874 +
						denovo17891 + denovo17940 + denovo17955 + denovo17958 + denovo17967 + denovo17969 + denovo17990 + denovo17991 +
						denovo18014 + denovo18038 + denovo18090 + denovo18093 + denovo18099 + denovo18105 + denovo18163 + denovo18176 +
						denovo18264 + denovo18296 + denovo18301 + denovo18302 + denovo18355 + denovo18368 + denovo18375 + denovo18416 +
						denovo18418 + denovo18452 + denovo18488 + denovo18498 + denovo18534 + denovo18539 + denovo18588 + denovo18589 +
						denovo18633 + denovo18674 + denovo18680 + denovo18716 + denovo18739 + denovo18760 + denovo18821 + denovo18848 +
						denovo18874 + denovo18906 + denovo18941 + denovo18975 + denovo18980 + denovo18994 + denovo19077 + denovo19096 +
						denovo19124 + denovo19146 + denovo19218 + denovo19265 + denovo19290 + denovo19392 + denovo19400 + denovo19416 +
						denovo19450 + denovo19499 + denovo19517 + denovo19560 + denovo19576 + denovo19619 + denovo19650 + denovo19664 +
						denovo19667 + denovo19715 + denovo19732 + denovo19773 + denovo19775 + denovo19829 + denovo19848 + denovo19892 +
						denovo19909 + denovo20052 + denovo20065 + denovo20075 + denovo20111 + denovo20165 + denovo20196 + denovo20215 +
						denovo20226 + denovo20234 + denovo20275 + denovo20282 + denovo20290 + denovo20345 + denovo20361 + denovo20371 +
						denovo20372 + denovo20374 + denovo20403 + denovo20451 + denovo20462 + denovo20466 + denovo20491 + denovo20522 +
						denovo20530 + denovo20586 + denovo20606 + denovo20682 + denovo20687 + denovo20693 + denovo20701 + denovo20741 +
						denovo20808 + denovo20813 + denovo20904 + denovo20940 + denovo20997 + denovo21001 + denovo21010 + denovo21049 +
						denovo21064 + denovo21071 + denovo21075 + denovo21085 + denovo21094 + denovo21096 + denovo21098 + denovo21150 +
						denovo21155 + denovo21169 + denovo21173 + denovo21178 + denovo21192 + denovo21200 + denovo21348 + denovo21353 +
						denovo21378 + denovo21389 + denovo21413 + denovo21425 + denovo21471 + denovo21491 + denovo21499 + denovo21525 +
						denovo21543 + denovo21568 + denovo21580 + denovo21587 + denovo21627 + denovo21634 + denovo21638 + denovo21642 +
						denovo21659 + denovo21700 + denovo21726 + denovo21739 + denovo21762 + denovo21767 + denovo21793 + denovo21843 +
						denovo21858 + denovo21874 + denovo21891 + denovo21917 + denovo21924 + denovo21953 + denovo21989 + denovo22004 +
						denovo22022 + denovo22029 + denovo22031 + denovo22065 + denovo22070 + denovo22110 + denovo22137 + denovo22139 +
						denovo22171 + denovo22205 + denovo22214 + denovo22233 + denovo22294 + denovo22309 + denovo22314 +
						denovo22330 + denovo22338 + denovo22344 + denovo22365 + denovo22423 + denovo22503 + denovo22539 + denovo22553 +
						denovo22561 + denovo22596 + denovo22629 + denovo22695 + denovo22721 + denovo22821 + denovo22850 + denovo22879 +
						denovo22881 + denovo22923 + denovo22934 + denovo22941 + denovo22972 + denovo22974 + denovo22975 + denovo23006 +
						denovo23064 + denovo23106 + denovo23118 + denovo23128 + denovo23143 + denovo23155 + denovo23161 + denovo23173 +
						denovo23184 + denovo23191 + denovo23192 + denovo23245 + denovo23266 + denovo23274 + denovo23280 + denovo23338 +
						denovo23345 + denovo23374 + denovo23377 + denovo23391 + denovo23424 + denovo23453 + denovo23509 +
						denovo23544 + denovo23552 + denovo23560 + denovo23562 + denovo23591 + denovo23625 + denovo23691 + denovo23702 +
						denovo23712 + denovo23734 + denovo23741 + denovo23769 + denovo23798 + denovo23816 + denovo23833 + denovo23893 +
						denovo23912 + denovo23954 + denovo23971 + denovo23986 + denovo23992 + denovo23999 + denovo24032 + denovo24045 +
						denovo24064 + denovo24080 + denovo24082 + denovo24117 + denovo24189 + denovo24262 + denovo24275 + denovo24338 +
						denovo24390 + denovo24397 + denovo24433 + denovo24439 + denovo24533 + denovo24547 + denovo24552 + denovo24608 +
						denovo24641 + denovo24647 + denovo24660 + denovo24668 + denovo24742 + denovo24751 + denovo24757 + denovo24758 +
						denovo24782 + denovo24869 + denovo24872 + denovo24881 + denovo24900 + denovo24939 + denovo24958 + denovo25001 +
						denovo25027 + denovo25105 + denovo25137 + denovo25154 + denovo25205 + denovo25228 + denovo25229 + denovo25247 +
						denovo25261 + denovo25280 + denovo25284 + denovo25306 + denovo25347 + denovo25367 + denovo25388 + denovo25393 +
						denovo25427 + denovo25429 + denovo25464 + denovo25493 + denovo25582 + denovo25590 + denovo25594 + denovo25662 +
						denovo25699 + denovo25700 + denovo25726 + denovo25761 + denovo25797 + denovo25890 + denovo25903 + denovo25927 +
						denovo25990 + denovo25995 + denovo26040 + denovo26076 + denovo26099 + denovo26102 + denovo26135 +
						denovo26184 + denovo26220 + denovo26268 + denovo26278 + denovo26310 + denovo26361 + denovo26424 + denovo26453 +
						denovo26471 + denovo26474 + denovo26511 + denovo26541 + denovo26561 + denovo26585 + denovo26591 + denovo26602 +
						denovo26610 + denovo26633 + denovo26683 + denovo26719 + denovo26752 + denovo26785 + denovo26893 + denovo26902 +
						denovo26927 + denovo26954 + denovo26979 + denovo27015 + denovo27024 + denovo27039 + denovo27045 + denovo27054 +
						denovo27066 + denovo27126 + denovo27139 + denovo27170 + denovo27218 + denovo27221 + denovo27254 + denovo27289 +
						denovo27354 + denovo27355 + denovo27412 + denovo27461 + denovo27504 + denovo27534 + denovo27535 + denovo27603 +
						denovo27648 + denovo27673 + denovo27682 + denovo27700 + denovo27702 + denovo27725 + denovo27869 + denovo27874 +
						denovo27878 + denovo27900 + denovo27928 + denovo27946 + denovo27968 + denovo28021 + denovo28113 + denovo28128 +
						denovo28131 + denovo28253 + denovo28266 + denovo28333 + denovo28387 + denovo28392 + denovo28451 + denovo28485 +
						denovo28538 + denovo28595 + denovo28605 + denovo28613 + denovo28614 + denovo28625 + denovo28633 +
						denovo28699 + denovo28732 + denovo28750 + denovo28787 + denovo28790 + denovo28797 + denovo28838 +
						denovo28848 + denovo28885 + denovo28979 + denovo29009 + denovo29016 + denovo29038 + denovo29062 +
						denovo29083 + denovo29163 + denovo29188 + denovo29294 + denovo29310 + denovo29327 + denovo29355 +
						denovo29362 + denovo29381 + denovo29384 + denovo29388 + denovo29399 + denovo29420 + denovo29437 + denovo29438 +
						denovo29442 + denovo29459 + denovo29567 + denovo29620 + denovo29639 + denovo29642 + denovo29724 + denovo29754 +
						denovo29762 + denovo29793 + denovo29812 + denovo29821 + denovo29841 + denovo29842 + denovo29844 + denovo29920 +
						denovo29976 + denovo30013 + denovo30050 + denovo30054 + denovo30058 + denovo30061 + denovo30083 +
						denovo30112 + denovo30146 + denovo30157 + denovo30167 + denovo30186 + denovo30218 + denovo30224 + denovo30283 +
						denovo30317 + denovo30366 + denovo30371 + denovo30383 + denovo30393 + denovo30420 + denovo30484 +
						denovo30573 + denovo30580 + denovo30592 + denovo30635 + denovo30693 + denovo30715 + denovo30728 + denovo30734 +
						denovo30756 + denovo30777 + denovo30782 + denovo30811 + denovo30875 + denovo30883 + denovo30919 + denovo30940 +
						denovo30941 + denovo31003 + denovo31036 + denovo31052 + denovo31070 + denovo31132 + denovo31177 + denovo31274 +
						denovo31322 + denovo31348 + denovo31376 + denovo31385 + denovo31403 + denovo31414 + denovo31420 + denovo31430 +
						denovo31449 + denovo31460 + denovo31476 + denovo31477 + denovo31486 + denovo31487 + denovo31506 + denovo31534 +
						denovo31546 + denovo31557 + denovo31574 + denovo31576 + denovo31669 + denovo31673 + denovo31688 + denovo31746 +
						denovo31773 + denovo31820 + denovo31843 + denovo31846 + denovo31904 + denovo31934 + denovo31944 + denovo31959 +
						denovo31968 + denovo31980 + denovo31993 + denovo31999 + denovo32001 + denovo32013 + denovo32049 + denovo32058 +
						denovo32062 + denovo32079 + denovo32084 + denovo32094 + denovo32099 + denovo32199 + denovo32203 + denovo32266 +
						denovo32269 + denovo32278 + denovo32304 + denovo32477 + denovo32498 + denovo32547 + denovo32564 + denovo32569 +
						denovo32601 + denovo32656 + denovo32679 + denovo32738 + denovo32824 + denovo32855 + denovo32901 + denovo32956 +
						denovo33007 + denovo33021 + denovo33022 + denovo33074 + denovo33136 + denovo33171 + denovo33245 +
						denovo33246 + denovo33316 + denovo33340 + denovo33347 + denovo33360 + denovo33376 + denovo33440 +
						denovo33472 + denovo33491 + denovo33512 + denovo33518 + denovo33524 + denovo33571 + denovo33607 +
						denovo33616 + denovo33632 + denovo33644 + denovo33702 + denovo33729 + denovo33733 + denovo33757 + denovo33778 +
						denovo33793 + denovo33800 + denovo33831 + denovo33834 + denovo33844 + denovo33877 + denovo33891 + denovo33925 +
						denovo33941 + denovo34003 + denovo34012 + denovo34019 + denovo34040 + denovo34041 + denovo34049 + denovo34061 +
						denovo34111 + denovo34149 + denovo34244 + denovo34260 + denovo34337 + denovo34351 + denovo34364 + denovo34369 +
						denovo34407 + denovo34455 + denovo34471 + denovo34492 + denovo34499 + denovo34513 + denovo34525 + denovo34537 +
						denovo34585 + denovo34616 + denovo34636 + denovo34680 + denovo34713 + denovo34718 + denovo34725 + denovo34751 +
						denovo34763 + denovo34776 + denovo34808 + denovo34816 + denovo34840 + denovo34852 + denovo34873 + denovo34881 +
						denovo34918 + denovo34937 + denovo34998 + denovo35012 + denovo35022 + denovo35036 + denovo35062 + denovo35083 +
						denovo35089 + denovo35121 + denovo35124 + denovo35129 + denovo35133 + denovo35140 + denovo35153 + denovo35174 +
						denovo35196 + denovo35207 + denovo35216 + denovo35227 + denovo35246 + denovo35348 + denovo35369 + denovo35381 +
						denovo35415 + denovo35426 + denovo35451 + denovo35470 + denovo35487 + denovo35540 + denovo35585 + denovo35609 +
						denovo35655 + denovo35689 + denovo35729 + denovo35734 + denovo35756 + denovo35763 + denovo35790 + denovo35861 +
						denovo35862 + denovo35910 + denovo35917 + denovo35969 + denovo36019 + denovo36051 + denovo36059 +
						denovo36156 + denovo36166 + denovo36224 + denovo36231 + denovo36273 + denovo36279 + denovo36295 + denovo36392 +
						denovo36397 + denovo36401 + denovo36409 + denovo36459 + denovo36484 + denovo36503 + denovo36505 + denovo36519 +
						denovo36522 + denovo36539 + denovo36541 + denovo36545 + denovo36556 + denovo36565 + denovo36602 + denovo36603 +
						denovo36617 + denovo36639 + denovo36661 + denovo36671 + denovo36676 + denovo36680 + denovo36684 + denovo36725 +
						denovo36788 + denovo36805 + denovo36855 + denovo36892 + denovo36939 + denovo37007 + denovo37009 + denovo37046 +
						denovo37143 + denovo37144 + denovo37176 + denovo37230 + denovo37232 + denovo37281 + denovo37289 + denovo37328 +
						denovo37335 + denovo37352 + denovo37354 + denovo37357 + denovo37361 + denovo37374 + denovo37394 + denovo37483 +
						denovo37554 + denovo37581 + denovo37584 + denovo37596 + denovo37605 + denovo37637 + denovo37656 + denovo37659 +
						denovo37680 + denovo37683 + denovo37685 + denovo37717 + denovo37740 + denovo37764 + denovo37786 +
						denovo37799 + denovo37800 + denovo37811 + denovo37819 + denovo37831 + denovo37850 + denovo37896 + denovo37933 +
						denovo37978 + denovo37983 + denovo37996 + denovo38175 + denovo38177 + denovo38192 + denovo38195 +
						denovo38198 + denovo38251 + denovo38259 + denovo38279 + denovo38297 + denovo38301 + denovo38303 + denovo38319 +
						denovo38347 + denovo38364 + denovo38434 + denovo38460 + denovo38479 + denovo38500 + denovo38517 + denovo38523 +
						denovo38524 + denovo38561 + denovo38576 + denovo38634 + denovo38640 + denovo38656 + denovo38692 + denovo38695 +
						denovo38726 + denovo38732 + denovo38737 + denovo38762 + denovo38815 + denovo38826 + denovo38868 + denovo38905 +
						denovo39013 + denovo39032 + denovo39039 + denovo39049 + denovo39090 + denovo39140 + denovo39143 + denovo39151 +
						denovo39164 + denovo39182 + denovo39216 + denovo39243 + denovo39281 + denovo39304 + denovo39307 + denovo39344 +
						denovo39432 + denovo39435 + denovo39464 + denovo39476 + denovo39485 + denovo39507 + denovo39569 + denovo39601 +
						denovo39623 + denovo39649 + denovo39674 + denovo39720 + denovo39729 + denovo39798 + denovo39817 + denovo39942 +
						denovo39955 + denovo39996 + denovo40008 + denovo40057 + denovo40065 + denovo40067 + denovo40071 + denovo40078 +
						denovo40121 + denovo40125 + denovo40136 + denovo40138 + denovo40143 + denovo40154 + denovo40171 + denovo40174 +
						denovo40212 + denovo40246 + denovo40248 + denovo40257 + denovo40290 + denovo40350 + denovo40395 +
						denovo40403 + denovo40437 + denovo40479 + denovo40480 + denovo40510 + denovo40526 + denovo40569 + denovo40570 +
						denovo40608 + denovo40637 + denovo40728 + denovo40744 + denovo40755 + denovo40774 + denovo40780 +
						denovo40843 + denovo40910 + denovo40932 + denovo40971 + denovo41003 + denovo41027 + denovo41029 +
						denovo41067 + denovo41108 + denovo41114 + denovo41164 + denovo41166 + denovo41169 + denovo41172 + denovo41182 +
						denovo41235 + denovo41238 + denovo41240 + denovo41287 + denovo41299 + denovo41306 + denovo41325 + denovo41363 +
						denovo41381 + denovo41393 + denovo41425 + denovo41433 + denovo41499 + denovo41548 + denovo41553 + denovo41626 +
						denovo41637 + denovo41671 + denovo41707 + denovo41719 + denovo41738 + denovo41751 + denovo41806 +
						denovo41845 + denovo41868 + denovo41888 + denovo41911 + denovo42002 + denovo42005 + denovo42058 + denovo42157 +
						denovo42159 + denovo42193 + denovo42253 + denovo42316 + denovo42327 + denovo42340 + denovo42364 + denovo42381 +
						denovo42405 + denovo42419 + denovo42436 + denovo42444 + denovo42456 + denovo42485 + denovo42521 + denovo42538 +
						denovo42592 + denovo42632 + denovo42711 + denovo42716 + denovo42905 + denovo42931 + denovo42936 + denovo42944 +
						denovo42957 + denovo42976 + denovo42979 + denovo42993 + denovo42997 + denovo43004 + denovo43005 + denovo43018 +
						denovo43065 + denovo43106 + denovo43124 + denovo43133 + denovo43190 + denovo43231 + denovo43232 + denovo43258 +
						denovo43259 + denovo43280 + denovo43304 + denovo43350 + denovo43357 + denovo43373 + denovo43425 + denovo43475 +
						denovo43561 + denovo43593 + denovo43601 + denovo43671 + denovo43711 + denovo43719 + denovo43806 + denovo43812 +
						denovo43854 + denovo43857 + denovo43970 + denovo43977 + denovo43984 + denovo43997 + denovo44016 + denovo44034 +
						denovo44038 + denovo44077 + denovo44079 + denovo44093 + denovo44095 + denovo44146 + denovo44163 +
						denovo44234 + denovo44260 + denovo44276 + denovo44289 + denovo44301 + denovo44340 + denovo44358 + denovo44414 +
						denovo44496 + denovo44499 + denovo44502 + denovo44510 + denovo44512 + denovo44527 + denovo44528 + denovo44532 +
						denovo44538 + denovo44545 + denovo44609 + denovo44662 + denovo44669 + denovo44678 + denovo44696 + denovo44702 +
						denovo44708 + denovo44715 + denovo44726 + denovo44745 + denovo44793 + denovo44816 + denovo44863 + denovo44867 +
						denovo44912 + denovo44940 + denovo44979 + denovo44994 + denovo45004 + denovo45013 + denovo45018 + denovo45024 +
						denovo45046 + denovo45063 + denovo45097 + denovo45108 + denovo45110 + denovo45140 + denovo45153 + denovo45175 +
						denovo45194 + denovo45195 + denovo45202 + denovo45230 + denovo45242 + denovo45298 + denovo45302 + denovo45316 +
						denovo45322 + denovo45419 + denovo45422 + denovo45437 + denovo45453 + denovo45463 + denovo45476 + denovo45516 +
						denovo45590 + denovo45637 + denovo45642 + denovo45680 + denovo45690 + denovo45701 + denovo45702 + denovo45707 +
						denovo45708 + denovo45770 + denovo45783 + denovo45876 + denovo45885 + denovo45898 + denovo45923 + denovo45928 +
						denovo45939 + denovo45967 + denovo45971 + denovo46052 + denovo46097 + denovo46109 + denovo46141 +
						denovo46165 + denovo46188 + denovo46203 + denovo46226 + denovo46235 + denovo46263 + denovo46268 + denovo46293 +
						denovo46322 + denovo46336 + denovo46344 + denovo46346 + denovo46373 + denovo46466 + denovo46472 + denovo46516 +
						denovo46542 + denovo46555 + denovo46579 + denovo46598 + denovo46654 + denovo46655 + denovo46656 + denovo46687 +
						denovo46702 + denovo46721 + denovo46816 + denovo46836 + denovo46870 + denovo46874 + denovo46879 + denovo46883 +
						denovo46906 + denovo46926 + denovo46934 + denovo46966 + denovo46984 + denovo47032 + denovo47041 + denovo47049 +
						denovo47052 + denovo47066 + denovo47109 + denovo47123 + denovo47259 + denovo47312 + denovo47382 + denovo47441 +
						denovo47510 + denovo47528 + denovo47557 + denovo47587 + denovo47595 + denovo47627 + denovo47661 + denovo47667 +
						denovo47682 + denovo47717 + denovo47731 + denovo47801 + denovo47803 + denovo47808 + denovo47816 + denovo47824 +
						denovo47826 + denovo47852 + denovo47880 + denovo47891 + denovo47919 + denovo47944 + denovo48005 + denovo48013 +
						denovo48031 + denovo48087 + denovo48101 + denovo48163 + denovo48177 + denovo48209 + denovo48210 + denovo48288 +
						denovo48298 + denovo48323 + denovo48353 + denovo48384 + denovo48434 + denovo48442 + denovo48449 + denovo48458 +
						denovo48488 + denovo48499 + denovo48520 + denovo48531 + denovo48543 + denovo48546 + denovo48554 + denovo48579 +
						denovo48582 + denovo48616 + denovo48660 + denovo48676 + denovo48708 + denovo48726 + denovo48747 + denovo48764 +
						denovo48766 + denovo48816 + denovo48846 + denovo48861 + denovo48863 + denovo48910 + denovo48915 + denovo48972 +
						denovo49012 + denovo49063 + denovo49070 + denovo49084 + denovo49086 + denovo49107 + denovo49111 + denovo49122 +
						denovo49127 + denovo49157 + denovo49176 + denovo49270 + denovo49280 + denovo49298 + denovo49364 + denovo49422 +
						denovo49447 + denovo49460 + denovo49518 + denovo49541 + denovo49572 + denovo49594 + denovo49602 +
						denovo49628 + denovo49679 + denovo49698 + denovo49735 + denovo49768 + denovo49782 + denovo49789 + denovo49803 +
						denovo49847 + denovo49872 + denovo49914 + denovo49942 + denovo49975 + denovo50005 + denovo50043 + denovo50053 +
						denovo50125 + denovo50128 + denovo50166 + denovo50185 + denovo50189 + denovo50210 + denovo50218 + denovo50222 +
						denovo50229 + denovo50243 + denovo50263 + denovo50284 + denovo50333 + denovo50368 + denovo50400 + denovo50434 +
						denovo50458 + denovo50542 + denovo50561 + denovo50626 + denovo50628 + denovo50698 + denovo50722 +
						denovo50727 + denovo50740 + denovo50741 + denovo50771 + denovo50802 + denovo50815 + denovo50818 + denovo50832 +
						denovo50880 + denovo50893 + denovo50894 + denovo50921 + denovo50951 + denovo50974 + denovo50985 +
						denovo51035 + denovo51045 + denovo51088 + denovo51117 + denovo51123 + denovo51153 + denovo51156 + denovo51159 +
						denovo51197 + denovo51213 + denovo51244 + denovo51273 + denovo51288 + denovo51323 + denovo51376 + denovo51405 +
						denovo51433 + denovo51477 + denovo51537 + denovo51540 + denovo51549 + denovo51562 + denovo51565 + denovo51599 +
						denovo51645 + denovo51650 + denovo51701 + denovo51721 + denovo51745 + denovo51756 + denovo51785 + denovo51787 +
						denovo51796 + denovo51826 + denovo51830 + denovo51843 + denovo51853 + denovo51878 + denovo51921 + denovo52036 +
						denovo52074 + denovo52146 + denovo52155 + denovo52226 + denovo52251 + denovo52256 + denovo52282 + denovo52295 +
						denovo52337 + denovo52345 + denovo52348 + denovo52363 + denovo52377 + denovo52433 + denovo52488 + denovo52528 +
						denovo52538 + denovo52542 + denovo52552 + denovo52562 + denovo52601 + denovo52605 + denovo52630 + denovo52631 +
						denovo52758 + denovo52780 + denovo52784 + denovo52836 + denovo52847 + denovo52851 + denovo52880 + denovo52891 +
						denovo52902 + denovo52934 + denovo52935 + denovo52967 + denovo53013 + denovo53024 + denovo53064 +
						denovo53073 + denovo53103 + denovo53110 + denovo53122 + denovo53139 + denovo53149 + denovo53216 + denovo53236 +
						denovo53250 + denovo53269 + denovo53278 + denovo53327 + denovo53332 + denovo53346 + denovo53356 +
						denovo53360 + denovo53364 + denovo53401 + denovo53439 + denovo53451 + denovo53468 + denovo53504 + denovo53550 +
						denovo53564 + denovo53626 + denovo53643 + denovo53665 + denovo53671 + denovo53772 + denovo53821 + denovo53826 +
						denovo53835 + denovo53852 + denovo53866 + denovo53913 + denovo53965 + denovo53967 + denovo54032 + denovo54068 +
						denovo54077 + denovo54159 + denovo54176 + denovo54181 + denovo54234 + denovo54251 + denovo54252 + denovo54266 +
						denovo54269 + denovo54296 + denovo54311 + denovo54321 + denovo54361 + denovo54433 + denovo54434 + denovo54440 +
						denovo54484 + denovo54527 + denovo54532 + denovo54582 + denovo54598 + denovo54605 + denovo54617 + denovo54635 +
						denovo54651 + denovo54652 + denovo54670 + denovo54697 + denovo54727 + denovo54792 + denovo54803 + denovo54806 +
						denovo54875 + denovo54883 + denovo54886 + denovo54900 + denovo54911 + denovo54917 + denovo54932 + denovo54998 +
						denovo55042 + denovo55134 + denovo55184 + denovo55189 + denovo55202 + denovo55260 + denovo55267 +
						denovo55269 + denovo55274 + denovo55316 + denovo55328 + denovo55341 + denovo55355 + denovo55371 + denovo55403 +
						denovo55425 + denovo55429 + denovo55498 + denovo55505 + denovo55522 + denovo55583 +
						denovo55607 + denovo55609 + denovo55767 + denovo55832 + denovo55845 + denovo55853 + denovo55891 + denovo55896 +
						denovo55930 + denovo55945 + denovo55952 + denovo55967 + denovo55995 + denovo56011 + denovo56050 + denovo56056 +
						denovo56087 + denovo56120 + denovo56147 + denovo56155 + denovo56156 + denovo56160 + denovo56174 + denovo56178 +
						denovo56215 + denovo56233 + denovo56251 + denovo56265 + denovo56344 + denovo56419 + denovo56444 + denovo56455 +
						denovo56457 + denovo56464 + denovo56467 + denovo56472 + denovo56492 + denovo56502 + denovo56504 + denovo56505 +
						denovo56521 + denovo56569 + denovo56587 + denovo56635 + denovo56639 + denovo56671 + denovo56709 + denovo56761 +
						denovo56763 + denovo56772 + denovo56778 + denovo56803 + denovo56830 + denovo56835 + denovo56853 + denovo56871 +
						denovo56894 + denovo56914 + denovo56920 + denovo56982 + denovo57017 + denovo57045 + denovo57138 +
						denovo57163 + denovo57166 + denovo57185 + denovo57207 + denovo57275 + denovo57300 + denovo57305 + denovo57353 +
						denovo57372 + denovo57421 + denovo57424 + denovo57437 + denovo57459 + denovo57476 + denovo57505 + denovo57514 +
						denovo57534 + denovo57570 + denovo57615 + denovo57699 + denovo57713 + denovo57759 + denovo57771 + denovo57838 +
						denovo57844 + denovo57868 + denovo57871 + denovo57875 + denovo57911 + denovo57914 + denovo57922 + denovo57972 +
						denovo58022 + denovo58034 + denovo58061 + denovo58083 + denovo58090 + denovo58111 + denovo58171 + denovo58245 +
						denovo58262 + denovo58328 + denovo58353 + denovo58378 + denovo58390 + denovo58520 + denovo58550 + denovo58552 +
						denovo58582 + denovo58600 + denovo58604 + denovo58635 + denovo58678 + denovo58681 + denovo58699 + denovo58709 +
						denovo58713 + denovo58717 + denovo58752 + denovo58784 + denovo58808 + denovo58829 + denovo58830 + denovo58848 +
						denovo58853 + denovo58856 + denovo58885 + denovo58948 + denovo58950 + denovo58976 + denovo58996 + denovo59030 +
						denovo59037 + denovo59060 + denovo59106 + denovo59209 + denovo59213 + denovo59238 + denovo59262 + denovo59276 +
						denovo59291 + denovo59294 + denovo59296 + denovo59302 + denovo59310 + denovo59335 + denovo59354 + denovo59370 +
						denovo59371 + denovo59372 + denovo59384 + denovo59392 + denovo59406 + denovo59415 + denovo59416 + denovo59430 +
						denovo59433 + denovo59437 + denovo59451 + denovo59460 + denovo59575 + denovo59578 + denovo59643 + denovo59658 +
						denovo59732 + denovo59744 + denovo59754 + denovo59755 + denovo59779 + denovo59817 + denovo59839 + denovo59844 +
						denovo59846 + denovo59872 + denovo59930 + denovo59960 + denovo59970 + denovo60009 + denovo60012 + denovo60068 +
						denovo60081 + denovo60094 + denovo60163 + denovo60199 + denovo60243 + denovo60309 + denovo60332 + denovo60358 +
						denovo60363 + denovo60379 + denovo60410 + denovo60469 + denovo60532 + denovo60574 + denovo60593 + denovo60650 +
						denovo60654 + denovo60692 + denovo60724 + denovo60739 + denovo60775 + denovo60794 + denovo60831 +
						denovo60878 + denovo60912 + denovo60917 + denovo60924 + denovo60946 + denovo60996 + denovo61066 + denovo61073 +
						denovo61094 + denovo61142 + denovo61145 + denovo61168 + denovo61218 + denovo61222 + denovo61249 + denovo61354 +
						denovo61369 + denovo61376 + denovo61378 + denovo61396 + denovo61402 + denovo61480 + denovo61498 + denovo61502 +
						denovo61535 + denovo61540 + denovo61563 + denovo61586 + denovo61635 + denovo61640 + denovo61657 + denovo61702 +
						denovo61722 + denovo61728 + denovo61733 + denovo61752 + denovo61787 + denovo61805 + denovo61820 + denovo61841 +
						denovo61842 + denovo61849 + denovo61912 + denovo61974 + denovo62011 + denovo62056 + denovo62127 + denovo62153 +
						denovo62166 + denovo62178 + denovo62238 + denovo62253 + denovo62282 + denovo62284 +
						denovo62303 + denovo62380 + denovo62433 + denovo62506 + denovo62540 + denovo62578 + denovo62579 + denovo62590 +
						denovo62604 + denovo62650 + denovo62670 + denovo62737 + denovo62783 + denovo62804 + denovo62883 + denovo62905 +
						denovo62943 + denovo62948 + denovo62961 + denovo63018 + denovo63037 + denovo63057 + denovo63082 + denovo63101 +
						denovo63126 + denovo63143 + denovo63152 + denovo63216 + denovo63228 + denovo63340 + denovo63341 + denovo63389 +
						denovo63415 + denovo63429 + denovo63434 + denovo63444 + denovo63482 + denovo63522 + denovo63539 + denovo63552 +
						denovo63560 + denovo63605 + denovo63698 + denovo63738 + denovo63757 + denovo63777 + denovo63857 +
						denovo63861 + denovo63863 + denovo63890 + denovo63899 + denovo63907 + denovo63913 + denovo63925 + denovo63973 +
						denovo63988 + denovo64010 + denovo64018 + denovo64026 + denovo64067 + denovo64081 + denovo64089 + denovo64093 +
						denovo64094 + denovo64121 + denovo64122 + denovo64131 + denovo64140 + denovo64141 + denovo64150 + denovo64208 +
						denovo64236 + denovo64237 + denovo64244 + denovo64260 + denovo64290 + denovo64309 + denovo64324 +
						denovo64359 + denovo64381 + denovo64430 + denovo64431 + denovo64455 + denovo64570 +
						denovo64608 + denovo64618 + denovo64631 + denovo64645 + denovo64747 + denovo64768 + denovo64769 + denovo64797 +
						denovo64804 + denovo64806 + denovo64839 + denovo64844 + denovo64852 + denovo64891 + denovo64906 + denovo64963 +
						denovo64989 + denovo65009 + denovo65010 + denovo65040 + denovo65065 + denovo65072 + denovo65088 + denovo65096 +
						denovo65113 + denovo65127 + denovo65151 + denovo65209 + denovo65233 + denovo65241 + denovo65258 + denovo65285 +
						denovo65327 + denovo65328 + denovo65338 + denovo65341 + denovo65354 + denovo65371 + denovo65398 + denovo65424 +
						denovo65476 + denovo65477 + denovo65563 + denovo65566 + denovo65638 + denovo65660 + denovo65708 + denovo65753 +
						denovo65759 + denovo65795 + denovo65830 + denovo65848 + denovo65855 + denovo65870 + denovo65873 + denovo65901 +
						denovo66008 + denovo66081 + denovo66089 + denovo66094 + denovo66129 + denovo66136 + denovo66138 + denovo66177 +
						denovo66227 + denovo66234 + denovo66236 + denovo66318 + denovo66356 + denovo66372 + denovo66380 +
						denovo66381 + denovo66418 + denovo66424 + denovo66457 + denovo66467 + denovo66483 + denovo66521 + denovo66524 +
						denovo66566 + denovo66615 + denovo66677 + denovo66684 + denovo66689 + denovo66690 + denovo66691 + denovo66698 +
						denovo66732 + denovo66741 + denovo66744 + denovo66746 + denovo66749 + denovo66757 + denovo66765 + denovo66776 +
						denovo66805 + denovo66806 + denovo66836 + denovo66876 + denovo66881 + denovo66904 + denovo66906 +
						denovo66947 + denovo66982 + denovo66987 + denovo66991 + denovo66996 + denovo67022 + denovo67074 + denovo67086 +
						denovo67112 + denovo67138 + denovo67157 + denovo67163 + denovo67181 + denovo67191 + denovo67266 + denovo67297 +
						denovo67318 + denovo67322 + denovo67378 + denovo67406 + denovo67418 + denovo67450 + denovo67470 +
						denovo67489 + denovo67492 + denovo67500 + denovo67581 + denovo67583 + denovo67599 + denovo67612 + denovo67657 +
						denovo67660 + denovo67662 + denovo67696 + denovo67725 + denovo67728 + denovo67729 + denovo67762 + denovo67769 +
						denovo67770 + denovo67807 + denovo67834 + denovo67871 + denovo68056 + denovo68058 + denovo68107 +
						denovo68129 + denovo68136 + denovo68153 + denovo68160 + denovo68170 + denovo68181 + denovo68220 + denovo68227 +
						denovo68256 + denovo68285 + denovo68301 + denovo68326 + denovo68332 + denovo68397 + denovo68408 + denovo68409 +
						denovo68422 + denovo68432 + denovo68458 + denovo68488 + denovo68494 + denovo68544 + denovo68556 + denovo68577 +
						denovo68598 + denovo68640 + denovo68646 + denovo68709 + denovo68724 + denovo68735 + denovo68767 + denovo68791 +
						denovo68834 + denovo68839 + denovo68867 + denovo68889 + denovo68914 + denovo68955 + denovo68962 + denovo68980 +
						denovo69002 + denovo69004 + denovo69018 + denovo69030 + denovo69112 + denovo69113 + denovo69145 + denovo69181 +
						denovo69213 + denovo69301 + denovo69339 + denovo69380 + denovo69414 + denovo69441 + denovo69444 +
						denovo69451 + denovo69453 + denovo69455 + denovo69474 + denovo69526 + denovo69559 + denovo69576 +
						denovo69602 + denovo69617 + denovo69620 + denovo69623 + denovo69626 + denovo69665 + denovo69705 +
						denovo69711 + denovo69717 + denovo69730 + denovo69742 + denovo69770 + denovo69772 + denovo69773 + denovo69781 +
						denovo69839 + denovo69845 + denovo69857 + denovo69866 + denovo69909 + denovo69945 + denovo69972 + denovo69980 +
						denovo70028 + denovo70083 + denovo70094 + denovo70110 + denovo70238 + denovo70241 + denovo70253 + denovo70260 +
						denovo70264 + denovo70290 + denovo70326 + denovo70346 + denovo70405 + denovo70420 + denovo70434 + denovo70491 +
						denovo70519 + denovo70615 + denovo70651 + denovo70654 + denovo70671 + denovo70684 + denovo70701 + denovo70782 +
						denovo70806 + denovo70816 + denovo70845 + denovo70854 + denovo70857 + denovo70861 + denovo70900 + denovo70906 +
						denovo70946 + denovo70950 + denovo70968 + denovo70990 + denovo71006 + denovo71008 + denovo71014 + denovo71030 +
						denovo71050 + denovo71102 + denovo71142 + denovo71236 + denovo71278 + denovo71315 + denovo71333 + denovo71348 +
						denovo71359 + denovo71366 + denovo71372 + denovo71403 + denovo71426 + denovo71436 + denovo71438 + denovo71440 +
						denovo71441 + denovo71446 + denovo71449 + denovo71482 + denovo71516 + denovo71525 + denovo71547 + denovo71576 +
						denovo71594 + denovo71598 + denovo71604 + denovo71649 + denovo71666 + denovo71672 + denovo71674 + denovo71677 +
						denovo71810 + denovo71840 + denovo71843 + denovo71913 + denovo71953 + denovo71958 + denovo72053 +
						denovo72109 + denovo72153 + denovo72161 + denovo72166 + denovo72174 + denovo72199 + denovo72204 +
						denovo72210 + denovo72233 + denovo72284 + denovo72320 + denovo72332 + denovo72457 + denovo72462 + denovo72466 +
						denovo72474 + denovo72485 + denovo72496 + denovo72517 + denovo72554 + denovo72619 + denovo72634 + denovo72695 +
						denovo72739 + denovo72752 + denovo72767 + denovo72820 + denovo72826 + denovo72870 + denovo72892 + denovo72940 +
						denovo72952 + denovo72971 + denovo73038 + denovo73045 + denovo73105 + denovo73190 + denovo73206 + denovo73214 +
						denovo73217 + denovo73347 + denovo73349 + denovo73400 + denovo73404 + denovo73412 + denovo73422 + denovo73448 +
						denovo73451 + denovo73457 + denovo73461 + denovo73581 + denovo73630 + denovo73674 + denovo73710 + denovo73713 +
						denovo73743 + denovo73798 + denovo73816 + denovo73853 + denovo73890 + denovo73902 + denovo73915 + denovo73924 +
						denovo73927 + denovo73955 + denovo73973 + denovo74010 + denovo74063 + denovo74064 + denovo74075 +
						denovo74082 + denovo74097 + denovo74129 + denovo74165 + denovo74194 + denovo74276 + denovo74285 + denovo74297 +
						denovo74301 + denovo74336 + denovo74435 + denovo74440 + denovo74444 + denovo74478 + denovo74493 + denovo74494 +
						denovo74516 + denovo74523 + denovo74525 + denovo74564 + denovo74585 + denovo74676 + denovo74692 + denovo74702 +
						denovo74703 + denovo74781 + denovo74829 + denovo74852 + denovo74929 + denovo74938 + denovo74948 + denovo74955 +
						denovo74957 + denovo74963 + denovo74966 + denovo75041 + denovo75077 + denovo75116 + denovo75125 + denovo75134 +
						denovo75141 + denovo75196 + denovo75213 + denovo75224 + denovo75239 + denovo75252 + denovo75273 + denovo75353 +
						denovo75373 + denovo75444 + denovo75545 + denovo75585 + denovo75648 + denovo75653 + denovo75659 +
						denovo75670 + denovo75694 + denovo75712 + denovo75717 + denovo75718 + denovo75779 + denovo75796 + denovo75851 +
						denovo75859 + denovo75891 + denovo75906 + denovo75910 + denovo75947 + denovo76024 + denovo76040 + denovo76046 +
						denovo76064 + denovo76088 + denovo76094 + denovo76127 + denovo76133 + denovo76135 + denovo76139 + denovo76188 +
						denovo76216 + denovo76330 + denovo76333 + denovo76340 + denovo76345 + denovo76399 + denovo76428 +
						denovo76437 + denovo76440 + denovo76441 + denovo76469 + denovo76479 + denovo76511 + denovo76580 + denovo76594 +
						denovo76603 + denovo76605 + denovo76735 + denovo76762 + denovo76812 + denovo76824 + denovo76837 + denovo76861 +
						denovo76924 + denovo76940 + denovo76966 + denovo76989 + denovo77043 + denovo77089 + denovo77115 + denovo77134 +
						denovo77166 + denovo77215 + denovo77353 + denovo77362 + denovo77377 + denovo77381 + denovo77392 + denovo77399 +
						denovo77415 + denovo77462 + denovo77504 + denovo77524 + denovo77540 + denovo77590 + denovo77626 +
						denovo77657 + denovo77668 + denovo77675 + denovo77728 + denovo77743 + denovo77786 + denovo77811 + denovo77834 +
						denovo77836 + denovo77858 + denovo77859 + denovo77873 + denovo77882 + denovo77905 + denovo77917 + denovo77938 +
						denovo77945 + denovo78076 + denovo78079 + denovo78113 + denovo78120 + denovo78132 + denovo78134 +
						denovo78194 + denovo78285 + denovo78304 + denovo78331 + denovo78347 + denovo78352 + denovo78354 + denovo78360 +
						denovo78372 + denovo78388 + denovo78391 + denovo78415 + denovo78446 + denovo78480 + denovo78488 + denovo78504 +
						denovo78536 + denovo78542 + denovo78557 + denovo78559 + denovo78572 + denovo78599 + denovo78671 + denovo78707 +
						denovo78708 + denovo78765 + denovo78794 + denovo78818 + denovo78873 + denovo78893 + denovo78895 + denovo78897 +
						denovo78963 + denovo78976 + denovo78985 + denovo79050 + denovo79057 + denovo79092 + denovo79121 + denovo79153 +
						denovo79198 + denovo79238 + denovo79254 + denovo79312 + denovo79320 + denovo79348 + denovo79390 + denovo79459 +
						denovo79494 + denovo79502 + denovo79525 + denovo79629 + denovo79677 + denovo79683 + denovo79693 + denovo79731 +
						denovo79749 + denovo79757 + denovo79766 + denovo79797 + denovo79836 + denovo79865 + denovo79899 +
						denovo79905 + denovo79969 + denovo80005 + denovo80009 + denovo80018 + denovo80024 + denovo80046 + denovo80055 +
						denovo80059 + denovo80061 + denovo80143 + denovo80178 + denovo80188 + denovo80197 + denovo80238 + denovo80257 +
						denovo80271 + denovo80301 + denovo80335 + denovo80358 + denovo80386 + denovo80398 + denovo80411 + denovo80430 +
						denovo80439 + denovo80492 + denovo80501 + denovo80535 + denovo80539 + denovo80582 + denovo80588 + denovo80592 +
						denovo80658 + denovo80660 + denovo80661 + denovo80729 + denovo80740 + denovo80820 + denovo80836 + denovo80839 +
						denovo80877 + denovo80889 + denovo80893 + denovo80967 + denovo80987 + denovo81100 + denovo81109 + denovo81123 +
						denovo81188 + denovo81191 + denovo81194 + denovo81213 + denovo81232 + denovo81312 + denovo81313 + denovo81316 +
						denovo81356 + denovo81406 + denovo81439 + denovo81567 + denovo81571 + denovo81603 + denovo81609 + denovo81613 +
						denovo81629 + denovo81666 + denovo81723 + denovo81754 + denovo81769 + denovo81843 + denovo81850 + denovo81870 +
						denovo81879 + denovo81920 + denovo81938 + denovo81974 + denovo81977 + denovo81982 + denovo82076 + denovo82092 +
						denovo82128 + denovo82139 + denovo82141 + denovo82154 + denovo82155 + denovo82156 + denovo82186 + denovo82209 +
						denovo82234 + denovo82317 + denovo82318 + denovo82322 + denovo82434 + denovo82462 + denovo82545 +
						denovo82600 + denovo82604 + denovo82633 + denovo82639 + denovo82649 + denovo82691 + denovo82694 + denovo82700 +
						denovo82747 + denovo82754 + denovo82775 + denovo82789 + denovo82791 + denovo82864 + denovo82875 + denovo82880 +
						denovo82884 + denovo82930 + denovo82993 + denovo82995 + denovo82997 + denovo82999 + denovo83002 + denovo83005 +
						denovo83021 + denovo83062 + denovo83119 + denovo83120 + denovo83131 + denovo83220 + denovo83240 + denovo83256 +
						denovo83284 + denovo83290 + denovo83314 + denovo83458 + denovo83479 + denovo83506 +
						denovo83513 + denovo83531 + denovo83545 + denovo83559 + denovo83574 + denovo83575 + denovo83588 + denovo83594 +
						denovo83605 + denovo83622 + denovo83625 + denovo83669 + denovo83694 + denovo83720 + denovo83725 + denovo83771 +
						denovo83799 + denovo83822 + denovo83855 + denovo83886 + denovo83913 + denovo83945 + denovo83991 + denovo84002 +
						denovo84020 + denovo84076 + denovo84125 + denovo84184 + denovo84224 + denovo84230 + denovo84304 + denovo84346 +
						denovo84364 + denovo84447 + denovo84463 + denovo84465 + denovo84489 + denovo84522 + denovo84545 + denovo84557 +
						denovo84620 + denovo84643 + denovo84660 + denovo84664 + denovo84665 + denovo84750 + denovo84768 + denovo84854 +
						denovo84860 + denovo84928 + denovo85016 + denovo85017 + denovo85043 + denovo85049 + denovo85080 + denovo85124 +
						denovo85145 + denovo85188 + denovo85191 + denovo85216 + denovo85235 + denovo85236 + denovo85241 + denovo85250 +
						denovo85253 + denovo85293 + denovo85297 + denovo85313 + denovo85326 + denovo85331 + denovo85404 + denovo85457 +
						denovo85481 + denovo85510 + denovo85514 + denovo85534 + denovo85537 + denovo85545 + denovo85567 + denovo85574 +
						denovo85616 + denovo85622 + denovo85650 + denovo85670 + denovo85722 + denovo85798 + denovo85851 + denovo85853 +
						denovo85883 + denovo85933 + denovo85937 + denovo85965 + denovo85966 + denovo85979 + denovo85993 + denovo85994 +
						denovo86079 + denovo86086 + denovo86135 + denovo86151 + denovo86182 + denovo86205 + denovo86212 + denovo86221 +
						denovo86254 + denovo86292 + denovo86319 + denovo86330 + denovo86337 + denovo86340 + denovo86350 + denovo86362 +
						denovo86419 + denovo86426 + denovo86531 + denovo86614 + denovo86616 + denovo86637 + denovo86655 + denovo86659 +
						denovo86680 + denovo86698 + denovo86712 + denovo86717 + denovo86733 + denovo86758 + denovo86767 + denovo86784 +
						denovo86797 + denovo86798 + denovo86853 + denovo86867 + denovo86893 + denovo86896 + denovo86911 + denovo86916 +
						denovo86964 + denovo86967 + denovo86991 + denovo87003 + denovo87007 + denovo87036 + denovo87064 + denovo87071 +
						denovo87190 + denovo87226 + denovo87257 + denovo87320 + denovo87382 + denovo87400 + denovo87447 + denovo87452 +
						denovo87483 + denovo87503 + denovo87516 + denovo87524 + denovo87533 + denovo87539 + denovo87551 + denovo87593 +
						denovo87671 + denovo87673 + denovo87675 + denovo87713 + denovo87727 + denovo87760 + denovo87761 + denovo87812 +
						denovo87895 + denovo87962 + denovo87981 + denovo88069 + denovo88077 + denovo88083 + denovo88101 + denovo88115 +
						denovo88163 + denovo88212 + denovo88216 + denovo88246 + denovo88285 + denovo88334 + denovo88346 +
						denovo88366 + denovo88404 + denovo88419 + denovo88449 + denovo88466 + denovo88491 + denovo88500 + denovo88507 +
						denovo88612 + denovo88624 + denovo88679 + denovo88688 + denovo88710 + denovo88730 + denovo88771 + denovo88799 +
						denovo88817 + denovo88837 + denovo88911 + denovo88934 + denovo89012 + denovo89048 + denovo89056 + denovo89069 +
						denovo89074 + denovo89077 + denovo89096 + denovo89117 + denovo89163 + denovo89168 + denovo89172 +
						denovo89187 + denovo89192 + denovo89198 + denovo89227 + denovo89231 + denovo89263 + denovo89272 + denovo89297 +
						denovo89309 + denovo89343 + denovo89370 + denovo89386 + denovo89417 + denovo89431 + denovo89444 + denovo89455 +
						denovo89457 + denovo89459 + denovo89513 + denovo89561 + denovo89574 + denovo89577 + denovo89598 + denovo89618 +
						denovo89648 + denovo89662 + denovo89671 + denovo89798 + denovo89803 + denovo89808 + denovo89828 + denovo89846 +
						denovo89870 + denovo89901 + denovo89905 + denovo89939 + denovo89975 + denovo90065 + denovo90084 + denovo90115 +
						denovo90159 + denovo90175 + denovo90194 + denovo90211 + denovo90264 + denovo90277 + denovo90319 + denovo90344 +
						denovo90353 + denovo90383 + denovo90404 + denovo90414 + denovo90475 + denovo90509 + denovo90530 + denovo90551 +
						denovo90553 + denovo90563 + denovo90564 + denovo90689 + denovo90693 + denovo90696 + denovo90707 + denovo90710 +
						denovo90713 + denovo90721 + denovo90734 + denovo90796 + denovo90806 + denovo90823 + denovo90853 + denovo90890 +
						denovo90930 + denovo91013 + denovo91033 + denovo91090 + denovo91095 + denovo91147 + denovo91161 + denovo91186 +
						denovo91205 + denovo91206 + denovo91214 + denovo91227 + denovo91228 + denovo91302 + denovo91346 +
						denovo91353 + denovo91357 + denovo91389 + denovo91393 + denovo91435 + denovo91466 + denovo91482 + denovo91498 +
						denovo91531 + denovo91568 + denovo91582 + denovo91608 + denovo91623 + denovo91671 + denovo91676 + denovo91677 +
						denovo91707 + denovo91729 + denovo91743 + denovo91747 + denovo91773 + denovo91779 + denovo91799 + denovo91808 +
						denovo91813 + denovo91826 + denovo91840 + denovo91945 + denovo91951 + denovo91956 + denovo91980 + denovo92003 +
						denovo92017 + denovo92024 + denovo92025 + denovo92044 + denovo92046 + denovo92077 + denovo92124 + denovo92132 +
						denovo92169 + denovo92172 + denovo92197 + denovo92214 + denovo92235 + denovo92237 + denovo92277 + denovo92284 +
						denovo92302 + denovo92327 + denovo92387 + denovo92401 + denovo92407 + denovo92478 + denovo92490 +
						denovo92507 + denovo92512 + denovo92588 + denovo92589 + denovo92607 + denovo92616 + denovo92622 + denovo92625 +
						denovo92634 + denovo92657 + denovo92677 + denovo92705 + denovo92731 + denovo92773 +
						denovo92809 + denovo92875 + denovo92885 + denovo92918 + denovo92978 + denovo92986 + denovo93012 + denovo93018 +
						denovo93057 + denovo93072 + denovo93106 + denovo93146 + denovo93156 + denovo93165 + denovo93183 + denovo93192 +
						denovo93221 + denovo93256 + denovo93292 + denovo93321 + denovo93346 + denovo93350 + denovo93355 + denovo93363 +
						denovo93376 + denovo93391 + denovo93405 + denovo93435 + denovo93446 + denovo93473 + denovo93480 +
						denovo93488 + denovo93730 + denovo93778 + denovo93786 + denovo93824 + denovo93862 + denovo93915 +
						denovo93918 + denovo93932 + denovo93935 + denovo93936 + denovo93957 + denovo93985 + denovo94010 + denovo94120 +
						denovo94141 + denovo94146 + denovo94151 + denovo94152 + denovo94170 + denovo94185 + denovo94209 + denovo94213 +
						denovo94217 + denovo94230 + denovo94306 + denovo94324 + denovo94332 + denovo94335 + denovo94336 + denovo94402 +
						denovo94530 + denovo94630 + denovo94632 + denovo94644 + denovo94715 + denovo94741 + denovo94761 + denovo94762 +
						denovo94784 + denovo94804 + denovo94840 + denovo94848 + denovo94862 + denovo94880 + denovo94899 + denovo94928 +
						denovo94978 + denovo94988 + denovo95005 + denovo95033 + denovo95057 + denovo95084 + denovo95087 + denovo95115 +
						denovo95131 + denovo95148 + denovo95157 + denovo95167 + denovo95172 + denovo95181 + denovo95233 + denovo95292 +
						denovo95302 + denovo95315 + denovo95370 + denovo95387 + denovo95404 + denovo95413 + denovo95438 +
						denovo95440 + denovo95463 + denovo95496 + denovo95521 + denovo95527 + denovo95588 + denovo95595 + denovo95596 +
						denovo95604 + denovo95614 + denovo95627 + denovo95648 + denovo95676 + denovo95763 + denovo95771 + denovo95820 +
						denovo95856 + denovo95900 + denovo95910 + denovo96003 + denovo96018 + denovo96022 + denovo96128 + denovo96134 +
						denovo96174 + denovo96188 + denovo96263 + denovo96302 + denovo96320 + denovo96338 + denovo96406 + denovo96433 +
						denovo96470 + denovo96497 + denovo96506 + denovo96532 + denovo96534 + denovo96654 + denovo96659 + denovo96677 +
						denovo96688 + denovo96702 + denovo96712 + denovo96719 + denovo96769 + denovo96816 + denovo96855 + denovo96940 +
						denovo96943 + denovo96948 + denovo96958 + denovo96966 + denovo96971 + denovo97027 + denovo97057 + denovo97076 +
						denovo97093 + denovo97106 + denovo97151 + denovo97154 + denovo97185 + denovo97204 + denovo97267 + denovo97334 +
						denovo97372 + denovo97396 + denovo97406 + denovo97410 + denovo97414 + denovo97432 + denovo97452 + denovo97510 +
						denovo97542 + denovo97577 + denovo97611 + denovo97618 + denovo97657 + denovo97686 + denovo97700 + denovo97710 +
						denovo97713 + denovo97721 + denovo97749 + denovo97803 + denovo97807 + denovo97844 + denovo97849 + denovo97874 +
						denovo97886 + denovo97913 + denovo97931 + denovo97954 + denovo98038 + denovo98083 + denovo98129 + denovo98155 +
						denovo98175 + denovo98183 + denovo98189 + denovo98235 + denovo98249 + denovo98257 + denovo98285 + denovo98308 +
						denovo98352 + denovo98363 + denovo98382 + denovo98428 + denovo98433 + denovo98554 + denovo98636 +
						denovo98661 + denovo98663 + denovo98698 + denovo98739 + denovo98749 + denovo98760 + denovo98762 +
						denovo98779 + denovo98809 + denovo98841 + denovo98886 + denovo98892 + denovo98898 + denovo98900 + denovo98907 +
						denovo98935 + denovo98942 + denovo98967 + denovo98992 + denovo99023 + denovo99074 + denovo99109 + denovo99118 +
						denovo99127 + denovo99133 + denovo99144 + denovo99173 + denovo99257 + denovo99281 + denovo99399 + denovo99450 +
						denovo99451 + denovo99452 + denovo99455 + denovo99559 + denovo99561 + denovo99563 + denovo99589 + denovo99606 +
						denovo99690 + denovo99736 + denovo99758 + denovo99780 + denovo99785 + denovo99854 + denovo99876 + denovo99896 +
						denovo99927 + denovo99929 + denovo99965 + denovo99978 + denovo100004 + denovo100052 + denovo100080 + denovo100132 +
						denovo100139 + denovo100160 + denovo100176 + denovo100188 + denovo100208 + denovo100289 + denovo100310 +
						denovo100334 + denovo100360 + denovo100408 + denovo100440 + denovo100446 + denovo100473 + denovo100534 + denovo100538 +
						denovo100550 + denovo100563 + denovo100576 + denovo100613 + denovo100617 + denovo100662 + denovo100783 +
						denovo100889 + denovo100892 + denovo100894 + denovo100906 + denovo100923 + denovo100953 + denovo100992 + denovo101027 +
						denovo101034 + denovo101049 + denovo101069 + denovo101094 + denovo101105 + denovo101112 + denovo101119 + denovo101164 +
						denovo101173 + denovo101225 + denovo101239 + denovo101252 + denovo101288 + denovo101296 + denovo101357 +
						denovo101382 + denovo101405 + denovo101423 + denovo101444 + denovo101458 + denovo101472 + denovo101503 + denovo101505 +
						denovo101517 + denovo101525 + denovo101598 + denovo101640 + denovo101644 + denovo101669 + denovo101744 + denovo101750 +
						denovo101784 + denovo101803 + denovo101815 + denovo101821 + denovo101830 + denovo101876 + denovo101896 +
						denovo101941 + denovo101953 + denovo101964 + denovo102006 + denovo102014 + denovo102021 + denovo102030 + denovo102108 +
						denovo102111 + denovo102168 + denovo102228 + denovo102239 + denovo102244 + denovo102260 + denovo102262 + denovo102288 +
						denovo102303 + denovo102329 + denovo102335 + denovo102366 + denovo102427 + denovo102433 + denovo102468 + denovo102480 +
						denovo102494 + denovo102600 + denovo102604 + denovo102609 + denovo102619 + denovo102620 +
						denovo102633 + denovo102681 + denovo102682 + denovo102699 + denovo102728 + denovo102737 + denovo102789 +
						denovo102865 + denovo102866 + denovo102887 + denovo102905 + denovo102982 + denovo103023 + denovo103036 + denovo103060 +
						denovo103083 + denovo103108 + denovo103133 + denovo103135 + denovo103138 + denovo103153 + denovo103165 + denovo103168 +
						denovo103172 + denovo103181 + denovo103223 + denovo103241 + denovo103277 + denovo103357 + denovo103358 + denovo103425 +
						denovo103444 + denovo103454 + denovo103470 + denovo103507 + denovo103509 + denovo103516 + denovo103544 + denovo103547 +
						denovo103555 + denovo103557 + denovo103569 + denovo103586 + denovo103611 + denovo103617 + denovo103638 + denovo103667 +
						denovo103720 + denovo103745 + denovo103781 + denovo103838 + denovo103864 + denovo103914 + denovo104042 + denovo104044 +
						denovo104051 + denovo104066 + denovo104079 + denovo104093 + denovo104100 + denovo104104 + denovo104107 + denovo104115 +
						denovo104127 + denovo104143 + denovo104177 + denovo104194 + denovo104200 + denovo104210 + denovo104214 + denovo104238 +
						denovo104252 + denovo104283 + denovo104328 + denovo104329 + denovo104364 + denovo104397 + denovo104419 + denovo104421 +
						denovo104422 + denovo104437 + denovo104494 + denovo104562 + denovo104591 + denovo104618 + denovo104633 + denovo104661 +
						denovo104667 + denovo104699 + denovo104759 + denovo104792 + denovo104801 + denovo104836 + denovo104873 + denovo104881 +
						denovo104885 + denovo104930 + denovo104936 + denovo104938 + denovo104986 + denovo105076 + denovo105086 + denovo105088 +
						denovo105146 + denovo105158 + denovo105231 + denovo105241 + denovo105341 + denovo105345 + denovo105354 + denovo105361 +
						denovo105367 + denovo105396 + denovo105417 + denovo105445 + denovo105492 + denovo105509 + denovo105545 + denovo105588 +
						denovo105623 + denovo105658 + denovo105661 + denovo105731 + denovo105732 + denovo105737 + denovo105755 + denovo105783 +
						denovo105799 + denovo105801 + denovo105804 + denovo105823 + denovo105948 + denovo105954 + denovo105956 + denovo105962 +
						denovo105965 + denovo105982 + denovo105985 + denovo105990 + denovo106067 + denovo106077 + denovo106093 + denovo106145 +
						denovo106160 + denovo106224 + denovo106281 + denovo106288 + denovo106289 + denovo106299 + denovo106381 + denovo106384 +
						denovo106392 + denovo106393 + denovo106431 + denovo106492 + denovo106496 + denovo106521 + denovo106532 + denovo106566 +
						denovo106578 + denovo106585 + denovo106614 + denovo106620 + denovo106627 + denovo106663 + denovo106671 + denovo106770 +
						denovo106801 + denovo106802 + denovo106804 + denovo106840 + denovo106857 + denovo106868 + denovo106875 + denovo106895 +
						denovo106904 + denovo106950 + denovo106999 + denovo107014 + denovo107017 + denovo107072 + denovo107105 + denovo107106 +
						denovo107202 + denovo107220 + denovo107255 + denovo107294 + denovo107340 + denovo107374 + denovo107388 + denovo107397 +
						denovo107409 + denovo107431 + denovo107464 + denovo107470 + denovo107499 + denovo107500 + denovo107508 +
						denovo107548 + denovo107571 + denovo107586 + denovo107612 + denovo107629 + denovo107631 + denovo107651 + denovo107669 +
						denovo107686 + denovo107732 + denovo107752 + denovo107758 + denovo107873 + denovo107882 + denovo107921 + denovo107930 +
						denovo107970 + denovo108021 + denovo108047 + denovo108074 + denovo108096 + denovo108107 + denovo108223 + denovo108249 +
						denovo108251 + denovo108291 + denovo108316 + denovo108355 + denovo108432 + denovo108433 + denovo108457 + denovo108483 +
						denovo108502 + denovo108527 + denovo108534 + denovo108552 + denovo108572 + denovo108588 + denovo108596 + denovo108601 +
						denovo108623 + denovo108632 + denovo108646 + denovo108671 + denovo108704 + denovo108723 + denovo108743 + denovo108832 +
						denovo108846 + denovo108912 + denovo108921 + denovo108932 + denovo108951 + denovo109023 + denovo109025 + denovo109034 +
						denovo109037 + denovo109061 + denovo109071 + denovo109094 + denovo109104 + denovo109109 + denovo109203 + denovo109266 +
						denovo109301 + denovo109324 + denovo109347 + denovo109398 + denovo109524 + denovo109556 + denovo109659 + denovo109682 +
						denovo109722 + denovo109734 + denovo109751 + denovo109754 + denovo109758 + denovo109851 + denovo109853 + denovo109863 +
						denovo109877 + denovo109880 + denovo109899 + denovo109921 + denovo109944 + denovo109977 + denovo109982 + denovo109992 +
						denovo110026 + denovo110034 + denovo110036 + denovo110079 + denovo110093 + denovo110121 + denovo110141 + denovo110157 +
						denovo110173 + denovo110180 + denovo110228 + denovo110234 + denovo110235 + denovo110248 + denovo110300 + denovo110312 +
						denovo110357 + denovo110371 + denovo110418 + denovo110483 + denovo110521 + denovo110533 + denovo110550 + denovo110586 +
						denovo110608 + denovo110610 + denovo110651 + denovo110699 + denovo110705 + denovo110713 + denovo110716 +
						denovo110721 + denovo110764 + denovo110798 + denovo110817 + denovo110826 + denovo110832 + denovo110881 + denovo110902 +
						denovo110932 + denovo110945 + denovo110950 + denovo110966 + denovo110969 + denovo110988 + denovo111002 + denovo111028 +
						denovo111048 + denovo111074 + denovo111078 + denovo111110 + denovo111126 + denovo111145 + denovo111146 +
						denovo111159 + denovo111166 + denovo111173 + denovo111197 + denovo111262 + denovo111274 + denovo111294 + denovo111300 +
						denovo111343 + denovo111347 + denovo111370 + denovo111376 + denovo111383 + denovo111411 + denovo111415 + denovo111497 +
						denovo111506 + denovo111508 + denovo111516 + denovo111528 + denovo111565 + denovo111596 + denovo111648 + denovo111664 +
						denovo111666 + denovo111669 + denovo111674 + denovo111728 + denovo111730 + denovo111764 + denovo111786 + denovo111817 +
						denovo111841 + denovo111851 + denovo111906 + denovo111939 + denovo111950 + denovo111960 + denovo111967 + denovo112075 +
						denovo112108 + denovo112121 + denovo112146 + denovo112211 + denovo112235 + denovo112262 + denovo112278 + denovo112317 +
						denovo112321 + denovo112322 + denovo112331 + denovo112345 + denovo112346 + denovo112431 + denovo112451 +
						denovo112490 + denovo112507 + denovo112552 + denovo112560 + denovo112590 + denovo112593 + denovo112615 +
						denovo112655 + denovo112706 + denovo112715 + denovo112728 + denovo112813 + denovo112818 +
						denovo112829 + denovo112837 + denovo112868 + denovo112884 + denovo112891 + denovo112893 + denovo112896 + denovo112903 +
						denovo112918 + denovo112970 + denovo112982 + denovo112987 + denovo113033 + denovo113077 + denovo113103 + denovo113109 +
						denovo113144 + denovo113176 + denovo113237 + denovo113257 + denovo113275 + denovo113279 + denovo113284 + denovo113304 +
						denovo113316 + denovo113345 + denovo113404 + denovo113413 + denovo113416 + denovo113425 + denovo113428 +
						denovo113452 + denovo113471 + denovo113502 + denovo113545 + denovo113567 + denovo113572 + denovo113592 +
						denovo113598 + denovo113616 + denovo113634 + denovo113638 + denovo113680 + denovo113685 + denovo113689 + denovo113698 +
						denovo113757 + denovo113760 + denovo113769 + denovo113770 + denovo113780 + denovo113786 + denovo113790 + denovo113791 +
						denovo113808 + denovo113815 + denovo113831 + denovo113859 + denovo113877 + denovo113881 + denovo113885 + denovo113925 +
						denovo113930 + denovo113961 + denovo113996 + denovo114061 + denovo114081 + denovo114191 + denovo114192 + denovo114215 +
						denovo114260 + denovo114268 + denovo114271 + denovo114290 + denovo114302 + denovo114316 + denovo114346 + denovo114370 +
						denovo114378 + denovo114380 + denovo114397 + denovo114413 + denovo114423 + denovo114461 + denovo114488 +
						denovo114508 + denovo114520 + denovo114527 + denovo114579 + denovo114593 + denovo114604 + denovo114678 +
						denovo114709 + denovo114712 + denovo114762 + denovo114780 + denovo114786 + denovo114821 + denovo114868 + denovo114869 +
						denovo114898 + denovo114908 + denovo114920 + denovo114926 + denovo114977 + denovo115038 + denovo115043 + denovo115090 +
						denovo115118 + denovo115175 + denovo115212 + denovo115219 + denovo115251 + denovo115298 + denovo115300 + denovo115302 +
						denovo115358 + denovo115362 + denovo115374 + denovo115377 + denovo115394 + denovo115399 + denovo115415 + denovo115421 +
						denovo115448 + denovo115451 + denovo115472 + denovo115474 + denovo115511 + denovo115518 + denovo115532 + denovo115536 +
						denovo115551 + denovo115562 + denovo115588 + denovo115624 + denovo115629 + denovo115665 + denovo115669 + denovo115695 +
						denovo115708 + denovo115744 + denovo115792 + denovo115803 + denovo115864 + denovo115973 + denovo115983 +
						denovo116020 + denovo116045 + denovo116056 + denovo116061 + denovo116071 + denovo116092 + denovo116096 +
						denovo116145 + denovo116156 + denovo116191 + denovo116193 + denovo116197 + denovo116223 + denovo116253 + denovo116277 +
						denovo116279 + denovo116317 + denovo116330 + denovo116341 + denovo116357 + denovo116374 + denovo116378 + denovo116454 +
						denovo116482 + denovo116517 + denovo116518 + denovo116544 + denovo116549 + denovo116578 + denovo116583 + denovo116619 +
						denovo116622 + denovo116642 + denovo116662 + denovo116670 + denovo116689 + denovo116718 + denovo116724 + denovo116728 +
						denovo116777 + denovo116784 + denovo116794 + denovo116796 + denovo116841 + denovo116843 + denovo116850 + denovo116866 +
						denovo116882 + denovo116891 + denovo116892 + denovo116938 + denovo116961 + denovo116966 + denovo116971 + denovo116989 +
						denovo116990 + denovo117001 + denovo117061 + denovo117112 + denovo117148 + denovo117157 + denovo117165 + denovo117201 +
						denovo117252 + denovo117256 + denovo117288 + denovo117320 + denovo117391 + denovo117394 + denovo117422 + denovo117462 +
						denovo117493 + denovo117517 + denovo117537 + denovo117558 + denovo117581 + denovo117599 + denovo117604 + denovo117633 +
						denovo117651 + denovo117652 + denovo117657 + denovo117660 + denovo117725 + denovo117742 + denovo117765 + denovo117777 +
						denovo117804 + denovo117812 + denovo117821 + denovo117831 + denovo117896 + denovo117904 + denovo117919 + denovo117923 +
						denovo117941 + denovo117962 + denovo117988 + denovo118061 + denovo118084 + denovo118146 + denovo118158 + denovo118174 +
						denovo118189 + denovo118205 + denovo118252 + denovo118271 + denovo118277 + denovo118303 + denovo118317 +
						denovo118321 + denovo118325 + denovo118355 + denovo118364 + denovo118412 + denovo118432 + denovo118460 + denovo118477 +
						denovo118486 + denovo118506 + denovo118525 + denovo118604 + denovo118636 + denovo118647 + denovo118670 + denovo118714 +
						denovo118741 + denovo118771 + denovo118799 + denovo118807 + denovo118820 + denovo118851 + denovo118877 + denovo118880 +
						denovo118886 + denovo118916 + denovo118923 + denovo118954 + denovo118957 + denovo118967 + denovo118968 + denovo118988 +
						denovo118989 + denovo118994 + denovo119047 + denovo119059 + denovo119068 + denovo119072 + denovo119086 + denovo119104 +
						denovo119115 + denovo119241 + denovo119250 + denovo119278 + denovo119320 + denovo119333 + denovo119341 + denovo119342 +
						denovo119366 + denovo119394 + denovo119466 + denovo119474 + denovo119546 + denovo119567 + denovo119590 + denovo119613 +
						denovo119640 + denovo119671 + denovo119678 + denovo119680 + denovo119700 + denovo119730 + denovo119742 + denovo119758 +
						denovo119769 + denovo119779 + denovo119793 + denovo119826 + denovo119869 + denovo119883 +
						denovo119958 + denovo119980 + denovo120005 + denovo120009 + denovo120013 + denovo120057 + denovo120068 + denovo120072 +
						denovo120076 + denovo120082 + denovo120109 + denovo120140 + denovo120152 + denovo120159 + denovo120210 + denovo120266 +
						denovo120269 + denovo120294 + denovo120309 + denovo120316 + denovo120358 + denovo120361 + denovo120363 + denovo120379 +
						denovo120401 + denovo120469 + denovo120495 + denovo120499 + denovo120503 + denovo120511 + denovo120531 + denovo120544 +
						denovo120590 + denovo120624 + denovo120636 + denovo120642 + denovo120653 + denovo120672 + denovo120685 + denovo120695 +
						denovo120706 + denovo120710 + denovo120779 + denovo120786 + denovo120790 + denovo120799 + denovo120828 + denovo120840 +
						denovo120845 + denovo120848 + denovo120858 + denovo120901 + denovo121027 + denovo121064 + denovo121110 + denovo121116 +
						denovo121120 + denovo121137 + denovo121219 + denovo121223 + denovo121269 + denovo121275 + denovo121294 + denovo121308 +
						denovo121312 + denovo121381 + denovo121384 + denovo121386 + denovo121401 + denovo121420 + denovo121433 + denovo121436 +
						denovo121441 + denovo121478 + denovo121489 + denovo121516 + denovo121519 + denovo121523 + denovo121551 + denovo121560 +
						denovo121570 + denovo121640 + denovo121675 + denovo121684 + denovo121721 + denovo121782 + denovo121785 + denovo121800 +
						denovo121801 + denovo121827 + denovo121860 + denovo121925 + denovo121970 + denovo121990 +
						denovo121999 + denovo122006 + denovo122072 + denovo122092 + denovo122096 + denovo122113 + denovo122132 + denovo122165 +
						denovo122168 + denovo122178 + denovo122213 + denovo122224 + denovo122259 + denovo122298 + denovo122316 + denovo122344 +
						denovo122378 + denovo122380 + denovo122414 + denovo122480 + denovo122490 + denovo122491 + denovo122513 + denovo122518 +
						denovo122520 + denovo122548 + denovo122553 + denovo122556 + denovo122565 + denovo122568 + denovo122573 + denovo122583 +
						denovo122592 + denovo122615 + denovo122639 + denovo122645 + denovo122666 + denovo122667 + denovo122746 + denovo122780 +
						denovo122797 + denovo122798 + denovo122845 + denovo122906 + denovo122921 + denovo122925 + denovo122930 + denovo122954 +
						denovo122960 + denovo122977 + denovo123001 + denovo123074 + denovo123077 + denovo123115 + denovo123118 + denovo123119 +
						denovo123145 + denovo123157 + denovo123162 + denovo123170 + denovo123214 + denovo123215 + denovo123261 + denovo123303 +
						denovo123321 + denovo123360 + denovo123447 + denovo123485 + denovo123509 + denovo123519 + denovo123533 +
						denovo123597 + denovo123620 + denovo123646 + denovo123648 + denovo123671 + denovo123675 + denovo123776 +
						denovo123806 + denovo123818 + denovo123827 + denovo123853 + denovo123874 + denovo123882 + denovo123884 + denovo123895 +
						denovo123914 + denovo123917 + denovo123967 + denovo123979 + denovo124066 + denovo124070 + denovo124078 + denovo124083 +
						denovo124096 + denovo124105 + denovo124106 + denovo124111 + denovo124135 + denovo124158 + denovo124251 + denovo124272 +
						denovo124332 + denovo124333 + denovo124360 + denovo124390 + denovo124413 + denovo124434 + denovo124468 + denovo124493 +
						denovo124512 + denovo124551 + denovo124552 + denovo124563 + denovo124573 + denovo124591 + denovo124611 + denovo124623 +
						denovo124626 + denovo124632 + denovo124639 + denovo124646 + denovo124648 + denovo124676 + denovo124682 +
						denovo124700 + denovo124701 + denovo124703 + denovo124711 + denovo124741 + denovo124750 + denovo124752 + denovo124772 +
						denovo124782 + denovo124817 + denovo124822 + denovo124839 + denovo124845 + denovo124858 + denovo124864 + denovo124877 +
						denovo124988 + denovo125002 + denovo125014 + denovo125074 + denovo125080 + denovo125130 + denovo125139 + denovo125140 +
						denovo125155 + denovo125170 + denovo125179 + denovo125207 + denovo125221 + denovo125244 + denovo125249 + denovo125251 +
						denovo125290 + denovo125335 + denovo125365 + denovo125458 + denovo125499 + denovo125550 + denovo125592 + denovo125594 +
						denovo125635 + denovo125689 + denovo125706 + denovo125732 + denovo125737 + denovo125744 + denovo125823 + denovo125854 +
						denovo125870 + denovo125883 + denovo125923 + denovo125977 + denovo125998 + denovo126003 + denovo126029 + denovo126072 +
						denovo126076 + denovo126091 + denovo126119 + denovo126148 + denovo126152 + denovo126178 + denovo126196 + denovo126200 +
						denovo126207 + denovo126209 + denovo126212 + denovo126226 + denovo126280 + denovo126357 + denovo126363 + denovo126378 +
						denovo126387 + denovo126396 + denovo126408 + denovo126415 + denovo126468 + denovo126472 + denovo126495 + denovo126510 +
						denovo126572 + denovo126600 + denovo126659 + denovo126671 + denovo126677 + denovo126699 + denovo126700 + denovo126726 +
						denovo126729 + denovo126741 + denovo126752 + denovo126781 + denovo126799 + denovo126852 + denovo126881 + denovo126891 +
						denovo126911 + denovo126980 + denovo127002 + denovo127054 + denovo127098 + denovo127118 + denovo127141 + denovo127143 +
						denovo127146 + denovo127164 + denovo127171 + denovo127187 + denovo127189 + denovo127198 + denovo127202 + denovo127219 +
						denovo127233 + denovo127250 + denovo127251 + denovo127357 + denovo127415 + denovo127440 + denovo127471 + denovo127477 +
						denovo127482 + denovo127499 + denovo127521 + denovo127568 + denovo127573 + denovo127575 + denovo127637 + denovo127663 +
						denovo127677 + denovo127692 + denovo127701 + denovo127721 + denovo127733 + denovo127758 + denovo127817 + denovo127851 +
						denovo127873 + denovo127900 + denovo127945 + denovo127955 + denovo127957 + denovo127966 + denovo127971 + denovo127972 +
						denovo127975 + denovo127980 + denovo127992 + denovo127998 + denovo128053 + denovo128071 + denovo128101 + denovo128107 +
						denovo128116 + denovo128125 + denovo128140 + denovo128175 + denovo128178 + denovo128221 + denovo128223 + denovo128250 +
						denovo128263 + denovo128293 + denovo128315 + denovo128325 + denovo128354 + denovo128371 + denovo128391 + denovo128398 +
						denovo128401 + denovo128462 + denovo128463 + denovo128476 + denovo128515 + denovo128557 + denovo128589 +
						denovo128592 + denovo128622 + denovo128623 + denovo128634 + denovo128644 + denovo128672 + denovo128675 + denovo128701 +
						denovo128736 + denovo128747 + denovo128762 + denovo128764 + denovo128783 + denovo128796 + denovo128814 + denovo128816 +
						denovo128829 + denovo128838 + denovo128845 + denovo128863 + denovo128893 + denovo128909 + denovo128917 + denovo128920 +
						denovo128936 + denovo128964 + denovo128968 + denovo128977 + denovo129007 + denovo129013 + denovo129141 + denovo129152 +
						denovo129177 + denovo129196 + denovo129198 + denovo129235 + denovo129238 + denovo129240 + denovo129254 + denovo129258 +
						denovo129277 + denovo129283 + denovo129284 + denovo129301 + denovo129310 + denovo129324 + denovo129333 + denovo129380 +
						denovo129384 + denovo129394 + denovo129399 + denovo129408 + denovo129428 + denovo129440 + denovo129441 + denovo129514 +
						denovo129543 + denovo129626 + denovo129692 + denovo129704 + denovo129782 + denovo129784 + denovo129803 + denovo129844 +
						denovo129848 + denovo129886 + denovo129904 + denovo129947 + denovo129952 + denovo129999 + denovo130002 + denovo130019 +
						denovo130022 + denovo130116 + denovo130120 + denovo130178 + denovo130234 + denovo130253 + denovo130254 +
						denovo130285 + denovo130302 + denovo130313 + denovo130341 + denovo130370 + denovo130377 + denovo130419 + denovo130437 +
						denovo130455 + denovo130482 + denovo130519 + denovo130524 + denovo130551 + denovo130565 + denovo130573 + denovo130614 +
						denovo130629 + denovo130658 + denovo130666 + denovo130722 + denovo130730 + denovo130742 + denovo130852 +
						denovo130858 + denovo130944 + denovo130962 + denovo130997 + denovo131009 + denovo131022 + denovo131053 + denovo131078 +
						denovo131087 + denovo131138 + denovo131144 + denovo131163 + denovo131166 + denovo131187 + denovo131208 + denovo131261 +
						denovo131302 + denovo131303 + denovo131309 + denovo131347 + denovo131348 + denovo131350 + denovo131358 + denovo131359 +
						denovo131404 + denovo131407 + denovo131411 + denovo131414 + denovo131426 + denovo131468 + denovo131531 + denovo131565 +
						denovo131570 + denovo131571 + denovo131574 + denovo131726 + denovo131779 + denovo131782 + denovo131820 + denovo131825 +
						denovo131832 + denovo131839 + denovo131843 + denovo131885 + denovo131932 + denovo131949 + denovo131965 + denovo131971 +
						denovo131973 + denovo132001 + denovo132004 + denovo132022 + denovo132024 + denovo132032 + denovo132095 + denovo132137 +
						denovo132140 + denovo132164 + denovo132196 + denovo132217 + denovo132263 + denovo132306 + denovo132322 + denovo132421 +
						denovo132425 + denovo132468 + denovo132469 + denovo132504 + denovo132511 + denovo132515 + denovo132535 + denovo132597 +
						denovo132714 + denovo132812 + denovo132839 + denovo132854 + denovo132897 + denovo132920 + denovo132929 + denovo132933 +
						denovo132986 + denovo132990 + denovo132993 + denovo133026 + denovo133027 + denovo133085 + denovo133221 + denovo133247 +
						denovo133293 + denovo133295 + denovo133319 + denovo133322 + denovo133331 + denovo133353 + denovo133358 +
						denovo133404 + denovo133446 + denovo133470 + denovo133520 + denovo133532 + denovo133542 + denovo133590 + denovo133599 +
						denovo133605 + denovo133632 + denovo133647 + denovo133665 + denovo133741 + denovo133746 + denovo133781 + denovo133792 +
						denovo133808 + denovo133858 + denovo133859 + denovo133872 + denovo133873 + denovo133891 + denovo133905 + denovo133916 +
						denovo133934 + denovo133939 + denovo133957 + denovo133971 + denovo134018 + denovo134052 + denovo134098 + denovo134109 +
						denovo134145 + denovo134274 + denovo134285 + denovo134292 + denovo134324 + denovo134332 + denovo134336 + denovo134409 +
						denovo134417 + denovo134436 + denovo134441 + denovo134444 + denovo134449 + denovo134469 + denovo134495 + denovo134517 +
						denovo134546 + denovo134579 + denovo134586 + denovo134629 + denovo134668 + denovo134672 + denovo134745 + denovo134782 +
						denovo134794 + denovo134813 + denovo134832 + denovo134893 + denovo134902 + denovo134907 + denovo134926 + denovo134972 +
						denovo135007 + denovo135013 + denovo135018 + denovo135063 + denovo135084 + denovo135105 + denovo135117 + denovo135133 +
						denovo135137 + denovo135158 + denovo135163 + denovo135166 + denovo135171 + denovo135221 + denovo135222 + denovo135226 +
						denovo135255 + denovo135261 + denovo135287 + denovo135298 + denovo135316 + denovo135330 + denovo135334 + denovo135335 +
						denovo135337 + denovo135388 + denovo135499 + denovo135508 + denovo135569 + denovo135601 + denovo135629 + denovo135693 +
						denovo135731 + denovo135790 + denovo135823 + denovo135837 + denovo135853 + denovo135855 + denovo135863 + denovo135884 +
						denovo135897 + denovo135943 + denovo135945 + denovo135968 + denovo135972 + denovo136019 + denovo136026 + denovo136067 +
						denovo136083 + denovo136087 + denovo136094 + denovo136097 + denovo136111 + denovo136164 + denovo136194 + denovo136196 +
						denovo136233 + denovo136274 + denovo136280 + denovo136286 + denovo136350 + denovo136375 +
						denovo136408 + denovo136446 + denovo136456 + denovo136545 + denovo136569 + denovo136588 + denovo136593 + denovo136598 +
						denovo136600 + denovo136623 + denovo136626 + denovo136633 + denovo136649 + denovo136655 + denovo136709 + denovo136826 +
						denovo136859 + denovo136861 + denovo136871 + denovo136920 + denovo136946 + denovo136990 + denovo137066 + denovo137149 +
						denovo137159 + denovo137221 + denovo137230 + denovo137292 + denovo137305 + denovo137318 + denovo137442 + denovo137473 +
						denovo137505 + denovo137514 + denovo137542 + denovo137549 + denovo137552 + denovo137553 + denovo137579 + denovo137590 +
						denovo137597 + denovo137627 + denovo137657 + denovo137664 + denovo137700 + denovo137738 + denovo137757 +
						denovo137793 + denovo137802 + denovo137845 + denovo137849 + denovo137881 + denovo137888 + denovo137912 + denovo137991 +
						denovo137998 + denovo138054 + denovo138099 + denovo138127 + denovo138173 + denovo138176 + denovo138225 + denovo138228 +
						denovo138233 + denovo138248 + denovo138266 + denovo138309 + denovo138313 + denovo138332 + denovo138360 + denovo138374 +
						denovo138375 + denovo138383 + denovo138431 + denovo138477 + denovo138486 + denovo138500 + denovo138501 + denovo138503 +
						denovo138504 + denovo138518 + denovo138524 + denovo138535 + denovo138557 + denovo138621 + denovo138641 + denovo138656 +
						denovo138662 + denovo138689 + denovo138741 + denovo138766 + denovo138774 + denovo138801 + denovo138815 + denovo138822 +
						denovo138850 + denovo138914 + denovo138947 + denovo138954 + denovo138957 + denovo138958 + denovo138967 + denovo138991 +
						denovo139052 + denovo139084 + denovo139136 + denovo139148 + denovo139156 + denovo139161 + denovo139186 + denovo139219 +
						denovo139237 + denovo139269 + denovo139275 + denovo139283 + denovo139344 + denovo139379 + denovo139386 + denovo139437 +
						denovo139438 + denovo139454 + denovo139466 + denovo139471 + denovo139477 + denovo139565 +
						denovo139581 + denovo139612 + denovo139638 + denovo139675 + denovo139683 + denovo139695 + denovo139722 + denovo139724 +
						denovo139740 + denovo139759 + denovo139772 + denovo139813 + denovo139852 + denovo139869 + denovo139880 + denovo139925 +
						denovo139994 + denovo140005 + denovo140026 + denovo140096 + denovo140101 + denovo140103 + denovo140116 +
						denovo140126 + denovo140169 + denovo140193 + denovo140220 + denovo140307 + denovo140332 + denovo140360 +
						denovo140363 + denovo140372 + denovo140384 + denovo140437 + denovo140455 + denovo140504 + denovo140538 +
						denovo140543 + denovo140584 + denovo140585 + denovo140625 + denovo140626 + denovo140653 + denovo140658 + denovo140685 +
						denovo140694 + denovo140749 + denovo140771 + denovo140787 + denovo140796 + denovo140870 + denovo140887 + denovo140895 +
						denovo140897 + denovo140970 + denovo140974 + denovo141042 + denovo141062 + denovo141085 + denovo141091 + denovo141149 +
						denovo141160 + denovo141196 + denovo141201 + denovo141206 + denovo141231 + denovo141256 + denovo141311 + denovo141331 +
						denovo141350 + denovo141361 + denovo141422 + denovo141435 + denovo141440 + denovo141472 + denovo141478 + denovo141506 +
						denovo141593 + denovo141648 + denovo141661 + denovo141678 + denovo141728 + denovo141732 + denovo141811 + denovo141821 +
						denovo141859 + denovo141870 + denovo141892 + denovo141899 + denovo141941 + denovo141947 + denovo142011 +
						denovo142035 + denovo142082 + denovo142106 + denovo142113 + denovo142166 + denovo142221 + denovo142243 + denovo142299 +
						denovo142302 + denovo142350 + denovo142365 + denovo142408 + denovo142446 + denovo142456 + denovo142479 + denovo142489 +
						denovo142514 + denovo142519 + denovo142588 + denovo142593 + denovo142635 + denovo142643 + denovo142645 + denovo142653 +
						denovo142654 + denovo142689 + denovo142690 + denovo142706 + denovo142712 + denovo142766 + denovo142821 + denovo142854 +
						denovo142920 + denovo142943 + denovo142956 + denovo142957 + denovo142964 + denovo142990 + denovo142991 + denovo143008 +
						denovo143014 + denovo143021 + denovo143034 + denovo143114 + denovo143135 + denovo143145 + denovo143160 + denovo143243 +
						denovo143269 + denovo143281 + denovo143301 + denovo143389 + denovo143420 + denovo143431 + denovo143435 + denovo143436 +
						denovo143440 + denovo143535 + denovo143541 + denovo143582 + denovo143593 + denovo143627 + denovo143666 + denovo143694 +
						denovo143734 + denovo143793 + denovo143847 + denovo143887 + denovo143888 + denovo143889 + denovo143902 +
						denovo143914 + denovo143919 + denovo143935 + denovo143946 + denovo143947 + denovo144004 + denovo144010 + denovo144038 +
						denovo144040 + denovo144046 + denovo144053 + denovo144073 + denovo144077 + denovo144119 + denovo144134 + denovo144191 +
						denovo144203 + denovo144206 + denovo144207 + denovo144230 + denovo144259 + denovo144275 + denovo144283 + denovo144294 +
						denovo144348 + denovo144356 + denovo144363 + denovo144442 + denovo144466 + denovo144516 + denovo144519 +
						denovo144525 + denovo144562 + denovo144584 + denovo144587 + denovo144604 + denovo144606 + denovo144610 + denovo144620 +
						denovo144657 + denovo144671 + denovo144697 + denovo144707 + denovo144712 + denovo144777 + denovo144819 + denovo144888 +
						denovo144900 + denovo144918 + denovo144983 + denovo145002 + denovo145075 + denovo145099 + denovo145102 + denovo145105 +
						denovo145118 + denovo145126 + denovo145137 + denovo145204 + denovo145281 + denovo145288 + denovo145289 + denovo145303 +
						denovo145319 + denovo145334 + denovo145392 + denovo145401 + denovo145413 + denovo145432 + denovo145441 + denovo145485 +
						denovo145492 + denovo145517 + denovo145520 + denovo145543 + denovo145562 + denovo145602 + denovo145623 + denovo145642 +
						denovo145673 + denovo145716 + denovo145859 + denovo145888 + denovo145939 + denovo145949 +
						denovo145999 + denovo146060 + denovo146078 + denovo146109 + denovo146129 + denovo146148 + denovo146159 +
						denovo146186 + denovo146190 + denovo146195 + denovo146212 + denovo146255 + denovo146276 + denovo146294 + denovo146306 +
						denovo146310 + denovo146361 + denovo146392 + denovo146399 + denovo146414 + denovo146424 + denovo146436 + denovo146479 +
						denovo146495 + denovo146516 + denovo146543 + denovo146596 + denovo146645 + denovo146697 + denovo146701 +
						denovo146772 + denovo146825 + denovo146840 + denovo146891 + denovo146942 + denovo146947 + denovo146954 + denovo146971 +
						denovo146982 + denovo146984 + denovo147036 + denovo147064 + denovo147069 + denovo147109 + denovo147133 + denovo147139 +
						denovo147158 + denovo147200 + denovo147243 + denovo147247 + denovo147252 + denovo147268 + denovo147281 + denovo147290 +
						denovo147291 + denovo147305 + denovo147314 + denovo147345 + denovo147349 + denovo147364 + denovo147375 + denovo147410 +
						denovo147424 + denovo147435 + denovo147496 + denovo147524 + denovo147559 + denovo147620 + denovo147627 + denovo147644 +
						denovo147669 + denovo147740 + denovo147741 + denovo147755 + denovo147767 + denovo147802 + denovo147804 + denovo147805 +
						denovo147816 + denovo147851 + denovo147944 + denovo147949 + denovo147961 + denovo147973 + denovo147994 + denovo147999 +
						denovo148044 + denovo148048 + denovo148074 + denovo148076 + denovo148085 + denovo148089 + denovo148124 + denovo148134 +
						denovo148171 + denovo148259 + denovo148270 + denovo148296 + denovo148354 + denovo148371 + denovo148374 +
						denovo148376 + denovo148381 + denovo148456 + denovo148459 + denovo148460 + denovo148463 + denovo148468 + denovo148508 +
						denovo148516 + denovo148541 + denovo148559 + denovo148616 + denovo148669 + denovo148691 + denovo148695 + denovo148790 +
						denovo148811 + denovo148821 + denovo148832 + denovo148838 + denovo148868 + denovo148871 + denovo148987 + denovo149045 +
						denovo149050 + denovo149095 + denovo149117 + denovo149163 + denovo149196 + denovo149230 + denovo149241 + denovo149328 +
						denovo149339 + denovo149382 + denovo149401 + denovo149405 + denovo149428 + denovo149436 + denovo149440 + denovo149449 +
						denovo149458 + denovo149462 + denovo149499 + denovo149528 + denovo149530 + denovo149547 + denovo149589 + denovo149594 +
						denovo149612 + denovo149664 + denovo149684 + denovo149731 + denovo149760 + denovo149761 + denovo149855 +
						denovo149875 + denovo149895 + denovo149913 + denovo149917 + denovo149962 + denovo150001 + denovo150020 +
						denovo150027 + denovo150028 + denovo150065 + denovo150067 + denovo150084 + denovo150085 + denovo150088 + denovo150090 +
						denovo150094 + denovo150111 + denovo150123 + denovo150129 + denovo150132 + denovo150180 + denovo150199 + denovo150242 +
						denovo150257 + denovo150291 + denovo150300 + denovo150390 + denovo150408 + denovo150411 + denovo150470 +
						denovo150541 + denovo150565 + denovo150602 + denovo150614 + denovo150619 + denovo150625 + denovo150627 + denovo150632 +
						denovo150679 + denovo150732 + denovo150788 + denovo150813 + denovo150853 + denovo150855 + denovo150868 + denovo150876 +
						denovo150910 + denovo150940 + denovo150950 + denovo150963 + denovo150967 + denovo150982 + denovo151017 + denovo151018 +
						denovo151030 + denovo151050 + denovo151055 + denovo151142 + denovo151150 + denovo151185 + denovo151215 + denovo151224 +
						denovo151242 + denovo151258 + denovo151267 + denovo151322 + denovo151323 + denovo151331 + denovo151402 + denovo151414 +
						denovo151433 + denovo151525 + denovo151568 + denovo151573 + denovo151582 + denovo151642 + denovo151685 + denovo151698 +
						denovo151733 + denovo151743 + denovo151784 + denovo151877 + denovo151886 + denovo151898 + denovo151922 + denovo151924 +
						denovo151944 + denovo151965 + denovo151971 + denovo152049 + denovo152075 + denovo152084 + denovo152115 + denovo152130 +
						denovo152220 + denovo152230 + denovo152283 + denovo152306 + denovo152322 + denovo152331 + denovo152345 +
						denovo152449 + denovo152453 + denovo152483 + denovo152512 + denovo152516 + denovo152519 + denovo152539 + denovo152545 +
						denovo152575 + denovo152620 + denovo152648 + denovo152680 + denovo152688 + denovo152750 + denovo152766 + denovo152781 +
						denovo152786 + denovo152857 + denovo152863 + denovo152888 + denovo152902 + denovo152928 + denovo152936 + denovo152953 +
						denovo152963 + denovo152968 + denovo152979 + denovo153009 + denovo153025 + denovo153031 + denovo153033 + denovo153075 +
						denovo153121 + denovo153138 + denovo153164 + denovo153207 + denovo153219 + denovo153242 + denovo153262 + denovo153332 +
						denovo153471 + denovo153472 + denovo153494 + denovo153526 + denovo153552 + denovo153580 + denovo153591 + denovo153601 +
						denovo153644 + denovo153662 + denovo153746 + denovo153756 + denovo153762 + denovo153764 + denovo153839 +
						denovo153886 + denovo153937 + denovo153999 + denovo154074 + denovo154104 + denovo154118 + denovo154144 +
						denovo154161 + denovo154203 + denovo154209 + denovo154234 + denovo154242 + denovo154283 + denovo154305 + denovo154315 +
						denovo154332 + denovo154347 + denovo154351 + denovo154362 + denovo154373 + denovo154388 + denovo154413 + denovo154504 +
						denovo154506 + denovo154526 + denovo154541 + denovo154548 + denovo154586 + denovo154616 + denovo154671 + denovo154686 +
						denovo154705 + denovo154759 + denovo154760 + denovo154876 + denovo154895 + denovo154917 + denovo154931 +
						denovo154935 + denovo154983 + denovo154988 + denovo155045 + denovo155104 + denovo155125 + denovo155137 + denovo155162 +
						denovo155194 + denovo155196 + denovo155203 + denovo155244 + denovo155264 + denovo155274 + denovo155334 + denovo155400 +
						denovo155469 + denovo155495 + denovo155501 + denovo155511 + denovo155527 + denovo155529 + denovo155563 + denovo155571 +
						denovo155573 + denovo155579 + denovo155596 + denovo155603 + denovo155616 + denovo155618 + denovo155636 +
						denovo155655 + denovo155699 + denovo155710 + denovo155743 + denovo155745 + denovo155842 + denovo155853 + denovo155870 +
						denovo155877 + denovo155929 + denovo155950 + denovo155962 + denovo155973 + denovo155982 + denovo155983 + denovo156065 +
						denovo156084 + denovo156087 + denovo156104 + denovo156108 + denovo156109 + denovo156113 + denovo156149 + denovo156156 +
						denovo156174 + denovo156186 + denovo156218 + denovo156219 + denovo156228 + denovo156281 + denovo156309 + denovo156313 +
						denovo156366 + denovo156377 + denovo156389 + denovo156422 + denovo156423 + denovo156500 + denovo156555 + denovo156557 +
						denovo156606 + denovo156742 + denovo156755 + denovo156756 + denovo156784 + denovo156814 + denovo156843 + denovo156861 +
						denovo156912 + denovo156931 + denovo156944 + denovo156949 + denovo156972 + denovo156975 + denovo157048 + denovo157085 +
						denovo157110 + denovo157135 + denovo157176 + denovo157232 + denovo157245 + denovo157289 + denovo157329 + denovo157347 +
						denovo157379 + denovo157433 + denovo157627 + denovo157657 + denovo157677 + denovo157690 + denovo157694 + denovo157749 +
						denovo157759 + denovo157764 + denovo157770 + denovo157777 + denovo157811 + denovo157827 + denovo157834 + denovo157849 +
						denovo157916 + denovo157945 + denovo157946 + denovo157995 + denovo158048 + denovo158080 + denovo158192 + denovo158202 +
						denovo158210 + denovo158241 + denovo158282 + denovo158349 + denovo158378 + denovo158400 + denovo158405 + denovo158407 +
						denovo158427 + denovo158483 + denovo158499 + denovo158505 + denovo158512 + denovo158513 + denovo158521 + denovo158528 +
						denovo158548 + denovo158571 + denovo158579 + denovo158615 + denovo158620 + denovo158626 + denovo158632 + denovo158653 +
						denovo158658 + denovo158659 + denovo158691 + denovo158692 + denovo158724 + denovo158725 + denovo158735 + denovo158757 +
						denovo158795 + denovo158825 + denovo158833 + denovo158907 + denovo158908 + denovo158913 + denovo159019 + denovo159033 +
						denovo159045 + denovo159057 + denovo159102 + denovo159109 + denovo159167 + denovo159171 + denovo159212 + denovo159258 +
						denovo159274 + denovo159327 + denovo159408 + denovo159420 + denovo159444 + denovo159450 + denovo159468 + denovo159503 +
						denovo159507 + denovo159573 + denovo159619 + denovo159622 + denovo159629 + denovo159639 + denovo159642 + denovo159686 +
						denovo159695 + denovo159718 + denovo159727 + denovo159745 + denovo159765 + denovo159849 + denovo159913 + denovo159918 +
						denovo160177 + denovo160247 + denovo160352 + denovo160386 + denovo160443 + denovo160493 + denovo160514 +
						denovo160530 + denovo160619 + denovo160734 + denovo160747 + denovo160773 + denovo160792 + denovo160797 +
						denovo160805 + denovo160833 + denovo160863 + denovo160887 + denovo160905 + denovo160995 + denovo161000 + denovo161022 +
						denovo161049 + denovo161058 + denovo161130 + denovo161185 + denovo161222 + denovo161229 + denovo161258 + denovo161264 +
						denovo161294 + denovo161305 + denovo161349 + denovo161365 + denovo161369 + denovo161405 + denovo161446 + denovo161454 +
						denovo161540 + denovo161594 + denovo161655 + denovo161676 + denovo161680 + denovo161694 + denovo161702 +
						denovo161707 + denovo161713 + denovo161734 + denovo161745 + denovo161789 + denovo161806 + denovo161809 + denovo161810 +
						denovo161870 + denovo161916 + denovo161942 + denovo161952 + denovo161968 + denovo162000 + denovo162014 + denovo162019 +
						denovo162065 + denovo162067 + denovo162098 + denovo162101 + denovo162151 + denovo162174 + denovo162176 +
						denovo162182 + denovo162215 + denovo162219 + denovo162253 + denovo162293 + denovo162340 + denovo162361 + denovo162390 +
						denovo162499 + denovo162506 + denovo162520 + denovo162524 + denovo162564 + denovo162600 + denovo162645 + denovo162697 +
						denovo162738 + denovo162743 + denovo162754 + denovo162755 + denovo162788 + denovo162799 + denovo162813 + denovo162893 +
						denovo162907 + denovo162922 + denovo162940 + denovo162974 + denovo163045 + denovo163052 + denovo163059 + denovo163072 +
						denovo163087 + denovo163128 + denovo163131 + denovo163159 + denovo163170 + denovo163211 + denovo163260 + denovo163296 +
						denovo163351 + denovo163409 + denovo163484 + denovo163499 + denovo163507 + denovo163530 + denovo163573 + denovo163575 +
						denovo163581 + denovo163662 + denovo163679 + denovo163688 + denovo163698 + denovo163725 + denovo163729 + denovo163810 +
						denovo163823 + denovo163857 + denovo163873 + denovo163927 + denovo163996 + denovo164006 + denovo164014 + denovo164117 +
						denovo164129 + denovo164130 + denovo164178 + denovo164212 + denovo164241 + denovo164263 + denovo164296 + denovo164298 +
						denovo164332 + denovo164377 + denovo164450 + denovo164464 + denovo164520 + denovo164608 + denovo164618 + denovo164625 +
						denovo164647 + denovo164656 + denovo164658 + denovo164669 + denovo164764 + denovo164770 + denovo164782 +
						denovo164808 + denovo164825 + denovo164851 + denovo164855 + denovo164883 + denovo164902 + denovo164936 + denovo164939 +
						denovo164959 + denovo164966 + denovo164975 + denovo164981 + denovo165053 + denovo165084 + denovo165093 + denovo165099 +
						denovo165110 + denovo165123 + denovo165139 + denovo165197 + denovo165214 + denovo165231 + denovo165340 + denovo165343 +
						denovo165348 + denovo165421 + denovo165465 + denovo165482 + denovo165506 + denovo165508 + denovo165548 + denovo165564 +
						denovo165568 + denovo165615 + denovo165621 + denovo165638 + denovo165662 + denovo165730 + denovo165751 + denovo165841 +
						denovo165845 + denovo165857 + denovo165860 + denovo165893 + denovo165899 + denovo165924 + denovo165927 + denovo165950 +
						denovo165951 + denovo165988 + denovo166030 + denovo166148 + denovo166168 + denovo166172 + denovo166188 + denovo166192 +
						denovo166243 + denovo166259 + denovo166310 + denovo166336 + denovo166363 + denovo166394 + denovo166400 + denovo166443 +
						denovo166561 + denovo166578 + denovo166592 + denovo166621 + denovo166662 + denovo166680 + denovo166730 + denovo166747 +
						denovo166783 + denovo166818 + denovo166848 + denovo166865 + denovo166882 + denovo166910 + denovo166936 + denovo166947 +
						denovo166994 + denovo166995 + denovo167037 + denovo167080 + denovo167136 + denovo167187 + denovo167228 +
						denovo167254 + denovo167319 + denovo167326 + denovo167377 + denovo167399 + denovo167402 + denovo167425 + denovo167432 +
						denovo167492 + denovo167505 + denovo167517 + denovo167527 + denovo167535 + denovo167551 + denovo167577 +
						denovo167585 + denovo167617 + denovo167665 + denovo167666 + denovo167718 + denovo167735 + denovo167738 + denovo167744 +
						denovo167751 + denovo167774 + denovo167788 + denovo167790 + denovo167837 + denovo167849 + denovo167863 + denovo167924 +
						denovo167925 + denovo167936 + denovo167973 + denovo167979 + denovo167993 + denovo168068 + denovo168150 +
						denovo168166 + denovo168201 + denovo168204 + denovo168218 + denovo168223 + denovo168225 + denovo168234 + denovo168250 +
						denovo168270 + denovo168294 + denovo168316 + denovo168393 + denovo168430 + denovo168464 + denovo168474 + denovo168497 +
						denovo168538 + denovo168552 + denovo168565 + denovo168576 + denovo168616 + denovo168642 + denovo168652 + denovo168674 +
						denovo168704 + denovo168725 + denovo168767 + denovo168775 + denovo168800 + denovo168847 + denovo168851 + denovo168912 +
						denovo168919 + denovo168986 + denovo168994 + denovo169055 + denovo169118 + denovo169121 + denovo169122 + denovo169146 +
						denovo169147 + denovo169222 + denovo169250 + denovo169311 + denovo169315 + denovo169372 + denovo169411 + denovo169417 +
						denovo169445 + denovo169448 + denovo169456 + denovo169464 + denovo169465 + denovo169473 + denovo169510 + denovo169689 +
						denovo169694 + denovo169695 + denovo169722 + denovo169732 + denovo169737 + denovo169756 + denovo169840 +
						denovo169850 + denovo169853 + denovo169866 + denovo169900 + denovo169922 + denovo169931 + denovo169953 + denovo169954 +
						denovo169956 + denovo169963 + denovo170099 + denovo170112 + denovo170180 + denovo170221 + denovo170253 + denovo170259 +
						denovo170334 + denovo170353 + denovo170385 + denovo170451 + denovo170455 + denovo170474 + denovo170533 + denovo170545 +
						denovo170589 + denovo170602 + denovo170652 + denovo170689 + denovo170775 + denovo170788 + denovo170799 + denovo170879 +
						denovo170932 + denovo170954 + denovo171006 + denovo171034 + denovo171040 + denovo171059 + denovo171063 + denovo171127 +
						denovo171138 + denovo171141 + denovo171143 + denovo171162 + denovo171194 + denovo171222 + denovo171231 + denovo171232 +
						denovo171239 + denovo171263 + denovo171287 + denovo171295 + denovo171321 + denovo171368 + denovo171372 +
						denovo171384 + denovo171449 + denovo171455 + denovo171459 + denovo171503 + denovo171567 + denovo171613 +
						denovo171639 + denovo171656 + denovo171720 + denovo171721 + denovo171733 + denovo171737 + denovo171746 + denovo171786 +
						denovo171824 + denovo171852 + denovo171855 + denovo171881 + denovo171914 + denovo171925 + denovo171998 + denovo172032 +
						denovo172173 + denovo172203 + denovo172255 + denovo172272 + denovo172274 + denovo172301 + denovo172336 + denovo172340 +
						denovo172422 + denovo172477 + denovo172479 + denovo172497 + denovo172541 + denovo172556 + denovo172581 + denovo172589 +
						denovo172607 + denovo172617 + denovo172630 + denovo172691 + denovo172744 + denovo172785 + denovo172796 + denovo172839 +
						denovo172862 + denovo172950 + denovo173049 + denovo173136 + denovo173178 + denovo173199 + denovo173240 +
						denovo173242 + denovo173272 + denovo173299 + denovo173328 + denovo173337 + denovo173339 + denovo173367 + denovo173397 +
						denovo173457 + denovo173475 + denovo173477 + denovo173491 + denovo173496 + denovo173508 + denovo173575 + denovo173591 +
						denovo173609 + denovo173634 + denovo173638 + denovo173691 + denovo173715 + denovo173738 + denovo173766 + denovo173772 +
						denovo173775 + denovo173828 + denovo174012 + denovo174014 + denovo174103 + denovo174123 + denovo174153 + denovo174157 +
						denovo174162 + denovo174198 + denovo174201 + denovo174212 + denovo174216 + denovo174217 + denovo174266 + denovo174271 +
						denovo174278 + denovo174293 + denovo174303 + denovo174385 + denovo174386 + denovo174430 + denovo174435 + denovo174459 +
						denovo174464 + denovo174470 + denovo174497 + denovo174518 + denovo174570 + denovo174613 + denovo174618 + denovo174655 +
						denovo174676 + denovo174691 + denovo174703 + denovo174765 + denovo174773 + denovo174780 + denovo174827 + denovo174846 +
						denovo174852 + denovo174857 + denovo174868 + denovo174886 + denovo174941 + denovo174974 + denovo175001 + denovo175012 +
						denovo175060 + denovo175095 + denovo175122 + denovo175143 + denovo175198 + denovo175200 + denovo175207 + denovo175208 +
						denovo175228 + denovo175317 + denovo175326 + denovo175332 + denovo175360 + denovo175367 + denovo175420 + denovo175461 +
						denovo175494 + denovo175503 + denovo175508 + denovo175538 + denovo175543 + denovo175558 + denovo175590 +
						denovo175645 + denovo175709 + denovo175713 + denovo175732 + denovo175733 + denovo175763 + denovo175777 + denovo175778 +
						denovo175939 + denovo175972 + denovo175974 + denovo175978 + denovo175983 + denovo175995 + denovo175996 + denovo176025 +
						denovo176046 + denovo176094 + denovo176100 + denovo176172 + denovo176206 + denovo176235 + denovo176238 + denovo176255 +
						denovo176337 + denovo176364 + denovo176412 + denovo176416 + denovo176418 + denovo176463 + denovo176467 + denovo176471 +
						denovo176503 + denovo176547 + denovo176561 + denovo176617 + denovo176634 + denovo176682 + denovo176685 +
						denovo176718 + denovo176757 + denovo176760 + denovo176765 + denovo176778 + denovo176800 + denovo176852 + denovo176861 +
						denovo176903 + denovo176966 + denovo176981 + denovo176988 + denovo176995 + denovo177020 + denovo177055 +
						denovo177096 + denovo177145 + denovo177198 + denovo177236 + denovo177260 + denovo177271 + denovo177279 + denovo177316 +
						denovo177339 + denovo177354 + denovo177366 + denovo177393 + denovo177406 + denovo177411 + denovo177500 + denovo177520 +
						denovo177557 + denovo177559 + denovo177590 + denovo177599 + denovo177707 + denovo177723 + denovo177769 +
						denovo177814 + denovo177830 + denovo177839 + denovo177843 + denovo177863 + denovo177898 + denovo177917 + denovo177929 +
						denovo177947 + denovo177970 + denovo178002 + denovo178062 + denovo178064 + denovo178070 + denovo178109 + denovo178137 +
						denovo178152 + denovo178161 + denovo178187 + denovo178202 + denovo178205 + denovo178212 + denovo178214 + denovo178295 +
						denovo178298 + denovo178300 + denovo178314 + denovo178326 + denovo178359 + denovo178397 + denovo178409 + denovo178423 +
						denovo178425 + denovo178501 + denovo178510 + denovo178514 + denovo178545 + denovo178546 + denovo178625 + denovo178646 +
						denovo178672 + denovo178716 + denovo178735 + denovo178829 + denovo178831 + denovo178903 + denovo178927 + denovo178990 +
						denovo179025 + denovo179027 + denovo179048 + denovo179086 + denovo179144 + denovo179159 + denovo179194 + denovo179210 +
						denovo179228 + denovo179243 + denovo179317 + denovo179330 + denovo179334 + denovo179341 + denovo179362 + denovo179384 +
						denovo179413 + denovo179416 + denovo179443 + denovo179481 + denovo179493 + denovo179499 + denovo179502 + denovo179513 +
						denovo179522 + denovo179523 + denovo179535 + denovo179538 + denovo179578 + denovo179583 + denovo179591 + denovo179593 +
						denovo179594 + denovo179613 + denovo179645 + denovo179693 + denovo179766 + denovo179785 + denovo179795 + denovo179798 +
						denovo179825 + denovo179835 + denovo179846 + denovo179873 + denovo179874 + denovo179907 + denovo179911 + denovo179915 +
						denovo179922 + denovo179925 + denovo179936 + denovo179999 + denovo180017 + denovo180072 + denovo180090 + denovo180154 +
						denovo180163 + denovo180171 + denovo180187 + denovo180191 + denovo180202 + denovo180218 + denovo180228 +
						denovo180244 + denovo180257 + denovo180276 + denovo180369 + denovo180385 + denovo180391 + denovo180397 + denovo180411 +
						denovo180417 + denovo180421 + denovo180444 + denovo180475 + denovo180487 + denovo180532 + denovo180618 + denovo180672 +
						denovo180767 + denovo180768 + denovo180788 + denovo180822 + denovo180853 + denovo180872 + denovo180911 +
						denovo180934 + denovo181058 + denovo181102 + denovo181118 + denovo181155 + denovo181192 + denovo181197 + denovo181206 +
						denovo181235 + denovo181244 + denovo181269 + denovo181288 + denovo181305 + denovo181340 + denovo181392 + denovo181455 +
						denovo181502 + denovo181530 + denovo181536 + denovo181545 + denovo181606 + denovo181636 + denovo181773 + denovo181783 +
						denovo181809 + denovo181834 + denovo181835 + denovo181874 + denovo181883 + denovo181892 + denovo181903 + denovo181909 +
						denovo181928 + denovo181998 + denovo182030 + denovo182057 + denovo182081 + denovo182086 + denovo182122 +
						denovo182172 + denovo182236 + denovo182267 + denovo182297 + denovo182322 + denovo182328 + denovo182351 + denovo182352 +
						denovo182356 + denovo182362 + denovo182366 + denovo182390 + denovo182403 + denovo182451 + denovo182496 + denovo182558 +
						denovo182593 + denovo182596 + denovo182638 + denovo182645 + denovo182743 + denovo182748 +
						denovo182756 + denovo182797 + denovo182820 + denovo182887 + denovo182895 + denovo182898 + denovo182906 + denovo182908 +
						denovo182909 + denovo182927 + denovo182959 + denovo182964 + denovo182971 + denovo182997 + denovo183026 + denovo183067 +
						denovo183072 + denovo183090 + denovo183095 + denovo183122 + denovo183159 + denovo183183 + denovo183190 + denovo183222 +
						denovo183226 + denovo183242 + denovo183248 + denovo183249 + denovo183261 + denovo183270 + denovo183305 + denovo183318 +
						denovo183323 + denovo183335 + denovo183355 + denovo183365 + denovo183385 + denovo183412 + denovo183440 + denovo183475 +
						denovo183490 + denovo183516 + denovo183562 + denovo183627 + denovo183650 + denovo183654 + denovo183699 + denovo183711 +
						denovo183723 + denovo183736 + denovo183748 + denovo183749 + denovo183822 + denovo183907 + denovo183916 +
						denovo183969 + denovo183975 + denovo184010 + denovo184034 + denovo184053 + denovo184063 + denovo184123 +
						denovo184124 + denovo184126 + denovo184238 + denovo184259 + denovo184261 + denovo184279 + denovo184296 + denovo184320 +
						denovo184371 + denovo184418 + denovo184460 + denovo184553 + denovo184560 + denovo184561 + denovo184568 +
						denovo184600 + denovo184607 + denovo184620 + denovo184628 + denovo184629 + denovo184637 + denovo184645 + denovo184674 +
						denovo184692 + denovo184700 + denovo184714 + denovo184733 + denovo184739 + denovo184751 + denovo184757 + denovo184789 +
						denovo184819 + denovo184827 + denovo184841 + denovo184861 + denovo184862 + denovo184865 + denovo184880 + denovo184924 +
						denovo184930 + denovo184999 + denovo185008 + denovo185025 + denovo185048 + denovo185056 + denovo185061 + denovo185071 +
						denovo185112 + denovo185164 + denovo185169 + denovo185217 + denovo185222 + denovo185240 + denovo185280 +
						denovo185331 + denovo185388 + denovo185422 + denovo185429 + denovo185432 + denovo185457 + denovo185473 + denovo185544 +
						denovo185545 + denovo185583 + denovo185608 + denovo185620 + denovo185635 + denovo185668 + denovo185808 + denovo185844 +
						denovo185845 + denovo185933 + denovo185938 + denovo185946 + denovo185965 + denovo185979 + denovo185980 +
						denovo186021 + denovo186022 + denovo186024 + denovo186051 + denovo186085 + denovo186110 + denovo186125 + denovo186145 +
						denovo186223 + denovo186225 + denovo186260 + denovo186285 + denovo186330 + denovo186353 + denovo186376 + denovo186388 +
						denovo186453 + denovo186462 + denovo186499 + denovo186502 + denovo186519 + denovo186534 + denovo186576 + denovo186652 +
						denovo186669 + denovo186689 + denovo186740 + denovo186756 + denovo186799 + denovo186815 + denovo186818 + denovo186855 +
						denovo186895 + denovo186916 + denovo186938 + denovo186974 + denovo186975 + denovo186994 + denovo186995 + denovo187013 +
						denovo187033 + denovo187068 + denovo187095 + denovo187154 + denovo187172 + denovo187252 + denovo187262 + denovo187278 +
						denovo187282 + denovo187365 + denovo187427 + denovo187440 + denovo187452 + denovo187490 + denovo187527 + denovo187528 +
						denovo187554 + denovo187603 + denovo187626 + denovo187685 + denovo187785 + denovo187789 + denovo187812 + denovo187914 +
						denovo187919 + denovo188002 + denovo188024 + denovo188060 + denovo188067 + denovo188091 + denovo188097 + denovo188101 +
						denovo188137 + denovo188159 + denovo188183 + denovo188234 + denovo188262 + denovo188300 + denovo188328 + denovo188376 +
						denovo188438 + denovo188471 + denovo188507 + denovo188526 + denovo188560 + denovo188576 + denovo188578 + denovo188587 +
						denovo188616 + denovo188639 + denovo188791 + denovo188840 + denovo188843 + denovo188854 + denovo188858 + denovo188895 +
						denovo188957 + denovo189047 + denovo189056 + denovo189091 + denovo189111 + denovo189154 + denovo189190 +
						denovo189192 + denovo189207 + denovo189218 + denovo189230 + denovo189272 + denovo189284 + denovo189354 + denovo189368 +
						denovo189388 + denovo189398 + denovo189430 + denovo189465 + denovo189475 + denovo189529 + denovo189594 + denovo189629 +
						denovo189654 + denovo189660 + denovo189665 + denovo189666 + denovo189671 + denovo189695 + denovo189762 + denovo189765 +
						denovo189833 + denovo189853 + denovo189854 + denovo189864 + denovo189893 + denovo189905 + denovo189959 + denovo189981 +
						denovo189993 + denovo190018 + denovo190047 + denovo190068 + denovo190071 + denovo190106 + denovo190157 + denovo190183 +
						denovo190205 + denovo190213 + denovo190269 + denovo190322 + denovo190324 + denovo190383 + denovo190400 +
						denovo190403 + denovo190409 + denovo190438 + denovo190461 + denovo190513 + denovo190568 + denovo190589 + denovo190599 +
						denovo190611 + denovo190614 + denovo190651 + denovo190668 + denovo190669 + denovo190699 + denovo190701 + denovo190714 +
						denovo190776 + denovo190785 + denovo190810 + denovo190831 + denovo190876 + denovo190883 + denovo190886 + denovo190897 +
						denovo190898 + denovo190943 + denovo190963 + denovo191057 + denovo191061 + denovo191084 + denovo191087 + denovo191146 +
						denovo191195 + denovo191196 + denovo191231 + denovo191237 + denovo191247 + denovo191282 + denovo191314 + denovo191336 +
						denovo191369 + denovo191388 + denovo191407 + denovo191424 + denovo191483 + denovo191537 + denovo191539 +
						denovo191547 + denovo191603 + denovo191659 + denovo191727 + denovo191731 + denovo191771 + denovo191825 + denovo191831 +
						denovo191848 + denovo191852 + denovo191874 + denovo191925 + denovo191948 + denovo191982 + denovo191995 + denovo192007 +
						denovo192020 + denovo192191 + denovo192231 + denovo192268 + denovo192304 + denovo192319 + denovo192326 + denovo192342 +
						denovo192427 + denovo192443 + denovo192452 + denovo192460 + denovo192482 + denovo192483 + denovo192580 +
						denovo192593 + denovo192615 + denovo192620 + denovo192637 + denovo192672 + denovo192718 + denovo192733 + denovo192739 +
						denovo192827 + denovo192858 + denovo192884 + denovo192965 + denovo192976 + denovo193072 + denovo193121 +
						denovo193191 + denovo193197 + denovo193205 + denovo193220 + denovo193269 + denovo193308 + denovo193319 + denovo193328 +
						denovo193333 + denovo193371 + denovo193392 + denovo193467 + denovo193471 + denovo193512 + denovo193542 + denovo193544 +
						denovo193581 + denovo193620 + denovo193627 + denovo193642 + denovo193664 + denovo193714 + denovo193723 + denovo193739 +
						denovo193740 + denovo193744 + denovo193791 + denovo193807 + denovo193817 + denovo193823 + denovo193853 + denovo193863 +
						denovo193918 + denovo193964 + denovo193971 + denovo193990 + denovo193999 + denovo194010 + denovo194014 +
						denovo194024 + denovo194035 + denovo194038 + denovo194046 + denovo194058 + denovo194077 + denovo194079 + denovo194096 +
						denovo194133 + denovo194171 + denovo194231 + denovo194298 + denovo194455 + denovo194524 + denovo194537 +
						denovo194549 + denovo194594 + denovo194628 + denovo194688 + denovo194709 + denovo194713 + denovo194724 + denovo194725 +
						denovo194764 + denovo194768 + denovo194769 + denovo194774 + denovo194776 + denovo194777 + denovo194808 + denovo194826 +
						denovo194888 + denovo194949 + denovo194970 + denovo195020 + denovo195087 + denovo195091 + denovo195103 + denovo195123 +
						denovo195174 + denovo195207 + denovo195222 + denovo195224 + denovo195251 + denovo195270 + denovo195272 + denovo195284 +
						denovo195385 + denovo195457 + denovo195470 + denovo195502 + denovo195526 + denovo195527 + denovo195562 + denovo195573 +
						denovo195632 + denovo195634 + denovo195639 + denovo195698 + denovo195700 + denovo195705 + denovo195731 + denovo195778 +
						denovo195791 + denovo195804 + denovo195819 + denovo195858 + denovo195871 + denovo195875 + denovo195885 + denovo195888 +
						denovo195918 + denovo195919 + denovo195944 + denovo195973 + denovo196007 + denovo196076 + denovo196089 + denovo196124 +
						denovo196129 + denovo196153 + denovo196161 + denovo196168 + denovo196207 + denovo196208 + denovo196237 + denovo196244 +
						denovo196247 + denovo196253 + denovo196265 + denovo196270 + denovo196272 + denovo196365 + denovo196372 + denovo196387 +
						denovo196466 + denovo196500 + denovo196575 + denovo196677 + denovo196693 + denovo196695 + denovo196696 +
						denovo196736 + denovo196737 + denovo196761 + denovo196791 + denovo196813 + denovo196824 + denovo196857 + denovo196863 +
						denovo196865 + denovo196899 + denovo196933 + denovo196950 + denovo196989 + denovo196998 + denovo196999 + denovo197012 +
						denovo197032 + denovo197033 + denovo197039 + denovo197074 + denovo197078 + denovo197101 + denovo197132 + denovo197136 +
						denovo197137 + denovo197189 + denovo197196 + denovo197199 + denovo197220 + denovo197241 + denovo197248 + denovo197250 +
						denovo197255 + denovo197371 + denovo197383 + denovo197400 + denovo197421 + denovo197460 +
						denovo197461 + denovo197524 + denovo197578 + denovo197599 + denovo197600 + denovo197612 + denovo197634 + denovo197649 +
						denovo197650 + denovo197681 + denovo197682 + denovo197715 + denovo197736 + denovo197738 + denovo197775 + denovo197788 +
						denovo197839 + denovo197845 + denovo197857 + denovo197877 + denovo197880 + denovo197904 + denovo197933 + denovo197964 +
						denovo197988 + denovo197994 + denovo198008 + denovo198060 + denovo198081 + denovo198098 + denovo198102 +
						denovo198124 + denovo198168 + denovo198213 + denovo198318 + denovo198392 + denovo198414 + denovo198471 + denovo198484 +
						denovo198494 + denovo198501 + denovo198563 + denovo198570 + denovo198587 + denovo198613 + denovo198614 + denovo198624 +
						denovo198648 + denovo198651 + denovo198663 + denovo198680 + denovo198685 + denovo198696 + denovo198726 + denovo198780 +
						denovo198781 + denovo198786 + denovo198790 + denovo198814 + denovo198858 + denovo198878 + denovo198897 + denovo198916 +
						denovo198928 + denovo198941 + denovo199007 + denovo199032 + denovo199052 + denovo199112 +
						denovo199116 + denovo199127 + denovo199144 + denovo199172 + denovo199184 + denovo199242 + denovo199287 + denovo199342 +
						denovo199344 + denovo199360 + denovo199373 + denovo199404 + denovo199452 + denovo199454 + denovo199463 +
						denovo199515 + denovo199521 + denovo199546 + denovo199605 + denovo199611 + denovo199613 + denovo199715 + denovo199771 +
						denovo199783 + denovo199788 + denovo199807 + denovo199828 + denovo199854 + denovo199864 + denovo199881 + denovo199901 +
						denovo199912 + denovo199944 + denovo199955 + denovo200029 + denovo200044 + denovo200047 + denovo200072 + denovo200080 +
						denovo200091 + denovo200094 + denovo200104 + denovo200115 + denovo200174 + denovo200183 + denovo200189 + denovo200222 +
						denovo200258 + denovo200262 + denovo200270 + denovo200280 + denovo200310 + denovo200343 + denovo200351 + denovo200408 +
						denovo200422 + denovo200423 + denovo200439 + denovo200445 + denovo200447 + denovo200495 + denovo200510 + denovo200531 +
						denovo200543 + denovo200597 + denovo200602 + denovo200608 + denovo200617 + denovo200635 + denovo200722 +
						denovo200767 + denovo200792 + denovo200839 + denovo200872 + denovo200893 + denovo200902 + denovo200905 + denovo200935 +
						denovo200947 + denovo200977 + denovo201015 + denovo201089 + denovo201109 + denovo201153 + denovo201162 + denovo201214 +
						denovo201266 + denovo201305 + denovo201328 + denovo201331 + denovo201370 + denovo201375 + denovo201376 + denovo201385 +
						denovo201419 + denovo201432 + denovo201492 + denovo201516 + denovo201534 + denovo201563 + denovo201571 + denovo201605 +
						denovo201622 + denovo201629 + denovo201748 + denovo201796 + denovo201803 + denovo202026 +
						denovo202032 + denovo202097 + denovo202150 + denovo202202 + denovo202219 + denovo202282 + denovo202396 + denovo202408 +
						denovo202444 + denovo202448 + denovo202479 + denovo202488 + denovo202524 + denovo202544 + denovo202551 +
						denovo202567 + denovo202613 + denovo202654 + denovo202679 + denovo202704 + denovo202723 + denovo202731 + denovo202754 +
						denovo202846 + denovo202867 + denovo202871 + denovo202960 + denovo203003 + denovo203007 + denovo203070 +
						denovo203079 + denovo203083 + denovo203117 + denovo203136 + denovo203151 + denovo203169 + denovo203172 + denovo203187 +
						denovo203209 + denovo203257 + denovo203273 + denovo203276 + denovo203284 + denovo203295 + denovo203305 +
						denovo203317 + denovo203320 + denovo203348 + denovo203373 + denovo203426 + denovo203431 + denovo203514 + denovo203590 +
						denovo203634 + denovo203671 + denovo203812 + denovo203891 + denovo203906 + denovo203909 + denovo203912 +
						denovo203965 + denovo204036 + denovo204045 + denovo204077 + denovo204089 + denovo204166 + denovo204189 + denovo204225 +
						denovo204237 + denovo204289 + denovo204292 + denovo204319 + denovo204396 + denovo204409 + denovo204419 + denovo204425 +
						denovo204457 + denovo204471 + denovo204475 + denovo204520 + denovo204571 + denovo204590 + denovo204610 + denovo204657 +
						denovo204694 + denovo204703 + denovo204771 + denovo204775 + denovo204802 + denovo204831 + denovo204836 +
						denovo204840 + denovo204856 + denovo204882 + denovo204923 + denovo204965 + denovo204989 + denovo205029 + denovo205062 +
						denovo205113 + denovo205127 + denovo205183 + denovo205192 + denovo205239 + denovo205242 +
						denovo205246 + denovo205252 + denovo205261 + denovo205263 + denovo205312 + denovo205354 + denovo205367 +
						denovo205491 + denovo205496 + denovo205519 + denovo205524 + denovo205530 + denovo205555 + denovo205564 + denovo205568 +
						denovo205580 + denovo205585 + denovo205594 + denovo205605 + denovo205653 + denovo205664 + denovo205688 + denovo205696 +
						denovo205703 + denovo205768 + denovo205769 + denovo205780 + denovo205796 + denovo205799 + denovo205836 + denovo205841 +
						denovo205868 + denovo205892 + denovo205894 + denovo205902 + denovo205924 + denovo205939 + denovo205946 + denovo205990 +
						denovo206025 + denovo206071 + denovo206078 + denovo206094 + denovo206100 + denovo206171 + denovo206216 + denovo206242 +
						denovo206252 + denovo206261 + denovo206315 + denovo206332 + denovo206345 + denovo206356 + denovo206360 + denovo206367 +
						denovo206454 + denovo206487 + denovo206494 + denovo206509 + denovo206542 + denovo206546 + denovo206554 + denovo206571 +
						denovo206576 + denovo206674 + denovo206686 + denovo206727 + denovo206743 + denovo206749 + denovo206816 + denovo206831 +
						denovo206896 + denovo206900 + denovo206942 + denovo206957 + denovo207002 + denovo207059 + denovo207076 +
						denovo207080 + denovo207103 + denovo207118 + denovo207137 + denovo207169 + denovo207176 + denovo207182 + denovo207234 +
						denovo207240 + denovo207257 + denovo207259 + denovo207392 + denovo207412 + denovo207422 + denovo207449 + denovo207501 +
						denovo207544 + denovo207599 + denovo207624 + denovo207661 + denovo207674 + denovo207678 + denovo207691 + denovo207699 +
						denovo207725 + denovo207736 + denovo207778 + denovo207793 + denovo207817 + denovo207855 + denovo207868 + denovo207878 +
						denovo207884 + denovo207909 + denovo207910 + denovo207917 + denovo207938 + denovo207976 + denovo207978 + denovo208008 +
						denovo208009 + denovo208023 + denovo208047 + denovo208096 + denovo208124 + denovo208125 + denovo208135 + denovo208224 +
						denovo208234 + denovo208357 + denovo208378 + denovo208405 + denovo208437 + denovo208449 + denovo208457 + denovo208459 +
						denovo208483 + denovo208500 + denovo208560 + denovo208588 + denovo208644 + denovo208651 + denovo208680 + denovo208700 +
						denovo208733 + denovo208743 + denovo208819 + denovo208846 + denovo208849 + denovo208867 + denovo208902 + denovo208956 +
						denovo208981 + denovo208995 + denovo209005 + denovo209013 + denovo209014 + denovo209054 + denovo209084 + denovo209106 +
						denovo209116 + denovo209124 + denovo209175 + denovo209191 + denovo209325 + denovo209331 + denovo209346 + denovo209353 +
						denovo209405 + denovo209447 + denovo209497 + denovo209501 + denovo209513 + denovo209526 + denovo209530 + denovo209592 +
						denovo209618 + denovo209632 + denovo209671 + denovo209714 + denovo209758 + denovo209763 + denovo209774 + denovo209781 +
						denovo209782 + denovo209841 + denovo209855 + denovo209874 + denovo209888 + denovo209929 + denovo210026 + denovo210029 +
						denovo210048 + denovo210051 + denovo210093 + denovo210166 + denovo210189 + denovo210190 + denovo210200 + denovo210206 +
						denovo210236 + denovo210261 + denovo210263 + denovo210308 + denovo210318 + denovo210352 + denovo210370 + denovo210375 +
						denovo210395 + denovo210413 + denovo210423 + denovo210440 + denovo210454 + denovo210456 + denovo210467 + denovo210473 +
						denovo210530 + denovo210536 + denovo210564 + denovo210575 + denovo210586 + denovo210587 + denovo210610 + denovo210630 +
						denovo210761 + denovo210779 + denovo210791 + denovo210855 + denovo210867 + denovo210883 + denovo210892 +
						denovo210936 + denovo210973 + denovo210989 + denovo211012 + denovo211014 + denovo211024 + denovo211032 + denovo211053 +
						denovo211116 + denovo211129 + denovo211137 + denovo211149 + denovo211157 + denovo211194 + denovo211208 + denovo211225 +
						denovo211263 + denovo211275 + denovo211287 + denovo211340 + denovo211393 + denovo211410 + denovo211426 +
						denovo211456 + denovo211459 + denovo211497 + denovo211500 + denovo211516 + denovo211519 + denovo211572 +
						denovo211599 + denovo211611 + denovo211650 + denovo211684 + denovo211745 + denovo211780 + denovo211839 + denovo211840 +
						denovo211862 + denovo211899 + denovo211900 + denovo211943 + denovo211955 + denovo211979 + denovo212022 + denovo212043 +
						denovo212062 + denovo212249 + denovo212253 + denovo212268 + denovo212281 + denovo212290 + denovo212299 + denovo212302 +
						denovo212321 + denovo212377 + denovo212394 + denovo212408 + denovo212419 + denovo212423 + denovo212425 +
						denovo212443 + denovo212467 + denovo212488 + denovo212524 + denovo212528 + denovo212532 + denovo212588 + denovo212607 +
						denovo212612 + denovo212656 + denovo212720 + denovo212728 + denovo212732 + denovo212768 + denovo212821 + denovo212830 +
						denovo212887 + denovo212888 + denovo212948 + denovo212951 + denovo212967 + denovo212980 + denovo212995 + denovo213006 +
						denovo213053 + denovo213188 + denovo213197 + denovo213205 + denovo213213 + denovo213222 + denovo213304 +
						denovo213311 + denovo213372 + denovo213380 + denovo213399 + denovo213416 + denovo213485 + denovo213508 + denovo213509 +
						denovo213537 + denovo213562 + denovo213576 + denovo213584 + denovo213592 + denovo213618 + denovo213625 + denovo213627 +
						denovo213708 + denovo213715 + denovo213719 + denovo213733 + denovo213817 + denovo213826 + denovo213850 + denovo213933 +
						denovo213955 + denovo213986 + denovo214061 + denovo214062 + denovo214135 + denovo214150 + denovo214153 +
						denovo214175 + denovo214197 + denovo214199 + denovo214222 + denovo214241 + denovo214299 + denovo214354 + denovo214373 +
						denovo214397 + denovo214466 + denovo214468 + denovo214473 + denovo214507 + denovo214508 + denovo214568 + denovo214579 +
						denovo214584 + denovo214660 + denovo214668 + denovo214685 + denovo214704 + denovo214719 + denovo214741 + denovo214778 +
						denovo214829 + denovo214837 + denovo214845 + denovo214853 + denovo214883 + denovo214885 + denovo214898 + denovo214912 +
						denovo214920 + denovo214962 + denovo214984 + denovo215052 + denovo215058 + denovo215072 + denovo215091 + denovo215095 +
						denovo215146 + denovo215213 + denovo215229 + denovo215257 + denovo215271 + denovo215297 + denovo215299 + denovo215306 +
						denovo215321 + denovo215336 + denovo215364 + denovo215371 + denovo215388 + denovo215454 + denovo215494 +
						denovo215513 + denovo215553 + denovo215602 + denovo215701 + denovo215712 + denovo215722 + denovo215731 + denovo215744 +
						denovo215755 + denovo215762 + denovo215772 + denovo215777 + denovo215797 + denovo215805 + denovo215832 + denovo215844 +
						denovo215863 + denovo215866 + denovo215903 + denovo215913 + denovo215937 + denovo215959 + denovo215983 + denovo216014 +
						denovo216015 + denovo216017 + denovo216030 + denovo216074 + denovo216075 + denovo216098 + denovo216185 +
						denovo216195 + denovo216238 + denovo216241 + denovo216249 + denovo216270 + denovo216273 + denovo216291 + denovo216301 +
						denovo216350 + denovo216366 + denovo216373 + denovo216385 + denovo216407 + denovo216461 + denovo216469 + denovo216471 +
						denovo216474 + denovo216475 + denovo216491 + denovo216498 + denovo216502 + denovo216508 + denovo216517 + denovo216528 +
						denovo216542 + denovo216582 + denovo216605 + denovo216626 + denovo216640 + denovo216699 + denovo216766 + denovo216788 +
						denovo216804 + denovo216869 + denovo216886 + denovo216902 + denovo216974 + denovo216975 + denovo216981 + denovo216998 +
						denovo217073 + denovo217075 + denovo217083 + denovo217142 + denovo217151 + denovo217175 + denovo217193 + denovo217265 +
						denovo217296 + denovo217312 + denovo217355 + denovo217356 + denovo217384 + denovo217403 + denovo217422 + denovo217443 +
						denovo217478 + denovo217483 + denovo217496 + denovo217506 + denovo217513 + denovo217589 + denovo217615 + denovo217672 +
						denovo217704 + denovo217711 + denovo217746 + denovo217778 + denovo217793 + denovo217811 + denovo217866 +
						denovo217906 + denovo217925 + denovo217995 + denovo218024 + denovo218046 + denovo218139 + denovo218141 + denovo218150 +
						denovo218159 + denovo218204 + denovo218229 + denovo218235 + denovo218282 + denovo218284 + denovo218304 +
						denovo218335 + denovo218370 + denovo218446 + denovo218473 + denovo218635 + denovo218671 + denovo218727 + denovo218758 +
						denovo218780 + denovo218790 + denovo218795 + denovo218838 + denovo218861 + denovo218875 + denovo218879 +
						denovo218895 + denovo218951 + denovo218996 + denovo219000 + denovo219004 + denovo219014 + denovo219015 +
						denovo219040 + denovo219051 + denovo219111 + denovo219219 + denovo219236 + denovo219253 + denovo219255 +
						denovo219272 + denovo219285 + denovo219341 + denovo219344 + denovo219378 + denovo219433 + denovo219453 + denovo219459 +
						denovo219463 + denovo219466 + denovo219507 + denovo219521 + denovo219537 + denovo219545 + denovo219571 + denovo219607 +
						denovo219636 + denovo219637 + denovo219704 + denovo219705 + denovo219721 + denovo219722 + denovo219781 + denovo219811 +
						denovo219823 + denovo219848 + denovo219851 + denovo219872 + denovo219965 + denovo220096 + denovo220126 +
						denovo220166 + denovo220195 + denovo220233 + denovo220258 + denovo220266 + denovo220325 + denovo220430 + denovo220476 +
						denovo220486 + denovo220535 + denovo220573 + denovo220585 + denovo220588 + denovo220625 + denovo220652 + denovo220654 +
						denovo220662 + denovo220667 + denovo220689 + denovo220768 + denovo220771 + denovo220810 + denovo220818 + denovo220822 +
						denovo220844 + denovo220855 + denovo220866 + denovo220879 + denovo220941 + denovo220956 + denovo220968 + denovo221034 +
						denovo221039 + denovo221058 + denovo221061 + denovo221068 + denovo221086 + denovo221099 + denovo221129 + denovo221140 +
						denovo221177 + denovo221182 + denovo221188 + denovo221216 + denovo221232 + denovo221239 + denovo221280 + denovo221293 +
						denovo221311 + denovo221337 + denovo221340 + denovo221342 + denovo221347 + denovo221371 + denovo221377 + denovo221405 +
						denovo221423 + denovo221486 + denovo221591 + denovo221653 + denovo221658 + denovo221668 + denovo221707 + denovo221725 +
						denovo221761 + denovo221788 + denovo221801 + denovo221813 + denovo221817 + denovo221870 + denovo221930 + denovo222007 +
						denovo222040 + denovo222049 + denovo222064 + denovo222117 + denovo222132 + denovo222143 + denovo222153 +
						denovo222160 + denovo222270 + denovo222357 + denovo222362 + denovo222366 + denovo222435 + denovo222447 +
						denovo222483 + denovo222528 + denovo222558 + denovo222560 + denovo222579 + denovo222623 + denovo222683 + denovo222703 +
						denovo222736 + denovo222739 + denovo222751 + denovo222783 + denovo222867 + denovo222882 + denovo222961 + denovo223014 +
						denovo223039 + denovo223069 + denovo223103 + denovo223120 + denovo223130 + denovo223135 + denovo223180 + denovo223183 +
						denovo223262 + denovo223321 + denovo223366 + denovo223381 + denovo223395 + denovo223401 + denovo223418 + denovo223477 +
						denovo223522 + denovo223528 + denovo223568 + denovo223606 + denovo223624 + denovo223644 + denovo223646 + denovo223686 +
						denovo223720 + denovo223723 + denovo223757 + denovo223778 + denovo223801 + denovo223815 + denovo223821 + denovo223837 +
						denovo223841 + denovo223842 + denovo223949 + denovo224012 + denovo224039 + denovo224088 + denovo224142 +
						denovo224166 + denovo224241 + denovo224242 + denovo224254 + denovo224261 + denovo224283 + denovo224291 + denovo224299 +
						denovo224308 + denovo224358 + denovo224369 + denovo224387 + denovo224404 + denovo224416 + denovo224420 + denovo224439 +
						denovo224484 + denovo224487 + denovo224552 + denovo224604 + denovo224609 + denovo224676 + denovo224681 +
						denovo224806 + denovo224812 + denovo224837 + denovo224839 + denovo224890 + denovo224894 + denovo224898 + denovo224942 +
						denovo224978 + denovo225028 + denovo225042 + denovo225131 + denovo225200 + denovo225227 + denovo225228 + denovo225327 +
						denovo225351 + denovo225358 + denovo225400 + denovo225431 + denovo225439 + denovo225444 + denovo225635 + denovo225652 +
						denovo225658 + denovo225673 + denovo225693 + denovo225711 + denovo225713 + denovo225720 + denovo225723 + denovo225725 +
						denovo225803 + denovo225804 + denovo225808 + denovo225855 + denovo225859 + denovo225890 + denovo225937 + denovo225956 +
						denovo225964 + denovo225990 + denovo226021 + denovo226026 + denovo226093 + denovo226098 + denovo226133 + denovo226135, data = presence_absence_2016_FFG, CV = TRUE)

DFA_predict <- table(presence_absence_2016_FFG$FFG, DFA_fit_JK$class)

# Measure Total Correctly Predicted Individuals by Functional Feeding Group
diag(prop.table(DFA_predict, 1))


## 10) Presence & Absence - Discriminant Function Analysis for Traditional Feeding Habit

# Discriminant Function Analysis Using Jacknifed Prediction
presence_absence_2016_feeding_habit <- read.csv("/data/shawn/2016_Multivariate_Analysis/presence_absence_2016_feeding_habit.csv")
presence_absence_2016_feeding_habit <- as.data.frame(presence_absence_2016_feeding_habit)
row.names(presence_absence_2016_feeding_habit) <- presence_absence_2016_feeding_habit$SampleID
presence_absence_2016_feeding_habit <- presence_absence_2016_feeding_habit[,-1]

DFA_fit_JK <- lda(FeedingHabit ~ denovo36 + denovo122 + denovo133 + denovo140 + denovo162 + denovo216 + denovo255 + 
								 denovo332 + denovo367 + denovo391 + denovo402 + denovo435 + denovo448 + denovo465 + denovo478 +
                                 denovo511 + denovo535 + denovo544 + denovo567 + denovo584 + denovo624 + denovo680 + denovo695 +
                                 denovo707 + denovo728 + denovo730 + denovo747 + denovo837 + denovo922 + denovo924 + denovo946 +
                                 denovo964 + denovo991 + denovo994 + denovo1013 + denovo1027 + denovo1041 + denovo1047 + denovo1048 +
                                 denovo1054 + denovo1064 + denovo1077 + denovo1106 + denovo1118 + denovo1142 + denovo1182 + denovo1211 +
                                 denovo1223 + denovo1229 + denovo1292 + denovo1332 + denovo1335 + denovo1339 + denovo1365 + denovo1386 +
                                 denovo1407 + denovo1408 + denovo1487 + denovo1524 + denovo1533 + denovo1540 + denovo1545 + denovo1553 +
                                 denovo1621 + denovo1625 + denovo1632 + denovo1642 + denovo1664 + denovo1669 + denovo1707 + denovo1767 +
                                 denovo1773 + denovo1836 + denovo1839 + denovo1887 + denovo1941 + denovo1944 + denovo1989 + denovo1992 +
                                 denovo2011 + denovo2029 + denovo2120 + denovo2147 + denovo2179 + denovo2189 + denovo2245 + denovo2269 +
                                 denovo2290 + denovo2297 + denovo2351 + denovo2359 + denovo2367 + denovo2418 + denovo2465 + denovo2479 +
                                 denovo2491 + denovo2505 + denovo2549 + denovo2554 + denovo2556 + denovo2567 + denovo2581 + denovo2606 +
                                 denovo2619 + denovo2677 + denovo2713 + denovo2800 + denovo2808 + denovo2824 + denovo2835 + denovo2839 +
                                 denovo2848 + denovo2862 + denovo2867 + denovo2952 + denovo3030 + denovo3065 + denovo3079 + denovo3151 +
                                 denovo3198 + denovo3246 + denovo3319 + denovo3335 + denovo3346 + denovo3363 + denovo3368 + denovo3373 +
                                 denovo3406 + denovo3426 + denovo3432 + denovo3449 + denovo3461 + denovo3480 + denovo3489 + denovo3536 +
                                 denovo3599 + denovo3643 + denovo3652 + denovo3654 + denovo3673 + denovo3679 + denovo3700 + denovo3745 +
                                 denovo3768 + denovo3777 + denovo3838 + denovo3853 + denovo3858 + denovo3882 + denovo3929 + denovo3959 +
                                 denovo4025 + denovo4040 + denovo4074 + denovo4079 + denovo4082 + denovo4106 + denovo4114 + denovo4154 +
                                 denovo4164 + denovo4205 + denovo4228 + denovo4236 + denovo4247 + denovo4283 + denovo4289 + denovo4312 +
                                 denovo4319 + denovo4325 + denovo4340 + denovo4352 + denovo4378 + denovo4408 + denovo4419 + denovo4424 +
                                 denovo4437 + denovo4440 + denovo4499 + denovo4501 + denovo4533 + denovo4575 + denovo4609 + denovo4631 +
                                 denovo4636 + denovo4711 + denovo4746 + denovo4759 + denovo4762 + denovo4770 + denovo4805 + denovo4834 +
                                 denovo4836 + denovo4859 + denovo4875 + denovo4900 + denovo4904 + denovo4984 + denovo4988 + denovo5067 +
                                 denovo5123 + denovo5133 + denovo5168 + denovo5198 + denovo5204 + denovo5226 + denovo5247 + denovo5263 +
                                 denovo5269 + denovo5288 + denovo5293 + denovo5324 + denovo5414 + denovo5448 + denovo5509 + denovo5550 +
                                 denovo5581 + denovo5643 + denovo5647 + denovo5690 + denovo5741 + denovo5776 + denovo5801 + denovo5859 +
                                 denovo5971 + denovo5975 + denovo5980 + denovo5986 + denovo6019 + denovo6122 + denovo6138 + denovo6148 +
                                 denovo6150 + denovo6209 + denovo6214 + denovo6226 + denovo6255 + denovo6316 + denovo6402 + denovo6411 +
                                 denovo6438 + denovo6439 + denovo6451 + denovo6453 + denovo6472 + denovo6473 + denovo6493 + denovo6496 +
                                 denovo6499 + denovo6504 + denovo6520 + denovo6521 + denovo6527 + denovo6550 + denovo6573 + denovo6577 +
                                 denovo6583 + denovo6608 + denovo6660 + denovo6689 + denovo6709 + denovo6735 + denovo6749 + denovo6811 +
                                 denovo6816 + denovo6833 + denovo6835 + denovo6856 + denovo6872 + denovo6878 + denovo6909 + denovo6918 +
                                 denovo7003 + denovo7005 + denovo7007 + denovo7053 + denovo7060 + denovo7074 + denovo7086 + denovo7124 +
                                 denovo7136 + denovo7141 + denovo7162 + denovo7171 + denovo7186 + denovo7187 + denovo7195 + denovo7197 +
                                 denovo7204 + denovo7206 + denovo7238 + denovo7239 + denovo7251 + denovo7257 + denovo7288 + denovo7303 +
                                 denovo7334 + denovo7336 + denovo7342 + denovo7347 + denovo7359 + denovo7361 + denovo7369 + denovo7414 +
                                 denovo7417 + denovo7432 + denovo7456 + denovo7459 + denovo7465 + denovo7470 + denovo7473 + denovo7495 +
                                 denovo7508 + denovo7522 + denovo7528 + denovo7541 + denovo7551 + denovo7559 + denovo7672 + denovo7682 +
                                 denovo7729 + denovo7741 + denovo7753 + denovo7771 + denovo7794 + denovo7800 + denovo7841 + denovo7866 +
                                 denovo7880 + denovo7895 + denovo7908 + denovo7939 + denovo7953 + denovo7958 + denovo7966 + denovo7995 +
                                 denovo8002 + denovo8017 + denovo8034 + denovo8057 + denovo8087 + denovo8121 + denovo8126 + denovo8155 +
                                 denovo8181 + denovo8186 + denovo8188 + denovo8230 + denovo8263 + denovo8288 + denovo8303 + denovo8350 +
                                 denovo8352 + denovo8357 + denovo8385 + denovo8386 + denovo8405 + denovo8428 + denovo8463 + denovo8485 +
                                 denovo8494 + denovo8502 + denovo8563 + denovo8582 + denovo8599 + denovo8602 + denovo8606 + denovo8696 +
                                 denovo8704 + denovo8714 + denovo8735 + denovo8741 + denovo8790 + denovo8817 + denovo8861 + denovo8879 +
                                 denovo8889 + denovo8906 + denovo8908 + denovo8929 + denovo8970 + denovo8978 + denovo8986 + denovo8991 +
                                 denovo9005 + denovo9017 + denovo9030 + denovo9075 + denovo9097 + denovo9115 + denovo9118 + denovo9181 +
                                 denovo9204 + denovo9206 + denovo9216 + denovo9227 + denovo9248 + denovo9265 + denovo9271 + denovo9276 +
                                 denovo9341 + denovo9350 + denovo9366 + denovo9409 + denovo9431 + denovo9457 + denovo9464 + denovo9476 +
                                 denovo9501 + denovo9502 + denovo9606 + denovo9615 + denovo9636 + denovo9646 + denovo9678 + denovo9684 +
                                 denovo9687 + denovo9711 + denovo9726 + denovo9758 + denovo9761 + denovo9771 + denovo9802 + denovo9809 +
                                 denovo9812 + denovo9894 + denovo9900 + denovo9914 + denovo9920 + denovo9932 + denovo9953 + denovo9957 +
                                 denovo9967 + denovo10061 + denovo10109 + denovo10127 + denovo10155 + denovo10157 + denovo10176 + denovo10181 +
                                 denovo10213 + denovo10226 + denovo10227 + denovo10241 + denovo10264 + denovo10282 + denovo10315 + denovo10329 +
                                 denovo10333 + denovo10352 + denovo10366 + denovo10376 + denovo10383 + denovo10395 + denovo10447 + denovo10452 +
                                 denovo10464 + denovo10487 + denovo10516 + denovo10559 + denovo10580 + denovo10581 + denovo10602 + denovo10618 +
                                 denovo10621 + denovo10634 + denovo10639 + denovo10646 + denovo10647 + denovo10648 + denovo10650 + denovo10675 +
                                 denovo10689 + denovo10718 + denovo10732 + denovo10752 + denovo10831 + denovo10842 + denovo10866 + denovo10889 +
                                 denovo10898 + denovo10992 + denovo11010 + denovo11038 + denovo11050 + denovo11051 + denovo11073 + denovo11095 +
                                 denovo11152 + denovo11191 + denovo11209 + denovo11229 + denovo11234 + denovo11235 + denovo11245 + denovo11247 +
                                 denovo11249 + denovo11286 + denovo11290 + denovo11329 + denovo11367 + denovo11392 + denovo11403 + denovo11410 +
                                 denovo11419 + denovo11430 + denovo11505 + denovo11556 + denovo11574 + denovo11588 + denovo11597 + denovo11623 +
                                 denovo11629 + denovo11646 + denovo11689 + denovo11730 + denovo11772 + denovo11778 + denovo11784 + denovo11789 +
                                 denovo11814 + denovo11817 + denovo11820 + denovo11844 + denovo11866 + denovo11873 + denovo11890 + denovo11958 +
                                 denovo11971 + denovo11980 + denovo12022 + denovo12027 + denovo12030 + denovo12052 + denovo12062 + denovo12087 +
                                 denovo12137 + denovo12185 + denovo12192 + denovo12221 + denovo12242 + denovo12248 + denovo12255 + denovo12265 +
                                 denovo12274 + denovo12275 + denovo12407 + denovo12420 + denovo12474 + denovo12497 + denovo12548 + denovo12609 +
                                 denovo12621 + denovo12641 + denovo12649 + denovo12658 + denovo12690 + denovo12693 + denovo12728 + denovo12736 +
                                 denovo12765 + denovo12777 + denovo12783 + denovo12808 + denovo12842 + denovo12855 + denovo12870 + denovo12873 +
                                 denovo12881 + denovo12890 + denovo12913 + denovo12943 + denovo12952 + denovo12960 + denovo12963 + denovo13006 +
                                 denovo13017 + denovo13068 + denovo13075 + denovo13081 + denovo13096 + denovo13120 + denovo13121 + denovo13148 +
                                 denovo13183 + denovo13197 + denovo13265 + denovo13269 + denovo13312 + denovo13334 + denovo13336 + denovo13348 +
                                 denovo13417 + denovo13466 + denovo13515 + denovo13517 + denovo13520 + denovo13524 + denovo13529 + denovo13536 +
                                 denovo13547 + denovo13569 + denovo13641 + denovo13669 + denovo13673 + denovo13708 + denovo13710 + denovo13767 +
                                 denovo13781 + denovo13800 + denovo13838 + denovo13877 + denovo13900 + denovo13928 + denovo13935 + denovo13948 +
                                 denovo13952 + denovo13996 + denovo14012 + denovo14019 + denovo14029 + denovo14044 + denovo14062 + denovo14089 +
                                 denovo14125 + denovo14153 + denovo14202 + denovo14214 + denovo14249 + denovo14250 + denovo14251 + denovo14255 +
                                 denovo14256 + denovo14273 + denovo14278 + denovo14294 + denovo14308 + denovo14329 + denovo14332 + denovo14346 +
                                 denovo14358 + denovo14400 + denovo14411 + denovo14544 + denovo14565 + denovo14640 + denovo14642 + denovo14666 +
                                 denovo14679 + denovo14703 + denovo14705 + denovo14713 + denovo14781 + denovo14802 + denovo14803 + denovo14809 +
                                 denovo14818 + denovo14820 + denovo14821 + denovo14960 + denovo14971 + denovo14978 + denovo14982 + denovo14991 +
                                 denovo15000 + denovo15015 + denovo15066 + denovo15076 + denovo15101 + denovo15120 + denovo15150 + denovo15151 +
                                 denovo15185 + denovo15212 + denovo15284 + denovo15307 + denovo15356 + denovo15359 + denovo15368 + denovo15371 +
                                 denovo15376 + denovo15430 + denovo15438 + denovo15449 + denovo15450 + denovo15509 + denovo15558 + denovo15581 +
                                 denovo15601 + denovo15711 + denovo15713 + denovo15714 + denovo15729 + denovo15764 + denovo15801 + denovo15811 +
                                 denovo15835 + denovo15849 + denovo15855 + denovo15895 + denovo15899 + denovo15968 + denovo15998 + denovo16021 +
                                 denovo16058 + denovo16064 + denovo16075 + denovo16081 + denovo16177 + denovo16182 + denovo16194 + denovo16353 +
                                 denovo16407 + denovo16437 + denovo16466 + denovo16470 + denovo16504 + denovo16516 + denovo16524 + denovo16528 +
                                 denovo16552 + denovo16577 + denovo16584 + denovo16603 + denovo16678 + denovo16722 + denovo16737 + denovo16746 +
                                 denovo16763 + denovo16782 + denovo16812 + denovo16819 + denovo16828 + denovo16857 + denovo16859 + denovo16861 +
                                 denovo16914 + denovo16915 + denovo16918 + denovo16923 + denovo16979 + denovo17022 + denovo17052 + denovo17078 +
                                 denovo17092 + denovo17190 + denovo17194 + denovo17232 + denovo17277 + denovo17292 + denovo17303 + denovo17304 +
                                 denovo17324 + denovo17336 + denovo17365 + denovo17373 + denovo17402 + denovo17436 + denovo17515 + denovo17562 +
                                 denovo17628 + denovo17633 + denovo17634 + denovo17645 + denovo17646 + denovo17683 + denovo17691 + denovo17713 +
                                 denovo17735 + denovo17762 + denovo17769 + denovo17836 + denovo17846 + denovo17856 + denovo17859 + denovo17874 +
                                 denovo17891 + denovo17940 + denovo17955 + denovo17958 + denovo17967 + denovo17969 + denovo17990 + denovo17991 +
                                 denovo18014 + denovo18038 + denovo18090 + denovo18093 + denovo18099 + denovo18105 + denovo18163 + denovo18176 +
                                 denovo18264 + denovo18296 + denovo18301 + denovo18302 + denovo18355 + denovo18368 + denovo18375 + denovo18416 +
                                 denovo18418 + denovo18452 + denovo18488 + denovo18498 + denovo18534 + denovo18539 + denovo18588 + denovo18589 +
                                 denovo18633 + denovo18674 + denovo18680 + denovo18716 + denovo18739 + denovo18760 + denovo18821 + denovo18848 +
                                 denovo18874 + denovo18906 + denovo18941 + denovo18975 + denovo18980 + denovo18994 + denovo19077 + denovo19096 +
                                 denovo19124 + denovo19146 + denovo19218 + denovo19265 + denovo19290 + denovo19392 + denovo19400 + denovo19416 +
                                 denovo19450 + denovo19499 + denovo19517 + denovo19560 + denovo19576 + denovo19619 + denovo19650 + denovo19664 +
                                 denovo19667 + denovo19715 + denovo19732 + denovo19773 + denovo19775 + denovo19829 + denovo19848 + denovo19892 +
                                 denovo19909 + denovo20052 + denovo20065 + denovo20075 + denovo20111 + denovo20165 + denovo20196 + denovo20215 +
                                 denovo20226 + denovo20234 + denovo20275 + denovo20282 + denovo20290 + denovo20345 + denovo20361 + denovo20371 +
                                 denovo20372 + denovo20374 + denovo20403 + denovo20451 + denovo20462 + denovo20466 + denovo20491 + denovo20522 +
                                 denovo20530 + denovo20586 + denovo20606 + denovo20682 + denovo20687 + denovo20693 + denovo20701 + denovo20741 +
                                 denovo20808 + denovo20813 + denovo20904 + denovo20940 + denovo20997 + denovo21001 + denovo21010 + denovo21049 +
                                 denovo21064 + denovo21071 + denovo21075 + denovo21085 + denovo21094 + denovo21096 + denovo21098 + denovo21150 +
                                 denovo21155 + denovo21169 + denovo21173 + denovo21178 + denovo21192 + denovo21200 + denovo21348 + denovo21353 +
                                 denovo21378 + denovo21389 + denovo21413 + denovo21425 + denovo21471 + denovo21491 + denovo21499 + denovo21525 +
                                 denovo21543 + denovo21568 + denovo21580 + denovo21587 + denovo21627 + denovo21634 + denovo21638 + denovo21642 +
                                 denovo21659 + denovo21700 + denovo21726 + denovo21739 + denovo21762 + denovo21767 + denovo21793 + denovo21843 +
                                 denovo21858 + denovo21874 + denovo21891 + denovo21917 + denovo21924 + denovo21953 + denovo21989 + denovo22004 +
                                 denovo22022 + denovo22029 + denovo22031 + denovo22065 + denovo22070 + denovo22110 + denovo22137 + denovo22139 +
                                 denovo22171 + denovo22205 + denovo22214 + denovo22233 + denovo22287 + denovo22294 + denovo22309 + denovo22314 +
                                 denovo22330 + denovo22338 + denovo22344 + denovo22365 + denovo22423 + denovo22503 + denovo22539 + denovo22553 +
                                 denovo22561 + denovo22596 + denovo22629 + denovo22695 + denovo22721 + denovo22821 + denovo22850 + denovo22879 +
                                 denovo22881 + denovo22923 + denovo22934 + denovo22941 + denovo22972 + denovo22974 + denovo22975 + denovo23006 +
                                 denovo23064 + denovo23106 + denovo23118 + denovo23128 + denovo23143 + denovo23155 + denovo23161 + denovo23173 +
                                 denovo23184 + denovo23191 + denovo23192 + denovo23245 + denovo23266 + denovo23274 + denovo23280 + denovo23338 +
                                 denovo23345 + denovo23374 + denovo23377 + denovo23391 + denovo23424 + denovo23448 + denovo23453 + denovo23509 +
                                 denovo23544 + denovo23552 + denovo23560 + denovo23562 + denovo23591 + denovo23625 + denovo23691 + denovo23702 +
                                 denovo23712 + denovo23734 + denovo23741 + denovo23769 + denovo23798 + denovo23816 + denovo23833 + denovo23893 +
                                 denovo23912 + denovo23954 + denovo23971 + denovo23986 + denovo23992 + denovo23999 + denovo24032 + denovo24045 +
                                 denovo24064 + denovo24080 + denovo24082 + denovo24117 + denovo24189 + denovo24262 + denovo24275 + denovo24338 +
                                 denovo24390 + denovo24397 + denovo24433 + denovo24439 + denovo24533 + denovo24547 + denovo24552 + denovo24608 +
                                 denovo24641 + denovo24647 + denovo24660 + denovo24668 + denovo24742 + denovo24751 + denovo24757 + denovo24758 +
                                 denovo24782 + denovo24869 + denovo24872 + denovo24881 + denovo24900 + denovo24939 + denovo24958 + denovo25001 +
                                 denovo25027 + denovo25105 + denovo25137 + denovo25154 + denovo25205 + denovo25228 + denovo25229 + denovo25247 +
                                 denovo25261 + denovo25280 + denovo25284 + denovo25306 + denovo25347 + denovo25367 + denovo25388 + denovo25393 +
                                 denovo25427 + denovo25429 + denovo25464 + denovo25493 + denovo25582 + denovo25590 + denovo25594 + denovo25662 +
                                 denovo25699 + denovo25700 + denovo25726 + denovo25761 + denovo25797 + denovo25890 + denovo25903 + denovo25927 +
                                 denovo25965 + denovo25990 + denovo25995 + denovo26040 + denovo26076 + denovo26099 + denovo26102 + denovo26135 +
                                 denovo26184 + denovo26220 + denovo26268 + denovo26278 + denovo26310 + denovo26361 + denovo26424 + denovo26453 +
                                 denovo26471 + denovo26474 + denovo26511 + denovo26541 + denovo26561 + denovo26585 + denovo26591 + denovo26602 +
                                 denovo26610 + denovo26633 + denovo26683 + denovo26719 + denovo26752 + denovo26785 + denovo26893 + denovo26902 +
                                 denovo26927 + denovo26954 + denovo26979 + denovo27015 + denovo27024 + denovo27039 + denovo27045 + denovo27054 +
                                 denovo27066 + denovo27126 + denovo27139 + denovo27170 + denovo27218 + denovo27221 + denovo27254 + denovo27289 +
                                 denovo27354 + denovo27355 + denovo27412 + denovo27461 + denovo27504 + denovo27534 + denovo27535 + denovo27603 +
                                 denovo27648 + denovo27673 + denovo27682 + denovo27700 + denovo27702 + denovo27725 + denovo27869 + denovo27874 +
                                 denovo27878 + denovo27900 + denovo27928 + denovo27946 + denovo27968 + denovo28021 + denovo28113 + denovo28128 +
                                 denovo28131 + denovo28253 + denovo28266 + denovo28333 + denovo28387 + denovo28392 + denovo28451 + denovo28485 +
                                 denovo28517 + denovo28538 + denovo28595 + denovo28605 + denovo28613 + denovo28614 + denovo28625 + denovo28633 +
                                 denovo28699 + denovo28732 + denovo28750 + denovo28787 + denovo28790 + denovo28797 + denovo28799 + denovo28838 +
                                 denovo28848 + denovo28872 + denovo28885 + denovo28979 + denovo29009 + denovo29016 + denovo29038 + denovo29062 +
                                 denovo29083 + denovo29163 + denovo29188 + denovo29270 + denovo29294 + denovo29310 + denovo29327 + denovo29355 +
                                 denovo29362 + denovo29381 + denovo29384 + denovo29388 + denovo29399 + denovo29420 + denovo29437 + denovo29438 +
                                 denovo29442 + denovo29459 + denovo29567 + denovo29620 + denovo29639 + denovo29642 + denovo29724 + denovo29754 +
                                 denovo29762 + denovo29793 + denovo29812 + denovo29821 + denovo29841 + denovo29842 + denovo29844 + denovo29920 +
                                 denovo29976 + denovo30001 + denovo30013 + denovo30050 + denovo30054 + denovo30058 + denovo30061 + denovo30083 +
                                 denovo30112 + denovo30146 + denovo30157 + denovo30167 + denovo30186 + denovo30218 + denovo30224 + denovo30283 +
                                 denovo30317 + denovo30366 + denovo30371 + denovo30383 + denovo30393 + denovo30420 + denovo30460 + denovo30484 +
                                 denovo30573 + denovo30580 + denovo30592 + denovo30635 + denovo30693 + denovo30715 + denovo30728 + denovo30734 +
                                 denovo30756 + denovo30777 + denovo30782 + denovo30811 + denovo30875 + denovo30883 + denovo30919 + denovo30940 +
                                 denovo30941 + denovo31003 + denovo31036 + denovo31052 + denovo31070 + denovo31132 + denovo31177 + denovo31274 +
                                 denovo31322 + denovo31348 + denovo31376 + denovo31385 + denovo31403 + denovo31414 + denovo31420 + denovo31430 +
                                 denovo31449 + denovo31460 + denovo31476 + denovo31477 + denovo31486 + denovo31487 + denovo31506 + denovo31534 +
                                 denovo31546 + denovo31557 + denovo31574 + denovo31576 + denovo31669 + denovo31673 + denovo31688 + denovo31746 +
                                 denovo31773 + denovo31820 + denovo31843 + denovo31846 + denovo31904 + denovo31934 + denovo31944 + denovo31959 +
                                 denovo31968 + denovo31980 + denovo31993 + denovo31999 + denovo32001 + denovo32013 + denovo32049 + denovo32058 +
                                 denovo32062 + denovo32079 + denovo32084 + denovo32094 + denovo32099 + denovo32199 + denovo32203 + denovo32266 +
                                 denovo32269 + denovo32278 + denovo32304 + denovo32477 + denovo32498 + denovo32547 + denovo32564 + denovo32569 +
                                 denovo32601 + denovo32656 + denovo32679 + denovo32738 + denovo32824 + denovo32855 + denovo32901 + denovo32956 +
                                 denovo32960 + denovo33007 + denovo33021 + denovo33022 + denovo33074 + denovo33136 + denovo33171 + denovo33245 +
                                 denovo33246 + denovo33276 + denovo33316 + denovo33340 + denovo33347 + denovo33360 + denovo33376 + denovo33440 +
                                 denovo33472 + denovo33491 + denovo33495 + denovo33512 + denovo33518 + denovo33524 + denovo33571 + denovo33607 +
                                 denovo33616 + denovo33632 + denovo33644 + denovo33702 + denovo33729 + denovo33733 + denovo33757 + denovo33778 +
                                 denovo33793 + denovo33800 + denovo33831 + denovo33834 + denovo33844 + denovo33877 + denovo33891 + denovo33925 +
                                 denovo33941 + denovo34003 + denovo34012 + denovo34019 + denovo34040 + denovo34041 + denovo34049 + denovo34061 +
                                 denovo34111 + denovo34149 + denovo34244 + denovo34260 + denovo34337 + denovo34351 + denovo34364 + denovo34369 +
                                 denovo34407 + denovo34455 + denovo34471 + denovo34492 + denovo34499 + denovo34513 + denovo34525 + denovo34537 +
                                 denovo34585 + denovo34616 + denovo34636 + denovo34680 + denovo34713 + denovo34718 + denovo34725 + denovo34751 +
                                 denovo34763 + denovo34776 + denovo34808 + denovo34816 + denovo34840 + denovo34852 + denovo34873 + denovo34881 +
                                 denovo34918 + denovo34937 + denovo34998 + denovo35012 + denovo35022 + denovo35036 + denovo35062 + denovo35083 +
                                 denovo35089 + denovo35121 + denovo35124 + denovo35129 + denovo35133 + denovo35140 + denovo35153 + denovo35174 +
                                 denovo35196 + denovo35207 + denovo35216 + denovo35227 + denovo35246 + denovo35348 + denovo35369 + denovo35381 +
                                 denovo35415 + denovo35426 + denovo35451 + denovo35470 + denovo35487 + denovo35540 + denovo35585 + denovo35609 +
                                 denovo35655 + denovo35689 + denovo35729 + denovo35734 + denovo35756 + denovo35763 + denovo35790 + denovo35861 +
                                 denovo35862 + denovo35899 + denovo35910 + denovo35917 + denovo35969 + denovo36019 + denovo36051 + denovo36059 +
                                 denovo36156 + denovo36166 + denovo36224 + denovo36231 + denovo36273 + denovo36279 + denovo36295 + denovo36392 +
                                 denovo36397 + denovo36401 + denovo36409 + denovo36459 + denovo36484 + denovo36503 + denovo36505 + denovo36519 +
                                 denovo36522 + denovo36539 + denovo36541 + denovo36545 + denovo36556 + denovo36565 + denovo36602 + denovo36603 +
                                 denovo36617 + denovo36639 + denovo36661 + denovo36671 + denovo36676 + denovo36680 + denovo36684 + denovo36725 +
                                 denovo36788 + denovo36805 + denovo36855 + denovo36892 + denovo36939 + denovo37007 + denovo37009 + denovo37046 +
                                 denovo37143 + denovo37144 + denovo37176 + denovo37230 + denovo37232 + denovo37281 + denovo37289 + denovo37328 +
                                 denovo37335 + denovo37352 + denovo37354 + denovo37357 + denovo37361 + denovo37374 + denovo37394 + denovo37483 +
                                 denovo37554 + denovo37581 + denovo37584 + denovo37596 + denovo37605 + denovo37637 + denovo37656 + denovo37659 +
                                 denovo37680 + denovo37683 + denovo37685 + denovo37692 + denovo37717 + denovo37740 + denovo37764 + denovo37786 +
                                 denovo37799 + denovo37800 + denovo37811 + denovo37819 + denovo37831 + denovo37850 + denovo37896 + denovo37933 +
                                 denovo37978 + denovo37983 + denovo37996 + denovo38060 + denovo38175 + denovo38177 + denovo38192 + denovo38195 +
                                 denovo38198 + denovo38251 + denovo38259 + denovo38279 + denovo38297 + denovo38301 + denovo38303 + denovo38319 +
                                 denovo38347 + denovo38364 + denovo38434 + denovo38460 + denovo38479 + denovo38500 + denovo38517 + denovo38523 +
                                 denovo38524 + denovo38561 + denovo38576 + denovo38634 + denovo38640 + denovo38656 + denovo38692 + denovo38695 +
                                 denovo38726 + denovo38732 + denovo38737 + denovo38762 + denovo38815 + denovo38826 + denovo38868 + denovo38905 +
                                 denovo39013 + denovo39032 + denovo39039 + denovo39049 + denovo39090 + denovo39140 + denovo39143 + denovo39151 +
                                 denovo39164 + denovo39182 + denovo39216 + denovo39243 + denovo39281 + denovo39304 + denovo39307 + denovo39344 +
                                 denovo39432 + denovo39435 + denovo39464 + denovo39476 + denovo39485 + denovo39507 + denovo39569 + denovo39601 +
                                 denovo39623 + denovo39649 + denovo39674 + denovo39720 + denovo39729 + denovo39798 + denovo39817 + denovo39942 +
                                 denovo39955 + denovo39996 + denovo40008 + denovo40057 + denovo40065 + denovo40067 + denovo40071 + denovo40078 +
                                 denovo40121 + denovo40125 + denovo40136 + denovo40138 + denovo40143 + denovo40154 + denovo40171 + denovo40174 +
                                 denovo40212 + denovo40246 + denovo40248 + denovo40257 + denovo40290 + denovo40345 + denovo40350 + denovo40395 +
                                 denovo40403 + denovo40437 + denovo40479 + denovo40480 + denovo40510 + denovo40526 + denovo40569 + denovo40570 +
                                 denovo40608 + denovo40637 + denovo40661 + denovo40728 + denovo40744 + denovo40755 + denovo40774 + denovo40780 +
                                 denovo40843 + denovo40910 + denovo40932 + denovo40971 + denovo41003 + denovo41023 + denovo41027 + denovo41029 +
                                 denovo41067 + denovo41108 + denovo41114 + denovo41164 + denovo41166 + denovo41169 + denovo41172 + denovo41182 +
                                 denovo41235 + denovo41238 + denovo41240 + denovo41287 + denovo41299 + denovo41306 + denovo41325 + denovo41363 +
                                 denovo41381 + denovo41393 + denovo41425 + denovo41433 + denovo41499 + denovo41548 + denovo41553 + denovo41626 +
                                 denovo41637 + denovo41671 + denovo41698 + denovo41707 + denovo41719 + denovo41738 + denovo41751 + denovo41806 +
                                 denovo41845 + denovo41868 + denovo41888 + denovo41911 + denovo42002 + denovo42005 + denovo42058 + denovo42157 +
                                 denovo42159 + denovo42193 + denovo42253 + denovo42316 + denovo42327 + denovo42340 + denovo42364 + denovo42381 +
                                 denovo42405 + denovo42419 + denovo42436 + denovo42444 + denovo42456 + denovo42485 + denovo42521 + denovo42538 +
                                 denovo42592 + denovo42632 + denovo42711 + denovo42716 + denovo42905 + denovo42931 + denovo42936 + denovo42944 +
                                 denovo42957 + denovo42976 + denovo42979 + denovo42993 + denovo42997 + denovo43004 + denovo43005 + denovo43018 +
                                 denovo43065 + denovo43106 + denovo43124 + denovo43133 + denovo43190 + denovo43231 + denovo43232 + denovo43258 +
                                 denovo43259 + denovo43280 + denovo43304 + denovo43350 + denovo43357 + denovo43373 + denovo43425 + denovo43475 +
                                 denovo43561 + denovo43593 + denovo43601 + denovo43671 + denovo43711 + denovo43719 + denovo43806 + denovo43812 +
                                 denovo43854 + denovo43857 + denovo43970 + denovo43977 + denovo43984 + denovo43997 + denovo44016 + denovo44034 +
                                 denovo44038 + denovo44058 + denovo44077 + denovo44079 + denovo44093 + denovo44095 + denovo44146 + denovo44163 +
                                 denovo44234 + denovo44260 + denovo44276 + denovo44289 + denovo44301 + denovo44340 + denovo44358 + denovo44414 +
                                 denovo44496 + denovo44499 + denovo44502 + denovo44510 + denovo44512 + denovo44527 + denovo44528 + denovo44532 +
                                 denovo44538 + denovo44545 + denovo44609 + denovo44662 + denovo44669 + denovo44678 + denovo44696 + denovo44702 +
                                 denovo44708 + denovo44715 + denovo44726 + denovo44745 + denovo44793 + denovo44816 + denovo44863 + denovo44867 +
                                 denovo44912 + denovo44940 + denovo44979 + denovo44994 + denovo45004 + denovo45013 + denovo45018 + denovo45024 +
                                 denovo45046 + denovo45063 + denovo45097 + denovo45108 + denovo45110 + denovo45140 + denovo45153 + denovo45175 +
                                 denovo45194 + denovo45195 + denovo45202 + denovo45230 + denovo45242 + denovo45298 + denovo45302 + denovo45316 +
                                 denovo45322 + denovo45419 + denovo45422 + denovo45437 + denovo45453 + denovo45463 + denovo45476 + denovo45516 +
                                 denovo45590 + denovo45637 + denovo45642 + denovo45680 + denovo45690 + denovo45701 + denovo45702 + denovo45707 +
                                 denovo45708 + denovo45770 + denovo45783 + denovo45876 + denovo45885 + denovo45898 + denovo45923 + denovo45928 +
                                 denovo45939 + denovo45967 + denovo45971 + denovo46052 + denovo46093 + denovo46097 + denovo46109 + denovo46141 +
                                 denovo46165 + denovo46188 + denovo46203 + denovo46226 + denovo46235 + denovo46263 + denovo46268 + denovo46293 +
                                 denovo46322 + denovo46336 + denovo46344 + denovo46346 + denovo46373 + denovo46466 + denovo46472 + denovo46516 +
                                 denovo46542 + denovo46555 + denovo46579 + denovo46598 + denovo46654 + denovo46655 + denovo46656 + denovo46687 +
                                 denovo46702 + denovo46721 + denovo46816 + denovo46836 + denovo46870 + denovo46874 + denovo46879 + denovo46883 +
                                 denovo46906 + denovo46926 + denovo46934 + denovo46966 + denovo46984 + denovo47032 + denovo47041 + denovo47049 +
                                 denovo47052 + denovo47066 + denovo47109 + denovo47123 + denovo47259 + denovo47312 + denovo47382 + denovo47441 +
                                 denovo47510 + denovo47528 + denovo47557 + denovo47587 + denovo47595 + denovo47627 + denovo47661 + denovo47667 +
                                 denovo47682 + denovo47717 + denovo47731 + denovo47801 + denovo47803 + denovo47808 + denovo47816 + denovo47824 +
                                 denovo47826 + denovo47852 + denovo47880 + denovo47891 + denovo47919 + denovo47944 + denovo48005 + denovo48013 +
                                 denovo48031 + denovo48087 + denovo48101 + denovo48163 + denovo48177 + denovo48209 + denovo48210 + denovo48288 +
                                 denovo48298 + denovo48323 + denovo48353 + denovo48384 + denovo48434 + denovo48442 + denovo48449 + denovo48458 +
                                 denovo48488 + denovo48499 + denovo48520 + denovo48531 + denovo48543 + denovo48546 + denovo48554 + denovo48579 +
                                 denovo48582 + denovo48616 + denovo48660 + denovo48676 + denovo48708 + denovo48726 + denovo48747 + denovo48764 +
                                 denovo48766 + denovo48816 + denovo48846 + denovo48861 + denovo48863 + denovo48910 + denovo48915 + denovo48972 +
                                 denovo49012 + denovo49063 + denovo49070 + denovo49084 + denovo49086 + denovo49107 + denovo49111 + denovo49122 +
                                 denovo49127 + denovo49157 + denovo49176 + denovo49270 + denovo49280 + denovo49298 + denovo49364 + denovo49422 +
                                 denovo49447 + denovo49460 + denovo49518 + denovo49539 + denovo49541 + denovo49572 + denovo49594 + denovo49602 +
                                 denovo49628 + denovo49679 + denovo49698 + denovo49735 + denovo49768 + denovo49782 + denovo49789 + denovo49803 +
                                 denovo49847 + denovo49872 + denovo49914 + denovo49942 + denovo49975 + denovo50005 + denovo50043 + denovo50053 +
                                 denovo50125 + denovo50128 + denovo50166 + denovo50185 + denovo50189 + denovo50210 + denovo50218 + denovo50222 +
                                 denovo50229 + denovo50243 + denovo50263 + denovo50284 + denovo50333 + denovo50368 + denovo50400 + denovo50434 +
                                 denovo50458 + denovo50542 + denovo50561 + denovo50626 + denovo50628 + denovo50698 + denovo50722 + denovo50725 +
                                 denovo50727 + denovo50740 + denovo50741 + denovo50771 + denovo50802 + denovo50815 + denovo50818 + denovo50832 +
                                 denovo50880 + denovo50893 + denovo50894 + denovo50921 + denovo50929 + denovo50951 + denovo50974 + denovo50985 +
                                 denovo51035 + denovo51045 + denovo51088 + denovo51117 + denovo51123 + denovo51153 + denovo51156 + denovo51159 +
                                 denovo51197 + denovo51213 + denovo51244 + denovo51273 + denovo51288 + denovo51323 + denovo51376 + denovo51405 +
                                 denovo51433 + denovo51477 + denovo51537 + denovo51540 + denovo51549 + denovo51562 + denovo51565 + denovo51599 +
                                 denovo51645 + denovo51650 + denovo51701 + denovo51721 + denovo51745 + denovo51756 + denovo51785 + denovo51787 +
                                 denovo51796 + denovo51826 + denovo51830 + denovo51843 + denovo51853 + denovo51878 + denovo51921 + denovo52036 +
                                 denovo52074 + denovo52146 + denovo52155 + denovo52226 + denovo52251 + denovo52256 + denovo52282 + denovo52295 +
                                 denovo52337 + denovo52345 + denovo52348 + denovo52363 + denovo52377 + denovo52433 + denovo52488 + denovo52528 +
                                 denovo52538 + denovo52542 + denovo52552 + denovo52562 + denovo52601 + denovo52605 + denovo52630 + denovo52631 +
                                 denovo52758 + denovo52780 + denovo52784 + denovo52836 + denovo52847 + denovo52851 + denovo52880 + denovo52891 +
                                 denovo52902 + denovo52934 + denovo52935 + denovo52967 + denovo52968 + denovo53013 + denovo53024 + denovo53064 +
                                 denovo53073 + denovo53103 + denovo53110 + denovo53122 + denovo53139 + denovo53149 + denovo53216 + denovo53236 +
                                 denovo53250 + denovo53269 + denovo53278 + denovo53320 + denovo53327 + denovo53332 + denovo53346 + denovo53356 +
                                 denovo53360 + denovo53364 + denovo53401 + denovo53439 + denovo53451 + denovo53468 + denovo53504 + denovo53550 +
                                 denovo53564 + denovo53626 + denovo53643 + denovo53665 + denovo53671 + denovo53772 + denovo53821 + denovo53826 +
                                 denovo53835 + denovo53852 + denovo53866 + denovo53913 + denovo53965 + denovo53967 + denovo54032 + denovo54068 +
                                 denovo54077 + denovo54159 + denovo54176 + denovo54181 + denovo54234 + denovo54251 + denovo54252 + denovo54266 +
                                 denovo54269 + denovo54296 + denovo54311 + denovo54321 + denovo54361 + denovo54433 + denovo54434 + denovo54440 +
                                 denovo54484 + denovo54527 + denovo54532 + denovo54582 + denovo54598 + denovo54605 + denovo54617 + denovo54635 +
                                 denovo54651 + denovo54652 + denovo54670 + denovo54697 + denovo54727 + denovo54792 + denovo54803 + denovo54806 +
                                 denovo54875 + denovo54883 + denovo54886 + denovo54900 + denovo54911 + denovo54917 + denovo54932 + denovo54998 +
                                 denovo55042 + denovo55134 + denovo55155 + denovo55184 + denovo55189 + denovo55202 + denovo55260 + denovo55267 +
                                 denovo55269 + denovo55274 + denovo55316 + denovo55328 + denovo55341 + denovo55355 + denovo55371 + denovo55403 +
                                 denovo55425 + denovo55429 + denovo55498 + denovo55505 + denovo55522 + denovo55583 + denovo55599 + denovo55601 +
                                 denovo55607 + denovo55609 + denovo55767 + denovo55832 + denovo55845 + denovo55853 + denovo55891 + denovo55896 +
                                 denovo55930 + denovo55945 + denovo55952 + denovo55967 + denovo55995 + denovo56011 + denovo56050 + denovo56056 +
                                 denovo56087 + denovo56120 + denovo56147 + denovo56155 + denovo56156 + denovo56160 + denovo56174 + denovo56178 +
                                 denovo56215 + denovo56233 + denovo56251 + denovo56265 + denovo56344 + denovo56419 + denovo56444 + denovo56455 +
                                 denovo56457 + denovo56464 + denovo56467 + denovo56472 + denovo56492 + denovo56502 + denovo56504 + denovo56505 +
                                 denovo56521 + denovo56569 + denovo56587 + denovo56635 + denovo56639 + denovo56671 + denovo56709 + denovo56761 +
                                 denovo56763 + denovo56772 + denovo56778 + denovo56803 + denovo56830 + denovo56835 + denovo56853 + denovo56871 +
                                 denovo56894 + denovo56914 + denovo56920 + denovo56965 + denovo56982 + denovo57017 + denovo57045 + denovo57138 +
                                 denovo57163 + denovo57166 + denovo57185 + denovo57207 + denovo57275 + denovo57300 + denovo57305 + denovo57353 +
                                 denovo57372 + denovo57421 + denovo57424 + denovo57437 + denovo57459 + denovo57476 + denovo57505 + denovo57514 +
                                 denovo57534 + denovo57570 + denovo57615 + denovo57699 + denovo57713 + denovo57759 + denovo57771 + denovo57838 +
                                 denovo57844 + denovo57868 + denovo57871 + denovo57875 + denovo57911 + denovo57914 + denovo57922 + denovo57972 +
                                 denovo58022 + denovo58034 + denovo58061 + denovo58083 + denovo58090 + denovo58111 + denovo58171 + denovo58245 +
                                 denovo58262 + denovo58328 + denovo58353 + denovo58378 + denovo58390 + denovo58520 + denovo58550 + denovo58552 +
                                 denovo58582 + denovo58600 + denovo58604 + denovo58635 + denovo58678 + denovo58681 + denovo58699 + denovo58709 +
                                 denovo58713 + denovo58717 + denovo58752 + denovo58784 + denovo58808 + denovo58829 + denovo58830 + denovo58848 +
                                 denovo58853 + denovo58856 + denovo58885 + denovo58948 + denovo58950 + denovo58976 + denovo58996 + denovo59030 +
                                 denovo59037 + denovo59060 + denovo59106 + denovo59209 + denovo59213 + denovo59238 + denovo59262 + denovo59276 +
                                 denovo59291 + denovo59294 + denovo59296 + denovo59302 + denovo59310 + denovo59335 + denovo59354 + denovo59370 +
                                 denovo59371 + denovo59372 + denovo59384 + denovo59392 + denovo59406 + denovo59415 + denovo59416 + denovo59430 +
                                 denovo59433 + denovo59437 + denovo59451 + denovo59460 + denovo59575 + denovo59578 + denovo59643 + denovo59658 +
                                 denovo59732 + denovo59744 + denovo59754 + denovo59755 + denovo59779 + denovo59817 + denovo59839 + denovo59844 +
                                 denovo59846 + denovo59872 + denovo59930 + denovo59960 + denovo59970 + denovo60009 + denovo60012 + denovo60068 +
                                 denovo60081 + denovo60094 + denovo60163 + denovo60199 + denovo60243 + denovo60309 + denovo60332 + denovo60358 +
                                 denovo60363 + denovo60379 + denovo60410 + denovo60469 + denovo60532 + denovo60574 + denovo60593 + denovo60650 +
                                 denovo60654 + denovo60692 + denovo60724 + denovo60739 + denovo60775 + denovo60794 + denovo60812 + denovo60831 +
                                 denovo60878 + denovo60912 + denovo60917 + denovo60924 + denovo60946 + denovo60996 + denovo61066 + denovo61073 +
                                 denovo61094 + denovo61142 + denovo61145 + denovo61168 + denovo61218 + denovo61222 + denovo61249 + denovo61354 +
                                 denovo61369 + denovo61376 + denovo61378 + denovo61396 + denovo61402 + denovo61480 + denovo61498 + denovo61502 +
                                 denovo61535 + denovo61540 + denovo61563 + denovo61586 + denovo61635 + denovo61640 + denovo61657 + denovo61702 +
                                 denovo61722 + denovo61728 + denovo61733 + denovo61752 + denovo61787 + denovo61805 + denovo61820 + denovo61841 +
                                 denovo61842 + denovo61849 + denovo61912 + denovo61974 + denovo62011 + denovo62056 + denovo62127 + denovo62153 +
                                 denovo62166 + denovo62178 + denovo62188 + denovo62238 + denovo62253 + denovo62257 + denovo62282 + denovo62284 +
                                 denovo62303 + denovo62380 + denovo62433 + denovo62506 + denovo62540 + denovo62578 + denovo62579 + denovo62590 +
                                 denovo62604 + denovo62650 + denovo62670 + denovo62737 + denovo62783 + denovo62804 + denovo62883 + denovo62905 +
                                 denovo62943 + denovo62948 + denovo62961 + denovo63018 + denovo63037 + denovo63057 + denovo63082 + denovo63101 +
                                 denovo63126 + denovo63143 + denovo63152 + denovo63216 + denovo63228 + denovo63340 + denovo63341 + denovo63389 +
                                 denovo63415 + denovo63429 + denovo63434 + denovo63444 + denovo63482 + denovo63522 + denovo63539 + denovo63552 +
                                 denovo63560 + denovo63586 + denovo63605 + denovo63698 + denovo63738 + denovo63757 + denovo63777 + denovo63857 +
                                 denovo63861 + denovo63863 + denovo63890 + denovo63899 + denovo63907 + denovo63913 + denovo63925 + denovo63973 +
                                 denovo63988 + denovo64010 + denovo64018 + denovo64026 + denovo64067 + denovo64081 + denovo64089 + denovo64093 +
                                 denovo64094 + denovo64121 + denovo64122 + denovo64131 + denovo64140 + denovo64141 + denovo64150 + denovo64208 +
                                 denovo64209 + denovo64236 + denovo64237 + denovo64244 + denovo64260 + denovo64290 + denovo64309 + denovo64324 +
                                 denovo64359 + denovo64381 + denovo64430 + denovo64431 + denovo64455 + denovo64503 + denovo64570 + denovo64592 +
                                 denovo64608 + denovo64618 + denovo64631 + denovo64645 + denovo64747 + denovo64768 + denovo64769 + denovo64797 +
                                 denovo64804 + denovo64806 + denovo64839 + denovo64844 + denovo64852 + denovo64891 + denovo64906 + denovo64963 +
                                 denovo64989 + denovo65009 + denovo65010 + denovo65040 + denovo65065 + denovo65072 + denovo65088 + denovo65096 +
                                 denovo65113 + denovo65127 + denovo65151 + denovo65209 + denovo65233 + denovo65241 + denovo65258 + denovo65285 +
                                 denovo65327 + denovo65328 + denovo65338 + denovo65341 + denovo65354 + denovo65371 + denovo65398 + denovo65424 +
                                 denovo65476 + denovo65477 + denovo65563 + denovo65566 + denovo65638 + denovo65660 + denovo65708 + denovo65753 +
                                 denovo65759 + denovo65795 + denovo65830 + denovo65848 + denovo65855 + denovo65870 + denovo65873 + denovo65901 +
                                 denovo66008 + denovo66081 + denovo66089 + denovo66094 + denovo66129 + denovo66136 + denovo66138 + denovo66177 +
                                 denovo66227 + denovo66234 + denovo66236 + denovo66318 + denovo66337 + denovo66356 + denovo66372 + denovo66380 +
                                 denovo66381 + denovo66418 + denovo66424 + denovo66457 + denovo66467 + denovo66483 + denovo66521 + denovo66524 +
                                 denovo66566 + denovo66615 + denovo66677 + denovo66684 + denovo66689 + denovo66690 + denovo66691 + denovo66698 +
                                 denovo66732 + denovo66741 + denovo66744 + denovo66746 + denovo66749 + denovo66757 + denovo66765 + denovo66776 +
                                 denovo66805 + denovo66806 + denovo66836 + denovo66842 + denovo66876 + denovo66881 + denovo66904 + denovo66906 +
                                 denovo66947 + denovo66982 + denovo66987 + denovo66991 + denovo66996 + denovo67022 + denovo67074 + denovo67086 +
                                 denovo67112 + denovo67138 + denovo67157 + denovo67163 + denovo67181 + denovo67191 + denovo67266 + denovo67297 +
                                 denovo67318 + denovo67320 + denovo67322 + denovo67378 + denovo67406 + denovo67418 + denovo67450 + denovo67470 +
                                 denovo67489 + denovo67492 + denovo67500 + denovo67581 + denovo67583 + denovo67599 + denovo67612 + denovo67657 +
                                 denovo67660 + denovo67662 + denovo67696 + denovo67725 + denovo67728 + denovo67729 + denovo67762 + denovo67769 +
                                 denovo67770 + denovo67807 + denovo67834 + denovo67853 + denovo67871 + denovo68056 + denovo68058 + denovo68107 +
                                 denovo68129 + denovo68136 + denovo68153 + denovo68160 + denovo68170 + denovo68181 + denovo68220 + denovo68227 +
                                 denovo68256 + denovo68285 + denovo68301 + denovo68326 + denovo68332 + denovo68397 + denovo68408 + denovo68409 +
                                 denovo68422 + denovo68432 + denovo68458 + denovo68488 + denovo68494 + denovo68544 + denovo68556 + denovo68577 +
                                 denovo68598 + denovo68640 + denovo68646 + denovo68709 + denovo68724 + denovo68735 + denovo68767 + denovo68791 +
                                 denovo68834 + denovo68839 + denovo68867 + denovo68889 + denovo68914 + denovo68955 + denovo68962 + denovo68980 +
                                 denovo69002 + denovo69004 + denovo69018 + denovo69030 + denovo69112 + denovo69113 + denovo69145 + denovo69181 +
                                 denovo69213 + denovo69281 + denovo69301 + denovo69339 + denovo69380 + denovo69414 + denovo69441 + denovo69444 +
                                 denovo69451 + denovo69453 + denovo69455 + denovo69474 + denovo69488 + denovo69526 + denovo69559 + denovo69576 +
                                 denovo69602 + denovo69617 + denovo69620 + denovo69623 + denovo69626 + denovo69653 + denovo69665 + denovo69705 +
                                 denovo69711 + denovo69717 + denovo69730 + denovo69742 + denovo69770 + denovo69772 + denovo69773 + denovo69781 +
                                 denovo69839 + denovo69845 + denovo69857 + denovo69866 + denovo69909 + denovo69945 + denovo69972 + denovo69980 +
                                 denovo70028 + denovo70083 + denovo70094 + denovo70110 + denovo70238 + denovo70241 + denovo70253 + denovo70260 +
                                 denovo70264 + denovo70290 + denovo70326 + denovo70346 + denovo70405 + denovo70420 + denovo70434 + denovo70491 +
                                 denovo70519 + denovo70615 + denovo70651 + denovo70654 + denovo70671 + denovo70684 + denovo70701 + denovo70782 +
                                 denovo70806 + denovo70816 + denovo70845 + denovo70854 + denovo70857 + denovo70861 + denovo70900 + denovo70906 +
                                 denovo70946 + denovo70950 + denovo70968 + denovo70990 + denovo71006 + denovo71008 + denovo71014 + denovo71030 +
                                 denovo71050 + denovo71102 + denovo71142 + denovo71236 + denovo71278 + denovo71315 + denovo71333 + denovo71348 +
                                 denovo71359 + denovo71366 + denovo71372 + denovo71403 + denovo71426 + denovo71436 + denovo71438 + denovo71440 +
                                 denovo71441 + denovo71446 + denovo71449 + denovo71482 + denovo71516 + denovo71525 + denovo71547 + denovo71576 +
                                 denovo71594 + denovo71598 + denovo71604 + denovo71649 + denovo71666 + denovo71672 + denovo71674 + denovo71677 +
                                 denovo71810 + denovo71840 + denovo71843 + denovo71913 + denovo71953 + denovo71958 + denovo72000 + denovo72053 +
                                 denovo72062 + denovo72109 + denovo72153 + denovo72161 + denovo72166 + denovo72174 + denovo72199 + denovo72204 +
                                 denovo72210 + denovo72233 + denovo72284 + denovo72320 + denovo72332 + denovo72457 + denovo72462 + denovo72466 +
                                 denovo72474 + denovo72485 + denovo72496 + denovo72517 + denovo72554 + denovo72619 + denovo72634 + denovo72695 +
                                 denovo72739 + denovo72752 + denovo72767 + denovo72820 + denovo72826 + denovo72870 + denovo72892 + denovo72940 +
                                 denovo72952 + denovo72971 + denovo73038 + denovo73045 + denovo73105 + denovo73190 + denovo73206 + denovo73214 +
                                 denovo73217 + denovo73347 + denovo73349 + denovo73400 + denovo73404 + denovo73412 + denovo73422 + denovo73448 +
                                 denovo73451 + denovo73457 + denovo73461 + denovo73581 + denovo73630 + denovo73674 + denovo73710 + denovo73713 +
                                 denovo73743 + denovo73798 + denovo73816 + denovo73853 + denovo73890 + denovo73902 + denovo73915 + denovo73924 +
                                 denovo73927 + denovo73930 + denovo73955 + denovo73973 + denovo74010 + denovo74063 + denovo74064 + denovo74075 +
                                 denovo74082 + denovo74097 + denovo74129 + denovo74165 + denovo74194 + denovo74276 + denovo74285 + denovo74297 +
                                 denovo74301 + denovo74336 + denovo74435 + denovo74440 + denovo74444 + denovo74478 + denovo74493 + denovo74494 +
                                 denovo74516 + denovo74523 + denovo74525 + denovo74564 + denovo74585 + denovo74676 + denovo74692 + denovo74702 +
                                 denovo74703 + denovo74781 + denovo74829 + denovo74852 + denovo74929 + denovo74938 + denovo74948 + denovo74955 +
                                 denovo74957 + denovo74963 + denovo74966 + denovo75041 + denovo75077 + denovo75116 + denovo75125 + denovo75134 +
                                 denovo75141 + denovo75196 + denovo75213 + denovo75224 + denovo75239 + denovo75252 + denovo75273 + denovo75353 +
                                 denovo75373 + denovo75420 + denovo75444 + denovo75545 + denovo75585 + denovo75648 + denovo75653 + denovo75659 +
                                 denovo75670 + denovo75694 + denovo75712 + denovo75717 + denovo75718 + denovo75779 + denovo75796 + denovo75851 +
                                 denovo75859 + denovo75891 + denovo75906 + denovo75910 + denovo75947 + denovo76024 + denovo76040 + denovo76046 +
                                 denovo76064 + denovo76088 + denovo76094 + denovo76127 + denovo76133 + denovo76135 + denovo76139 + denovo76188 +
                                 denovo76216 + denovo76258 + denovo76330 + denovo76333 + denovo76340 + denovo76345 + denovo76399 + denovo76428 +
                                 denovo76437 + denovo76440 + denovo76441 + denovo76469 + denovo76479 + denovo76511 + denovo76580 + denovo76594 +
                                 denovo76603 + denovo76605 + denovo76735 + denovo76762 + denovo76812 + denovo76824 + denovo76837 + denovo76861 +
                                 denovo76924 + denovo76940 + denovo76966 + denovo76989 + denovo77043 + denovo77089 + denovo77115 + denovo77134 +
                                 denovo77166 + denovo77215 + denovo77353 + denovo77362 + denovo77377 + denovo77381 + denovo77392 + denovo77399 +
                                 denovo77415 + denovo77462 + denovo77504 + denovo77512 + denovo77524 + denovo77540 + denovo77590 + denovo77626 +
                                 denovo77657 + denovo77668 + denovo77675 + denovo77728 + denovo77743 + denovo77786 + denovo77811 + denovo77834 +
                                 denovo77836 + denovo77858 + denovo77859 + denovo77873 + denovo77882 + denovo77905 + denovo77917 + denovo77938 +
                                 denovo77945 + denovo78076 + denovo78079 + denovo78113 + denovo78120 + denovo78132 + denovo78134 + denovo78175 +
                                 denovo78194 + denovo78285 + denovo78304 + denovo78331 + denovo78347 + denovo78352 + denovo78354 + denovo78360 +
                                 denovo78372 + denovo78388 + denovo78391 + denovo78415 + denovo78446 + denovo78480 + denovo78488 + denovo78504 +
                                 denovo78536 + denovo78542 + denovo78557 + denovo78559 + denovo78572 + denovo78599 + denovo78671 + denovo78707 +
                                 denovo78708 + denovo78765 + denovo78794 + denovo78818 + denovo78873 + denovo78893 + denovo78895 + denovo78897 +
                                 denovo78963 + denovo78976 + denovo78985 + denovo79050 + denovo79057 + denovo79092 + denovo79121 + denovo79153 +
                                 denovo79198 + denovo79238 + denovo79254 + denovo79312 + denovo79320 + denovo79348 + denovo79390 + denovo79459 +
                                 denovo79494 + denovo79502 + denovo79525 + denovo79629 + denovo79677 + denovo79683 + denovo79693 + denovo79731 +
                                 denovo79740 + denovo79749 + denovo79757 + denovo79766 + denovo79797 + denovo79836 + denovo79865 + denovo79899 +
                                 denovo79905 + denovo79969 + denovo80005 + denovo80009 + denovo80018 + denovo80024 + denovo80046 + denovo80055 +
                                 denovo80059 + denovo80061 + denovo80143 + denovo80178 + denovo80188 + denovo80197 + denovo80238 + denovo80257 +
                                 denovo80271 + denovo80301 + denovo80335 + denovo80358 + denovo80386 + denovo80398 + denovo80411 + denovo80430 +
                                 denovo80439 + denovo80492 + denovo80501 + denovo80535 + denovo80539 + denovo80582 + denovo80588 + denovo80592 +
                                 denovo80658 + denovo80660 + denovo80661 + denovo80729 + denovo80740 + denovo80820 + denovo80836 + denovo80839 +
                                 denovo80877 + denovo80889 + denovo80893 + denovo80967 + denovo80987 + denovo81100 + denovo81109 + denovo81123 +
                                 denovo81188 + denovo81191 + denovo81194 + denovo81213 + denovo81232 + denovo81312 + denovo81313 + denovo81316 +
                                 denovo81356 + denovo81406 + denovo81439 + denovo81567 + denovo81571 + denovo81603 + denovo81609 + denovo81613 +
                                 denovo81629 + denovo81666 + denovo81723 + denovo81754 + denovo81769 + denovo81843 + denovo81850 + denovo81870 +
                                 denovo81879 + denovo81920 + denovo81938 + denovo81974 + denovo81977 + denovo81982 + denovo82076 + denovo82092 +
                                 denovo82128 + denovo82139 + denovo82141 + denovo82154 + denovo82155 + denovo82156 + denovo82186 + denovo82209 +
                                 denovo82234 + denovo82317 + denovo82318 + denovo82322 + denovo82434 + denovo82455 + denovo82462 + denovo82545 +
                                 denovo82600 + denovo82604 + denovo82633 + denovo82639 + denovo82649 + denovo82691 + denovo82694 + denovo82700 +
                                 denovo82747 + denovo82754 + denovo82775 + denovo82789 + denovo82791 + denovo82864 + denovo82875 + denovo82880 +
                                 denovo82884 + denovo82930 + denovo82993 + denovo82995 + denovo82997 + denovo82999 + denovo83002 + denovo83005 +
                                 denovo83021 + denovo83062 + denovo83119 + denovo83120 + denovo83131 + denovo83220 + denovo83240 + denovo83256 +
                                 denovo83284 + denovo83290 + denovo83307 + denovo83314 + denovo83448 + denovo83458 + denovo83479 + denovo83506 +
                                 denovo83513 + denovo83531 + denovo83545 + denovo83559 + denovo83574 + denovo83575 + denovo83588 + denovo83594 +
                                 denovo83605 + denovo83622 + denovo83625 + denovo83669 + denovo83694 + denovo83720 + denovo83725 + denovo83771 +
                                 denovo83799 + denovo83822 + denovo83855 + denovo83886 + denovo83913 + denovo83945 + denovo83991 + denovo84002 +
                                 denovo84020 + denovo84076 + denovo84125 + denovo84184 + denovo84224 + denovo84230 + denovo84304 + denovo84346 +
                                 denovo84364 + denovo84447 + denovo84463 + denovo84465 + denovo84489 + denovo84522 + denovo84545 + denovo84557 +
                                 denovo84620 + denovo84643 + denovo84660 + denovo84664 + denovo84665 + denovo84750 + denovo84768 + denovo84854 +
                                 denovo84860 + denovo84928 + denovo85016 + denovo85017 + denovo85043 + denovo85049 + denovo85080 + denovo85124 +
                                 denovo85145 + denovo85188 + denovo85191 + denovo85216 + denovo85235 + denovo85236 + denovo85241 + denovo85250 +
                                 denovo85253 + denovo85293 + denovo85297 + denovo85313 + denovo85326 + denovo85331 + denovo85404 + denovo85457 +
                                 denovo85481 + denovo85510 + denovo85514 + denovo85534 + denovo85537 + denovo85545 + denovo85567 + denovo85574 +
                                 denovo85616 + denovo85622 + denovo85650 + denovo85670 + denovo85722 + denovo85798 + denovo85851 + denovo85853 +
                                 denovo85883 + denovo85933 + denovo85937 + denovo85965 + denovo85966 + denovo85979 + denovo85993 + denovo85994 +
                                 denovo86079 + denovo86086 + denovo86135 + denovo86151 + denovo86182 + denovo86205 + denovo86212 + denovo86221 +
                                 denovo86254 + denovo86292 + denovo86319 + denovo86330 + denovo86337 + denovo86340 + denovo86350 + denovo86362 +
                                 denovo86419 + denovo86426 + denovo86531 + denovo86614 + denovo86616 + denovo86637 + denovo86655 + denovo86659 +
                                 denovo86680 + denovo86698 + denovo86712 + denovo86717 + denovo86733 + denovo86758 + denovo86767 + denovo86784 +
                                 denovo86797 + denovo86798 + denovo86853 + denovo86867 + denovo86893 + denovo86896 + denovo86911 + denovo86916 +
                                 denovo86964 + denovo86967 + denovo86991 + denovo87003 + denovo87007 + denovo87036 + denovo87064 + denovo87071 +
                                 denovo87190 + denovo87226 + denovo87257 + denovo87320 + denovo87382 + denovo87400 + denovo87447 + denovo87452 +
                                 denovo87483 + denovo87503 + denovo87516 + denovo87524 + denovo87533 + denovo87539 + denovo87551 + denovo87593 +
                                 denovo87671 + denovo87673 + denovo87675 + denovo87713 + denovo87727 + denovo87760 + denovo87761 + denovo87812 +
                                 denovo87895 + denovo87962 + denovo87981 + denovo88069 + denovo88077 + denovo88083 + denovo88101 + denovo88115 +
                                 denovo88163 + denovo88212 + denovo88216 + denovo88234 + denovo88246 + denovo88285 + denovo88334 + denovo88346 +
                                 denovo88366 + denovo88404 + denovo88419 + denovo88449 + denovo88466 + denovo88491 + denovo88500 + denovo88507 +
                                 denovo88612 + denovo88624 + denovo88679 + denovo88688 + denovo88710 + denovo88730 + denovo88771 + denovo88799 +
                                 denovo88817 + denovo88837 + denovo88911 + denovo88934 + denovo89012 + denovo89048 + denovo89056 + denovo89069 +
                                 denovo89074 + denovo89077 + denovo89091 + denovo89096 + denovo89117 + denovo89163 + denovo89168 + denovo89172 +
                                 denovo89187 + denovo89192 + denovo89198 + denovo89227 + denovo89231 + denovo89263 + denovo89272 + denovo89297 +
                                 denovo89309 + denovo89343 + denovo89370 + denovo89386 + denovo89417 + denovo89431 + denovo89444 + denovo89455 +
                                 denovo89457 + denovo89459 + denovo89513 + denovo89561 + denovo89574 + denovo89577 + denovo89598 + denovo89618 +
                                 denovo89648 + denovo89662 + denovo89671 + denovo89798 + denovo89803 + denovo89808 + denovo89828 + denovo89846 +
                                 denovo89870 + denovo89901 + denovo89905 + denovo89939 + denovo89975 + denovo90065 + denovo90084 + denovo90115 +
                                 denovo90159 + denovo90175 + denovo90194 + denovo90211 + denovo90264 + denovo90277 + denovo90319 + denovo90344 +
                                 denovo90353 + denovo90383 + denovo90404 + denovo90414 + denovo90475 + denovo90509 + denovo90530 + denovo90551 +
                                 denovo90553 + denovo90563 + denovo90564 + denovo90689 + denovo90693 + denovo90696 + denovo90707 + denovo90710 +
                                 denovo90713 + denovo90721 + denovo90734 + denovo90796 + denovo90806 + denovo90823 + denovo90853 + denovo90890 +
                                 denovo90930 + denovo91013 + denovo91033 + denovo91090 + denovo91095 + denovo91147 + denovo91161 + denovo91186 +
                                 denovo91187 + denovo91205 + denovo91206 + denovo91214 + denovo91227 + denovo91228 + denovo91302 + denovo91346 +
                                 denovo91353 + denovo91357 + denovo91389 + denovo91393 + denovo91435 + denovo91466 + denovo91482 + denovo91498 +
                                 denovo91531 + denovo91568 + denovo91582 + denovo91608 + denovo91623 + denovo91671 + denovo91676 + denovo91677 +
                                 denovo91707 + denovo91729 + denovo91743 + denovo91747 + denovo91773 + denovo91779 + denovo91799 + denovo91808 +
                                 denovo91813 + denovo91826 + denovo91840 + denovo91945 + denovo91951 + denovo91956 + denovo91980 + denovo92003 +
                                 denovo92017 + denovo92024 + denovo92025 + denovo92044 + denovo92046 + denovo92077 + denovo92124 + denovo92132 +
                                 denovo92169 + denovo92172 + denovo92197 + denovo92214 + denovo92235 + denovo92237 + denovo92277 + denovo92284 +
                                 denovo92302 + denovo92327 + denovo92387 + denovo92401 + denovo92407 + denovo92478 + denovo92490 + denovo92505 +
                                 denovo92507 + denovo92512 + denovo92588 + denovo92589 + denovo92607 + denovo92616 + denovo92622 + denovo92625 +
                                 denovo92634 + denovo92657 + denovo92677 + denovo92705 + denovo92731 + denovo92773 + denovo92786 + denovo92796 +
                                 denovo92809 + denovo92875 + denovo92885 + denovo92918 + denovo92978 + denovo92986 + denovo93012 + denovo93018 +
                                 denovo93057 + denovo93072 + denovo93106 + denovo93146 + denovo93156 + denovo93165 + denovo93183 + denovo93192 +
                                 denovo93221 + denovo93256 + denovo93292 + denovo93321 + denovo93346 + denovo93350 + denovo93355 + denovo93363 +
                                 denovo93376 + denovo93391 + denovo93405 + denovo93435 + denovo93446 + denovo93457 + denovo93473 + denovo93480 +
                                 denovo93488 + denovo93730 + denovo93778 + denovo93786 + denovo93824 + denovo93859 + denovo93862 + denovo93915 +
                                 denovo93918 + denovo93932 + denovo93935 + denovo93936 + denovo93957 + denovo93985 + denovo94010 + denovo94120 +
                                 denovo94141 + denovo94146 + denovo94151 + denovo94152 + denovo94170 + denovo94185 + denovo94209 + denovo94213 +
                                 denovo94217 + denovo94230 + denovo94306 + denovo94324 + denovo94332 + denovo94335 + denovo94336 + denovo94402 +
                                 denovo94530 + denovo94630 + denovo94632 + denovo94644 + denovo94715 + denovo94741 + denovo94761 + denovo94762 +
                                 denovo94784 + denovo94804 + denovo94840 + denovo94848 + denovo94862 + denovo94880 + denovo94899 + denovo94928 +
                                 denovo94978 + denovo94988 + denovo95005 + denovo95033 + denovo95057 + denovo95084 + denovo95087 + denovo95115 +
                                 denovo95131 + denovo95148 + denovo95157 + denovo95167 + denovo95172 + denovo95181 + denovo95233 + denovo95292 +
                                 denovo95302 + denovo95315 + denovo95370 + denovo95387 + denovo95404 + denovo95413 + denovo95426 + denovo95438 +
                                 denovo95440 + denovo95463 + denovo95496 + denovo95521 + denovo95527 + denovo95588 + denovo95595 + denovo95596 +
                                 denovo95604 + denovo95614 + denovo95627 + denovo95648 + denovo95676 + denovo95763 + denovo95771 + denovo95820 +
                                 denovo95856 + denovo95900 + denovo95910 + denovo96003 + denovo96018 + denovo96022 + denovo96128 + denovo96134 +
                                 denovo96174 + denovo96188 + denovo96263 + denovo96302 + denovo96320 + denovo96338 + denovo96406 + denovo96433 +
                                 denovo96470 + denovo96497 + denovo96506 + denovo96532 + denovo96534 + denovo96654 + denovo96659 + denovo96677 +
                                 denovo96688 + denovo96702 + denovo96712 + denovo96719 + denovo96769 + denovo96816 + denovo96855 + denovo96940 +
                                 denovo96943 + denovo96948 + denovo96958 + denovo96966 + denovo96971 + denovo97027 + denovo97057 + denovo97076 +
                                 denovo97093 + denovo97106 + denovo97151 + denovo97154 + denovo97185 + denovo97204 + denovo97267 + denovo97334 +
                                 denovo97372 + denovo97396 + denovo97406 + denovo97410 + denovo97414 + denovo97432 + denovo97452 + denovo97510 +
                                 denovo97542 + denovo97577 + denovo97611 + denovo97618 + denovo97657 + denovo97686 + denovo97700 + denovo97710 +
                                 denovo97713 + denovo97721 + denovo97749 + denovo97803 + denovo97807 + denovo97844 + denovo97849 + denovo97874 +
                                 denovo97886 + denovo97913 + denovo97931 + denovo97954 + denovo98038 + denovo98083 + denovo98129 + denovo98155 +
                                 denovo98175 + denovo98183 + denovo98189 + denovo98235 + denovo98249 + denovo98257 + denovo98285 + denovo98308 +
                                 denovo98352 + denovo98363 + denovo98382 + denovo98428 + denovo98433 + denovo98554 + denovo98631 + denovo98636 +
                                 denovo98661 + denovo98663 + denovo98698 + denovo98709 + denovo98739 + denovo98749 + denovo98760 + denovo98762 +
                                 denovo98779 + denovo98809 + denovo98841 + denovo98886 + denovo98892 + denovo98898 + denovo98900 + denovo98907 +
                                 denovo98935 + denovo98942 + denovo98967 + denovo98992 + denovo99023 + denovo99074 + denovo99109 + denovo99118 +
                                 denovo99127 + denovo99133 + denovo99144 + denovo99173 + denovo99257 + denovo99281 + denovo99399 + denovo99450 +
                                 denovo99451 + denovo99452 + denovo99455 + denovo99559 + denovo99561 + denovo99563 + denovo99589 + denovo99606 +
                                 denovo99690 + denovo99736 + denovo99758 + denovo99780 + denovo99785 + denovo99854 + denovo99876 + denovo99896 +
                                 denovo99927 + denovo99929 + denovo99965 + denovo99978 + denovo100004 + denovo100052 + denovo100080 + denovo100132 +
                                 denovo100139 + denovo100160 + denovo100176 + denovo100188 + denovo100208 + denovo100251 + denovo100289 + denovo100310 +
                                 denovo100334 + denovo100360 + denovo100408 + denovo100440 + denovo100446 + denovo100473 + denovo100534 + denovo100538 +
                                 denovo100550 + denovo100563 + denovo100576 + denovo100594 + denovo100613 + denovo100617 + denovo100662 + denovo100783 +
                                 denovo100889 + denovo100892 + denovo100894 + denovo100906 + denovo100923 + denovo100953 + denovo100992 + denovo101027 +
                                 denovo101034 + denovo101049 + denovo101069 + denovo101094 + denovo101105 + denovo101112 + denovo101119 + denovo101164 +
                                 denovo101173 + denovo101225 + denovo101239 + denovo101246 + denovo101252 + denovo101288 + denovo101296 + denovo101357 +
                                 denovo101382 + denovo101405 + denovo101423 + denovo101444 + denovo101458 + denovo101472 + denovo101503 + denovo101505 +
                                 denovo101517 + denovo101525 + denovo101598 + denovo101640 + denovo101644 + denovo101669 + denovo101744 + denovo101750 +
                                 denovo101784 + denovo101803 + denovo101815 + denovo101821 + denovo101830 + denovo101876 + denovo101880 + denovo101896 +
                                 denovo101941 + denovo101953 + denovo101964 + denovo102006 + denovo102014 + denovo102021 + denovo102030 + denovo102108 +
                                 denovo102111 + denovo102168 + denovo102228 + denovo102239 + denovo102244 + denovo102260 + denovo102262 + denovo102288 +
                                 denovo102303 + denovo102329 + denovo102335 + denovo102366 + denovo102427 + denovo102433 + denovo102468 + denovo102480 +
                                 denovo102494 + denovo102562 + denovo102572 + denovo102600 + denovo102604 + denovo102609 + denovo102619 + denovo102620 +
                                 denovo102633 + denovo102681 + denovo102682 + denovo102699 + denovo102728 + denovo102737 + denovo102789 + denovo102834 +
                                 denovo102865 + denovo102866 + denovo102887 + denovo102905 + denovo102982 + denovo103023 + denovo103036 + denovo103060 +
                                 denovo103083 + denovo103108 + denovo103133 + denovo103135 + denovo103138 + denovo103153 + denovo103165 + denovo103168 +
                                 denovo103172 + denovo103181 + denovo103223 + denovo103241 + denovo103277 + denovo103357 + denovo103358 + denovo103425 +
                                 denovo103444 + denovo103454 + denovo103470 + denovo103507 + denovo103509 + denovo103516 + denovo103544 + denovo103547 +
                                 denovo103555 + denovo103557 + denovo103569 + denovo103586 + denovo103611 + denovo103617 + denovo103638 + denovo103667 +
                                 denovo103720 + denovo103745 + denovo103781 + denovo103838 + denovo103864 + denovo103914 + denovo104042 + denovo104044 +
                                 denovo104051 + denovo104066 + denovo104079 + denovo104093 + denovo104100 + denovo104104 + denovo104107 + denovo104115 +
                                 denovo104127 + denovo104143 + denovo104177 + denovo104194 + denovo104200 + denovo104210 + denovo104214 + denovo104238 +
                                 denovo104252 + denovo104283 + denovo104328 + denovo104329 + denovo104364 + denovo104397 + denovo104419 + denovo104421 +
                                 denovo104422 + denovo104437 + denovo104494 + denovo104562 + denovo104591 + denovo104618 + denovo104633 + denovo104661 +
                                 denovo104667 + denovo104699 + denovo104759 + denovo104792 + denovo104801 + denovo104836 + denovo104873 + denovo104881 +
                                 denovo104885 + denovo104930 + denovo104936 + denovo104938 + denovo104986 + denovo105076 + denovo105086 + denovo105088 +
                                 denovo105146 + denovo105158 + denovo105231 + denovo105241 + denovo105341 + denovo105345 + denovo105354 + denovo105361 +
                                 denovo105367 + denovo105396 + denovo105417 + denovo105445 + denovo105492 + denovo105509 + denovo105545 + denovo105588 +
                                 denovo105623 + denovo105658 + denovo105661 + denovo105731 + denovo105732 + denovo105737 + denovo105755 + denovo105783 +
                                 denovo105799 + denovo105801 + denovo105804 + denovo105823 + denovo105948 + denovo105954 + denovo105956 + denovo105962 +
                                 denovo105965 + denovo105982 + denovo105985 + denovo105990 + denovo106067 + denovo106077 + denovo106093 + denovo106145 +
                                 denovo106160 + denovo106224 + denovo106281 + denovo106288 + denovo106289 + denovo106299 + denovo106381 + denovo106384 +
                                 denovo106392 + denovo106393 + denovo106431 + denovo106492 + denovo106496 + denovo106521 + denovo106532 + denovo106566 +
                                 denovo106578 + denovo106585 + denovo106614 + denovo106620 + denovo106627 + denovo106663 + denovo106671 + denovo106770 +
                                 denovo106801 + denovo106802 + denovo106804 + denovo106840 + denovo106857 + denovo106868 + denovo106875 + denovo106895 +
                                 denovo106904 + denovo106950 + denovo106999 + denovo107014 + denovo107017 + denovo107072 + denovo107105 + denovo107106 +
                                 denovo107202 + denovo107220 + denovo107255 + denovo107294 + denovo107340 + denovo107374 + denovo107388 + denovo107397 +
                                 denovo107403 + denovo107409 + denovo107431 + denovo107464 + denovo107470 + denovo107499 + denovo107500 + denovo107508 +
                                 denovo107548 + denovo107571 + denovo107586 + denovo107612 + denovo107629 + denovo107631 + denovo107651 + denovo107669 +
                                 denovo107686 + denovo107732 + denovo107752 + denovo107758 + denovo107873 + denovo107882 + denovo107921 + denovo107930 +
                                 denovo107970 + denovo108021 + denovo108047 + denovo108074 + denovo108096 + denovo108107 + denovo108223 + denovo108249 +
                                 denovo108251 + denovo108291 + denovo108316 + denovo108355 + denovo108432 + denovo108433 + denovo108457 + denovo108483 +
                                 denovo108502 + denovo108527 + denovo108534 + denovo108552 + denovo108572 + denovo108588 + denovo108596 + denovo108601 +
                                 denovo108623 + denovo108632 + denovo108646 + denovo108671 + denovo108704 + denovo108723 + denovo108743 + denovo108832 +
                                 denovo108846 + denovo108912 + denovo108921 + denovo108932 + denovo108951 + denovo109023 + denovo109025 + denovo109034 +
                                 denovo109037 + denovo109061 + denovo109071 + denovo109094 + denovo109104 + denovo109109 + denovo109203 + denovo109266 +
                                 denovo109301 + denovo109324 + denovo109347 + denovo109398 + denovo109524 + denovo109556 + denovo109659 + denovo109682 +
                                 denovo109722 + denovo109734 + denovo109751 + denovo109754 + denovo109758 + denovo109851 + denovo109853 + denovo109863 +
                                 denovo109877 + denovo109880 + denovo109899 + denovo109921 + denovo109944 + denovo109977 + denovo109982 + denovo109992 +
                                 denovo110026 + denovo110034 + denovo110036 + denovo110079 + denovo110093 + denovo110121 + denovo110141 + denovo110157 +
                                 denovo110173 + denovo110180 + denovo110228 + denovo110234 + denovo110235 + denovo110248 + denovo110300 + denovo110312 +
                                 denovo110357 + denovo110371 + denovo110418 + denovo110483 + denovo110521 + denovo110533 + denovo110550 + denovo110586 +
                                 denovo110608 + denovo110610 + denovo110651 + denovo110659 + denovo110699 + denovo110705 + denovo110713 + denovo110716 +
                                 denovo110721 + denovo110764 + denovo110798 + denovo110817 + denovo110826 + denovo110832 + denovo110881 + denovo110902 +
                                 denovo110932 + denovo110945 + denovo110950 + denovo110966 + denovo110969 + denovo110988 + denovo111002 + denovo111028 +
                                 denovo111048 + denovo111074 + denovo111078 + denovo111095 + denovo111110 + denovo111126 + denovo111145 + denovo111146 +
                                 denovo111159 + denovo111166 + denovo111173 + denovo111197 + denovo111262 + denovo111274 + denovo111294 + denovo111300 +
                                 denovo111343 + denovo111347 + denovo111370 + denovo111376 + denovo111383 + denovo111411 + denovo111415 + denovo111497 +
                                 denovo111506 + denovo111508 + denovo111516 + denovo111528 + denovo111565 + denovo111596 + denovo111648 + denovo111664 +
                                 denovo111666 + denovo111669 + denovo111674 + denovo111728 + denovo111730 + denovo111764 + denovo111786 + denovo111817 +
                                 denovo111841 + denovo111851 + denovo111906 + denovo111939 + denovo111950 + denovo111960 + denovo111967 + denovo112075 +
                                 denovo112108 + denovo112121 + denovo112146 + denovo112211 + denovo112235 + denovo112262 + denovo112278 + denovo112317 +
                                 denovo112321 + denovo112322 + denovo112331 + denovo112345 + denovo112346 + denovo112381 + denovo112431 + denovo112451 +
                                 denovo112490 + denovo112507 + denovo112552 + denovo112560 + denovo112581 + denovo112590 + denovo112593 + denovo112615 +
                                 denovo112632 + denovo112655 + denovo112706 + denovo112715 + denovo112728 + denovo112771 + denovo112813 + denovo112818 +
                                 denovo112829 + denovo112837 + denovo112868 + denovo112884 + denovo112891 + denovo112893 + denovo112896 + denovo112903 +
                                 denovo112918 + denovo112970 + denovo112982 + denovo112987 + denovo113033 + denovo113077 + denovo113103 + denovo113109 +
                                 denovo113144 + denovo113176 + denovo113237 + denovo113257 + denovo113275 + denovo113279 + denovo113284 + denovo113304 +
                                 denovo113314 + denovo113316 + denovo113345 + denovo113404 + denovo113413 + denovo113416 + denovo113425 + denovo113428 +
                                 denovo113452 + denovo113471 + denovo113502 + denovo113545 + denovo113552 + denovo113567 + denovo113572 + denovo113592 +
                                 denovo113598 + denovo113616 + denovo113634 + denovo113638 + denovo113680 + denovo113685 + denovo113689 + denovo113698 +
                                 denovo113757 + denovo113760 + denovo113769 + denovo113770 + denovo113780 + denovo113786 + denovo113790 + denovo113791 +
                                 denovo113808 + denovo113815 + denovo113831 + denovo113859 + denovo113877 + denovo113881 + denovo113885 + denovo113925 +
                                 denovo113930 + denovo113961 + denovo113996 + denovo114061 + denovo114081 + denovo114191 + denovo114192 + denovo114215 +
                                 denovo114260 + denovo114268 + denovo114271 + denovo114290 + denovo114302 + denovo114316 + denovo114346 + denovo114370 +
                                 denovo114378 + denovo114380 + denovo114397 + denovo114413 + denovo114423 + denovo114461 + denovo114483 + denovo114488 +
                                 denovo114508 + denovo114520 + denovo114527 + denovo114579 + denovo114582 + denovo114593 + denovo114604 + denovo114678 +
                                 denovo114709 + denovo114712 + denovo114762 + denovo114780 + denovo114786 + denovo114821 + denovo114868 + denovo114869 +
                                 denovo114898 + denovo114908 + denovo114920 + denovo114926 + denovo114977 + denovo115038 + denovo115043 + denovo115090 +
                                 denovo115118 + denovo115175 + denovo115212 + denovo115219 + denovo115251 + denovo115298 + denovo115300 + denovo115302 +
                                 denovo115358 + denovo115362 + denovo115374 + denovo115377 + denovo115394 + denovo115399 + denovo115415 + denovo115421 +
                                 denovo115448 + denovo115451 + denovo115472 + denovo115474 + denovo115511 + denovo115518 + denovo115532 + denovo115536 +
                                 denovo115551 + denovo115562 + denovo115588 + denovo115624 + denovo115629 + denovo115665 + denovo115669 + denovo115695 +
                                 denovo115708 + denovo115722 + denovo115744 + denovo115792 + denovo115803 + denovo115864 + denovo115973 + denovo115983 +
                                 denovo116020 + denovo116045 + denovo116056 + denovo116061 + denovo116071 + denovo116092 + denovo116096 + denovo116137 +
                                 denovo116145 + denovo116156 + denovo116191 + denovo116193 + denovo116197 + denovo116223 + denovo116253 + denovo116277 +
                                 denovo116279 + denovo116317 + denovo116330 + denovo116341 + denovo116357 + denovo116374 + denovo116378 + denovo116454 +
                                 denovo116482 + denovo116517 + denovo116518 + denovo116544 + denovo116549 + denovo116578 + denovo116583 + denovo116619 +
                                 denovo116622 + denovo116642 + denovo116662 + denovo116670 + denovo116689 + denovo116718 + denovo116724 + denovo116728 +
                                 denovo116777 + denovo116784 + denovo116794 + denovo116796 + denovo116841 + denovo116843 + denovo116850 + denovo116866 +
                                 denovo116882 + denovo116891 + denovo116892 + denovo116938 + denovo116961 + denovo116966 + denovo116971 + denovo116989 +
                                 denovo116990 + denovo117001 + denovo117061 + denovo117112 + denovo117148 + denovo117157 + denovo117165 + denovo117201 +
                                 denovo117252 + denovo117256 + denovo117288 + denovo117320 + denovo117391 + denovo117394 + denovo117422 + denovo117462 +
                                 denovo117493 + denovo117517 + denovo117537 + denovo117558 + denovo117581 + denovo117599 + denovo117604 + denovo117633 +
                                 denovo117651 + denovo117652 + denovo117657 + denovo117660 + denovo117725 + denovo117742 + denovo117765 + denovo117777 +
                                 denovo117804 + denovo117812 + denovo117821 + denovo117831 + denovo117896 + denovo117904 + denovo117919 + denovo117923 +
                                 denovo117941 + denovo117962 + denovo117988 + denovo118061 + denovo118084 + denovo118146 + denovo118158 + denovo118174 +
                                 denovo118189 + denovo118205 + denovo118252 + denovo118271 + denovo118272 + denovo118277 + denovo118303 + denovo118317 +
                                 denovo118321 + denovo118325 + denovo118355 + denovo118364 + denovo118412 + denovo118432 + denovo118460 + denovo118477 +
                                 denovo118486 + denovo118506 + denovo118525 + denovo118604 + denovo118636 + denovo118647 + denovo118670 + denovo118714 +
                                 denovo118741 + denovo118771 + denovo118799 + denovo118807 + denovo118820 + denovo118851 + denovo118877 + denovo118880 +
                                 denovo118886 + denovo118916 + denovo118923 + denovo118954 + denovo118957 + denovo118967 + denovo118968 + denovo118988 +
                                 denovo118989 + denovo118994 + denovo119047 + denovo119059 + denovo119068 + denovo119072 + denovo119086 + denovo119104 +
                                 denovo119115 + denovo119241 + denovo119250 + denovo119278 + denovo119320 + denovo119333 + denovo119341 + denovo119342 +
                                 denovo119366 + denovo119394 + denovo119466 + denovo119474 + denovo119546 + denovo119567 + denovo119590 + denovo119613 +
                                 denovo119640 + denovo119671 + denovo119678 + denovo119680 + denovo119700 + denovo119730 + denovo119742 + denovo119758 +
                                 denovo119769 + denovo119770 + denovo119779 + denovo119793 + denovo119826 + denovo119832 + denovo119869 + denovo119883 +
                                 denovo119958 + denovo119980 + denovo120005 + denovo120009 + denovo120013 + denovo120057 + denovo120068 + denovo120072 +
                                 denovo120076 + denovo120082 + denovo120109 + denovo120140 + denovo120152 + denovo120159 + denovo120210 + denovo120266 +
                                 denovo120269 + denovo120294 + denovo120309 + denovo120316 + denovo120358 + denovo120361 + denovo120363 + denovo120379 +
                                 denovo120401 + denovo120469 + denovo120495 + denovo120499 + denovo120503 + denovo120511 + denovo120531 + denovo120544 +
                                 denovo120590 + denovo120624 + denovo120636 + denovo120642 + denovo120653 + denovo120672 + denovo120685 + denovo120695 +
                                 denovo120706 + denovo120710 + denovo120779 + denovo120786 + denovo120790 + denovo120799 + denovo120828 + denovo120840 +
                                 denovo120845 + denovo120848 + denovo120858 + denovo120901 + denovo121027 + denovo121064 + denovo121110 + denovo121116 +
                                 denovo121120 + denovo121137 + denovo121219 + denovo121223 + denovo121269 + denovo121275 + denovo121294 + denovo121308 +
                                 denovo121312 + denovo121381 + denovo121384 + denovo121386 + denovo121401 + denovo121420 + denovo121433 + denovo121436 +
                                 denovo121441 + denovo121478 + denovo121489 + denovo121516 + denovo121519 + denovo121523 + denovo121551 + denovo121560 +
                                 denovo121570 + denovo121640 + denovo121675 + denovo121684 + denovo121721 + denovo121782 + denovo121785 + denovo121800 +
                                 denovo121801 + denovo121827 + denovo121860 + denovo121867 + denovo121914 + denovo121925 + denovo121970 + denovo121990 +
                                 denovo121999 + denovo122006 + denovo122072 + denovo122092 + denovo122096 + denovo122113 + denovo122132 + denovo122165 +
                                 denovo122168 + denovo122178 + denovo122213 + denovo122224 + denovo122259 + denovo122298 + denovo122316 + denovo122344 +
                                 denovo122378 + denovo122380 + denovo122414 + denovo122480 + denovo122490 + denovo122491 + denovo122513 + denovo122518 +
                                 denovo122520 + denovo122548 + denovo122553 + denovo122556 + denovo122565 + denovo122568 + denovo122573 + denovo122583 +
                                 denovo122592 + denovo122615 + denovo122639 + denovo122645 + denovo122666 + denovo122667 + denovo122746 + denovo122780 +
                                 denovo122797 + denovo122798 + denovo122845 + denovo122906 + denovo122921 + denovo122925 + denovo122930 + denovo122954 +
                                 denovo122960 + denovo122977 + denovo123001 + denovo123074 + denovo123077 + denovo123115 + denovo123118 + denovo123119 +
                                 denovo123145 + denovo123157 + denovo123162 + denovo123170 + denovo123214 + denovo123215 + denovo123261 + denovo123303 +
                                 denovo123321 + denovo123360 + denovo123447 + denovo123485 + denovo123509 + denovo123519 + denovo123533 + denovo123583 +
                                 denovo123597 + denovo123620 + denovo123646 + denovo123648 + denovo123671 + denovo123675 + denovo123717 + denovo123776 +
                                 denovo123806 + denovo123818 + denovo123827 + denovo123853 + denovo123874 + denovo123882 + denovo123884 + denovo123895 +
                                 denovo123914 + denovo123917 + denovo123967 + denovo123979 + denovo124066 + denovo124070 + denovo124078 + denovo124083 +
                                 denovo124096 + denovo124105 + denovo124106 + denovo124111 + denovo124135 + denovo124158 + denovo124251 + denovo124272 +
                                 denovo124332 + denovo124333 + denovo124360 + denovo124390 + denovo124413 + denovo124434 + denovo124468 + denovo124493 +
                                 denovo124512 + denovo124551 + denovo124552 + denovo124563 + denovo124573 + denovo124591 + denovo124611 + denovo124623 +
                                 denovo124626 + denovo124632 + denovo124636 + denovo124639 + denovo124646 + denovo124648 + denovo124676 + denovo124682 +
                                 denovo124700 + denovo124701 + denovo124703 + denovo124711 + denovo124741 + denovo124750 + denovo124752 + denovo124772 +
                                 denovo124782 + denovo124817 + denovo124822 + denovo124839 + denovo124845 + denovo124858 + denovo124864 + denovo124877 +
                                 denovo124988 + denovo125002 + denovo125014 + denovo125074 + denovo125080 + denovo125130 + denovo125139 + denovo125140 +
                                 denovo125155 + denovo125170 + denovo125179 + denovo125207 + denovo125221 + denovo125244 + denovo125249 + denovo125251 +
                                 denovo125290 + denovo125335 + denovo125365 + denovo125458 + denovo125499 + denovo125550 + denovo125592 + denovo125594 +
                                 denovo125635 + denovo125689 + denovo125706 + denovo125732 + denovo125737 + denovo125744 + denovo125823 + denovo125854 +
                                 denovo125870 + denovo125883 + denovo125923 + denovo125977 + denovo125998 + denovo126003 + denovo126029 + denovo126072 +
                                 denovo126076 + denovo126091 + denovo126119 + denovo126148 + denovo126152 + denovo126178 + denovo126196 + denovo126200 +
                                 denovo126207 + denovo126209 + denovo126212 + denovo126226 + denovo126280 + denovo126357 + denovo126363 + denovo126378 +
                                 denovo126387 + denovo126396 + denovo126408 + denovo126415 + denovo126468 + denovo126472 + denovo126495 + denovo126510 +
                                 denovo126572 + denovo126600 + denovo126659 + denovo126671 + denovo126677 + denovo126699 + denovo126700 + denovo126726 +
                                 denovo126729 + denovo126741 + denovo126752 + denovo126781 + denovo126799 + denovo126852 + denovo126881 + denovo126891 +
                                 denovo126911 + denovo126980 + denovo127002 + denovo127054 + denovo127098 + denovo127118 + denovo127141 + denovo127143 +
                                 denovo127146 + denovo127164 + denovo127171 + denovo127187 + denovo127189 + denovo127198 + denovo127202 + denovo127219 +
                                 denovo127233 + denovo127250 + denovo127251 + denovo127357 + denovo127415 + denovo127440 + denovo127471 + denovo127477 +
                                 denovo127482 + denovo127499 + denovo127521 + denovo127568 + denovo127573 + denovo127575 + denovo127637 + denovo127663 +
                                 denovo127677 + denovo127692 + denovo127701 + denovo127721 + denovo127733 + denovo127758 + denovo127817 + denovo127851 +
                                 denovo127873 + denovo127900 + denovo127945 + denovo127955 + denovo127957 + denovo127966 + denovo127971 + denovo127972 +
                                 denovo127975 + denovo127980 + denovo127992 + denovo127998 + denovo128053 + denovo128071 + denovo128101 + denovo128107 +
                                 denovo128116 + denovo128125 + denovo128140 + denovo128175 + denovo128178 + denovo128221 + denovo128223 + denovo128250 +
                                 denovo128263 + denovo128293 + denovo128315 + denovo128325 + denovo128354 + denovo128371 + denovo128391 + denovo128398 +
                                 denovo128401 + denovo128462 + denovo128463 + denovo128476 + denovo128515 + denovo128523 + denovo128557 + denovo128589 +
                                 denovo128592 + denovo128622 + denovo128623 + denovo128634 + denovo128644 + denovo128672 + denovo128675 + denovo128701 +
                                 denovo128736 + denovo128747 + denovo128762 + denovo128764 + denovo128783 + denovo128796 + denovo128814 + denovo128816 +
                                 denovo128829 + denovo128838 + denovo128845 + denovo128863 + denovo128893 + denovo128909 + denovo128917 + denovo128920 +
                                 denovo128936 + denovo128964 + denovo128968 + denovo128977 + denovo129007 + denovo129013 + denovo129141 + denovo129152 +
                                 denovo129177 + denovo129196 + denovo129198 + denovo129235 + denovo129238 + denovo129240 + denovo129254 + denovo129258 +
                                 denovo129277 + denovo129283 + denovo129284 + denovo129301 + denovo129310 + denovo129324 + denovo129333 + denovo129380 +
                                 denovo129384 + denovo129394 + denovo129399 + denovo129408 + denovo129428 + denovo129440 + denovo129441 + denovo129514 +
                                 denovo129543 + denovo129626 + denovo129692 + denovo129704 + denovo129782 + denovo129784 + denovo129803 + denovo129844 +
                                 denovo129848 + denovo129886 + denovo129904 + denovo129947 + denovo129952 + denovo129999 + denovo130002 + denovo130019 +
                                 denovo130022 + denovo130116 + denovo130120 + denovo130178 + denovo130234 + denovo130249 + denovo130253 + denovo130254 +
                                 denovo130285 + denovo130302 + denovo130313 + denovo130341 + denovo130370 + denovo130377 + denovo130419 + denovo130437 +
                                 denovo130455 + denovo130482 + denovo130519 + denovo130524 + denovo130551 + denovo130565 + denovo130573 + denovo130614 +
                                 denovo130629 + denovo130658 + denovo130666 + denovo130722 + denovo130730 + denovo130742 + denovo130759 + denovo130852 +
                                 denovo130858 + denovo130944 + denovo130962 + denovo130997 + denovo131009 + denovo131022 + denovo131053 + denovo131078 +
                                 denovo131087 + denovo131138 + denovo131144 + denovo131163 + denovo131166 + denovo131187 + denovo131208 + denovo131261 +
                                 denovo131302 + denovo131303 + denovo131309 + denovo131347 + denovo131348 + denovo131350 + denovo131358 + denovo131359 +
                                 denovo131404 + denovo131407 + denovo131411 + denovo131414 + denovo131426 + denovo131468 + denovo131531 + denovo131565 +
                                 denovo131570 + denovo131571 + denovo131574 + denovo131726 + denovo131779 + denovo131782 + denovo131820 + denovo131825 +
                                 denovo131832 + denovo131839 + denovo131843 + denovo131885 + denovo131932 + denovo131949 + denovo131965 + denovo131971 +
                                 denovo131973 + denovo132001 + denovo132004 + denovo132022 + denovo132024 + denovo132032 + denovo132095 + denovo132137 +
                                 denovo132140 + denovo132164 + denovo132196 + denovo132217 + denovo132263 + denovo132306 + denovo132322 + denovo132421 +
                                 denovo132425 + denovo132468 + denovo132469 + denovo132504 + denovo132511 + denovo132515 + denovo132535 + denovo132597 +
                                 denovo132714 + denovo132812 + denovo132839 + denovo132854 + denovo132897 + denovo132920 + denovo132929 + denovo132933 +
                                 denovo132986 + denovo132990 + denovo132993 + denovo133026 + denovo133027 + denovo133085 + denovo133221 + denovo133247 +
                                 denovo133293 + denovo133295 + denovo133319 + denovo133322 + denovo133331 + denovo133353 + denovo133358 + denovo133403 +
                                 denovo133404 + denovo133446 + denovo133470 + denovo133520 + denovo133532 + denovo133542 + denovo133590 + denovo133599 +
                                 denovo133605 + denovo133632 + denovo133647 + denovo133665 + denovo133741 + denovo133746 + denovo133781 + denovo133792 +
                                 denovo133808 + denovo133858 + denovo133859 + denovo133872 + denovo133873 + denovo133891 + denovo133905 + denovo133916 +
                                 denovo133934 + denovo133939 + denovo133957 + denovo133971 + denovo134018 + denovo134052 + denovo134098 + denovo134109 +
                                 denovo134145 + denovo134274 + denovo134285 + denovo134292 + denovo134324 + denovo134332 + denovo134336 + denovo134409 +
                                 denovo134417 + denovo134436 + denovo134441 + denovo134444 + denovo134449 + denovo134469 + denovo134495 + denovo134517 +
                                 denovo134546 + denovo134579 + denovo134586 + denovo134629 + denovo134668 + denovo134672 + denovo134745 + denovo134782 +
                                 denovo134794 + denovo134813 + denovo134832 + denovo134893 + denovo134902 + denovo134907 + denovo134926 + denovo134972 +
                                 denovo135007 + denovo135013 + denovo135018 + denovo135063 + denovo135084 + denovo135105 + denovo135117 + denovo135133 +
                                 denovo135137 + denovo135158 + denovo135163 + denovo135166 + denovo135171 + denovo135221 + denovo135222 + denovo135226 +
                                 denovo135255 + denovo135261 + denovo135287 + denovo135298 + denovo135316 + denovo135330 + denovo135334 + denovo135335 +
                                 denovo135337 + denovo135388 + denovo135499 + denovo135508 + denovo135569 + denovo135601 + denovo135629 + denovo135693 +
                                 denovo135731 + denovo135790 + denovo135823 + denovo135837 + denovo135853 + denovo135855 + denovo135863 + denovo135884 +
                                 denovo135897 + denovo135943 + denovo135945 + denovo135968 + denovo135972 + denovo136019 + denovo136026 + denovo136067 +
                                 denovo136083 + denovo136087 + denovo136094 + denovo136097 + denovo136111 + denovo136164 + denovo136194 + denovo136196 +
                                 denovo136233 + denovo136255 + denovo136261 + denovo136274 + denovo136280 + denovo136286 + denovo136350 + denovo136375 +
                                 denovo136408 + denovo136446 + denovo136456 + denovo136545 + denovo136569 + denovo136588 + denovo136593 + denovo136598 +
                                 denovo136600 + denovo136623 + denovo136626 + denovo136633 + denovo136649 + denovo136655 + denovo136709 + denovo136826 +
                                 denovo136859 + denovo136861 + denovo136871 + denovo136920 + denovo136946 + denovo136990 + denovo137066 + denovo137149 +
                                 denovo137159 + denovo137221 + denovo137230 + denovo137292 + denovo137305 + denovo137318 + denovo137442 + denovo137473 +
                                 denovo137505 + denovo137514 + denovo137542 + denovo137549 + denovo137552 + denovo137553 + denovo137579 + denovo137590 +
                                 denovo137597 + denovo137627 + denovo137657 + denovo137664 + denovo137700 + denovo137738 + denovo137757 + denovo137765 +
                                 denovo137793 + denovo137802 + denovo137845 + denovo137849 + denovo137881 + denovo137888 + denovo137912 + denovo137991 +
                                 denovo137998 + denovo138054 + denovo138099 + denovo138127 + denovo138173 + denovo138176 + denovo138225 + denovo138228 +
                                 denovo138233 + denovo138248 + denovo138266 + denovo138309 + denovo138313 + denovo138332 + denovo138360 + denovo138374 +
                                 denovo138375 + denovo138383 + denovo138431 + denovo138477 + denovo138486 + denovo138500 + denovo138501 + denovo138503 +
                                 denovo138504 + denovo138518 + denovo138524 + denovo138535 + denovo138557 + denovo138621 + denovo138641 + denovo138656 +
                                 denovo138662 + denovo138689 + denovo138741 + denovo138766 + denovo138774 + denovo138801 + denovo138815 + denovo138822 +
                                 denovo138850 + denovo138914 + denovo138947 + denovo138954 + denovo138957 + denovo138958 + denovo138967 + denovo138991 +
                                 denovo139052 + denovo139084 + denovo139136 + denovo139148 + denovo139156 + denovo139161 + denovo139186 + denovo139219 +
                                 denovo139237 + denovo139269 + denovo139275 + denovo139283 + denovo139344 + denovo139379 + denovo139386 + denovo139437 +
                                 denovo139438 + denovo139454 + denovo139466 + denovo139471 + denovo139477 + denovo139520 + denovo139530 + denovo139565 +
                                 denovo139581 + denovo139612 + denovo139638 + denovo139675 + denovo139683 + denovo139695 + denovo139722 + denovo139724 +
                                 denovo139740 + denovo139759 + denovo139772 + denovo139813 + denovo139852 + denovo139869 + denovo139880 + denovo139925 +
                                 denovo139994 + denovo140005 + denovo140026 + denovo140096 + denovo140101 + denovo140103 + denovo140105 + denovo140116 +
                                 denovo140126 + denovo140169 + denovo140193 + denovo140220 + denovo140261 + denovo140307 + denovo140332 + denovo140360 +
                                 denovo140363 + denovo140372 + denovo140384 + denovo140437 + denovo140455 + denovo140504 + denovo140508 + denovo140538 +
                                 denovo140543 + denovo140584 + denovo140585 + denovo140625 + denovo140626 + denovo140653 + denovo140658 + denovo140685 +
                                 denovo140694 + denovo140749 + denovo140771 + denovo140787 + denovo140796 + denovo140870 + denovo140887 + denovo140895 +
                                 denovo140897 + denovo140970 + denovo140974 + denovo141042 + denovo141062 + denovo141085 + denovo141091 + denovo141149 +
                                 denovo141160 + denovo141196 + denovo141201 + denovo141206 + denovo141231 + denovo141256 + denovo141311 + denovo141331 +
                                 denovo141350 + denovo141361 + denovo141422 + denovo141435 + denovo141440 + denovo141472 + denovo141478 + denovo141506 +
                                 denovo141593 + denovo141648 + denovo141661 + denovo141678 + denovo141728 + denovo141732 + denovo141811 + denovo141821 +
                                 denovo141859 + denovo141870 + denovo141892 + denovo141899 + denovo141938 + denovo141941 + denovo141947 + denovo142011 +
                                 denovo142035 + denovo142082 + denovo142106 + denovo142113 + denovo142166 + denovo142221 + denovo142243 + denovo142299 +
                                 denovo142302 + denovo142350 + denovo142365 + denovo142408 + denovo142446 + denovo142456 + denovo142479 + denovo142489 +
                                 denovo142514 + denovo142519 + denovo142588 + denovo142593 + denovo142635 + denovo142643 + denovo142645 + denovo142653 +
                                 denovo142654 + denovo142689 + denovo142690 + denovo142706 + denovo142712 + denovo142766 + denovo142821 + denovo142854 +
                                 denovo142920 + denovo142943 + denovo142956 + denovo142957 + denovo142964 + denovo142990 + denovo142991 + denovo143008 +
                                 denovo143014 + denovo143021 + denovo143034 + denovo143114 + denovo143135 + denovo143145 + denovo143160 + denovo143243 +
                                 denovo143269 + denovo143281 + denovo143301 + denovo143389 + denovo143420 + denovo143431 + denovo143435 + denovo143436 +
                                 denovo143440 + denovo143535 + denovo143541 + denovo143582 + denovo143593 + denovo143627 + denovo143666 + denovo143694 +
                                 denovo143734 + denovo143738 + denovo143793 + denovo143847 + denovo143887 + denovo143888 + denovo143889 + denovo143902 +
                                 denovo143914 + denovo143919 + denovo143935 + denovo143946 + denovo143947 + denovo144004 + denovo144010 + denovo144038 +
                                 denovo144040 + denovo144046 + denovo144053 + denovo144073 + denovo144077 + denovo144119 + denovo144134 + denovo144191 +
                                 denovo144203 + denovo144206 + denovo144207 + denovo144230 + denovo144259 + denovo144275 + denovo144283 + denovo144294 +
                                 denovo144323 + denovo144348 + denovo144356 + denovo144363 + denovo144442 + denovo144466 + denovo144516 + denovo144519 +
                                 denovo144525 + denovo144562 + denovo144584 + denovo144587 + denovo144604 + denovo144606 + denovo144610 + denovo144620 +
                                 denovo144657 + denovo144671 + denovo144697 + denovo144707 + denovo144712 + denovo144777 + denovo144819 + denovo144888 +
                                 denovo144900 + denovo144918 + denovo144983 + denovo145002 + denovo145075 + denovo145099 + denovo145102 + denovo145105 +
                                 denovo145118 + denovo145126 + denovo145137 + denovo145204 + denovo145281 + denovo145288 + denovo145289 + denovo145303 +
                                 denovo145319 + denovo145334 + denovo145392 + denovo145401 + denovo145413 + denovo145432 + denovo145441 + denovo145485 +
                                 denovo145492 + denovo145517 + denovo145520 + denovo145543 + denovo145562 + denovo145602 + denovo145623 + denovo145642 +
                                 denovo145673 + denovo145716 + denovo145735 + denovo145740 + denovo145859 + denovo145888 + denovo145939 + denovo145949 +
                                 denovo145999 + denovo146060 + denovo146078 + denovo146109 + denovo146129 + denovo146140 + denovo146148 + denovo146159 +
                                 denovo146186 + denovo146190 + denovo146195 + denovo146212 + denovo146255 + denovo146276 + denovo146294 + denovo146306 +
                                 denovo146310 + denovo146361 + denovo146392 + denovo146399 + denovo146414 + denovo146424 + denovo146436 + denovo146479 +
                                 denovo146495 + denovo146516 + denovo146543 + denovo146596 + denovo146645 + denovo146669 + denovo146697 + denovo146701 +
                                 denovo146772 + denovo146825 + denovo146840 + denovo146891 + denovo146942 + denovo146947 + denovo146954 + denovo146971 +
                                 denovo146982 + denovo146984 + denovo147036 + denovo147064 + denovo147069 + denovo147109 + denovo147133 + denovo147139 +
                                 denovo147158 + denovo147200 + denovo147243 + denovo147247 + denovo147252 + denovo147268 + denovo147281 + denovo147290 +
                                 denovo147291 + denovo147305 + denovo147314 + denovo147345 + denovo147349 + denovo147364 + denovo147375 + denovo147410 +
                                 denovo147424 + denovo147435 + denovo147496 + denovo147524 + denovo147559 + denovo147620 + denovo147627 + denovo147644 +
                                 denovo147669 + denovo147740 + denovo147741 + denovo147755 + denovo147767 + denovo147802 + denovo147804 + denovo147805 +
                                 denovo147816 + denovo147851 + denovo147944 + denovo147949 + denovo147961 + denovo147973 + denovo147994 + denovo147999 +
                                 denovo148044 + denovo148048 + denovo148074 + denovo148076 + denovo148085 + denovo148089 + denovo148124 + denovo148134 +
                                 denovo148171 + denovo148245 + denovo148259 + denovo148270 + denovo148296 + denovo148354 + denovo148371 + denovo148374 +
                                 denovo148376 + denovo148381 + denovo148456 + denovo148459 + denovo148460 + denovo148463 + denovo148468 + denovo148508 +
                                 denovo148516 + denovo148541 + denovo148559 + denovo148616 + denovo148669 + denovo148691 + denovo148695 + denovo148790 +
                                 denovo148811 + denovo148821 + denovo148832 + denovo148838 + denovo148868 + denovo148871 + denovo148987 + denovo149045 +
                                 denovo149050 + denovo149095 + denovo149117 + denovo149163 + denovo149196 + denovo149230 + denovo149241 + denovo149328 +
                                 denovo149339 + denovo149382 + denovo149401 + denovo149405 + denovo149428 + denovo149436 + denovo149440 + denovo149449 +
                                 denovo149458 + denovo149462 + denovo149499 + denovo149528 + denovo149530 + denovo149547 + denovo149589 + denovo149594 +
                                 denovo149612 + denovo149664 + denovo149684 + denovo149731 + denovo149760 + denovo149761 + denovo149796 + denovo149855 +
                                 denovo149875 + denovo149895 + denovo149913 + denovo149917 + denovo149962 + denovo150001 + denovo150009 + denovo150020 +
                                 denovo150027 + denovo150028 + denovo150065 + denovo150067 + denovo150084 + denovo150085 + denovo150088 + denovo150090 +
                                 denovo150094 + denovo150111 + denovo150123 + denovo150129 + denovo150132 + denovo150180 + denovo150199 + denovo150242 +
                                 denovo150257 + denovo150291 + denovo150300 + denovo150390 + denovo150408 + denovo150411 + denovo150470 + denovo150512 +
                                 denovo150541 + denovo150565 + denovo150602 + denovo150614 + denovo150619 + denovo150625 + denovo150627 + denovo150632 +
                                 denovo150679 + denovo150732 + denovo150788 + denovo150813 + denovo150853 + denovo150855 + denovo150868 + denovo150876 +
                                 denovo150910 + denovo150940 + denovo150950 + denovo150963 + denovo150967 + denovo150982 + denovo151017 + denovo151018 +
                                 denovo151030 + denovo151050 + denovo151055 + denovo151142 + denovo151150 + denovo151185 + denovo151215 + denovo151224 +
                                 denovo151242 + denovo151258 + denovo151267 + denovo151322 + denovo151323 + denovo151331 + denovo151402 + denovo151414 +
                                 denovo151433 + denovo151525 + denovo151568 + denovo151573 + denovo151582 + denovo151642 + denovo151685 + denovo151698 +
                                 denovo151733 + denovo151743 + denovo151784 + denovo151877 + denovo151886 + denovo151898 + denovo151922 + denovo151924 +
                                 denovo151944 + denovo151965 + denovo151971 + denovo152049 + denovo152075 + denovo152084 + denovo152115 + denovo152130 +
                                 denovo152220 + denovo152230 + denovo152283 + denovo152306 + denovo152322 + denovo152331 + denovo152345 + denovo152385 +
                                 denovo152449 + denovo152453 + denovo152483 + denovo152512 + denovo152516 + denovo152519 + denovo152539 + denovo152545 +
                                 denovo152575 + denovo152620 + denovo152648 + denovo152680 + denovo152688 + denovo152750 + denovo152766 + denovo152781 +
                                 denovo152786 + denovo152857 + denovo152863 + denovo152888 + denovo152902 + denovo152928 + denovo152936 + denovo152953 +
                                 denovo152963 + denovo152968 + denovo152979 + denovo153009 + denovo153025 + denovo153031 + denovo153033 + denovo153075 +
                                 denovo153121 + denovo153138 + denovo153164 + denovo153207 + denovo153219 + denovo153242 + denovo153262 + denovo153332 +
                                 denovo153471 + denovo153472 + denovo153494 + denovo153526 + denovo153552 + denovo153580 + denovo153591 + denovo153601 +
                                 denovo153644 + denovo153662 + denovo153746 + denovo153756 + denovo153762 + denovo153764 + denovo153839 + denovo153868 +
                                 denovo153886 + denovo153916 + denovo153937 + denovo153999 + denovo154074 + denovo154104 + denovo154118 + denovo154144 +
                                 denovo154161 + denovo154203 + denovo154209 + denovo154234 + denovo154242 + denovo154283 + denovo154305 + denovo154315 +
                                 denovo154332 + denovo154347 + denovo154351 + denovo154362 + denovo154373 + denovo154388 + denovo154413 + denovo154504 +
                                 denovo154506 + denovo154526 + denovo154541 + denovo154548 + denovo154586 + denovo154616 + denovo154671 + denovo154686 +
                                 denovo154705 + denovo154759 + denovo154760 + denovo154875 + denovo154876 + denovo154895 + denovo154917 + denovo154931 +
                                 denovo154935 + denovo154983 + denovo154988 + denovo155045 + denovo155104 + denovo155125 + denovo155137 + denovo155162 +
                                 denovo155194 + denovo155196 + denovo155203 + denovo155244 + denovo155264 + denovo155274 + denovo155334 + denovo155400 +
                                 denovo155469 + denovo155495 + denovo155501 + denovo155511 + denovo155527 + denovo155529 + denovo155563 + denovo155571 +
                                 denovo155573 + denovo155579 + denovo155596 + denovo155603 + denovo155616 + denovo155618 + denovo155636 + denovo155654 +
                                 denovo155655 + denovo155699 + denovo155710 + denovo155743 + denovo155745 + denovo155842 + denovo155853 + denovo155870 +
                                 denovo155877 + denovo155929 + denovo155950 + denovo155962 + denovo155973 + denovo155982 + denovo155983 + denovo156065 +
                                 denovo156084 + denovo156087 + denovo156104 + denovo156108 + denovo156109 + denovo156113 + denovo156149 + denovo156156 +
                                 denovo156174 + denovo156186 + denovo156218 + denovo156219 + denovo156228 + denovo156281 + denovo156309 + denovo156313 +
                                 denovo156366 + denovo156377 + denovo156389 + denovo156422 + denovo156423 + denovo156500 + denovo156555 + denovo156557 +
                                 denovo156606 + denovo156742 + denovo156755 + denovo156756 + denovo156784 + denovo156814 + denovo156843 + denovo156861 +
                                 denovo156912 + denovo156931 + denovo156944 + denovo156949 + denovo156972 + denovo156975 + denovo157048 + denovo157085 +
                                 denovo157110 + denovo157135 + denovo157176 + denovo157232 + denovo157245 + denovo157289 + denovo157329 + denovo157347 +
                                 denovo157379 + denovo157433 + denovo157627 + denovo157657 + denovo157677 + denovo157690 + denovo157694 + denovo157749 +
                                 denovo157759 + denovo157764 + denovo157770 + denovo157777 + denovo157811 + denovo157827 + denovo157834 + denovo157849 +
                                 denovo157916 + denovo157945 + denovo157946 + denovo157995 + denovo158048 + denovo158080 + denovo158192 + denovo158202 +
                                 denovo158210 + denovo158241 + denovo158282 + denovo158349 + denovo158378 + denovo158400 + denovo158405 + denovo158407 +
                                 denovo158427 + denovo158483 + denovo158499 + denovo158505 + denovo158512 + denovo158513 + denovo158521 + denovo158528 +
                                 denovo158548 + denovo158571 + denovo158579 + denovo158615 + denovo158620 + denovo158626 + denovo158632 + denovo158653 +
                                 denovo158658 + denovo158659 + denovo158691 + denovo158692 + denovo158724 + denovo158725 + denovo158735 + denovo158757 +
                                 denovo158795 + denovo158825 + denovo158833 + denovo158907 + denovo158908 + denovo158913 + denovo159019 + denovo159033 +
                                 denovo159045 + denovo159057 + denovo159102 + denovo159109 + denovo159167 + denovo159171 + denovo159212 + denovo159258 +
                                 denovo159274 + denovo159327 + denovo159408 + denovo159420 + denovo159444 + denovo159450 + denovo159468 + denovo159503 +
                                 denovo159507 + denovo159573 + denovo159619 + denovo159622 + denovo159629 + denovo159639 + denovo159642 + denovo159686 +
                                 denovo159695 + denovo159718 + denovo159727 + denovo159745 + denovo159765 + denovo159849 + denovo159913 + denovo159918 +
                                 denovo159927 + denovo159947 + denovo159969 + denovo160035 + denovo160046 + denovo160060 + denovo160072 + denovo160128 +
                                 denovo160177 + denovo160247 + denovo160352 + denovo160386 + denovo160418 + denovo160443 + denovo160493 + denovo160514 +
                                 denovo160530 + denovo160619 + denovo160734 + denovo160747 + denovo160765 + denovo160773 + denovo160792 + denovo160797 +
                                 denovo160805 + denovo160833 + denovo160863 + denovo160887 + denovo160905 + denovo160995 + denovo161000 + denovo161022 +
                                 denovo161049 + denovo161058 + denovo161130 + denovo161185 + denovo161222 + denovo161229 + denovo161258 + denovo161264 +
                                 denovo161294 + denovo161305 + denovo161349 + denovo161365 + denovo161369 + denovo161405 + denovo161446 + denovo161454 +
                                 denovo161503 + denovo161540 + denovo161594 + denovo161655 + denovo161676 + denovo161680 + denovo161694 + denovo161702 +
                                 denovo161707 + denovo161713 + denovo161734 + denovo161745 + denovo161789 + denovo161806 + denovo161809 + denovo161810 +
                                 denovo161870 + denovo161916 + denovo161942 + denovo161952 + denovo161968 + denovo162000 + denovo162014 + denovo162019 +
                                 denovo162065 + denovo162067 + denovo162098 + denovo162101 + denovo162151 + denovo162155 + denovo162174 + denovo162176 +
                                 denovo162182 + denovo162215 + denovo162219 + denovo162253 + denovo162293 + denovo162340 + denovo162361 + denovo162390 +
                                 denovo162499 + denovo162506 + denovo162520 + denovo162524 + denovo162564 + denovo162600 + denovo162645 + denovo162697 +
                                 denovo162738 + denovo162743 + denovo162754 + denovo162755 + denovo162788 + denovo162799 + denovo162813 + denovo162893 +
                                 denovo162907 + denovo162922 + denovo162940 + denovo162974 + denovo163045 + denovo163052 + denovo163059 + denovo163072 +
                                 denovo163087 + denovo163128 + denovo163131 + denovo163159 + denovo163170 + denovo163211 + denovo163260 + denovo163296 +
                                 denovo163351 + denovo163409 + denovo163484 + denovo163499 + denovo163507 + denovo163530 + denovo163573 + denovo163575 +
                                 denovo163581 + denovo163662 + denovo163679 + denovo163688 + denovo163698 + denovo163725 + denovo163729 + denovo163810 +
                                 denovo163823 + denovo163857 + denovo163873 + denovo163927 + denovo163996 + denovo164006 + denovo164014 + denovo164117 +
                                 denovo164129 + denovo164130 + denovo164178 + denovo164212 + denovo164241 + denovo164263 + denovo164296 + denovo164298 +
                                 denovo164332 + denovo164377 + denovo164450 + denovo164464 + denovo164520 + denovo164608 + denovo164618 + denovo164625 +
                                 denovo164647 + denovo164656 + denovo164658 + denovo164669 + denovo164764 + denovo164770 + denovo164782 + denovo164800 +
                                 denovo164808 + denovo164825 + denovo164851 + denovo164855 + denovo164883 + denovo164902 + denovo164936 + denovo164939 +
                                 denovo164959 + denovo164966 + denovo164975 + denovo164981 + denovo165053 + denovo165084 + denovo165093 + denovo165099 +
                                 denovo165110 + denovo165123 + denovo165139 + denovo165197 + denovo165214 + denovo165231 + denovo165340 + denovo165343 +
                                 denovo165348 + denovo165421 + denovo165465 + denovo165482 + denovo165506 + denovo165508 + denovo165548 + denovo165564 +
                                 denovo165568 + denovo165615 + denovo165621 + denovo165638 + denovo165662 + denovo165730 + denovo165751 + denovo165841 +
                                 denovo165845 + denovo165857 + denovo165860 + denovo165893 + denovo165899 + denovo165924 + denovo165927 + denovo165950 +
                                 denovo165951 + denovo165988 + denovo166030 + denovo166148 + denovo166168 + denovo166172 + denovo166188 + denovo166192 +
                                 denovo166243 + denovo166259 + denovo166310 + denovo166336 + denovo166363 + denovo166394 + denovo166400 + denovo166443 +
                                 denovo166561 + denovo166578 + denovo166592 + denovo166621 + denovo166662 + denovo166680 + denovo166730 + denovo166747 +
                                 denovo166783 + denovo166818 + denovo166848 + denovo166865 + denovo166882 + denovo166910 + denovo166936 + denovo166947 +
                                 denovo166994 + denovo166995 + denovo167037 + denovo167080 + denovo167136 + denovo167187 + denovo167228 + denovo167246 +
                                 denovo167254 + denovo167319 + denovo167326 + denovo167377 + denovo167399 + denovo167402 + denovo167425 + denovo167432 +
                                 denovo167435 + denovo167492 + denovo167505 + denovo167517 + denovo167527 + denovo167535 + denovo167551 + denovo167577 +
                                 denovo167585 + denovo167617 + denovo167665 + denovo167666 + denovo167718 + denovo167735 + denovo167738 + denovo167744 +
                                 denovo167751 + denovo167774 + denovo167788 + denovo167790 + denovo167837 + denovo167849 + denovo167863 + denovo167924 +
                                 denovo167925 + denovo167936 + denovo167973 + denovo167979 + denovo167993 + denovo168068 + denovo168118 + denovo168150 +
                                 denovo168166 + denovo168201 + denovo168204 + denovo168218 + denovo168223 + denovo168225 + denovo168234 + denovo168250 +
                                 denovo168270 + denovo168294 + denovo168316 + denovo168393 + denovo168430 + denovo168464 + denovo168474 + denovo168497 +
                                 denovo168538 + denovo168552 + denovo168565 + denovo168576 + denovo168616 + denovo168642 + denovo168652 + denovo168674 +
                                 denovo168704 + denovo168725 + denovo168767 + denovo168775 + denovo168800 + denovo168847 + denovo168851 + denovo168912 +
                                 denovo168919 + denovo168986 + denovo168994 + denovo169055 + denovo169118 + denovo169121 + denovo169122 + denovo169146 +
                                 denovo169147 + denovo169222 + denovo169250 + denovo169311 + denovo169315 + denovo169372 + denovo169411 + denovo169417 +
                                 denovo169445 + denovo169448 + denovo169456 + denovo169464 + denovo169465 + denovo169473 + denovo169510 + denovo169689 +
                                 denovo169694 + denovo169695 + denovo169707 + denovo169722 + denovo169732 + denovo169737 + denovo169756 + denovo169840 +
                                 denovo169850 + denovo169853 + denovo169866 + denovo169900 + denovo169922 + denovo169931 + denovo169953 + denovo169954 +
                                 denovo169956 + denovo169963 + denovo170099 + denovo170112 + denovo170180 + denovo170221 + denovo170253 + denovo170259 +
                                 denovo170334 + denovo170353 + denovo170385 + denovo170451 + denovo170455 + denovo170474 + denovo170533 + denovo170545 +
                                 denovo170589 + denovo170602 + denovo170652 + denovo170689 + denovo170775 + denovo170788 + denovo170799 + denovo170879 +
                                 denovo170932 + denovo170954 + denovo171006 + denovo171034 + denovo171040 + denovo171059 + denovo171063 + denovo171127 +
                                 denovo171138 + denovo171141 + denovo171143 + denovo171162 + denovo171194 + denovo171222 + denovo171231 + denovo171232 +
                                 denovo171239 + denovo171250 + denovo171263 + denovo171287 + denovo171295 + denovo171321 + denovo171368 + denovo171372 +
                                 denovo171384 + denovo171449 + denovo171455 + denovo171459 + denovo171503 + denovo171567 + denovo171599 + denovo171613 +
                                 denovo171639 + denovo171656 + denovo171720 + denovo171721 + denovo171733 + denovo171737 + denovo171746 + denovo171786 +
                                 denovo171824 + denovo171852 + denovo171855 + denovo171881 + denovo171914 + denovo171925 + denovo171998 + denovo172032 +
                                 denovo172173 + denovo172203 + denovo172255 + denovo172272 + denovo172274 + denovo172301 + denovo172336 + denovo172340 +
                                 denovo172422 + denovo172477 + denovo172479 + denovo172497 + denovo172541 + denovo172556 + denovo172581 + denovo172589 +
                                 denovo172607 + denovo172617 + denovo172630 + denovo172691 + denovo172744 + denovo172785 + denovo172796 + denovo172839 +
                                 denovo172862 + denovo172941 + denovo172950 + denovo173049 + denovo173136 + denovo173178 + denovo173199 + denovo173240 +
                                 denovo173242 + denovo173272 + denovo173299 + denovo173328 + denovo173337 + denovo173339 + denovo173367 + denovo173397 +
                                 denovo173457 + denovo173475 + denovo173477 + denovo173491 + denovo173496 + denovo173508 + denovo173575 + denovo173591 +
                                 denovo173609 + denovo173634 + denovo173638 + denovo173691 + denovo173715 + denovo173738 + denovo173766 + denovo173772 +
                                 denovo173775 + denovo173828 + denovo174012 + denovo174014 + denovo174103 + denovo174123 + denovo174153 + denovo174157 +
                                 denovo174162 + denovo174198 + denovo174201 + denovo174212 + denovo174216 + denovo174217 + denovo174266 + denovo174271 +
                                 denovo174278 + denovo174293 + denovo174303 + denovo174385 + denovo174386 + denovo174430 + denovo174435 + denovo174459 +
                                 denovo174464 + denovo174470 + denovo174497 + denovo174518 + denovo174570 + denovo174613 + denovo174618 + denovo174655 +
                                 denovo174676 + denovo174691 + denovo174703 + denovo174765 + denovo174773 + denovo174780 + denovo174827 + denovo174846 +
                                 denovo174852 + denovo174857 + denovo174868 + denovo174886 + denovo174941 + denovo174974 + denovo175001 + denovo175012 +
                                 denovo175060 + denovo175095 + denovo175122 + denovo175143 + denovo175198 + denovo175200 + denovo175207 + denovo175208 +
                                 denovo175228 + denovo175317 + denovo175326 + denovo175332 + denovo175360 + denovo175367 + denovo175420 + denovo175461 +
                                 denovo175494 + denovo175503 + denovo175508 + denovo175514 + denovo175538 + denovo175543 + denovo175558 + denovo175590 +
                                 denovo175645 + denovo175709 + denovo175713 + denovo175732 + denovo175733 + denovo175763 + denovo175777 + denovo175778 +
                                 denovo175939 + denovo175972 + denovo175974 + denovo175978 + denovo175983 + denovo175995 + denovo175996 + denovo176025 +
                                 denovo176046 + denovo176094 + denovo176100 + denovo176172 + denovo176206 + denovo176235 + denovo176238 + denovo176255 +
                                 denovo176337 + denovo176364 + denovo176412 + denovo176416 + denovo176418 + denovo176463 + denovo176467 + denovo176471 +
                                 denovo176473 + denovo176503 + denovo176547 + denovo176561 + denovo176617 + denovo176634 + denovo176682 + denovo176685 +
                                 denovo176718 + denovo176757 + denovo176760 + denovo176765 + denovo176778 + denovo176800 + denovo176852 + denovo176861 +
                                 denovo176903 + denovo176966 + denovo176981 + denovo176988 + denovo176995 + denovo177020 + denovo177055 + denovo177087 +
                                 denovo177096 + denovo177145 + denovo177198 + denovo177236 + denovo177260 + denovo177271 + denovo177279 + denovo177316 +
                                 denovo177339 + denovo177354 + denovo177366 + denovo177393 + denovo177406 + denovo177411 + denovo177500 + denovo177520 +
                                 denovo177557 + denovo177559 + denovo177590 + denovo177599 + denovo177689 + denovo177707 + denovo177723 + denovo177769 +
                                 denovo177814 + denovo177830 + denovo177839 + denovo177843 + denovo177863 + denovo177898 + denovo177917 + denovo177929 +
                                 denovo177947 + denovo177970 + denovo178002 + denovo178062 + denovo178064 + denovo178070 + denovo178109 + denovo178137 +
                                 denovo178152 + denovo178161 + denovo178187 + denovo178202 + denovo178205 + denovo178212 + denovo178214 + denovo178295 +
                                 denovo178298 + denovo178300 + denovo178314 + denovo178326 + denovo178359 + denovo178397 + denovo178409 + denovo178423 +
                                 denovo178425 + denovo178501 + denovo178510 + denovo178514 + denovo178545 + denovo178546 + denovo178625 + denovo178646 +
                                 denovo178672 + denovo178716 + denovo178735 + denovo178829 + denovo178831 + denovo178903 + denovo178927 + denovo178990 +
                                 denovo179025 + denovo179027 + denovo179048 + denovo179086 + denovo179144 + denovo179159 + denovo179194 + denovo179210 +
                                 denovo179228 + denovo179243 + denovo179317 + denovo179330 + denovo179334 + denovo179341 + denovo179362 + denovo179384 +
                                 denovo179413 + denovo179416 + denovo179443 + denovo179481 + denovo179493 + denovo179499 + denovo179502 + denovo179513 +
                                 denovo179522 + denovo179523 + denovo179535 + denovo179538 + denovo179578 + denovo179583 + denovo179591 + denovo179593 +
                                 denovo179594 + denovo179613 + denovo179645 + denovo179693 + denovo179766 + denovo179785 + denovo179795 + denovo179798 +
                                 denovo179825 + denovo179835 + denovo179846 + denovo179873 + denovo179874 + denovo179907 + denovo179911 + denovo179915 +
                                 denovo179922 + denovo179925 + denovo179936 + denovo179999 + denovo180017 + denovo180072 + denovo180090 + denovo180154 +
                                 denovo180163 + denovo180164 + denovo180171 + denovo180187 + denovo180191 + denovo180202 + denovo180218 + denovo180228 +
                                 denovo180244 + denovo180257 + denovo180276 + denovo180369 + denovo180385 + denovo180391 + denovo180397 + denovo180411 +
                                 denovo180417 + denovo180421 + denovo180444 + denovo180475 + denovo180487 + denovo180532 + denovo180618 + denovo180672 +
                                 denovo180767 + denovo180768 + denovo180788 + denovo180822 + denovo180853 + denovo180872 + denovo180876 + denovo180911 +
                                 denovo180934 + denovo181058 + denovo181102 + denovo181118 + denovo181155 + denovo181192 + denovo181197 + denovo181206 +
                                 denovo181235 + denovo181244 + denovo181269 + denovo181288 + denovo181305 + denovo181340 + denovo181392 + denovo181455 +
                                 denovo181502 + denovo181530 + denovo181536 + denovo181545 + denovo181606 + denovo181636 + denovo181773 + denovo181783 +
                                 denovo181809 + denovo181834 + denovo181835 + denovo181874 + denovo181883 + denovo181892 + denovo181903 + denovo181909 +
                                 denovo181928 + denovo181998 + denovo182030 + denovo182052 + denovo182057 + denovo182081 + denovo182086 + denovo182122 +
                                 denovo182172 + denovo182236 + denovo182267 + denovo182297 + denovo182322 + denovo182328 + denovo182351 + denovo182352 +
                                 denovo182356 + denovo182362 + denovo182366 + denovo182390 + denovo182403 + denovo182451 + denovo182496 + denovo182558 +
                                 denovo182593 + denovo182596 + denovo182600 + denovo182638 + denovo182645 + denovo182650 + denovo182743 + denovo182748 +
                                 denovo182756 + denovo182797 + denovo182820 + denovo182887 + denovo182895 + denovo182898 + denovo182906 + denovo182908 +
                                 denovo182909 + denovo182927 + denovo182959 + denovo182964 + denovo182971 + denovo182997 + denovo183026 + denovo183067 +
                                 denovo183072 + denovo183090 + denovo183095 + denovo183122 + denovo183159 + denovo183183 + denovo183190 + denovo183222 +
                                 denovo183226 + denovo183242 + denovo183248 + denovo183249 + denovo183261 + denovo183270 + denovo183305 + denovo183318 +
                                 denovo183323 + denovo183335 + denovo183355 + denovo183365 + denovo183385 + denovo183412 + denovo183440 + denovo183475 +
                                 denovo183490 + denovo183516 + denovo183562 + denovo183627 + denovo183650 + denovo183654 + denovo183699 + denovo183711 +
                                 denovo183723 + denovo183734 + denovo183736 + denovo183748 + denovo183749 + denovo183822 + denovo183907 + denovo183916 +
                                 denovo183969 + denovo183975 + denovo184010 + denovo184034 + denovo184041 + denovo184053 + denovo184063 + denovo184123 +
                                 denovo184124 + denovo184126 + denovo184238 + denovo184259 + denovo184261 + denovo184279 + denovo184296 + denovo184320 +
                                 denovo184371 + denovo184397 + denovo184418 + denovo184460 + denovo184553 + denovo184560 + denovo184561 + denovo184568 +
                                 denovo184600 + denovo184607 + denovo184620 + denovo184628 + denovo184629 + denovo184637 + denovo184645 + denovo184674 +
                                 denovo184692 + denovo184700 + denovo184714 + denovo184733 + denovo184739 + denovo184751 + denovo184757 + denovo184789 +
                                 denovo184819 + denovo184827 + denovo184841 + denovo184861 + denovo184862 + denovo184865 + denovo184880 + denovo184924 +
                                 denovo184930 + denovo184999 + denovo185008 + denovo185025 + denovo185048 + denovo185056 + denovo185061 + denovo185071 +
                                 denovo185112 + denovo185146 + denovo185164 + denovo185169 + denovo185217 + denovo185222 + denovo185240 + denovo185280 +
                                 denovo185331 + denovo185388 + denovo185422 + denovo185429 + denovo185432 + denovo185457 + denovo185473 + denovo185544 +
                                 denovo185545 + denovo185583 + denovo185608 + denovo185620 + denovo185635 + denovo185668 + denovo185808 + denovo185844 +
                                 denovo185845 + denovo185930 + denovo185933 + denovo185938 + denovo185946 + denovo185965 + denovo185979 + denovo185980 +
                                 denovo186021 + denovo186022 + denovo186024 + denovo186051 + denovo186085 + denovo186110 + denovo186125 + denovo186145 +
                                 denovo186223 + denovo186225 + denovo186260 + denovo186285 + denovo186330 + denovo186353 + denovo186376 + denovo186388 +
                                 denovo186453 + denovo186462 + denovo186499 + denovo186502 + denovo186519 + denovo186534 + denovo186576 + denovo186652 +
                                 denovo186669 + denovo186689 + denovo186740 + denovo186756 + denovo186799 + denovo186815 + denovo186818 + denovo186855 +
                                 denovo186895 + denovo186916 + denovo186938 + denovo186974 + denovo186975 + denovo186994 + denovo186995 + denovo187013 +
                                 denovo187033 + denovo187068 + denovo187095 + denovo187154 + denovo187172 + denovo187252 + denovo187262 + denovo187278 +
                                 denovo187282 + denovo187365 + denovo187427 + denovo187440 + denovo187452 + denovo187490 + denovo187527 + denovo187528 +
                                 denovo187554 + denovo187603 + denovo187626 + denovo187685 + denovo187785 + denovo187789 + denovo187812 + denovo187914 +
                                 denovo187919 + denovo188002 + denovo188024 + denovo188060 + denovo188067 + denovo188091 + denovo188097 + denovo188101 +
                                 denovo188137 + denovo188159 + denovo188183 + denovo188234 + denovo188262 + denovo188300 + denovo188328 + denovo188376 +
                                 denovo188438 + denovo188471 + denovo188507 + denovo188526 + denovo188560 + denovo188576 + denovo188578 + denovo188587 +
                                 denovo188616 + denovo188639 + denovo188791 + denovo188840 + denovo188843 + denovo188854 + denovo188858 + denovo188895 +
                                 denovo188957 + denovo189047 + denovo189056 + denovo189091 + denovo189111 + denovo189154 + denovo189158 + denovo189190 +
                                 denovo189192 + denovo189207 + denovo189218 + denovo189230 + denovo189272 + denovo189284 + denovo189354 + denovo189368 +
                                 denovo189388 + denovo189398 + denovo189430 + denovo189465 + denovo189475 + denovo189529 + denovo189594 + denovo189629 +
                                 denovo189654 + denovo189660 + denovo189665 + denovo189666 + denovo189671 + denovo189695 + denovo189762 + denovo189765 +
                                 denovo189833 + denovo189853 + denovo189854 + denovo189864 + denovo189893 + denovo189905 + denovo189959 + denovo189981 +
                                 denovo189993 + denovo190018 + denovo190047 + denovo190068 + denovo190071 + denovo190106 + denovo190157 + denovo190183 +
                                 denovo190205 + denovo190213 + denovo190269 + denovo190322 + denovo190324 + denovo190364 + denovo190383 + denovo190400 +
                                 denovo190403 + denovo190409 + denovo190438 + denovo190461 + denovo190513 + denovo190568 + denovo190589 + denovo190599 +
                                 denovo190611 + denovo190614 + denovo190651 + denovo190668 + denovo190669 + denovo190699 + denovo190701 + denovo190714 +
                                 denovo190776 + denovo190785 + denovo190810 + denovo190831 + denovo190876 + denovo190883 + denovo190886 + denovo190897 +
                                 denovo190898 + denovo190943 + denovo190963 + denovo191057 + denovo191061 + denovo191084 + denovo191087 + denovo191146 +
                                 denovo191195 + denovo191196 + denovo191231 + denovo191237 + denovo191247 + denovo191282 + denovo191314 + denovo191336 +
                                 denovo191369 + denovo191388 + denovo191407 + denovo191424 + denovo191474 + denovo191483 + denovo191537 + denovo191539 +
                                 denovo191547 + denovo191603 + denovo191659 + denovo191727 + denovo191731 + denovo191771 + denovo191825 + denovo191831 +
                                 denovo191848 + denovo191852 + denovo191874 + denovo191925 + denovo191948 + denovo191982 + denovo191995 + denovo192007 +
                                 denovo192020 + denovo192191 + denovo192231 + denovo192268 + denovo192304 + denovo192319 + denovo192326 + denovo192342 +
                                 denovo192375 + denovo192427 + denovo192443 + denovo192452 + denovo192460 + denovo192482 + denovo192483 + denovo192580 +
                                 denovo192593 + denovo192615 + denovo192620 + denovo192637 + denovo192672 + denovo192718 + denovo192733 + denovo192739 +
                                 denovo192827 + denovo192858 + denovo192884 + denovo192965 + denovo192972 + denovo192976 + denovo193072 + denovo193121 +
                                 denovo193191 + denovo193197 + denovo193205 + denovo193220 + denovo193269 + denovo193308 + denovo193319 + denovo193328 +
                                 denovo193333 + denovo193371 + denovo193392 + denovo193467 + denovo193471 + denovo193512 + denovo193542 + denovo193544 +
                                 denovo193581 + denovo193620 + denovo193627 + denovo193642 + denovo193664 + denovo193714 + denovo193723 + denovo193739 +
                                 denovo193740 + denovo193744 + denovo193791 + denovo193807 + denovo193817 + denovo193823 + denovo193853 + denovo193863 +
                                 denovo193918 + denovo193951 + denovo193964 + denovo193971 + denovo193990 + denovo193999 + denovo194010 + denovo194014 +
                                 denovo194024 + denovo194035 + denovo194038 + denovo194046 + denovo194058 + denovo194077 + denovo194079 + denovo194096 +
                                 denovo194133 + denovo194171 + denovo194231 + denovo194298 + denovo194326 + denovo194455 + denovo194524 + denovo194537 +
                                 denovo194549 + denovo194594 + denovo194628 + denovo194688 + denovo194709 + denovo194713 + denovo194724 + denovo194725 +
                                 denovo194764 + denovo194768 + denovo194769 + denovo194774 + denovo194776 + denovo194777 + denovo194808 + denovo194826 +
                                 denovo194888 + denovo194949 + denovo194970 + denovo195020 + denovo195087 + denovo195091 + denovo195103 + denovo195123 +
                                 denovo195174 + denovo195207 + denovo195222 + denovo195224 + denovo195251 + denovo195270 + denovo195272 + denovo195284 +
                                 denovo195385 + denovo195457 + denovo195470 + denovo195502 + denovo195526 + denovo195527 + denovo195562 + denovo195573 +
                                 denovo195632 + denovo195634 + denovo195639 + denovo195698 + denovo195700 + denovo195705 + denovo195731 + denovo195778 +
                                 denovo195791 + denovo195804 + denovo195819 + denovo195858 + denovo195871 + denovo195875 + denovo195885 + denovo195888 +
                                 denovo195918 + denovo195919 + denovo195944 + denovo195973 + denovo196007 + denovo196076 + denovo196089 + denovo196124 +
                                 denovo196129 + denovo196153 + denovo196161 + denovo196168 + denovo196207 + denovo196208 + denovo196237 + denovo196244 +
                                 denovo196247 + denovo196253 + denovo196265 + denovo196270 + denovo196272 + denovo196365 + denovo196372 + denovo196387 +
                                 denovo196411 + denovo196466 + denovo196500 + denovo196575 + denovo196677 + denovo196693 + denovo196695 + denovo196696 +
                                 denovo196736 + denovo196737 + denovo196761 + denovo196791 + denovo196813 + denovo196824 + denovo196857 + denovo196863 +
                                 denovo196865 + denovo196899 + denovo196933 + denovo196950 + denovo196989 + denovo196998 + denovo196999 + denovo197012 +
                                 denovo197032 + denovo197033 + denovo197039 + denovo197074 + denovo197078 + denovo197101 + denovo197132 + denovo197136 +
                                 denovo197137 + denovo197189 + denovo197196 + denovo197199 + denovo197220 + denovo197241 + denovo197248 + denovo197250 +
                                 denovo197255 + denovo197276 + denovo197371 + denovo197383 + denovo197400 + denovo197421 + denovo197429 + denovo197460 +
                                 denovo197461 + denovo197524 + denovo197578 + denovo197599 + denovo197600 + denovo197612 + denovo197634 + denovo197649 +
                                 denovo197650 + denovo197681 + denovo197682 + denovo197715 + denovo197736 + denovo197738 + denovo197775 + denovo197788 +
                                 denovo197839 + denovo197845 + denovo197857 + denovo197877 + denovo197880 + denovo197904 + denovo197933 + denovo197964 +
                                 denovo197988 + denovo197992 + denovo197994 + denovo198008 + denovo198060 + denovo198081 + denovo198098 + denovo198102 +
                                 denovo198124 + denovo198168 + denovo198213 + denovo198318 + denovo198392 + denovo198414 + denovo198471 + denovo198484 +
                                 denovo198494 + denovo198501 + denovo198563 + denovo198570 + denovo198587 + denovo198613 + denovo198614 + denovo198624 +
                                 denovo198648 + denovo198651 + denovo198663 + denovo198680 + denovo198685 + denovo198696 + denovo198726 + denovo198780 +
                                 denovo198781 + denovo198786 + denovo198790 + denovo198814 + denovo198858 + denovo198878 + denovo198897 + denovo198916 +
                                 denovo198928 + denovo198941 + denovo198962 + denovo198979 + denovo199007 + denovo199032 + denovo199052 + denovo199112 +
                                 denovo199116 + denovo199127 + denovo199144 + denovo199172 + denovo199184 + denovo199242 + denovo199287 + denovo199342 +
                                 denovo199344 + denovo199360 + denovo199373 + denovo199397 + denovo199404 + denovo199452 + denovo199454 + denovo199463 +
                                 denovo199515 + denovo199521 + denovo199546 + denovo199605 + denovo199611 + denovo199613 + denovo199715 + denovo199771 +
                                 denovo199783 + denovo199788 + denovo199807 + denovo199828 + denovo199854 + denovo199864 + denovo199881 + denovo199901 +
                                 denovo199912 + denovo199944 + denovo199955 + denovo200029 + denovo200044 + denovo200047 + denovo200072 + denovo200080 +
                                 denovo200091 + denovo200094 + denovo200104 + denovo200115 + denovo200174 + denovo200183 + denovo200189 + denovo200222 +
                                 denovo200258 + denovo200262 + denovo200270 + denovo200280 + denovo200310 + denovo200343 + denovo200351 + denovo200408 +
                                 denovo200422 + denovo200423 + denovo200439 + denovo200445 + denovo200447 + denovo200495 + denovo200510 + denovo200531 +
                                 denovo200543 + denovo200597 + denovo200602 + denovo200608 + denovo200617 + denovo200635 + denovo200705 + denovo200722 +
                                 denovo200767 + denovo200792 + denovo200839 + denovo200872 + denovo200893 + denovo200902 + denovo200905 + denovo200935 +
                                 denovo200947 + denovo200977 + denovo201015 + denovo201089 + denovo201109 + denovo201153 + denovo201162 + denovo201214 +
                                 denovo201266 + denovo201305 + denovo201328 + denovo201331 + denovo201370 + denovo201375 + denovo201376 + denovo201385 +
                                 denovo201419 + denovo201432 + denovo201492 + denovo201516 + denovo201534 + denovo201563 + denovo201571 + denovo201605 +
                                 denovo201622 + denovo201629 + denovo201748 + denovo201796 + denovo201803 + denovo201880 + denovo201943 + denovo202026 +
                                 denovo202032 + denovo202097 + denovo202150 + denovo202202 + denovo202219 + denovo202282 + denovo202396 + denovo202408 +
                                 denovo202444 + denovo202448 + denovo202479 + denovo202488 + denovo202504 + denovo202524 + denovo202544 + denovo202551 +
                                 denovo202567 + denovo202613 + denovo202654 + denovo202679 + denovo202704 + denovo202723 + denovo202731 + denovo202754 +
                                 denovo202846 + denovo202867 + denovo202871 + denovo202960 + denovo203003 + denovo203007 + denovo203065 + denovo203070 +
                                 denovo203079 + denovo203083 + denovo203117 + denovo203136 + denovo203151 + denovo203169 + denovo203172 + denovo203187 +
                                 denovo203209 + denovo203257 + denovo203273 + denovo203276 + denovo203284 + denovo203293 + denovo203295 + denovo203305 +
                                 denovo203317 + denovo203320 + denovo203348 + denovo203373 + denovo203426 + denovo203431 + denovo203514 + denovo203590 +
                                 denovo203634 + denovo203671 + denovo203812 + denovo203838 + denovo203891 + denovo203906 + denovo203909 + denovo203912 +
                                 denovo203965 + denovo204036 + denovo204045 + denovo204077 + denovo204089 + denovo204166 + denovo204189 + denovo204225 +
                                 denovo204237 + denovo204289 + denovo204292 + denovo204319 + denovo204396 + denovo204409 + denovo204419 + denovo204425 +
                                 denovo204457 + denovo204471 + denovo204475 + denovo204520 + denovo204571 + denovo204590 + denovo204610 + denovo204657 +
                                 denovo204694 + denovo204703 + denovo204771 + denovo204775 + denovo204802 + denovo204807 + denovo204831 + denovo204836 +
                                 denovo204840 + denovo204856 + denovo204882 + denovo204923 + denovo204965 + denovo204989 + denovo205029 + denovo205062 +
                                 denovo205113 + denovo205127 + denovo205183 + denovo205192 + denovo205217 + denovo205231 + denovo205239 + denovo205242 +
                                 denovo205246 + denovo205252 + denovo205261 + denovo205263 + denovo205312 + denovo205354 + denovo205367 + denovo205453 +
                                 denovo205491 + denovo205496 + denovo205519 + denovo205524 + denovo205530 + denovo205555 + denovo205564 + denovo205568 +
                                 denovo205580 + denovo205585 + denovo205594 + denovo205605 + denovo205653 + denovo205664 + denovo205688 + denovo205696 +
                                 denovo205703 + denovo205768 + denovo205769 + denovo205780 + denovo205796 + denovo205799 + denovo205836 + denovo205841 +
                                 denovo205868 + denovo205892 + denovo205894 + denovo205902 + denovo205924 + denovo205939 + denovo205946 + denovo205990 +
                                 denovo206025 + denovo206071 + denovo206078 + denovo206094 + denovo206100 + denovo206171 + denovo206216 + denovo206242 +
                                 denovo206252 + denovo206261 + denovo206315 + denovo206332 + denovo206345 + denovo206356 + denovo206360 + denovo206367 +
                                 denovo206454 + denovo206487 + denovo206494 + denovo206509 + denovo206542 + denovo206546 + denovo206554 + denovo206571 +
                                 denovo206576 + denovo206674 + denovo206686 + denovo206727 + denovo206743 + denovo206749 + denovo206816 + denovo206831 +
                                 denovo206896 + denovo206900 + denovo206942 + denovo206953 + denovo206957 + denovo207002 + denovo207059 + denovo207076 +
                                 denovo207080 + denovo207103 + denovo207118 + denovo207137 + denovo207169 + denovo207176 + denovo207182 + denovo207234 +
                                 denovo207240 + denovo207257 + denovo207259 + denovo207392 + denovo207412 + denovo207422 + denovo207449 + denovo207501 +
                                 denovo207544 + denovo207599 + denovo207624 + denovo207661 + denovo207674 + denovo207678 + denovo207691 + denovo207699 +
                                 denovo207725 + denovo207736 + denovo207778 + denovo207793 + denovo207817 + denovo207855 + denovo207868 + denovo207878 +
                                 denovo207884 + denovo207909 + denovo207910 + denovo207917 + denovo207938 + denovo207976 + denovo207978 + denovo208008 +
                                 denovo208009 + denovo208023 + denovo208047 + denovo208096 + denovo208124 + denovo208125 + denovo208135 + denovo208224 +
                                 denovo208234 + denovo208357 + denovo208378 + denovo208405 + denovo208437 + denovo208449 + denovo208457 + denovo208459 +
                                 denovo208483 + denovo208500 + denovo208560 + denovo208588 + denovo208644 + denovo208651 + denovo208680 + denovo208700 +
                                 denovo208733 + denovo208743 + denovo208819 + denovo208846 + denovo208849 + denovo208867 + denovo208902 + denovo208956 +
                                 denovo208981 + denovo208995 + denovo209005 + denovo209013 + denovo209014 + denovo209054 + denovo209084 + denovo209106 +
                                 denovo209116 + denovo209124 + denovo209175 + denovo209191 + denovo209325 + denovo209331 + denovo209346 + denovo209353 +
                                 denovo209405 + denovo209447 + denovo209497 + denovo209501 + denovo209513 + denovo209526 + denovo209530 + denovo209592 +
                                 denovo209618 + denovo209632 + denovo209671 + denovo209714 + denovo209758 + denovo209763 + denovo209774 + denovo209781 +
                                 denovo209782 + denovo209841 + denovo209855 + denovo209874 + denovo209888 + denovo209929 + denovo210026 + denovo210029 +
                                 denovo210048 + denovo210051 + denovo210093 + denovo210166 + denovo210189 + denovo210190 + denovo210200 + denovo210206 +
                                 denovo210236 + denovo210261 + denovo210263 + denovo210308 + denovo210318 + denovo210352 + denovo210370 + denovo210375 +
                                 denovo210395 + denovo210413 + denovo210423 + denovo210440 + denovo210454 + denovo210456 + denovo210467 + denovo210473 +
                                 denovo210530 + denovo210536 + denovo210564 + denovo210575 + denovo210586 + denovo210587 + denovo210610 + denovo210630 +
                                 denovo210668 + denovo210761 + denovo210779 + denovo210791 + denovo210855 + denovo210867 + denovo210883 + denovo210892 +
                                 denovo210936 + denovo210973 + denovo210989 + denovo211012 + denovo211014 + denovo211024 + denovo211032 + denovo211053 +
                                 denovo211116 + denovo211129 + denovo211137 + denovo211149 + denovo211157 + denovo211194 + denovo211208 + denovo211225 +
                                 denovo211263 + denovo211275 + denovo211287 + denovo211316 + denovo211340 + denovo211393 + denovo211410 + denovo211426 +
                                 denovo211456 + denovo211459 + denovo211497 + denovo211500 + denovo211516 + denovo211519 + denovo211565 + denovo211572 +
                                 denovo211599 + denovo211611 + denovo211650 + denovo211684 + denovo211745 + denovo211780 + denovo211839 + denovo211840 +
                                 denovo211862 + denovo211899 + denovo211900 + denovo211943 + denovo211955 + denovo211979 + denovo212022 + denovo212043 +
                                 denovo212062 + denovo212249 + denovo212253 + denovo212268 + denovo212281 + denovo212290 + denovo212299 + denovo212302 +
                                 denovo212321 + denovo212340 + denovo212377 + denovo212394 + denovo212408 + denovo212419 + denovo212423 + denovo212425 +
                                 denovo212443 + denovo212467 + denovo212488 + denovo212524 + denovo212528 + denovo212532 + denovo212588 + denovo212607 +
                                 denovo212612 + denovo212656 + denovo212720 + denovo212728 + denovo212732 + denovo212768 + denovo212821 + denovo212830 +
                                 denovo212887 + denovo212888 + denovo212948 + denovo212951 + denovo212967 + denovo212980 + denovo212995 + denovo213006 +
                                 denovo213053 + denovo213188 + denovo213197 + denovo213205 + denovo213213 + denovo213222 + denovo213285 + denovo213304 +
                                 denovo213311 + denovo213372 + denovo213380 + denovo213399 + denovo213416 + denovo213485 + denovo213508 + denovo213509 +
                                 denovo213537 + denovo213562 + denovo213576 + denovo213584 + denovo213592 + denovo213618 + denovo213625 + denovo213627 +
                                 denovo213708 + denovo213715 + denovo213719 + denovo213733 + denovo213817 + denovo213826 + denovo213850 + denovo213933 +
                                 denovo213955 + denovo213986 + denovo214061 + denovo214062 + denovo214064 + denovo214135 + denovo214150 + denovo214153 +
                                 denovo214175 + denovo214197 + denovo214199 + denovo214222 + denovo214241 + denovo214299 + denovo214354 + denovo214373 +
                                 denovo214397 + denovo214466 + denovo214468 + denovo214473 + denovo214507 + denovo214508 + denovo214568 + denovo214579 +
                                 denovo214584 + denovo214660 + denovo214668 + denovo214685 + denovo214704 + denovo214719 + denovo214741 + denovo214778 +
                                 denovo214829 + denovo214837 + denovo214845 + denovo214853 + denovo214883 + denovo214885 + denovo214898 + denovo214912 +
                                 denovo214920 + denovo214962 + denovo214984 + denovo215052 + denovo215058 + denovo215072 + denovo215091 + denovo215095 +
                                 denovo215146 + denovo215213 + denovo215229 + denovo215257 + denovo215271 + denovo215297 + denovo215299 + denovo215306 +
                                 denovo215321 + denovo215336 + denovo215364 + denovo215371 + denovo215388 + denovo215451 + denovo215454 + denovo215494 +
                                 denovo215513 + denovo215553 + denovo215602 + denovo215701 + denovo215712 + denovo215722 + denovo215731 + denovo215744 +
                                 denovo215755 + denovo215762 + denovo215772 + denovo215777 + denovo215797 + denovo215805 + denovo215832 + denovo215844 +
                                 denovo215863 + denovo215866 + denovo215903 + denovo215913 + denovo215937 + denovo215959 + denovo215983 + denovo216014 +
                                 denovo216015 + denovo216017 + denovo216030 + denovo216074 + denovo216075 + denovo216098 + denovo216131 + denovo216185 +
                                 denovo216195 + denovo216238 + denovo216241 + denovo216249 + denovo216270 + denovo216273 + denovo216291 + denovo216301 +
                                 denovo216350 + denovo216366 + denovo216373 + denovo216385 + denovo216407 + denovo216461 + denovo216469 + denovo216471 +
                                 denovo216474 + denovo216475 + denovo216491 + denovo216498 + denovo216502 + denovo216508 + denovo216517 + denovo216528 +
                                 denovo216542 + denovo216582 + denovo216605 + denovo216626 + denovo216640 + denovo216699 + denovo216766 + denovo216788 +
                                 denovo216804 + denovo216869 + denovo216886 + denovo216902 + denovo216974 + denovo216975 + denovo216981 + denovo216998 +
                                 denovo217073 + denovo217075 + denovo217083 + denovo217142 + denovo217151 + denovo217175 + denovo217193 + denovo217265 +
                                 denovo217296 + denovo217312 + denovo217355 + denovo217356 + denovo217384 + denovo217403 + denovo217422 + denovo217443 +
                                 denovo217478 + denovo217483 + denovo217496 + denovo217506 + denovo217513 + denovo217589 + denovo217615 + denovo217672 +
                                 denovo217704 + denovo217711 + denovo217746 + denovo217778 + denovo217793 + denovo217811 + denovo217866 + denovo217868 +
                                 denovo217906 + denovo217925 + denovo217995 + denovo218024 + denovo218046 + denovo218139 + denovo218141 + denovo218150 +
                                 denovo218159 + denovo218204 + denovo218217 + denovo218229 + denovo218235 + denovo218282 + denovo218284 + denovo218304 +
                                 denovo218335 + denovo218370 + denovo218446 + denovo218473 + denovo218635 + denovo218671 + denovo218727 + denovo218758 +
                                 denovo218780 + denovo218790 + denovo218791 + denovo218795 + denovo218838 + denovo218861 + denovo218875 + denovo218879 +
                                 denovo218895 + denovo218951 + denovo218996 + denovo219000 + denovo219004 + denovo219014 + denovo219015 + denovo219036 +
                                 denovo219040 + denovo219051 + denovo219108 + denovo219111 + denovo219219 + denovo219236 + denovo219253 + denovo219255 +
                                 denovo219272 + denovo219285 + denovo219341 + denovo219344 + denovo219378 + denovo219433 + denovo219453 + denovo219459 +
                                 denovo219463 + denovo219466 + denovo219507 + denovo219521 + denovo219537 + denovo219545 + denovo219571 + denovo219607 +
                                 denovo219636 + denovo219637 + denovo219704 + denovo219705 + denovo219721 + denovo219722 + denovo219781 + denovo219811 +
                                 denovo219823 + denovo219848 + denovo219851 + denovo219872 + denovo219965 + denovo219996 + denovo220096 + denovo220126 +
                                 denovo220166 + denovo220195 + denovo220233 + denovo220258 + denovo220266 + denovo220325 + denovo220430 + denovo220476 +
                                 denovo220486 + denovo220535 + denovo220573 + denovo220585 + denovo220588 + denovo220625 + denovo220652 + denovo220654 +
                                 denovo220662 + denovo220667 + denovo220689 + denovo220768 + denovo220771 + denovo220810 + denovo220818 + denovo220822 +
                                 denovo220844 + denovo220855 + denovo220866 + denovo220879 + denovo220941 + denovo220956 + denovo220968 + denovo221034 +
                                 denovo221039 + denovo221058 + denovo221061 + denovo221068 + denovo221086 + denovo221099 + denovo221129 + denovo221140 +
                                 denovo221177 + denovo221182 + denovo221188 + denovo221216 + denovo221232 + denovo221239 + denovo221280 + denovo221293 +
                                 denovo221311 + denovo221337 + denovo221340 + denovo221342 + denovo221347 + denovo221371 + denovo221377 + denovo221405 +
                                 denovo221423 + denovo221486 + denovo221591 + denovo221653 + denovo221658 + denovo221668 + denovo221707 + denovo221725 +
                                 denovo221761 + denovo221788 + denovo221801 + denovo221813 + denovo221817 + denovo221870 + denovo221930 + denovo222007 +
                                 denovo222040 + denovo222049 + denovo222064 + denovo222115 + denovo222117 + denovo222132 + denovo222143 + denovo222153 +
                                 denovo222160 + denovo222240 + denovo222270 + denovo222357 + denovo222362 + denovo222366 + denovo222435 + denovo222447 +
                                 denovo222483 + denovo222528 + denovo222558 + denovo222560 + denovo222579 + denovo222623 + denovo222683 + denovo222703 +
                                 denovo222736 + denovo222739 + denovo222751 + denovo222783 + denovo222867 + denovo222882 + denovo222961 + denovo223014 +
                                 denovo223039 + denovo223069 + denovo223103 + denovo223120 + denovo223130 + denovo223135 + denovo223180 + denovo223183 +
                                 denovo223262 + denovo223321 + denovo223366 + denovo223381 + denovo223395 + denovo223401 + denovo223418 + denovo223477 +
                                 denovo223522 + denovo223528 + denovo223568 + denovo223606 + denovo223624 + denovo223644 + denovo223646 + denovo223686 +
                                 denovo223720 + denovo223723 + denovo223757 + denovo223778 + denovo223801 + denovo223815 + denovo223821 + denovo223837 +
                                 denovo223841 + denovo223842 + denovo223949 + denovo224012 + denovo224037 + denovo224039 + denovo224088 + denovo224142 +
                                 denovo224166 + denovo224241 + denovo224242 + denovo224254 + denovo224261 + denovo224283 + denovo224291 + denovo224299 +
                                 denovo224308 + denovo224358 + denovo224369 + denovo224387 + denovo224404 + denovo224416 + denovo224420 + denovo224439 +
                                 denovo224484 + denovo224487 + denovo224552 + denovo224580 + denovo224604 + denovo224609 + denovo224676 + denovo224681 +
                                 denovo224806 + denovo224812 + denovo224837 + denovo224839 + denovo224890 + denovo224894 + denovo224898 + denovo224942 +
                                 denovo224978 + denovo225028 + denovo225042 + denovo225131 + denovo225200 + denovo225227 + denovo225228 + denovo225327 +
                                 denovo225351 + denovo225358 + denovo225400 + denovo225431 + denovo225439 + denovo225444 + denovo225635 + denovo225652 +
                                 denovo225658 + denovo225673 + denovo225693 + denovo225711 + denovo225713 + denovo225720 + denovo225723 + denovo225725 +
                                 denovo225803 + denovo225804 + denovo225808 + denovo225855 + denovo225859 + denovo225890 + denovo225937 + denovo225956 +
                                 denovo225964 + denovo225990 + denovo226021 + denovo226026 + denovo226093 + denovo226098 + denovo226133 + denovo226135, data = presence_absence_2016_feeding_habit, CV = TRUE)

DFA_predict <- table(presence_absence_2016_feeding_habit$FeedingHabit, DFA_fit_JK$class)

# Measure Total Correctly Predicted Individuals by Traditional Feeding Habit
diag(prop.table(DFA_predict, 1))


## 11) Relative Abundance - Discriminant Function Analysis for Water Velocity

# Discriminant Function Analysis Using Jacknifed Prediction
presence_absence_2016_water_velocity <- read.csv("/data/shawn/2016_Multivariate_Analysis/presence_absence_2016_water_velocity.csv")
presence_absence_2016_water_velocity <- as.data.frame(presence_absence_2016_water_velocity)
row.names(presence_absence_2016_water_velocity) <- presence_absence_2016_water_velocity$SampleID
presence_absence_2016_water_velocity <- presence_absence_2016_water_velocity[,-1]

DFA_fit_JK <- lda(WaterVelocity ~ denovo36 + denovo122 + denovo133 + denovo140 + denovo162 + denovo216 + denovo255 + 
                                 denovo332 + denovo367 + denovo391 + denovo402 + denovo435 + denovo448 + denovo465 + denovo478 +
                                 denovo511 + denovo535 + denovo544 + denovo567 + denovo584 + denovo624 + denovo680 +
                                 denovo707 + denovo728 + denovo730 + denovo747 + denovo837 + denovo922 + denovo924 + denovo946 +
                                 denovo964 + denovo991 + denovo994 + denovo1013 + denovo1027 + denovo1041 + denovo1047 + denovo1048 +
                                 denovo1054 + denovo1077 + denovo1106 + denovo1118 + denovo1142 + denovo1182 + denovo1211 +
                                 denovo1223 + denovo1229 + denovo1292 + denovo1332 + denovo1335 + denovo1339 + denovo1365 +
                                 denovo1407 + denovo1408 + denovo1487 + denovo1524 + denovo1533 + denovo1540 + denovo1545 + denovo1553 +
                                 denovo1621 + denovo1625 + denovo1632 + denovo1642 + denovo1664 + denovo1669 + denovo1707 + denovo1767 +
                                 denovo1773 + denovo1836 + denovo1839 + denovo1887 + denovo1941 + denovo1944 + denovo1989 + denovo1992 +
                                 denovo2011 + denovo2029 + denovo2120 + denovo2147 + denovo2179 + denovo2189 + denovo2245 + denovo2269 +
                                 denovo2290 + denovo2297 + denovo2351 + denovo2359 + denovo2367 + denovo2418 + denovo2465 + denovo2479 +
                                 denovo2491 + denovo2505 + denovo2549 + denovo2554 + denovo2556 + denovo2567 + denovo2581 + denovo2606 +
                                 denovo2619 + denovo2677 + denovo2713 + denovo2800 + denovo2808 + denovo2824 + denovo2835 + denovo2839 +
                                 denovo2848 + denovo2867 + denovo2952 + denovo3030 + denovo3065 + denovo3079 + denovo3151 +
                                 denovo3198 + denovo3246 + denovo3319 + denovo3335 + denovo3346 + denovo3363 + denovo3368 +
                                 denovo3406 + denovo3426 + denovo3432 + denovo3449 + denovo3461 + denovo3480 + denovo3489 + denovo3536 +
                                 denovo3599 + denovo3643 + denovo3652 + denovo3654 + denovo3673 + denovo3679 + denovo3700 + denovo3745 +
                                 denovo3768 + denovo3777 + denovo3838 + denovo3853 + denovo3858 + denovo3882 + denovo3929 + denovo3959 +
                                 denovo4025 + denovo4040 + denovo4074 + denovo4079 + denovo4082 + denovo4114 + denovo4154 +
                                 denovo4164 + denovo4205 + denovo4228 + denovo4236 + denovo4247 + denovo4283 + denovo4289 + denovo4312 +
                                 denovo4319 + denovo4325 + denovo4340 + denovo4352 + denovo4378 + denovo4408 + denovo4419 + denovo4424 +
                                 denovo4437 + denovo4440 + denovo4499 + denovo4501 + denovo4533 + denovo4575 + denovo4609 + denovo4631 +
                                 denovo4636 + denovo4711 + denovo4746 + denovo4759 + denovo4762 + denovo4770 + denovo4805 + denovo4834 +
                                 denovo4836 + denovo4859 + denovo4875 + denovo4900 + denovo4904 + denovo4984 + denovo4988 + denovo5067 +
                                 denovo5133 + denovo5168 + denovo5198 + denovo5204 + denovo5226 + denovo5247 + denovo5263 +
                                 denovo5269 + denovo5288 + denovo5293 + denovo5324 + denovo5414 + denovo5448 + denovo5509 + denovo5550 +
                                 denovo5581 + denovo5647 + denovo5690 + denovo5741 + denovo5776 + denovo5801 + denovo5859 +
                                 denovo5971 + denovo5975 + denovo5980 + denovo5986 + denovo6019 + denovo6138 + denovo6148 +
                                 denovo6150 + denovo6209 + denovo6214 + denovo6226 + denovo6255 + denovo6316 + denovo6402 + denovo6411 +
                                 denovo6438 + denovo6439 + denovo6451 + denovo6453 + denovo6472 + denovo6473 + denovo6493 + denovo6496 +
                                 denovo6499 + denovo6504 + denovo6520 + denovo6521 + denovo6527 + denovo6550 + denovo6573 + denovo6577 +
                                 denovo6583 + denovo6608 + denovo6660 + denovo6689 + denovo6709 + denovo6735 + denovo6749 + denovo6811 +
                                 denovo6816 + denovo6833 + denovo6835 + denovo6856 + denovo6872 + denovo6878 + denovo6909 + denovo6918 +
                                 denovo7003 + denovo7005 + denovo7007 + denovo7053 + denovo7060 + denovo7074 + denovo7086 + denovo7124 +
                                 denovo7136 + denovo7141 + denovo7162 + denovo7171 + denovo7186 + denovo7187 + denovo7195 + denovo7197 +
                                 denovo7204 + denovo7206 + denovo7238 + denovo7239 + denovo7251 + denovo7257 + denovo7288 + denovo7303 +
                                 denovo7334 + denovo7336 + denovo7342 + denovo7347 + denovo7359 + denovo7361 + denovo7369 + denovo7414 +
                                 denovo7417 + denovo7432 + denovo7456 + denovo7459 + denovo7465 + denovo7470 + denovo7473 + denovo7495 +
                                 denovo7508 + denovo7522 + denovo7528 + denovo7541 + denovo7551 + denovo7559 + denovo7672 + denovo7682 +
                                 denovo7729 + denovo7741 + denovo7753 + denovo7771 + denovo7794 + denovo7800 + denovo7841 + denovo7866 +
                                 denovo7880 + denovo7895 + denovo7908 + denovo7939 + denovo7953 + denovo7958 + denovo7966 + denovo7995 +
                                 denovo8002 + denovo8017 + denovo8034 + denovo8057 + denovo8087 + denovo8121 + denovo8126 + denovo8155 +
                                 denovo8181 + denovo8186 + denovo8188 + denovo8230 + denovo8263 + denovo8288 + denovo8303 + denovo8350 +
                                 denovo8352 + denovo8357 + denovo8385 + denovo8386 + denovo8405 + denovo8428 + denovo8463 + denovo8485 +
                                 denovo8502 + denovo8563 + denovo8582 + denovo8599 + denovo8602 + denovo8606 + denovo8696 +
                                 denovo8704 + denovo8714 + denovo8735 + denovo8741 + denovo8790 + denovo8817 + denovo8861 + denovo8879 +
                                 denovo8889 + denovo8906 + denovo8908 + denovo8929 + denovo8970 + denovo8978 + denovo8986 + denovo8991 +
                                 denovo9005 + denovo9017 + denovo9030 + denovo9075 + denovo9097 + denovo9115 + denovo9118 + denovo9181 +
                                 denovo9204 + denovo9206 + denovo9216 + denovo9227 + denovo9248 + denovo9265 + denovo9271 + denovo9276 +
                                 denovo9341 + denovo9350 + denovo9366 + denovo9409 + denovo9431 + denovo9457 + denovo9464 + denovo9476 +
                                 denovo9501 + denovo9502 + denovo9606 + denovo9615 + denovo9636 + denovo9678 + denovo9684 +
                                 denovo9687 + denovo9711 + denovo9726 + denovo9758 + denovo9761 + denovo9771 + denovo9802 + denovo9809 +
                                 denovo9812 + denovo9900 + denovo9914 + denovo9920 + denovo9932 + denovo9953 + denovo9957 +
                                 denovo9967 + denovo10061 + denovo10127 + denovo10155 + denovo10157 + denovo10176 + denovo10181 +
                                 denovo10213 + denovo10226 + denovo10227 + denovo10241 + denovo10264 + denovo10282 + denovo10315 + denovo10329 +
                                 denovo10333 + denovo10352 + denovo10366 + denovo10376 + denovo10383 + denovo10395 + denovo10447 + denovo10452 +
                                 denovo10464 + denovo10487 + denovo10516 + denovo10559 + denovo10580 + denovo10581 + denovo10602 + denovo10618 +
                                 denovo10621 + denovo10634 + denovo10639 + denovo10646 + denovo10647 + denovo10648 + denovo10650 + denovo10675 +
                                 denovo10689 + denovo10718 + denovo10732 + denovo10752 + denovo10831 + denovo10842 + denovo10866 + denovo10889 +
                                 denovo10898 + denovo10992 + denovo11010 + denovo11038 + denovo11050 + denovo11051 + denovo11073 + denovo11095 +
                                 denovo11152 + denovo11191 + denovo11209 + denovo11229 + denovo11234 + denovo11235 + denovo11245 + denovo11247 +
                                 denovo11249 + denovo11286 + denovo11290 + denovo11329 + denovo11367 + denovo11403 + denovo11410 +
                                 denovo11419 + denovo11430 + denovo11505 + denovo11556 + denovo11574 + denovo11588 + denovo11597 + denovo11623 +
                                 denovo11629 + denovo11646 + denovo11689 + denovo11730 + denovo11772 + denovo11778 + denovo11784 + denovo11789 +
                                 denovo11814 + denovo11817 + denovo11820 + denovo11844 + denovo11866 + denovo11873 + denovo11890 + denovo11958 +
                                 denovo11971 + denovo11980 + denovo12022 + denovo12027 + denovo12030 + denovo12052 + denovo12062 + denovo12087 +
                                 denovo12137 + denovo12185 + denovo12192 + denovo12221 + denovo12242 + denovo12248 + denovo12255 + denovo12265 +
                                 denovo12274 + denovo12275 + denovo12407 + denovo12420 + denovo12474 + denovo12497 + denovo12548 + denovo12609 +
                                 denovo12621 + denovo12641 + denovo12649 + denovo12658 + denovo12690 + denovo12693 + denovo12728 + denovo12736 +
                                 denovo12765 + denovo12777 + denovo12783 + denovo12808 + denovo12842 + denovo12870 + denovo12873 +
                                 denovo12881 + denovo12890 + denovo12913 + denovo12943 + denovo12952 + denovo12960 + denovo12963 + denovo13006 +
                                 denovo13017 + denovo13068 + denovo13075 + denovo13081 + denovo13096 + denovo13120 + denovo13121 + denovo13148 +
                                 denovo13183 + denovo13197 + denovo13265 + denovo13269 + denovo13312 + denovo13334 + denovo13336 + denovo13348 +
                                 denovo13417 + denovo13466 + denovo13515 + denovo13517 + denovo13524 + denovo13529 + denovo13536 +
                                 denovo13547 + denovo13569 + denovo13641 + denovo13669 + denovo13673 + denovo13708 + denovo13710 + denovo13767 +
                                 denovo13781 + denovo13800 + denovo13838 + denovo13877 + denovo13900 + denovo13928 + denovo13935 + denovo13948 +
                                 denovo13952 + denovo13996 + denovo14012 + denovo14019 + denovo14029 + denovo14044 + denovo14089 +
                                 denovo14125 + denovo14153 + denovo14202 + denovo14214 + denovo14249 + denovo14250 + denovo14251 + denovo14255 +
                                 denovo14273 + denovo14278 + denovo14294 + denovo14308 + denovo14329 + denovo14332 + denovo14346 +
                                 denovo14358 + denovo14400 + denovo14411 + denovo14544 + denovo14565 + denovo14640 + denovo14642 + denovo14666 +
                                 denovo14679 + denovo14703 + denovo14705 + denovo14713 + denovo14781 + denovo14802 + denovo14803 + denovo14809 +
                                 denovo14818 + denovo14820 + denovo14821 + denovo14971 + denovo14978 + denovo14982 + denovo14991 +
                                 denovo15000 + denovo15015 + denovo15066 + denovo15076 + denovo15101 + denovo15150 + denovo15151 +
                                 denovo15185 + denovo15212 + denovo15284 + denovo15307 + denovo15356 + denovo15359 + denovo15368 + denovo15371 +
                                 denovo15376 + denovo15430 + denovo15438 + denovo15449 + denovo15450 + denovo15509 + denovo15558 + denovo15581 +
                                 denovo15601 + denovo15711 + denovo15713 + denovo15714 + denovo15729 + denovo15764 + denovo15801 + denovo15811 +
                                 denovo15835 + denovo15849 + denovo15855 + denovo15895 + denovo15899 + denovo15968 + denovo15998 + denovo16021 +
                                 denovo16058 + denovo16064 + denovo16075 + denovo16081 + denovo16177 + denovo16182 + denovo16194 + denovo16353 +
                                 denovo16407 + denovo16437 + denovo16466 + denovo16470 + denovo16504 + denovo16516 + denovo16524 + denovo16528 +
                                 denovo16552 + denovo16577 + denovo16584 + denovo16603 + denovo16678 + denovo16722 + denovo16737 + denovo16746 +
                                 denovo16763 + denovo16782 + denovo16812 + denovo16819 + denovo16828 + denovo16857 + denovo16859 + denovo16861 +
                                 denovo16914 + denovo16915 + denovo16918 + denovo16923 + denovo16979 + denovo17022 + denovo17052 + denovo17078 +
                                 denovo17092 + denovo17190 + denovo17194 + denovo17232 + denovo17277 + denovo17292 + denovo17303 + denovo17304 +
                                 denovo17324 + denovo17336 + denovo17365 + denovo17373 + denovo17402 + denovo17436 + denovo17515 + denovo17562 +
                                 denovo17628 + denovo17633 + denovo17634 + denovo17645 + denovo17646 + denovo17683 + denovo17691 + denovo17713 +
                                 denovo17735 + denovo17762 + denovo17769 + denovo17836 + denovo17846 + denovo17856 + denovo17859 + denovo17874 +
                                 denovo17891 + denovo17940 + denovo17955 + denovo17958 + denovo17967 + denovo17969 + denovo17990 + denovo17991 +
                                 denovo18014 + denovo18038 + denovo18090 + denovo18093 + denovo18099 + denovo18105 + denovo18163 + denovo18176 +
                                 denovo18264 + denovo18296 + denovo18301 + denovo18302 + denovo18355 + denovo18368 + denovo18375 + denovo18416 +
                                 denovo18418 + denovo18452 + denovo18488 + denovo18498 + denovo18534 + denovo18539 + denovo18588 + denovo18589 +
                                 denovo18633 + denovo18674 + denovo18680 + denovo18716 + denovo18739 + denovo18760 + denovo18821 + denovo18848 +
                                 denovo18874 + denovo18906 + denovo18941 + denovo18975 + denovo18980 + denovo18994 + denovo19077 + denovo19096 +
                                 denovo19124 + denovo19146 + denovo19218 + denovo19265 + denovo19290 + denovo19392 + denovo19400 + denovo19416 +
                                 denovo19450 + denovo19499 + denovo19517 + denovo19560 + denovo19576 + denovo19619 + denovo19650 + denovo19664 +
                                 denovo19667 + denovo19715 + denovo19732 + denovo19773 + denovo19775 + denovo19829 + denovo19848 + denovo19892 +
                                 denovo19909 + denovo20052 + denovo20065 + denovo20075 + denovo20111 + denovo20165 + denovo20196 + denovo20215 +
                                 denovo20226 + denovo20234 + denovo20275 + denovo20282 + denovo20290 + denovo20345 + denovo20361 + denovo20371 +
                                 denovo20372 + denovo20374 + denovo20403 + denovo20451 + denovo20462 + denovo20466 + denovo20491 + denovo20522 +
                                 denovo20530 + denovo20586 + denovo20606 + denovo20682 + denovo20687 + denovo20693 + denovo20701 + denovo20741 +
                                 denovo20808 + denovo20813 + denovo20904 + denovo20940 + denovo20997 + denovo21001 + denovo21010 + denovo21049 +
                                 denovo21064 + denovo21071 + denovo21075 + denovo21085 + denovo21094 + denovo21096 + denovo21098 + denovo21150 +
                                 denovo21155 + denovo21169 + denovo21173 + denovo21178 + denovo21192 + denovo21200 + denovo21348 + denovo21353 +
                                 denovo21378 + denovo21389 + denovo21413 + denovo21425 + denovo21471 + denovo21491 + denovo21499 + denovo21525 +
                                 denovo21543 + denovo21568 + denovo21580 + denovo21587 + denovo21627 + denovo21634 + denovo21638 + denovo21642 +
                                 denovo21659 + denovo21700 + denovo21726 + denovo21739 + denovo21762 + denovo21767 + denovo21793 + denovo21843 +
                                 denovo21858 + denovo21874 + denovo21891 + denovo21917 + denovo21924 + denovo21953 + denovo21989 + denovo22004 +
                                 denovo22022 + denovo22029 + denovo22031 + denovo22065 + denovo22070 + denovo22110 + denovo22137 + denovo22139 +
                                 denovo22171 + denovo22205 + denovo22214 + denovo22233 + denovo22294 + denovo22309 + denovo22314 +
                                 denovo22330 + denovo22338 + denovo22344 + denovo22365 + denovo22423 + denovo22503 + denovo22539 + denovo22553 +
                                 denovo22561 + denovo22596 + denovo22629 + denovo22695 + denovo22721 + denovo22821 + denovo22850 + denovo22879 +
                                 denovo22881 + denovo22923 + denovo22934 + denovo22941 + denovo22972 + denovo22974 + denovo22975 + denovo23006 +
                                 denovo23064 + denovo23106 + denovo23118 + denovo23128 + denovo23143 + denovo23155 + denovo23161 + denovo23173 +
                                 denovo23184 + denovo23191 + denovo23192 + denovo23245 + denovo23266 + denovo23274 + denovo23280 + denovo23338 +
                                 denovo23345 + denovo23374 + denovo23377 + denovo23391 + denovo23424 + denovo23453 + denovo23509 +
                                 denovo23544 + denovo23552 + denovo23560 + denovo23562 + denovo23591 + denovo23625 + denovo23691 + denovo23702 +
                                 denovo23712 + denovo23734 + denovo23741 + denovo23769 + denovo23798 + denovo23816 + denovo23833 + denovo23893 +
                                 denovo23912 + denovo23954 + denovo23971 + denovo23986 + denovo23992 + denovo23999 + denovo24032 + denovo24045 +
                                 denovo24064 + denovo24080 + denovo24082 + denovo24117 + denovo24189 + denovo24262 + denovo24275 + denovo24338 +
                                 denovo24390 + denovo24397 + denovo24433 + denovo24439 + denovo24533 + denovo24547 + denovo24552 + denovo24608 +
                                 denovo24641 + denovo24647 + denovo24660 + denovo24668 + denovo24742 + denovo24751 + denovo24757 + denovo24758 +
                                 denovo24782 + denovo24869 + denovo24872 + denovo24881 + denovo24900 + denovo24939 + denovo24958 + denovo25001 +
                                 denovo25027 + denovo25105 + denovo25137 + denovo25154 + denovo25205 + denovo25228 + denovo25229 + denovo25247 +
                                 denovo25261 + denovo25280 + denovo25284 + denovo25306 + denovo25347 + denovo25367 + denovo25388 + denovo25393 +
                                 denovo25427 + denovo25429 + denovo25464 + denovo25493 + denovo25582 + denovo25590 + denovo25594 + denovo25662 +
                                 denovo25699 + denovo25700 + denovo25726 + denovo25761 + denovo25797 + denovo25890 + denovo25903 + denovo25927 +
                                 denovo25990 + denovo25995 + denovo26040 + denovo26076 + denovo26099 + denovo26102 + denovo26135 +
                                 denovo26184 + denovo26220 + denovo26268 + denovo26278 + denovo26310 + denovo26361 + denovo26424 + denovo26453 +
                                 denovo26471 + denovo26474 + denovo26511 + denovo26541 + denovo26561 + denovo26585 + denovo26591 + denovo26602 +
                                 denovo26610 + denovo26633 + denovo26683 + denovo26719 + denovo26752 + denovo26785 + denovo26893 + denovo26902 +
                                 denovo26927 + denovo26954 + denovo26979 + denovo27015 + denovo27024 + denovo27039 + denovo27045 + denovo27054 +
                                 denovo27066 + denovo27126 + denovo27139 + denovo27170 + denovo27218 + denovo27221 + denovo27254 + denovo27289 +
                                 denovo27354 + denovo27355 + denovo27412 + denovo27461 + denovo27504 + denovo27534 + denovo27535 + denovo27603 +
                                 denovo27648 + denovo27673 + denovo27682 + denovo27700 + denovo27702 + denovo27725 + denovo27869 + denovo27874 +
                                 denovo27878 + denovo27900 + denovo27928 + denovo27946 + denovo27968 + denovo28021 + denovo28113 + denovo28128 +
                                 denovo28131 + denovo28253 + denovo28266 + denovo28333 + denovo28387 + denovo28392 + denovo28451 + denovo28485 +
                                 denovo28538 + denovo28595 + denovo28605 + denovo28613 + denovo28614 + denovo28625 + denovo28633 +
                                 denovo28699 + denovo28732 + denovo28750 + denovo28787 + denovo28790 + denovo28797 + denovo28838 +
                                 denovo28848 + denovo28885 + denovo28979 + denovo29009 + denovo29016 + denovo29038 + denovo29062 +
                                 denovo29083 + denovo29163 + denovo29188 + denovo29294 + denovo29310 + denovo29327 + denovo29355 +
                                 denovo29362 + denovo29381 + denovo29384 + denovo29388 + denovo29399 + denovo29420 + denovo29437 + denovo29438 +
                                 denovo29442 + denovo29459 + denovo29567 + denovo29620 + denovo29639 + denovo29642 + denovo29724 + denovo29754 +
                                 denovo29762 + denovo29793 + denovo29812 + denovo29821 + denovo29841 + denovo29842 + denovo29844 + denovo29920 +
                                 denovo29976 + denovo30013 + denovo30050 + denovo30054 + denovo30058 + denovo30061 + denovo30083 +
                                 denovo30112 + denovo30146 + denovo30157 + denovo30167 + denovo30186 + denovo30218 + denovo30224 + denovo30283 +
                                 denovo30317 + denovo30366 + denovo30371 + denovo30383 + denovo30393 + denovo30420 + denovo30484 +
                                 denovo30573 + denovo30580 + denovo30592 + denovo30635 + denovo30693 + denovo30715 + denovo30728 + denovo30734 +
                                 denovo30756 + denovo30777 + denovo30782 + denovo30811 + denovo30875 + denovo30883 + denovo30919 + denovo30940 +
                                 denovo30941 + denovo31003 + denovo31036 + denovo31052 + denovo31070 + denovo31132 + denovo31177 + denovo31274 +
                                 denovo31322 + denovo31348 + denovo31376 + denovo31385 + denovo31403 + denovo31414 + denovo31420 + denovo31430 +
                                 denovo31449 + denovo31460 + denovo31476 + denovo31477 + denovo31486 + denovo31487 + denovo31506 + denovo31534 +
                                 denovo31546 + denovo31557 + denovo31574 + denovo31576 + denovo31669 + denovo31673 + denovo31688 + denovo31746 +
                                 denovo31773 + denovo31820 + denovo31843 + denovo31846 + denovo31904 + denovo31934 + denovo31944 + denovo31959 +
                                 denovo31968 + denovo31980 + denovo31993 + denovo31999 + denovo32001 + denovo32013 + denovo32049 + denovo32058 +
                                 denovo32062 + denovo32079 + denovo32084 + denovo32094 + denovo32099 + denovo32199 + denovo32203 + denovo32266 +
                                 denovo32269 + denovo32278 + denovo32304 + denovo32477 + denovo32498 + denovo32547 + denovo32564 + denovo32569 +
                                 denovo32601 + denovo32656 + denovo32679 + denovo32738 + denovo32824 + denovo32855 + denovo32901 + denovo32956 +
                                 denovo33007 + denovo33021 + denovo33022 + denovo33074 + denovo33136 + denovo33171 + denovo33245 +
                                 denovo33246 + denovo33316 + denovo33340 + denovo33347 + denovo33360 + denovo33376 + denovo33440 +
                                 denovo33472 + denovo33491 + denovo33512 + denovo33518 + denovo33524 + denovo33571 + denovo33607 +
                                 denovo33616 + denovo33632 + denovo33644 + denovo33702 + denovo33729 + denovo33733 + denovo33757 + denovo33778 +
                                 denovo33793 + denovo33800 + denovo33831 + denovo33834 + denovo33844 + denovo33877 + denovo33891 + denovo33925 +
                                 denovo33941 + denovo34003 + denovo34012 + denovo34019 + denovo34040 + denovo34041 + denovo34049 + denovo34061 +
                                 denovo34111 + denovo34149 + denovo34244 + denovo34260 + denovo34337 + denovo34351 + denovo34364 + denovo34369 +
                                 denovo34407 + denovo34455 + denovo34471 + denovo34492 + denovo34499 + denovo34513 + denovo34525 + denovo34537 +
                                 denovo34585 + denovo34616 + denovo34636 + denovo34680 + denovo34713 + denovo34718 + denovo34725 + denovo34751 +
                                 denovo34763 + denovo34776 + denovo34808 + denovo34816 + denovo34840 + denovo34852 + denovo34873 + denovo34881 +
                                 denovo34918 + denovo34937 + denovo34998 + denovo35012 + denovo35022 + denovo35036 + denovo35062 + denovo35083 +
                                 denovo35089 + denovo35121 + denovo35124 + denovo35129 + denovo35133 + denovo35140 + denovo35153 + denovo35174 +
                                 denovo35196 + denovo35207 + denovo35216 + denovo35227 + denovo35246 + denovo35348 + denovo35369 + denovo35381 +
                                 denovo35415 + denovo35426 + denovo35451 + denovo35470 + denovo35487 + denovo35540 + denovo35585 + denovo35609 +
                                 denovo35655 + denovo35689 + denovo35729 + denovo35734 + denovo35756 + denovo35763 + denovo35790 + denovo35861 +
                                 denovo35862 + denovo35910 + denovo35917 + denovo35969 + denovo36019 + denovo36051 + denovo36059 +
                                 denovo36156 + denovo36166 + denovo36224 + denovo36231 + denovo36273 + denovo36279 + denovo36295 + denovo36392 +
                                 denovo36397 + denovo36401 + denovo36409 + denovo36459 + denovo36484 + denovo36503 + denovo36505 + denovo36519 +
                                 denovo36522 + denovo36539 + denovo36541 + denovo36545 + denovo36556 + denovo36565 + denovo36602 + denovo36603 +
                                 denovo36617 + denovo36639 + denovo36661 + denovo36671 + denovo36676 + denovo36680 + denovo36684 + denovo36725 +
                                 denovo36788 + denovo36805 + denovo36855 + denovo36892 + denovo36939 + denovo37007 + denovo37009 + denovo37046 +
                                 denovo37143 + denovo37144 + denovo37176 + denovo37230 + denovo37232 + denovo37281 + denovo37289 + denovo37328 +
                                 denovo37335 + denovo37352 + denovo37354 + denovo37357 + denovo37361 + denovo37374 + denovo37394 + denovo37483 +
                                 denovo37554 + denovo37581 + denovo37584 + denovo37596 + denovo37605 + denovo37637 + denovo37656 + denovo37659 +
                                 denovo37680 + denovo37683 + denovo37685 + denovo37717 + denovo37740 + denovo37764 + denovo37786 +
                                 denovo37799 + denovo37800 + denovo37811 + denovo37819 + denovo37831 + denovo37850 + denovo37896 + denovo37933 +
                                 denovo37978 + denovo37983 + denovo37996 + denovo38175 + denovo38177 + denovo38192 + denovo38195 +
                                 denovo38198 + denovo38251 + denovo38259 + denovo38279 + denovo38297 + denovo38301 + denovo38303 + denovo38319 +
                                 denovo38347 + denovo38364 + denovo38434 + denovo38460 + denovo38479 + denovo38500 + denovo38517 + denovo38523 +
                                 denovo38524 + denovo38561 + denovo38576 + denovo38634 + denovo38640 + denovo38656 + denovo38692 + denovo38695 +
                                 denovo38726 + denovo38732 + denovo38737 + denovo38762 + denovo38815 + denovo38826 + denovo38868 + denovo38905 +
                                 denovo39013 + denovo39032 + denovo39039 + denovo39049 + denovo39090 + denovo39140 + denovo39143 + denovo39151 +
                                 denovo39164 + denovo39182 + denovo39216 + denovo39243 + denovo39281 + denovo39304 + denovo39307 + denovo39344 +
                                 denovo39432 + denovo39435 + denovo39464 + denovo39476 + denovo39485 + denovo39507 + denovo39569 + denovo39601 +
                                 denovo39623 + denovo39649 + denovo39674 + denovo39720 + denovo39729 + denovo39798 + denovo39817 + denovo39942 +
                                 denovo39955 + denovo39996 + denovo40008 + denovo40057 + denovo40065 + denovo40067 + denovo40071 + denovo40078 +
                                 denovo40121 + denovo40125 + denovo40136 + denovo40138 + denovo40143 + denovo40154 + denovo40171 + denovo40174 +
                                 denovo40212 + denovo40246 + denovo40248 + denovo40257 + denovo40290 + denovo40350 + denovo40395 +
                                 denovo40403 + denovo40437 + denovo40479 + denovo40480 + denovo40510 + denovo40526 + denovo40569 + denovo40570 +
                                 denovo40608 + denovo40637 + denovo40728 + denovo40744 + denovo40755 + denovo40774 + denovo40780 +
                                 denovo40843 + denovo40910 + denovo40932 + denovo40971 + denovo41003 + denovo41027 + denovo41029 +
                                 denovo41067 + denovo41108 + denovo41114 + denovo41164 + denovo41166 + denovo41169 + denovo41172 + denovo41182 +
                                 denovo41235 + denovo41238 + denovo41240 + denovo41287 + denovo41299 + denovo41306 + denovo41325 + denovo41363 +
                                 denovo41381 + denovo41393 + denovo41425 + denovo41433 + denovo41499 + denovo41548 + denovo41553 + denovo41626 +
                                 denovo41637 + denovo41671 + denovo41707 + denovo41719 + denovo41738 + denovo41751 + denovo41806 +
                                 denovo41845 + denovo41868 + denovo41888 + denovo41911 + denovo42002 + denovo42005 + denovo42058 + denovo42157 +
                                 denovo42159 + denovo42193 + denovo42253 + denovo42316 + denovo42327 + denovo42340 + denovo42364 + denovo42381 +
                                 denovo42405 + denovo42419 + denovo42436 + denovo42444 + denovo42456 + denovo42485 + denovo42521 + denovo42538 +
                                 denovo42592 + denovo42632 + denovo42711 + denovo42716 + denovo42905 + denovo42931 + denovo42936 + denovo42944 +
                                 denovo42957 + denovo42976 + denovo42979 + denovo42993 + denovo42997 + denovo43004 + denovo43005 + denovo43018 +
                                 denovo43065 + denovo43106 + denovo43124 + denovo43133 + denovo43190 + denovo43231 + denovo43232 + denovo43258 +
                                 denovo43259 + denovo43280 + denovo43304 + denovo43350 + denovo43357 + denovo43373 + denovo43425 + denovo43475 +
                                 denovo43561 + denovo43593 + denovo43601 + denovo43671 + denovo43711 + denovo43719 + denovo43806 + denovo43812 +
                                 denovo43854 + denovo43857 + denovo43970 + denovo43977 + denovo43984 + denovo43997 + denovo44016 + denovo44034 +
                                 denovo44038 + denovo44077 + denovo44079 + denovo44093 + denovo44095 + denovo44146 + denovo44163 +
                                 denovo44234 + denovo44260 + denovo44276 + denovo44289 + denovo44301 + denovo44340 + denovo44358 + denovo44414 +
                                 denovo44496 + denovo44499 + denovo44502 + denovo44510 + denovo44512 + denovo44527 + denovo44528 + denovo44532 +
                                 denovo44538 + denovo44545 + denovo44609 + denovo44662 + denovo44669 + denovo44678 + denovo44696 + denovo44702 +
                                 denovo44708 + denovo44715 + denovo44726 + denovo44745 + denovo44793 + denovo44816 + denovo44863 + denovo44867 +
                                 denovo44912 + denovo44940 + denovo44979 + denovo44994 + denovo45004 + denovo45013 + denovo45018 + denovo45024 +
                                 denovo45046 + denovo45063 + denovo45097 + denovo45108 + denovo45110 + denovo45140 + denovo45153 + denovo45175 +
                                 denovo45194 + denovo45195 + denovo45202 + denovo45230 + denovo45242 + denovo45298 + denovo45302 + denovo45316 +
                                 denovo45322 + denovo45419 + denovo45422 + denovo45437 + denovo45453 + denovo45463 + denovo45476 + denovo45516 +
                                 denovo45590 + denovo45637 + denovo45642 + denovo45680 + denovo45690 + denovo45701 + denovo45702 + denovo45707 +
                                 denovo45708 + denovo45770 + denovo45783 + denovo45876 + denovo45885 + denovo45898 + denovo45923 + denovo45928 +
                                 denovo45939 + denovo45967 + denovo45971 + denovo46052 + denovo46097 + denovo46109 + denovo46141 +
                                 denovo46165 + denovo46188 + denovo46203 + denovo46226 + denovo46235 + denovo46263 + denovo46268 + denovo46293 +
                                 denovo46322 + denovo46336 + denovo46344 + denovo46346 + denovo46373 + denovo46466 + denovo46472 + denovo46516 +
                                 denovo46542 + denovo46555 + denovo46579 + denovo46598 + denovo46654 + denovo46655 + denovo46656 + denovo46687 +
                                 denovo46702 + denovo46721 + denovo46816 + denovo46836 + denovo46870 + denovo46874 + denovo46879 + denovo46883 +
                                 denovo46906 + denovo46926 + denovo46934 + denovo46966 + denovo46984 + denovo47032 + denovo47041 + denovo47049 +
                                 denovo47052 + denovo47066 + denovo47109 + denovo47123 + denovo47259 + denovo47312 + denovo47382 + denovo47441 +
                                 denovo47510 + denovo47528 + denovo47557 + denovo47587 + denovo47595 + denovo47627 + denovo47661 + denovo47667 +
                                 denovo47682 + denovo47717 + denovo47731 + denovo47801 + denovo47803 + denovo47808 + denovo47816 + denovo47824 +
                                 denovo47826 + denovo47852 + denovo47880 + denovo47891 + denovo47919 + denovo47944 + denovo48005 + denovo48013 +
                                 denovo48031 + denovo48087 + denovo48101 + denovo48163 + denovo48177 + denovo48209 + denovo48210 + denovo48288 +
                                 denovo48298 + denovo48323 + denovo48353 + denovo48384 + denovo48434 + denovo48442 + denovo48449 + denovo48458 +
                                 denovo48488 + denovo48499 + denovo48520 + denovo48531 + denovo48543 + denovo48546 + denovo48554 + denovo48579 +
                                 denovo48582 + denovo48616 + denovo48660 + denovo48676 + denovo48708 + denovo48726 + denovo48747 + denovo48764 +
                                 denovo48766 + denovo48816 + denovo48846 + denovo48861 + denovo48863 + denovo48910 + denovo48915 + denovo48972 +
                                 denovo49012 + denovo49063 + denovo49070 + denovo49084 + denovo49086 + denovo49107 + denovo49111 + denovo49122 +
                                 denovo49127 + denovo49157 + denovo49176 + denovo49270 + denovo49280 + denovo49298 + denovo49364 + denovo49422 +
                                 denovo49447 + denovo49460 + denovo49518 + denovo49541 + denovo49572 + denovo49594 + denovo49602 +
                                 denovo49628 + denovo49679 + denovo49698 + denovo49735 + denovo49768 + denovo49782 + denovo49789 + denovo49803 +
                                 denovo49847 + denovo49872 + denovo49914 + denovo49942 + denovo49975 + denovo50005 + denovo50043 + denovo50053 +
                                 denovo50125 + denovo50128 + denovo50166 + denovo50185 + denovo50189 + denovo50210 + denovo50218 + denovo50222 +
                                 denovo50229 + denovo50243 + denovo50263 + denovo50284 + denovo50333 + denovo50368 + denovo50400 + denovo50434 +
                                 denovo50458 + denovo50542 + denovo50561 + denovo50626 + denovo50628 + denovo50698 + denovo50722 +
                                 denovo50727 + denovo50740 + denovo50741 + denovo50771 + denovo50802 + denovo50815 + denovo50818 + denovo50832 +
                                 denovo50880 + denovo50893 + denovo50894 + denovo50921 + denovo50951 + denovo50974 + denovo50985 +
                                 denovo51035 + denovo51045 + denovo51088 + denovo51117 + denovo51123 + denovo51153 + denovo51156 + denovo51159 +
                                 denovo51197 + denovo51213 + denovo51244 + denovo51273 + denovo51288 + denovo51323 + denovo51376 + denovo51405 +
                                 denovo51433 + denovo51477 + denovo51537 + denovo51540 + denovo51549 + denovo51562 + denovo51565 + denovo51599 +
                                 denovo51645 + denovo51650 + denovo51701 + denovo51721 + denovo51745 + denovo51756 + denovo51785 + denovo51787 +
                                 denovo51796 + denovo51826 + denovo51830 + denovo51843 + denovo51853 + denovo51878 + denovo51921 + denovo52036 +
                                 denovo52074 + denovo52146 + denovo52155 + denovo52226 + denovo52251 + denovo52256 + denovo52282 + denovo52295 +
                                 denovo52337 + denovo52345 + denovo52348 + denovo52363 + denovo52377 + denovo52433 + denovo52488 + denovo52528 +
                                 denovo52538 + denovo52542 + denovo52552 + denovo52562 + denovo52601 + denovo52605 + denovo52630 + denovo52631 +
                                 denovo52758 + denovo52780 + denovo52784 + denovo52836 + denovo52847 + denovo52851 + denovo52880 + denovo52891 +
                                 denovo52902 + denovo52934 + denovo52935 + denovo52967 + denovo53013 + denovo53024 + denovo53064 +
                                 denovo53073 + denovo53103 + denovo53110 + denovo53122 + denovo53139 + denovo53149 + denovo53216 + denovo53236 +
                                 denovo53250 + denovo53269 + denovo53278 + denovo53327 + denovo53332 + denovo53346 + denovo53356 +
                                 denovo53360 + denovo53364 + denovo53401 + denovo53439 + denovo53451 + denovo53468 + denovo53504 + denovo53550 +
                                 denovo53564 + denovo53626 + denovo53643 + denovo53665 + denovo53671 + denovo53772 + denovo53821 + denovo53826 +
                                 denovo53835 + denovo53852 + denovo53866 + denovo53913 + denovo53965 + denovo53967 + denovo54032 + denovo54068 +
                                 denovo54077 + denovo54159 + denovo54176 + denovo54181 + denovo54234 + denovo54251 + denovo54252 + denovo54266 +
                                 denovo54269 + denovo54296 + denovo54311 + denovo54321 + denovo54361 + denovo54433 + denovo54434 + denovo54440 +
                                 denovo54484 + denovo54527 + denovo54532 + denovo54582 + denovo54598 + denovo54605 + denovo54617 + denovo54635 +
                                 denovo54651 + denovo54652 + denovo54670 + denovo54697 + denovo54727 + denovo54792 + denovo54803 + denovo54806 +
                                 denovo54875 + denovo54883 + denovo54886 + denovo54900 + denovo54911 + denovo54917 + denovo54932 + denovo54998 +
                                 denovo55042 + denovo55134 + denovo55184 + denovo55189 + denovo55202 + denovo55260 + denovo55267 +
                                 denovo55269 + denovo55274 + denovo55316 + denovo55328 + denovo55341 + denovo55355 + denovo55371 + denovo55403 +
                                 denovo55425 + denovo55429 + denovo55498 + denovo55505 + denovo55522 + denovo55583 +
                                 denovo55607 + denovo55609 + denovo55767 + denovo55832 + denovo55845 + denovo55853 + denovo55891 + denovo55896 +
                                 denovo55930 + denovo55945 + denovo55952 + denovo55967 + denovo55995 + denovo56011 + denovo56050 + denovo56056 +
                                 denovo56087 + denovo56120 + denovo56147 + denovo56155 + denovo56156 + denovo56160 + denovo56174 + denovo56178 +
                                 denovo56215 + denovo56233 + denovo56251 + denovo56265 + denovo56344 + denovo56419 + denovo56444 + denovo56455 +
                                 denovo56457 + denovo56464 + denovo56467 + denovo56472 + denovo56492 + denovo56502 + denovo56504 + denovo56505 +
                                 denovo56521 + denovo56569 + denovo56587 + denovo56635 + denovo56639 + denovo56671 + denovo56709 + denovo56761 +
                                 denovo56763 + denovo56772 + denovo56778 + denovo56803 + denovo56830 + denovo56835 + denovo56853 + denovo56871 +
                                 denovo56894 + denovo56914 + denovo56920 + denovo56982 + denovo57017 + denovo57045 + denovo57138 +
                                 denovo57163 + denovo57166 + denovo57185 + denovo57207 + denovo57275 + denovo57300 + denovo57305 + denovo57353 +
                                 denovo57372 + denovo57421 + denovo57424 + denovo57437 + denovo57459 + denovo57476 + denovo57505 + denovo57514 +
                                 denovo57534 + denovo57570 + denovo57615 + denovo57699 + denovo57713 + denovo57759 + denovo57771 + denovo57838 +
                                 denovo57844 + denovo57868 + denovo57871 + denovo57875 + denovo57911 + denovo57914 + denovo57922 + denovo57972 +
                                 denovo58022 + denovo58034 + denovo58061 + denovo58083 + denovo58090 + denovo58111 + denovo58171 + denovo58245 +
                                 denovo58262 + denovo58328 + denovo58353 + denovo58378 + denovo58390 + denovo58520 + denovo58550 + denovo58552 +
                                 denovo58582 + denovo58600 + denovo58604 + denovo58635 + denovo58678 + denovo58681 + denovo58699 + denovo58709 +
                                 denovo58713 + denovo58717 + denovo58752 + denovo58784 + denovo58808 + denovo58829 + denovo58830 + denovo58848 +
                                 denovo58853 + denovo58856 + denovo58885 + denovo58948 + denovo58950 + denovo58976 + denovo58996 + denovo59030 +
                                 denovo59037 + denovo59060 + denovo59106 + denovo59209 + denovo59213 + denovo59238 + denovo59262 + denovo59276 +
                                 denovo59291 + denovo59294 + denovo59296 + denovo59302 + denovo59310 + denovo59335 + denovo59354 + denovo59370 +
                                 denovo59371 + denovo59372 + denovo59384 + denovo59392 + denovo59406 + denovo59415 + denovo59416 + denovo59430 +
                                 denovo59433 + denovo59437 + denovo59451 + denovo59460 + denovo59575 + denovo59578 + denovo59643 + denovo59658 +
                                 denovo59732 + denovo59744 + denovo59754 + denovo59755 + denovo59779 + denovo59817 + denovo59839 + denovo59844 +
                                 denovo59846 + denovo59872 + denovo59930 + denovo59960 + denovo59970 + denovo60009 + denovo60012 + denovo60068 +
                                 denovo60081 + denovo60094 + denovo60163 + denovo60199 + denovo60243 + denovo60309 + denovo60332 + denovo60358 +
                                 denovo60363 + denovo60379 + denovo60410 + denovo60469 + denovo60532 + denovo60574 + denovo60593 + denovo60650 +
                                 denovo60654 + denovo60692 + denovo60724 + denovo60739 + denovo60775 + denovo60794 + denovo60831 +
                                 denovo60878 + denovo60912 + denovo60917 + denovo60924 + denovo60946 + denovo60996 + denovo61066 + denovo61073 +
                                 denovo61094 + denovo61142 + denovo61145 + denovo61168 + denovo61218 + denovo61222 + denovo61249 + denovo61354 +
                                 denovo61369 + denovo61376 + denovo61378 + denovo61396 + denovo61402 + denovo61480 + denovo61498 + denovo61502 +
                                 denovo61535 + denovo61540 + denovo61563 + denovo61586 + denovo61635 + denovo61640 + denovo61657 + denovo61702 +
                                 denovo61722 + denovo61728 + denovo61733 + denovo61752 + denovo61787 + denovo61805 + denovo61820 + denovo61841 +
                                 denovo61842 + denovo61849 + denovo61912 + denovo61974 + denovo62011 + denovo62056 + denovo62127 + denovo62153 +
                                 denovo62166 + denovo62178 + denovo62238 + denovo62253 + denovo62282 + denovo62284 +
                                 denovo62303 + denovo62380 + denovo62433 + denovo62506 + denovo62540 + denovo62578 + denovo62579 + denovo62590 +
                                 denovo62604 + denovo62650 + denovo62670 + denovo62737 + denovo62783 + denovo62804 + denovo62883 + denovo62905 +
                                 denovo62943 + denovo62948 + denovo62961 + denovo63018 + denovo63037 + denovo63057 + denovo63082 + denovo63101 +
                                 denovo63126 + denovo63143 + denovo63152 + denovo63216 + denovo63228 + denovo63340 + denovo63341 + denovo63389 +
                                 denovo63415 + denovo63429 + denovo63434 + denovo63444 + denovo63482 + denovo63522 + denovo63539 + denovo63552 +
                                 denovo63560 + denovo63605 + denovo63698 + denovo63738 + denovo63757 + denovo63777 + denovo63857 +
                                 denovo63861 + denovo63863 + denovo63890 + denovo63899 + denovo63907 + denovo63913 + denovo63925 + denovo63973 +
                                 denovo63988 + denovo64010 + denovo64018 + denovo64026 + denovo64067 + denovo64081 + denovo64089 + denovo64093 +
                                 denovo64094 + denovo64121 + denovo64122 + denovo64131 + denovo64140 + denovo64141 + denovo64150 + denovo64208 +
                                 denovo64236 + denovo64237 + denovo64244 + denovo64260 + denovo64290 + denovo64309 + denovo64324 +
                                 denovo64359 + denovo64381 + denovo64430 + denovo64431 + denovo64455 + denovo64570 +
                                 denovo64608 + denovo64618 + denovo64631 + denovo64645 + denovo64747 + denovo64768 + denovo64769 + denovo64797 +
                                 denovo64804 + denovo64806 + denovo64839 + denovo64844 + denovo64852 + denovo64891 + denovo64906 + denovo64963 +
                                 denovo64989 + denovo65009 + denovo65010 + denovo65040 + denovo65065 + denovo65072 + denovo65088 + denovo65096 +
                                 denovo65113 + denovo65127 + denovo65151 + denovo65209 + denovo65233 + denovo65241 + denovo65258 + denovo65285 +
                                 denovo65327 + denovo65328 + denovo65338 + denovo65341 + denovo65354 + denovo65371 + denovo65398 + denovo65424 +
                                 denovo65476 + denovo65477 + denovo65563 + denovo65566 + denovo65638 + denovo65660 + denovo65708 + denovo65753 +
                                 denovo65759 + denovo65795 + denovo65830 + denovo65848 + denovo65855 + denovo65870 + denovo65873 + denovo65901 +
                                 denovo66008 + denovo66081 + denovo66089 + denovo66094 + denovo66129 + denovo66136 + denovo66138 + denovo66177 +
                                 denovo66227 + denovo66234 + denovo66236 + denovo66318 + denovo66356 + denovo66372 + denovo66380 +
                                 denovo66381 + denovo66418 + denovo66424 + denovo66457 + denovo66467 + denovo66483 + denovo66521 + denovo66524 +
                                 denovo66566 + denovo66615 + denovo66677 + denovo66684 + denovo66689 + denovo66690 + denovo66691 + denovo66698 +
                                 denovo66732 + denovo66741 + denovo66744 + denovo66746 + denovo66749 + denovo66757 + denovo66765 + denovo66776 +
                                 denovo66805 + denovo66806 + denovo66836 + denovo66876 + denovo66881 + denovo66904 + denovo66906 +
                                 denovo66947 + denovo66982 + denovo66987 + denovo66991 + denovo66996 + denovo67022 + denovo67074 + denovo67086 +
                                 denovo67112 + denovo67138 + denovo67157 + denovo67163 + denovo67181 + denovo67191 + denovo67266 + denovo67297 +
                                 denovo67318 + denovo67322 + denovo67378 + denovo67406 + denovo67418 + denovo67450 + denovo67470 +
                                 denovo67489 + denovo67492 + denovo67500 + denovo67581 + denovo67583 + denovo67599 + denovo67612 + denovo67657 +
                                 denovo67660 + denovo67662 + denovo67696 + denovo67725 + denovo67728 + denovo67729 + denovo67762 + denovo67769 +
                                 denovo67770 + denovo67807 + denovo67834 + denovo67871 + denovo68056 + denovo68058 + denovo68107 +
                                 denovo68129 + denovo68136 + denovo68153 + denovo68160 + denovo68170 + denovo68181 + denovo68220 + denovo68227 +
                                 denovo68256 + denovo68285 + denovo68301 + denovo68326 + denovo68332 + denovo68397 + denovo68408 + denovo68409 +
                                 denovo68422 + denovo68432 + denovo68458 + denovo68488 + denovo68494 + denovo68544 + denovo68556 + denovo68577 +
                                 denovo68598 + denovo68640 + denovo68646 + denovo68709 + denovo68724 + denovo68735 + denovo68767 + denovo68791 +
                                 denovo68834 + denovo68839 + denovo68867 + denovo68889 + denovo68914 + denovo68955 + denovo68962 + denovo68980 +
                                 denovo69002 + denovo69004 + denovo69018 + denovo69030 + denovo69112 + denovo69113 + denovo69145 + denovo69181 +
                                 denovo69213 + denovo69301 + denovo69339 + denovo69380 + denovo69414 + denovo69441 + denovo69444 +
                                 denovo69451 + denovo69453 + denovo69455 + denovo69474 + denovo69526 + denovo69559 + denovo69576 +
                                 denovo69602 + denovo69617 + denovo69620 + denovo69623 + denovo69626 + denovo69665 + denovo69705 +
                                 denovo69711 + denovo69717 + denovo69730 + denovo69742 + denovo69770 + denovo69772 + denovo69773 + denovo69781 +
                                 denovo69839 + denovo69845 + denovo69857 + denovo69866 + denovo69909 + denovo69945 + denovo69972 + denovo69980 +
                                 denovo70028 + denovo70083 + denovo70094 + denovo70110 + denovo70238 + denovo70241 + denovo70253 + denovo70260 +
                                 denovo70264 + denovo70290 + denovo70326 + denovo70346 + denovo70405 + denovo70420 + denovo70434 + denovo70491 +
                                 denovo70519 + denovo70615 + denovo70651 + denovo70654 + denovo70671 + denovo70684 + denovo70701 + denovo70782 +
                                 denovo70806 + denovo70816 + denovo70845 + denovo70854 + denovo70857 + denovo70861 + denovo70900 + denovo70906 +
                                 denovo70946 + denovo70950 + denovo70968 + denovo70990 + denovo71006 + denovo71008 + denovo71014 + denovo71030 +
                                 denovo71050 + denovo71102 + denovo71142 + denovo71236 + denovo71278 + denovo71315 + denovo71333 + denovo71348 +
                                 denovo71359 + denovo71366 + denovo71372 + denovo71403 + denovo71426 + denovo71436 + denovo71438 + denovo71440 +
                                 denovo71441 + denovo71446 + denovo71449 + denovo71482 + denovo71516 + denovo71525 + denovo71547 + denovo71576 +
                                 denovo71594 + denovo71598 + denovo71604 + denovo71649 + denovo71666 + denovo71672 + denovo71674 + denovo71677 +
                                 denovo71810 + denovo71840 + denovo71843 + denovo71913 + denovo71953 + denovo71958 + denovo72053 +
                                 denovo72109 + denovo72153 + denovo72161 + denovo72166 + denovo72174 + denovo72199 + denovo72204 +
                                 denovo72210 + denovo72233 + denovo72284 + denovo72320 + denovo72332 + denovo72457 + denovo72462 + denovo72466 +
                                 denovo72474 + denovo72485 + denovo72496 + denovo72517 + denovo72554 + denovo72619 + denovo72634 + denovo72695 +
                                 denovo72739 + denovo72752 + denovo72767 + denovo72820 + denovo72826 + denovo72870 + denovo72892 + denovo72940 +
                                 denovo72952 + denovo72971 + denovo73038 + denovo73045 + denovo73105 + denovo73190 + denovo73206 + denovo73214 +
                                 denovo73217 + denovo73347 + denovo73349 + denovo73400 + denovo73404 + denovo73412 + denovo73422 + denovo73448 +
                                 denovo73451 + denovo73457 + denovo73461 + denovo73581 + denovo73630 + denovo73674 + denovo73710 + denovo73713 +
                                 denovo73743 + denovo73798 + denovo73816 + denovo73853 + denovo73890 + denovo73902 + denovo73915 + denovo73924 +
                                 denovo73927 + denovo73955 + denovo73973 + denovo74010 + denovo74063 + denovo74064 + denovo74075 +
                                 denovo74082 + denovo74097 + denovo74129 + denovo74165 + denovo74194 + denovo74276 + denovo74285 + denovo74297 +
                                 denovo74301 + denovo74336 + denovo74435 + denovo74440 + denovo74444 + denovo74478 + denovo74493 + denovo74494 +
                                 denovo74516 + denovo74523 + denovo74525 + denovo74564 + denovo74585 + denovo74676 + denovo74692 + denovo74702 +
                                 denovo74703 + denovo74781 + denovo74829 + denovo74852 + denovo74929 + denovo74938 + denovo74948 + denovo74955 +
                                 denovo74957 + denovo74963 + denovo74966 + denovo75041 + denovo75077 + denovo75116 + denovo75125 + denovo75134 +
                                 denovo75141 + denovo75196 + denovo75213 + denovo75224 + denovo75239 + denovo75252 + denovo75273 + denovo75353 +
                                 denovo75373 + denovo75444 + denovo75545 + denovo75585 + denovo75648 + denovo75653 + denovo75659 +
                                 denovo75670 + denovo75694 + denovo75712 + denovo75717 + denovo75718 + denovo75779 + denovo75796 + denovo75851 +
                                 denovo75859 + denovo75891 + denovo75906 + denovo75910 + denovo75947 + denovo76024 + denovo76040 + denovo76046 +
                                 denovo76064 + denovo76088 + denovo76094 + denovo76127 + denovo76133 + denovo76135 + denovo76139 + denovo76188 +
                                 denovo76216 + denovo76330 + denovo76333 + denovo76340 + denovo76345 + denovo76399 + denovo76428 +
                                 denovo76437 + denovo76440 + denovo76441 + denovo76469 + denovo76479 + denovo76511 + denovo76580 + denovo76594 +
                                 denovo76603 + denovo76605 + denovo76735 + denovo76762 + denovo76812 + denovo76824 + denovo76837 + denovo76861 +
                                 denovo76924 + denovo76940 + denovo76966 + denovo76989 + denovo77043 + denovo77089 + denovo77115 + denovo77134 +
                                 denovo77166 + denovo77215 + denovo77353 + denovo77362 + denovo77377 + denovo77381 + denovo77392 + denovo77399 +
                                 denovo77415 + denovo77462 + denovo77504 + denovo77524 + denovo77540 + denovo77590 + denovo77626 +
                                 denovo77657 + denovo77668 + denovo77675 + denovo77728 + denovo77743 + denovo77786 + denovo77811 + denovo77834 +
                                 denovo77836 + denovo77858 + denovo77859 + denovo77873 + denovo77882 + denovo77905 + denovo77917 + denovo77938 +
                                 denovo77945 + denovo78076 + denovo78079 + denovo78113 + denovo78120 + denovo78132 + denovo78134 +
                                 denovo78194 + denovo78285 + denovo78304 + denovo78331 + denovo78347 + denovo78352 + denovo78354 + denovo78360 +
                                 denovo78372 + denovo78388 + denovo78391 + denovo78415 + denovo78446 + denovo78480 + denovo78488 + denovo78504 +
                                 denovo78536 + denovo78542 + denovo78557 + denovo78559 + denovo78572 + denovo78599 + denovo78671 + denovo78707 +
                                 denovo78708 + denovo78765 + denovo78794 + denovo78818 + denovo78873 + denovo78893 + denovo78895 + denovo78897 +
                                 denovo78963 + denovo78976 + denovo78985 + denovo79050 + denovo79057 + denovo79092 + denovo79121 + denovo79153 +
                                 denovo79198 + denovo79238 + denovo79254 + denovo79312 + denovo79320 + denovo79348 + denovo79390 + denovo79459 +
                                 denovo79494 + denovo79502 + denovo79525 + denovo79629 + denovo79677 + denovo79683 + denovo79693 + denovo79731 +
                                 denovo79749 + denovo79757 + denovo79766 + denovo79797 + denovo79836 + denovo79865 + denovo79899 +
                                 denovo79905 + denovo79969 + denovo80005 + denovo80009 + denovo80018 + denovo80024 + denovo80046 + denovo80055 +
                                 denovo80059 + denovo80061 + denovo80143 + denovo80178 + denovo80188 + denovo80197 + denovo80238 + denovo80257 +
                                 denovo80271 + denovo80301 + denovo80335 + denovo80358 + denovo80386 + denovo80398 + denovo80411 + denovo80430 +
                                 denovo80439 + denovo80492 + denovo80501 + denovo80535 + denovo80539 + denovo80582 + denovo80588 + denovo80592 +
                                 denovo80658 + denovo80660 + denovo80661 + denovo80729 + denovo80740 + denovo80820 + denovo80836 + denovo80839 +
                                 denovo80877 + denovo80889 + denovo80893 + denovo80967 + denovo80987 + denovo81100 + denovo81109 + denovo81123 +
                                 denovo81188 + denovo81191 + denovo81194 + denovo81213 + denovo81232 + denovo81312 + denovo81313 + denovo81316 +
								 denovo81356 + denovo81406 + denovo81439 + denovo81567 + denovo81571 + denovo81603 + denovo81609 + denovo81613 +
                                 denovo81629 + denovo81666 + denovo81723 + denovo81754 + denovo81769 + denovo81843 + denovo81850 + denovo81870 +
                                 denovo81879 + denovo81920 + denovo81938 + denovo81974 + denovo81977 + denovo81982 + denovo82076 + denovo82092 +
                                 denovo82128 + denovo82139 + denovo82141 + denovo82154 + denovo82155 + denovo82156 + denovo82186 + denovo82209 +
                                 denovo82234 + denovo82317 + denovo82318 + denovo82322 + denovo82434 + denovo82462 + denovo82545 +
                                 denovo82600 + denovo82604 + denovo82633 + denovo82639 + denovo82649 + denovo82691 + denovo82694 + denovo82700 +
                                 denovo82747 + denovo82754 + denovo82775 + denovo82789 + denovo82791 + denovo82864 + denovo82875 + denovo82880 +
                                 denovo82884 + denovo82930 + denovo82993 + denovo82995 + denovo82997 + denovo82999 + denovo83002 + denovo83005 +
                                 denovo83021 + denovo83062 + denovo83119 + denovo83120 + denovo83131 + denovo83220 + denovo83240 + denovo83256 +
                                 denovo83284 + denovo83290 + denovo83314 + denovo83458 + denovo83479 + denovo83506 +
                                 denovo83513 + denovo83531 + denovo83545 + denovo83559 + denovo83574 + denovo83575 + denovo83588 + denovo83594 +
                                 denovo83605 + denovo83622 + denovo83625 + denovo83669 + denovo83694 + denovo83720 + denovo83725 + denovo83771 +
                                 denovo83799 + denovo83822 + denovo83855 + denovo83886 + denovo83913 + denovo83945 + denovo83991 + denovo84002 +
                                 denovo84020 + denovo84076 + denovo84125 + denovo84184 + denovo84224 + denovo84230 + denovo84304 + denovo84346 +
                                 denovo84364 + denovo84447 + denovo84463 + denovo84465 + denovo84489 + denovo84522 + denovo84545 + denovo84557 +
                                 denovo84620 + denovo84643 + denovo84660 + denovo84664 + denovo84665 + denovo84750 + denovo84768 + denovo84854 +
                                 denovo84860 + denovo84928 + denovo85016 + denovo85017 + denovo85043 + denovo85049 + denovo85080 + denovo85124 +
                                 denovo85145 + denovo85188 + denovo85191 + denovo85216 + denovo85235 + denovo85236 + denovo85241 + denovo85250 +
                                 denovo85253 + denovo85293 + denovo85297 + denovo85313 + denovo85326 + denovo85331 + denovo85404 + denovo85457 +
                                 denovo85481 + denovo85510 + denovo85514 + denovo85534 + denovo85537 + denovo85545 + denovo85567 + denovo85574 +
                                 denovo85616 + denovo85622 + denovo85650 + denovo85670 + denovo85722 + denovo85798 + denovo85851 + denovo85853 +
                                 denovo85883 + denovo85933 + denovo85937 + denovo85965 + denovo85966 + denovo85979 + denovo85993 + denovo85994 +
                                 denovo86079 + denovo86086 + denovo86135 + denovo86151 + denovo86182 + denovo86205 + denovo86212 + denovo86221 +
                                 denovo86254 + denovo86292 + denovo86319 + denovo86330 + denovo86337 + denovo86340 + denovo86350 + denovo86362 +
                                 denovo86419 + denovo86426 + denovo86531 + denovo86614 + denovo86616 + denovo86637 + denovo86655 + denovo86659 +
                                 denovo86680 + denovo86698 + denovo86712 + denovo86717 + denovo86733 + denovo86758 + denovo86767 + denovo86784 +
                                 denovo86797 + denovo86798 + denovo86853 + denovo86867 + denovo86893 + denovo86896 + denovo86911 + denovo86916 +
                                 denovo86964 + denovo86967 + denovo86991 + denovo87003 + denovo87007 + denovo87036 + denovo87064 + denovo87071 +
                                 denovo87190 + denovo87226 + denovo87257 + denovo87320 + denovo87382 + denovo87400 + denovo87447 + denovo87452 +
                                 denovo87483 + denovo87503 + denovo87516 + denovo87524 + denovo87533 + denovo87539 + denovo87551 + denovo87593 +
                                 denovo87671 + denovo87673 + denovo87675 + denovo87713 + denovo87727 + denovo87760 + denovo87761 + denovo87812 +
                                 denovo87895 + denovo87962 + denovo87981 + denovo88069 + denovo88077 + denovo88083 + denovo88101 + denovo88115 +
                                 denovo88163 + denovo88212 + denovo88216 + denovo88246 + denovo88285 + denovo88334 + denovo88346 +
                                 denovo88366 + denovo88404 + denovo88419 + denovo88449 + denovo88466 + denovo88491 + denovo88500 + denovo88507 +
                                 denovo88612 + denovo88624 + denovo88679 + denovo88688 + denovo88710 + denovo88730 + denovo88771 + denovo88799 +
                                 denovo88817 + denovo88837 + denovo88911 + denovo88934 + denovo89012 + denovo89048 + denovo89056 + denovo89069 +
                                 denovo89074 + denovo89077 + denovo89096 + denovo89117 + denovo89163 + denovo89168 + denovo89172 +
                                 denovo89187 + denovo89192 + denovo89198 + denovo89227 + denovo89231 + denovo89263 + denovo89272 + denovo89297 +
                                 denovo89309 + denovo89343 + denovo89370 + denovo89386 + denovo89417 + denovo89431 + denovo89444 + denovo89455 +
                                 denovo89457 + denovo89459 + denovo89513 + denovo89561 + denovo89574 + denovo89577 + denovo89598 + denovo89618 +
                                 denovo89648 + denovo89662 + denovo89671 + denovo89798 + denovo89803 + denovo89808 + denovo89828 + denovo89846 +
                                 denovo89870 + denovo89901 + denovo89905 + denovo89939 + denovo89975 + denovo90065 + denovo90084 + denovo90115 +
                                 denovo90159 + denovo90175 + denovo90194 + denovo90211 + denovo90264 + denovo90277 + denovo90319 + denovo90344 +
                                 denovo90353 + denovo90383 + denovo90404 + denovo90414 + denovo90475 + denovo90509 + denovo90530 + denovo90551 +
                                 denovo90553 + denovo90563 + denovo90564 + denovo90689 + denovo90693 + denovo90696 + denovo90707 + denovo90710 +
                                 denovo90713 + denovo90721 + denovo90734 + denovo90796 + denovo90806 + denovo90823 + denovo90853 + denovo90890 +
                                 denovo90930 + denovo91013 + denovo91033 + denovo91090 + denovo91095 + denovo91147 + denovo91161 + denovo91186 +
                                 denovo91205 + denovo91206 + denovo91214 + denovo91227 + denovo91228 + denovo91302 + denovo91346 +
                                 denovo91353 + denovo91357 + denovo91389 + denovo91393 + denovo91435 + denovo91466 + denovo91482 + denovo91498 +
                                 denovo91531 + denovo91568 + denovo91582 + denovo91608 + denovo91623 + denovo91671 + denovo91676 + denovo91677 +
                                 denovo91707 + denovo91729 + denovo91743 + denovo91747 + denovo91773 + denovo91779 + denovo91799 + denovo91808 +
                                 denovo91813 + denovo91826 + denovo91840 + denovo91945 + denovo91951 + denovo91956 + denovo91980 + denovo92003 +
                                 denovo92017 + denovo92024 + denovo92025 + denovo92044 + denovo92046 + denovo92077 + denovo92124 + denovo92132 +
                                 denovo92169 + denovo92172 + denovo92197 + denovo92214 + denovo92235 + denovo92237 + denovo92277 + denovo92284 +
                                 denovo92302 + denovo92327 + denovo92387 + denovo92401 + denovo92407 + denovo92478 + denovo92490 +
                                 denovo92507 + denovo92512 + denovo92588 + denovo92589 + denovo92607 + denovo92616 + denovo92622 + denovo92625 +
                                 denovo92634 + denovo92657 + denovo92677 + denovo92705 + denovo92731 + denovo92773 +
                                 denovo92809 + denovo92875 + denovo92885 + denovo92918 + denovo92978 + denovo92986 + denovo93012 + denovo93018 +
                                 denovo93057 + denovo93072 + denovo93106 + denovo93146 + denovo93156 + denovo93165 + denovo93183 + denovo93192 +
                                 denovo93221 + denovo93256 + denovo93292 + denovo93321 + denovo93346 + denovo93350 + denovo93355 + denovo93363 +
                                 denovo93376 + denovo93391 + denovo93405 + denovo93435 + denovo93446 + denovo93473 + denovo93480 +
                                 denovo93488 + denovo93730 + denovo93778 + denovo93786 + denovo93824 + denovo93862 + denovo93915 +
                                 denovo93918 + denovo93932 + denovo93935 + denovo93936 + denovo93957 + denovo93985 + denovo94010 + denovo94120 +
                                 denovo94141 + denovo94146 + denovo94151 + denovo94152 + denovo94170 + denovo94185 + denovo94209 + denovo94213 +
                                 denovo94217 + denovo94230 + denovo94306 + denovo94324 + denovo94332 + denovo94335 + denovo94336 + denovo94402 +
                                 denovo94530 + denovo94630 + denovo94632 + denovo94644 + denovo94715 + denovo94741 + denovo94761 + denovo94762 +
                                 denovo94784 + denovo94804 + denovo94840 + denovo94848 + denovo94862 + denovo94880 + denovo94899 + denovo94928 +
                                 denovo94978 + denovo94988 + denovo95005 + denovo95033 + denovo95057 + denovo95084 + denovo95087 + denovo95115 +
                                 denovo95131 + denovo95148 + denovo95157 + denovo95167 + denovo95172 + denovo95181 + denovo95233 + denovo95292 +
                                 denovo95302 + denovo95315 + denovo95370 + denovo95387 + denovo95404 + denovo95413 + denovo95438 +
                                 denovo95440 + denovo95463 + denovo95496 + denovo95521 + denovo95527 + denovo95588 + denovo95595 + denovo95596 +
                                 denovo95604 + denovo95614 + denovo95627 + denovo95648 + denovo95676 + denovo95763 + denovo95771 + denovo95820 +
                                 denovo95856 + denovo95900 + denovo95910 + denovo96003 + denovo96018 + denovo96022 + denovo96128 + denovo96134 +
                                 denovo96174 + denovo96188 + denovo96263 + denovo96302 + denovo96320 + denovo96338 + denovo96406 + denovo96433 +
                                 denovo96470 + denovo96497 + denovo96506 + denovo96532 + denovo96534 + denovo96654 + denovo96659 + denovo96677 +
                                 denovo96688 + denovo96702 + denovo96712 + denovo96719 + denovo96769 + denovo96816 + denovo96855 + denovo96940 +
                                 denovo96943 + denovo96948 + denovo96958 + denovo96966 + denovo96971 + denovo97027 + denovo97057 + denovo97076 +
                                 denovo97093 + denovo97106 + denovo97151 + denovo97154 + denovo97185 + denovo97204 + denovo97267 + denovo97334 +
                                 denovo97372 + denovo97396 + denovo97406 + denovo97410 + denovo97414 + denovo97432 + denovo97452 + denovo97510 +
                                 denovo97542 + denovo97577 + denovo97611 + denovo97618 + denovo97657 + denovo97686 + denovo97700 + denovo97710 +
                                 denovo97713 + denovo97721 + denovo97749 + denovo97803 + denovo97807 + denovo97844 + denovo97849 + denovo97874 +
                                 denovo97886 + denovo97913 + denovo97931 + denovo97954 + denovo98038 + denovo98083 + denovo98129 + denovo98155 +
                                 denovo98175 + denovo98183 + denovo98189 + denovo98235 + denovo98249 + denovo98257 + denovo98285 + denovo98308 +
                                 denovo98352 + denovo98363 + denovo98382 + denovo98428 + denovo98433 + denovo98554 + denovo98636 +
                                 denovo98661 + denovo98663 + denovo98698 + denovo98739 + denovo98749 + denovo98760 + denovo98762 +
                                 denovo98779 + denovo98809 + denovo98841 + denovo98886 + denovo98892 + denovo98898 + denovo98900 + denovo98907 +
                                 denovo98935 + denovo98942 + denovo98967 + denovo98992 + denovo99023 + denovo99074 + denovo99109 + denovo99118 +
                                 denovo99127 + denovo99133 + denovo99144 + denovo99173 + denovo99257 + denovo99281 + denovo99399 + denovo99450 +
                                 denovo99451 + denovo99452 + denovo99455 + denovo99559 + denovo99561 + denovo99563 + denovo99589 + denovo99606 +
                                 denovo99690 + denovo99736 + denovo99758 + denovo99780 + denovo99785 + denovo99854 + denovo99876 + denovo99896 +
                                 denovo99927 + denovo99929 + denovo99965 + denovo99978 + denovo100004 + denovo100052 + denovo100080 + denovo100132 +
                                 denovo100139 + denovo100160 + denovo100176 + denovo100188 + denovo100208 + denovo100289 + denovo100310 +
                                 denovo100334 + denovo100360 + denovo100408 + denovo100440 + denovo100446 + denovo100473 + denovo100534 + denovo100538 +
                                 denovo100550 + denovo100563 + denovo100576 + denovo100613 + denovo100617 + denovo100662 + denovo100783 +
                                 denovo100889 + denovo100892 + denovo100894 + denovo100906 + denovo100923 + denovo100953 + denovo100992 + denovo101027 +
                                 denovo101034 + denovo101049 + denovo101069 + denovo101094 + denovo101105 + denovo101112 + denovo101119 + denovo101164 +
                                 denovo101173 + denovo101225 + denovo101239 + denovo101252 + denovo101288 + denovo101296 + denovo101357 +
                                 denovo101382 + denovo101405 + denovo101423 + denovo101444 + denovo101458 + denovo101472 + denovo101503 + denovo101505 +
                                 denovo101517 + denovo101525 + denovo101598 + denovo101640 + denovo101644 + denovo101669 + denovo101744 + denovo101750 +
                                 denovo101784 + denovo101803 + denovo101815 + denovo101821 + denovo101830 + denovo101876 + denovo101896 +
                                 denovo101941 + denovo101953 + denovo101964 + denovo102006 + denovo102014 + denovo102021 + denovo102030 + denovo102108 +
                                 denovo102111 + denovo102168 + denovo102228 + denovo102239 + denovo102244 + denovo102260 + denovo102262 + denovo102288 +
                                 denovo102303 + denovo102329 + denovo102335 + denovo102366 + denovo102427 + denovo102433 + denovo102468 + denovo102480 +
                                 denovo102494 + denovo102600 + denovo102604 + denovo102609 + denovo102619 + denovo102620 +
                                 denovo102633 + denovo102681 + denovo102682 + denovo102699 + denovo102728 + denovo102737 + denovo102789 +
                                 denovo102865 + denovo102866 + denovo102887 + denovo102905 + denovo102982 + denovo103023 + denovo103036 + denovo103060 +
                                 denovo103083 + denovo103108 + denovo103133 + denovo103135 + denovo103138 + denovo103153 + denovo103165 + denovo103168 +
                                 denovo103172 + denovo103181 + denovo103223 + denovo103241 + denovo103277 + denovo103357 + denovo103358 + denovo103425 +
                                 denovo103444 + denovo103454 + denovo103470 + denovo103507 + denovo103509 + denovo103516 + denovo103544 + denovo103547 +
                                 denovo103555 + denovo103557 + denovo103569 + denovo103586 + denovo103611 + denovo103617 + denovo103638 + denovo103667 +
                                 denovo103720 + denovo103745 + denovo103781 + denovo103838 + denovo103864 + denovo103914 + denovo104042 + denovo104044 +
                                 denovo104051 + denovo104066 + denovo104079 + denovo104093 + denovo104100 + denovo104104 + denovo104107 + denovo104115 +
                                 denovo104127 + denovo104143 + denovo104177 + denovo104194 + denovo104200 + denovo104210 + denovo104214 + denovo104238 +
                                 denovo104252 + denovo104283 + denovo104328 + denovo104329 + denovo104364 + denovo104397 + denovo104419 + denovo104421 +
                                 denovo104422 + denovo104437 + denovo104494 + denovo104562 + denovo104591 + denovo104618 + denovo104633 + denovo104661 +
                                 denovo104667 + denovo104699 + denovo104759 + denovo104792 + denovo104801 + denovo104836 + denovo104873 + denovo104881 +
                                 denovo104885 + denovo104930 + denovo104936 + denovo104938 + denovo104986 + denovo105076 + denovo105086 + denovo105088 +
                                 denovo105146 + denovo105158 + denovo105231 + denovo105241 + denovo105341 + denovo105345 + denovo105354 + denovo105361 +
                                 denovo105367 + denovo105396 + denovo105417 + denovo105445 + denovo105492 + denovo105509 + denovo105545 + denovo105588 +
                                 denovo105623 + denovo105658 + denovo105661 + denovo105731 + denovo105732 + denovo105737 + denovo105755 + denovo105783 +
                                 denovo105799 + denovo105801 + denovo105804 + denovo105823 + denovo105948 + denovo105954 + denovo105956 + denovo105962 +
                                 denovo105965 + denovo105982 + denovo105985 + denovo105990 + denovo106067 + denovo106077 + denovo106093 + denovo106145 +
                                 denovo106160 + denovo106224 + denovo106281 + denovo106288 + denovo106289 + denovo106299 + denovo106381 + denovo106384 +
                                 denovo106392 + denovo106393 + denovo106431 + denovo106492 + denovo106496 + denovo106521 + denovo106532 + denovo106566 +
                                 denovo106578 + denovo106585 + denovo106614 + denovo106620 + denovo106627 + denovo106663 + denovo106671 + denovo106770 +
                                 denovo106801 + denovo106802 + denovo106804 + denovo106840 + denovo106857 + denovo106868 + denovo106875 + denovo106895 +
                                 denovo106904 + denovo106950 + denovo106999 + denovo107014 + denovo107017 + denovo107072 + denovo107105 + denovo107106 +
                                 denovo107202 + denovo107220 + denovo107255 + denovo107294 + denovo107340 + denovo107374 + denovo107388 + denovo107397 +
                                 denovo107409 + denovo107431 + denovo107464 + denovo107470 + denovo107499 + denovo107500 + denovo107508 +
                                 denovo107548 + denovo107571 + denovo107586 + denovo107612 + denovo107629 + denovo107631 + denovo107651 + denovo107669 +
                                 denovo107686 + denovo107732 + denovo107752 + denovo107758 + denovo107873 + denovo107882 + denovo107921 + denovo107930 +
                                 denovo107970 + denovo108021 + denovo108047 + denovo108074 + denovo108096 + denovo108107 + denovo108223 + denovo108249 +
                                 denovo108251 + denovo108291 + denovo108316 + denovo108355 + denovo108432 + denovo108433 + denovo108457 + denovo108483 +
                                 denovo108502 + denovo108527 + denovo108534 + denovo108552 + denovo108572 + denovo108588 + denovo108596 + denovo108601 +
                                 denovo108623 + denovo108632 + denovo108646 + denovo108671 + denovo108704 + denovo108723 + denovo108743 + denovo108832 +
                                 denovo108846 + denovo108912 + denovo108921 + denovo108932 + denovo108951 + denovo109023 + denovo109025 + denovo109034 +
                                 denovo109037 + denovo109061 + denovo109071 + denovo109094 + denovo109104 + denovo109109 + denovo109203 + denovo109266 +
                                 denovo109301 + denovo109324 + denovo109347 + denovo109398 + denovo109524 + denovo109556 + denovo109659 + denovo109682 +
                                 denovo109722 + denovo109734 + denovo109751 + denovo109754 + denovo109758 + denovo109851 + denovo109853 + denovo109863 +
                                 denovo109877 + denovo109880 + denovo109899 + denovo109921 + denovo109944 + denovo109977 + denovo109982 + denovo109992 +
                                 denovo110026 + denovo110034 + denovo110036 + denovo110079 + denovo110093 + denovo110121 + denovo110141 + denovo110157 +
                                 denovo110173 + denovo110180 + denovo110228 + denovo110234 + denovo110235 + denovo110248 + denovo110300 + denovo110312 +
                                 denovo110357 + denovo110371 + denovo110418 + denovo110483 + denovo110521 + denovo110533 + denovo110550 + denovo110586 +
                                 denovo110608 + denovo110610 + denovo110651 + denovo110699 + denovo110705 + denovo110713 + denovo110716 +
                                 denovo110721 + denovo110764 + denovo110798 + denovo110817 + denovo110826 + denovo110832 + denovo110881 + denovo110902 +
                                 denovo110932 + denovo110945 + denovo110950 + denovo110966 + denovo110969 + denovo110988 + denovo111002 + denovo111028 +
                                 denovo111048 + denovo111074 + denovo111078 + denovo111110 + denovo111126 + denovo111145 + denovo111146 +
                                 denovo111159 + denovo111166 + denovo111173 + denovo111197 + denovo111262 + denovo111274 + denovo111294 + denovo111300 +
                                 denovo111343 + denovo111347 + denovo111370 + denovo111376 + denovo111383 + denovo111411 + denovo111415 + denovo111497 +
                                 denovo111506 + denovo111508 + denovo111516 + denovo111528 + denovo111565 + denovo111596 + denovo111648 + denovo111664 +
                                 denovo111666 + denovo111669 + denovo111674 + denovo111728 + denovo111730 + denovo111764 + denovo111786 + denovo111817 +
                                 denovo111841 + denovo111851 + denovo111906 + denovo111939 + denovo111950 + denovo111960 + denovo111967 + denovo112075 +
                                 denovo112108 + denovo112121 + denovo112146 + denovo112211 + denovo112235 + denovo112262 + denovo112278 + denovo112317 +
                                 denovo112321 + denovo112322 + denovo112331 + denovo112345 + denovo112346 + denovo112431 + denovo112451 +
                                 denovo112490 + denovo112507 + denovo112552 + denovo112560 + denovo112590 + denovo112593 + denovo112615 +
                                 denovo112655 + denovo112706 + denovo112715 + denovo112728 + denovo112813 + denovo112818 +
                                 denovo112829 + denovo112837 + denovo112868 + denovo112884 + denovo112891 + denovo112893 + denovo112896 + denovo112903 +
                                 denovo112918 + denovo112970 + denovo112982 + denovo112987 + denovo113033 + denovo113077 + denovo113103 + denovo113109 +
                                 denovo113144 + denovo113176 + denovo113237 + denovo113257 + denovo113275 + denovo113279 + denovo113284 + denovo113304 +
                                 denovo113316 + denovo113345 + denovo113404 + denovo113413 + denovo113416 + denovo113425 + denovo113428 +
                                 denovo113452 + denovo113471 + denovo113502 + denovo113545 + denovo113567 + denovo113572 + denovo113592 +
                                 denovo113598 + denovo113616 + denovo113634 + denovo113638 + denovo113680 + denovo113685 + denovo113689 + denovo113698 +
                                 denovo113757 + denovo113760 + denovo113769 + denovo113770 + denovo113780 + denovo113786 + denovo113790 + denovo113791 +
                                 denovo113808 + denovo113815 + denovo113831 + denovo113859 + denovo113877 + denovo113881 + denovo113885 + denovo113925 +
                                 denovo113930 + denovo113961 + denovo113996 + denovo114061 + denovo114081 + denovo114191 + denovo114192 + denovo114215 +
                                 denovo114260 + denovo114268 + denovo114271 + denovo114290 + denovo114302 + denovo114316 + denovo114346 + denovo114370 +
                                 denovo114378 + denovo114380 + denovo114397 + denovo114413 + denovo114423 + denovo114461 + denovo114488 +
                                 denovo114508 + denovo114520 + denovo114527 + denovo114579 + denovo114593 + denovo114604 + denovo114678 +
                                 denovo114709 + denovo114712 + denovo114762 + denovo114780 + denovo114786 + denovo114821 + denovo114868 + denovo114869 +
                                 denovo114898 + denovo114908 + denovo114920 + denovo114926 + denovo114977 + denovo115038 + denovo115043 + denovo115090 +
                                 denovo115118 + denovo115175 + denovo115212 + denovo115219 + denovo115251 + denovo115298 + denovo115300 + denovo115302 +
                                 denovo115358 + denovo115362 + denovo115374 + denovo115377 + denovo115394 + denovo115399 + denovo115415 + denovo115421 +
                                 denovo115448 + denovo115451 + denovo115472 + denovo115474 + denovo115511 + denovo115518 + denovo115532 + denovo115536 +
                                 denovo115551 + denovo115562 + denovo115588 + denovo115624 + denovo115629 + denovo115665 + denovo115669 + denovo115695 +
                                 denovo115708 + denovo115744 + denovo115792 + denovo115803 + denovo115864 + denovo115973 + denovo115983 +
                                 denovo116020 + denovo116045 + denovo116056 + denovo116061 + denovo116071 + denovo116092 + denovo116096 +
                                 denovo116145 + denovo116156 + denovo116191 + denovo116193 + denovo116197 + denovo116223 + denovo116253 + denovo116277 +
                                 denovo116279 + denovo116317 + denovo116330 + denovo116341 + denovo116357 + denovo116374 + denovo116378 + denovo116454 +
                                 denovo116482 + denovo116517 + denovo116518 + denovo116544 + denovo116549 + denovo116578 + denovo116583 + denovo116619 +
                                 denovo116622 + denovo116642 + denovo116662 + denovo116670 + denovo116689 + denovo116718 + denovo116724 + denovo116728 +
                                 denovo116777 + denovo116784 + denovo116794 + denovo116796 + denovo116841 + denovo116843 + denovo116850 + denovo116866 +
                                 denovo116882 + denovo116891 + denovo116892 + denovo116938 + denovo116961 + denovo116966 + denovo116971 + denovo116989 +
                                 denovo116990 + denovo117001 + denovo117061 + denovo117112 + denovo117148 + denovo117157 + denovo117165 + denovo117201 +
                                 denovo117252 + denovo117256 + denovo117288 + denovo117320 + denovo117391 + denovo117394 + denovo117422 + denovo117462 +
                                 denovo117493 + denovo117517 + denovo117537 + denovo117558 + denovo117581 + denovo117599 + denovo117604 + denovo117633 +
                                 denovo117651 + denovo117652 + denovo117657 + denovo117660 + denovo117725 + denovo117742 + denovo117765 + denovo117777 +
                                 denovo117804 + denovo117812 + denovo117821 + denovo117831 + denovo117896 + denovo117904 + denovo117919 + denovo117923 +
                                 denovo117941 + denovo117962 + denovo117988 + denovo118061 + denovo118084 + denovo118146 + denovo118158 + denovo118174 +
                                 denovo118189 + denovo118205 + denovo118252 + denovo118271 + denovo118277 + denovo118303 + denovo118317 +
                                 denovo118321 + denovo118325 + denovo118355 + denovo118364 + denovo118412 + denovo118432 + denovo118460 + denovo118477 +
                                 denovo118486 + denovo118506 + denovo118525 + denovo118604 + denovo118636 + denovo118647 + denovo118670 + denovo118714 +
                                 denovo118741 + denovo118771 + denovo118799 + denovo118807 + denovo118820 + denovo118851 + denovo118877 + denovo118880 +
                                 denovo118886 + denovo118916 + denovo118923 + denovo118954 + denovo118957 + denovo118967 + denovo118968 + denovo118988 +
                                 denovo118989 + denovo118994 + denovo119047 + denovo119059 + denovo119068 + denovo119072 + denovo119086 + denovo119104 +
                                 denovo119115 + denovo119241 + denovo119250 + denovo119278 + denovo119320 + denovo119333 + denovo119341 + denovo119342 +
                                 denovo119366 + denovo119394 + denovo119466 + denovo119474 + denovo119546 + denovo119567 + denovo119590 + denovo119613 +
                                 denovo119640 + denovo119671 + denovo119678 + denovo119680 + denovo119700 + denovo119730 + denovo119742 + denovo119758 +
                                 denovo119769 + denovo119779 + denovo119793 + denovo119826 + denovo119869 + denovo119883 +
                                 denovo119958 + denovo119980 + denovo120005 + denovo120009 + denovo120013 + denovo120057 + denovo120068 + denovo120072 +
                                 denovo120076 + denovo120082 + denovo120109 + denovo120140 + denovo120152 + denovo120159 + denovo120210 + denovo120266 +
                                 denovo120269 + denovo120294 + denovo120309 + denovo120316 + denovo120358 + denovo120361 + denovo120363 + denovo120379 +
                                 denovo120401 + denovo120469 + denovo120495 + denovo120499 + denovo120503 + denovo120511 + denovo120531 + denovo120544 +
                                 denovo120590 + denovo120624 + denovo120636 + denovo120642 + denovo120653 + denovo120672 + denovo120685 + denovo120695 +
                                 denovo120706 + denovo120710 + denovo120779 + denovo120786 + denovo120790 + denovo120799 + denovo120828 + denovo120840 +
                                 denovo120845 + denovo120848 + denovo120858 + denovo120901 + denovo121027 + denovo121064 + denovo121110 + denovo121116 +
                                 denovo121120 + denovo121137 + denovo121219 + denovo121223 + denovo121269 + denovo121275 + denovo121294 + denovo121308 +
                                 denovo121312 + denovo121381 + denovo121384 + denovo121386 + denovo121401 + denovo121420 + denovo121433 + denovo121436 +
                                 denovo121441 + denovo121478 + denovo121489 + denovo121516 + denovo121519 + denovo121523 + denovo121551 + denovo121560 +
                                 denovo121570 + denovo121640 + denovo121675 + denovo121684 + denovo121721 + denovo121782 + denovo121785 + denovo121800 +
                                 denovo121801 + denovo121827 + denovo121860 + denovo121925 + denovo121970 + denovo121990 +
                                 denovo121999 + denovo122006 + denovo122072 + denovo122092 + denovo122096 + denovo122113 + denovo122132 + denovo122165 +
                                 denovo122168 + denovo122178 + denovo122213 + denovo122224 + denovo122259 + denovo122298 + denovo122316 + denovo122344 +
                                 denovo122378 + denovo122380 + denovo122414 + denovo122480 + denovo122490 + denovo122491 + denovo122513 + denovo122518 +
                                 denovo122520 + denovo122548 + denovo122553 + denovo122556 + denovo122565 + denovo122568 + denovo122573 + denovo122583 +
                                 denovo122592 + denovo122615 + denovo122639 + denovo122645 + denovo122666 + denovo122667 + denovo122746 + denovo122780 +
                                 denovo122797 + denovo122798 + denovo122845 + denovo122906 + denovo122921 + denovo122925 + denovo122930 + denovo122954 +
                                 denovo122960 + denovo122977 + denovo123001 + denovo123074 + denovo123077 + denovo123115 + denovo123118 + denovo123119 +
                                 denovo123145 + denovo123157 + denovo123162 + denovo123170 + denovo123214 + denovo123215 + denovo123261 + denovo123303 +
                                 denovo123321 + denovo123360 + denovo123447 + denovo123485 + denovo123509 + denovo123519 + denovo123533 +
                                 denovo123597 + denovo123620 + denovo123646 + denovo123648 + denovo123671 + denovo123675 + denovo123776 +
                                 denovo123806 + denovo123818 + denovo123827 + denovo123853 + denovo123874 + denovo123882 + denovo123884 + denovo123895 +
                                 denovo123914 + denovo123917 + denovo123967 + denovo123979 + denovo124066 + denovo124070 + denovo124078 + denovo124083 +
                                 denovo124096 + denovo124105 + denovo124106 + denovo124111 + denovo124135 + denovo124158 + denovo124251 + denovo124272 +
                                 denovo124332 + denovo124333 + denovo124360 + denovo124390 + denovo124413 + denovo124434 + denovo124468 + denovo124493 +
                                 denovo124512 + denovo124551 + denovo124552 + denovo124563 + denovo124573 + denovo124591 + denovo124611 + denovo124623 +
                                 denovo124626 + denovo124632 + denovo124639 + denovo124646 + denovo124648 + denovo124676 + denovo124682 +
                                 denovo124700 + denovo124701 + denovo124703 + denovo124711 + denovo124741 + denovo124750 + denovo124752 + denovo124772 +
                                 denovo124782 + denovo124817 + denovo124822 + denovo124839 + denovo124845 + denovo124858 + denovo124864 + denovo124877 +
                                 denovo124988 + denovo125002 + denovo125014 + denovo125074 + denovo125080 + denovo125130 + denovo125139 + denovo125140 +
                                 denovo125155 + denovo125170 + denovo125179 + denovo125207 + denovo125221 + denovo125244 + denovo125249 + denovo125251 +
                                 denovo125290 + denovo125335 + denovo125365 + denovo125458 + denovo125499 + denovo125550 + denovo125592 + denovo125594 +
                                 denovo125635 + denovo125689 + denovo125706 + denovo125732 + denovo125737 + denovo125744 + denovo125823 + denovo125854 +
                                 denovo125870 + denovo125883 + denovo125923 + denovo125977 + denovo125998 + denovo126003 + denovo126029 + denovo126072 +
                                 denovo126076 + denovo126091 + denovo126119 + denovo126148 + denovo126152 + denovo126178 + denovo126196 + denovo126200 +
                                 denovo126207 + denovo126209 + denovo126212 + denovo126226 + denovo126280 + denovo126357 + denovo126363 + denovo126378 +
                                 denovo126387 + denovo126396 + denovo126408 + denovo126415 + denovo126468 + denovo126472 + denovo126495 + denovo126510 +
                                 denovo126572 + denovo126600 + denovo126659 + denovo126671 + denovo126677 + denovo126699 + denovo126700 + denovo126726 +
                                 denovo126729 + denovo126741 + denovo126752 + denovo126781 + denovo126799 + denovo126852 + denovo126881 + denovo126891 +
                                 denovo126911 + denovo126980 + denovo127002 + denovo127054 + denovo127098 + denovo127118 + denovo127141 + denovo127143 +
                                 denovo127146 + denovo127164 + denovo127171 + denovo127187 + denovo127189 + denovo127198 + denovo127202 + denovo127219 +
                                 denovo127233 + denovo127250 + denovo127251 + denovo127357 + denovo127415 + denovo127440 + denovo127471 + denovo127477 +
                                 denovo127482 + denovo127499 + denovo127521 + denovo127568 + denovo127573 + denovo127575 + denovo127637 + denovo127663 +
                                 denovo127677 + denovo127692 + denovo127701 + denovo127721 + denovo127733 + denovo127758 + denovo127817 + denovo127851 +
                                 denovo127873 + denovo127900 + denovo127945 + denovo127955 + denovo127957 + denovo127966 + denovo127971 + denovo127972 +
                                 denovo127975 + denovo127980 + denovo127992 + denovo127998 + denovo128053 + denovo128071 + denovo128101 + denovo128107 +
                                 denovo128116 + denovo128125 + denovo128140 + denovo128175 + denovo128178 + denovo128221 + denovo128223 + denovo128250 +
                                 denovo128263 + denovo128293 + denovo128315 + denovo128325 + denovo128354 + denovo128371 + denovo128391 + denovo128398 +
                                 denovo128401 + denovo128462 + denovo128463 + denovo128476 + denovo128515 + denovo128557 + denovo128589 +
                                 denovo128592 + denovo128622 + denovo128623 + denovo128634 + denovo128644 + denovo128672 + denovo128675 + denovo128701 +
                                 denovo128736 + denovo128747 + denovo128762 + denovo128764 + denovo128783 + denovo128796 + denovo128814 + denovo128816 +
                                 denovo128829 + denovo128838 + denovo128845 + denovo128863 + denovo128893 + denovo128909 + denovo128917 + denovo128920 +
                                 denovo128936 + denovo128964 + denovo128968 + denovo128977 + denovo129007 + denovo129013 + denovo129141 + denovo129152 +
                                 denovo129177 + denovo129196 + denovo129198 + denovo129235 + denovo129238 + denovo129240 + denovo129254 + denovo129258 +
                                 denovo129277 + denovo129283 + denovo129284 + denovo129301 + denovo129310 + denovo129324 + denovo129333 + denovo129380 +
                                 denovo129384 + denovo129394 + denovo129399 + denovo129408 + denovo129428 + denovo129440 + denovo129441 + denovo129514 +
                                 denovo129543 + denovo129626 + denovo129692 + denovo129704 + denovo129782 + denovo129784 + denovo129803 + denovo129844 +
                                 denovo129848 + denovo129886 + denovo129904 + denovo129947 + denovo129952 + denovo129999 + denovo130002 + denovo130019 +
                                 denovo130022 + denovo130116 + denovo130120 + denovo130178 + denovo130234 + denovo130253 + denovo130254 +
                                 denovo130285 + denovo130302 + denovo130313 + denovo130341 + denovo130370 + denovo130377 + denovo130419 + denovo130437 +
                                 denovo130455 + denovo130482 + denovo130519 + denovo130524 + denovo130551 + denovo130565 + denovo130573 + denovo130614 +
                                 denovo130629 + denovo130658 + denovo130666 + denovo130722 + denovo130730 + denovo130742 + denovo130852 +
                                 denovo130858 + denovo130944 + denovo130962 + denovo130997 + denovo131009 + denovo131022 + denovo131053 + denovo131078 +
                                 denovo131087 + denovo131138 + denovo131144 + denovo131163 + denovo131166 + denovo131187 + denovo131208 + denovo131261 +
                                 denovo131302 + denovo131303 + denovo131309 + denovo131347 + denovo131348 + denovo131350 + denovo131358 + denovo131359 +
                                 denovo131404 + denovo131407 + denovo131411 + denovo131414 + denovo131426 + denovo131468 + denovo131531 + denovo131565 +
                                 denovo131570 + denovo131571 + denovo131574 + denovo131726 + denovo131779 + denovo131782 + denovo131820 + denovo131825 +
                                 denovo131832 + denovo131839 + denovo131843 + denovo131885 + denovo131932 + denovo131949 + denovo131965 + denovo131971 +
                                 denovo131973 + denovo132001 + denovo132004 + denovo132022 + denovo132024 + denovo132032 + denovo132095 + denovo132137 +
                                 denovo132140 + denovo132164 + denovo132196 + denovo132217 + denovo132263 + denovo132306 + denovo132322 + denovo132421 +
                                 denovo132425 + denovo132468 + denovo132469 + denovo132504 + denovo132511 + denovo132515 + denovo132535 + denovo132597 +
                                 denovo132714 + denovo132812 + denovo132839 + denovo132854 + denovo132897 + denovo132920 + denovo132929 + denovo132933 +
                                 denovo132986 + denovo132990 + denovo132993 + denovo133026 + denovo133027 + denovo133085 + denovo133221 + denovo133247 +
                                 denovo133293 + denovo133295 + denovo133319 + denovo133322 + denovo133331 + denovo133353 + denovo133358 +
                                 denovo133404 + denovo133446 + denovo133470 + denovo133520 + denovo133532 + denovo133542 + denovo133590 + denovo133599 +
                                 denovo133605 + denovo133632 + denovo133647 + denovo133665 + denovo133741 + denovo133746 + denovo133781 + denovo133792 +
                                 denovo133808 + denovo133858 + denovo133859 + denovo133872 + denovo133873 + denovo133891 + denovo133905 + denovo133916 +
                                 denovo133934 + denovo133939 + denovo133957 + denovo133971 + denovo134018 + denovo134052 + denovo134098 + denovo134109 +
                                 denovo134145 + denovo134274 + denovo134285 + denovo134292 + denovo134324 + denovo134332 + denovo134336 + denovo134409 +
                                 denovo134417 + denovo134436 + denovo134441 + denovo134444 + denovo134449 + denovo134469 + denovo134495 + denovo134517 +
                                 denovo134546 + denovo134579 + denovo134586 + denovo134629 + denovo134668 + denovo134672 + denovo134745 + denovo134782 +
                                 denovo134794 + denovo134813 + denovo134832 + denovo134893 + denovo134902 + denovo134907 + denovo134926 + denovo134972 +
                                 denovo135007 + denovo135013 + denovo135018 + denovo135063 + denovo135084 + denovo135105 + denovo135117 + denovo135133 +
                                 denovo135137 + denovo135158 + denovo135163 + denovo135166 + denovo135171 + denovo135221 + denovo135222 + denovo135226 +
                                 denovo135255 + denovo135261 + denovo135287 + denovo135298 + denovo135316 + denovo135330 + denovo135334 + denovo135335 +
                                 denovo135337 + denovo135388 + denovo135499 + denovo135508 + denovo135569 + denovo135601 + denovo135629 + denovo135693 +
                                 denovo135731 + denovo135790 + denovo135823 + denovo135837 + denovo135853 + denovo135855 + denovo135863 + denovo135884 +
                                 denovo135897 + denovo135943 + denovo135945 + denovo135968 + denovo135972 + denovo136019 + denovo136026 + denovo136067 +
                                 denovo136083 + denovo136087 + denovo136094 + denovo136097 + denovo136111 + denovo136164 + denovo136194 + denovo136196 +
                                 denovo136233 + denovo136274 + denovo136280 + denovo136286 + denovo136350 + denovo136375 +
                                 denovo136408 + denovo136446 + denovo136456 + denovo136545 + denovo136569 + denovo136588 + denovo136593 + denovo136598 +
                                 denovo136600 + denovo136623 + denovo136626 + denovo136633 + denovo136649 + denovo136655 + denovo136709 + denovo136826 +
                                 denovo136859 + denovo136861 + denovo136871 + denovo136920 + denovo136946 + denovo136990 + denovo137066 + denovo137149 +
                                 denovo137159 + denovo137221 + denovo137230 + denovo137292 + denovo137305 + denovo137318 + denovo137442 + denovo137473 +
                                 denovo137505 + denovo137514 + denovo137542 + denovo137549 + denovo137552 + denovo137553 + denovo137579 + denovo137590 +
                                 denovo137597 + denovo137627 + denovo137657 + denovo137664 + denovo137700 + denovo137738 + denovo137757 +
                                 denovo137793 + denovo137802 + denovo137845 + denovo137849 + denovo137881 + denovo137888 + denovo137912 + denovo137991 +
                                 denovo137998 + denovo138054 + denovo138099 + denovo138127 + denovo138173 + denovo138176 + denovo138225 + denovo138228 +
                                 denovo138233 + denovo138248 + denovo138266 + denovo138309 + denovo138313 + denovo138332 + denovo138360 + denovo138374 +
                                 denovo138375 + denovo138383 + denovo138431 + denovo138477 + denovo138486 + denovo138500 + denovo138501 + denovo138503 +
                                 denovo138504 + denovo138518 + denovo138524 + denovo138535 + denovo138557 + denovo138621 + denovo138641 + denovo138656 +
                                 denovo138662 + denovo138689 + denovo138741 + denovo138766 + denovo138774 + denovo138801 + denovo138815 + denovo138822 +
                                 denovo138850 + denovo138914 + denovo138947 + denovo138954 + denovo138957 + denovo138958 + denovo138967 + denovo138991 +
                                 denovo139052 + denovo139084 + denovo139136 + denovo139148 + denovo139156 + denovo139161 + denovo139186 + denovo139219 +
                                 denovo139237 + denovo139269 + denovo139275 + denovo139283 + denovo139344 + denovo139379 + denovo139386 + denovo139437 +
                                 denovo139438 + denovo139454 + denovo139466 + denovo139471 + denovo139477 + denovo139565 +
                                 denovo139581 + denovo139612 + denovo139638 + denovo139675 + denovo139683 + denovo139695 + denovo139722 + denovo139724 +
                                 denovo139740 + denovo139759 + denovo139772 + denovo139813 + denovo139852 + denovo139869 + denovo139880 + denovo139925 +
                                 denovo139994 + denovo140005 + denovo140026 + denovo140096 + denovo140101 + denovo140103 + denovo140116 +
                                 denovo140126 + denovo140169 + denovo140193 + denovo140220 + denovo140307 + denovo140332 + denovo140360 +
                                 denovo140363 + denovo140372 + denovo140384 + denovo140437 + denovo140455 + denovo140504 + denovo140538 +
                                 denovo140543 + denovo140584 + denovo140585 + denovo140625 + denovo140626 + denovo140653 + denovo140658 + denovo140685 +
                                 denovo140694 + denovo140749 + denovo140771 + denovo140787 + denovo140796 + denovo140870 + denovo140887 + denovo140895 +
                                 denovo140897 + denovo140970 + denovo140974 + denovo141042 + denovo141062 + denovo141085 + denovo141091 + denovo141149 +
                                 denovo141160 + denovo141196 + denovo141201 + denovo141206 + denovo141231 + denovo141256 + denovo141311 + denovo141331 +
                                 denovo141350 + denovo141361 + denovo141422 + denovo141435 + denovo141440 + denovo141472 + denovo141478 + denovo141506 +
                                 denovo141593 + denovo141648 + denovo141661 + denovo141678 + denovo141728 + denovo141732 + denovo141811 + denovo141821 +
                                 denovo141859 + denovo141870 + denovo141892 + denovo141899 + denovo141941 + denovo141947 + denovo142011 +
                                 denovo142035 + denovo142082 + denovo142106 + denovo142113 + denovo142166 + denovo142221 + denovo142243 + denovo142299 +
                                 denovo142302 + denovo142350 + denovo142365 + denovo142408 + denovo142446 + denovo142456 + denovo142479 + denovo142489 +
                                 denovo142514 + denovo142519 + denovo142588 + denovo142593 + denovo142635 + denovo142643 + denovo142645 + denovo142653 +
                                 denovo142654 + denovo142689 + denovo142690 + denovo142706 + denovo142712 + denovo142766 + denovo142821 + denovo142854 +
                                 denovo142920 + denovo142943 + denovo142956 + denovo142957 + denovo142964 + denovo142990 + denovo142991 + denovo143008 +
                                 denovo143014 + denovo143021 + denovo143034 + denovo143114 + denovo143135 + denovo143145 + denovo143160 + denovo143243 +
                                 denovo143269 + denovo143281 + denovo143301 + denovo143389 + denovo143420 + denovo143431 + denovo143435 + denovo143436 +
                                 denovo143440 + denovo143535 + denovo143541 + denovo143582 + denovo143593 + denovo143627 + denovo143666 + denovo143694 +
                                 denovo143734 + denovo143793 + denovo143847 + denovo143887 + denovo143888 + denovo143889 + denovo143902 +
                                 denovo143914 + denovo143919 + denovo143935 + denovo143946 + denovo143947 + denovo144004 + denovo144010 + denovo144038 +
                                 denovo144040 + denovo144046 + denovo144053 + denovo144073 + denovo144077 + denovo144119 + denovo144134 + denovo144191 +
                                 denovo144203 + denovo144206 + denovo144207 + denovo144230 + denovo144259 + denovo144275 + denovo144283 + denovo144294 +
                                 denovo144348 + denovo144356 + denovo144363 + denovo144442 + denovo144466 + denovo144516 + denovo144519 +
                                 denovo144525 + denovo144562 + denovo144584 + denovo144587 + denovo144604 + denovo144606 + denovo144610 + denovo144620 +
                                 denovo144657 + denovo144671 + denovo144697 + denovo144707 + denovo144712 + denovo144777 + denovo144819 + denovo144888 +
                                 denovo144900 + denovo144918 + denovo144983 + denovo145002 + denovo145075 + denovo145099 + denovo145102 + denovo145105 +
                                 denovo145118 + denovo145126 + denovo145137 + denovo145204 + denovo145281 + denovo145288 + denovo145289 + denovo145303 +
                                 denovo145319 + denovo145334 + denovo145392 + denovo145401 + denovo145413 + denovo145432 + denovo145441 + denovo145485 +
                                 denovo145492 + denovo145517 + denovo145520 + denovo145543 + denovo145562 + denovo145602 + denovo145623 + denovo145642 +
                                 denovo145673 + denovo145716 + denovo145859 + denovo145888 + denovo145939 + denovo145949 +
                                 denovo145999 + denovo146060 + denovo146078 + denovo146109 + denovo146129 + denovo146148 + denovo146159 +
                                 denovo146186 + denovo146190 + denovo146195 + denovo146212 + denovo146255 + denovo146276 + denovo146294 + denovo146306 +
                                 denovo146310 + denovo146361 + denovo146392 + denovo146399 + denovo146414 + denovo146424 + denovo146436 + denovo146479 +
                                 denovo146495 + denovo146516 + denovo146543 + denovo146596 + denovo146645 + denovo146697 + denovo146701 +
                                 denovo146772 + denovo146825 + denovo146840 + denovo146891 + denovo146942 + denovo146947 + denovo146954 + denovo146971 +
                                 denovo146982 + denovo146984 + denovo147036 + denovo147064 + denovo147069 + denovo147109 + denovo147133 + denovo147139 +
                                 denovo147158 + denovo147200 + denovo147243 + denovo147247 + denovo147252 + denovo147268 + denovo147281 + denovo147290 +
                                 denovo147291 + denovo147305 + denovo147314 + denovo147345 + denovo147349 + denovo147364 + denovo147375 + denovo147410 +
                                 denovo147424 + denovo147435 + denovo147496 + denovo147524 + denovo147559 + denovo147620 + denovo147627 + denovo147644 +
                                 denovo147669 + denovo147740 + denovo147741 + denovo147755 + denovo147767 + denovo147802 + denovo147804 + denovo147805 +
                                 denovo147816 + denovo147851 + denovo147944 + denovo147949 + denovo147961 + denovo147973 + denovo147994 + denovo147999 +
                                 denovo148044 + denovo148048 + denovo148074 + denovo148076 + denovo148085 + denovo148089 + denovo148124 + denovo148134 +
                                 denovo148171 + denovo148259 + denovo148270 + denovo148296 + denovo148354 + denovo148371 + denovo148374 +
                                 denovo148376 + denovo148381 + denovo148456 + denovo148459 + denovo148460 + denovo148463 + denovo148468 + denovo148508 +
                                 denovo148516 + denovo148541 + denovo148559 + denovo148616 + denovo148669 + denovo148691 + denovo148695 + denovo148790 +
                                 denovo148811 + denovo148821 + denovo148832 + denovo148838 + denovo148868 + denovo148871 + denovo148987 + denovo149045 +
                                 denovo149050 + denovo149095 + denovo149117 + denovo149163 + denovo149196 + denovo149230 + denovo149241 + denovo149328 +
                                 denovo149339 + denovo149382 + denovo149401 + denovo149405 + denovo149428 + denovo149436 + denovo149440 + denovo149449 +
                                 denovo149458 + denovo149462 + denovo149499 + denovo149528 + denovo149530 + denovo149547 + denovo149589 + denovo149594 +
                                 denovo149612 + denovo149664 + denovo149684 + denovo149731 + denovo149760 + denovo149761 + denovo149855 +
                                 denovo149875 + denovo149895 + denovo149913 + denovo149917 + denovo149962 + denovo150001 + denovo150020 +
                                 denovo150027 + denovo150028 + denovo150065 + denovo150067 + denovo150084 + denovo150085 + denovo150088 + denovo150090 +
                                 denovo150094 + denovo150111 + denovo150123 + denovo150129 + denovo150132 + denovo150180 + denovo150199 + denovo150242 +
                                 denovo150257 + denovo150291 + denovo150300 + denovo150390 + denovo150408 + denovo150411 + denovo150470 +
                                 denovo150541 + denovo150565 + denovo150602 + denovo150614 + denovo150619 + denovo150625 + denovo150627 + denovo150632 +
                                 denovo150679 + denovo150732 + denovo150788 + denovo150813 + denovo150853 + denovo150855 + denovo150868 + denovo150876 +
                                 denovo150910 + denovo150940 + denovo150950 + denovo150963 + denovo150967 + denovo150982 + denovo151017 + denovo151018 +
                                 denovo151030 + denovo151050 + denovo151055 + denovo151142 + denovo151150 + denovo151185 + denovo151215 + denovo151224 +
                                 denovo151242 + denovo151258 + denovo151267 + denovo151322 + denovo151323 + denovo151331 + denovo151402 + denovo151414 +
                                 denovo151433 + denovo151525 + denovo151568 + denovo151573 + denovo151582 + denovo151642 + denovo151685 + denovo151698 +
                                 denovo151733 + denovo151743 + denovo151784 + denovo151877 + denovo151886 + denovo151898 + denovo151922 + denovo151924 +
                                 denovo151944 + denovo151965 + denovo151971 + denovo152049 + denovo152075 + denovo152084 + denovo152115 + denovo152130 +
                                 denovo152220 + denovo152230 + denovo152283 + denovo152306 + denovo152322 + denovo152331 + denovo152345 +
                                 denovo152449 + denovo152453 + denovo152483 + denovo152512 + denovo152516 + denovo152519 + denovo152539 + denovo152545 +
                                 denovo152575 + denovo152620 + denovo152648 + denovo152680 + denovo152688 + denovo152750 + denovo152766 + denovo152781 +
                                 denovo152786 + denovo152857 + denovo152863 + denovo152888 + denovo152902 + denovo152928 + denovo152936 + denovo152953 +
                                 denovo152963 + denovo152968 + denovo152979 + denovo153009 + denovo153025 + denovo153031 + denovo153033 + denovo153075 +
                                 denovo153121 + denovo153138 + denovo153164 + denovo153207 + denovo153219 + denovo153242 + denovo153262 + denovo153332 +
                                 denovo153471 + denovo153472 + denovo153494 + denovo153526 + denovo153552 + denovo153580 + denovo153591 + denovo153601 +
                                 denovo153644 + denovo153662 + denovo153746 + denovo153756 + denovo153762 + denovo153764 + denovo153839 +
                                 denovo153886 + denovo153937 + denovo153999 + denovo154074 + denovo154104 + denovo154118 + denovo154144 +
                                 denovo154161 + denovo154203 + denovo154209 + denovo154234 + denovo154242 + denovo154283 + denovo154305 + denovo154315 +
                                 denovo154332 + denovo154347 + denovo154351 + denovo154362 + denovo154373 + denovo154388 + denovo154413 + denovo154504 +
                                 denovo154506 + denovo154526 + denovo154541 + denovo154548 + denovo154586 + denovo154616 + denovo154671 + denovo154686 +
                                 denovo154705 + denovo154759 + denovo154760 + denovo154876 + denovo154895 + denovo154917 + denovo154931 +
                                 denovo154935 + denovo154983 + denovo154988 + denovo155045 + denovo155104 + denovo155125 + denovo155137 + denovo155162 +
                                 denovo155194 + denovo155196 + denovo155203 + denovo155244 + denovo155264 + denovo155274 + denovo155334 + denovo155400 +
                                 denovo155469 + denovo155495 + denovo155501 + denovo155511 + denovo155527 + denovo155529 + denovo155563 + denovo155571 +
                                 denovo155573 + denovo155579 + denovo155596 + denovo155603 + denovo155616 + denovo155618 + denovo155636 +
                                 denovo155655 + denovo155699 + denovo155710 + denovo155743 + denovo155745 + denovo155842 + denovo155853 + denovo155870 +
                                 denovo155877 + denovo155929 + denovo155950 + denovo155962 + denovo155973 + denovo155982 + denovo155983 + denovo156065 +
                                 denovo156084 + denovo156087 + denovo156104 + denovo156108 + denovo156109 + denovo156113 + denovo156149 + denovo156156 +
                                 denovo156174 + denovo156186 + denovo156218 + denovo156219 + denovo156228 + denovo156281 + denovo156309 + denovo156313 +
                                 denovo156366 + denovo156377 + denovo156389 + denovo156422 + denovo156423 + denovo156500 + denovo156555 + denovo156557 +
                                 denovo156606 + denovo156742 + denovo156755 + denovo156756 + denovo156784 + denovo156814 + denovo156843 + denovo156861 +
                                 denovo156912 + denovo156931 + denovo156944 + denovo156949 + denovo156972 + denovo156975 + denovo157048 + denovo157085 +
                                 denovo157110 + denovo157135 + denovo157176 + denovo157232 + denovo157245 + denovo157289 + denovo157329 + denovo157347 +
                                 denovo157379 + denovo157433 + denovo157627 + denovo157657 + denovo157677 + denovo157690 + denovo157694 + denovo157749 +
                                 denovo157759 + denovo157764 + denovo157770 + denovo157777 + denovo157811 + denovo157827 + denovo157834 + denovo157849 +
                                 denovo157916 + denovo157945 + denovo157946 + denovo157995 + denovo158048 + denovo158080 + denovo158192 + denovo158202 +
                                 denovo158210 + denovo158241 + denovo158282 + denovo158349 + denovo158378 + denovo158400 + denovo158405 + denovo158407 +
                                 denovo158427 + denovo158483 + denovo158499 + denovo158505 + denovo158512 + denovo158513 + denovo158521 + denovo158528 +
                                 denovo158548 + denovo158571 + denovo158579 + denovo158615 + denovo158620 + denovo158626 + denovo158632 + denovo158653 +
                                 denovo158658 + denovo158659 + denovo158691 + denovo158692 + denovo158724 + denovo158725 + denovo158735 + denovo158757 +
                                 denovo158795 + denovo158825 + denovo158833 + denovo158907 + denovo158908 + denovo158913 + denovo159019 + denovo159033 +
                                 denovo159045 + denovo159057 + denovo159102 + denovo159109 + denovo159167 + denovo159171 + denovo159212 + denovo159258 +
                                 denovo159274 + denovo159327 + denovo159408 + denovo159420 + denovo159444 + denovo159450 + denovo159468 + denovo159503 +
                                 denovo159507 + denovo159573 + denovo159619 + denovo159622 + denovo159629 + denovo159639 + denovo159642 + denovo159686 +
                                 denovo159695 + denovo159718 + denovo159727 + denovo159745 + denovo159765 + denovo159849 + denovo159913 + denovo159918 +
                                 denovo160177 + denovo160247 + denovo160352 + denovo160386 + denovo160443 + denovo160493 + denovo160514 +
                                 denovo160530 + denovo160619 + denovo160734 + denovo160747 + denovo160773 + denovo160792 + denovo160797 +
                                 denovo160805 + denovo160833 + denovo160863 + denovo160887 + denovo160905 + denovo160995 + denovo161000 + denovo161022 +
                                 denovo161049 + denovo161058 + denovo161130 + denovo161185 + denovo161222 + denovo161229 + denovo161258 + denovo161264 +
                                 denovo161294 + denovo161305 + denovo161349 + denovo161365 + denovo161369 + denovo161405 + denovo161446 + denovo161454 +
                                 denovo161540 + denovo161594 + denovo161655 + denovo161676 + denovo161680 + denovo161694 + denovo161702 +
                                 denovo161707 + denovo161713 + denovo161734 + denovo161745 + denovo161789 + denovo161806 + denovo161809 + denovo161810 +
                                 denovo161870 + denovo161916 + denovo161942 + denovo161952 + denovo161968 + denovo162000 + denovo162014 + denovo162019 +
                                 denovo162065 + denovo162067 + denovo162098 + denovo162101 + denovo162151 + denovo162174 + denovo162176 +
                                 denovo162182 + denovo162215 + denovo162219 + denovo162253 + denovo162293 + denovo162340 + denovo162361 + denovo162390 +
                                 denovo162499 + denovo162506 + denovo162520 + denovo162524 + denovo162564 + denovo162600 + denovo162645 + denovo162697 +
                                 denovo162738 + denovo162743 + denovo162754 + denovo162755 + denovo162788 + denovo162799 + denovo162813 + denovo162893 +
                                 denovo162907 + denovo162922 + denovo162940 + denovo162974 + denovo163045 + denovo163052 + denovo163059 + denovo163072 +
                                 denovo163087 + denovo163128 + denovo163131 + denovo163159 + denovo163170 + denovo163211 + denovo163260 + denovo163296 +
                                 denovo163351 + denovo163409 + denovo163484 + denovo163499 + denovo163507 + denovo163530 + denovo163573 + denovo163575 +
                                 denovo163581 + denovo163662 + denovo163679 + denovo163688 + denovo163698 + denovo163725 + denovo163729 + denovo163810 +
                                 denovo163823 + denovo163857 + denovo163873 + denovo163927 + denovo163996 + denovo164006 + denovo164014 + denovo164117 +
                                 denovo164129 + denovo164130 + denovo164178 + denovo164212 + denovo164241 + denovo164263 + denovo164296 + denovo164298 +
                                 denovo164332 + denovo164377 + denovo164450 + denovo164464 + denovo164520 + denovo164608 + denovo164618 + denovo164625 +
                                 denovo164647 + denovo164656 + denovo164658 + denovo164669 + denovo164764 + denovo164770 + denovo164782 +
                                 denovo164808 + denovo164825 + denovo164851 + denovo164855 + denovo164883 + denovo164902 + denovo164936 + denovo164939 +
                                 denovo164959 + denovo164966 + denovo164975 + denovo164981 + denovo165053 + denovo165084 + denovo165093 + denovo165099 +
                                 denovo165110 + denovo165123 + denovo165139 + denovo165197 + denovo165214 + denovo165231 + denovo165340 + denovo165343 +
                                 denovo165348 + denovo165421 + denovo165465 + denovo165482 + denovo165506 + denovo165508 + denovo165548 + denovo165564 +
                                 denovo165568 + denovo165615 + denovo165621 + denovo165638 + denovo165662 + denovo165730 + denovo165751 + denovo165841 +
                                 denovo165845 + denovo165857 + denovo165860 + denovo165893 + denovo165899 + denovo165924 + denovo165927 + denovo165950 +
                                 denovo165951 + denovo165988 + denovo166030 + denovo166148 + denovo166168 + denovo166172 + denovo166188 + denovo166192 +
                                 denovo166243 + denovo166259 + denovo166310 + denovo166336 + denovo166363 + denovo166394 + denovo166400 + denovo166443 +
                                 denovo166561 + denovo166578 + denovo166592 + denovo166621 + denovo166662 + denovo166680 + denovo166730 + denovo166747 +
                                 denovo166783 + denovo166818 + denovo166848 + denovo166865 + denovo166882 + denovo166910 + denovo166936 + denovo166947 +
                                 denovo166994 + denovo166995 + denovo167037 + denovo167080 + denovo167136 + denovo167187 + denovo167228 +
                                 denovo167254 + denovo167319 + denovo167326 + denovo167377 + denovo167399 + denovo167402 + denovo167425 + denovo167432 +
                                 denovo167492 + denovo167505 + denovo167517 + denovo167527 + denovo167535 + denovo167551 + denovo167577 +
                                 denovo167585 + denovo167617 + denovo167665 + denovo167666 + denovo167718 + denovo167735 + denovo167738 + denovo167744 +
                                 denovo167751 + denovo167774 + denovo167788 + denovo167790 + denovo167837 + denovo167849 + denovo167863 + denovo167924 +
                                 denovo167925 + denovo167936 + denovo167973 + denovo167979 + denovo167993 + denovo168068 + denovo168150 +
                                 denovo168166 + denovo168201 + denovo168204 + denovo168218 + denovo168223 + denovo168225 + denovo168234 + denovo168250 +
                                 denovo168270 + denovo168294 + denovo168316 + denovo168393 + denovo168430 + denovo168464 + denovo168474 + denovo168497 +
                                 denovo168538 + denovo168552 + denovo168565 + denovo168576 + denovo168616 + denovo168642 + denovo168652 + denovo168674 +
                                 denovo168704 + denovo168725 + denovo168767 + denovo168775 + denovo168800 + denovo168847 + denovo168851 + denovo168912 +
                                 denovo168919 + denovo168986 + denovo168994 + denovo169055 + denovo169118 + denovo169121 + denovo169122 + denovo169146 +
                                 denovo169147 + denovo169222 + denovo169250 + denovo169311 + denovo169315 + denovo169372 + denovo169411 + denovo169417 +
                                 denovo169445 + denovo169448 + denovo169456 + denovo169464 + denovo169465 + denovo169473 + denovo169510 + denovo169689 +
                                 denovo169694 + denovo169695 + denovo169722 + denovo169732 + denovo169737 + denovo169756 + denovo169840 +
                                 denovo169850 + denovo169853 + denovo169866 + denovo169900 + denovo169922 + denovo169931 + denovo169953 + denovo169954 +
                                 denovo169956 + denovo169963 + denovo170099 + denovo170112 + denovo170180 + denovo170221 + denovo170253 + denovo170259 +
                                 denovo170334 + denovo170353 + denovo170385 + denovo170451 + denovo170455 + denovo170474 + denovo170533 + denovo170545 +
                                 denovo170589 + denovo170602 + denovo170652 + denovo170689 + denovo170775 + denovo170788 + denovo170799 + denovo170879 +
                                 denovo170932 + denovo170954 + denovo171006 + denovo171034 + denovo171040 + denovo171059 + denovo171063 + denovo171127 +
                                 denovo171138 + denovo171141 + denovo171143 + denovo171162 + denovo171194 + denovo171222 + denovo171231 + denovo171232 +
                                 denovo171239 + denovo171263 + denovo171287 + denovo171295 + denovo171321 + denovo171368 + denovo171372 +
                                 denovo171384 + denovo171449 + denovo171455 + denovo171459 + denovo171503 + denovo171567 + denovo171613 +
                                 denovo171639 + denovo171656 + denovo171720 + denovo171721 + denovo171733 + denovo171737 + denovo171746 + denovo171786 +
                                 denovo171824 + denovo171852 + denovo171855 + denovo171881 + denovo171914 + denovo171925 + denovo171998 + denovo172032 +
                                 denovo172173 + denovo172203 + denovo172255 + denovo172272 + denovo172274 + denovo172301 + denovo172336 + denovo172340 +
                                 denovo172422 + denovo172477 + denovo172479 + denovo172497 + denovo172541 + denovo172556 + denovo172581 + denovo172589 +
                                 denovo172607 + denovo172617 + denovo172630 + denovo172691 + denovo172744 + denovo172785 + denovo172796 + denovo172839 +
                                 denovo172862 + denovo172950 + denovo173049 + denovo173136 + denovo173178 + denovo173199 + denovo173240 +
                                 denovo173242 + denovo173272 + denovo173299 + denovo173328 + denovo173337 + denovo173339 + denovo173367 + denovo173397 +
                                 denovo173457 + denovo173475 + denovo173477 + denovo173491 + denovo173496 + denovo173508 + denovo173575 + denovo173591 +
                                 denovo173609 + denovo173634 + denovo173638 + denovo173691 + denovo173715 + denovo173738 + denovo173766 + denovo173772 +
                                 denovo173775 + denovo173828 + denovo174012 + denovo174014 + denovo174103 + denovo174123 + denovo174153 + denovo174157 +
                                 denovo174162 + denovo174198 + denovo174201 + denovo174212 + denovo174216 + denovo174217 + denovo174266 + denovo174271 +
                                 denovo174278 + denovo174293 + denovo174303 + denovo174385 + denovo174386 + denovo174430 + denovo174435 + denovo174459 +
                                 denovo174464 + denovo174470 + denovo174497 + denovo174518 + denovo174570 + denovo174613 + denovo174618 + denovo174655 +
                                 denovo174676 + denovo174691 + denovo174703 + denovo174765 + denovo174773 + denovo174780 + denovo174827 + denovo174846 +
                                 denovo174852 + denovo174857 + denovo174868 + denovo174886 + denovo174941 + denovo174974 + denovo175001 + denovo175012 +
                                 denovo175060 + denovo175095 + denovo175122 + denovo175143 + denovo175198 + denovo175200 + denovo175207 + denovo175208 +
                                 denovo175228 + denovo175317 + denovo175326 + denovo175332 + denovo175360 + denovo175367 + denovo175420 + denovo175461 +
                                 denovo175494 + denovo175503 + denovo175508 + denovo175538 + denovo175543 + denovo175558 + denovo175590 +
                                 denovo175645 + denovo175709 + denovo175713 + denovo175732 + denovo175733 + denovo175763 + denovo175777 + denovo175778 +
                                 denovo175939 + denovo175972 + denovo175974 + denovo175978 + denovo175983 + denovo175995 + denovo175996 + denovo176025 +
                                 denovo176046 + denovo176094 + denovo176100 + denovo176172 + denovo176206 + denovo176235 + denovo176238 + denovo176255 +
                                 denovo176337 + denovo176364 + denovo176412 + denovo176416 + denovo176418 + denovo176463 + denovo176467 + denovo176471 +
                                 denovo176503 + denovo176547 + denovo176561 + denovo176617 + denovo176634 + denovo176682 + denovo176685 +
                                 denovo176718 + denovo176757 + denovo176760 + denovo176765 + denovo176778 + denovo176800 + denovo176852 + denovo176861 +
                                 denovo176903 + denovo176966 + denovo176981 + denovo176988 + denovo176995 + denovo177020 + denovo177055 +
                                 denovo177096 + denovo177145 + denovo177198 + denovo177236 + denovo177260 + denovo177271 + denovo177279 + denovo177316 +
                                 denovo177339 + denovo177354 + denovo177366 + denovo177393 + denovo177406 + denovo177411 + denovo177500 + denovo177520 +
                                 denovo177557 + denovo177559 + denovo177590 + denovo177599 + denovo177707 + denovo177723 + denovo177769 +
                                 denovo177814 + denovo177830 + denovo177839 + denovo177843 + denovo177863 + denovo177898 + denovo177917 + denovo177929 +
                                 denovo177947 + denovo177970 + denovo178002 + denovo178062 + denovo178064 + denovo178070 + denovo178109 + denovo178137 +
                                 denovo178152 + denovo178161 + denovo178187 + denovo178202 + denovo178205 + denovo178212 + denovo178214 + denovo178295 +
                                 denovo178298 + denovo178300 + denovo178314 + denovo178326 + denovo178359 + denovo178397 + denovo178409 + denovo178423 +
                                 denovo178425 + denovo178501 + denovo178510 + denovo178514 + denovo178545 + denovo178546 + denovo178625 + denovo178646 +
                                 denovo178672 + denovo178716 + denovo178735 + denovo178829 + denovo178831 + denovo178903 + denovo178927 + denovo178990 +
                                 denovo179025 + denovo179027 + denovo179048 + denovo179086 + denovo179144 + denovo179159 + denovo179194 + denovo179210 +
                                 denovo179228 + denovo179243 + denovo179317 + denovo179330 + denovo179334 + denovo179341 + denovo179362 + denovo179384 +
                                 denovo179413 + denovo179416 + denovo179443 + denovo179481 + denovo179493 + denovo179499 + denovo179502 + denovo179513 +
                                 denovo179522 + denovo179523 + denovo179535 + denovo179538 + denovo179578 + denovo179583 + denovo179591 + denovo179593 +
                                 denovo179594 + denovo179613 + denovo179645 + denovo179693 + denovo179766 + denovo179785 + denovo179795 + denovo179798 +
                                 denovo179825 + denovo179835 + denovo179846 + denovo179873 + denovo179874 + denovo179907 + denovo179911 + denovo179915 +
                                 denovo179922 + denovo179925 + denovo179936 + denovo179999 + denovo180017 + denovo180072 + denovo180090 + denovo180154 +
                                 denovo180163 + denovo180171 + denovo180187 + denovo180191 + denovo180202 + denovo180218 + denovo180228 +
                                 denovo180244 + denovo180257 + denovo180276 + denovo180369 + denovo180385 + denovo180391 + denovo180397 + denovo180411 +
                                 denovo180417 + denovo180421 + denovo180444 + denovo180475 + denovo180487 + denovo180532 + denovo180618 + denovo180672 +
                                 denovo180767 + denovo180768 + denovo180788 + denovo180822 + denovo180853 + denovo180872 + denovo180911 +
                                 denovo180934 + denovo181058 + denovo181102 + denovo181118 + denovo181155 + denovo181192 + denovo181197 + denovo181206 +
                                 denovo181235 + denovo181244 + denovo181269 + denovo181288 + denovo181305 + denovo181340 + denovo181392 + denovo181455 +
                                 denovo181502 + denovo181530 + denovo181536 + denovo181545 + denovo181606 + denovo181636 + denovo181773 + denovo181783 +
                                 denovo181809 + denovo181834 + denovo181835 + denovo181874 + denovo181883 + denovo181892 + denovo181903 + denovo181909 +
                                 denovo181928 + denovo181998 + denovo182030 + denovo182057 + denovo182081 + denovo182086 + denovo182122 +
                                 denovo182172 + denovo182236 + denovo182267 + denovo182297 + denovo182322 + denovo182328 + denovo182351 + denovo182352 +
                                 denovo182356 + denovo182362 + denovo182366 + denovo182390 + denovo182403 + denovo182451 + denovo182496 + denovo182558 +
                                 denovo182593 + denovo182596 + denovo182638 + denovo182645 + denovo182743 + denovo182748 +
                                 denovo182756 + denovo182797 + denovo182820 + denovo182887 + denovo182895 + denovo182898 + denovo182906 + denovo182908 +
                                 denovo182909 + denovo182927 + denovo182959 + denovo182964 + denovo182971 + denovo182997 + denovo183026 + denovo183067 +
                                 denovo183072 + denovo183090 + denovo183095 + denovo183122 + denovo183159 + denovo183183 + denovo183190 + denovo183222 +
                                 denovo183226 + denovo183242 + denovo183248 + denovo183249 + denovo183261 + denovo183270 + denovo183305 + denovo183318 +
                                 denovo183323 + denovo183335 + denovo183355 + denovo183365 + denovo183385 + denovo183412 + denovo183440 + denovo183475 +
                                 denovo183490 + denovo183516 + denovo183562 + denovo183627 + denovo183650 + denovo183654 + denovo183699 + denovo183711 +
                                 denovo183723 + denovo183736 + denovo183748 + denovo183749 + denovo183822 + denovo183907 + denovo183916 +
                                 denovo183969 + denovo183975 + denovo184010 + denovo184034 + denovo184053 + denovo184063 + denovo184123 +
                                 denovo184124 + denovo184126 + denovo184238 + denovo184259 + denovo184261 + denovo184279 + denovo184296 + denovo184320 +
                                 denovo184371 + denovo184418 + denovo184460 + denovo184553 + denovo184560 + denovo184561 + denovo184568 +
                                 denovo184600 + denovo184607 + denovo184620 + denovo184628 + denovo184629 + denovo184637 + denovo184645 + denovo184674 +
                                 denovo184692 + denovo184700 + denovo184714 + denovo184733 + denovo184739 + denovo184751 + denovo184757 + denovo184789 +
                                 denovo184819 + denovo184827 + denovo184841 + denovo184861 + denovo184862 + denovo184865 + denovo184880 + denovo184924 +
                                 denovo184930 + denovo184999 + denovo185008 + denovo185025 + denovo185048 + denovo185056 + denovo185061 + denovo185071 +
                                 denovo185112 + denovo185164 + denovo185169 + denovo185217 + denovo185222 + denovo185240 + denovo185280 +
                                 denovo185331 + denovo185388 + denovo185422 + denovo185429 + denovo185432 + denovo185457 + denovo185473 + denovo185544 +
                                 denovo185545 + denovo185583 + denovo185608 + denovo185620 + denovo185635 + denovo185668 + denovo185808 + denovo185844 +
                                 denovo185845 + denovo185933 + denovo185938 + denovo185946 + denovo185965 + denovo185979 + denovo185980 +
                                 denovo186021 + denovo186022 + denovo186024 + denovo186051 + denovo186085 + denovo186110 + denovo186125 + denovo186145 +
                                 denovo186223 + denovo186225 + denovo186260 + denovo186285 + denovo186330 + denovo186353 + denovo186376 + denovo186388 +
                                 denovo186453 + denovo186462 + denovo186499 + denovo186502 + denovo186519 + denovo186534 + denovo186576 + denovo186652 +
                                 denovo186669 + denovo186689 + denovo186740 + denovo186756 + denovo186799 + denovo186815 + denovo186818 + denovo186855 +
                                 denovo186895 + denovo186916 + denovo186938 + denovo186974 + denovo186975 + denovo186994 + denovo186995 + denovo187013 +
                                 denovo187033 + denovo187068 + denovo187095 + denovo187154 + denovo187172 + denovo187252 + denovo187262 + denovo187278 +
                                 denovo187282 + denovo187365 + denovo187427 + denovo187440 + denovo187452 + denovo187490 + denovo187527 + denovo187528 +
                                 denovo187554 + denovo187603 + denovo187626 + denovo187685 + denovo187785 + denovo187789 + denovo187812 + denovo187914 +
                                 denovo187919 + denovo188002 + denovo188024 + denovo188060 + denovo188067 + denovo188091 + denovo188097 + denovo188101 +
                                 denovo188137 + denovo188159 + denovo188183 + denovo188234 + denovo188262 + denovo188300 + denovo188328 + denovo188376 +
                                 denovo188438 + denovo188471 + denovo188507 + denovo188526 + denovo188560 + denovo188576 + denovo188578 + denovo188587 +
                                 denovo188616 + denovo188639 + denovo188791 + denovo188840 + denovo188843 + denovo188854 + denovo188858 + denovo188895 +
                                 denovo188957 + denovo189047 + denovo189056 + denovo189091 + denovo189111 + denovo189154 + denovo189190 +
                                 denovo189192 + denovo189207 + denovo189218 + denovo189230 + denovo189272 + denovo189284 + denovo189354 + denovo189368 +
                                 denovo189388 + denovo189398 + denovo189430 + denovo189465 + denovo189475 + denovo189529 + denovo189594 + denovo189629 +
                                 denovo189654 + denovo189660 + denovo189665 + denovo189666 + denovo189671 + denovo189695 + denovo189762 + denovo189765 +
                                 denovo189833 + denovo189853 + denovo189854 + denovo189864 + denovo189893 + denovo189905 + denovo189959 + denovo189981 +
                                 denovo189993 + denovo190018 + denovo190047 + denovo190068 + denovo190071 + denovo190106 + denovo190157 + denovo190183 +
                                 denovo190205 + denovo190213 + denovo190269 + denovo190322 + denovo190324 + denovo190383 + denovo190400 +
                                 denovo190403 + denovo190409 + denovo190438 + denovo190461 + denovo190513 + denovo190568 + denovo190589 + denovo190599 +
                                 denovo190611 + denovo190614 + denovo190651 + denovo190668 + denovo190669 + denovo190699 + denovo190701 + denovo190714 +
                                 denovo190776 + denovo190785 + denovo190810 + denovo190831 + denovo190876 + denovo190883 + denovo190886 + denovo190897 +
                                 denovo190898 + denovo190943 + denovo190963 + denovo191057 + denovo191061 + denovo191084 + denovo191087 + denovo191146 +
                                 denovo191195 + denovo191196 + denovo191231 + denovo191237 + denovo191247 + denovo191282 + denovo191314 + denovo191336 +
                                 denovo191369 + denovo191388 + denovo191407 + denovo191424 + denovo191483 + denovo191537 + denovo191539 +
                                 denovo191547 + denovo191603 + denovo191659 + denovo191727 + denovo191731 + denovo191771 + denovo191825 + denovo191831 +
                                 denovo191848 + denovo191852 + denovo191874 + denovo191925 + denovo191948 + denovo191982 + denovo191995 + denovo192007 +
                                 denovo192020 + denovo192191 + denovo192231 + denovo192268 + denovo192304 + denovo192319 + denovo192326 + denovo192342 +
                                 denovo192427 + denovo192443 + denovo192452 + denovo192460 + denovo192482 + denovo192483 + denovo192580 +
                                 denovo192593 + denovo192615 + denovo192620 + denovo192637 + denovo192672 + denovo192718 + denovo192733 + denovo192739 +
                                 denovo192827 + denovo192858 + denovo192884 + denovo192965 + denovo192976 + denovo193072 + denovo193121 +
                                 denovo193191 + denovo193197 + denovo193205 + denovo193220 + denovo193269 + denovo193308 + denovo193319 + denovo193328 +
                                 denovo193333 + denovo193371 + denovo193392 + denovo193467 + denovo193471 + denovo193512 + denovo193542 + denovo193544 +
                                 denovo193581 + denovo193620 + denovo193627 + denovo193642 + denovo193664 + denovo193714 + denovo193723 + denovo193739 +
                                 denovo193740 + denovo193744 + denovo193791 + denovo193807 + denovo193817 + denovo193823 + denovo193853 + denovo193863 +
                                 denovo193918 + denovo193964 + denovo193971 + denovo193990 + denovo193999 + denovo194010 + denovo194014 +
                                 denovo194024 + denovo194035 + denovo194038 + denovo194046 + denovo194058 + denovo194077 + denovo194079 + denovo194096 +
                                 denovo194133 + denovo194171 + denovo194231 + denovo194298 + denovo194455 + denovo194524 + denovo194537 +
                                 denovo194549 + denovo194594 + denovo194628 + denovo194688 + denovo194709 + denovo194713 + denovo194724 + denovo194725 +
                                 denovo194764 + denovo194768 + denovo194769 + denovo194774 + denovo194776 + denovo194777 + denovo194808 + denovo194826 +
                                 denovo194888 + denovo194949 + denovo194970 + denovo195020 + denovo195087 + denovo195091 + denovo195103 + denovo195123 +
                                 denovo195174 + denovo195207 + denovo195222 + denovo195224 + denovo195251 + denovo195270 + denovo195272 + denovo195284 +
                                 denovo195385 + denovo195457 + denovo195470 + denovo195502 + denovo195526 + denovo195527 + denovo195562 + denovo195573 +
                                 denovo195632 + denovo195634 + denovo195639 + denovo195698 + denovo195700 + denovo195705 + denovo195731 + denovo195778 +
                                 denovo195791 + denovo195804 + denovo195819 + denovo195858 + denovo195871 + denovo195875 + denovo195885 + denovo195888 +
                                 denovo195918 + denovo195919 + denovo195944 + denovo195973 + denovo196007 + denovo196076 + denovo196089 + denovo196124 +
                                 denovo196129 + denovo196153 + denovo196161 + denovo196168 + denovo196207 + denovo196208 + denovo196237 + denovo196244 +
                                 denovo196247 + denovo196253 + denovo196265 + denovo196270 + denovo196272 + denovo196365 + denovo196372 + denovo196387 +
                                 denovo196466 + denovo196500 + denovo196575 + denovo196677 + denovo196693 + denovo196695 + denovo196696 +
                                 denovo196736 + denovo196737 + denovo196761 + denovo196791 + denovo196813 + denovo196824 + denovo196857 + denovo196863 +
                                 denovo196865 + denovo196899 + denovo196933 + denovo196950 + denovo196989 + denovo196998 + denovo196999 + denovo197012 +
                                 denovo197032 + denovo197033 + denovo197039 + denovo197074 + denovo197078 + denovo197101 + denovo197132 + denovo197136 +
                                 denovo197137 + denovo197189 + denovo197196 + denovo197199 + denovo197220 + denovo197241 + denovo197248 + denovo197250 +
                                 denovo197255 + denovo197371 + denovo197383 + denovo197400 + denovo197421 + denovo197460 +
                                 denovo197461 + denovo197524 + denovo197578 + denovo197599 + denovo197600 + denovo197612 + denovo197634 + denovo197649 +
                                 denovo197650 + denovo197681 + denovo197682 + denovo197715 + denovo197736 + denovo197738 + denovo197775 + denovo197788 +
                                 denovo197839 + denovo197845 + denovo197857 + denovo197877 + denovo197880 + denovo197904 + denovo197933 + denovo197964 +
                                 denovo197988 + denovo197994 + denovo198008 + denovo198060 + denovo198081 + denovo198098 + denovo198102 +
                                 denovo198124 + denovo198168 + denovo198213 + denovo198318 + denovo198392 + denovo198414 + denovo198471 + denovo198484 +
                                 denovo198494 + denovo198501 + denovo198563 + denovo198570 + denovo198587 + denovo198613 + denovo198614 + denovo198624 +
                                 denovo198648 + denovo198651 + denovo198663 + denovo198680 + denovo198685 + denovo198696 + denovo198726 + denovo198780 +
                                 denovo198781 + denovo198786 + denovo198790 + denovo198814 + denovo198858 + denovo198878 + denovo198897 + denovo198916 +
                                 denovo198928 + denovo198941 + denovo199007 + denovo199032 + denovo199052 + denovo199112 +
                                 denovo199116 + denovo199127 + denovo199144 + denovo199172 + denovo199184 + denovo199242 + denovo199287 + denovo199342 +
                                 denovo199344 + denovo199360 + denovo199373 + denovo199404 + denovo199452 + denovo199454 + denovo199463 +
                                 denovo199515 + denovo199521 + denovo199546 + denovo199605 + denovo199611 + denovo199613 + denovo199715 + denovo199771 +
                                 denovo199783 + denovo199788 + denovo199807 + denovo199828 + denovo199854 + denovo199864 + denovo199881 + denovo199901 +
                                 denovo199912 + denovo199944 + denovo199955 + denovo200029 + denovo200044 + denovo200047 + denovo200072 + denovo200080 +
                                 denovo200091 + denovo200094 + denovo200104 + denovo200115 + denovo200174 + denovo200183 + denovo200189 + denovo200222 +
                                 denovo200258 + denovo200262 + denovo200270 + denovo200280 + denovo200310 + denovo200343 + denovo200351 + denovo200408 +
                                 denovo200422 + denovo200423 + denovo200439 + denovo200445 + denovo200447 + denovo200495 + denovo200510 + denovo200531 +
                                 denovo200543 + denovo200597 + denovo200602 + denovo200608 + denovo200617 + denovo200635 + denovo200722 +
                                 denovo200767 + denovo200792 + denovo200839 + denovo200872 + denovo200893 + denovo200902 + denovo200905 + denovo200935 +
                                 denovo200947 + denovo200977 + denovo201015 + denovo201089 + denovo201109 + denovo201153 + denovo201162 + denovo201214 +
                                 denovo201266 + denovo201305 + denovo201328 + denovo201331 + denovo201370 + denovo201375 + denovo201376 + denovo201385 +
                                 denovo201419 + denovo201432 + denovo201492 + denovo201516 + denovo201534 + denovo201563 + denovo201571 + denovo201605 +
                                 denovo201622 + denovo201629 + denovo201748 + denovo201796 + denovo201803 + denovo202026 +
                                 denovo202032 + denovo202097 + denovo202150 + denovo202202 + denovo202219 + denovo202282 + denovo202396 + denovo202408 +
                                 denovo202444 + denovo202448 + denovo202479 + denovo202488 + denovo202524 + denovo202544 + denovo202551 +
                                 denovo202567 + denovo202613 + denovo202654 + denovo202679 + denovo202704 + denovo202723 + denovo202731 + denovo202754 +
                                 denovo202846 + denovo202867 + denovo202871 + denovo202960 + denovo203003 + denovo203007 + denovo203070 +
                                 denovo203079 + denovo203083 + denovo203117 + denovo203136 + denovo203151 + denovo203169 + denovo203172 + denovo203187 +
                                 denovo203209 + denovo203257 + denovo203273 + denovo203276 + denovo203284 + denovo203295 + denovo203305 +
                                 denovo203317 + denovo203320 + denovo203348 + denovo203373 + denovo203426 + denovo203431 + denovo203514 + denovo203590 +
                                 denovo203634 + denovo203671 + denovo203812 + denovo203891 + denovo203906 + denovo203909 + denovo203912 +
                                 denovo203965 + denovo204036 + denovo204045 + denovo204077 + denovo204089 + denovo204166 + denovo204189 + denovo204225 +
                                 denovo204237 + denovo204289 + denovo204292 + denovo204319 + denovo204396 + denovo204409 + denovo204419 + denovo204425 +
                                 denovo204457 + denovo204471 + denovo204475 + denovo204520 + denovo204571 + denovo204590 + denovo204610 + denovo204657 +
                                 denovo204694 + denovo204703 + denovo204771 + denovo204775 + denovo204802 + denovo204831 + denovo204836 +
                                 denovo204840 + denovo204856 + denovo204882 + denovo204923 + denovo204965 + denovo204989 + denovo205029 + denovo205062 +
                                 denovo205113 + denovo205127 + denovo205183 + denovo205192 + denovo205239 + denovo205242 +
                                 denovo205246 + denovo205252 + denovo205261 + denovo205263 + denovo205312 + denovo205354 + denovo205367 +
                                 denovo205491 + denovo205496 + denovo205519 + denovo205524 + denovo205530 + denovo205555 + denovo205564 + denovo205568 +
                                 denovo205580 + denovo205585 + denovo205594 + denovo205605 + denovo205653 + denovo205664 + denovo205688 + denovo205696 +
                                 denovo205703 + denovo205768 + denovo205769 + denovo205780 + denovo205796 + denovo205799 + denovo205836 + denovo205841 +
                                 denovo205868 + denovo205892 + denovo205894 + denovo205902 + denovo205924 + denovo205939 + denovo205946 + denovo205990 +
                                 denovo206025 + denovo206071 + denovo206078 + denovo206094 + denovo206100 + denovo206171 + denovo206216 + denovo206242 +
                                 denovo206252 + denovo206261 + denovo206315 + denovo206332 + denovo206345 + denovo206356 + denovo206360 + denovo206367 +
                                 denovo206454 + denovo206487 + denovo206494 + denovo206509 + denovo206542 + denovo206546 + denovo206554 + denovo206571 +
                                 denovo206576 + denovo206674 + denovo206686 + denovo206727 + denovo206743 + denovo206749 + denovo206816 + denovo206831 +
                                 denovo206896 + denovo206900 + denovo206942 + denovo206957 + denovo207002 + denovo207059 + denovo207076 +
                                 denovo207080 + denovo207103 + denovo207118 + denovo207137 + denovo207169 + denovo207176 + denovo207182 + denovo207234 +
                                 denovo207240 + denovo207257 + denovo207259 + denovo207392 + denovo207412 + denovo207422 + denovo207449 + denovo207501 +
                                 denovo207544 + denovo207599 + denovo207624 + denovo207661 + denovo207674 + denovo207678 + denovo207691 + denovo207699 +
                                 denovo207725 + denovo207736 + denovo207778 + denovo207793 + denovo207817 + denovo207855 + denovo207868 + denovo207878 +
                                 denovo207884 + denovo207909 + denovo207910 + denovo207917 + denovo207938 + denovo207976 + denovo207978 + denovo208008 +
                                 denovo208009 + denovo208023 + denovo208047 + denovo208096 + denovo208124 + denovo208125 + denovo208135 + denovo208224 +
                                 denovo208234 + denovo208357 + denovo208378 + denovo208405 + denovo208437 + denovo208449 + denovo208457 + denovo208459 +
                                 denovo208483 + denovo208500 + denovo208560 + denovo208588 + denovo208644 + denovo208651 + denovo208680 + denovo208700 +
                                 denovo208733 + denovo208743 + denovo208819 + denovo208846 + denovo208849 + denovo208867 + denovo208902 + denovo208956 +
                                 denovo208981 + denovo208995 + denovo209005 + denovo209013 + denovo209014 + denovo209054 + denovo209084 + denovo209106 +
                                 denovo209116 + denovo209124 + denovo209175 + denovo209191 + denovo209325 + denovo209331 + denovo209346 + denovo209353 +
                                 denovo209405 + denovo209447 + denovo209497 + denovo209501 + denovo209513 + denovo209526 + denovo209530 + denovo209592 +
                                 denovo209618 + denovo209632 + denovo209671 + denovo209714 + denovo209758 + denovo209763 + denovo209774 + denovo209781 +
                                 denovo209782 + denovo209841 + denovo209855 + denovo209874 + denovo209888 + denovo209929 + denovo210026 + denovo210029 +
                                 denovo210048 + denovo210051 + denovo210093 + denovo210166 + denovo210189 + denovo210190 + denovo210200 + denovo210206 +
                                 denovo210236 + denovo210261 + denovo210263 + denovo210308 + denovo210318 + denovo210352 + denovo210370 + denovo210375 +
                                 denovo210395 + denovo210413 + denovo210423 + denovo210440 + denovo210454 + denovo210456 + denovo210467 + denovo210473 +
                                 denovo210530 + denovo210536 + denovo210564 + denovo210575 + denovo210586 + denovo210587 + denovo210610 + denovo210630 +
                                 denovo210761 + denovo210779 + denovo210791 + denovo210855 + denovo210867 + denovo210883 + denovo210892 +
                                 denovo210936 + denovo210973 + denovo210989 + denovo211012 + denovo211014 + denovo211024 + denovo211032 + denovo211053 +
                                 denovo211116 + denovo211129 + denovo211137 + denovo211149 + denovo211157 + denovo211194 + denovo211208 + denovo211225 +
                                 denovo211263 + denovo211275 + denovo211287 + denovo211340 + denovo211393 + denovo211410 + denovo211426 +
                                 denovo211456 + denovo211459 + denovo211497 + denovo211500 + denovo211516 + denovo211519 + denovo211572 +
                                 denovo211599 + denovo211611 + denovo211650 + denovo211684 + denovo211745 + denovo211780 + denovo211839 + denovo211840 +
                                 denovo211862 + denovo211899 + denovo211900 + denovo211943 + denovo211955 + denovo211979 + denovo212022 + denovo212043 +
                                 denovo212062 + denovo212249 + denovo212253 + denovo212268 + denovo212281 + denovo212290 + denovo212299 + denovo212302 +
                                 denovo212321 + denovo212377 + denovo212394 + denovo212408 + denovo212419 + denovo212423 + denovo212425 +
                                 denovo212443 + denovo212467 + denovo212488 + denovo212524 + denovo212528 + denovo212532 + denovo212588 + denovo212607 +
                                 denovo212612 + denovo212656 + denovo212720 + denovo212728 + denovo212732 + denovo212768 + denovo212821 + denovo212830 +
                                 denovo212887 + denovo212888 + denovo212948 + denovo212951 + denovo212967 + denovo212980 + denovo212995 + denovo213006 +
                                 denovo213053 + denovo213188 + denovo213197 + denovo213205 + denovo213213 + denovo213222 + denovo213304 +
                                 denovo213311 + denovo213372 + denovo213380 + denovo213399 + denovo213416 + denovo213485 + denovo213508 + denovo213509 +
                                 denovo213537 + denovo213562 + denovo213576 + denovo213584 + denovo213592 + denovo213618 + denovo213625 + denovo213627 +
                                 denovo213708 + denovo213715 + denovo213719 + denovo213733 + denovo213817 + denovo213826 + denovo213850 + denovo213933 +
                                 denovo213955 + denovo213986 + denovo214061 + denovo214062 + denovo214135 + denovo214150 + denovo214153 +
                                 denovo214175 + denovo214197 + denovo214199 + denovo214222 + denovo214241 + denovo214299 + denovo214354 + denovo214373 +
                                 denovo214397 + denovo214466 + denovo214468 + denovo214473 + denovo214507 + denovo214508 + denovo214568 + denovo214579 +
                                 denovo214584 + denovo214660 + denovo214668 + denovo214685 + denovo214704 + denovo214719 + denovo214741 + denovo214778 +
                                 denovo214829 + denovo214837 + denovo214845 + denovo214853 + denovo214883 + denovo214885 + denovo214898 + denovo214912 +
                                 denovo214920 + denovo214962 + denovo214984 + denovo215052 + denovo215058 + denovo215072 + denovo215091 + denovo215095 +
                                 denovo215146 + denovo215213 + denovo215229 + denovo215257 + denovo215271 + denovo215297 + denovo215299 + denovo215306 +
                                 denovo215321 + denovo215336 + denovo215364 + denovo215371 + denovo215388 + denovo215454 + denovo215494 +
                                 denovo215513 + denovo215553 + denovo215602 + denovo215701 + denovo215712 + denovo215722 + denovo215731 + denovo215744 +
                                 denovo215755 + denovo215762 + denovo215772 + denovo215777 + denovo215797 + denovo215805 + denovo215832 + denovo215844 +
                                 denovo215863 + denovo215866 + denovo215903 + denovo215913 + denovo215937 + denovo215959 + denovo215983 + denovo216014 +
                                 denovo216015 + denovo216017 + denovo216030 + denovo216074 + denovo216075 + denovo216098 + denovo216185 +
                                 denovo216195 + denovo216238 + denovo216241 + denovo216249 + denovo216270 + denovo216273 + denovo216291 + denovo216301 +
                                 denovo216350 + denovo216366 + denovo216373 + denovo216385 + denovo216407 + denovo216461 + denovo216469 + denovo216471 +
                                 denovo216474 + denovo216475 + denovo216491 + denovo216498 + denovo216502 + denovo216508 + denovo216517 + denovo216528 +
                                 denovo216542 + denovo216582 + denovo216605 + denovo216626 + denovo216640 + denovo216699 + denovo216766 + denovo216788 +
                                 denovo216804 + denovo216869 + denovo216886 + denovo216902 + denovo216974 + denovo216975 + denovo216981 + denovo216998 +
                                 denovo217073 + denovo217075 + denovo217083 + denovo217142 + denovo217151 + denovo217175 + denovo217193 + denovo217265 +
                                 denovo217296 + denovo217312 + denovo217355 + denovo217356 + denovo217384 + denovo217403 + denovo217422 + denovo217443 +
                                 denovo217478 + denovo217483 + denovo217496 + denovo217506 + denovo217513 + denovo217589 + denovo217615 + denovo217672 +
                                 denovo217704 + denovo217711 + denovo217746 + denovo217778 + denovo217793 + denovo217811 + denovo217866 +
                                 denovo217906 + denovo217925 + denovo217995 + denovo218024 + denovo218046 + denovo218139 + denovo218141 + denovo218150 +
                                 denovo218159 + denovo218204 + denovo218229 + denovo218235 + denovo218282 + denovo218284 + denovo218304 +
                                 denovo218335 + denovo218370 + denovo218446 + denovo218473 + denovo218635 + denovo218671 + denovo218727 + denovo218758 +
                                 denovo218780 + denovo218790 + denovo218795 + denovo218838 + denovo218861 + denovo218875 + denovo218879 +
                                 denovo218895 + denovo218951 + denovo218996 + denovo219000 + denovo219004 + denovo219014 + denovo219015 +
                                 denovo219040 + denovo219051 + denovo219111 + denovo219219 + denovo219236 + denovo219253 + denovo219255 +
                                 denovo219272 + denovo219285 + denovo219341 + denovo219344 + denovo219378 + denovo219433 + denovo219453 + denovo219459 +
                                 denovo219463 + denovo219466 + denovo219507 + denovo219521 + denovo219537 + denovo219545 + denovo219571 + denovo219607 +
                                 denovo219636 + denovo219637 + denovo219704 + denovo219705 + denovo219721 + denovo219722 + denovo219781 + denovo219811 +
                                 denovo219823 + denovo219848 + denovo219851 + denovo219872 + denovo219965 + denovo220096 + denovo220126 +
                                 denovo220166 + denovo220195 + denovo220233 + denovo220258 + denovo220266 + denovo220325 + denovo220430 + denovo220476 +
                                 denovo220486 + denovo220535 + denovo220573 + denovo220585 + denovo220588 + denovo220625 + denovo220652 + denovo220654 +
                                 denovo220662 + denovo220667 + denovo220689 + denovo220768 + denovo220771 + denovo220810 + denovo220818 + denovo220822 +
                                 denovo220844 + denovo220855 + denovo220866 + denovo220879 + denovo220941 + denovo220956 + denovo220968 + denovo221034 +
                                 denovo221039 + denovo221058 + denovo221061 + denovo221068 + denovo221086 + denovo221099 + denovo221129 + denovo221140 +
                                 denovo221177 + denovo221182 + denovo221188 + denovo221216 + denovo221232 + denovo221239 + denovo221280 + denovo221293 +
                                 denovo221311 + denovo221337 + denovo221340 + denovo221342 + denovo221347 + denovo221371 + denovo221377 + denovo221405 +
                                 denovo221423 + denovo221486 + denovo221591 + denovo221653 + denovo221658 + denovo221668 + denovo221707 + denovo221725 +
                                 denovo221761 + denovo221788 + denovo221801 + denovo221813 + denovo221817 + denovo221870 + denovo221930 + denovo222007 +
                                 denovo222040 + denovo222049 + denovo222064 + denovo222117 + denovo222132 + denovo222143 + denovo222153 +
                                 denovo222160 + denovo222270 + denovo222357 + denovo222362 + denovo222366 + denovo222435 + denovo222447 +
                                 denovo222483 + denovo222528 + denovo222558 + denovo222560 + denovo222579 + denovo222623 + denovo222683 + denovo222703 +
                                 denovo222736 + denovo222739 + denovo222751 + denovo222783 + denovo222867 + denovo222882 + denovo222961 + denovo223014 +
                                 denovo223039 + denovo223069 + denovo223103 + denovo223120 + denovo223130 + denovo223135 + denovo223180 + denovo223183 +
                                 denovo223262 + denovo223321 + denovo223366 + denovo223381 + denovo223395 + denovo223401 + denovo223418 + denovo223477 +
                                 denovo223522 + denovo223528 + denovo223568 + denovo223606 + denovo223624 + denovo223644 + denovo223646 + denovo223686 +
                                 denovo223720 + denovo223723 + denovo223757 + denovo223778 + denovo223801 + denovo223815 + denovo223821 + denovo223837 +
                                 denovo223841 + denovo223842 + denovo223949 + denovo224012 + denovo224039 + denovo224088 + denovo224142 +
                                 denovo224166 + denovo224241 + denovo224242 + denovo224254 + denovo224261 + denovo224283 + denovo224291 + denovo224299 +
                                 denovo224308 + denovo224358 + denovo224369 + denovo224387 + denovo224404 + denovo224416 + denovo224420 + denovo224439 +
                                 denovo224484 + denovo224487 + denovo224552 + denovo224604 + denovo224609 + denovo224676 + denovo224681 +
                                 denovo224806 + denovo224812 + denovo224837 + denovo224839 + denovo224890 + denovo224894 + denovo224898 + denovo224942 +
                                 denovo224978 + denovo225028 + denovo225042 + denovo225131 + denovo225200 + denovo225227 + denovo225228 + denovo225327 +
                                 denovo225351 + denovo225358 + denovo225400 + denovo225431 + denovo225439 + denovo225444 + denovo225635 + denovo225652 +
                                 denovo225658 + denovo225673 + denovo225693 + denovo225711 + denovo225713 + denovo225720 + denovo225723 + denovo225725 +
                                 denovo225803 + denovo225804 + denovo225808 + denovo225855 + denovo225859 + denovo225890 + denovo225937 + denovo225956 +
                                 denovo225964 + denovo225990 + denovo226021 + denovo226026 + denovo226093 + denovo226098 + denovo226133 + denovo226135, data = presence_absence_2016_water_velocity, CV = TRUE)

DFA_predict <- table(presence_absence_2016_water_velocity$WaterVelocity, DFA_fit_JK$class)

# Measure Total Correctly Predicted Individuals by Water Velocity
diag(prop.table(DFA_predict, 1))
