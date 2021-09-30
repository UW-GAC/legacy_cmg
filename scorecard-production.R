### Code to generate GSPCC 'Mendeliant Traits by the Numbers'
### Steve Buyske
### changed by Ben Heavner to remove API key from script
### Production version 2016-08-28
## Lacks error checking, however

## Run regularly with
## R CMD BATCH scorecard-production.R

###Preliminaries

# if the lines below give an error message it should be fixed be running
# install.packages(c("personograph", "RColorBrewer", "stringr", "grImport"))

library(personograph)
library(RColorBrewer)
library(stringr)
library(grImport)



cols <- c("#A0CADD", "#2399DA", "#F5C58F", "#F38018", "#D4B4BB", "#D35F75")

total.number.of.genes <- 19580  ### From the UW algorithm, but I don't know where they got it

# currently getting key from system environment to keep it out of the
# github repo or any R history or R scripts.
# can be set in .Rprofile file or in .Renviron file
# There are other ways to manage API keys, too.
omim.key <- Sys.getenv("OMIM_KEY")

square <- readPicture("square.ps.xml") # only reason package grImport is needed


clip.comment <- function(x) gsub("# ", "", x)

combine.f <- function(x){
    x <- as.data.frame(x, stringsAsFactors = FALSE)
    names(x) <- names(mim.df)  # hard-coded to mim.df
    x.list <- split(x, apply(x[, c("prefix", "phenotype.mapping.key", "isComplex")], 1, paste, collapse = ":"))
    do.call("rbind", lapply(x.list, combine.mini.f))
 }

combine.mini.f <- function(xx){
    xx <- as.data.frame(xx, stringsAsFactors = FALSE)
    names(xx) <- names(mim.df)  # hard-coded to mim.df
    result <- xx[1, ]
    result["preferred.title.symbol"] <- paste(unique(xx[, "preferred.title.symbol"]), collapse = " | ")
    result["phenoname"] <- paste(unique(xx[, "phenoname"]), collapse = " | ")
    result
}



### Get files from OMIM and massage them
# note: have to build URLs b/c they include the omim API key
mimTitlesURL <- paste0("http://data.omim.org/downloads/", omim.key, "/mimTitles.txt")
mimTitleswget <- paste0("wget -N ", mimTitlesURL)
system(mimTitleswget)

morbidmapURL <- paste0("http://data.omim.org/downloads/", omim.key, "/morbidmap.txt")
morbidmapwget <- paste0("wget -N ", morbidmapURL)
system(morbidmapwget)

system("wget -N http://omim.org/static/omim/data/mim2gene.txt")

mimtitles <- read.table("mimTitles.txt", sep = "\t", header = FALSE, as.is = TRUE, fill = TRUE, quote = "")
morbidmap <- read.table("morbidmap.txt", sep = "\t", header = FALSE, as.is = TRUE, fill = TRUE, quote = "")
mim2gene <- read.table("mim2gene.txt", sep = "\t", header = FALSE, as.is = TRUE, fill = TRUE, quote = "")


names(mimtitles) <- gsub("\\.+", ".", tolower(make.names(clip.comment(read.table("mimTitles.txt", sep = "\t", comment.char = "", skip = 2, nrows = 1, as.is = TRUE)[1,]))))
names(morbidmap) <- gsub("\\.+", ".", tolower(make.names(clip.comment(read.table("morbidmap.txt", sep = "\t", comment.char = "", skip = 3, nrows = 1, as.is = TRUE)[1,]))))
names(mim2gene) <- make.names(c("mim.number", "MIM.Entry.Type", "Entrez Gene ID", "Approved Gene Symbol", "Ensembl Gene ID"))

names(morbidmap)[names(morbidmap) == "mim.number"] <- "mim.number.orig"
names(morbidmap)[names(morbidmap) == "phenotype"] <- "phenoname"  # just to reduce confusion


## extract phenotype mapping key from phenoname
## map e.g., ... (3) ...  to 3 . Takes the last if there are more than one. If there are, the early ones appear to be part of the descriptions
morbidmap$phenotype.mapping.key <- as.integer(gsub("[()]", "", sapply(str_extract_all(morbidmap$phenoname, "\\(([0-9])\\)", simplify = FALSE), function(x) x[length(x)])))


## extract phenotype MIM number from phenoname
temp <- sapply(strsplit(morbidmap$phenoname, ","), function(x) x[length(x)])

## a few pheno mim numbers are followed by a comma
temp2 <- as.integer(sapply(strsplit(morbidmap$phenoname, ","), function(x) x[length(x) - 1]))


morbidmap$mim.number.pheno <- sapply(strsplit(temp, " +"), function(x) x[length(x) - 1])
morbidmap$mim.number.pheno[nchar(morbidmap$mim.number.pheno) < 6] <- NA # junk from descriptions
morbidmap$mim.number.pheno <- as.integer(morbidmap$mim.number.pheno)
morbidmap$mim.number.pheno[!is.na(temp2) & nchar(temp2) > 3] <- temp2[!is.na(temp2) & nchar(temp2) > 3]
morbidmap$mim.number.pheno[is.na(morbidmap$mim.number.pheno) & morbidmap$phenotype.mapping.key %in% 4] <- morbidmap$mim.number.orig[is.na(morbidmap$mim.number.pheno) & morbidmap$phenotype.mapping.key %in% 4]
# if a chromosomal del/dup syndrome, the gene MIM *is usually* the phenotype MIM



### Create master file of phenotypes

mim.df <- merge(mimtitles, mim2gene[, c("mim.number", "MIM.Entry.Type")], all = TRUE)
mim.df <- merge(mim.df, morbidmap[!is.na(morbidmap$mim.number.pheno),], by.x = "mim.number", by.y = "mim.number.pheno", all = TRUE)
mim.df[is.na(mim.df$phenoname), "phenoname"] <- morbidmap[match(mim.df[is.na(mim.df$phenoname), "mim.number"], morbidmap[, "mim.number.orig"]), "phenoname"]

temp <- mim.df$mim.number.orig
temp[is.na(temp)] <- mim.df$mim.number[is.na(temp)]
mim.df$mim.number.joint <- apply(cbind(mim.df[, "mim.number"], temp, mim.df$phenotype.mapping.key), 1, paste, collapse = "-")

## isComplex:
# if phenoname text has
# somatic, then isComplex <- "somatic"
# carcinoma, tumor, cancer, leukemia, sarcoma, blastoma, adenoma, cytoma, myelodysplastic, myelofibrosis, oma (at end of word) as well as somatic
# then isComplex <= "cancer"
# none of the above, but has 
# risk, quantitative trait locus, qtl, suscep* to, [, {
# then isComplex <- "yes"
# otherwise isComplex <- "no"

## Primary difference from Jessica Chong's UW-CMG code is that that code used only phenoname and not preferred.title.symbol; we use both

temp1 <- mim.df[, "preferred.title.symbol"]
temp2 <- mim.df[, "phenoname"]
temp1[is.na(temp1)] <- ""
temp2[is.na(temp2)] <- ""

temp <- apply(cbind(temp1, temp2), 1, paste, collapse = ";")

mim.df$isComplex <- "no"

mim.df$isComplex[unique(c(grep("risk", temp, ignore.case = TRUE), 
	intersect(intersect(grep("quantitative", temp, ignore.case = TRUE), grep("trait", temp, ignore.case = TRUE)), grep("locus", temp, ignore.case = TRUE)),
	grep("qtl", temp, ignore.case = TRUE),
	grep("suscep[a-z]+ to", temp, ignore.case = TRUE), # catch typos for susceptability without catching too much else
	grep("susceptability", temp, ignore.case = TRUE),
	grep("\\[", temp, ignore.case = TRUE),
	grep("\\{", temp, ignore.case = TRUE)
))] <- "yes"

mim.df$isComplex[intersect(
    which(mim.df$isComplex != "yes"),
	grep("somatic", temp, ignore.case = TRUE))] <- "somatic"
	


mim.df$isComplex[intersect(unique(c(
    grep("tumor", temp, ignore.case = TRUE),
	grep("carcinoma", temp, ignore.case = TRUE),
	grep("cancer", temp, ignore.case = TRUE),
	grep("leukemia", temp, ignore.case = TRUE),
	grep("sarcoma", temp, ignore.case = TRUE),
	grep("blastoma", temp, ignore.case = TRUE),
	grep("adenoma", temp, ignore.case = TRUE),
	grep("cytoma", temp, ignore.case = TRUE),
	grep("myelodysplastic", temp, ignore.case = TRUE),
	grep("Myelofibrosis", temp, ignore.case = TRUE),
	grep("oma[, ]", temp, ignore.case = TRUE))),
	grep("somatic", temp, ignore.case = TRUE))] <- "cancer"
	
mim.list <- split(mim.df, mim.df$mim.number.joint)



### mim2.df is a version of mim.df with redundancy removed    	
	
mim2.df <- do.call("rbind",lapply(mim.list, combine.f))
row.names(mim2.df) <- 1:nrow(mim2.df)


mim2.reduced.df <- mim2.df[with(mim2.df, (isComplex %in% c("no", "somatic")) & (is.na(phenotype.mapping.key) | phenotype.mapping.key != 4)), ]

### fig1

# Number sign and is.complex is no or missing; looking among unique mim.numbers
fig1.num <- with(unique(mim2.reduced.df[, c("mim.number.joint","prefix", "isComplex")]), sum(prefix %in% "Number Sign"))

# Number sign, NULL, or percent sign and is.complex is no or missing; looking among unique mim.numbers
fig1.denom <- with(unique(mim2.reduced.df[, c("mim.number.joint","prefix", "isComplex")]), sum(prefix %in% c("Number Sign", "NULL", "Percent", "Plus")))



### fig2

# Number of unique mim number for locus among entries with number sign and is.complex is no or missing; looking among all entries in mim2gene
fig2.num <- length(unique(mim2.reduced.df[mim2.reduced.df$prefix %in% "Number Sign", "mim.number.orig"]))
fig2.denom <- total.number.of.genes 


### fig3

# number of above where the number shows up more than once

fig3.num <- sum(table(mim2.reduced.df[mim2.reduced.df$prefix %in% "Number Sign", "mim.number.orig"])>1)
fig3.denom <- fig2.num

### Created figures and captions

png("fig1.png")
personograph(list(first = 1 - fig1.num/fig1.denom, second = fig1.num/fig1.denom), colors = list(first = cols[1], second = cols[2]), icon = square, draw.legend = FALSE)
dev.off()

png("fig2.png")
personograph(list(first = 1 - fig2.num/fig2.denom, second = fig2.num/fig2.denom), colors = list(first = cols[1 + 2], second = cols[2 + 2]), icon = square, draw.legend = FALSE)
dev.off()

png("fig3.png")
personograph(list(first = 1 - fig3.num/fig3.denom, second = fig3.num/fig3.denom), colors = list(first = cols[1 + 4], second = cols[2 + 4]), icon = square, draw.legend = FALSE)
dev.off()


cat(file = "head_text.txt", "About Mendelian Conditions\n
A phenotype is the collection of observable or measurable traits of an individual. Phenotypes that result from changes (i.e., variants) in a single gene (i.e., monogenic) and that can be transmitted from parents to offspring in Mendelian patterns, such as autosomal dominant, autosomal recessive, X-linked, are known as Mendelian phenotypes.\n\nWhat is known about the genomic basis of Mendelian phenotypes?"
)


cat(file = "caption_fig1.txt", 
    round(fig1.num/fig1.denom * 100), "% (", 
    prettyNum(fig1.num, big.mark = ","), "/",  prettyNum(fig1.denom, big.mark = ","), 
    ") is the fraction of Mendelian phenotypes for which the underlying gene is known."
, sep = "")


cat(file = "caption_fig2.txt", 
    round(fig2.num/fig2.denom * 100), "% (",
    prettyNum(fig2.num, big.mark = ","), "/", prettyNum(fig2.denom, big.mark = ","),
    ") of human genes are known to underlie a Mendelian phenotype."
, sep = "")


cat(file = "caption_fig3.txt", 
    round(fig3.num/fig3.denom * 100), "% (",
    prettyNum(fig3.num, big.mark = ","), "/", prettyNum(fig3.denom, big.mark = ","),
    ") of genes known to underlie a Mendelian phenotype cause two or more different phenotypes."
, sep = "")


today.date <- format(Sys.time(), "%Y-%m-%d")

cat(file = "foot_text.txt", "Counts are based on calculations by the Genome Sequencing Program Coordinating Center using data extracted from OMIM (http://www.omim.org) on ", today.date, ".", sep = "")


### Write file of historic data. Older files can be purged, but the code to do so isn't included.

temp.df <- data.frame(stringsAsFactors = FALSE, fig1.num, fig1.denom, fig2.num, fig2.denom, fig3.num, fig3.denom, calc.date = today.date)

historic.files <- system("ls historic-data*.txt", intern = TRUE)
if (length(historic.files) == 0) {
    historic.df <- temp.df } else {
    historic.df <- unique(rbind(read.table(file = historic.files[length(historic.files)], header = TRUE, as.is = TRUE, sep = "\t"), temp.df))
}

write.table(historic.df, sep = "\t", file = paste0("historic-data-", format(Sys.time(), "%Y-%m-%d"), ".txt"), row.names = FALSE, col.names = TRUE, quote = FALSE)

### Save files from OMIM in case we change algorithm and need historic data. Assumes directory "./omim.archive" exists
# system("mkdir omim.archive")
system(paste0("cp mimTitles.txt omim.archive/", today.date, "-mimTitles.txt"))
system(paste0("cp morbidmap.txt omim.archive/", today.date, "-morbidmap.txt"))
system(paste0("cp mim2gene.txt omim.archive/", today.date, "-mim2gene.txt"))


