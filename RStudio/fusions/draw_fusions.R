library(dplyr)
library(plyr)
library(Rsamtools)
library(maftools)
SURFME <- openxlsx::read.xlsx(file.path(data, "SURFME_v2023.xlsx"))
  
  
alignmentsFiles <- list()
alignmentsFiles <- list.files(file.path("A:","GBM", "results", "05_star"),
                             pattern = "sortedByCoord.out.bam")
alignmentsFiles <- alignmentsFiles[!grepl("2DH",alignmentsFiles)]

alignmentsFiles_593 <- alignmentsFiles[grep("593",alignmentsFiles)]

VI_3429_593_2DN_fusions <- rbind(fusion.df[["VI-3429-593-2DN-1-star-fusion"]], 
                                 fusion.df[["VI-3429-593-2DN-2-star-fusion"]], 
                                 fusion.df[["VI-3429-593-2DN-3-star-fusion"]]) %>% 
  dplyr::group_by(X.FusionName, LeftBreakpoint, RightBreakpoint) %>%
  dplyr::summarise(across(where(is.numeric), \(x) mean(x, na.rm = T)),
            across(where(is.character), \(x) toString(unique(x)))) %>% 
  dplyr::ungroup() %>% 
  dplyr::mutate(
    JunctionReadCount = ceiling(JunctionReadCount),
    SpanningFragCount = ceiling(SpanningFragCount)
  )

VI_3429_593_2DN_fusions <- prep_fusions(VI_3429_593_2DN_fusions)
VI_3429_593_2DN_fusions <- VI_3429_593_2DN_fusions %>% 
  dplyr::filter(transcript_id1 != "." & transcript_id2 != "." &
                  stringr::str_detect(condition, ",")) %>% 
  dplyr::rowwise() %>% 
  dplyr::filter(any(grepl(gene1, SURFME$Gene.names1)) | 
                  any(grepl(gene2, SURFME$Gene.names1)))

VI_3429_593_fusions <- do.call(rbind, fusion.df[1:7]) %>% 
  dplyr::group_by(X.FusionName, LeftBreakpoint, RightBreakpoint) %>%
  dplyr::summarise(across(where(is.numeric), \(x) mean(x, na.rm = T)),
            across(where(is.character), \(x) toString(unique(x)))) %>% 
  dplyr::ungroup() %>% 
  dplyr::mutate(
    JunctionReadCount = ceiling(JunctionReadCount),
    SpanningFragCount = ceiling(SpanningFragCount)
  )
VI_3429_593_fusions <- prep_fusions(VI_3429_593_fusions)
openxlsx::write.xlsx(VI_3429_593_fusions, file.path(date, results_dir, "VI_3429_593_fusions.xlsx"))

VI_3429_593_fusions <- VI_3429_593_fusions %>% 
  dplyr::filter(transcript_id1 != "." & transcript_id2 != "." &
                  stringr::str_detect(condition, ",")) %>% 
  dplyr::rowwise() %>% 
  dplyr::filter(any(grepl(gene1, SURFME$Gene.names1)) | 
                  any(grepl(gene2, SURFME$Gene.names1)))

VI_3429_673_fusions <- do.call(rbind, fusion.df[8:14]) %>% 
  dplyr::group_by(X.FusionName, LeftBreakpoint, RightBreakpoint) %>%
  dplyr::summarise(across(where(is.numeric), \(x) mean(x, na.rm = T)),
                   across(where(is.character), \(x) toString(unique(x)))) %>% 
  dplyr::ungroup() %>% 
  dplyr::mutate(
    JunctionReadCount = ceiling(JunctionReadCount),
    SpanningFragCount = ceiling(SpanningFragCount)
  )
VI_3429_673_fusions <- prep_fusions(VI_3429_673_fusions)
openxlsx::write.xlsx(VI_3429_673_fusions, file.path(date, results_dir, "VI_3429_673_fusions.xlsx"))

VI_3429_673_fusions <- VI_3429_673_fusions %>% 
  dplyr::filter(transcript_id1 != "." & transcript_id2 != "." &
                  stringr::str_detect(condition, ",")) %>% 
  dplyr::rowwise() %>% 
  dplyr::filter(any(grepl(gene1, SURFME$Gene.names1)) | 
                  any(grepl(gene2, SURFME$Gene.names1)))

# read cytoband annotation
install.packages("D3GB")
library(D3GB)
cytobands <- GRCh38.bands
cytobands <- cytobands[order(cytobands$chr, cytobands$start, cytobands$end),]%>%
  mutate_if(sapply(cytobands, is.character), as.factor) %>% 
  setNames(c("contig","start","end","band","giemsa"))

# read exon annotation
exons <- scan(file.path("..","..","refs","ref_annot.gtf"), 
              what=list(contig="",src="",type="",start=0,end=0,score="",strand="",frame="",attributes=""), 
              sep="\t", comment.char="#", quote='"', multi.line=F)
attr(exons, "row.names") <- .set_row_names(length(exons[[1]]))
class(exons) <- "data.frame"

exons <- exons[exons$type %in% c("exon","CDS"),c("contig","type","start","end","strand","attributes")]
exons$contig <- removeChr(exons$contig)
exons$geneID <- parseGtfAttribute("gene_id", exons)
exons$geneName <- parseGtfAttribute("gene_name", exons)
exons$geneName <- ifelse(exons$geneName == "", exons$geneID, exons$geneName)
exons$transcript <- parseGtfAttribute("transcript_id", exons)
exons$transcript <- gsub("\\.*", "", exons$transcript)
exonsnum <- exons %>% 
  dplyr::group_by(geneID) %>%
  dplyr::filter(type == "exon") %>%
  dplyr::summarise(exonNumber = n())
exons <- exons %>%
  dplyr::left_join(exonsnum, by = c("geneID"))

ProteinDomains <- list()
ProteinDomains$PCOLCE2 <- scan(file.path("data/uniprot/PCOLCE2.gff3"), 
                   what=list(id="",src="",type="",start=0,end=0,score="",strand="",frame="",attributes=""), 
                   sep="\t", comment.char="", quote="", multi.line=F)
ProteinDomains$ATR <- scan(file.path("data/uniprot/ATR.gff3"), 
                               what=list(id="",src="",type="",start=0,end=0,score="",strand="",frame="",attributes=""), 
                               sep="\t", comment.char="", quote="", multi.line=F)
ProteinDomains$SCAP <- scan(file.path("data/uniprot/SCAP.gff3"), 
                               what=list(id="",src="",type="",start=0,end=0,score="",strand="",frame="",attributes=""), 
                               sep="\t", comment.char="", quote="", multi.line=F)
ProteinDomains$SLC6A20 <- scan(file.path("data/uniprot/SLC6A20.gff3"), 
                               what=list(id="",src="",type="",start=0,end=0,score="",strand="",frame="",attributes=""), 
                               sep="\t", comment.char="", quote="", multi.line=F)

ProteinDomains <- lapply(ProteinDomains, function(x){
  return(as.data.frame(x))
})
ProteinDomains <- do.call(rbind.fill, ProteinDomains)
ProteinDomains$contig <- parseGtfAttribute("Chromosome", ProteinDomains)
ProteinDomains$start <- parseGtfAttribute("GenomeLocStart", ProteinDomains)
ProteinDomains$end <- parseGtfAttribute("GenomeLocEnd", ProteinDomains)

mafDomains <-  system.file("extdata", "protein_domains.RDs", package = "maftools")
mafDomains <- readRDS(file = mafDomains)
data.table::setDT(x = mafDomains)

protein_colours <- c("#4477AA", "#228833", "#AA3377", "#BBBBBB", "#66CCEE",
                    "#CCBB44", "#EE6677", "#0077BB", "#EE7733", "#33BBEE",
                    "#CC3311", "#009988", "#EE3377", "#999933",  "#FE6100",
                    "#785EF0", "#FFB000", "#63ACBE", "#648FFF")
ProteinDomains$color <- protein_colours[as.factor(ProteinDomains$type)]

readCoverage <- function(alignmentsFile, contig, coverageRegion) {
  coverageData <- tryCatch(
    {
      alignments <- readGAlignments(alignmentsFile, param=ScanBamParam(which=GRanges(contig, coverageRegion)))
      coverage(alignments)[[contig]]
    },
    error=function(e) {
      alignments <- readGAlignments(alignmentsFile, param=ScanBamParam(which=GRanges(addChr(contig), coverageRegion)))
      coverage(alignments)[[addChr(contig)]]
    }
  )
  if (exists("alignments")) rm(alignments)
  return(coverageData)
}

coverage <- list()
regions <- list()

library(GenomicAlignments)
for (i in 1:nrow(VI_3429_593_2DN_fusions)) {
  
  message(paste0("Drawing fusion #", i, ": ", VI_3429_593_2DN_fusions[i,"gene1"], ":", VI_3429_593_2DN_fusions[i,"gene2"]))
  
	# find exons near breakpoint for gene 1
	exon1 <- exons[exons$geneID %in% VI_3429_593_2DN_fusions[i,"gene_id1"], ] # find closest exon
	gene1 <- IRanges(min(exon1$start), max(exon1$end))
	region1 <- IRanges(min(start(gene1), as.numeric(VI_3429_593_2DN_fusions[i,"breakpoint1"]-showVicinity[1])), 
	                   max(end(gene1),  as.numeric(VI_3429_593_2DN_fusions[i,"breakpoint1"]+showVicinity[2])))
		
	# find exons near breakpoint for gene 1
	exon2 <- exons[exons$geneID %in% VI_3429_593_2DN_fusions[i,"gene_id2"], ] # find closest exon
	gene2 <- IRanges(min(exon2$start), max(exon2$end))
	region2 <- IRanges(min(start(gene2), as.numeric(VI_3429_593_2DN_fusions[i,"breakpoint2"]-showVicinity[1])),
	                   max(end(gene2),  as.numeric(VI_3429_593_2DN_fusions[i,"breakpoint2"]+showVicinity[2])))
		
	# get coverage track
  coverage1 <- readCoverage(alignmentsFile, VI_3429_593_2DN_fusions[i,"contig1"], region1)
  coverage2 <- readCoverage(alignmentsFile, VI_3429_593_2DN_fusions[i,"contig2"], region2)
  # shrink coverage range to chromosome boundaries to avoid subscript out of bounds errors
  region1 <- IRanges(max(start(region1), min(start(coverage1))), min(end(region1), max(end(coverage1))))
  region2 <- IRanges(max(start(region2), min(start(coverage2))), min(end(region2), max(end(coverage2))))
  
  coverage[[VI_3429_593_2DN_fusions[["X.FusionName"]][i]]] <- list(coverage1, coverage2)
  regions[[VI_3429_593_2DN_fusions[["X.FusionName"]][i]]] <- list(region1, region2)

  rm(list = c("exon1","gene1","region1", "coverage1",
              "exon2","gene2","region2", "coverage2"))
  }


exons1 <- list()
exons2 <- list()
for (i in 1:nrow(VI_3429_593_2DN_fusions)) {
  # find all exons belonging to the fused genes
  coverages = coverage[[VI_3429_593_2DN_fusions[["X.FusionName"]][i]]]
	e1 <- findExons(exons = exons, contig = VI_3429_593_2DN_fusions[["contig1"]][i],
	                    geneID = VI_3429_593_2DN_fusions[["gene_id1"]][i], 
	                    direction =  VI_3429_593_2DN_fusions[["direction1"]][i], 
	                    breakpoint = VI_3429_593_2DN_fusions[["breakpoint1"]][i], 
	                    coverage = coverages[[1]],
	                    transcriptId =  VI_3429_593_2DN_fusions[["transcript_id1"]][i])

	exons1[[VI_3429_593_2DN_fusions[["X.FusionName"]][i]]] <- e1
	
	e2 <- findExons(exons = exons, contig = VI_3429_593_2DN_fusions[["contig2"]][i],
	                geneID = VI_3429_593_2DN_fusions[["gene_id2"]][i], 
	                direction =  VI_3429_593_2DN_fusions[["direction2"]][i], 
	                breakpoint = VI_3429_593_2DN_fusions[["breakpoint2"]][i], 
	                coverage = coverages[[2]],
	                transcriptId =  VI_3429_593_2DN_fusions[["transcript_id2"]][i])
	
	exons2[[VI_3429_593_2DN_fusions[["X.FusionName"]][i]]] <- e2
	rm(list = c("e1","e2","coverages"))
}


for (i in 1:nrow(VI_3429_593_2DN_fusions)) {
	# normalize coverage
  coverageNormalization <- function(coverage, coverageRegion, exons) {
    max(1, 
        max(as.numeric(coverage[IRanges(sapply(exons$start,max,min(start(coverage))),
                                        sapply(exons$end,min,max(end(coverage))))]
                       )))
  }
  
	covnorm1 <- coverageNormalization(coverage[[i]][[1]], regions[[i]][[1]], exons1[[1]])
	covnorm2 <- coverageNormalization(coverage[[i]][[2]], regions[[i]][[2]], exons2[[1]])
	
	
	coverage[[i]][[1]] <- coverage[[i]][[1]]/covnorm1
	coverage[[i]]["normFactor1"] <- covnorm1
	
	coverage[[i]][[2]] <- coverage[[i]][[2]]/covnorm2
	coverage[[i]]["normFactor2"] <- covnorm2
}


# sort coding exons last, such that they are drawn over the border of non-coding exons
exons1 <- lapply(exons1, function(x) x[order(x$start, -rank(x$type)),])
exons2 <- lapply(exons2, function(x) x[order(x$start, -rank(x$type)),])

shift_coords <- function(exons, breakpoint){
  exons$left <- exons$start
  exons$right <- exons$end
  
  exons$right <- exons$right - min(exons$left)
  exons$breakpoint <- breakpoint - min(exons$left)
  exons$left <- exons$left - min(exons$left)
  return(exons)
}

exons1$`PCOLCE2--ATR` <- shift_coords(exons = exons1$`PCOLCE2--ATR`,
                                      breakpoint = VI_3429_593_2DN_fusions$breakpoint1[1])
exons2$`PCOLCE2--ATR` <- shift_coords(exons = exons2$`PCOLCE2--ATR`,
                                      breakpoint = VI_3429_593_2DN_fusions$breakpoint2[1])

exons1$`SCAP--SLC6A20` <- shift_coords(exons = exons1$`SCAP--SLC6A20`,
                                      breakpoint = VI_3429_593_2DN_fusions$breakpoint1[2])
exons2$`SCAP--SLC6A20` <- shift_coords(exons = exons2$`SCAP--SLC6A20`,
                                      breakpoint = VI_3429_593_2DN_fusions$breakpoint2[2])

# scale exon sizes to fit on page
scalingFactor <- list()
scalingFactor$`PCOLCE2--ATR` <- max(exons1$`PCOLCE2--ATR`$right) + max(exons2$`PCOLCE2--ATR`$right)
scalingFactor$`SCAP--SLC6A20` <- max(exons1$`SCAP--SLC6A20`$right) + max(exons2$`SCAP--SLC6A20`$right)


exons1$`PCOLCE2--ATR`$left <- exons1$`PCOLCE2--ATR`$left / scalingFactor$`PCOLCE2--ATR`
exons1$`PCOLCE2--ATR`$right <- exons1$`PCOLCE2--ATR`$right / scalingFactor$`PCOLCE2--ATR`
exons1$`PCOLCE2--ATR`$breakpoint <- exons1$`PCOLCE2--ATR`$breakpoint / scalingFactor$`PCOLCE2--ATR`

exons2$`PCOLCE2--ATR`$left <- exons2$`PCOLCE2--ATR`$left / scalingFactor$`PCOLCE2--ATR`
exons2$`PCOLCE2--ATR`$right <- exons2$`PCOLCE2--ATR`$right / scalingFactor$`PCOLCE2--ATR`
exons2$`PCOLCE2--ATR`$breakpoint <- exons2$`PCOLCE2--ATR`$breakpoint / scalingFactor$`PCOLCE2--ATR`


fusion_POLCE2_ATR <- list(
  region1 = regions$`PCOLCE2--ATR`[[1]],
  region2 = regions$`PCOLCE2--ATR`[[2]],
  coverage1 = coverage$`PCOLCE2--ATR`[[1]],
  coverage2 = coverage$`PCOLCE2--ATR`[[2]],
  exon1 = exons1$`PCOLCE2--ATR`,
  exon2 = exons2$`PCOLCE2--ATR`,
  normFactor1 = coverage$`PCOLCE2--ATR`$normFactor1,
  normFactor2 = coverage$`PCOLCE2--ATR`$normFactor2,
  scalingFactor = scalingFactor$`PCOLCE2--ATR`
)
fusion_SCAP_SLC6A20 <- list(
  region1 = regions$`SCAP--SLC6A20`[[1]],
  region2 = regions$`SCAP--SLC6A20`[[2]],
  coverage1 = coverage$`SCAP--SLC6A20`[[1]],
  coverage2 = coverage$`SCAP--SLC6A20`[[2]],
  exon1 = exons1$`SCAP--SLC6A20`,
  exon2 = exons2$`SCAP--SLC6A20`,
  normFactor1 = coverage$`SCAP--SLC6A20`$normFactor1,
  normFactor2 = coverage$`SCAP--SLC6A20`$normFactor2,
  scalingFactor = scalingFactor$`SCAP--SLC6A20`
)

pdfAttrs = list( "pdfWidth" = 11.692, "pdfHeight" = 8.267, 
                 "color1" = "#e5a5a5", "color2" = "#a7c4e5", 
                 "fontSize" =  1,"fontFamily" = "Helvetica")
pdfAttrs$darkColor1 <- getDarkColor(pdfAttrs$color1)
pdfAttrs$darkColor2 <- getDarkColor(pdfAttrs$color2)
circosColors <- c(translocation="#000000", duplication="#00bb00", deletion="#ff0000", inversion="#0000ff")

plot_fusion <- function(fusion_data, table, filename, title){
  
  pdf(file.path(date, plots_dir, "fusion_SCAP_SLC6A20.pdf"), onefile=T,
      width = pdfAttrs[["pdfWidth"]], height= pdfAttrs[["pdfHeigth"]],
      title="SCAP--SLC6A20 fusion")
  par(family=pdfAttrs[["fontFamily"]])
  
  # shift gene2 to the right border of the page
	gene2Offset <- 1 + 0.05 - max(fusion_data$exon2$right)
	
	# plot gene fusion
	breakpoint1 <- unique(fusion_data$exon1[["breakpoint"]])
	breakpoint2 <- unique(fusion_data$exon2[["breakpoint"]])
	
	# center fusion horizontally
	fusionOffset1 <- (max(fusion_data$exon1$right)+gene2Offset)/2 - ifelse(table[["direction1"]] == "downstream", fusion_data$exon1$breakpoint, max(fusion_data$exon1$right)-fusion_data$exon1$breakpoint)
	fusionOffset2 <- fusionOffset1 + ifelse(table[["direction1"]] == "downstream", fusion_data$exon1$breakpoint, max(fusion_data$exon1$right)-fusion_data$exon1$breakpoint)

	# layout: fusion on top, circos plot on bottom left, protein domains on bottom center, statistics on bottom right
	layout(matrix(c(1,1,1,2,4,5,3,4,5), 3, 3, byrow=TRUE), widths=c(1.1, 1.2, 0.7), heights=c(1.55, 1.2, 0.25))
	par(mar=c(0, 0, 0, 0))
	plot(0, 0, type="l", xlim=c(-0.12, 1.12), ylim=c(0.4, 1.1), bty="n", xaxt="n", yaxt="n", xlab="", ylab="")

	# vertical coordinates of layers
	ySampleName <- 1.04
	yIdeograms <- ifelse(alignmentsFile != "", 0.94, 0.84)
	yBreakpointLabels <- ifelse(alignmentsFile != "", 0.86, 0.76)
	yCoverage <- 0.72
	yExons <- 0.67
	yGeneNames <- 0.58
	yFusion <- 0.5
	yTranscript <- 0.45
	yScale <- 0.407
	yTrajectoryBreakpointLabels <- yBreakpointLabels - 0.035
	yTrajectoryExonTop <- yExons + 0.03
	yTrajectoryExonBottom <- yExons - 0.055
	yTrajectoryFusion <- yFusion + 0.03

	# print sample name (title of page)
	text(0.5, ySampleName, "SCAP--SLC6A20 fusion", font=2, cex=pdfAttrs[["fontSize"]]*1.5, adj=c(0.5,0))

	# draw ideograms
	if (!is.null(cytobands)) {
		drawIdeogram("left", min(fusion_data$exon1$left), max(fusion_data$exon1$right), 
		             yIdeograms, cytobands, table[["contig1"]], table[["breakpoint1"]])
		drawIdeogram("right", gene2Offset, gene2Offset+ max(fusion_data$exon2$right),
		             yIdeograms, cytobands, table[["contig2"]], table[["breakpoint2"]])
	}

	# draw gene & transcript names
	if (table[["gene1"]] != ".") {
	  text(max(fusion_data$exon1$right)/2, yGeneNames, table[["gene1"]],
	       font=2, cex=pdfAttrs[["fontSize"]], adj=c(0.5,0))
	}
	if (table[["site1"]] != "intergenic") {
	  text(max(fusion_data$exon1$right)/2, yGeneNames-0.01, head(fusion_data$exon1$transcript,1), 
	       cex=0.9*pdfAttrs[["fontSize"]], adj=c(0.5,1))
	}
	if (table[["gene2"]] != ".") {
	  text(gene2Offset+max(fusion_data$exon2$right)/2, yGeneNames, table[["gene2"]],
	       font=2, cex=pdfAttrs[["fontSize"]], adj=c(0.5,0))
	}
	if (table[["site2"]] != "intergenic") {
	  text(gene2Offset+max(fusion_data$exon2$right)/2, yGeneNames-0.01, 
	       head(fusion_data$exon2$transcript,1), 
	       cex=0.9*pdfAttrs[["fontSize"]], adj=c(0.5,1))
	}

	# label breakpoints
	text(breakpoint1+0.01, yBreakpointLabels-0.03, 
	     paste0("breakpoint1\n", table[["display_contig1"]], ":", table[["breakpoint1"]]),
	     adj=c(1,0), cex=pdfAttrs[["fontSize"]])
	text(gene2Offset+breakpoint2-0.01, yBreakpointLabels-0.03,
	     paste0("breakpoint2\n", table[["display_contig2"]], ":", table[["breakpoint2"]]),
	     adj=c(0,0), cex=pdfAttrs[["fontSize"]])

	# draw coverage axis
	if (alignmentsFile != "") {
		# left axis (gene1)
		lines(c(-0.02, -0.01, -0.01, -0.02), c(yCoverage, yCoverage, yCoverage+0.1, yCoverage+0.1))
		text(-0.025, yCoverage, "0", adj=c(1,0.5), cex=0.9*pdfAttrs[["fontSize"]])
		text(-0.025, yCoverage+0.1, fusion_data$normFactor1, adj=c(1,0.5), cex=0.9*pdfAttrs[["fontSize"]])
		text(-0.05, yCoverage+0.08, "Coverage", srt=90, cex=0.9*pdfAttrs[["fontSize"]], adj=c(1,0.5))

		rightCoverageAxisX <- gene2Offset+max(fusion_data$exon2$right)
		lines(c(rightCoverageAxisX+0.02, rightCoverageAxisX+0.01, rightCoverageAxisX+0.01, rightCoverageAxisX+0.02), c(yCoverage, yCoverage, yCoverage+0.1, yCoverage+0.1))
		text(rightCoverageAxisX+0.025, yCoverage, "0", adj=c(0,0.5), cex=0.9*pdfAttrs[["fontSize"]])
		text(rightCoverageAxisX+0.025, yCoverage+0.1, fusion_data$normFactor2, adj=c(0,0.5), cex=0.9*pdfAttrs[["fontSize"]])
		text(rightCoverageAxisX+0.05, yCoverage+0.08, "Coverage", srt=90, cex=0.9*pdfAttrs[["fontSize"]], adj=c(1,0.5))
		
		# plot coverage 1
		rect(min(fusion_data$exon1$left), yCoverage, max(fusion_data$exon1$right), 
		     yCoverage+0.1, col="#eeeeee", border=NA)
		
		drawCoverage(min(fusion_data$exon1$left), max(fusion_data$exon1$right),
		             yCoverage, fusion_data$coverage1, min(fusion_data$exon1$start), max(fusion_data$exon1$end),
		             pdfAttrs[["color1"]])
		

		# plot coverage 2
		rect(gene2Offset+min(fusion_data$exon2$left), yCoverage, 
		     gene2Offset+max(fusion_data$exon2$right), yCoverage+0.1, col="#eeeeee", border=NA)
		
		drawCoverage(gene2Offset+min(fusion_data$exon2$left), gene2Offset+max(fusion_data$exon2$right), 
		             yCoverage, fusion_data$coverage2, min(fusion_data$exon2$start), max(fusion_data$exon2$end),
		             pdfAttrs[["color2"]])
		
	}

	# plot gene 1
	lines(c(min(fusion_data$exon1$left), max(fusion_data$exon1$right)), c(yExons, yExons), 
	      col=pdfAttrs[["darkColor1"]])
	for (gene in unique(fusion_data$exon1$geneName)) {
		drawStrand(min(fusion_data$exon1[fusion_data$exon1$geneName == gene,"left"]), 
		           max(fusion_data$exon1[fusion_data$exon1$geneName == gene,"right"]), 
		           yExons, pdfAttrs[["darkColor1"]], 
		           head(fusion_data$exon1[fusion_data$exon1$geneName == gene,"strand"],1))
	}
	for (i in 1:nrow(fusion_data$exon1)) {
		drawExon(fusion_data$exon1[["left"]][i], fusion_data$exon1[["right"]][i],
		         yExons, pdfAttrs[["color1"]], fusion_data$exon1[["exonNumber.y"]][i],
		         fusion_data$exon1[["type"]][i])
	}

	# plot gene 2
	lines(c(gene2Offset, gene2Offset+max(fusion_data$exon2$right)), c(yExons, yExons), 
	      col=pdfAttrs[["darkColor2"]])
	for (gene in unique(fusion_data$exon2$geneName)) {
		drawStrand(gene2Offset+min(fusion_data$exon2[fusion_data$exon2$geneName == gene,"left"]),
		           gene2Offset+max(fusion_data$exon2[fusion_data$exon2$geneName == gene,"right"]), 
		           yExons, pdfAttrs[["darkColor2"]],
		           head(fusion_data$exon2[fusion_data$exon2$geneName == gene,"strand"],1))
	}
	for (i in 1:nrow(fusion_data$exon2)) {
		drawExon(gene2Offset+fusion_data$exon2[["left"]][i],
		         gene2Offset+fusion_data$exon2[["right"]][i],
		         yExons, pdfAttrs[["color2"]], 
		         fusion_data$exon2[["exonNumber.y"]][i], fusion_data$exon2[["type"]][i])
	}
	# plot gene1 of fusion
	if (table[["direction1"]] == "downstream") {
	  # plot strands
		lines(c(fusionOffset1, fusionOffset1+breakpoint1),
		      c(yFusion, yFusion), col=pdfAttrs[["darkColor2"]])
		for (gene in unique(fusion_data$exon1$geneName)) {
			exonsOfGene <- fusion_data$exon1[fusion_data$exon1$geneName == gene,]
			if (min(exonsOfGene$start) <= table[["breakpoint1"]]) {
				drawStrand(fusionOffset1+min(exonsOfGene$left), 
				           fusionOffset1+min(exonsOfGene$breakpoint, max(exonsOfGene$right)),
				           yFusion, col=pdfAttrs[["darkColor1"]], exonsOfGene$strand[1])}
		}
		# plot exons
		for (i in 1:nrow(fusion_data$exon1)) {
			if (fusion_data$exon1[["start"]][i] <= table[["breakpoint1"]]) {
			  drawExon(fusionOffset1 + fusion_data$exon1[["left"]],
				         fusionOffset1 + min(breakpoint1, fusion_data$exon1[["right"]]),
				         yFusion, pdfAttrs[["color1"]], fusion_data$exon1[["exonNumber.y"]],
				         fusion_data$exon1[["type"]]) }
		}
		# plot trajectories
		lines(c(0, 0, fusionOffset1), c(yTrajectoryExonTop, yTrajectoryExonBottom, yTrajectoryFusion), 
		      col="red", lty=2)
		lines(c(breakpoint1, breakpoint1, fusionOffset1+breakpoint1),
		      c(yTrajectoryBreakpointLabels, yTrajectoryExonBottom, yTrajectoryFusion), 
		      col="red", lty=2)
	} else if (table[["direction1"]] == "upstream") {
	  # plot strands
		lines(c(fusionOffset1, fusionOffset2), c(yFusion, yFusion), col=pdfAttrs[["darkColor1"]])
	  
		for (gene in unique(fusion_data$exon1$geneName)) {
			exonsOfGene <- fusion_data$exon1[fusion_data$exon1$geneName == gene,]
			if (max(exonsOfGene$end+1) >= table[["breakpoint1"]]) {
			  drawStrand(left = fusionOffset2-max(exonsOfGene$right)+breakpoint1, 
				           right = min(fusionOffset2, fusionOffset2-min(exonsOfGene$left)+breakpoint1),
				           y = yFusion, color = pdfAttrs[["darkColor1"]], strand = chartr("+-", "-+", exonsOfGene$strand[1]))
			}
		}
		# plot exons
		for (i in 1:nrow(fusion_data$exon1)) {
			if (fusion_data$exon1[["end"]][i]+1 >= table[["breakpoint1"]]) {
			  drawExon(left = fusionOffset1+max(fusion_data$exon1$right)-fusion_data$exon1[["right"]][i],
				         right = min(fusionOffset2, fusionOffset1+max(fusion_data$exon1$right)-fusion_data$exon1[["left"]][i]),
				         y = yFusion, color = pdfAttrs[["color1"]],title = fusion_data$exon1[["exonNumber.y"]][i],
				         type = fusion_data$exon1[["type"]][i])
			}
		}
		# plot trajectories
		lines(c(max(fusion_data$exon1$right), max(fusion_data$exon1$right), fusionOffset1),
		      c(yTrajectoryExonTop, yTrajectoryExonBottom, yTrajectoryFusion), col="red", lty=2)
		lines(c(breakpoint1, breakpoint1, fusionOffset1+max(fusion_data$exon1$right)-breakpoint1),
		      c(yTrajectoryBreakpointLabels, yTrajectoryExonBottom, yTrajectoryFusion), col="red", lty=2)
  }
	
	# plot gene2 of fusion
	if (table[["direction2"]] == "downstream") {
		# plot strands
		lines(c(fusionOffset2, fusionOffset2+breakpoint2), 
		      c(yFusion, yFusion), col=pdfAttrs[["darkColor2"]])
		for (gene in unique(fusion_data$exon2$geneName)) {
			exonsOfGene <- fusion_data$exon2[fusion_data$exon2$geneName == gene,]
			if (min(exonsOfGene$start) <= table[["breakpoint2"]])
				drawStrand(max(fusionOffset2, fusionOffset2+breakpoint2-max(exonsOfGene$right)),
				           fusionOffset2+breakpoint2-min(exonsOfGene$left), yFusion, 
				           col=pdfAttrs[["darkColor2"]], chartr("+-", "-+", exonsOfGene$strand[1]))
		}
		# plot exons
		for (i in 1:nrow(fusion_data$exon2)) {
			if (fusion_data$exon2[["start"]][i] <= table[["breakpoint2"]]) {
			  drawExon(max(fusionOffset2, fusionOffset2+breakpoint2-fusion_data$exon2[["right"]][i]),
			           fusionOffset2+breakpoint2-fusion_data$exon2[["left"]][i],
			           yFusion, pdfAttrs$color2, fusion_data$exon2[["exonNumber.y"]][i], 
			           fusion_data$exon2[["type"]][i])}
		}
		# plot trajectories
		lines(c(gene2Offset, gene2Offset, fusionOffset2+breakpoint2),
		      c(yTrajectoryExonTop, yTrajectoryExonBottom, yTrajectoryFusion), 
		      col="red", lty=2)
		lines(c(gene2Offset+breakpoint2, gene2Offset+breakpoint2, fusionOffset2),
		      c(yTrajectoryBreakpointLabels, yTrajectoryExonBottom, yTrajectoryFusion),
		      col="red", lty=2)
	} else if (table[["direction2"]] == "upstream") {
		# plot strands
		lines(c(fusionOffset2, fusionOffset2+max(fusion_data$exon2$right)-breakpoint2),
		      c(yFusion, yFusion), col=pdfAttrs[["darkColor2"]])
		for (gene in unique(fusion_data$exon2$geneName)) {
			exonsOfGene <- fusion_data$exon2[fusion_data$exon2$geneName == gene,]
			if (max(exonsOfGene$end+1) >= table[["breakpoint2"]])
				drawStrand(max(fusionOffset2, fusionOffset2+min(exonsOfGene$left)-breakpoint2),
				           fusionOffset2+max(exonsOfGene$right)-breakpoint2, 
				           yFusion, col=pdfAttrs[["darkColor2"]], exonsOfGene$strand[1])
		}
		# plot exons
		for (i in 1:nrow(fusion_data$exon2)) {
			if (fusion_data$exon2[["end"]][i]+1 >= table[["breakpoint2"]]) {
				drawExon(max(fusionOffset2, fusionOffset2+fusion_data$exon2[["left"]][i]-breakpoint2),
				         fusionOffset2+fusion_data$exon2[["right"]]-breakpoint2, 
				         yFusion, pdfAttrs[["color2"]], fusion_data$exon2[["exonNumber.y"]][i],
				         fusion_data$exon2[["type"]][i])}
		}
		# plot trajectories
		lines(c(gene2Offset+max(fusion_data$exon2$right),
		        gene2Offset+max(fusion_data$exon2$right),
		        fusionOffset2+max(fusion_data$exon2$right)-breakpoint2),
		      c(yTrajectoryExonTop, yTrajectoryExonBottom, yTrajectoryFusion),
		      col="red", lty=2)
		lines(c(gene2Offset+breakpoint2, gene2Offset+breakpoint2, fusionOffset2),
		      c(yTrajectoryBreakpointLabels, yTrajectoryExonBottom, yTrajectoryFusion),
		      col="red", lty=2)
	}
	
	if (table[["fusion_transcript"]] != ".") {
		# print fusion transcript colored by gene of origin
		fusion_transcript1 <- gsub("\\|.*", "", table[["fusion_transcript"]], perl=T)
		fusion_transcript1 <- substr(fusion_transcript1, max(1, nchar(fusion_transcript1)-30), nchar(fusion_transcript1))
		fusion_transcript2 <- gsub(".*\\|", "", table[["fusion_transcript"]], perl=T)
		fusion_transcript2 <- substr(fusion_transcript2, 1, min(nchar(fusion_transcript2), 30))
		# check for non-template bases
		non_template_bases <- gsub(".*\\|([^|]*)\\|.*", "\\1", table[["fusion_transcript"]], perl=T)
		if (non_template_bases == table[["fusion_transcript"]]) # no non-template bases found
			non_template_bases <- ""
		# divide non-template bases half-and-half for centered alignment
		non_template_bases1 <- substr(non_template_bases, 1, floor(nchar(non_template_bases)/2))
		non_template_bases2 <- substr(non_template_bases, ceiling(nchar(non_template_bases)/2+0.5), 
		                              nchar(non_template_bases))
		# transcript 1
		text(fusionOffset2, yTranscript, bquote(.(fusion_transcript1) * phantom(.(non_template_bases1))), 
		     col=pdfAttrs$darkColor1, adj=c(1,0.5), cex=pdfAttrs$fontSize)
		# transcript 2
		text(fusionOffset2, yTranscript, bquote(phantom(.(non_template_bases2)) * .(fusion_transcript2)), 
		     col=pdfAttrs$darkColor2, adj=c(0,0.5), cex=pdfAttrs$fontSize)
		# non-template bases
		text(fusionOffset2, yTranscript, non_template_bases1, adj=c(1,0.5), cex=pdfAttrs$fontSize)
		text(fusionOffset2, yTranscript, non_template_bases2, adj=c(0,0.5), cex=pdfAttrs$fontSize)
	}

	# draw scale
	realScale <- max(fusion_data$exon1$end - fusion_data$exon1$start,
	                 fusion_data$exon2$end - fusion_data$exon2$start)
	mapScale <- max(fusion_data$exon1$right - fusion_data$exon1$left, 
	                fusion_data$exon2$right - fusion_data$exon2$left)
	# choose scale which is closest to desired scale length
	desiredScaleSize <- 0.2
	realScale <- desiredScaleSize / mapScale * realScale
	mapScale <- desiredScaleSize
	realScaleOptimalFit <- signif(realScale, 1) # round to most significant digit
	mapScaleOptimalFit <- realScaleOptimalFit / realScale * mapScale
	# draw scale line
	lines(c(1-mapScaleOptimalFit, 1), c(yScale, yScale)) # scale line
	lines(c(1-mapScaleOptimalFit, 1-mapScaleOptimalFit), c(yScale-0.007, yScale+0.007)) # left whisker
	lines(c(1, 1), c(yScale-0.007, yScale+0.007)) # right whisker
	# draw units above scale line
	realScaleThousands <- max(0, min(3, floor(log10(realScaleOptimalFit)/3)))
	scaleUnits <- c("bp", "kbp", "Mbp", "Gbp")
	scaleLabel <- paste(realScaleOptimalFit/max(1,1000^realScaleThousands), scaleUnits[realScaleThousands+1])
	text(1-mapScaleOptimalFit/2, yScale+0.005, scaleLabel, adj=c(0.5,0), cex=pdfAttrs$fontSize*0.9)
	
	# draw circos plot
	drawCircos(fusions = table, cytobands = cytobands, 
		           minConfidenceForCircosPlot = "medium", circosColors = circosColors)
	par(mar=c(0,0,0,0), xpd=F)

	# draw protein domains
	plot(0, 0, type="l", xlim=c(-0.1, 1.1), ylim=c(0, 1), bty="n", xaxt="n", yaxt="n", xlab="", ylab="")
	par(xpd=NA)
	if (!is.null(ProteinDomains) && "GenomicRanges" %in% names(sessionInfo()$otherPkgs)) {
		drawProteinDomains(table[1,], fusion_data$exon1, fusion_data$exon2,
		                   ProteinDomains, pdfAttrs$color1, pdfAttrs$color2,
		                   mergeDomainsOverlappingBy = 0.9, optimizeDomainColors = T)
	}
	par(xpd=F)

	# print statistics about supporting alignments
	plot(0, 0, type="l", xlim=c(0, 1), ylim=c(0, 1), bty="n", xaxt="n", yaxt="n", xlab="", ylab="")
	text(0, 0.575, "SUPPORTING READ COUNT", font=2, adj=c(0,0), cex=pdfAttrs$fontSize)
	text(0, 0.525, paste0("Split reads = ", table[["split_reads"]], "\n", "Discordant mates = ", 
	                      table[["discordant_mates"]]), adj=c(0,1), cex=pdfAttrs$fontSize)
	
	
	devNull <- dev.off()
	message("Done")
}
plot_fusion(fusion_data = fusion_POLCE2_ATR, table = VI_3429_593_2DN_fusions[1,],
            filename = "fusion_PCOLCE2_ATR.pdf", title = "PCOLCE2--ATR fusion")

plot_fusion(fusion_data = fusion_SCAP_SLC6A20, table = VI_3429_593_2DN_fusions[2,],
            filename = "fusion_SCAP_SLC6A20.pdf", title = "SCAP--SLC6A20 fusion")
