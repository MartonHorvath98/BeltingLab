# ------------------- Data processing functions --------------------------------
# Recursively merges a list of data frames into a single data frame
merge.rec <- function(.list, ...){
  if(length(.list)==1) return(.list[[1]])
  Recall(c(list(merge(.list[[1]], .list[[2]], ...)), .list[-(1:2)]), ...)
}
# Read tables from the input directory; by default only tables where the number of supporting 
# reads is counted by the junction reads only (JC)
read_tables <- function(path, condition, pattern = ".MATS.JC.txt", sep = "\t") {
  files <- list.files(file.path(path, condition),
                      pattern = pattern, full.names = T)
  tables <- lapply(files, function(file) read.table(file, header = T, sep = sep))
  names(tables) <- gsub(pattern, "", basename(files))
  return(tables)
}

filter_tables <- function(df, conditions = list(cond1 = "C1", cond2 = "C2"),
                          thresholds = list("FDR" = 0.05, "count" = 20)){
  
  columns = c("ID","GeneID", "IJC_SAMPLE_1", "SJC_SAMPLE_1",
              "IJC_SAMPLE_2", "SJC_SAMPLE_2", "PValue", "FDR", "IncLevelDifference")
  
  include = c(paste("IJC", conditions$cond1, c(1:3), sep = "_"),
              paste("IJC", conditions$cond2, c(1:3), sep = "_"))
  skip = c(paste("SJC", conditions$cond1, c(1:3), sep = "_"),
           paste("SJC", conditions$cond2, c(1:3), sep = "_"))
  
  df = df %>% 
    dplyr::select(., columns)
  
  df = df %>%
    tidyr::separate(., col = IJC_SAMPLE_1, into = include[1:3], sep = ",") %>%
    tidyr::separate(., col = SJC_SAMPLE_1, into = skip[1:3], sep = ",") %>%
    tidyr::separate(., col = IJC_SAMPLE_2, into = include[4:6], sep = ",") %>%
    tidyr::separate(., col = SJC_SAMPLE_2, into = skip[4:6], sep = ",") %>% 
    dplyr::mutate(across(!c("ID","GeneID"), as.numeric)) %>% 
    dplyr::rowwise() %>%
    dplyr::mutate(count1 = sum(c_across(contains(conditions$cond1)), na.rm = TRUE),
                  count2 = sum(c_across(contains(conditions$cond2)), na.rm = TRUE)) %>%
    dplyr::ungroup()
  
  df = df %>% 
    dplyr::mutate(GeneID = gsub("\\..*","",GeneID)) %>% 
    dplyr::mutate(
      Symbol = AnnotationDbi::mapIds(org.Hs.eg.db, keys = GeneID, 
                                     keytype = "ENSEMBL", column = "SYMBOL"),
      EntrezID = AnnotationDbi::mapIds(org.Hs.eg.db, keys = GeneID, 
                                                 keytype = "ENSEMBL", column = "ENTREZID")) %>% 
    dplyr::relocate(EntrezID, Symbol, .after = GeneID) %>%
    dplyr::filter(FDR < thresholds$FDR & 
                    count1 >= thresholds$count & 
                    count2 >= thresholds$count)
    
  return(df)  
}

# A function to get the variants extracted and filtered through somatic variants
# and SURFME genes
filter_vcf <- function(vcf_obj, info_col, somatic_table, surfme){
  # create a list of data frames from the vcf objects
  df <- vcfR::vcfR2tidy(vcf_obj)
  df <- df$fix
  # Split the Info column
  split_info <- strsplit(df[[info_col]], "\\|")
  # Convert the list to a matrix
  info_matrix <- do.call(rbind, lapply(split_info, "[", 2:4))
  info_matrix <- data.frame(info_matrix)
  colnames(info_matrix) <- c("mutType","mutEffect", "geneSymbol")
  
  df <- cbind(df, info_matrix)
  
  som.df <- merge(df, somatic_table, by = c("CHROM", "POS", "REF","ALT"))
  
  surf.df <- som.df[which(som.df$geneSymbol %in% surfme),]
  
  return(list(df = df, somatic = som.df, surface = surf.df))
}

# Plot pie chart
plot_pie <- function(data, col, title){
  df <- data.frame(mut = data[[col]])
    
  df <- df %>% 
    dplyr::rowwise() %>%
    dplyr::mutate(mut = strsplit(mut, "&")) %>%
    unnest(mut) %>% 
    dplyr::mutate(mut = as.factor(mut))
  
  data <- df %>%
    dplyr::group_by(mut) %>% 
    dplyr::summarise(count = n())
  
  data <- data %>% 
    arrange(desc(mut)) %>%
    mutate(prop = count / sum(data$count) *100)
  
  data$label <- paste0(data$mut, ": ", round(data$prop, 2), "%")
  
  
  # Create the pie chart
  return(
    ggplot(data, aes(x = "", y = prop, fill = label)) +
      geom_bar(stat = "identity", width = 1, color = "white") +
      coord_polar(theta = "y") +
      labs(title = title, x = NULL, y = NULL, fill = "") +
      theme_void() +
      scale_fill_viridis_d(
        guide = guide_legend(override.aes=list(shape = 16))))
}


# define colors
changeColorBrightness <- function(color, delta) {
  rgb(
    min(255,max(0,col2rgb(color)["red",]+delta)),
    min(255,max(0,col2rgb(color)["green",]+delta)),
    min(255,max(0,col2rgb(color)["blue",]+delta)),
    maxColorValue=255
  )
}
getDarkColor <- function(color) { changeColorBrightness(color, -100) }
getBrightColor <- function(color) { changeColorBrightness(color, +190) }

# convenience functions to add/remove "chr" prefix
addChr <- function(contig) {
  ifelse(contig == "MT", "chrM", paste0("chr", contig))
}
removeChr <- function(contig) {
  sub("^chr", "", sub("^chrM", "MT", contig, perl=T), perl=T)
}

# convenience function to check if a value is between two others
between <- function(value, start, end) {
  value >= start & value <= end
}

prep_fusions <- function(fusions){
  fusions$gene1 <- sub("\\^.*", "", fusions$LeftGene, perl=T)
  fusions$gene2 <- sub("\\^.*", "", fusions$RightGene, perl=T)
  fusions$strand1 <- sub(".*:(.)$", "\\1/\\1", fusions$LeftBreakpoint, perl=T)
  fusions$strand2 <- sub(".*:(.)$", "\\1/\\1", fusions$RightBreakpoint, perl=T)
  fusions$display_contig1 <- sub(":[^:]*:[^:]*$", "", fusions$LeftBreakpoint, perl=T)
  fusions$display_contig2 <- sub(":[^:]*:[^:]*$", "", fusions$RightBreakpoint, perl=T)
  fusions$contig1 <- removeChr(fusions$display_contig1)
  fusions$contig2 <- removeChr(fusions$display_contig2)
  fusions$breakpoint1 <- as.numeric(sub(".*:([^:]*):[^:]*$", "\\1", fusions$LeftBreakpoint, perl=T))
  fusions$breakpoint2 <- as.numeric(sub(".*:([^:]*):[^:]*$", "\\1", fusions$RightBreakpoint, perl=T))
  fusions$direction1 <- ifelse(grepl(":\\+$", fusions$LeftBreakpoint, perl=T), "downstream", "upstream")
  fusions$direction2 <- ifelse(grepl(":\\+$", fusions$RightBreakpoint, perl=T), "upstream", "downstream")
  fusions$gene_id1 <- sub(".*\\^", "", fusions$LeftGene, perl=T)
  fusions$gene_id2 <- sub(".*\\^", "", fusions$RightGene, perl=T)
  fusions$transcript_id1 <- ifelse(rep(!("CDS_LEFT_ID" %in% colnames(fusions)), nrow(fusions)), ".", fusions$CDS_LEFT_ID)
  fusions$transcript_id2 <- ifelse(rep(!("CDS_RIGHT_ID" %in% colnames(fusions)), nrow(fusions)), ".", fusions$CDS_RIGHT_ID)
  fusions$fusion_transcript <- ifelse(rep(!("FUSION_CDS" %in% colnames(fusions)), nrow(fusions)), ".", toupper(sub("([a-z]*)", "\\1|", fusions$FUSION_CDS, perl=T)))
  fusions$reading_frame <- ifelse(rep(!("PROT_FUSION_TYPE" %in% colnames(fusions)), nrow(fusions)), ".", ifelse(fusions$PROT_FUSION_TYPE == "INFRAME", "in-frame", ifelse(fusions$PROT_FUSION_TYPE == "FRAMESHIFT", "out-of-frame", ".")))
  fusions$split_reads <- fusions$JunctionReadCount
  fusions$discordant_mates <- fusions$SpanningFragCount
  fusions$site1 <- rep("exon", nrow(fusions))
  fusions$site2 <- rep("exon", nrow(fusions))
  fusions$confidence <- rep("high", nrow(fusions))
  fusions$type <- ifelse(fusions$contig1 != fusions$contig2, "translocation", ifelse(fusions$direction1 == fusions$direction2, "inversion", ifelse((fusions$direction1 == "downstream") == (fusions$breakpoint1 < fusions$breakpoint2), "deletion", "duplication")))
  return(fusions)
}

parseGtfAttribute <- function(attribute, gtf) {
  parsed <- sub(paste0(".*", attribute, "[ =]([^;]+).*"), "\\1", gtf$attributes, perl=T)
  failedToParse <- parsed == gtf$attributes
  if (any(failedToParse)) {
    warning(paste0("Failed to parse '", attribute, "' attribute of ", sum(failedToParse), " record(s)."))
    parsed <- ifelse(failedToParse, "", parsed)
  }
  return(parsed)
}

drawVerticalGradient <- function(left, right, y, color, selection=NULL) {
  # check if gradient should only be drawn in part of the region
  if (!is.null(selection)) {
    y <- y[selection]
    left <- left[selection]
    right <- right[selection]
  }
  # draw gradient
  for (i in 1:length(y)) {
    polygon(
      c(left[1:i], right[1:i]),
      c(y[1:i], y[i:1]),
      border=NA,
      col=rgb(col2rgb(color)["red",], col2rgb(color)["green",], col2rgb(color)["blue",], col2rgb(color, alpha=T)["alpha",]*(1/length(y)), max=255)
    )
  }
}

drawCurlyBrace <- function(left, right, top, bottom, tip) {
  smoothness <- 20
  x <- cumsum(exp(-seq(-2.5, 2.5, len=smoothness)^2))
  x <- x/max(x)
  y <- seq(top, bottom, len=smoothness)
  lines(left+(tip-left)+x*(left-tip), y)
  lines(tip+x*(right-tip), y)
}

drawIdeogram <- function(adjust, left, right, y, cytobands, contig, breakpoint) {
  # define design of ideogram
  bandColors <- setNames(rgb(100:0, 100:0, 100:0, maxColorValue=100), paste0("gpos", 0:100))
  bandColors <- c(bandColors, gneg="#ffffff", acen="#ec4f4f", gvar = "#cc66ff", stalk="#0000ff")
  cytobands$color <- bandColors[as.character(cytobands$giemsa)]
  arcSteps <- 30 # defines roundness of arc
  curlyBraceHeight <- 0.03
  ideogramHeight <- 0.04
  ideogramWidth <- 0.4
  # extract bands of given contig
  bands <- cytobands[cytobands$contig==contig,]
  if (nrow(bands) == 0) {
    warning(paste("Ideogram of contig", contig, "cannot be drawn, because no Giemsa staining information is available."))
    return(NULL)
  }
  # scale width of ideogram to fit inside given region
  bands$left <- bands$start / max(cytobands$end) * ideogramWidth
  bands$right <- bands$end / max(cytobands$end) * ideogramWidth
  # left/right-align cytobands
  offset <- ifelse(adjust=="left", left, right - max(bands$right))
  bands$left <- bands$left + offset
  bands$right <- bands$right + offset
  # draw curly braces
  tip <- min(bands$left) + (max(bands$right)-min(bands$left)) / (max(bands$end)-min(bands$start)) * breakpoint
  drawCurlyBrace(left, right, y-0.05+curlyBraceHeight, y-0.05, tip)
  # draw title of chromosome
  text((max(bands$right)+min(bands$left))/2, y+0.07, paste("chromosome", contig),
       font=2, cex=pdfAttrs$fontSize, adj=c(0.5,0))
  # draw name of band
  bandName <- bands[which(between(breakpoint, bands$start, bands$end)), "name"]
  text(tip, y+0.03, bandName, cex=pdfAttrs$fontSize, adj=c(0.5,0))
  # draw start of chromosome
  leftArcX <- bands[1,"left"] + (1+cos(seq(pi/2,1.5*pi,len=arcSteps))) * (bands[1,"right"]-bands[1,"left"])
  leftArcY <- y + sin(seq(pi/2,1.5*pi,len=arcSteps)) * (ideogramHeight/2)
  polygon(leftArcX, leftArcY, col=bands[1,"color"])
  # draw bands
  centromereStart <- NULL
  centromereEnd <- NULL
  for (band in 2:(nrow(bands)-1)) {
    if (bands[band,"giemsa"] != "acen") {
      rect(bands[band,"left"], y-ideogramHeight/2, bands[band,"right"], y+ideogramHeight/2, col=bands[band,"color"])
    } else { # draw centromere
      if (is.null(centromereStart)) {
        polygon(c(bands[band,"left"], bands[band,"right"], bands[band,"left"]), c(y-ideogramHeight/2, y, y+ideogramHeight/2), col=bands[band,"color"])
        centromereStart <- bands[band,"left"]
      } else {
        polygon(c(bands[band,"right"], bands[band,"left"], bands[band,"right"]), c(y-ideogramHeight/2, y, y+ideogramHeight/2), col=bands[band,"color"])
        centromereEnd <- bands[band,"right"]
      }
    }
  }
  # draw end of chromosome
  band <- nrow(bands)
  rightArcX <- bands[band,"right"] - (1+cos(seq(1.5*pi,pi/2,len=arcSteps))) * (bands[band,"right"]-bands[band,"left"])
  rightArcY <- y + sin(seq(pi/2,1.5*pi,len=arcSteps)) * ideogramHeight/2
  polygon(rightArcX, rightArcY, col=bands[band,"color"])
  # if there is no centromere, make an artificial one with length zero
  if (is.null(centromereStart) || is.null(centromereEnd)) {
    centromereStart <- bands[1,"right"]
    centromereEnd <- bands[1,"right"]
  }
  # draw gradients for 3D effect
  drawVerticalGradient(leftArcX, rep(centromereStart, arcSteps), leftArcY, rgb(0,0,0,0.8), 1:round(arcSteps*0.4)) # black from top on p-arm
  drawVerticalGradient(leftArcX, rep(centromereStart, arcSteps), leftArcY, rgb(1,1,1,0.7), round(arcSteps*0.4):round(arcSteps*0.1)) # white to top on p-arm
  drawVerticalGradient(leftArcX, rep(centromereStart, arcSteps), leftArcY, rgb(1,1,1,0.7), round(arcSteps*0.4):round(arcSteps*0.6)) # white to bottom on p-arm
  drawVerticalGradient(leftArcX, rep(centromereStart, arcSteps), leftArcY, rgb(0,0,0,0.9), arcSteps:round(arcSteps*0.5)) # black from bottom on p-arm
  drawVerticalGradient(rightArcX, rep(centromereEnd, arcSteps), rightArcY, rgb(0,0,0,0.8), 1:round(arcSteps*0.4)) # black from top on q-arm
  drawVerticalGradient(rightArcX, rep(centromereEnd, arcSteps), rightArcY, rgb(1,1,1,0.7), round(arcSteps*0.4):round(arcSteps*0.1)) # white to top on q-arm
  drawVerticalGradient(rightArcX, rep(centromereEnd, arcSteps), rightArcY, rgb(1,1,1,0.7), round(arcSteps*0.4):round(arcSteps*0.6)) # white to bottom on q-arm
  drawVerticalGradient(rightArcX, rep(centromereEnd, arcSteps), rightArcY, rgb(0,0,0,0.9), arcSteps:round(arcSteps*0.5)) # black from bottom on q-arm
}

drawCoverage <- function(left, right, y, coverage, start, end, color) {
  maxResolution <- 5000 # max number of data points to draw coverage
  # draw coverage as bars
  if (!is.null(coverage)) {
    coverageData <- as.numeric(coverage[IRanges(sapply(start, max, min(start(coverage))), sapply(end, min, max(end(coverage))))])
    # downsample to maxResolution, if there are too many data points
    coverageData <- aggregate(coverageData, by=list(round(1:length(coverageData) * (right-left) * maxResolution/length(coverageData))), mean)$x
    polygon(c(left, seq(left, right, length.out=length(coverageData)), right), c(y, y+coverageData*0.1, y), col=color, border=NA)
  }
}

drawStrand <- function(left, right, y, color, strand) {
  if (strand %in% c("+", "-")) {
    # draw strand
    lines(c(left+0.001, right-0.001), c(y, y), col=color, lwd=2)
    lines(c(left+0.001, right-0.001), c(y, y), col=rgb(1,1,1,0.1), lwd=1)
    # indicate orientation
    if (right - left > 0.01)
      for (i in seq(left+0.005, right-0.005, by=sign(right-left-2*0.005)*0.01)) {
        arrows(i, y, i+0.001*ifelse(strand=="+", 1, -1), y, col=color, length=0.05, lwd=2, angle=60)
        arrows(i, y, i+0.001*ifelse(strand=="+", 1, -1), y, col=rgb(1,1,1,0.1), length=0.05, lwd=1, angle=60)
      }
  }
}

drawExon <- function(left, right, y, color, title, type) {
  gradientSteps <- 10 # defines smoothness of gradient
  exonHeight <- 0.03
  if (type == "CDS") {
    # draw coding regions as thicker bars
    rect(left, y+exonHeight, right, y+exonHeight/2-0.001, col=color, border=NA)
    rect(left, y-exonHeight, right, y-exonHeight/2+0.001, col=color, border=NA)
    # draw border
    lines(c(left, left, right, right), c(y+exonHeight/2, y+exonHeight, y+exonHeight, y+exonHeight/2), col=getDarkColor(color), lend=2)
    lines(c(left, left, right, right), c(y-exonHeight/2, y-exonHeight, y-exonHeight, y-exonHeight/2), col=getDarkColor(color), lend=2)
    # draw gradients for 3D effect
    drawVerticalGradient(rep(left, gradientSteps), rep(right, gradientSteps), seq(y+0.03, y+0.015, len=gradientSteps), rgb(0,0,0,0.2))
    drawVerticalGradient(rep(left, gradientSteps), rep(right, gradientSteps), seq(y-0.03, y-0.015, len=gradientSteps), rgb(0,0,0,0.3))
  } else if (type == "exon") {
    rect(left, y+exonHeight/2, right, y-exonHeight/2, col=color, border=getDarkColor(color))
    # draw gradients for 3D effect
    drawVerticalGradient(rep(left, gradientSteps), rep(right, gradientSteps), seq(y, y+exonHeight/2, len=gradientSteps), rgb(1,1,1,0.6))
    drawVerticalGradient(rep(left, gradientSteps), rep(right, gradientSteps), seq(y, y-exonHeight/2, len=gradientSteps), rgb(1,1,1,0.6))
    # add exon label
    text((left+right)/2, y, title, cex=0.9*pdfAttrs[["fontSize"]])
  }
}

drawCircos <- function(fusions, cytobands, minConfidenceForCircosPlot, circosColors) {
  # check if Giemsa staining information is available
  for (contig in unlist(fusions[c("contig1", "contig2")])) {
    if (!any(cytobands$contig==contig)) {
      warning(paste0("Circos plot cannot be drawn, because no Giemsa staining information is available for contig ", contig, "."))
      # draw empty plots as placeholder
      plot(0, 0, type="l", xlim=c(0, 1), ylim=c(0, 1), bty="n", xaxt="n", yaxt="n", xlab="", ylab="")
      plot(0, 0, type="l", xlim=c(0, 1), ylim=c(0, 1), bty="n", xaxt="n", yaxt="n", xlab="", ylab="")
      return(NULL)
    }
  }
  # initialize with empty circos plot
  circos.clear()
  circos.initializeWithIdeogram(cytoband=cytobands, plotType=NULL)
  # use gene names as labels or <contig>:<position> for intergenic breakpoints
  geneLabels <- data.frame(
    contig=c(fusions[["contig1"]], fusions[["contig2"]]),
    start=c(fusions[["breakpoint1"]], fusions[["breakpoint2"]])
  )
  geneLabels$end <- geneLabels$start + 1
  geneLabels$gene <- c(fusions[["gene1"]], fusions[["gene2"]])
  geneLabels$gene <- ifelse(c(fusions[["site1"]], fusions[["site2"]]) == "intergenic",
                            paste0(c(fusions[["display_contig1"]], fusions[["display_contig2"]]), ":", geneLabels$start),
                            geneLabels$gene)
  # draw gene labels
  circos.genomicLabels(geneLabels, labels.column=4, side="outside", cex=pdfAttrs$fontSize, labels_height=0.27)
  # draw chromosome labels in connector plot
  for (contig in unique(cytobands$contig)) {
    set.current.cell(track.index=2, sector.index=contig) # draw in gene label connector track (track.index=2)
    circos.text(CELL_META$xcenter, CELL_META$ycenter, contig, cex=0.85)
  }
  # draw ideograms
  circos.genomicIdeogram(cytoband=cytobands)
  # draw arcs
  confidenceRank <- c(low=0, medium=1, high=2)
  for (i in c(setdiff(1:nrow(fusions), fusion), fusion)) { # draw fusion of interest last, such that its arc is on top
    f <- fusions[i,]
    if (any(cytobands$contig==f$contig1) && any(cytobands$contig==f$contig2)) # ignore viral contigs, because we have no cytoband information for them
      if (minConfidenceForCircosPlot != "none" && confidenceRank[f$confidence] >= confidenceRank[minConfidenceForCircosPlot] || i==fusion)
        circos.link(
          f$contig1, f$breakpoint1,
          f$contig2, f$breakpoint2,
          lwd=2, col=ifelse(i==fusion, circosColors[f$type], getBrightColor(circosColors[f$type]))
        )
  }
  # draw legend
  plot(0, 0, type="l", xlim=c(0, 1), ylim=c(0, 1), bty="n", xaxt="n", yaxt="n", ylab="", xlab="")
  legend(x="top", legend=names(circosColors), col=sapply(circosColors, getBrightColor), lwd=3, ncol=2, box.lty=0)
}

drawProteinDomains <- function(fusion, exons1, exons2, ProteinDomains, color1, color2, mergeDomainsOverlappingBy, optimizeDomainColors) {
  
  exonHeight <- 0.2
  exonsY <- 0.5
  geneNamesY <- exonsY - exonHeight/2 - 0.05
  
  # find coding exons
  codingExons1 <- exons1[exons1$type == "CDS" & fusion$site1 != "intergenic",]
  codingExons2 <- exons2[exons2$type == "CDS" & fusion$site2 != "intergenic",]
  
  # cut off coding regions beyond breakpoint
  if (fusion$direction1 == "upstream") {
    codingExons1 <- codingExons1[codingExons1$end >= fusion$breakpoint1,]
    codingExons1$start <- ifelse(codingExons1$start < fusion$breakpoint1, fusion$breakpoint1, codingExons1$start)
  } else {
    codingExons1 <- codingExons1[codingExons1$start <= fusion$breakpoint1,]
    codingExons1$end <- ifelse(codingExons1$end > fusion$breakpoint1, fusion$breakpoint1, codingExons1$end)
  }
  if (fusion$direction2 == "upstream") {
    codingExons2 <- codingExons2[codingExons2$end >= fusion$breakpoint2,]
    codingExons2$start <- ifelse(codingExons2$start < fusion$breakpoint2, fusion$breakpoint2, codingExons2$start)
  } else {
    codingExons2 <- codingExons2[codingExons2$start <= fusion$breakpoint2,]
    codingExons2$end <- ifelse(codingExons2$end > fusion$breakpoint2, fusion$breakpoint2, codingExons2$end)
  }
  
  # find overlapping domains
  exonsGRanges1 <- GRanges(codingExons1$contig, IRanges(codingExons1$start, codingExons1$end), strand=codingExons1$strand)
  exonsGRanges2 <- GRanges(codingExons2$contig, IRanges(codingExons2$start, codingExons2$end), strand=codingExons2$strand)
  domainsGRanges <- GRanges(ProteinDomains$contig, IRanges(ProteinDomains$start, ProteinDomains$end), strand=ProteinDomains$strand)
  domainsGRanges$proteinDomainName <- ProteinDomains$proteinDomainName
  domainsGRanges$proteinDomainID <- ProteinDomains$proteinDomainID
  domainsGRanges$color <- ProteinDomains$color
  domainsGRanges <- domainsGRanges[suppressWarnings(unique(queryHits(findOverlaps(domainsGRanges, union(exonsGRanges1, exonsGRanges2)))))]
  
  # group overlapping domains by domain ID
  domainsGRangesList <- GRangesList(lapply(unique(domainsGRanges$proteinDomainID), function(x) { domainsGRanges[domainsGRanges$proteinDomainID == x] }))
  
  # trim protein domains to exon boundaries
  trimDomains <- function(domainsGRangesList, exonsGRanges) {
    do.call(
      "rbind",
      lapply(
        domainsGRangesList,
        function(x) {
          intersected <- as.data.frame(reduce(suppressWarnings(intersect(x, exonsGRanges))))
          if (nrow(intersected) > 0) {
            intersected$proteinDomainName <- head(x$proteinDomainName, 1)
            intersected$proteinDomainID <- head(x$proteinDomainID, 1)
            intersected$color <- head(x$color, 1)
          } else {
            intersected$proteinDomainName <- character()
            intersected$proteinDomainID <- character()
            intersected$color <- character()
          }
          return(intersected)
        }
      )
    )
  }
  retainedDomains1 <- trimDomains(domainsGRangesList, exonsGRanges1)
  retainedDomains2 <- trimDomains(domainsGRangesList, exonsGRanges2)
  
  # calculate length of coding exons
  codingExons1$length <- codingExons1$end - codingExons1$start + 1
  codingExons2$length <- codingExons2$end - codingExons2$start + 1
  
  # abort, if there are no coding regions
  if (sum(exons1$type == "CDS") + sum(exons2$type == "CDS") == 0) {
    text(0.5, 0.5, "Genes are not protein-coding.")
    return(NULL)
  }
  codingLength1 <- sum(codingExons1$length)
  codingLength2 <- sum(codingExons2$length)
  if (codingLength1 + codingLength2 == 0) {
    text(0.5, 0.5, "No coding regions retained in fusion transcript.")
    return(NULL)
  }
  if ((codingLength1 == 0 || grepl("\\.$", fusion$strand1)) && (codingLength2 == 0 || grepl("\\.$", fusion$strand2))) {
    text(0.5, 0.5, "Failed to determine retained protein domains due to lack of strand information.")
    return(NULL)
  }
  antisenseTranscription1 <- sub("/.*", "", fusion$strand1) != sub(".*/", "", fusion$strand1)
  antisenseTranscription2 <- sub("/.*", "", fusion$strand2) != sub(".*/", "", fusion$strand2)
  if ((codingLength1 == 0 || antisenseTranscription1) && (codingLength2 == 0 || antisenseTranscription2)) {
    text(0.5, 0.5, "No coding regions due to antisense transcription.")
    return(NULL)
  }
  
  # remove introns from protein domains
  removeIntronsFromProteinDomains <- function(codingExons, retainedDomains) {
    if (nrow(codingExons) == 0) return(NULL)
    cumulativeIntronLength <- 0
    previousExonEnd <- 0
    for (exon in 1:nrow(codingExons)) {
      if (codingExons[exon,"start"] > previousExonEnd)
        cumulativeIntronLength <- cumulativeIntronLength + codingExons[exon,"start"] - previousExonEnd
      domainsInExon <- which(between(retainedDomains$start, codingExons[exon,"start"], codingExons[exon,"end"]))
      retainedDomains[domainsInExon,"start"] <- retainedDomains[domainsInExon,"start"] - cumulativeIntronLength
      domainsInExon <- which(between(retainedDomains$end, codingExons[exon,"start"], codingExons[exon,"end"]))
      retainedDomains[domainsInExon,"end"] <- retainedDomains[domainsInExon,"end"] - cumulativeIntronLength
      previousExonEnd <- codingExons[exon,"end"]
    }
    # merge adjacent domains
    retainedDomains <- do.call(
      "rbind",
      lapply(
        unique(retainedDomains$proteinDomainID),
        function(x) {
          domain <- retainedDomains[retainedDomains$proteinDomainID == x,]
          merged <- reduce(GRanges(domain$seqnames, IRanges(domain$start, domain$end), strand=domain$strand))
          merged$proteinDomainName <- head(domain$proteinDomainName, 1)
          merged$proteinDomainID <- head(domain$proteinDomainID, 1)
          merged$color <- head(domain$color, 1)
          return(as.data.frame(merged))
        }
      )
    )
    return(retainedDomains)
  }
  retainedDomains1 <- removeIntronsFromProteinDomains(codingExons1, retainedDomains1)
  retainedDomains2 <- removeIntronsFromProteinDomains(codingExons2, retainedDomains2)
  
  # abort, if no domains are retained
  if (is.null(retainedDomains1) && is.null(retainedDomains2)) {
    text(0.5, 0.5, "No protein domains retained in fusion.")
    return(NULL)
  }
  
  # merge domains with similar coordinates
  mergeSimilarDomains <- function(domains, mergeDomainsOverlappingBy) {
    if (is.null(domains)) return(domains)
    merged <- domains[F,] # create empty data frame
    domains <- domains[order(domains$end - domains$start, decreasing=T),] # start with bigger domains => bigger domains are retained
    for (domain in rownames(domains)) {
      if (!any((abs(merged$start - domains[domain,"start"]) + abs(merged$end - domains[domain,"end"])) / (domains[domain,"end"] - domains[domain,"start"] + 1) <= 1-mergeDomainsOverlappingBy))
        merged <- rbind(merged, domains[domain,])
    }
    return(merged)
  }
  retainedDomains1 <- mergeSimilarDomains(retainedDomains1, mergeDomainsOverlappingBy)
  retainedDomains2 <- mergeSimilarDomains(retainedDomains2, mergeDomainsOverlappingBy)
  
  # if desired, reassign colors to protein domains to maximize contrast
  if (optimizeDomainColors) {
    uniqueDomains <- unique(c(retainedDomains1$proteinDomainID, retainedDomains2$proteinDomainID))
    # make rainbow of pretty pastell colors
    colors <- rainbow(length(uniqueDomains))
    colors <- apply(col2rgb(colors), 2, function(x) { 0.3 + y/255 * 0.7 }) # make pastell colors
    colors <- apply(colors, 2, function(x) {rgb(x["red"], x["green"], x["blue"])}) # convert back to rgb
    # reassign colors
    names(colors) <- uniqueDomains
    retainedDomains1$color <- colors[retainedDomains1$proteinDomainID]
    retainedDomains2$color <- colors[retainedDomains2$proteinDomainID]
  }
  
  # reverse exons and protein domains, if on the reverse strand
  if (any(codingExons1$strand == "-")) {
    codingExons1$length <- rev(codingExons1$length)
    temp <- retainedDomains1$end
    retainedDomains1$end <- codingLength1 - retainedDomains1$start
    retainedDomains1$start <- codingLength1 - temp
  }
  if (any(codingExons2$strand == "-")) {
    codingExons2$length <- rev(codingExons2$length)
    temp <- retainedDomains2$end
    retainedDomains2$end <- codingLength2 - retainedDomains2$start
    retainedDomains2$start <- codingLength2 - temp
  }
  
  # normalize length to 1
  codingExons1$length <- codingExons1$length / (codingLength1 + codingLength2)
  codingExons2$length <- codingExons2$length / (codingLength1 + codingLength2)
  retainedDomains1$start <- retainedDomains1$start / (codingLength1 + codingLength2)
  retainedDomains1$end <- retainedDomains1$end / (codingLength1 + codingLength2)
  retainedDomains2$start <- retainedDomains2$start / (codingLength1 + codingLength2)
  retainedDomains2$end <- retainedDomains2$end / (codingLength1 + codingLength2)
  
  # draw coding regions
  rect(0, exonsY-exonHeight/2, sum(codingExons1$length), exonsY+exonHeight/2, col=color1, border=NA)
  rect(sum(codingExons1$length), exonsY-exonHeight/2, sum(codingExons1$length) + sum(codingExons2$length), exonsY+exonHeight/2, col=color2, border=NA)
  
  # indicate exon boundaries as dotted lines
  exonBoundaries <- cumsum(c(codingExons1$length, codingExons2$length))
  if (length(exonBoundaries) > 1) {
    exonBoundaries <- exonBoundaries[1:(length(exonBoundaries)-1)]
    for (exonBoundary in exonBoundaries)
      lines(c(exonBoundary, exonBoundary), c(exonsY-exonHeight, exonsY+exonHeight), col="white", lty=3)
  }
  
  # find overlapping domains
  # nest if one is contained in another
  # stack if they overlap partially
  nestDomains <- function(domains) {
    if (length(unlist(domains)) == 0) return(domains)
    domains <- domains[order(domains$end - domains$start, decreasing=T),]
    rownames(domains) <- 1:nrow(domains)
    # find nested domains and make tree structure
    domains$parent <- 0
    for (domain in rownames(domains))
      domains[domains$start >= domains[domain,"start"] & domains$end <= domains[domain,"end"] & rownames(domains) != domain,"parent"] <- domain
    # find partially overlapping domains
    maxOverlappingDomains <- max(1, as.integer(coverage(IRanges(domains$start*10e6, domains$end*10e6))))
    padding <- 1 / maxOverlappingDomains * 0.4
    domains$y <- 0
    domains$height <- 0
    adjustPositionAndHeight <- function(parentDomain, y, height, padding, e) {
      for (domain in which(e$domains$parent == parentDomain)) {
        overlappingDomains <- which((between(e$domains$start, e$domains[domain,"start"], e$domains[domain,"end"]) |
                                       between(e$domains$end  , e$domains[domain,"start"], e$domains[domain,"end"])) &
                                      e$domains$parent == parentDomain)
        e$domains[domain,"height"] <- height/length(overlappingDomains) - padding * (length(overlappingDomains)-1) / length(overlappingDomains)
        e$domains[domain,"y"] <- y + (which(domain==overlappingDomains)-1) * (e$domains[domain,"height"] + padding)
        adjustPositionAndHeight(domain, e$domains[domain,"y"]+padding, e$domains[domain,"height"]-2*padding, padding, e)
      }
    }
    adjustPositionAndHeight(0, 0, 1, padding, environment())
    domains <- domains[order(domains$height, decreasing=T),] # draw nested domains last
    return(domains)
  }
  retainedDomains1 <- nestDomains(retainedDomains1)
  retainedDomains2 <- nestDomains(retainedDomains2)
  retainedDomains1$y <- exonsY - exonHeight/2 + 0.025 + (exonHeight-2*0.025) * retainedDomains1$y
  retainedDomains2$y <- exonsY - exonHeight/2 + 0.025 + (exonHeight-2*0.025) * retainedDomains2$y
  retainedDomains1$height <- retainedDomains1$height * (exonHeight-2*0.025)
  retainedDomains2$height <- retainedDomains2$height * (exonHeight-2*0.025)
  
  # draw domains
  drawProteinDomainRect <- function(left, bottom, right, top, color) {
    rect(left, bottom, right, top, col=color, border=getDarkColor(color))
    # draw gradients for 3D effect
    gradientSteps <- 20
    drawVerticalGradient(rep(left, gradientSteps), rep(right, gradientSteps), seq(top, bottom, len=gradientSteps), rgb(1,1,1,0.7))
    drawVerticalGradient(rep(left, gradientSteps), rep(right, gradientSteps), seq(bottom, bottom+(top-bottom)*0.4, len=gradientSteps), rgb(0,0,0,0.1))
  }
  if (length(unlist(retainedDomains1)) > 0)
    for (domain in 1:nrow(retainedDomains1))
      drawProteinDomainRect(retainedDomains1[domain,"start"], retainedDomains1[domain,"y"], retainedDomains1[domain,"end"], retainedDomains1[domain,"y"]+retainedDomains1[domain,"height"], retainedDomains1[domain,"color"])
  if (length(unlist(retainedDomains2)) > 0)
    for (domain in 1:nrow(retainedDomains2))
      drawProteinDomainRect(sum(codingExons1$length)+retainedDomains2[domain,"start"], retainedDomains2[domain,"y"], sum(codingExons1$length)+retainedDomains2[domain,"end"], retainedDomains2[domain,"y"]+retainedDomains2[domain,"height"], retainedDomains2[domain,"color"])
  
  # draw gene names, if there are coding exons
  if (codingLength1 > 0)
    text(sum(codingExons1$length)/2, geneNamesY, fusion$gene1, font=2, cex=fontSize)
  if (codingLength2 > 0)
    text(sum(codingExons1$length)+sum(codingExons2$length)/2, geneNamesY, fusion$gene2, font=2, cex=fontSize)
  
  # calculate how many non-adjacent unique domains there are
  # we need this info to know where to place labels vertically
  countUniqueDomains <- function(domains) {
    uniqueDomains <- 0
    if (length(unlist(domains)) > 0) {
      uniqueDomains <- 1
      if (nrow(domains) > 1) {
        previousDomain <- domains[1,"proteinDomainID"]
        for (domain in 2:nrow(domains)) {
          if (previousDomain != domains[domain,"proteinDomainID"])
            uniqueDomains <- uniqueDomains + 1
          previousDomain <- domains[domain,"proteinDomainID"]
        }
      }
    }
    return(uniqueDomains)
  }
  if (length(unlist(retainedDomains1)) > 0)
    retainedDomains1 <- retainedDomains1[order(retainedDomains1$start),]
  uniqueDomains1 <- countUniqueDomains(retainedDomains1)
  if (length(unlist(retainedDomains2)) > 0)
    retainedDomains2 <- retainedDomains2[order(retainedDomains2$end, decreasing=T),]
  uniqueDomains2 <- countUniqueDomains(retainedDomains2)
  
  # draw title of plot
  titleY <- exonsY + exonHeight/2 + (uniqueDomains1 + 2) * 0.05
  text(0.5, titleY+0.01, "RETAINED PROTEIN DOMAINS", adj=c(0.5, 0), font=2, cex=fontSize)
  text(0.5, titleY, ifelse(fusion$reading_frame %in% c("in-frame", "out-of-frame"), paste(fusion$reading_frame, "fusion"), ifelse(fusion$reading_frame == "stop-codon", "stop codon before fusion junction", "reading frame unclear")), adj=c(0.5, 1), cex=fontSize)
  
  # draw domain labels for gene1
  if (length(unlist(retainedDomains1)) > 0) {
    previousConnectorX <- -1
    previousLabelX <- -1
    labelY <- exonsY + exonHeight/2 + uniqueDomains1 * 0.05
    for (domain in 1:nrow(retainedDomains1)) {
      # if possible avoid overlapping lines of labels
      connectorX <- min(retainedDomains1[domain,"start"] + 0.01, (retainedDomains1[domain,"start"] + retainedDomains1[domain,"end"])/2)
      if (connectorX - previousConnectorX < 0.01 && retainedDomains1[domain,"end"] > previousConnectorX + 0.01)
        connectorX <- previousConnectorX + 0.01
      labelX <- max(connectorX, previousLabelX) + 0.02
      # use a signle label for adjacent domains of same type
      adjacentDomainsOfSameType <- domain + 1 <= nrow(retainedDomains1) && retainedDomains1[domain+1,"proteinDomainID"] == retainedDomains1[domain,"proteinDomainID"]
      if (adjacentDomainsOfSameType) {
        labelX <- retainedDomains1[domain+1,"start"] + 0.015
      } else {
        text(labelX, labelY, retainedDomains1[domain,"proteinDomainName"], adj=c(0,0.5), col=getDarkColor(retainedDomains1[domain,"color"]), cex=fontSize)
      }
      lines(c(labelX-0.005, connectorX, connectorX), c(labelY, labelY, retainedDomains1[domain,"y"]+retainedDomains1[domain,"height"]), col=getDarkColor(retainedDomains1[domain,"color"]))
      if (!adjacentDomainsOfSameType)
        labelY <- labelY - 0.05
      previousConnectorX <- connectorX
      previousLabelX <- labelX
    }
  }
  
  # draw domain labels for gene2
  if (length(unlist(retainedDomains2)) > 0) {
    previousConnectorX <- 100
    previousLabelX <- 100
    labelY <- exonsY - exonHeight/2 - (uniqueDomains2+1) * 0.05
    for (domain in 1:nrow(retainedDomains2)) {
      # if possible avoid overlapping connector lines of labels
      connectorX <- sum(codingExons1$length) + max(retainedDomains2[domain,"end"] - 0.01, (retainedDomains2[domain,"start"] + retainedDomains2[domain,"end"])/2)
      if (previousConnectorX - connectorX < 0.01 && sum(codingExons1$length) + retainedDomains2[domain,"start"] < previousConnectorX - 0.01)
        connectorX <- previousConnectorX - 0.01
      labelX <- min(connectorX, previousLabelX) - 0.02
      # use a signle label for adjacent domains of same type
      adjacentDomainsOfSameType <- domain + 1 <= nrow(retainedDomains2) && retainedDomains2[domain+1,"proteinDomainID"] == retainedDomains2[domain,"proteinDomainID"]
      if (adjacentDomainsOfSameType) {
        labelX <- sum(codingExons1$length) + retainedDomains2[domain+1,"end"] - 0.015
      } else {
        text(labelX, labelY, retainedDomains2[domain,"proteinDomainName"], adj=c(1,0.5), col=getDarkColor(retainedDomains2[domain,"color"]), cex=fontSize)
      }
      lines(c(labelX+0.005, connectorX, connectorX), c(labelY, labelY, retainedDomains2[domain,"y"]), col=getDarkColor(retainedDomains2[domain,"color"]))
      if (!adjacentDomainsOfSameType)
        labelY <- labelY + 0.05
      previousConnectorX <- connectorX
      previousLabelX <- labelX
    }
  }
  
}

findExons <- function(exons, contig, geneID, direction, 
                      breakpoint, coverage, transcriptId) {
  message(paste0("Extract break points for ", geneID))
  # look for exon with breakpoint as splice site
  candidateExons <-  exons[exons$geneID == geneID & exons$contig == contig,]
  # if the gene has multiple transcripts, search for transcripts which
  # encompass the breakpoint
  if (length(unique(candidateExons$transcript)) > 1) {
    transcriptStart <- aggregate(candidateExons$start,
                                 by=list(candidateExons$transcript), min)
    rownames(transcriptStart) <- transcriptStart[,1]
    transcriptEnd <- aggregate(candidateExons$end, 
                               by=list(candidateExons$transcript), max)
    rownames(transcriptEnd) <- transcriptEnd[,1]
    encompassingExons <- between(breakpoint, 
                                 transcriptStart[candidateExons$transcript,2],
                                 transcriptEnd[candidateExons$transcript,2])
    
    if (any(encompassingExons))
      candidateExons <- candidateExons[encompassingExons,]
    }
  
  
  # find the consensus transcript, if there are multiple hits
  if (length(unique(candidateExons$transcript)) > 1) {
    consensusTranscript <-
      ifelse(grepl("appris_principal_1", candidateExons$attributes), 12,
             ifelse(grepl("appris_principal_2", candidateExons$attributes), 11,
                    ifelse(grepl("appris_principal_3", candidateExons$attributes), 10,
                           ifelse(grepl("appris_principal_4", candidateExons$attributes), 9,
                                  ifelse(grepl("appris_principal_5", candidateExons$attributes), 8,
                                         ifelse(grepl("appris_principal", candidateExons$attributes), 7,
                                                ifelse(grepl("appris_candidate_longest", candidateExons$attributes), 6,
                                                       ifelse(grepl("appris_candidate", candidateExons$attributes), 5,
                                                              ifelse(grepl("appris_alternative_1", candidateExons$attributes), 4,
                                                                     ifelse(grepl("appris_alternative_2", candidateExons$attributes), 3,
                                                                            ifelse(grepl("appris_alternative", candidateExons$attributes), 2,
                                                                                   ifelse(grepl("CCDS", candidateExons$attributes), 1,
                                                                                          0
                                                                                   ))))))))))))
    candidateExons <- candidateExons[consensusTranscript == max(consensusTranscript),]
  }
  # use the transcript with the longest coding sequence, if there are still multiple hits
  if (length(unique(candidateExons$transcript)) > 1) {
    codingSequenceLength <- ifelse(candidateExons$type == "CDS", candidateExons$end - candidateExons$start, 0)
    totalCodingSequenceLength <- aggregate(codingSequenceLength, by=list(candidateExons$transcript), sum)
    rownames(totalCodingSequenceLength) <- totalCodingSequenceLength[,1]
    candidateExons <- candidateExons[totalCodingSequenceLength[candidateExons$transcript,2] == max(totalCodingSequenceLength[,2]),]
  }
  # use the transcript with the longest overall sequence, if there are still multiple hits
  if (length(unique(candidateExons$transcript)) > 1) {
    exonLength <- candidateExons$end - candidateExons$start
    totalExonLength <- aggregate(exonLength, by=list(candidateExons$transcript), sum)
    rownames(totalExonLength) <- totalExonLength[,1]
    candidateExons <- candidateExons[totalExonLength[candidateExons$transcript,2] == max(totalExonLength[,2]),]
  }
  # if there are still multiple hits, select the first one
  candidateExons <- unique(candidateExons[candidateExons$transcript == head(unique(candidateExons$transcript), 1),])
  return(candidateExons)
}

coverageNormalization <- function(coverage, coverageRegion, exons) {
  max(1, ifelse(
    squishIntrons, # => ignore intronic coverage
    max(as.numeric(coverage[IRanges(sapply(exons$start,max,min(start(coverage))),sapply(exons$end,min,max(end(coverage))))])),
    round(quantile(coverage[coverageRegion], 0.9999)) # ignore coverage spikes from read-attracting regions
  ))
}

findClosestGene <- function(exons, contig, breakpoint, extraConditions) {
  
  # find exons near breakpoint (extraConditions must define what is considered "near")
  closestExons <- exons[exons$contig == contig & extraConditions,] # find closest exon
  closestExons <- exons[exons$contig == contig & exons$geneID %in% closestExons$geneID,] # select all exons of closest gene
  
  # when more than one gene found with the given name, use the closest one
  if (length(unique(closestExons$geneID)) > 1) { # more than one gene found with the given name => use the closest one
    distanceToBreakpoint <- aggregate(1:nrow(closestExons), by=list(closestExons$geneID), function(x) { min(abs(closestExons[x,"start"]-breakpoint), abs(closestExons[x,"end"]-breakpoint)) })
    closestGene <- head(distanceToBreakpoint[distanceToBreakpoint[,2] == min(distanceToBreakpoint[,2]),1], 1)
    closestExons <- closestExons[closestExons$geneID == closestGene,]
  }
  
  # when no gene was found, return default values
  if (nrow(closestExons) == 0) {
    return(IRanges(max(1, breakpoint-1000), breakpoint+1000))
  } else {
    return(IRanges(min(closestExons$start), max(closestExons$end)))
  }
}
