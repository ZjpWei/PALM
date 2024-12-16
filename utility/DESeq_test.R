test.DESeq <- function(data, label){
  # based on vignette in phyloseq
  # https://www.bioconductor.org/packages/release/bioc/vignettes/phyloseq/inst/doc/phyloseq-mixture-models.html
  stopifnot(ncol(data) == length(label))
  stopifnot(all(sort(unique(label)) == c(0, 1)))

  # convert to phyloseq object
  temp_otu <- phyloseq::otu_table(data, taxa_are_rows = TRUE)
  temp_sample <- as.matrix(label)
  rownames(temp_sample) <- colnames(data)
  colnames(temp_sample) <- 'label'
  temp_sample <- data.frame(temp_sample)
  temp_sample$label <- factor(temp_sample$label)
  physeq <- phyloseq::phyloseq(phyloseq::otu_table(temp_otu),
                               phyloseq::sample_data(temp_sample))
  ####
  # taken from phyloseq vignette
  diagdds = phyloseq::phyloseq_to_deseq2(physeq, ~label)
  # calculate geometric means prior to estimate size factors
  gm_mean = function(x, na.rm=TRUE){
    exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
  }
  geoMeans = apply(DESeq2::counts(diagdds), 1, gm_mean)
  diagdds = DESeq2::estimateSizeFactors(diagdds, geoMeans = geoMeans)
  diagdds = DESeq2::DESeq(diagdds, fitType="local")
  res = DESeq2::results(diagdds)
  return(res)
}
