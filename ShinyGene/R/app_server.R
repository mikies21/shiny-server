#' The application server-side
#'
#' @param input,output,session Internal parameters for {shiny}.
#'     DO NOT REMOVE.
#' @import shiny
#' @noRd
app_server <- function(input, output, session) {
  # Your application server logic
  GeneRaw <- shiny::reactive({
    if (input$SelectGeneRaw == "example data") {
      RAvHCcpm
    } else {
      file <- input$UploadData
      ext <- tools::file_ext(file$datapath)

      req(file)
      validate(need(ext == "csv", "Please upload a csv file"))

      read.csv(file$datapath, header = T, row.names = 1)
    }

  })

  ConditionData <- shiny::reactive({
    if (input$SelectGeneRaw == "example data") {
      ConditionHCvRA

    } else {
      file <- input$UploadDataCondition
      ext <- tools::file_ext(file$datapath)

      req(file)
      validate(need(ext == "csv", "Please upload a csv file"))

      read.csv(file$datapath, header = T, row.names = 1)
    }
  })

  dds <- eventReactive(input$RunDEseq2, {

    dds <- DESeq2::DESeqDataSetFromMatrix(countData = round(GeneRaw()),
                                  colData = ConditionData(),
                                  design = ~condition)
    keep <- rowSums(DESeq2::counts(dds)) >= 10
    dds <- dds[keep,]
    dds <- DESeq2::DESeq(dds)
    dds
  })

  ddsRes <- reactive({
    res <- DESeq2::results(dds())
    res
  })

  GeneSet <- reactive({
    set <- subset(x = ReactomeGeneSet, subset = gs_description %in% input$GeneSetFilter)$gene_symbol
    symbol <- subset(x = ReactomeGeneSet, subset = gene_symbol %in% input$GeneSymbolsFilter)$gene_symbol
    GeneSet <- unique(c(symbol, set))
  })


  output$cond <- DT::renderDT({
    as.data.frame(ddsRes())
  })


  output$VulacanoPlot <- shiny::renderPlot({
    dat <- as.data.frame(ddsRes()) |>
      tibble::rownames_to_column("symbol") |>
      dplyr::mutate(sig = dplyr::case_when(padj < 0.05 & log2FoldChange > log2(1.5) ~ "UP regulated",
                                           padj < 0.05 & log2FoldChange < -log2(1.5) ~ "DOWN regulated",
                                           T ~ "not significant"))

    p <- ggplot2::ggplot()+
      ggplot2::geom_point(data = dat,
                          ggplot2::aes(x = log2FoldChange, y = -log10(pvalue), colour = sig),
                          alpha = 0.4)+
      ggplot2::geom_hline(yintercept = -log10(0.05), linetype = "dotted")+
      ggplot2::geom_vline(xintercept = c(log2(1.5), -log2(1.5)), linetype = "dotted")

    if (length(input$GeneSetFilter) > 0 |
        length(input$GeneSymbolsFilter) > 0) {

      SelectedGenes <- dat |>
        dplyr::filter(symbol %in% GeneSet())

      p <- p +
        ggplot2::geom_point(data = SelectedGenes,
                            ggplot2::aes(x = log2FoldChange, y = -log10(pvalue)),
                            colour = "black")+
        ggrepel::geom_text_repel(data = SelectedGenes,
                                 ggplot2::aes(x = log2FoldChange, y = -log10(pvalue), label = symbol),
                                 colour = "black")

    }

    p +
      ggplot2::scale_colour_manual(values = c("UP regulated" = "darkred",
                                              "DOWN regulated" = "blue",
                                              "not significant" = "grey"))+
      ggplot2::labs(x = "Log2 fold change", y = "-Log10(p-value)") +
      ggprism::theme_prism()

  })

  output$nDiffExpGenes <- renderText({
    diffExpGenes <- as.data.frame(ddsRes()) |>
      dplyr::filter(padj < input$alphaFDR) |>
      nrow()

    paste(sprintf("%s differentially expressed genes at padj %s,\n
                it corresponds to respectively %s per cent of the whole number transcriptome (total number of mRNA is %s)",
                diffExpGenes, input$alphaFDR, round(ncol(as.data.frame(ddsRes()))/diffExpGenes, digits = 2), nrow(as.data.frame(ddsRes()))))
  })

  output$pValDistribution <- renderPlot({
    padj <- ggplot2::ggplot(as.data.frame(ddsRes()))+
      ggplot2::geom_histogram(ggplot2::aes(x = padj), fill = "lightblue", colour = "black")+
      ggprism::theme_prism()+
      ggplot2::labs(x = "adjusted pvalue", y = "Frequency", title = "Adjusted p-value distribution")

    p <- ggplot2::ggplot(as.data.frame(ddsRes()))+
      ggplot2::geom_histogram(ggplot2::aes(x = pvalue), fill = "grey", colour = "black")+
      ggprism::theme_prism()+
      ggplot2::labs(x = "adjusted pvalue", y = "Frequency", title = "Non-adjusted p-value distribution")

    cowplot::plot_grid(padj, p)
  })

  output$Genes4Boxplot_UI <- renderUI({
    shinyWidgets::pickerInput(
      inputId = "Genes4Boxplot",
      label = "Gene symbols",
      choices = GeneSet(),
      multiple = T,
      options = list(
        `live-search` = TRUE,
        `actions-box` = TRUE)
    )
  })

  output$Genes2PlotText <- renderText({
    paste0(input$Genes4Boxplot, collapse = ", ")
  })



  output$BoxPlotGenes <- renderPlot({
    dat <- merge(ConditionData(), t(subset(x = GeneRaw(), subset = rownames(GeneRaw()) %in% input$Genes4Boxplot)), by = 0) |>
      dplyr::rename(sampleID = Row.names) |>
      tidyr::pivot_longer(!c(sampleID, condition), names_to = "symbol", values_to = "value")

    pvals <- dat |>
      dplyr::group_by(symbol) |>
      rstatix::t_test(value ~ condition) |>
      rstatix::adjust_pvalue(method = "BH") |>
      rstatix::add_y_position(fun = "max")|>
      dplyr::arrange(dplyr::desc(symbol))

    pvals_filtered <- as.data.frame(ddsRes()) |>
      tibble::rownames_to_column("symbol") |>
      dplyr::filter(symbol %in% input$Genes4Boxplot) |>
      dplyr::arrange(dplyr::desc(symbol)) |>
      tibble::tibble()

    pvals$p <- pvals_filtered$pvalue
    pvals$p.adj <- pvals_filtered$padj

    pvals <- pvals |>
      rstatix::add_significance()

    ggpubr::ggboxplot(data = dat, x = "condition", y = "value", color = "condition", palette = "aaas", facet.by = "symbol", scales = "free_y", outlier.shape = NA, add = "jitter")+
      ggpubr::stat_pvalue_manual(pvals)+
      ggplot2::facet_wrap(~symbol, scales = "free_y")+
      ggprism::theme_prism()
      #ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.05, 0.1)))
  })

}
