#' The application User-Interface
#'
#' @param request Internal parameter for `{shiny}`.
#'     DO NOT REMOVE.
#' @import shiny
#' @noRd
app_ui <- function(request) {
  tagList(
    # Leave this function for adding external resources
    golem_add_external_resources(),
    # Your application UI logic
    fluidPage(
      fluidRow(
        h1("ShinyGene"),
        shiny::column(
          width = 3,
          shiny::radioButtons(
            inputId = "SelectGeneRaw",
            label = "select gene data",
            choices = c("example data", "Upload"),
            selected = "example data"
          ),
          shiny::conditionalPanel(
            condition = "input.SelectGeneRaw == 'Upload'",
            shiny::wellPanel(
              shiny::fileInput(
                inputId = "UploadData",
                label = "upload CSV",
                multiple = F,
                accept = ".csv"),
              shiny::fileInput(
                inputId = "UploadDataCondition",
                label = "upload metadata CSV",
                multiple = F,
                accept = ".csv")
              )
            ),
          shiny::actionButton(
            inputId = "RunDEseq2",
            label = "Run Analysis"
            ),
          shinyWidgets::multiInput(
            inputId = "GeneSetFilter",
            label = "Gene Set :",
            choices = NULL,
            choiceNames = MetabolicPatways,
            choiceValues = MetabolicPatways
          ),
          shinyWidgets::multiInput(
            inputId = "GeneSymbolsFilter",
            label = "Gene Set :",
            choices = NULL,
            choiceNames = GeneSymbols,
            choiceValues = GeneSymbols
          )
        ),
        shiny::column(
          width = 9,
          shiny::tabsetPanel(
            type = "tabs",
            tabPanel(title = "Vulcano",
                     plotOutput(outputId = "VulacanoPlot"),
                     shinydashboard::box(collapsible = T,
                                  collapsed = T,
                                  width = 6,
                                  title = "FDR",
                                  "When you perform thousands of statistical tests (one for each gene),
                                  you will by chance call genes differentially expressed while they are not (false positives).
                                  You can control for this by applying certain statistical procedures called multiple hypothesis test correction.", br(),
                                  "We can count the number of genes that are differentially regulated at a certain alpha level.", br(),
                                  numericInput(inputId = "alphaFDR", label = "alpha level", value = 0.05, min = 0, max = 1, step = 0.0001, width = "20%"),
                                  textOutput(outputId = "nDiffExpGenes"),
                                  plotOutput(outputId = "pValDistribution")
                                  )
                     ),
            tabPanel(title = "Boxplots",
                     fluidRow(
                       column(
                         width = 3,
                         uiOutput(outputId = "Genes4Boxplot_UI"))
                       ),
                     fluidRow(
                       column(
                         width = 9,
                         plotOutput(outputId = "BoxPlotGenes")),
                       column(
                         width = 3,
                         textOutput(outputId = "Genes2PlotText"))
                     )
                     ),
            tabPanel(title = "Results", DT::DTOutput(outputId = "cond"))
          )
        )
      )
    )
  )
}

#' Add external Resources to the Application
#'
#' This function is internally used to add external
#' resources inside the Shiny application.
#'
#' @import shiny
#' @importFrom golem add_resource_path activate_js favicon bundle_resources
#' @noRd
golem_add_external_resources <- function() {
  add_resource_path(
    "www",
    app_sys("app/www")
  )

  tags$head(
    favicon(),
    bundle_resources(
      path = app_sys("app/www"),
      app_title = "ShinyGene"
    ),
    # Add here other external resources
    # for example, you can add shinyalert::useShinyalert()
    shinyWidgets::useShinydashboard()
  )
}
