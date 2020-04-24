
shinyUI(
        navbarPage("SCAFFoLD",
            tabPanel("Map exploration", uiOutput("graphUI")),
            tabPanel("Run SCAFFoLD Analysis", uiOutput("analysisUI")),
            tabPanel("Run clustering", uiOutput("clusteringUI")),
            tabPanel("Add frequency statistics", uiOutput("freqstatsUI")),
            tabPanel("Add expression statistics", uiOutput("exprstatsUI")),
            tabPanel("Add correlation statistics", uiOutput("coorelationUI")),
            tabPanel("Map dataset", uiOutput("mappingUI"))
    )
)





