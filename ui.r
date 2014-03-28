library(shiny)
widget_style <- "display: inline-block; vertical-align: text-top; padding: 0px; border: solid;
                  border-width: 0px; border-radius: 0px; border-color: #CCC;"
shinyUI(pageWithSidebar(
    headerPanel("ENCODE RNAseq Evaluation"),
    sidebarPanel(
        wellPanel(
            div(style = widget_style,selectInput("protocol","Choose a protocal:",
                    list("PolyA_dUTP","PolyA_TruSeq","PolyA_SMARTseq", "Total_dUTP",
                         "Total_TruSeq","Total_SMARTseq"))),
            div(style = widget_style,selectInput("genetype","Choose a gene type:",
                    list("protein_coding", "pseudogene","lincRNA","antisense"))),
            div(style = widget_style,numericInput("cutD","Choose a cutoff D: (At FPKM 1: D= -4.67 for PolyA_dUTP; -4.68 for PolyA_TruSeq; -3.21 for PolyA_SMARTseq; -4.09 for Total_dUTP; -3.87 for Total_TruSeq; -2.95 for Total_SMARTseq)",-4.67,-10,15))
            ),
        wellPanel(
            div(style = widget_style,numericInput("xstart1", "x-axis start in SD plot:", -5,-30,10)),
            div(style = widget_style,numericInput("xend1", "x-axis end in SD plot:", 10,-10,20)),
            div(style = widget_style,numericInput("ystart1", "y-axis start in SD plot:", 0,0,10)),
            div(style = widget_style,numericInput("yend1", "y-axis end in SD plot:", 1,0.1,20))
            ),
        wellPanel(
            div(style = widget_style,numericInput("xstart2", "x-axis start in P0 plot:", -5,-30,-5)),
            div(style = widget_style,numericInput("xend2", "x-axis end in P0 plot:", 0,-5,20)),
            div(style = widget_style,numericInput("ystart2", "y-axis start in P0 plot:", 0,0,0.1)),
            div(style = widget_style,numericInput("yend2", "y-axis end in P0 plot:", 0.05,0.01,0.5))
            ),
        wellPanel(
            div(style = widget_style,numericInput("xstart3", "x-axis start in CAT plot:", 20,0,500)),
            div(style = widget_style,numericInput("xend3", "x-axis end in CAT plot:", 2000,100,10000)),
            div(style = widget_style,numericInput("ystart3", "y-axis start in CAT plot:", 0,0,1)),
            div(style = widget_style,numericInput("yend3", "y-axis end in CAT plot:", 1,0,1))
            ),
        wellPanel(
            div(style = widget_style,numericInput("xstart4", "x-axis start in CAT plot with array:", 20,0,500)),
            div(style = widget_style,numericInput("xend4", "x-axis end in CAT plot with array:", 2000,100,10000)),
            div(style = widget_style,numericInput("ystart4", "y-axis start in CAT plot with array:", 0,0,1)),
            div(style = widget_style,numericInput("yend4", "y-axis end in CAT plot with array:", 1,0,1))
            ),
        submitButton("Update View")
        ),
    mainPanel(
        h3(textOutput("caption")),
        tabsetPanel(
            tabPanel("SD_plot",plotOutput("sdplot",width="600px", height="600px")),
            tabPanel("P0_plot",plotOutput("p0plot",width="600px", height="600px")),
            tabPanel("CAT_plot",wellPanel(
                div(style = widget_style,plotOutput("catplot1",width="600px", height="600px")),
                div(style = widget_style,plotOutput("catplot2",width="600px", height="600px"))
                )),
            tabPanel("CAT_plot_array",
                     wellPanel(
                         div(style = widget_style,plotOutput("catplotarray1",width="600px", height="600px")),
                         div(style = widget_style,plotOutput("catplotarray2",width="600px", height="600px"))
                         ))
            ))))
