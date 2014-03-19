library(shiny)
widget_style1 <-  "display: inline-block; vertical-align: text-top; padding: 7px; border: solid;
                  border-width: 1px; border-radius: 4px; border-color: #CCC;"
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
                    list("protein_coding", "pseudogene","lincRNA","antisense")))
            ),
        wellPanel(
            div(style = widget_style,textInput("rsemtpm","RSEM TPM label:","rsem")),
            div(style = widget_style,numericInput("cutoff1", "Threshold for RSEM TPM:", 0,0,100)),
            div(style = widget_style,textInput("rsempmetpm","RSEM PMETPM label:","rsem_pme")),
            div(style = widget_style,numericInput("cutoff2", "Threshold for RSEM PMETPM:", 0,0,100)),
            div(style = widget_style,textInput("fluxcapacitor","Flux Capacitor label:","flux")),
            div(style =widget_style,numericInput("cutoff3","Threshold for Flux Capacitor:", 0,0,100)),
            div(style = widget_style,textInput("cufflinkss","Cufflinks STAR label:","cuff_s")),
            div(style =widget_style,numericInput("cutoff4","Threshold for Cufflinks STAR:", 0,0,100)),
            div(style = widget_style,textInput("cufflinkst","Cufflinks TopHat label:","cuff_t")),
            div(style=widget_style,numericInput("cutoff5","Threshold for Cufflinks TopHat:",0,0,100)),
            div(style = widget_style,textInput("sailfish","Sailfish label:","sailfish")),
            div(style = widget_style,numericInput("cutoff6","Threshold for Sailfish:",0,0,100)),
            div(style = widget_style,textInput("express","eXpress label:","express")),
            div(style = widget_style,numericInput("cutoff7","Threshold for eXpress:",0,0,100)),
            div(style = widget_style,textInput("naive","Naive label:","naive")),
            div(style = widget_style,numericInput("cutoff8","Threshold for Naive:",0,0,100))
            ),
        wellPanel(
            div(style = widget_style,numericInput("xstart1", "x-axis start in SD plot:", -15,-30,10)),
            div(style = widget_style,numericInput("xend1", "x-axis end in SD plot:", 10,-10,20)),
            div(style = widget_style,numericInput("ystart1", "y-axis start in SD plot:", 0,0,10)),
            div(style = widget_style,numericInput("yend1", "y-axis end in SD plot:", 2,0.1,20))
            ),
        wellPanel(
            div(style = widget_style,numericInput("xstart2", "x-axis start in P0 plot:", -15,-30,-5)),
            div(style = widget_style,numericInput("xend2", "x-axis end in P0 plot:", 0,-5,20)),
            div(style = widget_style,numericInput("ystart2", "y-axis start in P0 plot:", 0,0,0.1)),
            div(style = widget_style,numericInput("yend2", "y-axis end in P0 plot:", 0.07,0.02,0.5))
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
