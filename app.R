library(shiny)
library(ggplot2)
#Declaring variables for loops
x<-NULL
y<-NULL
groups<-NULL
weightsFac<-NULL

getDat <- function(nGroups,groupN,xScatter,popInt,popSlopeL1,popSlopeL2,popInter,intVar,slopeVar,noise){
    #Generating (# of groups * size of groups) x values
    x <- rnorm(nGroups*groupN,0,1)
    
    #Generating group values
    rSlopes <- rnorm(nGroups,popSlopeL1,slopeVar) #Generating a slope for each group, centered around the level 1 population level slope
    rInts <- rnorm(nGroups,popInt,intVar) #Generating an intercept for each group, centered around the population intercept
    xOffset<-rnorm(nGroups,0,xScatter) #Creating a value to offset each group's x variable - needed to create distinct patches of data
    gender<-c(0,1)
    weights<-round(rnorm(nGroups,0,1),0) #our Level 2 variable for each group
    
    for (i in 1:nGroups) {
        #Figuring out the boundaries for the group
        iMax = groupN*i 
        iMin = iMax-(groupN)+1
        
        #Applying our x offset
        x[iMin:iMax]<-x[iMin:iMax]+xOffset[i]
        
        #Creating y values from input values
        y[iMin:iMax]<-((x[iMin:iMax]*rSlopes[i])+(rInts[i])+(weights[i]*popSlopeL2)+(noise*rnorm(groupN,0,1)+(weights[i]*x[iMin:iMax]*popInter)))
        
        #Associating a group and weight with each data point
        groups[iMin:iMax]<-rep(i,groupN)
        weightsFac[iMin:iMax]<-rep(weights[i],groupN)
    }
    
    #Adding generated values to a data frame and returning it
    df <- data.frame("Level1"=x,y=y,"groups"=as.factor(groups),"Level2"=weightsFac)
    return(df)
}
ui <- fluidPage(
    titlePanel("Clustering Algorithm Tool"),
    sidebarLayout(
        sidebarPanel(
            actionButton("save", 
                         "Save Data as .CSV",
                         width="100%"),
            sliderInput("nGroups",
                        "Number of groups:",
                        min = 1,
                        max = 50,
                        value = 10),
            sliderInput("groupN",
                        "Size of groups:",
                        min = 1,
                        max = 1000,
                        value = 200),
            sliderInput("xScatter",
                        "Group X SD:",
                        min = 0,
                        max = 25,
                        value = 5),
            sliderInput("popInt",
                        "Population Intercept:",
                        min = 0,
                        max = 100,
                        value = 0),
            sliderInput("popSlopeL1",
                        "Population Slope (L1):",
                        min = -5,
                        max = 5,
                        value = 1.0,
                        step = 0.05),
            sliderInput("popSlopeL2",
                        "Population Slope (L2):",
                        min = -5,
                        max = 5,
                        value = 0.05,
                        step = 0.05),
            sliderInput("popInter",
                        "L1xL2 Interaction",
                        min = -5,
                        max = 5,
                        value = 0.05,
                        step = 0.05),
            sliderInput("intVar",
                        "Intercept SD:",
                        min = 0,
                        max = 100,
                        value = 10),
            sliderInput("slopeVar",
                        "Slope SD:",
                        min = 0,
                        max = 5,
                        value = 0.05,
                        step = 0.05),
            sliderInput("noise",
                        "L1 Error:",
                        min = 0,
                        max = 50,
                        value = 5,
                        step = 1)),
        mainPanel(
            tabsetPanel(
                tabPanel("Plots",
                         fluidRow(
                             column(6,
                                    plotOutput("grp"),
                                    plotOutput("have"),
                                    plotOutput("hsing")
                                    ),
                             column(6,
                                    plotOutput("kmeans"),
                                    plotOutput("hcomp"),
                                    plotOutput("hcent")))),
                tabPanel("Hierarchical Clustering Trees",
                         fluidRow(
                             column(6,
                                    plotOutput("Ha"),
                                    plotOutput("Hsi")
                                    ),
                             column(6,
                                    plotOutput("Hco"),
                                    plotOutput("Hce"))))
                )
            )
        )
    )
server <- function(input, output, session) {
    observeEvent(input$save, {
        write.csv(makeDF(),"data.csv")
        session$sendCustomMessage(type = 'Message',
                                  message = paste("File written to",getwd()))
    })
    makeDF <- reactive({
        getDat(input$nGroups,input$groupN,input$xScatter, input$popInt,input$popSlopeL1,input$popSlopeL2, input$popInter, input$intVar,input$slopeVar,input$noise)
    })
    output$grp <- renderPlot({
        plot <- ggplot(data=makeDF(),aes(x=Level1,y=y,color=groups))+
            geom_point()+         
            geom_hline(yintercept=0,lty=2)+
            geom_vline(xintercept=0,lty=2)+
            theme_classic()+
            ggtitle("Actual Groups")
        plot
    })
    output$kmeans <- renderPlot({
        dat<-makeDF()
        dat.k<-kmeans(dat,input$nGroups,100,25)
        plot <- ggplot(data=dat,aes(x=Level1,y=y,color=as.factor(dat.k$cluster)))+
            geom_point()+         
            geom_hline(yintercept=0,lty=2)+
            geom_vline(xintercept=0,lty=2)+
            theme_classic()+
            ggtitle("Kmeans Groups")
        plot
    })
    output$hcomp <- renderPlot({
        dat<-makeDF()
        dat.hc<-cutree(hclust(dist(dat),method="complete"),input$nGroups)
        plot <- ggplot(data=makeDF(),aes(x=Level1,y=y,color=as.factor(dat.hc)))+
            geom_point()+         
            geom_hline(yintercept=0,lty=2)+
            geom_vline(xintercept=0,lty=2)+
            theme_classic()+
            ggtitle("Hierarchical Clustering (complete)")
        plot
    })
    output$have <- renderPlot({
        dat<-makeDF()
        dat.hc<-cutree(hclust(dist(dat),method="average"),input$nGroups)
        plot <- ggplot(data=makeDF(),aes(x=Level1,y=y,color=as.factor(dat.hc)))+
            geom_point()+         
            geom_hline(yintercept=0,lty=2)+
            geom_vline(xintercept=0,lty=2)+
            theme_classic()+
            ggtitle("Hierarchical Clustering (average)")
        plot
    })
    output$hsing <- renderPlot({
        dat<-makeDF()
        dat.hc<-cutree(hclust(dist(dat),method="single"),input$nGroups)
        plot <- ggplot(data=makeDF(),aes(x=Level1,y=y,color=as.factor(dat.hc)))+
            geom_point()+         
            geom_hline(yintercept=0,lty=2)+
            geom_vline(xintercept=0,lty=2)+
            theme_classic()+
            ggtitle("Hierarchical Clustering (single)")
        plot
    })
    output$hcent <- renderPlot({
        dat<-makeDF()
        dat.hc<-cutree(hclust(dist(dat),method="centroid"),input$nGroups)
        plot <- ggplot(data=makeDF(),aes(x=Level1,y=y,color=as.factor(dat.hc)))+
            geom_point()+         
            geom_hline(yintercept=0,lty=2)+
            geom_vline(xintercept=0,lty=2)+
            theme_classic()+
            ggtitle("Hierarchical Clustering (centroid)")
        plot
    })
    output$Ha <- renderPlot({
        dat<-makeDF()
        dat.ha<-hclust(dist(dat),method="average")
        plot(dat.ha)
    })
    
    output$Hco <- renderPlot({
        dat<-makeDF()
        dat.hco<-hclust(dist(dat),method="complete")
        plot(dat.hco)
    })
    
    output$Hsi <- renderPlot({
        dat<-makeDF()
        dat.hsi<-hclust(dist(dat),method="single")
        plot(dat.hsi)
    })
    
    output$Hce <- renderPlot({
        dat<-makeDF()
        dat.hce<-hclust(dist(dat),method="centroid")
        plot(dat.hce)
    })
}
shinyApp(ui = ui, server = server)
