



med_dat <- read.delim('C:/Users/abig4/OneDrive/Documents/GitHub/Statistical-Learning_Lab_2/gtex_Kmeans/gtex.gct',
                      skip = 2 ,row.names=c(1) , header = TRUE)
#med_dat <- read.delim("C:/Users/Shahar/Documents/GitHub/Statistical-Learning_Lab_2/gtex_Kmeans/gtex.gct",
#                      skip = 2 ,row.names=c(1) , header = TRUE)
gen_names <- med_dat[, 1]
med_dat <- med_dat[,-1]
ti_name <- colnames(med_dat)
med_dat <- transpose(med_dat)
med_dat <-med_dat[,colMeans(med_dat) > 0]
v <- apply(med_dat,2, var)
med_dat <-med_dat[, v != 0]
med_dat <- as.data.frame(scale(med_dat,T,T))
m <- prcomp(med_dat,scale. = T,center = T)
pca_p <- m$x
sv <-m$sdev^2/sum(m$sdev^2)
sv <- cumsum(sv)
sv <- length(sv[sv < 0.95])
m <- m$x[,1:sv]

m<-med_dat

k_means_shiny <- function(k){
    data_prep <- m
    m_new <- as.data.frame(data_prep[sample(1:53,k),])
    clus <- apply(data_prep,1 ,FUN = c_fun, m=m_new)
    m_old <- as.data.frame(data_prep[sample(1:53,k),])
    data_prep <- cbind(data_prep, clus)
    iter <- 1
    while(any(abs(m_new - m_old)) > 0.1 & iter < 1){
        data_prep <- as.data.frame(data_prep)
        m_old <- m_new
        m_new <- aggregate(. ~ clus, data_prep, FUN = mean)
        m_new <- m_new[,-1]
        data_prep$clus<- data_prep[,-1]  %>% apply(1,FUN = c_fun,m = m_new)
        print(m_new)
        print(m_old)
        iter <- iter +1
        
    }
    return(data_prep)
}

c_fun <- function(d,m){
    p <- NULL
    for (i in 1:length(m[,1])) {
        x <- rbind(d,m[i,])
        di <- dist(x, method = 'euclidean')
        p <- cbind(p,di)
    }
    p <- which.min(p)
    return(p)
}



library(shiny)
library(ggplot2)

# Define UI for application that draws a histogram
ui <- fluidPage(

    # Application title
    titlePanel("Exploratory analysis of RNA seq data"),

    # Sidebar with a slider input for number of bins 
    sidebarLayout(
        sidebarPanel(
            sliderInput("bins",
                        "Number of clusters:",
                        min = 1,
                        max = 53,
                        value = 7)
        ),

        
        mainPanel(
           plotOutput("scatterplot")
        )
    )
)

# Define server logic required to draw a histogram
server <- function(input, output) {

    output$scatterplot <- renderPlot({
        # generate bins based on input$bins from ui.R
        x <- k_means_shiny(input$bins)
        x <-cbind(pca_p,x["clus"])
        ggplot(data = x, aes(x =PC1, y=PC2, color = factor(clus))) + geom_point()

    })
}

# Run the application 
shinyApp(ui = ui, server = server)
