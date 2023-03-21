####### WORKING Version
#### App for proteomics_aged_c57bl6_shiny
#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#
source("./code/Proteomics_Functions.R")

library(shiny)
library(plotly)
library(DT)
library(dplyr)
library(stringi)
library(ggrepel)
library(gplots)
library(shinyalert)
#library(ggplot2)


# Define UI 
ui <- navbarPage(strong("Aging B6 Proteomics"),		# Application title
                 
                 # Sidebar with a slider input for number of bins 
                 tabPanel( "Effects on proteins (boxplots)",
                           sidebarLayout(
                             sidebarPanel(width = 2,
                                          
                                          h4(strong("Protein(s) Analysis")), 
                                          textInput("Pname","Gene Name",""),
                                          radioButtons(inputId = "checkbox1",label = NULL, choices = c("Ensembl ID" = 1, "Gene Symbol" = 2), selected = 2),
                                          
                                          selectInput("boxT", "Boxplot Type:",
                                                      c("Sex by Age" = "SA",
                                                        "Sex only" = "sex",
                                                        "Age only" = "age"), selected = "SA"),
                                          
                                          actionButton("geneeffect", "Query"),
                                          h6("Note: Click for sample details in the plot",style = "color:blue"),
										  hr()
                             ),
                             
                             # Show plots from analysis    character(0)
                             
                             mainPanel(
                               h3("Box plots of Sex-Age effects for protein(s) in query"),
                               plotOutput("PBoxPlot",click = "plot_click",width = "120%",height = "auto"),
                               #plotlyOutput("PBoxPlot",width = "auto",height = "auto"),
                               fluidRow(
                                 column(dataTableOutput("click_info"), width = 12)
                               ),
                               
                               #########################################							  
                               hr()
                             )
                           )
                 ),
                 
                 
                 
                 tabPanel( "Effects across proteins (volcano plots)",
                           sidebarLayout(
                             sidebarPanel(width = 2,
                                          
                                          h4(strong("Age/Sex Effect Analysis")),
                                          selectInput("group1", "Factor:",
                                                      c("Age" = "Age", "Sex" = "Sex"),selected = "Age"),
                                          numericInput("qval","FDR Threshold:", 0.1, min = 0,max = 1,step = 0.05),
                                          h5(strong("Highlight Option:")),
                                          h5(strong("Enter 'Multiple Proteins' (comma sparated, no space; example: Igkc,Bcat1,Vcam1) or Select 'Annotated Sets'")),
                                          
                                          textInput("Pset1","1) Multiple Proteins:",""),
                                          selectizeInput('Pset2','2) Annotated Sets:',choices = NULL,selected = character(0)),
                                          #selectizeInput('Pset2','gs_name', choices = NULL, selected = NULL),
                                          #selectizeInput('Pset2', 'AND', choices = c("choose" = "", levels(gene_sets_dat$gs_name))),
                                          
                                          #radioButtons(inputId = "checkbox2",label = NULL, choices = c("Multiple Proteins" = 1, "Annotated Sets " = 2), selected = 2),
                                          actionButton("factoreffect", "Run"),
                                          h6("Note: Click for details or drag to zoom-in in the plot",style = "color:blue")
                             ),
                             
                             # Show plots from analysis    character(0)
                             
                             mainPanel(
                               h3("Volcano plots for Sex/Age effects on proteins across tissues"),
                               
                               #plotOutput("PvolPlot1",height = "900px"),
                               plotlyOutput("PvolPlot", height = "800px"),
                               
                               fluidRow(
                                 textOutput("selected_title"),
                                 plotOutput("MSPBoxPlot",width = "110%"),
                                 tags$head(tags$style("#selected_title{color: black;
                                    font-size: 15px; font-style: italic;}"))
                                 
                               )
                               
                             )
                           )
                 ),
                 
                 
                 
                 tabPanel( "Effects across tissues (heatmap)",
                           sidebarLayout(
                             sidebarPanel(width = 2,
                                          
                                          h4(strong("Age/Sex Effect Analysis")),
                                          selectInput("group2", "Factor:",
                                                      c("Age" = "Age", "Sex" = "Sex"),selected = "Age"),
                                          numericInput("Hqval","FDR Threshold to filter proteins:", 0.1, min = 0,max = 1,step = 0.05),
                                          
                                          actionButton("cluster", "Run Analysis"),
                                          h6("Note: Click and drag to zoom-in in the plot",style = "color:blue")
                                          
                             ),
                             
                             # Show a plot of the generated distribution
                             mainPanel(
                               h3("Heatmap for Sex/Age effects on proteins across tissues"), 
                               plotlyOutput("HeatmapPlot")
                             )
                           )
                 ),
                 
                 tabPanel( "Effect comparison between tissues (scatter plot)",
                           sidebarLayout(
                             sidebarPanel(width = 2,
                                          
                                          h4(strong("Age/Sex effect comparison between two tissues")),
                                          
                                          selectInput("group3", "Factor for comparison:",
                                                      c("Age" = "Age", "Sex" = "Sex"),selected = "Age"),
                                          
                                          
                                          selectizeInput(
                                            "group4",
                                            label = "Select Two Tissues:",
                                            choices = c(" ", "kidney","hippocampus","heart","fat","cerebellum","striatum","spleen","skeletalmuscle","lung","liver"),
                                            multiple = TRUE,
                                            options = list(maxItems = 2)
                                          ),
                                          ## select pathway
                                          uiOutput("pathway_select"),
                                          #numericInput("Hqval","FDR Threshold to filter proteins:", 0.1, min = 0,max = 1,step = 0.05),
                                          
                                          actionButton("Comp", "Compare"),
                                          h6("Note: Click and drag to zoom-in in the plot",style = "color:blue")
                                          
                             ),
                             
                             
                             # Show a plot of the generated distribution
                             mainPanel(
                               #h3("Heatmap for Sex/Age effects on proteins across tissues"), 
                               plotlyOutput("ComPlot", height = "800px")
                               #plotOutput("ComPlot")
                               #fluidRow(
                                 
                                 #column(dataTableOutput("click_info"), width = 12)
                               #)
                             )
                           )
                 ),
				 
				 tabPanel(downloadButton("downloadData", "Download All Data", style = "color: #fff; background-color: #27ae60; border-color: #fff;margin-top:-30px;margin-bottom: 0px"))
                 
                 
)

# Define server
server <- function(input, output, session) {
  
  # Downloadable csv of selected dataset ----
  output$downloadData <- downloadHandler(
    filename = function() {
      paste("protein_all_data", ".zip", sep = "")
    },
    content = function(file) {
      #write.csv(protein_all_tissues_plotting_dat, file, row.names = FALSE)
	  file.copy("./data/protein_all_data.zip", file)
    }
  )
    
  updateSelectizeInput(session, 'Pset2',
                       choices = c(unique(gene_sets_dat$gs_name)),
                       selected = character(0),
                       server = TRUE
  )
  
  PD_S <- reactive({
    
    Pname <- (unlist(strsplit(input$Pname, ","))) 
    if (input$checkbox1 == 1) {PD = subset(protein_all_tissues_plotting_dat,gene.id %in% stri_trans_toupper(Pname))}				### type==1: Ensembl ID
    if (input$checkbox1 == 2) {PD = subset(protein_all_tissues_plotting_dat,symbol %in% stri_trans_totitle(stri_trans_tolower(Pname)))}	
    if(nrow(PD)==0){
      shinyalert("Oops!", "No such gene", type = "error")
      PD <- NULL
    }	
	return(PD)
  })
  
  D = reactive({
    req(PD_S())
    PD <- PD_S()	
    if (is.null(input$plot_click)) return()
    PD$symbol <- toupper(PD$symbol)
    if(input$boxT=='SA'){PDT = subset(PD,tissue == input$plot_click$panelvar2 & symbol == input$plot_click$panelvar3)}
	if(input$boxT!='SA'){PDT = subset(PD,tissue == input$plot_click$panelvar1 & symbol == input$plot_click$panelvar2)}
    
    DTE <- nearPoints(PDT, input$plot_click,allRows = T,addDist = TRUE)  ##nearCountry <- nearPoints(plotData(), input$plot_click, x='Sex', y='Intensity',maxpoints = 1)
    DTE = subset(DTE, dist_ < 30)
    DTE <- DTE[order(DTE$dist_),]
    return(DTE)			
  })
  # ,"Sex","Intensity"
  
  observeEvent(input$geneeffect, {
    req(PD_S())
	PD <- PD_S()
    dd <- Sex_Age_Effect_2(PD,tissue_specific_results_long,input$boxT)
    #dd <- Sex_Age_Effect(protein_all_tissues_plotting_dat, input$Pname, input$checkbox1)
    #output$PBoxPlot <- renderPlotly(dd$S_A_plot) 
    number_of_plot <- length(unique(PD$tissue))
    number_of_gene <- length(unique(PD$symbol))
    #print(number_of_gene)
    TW = "auto"
    TH = 400#
    if ( number_of_plot < 6) {TW = 300*number_of_plot}
    if ( number_of_gene > 1) {TH = 250*number_of_gene}
    output$PBoxPlot <- renderPlot({dd$S_A_plot},width = TW, height = TH,res = 105) 
    output$click_info <- DT::renderDataTable(data.frame())
    #output$test_A <- renderDataTable(dd$result)		 
  })
  
  
  observeEvent(input$plot_click, 
               {
                 req(PD_S())
				 PD <- PD_S()
                 DD <- D()
                 PD$Psymbol <- toupper(PD$symbol) 
                 output$click_info <- DT::renderDataTable(
                   subset(PD, mouse.id == DD$mouse.id[1] & Psymbol == DD$symbol[1])[,-ncol(PD)], filter = 'top', server = FALSE, ### subset(PD, mouse.id == DD$mouse.id[1])
                   options = list(pageLength = 10, autoWidth = TRUE), rownames = FALSE
                 )
               })
  #output$click_info <- renderDataTable(subset(D(),mouse.id == D()$mouse.id))			
  #output$click_info <- renderDataTable(D())
  
  #	
  
  
  observeEvent(input$factoreffect, {
    
    showModal(modalDialog("Working On It! Just one minute", footer = NULL))
    
    #if (input$group1 == "Age") {dh1 <- Age_Effect(tissue_specific_results,tissue_specific_results_long)}
    #if (input$group1 == "Sex") {dh1 <- Sex_Effect_C(tissue_specific_results,tissue_specific_results_long)}		 
    #output$PvolPlot1 <- renderPlot({dh1$SA_plot}) 
    
    Sqval = input$qval
    if (input$qval > 1) {Sqval = 1}
    if (input$qval < 0) {Sqval = 0}
    
    HP_1 <- unlist(strsplit(input$Pset1, ","));HP_1 <- stri_trans_totitle(stri_trans_tolower(HP_1))
    HP_2 <- NULL
    if (length(input$Pset2) != 0) {HP_2 <- subset(gene_sets_dat,gs_name == input$Pset2) %>% pull(gene_symbol) %>% unique()}
    HP_c <- c(HP_1,HP_2)
    if (length(HP_c) == 0) {HP_c = NULL}
    
    proteasme_sets <- toupper(unique(HP_c))  ##### proteasme_sets is symbol or NULL
    
    
    
    if (input$group1 == "Sex") {
      sex_effect_limits <- c(lapply(tissue_specific_results, function(x) min(x$sex_effect)) %>%
                               unlist %>% min %>% floor,
                             lapply(tissue_specific_results, function(x) max(x$sex_effect)) %>%
                               unlist %>% max %>% ceiling)
      sex_logp_max <- lapply(tissue_specific_results, function(x) max(-log10(x$sex_pval))) %>%
        unlist %>% max %>% ceiling
      
      #  sex effects for tissues
      
      dda = tissue_specific_results_long %>% dplyr::mutate(sig = ifelse(sex_qval < Sqval, paste0('FDR<',Sqval), paste0('FDR>=', Sqval)))
      dda = dda[!is.na(dda$symbol),]
      dda$symbol = toupper(dda$symbol)
      dda = dda %>% dplyr::mutate(high = ifelse(symbol %in% proteasme_sets, 'Highlighted', 'Not-Highlighted'));#View(dda)
      dda$high = factor(dda$high, levels = c('Not-Highlighted','Highlighted'))
      dda$sig = factor(dda$sig, levels = unique(dda$sig))
      highlight_dat <- dda %>%  filter(symbol %in% proteasme_sets)

      A_plot = highlight_key(dda, ~ symbol) %>% ggplot(.,mapping = aes(x = sex_effect, y = -log10(sex_pval), group = symbol)) +  # group = symbol,
        geom_point(aes(fill = sig, shape = high, show.legend = FALSE)) +
       # geom_point(data = highlight_dat,
        #           col = 'red') +   color = sig,
        geom_vline(xintercept = 0, linetype = "dashed") +
        #scale_color_manual(values = c("orchid","gray")) +
        scale_fill_manual(values = c("gray","orchid")) +
        scale_shape_manual(values = c(1, 19)) +
        #geom_text(data = dda %>%  filter(high =='Highlighted' & -log10(sex_pval) > -log10(0.15)), mapping = aes(label = symbol), hjust = 0,colour = "black",  
          #nudge_x = 0.25, nudge_y = 0.25, 
         # size = 2.8, nudge_y = 0.5,
          #check_overlap = T ) +
        xlab("Estimated sex difference (abundance of females compared to males)") + ylab("-log10(p-value)") +
        ggtitle("Proteins with sex differences") +
        xlim(sex_effect_limits) + ylim(0, sex_logp_max) +
        guides(color = guide_legend(title = "",override.aes = aes(label = "")),fill = guide_legend(title = "",override.aes = aes(label = "")),shape = guide_legend(title = "",override.aes = aes(label = ""))) +
        facet_wrap(~ tissue,ncol = 5) +
        theme_bw()
      
      
      #Sex_plot =  ggplotly(Sex_plot,tooltip = "symbol")
      #SA_plot = highlight(Sex_plot,on = "plotly_hover", off = "plotly_doubleclick")
    }
    
    if (input$group1 == "Age") {
      
      ## Volcano plots
      age_effect_limits <- c(lapply(tissue_specific_results, function(x) min(x$age_effect)) %>%
                               unlist %>% min %>% floor,
                             lapply(tissue_specific_results, function(x) max(x$age_effect)) %>%
                               unlist %>% max %>% ceiling)
      age_logp_max <- lapply(tissue_specific_results, function(x) max(-log10(x$age_pval))) %>%
        unlist %>% max %>% ceiling
      
      
      # age effects for tissues
      dda = tissue_specific_results_long %>% dplyr::mutate(sig = ifelse(age_qval < Sqval, paste0('FDR<',Sqval), paste0('FDR>=', Sqval)))
      dda = dda[!is.na(dda$symbol),]	  
      dda$symbol = toupper(dda$symbol)
      dda = dda %>% dplyr::mutate(high = ifelse(symbol %in% proteasme_sets, 'Highlighted', 'Not-Highlighted'));#View(dda)
      dda$high = factor(dda$high, levels = c('Not-Highlighted','Highlighted'))
      dda$sig = factor(dda$sig, levels = unique(dda$sig))
  	  
      A_plot = highlight_key(dda, ~ symbol) %>% ggplot(.,mapping = aes(x = age_effect, y = -log10(age_pval), group = symbol)) +
        geom_point(aes(fill = sig, shape = high, show.legend = FALSE)) +
        geom_vline(xintercept = 0, linetype = "dashed") +
        #scale_color_manual(values = c("lightgreen","gray")) +
        # color = sig,
        scale_fill_manual(values = c("gray","lightgreen")) +
        scale_shape_manual(values = c(1, 19)) +
        #geom_text(data = dda %>%  filter(high =='Highlighted' & -log10(age_pval) > -log10(0.15)), mapping = aes(label = symbol), hjust = 0,colour = "black",  
        #nudge_x = 0.25, nudge_y = 0.25, 
        #size = 2.8, nudge_y = 0.5,
        #check_overlap = T ) +
        xlab("Estimated age difference (abundance of 18 month-olds compared to 8 month-olds)") + ylab("-log10(p-value)") +
        ggtitle("Proteins with age differences") +
        xlim(age_effect_limits) + ylim(0, age_logp_max) +
        guides(color = guide_legend(title = "",override.aes = aes(label = "")),fill = guide_legend(title = "",override.aes = aes(label = "")),shape = guide_legend(title = "",override.aes = aes(label = ""))) +
        facet_wrap(~ tissue,ncol = 5) +
        theme_bw()
      
      
      #Age_plot = ggplotly(Age_plot,tooltip = "symbol", source = "select")
      #SA_plot = highlight(Age_plot,on = "plotly_hover", off = "plotly_doubleclick")
      
    }
    
    
    #TSA_plot <- ggplotly(SA_plot) %>% event_register("plotly_click")
    A_plot = ggplotly(A_plot,tooltip = "symbol", source = "select") %>% style(textposition = "top") %>% event_register(c("plotly_click"))
    TSA_plot = highlight(A_plot,on = "plotly_hover",opacityDim = getOption("opacityDim", 1), color = 'blue', selected = attrs_selected(showlegend = FALSE)) %>% toWebGL()
    output$PvolPlot <- renderPlotly({TSA_plot})
    #dh1$SA_plot  ,"plotly_hover","plotly_doubleclick"
    #output$PvolPlot <- renderPlot({dh1$SA_plot}) 
    output$MSPBoxPlot <- renderPlot({})
    removeModal() 
  })
  # plotly_doubleclick
  ## partial_bundle()    off = "plotly_doubleclick",
  observeEvent(event_data(event = "plotly_click", source = "select"), {    
    geneD <- event_data(event = "plotly_click", source = "select")
    #print(geneD)
    tmp = protein_all_tissues_plotting_dat
    tmp$symbol = toupper(tmp$symbol)
    PSD = subset(tmp,symbol %in% geneD$key)
    dd <- Sex_Age_Effect_1(PSD)
    number_of_plot <- length(unique(PSD$tissue))
    TW = "auto"
    if ( number_of_plot < 6) {TW = 300*number_of_plot}
    output$MSPBoxPlot <- renderPlot({dd$S_A_plot},width = TW,res = 105)
    output$selected_title <- renderText({"Sex-Age effects for selected protein across tissues"})
    
  })
  
  ### plot heatmap 
  observeEvent(input$cluster, {
    
    showModal(modalDialog("Working On It! Just one minute", footer = NULL))
    
    #if (input$group1 == "Age") {dh1 <- Age_Effect(tissue_specific_results,tissue_specific_results_long)}
    #if (input$group1 == "Sex") {dh1 <- Sex_Effect_C(tissue_specific_results,tissue_specific_results_long)}		 
    #output$PvolPlot1 <- renderPlot({dh1$SA_plot}) 
    
    Cqval = input$Hqval
    if (input$Hqval > 1) {Cqval = 1}
    if (input$Hqval < 0) {Cqval = 0}
    
    ## Collapse list if effects data
    raw_tissue_specific_results_dat <- do.call("rbind", tissue_specific_results) %>%
      group_by(tissue) %>%
      mutate(age_qval = p.adjust(age_pval, method = "BH")) %>%
      mutate(sex_qval = p.adjust(sex_pval, method = "BH")) %>%
      mutate(sex_by_age_qval = p.adjust(sex_by_age_pval, method = "BH")) %>%
      ungroup %>%
      mutate(age_sig = case_when(age_qval < Cqval & age_effect < 0 ~ "down",
                                 age_qval < Cqval & age_effect > 0 ~ "up",
                                 age_qval > 0 ~ "ns")) %>%
      mutate(age_sig = factor(age_sig, levels = c("up", "down", "ns"))) %>%
      mutate(sex_sig = case_when(sex_qval < Cqval & sex_effect < 0 ~ "down",
                                 sex_qval < Cqval & sex_effect > 0 ~ "up",
                                 sex_qval > 0 ~ "ns")) %>%
      mutate(sex_sig = factor(sex_sig, levels = c("up", "down", "ns")))

    tissue_specific_results_dat <- raw_tissue_specific_results_dat %>%
      mutate(tissue = ifelse(tissue == "skeletalmuscle", "skeletal muscle", tissue)) %>%
      mutate(tissue = simple_cap(tissue)) %>%
      mutate(tissue = factor(tissue, levels = simple_cap(c("kidney", "liver", "fat", "spleen", "lung", 
                                                           "heart", "skeletal muscle", "striatum", "cerebellum", "hippocampus"))))
    #tissue_specific_results_dat = tissue_specific_results_dat[!is.na(tissue_specific_results_dat$symbol),]
    #tissue_specific_results_dat$symbol = toupper(tissue_specific_results_dat$symbol)  ### revise 1: symbol  
	if (input$group2 == "Sex") {
      # Sex effect matrix     revise 2: protein.id to symbol
      full_sex_effect_mat <- tissue_specific_results_dat %>%
        dplyr::select(tissue, protein.id, sex_effect) %>%
        pivot_wider(names_from = "tissue", values_from = "sex_effect") %>%
        column_to_rownames("protein.id") %>%
        as.matrix %>%
        t
      
      # Proteins with age effect   revise 3: protein.id to symbol
      sex_effect_proteins <- tissue_specific_results_dat %>%
        filter(sex_qval < Cqval) %>%
        pull(protein.id) %>%
        unique

      ## Process age effects
      reduced_sex_effect_mat <- full_sex_effect_mat[,sex_effect_proteins]
      # Set NAs to 0, necessary for the clustering to work
      reduced_sex_effect_mat[is.na(reduced_sex_effect_mat)] <- 0
      my_palette <- colorRampPalette(c("dodgerblue3", "white", "firebrick3")) (n = 20)
      
      sex_effects_heatmap <- heatmap.2(reduced_sex_effect_mat, trace = "none", na.color = "black", scale = "none", 
                                       col = my_palette, labCol = NA)
      tissue_order <- sex_effects_heatmap$rowInd
      protein_order <- sex_effects_heatmap$colInd
      
      ## Using ggplot
      ig_set1 <- colnames(reduced_sex_effect_mat[,protein_order])  
      
      Ctem = ggplot(data = tissue_specific_results_dat %>%
                      filter(protein.id %in% ig_set1) %>%
                      mutate(protein.id = factor(protein.id, levels = colnames(reduced_sex_effect_mat)[protein_order])) %>%
                      mutate(tissue = factor(tissue, levels = rownames(reduced_sex_effect_mat)[tissue_order])),
                    aes(y = tissue, x = protein.id, fill = sex_effect,text = paste0('Symbol: ',toupper(symbol)))) +
        geom_tile() +
        scale_fill_gradient2(low = "dodgerblue3", mid = "white", high = "firebrick3", name = "Sex effect\n(Female vs Male)") +
        xlab("Proteins") + ylab("") +
        theme(axis.text.y = element_text(size = 12),
              axis.text.x = element_blank(),
              axis.ticks.x = element_blank(),
              axis.ticks.y = element_blank())

    }
    
    if (input$group2 == "Age") {
      # Age effect matrix
      full_age_effect_mat <- tissue_specific_results_dat %>%
        dplyr::select(tissue, protein.id, age_effect) %>%
        pivot_wider(names_from = "tissue", values_from = "age_effect") %>%
        column_to_rownames("protein.id") %>%
        as.matrix %>%
        t
      
      # Proteins with age effect
      age_effect_proteins <- tissue_specific_results_dat %>%
        filter(age_qval < Cqval) %>%
        pull(protein.id) %>%
        unique
      
      ## Process age effects
      reduced_age_effect_mat <- full_age_effect_mat[,age_effect_proteins]
      # Set NAs to 0, necessary for the clustering to work
      reduced_age_effect_mat[is.na(reduced_age_effect_mat)] <- 0
      my_palette <- colorRampPalette(c("dodgerblue3", "white", "firebrick3")) (n = 20)
      
      age_effects_heatmap <- heatmap.2(reduced_age_effect_mat, trace = "none", na.color = "black", scale = "none", 
                                       col = my_palette, labCol = NA)
      tissue_order <- age_effects_heatmap$rowInd
      protein_order <- age_effects_heatmap$colInd
      
      ## Using ggplot

      ig_set1 <- colnames(reduced_age_effect_mat[,protein_order])  

      Ctem = ggplot(data = tissue_specific_results_dat %>%
               filter(protein.id %in% ig_set1) %>%
               mutate(protein.id = factor(protein.id, levels = colnames(reduced_age_effect_mat)[protein_order])) %>%
               mutate(tissue = factor(tissue, levels = rownames(reduced_age_effect_mat)[tissue_order])),
             aes(y = tissue, x = protein.id, fill = age_effect,text = paste0('Symbol: ',toupper(symbol)))) +
        geom_tile() +
        scale_fill_gradient2(low = "dodgerblue3", mid = "white", high = "firebrick3", name = "Age effect\n(18 vs 8 month)") +
        xlab("Proteins") + ylab("") +
        theme(axis.text.y = element_text(size = 12),
              axis.text.x = element_blank(),
              axis.ticks.x = element_blank(),
              axis.ticks.y = element_blank())

    }
    
    output$HeatmapPlot <- renderPlotly({ggplotly(Ctem)})
    removeModal() 
  })
  
  
  ### Age/Sex effect comparison
  ### select pathway danmicly 
  
  
  PSD_1 <- reactive({
  
  if (input$group3 == "Age") {NPD = fgsea_age_tables }				### type==1: Ensembl ID
  if (input$group3 == "Sex") {NPD = fgsea_sex_tables }	

  return(NPD)
  })
  
  
  PSD_2 <- reactive({
    
    pathway_dat = NULL
    if ( length(input$group4) == 2) {
    pathway_dat <- tissue_fgsea_compare_pathwaylevel(tissue1 = input$group4[1],
                                                     tissue2 = input$group4[2],
                                                     fgsea_list = PSD_1())
    }
    return(pathway_dat)
  })


  observeEvent(input$group4, {
    output$pathway_select <- renderUI({
      
      selectInput("path", "Pathway for comparison:", 
                  choices = PSD_2()$pathway,selected = NULL)
    })
  })
  
  
  observeEvent(input$Comp, {
    
    showModal(modalDialog("Working On It!", footer = NULL))
    
    dataC = tissue_fgsea_data(tissue1 = input$group4[1],
                              tissue2 = input$group4[2],
                              effects_list = tissue_specific_results,
							                fdr = 0.1,
							                gene_sets_dat = gene_sets_dat,
							                effect_type = tolower(input$group3),
                              fgsea_list = PSD_1(), 
                              highlight_pathway = input$path)
    
    genes_fgsea_dat = dataC
    tissue1 = input$group4[1];tissue2 = input$group4[2]
    effect_type = input$group3
    highlight_pathway = input$path
    fdr = 0.1; include_thresh = NULL; highlight_col = "dodgerblue3"
    
    include_thresh <- ifelse(is.null(include_thresh), fdr, include_thresh)
    
    p <- ggplot(data = genes_fgsea_dat %>%
                  arrange(highlight) %>%
                  filter(tissue1_padj < include_thresh | tissue2_padj < include_thresh),
                aes(y = tissue1_effect/tissue1_effect_se,
                    x = tissue2_effect/tissue2_effect_se,
                    size = fdr_status,
                    col = highlight,text = symbol)) +
      geom_point() +
      scale_color_manual(values = c("gray",highlight_col), drop = FALSE) +
      scale_size_manual(values = c(1, 2, 4), drop = FALSE) +
      geom_hline(yintercept = 0, linetype = "dashed") + geom_vline(xintercept = 0, linetype = "dashed") +
      ylab(paste("Standardized", effect_type, "effect in", tissue1)) + xlab(paste("Standardized", effect_type, "effect in", tissue2)) +
      ggtitle(paste("Comparison of", effect_type, "effects on proteins between", tissue1, "and", tissue2),
              subtitle = paste("GO category:", highlight_pathway)) +
      theme_bw() +
      theme(#legend.position = c(0.2, 0.85),
            legend.text = element_text(size = 8)) +
      guides(size = guide_legend(title = ""),
             col = guide_legend(title = "")) 
    
    C_plot = ggplotly(p,tooltip = "symbol") %>% 
            layout(title = list(text = paste0(paste("Comparison of", effect_type, "effects on proteins between", tissue1, "and", tissue2),
                                        '<br>',
                                        '<sup>',
                                        paste("GO category:", highlight_pathway),
                                        '</sup>'))) 

    output$ComPlot <- renderPlotly({C_plot})
    
    removeModal() 
  })
  

  
}

# Run the application 
shinyApp(ui = ui, server = server)
