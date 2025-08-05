# plot card overlay module for histology overlay app


#required libs
if (!require("shiny", quietly = TRUE)) install.packages("shiny")
if (!require("Cardinal", quietly = TRUE)) {
  if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
  BiocManager::install("Cardinal")
}
if (!require("DT", quietly = TRUE)) install.packages("DT")
if (!require("grid", quietly = TRUE)) install.packages("grid")
if (!require("abind", quietly = TRUE)) install.packages("abind")
if (!require("tools", quietly = TRUE)) install.packages("tools")
if (!require("png", quietly = TRUE)) install.packages("png")
if (!require("jpeg", quietly = TRUE)) install.packages("jpeg")
if (!require("ggplot2", quietly = TRUE)) install.packages("ggplot2")
if (!require("scales", quietly = TRUE)) install.packages("scales")
if (!require("viridis", quietly = TRUE)) install.packages("viridis")

library(shiny)
library(Cardinal)
library(DT)
library(grid)
library(abind)
library(tools)
library(png)
library(jpeg)
library(ggplot2)
library(scales)
library(viridis)

#color palette function - handles various color schemes
cpal <- function(name) {
  tryCatch({
    switch(name,
      "Spectral" = rainbow(256),
      "Viridis" = {
        if (requireNamespace("viridis", quietly = TRUE)) {
          viridis::viridis(256)
        } else {
          rainbow(256)
        }
      },
      "Plasma" = {
        if (requireNamespace("viridis", quietly = TRUE)) {
          viridis::plasma(256)
        } else {
          heat.colors(256)
        }
      },
      "Inferno" = {
        if (requireNamespace("viridis", quietly = TRUE)) {
          viridis::inferno(256)
        } else {
          heat.colors(256)
        }
      },
      "Cividis" = {
        if (requireNamespace("viridis", quietly = TRUE)) {
          viridis::cividis(256)
        } else {
          rainbow(256)
        }
      },
      #default case - try hcl.colors, fallback to rainbow
      tryCatch(hcl.colors(256, name), error = function(e) rainbow(256))
    )
  }, error = function(e) {
    rainbow(256)  #ultimate fallback
  })
}

plot_card_UI <- function(id) {
    
    ns <- NS(id)
    
    col_choices<-c(hcl.pals())
    initial_cols<-c("Spectral", "Cividis", "Viridis", "Inferno", "Plasma", 
                    "Zissou 1", "Purple-Green", "Berlin", "PiYG", "Grays", 
                    "Batlow", "turku", "YlOrRd", "Terrain", 
                    "PrGn", "Green-Brown", "Hawaii", "Cork", "Rocket", "RdYlBu")
    
    col_choices<-c(initial_cols, setdiff(col_choices, initial_cols))
  
    tagList(  
      
        fluidRow(
          tags$head(
            #tags$style("label{font-family: BentonSans Book;}")
            #tags$style("label{font-size: 11px;} ")
          ),
          
          column(3, offset = 0,
                 radioButtons(ns("ion_viz3"), "Image visualization ions",
                                 c("First ion" = "viz_first", 
                                   "All ions (TIC)"="viz_all", 
                                    "Custom single / multiple"="custom")),
                 #numericInput(ns("mz_viz3"), "mz value for visualization",255.2),
                 uiOutput(ns("mz_viz3a")),
                 
                 fluidRow(
                   column(6, selectInput(ns("contrast3"), "Contrast enhancement", c( "none", "histogram",  "adaptive"))),
                   column(6, selectInput(ns("color3"), "Colorscale", col_choices, selected="Spectral"))
                 ),
                 fluidRow(
                   column(6, selectInput(ns("smooth3"), "Smoothing options", c("none", "mean", "gaussian", "bilateral", "adaptive", "diffusion", "guided"))),
                   column(6, checkboxInput(ns("normalize3"), "Scale multiple images?", value = TRUE))
                   ),
                 
                 fluidRow(
                   column(6, checkboxInput(ns("colorkey3"), "Draw color key?", value = TRUE)),
                   column(6, checkboxInput(ns("dark_bg"), "Dark background?", value = FALSE))
                   ),
                   
                 
                 
                 fluidRow(
                   column(6, numericInput(ns("width_im"), "Image plot width (px)", value = 800, step = 50)),
                   column(6, numericInput(ns("height_im"), "Image plot height (px)", value=600, step = 50))
                 ),
                 checkboxInput(ns("plot_pdata"), "Plot Phenotype data?", value=FALSE),
                 uiOutput(ns("plotpdata")),
                 checkboxInput(ns("expand_fonts"), "Extended font options?", value=FALSE),
                 uiOutput(ns("fonts")),                 
                 checkboxInput(ns("expand_runs"), "Select individual runs for plotting only?", value=FALSE),
                 uiOutput(ns("select_runs")),
                 p("___________________________________________________________________")
                 

                 
          ),
          column(9, uiOutput(ns("plot.window"))
          )  
        ),
        
        uiOutput(ns("spectrum"))
    )
          
  }
  

  plot_card_server <- function(id, overview_peaks_sel, spatialOnly=FALSE, allInputs=NULL) {

    moduleServer(id, function(input, output, session){
      
      #for dynamic UI
      ns = session$ns
      
      graphics.off()
      
      cat("plot_card_server started\n")
      cat("Data class:", class(overview_peaks_sel), "\n")
      
      #check if this is PNG MSI data
      if (is.list(overview_peaks_sel) && !is.null(overview_peaks_sel$type) && overview_peaks_sel$type == "png") {
        cat("PNG MSI data detected\n")
        cat("PNG dimensions:", dim(overview_peaks_sel$image), "\n")
        
        #handle PNG MSI data separately
        observe({
          output$plot.window <- renderUI({
            fluidRow(
              column(12, imageOutput(ns("plot3_pk")))
            )
          })
        })
        
        output$plot3_pk <- renderImage({
          #temp file to save the output
          outfile <- tempfile(fileext = '.png')
          
          #get transformation values from allInputs
          alpha_val <- if (!is.null(allInputs$alpha)) allInputs$alpha else 0.5
          scalex <- if (!is.null(allInputs$scalex)) allInputs$scalex else 1
          scaley <- if (!is.null(allInputs$scaley)) allInputs$scaley else 1
          rotate <- if (!is.null(allInputs$rotate)) allInputs$rotate else 0
          translate_x <- if (!is.null(allInputs$translate_x)) allInputs$translate_x/100 else 0  #scale to 0-1
          translate_y <- if (!is.null(allInputs$translate_y)) allInputs$translate_y/100 else 0  #scale to 0-1
          
          #layer visibility controls
          show_msi <- if (!is.null(allInputs$show_msi_layer)) allInputs$show_msi_layer else TRUE
          show_histology <- if (!is.null(allInputs$show_histology_layer)) allInputs$show_histology_layer else TRUE
          
          png(outfile, width = 800, height = 600)
          
          #set up plot area
          par(mar = c(0, 0, 0, 0))
          plot(c(0, 1), c(0, 1), type = "n", xlab = "", ylab = "", axes = FALSE)
          
          #display the MSI PNG image first (as backgroun) if enabled
          if (show_msi) {
            img <- overview_peaks_sel$image
            rasterImage(img, 0, 0, 1, 1)
          }
          
          #apply histology overlay if available and enabled
          if (show_histology && !is.null(allInputs$histology_upload)) {
            tryCatch({
              #read histology image
              file_ext <- tools::file_ext(allInputs$histology_upload$datapath)
              if (file_ext %in% c("png", "PNG")) {
                hist_img <- png::readPNG(allInputs$histology_upload$datapath)
              } else if (file_ext %in% c("jpg", "jpeg", "JPG", "JPEG")) {
                hist_img <- jpeg::readJPEG(allInputs$histology_upload$datapath)
              } else {
                return()
              }
              
              #apply transparency by modifyin the alpha channel
              if (length(dim(hist_img)) == 3) {
                #if no alpha channel, add one
                if (dim(hist_img)[3] == 3) {
                  hist_img_alpha <- array(alpha_val, dim = c(dim(hist_img)[1:2], 1))
                  hist_img <- abind::abind(hist_img, hist_img_alpha, along = 3)
                } else if (dim(hist_img)[3] == 4) {
                  #if alpha channel exists, modify it
                  hist_img[,,4] <- hist_img[,,4] * alpha_val
                }
              }
              
              #TODO: add rotation support here
              # apply transformations to histology overlay
              # calculate position with scaling and translation
              center_x <- 0.5 + translate_x
              center_y <- 0.5 + translate_y
              
              #calculate scaled dimensions
              width_scaled <- scalex
              height_scaled <- scaley
              
              #calculate corners for transformed image
              x_left <- center_x - width_scaled/2
              x_right <- center_x + width_scaled/2
              y_bottom <- center_y - height_scaled/2
              y_top <- center_y + height_scaled/2
              
              #apply overlay (without alpha parameter since rasterImage doesn't support it)
              rasterImage(hist_img, x_left, y_bottom, x_right, y_top)
              
            }, error = function(e) {
              cat("Error overlaying histology:", e$message, "\n")
            })
          }
          
          dev.off()
          
          list(src = outfile, alt = "PNG MSI Image with Histology Overlay")
        }, deleteFile = TRUE)
        
        return()  #exit early for PNG data
      }
      
      cat("Data dimensions:", dim(overview_peaks_sel), "\n")

      #create new overview_peaks_sel object with mean values
      if(is.null(fData(overview_peaks_sel)$mean)) {
        cat("Computing feature summaries...\n")
        tryCatch({
          overview_peaks_sel<-Cardinal::summarizeFeatures(overview_peaks_sel, verbose=F)
          cat("Feature summaries computed successfully\n")
        }, error = function(e) {
          cat("Error in summarizeFeatures:", e$message, "\n")
          showNotification(paste("Error computing feature summaries:", e$message), type="error")
          return()
        })
        
        if(class(overview_peaks_sel)=="try-error") {
          showNotification("No data available, please check your parameters or dataset", type="error")
          return()
        }
      }
      # browser()
      # #create new overview_peaks_sel object with media values
      # if(is.null(fData(overview_peaks_sel)$non_zero)) {
      #   overview_peaks_sel$non_zero<-Cardinal::summarizeFeatures(overview_peaks_sel, "nnzero")
      #   if(class(overview_peaks_sel)=="try-error") {
      #     showNotification("No data available, please check your parameters or dataset", type="error")
      #     return()
      #   }
      # }
      # 
 
      observe({
        
        output$plot.window <- renderUI({
          
          if(input$ion_viz3=="custom") {
            fluidRow(
              column(12, imageOutput(ns("plot3_pk"))),
              fluidRow(
                column(12, br()),  # Adding space between the image and the table
                column(12, DT::dataTableOutput(ns('tbl')))
              )
            )
          } else {
            fluidRow(
              column(12, imageOutput(ns("plot3_pk")))
            )
          }
          
        })
        
      })
      
      

      
         output$mz_viz3a <- renderUI ({
           req(overview_peaks_sel)
           
           if(input$ion_viz3=="custom") {
             
            
             mz_list<- mz(overview_peaks_sel)
             tbl<-as.data.frame(fData(overview_peaks_sel))
             if(!is.null(tbl$freq)) {
               tbl$freq<-round(tbl$freq, 2)
             }
             if(!is.null(tbl$mean)) {
               tbl$mean<-round(tbl$mean, 1)
             }
             tbl$mz<-round(tbl$mz, 4)
             
             
             updateNumericInput(session, ("width_im"), value=800, step=50)
             updateNumericInput(session, ("height_im"), value=450, step=50)
             
        
             
             #add something here at some point to only select mz, freq, count, and mean columns
             
             
             output$tbl <-DT::renderDataTable({
               tbl
               
             }, selection = "multiple")
  
             list(
               checkboxInput(ns("superpose"), "Superpose images?", value = FALSE),
               selectInput(ns("display_mode"), "Ion math?", 
                          c("none", "sum", "ratio", "subtract", "min", "max", 
                            "mean", "sd", "var", "multiply")),
               numericInput(ns("plusminus_viz3"), "+/- m/z for visualization", 0.05)
             )
           }
         
          })
         
        
    
      
         mz_viz3 <- reactive({
           req(overview_peaks_sel)
           if (is.null(input$tbl_rows_selected) || length(input$tbl_rows_selected) == 0) {
             return(NULL)
           }
           tbl <- as.data.frame(fData(overview_peaks_sel))
           return(tbl[input$tbl_rows_selected, "mz"])
         })

      

      observe({
        
        output$select_runs <- renderUI ({
          req(overview_peaks_sel)
          
          if(input$expand_runs) {
            
            run_list<- unique(run((overview_peaks_sel)))
            
            list(
              selectInput(ns("select_runs"),
                          label= "run selection (plot only)",
                          multiple=TRUE,
                          choices = run_list,
                          selected = run_list)
            )
          }
        })
        
        
      })
      
      observe({
        output$fonts <- renderUI ({
          
          if(input$expand_fonts){
            list(
              numericInput(ns("axis_size"), label = "Axis font scaling (%)", min = 0, value=100),
              numericInput(ns("title_size"), label = "Title font scaling (%)", min = 0, value=100),
              numericInput(ns("label_size"), label = "Label font scaling (%)", min = 0, value=100),
              numericInput(ns("subtitle_size"), label = "Subtitle font scaling (%)", min = 0, value=100)
            )
          }
          
          
        })
      })
      
      observe({
        output$plotpdata <- renderUI ({
          
          if(input$plot_pdata){
            list(
              selectInput(ns("pdata_var_plot"), label="pData variable to plot", 
                          choices = colnames(as.data.frame(pData(overview_peaks_sel)))[-c(1:2)]
                          )
            )
          }
          
          
        })
      })

      #add spectrum or not
      output$spectrum<-renderUI({
        
        if(input$ion_viz3=="custom") {
          DT::dataTableOutput(ns('tbl'))
        }

        
        if(spatialOnly==FALSE){
          tagList(
            fluidRow(
              column(6,uiOutput(ns("plot_ranges"))),
              column(6,uiOutput(ns("plot_ranges2"))),
            ),
           
            fluidRow(
              
              column(2,  uiOutput(ns("x_target"))),
              column(2,  numericInput(ns("x_tol"), "+/- m/z window", value=5)),
              column(2,  uiOutput(ns("slider_y_max"))),
              column(3,numericInput(ns("width_ms"), "Plot width (px)", value = 800)),
              column(3,numericInput(ns("height_ms"), "Plot height (px)", value=300))
            ),
            fluidRow(
              column(10, imageOutput(ns("plot4")), style = "margin-bottom: 0px; padding-bottom: 0px;")
              ),
            fluidRow(
              column(3, checkboxInput(ns("calc_ppm"), label = "Show ppm error?", width = '100%')),
              column(3, checkboxInput(ns("show_int"), label = "Show intensity?", width = '100%')),
              column(3, numericInput(ns("show_mz"), label = "# mz values to show?", width = '100%', value=0))
            ),
            checkboxInput(ns("spectrum_expand_fonts"), "Extended font options?", value=FALSE),
            uiOutput(ns("spectrum_fonts")),
          )
        } 
      })
      
      observe({
        output$spectrum_fonts <- renderUI ({
          
          if(input$spectrum_expand_fonts) {
            fluidRow(
              column(
                3,
                numericInput(
                  ns("spectrum_axis_size"),
                  label = "Axis font scaling (%)",
                  min = 0,
                  value = 100
                )
              ),
              # column(
              #   2,
              #   numericInput(
              #     ns("spectrum_title_size"),
              #     label = "Title font scaling (%)",
              #     min = 0,
              #     value = 100
              #   )
              # ),
              column(
                3,
                numericInput(
                  ns("spectrum_label_size"),
                  label = "Label font scaling (%)",
                  min = 0,
                  value = 100
                )
              ),
              # column(
              #   2,
              #   numericInput(
              #     ns("spectrum_subtitle_size"),
              #     label = "Subtitle font scaling (%)",
              #     min = 0,
              #     value = 100
              #   )
              # ),
              column(
                3,
                numericInput(
                  ns("linewidth"),
                  label="Linewidth",
                  min=0,
                  value=100
                )
              )
            )
            
          }
          
          
        })
      })
          
       observe({
        output$plot3_pk <- renderImage( {  #plot image in overview after peakpicking / reading file
          
          #req(overview_peaks_sel)
          
          #add dependency on alpha slider to trigger re-rendering
          alpha_val <- allInputs$alpha
          
          # temp file to save the output.
          # This file will be removed later by renderImage
          outfile <- tempfile(fileext = '.png')
          
          png(outfile, width = input$width_im, height = input$height_im)
          
          #ion <- switch(input$mode,
          #              "p"=786,
          #             "n"=255.2)
          
          
          if(!is.null(input$select_runs)) {
            overview_peaks_sel <- subsetPixels(overview_peaks_sel, run %in% input$select_runs)
          }
          
          # Check for clean rendering mode
          if (!is.null(allInputs$render_mode) && allInputs$render_mode == "clean") {
            # Extract MSI ion image using calibrated coordinates
            tryCatch({
              # Get the selected m/z values for visualization
              if (input$ion_viz3 == "custom" && !is.null(mz_viz3()) && length(mz_viz3()) > 0) {
                selected_mz <- mz_viz3()[1]  # Use first selected m/z
              } else if (input$ion_viz3 == "viz_first") {
                selected_mz <- mz(overview_peaks_sel)[1]  # Use first m/z
              } else {
                # For TIC or other modes, use first m/z as fallback
                selected_mz <- mz(overview_peaks_sel)[1]
              }
              
              #ccreate clean MSI image using cardinal
              test_image <- Cardinal::image(overview_peaks_sel, 
                                          mz = selected_mz,
                                          tolerance = if(!is.null(input$plusminus_viz3)) input$plusminus_viz3 else 0.05,
                                          enhance = if(input$contrast3 != "none") input$contrast3 else NULL,
                                          smooth = if(input$smooth3 != "none") input$smooth3 else NULL,
                                          col = cpal(input$color3),
                                          scale = input$normalize3)
              
              #extract pixel data from Cardinal plot
              test_pixels <- test_image$plots[[1]]$marks$pixels$encoding
              
              ##cnstruct clean pixel level data frame
              static_coords <- data.frame(
                x = test_pixels$x, 
                y = test_pixels$y, 
                color = test_pixels$color
              )
              
              # Normalize color range
              static_coords$color_normalized <- scales::rescale(static_coords$color)
              
              # Create clean ggplot without borders/axes
              p_clean <- ggplot(static_coords, aes(x = x, y = y, fill = color_normalized)) +
                geom_tile(alpha = alpha_val) +  # Apply transparency from alpha slider
                scale_fill_viridis_c(option = "inferno") +
                scale_y_reverse() +  # matches imaging convention
                theme_void() +       # removes all axes, ticks, and borders
                theme(plot.margin = unit(c(0,0,0,0), "cm"))  # remove all margins
              
              # Save the clean plot with transparent background
              ggsave(outfile, plot = p_clean, width = input$width_im/100, height = input$height_im/100, 
                     units = "in", dpi = 100, bg = "transparent")
              
              # If histology image is available, create overlay using traditional plotting
              if (!is.null(allInputs$histology_upload)) {
                tryCatch({
                  # Reopene the plot device to add histology overlay
                  png(outfile, width = input$width_im, height = input$height_im, bg = "white")
                  
                  # Load and display histology image first (as background)
                  histology_path <- allInputs$histology_upload$datapath
                  file_type <- tools::file_ext(histology_path)
                  histology_image <- switch(file_type,
                                            "png" = png::readPNG(histology_path),
                                            "jpg" = jpeg::readJPEG(histology_path),
                                            "jpeg" = jpeg::readJPEG(histology_path),
                                            stop("Unsupported file type"))
                  
                  # Set up plot area
                  par(mar = c(0, 0, 0, 0))
                  plot(range(static_coords$x), range(static_coords$y), type = "n", 
                       xlab = "", ylab = "", axes = FALSE)
                  
                  # Display histology image as background
                  
                  # Apply scaling to histology image based on allInputs values
                  scalex <- if (!is.null(allInputs$scalex)) allInputs$scalex else 1
                  scaley <- if (!is.null(allInputs$scaley)) allInputs$scaley else 1
                  rotate <- if (!is.null(allInputs$rotate)) allInputs$rotate else 0
                  translate_x <- if (!is.null(allInputs$translate_x)) allInputs$translate_x/100 else 0
                  translate_y <- if (!is.null(allInputs$translate_y)) allInputs$translate_y/100 else 0
                  
                  # Calculate center position
                  center_x <- (min(static_coords$x) + max(static_coords$x))/2 + translate_x
                  center_y <- (min(static_coords$y) + max(static_coords$y))/2 + translate_y
                  
                  # Calculate dimensions with scaling
                  width_full <- max(static_coords$x) - min(static_coords$x)
                  height_full <- max(static_coords$y) - min(static_coords$y)
                  width_scaled <- width_full * scalex
                  height_scaled <- height_full * scaley
                  
                  # Calculate corners for transformed image
                  x_left <- center_x - width_scaled/2
                  x_right <- center_x + width_scaled/2
                  y_bottom <- center_y - height_scaled/2
                  y_top <- center_y + height_scaled/2
                  
                  # Apply transformations with rotation if needed
                  if (!is.null(rotate) && rotate != 0) {
                    # Apply transformations with rotation angle
                    rasterImage(histology_image, 
                                x_left, y_bottom, x_right, y_top,
                                angle = rotate)
                  } else {
                    # No rotation needed
                    rasterImage(histology_image, 
                                x_left, y_bottom, x_right, y_top)
                  }
                  
                  # Overlay the MSI data with transparency
                  # Create a color matrix from the MSI data
                  x_coords <- sort(unique(static_coords$x))
                  y_coords <- sort(unique(static_coords$y))
                  
                  # Create image matrix
                  img_matrix <- array(0, dim = c(length(y_coords), length(x_coords), 4))
                  
                  # Fill in the MSI data with colors and alpha
                  for (i in seq_len(nrow(static_coords))) {
                    x_idx <- which(x_coords == static_coords$x[i])
                    y_idx <- which(y_coords == static_coords$y[i])
                    
                    # Get color from viridis palette
                    color_val <- static_coords$color_normalized[i]
                    viridis_color <- viridis::inferno(256)[round(color_val * 255) + 1]
                    rgb_vals <- col2rgb(viridis_color) / 255
                    
                    img_matrix[y_idx, x_idx, 1] <- rgb_vals[1]  # R
                    img_matrix[y_idx, x_idx, 2] <- rgb_vals[2]  # G
                    img_matrix[y_idx, x_idx, 3] <- rgb_vals[3]  # B
                    img_matrix[y_idx, x_idx, 4] <- alpha_val    # Alpha
                  }
                  
                  # Display the MSI overlay
                  # Note: MSI data doesn't need scaling as it's already the reference frame
                  rasterImage(img_matrix, 
                              min(static_coords$x), min(static_coords$y),
                              max(static_coords$x), max(static_coords$y))
                  
                }, error = function(e) {
                  cat("Error creating histology overlay in clean mode:", e$message, "\n")
                })
              }
              
              dev.off()
              
              # Return the clean image
              return(list(src = outfile, alt = "Clean MSI Image"))
              
            }, error = function(e) {
              cat("Error in clean rendering mode:", e$message, "\n")
              # Fall back to regular rendering
            })
          }
          
          vp_orig<-vizi_par()
          
          #set sizes
          if(input$expand_fonts) {
            req(input$axis_size)
            
            vizi_par(
              cex.axis=input$axis_size/100,
              cex.lab=input$label_size/100,
              cex.main=input$title_size/100,
              cex.sub=input$subtitle_size/100
            )
            #get margins
            # cur_mar<-par()$mar
            # 
            # new_mar<-c(cur_mar[1]+cex.labp/8, cur_mar[2]+cex.labp/2, cur_mar[3]+cex.mainp/2, cur_mar[4])
            # 
            # #if drawing colorkey, add a little on the left
            # if(input$colorkey3) {
            #   new_mar[4]=new_mar[4]+cex.axisp
            # }
            # 
            # cur_mgp<-par()$mgp
            # 
            # #new_mgp<-c(cur_mgp[1]+cex.labp/8, cur_mgp[2], 0)
            # new_mgp<-c(3+max(new_mar[1:2])/20, 1,0)
            # 
            
            
            
            
           } else {
              vizi_par(
                cex.axis=1,
                cex.lab=1,
                cex.main=1,
                cex.sub=1
              )

           }
           
          if(input$plot_pdata){
            req(input$pdata_var_plot)
            
            #create list of arguments for image
            arg_list<-list(overview_peaks_sel, 
                       input$pdata_var_plot,
                        key=(input$colorkey3),
                        col=pals::alphabet())
                        
            if(input$dark_bg) {
              arg_list$style <- "dark"
            }
            
            plt_tmp<-do.call(Cardinal::image, arg_list)
            
            
            # Cardinal::image(overview_peaks_sel,
            #                         input$pdata_var_plot,
            #                         key=(input$colorkey3),
            #                         #superpose=input$superpose,
            #                         col=pals::alphabet())
            print(plt_tmp,
                                  #cex.axis=req(cex.axisp),
                                  #cex.lab=cex.labp,
                                  #cex.main=cex.mainp,
                                  #cex.sub=cex.subp,
                                  #mar=new_mar,
                                  #mgp=new_mgp
                  )
            
            vizi_par(vp_orig)
          } else if (input$ion_viz3=="viz_all") {
            
            
            mz_range=range(mz(overview_peaks_sel))
            
            #find closest mz value to middle of range
            test_value <- mean(mz_range)
            
            
            #calculate the absolute differences
            differences <- abs(mz(overview_peaks_sel) - test_value)
            
            #find the index of the minimum difference
            closest_index <- which.min(differences)
            
            mz_set=mz(overview_peaks_sel)[closest_index]
            tol=max(differences) + differences[closest_index]+1
            
            plusminus=tol
            
            #TODO: this old way might be more memory efficient, check later
            # smoothing_option <- if (input$smooth3 != "none") paste0(", smooth ='", input$smooth3,"'") else ""
            # enhance_option <- if (input$contrast3 != "none") paste0(", enhance ='", input$contrast3,"'") else ""
            # 
            # image_command <- paste("Cardinal::image(overview_peaks_sel, 
            #                       mz=mz_set,
            #                       tolerance=round(plusminus,3), 
            #                       units='mz',
            #                       col=cpal(input$color3)",
            #                       enhance_option,
            #                       smoothing_option,",
            #                       scale=input$normalize3,
            #                       #superpose=input$superpose,
            #                       key=(input$colorkey3),
            #                       #cex.axis=req(cex.axisp),
            #                       #cex.lab=cex.labp,
            #                       #cex.main=cex.mainp,
            #                       #cex.sub=cex.subp,
            #                       #mar=new_mar,
            #                       #mgp=new_mgp
            # )")
            # 
            
            #print(eval(parse(text = image_command)))
            
            
            
            smoothing_option <- if (input$smooth3 != "none")  input$smooth3 else NULL
            enhance_option <- if (input$contrast3 != "none")  input$contrast3 else NULL
            
            
            arg_list<-list(overview_peaks_sel, 
                           mz=mz_set,
                           tolerance=round(plusminus,3), 
                           units='mz',
                           col=cpal(input$color3),
                            enhance=enhance_option,
                           smooth=smoothing_option,
                           scale=input$normalize3,
                           #superpose=input$superpose,
                           key=(input$colorkey3))
            
            if(input$dark_bg) {
              arg_list$style <- "dark"
            }
            
            print(do.call(Cardinal::image, arg_list))
            
            
            
            

            vizi_par(vp_orig)
            
            
          } else if (input$ion_viz3=="custom"){
            
            if(is.null(mz_viz3())){
              ion=mz(overview_peaks_sel[1,])
            } else {
              
              # observeEvent(input$mz_viz3,{
                ion=as.numeric(mz_viz3())
              # })
            }
            
            
            if(!is.null(input$display_mode) && input$display_mode!="none"){
              if(input$display_mode%in%c("min", "max", "min", "mean", "sum", "sd", "var")) { 
          
                
                
                select_vec<-as.character(mz(overview_peaks_sel)) %in% as.character(ion)
                #test to make sure there are 2 or more elements
                if(sum(select_vec)<2){
                  showNotification("At least two ions required for this calculation", type="error")
                  message("At least two ions required for this calculation")
                  return()
                }
                
                sm<-summarizePixels(overview_peaks_sel[select_vec,], stat=c(xic=input$display_mode), as="DataFrame")
                pData(overview_peaks_sel)$xic<-sm$xic
                
                label_txt=paste(input$display_mode, "mz(s)=", paste(ion, collapse=", "))
            
              }else if(input$display_mode=="ratio"){
                if(length(ion)!=2){
                  showNotification("Exactly two ions required for ratio (mz1/mz2)", type="error")
                  message("Exactly two ions required for ratio (mz1/mz2)")
                  return()
                } else {
                  mz1 <- spectra(subsetFeatures(overview_peaks_sel, mz=ion[1]))[1,]
                  mz2 <- spectra(subsetFeatures(overview_peaks_sel, mz=ion[2]))[1,]
                  overview_peaks_sel$xic <- (1 + mz1) / (1 + mz2)
                  
                  ion=round(ion, 4)
                  label_txt=paste("ratio of",ion[1],"/",ion[2])
                }
                
              } else if(input$display_mode=="subtract") {
                if(length(ion)!=2){
                  showNotification("Exactly two ions required for subtraction (mz1-mz2)", type="error")
                  message("Exactly two ions required for subtraction (mz1-mz2)")
                  return()
                } else {
                  mz1 <- spectra(subsetFeatures(overview_peaks_sel, mz=ion[1]))[1,]
                  mz2 <- spectra(subsetFeatures(overview_peaks_sel, mz=ion[2]))[1,]
                  overview_peaks_sel$xic <- (mz1) - (mz2)
                  ion=round(ion, 4)
                  label_txt=paste("difference of",ion[1],"-",ion[2])
                  
                }
                
                
              }else if(input$display_mode=="multiply"){
                nelements=length(ion)
                xic <- 1+spectra(subsetFeatures(overview_peaks_sel, mz=ion[1]))[1,]
                
                for(i in 2:nelements){
                  mz2= 1+spectra(subsetFeatures(overview_peaks_sel, mz=ion[i]))[1,]
                  overview_peaks_sel$xic=(xic) * (mz2)
                  ion=round(ion, 4)
                  label_txt=paste(ion[1],"*",ion[2])
                  }
                }
            
  
              
              
              if(sum(is.na(pData(overview_peaks_sel)$xic))==length(overview_peaks_sel)) {
                showNotification("This calculation does not work!")
                return()
              }
              
              
              
              # smoothing_option <- if (input$smooth3 != "none") paste0(", smooth ='", input$smooth3,"'") else ""
              # enhance_option <- if (input$contrast3 != "none") paste0(", enhance ='", input$contrast3,"'") else ""
              # 
              plusminus=input$plusminus_viz3
              
            #   image_command <- paste("Cardinal::image(overview_peaks_sel, 'xic',
            #                       tolerance=round(plusminus,3), 
            #                       units='mz',
            #                       col=cpal(input$color3)",
            #                          enhance_option,
            #                          smoothing_option,",
            #                       scale=input$normalize3,
            #                       #superpose=input$superpose,
            #                       key=(input$colorkey3),
            #                       #cex.axis=req(cex.axisp),
            #                       #cex.lab=cex.labp,
            #                       #cex.main=cex.mainp,
            #                       #cex.sub=cex.subp,
            #                       #mar=new_mar,
            #                       #mgp=new_mgp
            # )")
            #   
              smoothing_option <- if (input$smooth3 != "none")  input$smooth3 else NULL
              enhance_option <- if (input$contrast3 != "none")  input$contrast3 else NULL
              
              
              arg_list<-list(overview_peaks_sel,
                             'xic',
                             tolerance=round(plusminus,3), 
                             units='mz',
                             col=cpal(input$color3),
                             enhance=enhance_option,
                             smooth=smoothing_option,
                             scale=input$normalize3,
                             #superpose=input$superpose,
                             key=(input$colorkey3))
              
              if(input$dark_bg) {
                arg_list$style <- "dark"
              }
              
              print(matter::as_facets(do.call(Cardinal::image, arg_list), labels=label_txt))
              
              vizi_par(vp_orig)

            } else {
              
              
              
              mz_set=ion
              
              mz_range=range(mz(overview_peaks_sel))
              
              
              # smoothing_option <- if (input$smooth3 != "none") paste0(", smoothing ='", input$smooth3,"'") else ""
              # enhance_option <- if (input$contrast3 != "none") paste0(", enhance ='", input$contrast3,"'") else ""
              # 
              plusminus=input$plusminus_viz3
              
            #   image_command <- paste("Cardinal::image(overview_peaks_sel, 
            #                       mz=mz_set,
            #                       tolerance=round(plusminus,3), 
            #                       units='mz',
            #                       col=cpal(input$color3)",
            #                          enhance_option,
            #                          smoothing_option,",
            #                       scale=input$normalize3,
            #                       superpose=input$superpose,
            #                       key=(input$colorkey3),
            #                       #cex.axis=req(cex.axisp),
            #                       #cex.lab=cex.labp,
            #                       #cex.main=cex.mainp,
            #                       #cex.sub=cex.subp,
            #                       #mar=new_mar,
            #                       #mgp=new_mgp
            # )")
            #   
              
              
              smoothing_option <- if (input$smooth3 != "none")  input$smooth3 else NULL
              enhance_option <- if (input$contrast3 != "none")  input$contrast3 else NULL
              
              
              arg_list<-list(overview_peaks_sel,
                             #'xic',
                             mz=mz_set,
                             tolerance=round(plusminus,3), 
                             units='mz',
                             col=cpal(input$color3),
                             enhance=enhance_option,
                             smooth=smoothing_option,
                             scale=input$normalize3,
                             superpose=input$superpose,
                             key=(input$colorkey3))
              
              if(input$dark_bg) {
                arg_list$style <- "dark"
              }
              
              
              
              
              print(do.call(Cardinal::image, arg_list))
              
              vizi_par(vp_orig)

              
            }
          } else if (input$ion_viz3=="viz_first") {
            
            
            tol=0.05
            #browser()
            image_command <-Cardinal::image(overview_peaks_sel, 
                                  col=cpal(input$color3),
                                  #enhance=input$contrast3,
                                  #smooth=input$smooth3,
                                  scale=input$normalize3,
                                  #superpose=input$superpose,
                                  key=(input$colorkey3),
                                  #cex.axis=req(cex.axisp),
                                  #cex.lab=cex.labp,
                                  #cex.main=cex.mainp,
                                  #cex.sub=cex.subp,
                                  #mar=new_mar,
                                  #mgp=new_mgp
            )
            
            print(image_command)
            
            vizi_par(vp_orig)
            
          }
          
          #apply histology overlay if images are uploaded
          # This must be done BEFORE dev.off()
          if (!is.null(allInputs$histology_upload) && 
              !is.null(allInputs$msi_upload)) {
            
            tryCatch({
              ensure_3d <- function(image) {
                if (length(dim(image)) == 2) {
                  image <- abind::abind(image, image, image, along = 3)
                }
                else if (length(dim(image)) > 3) {
                  image <- image[,,,1]
                }
                return(image)
              }
              
              cat("Applying histology overlay with alpha =", allInputs$alpha, "\n")
              
              #load histology image
              histology_path <- allInputs$histology_upload$datapath
              file_type <- tools::file_ext(histology_path)
              histology_image <- switch(file_type,
                                        "png" = png::readPNG(histology_path),
                                        "jpg" = jpeg::readJPEG(histology_path),
                                        "jpeg" = jpeg::readJPEG(histology_path),
                                        stop("Unsupported file type"))
              
              histology_image <- ensure_3d(histology_image)
              
              #always add alpha channel based on slider value
              # remove existing alpha channel if present
              if (dim(histology_image)[3] == 4) {
                histology_image <- histology_image[,,,1:3]
              }
              
              #add new alpha channel based on slider
              alpha_channel <- array(allInputs$alpha, 
                                     dim = c(dim(histology_image)[1], 
                                             dim(histology_image)[2]))
              histology_image <- abind::abind(histology_image, alpha_channel, along = 3)
              
              #create and apply histology grob
              histology_grob <- rasterGrob(histology_image, interpolate = TRUE)
              
              histology_grob <- editGrob(
                histology_grob,
                vp = viewport(
                  x = unit(0.5, "npc") + unit(allInputs$translate_x, "mm"),
                  y = unit(0.5, "npc") + unit(allInputs$translate_y, "mm"),
                  angle = allInputs$rotate,
                  width = unit(1, "npc") * allInputs$scalex,
                  height = unit(1, "npc") * allInputs$scaley,
                  just = c("center", "center")
                )
              )
              
              #draw the overlay on the current plot
              grid.draw(histology_grob)
              
            }, error = function(e) {
              cat("Error in histology overlay:", e$message, "\n")
            })
          }
          
          dev.off()
          
          
          #return a list containing the filename
          list(src = outfile,
               contentType = 'image/png',
               width = input$width_im,
               height = input$height_im,
               alt = "This is alternate text")
        }, deleteFile = TRUE)
       }) 
       
       
       if(spatialOnly==FALSE) {
         observe({  
          output$plot_ranges<- renderUI( {
            
            req(overview_peaks_sel)
            
            a<-Cardinal::plot(overview_peaks_sel)
            
  
                  sliderInput(ns("mass_range_plot"), 
                              label = p("m/z range for MS plot (X)"), 
                              min = round(a$channels$x$limits[1]),
                              max = round(a$channels$x$limits[2]), 
                              value = round(a$channels$x$limits), 
                              step = NULL,
                              round = TRUE)
  
            
          })
         })
         
         observe({
          output$plot_ranges2<- renderUI( {
            req(overview_peaks_sel)
            #browser()
            #overview_peaks_sel<-Cardinal::summarizeFeatures(overview_peaks_sel)
            a<-Cardinal::plot(overview_peaks_sel, "mean")
            
            
           
            
            sliderInput (ns("int_range_plot"), 
                        label = p("intensity range for MS plot (Y)"), 
                        min = round(a$channels$y$limits[1],0),
                        max = round(a$channels$y$limits[2],0)*1.05, 
                        value = req(input$param_numeric), #round(a$par$ylim,0), 
                        step = a$channels$y$limits[2]/20,
                        round = TRUE
            )
            
                        
                        
            
          })
         })
         
         
         #set center of observed spectrum if using custom ion visualization
         observe({
           output$x_target <- renderUI({
             
             if(input$ion_viz3!="custom") {
               numericInput(ns("x_target"), "Center m/z value", value = NULL)
             } else {
               numericInput(ns("x_target"), "Center m/z value", value = mz_viz3())
             }
             
           })
         })
         
        
         
         
         
         observe( {
           
           output$slider_y_max <- renderUI({
             
             req(overview_peaks_sel)
            
             
             a<-Cardinal::plot(overview_peaks_sel, "mean")
             
             
             numericInput(ns("param_numeric"),
                          "Manual y-axis intensity max value",
                          min = round(a$channels$y$limits[1],0),
                          max = round(a$channels$y$limits[2],0),
                          value = round(a$channels$y$limits[2],0)
             )
           })
         })
               
          
         if(!is.null(overview_peaks_sel)) {
           
           updateSliderInput(session, ns("int_range_plot"), value = c(round(Cardinal::plot(req(overview_peaks_sel))$channels$y$limits,0) 
                                                                      ))
           updateNumericInput(session, ns("param_numeric"), value = round(Cardinal::plot(req(overview_peaks_sel))$channels$y$limits[2],0))
          }
        
        
        
         
         observe({
           output$plot4 <- renderImage( {
             
             
             req(overview_peaks_sel)
             req(input$mass_range_plot)
             
             # A temp file to save the output.
             # This file will be removed later by renderImage
             outfile <- tempfile(fileext = '.png')
             
             png(outfile, width = input$width_ms, height = input$height_ms)
             
            
             
             if(length(input$int_range_plot)==1) {
               ylim=c(0, input$int_range_plot)
             } else {
               ylim = input$int_range_plot
             }
             
             #change xlimits based on custom ion or not
             if(is.null(input$x_target) || is.na(input$x_target)){
               xlim=input$mass_range_plot
               
               overview_peaks_sel_plot<-overview_peaks_sel
               
             } else {
               xlim=c(input$x_target-input$x_tol, input$x_target+input$x_tol)
               
               #subsetFeatures to only include mz values within range
               overview_peaks_sel_plot<-subsetFeatures(overview_peaks_sel, mz > xlim[1] & mz < xlim[2])
               
             }
             
             if(!is.finite(xlim[1])){
               xlim=input$mass_range_plot
             }
             
             vp_orig<-vizi_par()
             if(input$spectrum_expand_fonts) {
               req(input$spectrum_axis_size)
               
               vizi_par(
                 cex.axis = req(input$spectrum_axis_size)/100,
                 cex.lab = input$spectrum_label_size/100,
                 cex.main = input$spectrum_label_size/100,
                 lwd = input$linewidth/100,
                 mar = c(0, 0, 1, 1)
               )
               
               
               #get margins
               #cur_mar<-par()$mar
               
               #new_mar<-c(cur_mar[1]+cex.labp/8, cur_mar[2]+cex.labp/2, cur_mar[3], cur_mar[4])
               
               #cur_mgp<-par()$mgp
               
               #new_mgp<-c(cur_mgp[1]+cex.labp/8, cur_mgp[2], 0)
               #new_mgp<-c(3+max(new_mar[1:2])/20, 1,0)
               
               #lwd=lwdp
                          
               
               
               
               
             } else {
               vizi_par(
                 cex.axis=1,
                 cex.lab=1,
                 cex.main=1,
                 cex.sub=1
               )
             }
             
             #browser()
             #overview_peaks_sel<-Cardinal::summarizeFeatures(overview_peaks_sel)
             
             p1<-Cardinal::plot(overview_peaks_sel_plot,
                                xlim=xlim,
                                ylim =ylim,
                                #cex.axis=req(cex.axisp),
                                #cex.lab=cex.labp,
                                #cex.main=cex.mainp,
                                #cex.sub=cex.subp,
                                #lwd=lwdp,
                                # mar=new_mar, 
                                # mgp=new_mgp, 
                                "mean",
                                annPeaks=input$show_mz,
                                free="y")
             print(p1)
             vizi_par(vp_orig)
             
             #check for ppm calc
             if(input$calc_ppm) {
               
               dat=overview_peaks_sel_plot
               x=mz(dat)
               targ_mz<-req(input$x_target)
               x_sel<-subset(x, x>=xlim[1] & x<= xlim[2])
               
               ppm_error<- round(1e6*(x_sel-targ_mz)/targ_mz, 2)
               
               p1_coord<-p1[[1]][[1]]$marks$peaks$encoding
               
               y_labs<-p1_coord$y[p1_coord$x %in% x_sel]
               
               if(length(ppm_error)==0){
                 showNotification("No ppm error calculated, are there any peaks?", type="warning")
                 return()
               } else {
                print(text(x=x_sel, y=y_labs+ylim[2]*.25, req(ppm_error)))
               }
               
             }
             
             
             
             if(input$show_int) {
              
               
               ###
               p1_coord<-p1[[1]][[1]]$marks$peaks$encoding
               
               
               
               dat=overview_peaks_sel_plot
               x=mz(dat)
               
               x_sel<-subset(x, x>=xlim[1] & x<= xlim[2])
               
               
               
               y_labs<-p1_coord$y[p1_coord$x %in% x_sel]
               
               if(length(y_labs)==0){
                 showNotification("No intensities found, are there any peaks?", type="warning")
                 return()
               } else {
                 print(text(x=x_sel, y=y_labs+ylim[2]*.15, req(round(y_labs, 0))))
               }
               
             }
             
             

             
             dev.off()
             
             # Return a list containing the filename
             list(src = outfile,
                  contentType = 'image/png',
                  width = input$width,
                  height = input$height,
                  alt = "This is alternate text")
           }, deleteFile = TRUE)
           
         })
       }
  })
        
  
}

