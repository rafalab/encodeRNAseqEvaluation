library(shiny)
widget_style <- "display: inline-block; vertical-align: text-top; padding: 0px; border: solid;
                  border-width: 0px; border-radius: 0px; border-color: #CCC;"
shinyUI(pageWithSidebar(
  headerPanel("ENCODE RNAseq Evaluation"),
  sidebarPanel(
    style="min-width:235px;max-width:320px", 
    submitButton("Update"),
    br(),
    HTML("<style>
          div.pos1{
               width:80px;
               }
          div.pos2{
               width:80px;
               }
          .container div {
              display: inline-block;
              width:180px;
              }
          </style>
          <div class='container'>
              <div class='pos1'>
                 <label>Protocol</label>
                 <select name='protocol' style='width:120px'>
                    <option value='PolyA_dUTP' selected='selected'>PolyA_dUTP</option>
                    <option value='PolyA_TruSeq'>PolyA_TruSeq</option>
                    <option value='PolyA_SMARTseq'>PolyA_SMARTseq</option>
                    <option value='Total_dUTP'>Total_dUTP</option>
                    <option value='Total_TruSeq'>Total_TruSeq</option>
                    <option value='Total_SMARTseq'>Total_SMARTseq</option>
                 </select>
              </div>
              <div class='pos2'>
                 <label>Gene Type</label>
                 <select name='genetype' style='width:120px;'>
                    <option value='protein_coding' selected='selected'>protein coding gene</option>
                    <option value='pseudogene'>pseudogene</option>
                    <option value='lincRNA'>lincRNA</option>
                    <option value='antisense'>antisense</option>
                 </select>
              </div>
          </div>
          <label>Cutoff by FPKM (Default is 1)</label> 
          <div class='container'>   
              <div class='pos1'>
                 <input type='number' name='cutFPKM' style='width:105px' value='1' min='0' max='1000'/>
              </div>
          </div>

          <br>
          <h5>Axes of SD plot:</h5>
          <div class='container'>
              <div class='pos1'>
                 <label>x-min</label>
                 <input type='number' name='xstart1' style='width:105px' value='-5' min='-30' max='10'/>
              </div>
              <div class='pos2'>
                 <label>x-max</label>
                 <input type='number' name='xend1' style='width:105px' value='10' min='-10' max='20'/>
              </div>
          </div>
          <div class='container'>
              <div class='pos1'>
                 <label>y-min</label>
                 <input type='number' name='ystart1' style='width:105px' value='0' min='0' max='10'/>
              </div>
              <div class='pos2'>
                 <label>y-max</label>
                 <input type='number' name='yend1' style='width:105px' value='1' min='0.01' max='20'/>
              </div>
          </div>

          <br>
          <h5>Axes of P0 plot:</h5>
          <div class='container'>
              <div class='pos1'>
                 <label>x-min</label>
                 <input type='number' name='xstart2' style='width:105px' value='-5' min='-30' max='-3'/>
              </div>
              <div class='pos2'>
                 <label>x-max</label>
                 <input type='number' name='xend2' style='width:105px' value='0' min='-5' max='20'/>
              </div>
          </div>
          <div class='container'>
              <div class='pos1'>
                 <label>y-min</label>
                 <input type='number' name='ystart2' style='width:105px' value='0' min='0' max='0.1'/>
              </div>
              <div class='pos2'>
                 <label>y-max</label>
                 <input type='number' name='yend2' style='width:105px' value='0.05' min='0.01' max='0.5'/>
              </div>
          </div>

          <br>
          <h5>CAT plot between RNA-seq replicates:</h5>
          <div class='container'>
              <div class='pos1'>
                 <label>x-min</label>
                 <input type='number' name='xstart3' style='width:105px' value='20' min='1' max='500'/>
              </div>
              <div class='pos2'>
                 <label>x-max</label>
                 <input type='number' name='xend3' style='width:105px' value='2000' min='50' max='10000'/>
              </div>
          </div>
          <div class='container'>
              <div class='pos1'>
                 <label>y-min</label>
                 <input type='number' name='ystart3' style='width:105px' value='0' min='0' max='1'/>
              </div>
              <div class='pos2'>
                 <label>y-max</label>
                 <input type='number' name='yend3' style='width:105px' value='1' min='0' max='1'/>
              </div>
          </div>
          <label>Constant to be added for the second plot</label>
          <label>(unit: FPKM; should be <= cutoff)</label> 
          <div class='container'>
              <div class='pos1'>
                 <input type='number' name='constant1' style='width:105px' value='1' min='0.00001' max='1000'/>
              </div>
          </div>


          <br>
          <h5>CAT plot between RNA-seq and microarray:</h5>
          <div class='container'>
              <div class='pos1'>
                 <label>x-min</label>
                 <input type='number' name='xstart4' style='width:105px' value='20' min='1' max='500'/>
              </div>
              <div class='pos2'>
                 <label>x-max</label>
                 <input type='number' name='xend4' style='width:105px' value='2000' min='50' max='10000'/>
              </div>
          </div>
          <div class='container'>
              <div class='pos1'>
                 <label>y-min</label>
                 <input type='number' name='ystart4' style='width:105px' value='0' min='0' max='1'/>
              </div>
              <div class='pos2'>
                 <label>y-max</label>
                 <input type='number' name='yend4' style='width:105px' value='1' min='0' max='1'/>
              </div>
          </div>
          <label>Constant to be added for the second plot</label>
          <label>(unit: FPKM; should be <= cutoff)</label> 
          <div class='container'>
              <div class='pos1'>
                 <input type='number' name='constant2' style='width:105px' value='1' min='0.00001' max='1000'/>
              </div>
          </div>
"
         )
    ),

      mainPanel(
        h3(textOutput("caption")),
        tabsetPanel(
            tabPanel("SD versus detrended log signal",plotOutput("sdplot",width="600px", height="600px")),
            tabPanel("Proportion of zeros",wellPanel(
              div(style = widget_style,plotOutput("p0plot",width="600px", height="600px")),
              div(style = widget_style,tableOutput("p0stbl"))
              )),
            tabPanel("CAT plot for replicates",wellPanel(
              #h5("First is plot based on all non-zero genes,with fold changes of genes at least have one zero in all replicates being set as 0."),
              #h5("Second is plot based on all genes, with constant added for zeros."),
                div(style = widget_style,plotOutput("catplot1",width="600px", height="600px")),
                div(style = widget_style,plotOutput("catplot2",width="600px", height="600px"))
                )),
            tabPanel("CAT plot comparing to microarrays", wellPanel(
              #h5("Only genes both in RNA-seq and microarray are displayed."),
              #h5("First is plot based on all non-zero genes,with fold changes of genes at least have one zero in all replicates being set as 0."),
              #h5("Second is plot based on all genes, with constant added for zeros."),
                         div(style = widget_style,plotOutput("catplotarray1",width="600px", height="600px")),
                         div(style = widget_style,plotOutput("catplotarray2",width="600px", height="600px"))
                         ))
            ))))
