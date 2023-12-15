## Title: <!--project dashboard title--> Differential abundance Operon Viewer (ASF519 Extracellular Vesicle)

<!--put a short project description here--> 
This Dashboard Allows Zooming in and scrolling through the Bacterial genes, ordered in its relative genomic position of the ASF519 reference genome, to show the relationship between the log2fold change, significance level, etc and the operons units.

## Motivation & Background

<!--write a paragraph or two about your motivation for creating this dashboard 
and provide some background information on the dashboard content if needed--> 
The Genetic components, mechanism and the function of Bacterial Extracellular Vesicles (EVs) are understudied. My current research is focused on comparing the genomic elements of EVs and its donor bacteria. I hypothesize that the DNA components of EVs are strategically sorted by the Donor Bacteria then pacakaged into EV cargos. Since Bacterial Genomes consists of genes that are grouped into small neighboring units of operons, which serves similar function, showing the relationship between the differential abundance (EVs vs. Donor Bacteria) pattern will partly provice evidence that EV's DNA components are biased towards serving specific functions, and that donors may selectively (and not randomly) choose which DNA are passed down to EVs. Since there are around 5500 locus in ASF519 model bacteria and around 1200 operon units, it is an extremely difficult visualization task to show operon units, LFCs, P-values, and its pattern all together in a comprehensive and intuitive way. As far as my knowledge goes, this will be the very first attempt in the field of microbiology to attempt this task. By incorporating the interactive aspect of R Shiny, I believe that we could deliver our findings in a more flexible way compared to the limitations of a static image.

# Operon Data Visualization Application: Dashboard Description & Usage

<!-- 
In the first phase, list all the inputs / outputs that will be available in your
dashboard. Your app must have multiple tabs, so make sure to list them accordingly

example:

- <input/output seen at all times> 
- Tab 1
  - <input / output of tab 1> File Upload button [Sorted.Counts.txt files from mapping EV reads to Reference Genome] / Flag for correct filetype [Message saying whether the filetype is correct or not]
  - <input / output of tab 1> File Upload button [DNA copy counts annotation files] / Flag for correct filetype [Message saying whether the filetype is correct or not
  - <input / output of tab 1> File Upload button [Operon annotation file] / Flag for correct filetype [Message saying whether the filetype is correct or not]
  - <input / output of tab 1> Submit button [Submitting visualization job] / Activated Tab 2 
- Tab 2
  - <input / output of tab 2> checkbox for coloring LFC bars by COG functional categories / changed barplot colorscheme from black to colorscheme
  - <input / output of tab 2> checkbox for p-value viewing / appearing or disappearing p-value annotation from visualization
  - <input / output of tab 2> zoom +/- buttons / zooming in or out of visualization
  - <input / output of tab 2> horizontal scroll bars / visualization of differnt part of the relative genomic position.
...

We encourage you to make a sketch of your overall app layout (and for each of 
the tabs), which you are welcome to insert here as images or add to the pull
request so we can get a better idea of what you have in mind.
-->


## Overview
This Shiny application provides an interactive visualization of operon-unit masked log2fold normalized gene counts change, ordered by relative genomics position for microbial . It integrates various R packages to enable robust data processing and interactive plotting functionalities.

## Features
- **Data Upload**: Users can upload their own `log2getmm_deseq_full` and `operon_data` CSV files.
- **Interactive Visualization**: The application utilizes `ggplot2` and `plotly` for dynamic and interactive data visualization.
- **Data Filtering**: Users can select ranges of relative positions and customize the visualization according to their needs.
- **Color Customization**: The application allows for the customization of colors for different categories.

## Installation

Before running the application, ensure that you have R installed on your system. You can download R from [The Comprehensive R Archive Network (CRAN)](https://cran.r-project.org/).

Once R is installed, you need to install the required packages. Run the following commands in R to install them:

```R
install.packages("shiny")
install.packages("ggplot2")
install.packages("DT")
install.packages("plotly")
install.packages("dplyr")
install.packages("zoo")
install.packages("stringr")
install.packages("egg")
install.packages("tidyr")
install.packages("shinyWidgets")
install.packages("colourpicker")
```

## Running the Application

To run the application, clone this repository or download the source code to your local machine. Open R and set your working directory to the folder containing the application files. Then, run the application with the following command:

```R
shiny::runApp()
```

## Usage

### Data Upload
- **log2getmm_deseq_full CSV**: Upload your `log2getmm_deseq_full` file in CSV format.
- **operon_data CSV**: Upload your `operon_data` file in CSV format.

### Interactive Controls
- **Position Range Slider**: Select the range of relative positions for which you want to visualize the data.
- **Color Pickers**: Customize the colors for different categories as per your preference.
- **View Comparative Microbial Operon Plot**: Click this button to render the plot based on the selected parameters.

### Visualization
- The main panel displays the interactive plot. Hover over elements in the plot to see more details.

## Citation
To reference this tool in your research or publications, please cite the following:  
Byeongyeon Cho, Grace Moore, Loc-Duyen Pham et al., "DNA characterization reveals potential operon-unit packaging of extracellular vesicle cargo from a gut bacterial symbiont." Preprint published on December 4, 2023, Version 1. Available at Research Square: [DOI: 10.21203/rs.3.rs-3689023/v1](https://doi.org/10.21203/rs.3.rs-3689023/v1).

## License
This project is open source and available under the MIT License. For more details, please refer to the [MIT License documentation](https://opensource.org/licenses/MIT).

## Contact
For inquiries, contributions, or further information, please feel free to reach out to Byeongyeon Cho. Contact: [Byeongyeon Cho](mailto:byeongyeon_cho@hms.harvard.edu).


***
