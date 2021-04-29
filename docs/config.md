<!--
Add here global page variables to use throughout your website.
-->
+++
prepath = "heterogeneous-graphlet-counting"
author = "Joel Robertson"
mintoclevel = 2
menu1_name = "GSE68559"
menu1_desc =  "Homo sapiens smoking data"
menu1_cut = 25
menu1_normalisation_method = "upper_quartile"
menu1_threshold_method = "empirical_dist"
menu1_variance_cut = 0.01
menu1_coexpression_measure = "pidc"

menu2_desc = "Chickpea salinity tolerance data"
menu2_name = "Mayank_de_novo"
menu2_cut = 25
menu2_normalisation_method = "upper_quartile"
menu2_threshold_method = "empirical_dist"
menu2_variance_cut = 0.01
menu2_coexpression_measure = "pcit"
# Add here files or directories that should be ignored by Franklin, otherwise
# these files might be copied and, if markdown, processed by Franklin which
# you might not want. Indicate directories by ending the name with a `/`.
# Base files such as LICENSE.md and README.md are ignored by default.
ignore = ["node_modules/"]

# RSS (the website_{title, descr, url} must be defined to get RSS)
generate_rss = true
website_title = "Heterogeneous graphlet counting"
website_descr = "outlining results of network construction and analysis."
website_url   = "https://tlienart.github.io/FranklinTemplates.jl/"
+++

<!--
Add here global latex commands to use throughout your pages.
-->
\newcommand{\R}{\mathbb R}
\newcommand{\scal}[1]{\langle #1 \rangle}
