# Welcome at GENOVA!
### GENome Organisation Visual Analytics
<img src="https://github.com/robinweide/HiSee/raw/master/LOGO.jpg" width="250">

# Data structures
As input for the loop- and TAD-specific analyses, we request a specific input. For loop-analyses, the format is a data.frame containing:

| Chromosome of upstream loop-anchor | Start of upstream loop-anchor | Stop of upstream loop-anchor| Chromosome of downstream loop-anchor | Start of downstream loop-anchor | Stop of downstream loop-anchor  |
| ------------- | ------------- |------------- |------------- |------------- |------------- |
| Chr1  | 100000  |  120000  |  Chr1  | 200000  |  220000  | 

For TAD-based analyses, we request a data.frame with:

| Chromosome of TAD | Start of TAD | Stop of TAD| 
| ------------- | ------------- |------------- |
| Chr1  | 100000  |  220000  | 
