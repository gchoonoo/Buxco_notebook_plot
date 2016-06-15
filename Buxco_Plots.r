source('utilities_MM.R')

# Heatmap

# Read in database
August2013_database.db <- makeBuxcoDB(db.name=file.path("August2013_database.db"))

# Output results for experiment rows only (EXP) and stratify by virus 
# Make sure that the names in outer.cols match the virus labels in the raw data (e.g. SARS vs. sars)

mvtsplot(August2013_database.db, outer.group.name='Virus', outer.cols=c(FLU="red", SARS="blue", Mock="springgreen4"), Break_type_label='EXP')

# Boxplot

# Observe variables
variables(August2013_database.db)

# Choose the variable and category to visualize
exp.penh <- retrieveData(August2013_database.db, variables="Penh", Break_type_label = 'EXP')

# Get table of categories
with(exp.penh, table(Days, Break_type_label))

# Create Boxplot
boxplot(Value~Sample_Name, data=exp.penh)

# Time Series
source("Buxco_Plot_Functions.r")

read.delim(file="./full_buxco_annotation_auc_v2.txt",header=T,sep="\t",colClasses=c('Sex'='character')) -> buxco_annot

# Set data1 to full annotation data, select Batch Date, Mating, Variable (Penh, EF50, Rpef, etc.), and Virus
buxco_plot(data1=buxco_annot, batch_date="November_2013", mating="13067x16912", var_data="Log(Penh)", virus=c("FLU","SARS","Mock"))

# Dot plot
buxco_dot_plot_data = dot_plot_data(var="Log(Penh)", virus=c("SARS"), lines=unique(buxco_annot$Mating), xlab=NULL, day=4, day_summary=4, vert_line=0)

dot_plot(buxco_dot_plot_data)
