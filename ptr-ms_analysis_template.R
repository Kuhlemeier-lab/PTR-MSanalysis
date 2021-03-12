# Analyse PTR-MS results
# Kuhlemeier group
# Marta Binaghi <marta.binaghi at ips.unibe.ch>
# Created January 6th 2021
# Last modified March 12th 2021
# License: GNU GPL-3

# This script is to analyse the PTR-MS data.
#
# Inputs:
# * raw data files
# * a sample file
# * mass settings of PTR-MS
#
# RAW DATA
# should be in a subfolder of the working directory.
# should be as exported from the PTR-MS control program to a text file
# should have the date written in the file name (better in this form:
# YYYYMMDD.txt), to associate the data file to the samples.
#
# SAMPLE FILE
# should have the following columns:
# "date"
# "sample"
# "cycle.start"
# "cycle.end"
# The table is comma separated and has a header as above.
#
# MASSES FILE
# This contains the mass of each measured compound as set in the machine
# (eg benzaldehyde true mass in 106, PTR-MS mass is 107), and their name.
# The PTR-MS mass gets rounded to the closest integer.
# The columns "compound" and "mass_ptrms" should be in the file.
# Note that is a mass is measured but not present in the masses file, it will
# be excluded from the processing!
#
#
# Outputs:
#
# The script creates an output directory in the working directory. In the
# directory we will have:
# * output_rawdata_onlysamples.csv
# * output_plantMean_selectedCompounds.csv
# * some plots


# libraries ---------------------------------------------------------------
library(plyr)
library(ggplot2)

# output folder -----------------------------------------------------------
outdir <- paste0("output_", format(Sys.time(), "%Y%m%d-%H%M"))

if ( dir.exists(paths = outdir) ) {
  stop("An output folder with same date and time already exist.
Rename ther existing folder or wait for one minute.")
} else {
  dir.create(path = outdir)
}

# masses file -------------------------------------------------------------

# The compund masses file is supposed to be in the working directory.
# The PTR-MS data files (as .txt) and the sample table file should be in a
# subfolder:

# Raw data subfolder:
rawdata_fd <- "rawdata"

if ( !dir.exists(paths = rawdata_fd) ) {
  stop("The specified folder is not found.")
}


## Read in the compound masses

# If you are using a different setting, you have to supply an appropriate
# file that contains the mass of each measured compound as set in the PTR-MS
# setting. The table has to have at least two columns:
# "compound" holding the name of each measured compound
# "mass_ptrms" holding the mass as set in the PTR-MS. Note that the mass_ptrms
# has to match the header of the measured data file (.txt file).

masses_file <- "ptrms_compound_mass.csv"

masses <- read.table(masses_file,
                     sep = ",",
                     header = TRUE,
                     stringsAsFactors = FALSE)
colnames(masses) <- tolower(colnames(masses))

if ( sum( c(grepl( pattern = "compound",
                 x = colnames(masses)),
            grepl( pattern = "mass_ptrms",
                   x = colnames(masses))
)) != 2 ) {
  stop("Your mass file is missing the 'compound' or 'mass_ptrms' column.")
} else {
  if ( is.numeric(masses$mass_ptrms) ) {
    print("Masses file seems good!")
  }
  else {
    stop("The content of the 'compound' column is not numeric.")
  }
}

# round PTR-MS masses to closest integer
masses$mass_ptrms <- round(masses$mass_ptrms, digits = 0)


# sample file -------------------------------------------------------------

## Read the sample table

# We import the table with the samples and the cycle start and end.
# Here we have a table that includes:
# "date"
# "sample"
# "cycle.start"
# "cycle.end"
# The table is comma separated and has a header as above.

sample_table <- "axw115_sample_cycles.csv"

sample_df <- read.table(paste0(rawdata_fd, "/", sample_table),
                        sep = ",",
                        header = TRUE,
                        stringsAsFactors = FALSE)
colnames(sample_df) <- tolower(colnames(sample_df))

# Test if all cycle starts are lower then cycle end
if ( sum(sample_df$cycle.end < sample_df$cycle.start) ) {
  broken_rows <- which(sample_df$cycle.end < sample_df$cycle.start)
  stop(paste0("In row(s) ",
              paste(broken_rows, collapse = ", "),
               " cycle.end number is smaller than cycle.start number. Please, fix the input file."))
}
# add a date column formatted as in the file names
sample_df$date_formatted <- gsub(pattern = "\\s*(\\d\\d)/(\\d\\d)/(\\d\\d\\d\\d)\\s*",
              replacement = "\\3\\2\\1",
              sample_df$date,
              perl = TRUE)
# sort by date
sample_df <- sample_df[order(sample_df$date_formatted), ]
# add a replicate id for each sample every time it is measured (allows to distinguish different
# flowers from same plant measured on the same day)
sample_df$replicate <- rep(0, times = nrow(sample_df))
for (sample in unique(sample_df$sample)) {
  sample_df$replicate[sample_df$sample == sample] <- 1:sum(sample_df$sample == sample)
}


# read data ---------------------------------------------------------------

## Read the PTR-MS data files

# In this case we read a text file as produced by manually exporting the
# data from the PTR-MS control program.

# To import all the files at once, I use the date as stored in the sample table,
# and look for text files with a matching name (filename contains the date, even
# if differently formatted).
# The date in the sample file is stored as dd/mm/yyyy
# The date in the filenames is stored as yyyymmdd

dates <- unique(sample_df$date_formatted)
print(paste0(length(dates),
       " dates have been found. Will look for corresponding data files."))

# list text files in the directory
alltxts <- list.files(path = rawdata_fd,
                      pattern = ".txt$")

# keep only the txt files that have date in the name matching our list of dates
rawdata_files <- character()
for ( dt in dates ) {
  matching_files <- alltxts[grep(pattern = dt,
                                 x = alltxts)]
  if ( length(matching_files) < 1 ) {
    stop(paste0("No file matching date ",
                dt,
                " was found. Please check."))
  } else if ( length(matching_files) > 1 ) {
    stop(paste0("More than one file was found matching date ",
                dt,
                ". Please check."))
  } else {
    rawdata_files <- c(rawdata_files, matching_files)
  }
  if (dt == dates[length(dates)]) {
    print("One data file was found for each date provided.")
  }
}
rm(dt, matching_files, alltxts)


# initialise a df to store the raw data from every file together
main_data <- data.frame()
# make a vector of masses as they should be in the header of the raw data
masses.tmp <- paste0("mz", masses$mass_ptrms)
# Read in data files
for ( rawfile in rawdata_files ) {
  rawdf.tmp <- read.table(paste0(rawdata_fd, "/",
                                 rawfile),
                          sep = "\t",
                          dec = ",",
                          na.strings = c("NA", "NaN"),
                          header = T,
                          stringsAsFactors = FALSE)
  colnames(rawdf.tmp) <- tolower(colnames(rawdf.tmp))
  # save cycle number in a new column, remove useless repeated cols of cycle number
  cycles.tmp <- rawdf.tmp[ , grep("cycle", colnames(rawdf.tmp))[1]]
  rawdf.tmp <- rawdf.tmp[ , !grepl("(cycle|x).*", colnames(rawdf.tmp))]
  # clean col names to have only mass rounded to closest integer
  colnames.rawdf.tmp <- round(as.numeric(gsub(pattern = "amplitude[.]{3}m[.]z[.](\\d+[.]\\d{2})[.]ch\\d+",
                                              replacement = "\\1",
                                              x = colnames(rawdf.tmp),
                                              perl = TRUE)),
                              digits = 0)
  colnames(rawdf.tmp) <- paste0("mz", colnames.rawdf.tmp)
  # check all the masses in the data are found in the masses reference df
  if ( sum( !colnames.rawdf.tmp %in% masses$mass_ptrms ) > 0 ) {
    unknownmasses.tmp <- colnames.rawdf.tmp[!colnames.rawdf.tmp %in% masses$mass_ptrms]
    warning(paste0("File ",
                   rawfile,"
                   Some masses that are present in the data measured are
                   not listed in the masses file. These masses measured will not be
                   included in the final dataframe: ",
                   paste0(unknownmasses.tmp, collapse = ", ") ))
    # remove columns with undefined masses
    rawdf.tmp <- rawdf.tmp[ , colnames(rawdf.tmp) %in% masses.tmp]
  }
  rm(unknownmasses.tmp, colnames.rawdf.tmp)
  # add a column with the file name and a column with the cycle
  rawdf.tmp <- cbind(filename = rep(rawfile, length.out = dim(rawdf.tmp)[1]),
                     cycle = cycles.tmp,
                     rawdf.tmp,
                     stringsAsFactors = FALSE)
  # append to df
  main_data <- rbind(main_data, rawdf.tmp)
  rm(rawdf.tmp, cycles.tmp)
}

## Add the sample information to the main_data df

# add formatted date to main_data to get correspondence with sample_df
main_data$date_formatted <- gsub(pattern = ".*(\\d{8}).*",
     replacement = "\\1",
     x = main_data$filename)
# add a column for the sample
main_data$sample <- as.character(NA)
# and one for the replicate
main_data$replicate <- as.character(NA)
# add sample info to main_data
for ( entry in 1:nrow(sample_df) ) {
  main_data$sample[main_data$date_formatted == sample_df$date_formatted[entry] &
                     main_data$cycle >= sample_df$cycle.start[entry] &
                     main_data$cycle <= sample_df$cycle.end[entry]
                   ] <- sample_df$sample[entry]
  main_data$replicate[main_data$date_formatted == sample_df$date_formatted[entry] &
                     main_data$cycle >= sample_df$cycle.start[entry] &
                     main_data$cycle <= sample_df$cycle.end[entry]
                   ] <- sample_df$replicate[entry]
}

# keep only rows that have a sample
clean_data <- main_data[!is.na(main_data$sample), ]

# write the clean table to a file
write.csv(x = clean_data, file = paste0(outdir, "/output_rawdata_onlysamples.csv"), quote = FALSE,
          row.names = FALSE)


# summarise data by flower ------------------------------------------------

# we usually measure a flower for 20 cycles. We then keep the last five cycles
# and make a mean of those. The mean represents the flower measurement.

# sort by date and cycle
clean_data <- clean_data[order(clean_data$date_formatted, clean_data$cycle), ]

# calculate the column mean of the last five rows in the dataframe (ie the last
# five cycles values of each  compound)
meanLastFive <- function(x) {
  return(apply(X = x[(nrow(x)-4):nrow(x), grepl("mz", colnames(x))],
               MARGIN = 2,
               FUN = mean, na.rm = TRUE))
}
# calculate similarly the SD
sdLastFive <- function(x) {
  return(apply(X = x[(nrow(x)-4):nrow(x), grepl("mz", colnames(x))],
               MARGIN = 2,
               FUN = sd, na.rm = TRUE))
}
# apply the meanLastFive to each sample in each day (group by date_formatted,
# sample and replicate)
clean_data_sampleMean <- ddply(clean_data[ , 3:33], .(date_formatted, sample, replicate), meanLastFive)
# and apply the sd (standard deviation) similarly. The SD dataframe is temporary
# and is used only to explore the variation in each measurement.
clean_data_sampleSD <- ddply(clean_data[ , 3:33], .(date_formatted, sample, replicate), sdLastFive)



# explore stability across cycles -----------------------------------------

# exploratory plot for oxygen levels (if > 0.3 there was a problem)
ggplot(data = clean_data_sampleMean,
       aes(x = 1:nrow(clean_data_sampleMean),
           y = mz32,
           col = date_formatted)) +
  geom_point() +
  geom_errorbar(aes(x = 1:nrow(clean_data_sampleMean),
                    ymin = mz32 - clean_data_sampleSD$mz32,
                    ymax = mz32 + clean_data_sampleSD$mz32)) +
  ylim(0, 0.5) +
  geom_hline(yintercept = 0.3, col = "red") +
  ggtitle("Sample mean and SD of oxygen content") +
  xlab("Sample") +
  ylab("Oxygen content +/- SD") +
  scale_color_discrete(name = "Date") +
  theme_classic()
ggsave(filename = paste0(outdir, "/plot_exploratory_oxygen_level.png"),
       width = 6, height = 4)

# Check if the last five cycles of each measurement are stable.
# To do this we look at the SD of each measurement, for one compound.
# Here we set the compound to be methylbenzoate (mz137)
compound <- "mz137"

ggplot(data = clean_data_sampleMean,
       aes(x = date_formatted,
           y = get(compound),
           col = sample)) +
  geom_point(position = position_dodge(0.7)) +
  geom_errorbar(aes(x = date_formatted,
                    ymin = get(compound) - clean_data_sampleSD[ , compound],
                    ymax = get(compound) + clean_data_sampleSD[ , compound]),
                position = position_dodge(0.7),
                width = 0.3) +
  ggtitle(paste0("Sample mean and SD of ", compound)) +
  xlab("Date") +
  ylab("Content mean +/- SD") +
  scale_color_discrete(name = "Sample") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(filename = paste0(outdir, "/plot_exploratory_", compound, "_level.png"),
       width = 6, height = 4)

# if you notice some too wide SD, you can get the raw data that produced it by
# looking at the clean_data df, and selecting the date and sample:
clean_data[clean_data$date_formatted == "20210104" &
             clean_data$sample == "1", ]
# we can also select the replicate:
clean_data[clean_data$date_formatted == "20210104" &
             clean_data$sample == "1" &
             clean_data$replicate == "4", ]

# in my case I want to exclude sample 1 rep 4 and sample 4 rep 4 because they
# show an irregular behaviour during the last 5 cycles analysed.
# To remove certain samples from the dataframe we can do:
clean_data_nooutlier <- clean_data[!(clean_data$sample == "1" &
                                      clean_data$replicate == "4"), ]
clean_data_nooutlier <- clean_data_nooutlier[!(clean_data_nooutlier$sample == "4" &
                                       clean_data_nooutlier$replicate == "4"), ]
# remember that if we remove some samples, we will have to recalculate the mean
# and SD (the mean and SD dataframes still include those outlier samples).
clean_data_sampleMean <- ddply(clean_data_nooutlier[ , 3:33], .(date_formatted, sample, replicate), meanLastFive)
# and apply the sd (standard deviation) similarly. The SD dataframe is temporary
# and is used only to explore the variation in each measurement.
clean_data_sampleSD <- ddply(clean_data_nooutlier[ , 3:33], .(date_formatted, sample, replicate), sdLastFive)

# If we make the diagnostic plot again, we should now see no big SD:
ggplot(data = clean_data_sampleMean,
       aes(x = date_formatted,
           y = get(compound),
           col = sample)) +
  geom_point(position = position_dodge(0.7)) +
  geom_errorbar(aes(x = date_formatted,
                    ymin = get(compound) - clean_data_sampleSD[ , compound],
                    ymax = get(compound) + clean_data_sampleSD[ , compound]),
                position = position_dodge(0.7),
                width = 0.3) +
  ggtitle(paste0("Sample mean and SD of ", compound)) +
  xlab("Date") +
  ylab("Content mean +/- SD") +
  scale_color_discrete(name = "Sample") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))



# plant mean across replicates --------------------------------------------

# We can now use the mean across the last five cycles in each replicate to
# calculate a mean value per plant (different flowers' mean).

meanPlant <- function(x) {
  return(apply(X = x[ , grepl("mz", colnames(x))],
               MARGIN = 2,
               FUN = mean, na.rm = TRUE))
}
sdPlant <- function(x) {
  return(apply(X = x[ , grepl("mz", colnames(x))],
               MARGIN = 2,
               FUN = sd, na.rm = TRUE))
}

plantMean <- ddply(clean_data_sampleMean[ , c(2, 4:31)], .(sample), meanPlant)
plantSD <- ddply(clean_data_sampleMean[ , c(2, 4:31)], .(sample), sdPlant)


# plot plant mean and SD --------------------------------------------------

# for one compound:
compound <- "mz137"
ggplot(data = plantMean,
       aes(x = sample,
           y = get(compound),
           col = sample)) +
  geom_point() +
  geom_errorbar(aes(x = sample,
                    ymin = get(compound) - plantSD[ , compound],
                    ymax = get(compound) + plantSD[ , compound])) +
  ggtitle(paste0("Plant mean and SD of ", compound)) +
  xlab("Sample") +
  ylab("Content mean +/- SD") +
  scale_color_discrete(name = "Sample") +
  theme_classic()

# we can add a label to distinguish conditions on the samples, eg genotypes
# create a correspondence sample - label
plant_gt <- c("1" = "W115",
              "2" = "W115",
              "3" = "W115",
              "4" = "W115",
              "5" = "W115",
              "6" = "P.axN",
              "7" = "P.axN",
              "8" = "P.axN",
              "9" = "P.axN",
              "10" = "P.axN")
# add a column with the label
plantMean$gt <- plant_gt[match(plantMean$sample, names(plant_gt))]
plantSD$gt <- plant_gt[match(plantSD$sample, names(plant_gt))]

# plot with colour based on the gt label:
ggplot(data = plantMean,
       aes(x = sample,
           y = get(compound),
           col = gt)) +
  geom_point() +
  geom_errorbar(aes(x = sample,
                    ymin = get(compound) - plantSD[ , compound],
                    ymax = get(compound) + plantSD[ , compound])) +
  ggtitle(paste0("Plant mean and SD of ", compound)) +
  xlab("Sample") +
  ylab("Content mean +/- SD") +
  scale_color_discrete(name = "Sample") +
  theme_classic()

# we can also plot the mean per label
ggplot(data = plantMean,
       aes(x = gt,
           y = get(compound),
           fill = gt)) +
  stat_summary_bin(fun = "mean",
                   geom = "bar") +
  geom_point(alpha = 0.6) +
  ggtitle(paste0("Genotype mean and SD of ", compound)) +
  xlab("Genotype") +
  ylab("Content mean +/- SD") +
  theme_classic() +
  theme(legend.position = "none")

# and finally we can plot a set of several compounds (listed in compounds)
compounds <- c("107", "137", "121", "123", "151", "153", "165", "213", "229",
               "109", "167")
# we have to reformat the data so we can plot all compounds in the same figure.
# We select only the compounds in compounds
plantMean_selected <- plantMean[ , colnames(plantMean) %in% c("sample", "gt",
                                                              paste0("mz", compounds))]
plantSD_selected <- plantSD[ , colnames(plantSD) %in% c("sample", "gt",
                                                            paste0("mz", compounds))]
# then we transform the tables into the long format to allow plotting
mean.tmp <- reshape(plantMean_selected,
        idvar = c("sample", "gt"),
        varying = list(names(plantMean_selected)[2:12]),
        v.names = "mean",
        direction = "long",
        times = names(plantMean_selected)[2:12],
        timevar = "mass")
# same with SD
SD.tmp <- reshape(plantSD_selected,
                    idvar = c("sample", "gt"),
                    varying = list(names(plantSD_selected)[2:12]),
                    v.names = "SD",
                    direction = "long",
                    times = names(plantSD_selected)[2:12],
                    timevar = "mass")
# we then merge them together
plant_df <- merge.data.frame(mean.tmp,
                             SD.tmp)
rm(mean.tmp, SD.tmp)

# and save them to a file
write.csv(x = plant_df, file = paste0(outdir, "/output_plantMean_selectedCompounds.csv"),
          quote = FALSE,
          row.names = FALSE)

# to have meaningful names of compounds instead of the mass, I make a named vector
# for the correspondence between the mass and the compound name
# the name of the compound is shortened to max 20 characters
label_compounds <- substr(masses$compound, 1, 20)
names(label_compounds) <- paste0("mz", masses$mass_ptrms)
# plot
ggplot(plant_df,
       aes(x = gt,
           y = mean,
           col = gt)) +
  geom_boxplot() +
  facet_wrap(~mass, scales = "free_y",
             labeller = as_labeller(label_compounds)) +
  scale_y_continuous(expand = expansion(mult=c(0,0.1)), limits = c(0, NA)) +
  ggtitle("Genotype mean and SD of several compounds") +
  xlab("Genotype") +
  ylab("Content mean +/- SD") +
  theme_classic() +
  theme(legend.position = "none")
ggsave(filename = paste0(outdir, "/plot_genotypeMean.png"),
       width = 8, height = 8)
